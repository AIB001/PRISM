#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Module Communication System

Provides standardized communication mechanisms between PRISM modules,
including data exchange, event notifications, and resource coordination.
"""

import json
import queue
import threading
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, List, Callable, Union
from dataclasses import dataclass, asdict
from enum import Enum
import pickle
import tempfile

from .logging_system import PrismLogger, LogLevel, EventType


class MessageType(Enum):
    """Types of inter-module messages"""
    DATA_REQUEST = "data_request"
    DATA_RESPONSE = "data_response" 
    STATUS_UPDATE = "status_update"
    ERROR_NOTIFICATION = "error_notification"
    COMPLETION_NOTIFICATION = "completion_notification"
    RESOURCE_REQUEST = "resource_request"
    RESOURCE_RESPONSE = "resource_response"
    CONFIGURATION_UPDATE = "configuration_update"


class MessagePriority(Enum):
    """Message priority levels"""
    LOW = 1
    NORMAL = 2
    HIGH = 3
    CRITICAL = 4


@dataclass
class ModuleMessage:
    """Standard message structure for inter-module communication"""
    sender: str
    receiver: str
    message_type: MessageType
    priority: MessagePriority
    timestamp: datetime
    message_id: str
    data: Dict[str, Any]
    requires_response: bool = False
    response_timeout: Optional[float] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'sender': self.sender,
            'receiver': self.receiver,
            'message_type': self.message_type.value,
            'priority': self.priority.value,
            'timestamp': self.timestamp.isoformat(),
            'message_id': self.message_id,
            'data': self.data,
            'requires_response': self.requires_response,
            'response_timeout': self.response_timeout
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ModuleMessage':
        """Create from dictionary"""
        return cls(
            sender=data['sender'],
            receiver=data['receiver'],
            message_type=MessageType(data['message_type']),
            priority=MessagePriority(data['priority']),
            timestamp=datetime.fromisoformat(data['timestamp']),
            message_id=data['message_id'],
            data=data['data'],
            requires_response=data.get('requires_response', False),
            response_timeout=data.get('response_timeout')
        )


class DataTransferProtocol:
    """Handles different types of data transfer between modules"""
    
    @staticmethod
    def transfer_small_data(data: Dict[str, Any]) -> Dict[str, Any]:
        """Transfer small data directly in memory"""
        return data.copy()
    
    @staticmethod
    def transfer_large_data(data: Dict[str, Any], temp_dir: Optional[Path] = None) -> Dict[str, Any]:
        """Transfer large data via temporary files"""
        if temp_dir is None:
            temp_dir = Path(tempfile.gettempdir()) / "prism_data_transfer"
        temp_dir.mkdir(exist_ok=True)
        
        # Create temporary file
        temp_file = temp_dir / f"data_transfer_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}.pkl"
        
        # Serialize data
        with open(temp_file, 'wb') as f:
            pickle.dump(data, f)
        
        return {
            'transfer_method': 'file',
            'file_path': str(temp_file),
            'data_size': temp_file.stat().st_size,
            'timestamp': datetime.now().isoformat()
        }
    
    @staticmethod
    def receive_large_data(transfer_info: Dict[str, Any]) -> Dict[str, Any]:
        """Receive large data from temporary file"""
        if transfer_info.get('transfer_method') != 'file':
            raise ValueError("Invalid transfer method for large data")
        
        file_path = Path(transfer_info['file_path'])
        if not file_path.exists():
            raise FileNotFoundError(f"Transfer file not found: {file_path}")
        
        # Deserialize data
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        
        # Cleanup temporary file
        try:
            file_path.unlink()
        except Exception:
            pass  # Ignore cleanup errors
        
        return data
    
    @staticmethod
    def estimate_data_size(data: Dict[str, Any]) -> int:
        """Estimate serialized size of data"""
        try:
            return len(pickle.dumps(data))
        except Exception:
            return len(str(data).encode('utf-8'))


class MessageBroker:
    """Central message broker for module communication"""
    
    def __init__(self, max_queue_size: int = 1000, temp_dir: Optional[Path] = None):
        self.message_queues: Dict[str, queue.PriorityQueue] = {}
        self.subscribers: Dict[MessageType, List[str]] = {}
        self.message_handlers: Dict[str, Dict[MessageType, Callable]] = {}
        self.active_modules: Dict[str, bool] = {}
        self.max_queue_size = max_queue_size
        self.temp_dir = temp_dir or Path(tempfile.gettempdir()) / "prism_communication"
        self.temp_dir.mkdir(exist_ok=True)
        
        # Statistics
        self.message_stats: Dict[str, int] = {
            'sent': 0,
            'delivered': 0,
            'failed': 0,
            'queued': 0
        }
        
        # Thread safety
        self._lock = threading.RLock()
        self._shutdown = threading.Event()
        
        # Background processing
        self._background_thread = threading.Thread(target=self._process_messages, daemon=True)
        self._background_thread.start()
        
        self.logger = PrismLogger("prism.message_broker")
        self.logger.info("Message broker initialized")
    
    def register_module(self, module_name: str) -> None:
        """Register a module with the broker"""
        with self._lock:
            if module_name not in self.message_queues:
                self.message_queues[module_name] = queue.PriorityQueue(self.max_queue_size)
                self.message_handlers[module_name] = {}
                self.active_modules[module_name] = True
                self.logger.info(f"Registered module: {module_name}")
    
    def unregister_module(self, module_name: str) -> None:
        """Unregister a module from the broker"""
        with self._lock:
            if module_name in self.active_modules:
                self.active_modules[module_name] = False
                self.logger.info(f"Unregistered module: {module_name}")
    
    def subscribe(self, module_name: str, message_type: MessageType, 
                 handler: Callable[[ModuleMessage], None]) -> None:
        """Subscribe a module to specific message types"""
        with self._lock:
            if message_type not in self.subscribers:
                self.subscribers[message_type] = []
            
            if module_name not in self.subscribers[message_type]:
                self.subscribers[message_type].append(module_name)
            
            self.message_handlers[module_name][message_type] = handler
            self.logger.info(f"Module {module_name} subscribed to {message_type.value}")
    
    def send_message(self, message: ModuleMessage) -> bool:
        """Send a message to target module"""
        try:
            # Handle large data transfers
            if 'data' in message.data and DataTransferProtocol.estimate_data_size(message.data) > 10240:  # 10KB threshold
                transfer_info = DataTransferProtocol.transfer_large_data(message.data, self.temp_dir)
                message.data = {'large_data_transfer': transfer_info}
            
            # Add to queue with priority
            priority_value = 5 - message.priority.value  # Lower number = higher priority in queue
            queue_item = (priority_value, time.time(), message)
            
            with self._lock:
                if message.receiver in self.message_queues:
                    try:
                        self.message_queues[message.receiver].put_nowait(queue_item)
                        self.message_stats['sent'] += 1
                        self.message_stats['queued'] += 1
                        return True
                    except queue.Full:
                        self.logger.warning(f"Message queue full for {message.receiver}")
                        self.message_stats['failed'] += 1
                        return False
                else:
                    self.logger.error(f"Module {message.receiver} not registered")
                    self.message_stats['failed'] += 1
                    return False
        
        except Exception as e:
            self.logger.error(f"Failed to send message: {e}")
            self.message_stats['failed'] += 1
            return False
    
    def get_message(self, module_name: str, timeout: Optional[float] = None) -> Optional[ModuleMessage]:
        """Get next message for a module"""
        if module_name not in self.message_queues:
            return None
        
        try:
            priority, timestamp, message = self.message_queues[module_name].get(timeout=timeout)
            self.message_stats['delivered'] += 1
            self.message_stats['queued'] -= 1
            
            # Handle large data transfers
            if 'large_data_transfer' in message.data:
                transfer_info = message.data['large_data_transfer']
                message.data = DataTransferProtocol.receive_large_data(transfer_info)
            
            return message
            
        except queue.Empty:
            return None
        except Exception as e:
            self.logger.error(f"Error getting message for {module_name}: {e}")
            return None
    
    def broadcast_message(self, message: ModuleMessage) -> int:
        """Broadcast message to all subscribers of the message type"""
        if message.message_type not in self.subscribers:
            return 0
        
        delivered_count = 0
        for subscriber in self.subscribers[message.message_type]:
            if self.active_modules.get(subscriber, False):
                # Create copy for each subscriber
                subscriber_message = ModuleMessage(
                    sender=message.sender,
                    receiver=subscriber,
                    message_type=message.message_type,
                    priority=message.priority,
                    timestamp=message.timestamp,
                    message_id=f"{message.message_id}_{subscriber}",
                    data=message.data.copy(),
                    requires_response=message.requires_response,
                    response_timeout=message.response_timeout
                )
                
                if self.send_message(subscriber_message):
                    delivered_count += 1
        
        return delivered_count
    
    def _process_messages(self) -> None:
        """Background thread for processing messages"""
        while not self._shutdown.is_set():
            try:
                # Process messages for all active modules
                for module_name in list(self.active_modules.keys()):
                    if not self.active_modules.get(module_name, False):
                        continue
                    
                    message = self.get_message(module_name, timeout=0.1)
                    if message and message.message_type in self.message_handlers.get(module_name, {}):
                        handler = self.message_handlers[module_name][message.message_type]
                        try:
                            handler(message)
                        except Exception as e:
                            self.logger.error(f"Message handler error in {module_name}: {e}")
                
                time.sleep(0.1)  # Small delay to prevent excessive CPU usage
                
            except Exception as e:
                self.logger.error(f"Error in message processing: {e}")
                time.sleep(1)
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get communication statistics"""
        with self._lock:
            return {
                'message_stats': self.message_stats.copy(),
                'active_modules': sum(self.active_modules.values()),
                'total_registered': len(self.message_queues),
                'queue_sizes': {name: q.qsize() for name, q in self.message_queues.items()},
                'subscriber_counts': {msg_type.value: len(subs) 
                                    for msg_type, subs in self.subscribers.items()}
            }
    
    def shutdown(self) -> None:
        """Shutdown the message broker"""
        self._shutdown.set()
        if self._background_thread.is_alive():
            self._background_thread.join(timeout=5)
        self.logger.info("Message broker shutdown")


class ModuleCommunicator:
    """High-level communication interface for modules"""
    
    def __init__(self, module_name: str, broker: MessageBroker):
        self.module_name = module_name
        self.broker = broker
        self.logger = PrismLogger(f"prism.{module_name}.communicator")
        self._message_id_counter = 0
        self._pending_responses: Dict[str, threading.Event] = {}
        self._response_data: Dict[str, Any] = {}
        
        # Register with broker
        self.broker.register_module(module_name)
        
        # Subscribe to response messages
        self.broker.subscribe(module_name, MessageType.DATA_RESPONSE, self._handle_response)
    
    def send_data(self, target_module: str, data: Dict[str, Any], 
                 priority: MessagePriority = MessagePriority.NORMAL,
                 requires_response: bool = False, timeout: Optional[float] = None) -> Optional[Dict[str, Any]]:
        """Send data to another module"""
        message_id = self._generate_message_id()
        
        message = ModuleMessage(
            sender=self.module_name,
            receiver=target_module,
            message_type=MessageType.DATA_REQUEST,
            priority=priority,
            timestamp=datetime.now(),
            message_id=message_id,
            data=data,
            requires_response=requires_response,
            response_timeout=timeout
        )
        
        if requires_response:
            # Setup response waiting
            response_event = threading.Event()
            self._pending_responses[message_id] = response_event
        
        # Send message
        success = self.broker.send_message(message)
        
        if success:
            self.logger.info(f"Sent data to {target_module}, message_id: {message_id}")
            
            if requires_response:
                # Wait for response
                if response_event.wait(timeout=timeout):
                    response_data = self._response_data.pop(message_id, None)
                    self._pending_responses.pop(message_id, None)
                    return response_data
                else:
                    self.logger.warning(f"Response timeout for message {message_id}")
                    self._pending_responses.pop(message_id, None)
                    return None
        else:
            self.logger.error(f"Failed to send data to {target_module}")
            if requires_response:
                self._pending_responses.pop(message_id, None)
        
        return None
    
    def broadcast_status(self, status_data: Dict[str, Any], 
                        priority: MessagePriority = MessagePriority.NORMAL) -> int:
        """Broadcast status update to all interested modules"""
        message_id = self._generate_message_id()
        
        message = ModuleMessage(
            sender=self.module_name,
            receiver="*",  # Broadcast indicator
            message_type=MessageType.STATUS_UPDATE,
            priority=priority,
            timestamp=datetime.now(),
            message_id=message_id,
            data=status_data
        )
        
        delivered_count = self.broker.broadcast_message(message)
        self.logger.info(f"Broadcasted status update to {delivered_count} modules")
        return delivered_count
    
    def notify_completion(self, results: Dict[str, Any]) -> int:
        """Notify completion to all interested modules"""
        message_id = self._generate_message_id()
        
        message = ModuleMessage(
            sender=self.module_name,
            receiver="*",
            message_type=MessageType.COMPLETION_NOTIFICATION,
            priority=MessagePriority.HIGH,
            timestamp=datetime.now(),
            message_id=message_id,
            data=results
        )
        
        delivered_count = self.broker.broadcast_message(message)
        self.logger.info(f"Notified completion to {delivered_count} modules")
        return delivered_count
    
    def notify_error(self, error_info: Dict[str, Any]) -> int:
        """Notify error to all interested modules"""
        message_id = self._generate_message_id()
        
        message = ModuleMessage(
            sender=self.module_name,
            receiver="*",
            message_type=MessageType.ERROR_NOTIFICATION,
            priority=MessagePriority.CRITICAL,
            timestamp=datetime.now(),
            message_id=message_id,
            data=error_info
        )
        
        delivered_count = self.broker.broadcast_message(message)
        self.logger.error(f"Notified error to {delivered_count} modules")
        return delivered_count
    
    def _handle_response(self, message: ModuleMessage) -> None:
        """Handle response messages"""
        response_message_id = message.data.get('response_to')
        if response_message_id in self._pending_responses:
            self._response_data[response_message_id] = message.data
            self._pending_responses[response_message_id].set()
    
    def _generate_message_id(self) -> str:
        """Generate unique message ID"""
        self._message_id_counter += 1
        return f"{self.module_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{self._message_id_counter}"
    
    def cleanup(self) -> None:
        """Cleanup communicator resources"""
        self.broker.unregister_module(self.module_name)
        self._pending_responses.clear()
        self._response_data.clear()


# Convenience functions for setting up communication
def create_message_broker(config: Optional[Dict[str, Any]] = None) -> MessageBroker:
    """Create and configure message broker"""
    config = config or {}
    
    broker = MessageBroker(
        max_queue_size=config.get('max_queue_size', 1000),
        temp_dir=Path(config.get('temp_dir', tempfile.gettempdir())) / "prism_communication"
    )
    
    return broker


def create_module_communicator(module_name: str, broker: MessageBroker) -> ModuleCommunicator:
    """Create module communicator"""
    return ModuleCommunicator(module_name, broker)