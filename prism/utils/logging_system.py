#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Advanced Logging System

Provides comprehensive logging, monitoring, and event tracking capabilities
for PRISM molecular dynamics simulations with structured logging, performance
tracking, and real-time monitoring.
"""

import os
import json
import time
import logging
import threading
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, List, Union, Callable
from contextlib import contextmanager
from dataclasses import dataclass, asdict
from enum import Enum
import subprocess
import psutil


class LogLevel(Enum):
    """Extended log levels for PRISM"""
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"
    PERFORMANCE = "PERFORMANCE"
    WORKFLOW = "WORKFLOW"
    SCIENTIFIC = "SCIENTIFIC"


class EventType(Enum):
    """Event types for structured logging"""
    SIMULATION_START = "simulation_start"
    SIMULATION_END = "simulation_end"
    STEP_START = "step_start"
    STEP_END = "step_end"
    ERROR_OCCURRED = "error_occurred"
    WARNING_ISSUED = "warning_issued"
    PERFORMANCE_METRIC = "performance_metric"
    RESOURCE_USAGE = "resource_usage"
    VALIDATION_RESULT = "validation_result"
    CONFIGURATION_CHANGE = "configuration_change"
    FILE_OPERATION = "file_operation"
    EXTERNAL_TOOL = "external_tool"


@dataclass
class LogEvent:
    """Structured log event"""
    timestamp: datetime
    event_type: EventType
    level: LogLevel
    message: str
    module: str
    function: str
    context: Dict[str, Any]
    duration: Optional[float] = None
    session_id: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        data = asdict(self)
        data['timestamp'] = self.timestamp.isoformat()
        data['event_type'] = self.event_type.value
        data['level'] = self.level.value
        return data


@dataclass
class PerformanceMetrics:
    """Performance metrics for monitoring"""
    cpu_percent: float
    memory_mb: float
    disk_io_read_mb: float
    disk_io_write_mb: float
    network_io_kb: Optional[float] = None
    gpu_percent: Optional[float] = None
    gpu_memory_mb: Optional[float] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class StructuredFormatter(logging.Formatter):
    """Custom formatter for structured logging"""
    
    def format(self, record):
        log_data = {
            'timestamp': datetime.fromtimestamp(record.created).isoformat(),
            'level': record.levelname,
            'logger': record.name,
            'message': record.getMessage(),
            'module': record.module,
            'function': record.funcName,
            'line': record.lineno
        }
        
        # Add extra fields if present
        if hasattr(record, 'context'):
            log_data['context'] = record.context
        if hasattr(record, 'event_type'):
            log_data['event_type'] = record.event_type
        if hasattr(record, 'session_id'):
            log_data['session_id'] = record.session_id
        if hasattr(record, 'duration'):
            log_data['duration'] = record.duration
        
        return json.dumps(log_data, default=str, ensure_ascii=False)


class PrismLogger:
    """Advanced logger for PRISM with monitoring capabilities"""
    
    def __init__(self, name: str, session_id: Optional[str] = None):
        self.name = name
        self.session_id = session_id or self._generate_session_id()
        self.logger = logging.getLogger(name)
        self.events: List[LogEvent] = []
        self.performance_history: List[PerformanceMetrics] = []
        self.start_time = time.time()
        self._setup_logger()
        
    def _generate_session_id(self) -> str:
        """Generate unique session ID"""
        return f"prism_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
    
    def _setup_logger(self):
        """Setup logger with appropriate handlers"""
        self.logger.setLevel(logging.DEBUG)
        
        # Remove existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)
        
        # Add console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
    
    def add_file_handler(self, log_file: Union[str, Path], 
                        structured: bool = True, level: str = "DEBUG"):
        """Add file handler for logging"""
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(getattr(logging, level))
        
        if structured:
            file_handler.setFormatter(StructuredFormatter())
        else:
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s'
            ))
        
        self.logger.addHandler(file_handler)
        self.info(f"Added file handler: {log_file}")
    
    def log_event(self, event_type: EventType, level: LogLevel, message: str,
                  context: Optional[Dict[str, Any]] = None, duration: Optional[float] = None):
        """Log structured event"""
        import inspect
        
        # Get caller information
        frame = inspect.currentframe().f_back
        module = frame.f_globals.get('__name__', 'unknown')
        function = frame.f_code.co_name
        
        event = LogEvent(
            timestamp=datetime.now(),
            event_type=event_type,
            level=level,
            message=message,
            module=module,
            function=function,
            context=context or {},
            duration=duration,
            session_id=self.session_id
        )
        
        self.events.append(event)
        
        # Log to standard logger
        log_level = getattr(logging, level.value)
        extra = {
            'context': context,
            'event_type': event_type.value,
            'session_id': self.session_id,
            'duration': duration
        }
        self.logger.log(log_level, message, extra=extra)
    
    def debug(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Debug level logging"""
        self.log_event(EventType.STEP_START, LogLevel.DEBUG, message, context)
    
    def info(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Info level logging"""
        self.log_event(EventType.STEP_START, LogLevel.INFO, message, context)
    
    def warning(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Warning level logging"""
        self.log_event(EventType.WARNING_ISSUED, LogLevel.WARNING, message, context)
    
    def error(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Error level logging"""
        self.log_event(EventType.ERROR_OCCURRED, LogLevel.ERROR, message, context)
    
    def critical(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Critical level logging"""
        self.log_event(EventType.ERROR_OCCURRED, LogLevel.CRITICAL, message, context)
    
    def performance(self, message: str, duration: float, context: Optional[Dict[str, Any]] = None):
        """Performance level logging"""
        self.log_event(EventType.PERFORMANCE_METRIC, LogLevel.PERFORMANCE, 
                      message, context, duration)
    
    def workflow(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Workflow level logging"""
        self.log_event(EventType.STEP_START, LogLevel.WORKFLOW, message, context)
    
    def scientific(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Scientific level logging"""
        self.log_event(EventType.VALIDATION_RESULT, LogLevel.SCIENTIFIC, message, context)
    
    @contextmanager
    def timed_operation(self, operation_name: str, context: Optional[Dict[str, Any]] = None):
        """Context manager for timing operations"""
        start_time = time.time()
        self.log_event(EventType.STEP_START, LogLevel.INFO, 
                      f"Starting {operation_name}", context)
        
        try:
            yield
            duration = time.time() - start_time
            self.log_event(EventType.STEP_END, LogLevel.INFO, 
                          f"Completed {operation_name}", context, duration)
        except Exception as e:
            duration = time.time() - start_time
            error_context = (context or {}).copy()
            error_context['error'] = str(e)
            error_context['duration'] = duration
            self.log_event(EventType.ERROR_OCCURRED, LogLevel.ERROR, 
                          f"Failed {operation_name}: {e}", error_context, duration)
            raise
    
    def log_external_tool(self, tool_name: str, command: str, exit_code: int,
                         stdout: str = "", stderr: str = "", duration: float = 0):
        """Log external tool execution"""
        context = {
            'tool_name': tool_name,
            'command': command,
            'exit_code': exit_code,
            'stdout_length': len(stdout),
            'stderr_length': len(stderr)
        }
        
        level = LogLevel.INFO if exit_code == 0 else LogLevel.ERROR
        event_type = EventType.EXTERNAL_TOOL
        
        message = f"{tool_name} {'succeeded' if exit_code == 0 else 'failed'} (exit code: {exit_code})"
        self.log_event(event_type, level, message, context, duration)
    
    def collect_performance_metrics(self) -> PerformanceMetrics:
        """Collect current system performance metrics"""
        try:
            process = psutil.Process()
            
            # CPU and memory
            cpu_percent = process.cpu_percent()
            memory_info = process.memory_info()
            memory_mb = memory_info.rss / 1024 / 1024
            
            # I/O
            io_counters = process.io_counters()
            disk_io_read_mb = io_counters.read_bytes / 1024 / 1024
            disk_io_write_mb = io_counters.write_bytes / 1024 / 1024
            
            # GPU metrics (if available)
            gpu_percent = None
            gpu_memory_mb = None
            try:
                import GPUtil
                gpus = GPUtil.getGPUs()
                if gpus:
                    gpu = gpus[0]  # Use first GPU
                    gpu_percent = gpu.load * 100
                    gpu_memory_mb = gpu.memoryUsed
            except ImportError:
                pass
            
            metrics = PerformanceMetrics(
                cpu_percent=cpu_percent,
                memory_mb=memory_mb,
                disk_io_read_mb=disk_io_read_mb,
                disk_io_write_mb=disk_io_write_mb,
                gpu_percent=gpu_percent,
                gpu_memory_mb=gpu_memory_mb
            )
            
            self.performance_history.append(metrics)
            
            # Log if metrics are concerning
            if cpu_percent > 90:
                self.warning(f"High CPU usage: {cpu_percent:.1f}%")
            if memory_mb > 8000:  # 8GB threshold
                self.warning(f"High memory usage: {memory_mb:.1f} MB")
            
            return metrics
            
        except Exception as e:
            self.error(f"Failed to collect performance metrics: {e}")
            return PerformanceMetrics(0, 0, 0, 0)
    
    def get_session_summary(self) -> Dict[str, Any]:
        """Get session summary statistics"""
        now = time.time()
        session_duration = now - self.start_time
        
        # Event statistics
        event_counts = {}
        for event in self.events:
            event_type = event.event_type.value
            event_counts[event_type] = event_counts.get(event_type, 0) + 1
        
        # Level statistics
        level_counts = {}
        for event in self.events:
            level = event.level.value
            level_counts[level] = level_counts.get(level, 0) + 1
        
        # Performance summary
        performance_summary = {}
        if self.performance_history:
            cpu_values = [m.cpu_percent for m in self.performance_history]
            memory_values = [m.memory_mb for m in self.performance_history]
            
            performance_summary = {
                'avg_cpu_percent': sum(cpu_values) / len(cpu_values),
                'max_cpu_percent': max(cpu_values),
                'avg_memory_mb': sum(memory_values) / len(memory_values),
                'max_memory_mb': max(memory_values),
                'metrics_collected': len(self.performance_history)
            }
        
        return {
            'session_id': self.session_id,
            'logger_name': self.name,
            'session_duration_seconds': session_duration,
            'total_events': len(self.events),
            'event_counts_by_type': event_counts,
            'event_counts_by_level': level_counts,
            'performance_summary': performance_summary,
            'start_time': datetime.fromtimestamp(self.start_time).isoformat(),
            'end_time': datetime.now().isoformat()
        }
    
    def export_events(self, output_file: Union[str, Path], format: str = "json"):
        """Export events to file"""
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        if format == "json":
            events_data = [event.to_dict() for event in self.events]
            with open(output_file, 'w') as f:
                json.dump(events_data, f, indent=2, default=str)
        elif format == "csv":
            import csv
            with open(output_file, 'w', newline='') as f:
                if self.events:
                    writer = csv.DictWriter(f, fieldnames=self.events[0].to_dict().keys())
                    writer.writeheader()
                    for event in self.events:
                        writer.writerow(event.to_dict())
        
        self.info(f"Exported {len(self.events)} events to {output_file}")


class MonitoringManager:
    """Real-time monitoring manager for PRISM simulations"""
    
    def __init__(self, logger: PrismLogger, monitoring_interval: float = 30.0):
        self.logger = logger
        self.monitoring_interval = monitoring_interval
        self.monitoring_active = False
        self.monitoring_thread = None
        self.callbacks: List[Callable[[PerformanceMetrics], None]] = []
        
    def add_callback(self, callback: Callable[[PerformanceMetrics], None]):
        """Add callback for performance metrics"""
        self.callbacks.append(callback)
    
    def start_monitoring(self):
        """Start background monitoring"""
        if self.monitoring_active:
            return
        
        self.monitoring_active = True
        self.monitoring_thread = threading.Thread(target=self._monitoring_loop, daemon=True)
        self.monitoring_thread.start()
        self.logger.info("Started performance monitoring")
    
    def stop_monitoring(self):
        """Stop background monitoring"""
        self.monitoring_active = False
        if self.monitoring_thread:
            self.monitoring_thread.join(timeout=5)
        self.logger.info("Stopped performance monitoring")
    
    def _monitoring_loop(self):
        """Background monitoring loop"""
        while self.monitoring_active:
            try:
                metrics = self.logger.collect_performance_metrics()
                
                # Log performance metrics
                self.logger.log_event(
                    EventType.RESOURCE_USAGE,
                    LogLevel.PERFORMANCE,
                    "System resource usage",
                    metrics.to_dict()
                )
                
                # Call callbacks
                for callback in self.callbacks:
                    try:
                        callback(metrics)
                    except Exception as e:
                        self.logger.error(f"Callback error: {e}")
                
                time.sleep(self.monitoring_interval)
                
            except Exception as e:
                self.logger.error(f"Monitoring error: {e}")
                time.sleep(self.monitoring_interval)


class LogAnalyzer:
    """Analyzer for PRISM log data"""
    
    def __init__(self, log_file: Union[str, Path]):
        self.log_file = Path(log_file)
        self.events = self._load_events()
    
    def _load_events(self) -> List[Dict[str, Any]]:
        """Load events from log file"""
        events = []
        if not self.log_file.exists():
            return events
        
        try:
            with open(self.log_file, 'r') as f:
                for line in f:
                    try:
                        event = json.loads(line.strip())
                        events.append(event)
                    except json.JSONDecodeError:
                        continue
            return events
        except Exception:
            return events
    
    def get_performance_report(self) -> Dict[str, Any]:
        """Generate performance analysis report"""
        if not self.events:
            return {"error": "No events found"}
        
        # Filter performance events
        perf_events = [e for e in self.events if e.get('event_type') == 'performance_metric']
        
        if not perf_events:
            return {"error": "No performance events found"}
        
        # Analyze durations
        durations = [e.get('duration', 0) for e in perf_events if e.get('duration')]
        
        if durations:
            avg_duration = sum(durations) / len(durations)
            max_duration = max(durations)
            min_duration = min(durations)
        else:
            avg_duration = max_duration = min_duration = 0
        
        return {
            'total_performance_events': len(perf_events),
            'average_duration': avg_duration,
            'max_duration': max_duration,
            'min_duration': min_duration,
            'total_events_analyzed': len(self.events)
        }
    
    def get_error_summary(self) -> Dict[str, Any]:
        """Get error summary from logs"""
        error_events = [e for e in self.events if e.get('level') in ['ERROR', 'CRITICAL']]
        
        # Group errors by message pattern
        error_patterns = {}
        for event in error_events:
            message = event.get('message', 'Unknown error')
            # Simplify message for grouping
            pattern = message.split(':')[0] if ':' in message else message
            error_patterns[pattern] = error_patterns.get(pattern, 0) + 1
        
        return {
            'total_errors': len(error_events),
            'error_patterns': error_patterns,
            'recent_errors': error_events[-5:] if error_events else []
        }


# Convenience functions
def get_prism_logger(name: str, log_dir: Optional[Union[str, Path]] = None,
                    enable_monitoring: bool = True) -> PrismLogger:
    """Get configured PRISM logger"""
    logger = PrismLogger(name)
    
    if log_dir:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        
        # Add structured log file
        structured_log = log_dir / f"{name}_structured.jsonl"
        logger.add_file_handler(structured_log, structured=True)
        
        # Add human-readable log file
        readable_log = log_dir / f"{name}.log"
        logger.add_file_handler(readable_log, structured=False)
    
    if enable_monitoring:
        monitor = MonitoringManager(logger)
        monitor.start_monitoring()
        # Store monitor reference to prevent garbage collection
        logger._monitor = monitor
    
    return logger


def setup_prism_logging(log_dir: Union[str, Path], log_level: str = "INFO") -> PrismLogger:
    """Setup PRISM logging system"""
    log_dir = Path(log_dir)
    
    # Configure root logger
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Create main PRISM logger
    main_logger = get_prism_logger("prism", log_dir, enable_monitoring=True)
    main_logger.info("PRISM logging system initialized")
    
    return main_logger


@contextmanager
def log_simulation_session(logger: PrismLogger, simulation_type: str, 
                          config: Optional[Dict[str, Any]] = None):
    """Context manager for logging simulation sessions"""
    start_time = time.time()
    
    context = {
        'simulation_type': simulation_type,
        'session_id': logger.session_id
    }
    if config:
        context['config_summary'] = {
            'sections': list(config.keys()),
            'total_parameters': sum(len(v) if isinstance(v, dict) else 1 for v in config.values())
        }
    
    logger.log_event(EventType.SIMULATION_START, LogLevel.WORKFLOW, 
                    f"Starting {simulation_type} simulation", context)
    
    try:
        yield logger
        
        duration = time.time() - start_time
        logger.log_event(EventType.SIMULATION_END, LogLevel.WORKFLOW, 
                        f"Completed {simulation_type} simulation", context, duration)
        
        # Generate session summary
        summary = logger.get_session_summary()
        logger.scientific(f"Session completed successfully", summary)
        
    except Exception as e:
        duration = time.time() - start_time
        error_context = context.copy()
        error_context['error'] = str(e)
        error_context['duration'] = duration
        
        logger.log_event(EventType.ERROR_OCCURRED, LogLevel.CRITICAL, 
                        f"Simulation failed: {e}", error_context, duration)
        raise