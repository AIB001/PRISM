#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Workflow Orchestrator System

Manages distributed workflow execution across multiple nodes and clusters,
providing load balancing, fault tolerance, and resource optimization.
"""

import time
import json
import socket
import threading
import concurrent.futures
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Callable, Tuple, Set
from dataclasses import dataclass, field
from enum import Enum
from abc import ABC, abstractmethod
from collections import defaultdict, deque
import heapq

from ..utils.logging_system import PrismLogger
from .scheduler import TaskScheduler, ResourceManager, ScheduledTask, TaskStatus
from .workflow import Workflow, WorkflowStatus, Task
from .executor import WorkflowExecutor, ExecutionContext, ExecutionResult


class NodeStatus(Enum):
    """Node status in the cluster"""
    ACTIVE = "active"
    BUSY = "busy"
    OFFLINE = "offline"
    MAINTENANCE = "maintenance"
    FAILED = "failed"


class LoadBalancingStrategy(Enum):
    """Load balancing strategies"""
    ROUND_ROBIN = "round_robin"
    LEAST_LOADED = "least_loaded"
    RESOURCE_AWARE = "resource_aware"
    LATENCY_AWARE = "latency_aware"
    CUSTOM = "custom"


@dataclass
class ComputeNode:
    """Represents a compute node in the cluster"""
    node_id: str
    hostname: str
    port: int = 8080
    
    # Node capabilities
    total_resources: Dict[str, float] = field(default_factory=dict)
    available_resources: Dict[str, float] = field(default_factory=dict)
    
    # Node state
    status: NodeStatus = NodeStatus.OFFLINE
    last_heartbeat: float = 0.0
    load_factor: float = 0.0
    
    # Execution tracking
    active_workflows: Set[str] = field(default_factory=set)
    active_tasks: Set[str] = field(default_factory=set)
    completed_workflows: int = 0
    failed_workflows: int = 0
    
    # Performance metrics
    avg_task_duration: float = 0.0
    network_latency: float = 0.0
    success_rate: float = 1.0
    
    @property
    def address(self) -> str:
        return f"{self.hostname}:{self.port}"
    
    @property
    def resource_utilization(self) -> float:
        """Calculate overall resource utilization"""
        if not self.total_resources:
            return 0.0
        
        utilization = 0.0
        count = 0
        
        for resource_type, total in self.total_resources.items():
            if total > 0:
                available = self.available_resources.get(resource_type, total)
                utilization += 1.0 - (available / total)
                count += 1
        
        return utilization / count if count > 0 else 0.0
    
    def can_accept_workflow(self, workflow: Workflow) -> bool:
        """Check if node can accept a workflow"""
        if self.status != NodeStatus.ACTIVE:
            return False
        
        # Check resource requirements
        peak_resources = workflow.get_resource_requirements()
        
        return (
            self.available_resources.get('cpu_cores', 0) >= peak_resources.cpu_cores and
            self.available_resources.get('memory_gb', 0) >= peak_resources.memory_gb and
            self.resource_utilization < 0.9  # Don't overload
        )
    
    def update_resources(self, allocated_resources: Dict[str, float]):
        """Update available resources after allocation"""
        for resource_type, amount in allocated_resources.items():
            current = self.available_resources.get(resource_type, 0)
            self.available_resources[resource_type] = max(0, current - amount)
    
    def release_resources(self, released_resources: Dict[str, float]):
        """Release resources after task completion"""
        for resource_type, amount in released_resources.items():
            current = self.available_resources.get(resource_type, 0)
            total = self.total_resources.get(resource_type, current + amount)
            self.available_resources[resource_type] = min(total, current + amount)


class ClusterManager:
    """Manages a cluster of compute nodes"""
    
    def __init__(self, cluster_id: str = "default"):
        self.cluster_id = cluster_id
        self.logger = PrismLogger(f"cluster_manager.{cluster_id}")
        
        # Node management
        self.nodes: Dict[str, ComputeNode] = {}
        self.load_balancer: Optional['LoadBalancer'] = None
        
        # Cluster state
        self.heartbeat_interval = 30.0  # seconds
        self.node_timeout = 120.0  # seconds
        
        # Monitoring and statistics
        self.cluster_stats = {
            'total_nodes': 0,
            'active_nodes': 0,
            'total_workflows': 0,
            'completed_workflows': 0,
            'failed_workflows': 0,
            'avg_node_utilization': 0.0
        }
        
        # Threading
        self.lock = threading.RLock()
        self.monitor_thread: Optional[threading.Thread] = None
        self.running = False
    
    def add_node(self, node: ComputeNode) -> bool:
        """Add a node to the cluster"""
        with self.lock:
            if node.node_id in self.nodes:
                self.logger.warning(f"Node {node.node_id} already exists")
                return False
            
            self.nodes[node.node_id] = node
            node.status = NodeStatus.ACTIVE
            node.last_heartbeat = time.time()
            
            self._update_cluster_stats()
            self.logger.info(f"Added node {node.node_id} ({node.address})")
            return True
    
    def remove_node(self, node_id: str) -> bool:
        """Remove a node from the cluster"""
        with self.lock:
            if node_id not in self.nodes:
                return False
            
            node = self.nodes[node_id]
            node.status = NodeStatus.OFFLINE
            
            # Reassign active workflows if any
            if node.active_workflows:
                self.logger.warning(f"Node {node_id} removed with {len(node.active_workflows)} active workflows")
                # In a real system, you would reassign these workflows
            
            del self.nodes[node_id]
            self._update_cluster_stats()
            self.logger.info(f"Removed node {node_id}")
            return True
    
    def get_available_nodes(self) -> List[ComputeNode]:
        """Get list of available nodes"""
        with self.lock:
            return [node for node in self.nodes.values() 
                   if node.status == NodeStatus.ACTIVE]
    
    def select_node_for_workflow(self, workflow: Workflow, 
                                strategy: LoadBalancingStrategy = LoadBalancingStrategy.RESOURCE_AWARE) -> Optional[ComputeNode]:
        """Select the best node for executing a workflow"""
        available_nodes = self.get_available_nodes()
        
        if not available_nodes:
            return None
        
        # Filter nodes that can handle the workflow
        capable_nodes = [node for node in available_nodes 
                        if node.can_accept_workflow(workflow)]
        
        if not capable_nodes:
            return None
        
        # Select based on strategy
        if strategy == LoadBalancingStrategy.ROUND_ROBIN:
            return self._round_robin_selection(capable_nodes)
        elif strategy == LoadBalancingStrategy.LEAST_LOADED:
            return self._least_loaded_selection(capable_nodes)
        elif strategy == LoadBalancingStrategy.RESOURCE_AWARE:
            return self._resource_aware_selection(capable_nodes, workflow)
        elif strategy == LoadBalancingStrategy.LATENCY_AWARE:
            return self._latency_aware_selection(capable_nodes)
        else:
            return capable_nodes[0]  # Default to first available
    
    def update_node_heartbeat(self, node_id: str, metrics: Optional[Dict[str, Any]] = None):
        """Update node heartbeat and metrics"""
        with self.lock:
            if node_id not in self.nodes:
                return
            
            node = self.nodes[node_id]
            node.last_heartbeat = time.time()
            
            if metrics:
                # Update node metrics
                node.load_factor = metrics.get('load_factor', node.load_factor)
                node.available_resources.update(metrics.get('available_resources', {}))
                node.avg_task_duration = metrics.get('avg_task_duration', node.avg_task_duration)
                node.network_latency = metrics.get('network_latency', node.network_latency)
    
    def start_monitoring(self):
        """Start cluster monitoring"""
        if self.running:
            return
        
        self.running = True
        
        def monitor_cluster():
            while self.running:
                try:
                    self._check_node_health()
                    self._update_cluster_stats()
                    time.sleep(self.heartbeat_interval)
                except Exception as e:
                    self.logger.error(f"Cluster monitoring error: {e}")
        
        self.monitor_thread = threading.Thread(target=monitor_cluster, daemon=True)
        self.monitor_thread.start()
        self.logger.info("Started cluster monitoring")
    
    def stop_monitoring(self):
        """Stop cluster monitoring"""
        self.running = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=10.0)
    
    def _round_robin_selection(self, nodes: List[ComputeNode]) -> ComputeNode:
        """Round-robin node selection"""
        # Simple round-robin based on completed workflows
        return min(nodes, key=lambda n: n.completed_workflows)
    
    def _least_loaded_selection(self, nodes: List[ComputeNode]) -> ComputeNode:
        """Select least loaded node"""
        return min(nodes, key=lambda n: n.resource_utilization)
    
    def _resource_aware_selection(self, nodes: List[ComputeNode], workflow: Workflow) -> ComputeNode:
        """Select node based on resource fit"""
        peak_resources = workflow.get_resource_requirements()
        
        def resource_score(node: ComputeNode) -> float:
            # Higher score for better resource fit
            cpu_fit = node.available_resources.get('cpu_cores', 0) / max(1, peak_resources.cpu_cores)
            mem_fit = node.available_resources.get('memory_gb', 0) / max(1, peak_resources.memory_gb)
            
            # Prefer nodes that are not overloaded but can handle the workload
            utilization_penalty = node.resource_utilization
            
            return (cpu_fit + mem_fit) / 2.0 - utilization_penalty
        
        return max(nodes, key=resource_score)
    
    def _latency_aware_selection(self, nodes: List[ComputeNode]) -> ComputeNode:
        """Select node with lowest latency"""
        return min(nodes, key=lambda n: n.network_latency)
    
    def _check_node_health(self):
        """Check health of all nodes"""
        current_time = time.time()
        
        with self.lock:
            for node in list(self.nodes.values()):
                time_since_heartbeat = current_time - node.last_heartbeat
                
                if time_since_heartbeat > self.node_timeout:
                    if node.status == NodeStatus.ACTIVE:
                        node.status = NodeStatus.OFFLINE
                        self.logger.warning(f"Node {node.node_id} marked offline (no heartbeat)")
    
    def _update_cluster_stats(self):
        """Update cluster statistics"""
        with self.lock:
            active_nodes = [n for n in self.nodes.values() if n.status == NodeStatus.ACTIVE]
            
            self.cluster_stats.update({
                'total_nodes': len(self.nodes),
                'active_nodes': len(active_nodes),
                'total_workflows': sum(n.completed_workflows + n.failed_workflows for n in self.nodes.values()),
                'completed_workflows': sum(n.completed_workflows for n in self.nodes.values()),
                'failed_workflows': sum(n.failed_workflows for n in self.nodes.values()),
                'avg_node_utilization': sum(n.resource_utilization for n in active_nodes) / len(active_nodes) if active_nodes else 0.0
            })
    
    def get_cluster_status(self) -> Dict[str, Any]:
        """Get cluster status and statistics"""
        with self.lock:
            return {
                'cluster_id': self.cluster_id,
                'stats': self.cluster_stats.copy(),
                'nodes': {
                    node_id: {
                        'node_id': node.node_id,
                        'address': node.address,
                        'status': node.status.value,
                        'resource_utilization': node.resource_utilization,
                        'active_workflows': len(node.active_workflows),
                        'active_tasks': len(node.active_tasks),
                        'success_rate': node.success_rate,
                        'last_heartbeat': node.last_heartbeat
                    }
                    for node_id, node in self.nodes.items()
                }
            }


class WorkflowOrchestrator:
    """Central orchestrator for managing distributed workflow execution"""
    
    def __init__(self, cluster_manager: ClusterManager):
        self.cluster_manager = cluster_manager
        self.logger = PrismLogger("workflow_orchestrator")
        
        # Workflow management
        self.active_workflows: Dict[str, Dict[str, Any]] = {}
        self.workflow_history: List[Dict[str, Any]] = []
        
        # Load balancing
        self.load_balancing_strategy = LoadBalancingStrategy.RESOURCE_AWARE
        
        # Threading
        self.lock = threading.RLock()
        self.executor_pool = concurrent.futures.ThreadPoolExecutor(max_workers=10)
    
    def submit_workflow(self, workflow: Workflow, 
                       execution_context: Optional[ExecutionContext] = None,
                       callback: Optional[Callable] = None) -> Optional[str]:
        """Submit a workflow for distributed execution"""
        
        # Select a node for execution
        selected_node = self.cluster_manager.select_node_for_workflow(
            workflow, self.load_balancing_strategy
        )
        
        if not selected_node:
            self.logger.error(f"No available node for workflow {workflow.workflow_id}")
            return None
        
        # Create execution context if not provided
        if execution_context is None:
            execution_context = ExecutionContext(
                workflow_id=workflow.workflow_id,
                execution_id=f"{workflow.workflow_id}_{int(time.time())}",
                working_directory=Path.cwd()
            )
        
        # Track workflow
        workflow_info = {
            'workflow': workflow,
            'execution_context': execution_context,
            'assigned_node': selected_node.node_id,
            'status': 'submitted',
            'start_time': time.time(),
            'callback': callback
        }
        
        with self.lock:
            self.active_workflows[execution_context.execution_id] = workflow_info
            selected_node.active_workflows.add(execution_context.execution_id)
        
        # Submit for execution
        future = self.executor_pool.submit(
            self._execute_workflow_on_node, 
            workflow, execution_context, selected_node
        )
        
        if callback:
            future.add_done_callback(lambda f: callback(f.result()))
        
        self.logger.info(f"Submitted workflow {workflow.workflow_id} to node {selected_node.node_id}")
        return execution_context.execution_id
    
    def _execute_workflow_on_node(self, workflow: Workflow, 
                                 context: ExecutionContext,
                                 node: ComputeNode) -> ExecutionResult:
        """Execute workflow on a specific node"""
        
        try:
            # In a real distributed system, this would send the workflow to the remote node
            # For now, we'll simulate local execution
            
            # Update workflow status
            with self.lock:
                if context.execution_id in self.active_workflows:
                    self.active_workflows[context.execution_id]['status'] = 'running'
            
            # Allocate resources on the node
            peak_resources = workflow.get_resource_requirements()
            resource_dict = {
                'cpu_cores': peak_resources.cpu_cores,
                'memory_gb': peak_resources.memory_gb,
                'gpu_devices': peak_resources.gpu_devices
            }
            node.update_resources(resource_dict)
            
            # Create local executor for simulation
            local_executor = WorkflowExecutor()
            result = local_executor._execute_workflow_sync(workflow, context)
            
            # Update node statistics
            node.completed_workflows += 1 if result.status == WorkflowStatus.COMPLETED else 0
            node.failed_workflows += 1 if result.status == WorkflowStatus.FAILED else 0
            
            # Update success rate
            total = node.completed_workflows + node.failed_workflows
            node.success_rate = node.completed_workflows / total if total > 0 else 1.0
            
            # Release resources
            node.release_resources(resource_dict)
            
            # Remove from active workflows
            with self.lock:
                if context.execution_id in self.active_workflows:
                    workflow_info = self.active_workflows[context.execution_id]
                    workflow_info['status'] = 'completed'
                    workflow_info['end_time'] = time.time()
                    workflow_info['result'] = result
                    
                    # Move to history
                    self.workflow_history.append(workflow_info)
                    del self.active_workflows[context.execution_id]
                
                node.active_workflows.discard(context.execution_id)
            
            self.logger.info(f"Workflow {workflow.workflow_id} completed on node {node.node_id}")
            return result
            
        except Exception as e:
            self.logger.error(f"Workflow execution failed on node {node.node_id}: {e}")
            
            # Handle failure
            node.failed_workflows += 1
            node.release_resources(resource_dict)
            
            with self.lock:
                if context.execution_id in self.active_workflows:
                    workflow_info = self.active_workflows[context.execution_id]
                    workflow_info['status'] = 'failed'
                    workflow_info['error'] = str(e)
                    workflow_info['end_time'] = time.time()
                    
                    self.workflow_history.append(workflow_info)
                    del self.active_workflows[context.execution_id]
                
                node.active_workflows.discard(context.execution_id)
            
            # Return failed result
            return ExecutionResult(
                workflow_id=workflow.workflow_id,
                execution_id=context.execution_id,
                status=WorkflowStatus.FAILED,
                start_time=context.start_time or time.time(),
                errors=[str(e)]
            )
    
    def get_execution_status(self, execution_id: str) -> Optional[Dict[str, Any]]:
        """Get status of workflow execution"""
        with self.lock:
            if execution_id in self.active_workflows:
                workflow_info = self.active_workflows[execution_id]
                return {
                    'execution_id': execution_id,
                    'status': workflow_info['status'],
                    'assigned_node': workflow_info['assigned_node'],
                    'start_time': workflow_info['start_time'],
                    'workflow_id': workflow_info['workflow'].workflow_id
                }
            
            # Check history
            for workflow_info in reversed(self.workflow_history):
                if workflow_info['execution_context'].execution_id == execution_id:
                    return {
                        'execution_id': execution_id,
                        'status': workflow_info['status'],
                        'assigned_node': workflow_info['assigned_node'],
                        'start_time': workflow_info['start_time'],
                        'end_time': workflow_info.get('end_time'),
                        'workflow_id': workflow_info['workflow'].workflow_id,
                        'result': workflow_info.get('result')
                    }
        
        return None
    
    def cancel_workflow(self, execution_id: str) -> bool:
        """Cancel a workflow execution"""
        with self.lock:
            if execution_id not in self.active_workflows:
                return False
            
            workflow_info = self.active_workflows[execution_id]
            workflow_info['status'] = 'cancelled'
            
            # Remove from node tracking
            node_id = workflow_info['assigned_node']
            if node_id in self.cluster_manager.nodes:
                node = self.cluster_manager.nodes[node_id]
                node.active_workflows.discard(execution_id)
            
            # Move to history
            workflow_info['end_time'] = time.time()
            self.workflow_history.append(workflow_info)
            del self.active_workflows[execution_id]
            
            self.logger.info(f"Cancelled workflow execution {execution_id}")
            return True
    
    def get_orchestrator_stats(self) -> Dict[str, Any]:
        """Get orchestrator statistics"""
        with self.lock:
            total_workflows = len(self.workflow_history)
            completed_workflows = sum(1 for w in self.workflow_history if w['status'] == 'completed')
            failed_workflows = sum(1 for w in self.workflow_history if w['status'] == 'failed')
            
            return {
                'active_workflows': len(self.active_workflows),
                'total_workflows_submitted': total_workflows,
                'completed_workflows': completed_workflows,
                'failed_workflows': failed_workflows,
                'success_rate': completed_workflows / total_workflows if total_workflows > 0 else 0.0,
                'load_balancing_strategy': self.load_balancing_strategy.value,
                'cluster_stats': self.cluster_manager.get_cluster_status()
            }
    
    def shutdown(self):
        """Shutdown the orchestrator"""
        self.logger.info("Shutting down workflow orchestrator")
        self.executor_pool.shutdown(wait=True, timeout=30.0)


class DistributedOrchestrator(WorkflowOrchestrator):
    """Extended orchestrator with advanced distributed features"""
    
    def __init__(self, cluster_manager: ClusterManager):
        super().__init__(cluster_manager)
        
        # Advanced features
        self.workflow_queues: Dict[str, deque] = defaultdict(deque)  # Priority queues
        self.resource_predictor: Optional['ResourcePredictor'] = None
        self.fault_tolerance_enabled = True
        
    def enable_auto_scaling(self, min_nodes: int = 1, max_nodes: int = 10):
        """Enable auto-scaling of cluster nodes"""
        # This would integrate with cloud APIs to automatically scale nodes
        self.logger.info(f"Auto-scaling enabled: {min_nodes}-{max_nodes} nodes")
    
    def enable_workflow_checkpointing(self, checkpoint_interval: float = 300.0):
        """Enable workflow checkpointing for fault tolerance"""
        self.checkpoint_interval = checkpoint_interval
        self.logger.info(f"Workflow checkpointing enabled (interval: {checkpoint_interval}s)")
    
    def submit_batch_workflows(self, workflows: List[Workflow], 
                              priority: str = "normal") -> List[str]:
        """Submit multiple workflows as a batch"""
        execution_ids = []
        
        for workflow in workflows:
            execution_id = self.submit_workflow(workflow)
            if execution_id:
                execution_ids.append(execution_id)
        
        self.logger.info(f"Submitted batch of {len(execution_ids)} workflows")
        return execution_ids


# Convenience functions
def create_single_node_orchestrator() -> WorkflowOrchestrator:
    """Create orchestrator with single local node"""
    cluster_manager = ClusterManager("local")
    
    # Add local node
    local_node = ComputeNode(
        node_id="local_node",
        hostname="localhost",
        total_resources={
            'cpu_cores': 8.0,
            'memory_gb': 16.0,
            'gpu_devices': 0.0,
            'disk_space_gb': 100.0
        },
        available_resources={
            'cpu_cores': 8.0,
            'memory_gb': 16.0,
            'gpu_devices': 0.0,
            'disk_space_gb': 100.0
        }
    )
    
    cluster_manager.add_node(local_node)
    cluster_manager.start_monitoring()
    
    return WorkflowOrchestrator(cluster_manager)


def create_distributed_orchestrator(nodes: List[Dict[str, Any]]) -> DistributedOrchestrator:
    """Create distributed orchestrator with multiple nodes"""
    cluster_manager = ClusterManager("distributed")
    
    for node_config in nodes:
        node = ComputeNode(
            node_id=node_config['node_id'],
            hostname=node_config['hostname'],
            port=node_config.get('port', 8080),
            total_resources=node_config.get('resources', {}),
            available_resources=node_config.get('resources', {}).copy()
        )
        cluster_manager.add_node(node)
    
    cluster_manager.start_monitoring()
    return DistributedOrchestrator(cluster_manager)