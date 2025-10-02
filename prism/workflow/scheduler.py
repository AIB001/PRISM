#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Task Scheduling System

Intelligent task scheduler with resource management, priority queuing,
and load balancing for molecular dynamics workflows.
"""

import os
import sys
import time
import threading
import multiprocessing as mp
import psutil
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Callable, Tuple, Set
from dataclasses import dataclass, field, asdict
from enum import Enum, IntEnum
from queue import PriorityQueue, Queue, Empty
import heapq
from collections import defaultdict, deque
from abc import ABC, abstractmethod

from ..utils.logging_system import PrismLogger


class TaskPriority(IntEnum):
    """Task priority levels (lower value = higher priority)"""
    CRITICAL = 0
    HIGH = 1
    NORMAL = 2
    LOW = 3
    BACKGROUND = 4


class TaskStatus(Enum):
    """Task execution status"""
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    TIMEOUT = "timeout"


class ResourceType(Enum):
    """Available resource types"""
    CPU_CORES = "cpu_cores"
    MEMORY_GB = "memory_gb"
    GPU_DEVICES = "gpu_devices"
    DISK_SPACE_GB = "disk_space_gb"
    NETWORK_BANDWIDTH = "network_bandwidth"


class SchedulingPolicy(Enum):
    """Task scheduling policies"""
    FIFO = "fifo"                    # First In, First Out
    PRIORITY = "priority"            # Priority-based scheduling
    ROUND_ROBIN = "round_robin"      # Round-robin scheduling
    SHORTEST_JOB_FIRST = "sjf"       # Shortest Job First
    FAIR_SHARE = "fair_share"        # Fair resource sharing
    RESOURCE_AWARE = "resource_aware" # Resource-aware scheduling


@dataclass
class ResourceRequirements:
    """Resource requirements for a task"""
    cpu_cores: int = 1
    memory_gb: float = 1.0
    gpu_devices: int = 0
    disk_space_gb: float = 1.0
    network_bandwidth: float = 0.0
    estimated_runtime_minutes: float = 60.0
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
    
    def __le__(self, other: 'ResourceRequirements') -> bool:
        """Check if this resource requirement is less than or equal to another"""
        return (
            self.cpu_cores <= other.cpu_cores and
            self.memory_gb <= other.memory_gb and
            self.gpu_devices <= other.gpu_devices and
            self.disk_space_gb <= other.disk_space_gb
        )
    
    def __add__(self, other: 'ResourceRequirements') -> 'ResourceRequirements':
        """Add two resource requirements"""
        return ResourceRequirements(
            cpu_cores=self.cpu_cores + other.cpu_cores,
            memory_gb=self.memory_gb + other.memory_gb,
            gpu_devices=self.gpu_devices + other.gpu_devices,
            disk_space_gb=self.disk_space_gb + other.disk_space_gb,
            network_bandwidth=self.network_bandwidth + other.network_bandwidth,
            estimated_runtime_minutes=max(self.estimated_runtime_minutes, 
                                        other.estimated_runtime_minutes)
        )
    
    def __sub__(self, other: 'ResourceRequirements') -> 'ResourceRequirements':
        """Subtract resource requirements"""
        return ResourceRequirements(
            cpu_cores=max(0, self.cpu_cores - other.cpu_cores),
            memory_gb=max(0.0, self.memory_gb - other.memory_gb),
            gpu_devices=max(0, self.gpu_devices - other.gpu_devices),
            disk_space_gb=max(0.0, self.disk_space_gb - other.disk_space_gb),
            network_bandwidth=max(0.0, self.network_bandwidth - other.network_bandwidth),
            estimated_runtime_minutes=self.estimated_runtime_minutes
        )


@dataclass
class ScheduledTask:
    """A task with scheduling metadata"""
    task_id: str
    command: str
    args: List[str] = field(default_factory=list)
    kwargs: Dict[str, Any] = field(default_factory=dict)
    priority: TaskPriority = TaskPriority.NORMAL
    resources: ResourceRequirements = field(default_factory=ResourceRequirements)
    dependencies: Set[str] = field(default_factory=set)
    max_retries: int = 3
    timeout_minutes: float = 120.0
    
    # Scheduling metadata
    status: TaskStatus = TaskStatus.PENDING
    submitted_time: float = field(default_factory=time.time)
    queued_time: Optional[float] = None
    started_time: Optional[float] = None
    completed_time: Optional[float] = None
    retry_count: int = 0
    assigned_worker: Optional[str] = None
    
    # Results
    result: Optional[Any] = None
    error_message: Optional[str] = None
    execution_metrics: Dict[str, Any] = field(default_factory=dict)
    
    def __lt__(self, other: 'ScheduledTask') -> bool:
        """Priority comparison for queue ordering"""
        if self.priority != other.priority:
            return self.priority < other.priority
        return self.submitted_time < other.submitted_time
    
    @property
    def wait_time(self) -> float:
        """Calculate wait time in queue"""
        if self.started_time and self.queued_time:
            return self.started_time - self.queued_time
        elif self.queued_time:
            return time.time() - self.queued_time
        return 0.0
    
    @property
    def execution_time(self) -> float:
        """Calculate execution time"""
        if self.completed_time and self.started_time:
            return self.completed_time - self.started_time
        elif self.started_time:
            return time.time() - self.started_time
        return 0.0


class ResourceManager:
    """Manages system resources and allocation"""
    
    def __init__(self):
        self.logger = PrismLogger("resource_manager")
        self.lock = threading.RLock()
        
        # System resources
        self.total_resources = self._detect_system_resources()
        self.available_resources = ResourceRequirements(**self.total_resources.to_dict())
        
        # Resource tracking
        self.allocated_resources: Dict[str, ResourceRequirements] = {}
        self.resource_history: List[Dict[str, Any]] = []
        
        # Resource limits
        self.max_cpu_utilization = 0.9
        self.max_memory_utilization = 0.85
        self.reservation_buffer = 0.1
        
    def _detect_system_resources(self) -> ResourceRequirements:
        """Detect available system resources"""
        try:
            cpu_cores = os.cpu_count() or 1
            memory_gb = psutil.virtual_memory().total / (1024**3)
            
            # Try to detect GPU devices
            gpu_devices = 0
            try:
                import GPUtil
                gpu_devices = len(GPUtil.getGPUs())
            except ImportError:
                pass
            
            # Disk space (current directory)
            disk_space_gb = psutil.disk_usage('.').free / (1024**3)
            
            resources = ResourceRequirements(
                cpu_cores=cpu_cores,
                memory_gb=memory_gb,
                gpu_devices=gpu_devices,
                disk_space_gb=disk_space_gb,
                network_bandwidth=1000.0  # Assume 1 Gbps
            )
            
            self.logger.info(f"Detected system resources: {resources}")
            return resources
            
        except Exception as e:
            self.logger.error(f"Failed to detect system resources: {e}")
            # Return conservative defaults
            return ResourceRequirements(cpu_cores=1, memory_gb=4.0)
    
    def can_allocate(self, requirements: ResourceRequirements) -> bool:
        """Check if resources can be allocated"""
        with self.lock:
            # Apply utilization limits
            max_cpu = int(self.total_resources.cpu_cores * self.max_cpu_utilization)
            max_memory = self.total_resources.memory_gb * self.max_memory_utilization
            
            return (
                requirements.cpu_cores <= max_cpu and
                requirements.memory_gb <= max_memory and
                requirements.gpu_devices <= self.available_resources.gpu_devices and
                requirements.disk_space_gb <= self.available_resources.disk_space_gb and
                requirements <= self.available_resources
            )
    
    def allocate_resources(self, task_id: str, requirements: ResourceRequirements) -> bool:
        """Allocate resources for a task"""
        with self.lock:
            if not self.can_allocate(requirements):
                return False
            
            self.allocated_resources[task_id] = requirements
            self.available_resources = self.available_resources - requirements
            
            self._record_resource_usage()
            
            self.logger.debug(f"Allocated resources for {task_id}: {requirements}")
            return True
    
    def release_resources(self, task_id: str) -> bool:
        """Release resources from a task"""
        with self.lock:
            if task_id not in self.allocated_resources:
                return False
            
            requirements = self.allocated_resources[task_id]
            del self.allocated_resources[task_id]
            
            self.available_resources = self.available_resources + requirements
            
            self._record_resource_usage()
            
            self.logger.debug(f"Released resources for {task_id}: {requirements}")
            return True
    
    def get_utilization_stats(self) -> Dict[str, Any]:
        """Get current resource utilization statistics"""
        with self.lock:
            total_allocated = ResourceRequirements()
            for req in self.allocated_resources.values():
                total_allocated = total_allocated + req
            
            return {
                'total_resources': self.total_resources.to_dict(),
                'available_resources': self.available_resources.to_dict(),
                'allocated_resources': total_allocated.to_dict(),
                'cpu_utilization': 1.0 - (self.available_resources.cpu_cores / self.total_resources.cpu_cores),
                'memory_utilization': 1.0 - (self.available_resources.memory_gb / self.total_resources.memory_gb),
                'active_tasks': len(self.allocated_resources)
            }
    
    def _record_resource_usage(self):
        """Record current resource usage for history"""
        stats = self.get_utilization_stats()
        stats['timestamp'] = time.time()
        self.resource_history.append(stats)
        
        # Keep only recent history (last 1000 records)
        if len(self.resource_history) > 1000:
            self.resource_history = self.resource_history[-1000:]


class TaskScheduler:
    """Intelligent task scheduler with multiple scheduling policies"""
    
    def __init__(self, resource_manager: ResourceManager, 
                 policy: SchedulingPolicy = SchedulingPolicy.PRIORITY):
        self.resource_manager = resource_manager
        self.policy = policy
        self.logger = PrismLogger("task_scheduler")
        
        # Task queues and tracking
        self.pending_queue = PriorityQueue()
        self.running_tasks: Dict[str, ScheduledTask] = {}
        self.completed_tasks: Dict[str, ScheduledTask] = {}
        self.failed_tasks: Dict[str, ScheduledTask] = {}
        
        # Scheduling state
        self.scheduler_thread: Optional[threading.Thread] = None
        self.running = False
        self.lock = threading.RLock()
        
        # Scheduling statistics
        self.stats = {
            'tasks_submitted': 0,
            'tasks_completed': 0,
            'tasks_failed': 0,
            'total_execution_time': 0.0,
            'total_wait_time': 0.0
        }
        
        # Fair share tracking (for fair share scheduling)
        self.user_usage: Dict[str, float] = defaultdict(float)
        self.user_shares: Dict[str, float] = defaultdict(lambda: 1.0)
    
    def submit_task(self, task: ScheduledTask) -> bool:
        """Submit a task for scheduling"""
        with self.lock:
            # Validate task
            if not self._validate_task(task):
                return False
            
            # Check dependencies
            if not self._check_dependencies(task):
                task.status = TaskStatus.PENDING
                self.logger.info(f"Task {task.task_id} waiting for dependencies")
                return True
            
            # Add to queue
            task.status = TaskStatus.QUEUED
            task.queued_time = time.time()
            self.pending_queue.put(task)
            
            self.stats['tasks_submitted'] += 1
            self.logger.info(f"Queued task {task.task_id} with priority {task.priority}")
            
            return True
    
    def start_scheduler(self, max_workers: Optional[int] = None):
        """Start the task scheduler"""
        if self.running:
            return
        
        self.running = True
        max_workers = max_workers or min(32, (os.cpu_count() or 1) * 2)
        
        def scheduler_worker():
            self.logger.info(f"Starting task scheduler with {max_workers} max workers")
            
            while self.running:
                try:
                    # Get next task to schedule
                    task = self._get_next_task()
                    if not task:
                        time.sleep(0.1)
                        continue
                    
                    # Try to allocate resources
                    if self.resource_manager.can_allocate(task.resources):
                        if self.resource_manager.allocate_resources(task.task_id, task.resources):
                            self._start_task_execution(task)
                        else:
                            # Put task back in queue
                            self.pending_queue.put(task)
                    else:
                        # Resource not available, wait
                        self.pending_queue.put(task)
                        time.sleep(1.0)
                    
                except Exception as e:
                    self.logger.error(f"Scheduler error: {e}")
                    time.sleep(1.0)
        
        self.scheduler_thread = threading.Thread(target=scheduler_worker, daemon=True)
        self.scheduler_thread.start()
    
    def stop_scheduler(self, timeout: float = 30.0):
        """Stop the task scheduler"""
        self.running = False
        
        if self.scheduler_thread:
            self.scheduler_thread.join(timeout=timeout)
        
        # Cancel pending tasks
        with self.lock:
            while not self.pending_queue.empty():
                try:
                    task = self.pending_queue.get_nowait()
                    task.status = TaskStatus.CANCELLED
                    self.completed_tasks[task.task_id] = task
                except Empty:
                    break
        
        self.logger.info("Task scheduler stopped")
    
    def cancel_task(self, task_id: str) -> bool:
        """Cancel a pending or running task"""
        with self.lock:
            # Check if task is running
            if task_id in self.running_tasks:
                task = self.running_tasks[task_id]
                task.status = TaskStatus.CANCELLED
                # Resource cleanup will happen in task completion
                self.logger.info(f"Cancelled running task {task_id}")
                return True
            
            # Check pending queue (this is inefficient but necessary)
            temp_tasks = []
            cancelled = False
            
            while not self.pending_queue.empty():
                try:
                    task = self.pending_queue.get_nowait()
                    if task.task_id == task_id:
                        task.status = TaskStatus.CANCELLED
                        self.completed_tasks[task_id] = task
                        cancelled = True
                        self.logger.info(f"Cancelled pending task {task_id}")
                    else:
                        temp_tasks.append(task)
                except Empty:
                    break
            
            # Put remaining tasks back
            for task in temp_tasks:
                self.pending_queue.put(task)
            
            return cancelled
    
    def get_task_status(self, task_id: str) -> Optional[TaskStatus]:
        """Get status of a task"""
        with self.lock:
            if task_id in self.running_tasks:
                return self.running_tasks[task_id].status
            elif task_id in self.completed_tasks:
                return self.completed_tasks[task_id].status
            elif task_id in self.failed_tasks:
                return self.failed_tasks[task_id].status
        return None
    
    def get_queue_stats(self) -> Dict[str, Any]:
        """Get scheduler statistics"""
        with self.lock:
            total_requests = self.stats['tasks_submitted']
            avg_wait_time = (self.stats['total_wait_time'] / total_requests 
                           if total_requests > 0 else 0.0)
            avg_execution_time = (self.stats['total_execution_time'] / 
                                self.stats['tasks_completed'] 
                                if self.stats['tasks_completed'] > 0 else 0.0)
            
            return {
                **self.stats,
                'pending_tasks': self.pending_queue.qsize(),
                'running_tasks': len(self.running_tasks),
                'completed_tasks': len(self.completed_tasks),
                'failed_tasks': len(self.failed_tasks),
                'average_wait_time': avg_wait_time,
                'average_execution_time': avg_execution_time,
                'success_rate': (self.stats['tasks_completed'] / total_requests 
                               if total_requests > 0 else 0.0)
            }
    
    def _validate_task(self, task: ScheduledTask) -> bool:
        """Validate task before scheduling"""
        if not task.task_id:
            self.logger.error("Task ID is required")
            return False
        
        if task.task_id in self.running_tasks or task.task_id in self.completed_tasks:
            self.logger.error(f"Task {task.task_id} already exists")
            return False
        
        return True
    
    def _check_dependencies(self, task: ScheduledTask) -> bool:
        """Check if task dependencies are satisfied"""
        for dep_id in task.dependencies:
            if dep_id not in self.completed_tasks:
                return False
            
            dep_task = self.completed_tasks[dep_id]
            if dep_task.status != TaskStatus.COMPLETED:
                return False
        
        return True
    
    def _get_next_task(self) -> Optional[ScheduledTask]:
        """Get next task to execute based on scheduling policy"""
        try:
            if self.policy == SchedulingPolicy.FIFO:
                return self._get_fifo_task()
            elif self.policy == SchedulingPolicy.PRIORITY:
                return self._get_priority_task()
            elif self.policy == SchedulingPolicy.SHORTEST_JOB_FIRST:
                return self._get_sjf_task()
            elif self.policy == SchedulingPolicy.FAIR_SHARE:
                return self._get_fair_share_task()
            elif self.policy == SchedulingPolicy.RESOURCE_AWARE:
                return self._get_resource_aware_task()
            else:
                return self._get_priority_task()  # Default
                
        except Empty:
            return None
    
    def _get_fifo_task(self) -> Optional[ScheduledTask]:
        """Get next task using FIFO policy"""
        return self.pending_queue.get_nowait()
    
    def _get_priority_task(self) -> Optional[ScheduledTask]:
        """Get next task using priority policy"""
        return self.pending_queue.get_nowait()
    
    def _get_sjf_task(self) -> Optional[ScheduledTask]:
        """Get shortest job first"""
        # This is inefficient but demonstrates the concept
        tasks = []
        while not self.pending_queue.empty():
            try:
                tasks.append(self.pending_queue.get_nowait())
            except Empty:
                break
        
        if not tasks:
            return None
        
        # Sort by estimated runtime
        tasks.sort(key=lambda t: t.resources.estimated_runtime_minutes)
        
        # Put all but the first back
        for task in tasks[1:]:
            self.pending_queue.put(task)
        
        return tasks[0]
    
    def _get_fair_share_task(self) -> Optional[ScheduledTask]:
        """Get task using fair share policy"""
        # Simplified fair share - would need user identification in practice
        return self.pending_queue.get_nowait()
    
    def _get_resource_aware_task(self) -> Optional[ScheduledTask]:
        """Get task that best fits available resources"""
        tasks = []
        while not self.pending_queue.empty():
            try:
                tasks.append(self.pending_queue.get_nowait())
            except Empty:
                break
        
        if not tasks:
            return None
        
        # Find task that best fits available resources
        available = self.resource_manager.available_resources
        
        def resource_fit_score(task: ScheduledTask) -> float:
            req = task.resources
            if req <= available:
                # Higher score for better resource utilization
                return (req.cpu_cores / available.cpu_cores + 
                       req.memory_gb / available.memory_gb) / 2.0
            else:
                return -1.0  # Cannot fit
        
        # Sort by fit score (descending)
        tasks.sort(key=resource_fit_score, reverse=True)
        
        # Put all but the first back (if it can fit)
        selected_task = None
        for i, task in enumerate(tasks):
            if i == 0 and resource_fit_score(task) > 0:
                selected_task = task
            else:
                self.pending_queue.put(task)
        
        return selected_task
    
    def _start_task_execution(self, task: ScheduledTask):
        """Start executing a task"""
        with self.lock:
            task.status = TaskStatus.RUNNING
            task.started_time = time.time()
            self.running_tasks[task.task_id] = task
            
            # Update wait time statistics
            if task.queued_time:
                self.stats['total_wait_time'] += task.wait_time
        
        # Start execution in separate thread
        def execute_task():
            try:
                # Simulate task execution (in real system, this would call actual execution)
                execution_time = task.resources.estimated_runtime_minutes * 60
                time.sleep(min(execution_time, 10))  # Cap at 10 seconds for demo
                
                # Mark as completed
                self._complete_task(task, success=True)
                
            except Exception as e:
                self._complete_task(task, success=False, error=str(e))
        
        execution_thread = threading.Thread(target=execute_task, daemon=True)
        execution_thread.start()
    
    def _complete_task(self, task: ScheduledTask, success: bool, error: Optional[str] = None):
        """Complete a task and update statistics"""
        with self.lock:
            task.completed_time = time.time()
            
            if success:
                task.status = TaskStatus.COMPLETED
                self.completed_tasks[task.task_id] = task
                self.stats['tasks_completed'] += 1
            else:
                task.status = TaskStatus.FAILED
                task.error_message = error
                
                # Retry logic
                if task.retry_count < task.max_retries:
                    task.retry_count += 1
                    task.status = TaskStatus.QUEUED
                    task.queued_time = time.time()
                    self.pending_queue.put(task)
                    self.logger.info(f"Retrying task {task.task_id} (attempt {task.retry_count})")
                else:
                    self.failed_tasks[task.task_id] = task
                    self.stats['tasks_failed'] += 1
                    self.logger.error(f"Task {task.task_id} failed after {task.max_retries} attempts")
            
            # Update statistics
            if task.started_time:
                self.stats['total_execution_time'] += task.execution_time
            
            # Release resources
            self.resource_manager.release_resources(task.task_id)
            
            # Remove from running tasks
            if task.task_id in self.running_tasks:
                del self.running_tasks[task.task_id]
            
            self.logger.info(f"Task {task.task_id} completed with status {task.status}")


# Convenience functions
def create_task_scheduler(policy: SchedulingPolicy = SchedulingPolicy.PRIORITY) -> TaskScheduler:
    """Create a task scheduler with resource manager"""
    resource_manager = ResourceManager()
    return TaskScheduler(resource_manager, policy)


def create_pmf_task(task_id: str, system_path: str, 
                   cpu_cores: int = 4, memory_gb: float = 8.0,
                   priority: TaskPriority = TaskPriority.NORMAL) -> ScheduledTask:
    """Create a PMF calculation task"""
    resources = ResourceRequirements(
        cpu_cores=cpu_cores,
        memory_gb=memory_gb,
        estimated_runtime_minutes=120.0
    )
    
    return ScheduledTask(
        task_id=task_id,
        command="run_pmf_calculation",
        args=[system_path],
        resources=resources,
        priority=priority,
        timeout_minutes=180.0
    )