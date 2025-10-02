#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Workflow Executor

Executes workflows by managing task scheduling, execution context,
and result handling with support for distributed execution.
"""

import time
import threading
import concurrent.futures
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Callable, Tuple
from dataclasses import dataclass, field
from enum import Enum
from collections import defaultdict

from ..utils.logging_system import PrismLogger
from .scheduler import TaskScheduler, ScheduledTask, TaskStatus, SchedulingPolicy
from .workflow import Workflow, Task, WorkflowStatus


class ExecutionStrategy(Enum):
    """Workflow execution strategies"""
    SEQUENTIAL = "sequential"      # Execute tasks one by one
    PARALLEL = "parallel"         # Execute independent tasks in parallel
    ADAPTIVE = "adaptive"         # Dynamically adjust based on resources
    DISTRIBUTED = "distributed"   # Execute across multiple nodes


@dataclass
class ExecutionContext:
    """Context for workflow execution"""
    workflow_id: str
    execution_id: str
    working_directory: Path
    environment_variables: Dict[str, str] = field(default_factory=dict)
    execution_strategy: ExecutionStrategy = ExecutionStrategy.PARALLEL
    max_concurrent_tasks: int = 4
    checkpoint_enabled: bool = True
    checkpoint_interval: float = 300.0  # 5 minutes
    
    # Execution state
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    current_phase: str = ""
    
    # Results and outputs
    task_results: Dict[str, Any] = field(default_factory=dict)
    execution_log: List[str] = field(default_factory=list)
    
    def log_event(self, message: str):
        """Log an execution event"""
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        log_entry = f"[{timestamp}] {message}"
        self.execution_log.append(log_entry)


@dataclass
class ExecutionResult:
    """Result of workflow execution"""
    workflow_id: str
    execution_id: str
    status: WorkflowStatus
    start_time: float
    end_time: Optional[float] = None
    
    # Task results
    completed_tasks: List[str] = field(default_factory=list)
    failed_tasks: List[str] = field(default_factory=list)
    skipped_tasks: List[str] = field(default_factory=list)
    
    # Execution metrics
    total_execution_time: float = 0.0
    total_cpu_hours: float = 0.0
    peak_memory_usage: float = 0.0
    
    # Results and errors
    results: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    
    @property
    def success_rate(self) -> float:
        """Calculate task success rate"""
        total_tasks = len(self.completed_tasks) + len(self.failed_tasks) + len(self.skipped_tasks)
        if total_tasks == 0:
            return 0.0
        return len(self.completed_tasks) / total_tasks
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'workflow_id': self.workflow_id,
            'execution_id': self.execution_id,
            'status': self.status.value,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'total_execution_time': self.total_execution_time,
            'completed_tasks': len(self.completed_tasks),
            'failed_tasks': len(self.failed_tasks),
            'skipped_tasks': len(self.skipped_tasks),
            'success_rate': self.success_rate,
            'total_cpu_hours': self.total_cpu_hours,
            'peak_memory_usage': self.peak_memory_usage,
            'results': self.results,
            'errors': self.errors,
            'warnings': self.warnings
        }


class TaskExecutor:
    """Executes individual tasks with proper context management"""
    
    def __init__(self):
        self.logger = PrismLogger("task_executor")
    
    def execute_task(self, task: Task, context: ExecutionContext) -> Tuple[bool, Any, Optional[str]]:
        """Execute a single task"""
        self.logger.info(f"Executing task {task.task_id}: {task.command}")
        
        try:
            # Set up execution environment
            self._setup_task_environment(task, context)
            
            # Execute the task based on command
            result = self._execute_command(task, context)
            
            # Store result in context
            context.task_results[task.task_id] = result
            context.log_event(f"Task {task.task_id} completed successfully")
            
            return True, result, None
            
        except Exception as e:
            error_msg = f"Task {task.task_id} failed: {str(e)}"
            self.logger.error(error_msg)
            context.log_event(error_msg)
            return False, None, error_msg
    
    def _setup_task_environment(self, task: Task, context: ExecutionContext):
        """Set up environment for task execution"""
        # Set environment variables
        import os
        for key, value in context.environment_variables.items():
            os.environ[key] = value
        
        # Change to working directory
        os.chdir(context.working_directory)
        
        # Set task-specific environment
        os.environ['PRISM_TASK_ID'] = task.task_id
        os.environ['PRISM_WORKFLOW_ID'] = context.workflow_id
    
    def _execute_command(self, task: Task, context: ExecutionContext) -> Any:
        """Execute the task command"""
        command = task.command
        
        # This is where you would implement actual command execution
        # For now, we'll simulate different types of tasks
        
        if command == "prepare_pmf_system":
            return self._execute_system_preparation(task, context)
        elif command == "run_equilibration":
            return self._execute_equilibration(task, context)
        elif command == "run_smd_pulling":
            return self._execute_smd_pulling(task, context)
        elif command == "run_umbrella_sampling":
            return self._execute_umbrella_sampling(task, context)
        elif command == "analyze_pmf":
            return self._execute_pmf_analysis(task, context)
        elif command == "validate_pmf_results":
            return self._execute_validation(task, context)
        else:
            # Generic command execution
            return self._execute_generic_command(task, context)
    
    def _execute_system_preparation(self, task: Task, context: ExecutionContext) -> Dict[str, Any]:
        """Execute system preparation task"""
        # Simulate system preparation
        time.sleep(0.5)  # Simulate work
        
        return {
            'status': 'completed',
            'output_files': ['prepared_system.gro', 'system.top'],
            'system_size': 50000,
            'preparation_time': 0.5
        }
    
    def _execute_equilibration(self, task: Task, context: ExecutionContext) -> Dict[str, Any]:
        """Execute equilibration task"""
        time.sleep(1.0)  # Simulate equilibration
        
        return {
            'status': 'completed',
            'final_energy': -450000.0,
            'temperature': 300.0,
            'pressure': 1.0,
            'equilibration_time': 1.0
        }
    
    def _execute_smd_pulling(self, task: Task, context: ExecutionContext) -> Dict[str, Any]:
        """Execute SMD pulling task"""
        time.sleep(2.0)  # Simulate SMD
        
        return {
            'status': 'completed',
            'max_force': 1500.0,
            'pull_distance': 3.5,
            'suggested_windows': 30,
            'trajectory_file': 'smd_pull.xtc'
        }
    
    def _execute_umbrella_sampling(self, task: Task, context: ExecutionContext) -> Dict[str, Any]:
        """Execute umbrella sampling task"""
        time.sleep(3.0)  # Simulate umbrella sampling
        
        # Get SMD results for window setup
        smd_result = context.task_results.get('smd_pulling', {})
        num_windows = smd_result.get('suggested_windows', 30)
        
        return {
            'status': 'completed',
            'windows_completed': num_windows,
            'total_simulation_time': num_windows * 10.0,  # ns
            'convergence_achieved': True
        }
    
    def _execute_pmf_analysis(self, task: Task, context: ExecutionContext) -> Dict[str, Any]:
        """Execute PMF analysis task"""
        time.sleep(1.5)  # Simulate analysis
        
        return {
            'status': 'completed',
            'binding_energy': -8.5,
            'pmf_profile': list(range(30)),  # Simplified
            'error_estimate': 0.5,
            'convergence': True,
            'method': 'WHAM'
        }
    
    def _execute_validation(self, task: Task, context: ExecutionContext) -> Dict[str, Any]:
        """Execute validation task"""
        time.sleep(0.3)  # Simulate validation
        
        # Check PMF results
        pmf_result = context.task_results.get('pmf_analysis', {})
        binding_energy = pmf_result.get('binding_energy', 0.0)
        
        validation_passed = abs(binding_energy) > 1.0  # Simple validation
        
        return {
            'status': 'completed' if validation_passed else 'warning',
            'validation_passed': validation_passed,
            'checks_performed': ['energy_range', 'convergence', 'profile_smoothness'],
            'warnings': [] if validation_passed else ['Binding energy suspiciously low']
        }
    
    def _execute_generic_command(self, task: Task, context: ExecutionContext) -> Dict[str, Any]:
        """Execute generic command"""
        # Simulate generic execution
        execution_time = task.resources.estimated_runtime_minutes / 60.0
        time.sleep(min(execution_time, 2.0))  # Cap simulation time
        
        return {
            'status': 'completed',
            'execution_time': execution_time,
            'command': task.command,
            'args': task.args
        }


class WorkflowExecutor:
    """Main workflow executor with scheduling and monitoring"""
    
    def __init__(self, scheduler: Optional[TaskScheduler] = None,
                 max_concurrent_workflows: int = 2):
        
        self.scheduler = scheduler or TaskScheduler.create_default()
        self.task_executor = TaskExecutor()
        self.max_concurrent_workflows = max_concurrent_workflows
        
        self.logger = PrismLogger("workflow_executor")
        
        # Execution tracking
        self.active_executions: Dict[str, ExecutionContext] = {}
        self.completed_executions: Dict[str, ExecutionResult] = {}
        
        # Threading
        self.executor_pool = concurrent.futures.ThreadPoolExecutor(
            max_workers=max_concurrent_workflows
        )
        self.lock = threading.RLock()
    
    def execute_workflow(self, workflow: Workflow, 
                        execution_context: Optional[ExecutionContext] = None,
                        callback: Optional[Callable[[ExecutionResult], None]] = None) -> str:
        """Execute a workflow asynchronously"""
        
        # Generate execution ID
        execution_id = f"{workflow.workflow_id}_{int(time.time())}"
        
        # Create execution context if not provided
        if execution_context is None:
            execution_context = ExecutionContext(
                workflow_id=workflow.workflow_id,
                execution_id=execution_id,
                working_directory=Path.cwd(),
                execution_strategy=ExecutionStrategy.PARALLEL
            )
        
        execution_context.start_time = time.time()
        
        with self.lock:
            self.active_executions[execution_id] = execution_context
        
        # Submit workflow for execution
        future = self.executor_pool.submit(self._execute_workflow_sync, workflow, execution_context)
        
        if callback:
            future.add_done_callback(lambda f: callback(f.result()))
        
        self.logger.info(f"Started workflow execution {execution_id}")
        return execution_id
    
    def _execute_workflow_sync(self, workflow: Workflow, context: ExecutionContext) -> ExecutionResult:
        """Execute workflow synchronously"""
        
        self.logger.info(f"Executing workflow {workflow.workflow_id}")
        
        # Initialize execution result
        result = ExecutionResult(
            workflow_id=workflow.workflow_id,
            execution_id=context.execution_id,
            status=WorkflowStatus.RUNNING,
            start_time=context.start_time
        )
        
        try:
            workflow.status = WorkflowStatus.RUNNING
            workflow.start_time = context.start_time
            
            context.log_event(f"Started workflow execution")
            
            # Execute based on strategy
            if context.execution_strategy == ExecutionStrategy.SEQUENTIAL:
                self._execute_sequential(workflow, context, result)
            elif context.execution_strategy == ExecutionStrategy.PARALLEL:
                self._execute_parallel(workflow, context, result)
            else:
                self._execute_adaptive(workflow, context, result)
            
            # Determine final status
            if result.failed_tasks:
                result.status = WorkflowStatus.FAILED
                workflow.status = WorkflowStatus.FAILED
            else:
                result.status = WorkflowStatus.COMPLETED
                workflow.status = WorkflowStatus.COMPLETED
            
            context.log_event(f"Workflow execution completed with status: {result.status.value}")
            
        except Exception as e:
            error_msg = f"Workflow execution failed: {str(e)}"
            self.logger.error(error_msg)
            result.status = WorkflowStatus.FAILED
            result.errors.append(error_msg)
            workflow.status = WorkflowStatus.FAILED
            context.log_event(error_msg)
        
        finally:
            # Finalize execution
            end_time = time.time()
            result.end_time = end_time
            result.total_execution_time = end_time - result.start_time
            workflow.end_time = end_time
            context.end_time = end_time
            
            # Move to completed executions
            with self.lock:
                if context.execution_id in self.active_executions:
                    del self.active_executions[context.execution_id]
                self.completed_executions[context.execution_id] = result
        
        self.logger.info(f"Workflow {workflow.workflow_id} execution completed: {result.status.value}")
        return result
    
    def _execute_sequential(self, workflow: Workflow, context: ExecutionContext, result: ExecutionResult):
        """Execute workflow sequentially"""
        execution_order = workflow.get_execution_order()
        
        for level in execution_order:
            for task_id in level:
                task = workflow.tasks[task_id]
                self._execute_single_task(task, context, result)
                
                if task_id in result.failed_tasks:
                    self.logger.error(f"Task {task_id} failed, stopping sequential execution")
                    return
    
    def _execute_parallel(self, workflow: Workflow, context: ExecutionContext, result: ExecutionResult):
        """Execute workflow with parallel task execution"""
        execution_order = workflow.get_execution_order()
        
        for level in execution_order:
            if len(level) == 1:
                # Single task, execute directly
                task = workflow.tasks[level[0]]
                self._execute_single_task(task, context, result)
            else:
                # Multiple tasks, execute in parallel
                self._execute_parallel_level(level, workflow, context, result)
            
            # Check for failures at each level
            failed_in_level = [task_id for task_id in level if task_id in result.failed_tasks]
            if failed_in_level:
                self.logger.warning(f"Tasks failed in level: {failed_in_level}")
                # Continue with remaining levels (could be made configurable)
    
    def _execute_adaptive(self, workflow: Workflow, context: ExecutionContext, result: ExecutionResult):
        """Execute workflow with adaptive strategy"""
        # For now, use parallel execution as adaptive
        self._execute_parallel(workflow, context, result)
    
    def _execute_parallel_level(self, task_ids: List[str], workflow: Workflow, 
                               context: ExecutionContext, result: ExecutionResult):
        """Execute a level of tasks in parallel"""
        
        max_workers = min(len(task_ids), context.max_concurrent_tasks)
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks in this level
            future_to_task = {
                executor.submit(self._execute_single_task, workflow.tasks[task_id], context, result): task_id
                for task_id in task_ids
            }
            
            # Wait for completion
            for future in concurrent.futures.as_completed(future_to_task):
                task_id = future_to_task[future]
                try:
                    future.result()  # This will raise exception if task failed
                except Exception as e:
                    self.logger.error(f"Parallel task {task_id} failed: {e}")
    
    def _execute_single_task(self, task: Task, context: ExecutionContext, result: ExecutionResult):
        """Execute a single task and update results"""
        task_start_time = time.time()
        
        try:
            # Execute the task
            success, task_result, error_msg = self.task_executor.execute_task(task, context)
            
            task_end_time = time.time()
            task_duration = task_end_time - task_start_time
            
            if success:
                result.completed_tasks.append(task.task_id)
                if task_result:
                    result.results[task.task_id] = task_result
                
                # Update task status
                task.status = "completed"
                task.result = task_result
                task.start_time = task_start_time
                task.end_time = task_end_time
                
                # Update execution metrics
                cpu_hours = task.resources.cpu_cores * (task_duration / 3600.0)
                result.total_cpu_hours += cpu_hours
                
                # Estimate peak memory (simplified)
                if task.resources.memory_gb > result.peak_memory_usage:
                    result.peak_memory_usage = task.resources.memory_gb
                
            else:
                result.failed_tasks.append(task.task_id)
                if error_msg:
                    result.errors.append(f"Task {task.task_id}: {error_msg}")
                
                task.status = "failed"
                task.error_message = error_msg
                task.start_time = task_start_time
                task.end_time = task_end_time
            
        except Exception as e:
            error_msg = f"Unexpected error in task {task.task_id}: {str(e)}"
            self.logger.error(error_msg)
            result.failed_tasks.append(task.task_id)
            result.errors.append(error_msg)
            task.status = "failed"
            task.error_message = error_msg
    
    def get_execution_status(self, execution_id: str) -> Optional[Dict[str, Any]]:
        """Get status of a workflow execution"""
        
        # Check active executions
        if execution_id in self.active_executions:
            context = self.active_executions[execution_id]
            return {
                'execution_id': execution_id,
                'status': 'running',
                'start_time': context.start_time,
                'current_phase': context.current_phase,
                'completed_tasks': len([t for t in context.task_results.keys()]),
                'total_tasks': len(context.task_results),
                'execution_time': time.time() - (context.start_time or time.time())
            }
        
        # Check completed executions
        if execution_id in self.completed_executions:
            result = self.completed_executions[execution_id]
            return result.to_dict()
        
        return None
    
    def cancel_execution(self, execution_id: str) -> bool:
        """Cancel a workflow execution"""
        if execution_id in self.active_executions:
            context = self.active_executions[execution_id]
            context.log_event("Execution cancelled by user")
            # In a real implementation, you would signal the execution thread to stop
            self.logger.info(f"Cancelled execution {execution_id}")
            return True
        return False
    
    def get_execution_statistics(self) -> Dict[str, Any]:
        """Get executor statistics"""
        with self.lock:
            total_executions = len(self.active_executions) + len(self.completed_executions)
            
            completed_results = list(self.completed_executions.values())
            if completed_results:
                avg_execution_time = sum(r.total_execution_time for r in completed_results) / len(completed_results)
                success_rate = sum(1 for r in completed_results if r.status == WorkflowStatus.COMPLETED) / len(completed_results)
                total_cpu_hours = sum(r.total_cpu_hours for r in completed_results)
            else:
                avg_execution_time = 0.0
                success_rate = 0.0
                total_cpu_hours = 0.0
            
            return {
                'total_executions': total_executions,
                'active_executions': len(self.active_executions),
                'completed_executions': len(self.completed_executions),
                'average_execution_time': avg_execution_time,
                'success_rate': success_rate,
                'total_cpu_hours': total_cpu_hours,
                'max_concurrent_workflows': self.max_concurrent_workflows
            }
    
    def shutdown(self, timeout: float = 30.0):
        """Shutdown the workflow executor"""
        self.logger.info("Shutting down workflow executor")
        self.executor_pool.shutdown(wait=True, timeout=timeout)


# Convenience functions
def create_workflow_executor(max_concurrent_workflows: int = 2) -> WorkflowExecutor:
    """Create a workflow executor with default configuration"""
    from .scheduler import create_task_scheduler
    
    scheduler = create_task_scheduler()
    scheduler.start_scheduler()
    
    return WorkflowExecutor(scheduler, max_concurrent_workflows)


def execute_pmf_workflow(system_path: str, working_dir: Optional[Path] = None) -> str:
    """Execute a standard PMF workflow"""
    from .workflow import WorkflowTemplates
    
    # Create workflow
    workflow_id = f"pmf_{int(time.time())}"
    workflow = WorkflowTemplates.create_pmf_workflow(workflow_id, system_path)
    
    # Create executor
    executor = create_workflow_executor()
    
    # Create execution context
    context = ExecutionContext(
        workflow_id=workflow_id,
        execution_id=f"{workflow_id}_exec",
        working_directory=working_dir or Path.cwd()
    )
    
    # Execute workflow
    return executor.execute_workflow(workflow, context)