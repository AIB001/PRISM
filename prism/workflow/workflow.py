#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Workflow Definition System

Defines workflows, tasks, and their dependencies for molecular dynamics
simulations with support for complex workflow patterns and execution strategies.
"""

import time
import json
import uuid
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Callable, Set, Tuple
from dataclasses import dataclass, field, asdict
from enum import Enum
from abc import ABC, abstractmethod
from collections import defaultdict, deque

from ..utils.logging_system import PrismLogger
from .scheduler import ScheduledTask, TaskPriority, ResourceRequirements


class WorkflowStatus(Enum):
    """Workflow execution status"""
    CREATED = "created"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    PAUSED = "paused"


class TaskType(Enum):
    """Types of tasks in workflows"""
    SYSTEM_PREP = "system_preparation"
    MD_SIMULATION = "md_simulation"
    SMD_PULLING = "smd_pulling"
    UMBRELLA_SAMPLING = "umbrella_sampling"
    PMF_ANALYSIS = "pmf_analysis"
    DATA_PROCESSING = "data_processing"
    VALIDATION = "validation"
    CLEANUP = "cleanup"
    CUSTOM = "custom"


class DependencyType(Enum):
    """Types of task dependencies"""
    SEQUENTIAL = "sequential"        # Task B starts after Task A completes
    DATA_DEPENDENCY = "data"        # Task B needs data from Task A
    RESOURCE_DEPENDENCY = "resource" # Task B needs resources released by Task A
    CONDITIONAL = "conditional"      # Task B starts if Task A meets condition


@dataclass
class TaskDependency:
    """Represents a dependency between tasks"""
    upstream_task_id: str
    downstream_task_id: str
    dependency_type: DependencyType
    condition: Optional[Callable[[Any], bool]] = None
    data_key: Optional[str] = None
    
    def is_satisfied(self, upstream_result: Any = None) -> bool:
        """Check if dependency is satisfied"""
        if self.dependency_type == DependencyType.CONDITIONAL:
            return self.condition(upstream_result) if self.condition else True
        return True  # Other types are satisfied when upstream task completes


@dataclass
class Task:
    """Represents a single task in a workflow"""
    task_id: str
    task_type: TaskType
    command: str
    args: List[str] = field(default_factory=list)
    kwargs: Dict[str, Any] = field(default_factory=dict)
    
    # Resource requirements
    resources: ResourceRequirements = field(default_factory=ResourceRequirements)
    
    # Execution parameters
    priority: TaskPriority = TaskPriority.NORMAL
    max_retries: int = 3
    timeout_minutes: float = 120.0
    
    # Dependencies
    dependencies: Set[str] = field(default_factory=set)
    
    # Task metadata
    description: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    # Execution state
    status: str = "pending"
    result: Optional[Any] = None
    error_message: Optional[str] = None
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    
    def to_scheduled_task(self) -> ScheduledTask:
        """Convert to ScheduledTask for execution"""
        return ScheduledTask(
            task_id=self.task_id,
            command=self.command,
            args=self.args,
            kwargs=self.kwargs,
            priority=self.priority,
            resources=self.resources,
            dependencies=self.dependencies,
            max_retries=self.max_retries,
            timeout_minutes=self.timeout_minutes
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert task to dictionary"""
        data = asdict(self)
        data['task_type'] = self.task_type.value
        data['priority'] = self.priority.value
        data['resources'] = self.resources.to_dict()
        data['dependencies'] = list(self.dependencies)
        return data


class Workflow:
    """Represents a complete workflow with tasks and dependencies"""
    
    def __init__(self, workflow_id: str, name: str, description: str = ""):
        self.workflow_id = workflow_id
        self.name = name
        self.description = description
        
        # Workflow components
        self.tasks: Dict[str, Task] = {}
        self.dependencies: List[TaskDependency] = []
        
        # Execution state
        self.status = WorkflowStatus.CREATED
        self.start_time: Optional[float] = None
        self.end_time: Optional[float] = None
        self.current_phase: str = ""
        
        # Results and metadata
        self.results: Dict[str, Any] = {}
        self.metadata: Dict[str, Any] = {}
        
        self.logger = PrismLogger(f"workflow.{workflow_id}")
    
    def add_task(self, task: Task):
        """Add a task to the workflow"""
        if task.task_id in self.tasks:
            raise ValueError(f"Task {task.task_id} already exists in workflow")
        
        self.tasks[task.task_id] = task
        self.logger.info(f"Added task {task.task_id} ({task.task_type.value})")
    
    def add_dependency(self, dependency: TaskDependency):
        """Add a dependency between tasks"""
        # Validate that both tasks exist
        if dependency.upstream_task_id not in self.tasks:
            raise ValueError(f"Upstream task {dependency.upstream_task_id} not found")
        
        if dependency.downstream_task_id not in self.tasks:
            raise ValueError(f"Downstream task {dependency.downstream_task_id} not found")
        
        self.dependencies.append(dependency)
        
        # Update task dependencies
        downstream_task = self.tasks[dependency.downstream_task_id]
        downstream_task.dependencies.add(dependency.upstream_task_id)
        
        self.logger.info(f"Added dependency: {dependency.upstream_task_id} -> {dependency.downstream_task_id}")
    
    def validate(self) -> Tuple[bool, List[str]]:
        """Validate workflow for correctness"""
        errors = []
        
        # Check for cycles in dependencies
        if self._has_cycles():
            errors.append("Workflow contains dependency cycles")
        
        # Check that all task dependencies exist
        for task in self.tasks.values():
            for dep_id in task.dependencies:
                if dep_id not in self.tasks:
                    errors.append(f"Task {task.task_id} depends on non-existent task {dep_id}")
        
        # Check resource requirements are reasonable
        for task in self.tasks.values():
            if task.resources.cpu_cores <= 0:
                errors.append(f"Task {task.task_id} has invalid CPU requirement")
            
            if task.resources.memory_gb <= 0:
                errors.append(f"Task {task.task_id} has invalid memory requirement")
        
        return len(errors) == 0, errors
    
    def get_ready_tasks(self) -> List[Task]:
        """Get tasks that are ready to execute (all dependencies satisfied)"""
        ready_tasks = []
        
        for task in self.tasks.values():
            if task.status != "pending":
                continue
            
            # Check if all dependencies are completed
            deps_satisfied = True
            for dep_id in task.dependencies:
                dep_task = self.tasks.get(dep_id)
                if not dep_task or dep_task.status != "completed":
                    deps_satisfied = False
                    break
            
            if deps_satisfied:
                ready_tasks.append(task)
        
        return ready_tasks
    
    def get_execution_order(self) -> List[List[str]]:
        """Get topological ordering of tasks for execution"""
        # Build adjacency list
        graph = defaultdict(list)
        in_degree = defaultdict(int)
        
        # Initialize all tasks
        for task_id in self.tasks:
            in_degree[task_id] = 0
        
        # Build graph from dependencies
        for dep in self.dependencies:
            graph[dep.upstream_task_id].append(dep.downstream_task_id)
            in_degree[dep.downstream_task_id] += 1
        
        # Topological sort with level grouping
        levels = []
        queue = deque([task_id for task_id in self.tasks if in_degree[task_id] == 0])
        
        while queue:
            # Process all tasks at current level
            current_level = []
            level_size = len(queue)
            
            for _ in range(level_size):
                task_id = queue.popleft()
                current_level.append(task_id)
                
                # Update neighbors
                for neighbor in graph[task_id]:
                    in_degree[neighbor] -= 1
                    if in_degree[neighbor] == 0:
                        queue.append(neighbor)
            
            if current_level:
                levels.append(current_level)
        
        return levels
    
    def estimate_execution_time(self) -> float:
        """Estimate total workflow execution time"""
        execution_order = self.get_execution_order()
        
        total_time = 0.0
        for level in execution_order:
            # For parallel execution, take the maximum time in each level
            level_time = max(
                self.tasks[task_id].resources.estimated_runtime_minutes 
                for task_id in level
            )
            total_time += level_time
        
        return total_time
    
    def get_resource_requirements(self) -> ResourceRequirements:
        """Get peak resource requirements for workflow"""
        execution_order = self.get_execution_order()
        
        peak_resources = ResourceRequirements()
        
        for level in execution_order:
            level_resources = ResourceRequirements()
            for task_id in level:
                task_resources = self.tasks[task_id].resources
                level_resources = level_resources + task_resources
            
            # Update peak if this level requires more resources
            if (level_resources.cpu_cores > peak_resources.cpu_cores or
                level_resources.memory_gb > peak_resources.memory_gb):
                peak_resources = level_resources
        
        return peak_resources
    
    def _has_cycles(self) -> bool:
        """Check for cycles using DFS"""
        WHITE, GRAY, BLACK = 0, 1, 2
        color = {task_id: WHITE for task_id in self.tasks}
        
        # Build adjacency list
        graph = defaultdict(list)
        for dep in self.dependencies:
            graph[dep.upstream_task_id].append(dep.downstream_task_id)
        
        def dfs(node):
            if color[node] == GRAY:
                return True  # Back edge found, cycle detected
            
            if color[node] == BLACK:
                return False
            
            color[node] = GRAY
            
            for neighbor in graph[node]:
                if dfs(neighbor):
                    return True
            
            color[node] = BLACK
            return False
        
        for task_id in self.tasks:
            if color[task_id] == WHITE:
                if dfs(task_id):
                    return True
        
        return False
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert workflow to dictionary"""
        return {
            'workflow_id': self.workflow_id,
            'name': self.name,
            'description': self.description,
            'status': self.status.value,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'current_phase': self.current_phase,
            'tasks': {task_id: task.to_dict() for task_id, task in self.tasks.items()},
            'dependencies': [
                {
                    'upstream': dep.upstream_task_id,
                    'downstream': dep.downstream_task_id,
                    'type': dep.dependency_type.value,
                    'data_key': dep.data_key
                }
                for dep in self.dependencies
            ],
            'results': self.results,
            'metadata': self.metadata,
            'estimated_execution_time': self.estimate_execution_time(),
            'peak_resources': self.get_resource_requirements().to_dict()
        }


class WorkflowBuilder:
    """Builder pattern for creating workflows"""
    
    def __init__(self, workflow_id: Optional[str] = None, name: str = "Unnamed Workflow"):
        if workflow_id is None:
            workflow_id = str(uuid.uuid4())
        
        self.workflow = Workflow(workflow_id, name)
        self.logger = PrismLogger(f"workflow_builder.{workflow_id}")
    
    def set_description(self, description: str) -> 'WorkflowBuilder':
        """Set workflow description"""
        self.workflow.description = description
        return self
    
    def add_task(self, task_id: str, task_type: TaskType, command: str,
                args: Optional[List[str]] = None,
                kwargs: Optional[Dict[str, Any]] = None,
                resources: Optional[ResourceRequirements] = None,
                priority: TaskPriority = TaskPriority.NORMAL,
                description: str = "") -> 'WorkflowBuilder':
        """Add a task to the workflow"""
        
        task = Task(
            task_id=task_id,
            task_type=task_type,
            command=command,
            args=args or [],
            kwargs=kwargs or {},
            resources=resources or ResourceRequirements(),
            priority=priority,
            description=description
        )
        
        self.workflow.add_task(task)
        return self
    
    def add_dependency(self, upstream_task_id: str, downstream_task_id: str,
                      dependency_type: DependencyType = DependencyType.SEQUENTIAL,
                      condition: Optional[Callable] = None,
                      data_key: Optional[str] = None) -> 'WorkflowBuilder':
        """Add a dependency between tasks"""
        
        dependency = TaskDependency(
            upstream_task_id=upstream_task_id,
            downstream_task_id=downstream_task_id,
            dependency_type=dependency_type,
            condition=condition,
            data_key=data_key
        )
        
        self.workflow.add_dependency(dependency)
        return self
    
    def add_sequential_tasks(self, task_configs: List[Dict[str, Any]]) -> 'WorkflowBuilder':
        """Add a sequence of dependent tasks"""
        
        previous_task_id = None
        
        for config in task_configs:
            task_id = config['task_id']
            
            # Add the task
            self.add_task(
                task_id=task_id,
                task_type=TaskType(config.get('task_type', 'custom')),
                command=config['command'],
                args=config.get('args', []),
                kwargs=config.get('kwargs', {}),
                resources=ResourceRequirements(**config.get('resources', {})),
                priority=TaskPriority(config.get('priority', TaskPriority.NORMAL.value)),
                description=config.get('description', '')
            )
            
            # Add dependency to previous task
            if previous_task_id:
                self.add_dependency(previous_task_id, task_id)
            
            previous_task_id = task_id
        
        return self
    
    def add_parallel_tasks(self, task_configs: List[Dict[str, Any]], 
                          upstream_task_id: Optional[str] = None,
                          downstream_task_id: Optional[str] = None) -> 'WorkflowBuilder':
        """Add parallel tasks with optional upstream/downstream dependencies"""
        
        parallel_task_ids = []
        
        # Add all parallel tasks
        for config in task_configs:
            task_id = config['task_id']
            
            self.add_task(
                task_id=task_id,
                task_type=TaskType(config.get('task_type', 'custom')),
                command=config['command'],
                args=config.get('args', []),
                kwargs=config.get('kwargs', {}),
                resources=ResourceRequirements(**config.get('resources', {})),
                priority=TaskPriority(config.get('priority', TaskPriority.NORMAL.value)),
                description=config.get('description', '')
            )
            
            parallel_task_ids.append(task_id)
            
            # Add dependency from upstream task
            if upstream_task_id:
                self.add_dependency(upstream_task_id, task_id)
        
        # Add dependencies to downstream task
        if downstream_task_id:
            for task_id in parallel_task_ids:
                self.add_dependency(task_id, downstream_task_id)
        
        return self
    
    def build(self) -> Workflow:
        """Build and validate the workflow"""
        is_valid, errors = self.workflow.validate()
        
        if not is_valid:
            error_msg = "Workflow validation failed: " + "; ".join(errors)
            self.logger.error(error_msg)
            raise ValueError(error_msg)
        
        self.logger.info(f"Built workflow {self.workflow.workflow_id} with {len(self.workflow.tasks)} tasks")
        return self.workflow


# Predefined workflow templates
class WorkflowTemplates:
    """Collection of predefined workflow templates"""
    
    @staticmethod
    def create_pmf_workflow(workflow_id: str, system_path: str,
                           ligand_residue: str = "LIG",
                           protein_residue: str = "Protein") -> Workflow:
        """Create a standard PMF calculation workflow"""
        
        builder = WorkflowBuilder(workflow_id, "PMF Calculation Workflow")
        builder.set_description(f"Complete PMF calculation for {system_path}")
        
        # Define resource requirements for each step
        prep_resources = ResourceRequirements(cpu_cores=2, memory_gb=4.0, estimated_runtime_minutes=30)
        equilibration_resources = ResourceRequirements(cpu_cores=4, memory_gb=8.0, estimated_runtime_minutes=60)
        smd_resources = ResourceRequirements(cpu_cores=8, memory_gb=16.0, estimated_runtime_minutes=120)
        umbrella_resources = ResourceRequirements(cpu_cores=16, memory_gb=32.0, estimated_runtime_minutes=240)
        analysis_resources = ResourceRequirements(cpu_cores=4, memory_gb=8.0, estimated_runtime_minutes=30)
        
        # Build workflow
        builder.add_task(
            task_id="system_preparation",
            task_type=TaskType.SYSTEM_PREP,
            command="prepare_pmf_system",
            args=[system_path, ligand_residue, protein_residue],
            resources=prep_resources,
            priority=TaskPriority.HIGH,
            description="Prepare system for PMF calculation"
        ).add_task(
            task_id="equilibration",
            task_type=TaskType.MD_SIMULATION,
            command="run_equilibration",
            args=[],
            resources=equilibration_resources,
            description="Equilibrate the system"
        ).add_task(
            task_id="smd_pulling",
            task_type=TaskType.SMD_PULLING,
            command="run_smd_pulling",
            args=[],
            resources=smd_resources,
            description="Run steered MD pulling simulation"
        ).add_task(
            task_id="umbrella_sampling",
            task_type=TaskType.UMBRELLA_SAMPLING,
            command="run_umbrella_sampling",
            args=[],
            resources=umbrella_resources,
            description="Run umbrella sampling simulations"
        ).add_task(
            task_id="pmf_analysis",
            task_type=TaskType.PMF_ANALYSIS,
            command="analyze_pmf",
            args=[],
            resources=analysis_resources,
            priority=TaskPriority.HIGH,
            description="Calculate PMF from umbrella sampling"
        ).add_task(
            task_id="validation",
            task_type=TaskType.VALIDATION,
            command="validate_pmf_results",
            args=[],
            resources=ResourceRequirements(cpu_cores=2, memory_gb=4.0, estimated_runtime_minutes=15),
            description="Validate PMF results"
        )
        
        # Add dependencies
        builder.add_dependency("system_preparation", "equilibration")
        builder.add_dependency("equilibration", "smd_pulling")
        builder.add_dependency("smd_pulling", "umbrella_sampling")
        builder.add_dependency("umbrella_sampling", "pmf_analysis")
        builder.add_dependency("pmf_analysis", "validation")
        
        return builder.build()
    
    @staticmethod
    def create_high_throughput_workflow(workflow_id: str, 
                                      system_paths: List[str]) -> Workflow:
        """Create a high-throughput screening workflow"""
        
        builder = WorkflowBuilder(workflow_id, "High-Throughput Screening")
        builder.set_description(f"Process {len(system_paths)} systems in parallel")
        
        # Add preparation phase
        builder.add_task(
            task_id="global_preparation",
            task_type=TaskType.SYSTEM_PREP,
            command="prepare_batch_systems",
            args=system_paths,
            resources=ResourceRequirements(cpu_cores=4, memory_gb=8.0, estimated_runtime_minutes=20),
            priority=TaskPriority.HIGH,
            description="Prepare all systems for processing"
        )
        
        # Add parallel processing tasks
        parallel_tasks = []
        for i, system_path in enumerate(system_paths):
            task_config = {
                'task_id': f"process_system_{i:03d}",
                'task_type': TaskType.MD_SIMULATION.value,
                'command': 'process_single_system',
                'args': [system_path],
                'resources': {
                    'cpu_cores': 4,
                    'memory_gb': 8.0,
                    'estimated_runtime_minutes': 90
                },
                'description': f"Process system {system_path}"
            }
            parallel_tasks.append(task_config)
        
        builder.add_parallel_tasks(
            parallel_tasks,
            upstream_task_id="global_preparation",
            downstream_task_id="final_analysis"
        )
        
        # Add final analysis
        builder.add_task(
            task_id="final_analysis",
            task_type=TaskType.DATA_PROCESSING,
            command="analyze_batch_results",
            args=[],
            resources=ResourceRequirements(cpu_cores=8, memory_gb=16.0, estimated_runtime_minutes=45),
            priority=TaskPriority.HIGH,
            description="Analyze results from all systems"
        )
        
        return builder.build()


# Convenience functions
def create_pmf_workflow(workflow_id: str, system_path: str) -> Workflow:
    """Create a standard PMF workflow"""
    return WorkflowTemplates.create_pmf_workflow(workflow_id, system_path)


def create_batch_workflow(workflow_id: str, system_paths: List[str]) -> Workflow:
    """Create a batch processing workflow"""
    return WorkflowTemplates.create_high_throughput_workflow(workflow_id, system_paths)


def load_workflow_from_json(json_path: Union[str, Path]) -> Workflow:
    """Load workflow from JSON file"""
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Reconstruct workflow from dictionary
    workflow = Workflow(
        data['workflow_id'],
        data['name'],
        data['description']
    )
    
    # Add tasks
    for task_id, task_data in data['tasks'].items():
        task = Task(
            task_id=task_id,
            task_type=TaskType(task_data['task_type']),
            command=task_data['command'],
            args=task_data['args'],
            kwargs=task_data['kwargs'],
            resources=ResourceRequirements(**task_data['resources']),
            priority=TaskPriority(task_data['priority']),
            max_retries=task_data.get('max_retries', 3),
            timeout_minutes=task_data.get('timeout_minutes', 120),
            description=task_data.get('description', '')
        )
        workflow.add_task(task)
    
    # Add dependencies
    for dep_data in data['dependencies']:
        dependency = TaskDependency(
            upstream_task_id=dep_data['upstream'],
            downstream_task_id=dep_data['downstream'],
            dependency_type=DependencyType(dep_data['type']),
            data_key=dep_data.get('data_key')
        )
        workflow.add_dependency(dependency)
    
    return workflow


def save_workflow_to_json(workflow: Workflow, json_path: Union[str, Path]):
    """Save workflow to JSON file"""
    with open(json_path, 'w') as f:
        json.dump(workflow.to_dict(), f, indent=2)