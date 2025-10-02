#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Workflow Monitoring System

Real-time monitoring and metrics collection for workflows, tasks, and resources
with alerting, performance tracking, and comprehensive reporting capabilities.
"""

import time
import json
import threading
import statistics
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Callable, Tuple
from dataclasses import dataclass, field
from enum import Enum
from abc import ABC, abstractmethod
from collections import defaultdict, deque
import heapq

from ..utils.logging_system import PrismLogger
from .scheduler import TaskScheduler, ScheduledTask, TaskStatus, ResourceManager
from .workflow import Workflow, WorkflowStatus, Task
from .executor import WorkflowExecutor, ExecutionResult
from .orchestrator import WorkflowOrchestrator, ComputeNode


class MetricType(Enum):
    """Types of metrics"""
    COUNTER = "counter"
    GAUGE = "gauge"
    HISTOGRAM = "histogram"
    SUMMARY = "summary"


class AlertSeverity(Enum):
    """Alert severity levels"""
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


@dataclass
class Metric:
    """Represents a single metric measurement"""
    name: str
    value: float
    timestamp: float
    labels: Dict[str, str] = field(default_factory=dict)
    metric_type: MetricType = MetricType.GAUGE
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'name': self.name,
            'value': self.value,
            'timestamp': self.timestamp,
            'labels': self.labels,
            'type': self.metric_type.value
        }


@dataclass
class Alert:
    """Represents an alert condition"""
    alert_id: str
    title: str
    description: str
    severity: AlertSeverity
    timestamp: float
    source: str
    labels: Dict[str, str] = field(default_factory=dict)
    resolved: bool = False
    resolved_timestamp: Optional[float] = None
    
    def resolve(self):
        """Mark alert as resolved"""
        self.resolved = True
        self.resolved_timestamp = time.time()
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'alert_id': self.alert_id,
            'title': self.title,
            'description': self.description,
            'severity': self.severity.value,
            'timestamp': self.timestamp,
            'source': self.source,
            'labels': self.labels,
            'resolved': self.resolved,
            'resolved_timestamp': self.resolved_timestamp
        }


class ExecutionMetrics:
    """Collects and manages execution metrics"""
    
    def __init__(self):
        self.metrics: List[Metric] = []
        self.metric_history: Dict[str, deque] = defaultdict(lambda: deque(maxlen=1000))
        self.lock = threading.RLock()
        
        # Performance counters
        self.counters = {
            'workflows_submitted': 0,
            'workflows_completed': 0,
            'workflows_failed': 0,
            'tasks_executed': 0,
            'tasks_failed': 0,
            'alerts_triggered': 0
        }
        
        # Time series data
        self.time_series = {
            'workflow_execution_times': deque(maxlen=100),
            'task_execution_times': deque(maxlen=1000),
            'resource_utilization': deque(maxlen=1000),
            'queue_lengths': deque(maxlen=1000)
        }
    
    def record_metric(self, name: str, value: float, 
                     labels: Optional[Dict[str, str]] = None,
                     metric_type: MetricType = MetricType.GAUGE):
        """Record a metric measurement"""
        with self.lock:
            metric = Metric(
                name=name,
                value=value,
                timestamp=time.time(),
                labels=labels or {},
                metric_type=metric_type
            )
            
            self.metrics.append(metric)
            self.metric_history[name].append(metric)
    
    def increment_counter(self, name: str, value: float = 1.0):
        """Increment a counter metric"""
        with self.lock:
            self.counters[name] = self.counters.get(name, 0) + value
            self.record_metric(name, self.counters[name], metric_type=MetricType.COUNTER)
    
    def record_execution_time(self, execution_type: str, duration: float, labels: Optional[Dict[str, str]] = None):
        """Record execution time metrics"""
        self.record_metric(f"{execution_type}_execution_time", duration, labels, MetricType.HISTOGRAM)
        
        if execution_type == "workflow":
            self.time_series['workflow_execution_times'].append(duration)
        elif execution_type == "task":
            self.time_series['task_execution_times'].append(duration)
    
    def record_resource_utilization(self, resource_type: str, utilization: float, node_id: Optional[str] = None):
        """Record resource utilization metrics"""
        labels = {'resource_type': resource_type}
        if node_id:
            labels['node_id'] = node_id
        
        self.record_metric("resource_utilization", utilization, labels)
        self.time_series['resource_utilization'].append((time.time(), resource_type, utilization))
    
    def get_metric_summary(self, metric_name: str, time_window: float = 3600.0) -> Dict[str, Any]:
        """Get summary statistics for a metric"""
        cutoff_time = time.time() - time_window
        
        with self.lock:
            recent_metrics = [
                m.value for m in self.metric_history[metric_name]
                if m.timestamp >= cutoff_time
            ]
        
        if not recent_metrics:
            return {'count': 0}
        
        return {
            'count': len(recent_metrics),
            'min': min(recent_metrics),
            'max': max(recent_metrics),
            'mean': statistics.mean(recent_metrics),
            'median': statistics.median(recent_metrics),
            'stddev': statistics.stdev(recent_metrics) if len(recent_metrics) > 1 else 0.0,
            'latest': recent_metrics[-1]
        }
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get overall performance summary"""
        with self.lock:
            workflow_times = list(self.time_series['workflow_execution_times'])
            task_times = list(self.time_series['task_execution_times'])
            
            return {
                'counters': self.counters.copy(),
                'workflow_performance': {
                    'count': len(workflow_times),
                    'avg_execution_time': statistics.mean(workflow_times) if workflow_times else 0.0,
                    'median_execution_time': statistics.median(workflow_times) if workflow_times else 0.0,
                    'success_rate': (self.counters['workflows_completed'] / 
                                   max(1, self.counters['workflows_submitted']))
                },
                'task_performance': {
                    'count': len(task_times),
                    'avg_execution_time': statistics.mean(task_times) if task_times else 0.0,
                    'median_execution_time': statistics.median(task_times) if task_times else 0.0,
                    'failure_rate': (self.counters['tasks_failed'] / 
                                   max(1, self.counters['tasks_executed']))
                },
                'current_timestamp': time.time()
            }


class AlertManager:
    """Manages alerts and notifications"""
    
    def __init__(self):
        self.logger = PrismLogger("alert_manager")
        self.alerts: Dict[str, Alert] = {}
        self.alert_rules: List[Dict[str, Any]] = []
        self.alert_handlers: List[Callable[[Alert], None]] = []
        
        # Alert configuration
        self.alert_cooldown = 300.0  # 5 minutes
        self.max_alerts = 1000
        
        self.lock = threading.RLock()
    
    def add_alert_rule(self, rule: Dict[str, Any]):
        """Add an alert rule"""
        required_fields = ['name', 'condition', 'severity', 'description']
        if not all(field in rule for field in required_fields):
            raise ValueError(f"Alert rule must contain: {required_fields}")
        
        self.alert_rules.append(rule)
        self.logger.info(f"Added alert rule: {rule['name']}")
    
    def add_alert_handler(self, handler: Callable[[Alert], None]):
        """Add an alert handler function"""
        self.alert_handlers.append(handler)
        self.logger.info(f"Added alert handler: {handler.__name__}")
    
    def trigger_alert(self, alert_id: str, title: str, description: str, 
                     severity: AlertSeverity, source: str, 
                     labels: Optional[Dict[str, str]] = None) -> bool:
        """Trigger an alert"""
        with self.lock:
            # Check if alert already exists and is not resolved
            if alert_id in self.alerts and not self.alerts[alert_id].resolved:
                return False  # Don't trigger duplicate alerts
            
            alert = Alert(
                alert_id=alert_id,
                title=title,
                description=description,
                severity=severity,
                timestamp=time.time(),
                source=source,
                labels=labels or {}
            )
            
            self.alerts[alert_id] = alert
            
            # Notify handlers
            for handler in self.alert_handlers:
                try:
                    handler(alert)
                except Exception as e:
                    self.logger.error(f"Alert handler failed: {e}")
            
            self.logger.warning(f"Alert triggered: {title} ({severity.value})")
            return True
    
    def resolve_alert(self, alert_id: str) -> bool:
        """Resolve an alert"""
        with self.lock:
            if alert_id not in self.alerts:
                return False
            
            alert = self.alerts[alert_id]
            if not alert.resolved:
                alert.resolve()
                self.logger.info(f"Alert resolved: {alert.title}")
                return True
            
            return False
    
    def check_alert_conditions(self, metrics: ExecutionMetrics):
        """Check all alert rules against current metrics"""
        for rule in self.alert_rules:
            try:
                if self._evaluate_condition(rule['condition'], metrics):
                    alert_id = f"{rule['name']}_{int(time.time())}"
                    self.trigger_alert(
                        alert_id=alert_id,
                        title=rule['name'],
                        description=rule['description'],
                        severity=AlertSeverity(rule['severity']),
                        source="alert_rule_engine"
                    )
            except Exception as e:
                self.logger.error(f"Error evaluating alert rule {rule['name']}: {e}")
    
    def _evaluate_condition(self, condition: str, metrics: ExecutionMetrics) -> bool:
        """Evaluate alert condition (simplified)"""
        # This would be a more sophisticated expression evaluator in practice
        try:
            # Example conditions:
            # "workflow_failure_rate > 0.1"
            # "resource_utilization > 0.9"
            # "queue_length > 100"
            
            if "workflow_failure_rate" in condition:
                total_workflows = metrics.counters.get('workflows_submitted', 0)
                failed_workflows = metrics.counters.get('workflows_failed', 0)
                failure_rate = failed_workflows / max(1, total_workflows)
                
                if "> 0.1" in condition:
                    return failure_rate > 0.1
            
            # Add more condition evaluations as needed
            return False
            
        except Exception:
            return False
    
    def get_active_alerts(self) -> List[Alert]:
        """Get all active (unresolved) alerts"""
        with self.lock:
            return [alert for alert in self.alerts.values() if not alert.resolved]
    
    def get_alert_summary(self) -> Dict[str, Any]:
        """Get alert summary statistics"""
        with self.lock:
            active_alerts = self.get_active_alerts()
            all_alerts = list(self.alerts.values())
            
            severity_counts = defaultdict(int)
            for alert in active_alerts:
                severity_counts[alert.severity.value] += 1
            
            return {
                'total_alerts': len(all_alerts),
                'active_alerts': len(active_alerts),
                'resolved_alerts': len(all_alerts) - len(active_alerts),
                'severity_breakdown': dict(severity_counts),
                'alert_rules_count': len(self.alert_rules)
            }


class TaskMonitor:
    """Monitors individual task execution"""
    
    def __init__(self, metrics: ExecutionMetrics, alert_manager: AlertManager):
        self.metrics = metrics
        self.alert_manager = alert_manager
        self.logger = PrismLogger("task_monitor")
        
        # Task tracking
        self.active_tasks: Dict[str, Dict[str, Any]] = {}
        self.task_history: List[Dict[str, Any]] = []
        
        # Monitoring thresholds
        self.slow_task_threshold = 600.0  # 10 minutes
        self.task_timeout_threshold = 3600.0  # 1 hour
        
        self.lock = threading.RLock()
    
    def start_monitoring_task(self, task: ScheduledTask):
        """Start monitoring a task"""
        with self.lock:
            task_info = {
                'task_id': task.task_id,
                'start_time': time.time(),
                'status': task.status.value,
                'estimated_duration': task.resources.estimated_runtime_minutes * 60,
                'timeout': task.timeout_minutes * 60
            }
            
            self.active_tasks[task.task_id] = task_info
            self.metrics.increment_counter('tasks_executed')
    
    def update_task_status(self, task_id: str, status: TaskStatus, 
                          result: Optional[Any] = None, error: Optional[str] = None):
        """Update task status and metrics"""
        with self.lock:
            if task_id not in self.active_tasks:
                return
            
            task_info = self.active_tasks[task_id]
            task_info['status'] = status.value
            current_time = time.time()
            execution_time = current_time - task_info['start_time']
            
            if status in [TaskStatus.COMPLETED, TaskStatus.FAILED, TaskStatus.CANCELLED]:
                # Task finished
                task_info['end_time'] = current_time
                task_info['execution_time'] = execution_time
                task_info['result'] = result
                task_info['error'] = error
                
                # Record metrics
                self.metrics.record_execution_time('task', execution_time)
                
                if status == TaskStatus.COMPLETED:
                    self.metrics.increment_counter('workflows_completed')
                elif status == TaskStatus.FAILED:
                    self.metrics.increment_counter('tasks_failed')
                    
                    # Trigger failure alert
                    self.alert_manager.trigger_alert(
                        alert_id=f"task_failure_{task_id}",
                        title="Task Execution Failed",
                        description=f"Task {task_id} failed: {error or 'Unknown error'}",
                        severity=AlertSeverity.WARNING,
                        source="task_monitor",
                        labels={'task_id': task_id}
                    )
                
                # Move to history
                self.task_history.append(task_info.copy())
                del self.active_tasks[task_id]
            
            # Check for slow tasks
            if (status == TaskStatus.RUNNING and 
                execution_time > self.slow_task_threshold):
                
                self.alert_manager.trigger_alert(
                    alert_id=f"slow_task_{task_id}",
                    title="Slow Task Execution",
                    description=f"Task {task_id} has been running for {execution_time:.1f} seconds",
                    severity=AlertSeverity.INFO,
                    source="task_monitor",
                    labels={'task_id': task_id}
                )
    
    def check_task_timeouts(self):
        """Check for task timeouts"""
        current_time = time.time()
        
        with self.lock:
            for task_id, task_info in list(self.active_tasks.items()):
                execution_time = current_time - task_info['start_time']
                
                if execution_time > task_info['timeout']:
                    self.alert_manager.trigger_alert(
                        alert_id=f"task_timeout_{task_id}",
                        title="Task Timeout",
                        description=f"Task {task_id} exceeded timeout ({task_info['timeout']} seconds)",
                        severity=AlertSeverity.ERROR,
                        source="task_monitor",
                        labels={'task_id': task_id}
                    )
    
    def get_task_statistics(self) -> Dict[str, Any]:
        """Get task execution statistics"""
        with self.lock:
            recent_tasks = [t for t in self.task_history if 
                          time.time() - t.get('end_time', 0) < 3600]  # Last hour
            
            if recent_tasks:
                execution_times = [t['execution_time'] for t in recent_tasks 
                                 if 'execution_time' in t]
                success_rate = sum(1 for t in recent_tasks if t['status'] == 'completed') / len(recent_tasks)
            else:
                execution_times = []
                success_rate = 0.0
            
            return {
                'active_tasks': len(self.active_tasks),
                'total_tasks_processed': len(self.task_history),
                'recent_tasks_count': len(recent_tasks),
                'success_rate': success_rate,
                'avg_execution_time': statistics.mean(execution_times) if execution_times else 0.0,
                'median_execution_time': statistics.median(execution_times) if execution_times else 0.0
            }


class ResourceMonitor:
    """Monitors resource utilization"""
    
    def __init__(self, metrics: ExecutionMetrics, alert_manager: AlertManager):
        self.metrics = metrics
        self.alert_manager = alert_manager
        self.logger = PrismLogger("resource_monitor")
        
        # Resource tracking
        self.resource_usage_history: Dict[str, deque] = defaultdict(lambda: deque(maxlen=1000))
        
        # Monitoring thresholds
        self.high_utilization_threshold = 0.9
        self.low_utilization_threshold = 0.1
        
    def record_resource_usage(self, node_id: str, resource_data: Dict[str, float]):
        """Record resource usage data"""
        timestamp = time.time()
        
        for resource_type, utilization in resource_data.items():
            # Record metric
            self.metrics.record_resource_utilization(resource_type, utilization, node_id)
            
            # Store in history
            self.resource_usage_history[f"{node_id}_{resource_type}"].append({
                'timestamp': timestamp,
                'utilization': utilization
            })
            
            # Check thresholds
            if utilization > self.high_utilization_threshold:
                self.alert_manager.trigger_alert(
                    alert_id=f"high_resource_usage_{node_id}_{resource_type}",
                    title="High Resource Utilization",
                    description=f"Node {node_id} {resource_type} utilization: {utilization:.1%}",
                    severity=AlertSeverity.WARNING,
                    source="resource_monitor",
                    labels={'node_id': node_id, 'resource_type': resource_type}
                )
    
    def get_resource_trends(self, node_id: str, resource_type: str, time_window: float = 3600.0) -> Dict[str, Any]:
        """Get resource utilization trends"""
        key = f"{node_id}_{resource_type}"
        history = list(self.resource_usage_history[key])
        
        cutoff_time = time.time() - time_window
        recent_data = [entry for entry in history if entry['timestamp'] >= cutoff_time]
        
        if not recent_data:
            return {'trend': 'no_data'}
        
        utilizations = [entry['utilization'] for entry in recent_data]
        
        return {
            'current': utilizations[-1] if utilizations else 0.0,
            'avg': statistics.mean(utilizations),
            'max': max(utilizations),
            'min': min(utilizations),
            'trend': 'increasing' if len(utilizations) > 1 and utilizations[-1] > utilizations[0] else 'stable'
        }


class WorkflowMonitor:
    """Main workflow monitoring coordinator"""
    
    def __init__(self):
        self.logger = PrismLogger("workflow_monitor")
        
        # Core monitoring components
        self.metrics = ExecutionMetrics()
        self.alert_manager = AlertManager()
        self.task_monitor = TaskMonitor(self.metrics, self.alert_manager)
        self.resource_monitor = ResourceMonitor(self.metrics, self.alert_manager)
        
        # Monitoring state
        self.monitoring_active = False
        self.monitor_thread: Optional[threading.Thread] = None
        self.monitoring_interval = 30.0  # seconds
        
        # Setup default alert rules
        self._setup_default_alert_rules()
        
        # Setup default alert handlers
        self._setup_default_alert_handlers()
    
    def _setup_default_alert_rules(self):
        """Setup default alert rules"""
        default_rules = [
            {
                'name': 'High Workflow Failure Rate',
                'condition': 'workflow_failure_rate > 0.2',
                'severity': 'warning',
                'description': 'Workflow failure rate exceeds 20%'
            },
            {
                'name': 'High Task Failure Rate',
                'condition': 'task_failure_rate > 0.1',
                'severity': 'warning',
                'description': 'Task failure rate exceeds 10%'
            }
        ]
        
        for rule in default_rules:
            self.alert_manager.add_alert_rule(rule)
    
    def _setup_default_alert_handlers(self):
        """Setup default alert handlers"""
        def log_alert_handler(alert: Alert):
            self.logger.warning(f"ALERT: {alert.title} - {alert.description}")
        
        self.alert_manager.add_alert_handler(log_alert_handler)
    
    def start_monitoring(self):
        """Start monitoring"""
        if self.monitoring_active:
            return
        
        self.monitoring_active = True
        
        def monitoring_loop():
            while self.monitoring_active:
                try:
                    # Check task timeouts
                    self.task_monitor.check_task_timeouts()
                    
                    # Check alert conditions
                    self.alert_manager.check_alert_conditions(self.metrics)
                    
                    # Sleep until next check
                    time.sleep(self.monitoring_interval)
                    
                except Exception as e:
                    self.logger.error(f"Monitoring error: {e}")
        
        self.monitor_thread = threading.Thread(target=monitoring_loop, daemon=True)
        self.monitor_thread.start()
        
        self.logger.info("Workflow monitoring started")
    
    def stop_monitoring(self):
        """Stop monitoring"""
        self.monitoring_active = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=10.0)
        self.logger.info("Workflow monitoring stopped")
    
    def record_workflow_event(self, event_type: str, workflow_id: str, **kwargs):
        """Record a workflow event"""
        self.metrics.record_metric(
            name=f"workflow_{event_type}",
            value=1.0,
            labels={'workflow_id': workflow_id, **kwargs},
            metric_type=MetricType.COUNTER
        )
        
        if event_type == "submitted":
            self.metrics.increment_counter('workflows_submitted')
        elif event_type == "completed":
            self.metrics.increment_counter('workflows_completed')
            if 'execution_time' in kwargs:
                self.metrics.record_execution_time('workflow', kwargs['execution_time'])
        elif event_type == "failed":
            self.metrics.increment_counter('workflows_failed')
    
    def get_monitoring_dashboard(self) -> Dict[str, Any]:
        """Get comprehensive monitoring dashboard data"""
        return {
            'timestamp': time.time(),
            'system_performance': self.metrics.get_performance_summary(),
            'task_statistics': self.task_monitor.get_task_statistics(),
            'alert_summary': self.alert_manager.get_alert_summary(),
            'active_alerts': [alert.to_dict() for alert in self.alert_manager.get_active_alerts()],
            'recent_metrics': {
                name: self.metrics.get_metric_summary(name)
                for name in ['workflow_execution_time', 'task_execution_time', 'resource_utilization']
            }
        }


# Convenience functions
def create_workflow_monitor() -> WorkflowMonitor:
    """Create a workflow monitor with default configuration"""
    monitor = WorkflowMonitor()
    monitor.start_monitoring()
    return monitor


def setup_monitoring_for_orchestrator(orchestrator: WorkflowOrchestrator) -> WorkflowMonitor:
    """Setup monitoring for an existing orchestrator"""
    monitor = create_workflow_monitor()
    
    # Add custom alert handler to integrate with orchestrator
    def orchestrator_alert_handler(alert: Alert):
        # In a real system, this could trigger orchestrator actions
        # like rescheduling workflows, scaling resources, etc.
        pass
    
    monitor.alert_manager.add_alert_handler(orchestrator_alert_handler)
    return monitor