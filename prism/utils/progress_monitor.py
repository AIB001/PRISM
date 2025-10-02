#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Progress Monitoring and Performance Tracking

Provides comprehensive progress tracking, performance monitoring, and
real-time status updates for PRISM molecular dynamics simulations.
"""

import time
import threading
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, List, Callable, Union
from dataclasses import dataclass, asdict
from contextlib import contextmanager
from enum import Enum
import json
import subprocess
import re


class WorkflowStage(Enum):
    """Workflow stages for progress tracking"""
    INITIALIZATION = "initialization"
    VALIDATION = "validation"
    SYSTEM_BUILDING = "system_building"
    EQUILIBRATION = "equilibration"
    SMD_PREPARATION = "smd_preparation"
    SMD_EXECUTION = "smd_execution"
    UMBRELLA_PREPARATION = "umbrella_preparation"
    UMBRELLA_EXECUTION = "umbrella_execution"
    ANALYSIS = "analysis"
    VISUALIZATION = "visualization"
    CLEANUP = "cleanup"
    COMPLETED = "completed"
    FAILED = "failed"


class StageStatus(Enum):
    """Status of workflow stages"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class StageProgress:
    """Progress information for a workflow stage"""
    stage: WorkflowStage
    status: StageStatus
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    progress_percent: float = 0.0
    current_step: str = ""
    total_steps: int = 0
    completed_steps: int = 0
    estimated_duration: Optional[float] = None
    actual_duration: Optional[float] = None
    error_message: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data['stage'] = self.stage.value
        data['status'] = self.status.value
        if self.start_time:
            data['start_time'] = self.start_time.isoformat()
        if self.end_time:
            data['end_time'] = self.end_time.isoformat()
        return data
    
    @property
    def duration(self) -> Optional[float]:
        """Get duration in seconds"""
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds()
        elif self.start_time:
            return (datetime.now() - self.start_time).total_seconds()
        return None


@dataclass
class SimulationProgress:
    """Overall simulation progress tracking"""
    simulation_id: str
    simulation_type: str
    start_time: datetime
    current_stage: WorkflowStage
    overall_progress: float = 0.0
    stages: Dict[WorkflowStage, StageProgress] = None
    estimated_completion: Optional[datetime] = None
    
    def __post_init__(self):
        if self.stages is None:
            self.stages = {}
            # Initialize all stages
            for stage in WorkflowStage:
                if stage not in [WorkflowStage.COMPLETED, WorkflowStage.FAILED]:
                    self.stages[stage] = StageProgress(stage, StageStatus.PENDING)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'simulation_id': self.simulation_id,
            'simulation_type': self.simulation_type,
            'start_time': self.start_time.isoformat(),
            'current_stage': self.current_stage.value,
            'overall_progress': self.overall_progress,
            'stages': {stage.value: progress.to_dict() 
                      for stage, progress in self.stages.items()},
            'estimated_completion': self.estimated_completion.isoformat() 
                                  if self.estimated_completion else None
        }


class ProgressTracker:
    """Progress tracker for PRISM simulations"""
    
    def __init__(self, simulation_id: str, simulation_type: str):
        self.simulation_id = simulation_id
        self.simulation_type = simulation_type
        self.progress = SimulationProgress(
            simulation_id=simulation_id,
            simulation_type=simulation_type,
            start_time=datetime.now(),
            current_stage=WorkflowStage.INITIALIZATION
        )
        self.callbacks: List[Callable[[SimulationProgress], None]] = []
        self.stage_weights = self._get_stage_weights()
    
    def _get_stage_weights(self) -> Dict[WorkflowStage, float]:
        """Get relative weights for different stages"""
        if self.simulation_type.lower() == "pmf":
            return {
                WorkflowStage.INITIALIZATION: 0.02,
                WorkflowStage.VALIDATION: 0.03,
                WorkflowStage.SYSTEM_BUILDING: 0.10,
                WorkflowStage.EQUILIBRATION: 0.15,
                WorkflowStage.SMD_PREPARATION: 0.05,
                WorkflowStage.SMD_EXECUTION: 0.20,
                WorkflowStage.UMBRELLA_PREPARATION: 0.05,
                WorkflowStage.UMBRELLA_EXECUTION: 0.30,
                WorkflowStage.ANALYSIS: 0.08,
                WorkflowStage.VISUALIZATION: 0.02
            }
        else:
            # Default weights for other simulation types
            return {stage: 1.0 / len(WorkflowStage) for stage in WorkflowStage 
                   if stage not in [WorkflowStage.COMPLETED, WorkflowStage.FAILED]}
    
    def add_callback(self, callback: Callable[[SimulationProgress], None]):
        """Add callback for progress updates"""
        self.callbacks.append(callback)
    
    def _notify_callbacks(self):
        """Notify all callbacks of progress update"""
        for callback in self.callbacks:
            try:
                callback(self.progress)
            except Exception as e:
                print(f"Callback error: {e}")
    
    def start_stage(self, stage: WorkflowStage, estimated_duration: Optional[float] = None,
                   total_steps: int = 0):
        """Start a new workflow stage"""
        if stage in self.progress.stages:
            stage_progress = self.progress.stages[stage]
            stage_progress.status = StageStatus.RUNNING
            stage_progress.start_time = datetime.now()
            stage_progress.estimated_duration = estimated_duration
            stage_progress.total_steps = total_steps
            stage_progress.completed_steps = 0
            stage_progress.progress_percent = 0.0
        
        self.progress.current_stage = stage
        self._update_overall_progress()
        self._notify_callbacks()
    
    def update_stage_progress(self, stage: WorkflowStage, progress_percent: float = None,
                            current_step: str = "", completed_steps: int = None):
        """Update progress for current stage"""
        if stage in self.progress.stages:
            stage_progress = self.progress.stages[stage]
            
            if progress_percent is not None:
                stage_progress.progress_percent = min(100.0, max(0.0, progress_percent))
            
            if current_step:
                stage_progress.current_step = current_step
            
            if completed_steps is not None:
                stage_progress.completed_steps = completed_steps
                if stage_progress.total_steps > 0:
                    stage_progress.progress_percent = (completed_steps / stage_progress.total_steps) * 100
        
        self._update_overall_progress()
        self._estimate_completion()
        self._notify_callbacks()
    
    def complete_stage(self, stage: WorkflowStage):
        """Mark stage as completed"""
        if stage in self.progress.stages:
            stage_progress = self.progress.stages[stage]
            stage_progress.status = StageStatus.COMPLETED
            stage_progress.end_time = datetime.now()
            stage_progress.progress_percent = 100.0
            stage_progress.actual_duration = stage_progress.duration
        
        self._update_overall_progress()
        self._notify_callbacks()
    
    def fail_stage(self, stage: WorkflowStage, error_message: str):
        """Mark stage as failed"""
        if stage in self.progress.stages:
            stage_progress = self.progress.stages[stage]
            stage_progress.status = StageStatus.FAILED
            stage_progress.end_time = datetime.now()
            stage_progress.error_message = error_message
        
        self.progress.current_stage = WorkflowStage.FAILED
        self._notify_callbacks()
    
    def skip_stage(self, stage: WorkflowStage, reason: str = ""):
        """Mark stage as skipped"""
        if stage in self.progress.stages:
            stage_progress = self.progress.stages[stage]
            stage_progress.status = StageStatus.SKIPPED
            stage_progress.error_message = reason
        
        self._update_overall_progress()
        self._notify_callbacks()
    
    def _update_overall_progress(self):
        """Update overall progress based on stage progress"""
        total_weight = 0.0
        weighted_progress = 0.0
        
        for stage, stage_progress in self.progress.stages.items():
            if stage in self.stage_weights:
                weight = self.stage_weights[stage]
                total_weight += weight
                
                if stage_progress.status == StageStatus.COMPLETED:
                    weighted_progress += weight * 100
                elif stage_progress.status == StageStatus.RUNNING:
                    weighted_progress += weight * stage_progress.progress_percent
                elif stage_progress.status == StageStatus.SKIPPED:
                    total_weight -= weight  # Don't count skipped stages
        
        if total_weight > 0:
            self.progress.overall_progress = weighted_progress / total_weight
        else:
            self.progress.overall_progress = 0.0
    
    def _estimate_completion(self):
        """Estimate completion time based on current progress"""
        if self.progress.overall_progress <= 0:
            return
        
        elapsed = (datetime.now() - self.progress.start_time).total_seconds()
        total_estimated = elapsed / (self.progress.overall_progress / 100.0)
        remaining = total_estimated - elapsed
        
        if remaining > 0:
            self.progress.estimated_completion = datetime.now() + timedelta(seconds=remaining)
        else:
            self.progress.estimated_completion = None
    
    def get_status_summary(self) -> Dict[str, Any]:
        """Get comprehensive status summary"""
        elapsed = (datetime.now() - self.progress.start_time).total_seconds()
        
        # Stage summaries
        stage_summaries = {}
        for stage, stage_progress in self.progress.stages.items():
            stage_summaries[stage.value] = {
                'status': stage_progress.status.value,
                'progress': stage_progress.progress_percent,
                'duration': stage_progress.duration,
                'current_step': stage_progress.current_step
            }
        
        return {
            'simulation_id': self.simulation_id,
            'simulation_type': self.simulation_type,
            'current_stage': self.progress.current_stage.value,
            'overall_progress': self.progress.overall_progress,
            'elapsed_time': elapsed,
            'estimated_completion': self.progress.estimated_completion.isoformat() 
                                  if self.progress.estimated_completion else None,
            'stages': stage_summaries
        }


class PerformanceMonitor:
    """Performance monitoring for GROMACS simulations"""
    
    def __init__(self, log_file: Optional[Path] = None):
        self.log_file = log_file
        self.monitoring_active = False
        self.monitoring_thread = None
        self.performance_data: List[Dict[str, Any]] = []
    
    def start_monitoring_md_log(self, md_log_file: Path, callback: Optional[Callable] = None):
        """Monitor MD simulation log file for progress"""
        if not md_log_file.exists():
            return
        
        self.monitoring_active = True
        self.monitoring_thread = threading.Thread(
            target=self._monitor_md_log, 
            args=(md_log_file, callback),
            daemon=True
        )
        self.monitoring_thread.start()
    
    def stop_monitoring(self):
        """Stop monitoring"""
        self.monitoring_active = False
        if self.monitoring_thread:
            self.monitoring_thread.join(timeout=5)
    
    def _monitor_md_log(self, log_file: Path, callback: Optional[Callable] = None):
        """Monitor GROMACS MD log file"""
        last_position = 0
        
        while self.monitoring_active:
            try:
                if log_file.exists():
                    with open(log_file, 'r') as f:
                        f.seek(last_position)
                        new_content = f.read()
                        if new_content:
                            self._parse_md_progress(new_content, callback)
                        last_position = f.tell()
                
                time.sleep(5)  # Check every 5 seconds
                
            except Exception as e:
                print(f"Error monitoring log: {e}")
                time.sleep(10)
    
    def _parse_md_progress(self, content: str, callback: Optional[Callable] = None):
        """Parse GROMACS log content for progress information"""
        lines = content.strip().split('\n')
        
        for line in lines:
            # Parse step information
            step_match = re.search(r'Step\s+(\d+)', line)
            if step_match:
                current_step = int(step_match.group(1))
                
                # Extract performance info if available
                perf_data = {'step': current_step, 'timestamp': datetime.now()}
                
                # Look for performance metrics in the same line
                if 'ns/day' in line:
                    ns_day_match = re.search(r'([\d.]+)\s+ns/day', line)
                    if ns_day_match:
                        perf_data['ns_per_day'] = float(ns_day_match.group(1))
                
                if 'hour/ns' in line:
                    hour_ns_match = re.search(r'([\d.]+)\s+hour/ns', line)
                    if hour_ns_match:
                        perf_data['hours_per_ns'] = float(hour_ns_match.group(1))
                
                self.performance_data.append(perf_data)
                
                if callback:
                    callback(perf_data)


@contextmanager
def track_progress(simulation_id: str, simulation_type: str = "MD",
                  progress_file: Optional[Path] = None):
    """Context manager for progress tracking"""
    tracker = ProgressTracker(simulation_id, simulation_type)
    
    # Save progress to file if specified
    if progress_file:
        def save_progress(progress: SimulationProgress):
            try:
                with open(progress_file, 'w') as f:
                    json.dump(progress.to_dict(), f, indent=2)
            except Exception as e:
                print(f"Error saving progress: {e}")
        
        tracker.add_callback(save_progress)
    
    try:
        tracker.start_stage(WorkflowStage.INITIALIZATION)
        yield tracker
        
        # Mark as completed if we reach this point
        tracker.progress.current_stage = WorkflowStage.COMPLETED
        tracker._notify_callbacks()
        
    except Exception as e:
        # Mark as failed
        tracker.fail_stage(tracker.progress.current_stage, str(e))
        raise


class ProgressReporter:
    """Generate progress reports and summaries"""
    
    @staticmethod
    def generate_text_report(progress: SimulationProgress) -> str:
        """Generate human-readable progress report"""
        lines = []
        lines.append("=" * 60)
        lines.append(f"PRISM SIMULATION PROGRESS REPORT")
        lines.append("=" * 60)
        lines.append(f"Simulation ID: {progress.simulation_id}")
        lines.append(f"Type: {progress.simulation_type}")
        lines.append(f"Started: {progress.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Current Stage: {progress.current_stage.value.replace('_', ' ').title()}")
        lines.append(f"Overall Progress: {progress.overall_progress:.1f}%")
        
        if progress.estimated_completion:
            lines.append(f"Estimated Completion: {progress.estimated_completion.strftime('%Y-%m-%d %H:%M:%S')}")
        
        lines.append("")
        lines.append("STAGE DETAILS:")
        lines.append("-" * 40)
        
        for stage, stage_progress in progress.stages.items():
            if stage_progress.status != StageStatus.PENDING:
                status_icon = {
                    StageStatus.RUNNING: "🔄",
                    StageStatus.COMPLETED: "✅",
                    StageStatus.FAILED: "❌",
                    StageStatus.SKIPPED: "⏭️"
                }.get(stage_progress.status, "⚪")
                
                stage_name = stage.value.replace('_', ' ').title()
                lines.append(f"{status_icon} {stage_name}")
                lines.append(f"   Progress: {stage_progress.progress_percent:.1f}%")
                
                if stage_progress.current_step:
                    lines.append(f"   Current: {stage_progress.current_step}")
                
                if stage_progress.duration:
                    lines.append(f"   Duration: {stage_progress.duration:.1f}s")
                
                if stage_progress.error_message:
                    lines.append(f"   Error: {stage_progress.error_message}")
                
                lines.append("")
        
        lines.append("=" * 60)
        return "\n".join(lines)
    
    @staticmethod
    def generate_html_report(progress: SimulationProgress) -> str:
        """Generate HTML progress report"""
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>PRISM Progress Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 10px; border-radius: 5px; }}
                .progress-bar {{ 
                    width: 100%; height: 20px; background-color: #ddd; 
                    border-radius: 10px; overflow: hidden; 
                }}
                .progress-fill {{ 
                    height: 100%; background-color: #4CAF50; 
                    width: {progress.overall_progress}%; 
                }}
                .stage {{ margin: 10px 0; padding: 10px; border-left: 4px solid #ccc; }}
                .stage.running {{ border-color: #2196F3; }}
                .stage.completed {{ border-color: #4CAF50; }}
                .stage.failed {{ border-color: #F44336; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>PRISM Simulation Progress</h1>
                <p><strong>ID:</strong> {progress.simulation_id}</p>
                <p><strong>Type:</strong> {progress.simulation_type}</p>
                <p><strong>Started:</strong> {progress.start_time.strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>Current Stage:</strong> {progress.current_stage.value.replace('_', ' ').title()}</p>
            </div>
            
            <h2>Overall Progress: {progress.overall_progress:.1f}%</h2>
            <div class="progress-bar">
                <div class="progress-fill"></div>
            </div>
            
            <h2>Stage Details</h2>
        """
        
        for stage, stage_progress in progress.stages.items():
            if stage_progress.status != StageStatus.PENDING:
                status_class = stage_progress.status.value
                stage_name = stage.value.replace('_', ' ').title()
                
                html += f"""
                <div class="stage {status_class}">
                    <h3>{stage_name} ({stage_progress.status.value.title()})</h3>
                    <p>Progress: {stage_progress.progress_percent:.1f}%</p>
                """
                
                if stage_progress.current_step:
                    html += f"<p>Current: {stage_progress.current_step}</p>"
                
                if stage_progress.duration:
                    html += f"<p>Duration: {stage_progress.duration:.1f}s</p>"
                
                html += "</div>"
        
        html += """
        </body>
        </html>
        """
        return html


# Convenience functions
def create_progress_tracker(simulation_id: str, simulation_type: str = "MD") -> ProgressTracker:
    """Create a progress tracker"""
    return ProgressTracker(simulation_id, simulation_type)


def monitor_gromacs_simulation(log_file: Path, progress_callback: Optional[Callable] = None) -> PerformanceMonitor:
    """Monitor GROMACS simulation progress"""
    monitor = PerformanceMonitor()
    monitor.start_monitoring_md_log(log_file, progress_callback)
    return monitor