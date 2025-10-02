#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Core API

Core API classes for PMF calculations, workflow management, and data analysis.
"""

import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Callable
from dataclasses import dataclass, asdict
from enum import Enum

from ..utils.logging_system import PrismLogger
from .exceptions import *


class CalculationStatus(Enum):
    """Enumeration of calculation status states"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class WorkflowPriority(Enum):
    """Enumeration of workflow priority levels"""
    LOW = "low"
    NORMAL = "normal"
    HIGH = "high"
    CRITICAL = "critical"


@dataclass
class PMFCalculationConfig:
    """Configuration for PMF calculation"""
    system_file: str
    ligand: str
    protein: str
    method: str = "umbrella_sampling"
    num_windows: int = 30
    window_spacing: float = 0.1
    force_constant: float = 1000.0
    equilibration_time: int = 1000  # ps
    sampling_time: int = 10000  # ps
    temperature: float = 298.15  # K
    pressure: float = 1.0  # bar
    output_directory: Optional[str] = None
    
    def validate(self):
        """Validate configuration parameters"""
        if not Path(self.system_file).exists():
            raise ValidationError(f"System file not found: {self.system_file}", "system_file")
        
        if self.num_windows <= 0:
            raise ValidationError("Number of windows must be positive", "num_windows", self.num_windows)
        
        if self.force_constant <= 0:
            raise ValidationError("Force constant must be positive", "force_constant", self.force_constant)
        
        if self.method not in ["umbrella_sampling", "smd", "adiabatic_bias"]:
            raise ValidationError(f"Unsupported method: {self.method}", "method", self.method)


@dataclass 
class PMFResult:
    """Results from PMF calculation"""
    calculation_id: str
    status: CalculationStatus
    config: PMFCalculationConfig
    start_time: float
    end_time: Optional[float] = None
    pmf_profile: Optional[List[Dict]] = None
    binding_energy: Optional[float] = None
    convergence_achieved: bool = False
    error_estimate: Optional[float] = None
    output_files: Optional[Dict[str, str]] = None
    error_message: Optional[str] = None
    
    @property
    def runtime(self) -> Optional[float]:
        """Calculate runtime in seconds"""
        if self.end_time:
            return self.end_time - self.start_time
        return None
    
    @property
    def is_complete(self) -> bool:
        """Check if calculation is complete"""
        return self.status in [CalculationStatus.COMPLETED, CalculationStatus.FAILED]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary"""
        result = asdict(self)
        result['status'] = self.status.value
        result['config'] = asdict(self.config)
        return result


class PMFCalculator:
    """Main class for PMF calculations"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger("pmf_calculator")
        self._active_calculations: Dict[str, PMFResult] = {}
        self._calculation_callbacks: Dict[str, List[Callable]] = {}
    
    def create_calculation(self, config: Union[PMFCalculationConfig, Dict[str, Any]]) -> str:
        """Create a new PMF calculation"""
        
        # Convert dict to config object if needed
        if isinstance(config, dict):
            config = PMFCalculationConfig(**config)
        
        # Validate configuration
        config.validate()
        
        # Generate unique calculation ID
        calc_id = f"pmf_{uuid.uuid4().hex[:8]}"
        
        # Create result object
        result = PMFResult(
            calculation_id=calc_id,
            status=CalculationStatus.PENDING,
            config=config,
            start_time=time.time()
        )
        
        self._active_calculations[calc_id] = result
        
        self.logger.info(f"Created PMF calculation: {calc_id}")
        return calc_id
    
    def submit_calculation(self, calculation_id: str, blocking: bool = False) -> PMFResult:
        """Submit PMF calculation for execution"""
        
        if calculation_id not in self._active_calculations:
            raise ValidationError(f"Calculation not found: {calculation_id}", "calculation_id")
        
        result = self._active_calculations[calculation_id]
        
        if result.status != CalculationStatus.PENDING:
            raise ValidationError(f"Calculation already submitted: {calculation_id}")
        
        self.logger.info(f"Submitting PMF calculation: {calculation_id}")
        
        try:
            # Update status
            result.status = CalculationStatus.RUNNING
            
            # Execute calculation (mock implementation)
            if blocking:
                self._execute_calculation_sync(result)
            else:
                self._execute_calculation_async(result)
            
            return result
            
        except Exception as e:
            result.status = CalculationStatus.FAILED
            result.error_message = str(e)
            result.end_time = time.time()
            self.logger.error(f"Calculation failed: {calculation_id} - {e}")
            raise CalculationError(f"PMF calculation failed: {e}", "execution")
    
    def get_result(self, calculation_id: str) -> PMFResult:
        """Get calculation result"""
        
        if calculation_id not in self._active_calculations:
            raise ValidationError(f"Calculation not found: {calculation_id}", "calculation_id")
        
        return self._active_calculations[calculation_id]
    
    def cancel_calculation(self, calculation_id: str) -> bool:
        """Cancel running calculation"""
        
        if calculation_id not in self._active_calculations:
            raise ValidationError(f"Calculation not found: {calculation_id}", "calculation_id")
        
        result = self._active_calculations[calculation_id]
        
        if result.status in [CalculationStatus.COMPLETED, CalculationStatus.FAILED]:
            return False
        
        result.status = CalculationStatus.CANCELLED
        result.end_time = time.time()
        
        self.logger.info(f"Cancelled calculation: {calculation_id}")
        return True
    
    def list_calculations(self, status_filter: Optional[CalculationStatus] = None) -> List[PMFResult]:
        """List all calculations with optional status filter"""
        
        calculations = list(self._active_calculations.values())
        
        if status_filter:
            calculations = [calc for calc in calculations if calc.status == status_filter]
        
        return sorted(calculations, key=lambda x: x.start_time, reverse=True)
    
    def add_callback(self, calculation_id: str, callback: Callable[[PMFResult], None]):
        """Add callback for calculation completion"""
        
        if calculation_id not in self._calculation_callbacks:
            self._calculation_callbacks[calculation_id] = []
        
        self._calculation_callbacks[calculation_id].append(callback)
    
    def run_pmf_calculation(self, system_file: str, ligand: str, protein: str, **kwargs) -> PMFResult:
        """Convenience method to run complete PMF calculation"""
        
        # Create configuration
        config = PMFCalculationConfig(
            system_file=system_file,
            ligand=ligand,
            protein=protein,
            **kwargs
        )
        
        # Create and submit calculation
        calc_id = self.create_calculation(config)
        result = self.submit_calculation(calc_id, blocking=True)
        
        return result
    
    def _execute_calculation_sync(self, result: PMFResult):
        """Execute calculation synchronously (mock implementation)"""
        
        config = result.config
        
        self.logger.info(f"Starting PMF calculation: {result.calculation_id}")
        
        # Simulate calculation steps
        steps = [
            ("System preparation", 2),
            ("Topology generation", 1),
            ("Energy minimization", 3),
            ("Equilibration", 5),
            ("PMF sampling", 10),
            ("Analysis and convergence check", 2)
        ]
        
        try:
            for step_name, duration in steps:
                self.logger.info(f"  {step_name}...")
                time.sleep(duration * 0.1)  # Simulate work
                
                # Call progress callbacks
                if result.calculation_id in self._calculation_callbacks:
                    for callback in self._calculation_callbacks[result.calculation_id]:
                        try:
                            callback(result)
                        except Exception as e:
                            self.logger.warning(f"Callback error: {e}")
            
            # Generate mock results
            result.pmf_profile = self._generate_mock_pmf_profile(config.num_windows)
            result.binding_energy = -8.5 + (hash(result.calculation_id) % 100) * 0.1
            result.convergence_achieved = True
            result.error_estimate = 0.5
            result.output_files = {
                'pmf_profile': f"{config.output_directory or '.'}/pmf_profile.dat",
                'umbrella_histograms': f"{config.output_directory or '.'}/histograms.dat",
                'analysis_log': f"{config.output_directory or '.'}/analysis.log"
            }
            
            result.status = CalculationStatus.COMPLETED
            result.end_time = time.time()
            
            self.logger.info(f"PMF calculation completed: {result.calculation_id}")
            
        except Exception as e:
            result.status = CalculationStatus.FAILED
            result.error_message = str(e)
            result.end_time = time.time()
            raise
    
    def _execute_calculation_async(self, result: PMFResult):
        """Execute calculation asynchronously (mock implementation)"""
        # For this demo, just mark as completed after short delay
        import threading
        
        def async_execution():
            time.sleep(2)
            self._execute_calculation_sync(result)
        
        thread = threading.Thread(target=async_execution)
        thread.daemon = True
        thread.start()
    
    def _generate_mock_pmf_profile(self, num_windows: int) -> List[Dict]:
        """Generate mock PMF profile data"""
        profile = []
        
        for i in range(num_windows):
            x = i * 0.1  # reaction coordinate
            # Gaussian-like potential with minimum at x=1.5
            energy = 20 * ((x - 1.5) ** 2) - 8.5
            error = 0.3 + 0.2 * abs(x - 1.5)  # Higher error at extremes
            
            profile.append({
                'reaction_coordinate': x,
                'pmf_energy': energy,
                'error_estimate': error,
                'sampling_count': 10000 - i * 100
            })
        
        return profile


class WorkflowManager:
    """Manager for PMF workflow orchestration"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger("workflow_manager")
        self._workflows: Dict[str, Dict[str, Any]] = {}
        self.pmf_calculator = PMFCalculator(logger)
    
    def create_workflow(self, config: Dict[str, Any], priority: WorkflowPriority = WorkflowPriority.NORMAL) -> str:
        """Create a new workflow"""
        
        workflow_id = f"workflow_{uuid.uuid4().hex[:8]}"
        
        workflow = {
            'id': workflow_id,
            'name': config.get('name', f'PMF Workflow {workflow_id}'),
            'priority': priority,
            'status': 'created',
            'config': config,
            'tasks': [],
            'created_time': time.time(),
            'submitted_time': None,
            'start_time': None,
            'end_time': None,
            'progress': 0.0,
            'results': {}
        }
        
        self._workflows[workflow_id] = workflow
        
        self.logger.info(f"Created workflow: {workflow_id}")
        return workflow_id
    
    def submit_workflow(self, workflow_id: str) -> bool:
        """Submit workflow for execution"""
        
        if workflow_id not in self._workflows:
            raise ValidationError(f"Workflow not found: {workflow_id}", "workflow_id")
        
        workflow = self._workflows[workflow_id]
        
        if workflow['status'] != 'created':
            raise ValidationError(f"Workflow already submitted: {workflow_id}")
        
        workflow['status'] = 'submitted'
        workflow['submitted_time'] = time.time()
        
        self.logger.info(f"Submitted workflow: {workflow_id}")
        return True
    
    def get_status(self, workflow_id: str) -> Dict[str, Any]:
        """Get workflow status"""
        
        if workflow_id not in self._workflows:
            raise ValidationError(f"Workflow not found: {workflow_id}", "workflow_id")
        
        return self._workflows[workflow_id].copy()
    
    def list_workflows(self, status_filter: Optional[str] = None) -> List[Dict[str, Any]]:
        """List all workflows with optional status filter"""
        
        workflows = list(self._workflows.values())
        
        if status_filter:
            workflows = [wf for wf in workflows if wf['status'] == status_filter]
        
        return sorted(workflows, key=lambda x: x['created_time'], reverse=True)
    
    def cancel_workflow(self, workflow_id: str) -> bool:
        """Cancel workflow execution"""
        
        if workflow_id not in self._workflows:
            raise ValidationError(f"Workflow not found: {workflow_id}", "workflow_id")
        
        workflow = self._workflows[workflow_id]
        
        if workflow['status'] in ['completed', 'failed', 'cancelled']:
            return False
        
        workflow['status'] = 'cancelled'
        workflow['end_time'] = time.time()
        
        self.logger.info(f"Cancelled workflow: {workflow_id}")
        return True


class DataAnalyzer:
    """Analyzer for PMF results and simulation data"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger("data_analyzer")
    
    def analyze_pmf_profile(self, pmf_data: List[Dict]) -> Dict[str, Any]:
        """Analyze PMF profile data"""
        
        if not pmf_data:
            raise DataError("Empty PMF profile data")
        
        energies = [point['pmf_energy'] for point in pmf_data]
        coordinates = [point['reaction_coordinate'] for point in pmf_data]
        
        analysis = {
            'num_points': len(pmf_data),
            'coordinate_range': (min(coordinates), max(coordinates)),
            'energy_range': (min(energies), max(energies)),
            'binding_energy': max(energies) - min(energies),
            'minimum_energy': min(energies),
            'maximum_energy': max(energies),
            'energy_barrier': max(energies) - energies[0] if energies else 0,
            'convergence_quality': self._assess_convergence(pmf_data)
        }
        
        return analysis
    
    def compare_profiles(self, profiles: List[List[Dict]], labels: Optional[List[str]] = None) -> Dict[str, Any]:
        """Compare multiple PMF profiles"""
        
        if not profiles:
            raise DataError("No profiles provided for comparison")
        
        labels = labels or [f"Profile {i+1}" for i in range(len(profiles))]
        
        analyses = []
        for i, profile in enumerate(profiles):
            analysis = self.analyze_pmf_profile(profile)
            analysis['label'] = labels[i]
            analyses.append(analysis)
        
        comparison = {
            'num_profiles': len(profiles),
            'individual_analyses': analyses,
            'binding_energy_comparison': [a['binding_energy'] for a in analyses],
            'energy_range_comparison': [a['energy_range'] for a in analyses],
            'statistical_summary': self._calculate_profile_statistics(analyses)
        }
        
        return comparison
    
    def calculate_statistics(self, data: List[float]) -> Dict[str, float]:
        """Calculate basic statistics for data array"""
        
        if not data:
            raise DataError("Empty data array for statistics")
        
        import statistics
        
        stats = {
            'count': len(data),
            'mean': statistics.mean(data),
            'median': statistics.median(data),
            'stdev': statistics.stdev(data) if len(data) > 1 else 0.0,
            'min': min(data),
            'max': max(data),
            'range': max(data) - min(data)
        }
        
        return stats
    
    def _assess_convergence(self, pmf_data: List[Dict]) -> str:
        """Assess convergence quality of PMF profile"""
        
        if len(pmf_data) < 10:
            return "insufficient_data"
        
        # Simple convergence assessment based on error estimates
        errors = [point.get('error_estimate', 1.0) for point in pmf_data]
        avg_error = sum(errors) / len(errors)
        
        if avg_error < 0.5:
            return "excellent"
        elif avg_error < 1.0:
            return "good"
        elif avg_error < 2.0:
            return "acceptable"
        else:
            return "poor"
    
    def _calculate_profile_statistics(self, analyses: List[Dict]) -> Dict[str, float]:
        """Calculate statistical summary of multiple profile analyses"""
        
        binding_energies = [a['binding_energy'] for a in analyses]
        
        if not binding_energies:
            return {}
        
        return self.calculate_statistics(binding_energies)