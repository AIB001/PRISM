#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Interface APIs

Interface classes for internal module communication and system integration.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Any, Union, Tuple
from pathlib import Path
from dataclasses import dataclass
from enum import Enum

from ..utils.logging_system import PrismLogger
from .exceptions import *


class SimulationStatus(Enum):
    """Status of molecular dynamics simulation"""
    INITIALIZED = "initialized"
    RUNNING = "running" 
    PAUSED = "paused"
    COMPLETED = "completed"
    FAILED = "failed"


@dataclass
class SystemConfig:
    """System configuration for MD simulations"""
    topology_file: str
    coordinate_file: str
    parameter_files: List[str]
    system_name: str
    temperature: float = 298.15
    pressure: float = 1.0
    solvent: str = "water"
    ion_concentration: float = 0.15


@dataclass
class SimulationParameters:
    """Parameters for MD simulation"""
    timestep: float = 0.002  # ps
    total_time: int = 1000  # ps
    output_frequency: int = 100
    trajectory_format: str = "xtc"
    energy_output: bool = True
    restart_frequency: int = 1000


class BaseInterface(ABC):
    """Base interface for all PRISM module interfaces"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger(self.__class__.__name__.lower())
        self._initialized = False
        self._status = "created"
    
    @abstractmethod
    def initialize(self, config: Dict[str, Any]) -> bool:
        """Initialize the interface"""
        pass
    
    @abstractmethod
    def shutdown(self) -> bool:
        """Shutdown the interface"""
        pass
    
    @property
    def is_initialized(self) -> bool:
        """Check if interface is initialized"""
        return self._initialized
    
    @property
    def status(self) -> str:
        """Get current interface status"""
        return self._status


class PMFInterface(BaseInterface):
    """Interface for PMF-specific calculations and analysis"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        super().__init__(logger)
        self._umbrella_windows = []
        self._force_constants = {}
        self._pmf_profile = None
    
    def initialize(self, config: Dict[str, Any]) -> bool:
        """Initialize PMF interface"""
        try:
            self.logger.info("Initializing PMF interface")
            
            # Validate required configuration
            required_keys = ['method', 'num_windows', 'force_constant']
            for key in required_keys:
                if key not in config:
                    raise ConfigurationError(f"Missing required PMF config: {key}", key)
            
            self._method = config['method']
            self._num_windows = config['num_windows']
            self._base_force_constant = config['force_constant']
            
            self._initialized = True
            self._status = "initialized"
            
            self.logger.info(f"PMF interface initialized: {self._method} with {self._num_windows} windows")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to initialize PMF interface: {e}")
            self._status = "error"
            return False
    
    def setup_umbrella_sampling(self, windows: List[float], force_constants: Optional[List[float]] = None) -> bool:
        """Setup umbrella sampling windows"""
        
        if not self._initialized:
            raise APIError("PMF interface not initialized")
        
        self._umbrella_windows = windows.copy()
        
        # Use provided force constants or default values
        if force_constants:
            if len(force_constants) != len(windows):
                raise ValidationError("Force constants count must match windows count")
            self._force_constants = dict(zip(windows, force_constants))
        else:
            # Use base force constant for all windows
            self._force_constants = {w: self._base_force_constant for w in windows}
        
        self.logger.info(f"Setup {len(windows)} umbrella sampling windows")
        return True
    
    def setup_smd_pulling(self, start_position: float, end_position: float, 
                         pulling_rate: float, spring_constant: float) -> bool:
        """Setup steered molecular dynamics pulling"""
        
        if not self._initialized:
            raise APIError("PMF interface not initialized")
        
        self._smd_config = {
            'start_position': start_position,
            'end_position': end_position,
            'pulling_rate': pulling_rate,
            'spring_constant': spring_constant,
            'total_distance': abs(end_position - start_position),
            'pulling_time': abs(end_position - start_position) / pulling_rate
        }
        
        self.logger.info(f"Setup SMD pulling: {start_position} -> {end_position} nm")
        return True
    
    def calculate_pmf_profile(self, histogram_data: Dict[float, List[int]]) -> List[Dict[str, float]]:
        """Calculate PMF profile from histogram data"""
        
        if not self._initialized:
            raise APIError("PMF interface not initialized")
        
        # Mock PMF calculation using weighted histogram analysis method (WHAM)
        profile = []
        
        # Simple mock calculation - in real implementation this would use WHAM or similar
        for window_pos in sorted(histogram_data.keys()):
            histogram = histogram_data[window_pos]
            
            # Calculate free energy (mock implementation)
            if histogram:
                # Use histogram statistics to estimate free energy
                total_counts = sum(histogram)
                if total_counts > 0:
                    # Simple approximation: -kT * ln(probability)
                    kT = 0.593 * 298.15  # kJ/mol at 298.15 K
                    max_counts = max(histogram)
                    probability = max_counts / total_counts
                    free_energy = -kT * (probability if probability > 0 else 1e-10)
                else:
                    free_energy = 0.0
            else:
                free_energy = 0.0
            
            profile.append({
                'reaction_coordinate': window_pos,
                'pmf_energy': free_energy,
                'error_estimate': 0.5,  # Mock error estimate
                'sampling_count': sum(histogram) if histogram else 0
            })
        
        # Normalize profile (set minimum to zero)
        if profile:
            min_energy = min(p['pmf_energy'] for p in profile)
            for point in profile:
                point['pmf_energy'] -= min_energy
        
        self._pmf_profile = profile
        self.logger.info(f"Calculated PMF profile with {len(profile)} points")
        
        return profile
    
    def analyze_convergence(self, time_series_data: List[Dict]) -> Dict[str, Any]:
        """Analyze PMF calculation convergence"""
        
        if not time_series_data:
            raise DataError("No time series data provided for convergence analysis")
        
        # Mock convergence analysis
        convergence_info = {
            'is_converged': True,
            'convergence_time': len(time_series_data) * 0.8,  # Mock: 80% of simulation time
            'final_rmsd': 0.3,  # Mock RMSD from final profile
            'quality_score': 0.85,  # Mock quality score
            'recommendations': [
                "PMF profile appears well converged",
                "Consider additional sampling at window positions with high error"
            ]
        }
        
        return convergence_info
    
    def get_pmf_profile(self) -> Optional[List[Dict[str, float]]]:
        """Get calculated PMF profile"""
        return self._pmf_profile.copy() if self._pmf_profile else None
    
    def shutdown(self) -> bool:
        """Shutdown PMF interface"""
        self.logger.info("Shutting down PMF interface")
        self._initialized = False
        self._status = "shutdown"
        return True


class MDInterface(BaseInterface):
    """Interface for molecular dynamics simulation integration"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        super().__init__(logger)
        self._system_config = None
        self._simulation_params = None
        self._current_simulation = None
    
    def initialize(self, config: Dict[str, Any]) -> bool:
        """Initialize MD interface"""
        try:
            self.logger.info("Initializing MD interface")
            
            # Set MD engine
            self._md_engine = config.get('md_engine', 'gromacs')
            self._working_directory = Path(config.get('working_directory', '.'))
            
            # Validate MD engine
            supported_engines = ['gromacs', 'amber', 'namd', 'openmm']
            if self._md_engine not in supported_engines:
                raise ConfigurationError(f"Unsupported MD engine: {self._md_engine}")
            
            self._initialized = True
            self._status = "initialized"
            
            self.logger.info(f"MD interface initialized with engine: {self._md_engine}")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to initialize MD interface: {e}")
            self._status = "error"
            return False
    
    def setup_system(self, system_config: Union[SystemConfig, Dict[str, Any]]) -> bool:
        """Setup molecular system for simulation"""
        
        if not self._initialized:
            raise APIError("MD interface not initialized")
        
        # Convert dict to SystemConfig if needed
        if isinstance(system_config, dict):
            system_config = SystemConfig(**system_config)
        
        # Validate system files exist
        if not Path(system_config.coordinate_file).exists():
            raise ValidationError(f"Coordinate file not found: {system_config.coordinate_file}")
        
        if not Path(system_config.topology_file).exists():
            raise ValidationError(f"Topology file not found: {system_config.topology_file}")
        
        self._system_config = system_config
        self.logger.info(f"System setup complete: {system_config.system_name}")
        
        return True
    
    def setup_simulation(self, parameters: Union[SimulationParameters, Dict[str, Any]]) -> bool:
        """Setup simulation parameters"""
        
        if not self._initialized:
            raise APIError("MD interface not initialized")
        
        if not self._system_config:
            raise APIError("System not configured - call setup_system() first")
        
        # Convert dict to SimulationParameters if needed
        if isinstance(parameters, dict):
            parameters = SimulationParameters(**parameters)
        
        self._simulation_params = parameters
        self.logger.info(f"Simulation parameters set: {parameters.total_time} ps")
        
        return True
    
    def run_energy_minimization(self, max_iterations: int = 1000, 
                              tolerance: float = 10.0) -> Dict[str, Any]:
        """Run energy minimization"""
        
        if not self._system_config:
            raise APIError("System not configured")
        
        self.logger.info("Running energy minimization")
        
        # Mock energy minimization
        result = {
            'status': 'completed',
            'iterations': min(max_iterations, 450),  # Mock: converged early
            'final_energy': -125430.5,  # kJ/mol
            'energy_change': -2350.2,  # kJ/mol
            'converged': True,
            'output_files': {
                'minimized_structure': str(self._working_directory / 'em.gro'),
                'energy_log': str(self._working_directory / 'em.log')
            }
        }
        
        self.logger.info("Energy minimization completed successfully")
        return result
    
    def run_equilibration(self, phases: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Run equilibration phases"""
        
        if not self._system_config or not self._simulation_params:
            raise APIError("System or simulation parameters not configured")
        
        self.logger.info(f"Running equilibration with {len(phases)} phases")
        
        # Mock equilibration
        results = []
        total_time = 0
        
        for i, phase in enumerate(phases):
            phase_time = phase.get('time', 100)  # ps
            phase_name = phase.get('name', f'Phase {i+1}')
            
            phase_result = {
                'phase': phase_name,
                'time': phase_time,
                'status': 'completed',
                'final_temperature': phase.get('target_temperature', 298.15),
                'final_pressure': phase.get('target_pressure', 1.0),
                'output_files': {
                    'trajectory': str(self._working_directory / f'eq{i+1}.xtc'),
                    'final_structure': str(self._working_directory / f'eq{i+1}.gro'),
                    'log': str(self._working_directory / f'eq{i+1}.log')
                }
            }
            
            results.append(phase_result)
            total_time += phase_time
        
        equilibration_result = {
            'status': 'completed',
            'total_time': total_time,
            'phases': results,
            'system_ready': True
        }
        
        self.logger.info(f"Equilibration completed: {total_time} ps total")
        return equilibration_result
    
    def run_production_simulation(self, output_prefix: str = "prod") -> Dict[str, Any]:
        """Run production molecular dynamics simulation"""
        
        if not self._system_config or not self._simulation_params:
            raise APIError("System or simulation parameters not configured")
        
        self.logger.info("Starting production MD simulation")
        
        # Mock production simulation
        result = {
            'status': 'completed',
            'simulation_time': self._simulation_params.total_time,
            'timestep': self._simulation_params.timestep,
            'frames_written': int(self._simulation_params.total_time / 
                                (self._simulation_params.timestep * self._simulation_params.output_frequency)),
            'average_temperature': self._system_config.temperature + 0.5,  # Mock small deviation
            'average_pressure': self._system_config.pressure + 0.02,  # Mock small deviation
            'output_files': {
                'trajectory': str(self._working_directory / f'{output_prefix}.xtc'),
                'structure': str(self._working_directory / f'{output_prefix}.gro'),
                'energy': str(self._working_directory / f'{output_prefix}.edr'),
                'log': str(self._working_directory / f'{output_prefix}.log')
            }
        }
        
        self._current_simulation = result
        self.logger.info(f"Production simulation completed: {result['simulation_time']} ps")
        
        return result
    
    def get_simulation_status(self) -> Dict[str, Any]:
        """Get current simulation status"""
        
        if self._current_simulation:
            return self._current_simulation.copy()
        else:
            return {'status': 'no_simulation', 'message': 'No active simulation'}
    
    def shutdown(self) -> bool:
        """Shutdown MD interface"""
        self.logger.info("Shutting down MD interface")
        self._initialized = False
        self._status = "shutdown"
        return True


class SystemInterface(BaseInterface):
    """Interface for system-level operations and resource management"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        super().__init__(logger)
        self._resource_limits = {}
        self._active_processes = {}
    
    def initialize(self, config: Dict[str, Any]) -> bool:
        """Initialize system interface"""
        try:
            self.logger.info("Initializing system interface")
            
            # Set resource limits
            self._resource_limits = {
                'max_cpu_cores': config.get('max_cpu_cores', 8),
                'max_memory_gb': config.get('max_memory_gb', 16.0),
                'max_disk_gb': config.get('max_disk_gb', 100.0),
                'max_processes': config.get('max_processes', 10)
            }
            
            # Set system paths
            self._working_directory = Path(config.get('working_directory', '.'))
            self._temp_directory = Path(config.get('temp_directory', '/tmp'))
            
            self._initialized = True
            self._status = "initialized"
            
            self.logger.info("System interface initialized")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to initialize system interface: {e}")
            self._status = "error"
            return False
    
    def check_resources(self, required_resources: Dict[str, Any]) -> Tuple[bool, Dict[str, Any]]:
        """Check if required resources are available"""
        
        if not self._initialized:
            raise APIError("System interface not initialized")
        
        availability = {}
        all_available = True
        
        # Check CPU cores
        required_cpu = required_resources.get('cpu_cores', 1)
        max_cpu = self._resource_limits['max_cpu_cores']
        cpu_available = required_cpu <= max_cpu
        availability['cpu_cores'] = {
            'required': required_cpu,
            'available': max_cpu,
            'sufficient': cpu_available
        }
        all_available = all_available and cpu_available
        
        # Check memory
        required_memory = required_resources.get('memory_gb', 1.0)
        max_memory = self._resource_limits['max_memory_gb']
        memory_available = required_memory <= max_memory
        availability['memory_gb'] = {
            'required': required_memory,
            'available': max_memory,
            'sufficient': memory_available
        }
        all_available = all_available and memory_available
        
        # Check disk space
        required_disk = required_resources.get('disk_gb', 1.0)
        max_disk = self._resource_limits['max_disk_gb']
        disk_available = required_disk <= max_disk
        availability['disk_gb'] = {
            'required': required_disk,
            'available': max_disk,
            'sufficient': disk_available
        }
        all_available = all_available and disk_available
        
        return all_available, availability
    
    def allocate_resources(self, process_id: str, resources: Dict[str, Any]) -> bool:
        """Allocate resources for a process"""
        
        if not self._initialized:
            raise APIError("System interface not initialized")
        
        # Check if resources are available
        available, availability_info = self.check_resources(resources)
        
        if not available:
            raise ResourceError("Insufficient resources available", 
                              resource_type="mixed", 
                              required_amount=resources)
        
        # Allocate resources (mock implementation)
        self._active_processes[process_id] = {
            'allocated_resources': resources,
            'allocation_time': time.time(),
            'status': 'allocated'
        }
        
        self.logger.info(f"Resources allocated for process {process_id}: {resources}")
        return True
    
    def deallocate_resources(self, process_id: str) -> bool:
        """Deallocate resources for a process"""
        
        if not self._initialized:
            raise APIError("System interface not initialized")
        
        if process_id in self._active_processes:
            resources = self._active_processes[process_id]['allocated_resources']
            del self._active_processes[process_id]
            
            self.logger.info(f"Resources deallocated for process {process_id}: {resources}")
            return True
        else:
            self.logger.warning(f"Process {process_id} not found for deallocation")
            return False
    
    def get_system_status(self) -> Dict[str, Any]:
        """Get current system status"""
        
        return {
            'initialized': self._initialized,
            'status': self._status,
            'resource_limits': self._resource_limits.copy(),
            'active_processes': len(self._active_processes),
            'allocated_resources': {
                pid: proc['allocated_resources'] 
                for pid, proc in self._active_processes.items()
            }
        }
    
    def shutdown(self) -> bool:
        """Shutdown system interface"""
        self.logger.info("Shutting down system interface")
        
        # Deallocate all resources
        for process_id in list(self._active_processes.keys()):
            self.deallocate_resources(process_id)
        
        self._initialized = False
        self._status = "shutdown"
        return True


import time