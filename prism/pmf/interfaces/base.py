#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Module Interfaces

Defines standardized interfaces for PMF module integration with the PRISM system.
This module provides base classes and protocols for modular design and extensibility.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, Union, Protocol, Callable
from pathlib import Path
from dataclasses import dataclass
from enum import Enum
import logging

from ...utils.logging_system import PrismLogger, LogLevel, EventType


class ModulePhase(Enum):
    """Standard phases for PRISM modules"""
    INITIALIZATION = "initialization"
    VALIDATION = "validation"
    PREPARATION = "preparation"
    EXECUTION = "execution"
    ANALYSIS = "analysis"
    CLEANUP = "cleanup"
    COMPLETED = "completed"
    FAILED = "failed"


class ModuleStatus(Enum):
    """Status of module execution"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class ModuleResult:
    """Standardized result structure for PRISM modules"""
    module_name: str
    phase: ModulePhase
    status: ModuleStatus
    output_data: Dict[str, Any]
    output_files: List[Path]
    metrics: Dict[str, float]
    errors: List[str]
    warnings: List[str]
    execution_time: float
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'module_name': self.module_name,
            'phase': self.phase.value,
            'status': self.status.value,
            'output_data': self.output_data,
            'output_files': [str(f) for f in self.output_files],
            'metrics': self.metrics,
            'errors': self.errors,
            'warnings': self.warnings,
            'execution_time': self.execution_time
        }


class ModuleInterface(ABC):
    """Base interface for all PRISM modules"""
    
    def __init__(self, name: str, config: Optional[Dict[str, Any]] = None, 
                 logger: Optional[PrismLogger] = None):
        self.name = name
        self.config = config or {}
        self.logger = logger or PrismLogger(f"prism.{name}")
        self.status = ModuleStatus.PENDING
        self.current_phase = ModulePhase.INITIALIZATION
        self.results: List[ModuleResult] = []
        self._callbacks: List[Callable[[ModuleResult], None]] = []
    
    @abstractmethod
    def validate_inputs(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Validate input parameters and return validation results"""
        pass
    
    @abstractmethod
    def prepare(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare module for execution"""
        pass
    
    @abstractmethod
    def execute(self, inputs: Dict[str, Any]) -> ModuleResult:
        """Execute the main module functionality"""
        pass
    
    @abstractmethod
    def analyze_results(self, execution_result: ModuleResult) -> ModuleResult:
        """Analyze execution results"""
        pass
    
    @abstractmethod
    def cleanup(self) -> None:
        """Cleanup temporary files and resources"""
        pass
    
    def add_callback(self, callback: Callable[[ModuleResult], None]) -> None:
        """Add callback for result notifications"""
        self._callbacks.append(callback)
    
    def _notify_callbacks(self, result: ModuleResult) -> None:
        """Notify all registered callbacks"""
        for callback in self._callbacks:
            try:
                callback(result)
            except Exception as e:
                self.logger.error(f"Callback error: {e}")
    
    def run_full_workflow(self, inputs: Dict[str, Any]) -> List[ModuleResult]:
        """Run complete module workflow"""
        results = []
        
        try:
            # Phase 1: Validation
            self._set_phase(ModulePhase.VALIDATION)
            validation_result = self.validate_inputs(inputs)
            self.logger.info(f"Validation completed for {self.name}")
            
            # Phase 2: Preparation
            self._set_phase(ModulePhase.PREPARATION)
            prep_result = self.prepare(inputs)
            self.logger.info(f"Preparation completed for {self.name}")
            
            # Phase 3: Execution
            self._set_phase(ModulePhase.EXECUTION)
            exec_result = self.execute(inputs)
            results.append(exec_result)
            self._notify_callbacks(exec_result)
            
            # Phase 4: Analysis
            self._set_phase(ModulePhase.ANALYSIS)
            analysis_result = self.analyze_results(exec_result)
            results.append(analysis_result)
            self._notify_callbacks(analysis_result)
            
            # Phase 5: Cleanup
            self._set_phase(ModulePhase.CLEANUP)
            self.cleanup()
            
            # Phase 6: Completion
            self._set_phase(ModulePhase.COMPLETED)
            self.status = ModuleStatus.COMPLETED
            self.logger.info(f"Module {self.name} completed successfully")
            
        except Exception as e:
            self._set_phase(ModulePhase.FAILED)
            self.status = ModuleStatus.FAILED
            error_result = ModuleResult(
                module_name=self.name,
                phase=self.current_phase,
                status=ModuleStatus.FAILED,
                output_data={},
                output_files=[],
                metrics={},
                errors=[str(e)],
                warnings=[],
                execution_time=0.0
            )
            results.append(error_result)
            self._notify_callbacks(error_result)
            self.logger.error(f"Module {self.name} failed: {e}")
            raise
        
        self.results = results
        return results
    
    def _set_phase(self, phase: ModulePhase) -> None:
        """Set current phase and log transition"""
        self.current_phase = phase
        self.logger.workflow(f"Entering phase: {phase.value}")
    
    def get_status_summary(self) -> Dict[str, Any]:
        """Get current status summary"""
        return {
            'name': self.name,
            'status': self.status.value,
            'current_phase': self.current_phase.value,
            'results_count': len(self.results),
            'has_errors': any(r.errors for r in self.results),
            'has_warnings': any(r.warnings for r in self.results)
        }


class PMFModuleInterface(ModuleInterface):
    """Specialized interface for PMF-related modules"""
    
    @abstractmethod
    def setup_md_system(self, system_path: Path) -> Dict[str, Any]:
        """Setup MD system for PMF calculation"""
        pass
    
    @abstractmethod
    def generate_configurations(self) -> List[Path]:
        """Generate necessary configuration files"""
        pass
    
    @abstractmethod
    def calculate_pmf(self) -> Dict[str, float]:
        """Calculate PMF values"""
        pass
    
    def validate_md_system(self, system_path: Path) -> bool:
        """Validate MD system structure"""
        required_files = [
            "system.gro",
            "system.top", 
            "production.mdp"
        ]
        
        system_dir = Path(system_path)
        if not system_dir.exists():
            self.logger.error(f"System directory not found: {system_path}")
            return False
        
        missing_files = []
        for required_file in required_files:
            if not (system_dir / required_file).exists():
                missing_files.append(required_file)
        
        if missing_files:
            self.logger.warning(f"Missing files in {system_path}: {missing_files}")
            return False
        
        return True


class ModuleCommunicationProtocol(Protocol):
    """Protocol for inter-module communication"""
    
    def send_data(self, target_module: str, data: Dict[str, Any]) -> bool:
        """Send data to another module"""
        ...
    
    def receive_data(self, source_module: str) -> Optional[Dict[str, Any]]:
        """Receive data from another module"""
        ...
    
    def broadcast_status(self, status: ModuleStatus) -> None:
        """Broadcast status update to all connected modules"""
        ...


class ModuleRegistry:
    """Registry for managing PRISM modules"""
    
    def __init__(self):
        self.modules: Dict[str, ModuleInterface] = {}
        self.dependencies: Dict[str, List[str]] = {}
        self.communication_channels: Dict[str, Dict[str, Any]] = {}
        self.logger = PrismLogger("prism.module_registry")
    
    def register_module(self, module: ModuleInterface, 
                       dependencies: Optional[List[str]] = None) -> None:
        """Register a module with optional dependencies"""
        self.modules[module.name] = module
        self.dependencies[module.name] = dependencies or []
        self.communication_channels[module.name] = {}
        
        self.logger.info(f"Registered module: {module.name}")
        if dependencies:
            self.logger.info(f"Module {module.name} depends on: {dependencies}")
    
    def get_module(self, name: str) -> Optional[ModuleInterface]:
        """Get module by name"""
        return self.modules.get(name)
    
    def get_execution_order(self) -> List[str]:
        """Get modules in dependency-resolved execution order"""
        # Simple topological sort
        visited = set()
        temp_visited = set()
        result = []
        
        def visit(module_name: str):
            if module_name in temp_visited:
                raise ValueError(f"Circular dependency detected involving {module_name}")
            if module_name in visited:
                return
            
            temp_visited.add(module_name)
            for dependency in self.dependencies.get(module_name, []):
                if dependency in self.modules:
                    visit(dependency)
            temp_visited.remove(module_name)
            visited.add(module_name)
            result.append(module_name)
        
        for module_name in self.modules:
            if module_name not in visited:
                visit(module_name)
        
        return result
    
    def setup_communication_channel(self, sender: str, receiver: str, 
                                  channel_config: Dict[str, Any]) -> None:
        """Setup communication channel between modules"""
        if sender not in self.modules or receiver not in self.modules:
            raise ValueError(f"Both sender ({sender}) and receiver ({receiver}) must be registered")
        
        channel_id = f"{sender}->{receiver}"
        self.communication_channels[sender][receiver] = channel_config
        self.logger.info(f"Setup communication channel: {channel_id}")
    
    def get_module_status_summary(self) -> Dict[str, Dict[str, Any]]:
        """Get status summary for all modules"""
        return {name: module.get_status_summary() 
                for name, module in self.modules.items()}


class WorkflowOrchestrator:
    """Orchestrates execution of multiple modules in a workflow"""
    
    def __init__(self, registry: ModuleRegistry):
        self.registry = registry
        self.logger = PrismLogger("prism.workflow_orchestrator")
        self.execution_results: Dict[str, List[ModuleResult]] = {}
    
    def execute_workflow(self, workflow_config: Dict[str, Any]) -> Dict[str, List[ModuleResult]]:
        """Execute complete workflow based on configuration"""
        execution_order = self.registry.get_execution_order()
        self.logger.info(f"Executing workflow with modules: {execution_order}")
        
        results = {}
        shared_data = {}
        
        for module_name in execution_order:
            module = self.registry.get_module(module_name)
            if not module:
                self.logger.error(f"Module not found: {module_name}")
                continue
            
            self.logger.info(f"Executing module: {module_name}")
            
            try:
                # Prepare inputs with shared data from previous modules
                module_inputs = workflow_config.get(module_name, {})
                module_inputs.update(shared_data)
                
                # Execute module
                module_results = module.run_full_workflow(module_inputs)
                results[module_name] = module_results
                
                # Extract shared data for next modules
                for result in module_results:
                    if result.status == ModuleStatus.COMPLETED:
                        shared_data.update(result.output_data)
                
                self.logger.info(f"Module {module_name} completed successfully")
                
            except Exception as e:
                self.logger.error(f"Module {module_name} failed: {e}")
                # Continue with next module or stop based on configuration
                if workflow_config.get('stop_on_error', True):
                    break
        
        self.execution_results = results
        return results
    
    def get_workflow_summary(self) -> Dict[str, Any]:
        """Get summary of workflow execution"""
        total_modules = len(self.execution_results)
        successful_modules = sum(1 for results in self.execution_results.values()
                               if any(r.status == ModuleStatus.COMPLETED for r in results))
        failed_modules = total_modules - successful_modules
        
        total_execution_time = sum(
            sum(r.execution_time for r in results)
            for results in self.execution_results.values()
        )
        
        return {
            'total_modules': total_modules,
            'successful_modules': successful_modules,
            'failed_modules': failed_modules,
            'success_rate': successful_modules / total_modules if total_modules > 0 else 0,
            'total_execution_time': total_execution_time,
            'module_results': {name: [r.to_dict() for r in results] 
                             for name, results in self.execution_results.items()}
        }


# Factory function for creating PMF modules
def create_pmf_module(module_type: str, config: Optional[Dict[str, Any]] = None,
                     logger: Optional[PrismLogger] = None) -> PMFModuleInterface:
    """Factory function to create PMF modules"""
    
    if module_type == "smd":
        from ..methods.smd import SMDManager as SMDModule
        return SMDModule(config=config, logger=logger)
    elif module_type == "umbrella":
        from ..methods.umbrella import UmbrellaManager as UmbrellaModule  
        return UmbrellaModule(config=config, logger=logger)
    elif module_type == "analyzer":
        from ..analysis.analyzer import PMFAnalyzer as PMFAnalyzerModule
        return PMFAnalyzerModule(config=config, logger=logger)
    else:
        raise ValueError(f"Unknown PMF module type: {module_type}")


# Convenience function for setting up PMF workflow
def setup_pmf_workflow(config: Dict[str, Any], 
                      logger: Optional[PrismLogger] = None) -> WorkflowOrchestrator:
    """Setup complete PMF workflow with all modules"""
    
    # Create registry
    registry = ModuleRegistry()
    
    # Create and register PMF modules
    smd_module = create_pmf_module("smd", config.get("smd", {}), logger)
    umbrella_module = create_pmf_module("umbrella", config.get("umbrella", {}), logger)
    analyzer_module = create_pmf_module("analyzer", config.get("analysis", {}), logger)
    
    # Register modules with dependencies
    registry.register_module(smd_module)
    registry.register_module(umbrella_module, dependencies=["smd"])
    registry.register_module(analyzer_module, dependencies=["umbrella"])
    
    # Setup communication channels
    registry.setup_communication_channel("smd", "umbrella", {"data_type": "trajectory"})
    registry.setup_communication_channel("umbrella", "analyzer", {"data_type": "windows"})
    
    # Create orchestrator
    orchestrator = WorkflowOrchestrator(registry)
    
    return orchestrator