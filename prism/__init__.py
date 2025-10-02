#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM - Potential of Mean Force (PMF) Calculation System

A comprehensive molecular dynamics PMF calculation framework with advanced
workflow management, analysis tools, and system integration capabilities.

This package provides:
- Core PMF calculation modules (smd, umbrella sampling, workflow)
- Advanced error handling and logging systems
- Configuration management and validation
- Data management with multi-format support
- CLI interface and API access
- Workflow orchestration and monitoring
- Performance optimization and testing framework
"""

# Version and metadata
__version__ = "2.0.0"
__author__ = "PRISM Development Team"
__email__ = "prism@research.org"
__license__ = "MIT"
__title__ = "PRISM PMF System"
__description__ = "Advanced Potential of Mean Force Calculation Framework"
__url__ = "https://github.com/prism/pmf-system"

# System requirements
import sys
from pathlib import Path

MINIMUM_PYTHON = (3, 7)
if sys.version_info < MINIMUM_PYTHON:
    raise RuntimeError(
        f"PRISM requires Python {MINIMUM_PYTHON[0]}.{MINIMUM_PYTHON[1]} or higher. "
        f"Current version: {sys.version_info.major}.{sys.version_info.minor}"
    )

# System paths and directories
PRISM_ROOT = Path(__file__).parent
PRISM_CONFIG_DIR = Path.home() / '.prism'
PRISM_DATA_DIR = Path.cwd() / 'prism_data'
PRISM_LOG_DIR = PRISM_CONFIG_DIR / 'logs'
PRISM_CACHE_DIR = PRISM_CONFIG_DIR / 'cache'
PRISM_PLUGINS_DIR = PRISM_CONFIG_DIR / 'plugins'

# Create necessary directories
for directory in [PRISM_CONFIG_DIR, PRISM_DATA_DIR, PRISM_LOG_DIR, PRISM_CACHE_DIR, PRISM_PLUGINS_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Core module imports with graceful error handling
_system_logger = None

try:
    # Initialize core utilities first
    from .utils.logging_system import PrismLogger
    from .utils.config_validator import ConfigurationValidator, ValidationLevel
    from .utils.error_handler import ErrorHandler, ErrorSeverity
    
    # Initialize system logger
    _system_logger = PrismLogger("prism_system")
    _system_logger.info(f"PRISM PMF System v{__version__} initializing...")
    
except ImportError as e:
    # Graceful degradation if core utilities are not available
    import logging
    logging.basicConfig(level=logging.INFO)
    _system_logger = logging.getLogger("prism_system")
    _system_logger.warning(f"Core utilities not fully available: {e}")

# Import existing PRISM components
try:
    from .builder import PRISMBuilder
    from .core import PRISMSystem  
    from .sim import model
    from .analysis import TrajAnalysis
    from . import modeling
    from . import pmf
    
    _has_legacy_components = True
    if _system_logger:
        _system_logger.info("Legacy PRISM components loaded")
        
except ImportError as e:
    _has_legacy_components = False
    if _system_logger:
        _system_logger.warning(f"Legacy components not available: {e}")

# Import new optimized systems
_available_systems = {
    'api': False,
    'cli': False,
    'workflow': False,
    'monitoring': False,
    'testing': False,
    'data_management': False
}

# API System
try:
    from .api import PMFCalculator, WorkflowManager, DataAnalyzer, create_client
    from .api.data import FormatConverter, DataExporter, DataImporter
    _available_systems['api'] = True
    if _system_logger:
        _system_logger.info("API system loaded")
except ImportError as e:
    if _system_logger:
        _system_logger.warning(f"API system not available: {e}")

# CLI System  
try:
    from .cli.main import PrismCLI
    _available_systems['cli'] = True
    if _system_logger:
        _system_logger.info("CLI system loaded")
except ImportError as e:
    if _system_logger:
        _system_logger.warning(f"CLI system not available: {e}")

# Workflow System
try:
    from .workflow.orchestrator import WorkflowOrchestrator
    from .workflow.scheduler import TaskScheduler
    _available_systems['workflow'] = True
    if _system_logger:
        _system_logger.info("Workflow system loaded")
except ImportError as e:
    if _system_logger:
        _system_logger.warning(f"Workflow system not available: {e}")

# Monitoring System
try:
    from .monitoring.performance import PerformanceMonitor
    from .monitoring.alerts import AlertManager
    _available_systems['monitoring'] = True
    if _system_logger:
        _system_logger.info("Monitoring system loaded")
except ImportError as e:
    if _system_logger:
        _system_logger.warning(f"Monitoring system not available: {e}")

# Testing Framework
try:
    from .testing.framework import PrismTestRunner
    _available_systems['testing'] = True
    if _system_logger:
        _system_logger.info("Testing framework loaded")
except ImportError as e:
    if _system_logger:
        _system_logger.warning(f"Testing framework not available: {e}")

# Data Management System
try:
    from .utils.data_management import DataManager
    from .utils.storage_manager import StorageManager
    _available_systems['data_management'] = True
    if _system_logger:
        _system_logger.info("Data management system loaded")
except ImportError as e:
    if _system_logger:
        _system_logger.warning(f"Data management system not available: {e}")

# High-level convenience functions for PRISM integration
def get_system_info():
    """
    Get comprehensive PRISM system information and available components.
    
    Returns:
    --------
    dict
        Complete system information including available modules and dependencies
    """
    info = {
        "version": __version__,
        "title": __title__,
        "description": __description__,
        "author": __author__,
        "license": __license__,
        "url": __url__,
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "system_paths": {
            "prism_root": str(PRISM_ROOT),
            "config_dir": str(PRISM_CONFIG_DIR),
            "data_dir": str(PRISM_DATA_DIR),
            "log_dir": str(PRISM_LOG_DIR),
            "cache_dir": str(PRISM_CACHE_DIR),
            "plugins_dir": str(PRISM_PLUGINS_DIR)
        },
        "legacy_components": _has_legacy_components,
        "available_systems": _available_systems.copy(),
        "dependencies": check_dependencies()
    }
    
    # Add component counts
    available_count = sum(1 for available in _available_systems.values() if available)
    info["optimization_systems"] = {
        "total": len(_available_systems),
        "available": available_count,
        "coverage": f"{available_count/len(_available_systems)*100:.1f}%"
    }
    
    return info

def get_available_apis():
    """
    Get list of available API interfaces and their status.
    
    Returns:
    --------
    dict
        Dictionary of available APIs and their import status
    """
    apis = {}
    
    # Check API system
    if _available_systems.get('api', False):
        try:
            from .api import PMFCalculator, WorkflowManager, DataAnalyzer
            apis['library_api'] = {
                'available': True,
                'classes': ['PMFCalculator', 'WorkflowManager', 'DataAnalyzer'],
                'description': 'Python Library API for direct access'
            }
        except ImportError:
            apis['library_api'] = {'available': False, 'error': 'Import failed'}
        
        try:
            from .api.client import PrismAPIClient
            apis['client_api'] = {
                'available': True,
                'classes': ['PrismAPIClient'],
                'description': 'Client API for external access'
            }
        except ImportError:
            apis['client_api'] = {'available': False, 'error': 'Import failed'}
            
        try:
            from .api.data import FormatConverter, DataExporter, DataImporter
            apis['data_api'] = {
                'available': True,
                'classes': ['FormatConverter', 'DataExporter', 'DataImporter'],
                'description': 'Data exchange and format conversion'
            }
        except ImportError:
            apis['data_api'] = {'available': False, 'error': 'Import failed'}
    
    # Check CLI system
    if _available_systems.get('cli', False):
        apis['cli_api'] = {
            'available': True,
            'classes': ['PrismCLI'],
            'description': 'Command-line interface and management tools'
        }
    
    return apis

def initialize_system(force_reload=False, verbose=False):
    """
    Initialize or reinitialize the PRISM system with all available components.
    
    Parameters:
    -----------
    force_reload : bool, optional
        Force reload of all modules (default: False)
    verbose : bool, optional
        Print detailed initialization information (default: False)
        
    Returns:
    --------
    dict
        Initialization results and status
    """
    global _system_logger, _available_systems, _has_legacy_components
    
    if verbose:
        print(f"🚀 Initializing PRISM PMF System v{__version__}")
        print(f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
        print(f"Root directory: {PRISM_ROOT}")
    
    initialization_results = {
        "version": __version__,
        "success": True,
        "systems_loaded": [],
        "systems_failed": [],
        "legacy_components": _has_legacy_components,
        "total_systems": len(_available_systems)
    }
    
    # Re-attempt loading systems if force_reload
    if force_reload:
        _available_systems = {key: False for key in _available_systems}
    
    # Check each optimization system
    for system_name in _available_systems:
        if not _available_systems[system_name] or force_reload:
            try:
                if system_name == 'api':
                    from .api import PMFCalculator
                    _available_systems['api'] = True
                    if verbose: print(f"  ✅ API system loaded")
                    
                elif system_name == 'cli':
                    from .cli.main import PrismCLI
                    _available_systems['cli'] = True
                    if verbose: print(f"  ✅ CLI system loaded")
                    
                elif system_name == 'workflow':
                    from .workflow.orchestrator import WorkflowOrchestrator
                    _available_systems['workflow'] = True
                    if verbose: print(f"  ✅ Workflow system loaded")
                    
                elif system_name == 'monitoring':
                    from .monitoring.performance import PerformanceMonitor
                    _available_systems['monitoring'] = True
                    if verbose: print(f"  ✅ Monitoring system loaded")
                    
                elif system_name == 'testing':
                    from .testing.framework import PrismTestRunner
                    _available_systems['testing'] = True
                    if verbose: print(f"  ✅ Testing framework loaded")
                    
                elif system_name == 'data_management':
                    from .utils.data_management import DataManager
                    _available_systems['data_management'] = True
                    if verbose: print(f"  ✅ Data management loaded")
                    
                initialization_results["systems_loaded"].append(system_name)
                
            except ImportError as e:
                initialization_results["systems_failed"].append({
                    "system": system_name,
                    "error": str(e)
                })
                if verbose: print(f"  ❌ {system_name} system failed: {e}")
    
    # Update system logger
    if _system_logger:
        loaded_count = len(initialization_results["systems_loaded"])
        total_count = len(_available_systems)
        _system_logger.info(f"System initialization complete: {loaded_count}/{total_count} systems loaded")
    
    # Calculate success metrics
    loaded_systems = sum(1 for available in _available_systems.values() if available)
    initialization_results["systems_available"] = loaded_systems
    initialization_results["success_rate"] = f"{loaded_systems/len(_available_systems)*100:.1f}%"
    
    if verbose:
        print(f"\n📊 Initialization Summary:")
        print(f"  Systems loaded: {loaded_systems}/{len(_available_systems)}")
        print(f"  Success rate: {initialization_results['success_rate']}")
        print(f"  Legacy components: {'✅' if _has_legacy_components else '❌'}")
    
    return initialization_results

def create_manager():
    """
    Create a PRISM manager that provides unified access to all available functionality.
    
    The manager automatically detects available systems and provides a consistent
    interface for accessing PMF calculations, data analysis, workflow management,
    and other PRISM features.
    
    Returns:
    --------
    PrismManager
        Manager instance with methods for all available systems
        
    Examples:
    ---------
    >>> import prism
    >>> manager = prism.create_manager()
    >>> status = manager.get_status()
    >>> result = manager.run_pmf_calculation(
    ...     system_file="complex.gro",
    ...     ligand="LIG",
    ...     protein="Protein"
    ... )
    >>> analysis = manager.analyze_data(result)
    """
    class PrismManager:
        """PRISM system unified manager
        
        Provides a single interface to access all PRISM functionality including
        PMF calculations, data analysis, workflow management, and system monitoring.
        Automatically adapts to available system components and provides graceful
        degradation when optional components are missing.
        """
        
        def __init__(self):
            self.version = __version__
            self.available_systems = _available_systems.copy()
            self._initialize_components()
        
        def _initialize_components(self):
            """Initialize available components based on system detection"""
            # API components
            if self.available_systems.get('api', False):
                try:
                    from .api import PMFCalculator, WorkflowManager, DataAnalyzer, create_client
                    self.pmf_calculator = PMFCalculator()
                    self.workflow_manager = WorkflowManager() 
                    self.data_analyzer = DataAnalyzer()
                    self.api_client = create_client()
                except ImportError:
                    pass
            
            # CLI components
            if self.available_systems.get('cli', False):
                try:
                    from .cli.main import PrismCLI
                    self.cli = PrismCLI()
                except ImportError:
                    pass
            
            # Workflow components
            if self.available_systems.get('workflow', False):
                try:
                    from .workflow.orchestrator import WorkflowOrchestrator
                    self.orchestrator = WorkflowOrchestrator()
                except ImportError:
                    pass
            
            # Data management components
            if self.available_systems.get('data_management', False):
                try:
                    from .utils.data_management import DataManager
                    self.data_manager = DataManager()
                except ImportError:
                    pass
            
            # Legacy components (backward compatibility)
            if _has_legacy_components:
                self.system = system
                self.build_system = build_system
                self.analyze_trajectory = analyze_trajectory
        
        def get_status(self):
            """Get manager and system status
            
            Returns:
            --------
            dict
                Dictionary containing version, available systems, and capabilities
            """
            return {
                "version": self.version,
                "available_systems": self.available_systems,
                "legacy_support": _has_legacy_components,
                "capabilities": self.get_capabilities()
            }
        
        def get_capabilities(self):
            """Get available capabilities
            
            Returns:
            --------
            dict
                Dictionary of available capabilities and their status
            """
            capabilities = {
                "pmf_calculation": hasattr(self, 'pmf_calculator'),
                "data_analysis": hasattr(self, 'data_analyzer'),
                "workflow_management": hasattr(self, 'workflow_manager'),
                "cli_interface": hasattr(self, 'cli'),
                "data_management": hasattr(self, 'data_manager'),
                "legacy_functions": _has_legacy_components
            }
            return capabilities
        
        def run_pmf_calculation(self, *args, **kwargs):
            """Run PMF calculation using the best available system
            
            Automatically selects between optimized API or legacy functions.
            
            Parameters:
            -----------
            *args, **kwargs
                Arguments passed to the PMF calculation function
                
            Returns:
            --------
            PMF calculation results
            """
            if hasattr(self, 'pmf_calculator'):
                return self.pmf_calculator.run_pmf_calculation(*args, **kwargs)
            elif _has_legacy_components:
                return run_pmf_workflow(*args, **kwargs)
            else:
                raise RuntimeError(
                    "No PMF calculation system available. "
                    "Please install required dependencies or check system configuration."
                )
        
        def analyze_data(self, *args, **kwargs):
            """Analyze data using the best available system
            
            Automatically selects between optimized API or legacy functions.
            
            Parameters:
            -----------
            *args, **kwargs
                Arguments passed to the data analysis function
                
            Returns:
            --------
            Data analysis results
            """
            if hasattr(self, 'data_analyzer'):
                return self.data_analyzer.analyze_pmf_profile(*args, **kwargs)
            elif _has_legacy_components:
                return analyze_trajectory(*args, **kwargs)
            else:
                raise RuntimeError(
                    "No data analysis system available. "
                    "Please install required dependencies or check system configuration."
                )
        
        def create_workflow(self, *args, **kwargs):
            """Create and manage workflows
            
            Parameters:
            -----------
            *args, **kwargs
                Arguments for workflow creation
                
            Returns:
            --------
            Workflow ID or workflow object
            """
            if hasattr(self, 'workflow_manager'):
                return self.workflow_manager.create_workflow(*args, **kwargs)
            elif hasattr(self, 'orchestrator'):
                return self.orchestrator.create_workflow(*args, **kwargs)
            else:
                raise RuntimeError("No workflow management system available")
        
        def manage_data(self, operation, *args, **kwargs):
            """Manage data storage and retrieval
            
            Parameters:
            -----------
            operation : str
                Data operation to perform
            *args, **kwargs
                Arguments for the data operation
                
            Returns:
            --------
            Operation results
            """
            if hasattr(self, 'data_manager'):
                return getattr(self.data_manager, operation)(*args, **kwargs)
            else:
                raise RuntimeError("No data management system available")
        
        def execute_cli_command(self, *args, **kwargs):
            """Execute CLI commands programmatically
            
            Parameters:
            -----------
            *args, **kwargs
                CLI command arguments
                
            Returns:
            --------
            Command execution results
            """
            if hasattr(self, 'cli'):
                return self.cli.execute_command(*args, **kwargs)
            else:
                raise RuntimeError("No CLI system available")
        
        def __repr__(self):
            """String representation of the manager"""
            available_count = sum(1 for available in self.available_systems.values() if available)
            total_count = len(self.available_systems)
            return f"PrismManager(version={self.version}, systems={available_count}/{total_count})"
    
    return PrismManager()

def create_unified_client():
    """
    Create a unified client (deprecated).
    
    .. deprecated:: 2.0.0
        Use :func:`create_manager` instead.
        
    Returns:
    --------
    PrismManager
        Manager instance (same as create_manager())
    """
    import warnings
    warnings.warn(
        "create_unified_client() is deprecated and will be removed in v3.0. "
        "Please use create_manager() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return create_manager()

# Legacy high-level API functions
def system(protein_path, ligand_path, config=None, **kwargs):
    """
    Create a protein-ligand system for MD simulation.

    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)
    config : str or dict, optional
        Configuration file path or configuration dictionary
    **kwargs : optional
        Additional parameters (output_dir, ligand_forcefield, forcefield, etc.)

    Returns:
    --------
    PRISMSystem
        A system object with methods to build and run simulations

    Examples:
    ---------
    >>> import prism as pm
    >>> system = pm.system("protein.pdb", "ligand.mol2")
    >>> system.build()
    >>> system.generate_mdp_files()

    >>> # With custom configuration
    >>> system = pm.system("protein.pdb", "ligand.sdf",
    ...                    config="config.yaml",
    ...                    ligand_forcefield="openff")
    >>> system.build()
    """
    return PRISMSystem(protein_path, ligand_path, config=config, **kwargs)

def build_system(protein_path, ligand_path, output_dir="prism_output", **kwargs):
    """
    Build a complete protein-ligand system (one-step function).

    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)
    output_dir : str
        Output directory for generated files
    **kwargs : optional
        Additional parameters (forcefield, water_model, ligand_forcefield, etc.)

    Returns:
    --------
    str
        Path to the output directory

    Examples:
    ---------
    >>> import prism as pm
    >>> output_path = pm.build_system("protein.pdb", "ligand.mol2")

    >>> # With custom force fields
    >>> output_path = pm.build_system("protein.pdb", "ligand.mol2",
    ...                               forcefield="amber99sb-ildn",
    ...                               water_model="tip4p",
    ...                               ligand_forcefield="openff")
    """
    system_obj = PRISMSystem(protein_path, ligand_path, output_dir=output_dir, **kwargs)
    return system_obj.build()

def visualize_trajectory(trajectory, topology, ligand, output="contact_analysis.html", **kwargs):
    """
    Generate interactive HTML visualization of protein-ligand contacts from MD trajectory.
    
    This function creates an interactive HTML file that visualizes the contact
    patterns between protein and ligand throughout the trajectory.
    
    Parameters:
    -----------
    trajectory : str
        Path to trajectory file (.xtc, .dcd, .trr, etc.)
    topology : str
        Path to topology/structure file (.pdb, .gro, etc.)
    ligand : str
        Path to ligand structure file (.sdf, .mol, .mol2)
    output : str, optional
        Output HTML file path (default: "contact_analysis.html")
    **kwargs : optional
        Additional parameters for visualization
        
    Returns:
    --------
    str
        Path to the generated HTML file
        
    Examples:
    ---------
    >>> import prism as pm
    >>> # Basic visualization
    >>> pm.visualize_trajectory("md.xtc", "system.gro", "ligand.sdf")
    
    >>> # Custom output file
    >>> pm.visualize_trajectory("trajectory.xtc", "topology.pdb", "ligand.mol2",
    ...                        output="my_analysis.html")
    
    Notes:
    ------
    The generated HTML file is self-contained and can be opened in any modern
    web browser. It includes interactive features like:
    - Drag to reposition residues
    - Zoom in/out
    - Toggle between 2D and 3D views
    - Export high-resolution images
    """
    try:
        from .analysis.visualization import generate_html
    except ImportError as e:
        raise ImportError(
            "Visualization module not available. "
            "Please ensure mdtraj is installed: pip install mdtraj"
        ) from e
    
    return generate_html(trajectory, topology, ligand, output, **kwargs)

def analyze_trajectory(topology, trajectory, ligand_resname="LIG", output_dir="analysis_results"):
    """
    Perform comprehensive trajectory analysis with visualization.
    
    This is a convenience function that combines trajectory analysis with
    HTML visualization generation.
    
    Parameters:
    -----------
    topology : str
        Path to topology/structure file (.pdb, .gro, etc.)
    trajectory : str
        Path to trajectory file (.xtc, .dcd, .trr, etc.)
    ligand_resname : str, optional
        Residue name of ligand (default: "LIG")
    output_dir : str, optional
        Output directory for analysis results (default: "analysis_results")
        
    Returns:
    --------
    TrajAnalysis
        Analysis object with results
        
    Examples:
    ---------
    >>> import prism as pm
    >>> analysis = pm.analyze_trajectory("system.gro", "md.xtc")
    >>> # Results are saved in analysis_results/ directory
    
    >>> # Custom ligand residue name
    >>> analysis = pm.analyze_trajectory("complex.pdb", "trajectory.dcd",
    ...                                  ligand_resname="MOL")
    """
    traj = TrajAnalysis(topology, trajectory, ligand_resname=ligand_resname)
    traj.analyze_all(output_dir=output_dir)
    return traj

def check_dependencies():
    """
    Check if all required dependencies are available.

    Returns:
    --------
    dict
        Dictionary showing availability of each dependency

    Example:
    --------
    >>> import prism as pm
    >>> deps = pm.check_dependencies()
    >>> print(deps)
    {'gromacs': True, 'pdbfixer': True, 'antechamber': True, 'openff': False, 'mdtraj': True}
    """
    from .utils.environment import GromacsEnvironment
    import subprocess

    dependencies = {
        'gromacs': False,
        'pdbfixer': False,
        'antechamber': False,
        'openff': False,
        'mdtraj': False,
        'rdkit': False
    }

    # Check GROMACS
    try:
        env = GromacsEnvironment()
        dependencies['gromacs'] = True
    except:
        pass

    # Check pdbfixer
    try:
        import pdbfixer
        dependencies['pdbfixer'] = True
    except:
        pass

    # Check for AmberTools (antechamber)
    try:
        subprocess.run(['antechamber', '-h'], capture_output=True, check=True)
        dependencies['antechamber'] = True
    except:
        pass

    # Check OpenFF
    try:
        import openff.toolkit
        dependencies['openff'] = True
    except:
        pass
    
    # Check MDTraj (for visualization)
    try:
        import mdtraj
        dependencies['mdtraj'] = True
    except:
        pass
    
    # Check RDKit (optional for enhanced visualization)
    try:
        import rdkit
        dependencies['rdkit'] = True
    except:
        pass

    return dependencies

def list_forcefields():
    """
    List available force fields from GROMACS installation.

    Returns:
    --------
    list
        List of available force field names

    Example:
    --------
    >>> import prism as pm
    >>> ffs = pm.list_forcefields()
    >>> for ff in ffs:
    ...     print(ff)
    amber99sb (index: 1)
    amber99sb-ildn (index: 2)
    amber14sb (index: 3)
    """
    from .utils.environment import GromacsEnvironment

    try:
        env = GromacsEnvironment()
        return env.list_force_fields()
    except Exception as e:
        print(f"Error detecting force fields: {e}")
        return []

def get_version():
    """Get PRISM version."""
    return __version__

def get_html_generator():
    """
    Get HTMLGenerator class for advanced visualization usage.
    
    Returns:
    --------
    HTMLGenerator
        The HTML generator class
        
    Example:
    --------
    >>> import prism as pm
    >>> HTMLGen = pm.get_html_generator()
    >>> generator = HTMLGen("trajectory.xtc", "topology.gro", "ligand.sdf")
    >>> generator.analyze()
    >>> generator.generate("output.html")
    """
    try:
        from .analysis.visualization import HTMLGenerator
        return HTMLGenerator
    except ImportError as e:
        raise ImportError(
            "HTMLGenerator not available. "
            "Please ensure mdtraj is installed: pip install mdtraj"
        ) from e

# PMF calculation functions
def run_pmf_workflow(md_system_dir, output_dir, config=None, **kwargs):
    """
    Run complete PMF calculation workflow.
    
    This is a high-level API function that automatically executes the complete
    workflow from MD results to PMF calculation, including SMD, umbrella sampling,
    and WHAM analysis.
    
    Parameters:
    -----------
    md_system_dir : str
        Path to MD results directory containing GMX_PROLIG_MD subdirectory
    output_dir : str
        Output directory for PMF calculation results
    config : str or dict, optional
        PMF configuration file path or configuration dictionary
    **kwargs : optional
        Additional configuration parameters
        
    Returns:
    --------
    dict
        Dictionary containing PMF calculation results including binding energy
        
    Examples:
    ---------
    >>> import prism as pm
    >>> # Basic usage
    >>> results = pm.run_pmf_workflow("./md_results", "./pmf_output")
    >>> print(f"Binding energy: {results['binding_energy']['value']:.2f} kcal/mol")
    
    >>> # Using custom configuration
    >>> config = {
    ...     'smd': {'pull_rate': 0.01, 'nsteps': 1000000},
    ...     'umbrella': {'production_time_ps': 15000}
    ... }
    >>> results = pm.run_pmf_workflow("./md_results", "./pmf_output", config=config)
    """
    return pmf.run_pmf_workflow(md_system_dir, output_dir, config, **kwargs)

def pmf_system(system_dir, output_dir, config=None, **kwargs):
    """
    Create PMF system object for manual control of PMF calculations.
    
    This function creates a PMF system object that allows users to control
    the PMF calculation process step by step, suitable for scenarios requiring
    fine-grained control or debugging.
    
    Parameters:
    -----------
    system_dir : str
        System directory path containing MD results or system to be calculated
    output_dir : str
        Output directory for PMF calculation results
    config : str or dict, optional
        PMF configuration file path or configuration dictionary
    **kwargs : optional
        Additional configuration parameters
        
    Returns:
    --------
    PMFSystem
        PMF system object with step-by-step execution methods
        
    Examples:
    ---------
    >>> import prism as pm
    >>> # Create PMF system
    >>> pmf_sys = pm.pmf_system("./gaff_model", "./pmf_results")
    >>> 
    >>> # Step-by-step execution
    >>> smd_results = pmf_sys.build(step='smd')
    >>> # Manual SMD run: bash ./gaff_model/pmf_smd/run_smd.sh
    >>> 
    >>> umbrella_results = pmf_sys.build_umbrella_step()
    >>> # Manual umbrella run: bash ./gaff_model/pmf_umbrella/run_all_umbrella.sh
    >>> 
    >>> analysis_results = pmf_sys.build_analysis_step()
    >>> print(f"Binding energy: {analysis_results['binding_energy']['value']:.2f} kcal/mol")
    
    >>> # Or automated execution
    >>> results = pmf_sys.run(mode='auto')
    """
    return pmf.pmf_system(system_dir, output_dir, config, **kwargs)

def run_pmf_step(md_system_dir, output_dir, step, config=None, **kwargs):
    """
    Run a single step of PMF calculation.
    
    Allows users to execute a specific step of PMF calculation independently,
    suitable for phased execution or error recovery scenarios.
    
    Parameters:
    -----------
    md_system_dir : str
        MD results directory path
    output_dir : str
        Output directory path
    step : str
        Step to execute: 'builder', 'smd', 'umbrella', 'analysis'
    config : str or dict, optional
        Configuration file path or configuration dictionary
    **kwargs : optional
        Additional configuration parameters
        
    Returns:
    --------
    dict
        Execution results for the specified step
        
    Examples:
    ---------
    >>> import prism as pm
    >>> # Prepare SMD only
    >>> smd_results = pm.run_pmf_step("./md_results", "./pmf_output", "smd")
    >>> 
    >>> # Run analysis only
    >>> analysis_results = pm.run_pmf_step("./md_results", "./pmf_output", "analysis")
    """
    return pmf.run_pmf_step(md_system_dir, output_dir, step, config, **kwargs)

def create_pmf_config(output_file, template='default'):
    """
    Create PMF configuration file template.
    
    Parameters:
    -----------
    output_file : str
        Output configuration file path
    template : str, optional
        Configuration template type: 'default', 'fast', 'accurate'
        
    Examples:
    ---------
    >>> import prism as pm
    >>> pm.create_pmf_config("my_pmf_config.yaml", template="accurate")
    """
    return pmf.create_pmf_config(output_file, template)

def get_pmf_info():
    """
    Get PMF module information.
    
    Returns:
    --------
    dict
        Detailed information about the PMF module
        
    Examples:
    ---------
    >>> import prism as pm
    >>> info = pm.get_pmf_info()
    >>> print(f"PMF module version: {info['version']}")
    >>> print(f"Supported forcefields: {info['supported_forcefields']}")
    """
    return pmf.get_pmf_info()

# System initialization and optimization systems check
def _log_system_status():
    """Log system initialization status"""
    if _system_logger:
        total_systems = len(_available_systems)
        available_systems = sum(1 for available in _available_systems.values() if available)
        
        _system_logger.info(f"PRISM PMF System v{__version__} ready")
        _system_logger.info(f"Optimization systems: {available_systems}/{total_systems} available")
        
        if _has_legacy_components:
            _system_logger.info("Legacy PRISM components: Available")
        
        for system_name, available in _available_systems.items():
            status = "✓" if available else "✗"
            _system_logger.debug(f"  {status} {system_name.replace('_', ' ').title()} System")

# Auto-initialize system on import
try:
    _log_system_status()
except Exception as e:
    if _system_logger:
        _system_logger.warning(f"System status logging failed: {e}")

# Export main classes and functions
__all__ = [
    # Core legacy classes
    "PRISMBuilder",
    "PRISMSystem", 
    "TrajAnalysis",
    
    # High-level API functions (legacy)
    "system",
    "build_system", 
    "model",
    
    # Analysis functions (legacy)
    "analyze_trajectory",
    "visualize_trajectory",
    "get_html_generator",
    
    # Utility functions
    "check_dependencies",
    "list_forcefields",
    "get_version",
    
    # Module imports
    "modeling",
    "pmf",
    
    # PMF module APIs (legacy)
    "run_pmf_workflow",
    "pmf_system",
    "run_pmf_step", 
    "create_pmf_config",
    "get_pmf_info",
    
    # New unified system functions
    "get_system_info",
    "get_available_apis",
    "initialize_system",
    "create_manager",
    "create_unified_client",  # Deprecated, use create_manager
    
    # Version and metadata
    "__version__",
    "__author__",
    "__license__",
    "__title__",
    "__description__",
    "__url__",
]

# Conditional imports for optimization systems
if _available_systems.get('api', False):
    try:
        # Make key API classes available at package level
        from .api import PMFCalculator, WorkflowManager, DataAnalyzer
        from .api import create_client as create_api_client
        from .api.data import FormatConverter, DataExporter, DataImporter
        
        __all__.extend([
            "PMFCalculator",
            "WorkflowManager", 
            "DataAnalyzer",
            "create_api_client",
            "FormatConverter",
            "DataExporter",
            "DataImporter"
        ])
        
        if _system_logger:
            _system_logger.debug("API classes imported at package level")
            
    except ImportError as e:
        if _system_logger:
            _system_logger.warning(f"Failed to import API classes at package level: {e}")

if _available_systems.get('cli', False):
    try:
        from .cli.main import PrismCLI
        __all__.append("PrismCLI")
        
        if _system_logger:
            _system_logger.debug("CLI classes imported at package level")
            
    except ImportError as e:
        if _system_logger:
            _system_logger.warning(f"Failed to import CLI classes at package level: {e}")

if _available_systems.get('workflow', False):
    try:
        from .workflow.orchestrator import WorkflowOrchestrator
        from .workflow.scheduler import TaskScheduler  
        __all__.extend(["WorkflowOrchestrator", "TaskScheduler"])
        
        if _system_logger:
            _system_logger.debug("Workflow classes imported at package level")
            
    except ImportError as e:
        if _system_logger:
            _system_logger.warning(f"Failed to import workflow classes at package level: {e}")

if _available_systems.get('monitoring', False):
    try:
        from .monitoring.performance import PerformanceMonitor
        from .monitoring.alerts import AlertManager
        __all__.extend(["PerformanceMonitor", "AlertManager"])
        
        if _system_logger:
            _system_logger.debug("Monitoring classes imported at package level")
            
    except ImportError as e:
        if _system_logger:
            _system_logger.warning(f"Failed to import monitoring classes at package level: {e}")

if _available_systems.get('testing', False):
    try:
        from .testing.framework import PrismTestRunner
        __all__.append("PrismTestRunner")
        
        if _system_logger:
            _system_logger.debug("Testing classes imported at package level")
            
    except ImportError as e:
        if _system_logger:
            _system_logger.warning(f"Failed to import testing classes at package level: {e}")

if _available_systems.get('data_management', False):
    try:
        from .utils.data_management import DataManager
        from .utils.storage_manager import StorageManager
        __all__.extend(["DataManager", "StorageManager"])
        
        if _system_logger:
            _system_logger.debug("Data management classes imported at package level")
            
    except ImportError as e:
        if _system_logger:
            _system_logger.warning(f"Failed to import data management classes at package level: {e}")

# Convenience function for quick system check
def quick_system_check():
    """
    Quick system status check for debugging and verification.
    
    Prints a formatted summary of system status to console.
    """
    print(f"\n🔬 PRISM PMF System v{__version__}")
    print(f"📁 Root: {PRISM_ROOT}")
    print(f"🐍 Python: {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
    
    print("\n📦 Optimization Systems:")
    for system_name, available in _available_systems.items():
        status = "✅" if available else "❌"
        print(f"  {status} {system_name.replace('_', ' ').title()}")
    
    print(f"\n🏗️ Legacy Components: {'✅ Available' if _has_legacy_components else '❌ Not available'}")
    
    available_count = sum(1 for available in _available_systems.values() if available)
    print(f"\n📊 System Status: {available_count}/{len(_available_systems)} optimization systems available")
    
    if available_count == len(_available_systems):
        print("🎉 All systems operational!")
    elif available_count > 0:
        print("✅ Partial functionality available")
    else:
        print("⚠️  Limited functionality - check installation")

__all__.append("quick_system_check")

# Module-level metadata for introspection
import datetime
_module_info = {
    "initialization_timestamp": datetime.datetime.now().isoformat(),
    "optimization_systems_loaded": list(_available_systems.keys()),
    "optimization_systems_available": [k for k, v in _available_systems.items() if v],
    "legacy_components_available": _has_legacy_components,
    "total_exported_symbols": len(__all__),
    "system_directories": {
        "config": str(PRISM_CONFIG_DIR),
        "data": str(PRISM_DATA_DIR), 
        "logs": str(PRISM_LOG_DIR),
        "cache": str(PRISM_CACHE_DIR),
        "plugins": str(PRISM_PLUGINS_DIR)
    }
}

# Make module info available  
__all__.append("_module_info")

# Final system initialization log
if _system_logger:
    available_count = sum(1 for available in _available_systems.values() if available)
    total_count = len(_available_systems)
    _system_logger.info(f"PRISM initialization complete - {available_count}/{total_count} optimization systems available")
    
    if available_count == total_count:
        _system_logger.info("🎉 All optimization systems successfully loaded!")
    elif available_count > total_count // 2:
        _system_logger.info(f"✅ Most optimization systems loaded ({available_count}/{total_count})")
    else:
        _system_logger.warning(f"⚠️  Limited optimization systems available ({available_count}/{total_count})")
    
    _system_logger.info(f"PRISM package initialization completed successfully at {_module_info['initialization_timestamp']}")
    _system_logger.debug(f"Exported {len(__all__)} symbols to package namespace")

# End of PRISM package initialization