#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Core - High-level PMF calculation interface

This module provides a simplified API for PMF calculations following PRISM architecture patterns.
"""

import os
import yaml
import logging
from pathlib import Path
from typing import Dict, Optional, Union

from .workflow import PMFWorkflow
from .pmf_builder import PMFBuilder
from ..utils.config import ConfigurationManager

logger = logging.getLogger(__name__)


class PMFSystem:
    """
    High-level interface for PMF calculations.
    
    This class provides a simplified API for calculating potential of mean force
    between protein-ligand complexes, following PRISM's architectural patterns.
    """
    
    def __init__(self, system_dir: str, output_dir: str = "pmf_output",
                 config: Optional[Union[str, Dict]] = None, rebuild_system: bool = False,
                 use_existing_files: bool = True, **kwargs):
        """
        Initialize a PMF calculation system.
        
        Parameters:
        -----------
        system_dir : str
            Path to the built system directory (containing GMX_PROLIG_MD)
        output_dir : str, optional
            Output directory for PMF results (default: "pmf_output")
        config : str or dict, optional
            PMF configuration file path or dictionary
        rebuild_system : bool, optional
            Whether to rebuild system with PMF-optimized geometry (default: False)
        use_existing_files : bool, optional
            Whether to use existing files directly without copying (default: True)
        **kwargs : optional
            Additional configuration parameters
        """
        self.system_dir = Path(system_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.rebuild_system = rebuild_system
        self.use_existing_files = use_existing_files
        
        # Validate system directory
        self._validate_system()
        
        # Process configuration
        self.config = self._process_config(config, **kwargs)
        
        # Initialize PMF builder if system rebuild is requested
        self.pmf_builder = None
        if rebuild_system:
            self.pmf_builder = PMFBuilder(
                md_results_dir=self.system_dir,
                output_dir=self.output_dir,
                use_existing_files=use_existing_files,
                config=self.config
            )
        
        # Initialize workflow manager using existing files directly
        # This ensures we work with original files without copying
        self.workflow = PMFWorkflow(
            system_dir=self.gmx_dir,
            output_dir=self.output_dir,
            config=self.config
        )
        
        logger.info(f"PMF system initialized: {self.system_dir} -> {self.output_dir}")
        if rebuild_system:
            logger.info("PMF builder enabled for system optimization")
    
    def _validate_system(self):
        """Validate system directory structure"""
        if not self.system_dir.exists():
            raise FileNotFoundError(f"System directory not found: {self.system_dir}")
        
        # Look for GMX_PROLIG_MD structure
        gmx_dir = None
        
        # Check if system_dir itself is GMX_PROLIG_MD
        if (self.system_dir / "solv_ions.gro").exists():
            gmx_dir = self.system_dir
        else:
            # Look for GMX_PROLIG_MD in subdirectories
            for candidate in self.system_dir.iterdir():
                if candidate.is_dir():
                    if (candidate / "GMX_PROLIG_MD").exists():
                        gmx_dir = candidate / "GMX_PROLIG_MD"
                        break
                    elif candidate.name == "GMX_PROLIG_MD":
                        gmx_dir = candidate
                        break
        
        if not gmx_dir:
            raise FileNotFoundError(
                f"GMX_PROLIG_MD directory not found in {self.system_dir}. "
                "Please run system.build() first."
            )
        
        self.gmx_dir = gmx_dir
        logger.debug(f"Found GMX directory: {self.gmx_dir}")
    
    def _process_config(self, config: Optional[Union[str, Dict]], **kwargs) -> Dict:
        """Process PMF configuration"""
        # Start with default configuration
        default_config = self._get_default_config()
        
        # Load from file if provided
        if isinstance(config, str) and Path(config).exists():
            with open(config, 'r') as f:
                file_config = yaml.safe_load(f)
            default_config.update(file_config)
        elif isinstance(config, dict):
            default_config.update(config)
        
        # Apply any kwargs
        if kwargs:
            # Flatten kwargs into config structure
            for key, value in kwargs.items():
                if key in ['reference_group', 'moving_group']:
                    default_config[key] = value
                elif key.startswith('smd_'):
                    default_config.setdefault('smd', {})[key[4:]] = value
                elif key.startswith('umbrella_'):
                    default_config.setdefault('umbrella', {})[key[9:]] = value
                elif key.startswith('analysis_'):
                    default_config.setdefault('analysis', {})[key[9:]] = value
                else:
                    default_config[key] = value
        
        return default_config
    
    def _get_default_config(self) -> Dict:
        """Get default PMF configuration following PRISM patterns"""
        return {
            'reference_group': 'Protein',
            'moving_group': 'LIG',
            'smd': {
                'pull_rate': 0.005,      # nm/ps
                'pull_k': 1000.0,        # kJ/mol/nm²
                'dt': 0.002,             # ps
            },
            'distance': {
                'start': 0.3,            # nm
                'end': 2.0,              # nm
            },
            'umbrella': {
                'sample_interval_near': 0.1,    # nm for close distances
                'sample_interval_far': 0.2,     # nm for far distances
                'cutoff_distance': 1.5,         # nm cutoff for adaptive sampling
                'production_time_ps': 22000,    # 22 ns per window
                'sampling_interval_ps': 10,     # Output every 10 ps
                'force_constant': 1000.0,       # kJ/mol/nm²
            },
            'simulation': {
                'temperature': 310.0,    # K
                'pressure': 1.0,         # bar
                'dt': 0.002,             # ps
            },
            'analysis': {
                'begin_time_ps': 2000,           # Skip first 2 ns
                'bootstrap_iterations': 50,      # Bootstrap samples
                'energy_unit': 'kCal',           # Output unit
            }
        }
    
    def rebuild_system(self, frame: Optional[int] = None) -> Dict:
        """
        Rebuild system with PMF-optimized geometry
        
        Parameters:
        -----------
        frame : int, optional
            Frame number to extract from MD trajectory (-1 for last)
            
        Returns:
        --------
        Dict : Rebuild results
        """
        if not self.pmf_builder:
            self.pmf_builder = PMFBuilder(
                md_results_dir=self.system_dir,
                output_dir=self.output_dir / "rebuilt_system", 
                config=self.config
            )
        
        logger.info("Rebuilding system with PMF-optimized geometry")
        rebuild_results = self.pmf_builder.build(frame)
        
        # Reinitialize workflow with rebuilt system
        rebuilt_system_dir = Path(rebuild_results['system_dir'])
        self.workflow = PMFWorkflow(
            system_dir=rebuilt_system_dir,
            output_dir=self.output_dir,
            config=self.config
        )
        
        logger.info("PMF workflow updated to use rebuilt system")
        return rebuild_results
    
    def build(self, step: Optional[str] = None) -> Dict:
        """
        Build PMF calculation (prepare all components).
        
        Parameters:
        -----------
        step : str, optional
            Specific step to prepare ('smd', 'umbrella', 'analysis')
            If None, prepares all steps
            
        Returns:
        --------
        Dict : Build results with paths and instructions
        """
        if step is None:
            # Prepare all steps
            results = {}
            results['smd'] = self.workflow.prepare_smd()
            results['umbrella'] = self.workflow.prepare_umbrella()
            results['analysis'] = self.workflow.prepare_analysis()
            
            # Provide execution instructions
            results['instructions'] = self._generate_build_instructions(results)
            
            logger.info("PMF calculation prepared successfully")
            return results
        
        elif step.lower() == 'smd':
            return self.workflow.prepare_smd()
        elif step.lower() == 'umbrella':
            return self.workflow.prepare_umbrella()
        elif step.lower() == 'analysis':
            return self.workflow.prepare_analysis()
        else:
            raise ValueError(f"Unknown step: {step}. Must be 'smd', 'umbrella', or 'analysis'")
    
    def build_step_by_step(self) -> Dict:
        """
        Build PMF calculation with step-by-step control.
        
        Returns:
        --------
        Dict : Results for each step with manual execution instructions
        """
        results = {
            'workflow': 'step_by_step',
            'steps': {}
        }
        
        # Step 1: SMD Preparation
        logger.info("=== Step 1: SMD Preparation ===")
        smd_results = self.workflow.prepare_smd()
        results['steps']['smd'] = smd_results
        
        # Provide clear instructions for manual execution
        results['instructions'] = {
            'step1': {
                'description': 'Run SMD simulation',
                'command': f"cd {smd_results['smd_dir']} && bash run_smd.sh",
                'expected_output': ['results/smd_pullf.xvg', 'results/smd_pullx.xvg'],
                'next_step': 'Call build_umbrella_step() after SMD completion'
            }
        }
        
        return results
    
    def build_umbrella_step(self) -> Dict:
        """Build umbrella sampling step (after SMD completion)"""
        logger.info("=== Step 2: Umbrella Sampling Preparation ===")
        umbrella_results = self.workflow.prepare_umbrella()
        
        return {
            'umbrella': umbrella_results,
            'instructions': {
                'step2': {
                    'description': 'Run umbrella sampling',
                    'command': f"cd {umbrella_results['umbrella_dir']} && bash run_all_umbrella.sh parallel",
                    'expected_output': 'window_*/results/umbrella_pullf.xvg files',
                    'next_step': 'Call build_analysis_step() after umbrella completion'
                }
            }
        }
    
    def build_analysis_step(self) -> Dict:
        """Build WHAM analysis step (after umbrella completion)"""
        logger.info("=== Step 3: WHAM Analysis ===")
        analysis_results = self.workflow.run_analysis()
        
        return {
            'analysis': analysis_results,
            'instructions': {
                'step3': {
                    'description': 'PMF analysis completed',
                    'results': analysis_results,
                    'binding_energy': analysis_results.get('binding_energy', {}),
                    'plots': analysis_results.get('plots', {})
                }
            }
        }
    
    def run(self, mode: str = 'auto') -> Dict:
        """
        Run PMF calculation.
        
        Parameters:
        -----------
        mode : str
            Execution mode ('auto' for automated, 'manual' for step-by-step)
            
        Returns:
        --------
        Dict : Complete PMF results
        """
        if mode.lower() == 'auto':
            return self.workflow.run_complete()
        elif mode.lower() == 'manual':
            return self.build_step_by_step()
        else:
            raise ValueError(f"Unknown mode: {mode}. Must be 'auto' or 'manual'")
    
    def get_status(self) -> Dict:
        """Get current PMF calculation status"""
        return self.workflow.get_status()
    
    def analyze_smd(self) -> Dict:
        """Analyze SMD results and generate plots"""
        return self.workflow.analyze_smd()
    
    def analyze_pmf(self) -> Dict:
        """Analyze final PMF results with detailed plots"""
        return self.workflow.analyze_pmf()
    
    def clean(self, components: Optional[list] = None) -> None:
        """
        Clean PMF calculation results.
        
        Parameters:
        -----------
        components : list, optional
            Components to clean ['smd', 'umbrella', 'analysis']
            If None, cleans all components
        """
        self.workflow.clean(components)
    
    def export_config(self, output_path: str) -> str:
        """
        Export current configuration to file.
        
        Parameters:
        -----------
        output_path : str
            Path for configuration file
            
        Returns:
        --------
        str : Path to exported configuration file
        """
        config_path = Path(output_path)
        with open(config_path, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False, indent=2)
        
        logger.info(f"Configuration exported to: {config_path}")
        return str(config_path)
    
    def _generate_build_instructions(self, results: Dict) -> Dict:
        """Generate comprehensive build instructions"""
        return {
            'workflow_type': 'complete_build',
            'execution_order': [
                {
                    'step': 1,
                    'name': 'SMD Simulation',
                    'command': f"cd {results['smd']['smd_dir']} && bash run_smd.sh",
                    'description': 'Run steered molecular dynamics simulation',
                    'estimated_time': '1-6 hours'
                },
                {
                    'step': 2,
                    'name': 'Umbrella Sampling',
                    'command': f"cd {results['umbrella']['umbrella_dir']} && bash run_all_umbrella.sh parallel",
                    'description': 'Run umbrella sampling windows',
                    'estimated_time': 'Several hours to days (depending on parallelization)'
                },
                {
                    'step': 3,
                    'name': 'WHAM Analysis',
                    'description': 'Analysis will be performed automatically after umbrella completion',
                    'note': 'Call system.build_analysis_step() to run WHAM analysis'
                }
            ],
            'output_directory': str(self.output_dir),
            'final_results': str(self.output_dir / "analysis" / "pmf_analysis_report.txt")
        }


# High-level convenience function following PRISM patterns
def pmf_system(system_dir: str, output_dir: str = "pmf_output", 
               config: Optional[Union[str, Dict]] = None, rebuild_system: bool = False,
               use_existing_files: bool = True, **kwargs) -> PMFSystem:
    """
    Create a PMF calculation system (convenience function).
    
    Parameters:
    -----------
    system_dir : str
        Path to built system directory
    output_dir : str, optional
        PMF output directory
    config : str or dict, optional
        Configuration file or dictionary
    rebuild_system : bool, optional
        Whether to rebuild system with PMF-optimized geometry
    use_existing_files : bool, optional
        Whether to use existing files directly without copying (default: True)
    **kwargs : optional
        Additional configuration parameters
        
    Returns:
    --------
    PMFSystem : PMF system instance
    
    Examples:
    --------
    >>> import prism
    >>> # Use existing files mode (recommended) - avoids file copying
    >>> pmf = prism.pmf.pmf_system("./gaff_model", "./pmf_results", use_existing_files=True)
    >>> results = pmf.run(mode='manual')
    
    >>> # Traditional mode - copies files to separate working directory
    >>> pmf = prism.pmf.pmf_system("./gaff_model", "./pmf_results", use_existing_files=False)
    >>> results = pmf.run(mode='auto')
    
    >>> # Rebuild system with PMF optimization
    >>> pmf = prism.pmf.pmf_system("./gaff_model", "./pmf_results", 
    ...                           rebuild_system=True, use_existing_files=True)
    >>> rebuild_results = pmf.rebuild_system()  # Optimize geometry
    >>> results = pmf.run(mode='auto')          # Run PMF calculation
    """
    return PMFSystem(system_dir, output_dir, config, rebuild_system, use_existing_files, **kwargs)