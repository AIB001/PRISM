#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Runner - Unified API for PMF calculations

This module provides a clean, unified API for running complete PMF workflows
with configurable parameters, replacing the need for shell scripts.
"""

import os
import yaml
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Union

from .core import PMFSystem, pmf_system
from .pmf_builder import PMFBuilder
from ..utils.config import ConfigurationManager

logger = logging.getLogger(__name__)


class PMFRunner:
    """
    Unified PMF Runner for complete workflows
    
    Provides a clean API to execute complete PMF calculations with
    configurable parameters, replacing shell scripts with Python API.
    """
    
    def __init__(self, config: Optional[Union[str, Dict]] = None):
        """
        Initialize PMF Runner
        
        Parameters:
        -----------
        config : str, dict, or None
            Configuration file path, configuration dictionary, or None for defaults
        """
        self.config = self._load_configuration(config)
        self.results = {}
        self.start_time = None
        
        # Configure logging
        log_level = self.config.get('output', {}).get('log_level', 'INFO')
        logging.getLogger().setLevel(getattr(logging, log_level))
        
        logger.info("PMF Runner initialized")
    
    def _load_configuration(self, config: Optional[Union[str, Dict]]) -> Dict:
        """Load and validate configuration"""
        if config is None:
            # Use default config
            default_config_path = Path(__file__).parent.parent / "configs" / "pmf_config.yaml"
            config = str(default_config_path)
        
        if isinstance(config, str):
            # Load from file
            with open(config, 'r') as f:
                config_dict = yaml.safe_load(f)
        else:
            # Use provided dictionary
            config_dict = config
        
        return config_dict
    
    def run_complete_workflow(self, md_system_dir: str, output_dir: str,
                            steps: Optional[List[str]] = None) -> Dict:
        """
        Run complete PMF workflow
        
        Parameters:
        -----------
        md_system_dir : str
            Path to MD system directory
        output_dir : str  
            Output directory for PMF calculations
        steps : list, optional
            List of steps to execute. If None, runs all steps:
            ['builder', 'smd', 'umbrella', 'analysis']
            
        Returns:
        --------
        Dict : Complete workflow results
        """
        if steps is None:
            steps = ['builder', 'smd', 'umbrella', 'analysis']
        
        logger.info("=== Starting Complete PMF Workflow ===")
        self.start_time = datetime.now()
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save configuration
        self._save_used_config(output_path)
        
        workflow_results = {
            'workflow_type': 'complete',
            'steps_executed': [],
            'start_time': self.start_time.isoformat(),
            'config_used': self.config,
            'output_directory': str(output_path)
        }
        
        try:
            # Step 1: System Builder (optional)
            if 'builder' in steps:
                builder_results = self._run_builder_step(md_system_dir, output_path)
                workflow_results['steps_executed'].append('builder')
                workflow_results['builder'] = builder_results
                
                # Use rebuilt system for subsequent steps
                system_dir = builder_results.get('final_system_dir', md_system_dir)
            else:
                system_dir = md_system_dir
            
            # Create PMF system for subsequent steps
            pmf = pmf_system(
                system_dir=system_dir,
                output_dir=str(output_path),
                config=self.config
            )
            
            # Step 2: SMD Simulation
            if 'smd' in steps:
                smd_results = self._run_smd_step(pmf)
                workflow_results['steps_executed'].append('smd')
                workflow_results['smd'] = smd_results
            
            # Step 3: Umbrella Sampling
            if 'umbrella' in steps:
                umbrella_results = self._run_umbrella_step(pmf)
                workflow_results['steps_executed'].append('umbrella')
                workflow_results['umbrella'] = umbrella_results
            
            # Step 4: WHAM Analysis
            if 'analysis' in steps:
                analysis_results = self._run_analysis_step(pmf)
                workflow_results['steps_executed'].append('analysis')
                workflow_results['analysis'] = analysis_results
                
                # Extract binding energy for summary
                if 'binding_energy' in analysis_results:
                    workflow_results['binding_energy'] = analysis_results['binding_energy']
            
            # Workflow completion
            end_time = datetime.now()
            workflow_results['end_time'] = end_time.isoformat()
            workflow_results['duration'] = str(end_time - self.start_time)
            workflow_results['status'] = 'completed'
            
            logger.info("=== PMF Workflow Completed Successfully ===")
            logger.info(f"Duration: {workflow_results['duration']}")
            
            if 'binding_energy' in workflow_results:
                be = workflow_results['binding_energy']
                unit = workflow_results.get('analysis', {}).get('energy_unit', 'kcal/mol')
                logger.info(f"Binding Energy: {be['value']:.2f} {unit}")
            
            return workflow_results
            
        except Exception as e:
            workflow_results['status'] = 'failed'
            workflow_results['error'] = str(e)
            workflow_results['end_time'] = datetime.now().isoformat()
            
            logger.error(f"PMF workflow failed: {e}")
            raise
    
    def _run_builder_step(self, md_system_dir: str, output_path: Path) -> Dict:
        """Execute system builder step"""
        logger.info("=== Step: System Builder ===")
        
        builder_config = self.config.get('builder', {})
        if not builder_config.get('rebuild_system', True):
            logger.info("System rebuilding disabled, skipping builder")
            return {'status': 'skipped', 'reason': 'rebuild_system = false'}
        
        # Create PMF builder
        builder_output = output_path / "rebuilt_system"
        builder = PMFBuilder(
            md_results_dir=md_system_dir,
            output_dir=builder_output,
            config=self.config
        )
        
        # Run system building with equilibration
        equilibrate = builder_config.get('run_equilibration', True)
        build_results = builder.build(equilibrate=equilibrate)
        
        logger.info("System builder completed")
        
        # Return results with final system directory
        results = {
            'status': 'completed',
            'build_results': build_results,
            'final_system_dir': build_results.get('system_dir'),
            'equilibrated': equilibrate
        }
        
        return results
    
    def _run_smd_step(self, pmf: PMFSystem) -> Dict:
        """Execute SMD simulation step"""
        logger.info("=== Step: SMD Simulation ===")
        
        smd_results = pmf.build(step='smd')
        
        logger.info("SMD preparation completed")
        logger.info(f"Manual execution: cd {smd_results.get('smd_dir', 'N/A')} && bash run_smd.sh")
        
        return {
            'status': 'prepared',
            'smd_results': smd_results,
            'manual_command': f"cd {smd_results.get('smd_dir', 'N/A')} && bash run_smd.sh"
        }
    
    def _run_umbrella_step(self, pmf: PMFSystem) -> Dict:
        """Execute umbrella sampling step"""
        logger.info("=== Step: Umbrella Sampling ===")
        
        umbrella_results = pmf.build(step='umbrella')
        
        logger.info("Umbrella sampling preparation completed")
        umbrella_dir = umbrella_results.get('umbrella_dir', 'N/A')
        logger.info(f"Manual execution: cd {umbrella_dir} && bash run_all_umbrella.sh parallel")
        
        return {
            'status': 'prepared',
            'umbrella_results': umbrella_results,
            'manual_command': f"cd {umbrella_dir} && bash run_all_umbrella.sh parallel"
        }
    
    def _run_analysis_step(self, pmf: PMFSystem) -> Dict:
        """Execute WHAM analysis step"""
        logger.info("=== Step: WHAM Analysis ===")
        
        analysis_results = pmf.build(step='analysis')
        
        logger.info("WHAM analysis completed")
        
        return {
            'status': 'completed',
            'analysis_results': analysis_results
        }
    
    def _save_used_config(self, output_path: Path):
        """Save the configuration used for this run"""
        config_file = output_path / "pmf_config_used.yaml"
        with open(config_file, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False, indent=2)
        logger.info(f"Configuration saved: {config_file}")


# =============================================================================
# High-level convenience functions
# =============================================================================

def run_pmf_workflow(md_system_dir: str, output_dir: str,
                    config: Optional[Union[str, Dict]] = None,
                    steps: Optional[List[str]] = None) -> Dict:
    """
    Run complete PMF workflow (convenience function)
    
    Parameters:
    -----------
    md_system_dir : str
        Path to MD system directory
    output_dir : str
        Output directory for PMF calculations
    config : str, dict, or None
        Configuration file path, dictionary, or None for defaults
    steps : list, optional
        Steps to execute: ['builder', 'smd', 'umbrella', 'analysis']
        
    Returns:
    --------
    Dict : Workflow results
    
    Example:
    --------
    >>> import prism.pmf as pmf
    >>> results = pmf.run_pmf_workflow("./md_system", "./pmf_results")
    >>> print(f"Binding energy: {results['binding_energy']['value']:.2f}")
    """
    runner = PMFRunner(config)
    return runner.run_complete_workflow(md_system_dir, output_dir, steps)


def run_pmf_step(md_system_dir: str, output_dir: str, step: str,
                config: Optional[Union[str, Dict]] = None) -> Dict:
    """
    Run individual PMF step
    
    Parameters:
    -----------
    md_system_dir : str
        Path to MD system directory  
    output_dir : str
        Output directory
    step : str
        Step to execute: 'builder', 'smd', 'umbrella', or 'analysis'
    config : str, dict, or None
        Configuration
        
    Returns:
    --------
    Dict : Step results
    
    Example:
    --------
    >>> import prism.pmf as pmf
    >>> results = pmf.run_pmf_step("./system", "./output", "smd")
    """
    return run_pmf_workflow(md_system_dir, output_dir, config, [step])


# =============================================================================
# Configuration utilities  
# =============================================================================

def create_pmf_config(output_file: str, template: str = "default") -> str:
    """
    Create PMF configuration file from template
    
    Parameters:
    -----------
    output_file : str
        Output configuration file path
    template : str
        Template name: "default", "fast", "accurate"
        
    Returns:
    --------
    str : Path to created configuration file
    """
    templates = {
        "default": _get_default_config(),
        "fast": _get_fast_config(), 
        "accurate": _get_accurate_config()
    }
    
    if template not in templates:
        raise ValueError(f"Unknown template: {template}")
    
    config = templates[template]
    
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, indent=2)
    
    logger.info(f"PMF configuration created: {output_file}")
    return output_file


def _get_default_config() -> Dict:
    """Get default PMF configuration"""
    return {
        'builder': {
            'rebuild_system': True,
            'z_extension': 2.5,
            'protein_forcefield': 'amber99sb',
            'ligand_forcefield': 'gaff'
        },
        'smd': {
            'pull_rate': 0.005,
            'pull_k': 1000.0,
            'nsteps': 2500000  # 5 ns
        },
        'umbrella': {
            'sample_interval_near': 0.1,
            'sample_interval_far': 0.2,
            'production_time_ps': 20000  # 20 ns per window
        }
    }


def _get_fast_config() -> Dict:
    """Get fast PMF configuration for testing"""
    config = _get_default_config()
    config['smd']['nsteps'] = 1250000  # 2.5 ns
    config['umbrella']['production_time_ps'] = 10000  # 10 ns per window
    return config


def _get_accurate_config() -> Dict:
    """Get accurate PMF configuration for production"""
    config = _get_default_config()
    config['smd']['pull_rate'] = 0.002  # Slower pulling
    config['umbrella']['sample_interval_near'] = 0.05  # Denser sampling
    config['umbrella']['production_time_ps'] = 50000  # 50 ns per window
    return config