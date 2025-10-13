#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF System - Main PMF workflow manager
"""

import os
import shutil
from pathlib import Path
import yaml

from .smd import SMDManager
from .umbrella import UmbrellaManager
from .analyzer import PMFAnalyzer
from ..utils import run_command


class PMFSystem:
    """
    PMF calculation system manager
    
    This class is created by PRISMSystem.build_pmf() and manages
    the entire PMF calculation workflow.
    """
    
    def __init__(self, system_dir, output_dir, config=None):
        """
        Initialize PMF system
        
        Parameters:
        -----------
        system_dir : str
            Path to GMX_PROLIG_MD directory
        output_dir : str
            Output directory for PMF calculations
        config : dict or str
            PMF configuration
        """
        self.system_dir = Path(system_dir).absolute()
        self.output_dir = Path(output_dir).absolute()
        
        # Validate system directory
        self.solv_ions_gro = self.system_dir / "solv_ions.gro"
        self.topol_top = self.system_dir / "topol.top"
        self.index_ndx = self.system_dir / "index.ndx"
        
        if not self.solv_ions_gro.exists():
            raise FileNotFoundError(
                f"System file not found: {self.solv_ions_gro}\n"
                "Please run system.build() first to create the system."
            )
        
        # Create output structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.smd_dir = self.output_dir / "smd"
        self.umbrella_dir = self.output_dir / "umbrella"
        self.analysis_dir = self.output_dir / "analysis"
        
        # Load configuration
        self.config = self._load_config(config)
        
        # Initialize managers
        self.smd_manager = SMDManager(self)
        self.umbrella_manager = UmbrellaManager(self)
        self.analyzer = PMFAnalyzer(self)
        
        # State tracking
        self.state = {
            'smd_completed': False,
            'umbrellas_setup': False,
            'umbrellas_completed': False,
            'analysis_completed': False
        }
        
        print(f"PMF system initialized")
        print(f"  System directory: {self.system_dir}")
        print(f"  Output directory: {self.output_dir}")
    
    def run_smd(self, **kwargs):
        """Run SMD simulation"""
        print("\n=== Running SMD Simulation ===")
        self.smd_manager.run(**kwargs)
        self.state['smd_completed'] = True
        return self.smd_dir
    
    def setup_umbrellas(self, interval=None, **kwargs):
        """Setup umbrella sampling windows"""
        if not self.state['smd_completed']:
            raise RuntimeError("SMD must be completed before setting up umbrellas")
        
        print("\n=== Setting Up Umbrella Windows ===")
        self.umbrella_manager.setup_windows(interval=interval, **kwargs)
        self.state['umbrellas_setup'] = True
        return self.umbrella_dir
    
    def run_umbrellas(self, **kwargs):
        """Run umbrella sampling simulations"""
        if not self.state['umbrellas_setup']:
            raise RuntimeError("Umbrellas must be setup before running")
        
        print("\n=== Running Umbrella Sampling ===")
        self.umbrella_manager.run(**kwargs)
        self.state['umbrellas_completed'] = True
        return self.umbrella_dir
    
    def analyze(self, **kwargs):
        """Analyze PMF with WHAM"""
        if not self.state['umbrellas_completed']:
            raise RuntimeError("Umbrella sampling must be completed before analysis")
        
        print("\n=== Analyzing PMF ===")
        results = self.analyzer.run_wham(**kwargs)
        self.state['analysis_completed'] = True
        return results
    
    def visualize(self, **kwargs):
        """Visualize PMF results"""
        if not self.state['analysis_completed']:
            raise RuntimeError("Analysis must be completed before visualization")
        
        print("\n=== Visualizing Results ===")
        self.analyzer.visualize(**kwargs)
    
    def run_complete_workflow(self, **kwargs):
        """Run complete PMF workflow"""
        self.run_smd()
        self.setup_umbrellas()
        self.run_umbrellas()
        self.analyze()
        self.visualize()
        print(f"\nPMF calculation completed! Results in: {self.analysis_dir}")
    
    def _load_config(self, config):
        """Load PMF configuration"""
        default_config = self._get_default_config()
        
        if isinstance(config, dict):
            return {**default_config, **config}
        elif isinstance(config, str) and Path(config).exists():
            with open(config, 'r') as f:
                user_config = yaml.safe_load(f)
            return {**default_config, **user_config}
        
        return default_config
    
    def _get_default_config(self):
        """Get default PMF configuration"""
        return {
            'smd': {
                'pull_groups': ['Protein', 'LIG'],
                'pull_coord1_dim': 'Y Y Y',
                'pull_coord1_rate': 0.005,  # nm/ps
                'pull_coord1_k': 1000,      # kJ/mol/nm²
                'nsteps': 2000000,          # 4 ns
                'dt': 0.002
            },
            'umbrella': {
                'spacing': 0.1,             # nm
                'force_constant': 1000,     # kJ/mol/nm²
                'simulation_time': 22,      # ns
                'equilibration_time': 2,    # ns
                'dt': 0.002
            },
            'analysis': {
                'temperature': 310,
                'bootstrap': 100,
                'tolerance': 1e-6
            }
        }