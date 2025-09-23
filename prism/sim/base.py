#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Base simulation model class for PRISM MD simulations
"""

import os
import shutil
from pathlib import Path
from .gmxsim import GMXSimulator
from .openmmsim import OpenMMSimulator
from .utils import validate_gmx_directory, find_forcefield_dir


class SimulationModel:
    """
    Base class for running MD simulations from GROMACS-prepared systems.
    
    This class provides a unified interface for running simulations using
    either GROMACS or OpenMM engines.
    """
    
    def __init__(self, gmx_dir):
        """
        Initialize simulation model.
        
        Parameters:
        -----------
        gmx_dir : str
            Path to the GMX_PROLIG_MD directory containing the prepared system
        """
        self.gmx_dir = os.path.abspath(gmx_dir)
        # The parent directory contains GMX_PROLIG_MD, mdps, and LIG.xxx2gmx
        self.output_dir = os.path.dirname(self.gmx_dir)
        
        # Validate directory structure
        self._validate_directory()
        
        # Find associated force field directory (at same level as GMX_PROLIG_MD)
        self.ff_dir = find_forcefield_dir(self.output_dir)
        
        # Get MDP directory (at same level as GMX_PROLIG_MD)
        self.mdp_dir = os.path.join(self.output_dir, "mdps")
        if not os.path.exists(self.mdp_dir):
            raise FileNotFoundError(
                f"MDP directory not found at {self.mdp_dir}. "
                f"Expected directory structure:\n"
                f"  {self.output_dir}/\n"
                f"    ├── GMX_PROLIG_MD/\n"
                f"    ├── mdps/\n"
                f"    └── LIG.xxx2gmx/"
            )
        
        # System files
        self.system_files = {
            'gro': os.path.join(self.gmx_dir, 'solv_ions.gro'),
            'top': os.path.join(self.gmx_dir, 'topol.top'),
            'em_mdp': os.path.join(self.mdp_dir, 'em.mdp'),
            'nvt_mdp': os.path.join(self.mdp_dir, 'nvt.mdp'),
            'npt_mdp': os.path.join(self.mdp_dir, 'npt.mdp'),
            'md_mdp': os.path.join(self.mdp_dir, 'md.mdp')
        }
        
        # Available engines
        self.engines = {
            'gmx': GMXSimulator,
            'gromacs': GMXSimulator,
            'openmm': OpenMMSimulator
        }
        
        print(f"Simulation model initialized:")
        print(f"  GMX directory: {self.gmx_dir}")
        print(f"  Force field directory: {self.ff_dir}")
        print(f"  MDP directory: {self.mdp_dir}")
    
    def _validate_directory(self):
        """Validate that the directory contains necessary GROMACS files"""
        if not os.path.exists(self.gmx_dir):
            raise FileNotFoundError(f"GMX directory not found: {self.gmx_dir}")
        
        required_files = ['solv_ions.gro', 'topol.top']
        missing = []
        
        for file in required_files:
            if not os.path.exists(os.path.join(self.gmx_dir, file)):
                missing.append(file)
        
        if missing:
            raise FileNotFoundError(
                f"Missing required files in {self.gmx_dir}: {', '.join(missing)}"
            )
    
    def run(self, engine='gmx', **kwargs):
        """
        Run MD simulation using specified engine.
        
        Parameters:
        -----------
        engine : str
            Simulation engine to use ('gmx' or 'openmm')
        **kwargs : dict
            Additional parameters for the simulation engine
            
        Returns:
        --------
        dict
            Dictionary containing paths to output files
            
        Examples:
        ---------
        >>> sim = SimulationModel("output/GMX_PROLIG_MD")
        >>> results = sim.run(engine="gmx")
        >>> results = sim.run(engine="openmm", platform="CUDA")
        """
        engine_name = engine.lower()
        
        if engine_name not in self.engines:
            raise ValueError(
                f"Unknown engine: {engine}. Available: {list(self.engines.keys())}"
            )
        
        print(f"\n{'='*60}")
        print(f"Running MD simulation with {engine_name.upper()}")
        print(f"{'='*60}")
        
        # Create simulator instance
        simulator_class = self.engines[engine_name]
        simulator = simulator_class(self.gmx_dir, self.system_files, self.ff_dir)
        
        # Run simulation
        results = simulator.run(**kwargs)
        
        print(f"\n{'='*60}")
        print(f"Simulation completed successfully!")
        print(f"{'='*60}")
        
        return results
    
    def info(self):
        """Print information about the simulation system"""
        print(f"\nSimulation System Information:")
        print(f"  Output directory: {self.output_dir}")
        print(f"  GMX directory: {os.path.basename(self.gmx_dir)}")
        print(f"  Force field directory: {os.path.basename(self.ff_dir)}")
        print(f"  MDP directory: {os.path.basename(self.mdp_dir)}")
        
        print(f"\n  Directory structure:")
        print(f"    {self.output_dir}/")
        print(f"      ├── {os.path.basename(self.gmx_dir)}/")
        print(f"      ├── {os.path.basename(self.mdp_dir)}/")
        print(f"      └── {os.path.basename(self.ff_dir)}/")
        
        print(f"\n  System files:")
        for name, path in self.system_files.items():
            exists = "✓" if os.path.exists(path) else "✗"
            if 'mdp' in name:
                rel_path = os.path.join(os.path.basename(self.mdp_dir), os.path.basename(path))
            else:
                rel_path = os.path.join(os.path.basename(self.gmx_dir), os.path.basename(path))
            print(f"    {exists} {name}: {rel_path}")
        
        # Check for existing simulation outputs
        print(f"\n  Existing simulation outputs in {os.path.basename(self.gmx_dir)}:")
        
        # Check subdirectories
        stages_output = {
            "em": ["em/em.gro", "em/em.edr", "em/em.log"],
            "nvt": ["nvt/nvt.gro", "nvt/nvt.cpt", "nvt/nvt.log"],
            "npt": ["npt/npt.gro", "npt/npt.cpt", "npt/npt.log"],
            "prod": ["prod/md.gro", "prod/md.xtc", "prod/md.cpt", "prod/md.log"]
        }
        
        for stage, files in stages_output.items():
            stage_exists = False
            for file in files:
                if os.path.exists(os.path.join(self.gmx_dir, file)):
                    stage_exists = True
                    break
            status = "✓" if stage_exists else "✗"
            print(f"    {status} {stage.upper()} stage")
    
    def __repr__(self):
        """String representation"""
        return f"SimulationModel(gmx_dir='{self.gmx_dir}')"