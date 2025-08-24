#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GROMACS simulator for PRISM MD simulations
"""

import os
import subprocess
import shutil
from pathlib import Path


class GMXSimulator:
    """
    GROMACS simulator class for running MD simulations.
    
    This class wraps GROMACS commands and uses the localrun.sh script
    for running standard MD simulation protocols.
    """
    
    def __init__(self, gmx_dir, system_files, ff_dir):
        """
        Initialize GROMACS simulator.
        
        Parameters:
        -----------
        gmx_dir : str
            Path to GMX_PROLIG_MD directory
        system_files : dict
            Dictionary of system file paths
        ff_dir : str
            Path to force field directory
        """
        self.gmx_dir = gmx_dir
        self.system_files = system_files
        self.ff_dir = ff_dir
        
        # Get localrun.sh script path
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.localrun_script = os.path.join(self.script_dir, "scripts", "localrun.sh")
        
        # Check if GROMACS is available
        self._check_gromacs()
    
    def _check_gromacs(self):
        """Check if GROMACS is available in the system"""
        try:
            result = subprocess.run(
                ["gmx", "--version"],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode != 0:
                raise RuntimeError("GROMACS not found. Please install GROMACS.")
            
            # Extract version info
            for line in result.stdout.split('\n'):
                if line.startswith("GROMACS version"):
                    print(f"Found {line}")
                    break
        except FileNotFoundError:
            raise RuntimeError("GROMACS not found. Please install GROMACS.")
    
    def run(self, stages=None, gpu_id=0, ntomp=10, ntmpi=1, continue_from=None, **kwargs):
        """
        Run GROMACS MD simulation.
        
        Parameters:
        -----------
        stages : list, optional
            List of stages to run ['em', 'nvt', 'npt', 'prod']
            If None, runs all stages
        gpu_id : int
            GPU device ID (default: 0)
        ntomp : int
            Number of OpenMP threads (default: 10)
        ntmpi : int
            Number of thread-MPI ranks (default: 1)
        continue_from : str, optional
            Continue from a specific stage ('em', 'nvt', 'npt')
        **kwargs : dict
            Additional parameters (ignored for compatibility)
            
        Returns:
        --------
        dict
            Dictionary with paths to output files
        """
        # Default stages
        if stages is None:
            stages = ['em', 'nvt', 'npt', 'prod']
        
        # Copy localrun.sh to GMX directory
        target_script = os.path.join(self.gmx_dir, "localrun.sh")
        if not os.path.exists(self.localrun_script):
            # Create localrun.sh from the content
            self._create_localrun_script(target_script)
        else:
            shutil.copy2(self.localrun_script, target_script)
        
        # Make script executable
        os.chmod(target_script, 0o755)
        
        # Modify script parameters if needed
        self._modify_script_parameters(target_script, gpu_id, ntomp, ntmpi)
        
        # Change to GMX directory
        original_dir = os.getcwd()
        os.chdir(self.gmx_dir)
        
        try:
            # Run the simulation
            print(f"Running GROMACS simulation in: {self.gmx_dir}")
            print(f"Stages: {', '.join(stages)}")
            print(f"GPU ID: {gpu_id}, Threads: {ntomp}")
            
            # Execute localrun.sh
            result = subprocess.run(
                ["bash", "localrun.sh"],
                check=True,
                text=True,
                capture_output=False  # Let output stream to console
            )
            
            # Collect output files
            outputs = self._collect_outputs()
            
            return outputs
            
        except subprocess.CalledProcessError as e:
            print(f"Error running GROMACS simulation: {e}")
            raise
        finally:
            # Return to original directory
            os.chdir(original_dir)
    
    def _create_localrun_script(self, target_path):
        """Create localrun.sh script from embedded content"""
        script_content = '''#!/bin/bash

######################################################
# SIMULATION PART
######################################################

# Energy Minimization (EM)
mkdir -p em
if [ -f ./em/em.gro ]; then
    echo "EM already completed, skipping..."
elif [ -f ./em/em.tpr ]; then
    echo "EM tpr file found, continuing from checkpoint..."
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v -cpi ./em/em.cpt
else
    echo "Starting EM from scratch..."
    gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o ./em/em.tpr -maxwarn 10
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v
fi

# NVT Equilibration
mkdir -p nvt
if [ -f ./nvt/nvt.gro ]; then
    echo "NVT already completed, skipping..."
elif [ -f ./nvt/nvt.tpr ]; then
    echo "NVT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v -cpi ./nvt/nvt.cpt
else
    echo "Starting NVT from scratch..."
    gmx grompp -f ../mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro -p topol.top -o ./nvt/nvt.tpr -maxwarn 10
    gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v
fi

# NPT Equilibration
mkdir -p npt
if [ -f ./npt/npt.gro ]; then
    echo "NPT already completed, skipping..."
elif [ -f ./npt/npt.tpr ]; then
    echo "NPT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v -cpi ./npt/npt.cpt
else
    echo "Starting NPT from scratch..."
    gmx grompp -f ../mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn 10
    gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v
fi

# Production MD
mkdir -p prod
if [ -f ./prod/md.gro ]; then
    echo "Production MD already completed, skipping..."
elif [ -f ./prod/md.tpr ]; then
    echo "Production MD tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v -cpi ./prod/md.cpt
else
    echo "Starting Production MD from scratch..."
    gmx grompp -f ../mdps/md.mdp -c ./npt/npt.gro -r ./npt/npt.gro -p topol.top -o ./prod/md.tpr -maxwarn 10
    gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v
fi'''
        
        with open(target_path, 'w') as f:
            f.write(script_content)
    
    def _modify_script_parameters(self, script_path, gpu_id, ntomp, ntmpi):
        """Modify script parameters for GPU and thread settings"""
        with open(script_path, 'r') as f:
            content = f.read()
        
        # Replace GPU ID
        content = content.replace('-gpu_id 0', f'-gpu_id {gpu_id}')
        
        # Replace thread settings
        content = content.replace('-ntomp 10', f'-ntomp {ntomp}')
        content = content.replace('-ntomp 15', f'-ntomp {ntomp}')
        content = content.replace('-ntmpi 1', f'-ntmpi {ntmpi}')
        
        with open(script_path, 'w') as f:
            f.write(content)
    
    def _collect_outputs(self):
        """Collect output files from simulation"""
        outputs = {}
        
        # Check for output files
        stage_outputs = {
            'em': {
                'gro': 'em/em.gro',
                'edr': 'em/em.edr',
                'log': 'em/em.log',
                'tpr': 'em/em.tpr'
            },
            'nvt': {
                'gro': 'nvt/nvt.gro',
                'edr': 'nvt/nvt.edr',
                'log': 'nvt/nvt.log',
                'tpr': 'nvt/nvt.tpr',
                'cpt': 'nvt/nvt.cpt'
            },
            'npt': {
                'gro': 'npt/npt.gro',
                'edr': 'npt/npt.edr',
                'log': 'npt/npt.log',
                'tpr': 'npt/npt.tpr',
                'cpt': 'npt/npt.cpt'
            },
            'prod': {
                'gro': 'prod/md.gro',
                'xtc': 'prod/md.xtc',
                'edr': 'prod/md.edr',
                'log': 'prod/md.log',
                'tpr': 'prod/md.tpr',
                'cpt': 'prod/md.cpt'
            }
        }
        
        for stage, files in stage_outputs.items():
            stage_dict = {}
            for file_type, filename in files.items():
                filepath = os.path.join(self.gmx_dir, filename)
                if os.path.exists(filepath):
                    stage_dict[file_type] = filepath
            if stage_dict:
                outputs[stage] = stage_dict
        
        return outputs