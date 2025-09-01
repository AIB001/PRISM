#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SMD (Steered Molecular Dynamics) module
"""

import os
import shutil
from pathlib import Path
import subprocess
import numpy as np

from .utils import create_mdp_file, extract_pull_groups


class SMDManager:
    """Manages SMD simulations"""
    
    def __init__(self, pmf_system):
        self.pmf_system = pmf_system
        self.smd_dir = pmf_system.smd_dir
        self.config = pmf_system.config['smd']
    
    def run(self, pull_groups=None, direction=None, rate=None, force=None):
        """
        Run SMD simulation
        
        Parameters:
        -----------
        pull_groups : list, optional
            Groups for pulling [reference, pulled]
        direction : list, optional
            Pull direction vector [x, y, z]
        rate : float, optional
            Pull rate in nm/ps
        force : float, optional
            Force constant in kJ/mol/nm²
        """
        # Create SMD directory
        self.smd_dir.mkdir(exist_ok=True)
        os.chdir(self.smd_dir)
        
        # Copy necessary files
        self._prepare_files()
        
        # Create index file if needed
        if not (self.smd_dir / "index.ndx").exists():
            self._create_index_file(pull_groups)
        
        # Generate SMD MDP file
        self._create_smd_mdp(direction, rate, force)
        
        # Run SMD
        self._run_smd_simulation()
        
        # Extract trajectory
        self._extract_trajectory()
        
        # Calculate distances
        self._calculate_distances()
        
        print(f"SMD completed. Results in: {self.smd_dir}")
    
    def _prepare_files(self):
        """Copy necessary files from system directory"""
        files_to_copy = ["solv_ions.gro", "topol.top"]
        for filename in files_to_copy:
            src = self.pmf_system.system_dir / filename
            dst = self.smd_dir / filename
            if src.exists() and not dst.exists():
                shutil.copy2(src, dst)
    
    def _create_index_file(self, pull_groups=None):
        """Create GROMACS index file"""
        if pull_groups is None:
            pull_groups = self.config['pull_groups']
        
        # Use gmx make_ndx to create index
        cmd = [
            "gmx", "make_ndx",
            "-f", "solv_ions.gro",
            "-o", "index.ndx"
        ]
        
        # Create groups for Protein and LIG
        input_text = "q\n"  # Quit, use default groups
        
        subprocess.run(cmd, input=input_text, text=True)
    
    def _create_smd_mdp(self, direction=None, rate=None, force=None):
        """Create SMD MDP file"""
        if direction is None:
            direction = self.config.get('pull_direction', [0, 0, 1])
        if rate is None:
            rate = self.config['pull_coord1_rate']
        if force is None:
            force = self.config['pull_coord1_k']
        
        template_path = Path(__file__).parent / "templates" / "mdp" / "pull_smd.mdp"
        
        if template_path.exists():
            # Use template
            with open(template_path, 'r') as f:
                mdp_content = f.read()
            
            # Replace placeholders
            mdp_content = mdp_content.replace("PULL_RATE", str(rate))
            mdp_content = mdp_content.replace("PULL_FORCE", str(force))
            mdp_content = mdp_content.replace("PULL_VEC", f"{direction[0]} {direction[1]} {direction[2]}")
        else:
            # Generate MDP content
            mdp_content = self._generate_smd_mdp_content(direction, rate, force)
        
        with open("smd.mdp", 'w') as f:
            f.write(mdp_content)
    
    def _run_smd_simulation(self):
        """Run the SMD simulation"""
        print("Preparing SMD simulation...")
        subprocess.run([
            "gmx", "grompp",
            "-f", "smd.mdp",
            "-c", "solv_ions.gro",
            "-p", "topol.top",
            "-n", "index.ndx",
            "-o", "smd.tpr",
            "-maxwarn", "2"
        ], check=True)
        
        print("Running SMD simulation...")
        subprocess.run([
            "gmx", "mdrun",
            "-deffnm", "smd",
            "-px", "pullx.xvg",
            "-pf", "pullf.xvg"
        ], check=True)
    
    def _extract_trajectory(self):
        """Extract configurations from trajectory"""
        print("Extracting configurations...")
        conf_dir = self.smd_dir / "conf"
        conf_dir.mkdir(exist_ok=True)
        
        # Extract frames
        subprocess.run([
            "gmx", "trjconv",
            "-s", "smd.tpr",
            "-f", "smd.xtc",
            "-o", str(conf_dir / "conf.gro"),
            "-sep"
        ], input="0\n", text=True, check=True)
    
    def _calculate_distances(self):
        """Calculate COM distances for each frame"""
        conf_dir = self.smd_dir / "conf"
        distances = []
        
        # Get list of conf files
        conf_files = sorted(conf_dir.glob("conf*.gro"))
        
        print(f"Calculating distances for {len(conf_files)} frames...")
        
        with open("distances.dat", 'w') as f:
            for i, conf_file in enumerate(conf_files):
                # Calculate distance using gmx distance
                result = subprocess.run([
                    "gmx", "distance",
                    "-s", "smd.tpr",
                    "-f", str(conf_file),
                    "-n", "index.ndx",
                    "-select", f"com of group {self.config['pull_groups'][0]} plus com of group {self.config['pull_groups'][1]}",
                    "-oall", f"dist{i}.xvg"
                ], capture_output=True, text=True)
                
                # Extract distance from output file
                if os.path.exists(f"dist{i}.xvg"):
                    with open(f"dist{i}.xvg", 'r') as dist_file:
                        for line in dist_file:
                            if not line.startswith(('#', '@')):
                                parts = line.split()
                                if len(parts) >= 2:
                                    distance = float(parts[1])
                                    distances.append((i, distance))
                                    f.write(f"{i} {distance}\n")
                                    break
                    os.remove(f"dist{i}.xvg")
        
        print(f"Distance calculation completed. Results saved to distances.dat")


class SMDBuilder:
    """Helper class to build custom SMD protocols"""
    
    @staticmethod
    def create_custom_mdp(output_file, **params):
        """Create custom SMD MDP file with specified parameters"""
        # Implementation for creating custom MDP files
        pass