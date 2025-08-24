#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Umbrella sampling module
"""

import os
import shutil
from pathlib import Path
import subprocess
import numpy as np


class UmbrellaManager:
    """Manages umbrella sampling simulations"""
    
    def __init__(self, pmf_system):
        self.pmf_system = pmf_system
        self.umbrella_dir = pmf_system.umbrella_dir
        self.config = pmf_system.config['umbrella']
        self.windows = []
    
    def setup_windows(self, interval=None, distance_file=None):
        """
        Setup umbrella sampling windows
        
        Based on setupUmbrella.py functionality
        """
        if interval is None:
            interval = self.config['spacing']
        
        if distance_file is None:
            distance_file = self.pmf_system.smd_dir / "distances.dat"
        
        # Create umbrella directory
        self.umbrella_dir.mkdir(exist_ok=True)
        
        # Read distances
        distances = self._read_distances(distance_file)
        
        # Select windows
        self.windows = self._select_windows(distances, interval)
        
        # Setup each window
        self._setup_window_directories()
        
        # Generate submission scripts
        self._generate_scripts()
        
        print(f"Setup {len(self.windows)} umbrella windows")
        return self.windows
    
    def run(self, parallel=False, max_parallel=10):
        """Run umbrella sampling simulations"""
        if parallel:
            self._run_parallel(max_parallel)
        else:
            self._run_sequential()
    
    def _read_distances(self, distance_file):
        """Read distance file"""
        distances = []
        with open(distance_file, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.split()
                    if len(parts) >= 2:
                        frame = int(parts[0])
                        dist = float(parts[1])
                        distances.append((frame, dist))
        return distances
    
    def _select_windows(self, distances, interval):
        """Select windows for umbrella sampling"""
        windows = []
        if not distances:
            return windows
        
        # Start with first frame
        windows.append(distances[0])
        current_dist = distances[0][1]
        
        # Select frames at interval spacing
        for frame, dist in distances[1:]:
            if abs(dist - current_dist) >= interval:
                windows.append((frame, dist))
                current_dist = dist
        
        return windows
    
    def _setup_window_directories(self):
        """Setup directories for each window"""
        for i, (frame, dist) in enumerate(self.windows):
            window_dir = self.umbrella_dir / f"window_{i:03d}"
            window_dir.mkdir(exist_ok=True)
            
            # Copy necessary files
            conf_file = self.pmf_system.smd_dir / "conf" / f"conf{frame}.gro"
            if conf_file.exists():
                shutil.copy2(conf_file, window_dir / "conf.gro")
            
            # Copy topology and index
            for filename in ["topol.top", "index.ndx"]:
                src = self.pmf_system.smd_dir / filename
                if src.exists():
                    shutil.copy2(src, window_dir / filename)
            
            # Create umbrella MDP
            self._create_umbrella_mdp(window_dir, dist)
    
    def _create_umbrella_mdp(self, window_dir, distance):
        """Create umbrella sampling MDP for a window"""
        template_path = Path(__file__).parent / "templates" / "mdp" / "umbrella.mdp"
        
        if template_path.exists():
            with open(template_path, 'r') as f:
                mdp_content = f.read()
            
            # Replace placeholders
            mdp_content = mdp_content.replace("UMBRELLA_POS", str(distance))
            mdp_content = mdp_content.replace("FORCE_CONSTANT", str(self.config['force_constant']))
            mdp_content = mdp_content.replace("NSTEPS", str(int(self.config['simulation_time'] * 1000 / self.config['dt'])))
        else:
            # Generate content
            mdp_content = self._generate_umbrella_mdp_content(distance)
        
        with open(window_dir / "umbrella.mdp", 'w') as f:
            f.write(mdp_content)
    
    def _run_sequential(self):
        """Run umbrella windows sequentially"""
        for i, (frame, dist) in enumerate(self.windows):
            window_dir = self.umbrella_dir / f"window_{i:03d}"
            os.chdir(window_dir)
            
            print(f"Running window {i+1}/{len(self.windows)} (distance: {dist:.3f} nm)")
            
            # grompp
            subprocess.run([
                "gmx", "grompp",
                "-f", "umbrella.mdp",
                "-c", "conf.gro",
                "-p", "topol.top",
                "-n", "index.ndx",
                "-o", "umbrella.tpr",
                "-maxwarn", "2"
            ], check=True)
            
            # mdrun
            subprocess.run([
                "gmx", "mdrun",
                "-deffnm", "umbrella",
                "-px", "pullx.xvg",
                "-pf", "pullf.xvg"
            ], check=True)


class WindowSelector:
    """Advanced window selection strategies"""
    
    @staticmethod
    def select_by_force(forces, threshold=100):
        """Select windows based on force threshold"""
        # Implementation
        pass
    
    @staticmethod
    def select_adaptive(distances, min_spacing=0.05, max_spacing=0.2):
        """Adaptive window selection based on PMF gradient"""
        # Implementation
        pass