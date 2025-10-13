#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced Umbrella sampling module compatible with step-by-step API
"""

import os
import shutil
from pathlib import Path
import subprocess
import numpy as np
import logging

from ..utils.exceptions import (
    UmbrellaIncompleteError, RequiredFileNotFoundError, PrerequisiteNotMetError,
    SMDExecutionError, PMFErrorCode
)
from ..utils.error_handling import (
    with_retry, error_context, safe_file_operation, validate_prerequisites,
    ErrorCollector, handle_external_tool_error
)

logger = logging.getLogger(__name__)


class UmbrellaManager:
    """Enhanced umbrella sampling simulation manager"""
    
    def __init__(self, pmf_system):
        self.pmf_system = pmf_system
        self.umbrella_dir = pmf_system.umbrella_dir
        self.config = pmf_system.config
        self.windows = []
    
    def setup_windows(self, interval=None, distance_file=None, **kwargs):
        """
        Setup umbrella sampling windows (enhanced with adaptive sampling)
        
        Parameters:
        -----------
        interval : float, optional
            Window spacing in nm
        distance_file : str, optional
            Path to distance file from SMD
            
        Raises:
        -------
        UmbrellaIncompleteError
            If umbrella window setup fails
        PrerequisiteNotMetError
            If SMD prerequisites are not met
        """
        with error_context("umbrella windows setup", {"interval": interval}):
            logger.info("=== Setting Up Umbrella Windows ===")
            
            # Validate prerequisites
            self._validate_umbrella_prerequisites()
            
            # Create umbrella directory safely
            self._create_umbrella_directory()
        
        # Extract trajectory frames if needed
        logger.info("Extracting trajectory frames...")
        frames_info = self._extract_trajectory_frames()
        
        # Calculate distances and generate windows
        logger.info("Generating umbrella windows...")
        self.windows = self._generate_umbrella_windows(frames_info)
        
        # Create umbrella directory structure and MDP files
        logger.info("Setting up umbrella directories...")
        self._setup_umbrella_directories()
        
        # Generate umbrella run scripts
        logger.info("Generating umbrella run scripts...")
        run_scripts = self._generate_umbrella_run_scripts()
        
        # Clean up temporary trajectory frames
        logger.info("Cleaning up temporary frames...")
        frames_dir = self.pmf_system.smd_dir / "trajectory_frames"
        if frames_dir.exists():
            shutil.rmtree(frames_dir)
            logger.info("Temporary frames cleaned up")
        
        results = {
            'umbrella_dir': str(self.umbrella_dir),
            'windows': [{'id': w['window_id'], 'distance': w['distance']} for w in self.windows],
            'n_windows': len(self.windows),
            'run_script': str(run_scripts['master_script']),
            'individual_scripts': run_scripts['individual_scripts'],
            'status': 'ready_for_umbrella'
        }
        
        logger.info("Umbrella sampling setup completed!")
        logger.info(f"Generated {len(self.windows)} umbrella windows")
        logger.info(f"Run umbrella sampling using: {run_scripts['master_script']}")
        
        return results
    
    def run(self, parallel=False, max_parallel=10):
        """Run umbrella sampling simulations (automated execution)"""
        if parallel:
            self._run_parallel(max_parallel)
        else:
            self._run_sequential()
    
    def _validate_umbrella_prerequisites(self) -> None:
        """Validate all prerequisites for umbrella sampling setup"""
        
        requirements = {
            'smd_completed': lambda: self._check_smd_completion(),
            'umbrella_directory_writable': lambda: self._check_directory_writable(self.umbrella_dir),
            'sufficient_disk_space': lambda: self._check_disk_space(self.umbrella_dir, required_gb=10.0)
        }
        
        validate_prerequisites(requirements, 'umbrella sampling setup')
    
    def _check_smd_completion(self) -> bool:
        """Check that SMD simulation has been completed with enhanced validation"""
        smd_dir = self.pmf_system.smd_dir
        
        required_files = [
            smd_dir / "results" / "smd.xtc", 
            smd_dir / "results" / "smd_pullf.xvg",
            smd_dir / "results" / "smd_pullx.xvg",
            smd_dir / "results" / "smd.tpr"
        ]
        
        missing_files = [str(f) for f in required_files if not f.exists()]
        
        if missing_files:
            raise RequiredFileNotFoundError(
                file_path=str(smd_dir / "results"),
                step="SMD simulation completion check",
                alternative_files=missing_files
            )
        
        # Validate file integrity
        for file_path in required_files:
            if file_path.stat().st_size == 0:
                from ..utils.exceptions import FileCorruptionError
                raise FileCorruptionError(
                    file_path=str(file_path),
                    expected_format="SMD output",
                    actual_issue="File is empty"
                )
        
        logger.info("SMD results validated")
        return True
    
    def _check_directory_writable(self, directory: Path) -> bool:
        """Check if directory is writable"""
        try:
            directory.mkdir(parents=True, exist_ok=True)
            test_file = directory / ".write_test"
            test_file.write_text("test")
            test_file.unlink()
            return True
        except (PermissionError, OSError):
            return False
    
    def _check_disk_space(self, directory: Path, required_gb: float) -> bool:
        """Check if sufficient disk space is available"""
        try:
            import shutil
            free_bytes = shutil.disk_usage(directory).free
            free_gb = free_bytes / (1024**3)
            
            if free_gb < required_gb:
                from ..utils.exceptions import InsufficientDiskSpaceError
                raise InsufficientDiskSpaceError(
                    required_space_gb=required_gb,
                    available_space_gb=free_gb,
                    operation="umbrella sampling setup"
                )
            
            return True
        except Exception as exc:
            logger.warning(f"Could not check disk space: {exc}")
            return True
    
    def _create_umbrella_directory(self) -> None:
        """Create umbrella directory safely"""
        try:
            self.umbrella_dir.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Created umbrella directory: {self.umbrella_dir}")
        except Exception as exc:
            raise UmbrellaIncompleteError(
                completed_windows=0,
                total_windows=0,
                failed_windows=[str(self.umbrella_dir)]
            ) from exc
    
    def _extract_trajectory_frames(self, n_frames=300):
        """Extract frames from SMD trajectory"""
        smd_dir = self.pmf_system.smd_dir
        logger.info(f"Extracting {n_frames} trajectory frames...")
        
        frames_dir = smd_dir / "trajectory_frames"
        frames_dir.mkdir(exist_ok=True)
        
        # Calculate frame extraction parameters
        distance_start = self.config.get('distance', {}).get('start', 0.3)
        distance_end = self.config.get('distance', {}).get('end', 2.0)
        pull_rate = self.config.get('smd', {}).get('pull_rate', 0.005)
        
        distance_range = distance_end - distance_start
        smd_time_ps = distance_range / pull_rate
        frame_interval_ps = smd_time_ps / n_frames
        
        # Extract frames using trjconv
        trjconv_cmd = [
            "gmx", "trjconv",
            "-f", str(smd_dir / "results" / "smd.xtc"),
            "-s", str(smd_dir / "results" / "smd.tpr"),
            "-o", str(frames_dir / "frame.gro"),
            "-sep",
            "-dt", str(frame_interval_ps)
        ]
        
        try:
            process = subprocess.Popen(
                trjconv_cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(smd_dir)
            )
            
            stdout, stderr = process.communicate(input="0\n")  # Select System
            
            if process.returncode != 0:
                handle_external_tool_error(
                    command=" ".join(trjconv_cmd),
                    exit_code=process.returncode,
                    stderr=stderr,
                    tool_name="GROMACS trjconv"
                )
            
            # Count extracted frames
            extracted_frames = len(list(frames_dir.glob("frame*.gro")))
            logger.info(f"Extracted {extracted_frames} trajectory frames")
            
            return {
                'frames_dir': str(frames_dir),
                'n_frames': extracted_frames,
                'frame_interval_ps': frame_interval_ps
            }
            
        except Exception as e:
            logger.error(f"Frame extraction failed: {e}")
            raise
    
    def _run_sequential(self):
        """Run umbrella windows sequentially"""
        for i, window in enumerate(self.windows):
            window_dir = window['directory']
            os.chdir(window_dir)
            
            logger.info(f"Running window {i+1}/{len(self.windows)} (distance: {window['distance']:.3f} nm)")
            
            # grompp
            subprocess.run([
                "gmx", "grompp",
                "-f", "umbrella.mdp",
                "-c", "start.gro", 
                "-p", "../topol.top",
                "-n", "../index.ndx",
                "-o", "umbrella.tpr",
                "-maxwarn", "10"
            ], check=True)
            
            # mdrun
            subprocess.run([
                "gmx", "mdrun",
                "-s", "umbrella.tpr",
                "-deffnm", "umbrella",
                "-ntmpi", "1",
                "-ntomp", "10",
                "-v"
            ], check=True)
            
            # Move results
            results_dir = window_dir / "results"
            results_dir.mkdir(exist_ok=True)
            for result_file in window_dir.glob("umbrella.*"):
                shutil.move(result_file, results_dir / result_file.name)
    
    def _run_parallel(self, max_parallel=4):
        """Run umbrella windows in parallel (placeholder)"""
        logger.warning("Parallel execution not implemented. Use sequential mode or manual scripts.")
        self._run_sequential()
    
    def _generate_umbrella_windows(self, frames_info):
        """Generate umbrella sampling windows based on trajectory frames"""
        frames_dir = Path(frames_info['frames_dir'])
        smd_dir = self.pmf_system.smd_dir
        
        # Calculate distances for all frames
        logger.info("Calculating distances for frames...")
        distance_data = self._calculate_frame_distances(smd_dir, frames_dir)
        
        # Generate adaptive windows
        sample_interval_near = self.config.get('umbrella', {}).get('sample_interval_near', 0.1)
        sample_interval_far = self.config.get('umbrella', {}).get('sample_interval_far', 0.2)
        cutoff_distance = self.config.get('umbrella', {}).get('cutoff_distance', 1.5)
        
        # Adaptive sampling logic
        windows = []
        selected_indices = []
        
        if not distance_data:
            raise UmbrellaIncompleteError(
                completed_windows=0,
                total_windows=0,
                failed_windows=["No distance data available for window generation"]
            )
        
        current_idx = 0
        selected_indices.append(current_idx)
        
        distances = [d[1] for d in distance_data]
        
        while current_idx < len(distances) - 1:
            current_distance = distances[current_idx]
            
            # Choose interval based on distance
            target_interval = sample_interval_near if current_distance < cutoff_distance else sample_interval_far
            target_distance = current_distance + target_interval
            
            # Find next frame closest to target distance
            best_idx = current_idx
            best_diff = float('inf')
            
            for idx in range(current_idx + 1, len(distances)):
                diff = abs(distances[idx] - target_distance)
                if diff < best_diff:
                    best_diff = diff
                    best_idx = idx
            
            if best_idx > current_idx:
                selected_indices.append(best_idx)
                current_idx = best_idx
            else:
                break
        
        # Create window information
        for i, idx in enumerate(selected_indices):
            frame, distance = distance_data[idx]
            windows.append({
                'window_id': i,
                'frame': frame,
                'distance': distance,
                'frame_file': frames_dir / f"frame{frame}.gro"
            })
        
        logger.info(f"Generated {len(windows)} umbrella windows")
        return windows
    
    def _calculate_frame_distances(self, smd_dir, frames_dir):
        """Calculate distances for trajectory frames"""
        reference_group = self.config.get('reference_group', 'Protein')
        moving_group = self.config.get('moving_group', 'LIG')
        
        frame_files = sorted(frames_dir.glob("frame*.gro"))
        distance_data = []
        
        temp_dir = smd_dir / "temp_distance"
        temp_dir.mkdir(exist_ok=True)
        
        try:
            for frame_file in frame_files:
                # Extract frame number from filename
                frame_num = int(frame_file.stem.replace('frame', ''))
                
                # Calculate distance using gmx distance
                distance_cmd = [
                    "gmx", "distance",
                    "-s", str(smd_dir / "results" / "smd.tpr"),
                    "-f", str(frame_file),
                    "-n", str(smd_dir / "index.ndx"),
                    "-select", f"com of group {reference_group} plus com of group {moving_group}",
                    "-oall", str(temp_dir / "dist_temp.xvg")
                ]
                
                try:
                    subprocess.run(distance_cmd, check=True, capture_output=True, text=True)
                    
                    # Read distance value
                    dist_file = temp_dir / "dist_temp.xvg"
                    with open(dist_file, 'r') as f:
                        lines = f.readlines()
                        for line in reversed(lines):
                            if not line.startswith('#') and not line.startswith('@') and line.strip():
                                distance = float(line.split()[1])
                                distance_data.append((frame_num, distance))
                                break
                
                except subprocess.CalledProcessError as e:
                    logger.warning(f"Could not calculate distance for frame {frame_num}: {e}")
                    continue
        
        finally:
            # Clean up temporary directory
            shutil.rmtree(temp_dir, ignore_errors=True)
        
        logger.info(f"Calculated distances for {len(distance_data)} frames")
        return distance_data
    
    def _setup_umbrella_directories(self):
        """Setup umbrella sampling directory structure"""
        # Generate umbrella MDP template
        umbrella_mdp_content = self._generate_umbrella_mdp_content()
        
        for window in self.windows:
            window_dir = self.umbrella_dir / f"window_{window['window_id']:03d}"
            window_dir.mkdir(exist_ok=True)
            
            # Create window-specific MDP file
            window_mdp_content = umbrella_mdp_content.replace(
                "DISTANCE_PLACEHOLDER", f"{window['distance']:.6f}"
            ).replace(
                "WINDOW_TITLE", f"Umbrella sampling at {window['distance']:.3f} nm (window {window['window_id']:03d})"
            )
            
            window_mdp = window_dir / "umbrella.mdp"
            with open(window_mdp, 'w') as f:
                f.write(window_mdp_content)
            
            # Copy starting structure
            if window['frame_file'].exists():
                dest_file = window_dir / "start.gro"
                shutil.copy2(window['frame_file'], dest_file)
            
            # Update window info with directory
            window['directory'] = window_dir
        
        logger.info(f"Setup {len(self.windows)} umbrella directories")
    
    def _generate_umbrella_mdp_content(self):
        """Generate umbrella sampling MDP template content"""
        dt = self.config.get('simulation', {}).get('dt', 0.002)
        production_time = self.config.get('umbrella', {}).get('production_time_ps', 22000)
        sampling_interval = self.config.get('umbrella', {}).get('sampling_interval_ps', 10)
        reference_group = self.config.get('reference_group', 'Protein')
        moving_group = self.config.get('moving_group', 'LIG')
        pull_k = self.config.get('smd', {}).get('pull_k', 1000.0)
        temp = self.config.get('simulation', {}).get('temperature', 310.0)
        pressure = self.config.get('simulation', {}).get('pressure', 1.0)
        
        nsteps = int(production_time / dt)
        output_interval = int(sampling_interval / dt)
        
        return f"""; WINDOW_TITLE
title               = WINDOW_TITLE
define              = -DPOSRE
integrator          = md
nsteps              = {nsteps}
dt                  = {dt}

; Output control
nstxtcout           = {output_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {output_interval}

; Bonds and constraints
continuation        = yes
constraint_algorithm = lincs
constraints         = h-bonds
lincs_iter          = 1
lincs_order         = 4

; Neighbor searching and short-range nonbonded interactions
cutoff-scheme       = Verlet
nstlist             = 10
rlist               = 1.0
rcoulomb            = 1.0
rvdw                = 1.0

; Electrostatics
coulombtype         = PME
pme_order           = 4
fourierspacing      = 0.16

; Temperature coupling
tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1     0.1
ref_t               = {temp}     {temp}

; Pressure coupling
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 1.0
ref_p               = {pressure}
compressibility     = 4.5e-5

; Periodic boundary conditions
pbc                 = xyz
DispCorr            = EnerPres

; Velocity generation
gen_vel             = no

; COM motion removal
refcoord_scaling    = com
comm-mode           = Linear
comm-grps           = Protein Non-Protein

; Pull code for umbrella sampling
pull                = yes
pull_ncoords        = 1
pull_ngroups        = 2
pull_group1_name    = {reference_group}
pull_group2_name    = {moving_group}
pull_coord1_type    = umbrella
pull_coord1_geometry = distance
pull_coord1_dim     = Y Y Y
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = 0.0
pull_coord1_k       = {pull_k}
pull_coord1_r0      = DISTANCE_PLACEHOLDER
pull-pbc-ref-prev-step-com = yes
"""
    
    def _generate_umbrella_run_scripts(self):
        """Generate shell scripts to run umbrella sampling"""
        smd_dir = self.pmf_system.smd_dir
        
        # Copy shared files to umbrella directory
        shared_files = ['topol.top', 'index.ndx']
        for filename in shared_files:
            source_file = smd_dir / filename
            dest_file = self.umbrella_dir / filename
            if source_file.exists():
                shutil.copy2(source_file, dest_file)
        
        # Copy LIG force field directory
        for lig_dir_name in ['LIG.amb2gmx', 'LIG.openff2gmx']:
            lig_source = smd_dir.parent / lig_dir_name
            lig_dest = self.umbrella_dir / lig_dir_name
            if lig_source.exists() and not lig_dest.exists():
                shutil.copytree(lig_source, lig_dest)
        
        # Individual window scripts
        individual_scripts = []
        for window in self.windows:
            window_dir = window['directory']
            
            script_content = f"""#!/bin/bash
######################################################
# UMBRELLA WINDOW {window['window_id']:03d} SCRIPT
# Distance: {window['distance']:.3f} nm
######################################################

set -e
echo "=== Running Umbrella Window {window['window_id']:03d} ==="
echo "Target distance: {window['distance']:.3f} nm"

# Create results directory
mkdir -p results

echo "Step 1: Generating TPR file..."
gmx grompp -f umbrella.mdp -c start.gro -n ../index.ndx -p ../topol.top -o umbrella.tpr -maxwarn 10

echo "Step 2: Running umbrella sampling..."
echo "This will take approximately 22 ns simulation time..."
gmx mdrun -s umbrella.tpr -deffnm umbrella -ntmpi 1 -ntomp 10 -v

echo "Step 3: Moving results..."
mv umbrella.* results/

echo "Step 4: Validating results..."
if [ -f "results/umbrella.gro" ] && [ -f "results/umbrella_pullf.xvg" ]; then
    echo "Window {window['window_id']:03d} completed successfully!"
else
    echo "✗ Window {window['window_id']:03d} failed!"
    exit 1
fi
"""

            script_file = window_dir / "run_window.sh"
            with open(script_file, 'w') as f:
                f.write(script_content)
            
            os.chmod(script_file, 0o755)
            individual_scripts.append(str(script_file))
        
        # Master run script
        master_script_content = f"""#!/bin/bash
######################################################
# MASTER UMBRELLA SAMPLING SCRIPT
# Total windows: {len(self.windows)}
######################################################

set -e
echo "=== Running All Umbrella Sampling Windows ==="
echo "Total windows: {len(self.windows)}"
echo "Estimated total time: Very long! Consider running in parallel."
echo ""

# Function to run a single window
run_window() {{
    local window_id=$1
    local window_dir="window_$(printf "%03d" $window_id)"
    
    echo "Starting window $window_id..."
    cd "$window_dir"
    
    if bash run_window.sh; then
        echo "Window $window_id completed"
        cd ..
        return 0
    else
        echo "✗ Window $window_id failed"
        cd ..
        return 1
    fi
}}

# Option 1: Sequential execution (very slow)
if [ "$1" = "sequential" ]; then
    echo "Running windows sequentially..."
    failed_windows=()
    
    for i in {{0..{len(self.windows)-1}}}; do
        if ! run_window $i; then
            failed_windows+=($i)
        fi
    done
    
    if [ ${{#failed_windows[@]}} -eq 0 ]; then
        echo "All windows completed successfully!"
    else
        echo "✗ Failed windows: ${{failed_windows[*]}}"
        exit 1
    fi

# Option 2: Parallel execution (recommended)
elif [ "$1" = "parallel" ]; then
    echo "Running windows in parallel..."
    echo "Make sure your system can handle multiple GROMACS jobs!"
    
    # Run windows in parallel (adjust -P based on your system)
    seq 0 {len(self.windows)-1} | xargs -n 1 -P 4 -I {{}} bash -c 'run_window {{}}'
    
    echo "Parallel execution completed!"

else
    echo "Usage: $0 [sequential|parallel]"
    echo ""
    echo "Options:"
    echo "  sequential - Run windows one by one (very slow)"
    echo "  parallel   - Run multiple windows simultaneously"
    echo ""
    echo "Alternatively, run individual windows manually:"
fi
"""
        
        master_script = self.umbrella_dir / "run_all_umbrella.sh"
        with open(master_script, 'w') as f:
            f.write(master_script_content)
        
        os.chmod(master_script, 0o755)
        
        logger.info("Umbrella run scripts generated")
        
        return {
            'master_script': master_script,
            'individual_scripts': individual_scripts
        }


class WindowSelector:
    """Advanced window selection strategies for umbrella sampling"""
    
    @staticmethod
    def select_by_force(force_data, threshold=100):
        """
        Select windows based on force threshold
        
        Parameters:
        -----------
        force_data : list
            List of (frame, force) tuples from SMD
        threshold : float
            Force threshold in kJ/mol/nm
            
        Returns:
        --------
        list : Selected frames where force exceeds threshold
        """
        selected = []
        for frame, force in force_data:
            if abs(force) > threshold:
                selected.append((frame, force))
        return selected
    
    @staticmethod 
    def select_adaptive(distance_data, min_spacing=0.05, max_spacing=0.2, gradient_threshold=0.1):
        """
        Adaptive window selection based on distance gradient
        
        Parameters:
        -----------
        distance_data : list
            List of (frame, distance) tuples
        min_spacing : float
            Minimum spacing between windows (nm)
        max_spacing : float
            Maximum spacing between windows (nm)
        gradient_threshold : float
            Gradient threshold for dense sampling
            
        Returns:
        --------
        list : Adaptively selected windows
        """
        if len(distance_data) < 2:
            return distance_data
        
        selected = [distance_data[0]]  # Always include first frame
        
        for i in range(1, len(distance_data)):
            current_frame, current_dist = distance_data[i]
            last_selected_frame, last_selected_dist = selected[-1]
            
            # Calculate spacing based on local gradient
            dist_diff = abs(current_dist - last_selected_dist)
            frame_diff = current_frame - last_selected_frame
            
            # Estimate local gradient
            if frame_diff > 0:
                gradient = dist_diff / frame_diff
                
                # Dense sampling in high gradient regions
                if gradient > gradient_threshold:
                    required_spacing = min_spacing
                else:
                    required_spacing = max_spacing
                
                # Add window if spacing requirement is met
                if dist_diff >= required_spacing:
                    selected.append((current_frame, current_dist))
        
        return selected
    
    @staticmethod
    def select_uniform(distance_data, n_windows=20):
        """
        Select uniformly spaced windows
        
        Parameters:
        -----------
        distance_data : list
            List of (frame, distance) tuples
        n_windows : int
            Number of windows to select
            
        Returns:
        --------
        list : Uniformly spaced windows
        """
        if len(distance_data) <= n_windows:
            return distance_data
        
        # Select every nth frame to get uniform spacing
        step = len(distance_data) // n_windows
        selected = [distance_data[i * step] for i in range(n_windows)]
        
        # Always include the last frame
        if selected[-1] != distance_data[-1]:
            selected.append(distance_data[-1])
        
        return selected