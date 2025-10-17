#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced SMD (Steered Molecular Dynamics) module compatible with step-by-step API
"""

import os
import shutil
from pathlib import Path
import subprocess
import numpy as np
import logging

from .utils import create_mdp_file, extract_pull_groups

logger = logging.getLogger(__name__)


class SMDManager:
    """Enhanced SMD simulation manager with step-by-step compatibility"""
    
    def __init__(self, pmf_system):
        self.pmf_system = pmf_system
        self.smd_dir = pmf_system.smd_dir
        self.config = pmf_system.config
        
        # Check if PMF system has been properly built with box remodeling
        self.pmf_system_built = self._check_pmf_system_built()
    
    def prepare(self, use_pmf_rebuilt_system: bool = True, **kwargs):
        """
        Prepare SMD simulation (step-by-step compatible)
        
        Parameters:
        -----------
        use_pmf_rebuilt_system : bool, optional
            Whether to use PMF-rebuilt system (default: True)
            If True, uses the PMF-optimized system with Z-axis alignment
            If False, uses the original MD results
        
        Returns:
        --------
        Dict : SMD preparation results with execution instructions
        """
        logger.info("Preparing SMD simulation...")
        
        if use_pmf_rebuilt_system and not self.pmf_system_built:
            logger.warning("PMF-rebuilt system not found. Consider running PMF builder first.")
            logger.warning("Falling back to original MD results...")
            use_pmf_rebuilt_system = False
        
        # Create SMD directory
        self.smd_dir.mkdir(exist_ok=True)
        
        # Copy necessary files (use PMF system if available)
        if use_pmf_rebuilt_system:
            self._prepare_files_from_pmf_system()
        else:
            self._prepare_files_from_md_results()
        
        # Create index file
        self._create_index_file()
        
        # Generate SMD MDP file (optimized for Z-axis pulling if using PMF system)
        self._create_smd_mdp(use_z_axis_optimized=use_pmf_rebuilt_system)
        
        # Generate run script
        run_script = self._generate_run_script(use_pmf_rebuilt_system)
        
        # Generate analysis script
        analysis_script = self._generate_analysis_script()
        
        results = {
            'smd_dir': str(self.smd_dir),
            'smd_mdp': str(self.smd_dir / "smd.mdp"),
            'run_script': str(run_script),
            'analysis_script': str(analysis_script),
            'system_type': 'pmf_rebuilt' if use_pmf_rebuilt_system else 'original_md',
            'z_axis_optimized': use_pmf_rebuilt_system,
            'status': 'prepared',
            'next_step': f'Run SMD simulation using: {run_script}'
        }
        
        logger.info(f"SMD preparation completed. System type: {'PMF-rebuilt (Z-optimized)' if use_pmf_rebuilt_system else 'Original MD'}")
        logger.info(f"Run SMD: {run_script}")
        return results
    
    def run(self, pull_groups=None, direction=None, rate=None, force=None):
        """
        Run SMD simulation (automated execution)
        
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
        # Prepare if not already done
        if not self.smd_dir.exists():
            self.prepare()
        
        os.chdir(self.smd_dir)
        
        # Run SMD
        self._run_smd_simulation()
        
        # Extract trajectory
        self._extract_trajectory()
        
        # Calculate distances
        self._calculate_distances()
        
        logger.info(f"SMD completed. Results in: {self.smd_dir}")
    
    def _check_pmf_system_built(self) -> bool:
        """Check if PMF system has been properly built"""
        pmf_system_dir = getattr(self.pmf_system, 'rebuilt_system_dir', None)
        
        if pmf_system_dir is None:
            return False
        
        pmf_system_path = Path(pmf_system_dir)
        required_files = [
            pmf_system_path / "solv_ions.gro",
            pmf_system_path / "topol.top"
        ]
        
        return all(f.exists() for f in required_files)
    
    def _prepare_files_from_pmf_system(self):
        """Prepare files from PMF-rebuilt system"""
        logger.info("Using PMF-rebuilt system (Z-axis optimized)")
        
        pmf_system_dir = Path(self.pmf_system.rebuilt_system_dir)
        
        # Check for equilibrated system first
        equilibrated_structure = pmf_system_dir.parent / "equilibration" / "npt" / "npt_final.gro"
        equilibrated_checkpoint = pmf_system_dir.parent / "equilibration" / "npt" / "npt_final.cpt"
        
        if equilibrated_structure.exists():
            logger.info("Using equilibrated PMF system")
            shutil.copy2(equilibrated_structure, self.smd_dir / "md.gro")
            if equilibrated_checkpoint.exists():
                shutil.copy2(equilibrated_checkpoint, self.smd_dir / "equilibrated.cpt")
        else:
            logger.info("Using unequilibrated PMF-rebuilt system")
            shutil.copy2(pmf_system_dir / "solv_ions.gro", self.smd_dir / "md.gro")
            
        # Copy topology
        shutil.copy2(pmf_system_dir / "topol.top", self.smd_dir / "topol.top")
        
        # Copy other necessary files
        for file_name in ["posre.itp"]:
            source_file = pmf_system_dir / file_name
            if source_file.exists():
                shutil.copy2(source_file, self.smd_dir / file_name)
        
        # Copy ligand force field directory
        self._copy_ligand_forcefield_dir(pmf_system_dir.parent)
        
        logger.info("PMF-rebuilt system files prepared")
    
    def _prepare_files_from_md_results(self):
        """Prepare files from original MD results (fallback)"""
        logger.info("Using original MD results (fallback mode)")
        
        # Get MD results directory
        md_dir = getattr(self.pmf_system, 'md_results_dir', self.pmf_system.system_dir)
        
        # Copy files to SMD directory
        self._prepare_files_copy_mode(md_dir)
    
    def _prepare_files_copy_mode(self, md_dir):
        """Prepare files by copying them to SMD directory with proper preprocessing"""
        logger.info("Preparing files for SMD simulation with preprocessing")
        
        # Step 1: Copy necessary files first
        self._copy_essential_files(md_dir)
        
        # Step 2: Preprocess structure file (remove PBC, center system)
        self._preprocess_structure_file(md_dir)
        
        # Step 3: Copy additional files
        self._copy_additional_files(md_dir)
    
    def _copy_essential_files(self, md_dir):
        """Copy essential files for PMF calculation"""
        logger.info("Copying essential files for PMF simulation")
        
        # Map of essential files to copy
        file_mapping = {
            'top': {
                'sources': ['topol.top'],
                'dest': 'topol.top',
                'description': 'Topology file'
            },
            'posre': {
                'sources': ['posre.itp'],
                'dest': 'posre.itp',
                'description': 'Position restraint file'
            },
            'tpr': {
                'sources': ['prod/md.tpr', 'md.tpr'],
                'dest': 'md.tpr',
                'description': 'MD binary file (required for trjconv)'
            }
        }
        
        for file_type, file_info in file_mapping.items():
            copied = False
            for source_name in file_info['sources']:
                source_file = md_dir / source_name
                if source_file.exists():
                    dest_file = self.smd_dir / file_info['dest']
                    if not dest_file.exists():
                        shutil.copy2(source_file, dest_file)
                        logger.info(f"Copied {file_info['description']}: {source_name}")
                    else:
                        logger.info(f"Using existing {file_info['description']}: {source_name}")
                    copied = True
                    break
            
            if not copied:
                # TPR file is optional for some operations, others are required
                if file_type == 'tpr':
                    logger.warning(f"Could not find {file_info['description']}. Some advanced features may not work.")
                    logger.warning(f"Tried: {file_info['sources']} in {md_dir}")
                else:
                    raise FileNotFoundError(
                        f"Could not find {file_info['description']}. "
                        f"Tried: {file_info['sources']} in {md_dir}"
                    )
    
    def _preprocess_structure_file(self, md_dir):
        """Preprocess structure file: remove PBC effects from MD final frame"""
        logger.info("Preprocessing structure file for PMF calculation")
        
        # Find structure file (MD final frame)
        structure_sources = ['prod/md.gro', 'md.gro']
        tpr_sources = ['prod/md.tpr', 'md.tpr']
        
        source_structure = None
        source_tpr = None
        
        # Find structure file
        for source_name in structure_sources:
            source_file = md_dir / source_name
            if source_file.exists():
                source_structure = source_file
                break
        
        if not source_structure:
            raise FileNotFoundError(
                f"Could not find MD structure file. "
                f"Tried: {structure_sources} in {md_dir}"
            )
        
        # Find TPR file (for PBC processing)
        for tpr_name in tpr_sources:
            tpr_file = md_dir / tpr_name
            if tpr_file.exists():
                source_tpr = tpr_file
                break
        
        dest_structure = self.smd_dir / "md.gro"
        
        # Process structure file to remove PBC
        self._process_structure_with_pbc_removal(source_structure, source_tpr, dest_structure)
    
    def _process_structure_with_pbc_removal(self, structure_file, tpr_file, output_file):
        """Process structure file to remove PBC effects"""
        logger.info(f"Processing structure file for PBC removal: {structure_file.name}")
        
        try:
            # Use gmx trjconv to remove PBC - unified approach for any structure file
            cmd = [
                "gmx", "trjconv",
                "-f", str(structure_file),  # Input can be .gro, .xtc, etc.
                "-s", str(tpr_file) if tpr_file else str(structure_file),
                "-o", str(output_file),
                "-pbc", "mol",
                "-quiet"
            ]
            
            process = subprocess.run(
                cmd,
                text=True,
                capture_output=True,
                cwd=str(self.smd_dir),
                check=True
            )
            
            logger.info(f"Successfully processed structure with PBC removal: {output_file.name}")
            
        except subprocess.CalledProcessError as e:
            logger.warning(f"PBC removal failed: {e.stderr}")
            logger.info("Using direct file copy as fallback")
            shutil.copy2(structure_file, output_file)
            logger.info(f"Copied structure file: {output_file.name}")
    
    def _copy_additional_files(self, md_dir):
        """Copy additional files after preprocessing"""
        logger.info("Copying additional force field files")
        self._copy_ligand_forcefield_dir(md_dir.parent)
    
    def _copy_ligand_forcefield_dir(self, parent_dir):
        """Copy ligand force field directory"""
        # Copy LIG force field directory
        for lig_dir_name in ['LIG.amb2gmx', 'LIG.openff2gmx']:
            lig_source = parent_dir / lig_dir_name
            if lig_source.exists():
                lig_dest = self.smd_dir.parent / lig_dir_name
                if not lig_dest.exists():
                    shutil.copytree(lig_source, lig_dest)
                    logger.info(f"Copied ligand directory: {lig_dir_name}")
                else:
                    logger.info(f"Using existing ligand directory: {lig_dir_name}")
                break
        else:
            logger.warning("No ligand force field directory found")
    
    def _create_index_file(self):
        """Create GROMACS index file automatically"""
        index_file = self.smd_dir / "index.ndx"
        
        if index_file.exists():
            logger.info("Index file already exists")
            return
        
        moving_group = self.config.get('moving_group', 'LIG')
        
        try:
            # Use gmx make_ndx to create index
            cmd = [
                "gmx", "make_ndx",
                "-f", "md.gro",
                "-o", "index.ndx"
            ]
            
            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(self.smd_dir)
            )
            
            # Create the moving group and quit
            input_commands = f"r {moving_group}\nq\n"
            stdout, stderr = process.communicate(input=input_commands)
            
            if process.returncode == 0:
                logger.info("Index file created successfully")
            else:
                # Fallback to basic index
                basic_cmd = ["gmx", "make_ndx", "-f", "md.gro", "-o", "index.ndx"]
                subprocess.run(basic_cmd, input="q\n", text=True, cwd=str(self.smd_dir), check=True)
                logger.info("Basic index file created")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create index file: {e}")
            raise RuntimeError(f"Index file creation failed: {e}")
    
    def _create_smd_mdp(self, use_z_axis_optimized: bool = True):
        """Create SMD MDP file with current configuration"""
        # Get configuration parameters
        reference_group = self.config.get('reference_group', 'Protein')
        moving_group = self.config.get('moving_group', 'LIG')
        pull_rate = self.config.get('smd', {}).get('pull_rate', 0.005)
        pull_k = self.config.get('smd', {}).get('pull_k', 1000.0)
        dt = self.config.get('simulation', {}).get('dt', 0.002)
        
        distance_start = self.config.get('distance', {}).get('start', 0.3)
        distance_end = self.config.get('distance', {}).get('end', 2.0)
        
        # Calculate simulation parameters
        distance_range = distance_end - distance_start
        smd_time_ps = distance_range / pull_rate
        nsteps = int(smd_time_ps / dt)
        output_interval = int(5.0 / dt)  # Output every 5 ps
        
        # Generate MDP content (Z-axis optimized if using PMF system)
        mdp_content = self._generate_smd_mdp_content(
            reference_group, moving_group, pull_rate, pull_k, 
            nsteps, dt, output_interval, distance_range, smd_time_ps,
            z_axis_optimized=use_z_axis_optimized
        )
        
        smd_mdp_file = self.smd_dir / "smd.mdp"
        with open(smd_mdp_file, 'w') as f:
            f.write(mdp_content)
        
        optimization_note = " (Z-axis optimized)" if use_z_axis_optimized else " (general geometry)"
        logger.info(f"SMD MDP file created{optimization_note}: {smd_mdp_file}")
    
    def _generate_smd_mdp_content(self, reference_group, moving_group, pull_rate, pull_k, 
                                 nsteps, dt, output_interval, distance_range, smd_time_ps,
                                 z_axis_optimized: bool = True):
        """Generate SMD MDP file content"""
        temp = self.config.get('simulation', {}).get('temperature', 310.0)
        pressure = self.config.get('simulation', {}).get('pressure', 1.0)
        
        # Pull direction: Z-axis for PMF-optimized systems, all directions for general
        pull_dim = "N N Y" if z_axis_optimized else "Y Y Y"
        system_note = "Z-axis optimized" if z_axis_optimized else "general geometry"
        
        # For PMF systems, we can use more aggressive settings
        if z_axis_optimized:
            # Optimized for PMF-rebuilt systems
            continuation = "no"  # Start fresh from equilibrated PMF system
            gen_vel = "yes"  # Generate new velocities for clean SMD start
        else:
            # Conservative settings for original MD systems
            continuation = "yes"
            gen_vel = "no"
        
        return f"""; SMD (Steered Molecular Dynamics) Parameters - {system_note}
title               = SMD simulation ({distance_range:.2f} nm in {smd_time_ps:.1f} ps)
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
continuation        = {continuation}
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
gen_vel             = {gen_vel}
{'gen_temp        = ' + str(temp) if gen_vel == 'yes' else ''}
{'gen_seed        = -1' if gen_vel == 'yes' else ''}

; COM motion removal
refcoord_scaling    = com
comm-mode           = Linear
comm-grps           = Protein Non-Protein

; Pull code for SMD - {'Z-axis pulling' if z_axis_optimized else 'distance pulling'}
pull                = yes
pull_ncoords        = 1
pull_ngroups        = 2
pull_group1_name    = {reference_group}
pull_group2_name    = {moving_group}
pull_coord1_type    = umbrella
pull_coord1_geometry = distance
pull_coord1_dim     = {pull_dim}{'  ; Z-axis optimized pulling' if z_axis_optimized else '  ; general distance pulling'}
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = {pull_rate}
pull_coord1_k       = {pull_k}
pull-pbc-ref-prev-step-com = yes
"""
    
    def _generate_run_script(self, use_pmf_system: bool = True):
        """Generate shell script to run SMD simulation"""
        distance_start = self.config.get('distance', {}).get('start', 0.3)
        distance_end = self.config.get('distance', {}).get('end', 2.0)
        pull_rate = self.config.get('smd', {}).get('pull_rate', 0.005)
        
        distance_range = distance_end - distance_start
        smd_time_ps = distance_range / pull_rate
        
        # Get file paths for script generation
        md_gro = 'md.gro'
        topol_top = 'topol.top'
        
        # Additional parameters for PMF systems
        system_type = "PMF-rebuilt (Z-optimized)" if use_pmf_system else "Original MD"
        checkpoint_option = "-t equilibrated.cpt" if use_pmf_system and (self.smd_dir / "equilibrated.cpt").exists() else ""
        
        script_content = f"""#!/bin/bash
######################################################
# SMD SIMULATION SCRIPT - {system_type}
# Distance: {distance_start:.2f} -> {distance_end:.2f} nm
# Pull rate: {pull_rate} nm/ps
# Simulation time: {smd_time_ps:.1f} ps
{'# System: Z-axis optimized for efficient pulling' if use_pmf_system else '# System: General MD results'}
######################################################

set -e
echo "=== SMD Simulation ({system_type}) ==="
{'echo "Using Z-axis optimized PMF system for efficient pulling"' if use_pmf_system else 'echo "Using original MD results"'}
echo ""

# Set GROMACS environment
export GMX_MAXBACKUP=-1  # Disable backup files
NCPU=10  # Adjust as needed

# Create results directory
mkdir -p results

echo "Step 1: Generating TPR file..."
# Check for equilibrated checkpoint
if [ -f "equilibrated.cpt" ]; then
    echo "Using equilibrated checkpoint for optimal SMD start"
    gmx grompp -f smd.mdp -c "{md_gro}" -n index.ndx -p "{topol_top}" -r "{md_gro}" {checkpoint_option} -o smd.tpr -maxwarn 10
else
    echo "Starting from structure file"
    gmx grompp -f smd.mdp -c "{md_gro}" -n index.ndx -p "{topol_top}" -r "{md_gro}" -o smd.tpr -maxwarn 10
fi

echo "Step 2: Running SMD simulation..."
echo "This may take a while depending on your system size and simulation time..."
echo "Pull direction: {'Z-axis (optimized)' if use_pmf_system else 'distance (general)'}"

# Try GPU first, fallback to CPU
if command -v nvidia-smi &> /dev/null && nvidia-smi &> /dev/null; then
    echo "Using GPU acceleration"
    gmx mdrun -s smd.tpr -deffnm smd -ntmpi 1 -ntomp $NCPU -nb gpu -pme gpu -v
else
    echo "Using CPU only"
    gmx mdrun -s smd.tpr -deffnm smd -ntmpi 1 -ntomp $NCPU -v
fi

echo "Step 3: Moving results to results directory..."
mv smd.* results/

echo "Step 4: Basic result check..."
if [ -f "results/smd.gro" ] && [ -f "results/smd.xtc" ] && [ -f "results/smd_pullf.xvg" ]; then
    echo "SMD simulation completed successfully!"
    echo "System type: {system_type}"
    echo "Generated files:"
    echo "  - results/smd.gro (final structure)"
    echo "  - results/smd.xtc (trajectory)"
    echo "  - results/smd_pullf.xvg (pull force)"
    echo "  - results/smd_pullx.xvg (pull distance)"
    echo ""
    echo "Next steps:"
    echo "1. Run analysis: bash analyze_smd.sh"
    echo "2. Then proceed to umbrella sampling preparation"
else
    echo "✗ SMD simulation failed or incomplete!"
    echo "Check the mdrun output for errors."
    exit 1
fi
"""
        
        run_script = self.smd_dir / "run_smd.sh"
        with open(run_script, 'w') as f:
            f.write(script_content)
        
        os.chmod(run_script, 0o755)
        logger.info(f"SMD run script created: {run_script}")
        
        return run_script
    
    def _generate_analysis_script(self):
        """Generate script to analyze SMD results"""
        script_content = """#!/bin/bash
######################################################
# SMD ANALYSIS SCRIPT
######################################################

set -e
echo "=== SMD Analysis ==="

# Check if results exist
if [ ! -f "results/smd_pullf.xvg" ] || [ ! -f "results/smd_pullx.xvg" ]; then
    echo "✗ SMD results not found! Run SMD simulation first."
    exit 1
fi

echo "Analyzing SMD results..."

# Create analysis directory
mkdir -p analysis

# Copy analysis files
cp results/smd_pullf.xvg analysis/
cp results/smd_pullx.xvg analysis/

echo "SMD analysis files prepared"
echo "Files available for plotting:"
echo "  - analysis/smd_pullf.xvg (force vs time)"
echo "  - analysis/smd_pullx.xvg (distance vs time)"
echo ""
echo "Use Python analysis functions for detailed plots"
"""
        
        analysis_script = self.smd_dir / "analyze_smd.sh"
        with open(analysis_script, 'w') as f:
            f.write(script_content)
        
        os.chmod(analysis_script, 0o755)
        logger.info(f"SMD analysis script created: {analysis_script}")
        
        return analysis_script
    
    def _run_smd_simulation(self):
        """Run the SMD simulation (automated execution)"""
        results_dir = self.smd_dir / "results"
        results_dir.mkdir(exist_ok=True)
        
        logger.info("Preparing SMD simulation...")
        subprocess.run([
            "gmx", "grompp",
            "-f", "smd.mdp",
            "-c", "md.gro",
            "-p", "topol.top",
            "-n", "index.ndx", 
            "-r", "md.gro",
            "-o", "smd.tpr",
            "-maxwarn", "10"
        ], check=True)
        
        logger.info("Running SMD simulation...")
        subprocess.run([
            "gmx", "mdrun",
            "-s", "smd.tpr",
            "-deffnm", "smd",
            "-ntmpi", "1",
            "-ntomp", "10",
            "-v"
        ], check=True)
        
        # Move results to results directory
        for result_file in self.smd_dir.glob("smd.*"):
            shutil.move(result_file, results_dir / result_file.name)
        
        logger.info(f"SMD results moved to: {results_dir}")
    
    def _extract_trajectory(self):
        """Extract configurations from trajectory for umbrella sampling"""
        logger.info("Extracting trajectory configurations...")
        
        results_dir = self.smd_dir / "results"
        conf_dir = self.smd_dir / "trajectory_frames"
        conf_dir.mkdir(exist_ok=True)
        
        # Extract frames at regular intervals
        n_frames = 300  # Extract 300 frames
        
        # Calculate frame extraction parameters  
        distance_start = self.config.get('distance', {}).get('start', 0.3)
        distance_end = self.config.get('distance', {}).get('end', 2.0)
        pull_rate = self.config.get('smd', {}).get('pull_rate', 0.005)
        
        distance_range = distance_end - distance_start
        smd_time_ps = distance_range / pull_rate
        frame_interval_ps = smd_time_ps / n_frames
        
        try:
            subprocess.run([
                "gmx", "trjconv",
                "-s", str(results_dir / "smd.tpr"),
                "-f", str(results_dir / "smd.xtc"),
                "-o", str(conf_dir / "frame.gro"),
                "-sep",
                "-dt", str(frame_interval_ps)
            ], input="0\n", text=True, check=True)
            
            # Count extracted frames
            extracted_frames = len(list(conf_dir.glob("frame*.gro")))
            logger.info(f"Extracted {extracted_frames} trajectory frames")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Frame extraction failed: {e}")
            raise RuntimeError(f"Failed to extract trajectory frames: {e}")
    
    def _calculate_distances(self):
        """Calculate COM distances for trajectory frames"""
        conf_dir = self.smd_dir / "trajectory_frames"
        results_dir = self.smd_dir / "results"
        
        if not conf_dir.exists():
            logger.warning("No trajectory frames found to calculate distances")
            return
        
        reference_group = self.config.get('reference_group', 'Protein')
        moving_group = self.config.get('moving_group', 'LIG')
        
        # Get list of frame files
        frame_files = sorted(conf_dir.glob("frame*.gro"))
        
        logger.info(f"Calculating distances for {len(frame_files)} frames...")
        
        distance_data = []
        temp_dir = self.smd_dir / "temp_distance"
        temp_dir.mkdir(exist_ok=True)
        
        try:
            for frame_file in frame_files:
                # Extract frame number from filename
                frame_num = int(frame_file.stem.replace('frame', ''))
                
                # Calculate distance using gmx distance
                distance_cmd = [
                    "gmx", "distance",
                    "-s", str(results_dir / "smd.tpr"),
                    "-f", str(frame_file),
                    "-n", "index.ndx",
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
            
            # Write distances file
            distances_file = self.smd_dir / "distances.dat"
            with open(distances_file, 'w') as f:
                f.write("# Frame Distance(nm)\n")
                for frame, distance in sorted(distance_data):
                    f.write(f"{frame} {distance:.6f}\n")
            
            logger.info(f"Distance calculation completed: {len(distance_data)} frames")
            logger.info(f"Results saved to: {distances_file}")
            
        finally:
            # Clean up temporary directory
            if temp_dir.exists():
                shutil.rmtree(temp_dir)


class SMDBuilder:
    """Helper class to build custom SMD protocols"""
    
    @staticmethod
    def create_custom_mdp(output_file, reference_group='Protein', moving_group='LIG', 
                         pull_rate=0.005, pull_k=1000.0, nsteps=2000000, dt=0.002, **params):
        """
        Create custom SMD MDP file with specified parameters
        
        Parameters:
        -----------
        output_file : str
            Output MDP file path
        reference_group : str
            Reference group name
        moving_group : str
            Moving group name  
        pull_rate : float
            Pull rate in nm/ps
        pull_k : float
            Force constant in kJ/mol/nm²
        nsteps : int
            Number of simulation steps
        dt : float
            Time step in ps
        **params : dict
            Additional MDP parameters
        """
        temp = params.get('temperature', 310.0)
        pressure = params.get('pressure', 1.0)
        output_interval = params.get('output_interval', int(5.0 / dt))
        
        mdp_content = f"""; Custom SMD Parameters
title               = Custom SMD simulation
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

; Pull code for SMD
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
pull_coord1_rate    = {pull_rate}
pull_coord1_k       = {pull_k}
pull-pbc-ref-prev-step-com = yes
"""
        
        # Add any additional parameters
        if 'additional_parameters' in params:
            mdp_content += "\n; Additional custom parameters\n"
            for key, value in params['additional_parameters'].items():
                mdp_content += f"{key}             = {value}\n"
        
        with open(output_file, 'w') as f:
            f.write(mdp_content)
        
        logger.info(f"Custom SMD MDP file created: {output_file}")
    
    @staticmethod
    def estimate_simulation_time(distance_start, distance_end, pull_rate):
        """
        Estimate SMD simulation time
        
        Parameters:
        -----------
        distance_start : float
            Starting distance (nm)
        distance_end : float
            Ending distance (nm)
        pull_rate : float
            Pull rate (nm/ps)
            
        Returns:
        --------
        float : Estimated simulation time in ps
        """
        distance_range = abs(distance_end - distance_start)
        return distance_range / pull_rate
    
    @staticmethod
    def validate_smd_parameters(config):
        """
        Validate SMD parameters
        
        Parameters:
        -----------
        config : dict
            SMD configuration dictionary
            
        Returns:
        --------
        dict : Validation results
        """
        issues = []
        warnings = []
        
        # Check pull rate
        pull_rate = config.get('smd', {}).get('pull_rate', 0.005)
        if pull_rate > 0.01:
            warnings.append(f"Pull rate {pull_rate} nm/ps is quite fast, may cause artifacts")
        elif pull_rate < 0.001:
            warnings.append(f"Pull rate {pull_rate} nm/ps is very slow, simulation will take long")
        
        # Check force constant
        pull_k = config.get('smd', {}).get('pull_k', 1000.0)
        if pull_k < 500:
            warnings.append(f"Force constant {pull_k} kJ/mol/nm² is low, may cause poor pulling")
        elif pull_k > 5000:
            warnings.append(f"Force constant {pull_k} kJ/mol/nm² is very high, may cause instability")
        
        # Check distance range
        distance_start = config.get('distance', {}).get('start', 0.3)
        distance_end = config.get('distance', {}).get('end', 2.0)
        distance_range = distance_end - distance_start
        
        if distance_range < 0.5:
            issues.append(f"Distance range {distance_range:.2f} nm is too small")
        elif distance_range > 3.0:
            warnings.append(f"Distance range {distance_range:.2f} nm is quite large")
        
        return {
            'valid': len(issues) == 0,
            'issues': issues,
            'warnings': warnings,
            'estimated_time_ps': SMDBuilder.estimate_simulation_time(distance_start, distance_end, pull_rate)
        }