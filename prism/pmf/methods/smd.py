#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unified SMD (Steered Molecular Dynamics) Module

This module provides both legacy compatibility and modern modular interface
for SMD calculations, eliminating code duplication while maintaining backward
compatibility.
"""

import os
import shutil
import time
from pathlib import Path
from typing import Dict, Any, Optional, List, Union
import subprocess
import numpy as np
import logging

# Import from optimization systems
from ..utils.exceptions import (
    SMDPreparationError, SMDExecutionError, RequiredFileNotFoundError,
    SystemNotFoundException, PrerequisiteNotMetError, PMFErrorCode
)
from ..utils.error_handling import (
    with_retry, error_context, safe_file_operation, validate_prerequisites,
    handle_external_tool_error
)

try:
    from ..utils.logging_system import PrismLogger, LogLevel, EventType
    from ..utils.module_communication import ModuleCommunicator, MessageBroker
except ImportError:
    # Fallback for when optimization systems are not available
    PrismLogger = logging.getLogger
    ModuleCommunicator = None
    MessageBroker = None

logger = logging.getLogger(__name__)


class SMDCore:
    """Core SMD functionality - shared between legacy and modular interfaces"""
    
    def __init__(self, pmf_system=None, config=None):
        self.pmf_system = pmf_system
        self.config = config or {}
        self.smd_dir = None
        
        # Initialize from pmf_system if available
        if pmf_system:
            self.smd_dir = getattr(pmf_system, 'smd_dir', None)
            if hasattr(pmf_system, 'config'):
                self.config = {**pmf_system.config, **self.config}
        
        # Check PMF system status
        self.pmf_system_built = self._check_pmf_system_built() if pmf_system else False
    
    def _check_pmf_system_built(self) -> bool:
        """Check if PMF system has been properly built"""
        if not self.pmf_system:
            return False
            
        pmf_system_dir = getattr(self.pmf_system, 'rebuilt_system_dir', None)
        if pmf_system_dir is None:
            pmf_system_dir = getattr(self.pmf_system, 'md_results_dir', None)
        
        if pmf_system_dir is None:
            return False
        
        pmf_system_path = Path(pmf_system_dir)
        required_files = [
            pmf_system_path / "solv_ions.gro",
            pmf_system_path / "topol.top"
        ]
        
        # Ensure this is a PMF system
        is_pmf_system = (
            "GMX_PMF_SYSTEM" in str(pmf_system_path) or
            pmf_system_path.name == "GMX_PMF_SYSTEM" or
            (pmf_system_path.parent / "GMX_PMF_SYSTEM").exists()
        )
        
        return all(f.exists() for f in required_files) and is_pmf_system
    
    def validate_smd_prerequisites(self, use_pmf_system: bool = True) -> None:
        """Validate all prerequisites for SMD preparation"""
        requirements = {
            'smd_directory_writable': lambda: self._check_directory_writable(self.smd_dir),
            'sufficient_disk_space': lambda: self._check_disk_space(self.smd_dir, required_gb=5.0)
        }
        
        if use_pmf_system:
            requirements.update({
                'pmf_system_built': lambda: self._check_pmf_system_built_validated(),
                'pmf_system_files': lambda: self._check_pmf_system_files()
            })
        
        validate_prerequisites(requirements, 'SMD preparation')
    
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
            free_bytes = shutil.disk_usage(directory).free
            free_gb = free_bytes / (1024**3)
            
            if free_gb < required_gb:
                from ..utils.exceptions import InsufficientDiskSpaceError
                raise InsufficientDiskSpaceError(
                    required_space_gb=required_gb,
                    available_space_gb=free_gb,
                    operation="SMD preparation"
                )
            return True
        except Exception as exc:
            logger.warning(f"Could not check disk space: {exc}")
            return True
    
    def _check_pmf_system_built_validated(self) -> bool:
        """Check if PMF system has been properly built with detailed validation"""
        if not self.pmf_system:
            return False
            
        try:
            pmf_system_dir = getattr(self.pmf_system, 'rebuilt_system_dir', None)
            if pmf_system_dir is None:
                pmf_system_dir = getattr(self.pmf_system, 'md_results_dir', None)
            
            if pmf_system_dir is None:
                raise SystemNotFoundException(
                    system_path="unknown",
                    searched_paths=["rebuilt_system_dir", "md_results_dir"]
                )
            
            pmf_system_path = Path(pmf_system_dir)
            required_files = [
                pmf_system_path / "solv_ions.gro",
                pmf_system_path / "topol.top"
            ]
            
            missing_files = [str(f) for f in required_files if not f.exists()]
            if missing_files:
                raise RequiredFileNotFoundError(
                    file_path=pmf_system_path,
                    step="PMF system validation",
                    alternative_files=missing_files
                )
            
            return True
            
        except Exception as exc:
            if isinstance(exc, (RequiredFileNotFoundError, SystemNotFoundException)):
                raise exc
            else:
                raise SMDPreparationError(
                    reason=f"PMF system validation failed: {exc}",
                    missing_components=["PMF system validation"]
                ) from exc
    
    def _check_pmf_system_files(self) -> bool:
        """Check that all required PMF system files are present and valid"""
        if not self.pmf_system:
            return False
            
        pmf_system_dir = Path(getattr(self.pmf_system, 'rebuilt_system_dir', None) or 
                             getattr(self.pmf_system, 'md_results_dir', None))
        
        if not pmf_system_dir:
            return False
        
        essential_files = {
            'structure': pmf_system_dir / "solv_ions.gro",
            'topology': pmf_system_dir / "topol.top"
        }
        
        for file_type, file_path in essential_files.items():
            if not file_path.exists():
                raise RequiredFileNotFoundError(
                    file_path=str(file_path),
                    step="PMF system files check"
                )
            
            if file_path.stat().st_size == 0:
                from ..utils.exceptions import FileCorruptionError
                raise FileCorruptionError(
                    file_path=str(file_path),
                    expected_format=file_type,
                    actual_issue="File is empty"
                )
        
        return True
    
    def generate_smd_mdp(self, use_z_axis_optimized: bool = True, **params) -> str:
        """Generate SMD MDP file content using user-defined parameters"""
        # Get configuration parameters - user parameters have priority
        reference_group = params.get('reference_group') or self.config.get('reference_group', 'Protein')
        moving_group = params.get('moving_group') or self.config.get('moving_group', 'LIG')

        # Use user-defined parameters or method parameters
        pull_rate = params.get('pull_rate') or self.pull_rate
        pull_k = params.get('pull_k') or self.pull_k
        dt = params.get('dt') or self.dt
        nsteps = params.get('nsteps') or self.nsteps

        temp = params.get('temperature') or self.temperature
        pressure = params.get('pressure') or self.pressure
        output_interval = params.get('output_interval') or int(5.0 / dt) if dt else 2500
        
        # Pull direction: Z-axis for PMF-optimized systems, all directions for general
        pull_dim = "N N Y" if use_z_axis_optimized else "Y Y Y"
        system_note = "Z-axis optimized" if use_z_axis_optimized else "general geometry"
        
        # For PMF systems, we can use more aggressive settings
        if use_z_axis_optimized:
            continuation = "no"
            gen_vel = "yes"
        else:
            continuation = "yes"
            gen_vel = "no"
        
        return f"""; SMD (Steered Molecular Dynamics) Parameters - {system_note}
title               = SMD simulation
define              = -DPOSRES_P
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

; Pull code for SMD - {'Z-axis pulling' if use_z_axis_optimized else 'distance pulling'}
pull                = yes
pull_ncoords        = 1
pull_ngroups        = 2
pull_group1_name    = {reference_group}
pull_group2_name    = {moving_group}
pull_coord1_type    = umbrella
pull_coord1_geometry = distance
pull_coord1_dim     = {pull_dim}{'  ; Z-axis optimized pulling' if use_z_axis_optimized else '  ; general distance pulling'}
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = {pull_rate}
pull_coord1_k       = {pull_k}
pull-pbc-ref-prev-step-com = yes
"""

    def _locate_equilibrated_system(self, pmf_system_dir: Path):
        """Locate equilibrated system files with intelligent fallback

        Search priority:
        1. Standard structure: equilibration/npt/npt.gro (current standard)
        2. Legacy structure: npt/npt.gro (backward compatibility)

        Returns:
            tuple: (structure_path, checkpoint_path) or (None, None)
        """
        # Priority 1: Standard structure (equilibration/npt/)
        equilibrated_structure = pmf_system_dir / "equilibration" / "npt" / "npt.gro"
        equilibrated_checkpoint = pmf_system_dir / "equilibration" / "npt" / "npt.cpt"

        if equilibrated_structure.exists():
            return equilibrated_structure, equilibrated_checkpoint

        # Priority 2: Legacy structure (backward compatibility)
        legacy_structure = pmf_system_dir / "npt" / "npt.gro"
        legacy_checkpoint = pmf_system_dir / "npt" / "npt.cpt"

        if legacy_structure.exists():
            return legacy_structure, legacy_checkpoint

        # No equilibrated system found
        return None, None

    def prepare_files(self, smd_dir: Path, use_pmf_system: bool = True) -> List[str]:
        """Prepare SMD files - unified preparation logic"""
        prepared_files = []
        
        smd_dir.mkdir(exist_ok=True)
        
        if use_pmf_system and self.pmf_system_built:
            prepared_files.extend(self._prepare_files_from_pmf_system(smd_dir))
        else:
            prepared_files.extend(self._prepare_files_from_md_results(smd_dir))
        
        # Generate MDP file
        mdp_content = self.generate_smd_mdp(use_z_axis_optimized=use_pmf_system)
        mdp_file = smd_dir / "smd.mdp"
        with open(mdp_file, 'w') as f:
            f.write(mdp_content)
        prepared_files.append(str(mdp_file))
        
        # Create index file if needed
        self._create_index_file(smd_dir)
        
        return prepared_files
    
    def _prepare_files_from_pmf_system(self, smd_dir: Path) -> List[str]:
        """Prepare files from PMF-rebuilt system"""
        logger.info("Using PMF-rebuilt system (Z-axis optimized)")
        prepared_files = []

        pmf_system_dir = Path(getattr(self.pmf_system, 'rebuilt_system_dir', None) or
                             getattr(self.pmf_system, 'md_results_dir', None))

        # Locate equilibrated system with smart fallback
        equilibrated_structure, equilibrated_checkpoint = self._locate_equilibrated_system(pmf_system_dir)

        target_structure = smd_dir / "md.gro"
        target_checkpoint = smd_dir / "equilibrated.cpt"
        target_topology = smd_dir / "topol.top"

        if equilibrated_structure and equilibrated_structure.exists():
            logger.info(f"Using equilibrated system: {equilibrated_structure.relative_to(pmf_system_dir)}")
            shutil.copy2(equilibrated_structure, target_structure)
            prepared_files.append(str(target_structure))

            if equilibrated_checkpoint and equilibrated_checkpoint.exists():
                shutil.copy2(equilibrated_checkpoint, target_checkpoint)
                prepared_files.append(str(target_checkpoint))
        else:
            logger.warning("No equilibrated system found, using unequilibrated solv_ions.gro")
            logger.warning("For better results, run equilibration first")
            base_structure = pmf_system_dir / "solv_ions.gro"
            shutil.copy2(base_structure, target_structure)
            prepared_files.append(str(target_structure))
        
        # Copy topology
        topology_source = pmf_system_dir / "topol.top"
        shutil.copy2(topology_source, target_topology)
        prepared_files.append(str(target_topology))
        
        return prepared_files
    
    def _prepare_files_from_md_results(self, smd_dir: Path) -> List[str]:
        """Prepare files from original MD results (fallback)"""
        logger.info("Using original MD results (fallback mode)")
        prepared_files = []
        
        md_dir = getattr(self.pmf_system, 'md_results_dir', None) if self.pmf_system else None
        if not md_dir:
            return prepared_files
            
        md_dir = Path(md_dir)
        
        # Copy essential files
        file_mapping = {
            'topol.top': 'topol.top',
            'md.gro': 'md.gro'
        }
        
        for source_name, dest_name in file_mapping.items():
            source_file = md_dir / source_name
            if source_file.exists():
                dest_file = smd_dir / dest_name
                shutil.copy2(source_file, dest_file)
                prepared_files.append(str(dest_file))
        
        return prepared_files
    
    def _create_index_file(self, smd_dir: Path):
        """Create GROMACS index file automatically"""
        index_file = smd_dir / "index.ndx"
        
        if index_file.exists():
            logger.info("Index file already exists")
            return
        
        moving_group = self.config.get('moving_group', 'LIG')
        
        try:
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
                cwd=str(smd_dir)
            )
            
            input_commands = f"r {moving_group}\nq\n"
            stdout, stderr = process.communicate(input=input_commands)
            
            if process.returncode == 0:
                logger.info("Index file created successfully")
            else:
                basic_cmd = ["gmx", "make_ndx", "-f", "md.gro", "-o", "index.ndx"]
                subprocess.run(basic_cmd, input="q\n", text=True, cwd=str(smd_dir), check=True)
                logger.info("Basic index file created")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create index file: {e}")
            raise RuntimeError(f"Index file creation failed: {e}")
    
    def run_smd_simulation(self, smd_dir: Path) -> Dict[str, Any]:
        """Run SMD simulation"""
        results_dir = smd_dir / "results"
        results_dir.mkdir(exist_ok=True)
        
        # GROMACS preprocessing
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
        ], cwd=smd_dir, check=True)
        
        # Run simulation
        logger.info("Running SMD simulation...")
        subprocess.run([
            "gmx", "mdrun",
            "-s", "smd.tpr",
            "-deffnm", "smd",
            "-ntmpi", "1",
            "-ntomp", "10",
            "-v"
        ], cwd=smd_dir, check=True)
        
        # Move results to results directory
        for result_file in smd_dir.glob("smd.*"):
            shutil.move(result_file, results_dir / result_file.name)
        
        logger.info(f"SMD results moved to: {results_dir}")
        
        return {
            'smd_dir': str(smd_dir),
            'results_dir': str(results_dir),
            'trajectory_file': str(results_dir / 'smd.xtc'),
            'pullx_file': str(results_dir / 'smd_pullx.xvg'),
            'pullf_file': str(results_dir / 'smd_pullf.xvg')
        }
    
    def parse_xvg_file(self, xvg_file: Path) -> List[float]:
        """Parse GROMACS XVG file and extract data"""
        data = []
        try:
            with open(xvg_file, 'r') as f:
                for line in f:
                    if not line.startswith('#') and not line.startswith('@'):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            data.append(float(parts[1]))
        except Exception as e:
            logger.warning(f"Error parsing {xvg_file}: {e}")
        
        return data


class SMDManager(SMDCore):
    """Legacy SMD Manager - maintains backward compatibility"""
    
    def __init__(self, pmf_system):
        super().__init__(pmf_system)
        logger.info("SMDManager initialized (legacy compatibility mode)")
    
    def prepare(self, use_pmf_rebuilt_system: bool = True, **kwargs):
        """Prepare SMD simulation (legacy interface)"""
        with error_context("SMD preparation", {"use_pmf_system": use_pmf_rebuilt_system}):
            logger.info("Preparing SMD simulation...")
            
            # Validate prerequisites
            self.validate_smd_prerequisites(use_pmf_rebuilt_system)
            
            if use_pmf_rebuilt_system and not self.pmf_system_built:
                logger.warning("PMF-rebuilt system not found. Falling back to original MD results...")
                use_pmf_rebuilt_system = False
            
            # Prepare files
            prepared_files = self.prepare_files(self.smd_dir, use_pmf_rebuilt_system)
            
            # Generate scripts
            run_script = self._generate_run_script(use_pmf_rebuilt_system)
            analysis_script = self._generate_analysis_script()
            
            results = {
                'smd_dir': str(self.smd_dir),
                'smd_mdp': str(self.smd_dir / "smd.mdp"),
                'run_script': str(run_script),
                'analysis_script': str(analysis_script),
                'system_type': 'pmf_rebuilt' if use_pmf_rebuilt_system else 'original_md',
                'z_axis_optimized': use_pmf_rebuilt_system,
                'status': 'prepared',
                'prepared_files': prepared_files,
                'next_step': f'Run SMD simulation using: {run_script}'
            }
            
            logger.info(f"SMD preparation completed. System type: {'PMF-rebuilt (Z-optimized)' if use_pmf_rebuilt_system else 'Original MD'}")
            return results
    
    def run(self, **kwargs):
        """Run SMD simulation (legacy interface)"""
        if not self.smd_dir.exists():
            self.prepare()
        
        results = self.run_smd_simulation(self.smd_dir)
        logger.info(f"SMD completed. Results in: {self.smd_dir}")
        return results
    
    def _generate_run_script(self, use_pmf_system: bool = True) -> Path:
        """Generate shell script to run SMD simulation"""
        system_type = "PMF-rebuilt (Z-optimized)" if use_pmf_system else "Original MD"
        checkpoint_option = "-t equilibrated.cpt" if use_pmf_system and (self.smd_dir / "equilibrated.cpt").exists() else ""
        
        script_content = f"""#!/bin/bash
######################################################
# SMD SIMULATION SCRIPT - {system_type}
######################################################

set -e
echo "=== SMD Simulation ({system_type}) ==="

export GMX_MAXBACKUP=-1
NCPU=10

mkdir -p results

echo "Step 1: Generating TPR file..."
gmx grompp -f smd.mdp -c md.gro -n index.ndx -p topol.top -r md.gro {checkpoint_option} -o smd.tpr -maxwarn 10

echo "Step 2: Running SMD simulation..."
gmx mdrun -s smd.tpr -deffnm smd -ntmpi 1 -ntomp $NCPU -v

echo "Step 3: Moving results..."
mv smd.* results/

echo "SMD simulation completed successfully!"
"""
        
        run_script = self.smd_dir / "run_smd.sh"
        with open(run_script, 'w') as f:
            f.write(script_content)
        
        os.chmod(run_script, 0o755)
        logger.info(f"SMD run script created: {run_script}")
        
        return run_script
    
    def _generate_analysis_script(self) -> Path:
        """Generate script to analyze SMD results"""
        script_content = """#!/bin/bash
######################################################
# SMD ANALYSIS SCRIPT
######################################################

set -e
echo "=== SMD Analysis ==="

if [ ! -f "results/smd_pullf.xvg" ] || [ ! -f "results/smd_pullx.xvg" ]; then
    echo "✗ SMD results not found! Run SMD simulation first."
    exit 1
fi

echo "Analyzing SMD results..."
mkdir -p analysis
cp results/smd_pullf.xvg analysis/
cp results/smd_pullx.xvg analysis/

echo "SMD analysis files prepared"
echo "Files available for plotting:"
echo "  - analysis/smd_pullf.xvg (force vs time)"
echo "  - analysis/smd_pullx.xvg (distance vs time)"
"""
        
        analysis_script = self.smd_dir / "analyze_smd.sh"
        with open(analysis_script, 'w') as f:
            f.write(script_content)
        
        os.chmod(analysis_script, 0o755)
        logger.info(f"SMD analysis script created: {analysis_script}")
        
        return analysis_script


# Legacy compatibility alias - SMDBuilder is now just an alias to SMDManager
SMDBuilder = SMDManager