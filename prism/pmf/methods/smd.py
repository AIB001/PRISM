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
from ..interfaces.base import PMFModuleInterface, ModuleResult, ModulePhase, ModuleStatus
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
        
        # Check for equilibrated system first
        equilibrated_structure = pmf_system_dir / "npt" / "npt.gro"
        equilibrated_checkpoint = pmf_system_dir / "npt" / "npt.cpt"
        
        target_structure = smd_dir / "md.gro"
        target_checkpoint = smd_dir / "equilibrated.cpt"
        target_topology = smd_dir / "topol.top"
        
        if equilibrated_structure.exists():
            logger.info("Using equilibrated PMF system")
            shutil.copy2(equilibrated_structure, target_structure)
            prepared_files.append(str(target_structure))
            
            if equilibrated_checkpoint.exists():
                shutil.copy2(equilibrated_checkpoint, target_checkpoint)
                prepared_files.append(str(target_checkpoint))
        else:
            logger.info("Using unequilibrated PMF-rebuilt system")
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


class SMDModule(SMDCore, PMFModuleInterface):
    """Modern modular SMD implementation - Fully user-customizable"""

    def __init__(self, config: Optional[Dict[str, Any]] = None,
                 logger: Optional[PrismLogger] = None,
                 communicator: Optional[ModuleCommunicator] = None,
                 **user_params):
        # Initialize PMFModuleInterface
        PMFModuleInterface.__init__(self, "smd", config, logger)
        # Initialize SMDCore
        SMDCore.__init__(self, config=config)

        self.communicator = communicator
        self.system_path = None
        self.output_files = []

        # SMD-specific configuration - FULLY USER CUSTOMIZABLE
        self.smd_config = self.config.get('smd', {})

        # Apply user parameters with highest priority
        self.smd_config.update(user_params)

        # User-defined parameters (no defaults, user must specify)
        self.pull_rate = self.smd_config.get('pull_rate')
        self.nsteps = self.smd_config.get('nsteps')
        self.dt = self.smd_config.get('dt')
        self.pull_k = self.smd_config.get('pull_k')

        # Optional parameters with sensible defaults only
        self.temperature = self.smd_config.get('temperature', 310.0)
        self.pressure = self.smd_config.get('pressure', 1.0)

        self.logger.info(f"SMD Module initialized - User-defined parameters:")
        self.logger.info(f"  pull_rate: {self.pull_rate}")
        self.logger.info(f"  nsteps: {self.nsteps}")
        self.logger.info(f"  dt: {self.dt}")
        self.logger.info(f"  pull_k: {self.pull_k}")
    
    def validate_inputs(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Validate SMD input parameters - User-defined validation"""
        validation_results = {
            'valid': True,
            'errors': [],
            'warnings': []
        }

        try:
            # Check required inputs
            required_keys = ['system_path', 'output_dir']
            for key in required_keys:
                if key not in inputs:
                    validation_results['errors'].append(f"Missing required input: {key}")
                    validation_results['valid'] = False

            if not validation_results['valid']:
                return validation_results

            # Validate system path
            system_path = Path(inputs['system_path'])
            if not self.validate_md_system(system_path):
                validation_results['errors'].append(f"Invalid MD system at {system_path}")
                validation_results['valid'] = False

            # Check that user has provided required SMD parameters
            required_smd_params = ['pull_rate', 'nsteps', 'dt', 'pull_k']
            missing_params = []
            for param in required_smd_params:
                if getattr(self, param) is None:
                    missing_params.append(param)

            if missing_params:
                validation_results['errors'].append(
                    f"Missing required SMD parameters: {', '.join(missing_params)}"
                )
                validation_results['valid'] = False

            # Basic sanity checks (no recommendations, just prevent obvious errors)
            if self.pull_rate is not None and self.pull_rate <= 0:
                validation_results['errors'].append("pull_rate must be positive")
                validation_results['valid'] = False

            if self.nsteps is not None and self.nsteps <= 0:
                validation_results['errors'].append("nsteps must be positive")
                validation_results['valid'] = False

            if self.dt is not None and self.dt <= 0:
                validation_results['errors'].append("dt must be positive")
                validation_results['valid'] = False

            if self.pull_k is not None and self.pull_k <= 0:
                validation_results['errors'].append("pull_k must be positive")
                validation_results['valid'] = False

            self.logger.info(f"SMD input validation completed: {len(validation_results['errors'])} errors")

        except Exception as e:
            validation_results['valid'] = False
            validation_results['errors'].append(f"Validation error: {str(e)}")
            self.logger.error(f"SMD validation failed: {e}")

        return validation_results
    
    def prepare(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare SMD simulation (modular interface)"""
        try:
            self.system_path = Path(inputs['system_path'])
            output_dir = Path(inputs['output_dir'])
            
            # Setup SMD directory
            smd_dir = output_dir / "smd"
            self.smd_dir = smd_dir
            
            self.logger.info(f"Preparing SMD simulation in {smd_dir}")
            
            # Use unified preparation logic
            prepared_files = self.prepare_files(smd_dir, use_pmf_system=inputs.get('use_pmf_system', True))
            
            self.logger.info(f"SMD preparation completed: {len(prepared_files)} files prepared")
            
            # Communicate preparation status to other modules
            if self.communicator:
                self.communicator.broadcast_status({
                    'module': 'smd',
                    'phase': 'preparation',
                    'status': 'completed',
                    'output_dir': str(smd_dir),
                    'files_prepared': len(prepared_files)
                })
            
            return {
                'smd_dir': str(smd_dir),
                'files_prepared': prepared_files,
                'pull_rate': self.pull_rate,
                'estimated_time': self._estimate_smd_time()
            }
            
        except Exception as e:
            self.logger.error(f"SMD preparation failed: {e}")
            raise SMDPreparationError(f"Failed to prepare SMD: {e}")
    
    def execute(self, inputs: Dict[str, Any]) -> ModuleResult:
        """Execute SMD simulation (modular interface)"""
        start_time = time.time()
        
        try:
            # Get SMD directory from preparation
            smd_dir = Path(inputs.get('smd_dir', inputs['output_dir']) + '/smd')
            self.smd_dir = smd_dir
            
            self.logger.info(f"Starting SMD execution in {smd_dir}")
            
            # Run SMD simulation using unified logic
            results = self.run_smd_simulation(smd_dir)
            
            # Calculate metrics
            metrics = self._calculate_smd_metrics(smd_dir)
            
            execution_time = time.time() - start_time
            
            result = ModuleResult(
                module_name=self.name,
                phase=ModulePhase.EXECUTION,
                status=ModuleStatus.COMPLETED,
                output_data=results,
                output_files=list(Path(f) for f in results.values() if isinstance(f, str) and Path(f).exists()),
                metrics=metrics,
                errors=[],
                warnings=[],
                execution_time=execution_time
            )
            
            self.logger.info(f"SMD execution completed successfully in {execution_time:.2f}s")
            
            # Communicate completion to other modules
            if self.communicator:
                self.communicator.notify_completion({
                    'module': 'smd',
                    'trajectory_file': results.get('trajectory_file'),
                    'pullx_file': results.get('pullx_file'),
                    'pullf_file': results.get('pullf_file'),
                    'suggested_windows': metrics.get('suggested_windows', 40)
                })
            
            return result
            
        except Exception as e:
            execution_time = time.time() - start_time
            error_msg = f"SMD execution failed: {str(e)}"
            self.logger.error(error_msg)
            
            return ModuleResult(
                module_name=self.name,
                phase=ModulePhase.EXECUTION,
                status=ModuleStatus.FAILED,
                output_data={},
                output_files=[],
                metrics={},
                errors=[error_msg],
                warnings=[],
                execution_time=execution_time
            )
    
    def analyze_results(self, execution_result: ModuleResult) -> ModuleResult:
        """Analyze SMD results (modular interface)"""
        start_time = time.time()
        
        try:
            if execution_result.status != ModuleStatus.COMPLETED:
                return execution_result
            
            smd_dir = Path(execution_result.output_data['smd_dir'])
            
            self.logger.info("Analyzing SMD results")
            
            # Analyze pull data
            analysis_results = self._analyze_pull_data(smd_dir)
            
            # Generate plots if matplotlib available
            plot_files = self._generate_analysis_plots(smd_dir, analysis_results)
            
            # Calculate suggested umbrella windows
            windows = self._calculate_umbrella_windows(analysis_results)
            
            execution_time = time.time() - start_time
            
            return ModuleResult(
                module_name=f"{self.name}_analysis",
                phase=ModulePhase.ANALYSIS,
                status=ModuleStatus.COMPLETED,
                output_data={
                    'analysis_results': analysis_results,
                    'suggested_windows': windows,
                    'max_force': analysis_results.get('max_force', 0),
                    'avg_force': analysis_results.get('avg_force', 0)
                },
                output_files=plot_files + execution_result.output_files,
                metrics={
                    'max_pull_force': analysis_results.get('max_force', 0),
                    'average_force': analysis_results.get('avg_force', 0),
                    'suggested_windows': len(windows),
                    'analysis_time': execution_time
                },
                errors=[],
                warnings=[],
                execution_time=execution_time
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            error_msg = f"SMD analysis failed: {str(e)}"
            self.logger.error(error_msg)
            
            return ModuleResult(
                module_name=f"{self.name}_analysis",
                phase=ModulePhase.ANALYSIS,
                status=ModuleStatus.FAILED,
                output_data={},
                output_files=[],
                metrics={},
                errors=[error_msg],
                warnings=[],
                execution_time=execution_time
            )
    
    def cleanup(self) -> None:
        """Cleanup temporary files"""
        try:
            if self.config.get('cleanup_temp_files', False):
                # Implementation for cleanup
                self.logger.info("SMD module cleanup completed")
        except Exception as e:
            self.logger.warning(f"SMD cleanup warning: {e}")
    
    def _estimate_smd_time(self) -> float:
        """Calculate actual SMD execution time based on user parameters"""
        if self.nsteps and self.dt:
            simulation_time_ps = self.nsteps * self.dt
            # Simple estimation: ~1 second real time per 1 ps simulation time (rough estimate)
            estimated_seconds = simulation_time_ps
            return estimated_seconds
        return 0.0  # Cannot estimate without user parameters
    
    def _calculate_smd_metrics(self, smd_dir: Path) -> Dict[str, float]:
        """Calculate SMD simulation metrics"""
        metrics = {}
        
        try:
            metrics['pull_rate'] = self.pull_rate
            metrics['total_steps'] = self.nsteps
            metrics['simulation_time_ns'] = self.nsteps * self.dt / 1000
            
            # Analyze pull files if available
            results_dir = smd_dir / "results"
            pullx_file = results_dir / 'smd_pullx.xvg'
            if pullx_file.exists():
                distances = self.parse_xvg_file(pullx_file)
                if distances:
                    metrics['max_distance'] = max(distances)
                    metrics['min_distance'] = min(distances)
                    metrics['distance_range'] = metrics['max_distance'] - metrics['min_distance']
                    metrics['suggested_windows'] = max(20, int(metrics['distance_range'] / 0.1))
            
            pullf_file = results_dir / 'smd_pullf.xvg'
            if pullf_file.exists():
                forces = self.parse_xvg_file(pullf_file)
                if forces:
                    metrics['max_force'] = max(forces)
                    metrics['avg_force'] = sum(forces) / len(forces)
        
        except Exception as e:
            self.logger.warning(f"Error calculating SMD metrics: {e}")
        
        return metrics
    
    def _analyze_pull_data(self, smd_dir: Path) -> Dict[str, Any]:
        """Analyze pull force and distance data"""
        analysis = {}
        results_dir = smd_dir / "results"
        
        # Parse distance data
        pullx_file = results_dir / 'smd_pullx.xvg'
        if pullx_file.exists():
            distances = self.parse_xvg_file(pullx_file)
            if distances:
                analysis['distances'] = distances
                analysis['total_distance'] = max(distances) - min(distances)
        
        # Parse force data
        pullf_file = results_dir / 'smd_pullf.xvg'
        if pullf_file.exists():
            forces = self.parse_xvg_file(pullf_file)
            if forces:
                analysis['forces'] = forces
                analysis['max_force'] = max(forces)
                analysis['avg_force'] = sum(forces) / len(forces)
        
        return analysis
    
    def _generate_analysis_plots(self, smd_dir: Path, analysis_results: Dict[str, Any]) -> List[Path]:
        """Generate analysis plots"""
        plot_files = []
        
        try:
            import matplotlib.pyplot as plt
            
            # Force vs time plot
            if 'forces' in analysis_results:
                plt.figure(figsize=(10, 6))
                plt.plot(analysis_results['forces'])
                plt.xlabel('Time Step')
                plt.ylabel('Pull Force (kJ/mol/nm)')
                plt.title('SMD Pull Force vs Time')
                plt.grid(True)
                
                force_plot = smd_dir / 'smd_force_analysis.png'
                plt.savefig(force_plot, dpi=300, bbox_inches='tight')
                plt.close()
                plot_files.append(force_plot)
            
            # Distance vs time plot
            if 'distances' in analysis_results:
                plt.figure(figsize=(10, 6))
                plt.plot(analysis_results['distances'])
                plt.xlabel('Time Step')
                plt.ylabel('Distance (nm)')
                plt.title('SMD Pull Distance vs Time')
                plt.grid(True)
                
                distance_plot = smd_dir / 'smd_distance_analysis.png'
                plt.savefig(distance_plot, dpi=300, bbox_inches='tight')
                plt.close()
                plot_files.append(distance_plot)
                
        except ImportError:
            self.logger.warning("Matplotlib not available for plotting")
        except Exception as e:
            self.logger.warning(f"Error generating plots: {e}")
        
        return plot_files
    
    def _calculate_umbrella_windows(self, analysis_results: Dict[str, Any]) -> List[float]:
        """Calculate suggested umbrella sampling windows"""
        windows = []
        
        if 'distances' in analysis_results:
            distances = analysis_results['distances']
            min_dist = min(distances)
            max_dist = max(distances)
            
            # Generate windows with 0.1 nm spacing
            window_spacing = 0.1
            num_windows = int((max_dist - min_dist) / window_spacing) + 1
            
            for i in range(num_windows):
                window_center = min_dist + i * window_spacing
                windows.append(round(window_center, 2))
        
        return windows[:50]  # Limit to 50 windows maximum


# Convenience functions for easy access
def create_smd_manager(pmf_system) -> SMDManager:
    """Create legacy SMD manager"""
    return SMDManager(pmf_system)


def create_smd_module(config: Optional[Dict[str, Any]] = None,
                     logger: Optional[PrismLogger] = None,
                     communicator: Optional[ModuleCommunicator] = None,
                     **user_params) -> SMDModule:
    """
    Create modern SMD module with full user customization

    Parameters:
    -----------
    config : dict, optional
        Configuration dictionary
    logger : PrismLogger, optional
        Logger instance
    communicator : ModuleCommunicator, optional
        Communication module
    **user_params :
        User-defined SMD parameters:
        - pull_rate : float (required) - Pulling rate in nm/ps
        - nsteps : int (required) - Number of simulation steps
        - dt : float (required) - Time step in ps
        - pull_k : float (required) - Force constant in kJ/mol/nm²
        - temperature : float (optional) - Temperature in K (default: 310.0)
        - pressure : float (optional) - Pressure in bar (default: 1.0)

    Examples:
    ---------
    >>> # Basic usage - user must provide all required parameters
    >>> smd = create_smd_module(
    ...     pull_rate=0.005,
    ...     nsteps=5000000,
    ...     dt=0.002,
    ...     pull_k=1000.0
    ... )
    >>>
    >>> # With custom temperature and pressure
    >>> smd = create_smd_module(
    ...     pull_rate=0.01,
    ...     nsteps=2000000,
    ...     dt=0.001,
    ...     pull_k=500.0,
    ...     temperature=300.0,
    ...     pressure=1.5
    ... )
    """
    return SMDModule(config, logger, communicator, **user_params)


# Legacy compatibility alias - SMDBuilder is now just an alias to SMDManager
SMDBuilder = SMDManager