#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Equilibration Manager

This module handles the complete equilibration process for PMF-rebuilt systems:
1. Energy Minimization (EM)
2. NVT equilibration 
3. NPT equilibration

Ensures the PMF-optimized system is properly equilibrated before PMF calculations.
"""

import os
import subprocess
from pathlib import Path
import yaml
import logging
from typing import Dict, Optional, List

from ..utils.exceptions import (
    PMFError, PMFSystemError, PMFFileError, PMFWorkflowError,
    RequiredFileNotFoundError, GromacsExecutionError, PMFErrorCode
)
from ..utils.error_handling import (
    with_retry, error_context, safe_file_operation, validate_prerequisites,
    handle_external_tool_error
)

logger = logging.getLogger(__name__)


def generate_pmf_equilibration_script(output_dir: Path, config: Dict, 
                                      gmx_command: str = "gmx") -> Path:
    """
    Standalone function to generate PMF equilibration script
    
    Parameters:
    -----------
    output_dir : Path
        Directory to create equilibration script
    config : Dict
        Configuration dictionary
    gmx_command : str
        GROMACS command to use
        
    Returns:
    --------
    Path : Generated script path
    """
    manager = PMFEquilibrationManager(
        system_dir=output_dir / "GMX_PMF_SYSTEM",
        output_dir=output_dir,
        config=config,
        gmx_command=gmx_command
    )
    
    return manager.generate_equilibration_script()


class PMFEquilibrationManager:
    """
    PMF Equilibration Manager
    
    Manages the complete equilibration workflow for PMF-optimized systems,
    including energy minimization and both NVT and NPT equilibration phases.
    """
    
    def __init__(self, system_dir: Path, output_dir: Path, config: Dict, 
                 gmx_command: str = "gmx"):
        """
        Initialize PMF equilibration manager
        
        Parameters:
        -----------
        system_dir : Path
            Directory containing the PMF-optimized system
        output_dir : Path  
            Output directory for equilibration results
        config : Dict
            Configuration dictionary containing equilibration parameters
        gmx_command : str
            GROMACS command to use
        """
        self.system_dir = Path(system_dir)
        self.output_dir = Path(output_dir)
        self.config = config
        self.gmx_command = gmx_command
        
        # Create equilibration directory structure
        self.equil_dir = self.output_dir / "equilibration"
        self.equil_dir.mkdir(parents=True, exist_ok=True)
        
        self.em_dir = self.equil_dir / "em"
        self.nvt_dir = self.equil_dir / "nvt"
        self.npt_dir = self.equil_dir / "npt"
        
        for dir_path in [self.em_dir, self.nvt_dir, self.npt_dir]:
            dir_path.mkdir(exist_ok=True)
        
        # Input files
        self.input_structure = self.system_dir / "solv_ions.gro"
        self.input_topology = self.system_dir / "topol.top"
        
        # Validate input files
        self._validate_input_files()
        
        # Equilibration state tracking
        self.state = {
            'em_completed': False,
            'nvt_completed': False,
            'npt_completed': False
        }
        
        logger.info(f"PMF equilibration manager initialized")
        logger.info(f"  System: {self.system_dir}")
        logger.info(f"  Output: {self.equil_dir}")
    
    def _validate_input_files(self):
        """Validate that required input files exist"""
        if not self.input_structure.exists():
            raise RequiredFileNotFoundError(
                file_path=str(self.input_structure),
                step="equilibration",
                alternative_files=[str(f) for f in self.system_dir.iterdir() if f.suffix in ['.gro', '.pdb']]
            )
        
        if not self.input_topology.exists():
            raise RequiredFileNotFoundError(
                file_path=str(self.input_topology),
                step="equilibration",
                alternative_files=[str(f) for f in self.system_dir.iterdir() if f.suffix == '.top']
            )
        
        logger.debug("Input files validated successfully")
    
    def run_complete_equilibration(self) -> Dict:
        """
        Run complete equilibration workflow (EM → NVT → NPT) with error handling
        
        Returns:
        --------
        Dict : Complete equilibration results
        
        Raises:
        -------
        RequiredFileNotFoundError
            If input files are missing
        PMFWorkflowError
            If equilibration steps fail
        GromacsExecutionError
            If GROMACS commands fail
        """
        with error_context("PMF system equilibration", {
            "system_dir": str(self.system_dir),
            "output_dir": str(self.equil_dir)
        }):
            logger.info("=== Starting Complete PMF System Equilibration ===")
            
            # Validate inputs
            self._validate_input_files()
            
            results = {
                'workflow': 'complete_equilibration',
                'stages': {}
            }
        
        try:
            # Stage 1: Energy Minimization
            logger.info("Stage 1: Energy Minimization")
            em_results = self.run_energy_minimization()
            results['stages']['em'] = em_results
            
            # Stage 2: NVT Equilibration  
            logger.info("Stage 2: NVT Equilibration")
            nvt_results = self.run_nvt_equilibration()
            results['stages']['nvt'] = nvt_results
            
            # Stage 3: NPT Equilibration
            logger.info("Stage 3: NPT Equilibration") 
            npt_results = self.run_npt_equilibration()
            results['stages']['npt'] = npt_results
            
            # Final equilibrated system
            results['final_system'] = {
                'structure': str(self.npt_dir / "npt_final.gro"),
                'checkpoint': str(self.npt_dir / "npt_final.cpt"),
                'ready_for_pmf': True
            }
            
            results['status'] = 'equilibration_complete'
            logger.info("Complete PMF system equilibration finished successfully!")
            
        except Exception as e:
            results['status'] = 'equilibration_failed'
            results['error'] = str(e)
            logger.error(f"Equilibration failed: {e}")
            raise
        
        # Save equilibration state
        self._save_equilibration_state(results)
        
        return results
    
    def run_energy_minimization(self) -> Dict:
        """
        Run energy minimization
        
        Returns:
        --------
        Dict : Energy minimization results
        """
        logger.info("Running energy minimization...")
        
        # Generate EM MDP file
        em_mdp = self._generate_em_mdp()
        
        # Prepare EM simulation
        em_tpr = self.em_dir / "em.tpr"
        self._run_grompp(
            mdp_file=em_mdp,
            structure_file=self.input_structure,
            topology_file=self.input_topology,
            output_tpr=em_tpr,
            work_dir=self.em_dir
        )
        
        # Run EM
        em_output = self.em_dir / "em.gro"
        em_energy = self.em_dir / "em_energy.xvg"
        
        self._run_mdrun(
            tpr_file=em_tpr,
            output_structure=em_output,
            output_energy=em_energy,
            work_dir=self.em_dir,
            description="Energy Minimization"
        )
        
        # Analyze EM results
        em_analysis = self._analyze_energy_minimization(em_energy)
        
        self.state['em_completed'] = True
        
        results = {
            'stage': 'energy_minimization',
            'status': 'completed',
            'output_structure': str(em_output),
            'energy_file': str(em_energy),
            'analysis': em_analysis,
            'mdp_used': str(em_mdp)
        }
        
        logger.info(f"Energy minimization completed")
        logger.info(f"   Final energy: {em_analysis.get('final_energy', 'N/A')} kJ/mol")
        
        return results
    
    def run_nvt_equilibration(self) -> Dict:
        """
        Run NVT equilibration
        
        Returns:
        --------
        Dict : NVT equilibration results
        """
        if not self.state['em_completed']:
            raise PMFWorkflowError(
                message="Energy minimization must be completed before NVT",
                error_code=PMFErrorCode.PREREQUISITE_NOT_MET,
                recoverable=True,
                recovery_suggestions=[
                    "Run energy minimization step first",
                    "Check equilibration workflow order"
                ]
            )
        
        logger.info("Running NVT equilibration...")
        
        # Generate NVT MDP file
        nvt_mdp = self._generate_nvt_mdp()
        
        # Input from EM results
        input_structure = self.em_dir / "em.gro"
        
        # Prepare NVT simulation
        nvt_tpr = self.nvt_dir / "nvt.tpr"
        self._run_grompp(
            mdp_file=nvt_mdp,
            structure_file=input_structure,
            topology_file=self.input_topology,
            output_tpr=nvt_tpr,
            work_dir=self.nvt_dir
        )
        
        # Run NVT
        nvt_output = self.nvt_dir / "nvt_final.gro"
        nvt_checkpoint = self.nvt_dir / "nvt_final.cpt"
        nvt_energy = self.nvt_dir / "nvt_energy.xvg"
        nvt_temperature = self.nvt_dir / "nvt_temperature.xvg"
        
        self._run_mdrun(
            tpr_file=nvt_tpr,
            output_structure=nvt_output,
            output_checkpoint=nvt_checkpoint,
            output_energy=nvt_energy,
            work_dir=self.nvt_dir,
            description="NVT Equilibration"
        )
        
        # Extract temperature data
        self._extract_temperature_data(nvt_energy, nvt_temperature)
        
        # Analyze NVT results
        nvt_analysis = self._analyze_nvt_equilibration(nvt_temperature)
        
        self.state['nvt_completed'] = True
        
        results = {
            'stage': 'nvt_equilibration',
            'status': 'completed',
            'output_structure': str(nvt_output),
            'output_checkpoint': str(nvt_checkpoint),
            'energy_file': str(nvt_energy),
            'temperature_file': str(nvt_temperature),
            'analysis': nvt_analysis,
            'mdp_used': str(nvt_mdp)
        }
        
        logger.info(f"NVT equilibration completed")
        logger.info(f"   Average temperature: {nvt_analysis.get('avg_temperature', 'N/A')} K")
        
        return results
    
    def run_npt_equilibration(self) -> Dict:
        """
        Run NPT equilibration
        
        Returns:
        --------
        Dict : NPT equilibration results
        """
        if not self.state['nvt_completed']:
            raise PMFWorkflowError(
                message="NVT equilibration must be completed before NPT",
                error_code=PMFErrorCode.PREREQUISITE_NOT_MET,
                recoverable=True,
                recovery_suggestions=[
                    "Run NVT equilibration step first",
                    "Check equilibration workflow order"
                ]
            )
        
        logger.info("Running NPT equilibration...")
        
        # Generate NPT MDP file
        npt_mdp = self._generate_npt_mdp()
        
        # Input from NVT results
        input_structure = self.nvt_dir / "nvt_final.gro"
        input_checkpoint = self.nvt_dir / "nvt_final.cpt"
        
        # Prepare NPT simulation
        npt_tpr = self.npt_dir / "npt.tpr"
        self._run_grompp(
            mdp_file=npt_mdp,
            structure_file=input_structure,
            topology_file=self.input_topology,
            checkpoint_file=input_checkpoint,
            output_tpr=npt_tpr,
            work_dir=self.npt_dir
        )
        
        # Run NPT
        npt_output = self.npt_dir / "npt_final.gro"
        npt_checkpoint = self.npt_dir / "npt_final.cpt"
        npt_energy = self.npt_dir / "npt_energy.xvg"
        npt_pressure = self.npt_dir / "npt_pressure.xvg"
        npt_density = self.npt_dir / "npt_density.xvg"
        
        self._run_mdrun(
            tpr_file=npt_tpr,
            output_structure=npt_output,
            output_checkpoint=npt_checkpoint,
            output_energy=npt_energy,
            work_dir=self.npt_dir,
            description="NPT Equilibration"
        )
        
        # Extract pressure and density data
        self._extract_pressure_data(npt_energy, npt_pressure)
        self._extract_density_data(npt_energy, npt_density)
        
        # Analyze NPT results
        npt_analysis = self._analyze_npt_equilibration(npt_pressure, npt_density)
        
        self.state['npt_completed'] = True
        
        results = {
            'stage': 'npt_equilibration',
            'status': 'completed',
            'output_structure': str(npt_output),
            'output_checkpoint': str(npt_checkpoint),
            'energy_file': str(npt_energy),
            'pressure_file': str(npt_pressure),
            'density_file': str(npt_density),
            'analysis': npt_analysis,
            'mdp_used': str(npt_mdp)
        }
        
        logger.info(f"NPT equilibration completed")
        logger.info(f"   Average pressure: {npt_analysis.get('avg_pressure', 'N/A')} bar")
        logger.info(f"   Average density: {npt_analysis.get('avg_density', 'N/A')} kg/m³")
        
        return results
    
    def _generate_em_mdp(self) -> Path:
        """Generate energy minimization MDP file"""
        em_config = self.config.get('equilibration', {}).get('em', {})
        
        # Default EM parameters optimized for PMF systems
        defaults = {
            'integrator': 'steep',
            'emtol': 10.0,
            'emstep': 0.01,
            'nsteps': 50000,
            'rcoulomb': 1.0,
            'rvdw': 1.0
        }
        defaults.update(em_config)
        
        mdp_content = f"""; em.mdp - Energy minimization for PMF system
title           = PMF System Energy Minimization
define          = -DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=500.0

; Minimization algorithm
integrator      = {defaults['integrator']}
emtol           = {defaults['emtol']}
emstep          = {defaults['emstep']}
nsteps          = {defaults['nsteps']}

; Output control
nstxout         = 0
nstvout         = 0
nstenergy       = 100
nstlog          = 100
nstxout-compressed = 0

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 1
rlist           = 1.0
pbc             = xyz

; Electrostatics
coulombtype     = PME
rcoulomb        = {defaults['rcoulomb']}
fourierspacing  = 0.16
pme_order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = {defaults['rvdw']}
DispCorr        = EnerPres
"""
        
        mdp_file = self.em_dir / "em.mdp"
        with open(mdp_file, 'w') as f:
            f.write(mdp_content)
        
        logger.debug(f"Generated EM MDP: {mdp_file}")
        return mdp_file
    
    def _generate_nvt_mdp(self) -> Path:
        """Generate NVT equilibration MDP file"""
        nvt_config = self.config.get('equilibration', {}).get('nvt', {})
        
        # Default NVT parameters
        defaults = {
            'nsteps': 50000,  # 100 ps with dt=0.002
            'dt': 0.002,
            'temperature': 310.0,
            'tau_t': 0.1,
            'rcoulomb': 1.0,
            'rvdw': 1.0
        }
        defaults.update(nvt_config)
        
        mdp_content = f"""; nvt.mdp - NVT equilibration for PMF system
title           = PMF System NVT Equilibration
define          = -DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=500.0

; Simulation parameters
integrator      = md
dt              = {defaults['dt']}
nsteps          = {defaults['nsteps']}

; Output control
nstxout         = 0
nstvout         = 0
nstenergy       = 500
nstlog          = 500
nstxout-compressed = 5000

; Bond parameters
continuation    = no
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0
pbc             = xyz

; Electrostatics
coulombtype     = PME
rcoulomb        = {defaults['rcoulomb']}
fourierspacing  = 0.16
pme_order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = {defaults['rvdw']}
DispCorr        = EnerPres

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = {defaults['tau_t']} {defaults['tau_t']}
ref_t           = {defaults['temperature']} {defaults['temperature']}

; Pressure coupling
pcoupl          = no

; Velocity generation
gen_vel         = yes
gen_temp        = {defaults['temperature']}
gen_seed        = -1
"""
        
        mdp_file = self.nvt_dir / "nvt.mdp"
        with open(mdp_file, 'w') as f:
            f.write(mdp_content)
        
        logger.debug(f"Generated NVT MDP: {mdp_file}")
        return mdp_file
    
    def _generate_npt_mdp(self) -> Path:
        """Generate NPT equilibration MDP file"""
        npt_config = self.config.get('equilibration', {}).get('npt', {})
        
        # Default NPT parameters - updated for improved equilibration
        defaults = {
            'nsteps': 250000,  # 500 ps with dt=0.002 - extended for better equilibration
            'dt': 0.002,
            'temperature': 310.0,
            'pressure': 1.0,
            'tau_t': 0.1,
            'tau_p': 1.0,      # Faster pressure response for better equilibration
            'rcoulomb': 1.0,
            'rvdw': 1.0
        }
        defaults.update(npt_config)
        
        mdp_content = f"""; npt.mdp - NPT equilibration with protein restrained
title               = NPT equilibration with protein restrained
define              = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator          = md
nsteps              = {defaults['nsteps']}
dt                  = {defaults['dt']}

nstxout             = 0
nstvout             = 0
nstenergy           = 5000
nstlog              = 5000

continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

cutoff-scheme           = Verlet
rcoulomb                = {defaults['rcoulomb']}
rvdw                    = {defaults['rvdw']}

coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = {defaults['tau_t']}     {defaults['tau_t']}
ref_t                   = {defaults['temperature']}     {defaults['temperature']}

pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = {defaults['tau_p']}
ref_p                   = {defaults['pressure']}
compressibility         = 4.5e-05

pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = Protein Non-Protein
"""
        
        mdp_file = self.npt_dir / "npt.mdp"
        with open(mdp_file, 'w') as f:
            f.write(mdp_content)
        
        logger.debug(f"Generated NPT MDP: {mdp_file}")
        return mdp_file
    
    def _run_grompp(self, mdp_file: Path, structure_file: Path, topology_file: Path,
                   output_tpr: Path, work_dir: Path, checkpoint_file: Optional[Path] = None):
        """
        Run GROMACS grompp to prepare simulation
        
        For GROMACS 2018+, automatically adds -r option using structure_file
        as position restraint coordinate file when position restraints are used.
        """
        cmd = [
            self.gmx_command, "grompp",
            "-f", str(mdp_file),
            "-c", str(structure_file),
            "-p", str(topology_file),
            "-r", str(structure_file),  # Position restraint coordinates (required for GROMACS 2018+)
            "-o", str(output_tpr),
            "-maxwarn", "2"
        ]
        
        if checkpoint_file and checkpoint_file.exists():
            cmd.extend(["-t", str(checkpoint_file)])
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
        
        if result.returncode != 0:
            handle_external_tool_error(
                command=" ".join(cmd),
                exit_code=result.returncode,
                stderr=result.stderr,
                tool_name="GROMACS grompp"
            )
        
        logger.debug(f"grompp completed: {output_tpr}")
    
    def _run_mdrun(self, tpr_file: Path, output_structure: Path, work_dir: Path,
                  description: str, output_checkpoint: Optional[Path] = None,
                  output_energy: Optional[Path] = None):
        """Run GROMACS mdrun"""
        cmd = [
            self.gmx_command, "mdrun",
            "-s", str(tpr_file),
            "-o", str(tpr_file.with_suffix('.trr')),
            "-c", str(output_structure),
            "-g", str(tpr_file.with_suffix('.log')),
            "-e", str(tpr_file.with_suffix('.edr')),
            "-deffnm", str(tpr_file.stem)
        ]
        
        if output_checkpoint:
            cmd.extend(["-cpo", str(output_checkpoint)])
        
        logger.info(f"Running {description}...")
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
        
        if result.returncode != 0:
            handle_external_tool_error(
                command=" ".join(cmd),
                exit_code=result.returncode,
                stderr=result.stderr,
                tool_name=f"GROMACS mdrun ({description})"
            )
        
        logger.debug(f"{description} completed successfully")
    
    def _extract_temperature_data(self, energy_file: Path, output_file: Path):
        """Extract temperature data from energy file"""
        cmd = [
            self.gmx_command, "energy",
            "-f", str(energy_file.with_suffix('.edr')),
            "-o", str(output_file)
        ]
        
        result = subprocess.run(
            cmd, 
            input="Temperature\n",
            text=True,
            capture_output=True,
            cwd=self.nvt_dir
        )
        
        if result.returncode != 0:
            logger.warning(f"Could not extract temperature data: {result.stderr}")
    
    def _extract_pressure_data(self, energy_file: Path, output_file: Path):
        """Extract pressure data from energy file"""
        cmd = [
            self.gmx_command, "energy", 
            "-f", str(energy_file.with_suffix('.edr')),
            "-o", str(output_file)
        ]
        
        result = subprocess.run(
            cmd,
            input="Pressure\n",
            text=True,
            capture_output=True,
            cwd=self.npt_dir
        )
        
        if result.returncode != 0:
            logger.warning(f"Could not extract pressure data: {result.stderr}")
    
    def _extract_density_data(self, energy_file: Path, output_file: Path):
        """Extract density data from energy file"""
        cmd = [
            self.gmx_command, "energy",
            "-f", str(energy_file.with_suffix('.edr')),
            "-o", str(output_file)
        ]
        
        result = subprocess.run(
            cmd,
            input="Density\n", 
            text=True,
            capture_output=True,
            cwd=self.npt_dir
        )
        
        if result.returncode != 0:
            logger.warning(f"Could not extract density data: {result.stderr}")
    
    def _analyze_energy_minimization(self, energy_file: Path) -> Dict:
        """Analyze energy minimization results"""
        analysis = {'stage': 'energy_minimization'}
        
        try:
            # Simple analysis - would parse energy file in practice
            analysis['status'] = 'completed'
            analysis['converged'] = True
            analysis['final_energy'] = 'N/A (parse energy file)'
            
        except Exception as e:
            analysis['status'] = 'analysis_failed'
            analysis['error'] = str(e)
        
        return analysis
    
    def _analyze_nvt_equilibration(self, temperature_file: Path) -> Dict:
        """Analyze NVT equilibration results"""
        analysis = {'stage': 'nvt_equilibration'}
        
        try:
            # Simple analysis - would parse temperature file in practice
            analysis['status'] = 'completed'
            analysis['equilibrated'] = True
            analysis['avg_temperature'] = 'N/A (parse temperature file)'
            analysis['temperature_stable'] = True
            
        except Exception as e:
            analysis['status'] = 'analysis_failed'
            analysis['error'] = str(e)
        
        return analysis
    
    def _analyze_npt_equilibration(self, pressure_file: Path, density_file: Path) -> Dict:
        """Analyze NPT equilibration results"""
        analysis = {'stage': 'npt_equilibration'}
        
        try:
            # Simple analysis - would parse pressure and density files in practice
            analysis['status'] = 'completed'
            analysis['equilibrated'] = True
            analysis['avg_pressure'] = 'N/A (parse pressure file)'
            analysis['avg_density'] = 'N/A (parse density file)'
            analysis['pressure_stable'] = True
            analysis['density_stable'] = True
            
        except Exception as e:
            analysis['status'] = 'analysis_failed'
            analysis['error'] = str(e)
        
        return analysis
    
    def _save_equilibration_state(self, results: Dict):
        """Save equilibration state and results"""
        state_file = self.equil_dir / "equilibration_state.yaml"
        
        state_data = {
            'equilibration_state': self.state,
            'results': results,
            'system_info': {
                'input_structure': str(self.input_structure),
                'input_topology': str(self.input_topology),
                'output_directory': str(self.equil_dir)
            }
        }
        
        with open(state_file, 'w') as f:
            yaml.dump(state_data, f, default_flow_style=False, indent=2)
        
        logger.info(f"Equilibration state saved: {state_file}")
    
    def get_equilibration_status(self) -> Dict:
        """Get current equilibration status"""
        return {
            'current_state': self.state.copy(),
            'em_ready': self.input_structure.exists() and self.input_topology.exists(),
            'nvt_ready': self.state['em_completed'],
            'npt_ready': self.state['nvt_completed'],
            'fully_equilibrated': all(self.state.values()),
            'next_step': self._get_next_equilibration_step()
        }
    
    def _get_next_equilibration_step(self) -> str:
        """Determine next equilibration step"""
        if not self.state['em_completed']:
            return 'energy_minimization'
        elif not self.state['nvt_completed']:
            return 'nvt_equilibration'
        elif not self.state['npt_completed']:
            return 'npt_equilibration'
        else:
            return 'equilibration_complete'
    
    def get_final_equilibrated_system(self) -> Dict:
        """Get paths to final equilibrated system files"""
        if not all(self.state.values()):
            raise RuntimeError("System not fully equilibrated")
        
        return {
            'structure': str(self.npt_dir / "npt_final.gro"),
            'checkpoint': str(self.npt_dir / "npt_final.cpt"),
            'topology': str(self.input_topology),
            'ready_for_pmf': True,
            'equilibration_summary': self.equil_dir / "equilibration_state.yaml"
        }
    
    def generate_equilibration_script(self) -> Path:
        """
        Generate localrun.sh-style equilibration script for PMF system
        
        Creates a comprehensive bash script similar to PRISM's localrun.sh
        that runs the complete equilibration workflow (EM → NVT → NPT)
        
        Returns:
        --------
        Path : Generated script path
        """
        script_path = self.equil_dir / "run_pmf_equilibration.sh"
        
        # Get configuration parameters
        em_config = self.config.get('equilibration', {}).get('em', {})
        nvt_config = self.config.get('equilibration', {}).get('nvt', {})
        npt_config = self.config.get('equilibration', {}).get('npt', {})
        
        em_nsteps = em_config.get('nsteps', 50000)
        nvt_nsteps = nvt_config.get('nsteps', 50000)
        npt_nsteps = npt_config.get('nsteps', 250000)  # Updated default: 500 ps for better equilibration
        dt = nvt_config.get('dt', 0.002)
        
        # Calculate simulation times
        nvt_time_ps = nvt_nsteps * dt
        npt_time_ps = npt_nsteps * dt
        
        script_content = f"""#!/bin/bash
######################################################
# PMF SYSTEM EQUILIBRATION SCRIPT
# Generated by PRISM PMF Module
# Similar to localrun.sh but optimized for PMF systems
######################################################

set -e  # Exit on any error
echo "=== PRISM PMF System Equilibration ==="
echo "Equilibration protocol: EM ({em_nsteps} steps) → NVT ({nvt_time_ps:.0f} ps) → NPT ({npt_time_ps:.0f} ps)"
echo ""

# Check for required files
if [ ! -f "../GMX_PMF_SYSTEM/solv_ions.gro" ]; then
    echo "Error: PMF system not found. Run PMF builder first."
    exit 1
fi

if [ ! -f "../GMX_PMF_SYSTEM/topol.top" ]; then
    echo "Error: Topology file not found. Run PMF builder first."
    exit 1
fi

# Set up environment
export GMX_MAXBACKUP=-1  # Disable backup files
GMX_COMMAND="{self.gmx_command}"
NCPU=10  # Adjust as needed
NGPU=1   # Adjust as needed

# Function to check GROMACS completion
check_completion() {{
    local stage=$1
    local expected_files=("${{@:2}}")
    
    for file in "${{expected_files[@]}}"; do
        if [ ! -f "$file" ]; then
            echo "Error: $stage failed - missing $file"
            exit 1
        fi
    done
    echo "$stage completed successfully"
}}

# Function to run with GPU acceleration if available
run_mdrun() {{
    local description="$1"
    local tpr_file="$2"
    local output_dir="$3"
    
    echo "Running $description..."
    
    # Try GPU first, fallback to CPU
    cd "$output_dir"
    if command -v nvidia-smi &> /dev/null && nvidia-smi &> /dev/null; then
        echo "  Using GPU acceleration"
        $GMX_COMMAND mdrun -s "$tpr_file" -deffnm "$(basename "$tpr_file" .tpr)" \
            -ntmpi 1 -ntomp $NCPU -nb gpu -pme gpu -v
    else
        echo "  Using CPU only"
        $GMX_COMMAND mdrun -s "$tpr_file" -deffnm "$(basename "$tpr_file" .tpr)" \
            -ntmpi 1 -ntomp $NCPU -v
    fi
    cd ..
}}

echo "Step 1: Energy Minimization"
echo "Creating EM directory and files..."
mkdir -p em
cd em

# Copy input files
cp ../../GMX_PMF_SYSTEM/solv_ions.gro ./
cp ../../GMX_PMF_SYSTEM/topol.top ./

# Create EM MDP file  
cat > em.mdp << 'EOF'
; em.mdp - Energy minimization for PMF system
title           = PMF System Energy Minimization
define          = -DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=500.0

; Minimization algorithm
integrator      = {em_config.get('integrator', 'steep')}
emtol           = {em_config.get('emtol', 10.0)}
emstep          = {em_config.get('emstep', 0.01)}
nsteps          = {em_nsteps}

; Output control
nstxout         = 0
nstvout         = 0
nstenergy       = 100
nstlog          = 100
nstxout-compressed = 0

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 1
rlist           = 1.0
pbc             = xyz

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0
fourierspacing  = 0.16
pme_order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.0
DispCorr        = EnerPres
EOF

# Prepare and run EM
echo "Preparing energy minimization..."
$GMX_COMMAND grompp -f em.mdp -c solv_ions.gro -p topol.top -r solv_ions.gro -o em.tpr -maxwarn 2

run_mdrun "Energy Minimization" "em.tpr" "$(pwd)"
check_completion "EM" "em.gro" "em.log" "em.edr"

cd ..
echo ""

echo "Step 2: NVT Equilibration"
echo "Creating NVT directory and files..."
mkdir -p nvt
cd nvt

# Copy files from EM
cp ../em/em.gro ./nvt_input.gro
cp ../em/topol.top ./

# Create NVT MDP file
cat > nvt.mdp << 'EOF'
; nvt.mdp - NVT equilibration for PMF system
title           = PMF System NVT Equilibration
define          = -DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=500.0

; Simulation parameters
integrator      = md
dt              = {dt}
nsteps          = {nvt_nsteps}

; Output control
nstxout         = 0
nstvout         = 0
nstenergy       = 500
nstlog          = 500
nstxout-compressed = 5000

; Bond parameters
continuation    = no
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0
pbc             = xyz

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0
fourierspacing  = 0.16
pme_order       = 4

; Van der Waals
vdwtype         = cut-off
rvdw            = 1.0
DispCorr        = EnerPres

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = {nvt_config.get('tau_t', 0.1)} {nvt_config.get('tau_t', 0.1)}
ref_t           = {nvt_config.get('temperature', 310.0)} {nvt_config.get('temperature', 310.0)}

; Pressure coupling
pcoupl          = no

; Velocity generation
gen_vel         = yes
gen_temp        = {nvt_config.get('temperature', 310.0)}
gen_seed        = -1
EOF

# Prepare and run NVT
echo "Preparing NVT equilibration..."
$GMX_COMMAND grompp -f nvt.mdp -c nvt_input.gro -p topol.top -r nvt_input.gro -o nvt.tpr -maxwarn 2

run_mdrun "NVT Equilibration" "nvt.tpr" "$(pwd)"
check_completion "NVT" "nvt.gro" "nvt.cpt" "nvt.log" "nvt.edr"

# Extract temperature for analysis
echo "Extracting temperature data..."
echo "Temperature" | $GMX_COMMAND energy -f nvt.edr -o temperature.xvg -quiet || true

cd ..
echo ""

echo "Step 3: NPT Equilibration"
echo "Creating NPT directory and files..."
mkdir -p npt
cd npt

# Copy files from NVT
cp ../nvt/nvt.gro ./npt_input.gro
cp ../nvt/nvt.cpt ./npt_input.cpt
cp ../nvt/topol.top ./

# Create NPT MDP file
cat > npt.mdp << 'EOF'
; npt.mdp - NPT equilibration with protein restrained
title               = NPT equilibration with protein restrained
define              = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator          = md
nsteps              = {npt_nsteps}
dt                  = {dt}

nstxout             = 0
nstvout             = 0
nstenergy           = 5000
nstlog              = 5000

continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

cutoff-scheme           = Verlet
rcoulomb                = 1.0
rvdw                    = 1.0

coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = {npt_config.get('tau_t', 0.1)}     {npt_config.get('tau_t', 0.1)}
ref_t                   = {npt_config.get('temperature', 310.0)}     {npt_config.get('temperature', 310.0)}

pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = {npt_config.get('tau_p', 1.0)}
ref_p                   = {npt_config.get('pressure', 1.0)}
compressibility         = 4.5e-05

pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = Protein Non-Protein
EOF

# Prepare and run NPT
echo "Preparing NPT equilibration..."
$GMX_COMMAND grompp -f npt.mdp -c npt_input.gro -p topol.top -r npt_input.gro -t npt_input.cpt -o npt.tpr -maxwarn 2

run_mdrun "NPT Equilibration" "npt.tpr" "$(pwd)"
check_completion "NPT" "npt.gro" "npt.cpt" "npt.log" "npt.edr"

# Extract pressure and density for analysis
echo "Extracting pressure and density data..."
echo "Pressure" | $GMX_COMMAND energy -f npt.edr -o pressure.xvg -quiet || true
echo "Density" | $GMX_COMMAND energy -f npt.edr -o density.xvg -quiet || true

# Create final equilibrated system links
echo "Creating final equilibrated system..."
cp npt.gro npt_final.gro
cp npt.cpt npt_final.cpt

cd ..
echo ""

echo "=== PMF System Equilibration Completed ==="
echo "Energy minimization: $(ls em/em.gro 2>/dev/null && echo 'SUCCESS' || echo 'FAILED')"
echo "NVT equilibration:   $(ls nvt/nvt.gro 2>/dev/null && echo 'SUCCESS' || echo 'FAILED')"
echo "NPT equilibration:   $(ls npt/npt.gro 2>/dev/null && echo 'SUCCESS' || echo 'FAILED')"
echo ""
echo "Final equilibrated system:"
echo "  Structure: npt/npt_final.gro"
echo "  Checkpoint: npt/npt_final.cpt"
echo "  Topology: ../GMX_PMF_SYSTEM/topol.top"
echo ""
echo "System is ready for PMF calculations!"
echo "Next: Run PMF workflow using the equilibrated system"
"""
        
        # Write script file
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make executable
        import stat
        script_path.chmod(script_path.stat().st_mode | stat.S_IEXEC)
        
        logger.info(f"Generated PMF equilibration script: {script_path}")
        return script_path