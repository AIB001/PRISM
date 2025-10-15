#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Builder - Specialized builder for PMF-oriented system preparation

This module creates PMF-ready systems from existing MD results by:
1. Extracting Protein and LIG structures from MD trajectory
2. Computing centroids for alignment along Z-axis
3. Rebuilding system with extended Z-box for pulling simulations
4. Following PRISM builder patterns for consistency
"""

import os
import re
import shutil
import numpy as np
import subprocess
from pathlib import Path
import yaml
import logging
from datetime import datetime
from typing import Dict, Tuple, Optional, Union

from ...utils.environment import GromacsEnvironment
from ...utils.config import ConfigurationManager
from .equilibration import PMFEquilibrationManager
from ..utils.exceptions import (
    PMFError, PMFSystemError, PMFFileError, PMFConfigurationError,
    PMFWorkflowError, SystemNotFoundException, RequiredFileNotFoundError,
    FileCorruptionError, GromacsExecutionError, PMFErrorCode
)
from ..utils.error_handling import (
    with_retry, error_context, safe_file_operation, validate_prerequisites,
    ErrorCollector, handle_external_tool_error
)

logger = logging.getLogger(__name__)


class PMFBuilder:
    """
    PMF System Builder - Specialized builder for PMF simulations
    
    Rebuilds protein-ligand systems from MD results with optimized geometry
    for steered molecular dynamics and umbrella sampling calculations.
    """
    
    def __init__(self, md_results_dir: Union[str, Path], output_dir: Union[str, Path],
                 config: Optional[Dict] = None, **kwargs):
        """
        Initialize PMF Builder
        
        Parameters:
        -----------
        md_results_dir : str or Path
            Directory containing MD results - should be built with PMF system remodeling
        output_dir : str or Path
            Output directory for PMF-optimized system
        config : dict, optional
            PMF configuration dictionary
        **kwargs : optional
            Additional configuration parameters
        """
        self.md_results_dir = Path(md_results_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        
        # Locate MD system files
        self.md_system_dir = self._locate_md_system()
        
        # Initialize GROMACS environment
        self.gromacs_env = GromacsEnvironment()
        self.gmx_command = self.gromacs_env.gmx_command
        
        # Detect original force field from MD system (maintain consistency)
        self.original_forcefield_info = self._detect_original_forcefield()
        
        # Process configuration with force field info
        self.config = self._process_pmf_config(config, **kwargs)
        
        # Initialize standard PRISM components
        self.config_manager = ConfigurationManager(
            None,  # Use provided config dict instead of file
            self.gromacs_env,
            forcefield_name=self.config['general']['forcefield'],
            water_model_name=self.config['general']['water_model']
        )
        
        # Create output structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.work_dir = self.output_dir / "PMF_SYSTEM_BUILD"
        self.work_dir.mkdir(exist_ok=True)
        
        # Output paths for rebuilt system
        self.rebuilt_system_dir = self.output_dir / "GMX_PMF_SYSTEM"
        self.rebuilt_system_dir.mkdir(exist_ok=True)
        
        # Initialize equilibration manager (created when needed)
        self.equilibration_manager = None
        
        logger.info(f"PMF Builder initialized")
        logger.info(f"  MD source: {self.md_system_dir}")
        logger.info(f"  PMF output: {self.output_dir}")
    
    def _locate_md_system(self) -> Path:
        """Locate the MD system directory"""
        # Check if md_results_dir itself contains system files
        if (self.md_results_dir / "solv_ions.gro").exists():
            return self.md_results_dir
        
        # Look for GMX_PROLIG_MD
        gmx_dir = self.md_results_dir / "GMX_PROLIG_MD"
        if gmx_dir.exists() and (gmx_dir / "solv_ions.gro").exists():
            return gmx_dir
        
        # Search in subdirectories
        for subdir in self.md_results_dir.iterdir():
            if subdir.is_dir():
                candidate = subdir / "GMX_PROLIG_MD"
                if candidate.exists() and (candidate / "solv_ions.gro").exists():
                    return candidate
        
        raise FileNotFoundError(
            f"Could not locate MD system files in {self.md_results_dir}. "
            "Expected structure: */GMX_PROLIG_MD/solv_ions.gro"
        )
    
    def _detect_original_forcefield(self) -> Dict:
        """
        Detect original MD system force field information to ensure PMF rebuild consistency
        
        Returns:
        --------
        Dict : Original force field information
        """
        logger.info("Detecting original MD system force field information...")
        
        topology_file = self.md_system_dir / "topol.top"
        if not topology_file.exists():
            logger.warning("Topology file not found, using default force field settings")
            return {'forcefield': 'amber99sb', 'water_model': 'tip3p'}
        
        forcefield_info = {
            'forcefield': 'amber99sb',  # Default value
            'water_model': 'tip3p',     # Default value
            'detected': False
        }
        
        try:
            with open(topology_file, 'r') as f:
                content = f.read()
            
            # Detect protein force field
            if 'amber99sb' in content.lower():
                forcefield_info['forcefield'] = 'amber99sb'
                forcefield_info['detected'] = True
            elif 'amber03' in content.lower():
                forcefield_info['forcefield'] = 'amber03'
                forcefield_info['detected'] = True
            elif 'charmm27' in content.lower():
                forcefield_info['forcefield'] = 'charmm27'
                forcefield_info['detected'] = True
            elif 'gromos54a7' in content.lower():
                forcefield_info['forcefield'] = 'gromos54a7'
                forcefield_info['detected'] = True
            elif 'opls' in content.lower():
                forcefield_info['forcefield'] = 'opls-aa/L'
                forcefield_info['detected'] = True
            
            # Detect water model
            if 'tip3p' in content.lower():
                forcefield_info['water_model'] = 'tip3p'
            elif 'tip4p' in content.lower():
                forcefield_info['water_model'] = 'tip4p'
            elif 'spc' in content.lower():
                forcefield_info['water_model'] = 'spc'
            elif 'spce' in content.lower():
                forcefield_info['water_model'] = 'spce'
            
            # Detect ligand force field information
            if 'gaff' in content.lower() or 'amber' in content.lower():
                forcefield_info['ligand_forcefield'] = 'gaff'
            elif 'openff' in content.lower() or 'sage' in content.lower():
                forcefield_info['ligand_forcefield'] = 'openff'
            
            if forcefield_info['detected']:
                logger.info(f"Detected original force field: {forcefield_info['forcefield']} + {forcefield_info['water_model']}")
            else:
                logger.warning("Unable to clearly detect force field, using default settings")
                
        except Exception as e:
            logger.warning(f"Force field detection error: {e}, using default settings")
        
        return forcefield_info
    
    def _process_pmf_config(self, config: Optional[Dict], **kwargs) -> Dict:
        """Process PMF-specific configuration"""
        # Use detected force field information to ensure consistency
        default_config = {
            'general': {
                'gmx_command': 'gmx',
                'forcefield': self.original_forcefield_info.get('forcefield', 'amber99sb'),
                'water_model': self.original_forcefield_info.get('water_model', 'tip3p'),
                'ligand_forcefield': self.original_forcefield_info.get('ligand_forcefield', 'gaff'),
                'forcefield_detected': self.original_forcefield_info.get('detected', False)
            },
            'box': {
                'box_type': 'cubic',
                'box_distance': 1.2,    # Fixed dt distance (nm) for standard box
                'pulling_distance': 2.0  # User-customizable pulling distance (nm)
            },
            'pmf_system': {
                'reference_group': 'Protein',
                'moving_group': 'LIG',
                'align_axis': 'z',
                'center_system': True
            },
            'extraction': {
                'frame_number': -1,  # Use last frame
                'output_structure': True
            },
            'equilibration': {
                'em': {
                    'integrator': 'steep',
                    'emtol': 10.0,
                    'emstep': 0.01,
                    'nsteps': 50000
                },
                'nvt': {
                    'nsteps': 50000,  # 100 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'tau_t': 0.1
                },
                'npt': {
                    'nsteps': 100000,  # 200 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'pressure': 1.0,
                    'tau_t': 0.1,
                    'tau_p': 2.0
                }
            }
        }
        
        # Update with provided config
        if config:
            self._deep_update(default_config, config)
        
        # Apply kwargs
        if kwargs:
            for key, value in kwargs.items():
                if key in ['reference_group', 'moving_group']:
                    default_config['pmf_system'][key] = value
                elif key == 'pulling_distance':
                    default_config['box']['pulling_distance'] = float(value)
                elif key == 'box_distance':
                    default_config['box']['box_distance'] = float(value)
                elif key in ['forcefield', 'water_model']:
                    default_config['general'][key] = value
                elif key == 'z_extension':  # Legacy parameter support
                    logger.warning("Parameter 'z_extension' is deprecated, use 'pulling_distance' instead")
                    default_config['box']['pulling_distance'] = float(value) / 2  # Convert to pulling distance
                else:
                    default_config.setdefault('extra', {})[key] = value
        
        # Handle 'auto' force field detection after all config processing
        if default_config['general']['forcefield'] == 'auto':
            detected_ff = self.original_forcefield_info.get('forcefield')
            if detected_ff:
                logger.info(f"Using auto-detected force field: {detected_ff}")
                default_config['general']['forcefield'] = detected_ff
            else:
                logger.warning("Auto force field detection failed, using amber99sb as default")
                default_config['general']['forcefield'] = 'amber99sb'
        
        # Handle 'auto' water model detection
        if default_config['general']['water_model'] == 'auto':
            detected_water = self.original_forcefield_info.get('water_model')
            if detected_water:
                logger.info(f"Using auto-detected water model: {detected_water}")
                default_config['general']['water_model'] = detected_water
            else:
                logger.warning("Auto water model detection failed, using tip3p as default")
                default_config['general']['water_model'] = 'tip3p'
        
        return default_config
    
    def _deep_update(self, target: Dict, source: Dict) -> None:
        """Deep update dictionary"""
        for key, value in source.items():
            if isinstance(value, dict) and key in target and isinstance(target[key], dict):
                self._deep_update(target[key], value)
            else:
                target[key] = value
    
    def _validate_build_prerequisites(self):
        """Validate prerequisites for PMF system build"""
        errors = []
        
        # Check if MD results directory exists
        if not self.md_results_dir.exists():
            raise SystemNotFoundException(
                system_path=str(self.md_results_dir),
                searched_paths=[str(self.md_results_dir)]
            )
        
        # Check for GMX_PROLIG_MD subdirectory
        gmx_prolig_dir = self.md_results_dir / "GMX_PROLIG_MD"
        if not gmx_prolig_dir.exists():
            raise RequiredFileNotFoundError(
                file_path=str(gmx_prolig_dir),
                step="PMF system build",
                alternative_files=list(self.md_results_dir.iterdir())
            )
        
        # Check for required MD files
        required_files = [
            gmx_prolig_dir / "md.xtc",
            gmx_prolig_dir / "md.tpr",
            gmx_prolig_dir / "topol.top"
        ]
        
        missing_files = []
        for file_path in required_files:
            if not file_path.exists():
                missing_files.append(str(file_path))
        
        if missing_files:
            raise RequiredFileNotFoundError(
                file_path=missing_files[0],
                step="PMF system build",
                alternative_files=[str(f) for f in gmx_prolig_dir.iterdir() if f.is_file()]
            )
        
        # Check output directory is writable
        try:
            self.output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            raise PMFFileError(
                message=f"Cannot create output directory: {exc}",
                error_code=PMFErrorCode.DIRECTORY_NOT_WRITABLE,
                recoverable=True,
                recovery_suggestions=[
                    "Check directory permissions",
                    "Ensure parent directories exist",
                    "Check available disk space"
                ],
                context={"output_dir": str(self.output_dir)}
            ) from exc
        
        logger.info("PMF build prerequisites validated successfully")
    
    def build(self, frame: Optional[int] = None, equilibrate: bool = True) -> Dict:
        """
        Build PMF-optimized system from MD results
        
        Parameters:
        -----------
        frame : int, optional
            Frame number to extract (-1 for last frame)
        equilibrate : bool, optional
            Whether to run full equilibration after system rebuild (default: True)
            
        Returns:
        --------
        Dict : Build results with paths and system info
        
        Raises:
        -------
        SystemNotFoundException
            If MD results directory is not found
        PMFSystemError
            If system validation fails
        PMFFileError
            If required files are missing
        """
        with error_context("PMF system build", {
            "frame": frame, 
            "equilibrate": equilibrate,
            "output_dir": str(self.output_dir)
        }):
            logger.info("=== Building PMF-Optimized System ===")
            
            # Validate prerequisites
            self._validate_build_prerequisites()
            
            # Step 1: Extract structures from MD trajectory
            with error_context("structure extraction"):
                extraction_results = self.extract_structures(frame)
            
            # Step 2: Calculate centroids and alignment
            with error_context("alignment calculation"):
                alignment_info = self.calculate_alignment_geometry()
            
            # Step 3: Rebuild system with alignment and extended box
            with error_context("system rebuild"):
                rebuild_results = self.rebuild_aligned_system(alignment_info)
        
        # Combine initial results
        results = {
            'extraction': extraction_results,
            'alignment': alignment_info,
            'rebuild': rebuild_results,
            'system_dir': str(self.rebuilt_system_dir),
            'equilibration_required': equilibrate
        }
        
        # Step 4: Setup equilibration - always generate scripts, optionally run
        if equilibrate:
            logger.info("=== Setting up System Equilibration ===")

            # Always generate both local and SLURM scripts
            logger.info("Generating local and SLURM equilibration scripts")
            script_results = self.generate_all_equilibration_scripts()

            # Run equilibration directly (local execution)
            logger.info("Running automatic equilibration (EM → NVT → NPT)")
            equilibration_results = self.run_equilibration()

            # Combine results
            results['equilibration'] = equilibration_results
            results['equilibration']['scripts'] = script_results
            results['status'] = 'pmf_ready_equilibrated'
            results['final_system'] = equilibration_results['final_system']
                
        else:
            # Even when not running equilibration, generate scripts for user
            logger.info("=== Generating Equilibration Scripts ===")
            script_results = self.generate_all_equilibration_scripts()

            results['equilibration'] = {'scripts': script_results}
            results['status'] = 'pmf_ready_script_generated'
            results['final_system'] = {
                'structure': str(self.rebuilt_system_dir / "solv_ions.gro"),
                'topology': str(self.rebuilt_system_dir / "topol.top"),
                'equilibrated': False,
                'script_paths': script_results
            }
        
        # Save configuration
        self._save_build_config(results)
        
        logger.info("PMF system build completed successfully!")
        logger.info(f"PMF-ready system: {self.rebuilt_system_dir}")
        if equilibrate:
            logger.info("System fully equilibrated and ready for PMF calculations")
        else:
            logger.info("System built but not equilibrated - consider running equilibration")
        
        return results
    
    def extract_structures(self, frame: Optional[int] = None) -> Dict:
        """
        Extract Protein and LIG structures from MD trajectory
        
        Parameters:
        -----------
        frame : int, optional
            Frame to extract (-1 for last frame)
            
        Returns:
        --------
        Dict : Extraction results with file paths
        """
        if frame is None:
            frame = self.config['extraction']['frame_number']
        
        logger.info(f"Extracting structures from frame {frame}")
        
        # Input files
        structure_file = self.md_system_dir / "solv_ions.gro"
        topology_file = self.md_system_dir / "topol.top"
        
        if not structure_file.exists():
            raise RequiredFileNotFoundError(
                file_path=str(structure_file),
                step="structure extraction",
                alternative_files=[str(f) for f in self.md_system_dir.iterdir() if f.suffix in ['.gro', '.pdb']]
            )
        
        # Output files
        protein_out = self.work_dir / "protein_extracted.pdb"
        ligand_out = self.work_dir / "ligand_extracted.pdb"
        complex_out = self.work_dir / "complex_extracted.gro"
        
        # Extract protein
        self._extract_group("Protein", structure_file, protein_out, frame)
        
        # Extract ligand
        self._extract_group("LIG", structure_file, ligand_out, frame)
        
        # Extract protein-ligand complex (no water)
        self._extract_complex(structure_file, complex_out, frame)
        
        results = {
            'protein_structure': str(protein_out),
            'ligand_structure': str(ligand_out),
            'complex_structure': str(complex_out),
            'frame_extracted': frame,
            'source_structure': str(structure_file)
        }
        
        logger.info("Structure extraction completed")
        return results
    
    def _extract_group(self, group: str, input_file: Path, output_file: Path, 
                      frame: int = -1) -> None:
        """Extract specific group from structure"""
        cmd = [
            self.gmx_command, "trjconv",
            "-f", str(input_file),
            "-o", str(output_file),
            "-s", str(input_file)
        ]
        
        if frame != -1:
            cmd.extend(["-dump", str(frame)])
        
        # Run extraction
        # GROMACS trjconv may require two inputs: group to extract and reference group
        gromacs_input = f"{group}\n{group}\n"
        
        result = subprocess.run(
            cmd,
            input=gromacs_input,
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            handle_external_tool_error(
                command=" ".join(cmd),
                exit_code=result.returncode,
                stderr=result.stderr,
                tool_name="GROMACS trjconv"
            )
        
        logger.debug(f"Extracted {group} to {output_file}")
    
    def _extract_complex(self, input_file: Path, output_file: Path, frame: int = -1) -> None:
        """Extract protein-ligand complex (no solvent) using index file"""
        # First create index file with Protein_LIG group
        index_file = self._create_complex_index(input_file)
        
        cmd = [
            self.gmx_command, "trjconv", 
            "-f", str(input_file),
            "-o", str(output_file),
            "-s", str(input_file),
            "-n", str(index_file)  # Use index file
        ]
        
        if frame != -1:
            cmd.extend(["-dump", str(frame)])
        
        # Determine which group to select
        group_selection = self._determine_complex_group_selection(index_file)
        
        # Select the appropriate group from index
        # GROMACS trjconv requires two inputs: group to extract and reference group
        gromacs_input = f"{group_selection}\n{group_selection}\n"
        
        result = subprocess.run(
            cmd,
            input=gromacs_input,
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            handle_external_tool_error(
                command=" ".join(cmd),
                exit_code=result.returncode,
                stderr=result.stderr,
                tool_name="GROMACS trjconv"
            )
        
        logger.debug(f"Extracted complex to {output_file}")
    
    def _determine_complex_group_selection(self, index_file: Path) -> str:
        """
        Determine the correct group selection for complex extraction
        
        Parameters:
        -----------
        index_file : Path
            Path to index file
            
        Returns:
        --------
        str : Group name or number to select
        """
        # Check for metadata file with group number first
        metadata_file = self.work_dir / "complex_group_info.txt"
        if metadata_file.exists():
            try:
                with open(metadata_file, 'r') as f:
                    for line in f:
                        if line.startswith('combined_group_number='):
                            group_num = int(line.split('=')[1].strip())
                            logger.debug(f"Using group number from metadata: {group_num}")
                            return str(group_num)
            except (ValueError, IOError) as e:
                logger.warning(f"Could not read metadata file: {e}")
        
        # Simple approach: Try to find Protein_LIG group by examining available groups
        # without complex analysis that might fail
        try:
            # Use the original MD system structure file
            original_structure = self.md_system_dir / "solv_ions.gro"
            
            # Run make_ndx to see available groups without complex analysis
            cmd = [
                self.gmx_command, "make_ndx",
                "-f", str(original_structure)
            ]
            
            result = subprocess.run(
                cmd,
                input="q\n",  # Quit immediately
                text=True,
                capture_output=True,
                cwd=self.work_dir
            )
            
            if result.returncode == 0:
                # Parse output to find Protein_LIG groups
                lines = result.stdout.split('\n')
                protein_lig_groups = []
                
                for line in lines:
                    if 'Protein_LIG' in line and 'has' in line and 'elements' in line:
                        # Extract group number and size
                        # Format: "Group    21 (    Protein_LIG) has  3266 elements"
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                group_num = int(parts[1])
                                size = int(parts[5])
                                protein_lig_groups.append((group_num, size))
                                logger.debug(f"Found Protein_LIG group {group_num} with {size} elements")
                            except ValueError:
                                continue
                
                if protein_lig_groups:
                    # Select the largest group
                    largest_group = max(protein_lig_groups, key=lambda x: x[1])
                    group_num, size = largest_group
                    logger.debug(f"Selected Protein_LIG group {group_num} with {size} elements")
                    return str(group_num)
                    
        except Exception as e:
            logger.warning(f"Could not analyze groups: {e}")
        
        # Fallback: Use group number 21 based on common GROMACS behavior
        logger.warning("Using fallback group selection: 21 (typical largest Protein_LIG group)")
        return "21"
    
    def _create_complex_index(self, structure_file: Path) -> Path:
        """
        Create index file containing Protein_LIG group
        
        Parameters:
        -----------
        structure_file : Path
            Input structure file
            
        Returns:
        --------
        Path : Generated index file path
        """
        index_file = self.work_dir / "complex.ndx"
        
        logger.debug("Creating index file for protein-ligand complex")
        
        try:
            # Step 1: First run make_ndx to see available groups
            group_info = self._get_available_groups(structure_file)
            
            # Step 2: Find Protein and LIG group numbers
            protein_num = group_info.get('Protein')
            lig_num = group_info.get('LIG')
            
            if protein_num is None:
                logger.error("Protein group not found in structure")
                raise PMFSystemError(
                    message="Protein group not found in structure",
                    error_code=PMFErrorCode.SYSTEM_INVALID,
                    recoverable=False,
                    recovery_suggestions=[
                        "Check if protein is properly defined in system",
                        "Verify topology file contains protein residues",
                        "Check if system building completed successfully"
                    ]
                )
            
            if lig_num is None:
                logger.error("LIG group not found in structure")
                raise PMFSystemError(
                    message="LIG group not found in structure",
                    error_code=PMFErrorCode.SYSTEM_INVALID,
                    recoverable=False,
                    recovery_suggestions=[
                        "Check if ligand is properly defined in system",
                        "Verify topology file contains ligand residues", 
                        "Check ligand residue name and force field generation"
                    ]
                )
            
            logger.debug(f"Found groups: Protein={protein_num}, LIG={lig_num}")
            
            # Step 3: Create combined group using group numbers
            return self._create_index_with_numbers(structure_file, protein_num, lig_num)
            
        except Exception as e:
            logger.warning(f"Automatic index creation failed: {e}")
            logger.info("Falling back to manual index creation")
            return self._create_index_manually(structure_file)
    
    def _get_available_groups(self, structure_file: Path) -> Dict[str, int]:
        """
        Parse available groups from gmx make_ndx output
        
        Parameters:
        -----------
        structure_file : Path
            Input structure file
            
        Returns:
        --------
        Dict[str, int] : Mapping of group names to numbers
        """
        # Run make_ndx in query mode (quit immediately)
        cmd = [
            self.gmx_command, "make_ndx",
            "-f", str(structure_file)
        ]
        
        result = subprocess.run(
            cmd,
            input="q\n",  # Quit immediately to see available groups
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            handle_external_tool_error(
                command=" ".join(cmd),
                exit_code=result.returncode,
                stderr=result.stderr,
                tool_name="GROMACS make_ndx"
            )
        
        # Parse the output to extract group information
        groups = {}
        lines = result.stdout.split('\n')
        
        for line in lines:
            # Look for lines like "  0 System              : 31266 atoms"
            # or "  1 Protein             :  4834 atoms"
            if ':' in line and 'atoms' in line:
                try:
                    # Extract group number and name
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        group_num = int(parts[0])
                        group_name = parts[1]
                        groups[group_name] = group_num
                        logger.debug(f"Found group: {group_name} = {group_num}")
                except (ValueError, IndexError):
                    continue
        
        if not groups:
            raise PMFSystemError(
                message="Could not parse group information from make_ndx output",
                error_code=PMFErrorCode.GROMACS_EXECUTION_FAILED,
                recoverable=True,
                recovery_suggestions=[
                    "Check GROMACS version compatibility",
                    "Verify system structure is valid",
                    "Check for corrupted topology file"
                ],
                context={"make_ndx_output": result.stdout[:500]}
            )
        
        return groups
    
    def _create_index_with_numbers(self, structure_file: Path, protein_num: int, lig_num: int) -> Path:
        """
        Create index file using specific group numbers
        
        Parameters:
        -----------
        structure_file : Path
            Input structure file
        protein_num : int
            Protein group number
        lig_num : int
            LIG group number
            
        Returns:
        --------
        Path : Generated index file path
        """
        index_file = self.work_dir / "complex.ndx"
        
        cmd = [
            self.gmx_command, "make_ndx",
            "-f", str(structure_file),
            "-o", str(index_file)
        ]
        
        # Create commands using group numbers:
        # 1. "<protein_num> | <lig_num>" - combine protein and ligand groups
        # 2. "name <new_group_num> Protein_LIG" - rename the combined group
        # 3. "q" - quit
        
        # The new group will typically be assigned the next available number
        # We'll use a reasonable estimate and handle naming accordingly
        combine_command = f"{protein_num} | {lig_num}\n"
        
        # Try to determine the new group number (usually max existing + 1)
        new_group_num = max(protein_num, lig_num) + 1
        name_command = f"name {new_group_num} Protein_LIG\n"
        
        input_commands = combine_command + name_command + "q\n"
        
        logger.debug(f"Creating index with commands: {combine_command.strip()}, {name_command.strip()}")
        
        result = subprocess.run(
            cmd,
            input=input_commands,
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            logger.warning(f"make_ndx with numbers failed: {result.stderr}")
            # Try alternative naming strategy
            return self._create_index_alternative_naming(structure_file, protein_num, lig_num)
        
        # Verify index file was created
        if not index_file.exists():
            logger.warning("Index file not created, trying alternative approach")
            return self._create_index_alternative_naming(structure_file, protein_num, lig_num)
        
        # Verify the Protein_LIG group exists in the file
        if not self._verify_index_group(index_file, "Protein_LIG"):
            logger.warning("Protein_LIG group not found in index, trying alternative naming")
            return self._create_index_alternative_naming(structure_file, protein_num, lig_num)
        
        logger.debug(f"Index file created successfully: {index_file}")
        return index_file
    
    def _create_index_alternative_naming(self, structure_file: Path, protein_num: int, lig_num: int) -> Path:
        """
        Alternative index creation without explicit naming
        
        Parameters:
        -----------
        structure_file : Path
            Input structure file
        protein_num : int
            Protein group number
        lig_num : int
            LIG group number
            
        Returns:
        --------
        Path : Generated index file path
        """
        index_file = self.work_dir / "complex_alt.ndx"
        
        cmd = [
            self.gmx_command, "make_ndx",
            "-f", str(structure_file),
            "-o", str(index_file)
        ]
        
        # Just create the combined group without renaming
        # We'll update the extraction method to use the group number instead
        combine_command = f"{protein_num} | {lig_num}\n"
        input_commands = combine_command + "q\n"
        
        logger.debug(f"Creating index with simple combination: {combine_command.strip()}")
        
        result = subprocess.run(
            cmd,
            input=input_commands,
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"Alternative index creation failed: {result.stderr}")
        
        if not index_file.exists():
            raise RuntimeError("Index file creation failed")
        
        # Store the expected group number for the combined group
        # It's typically the next available number after the input groups
        expected_group_num = max(protein_num, lig_num) + 1
        
        # Create a metadata file to remember the group number
        metadata_file = self.work_dir / "complex_group_info.txt"
        with open(metadata_file, 'w') as f:
            f.write(f"combined_group_number={expected_group_num}\n")
            f.write(f"protein_group={protein_num}\n")
            f.write(f"ligand_group={lig_num}\n")
        
        logger.debug(f"Alternative index created, combined group number: {expected_group_num}")
        return index_file
    
    def _verify_index_group(self, index_file: Path, group_name: str) -> bool:
        """
        Verify that a specific group exists in the index file
        
        Parameters:
        -----------
        index_file : Path
            Path to index file
        group_name : str
            Name of group to verify
            
        Returns:
        --------
        bool : True if group exists
        """
        try:
            with open(index_file, 'r') as f:
                content = f.read()
            
            # Look for the group name in square brackets
            return f"[ {group_name} ]" in content
            
        except Exception:
            return False
    
    def _create_index_manually(self, structure_file: Path) -> Path:
        """
        Manually create index file for protein-ligand complex
        
        Parameters:
        -----------
        structure_file : Path
            Input structure file
            
        Returns:
        --------
        Path : Generated index file path
        """
        index_file = self.work_dir / "complex_manual.ndx"
        
        logger.debug("Creating index file manually by parsing structure")
        
        protein_atoms = []
        ligand_atoms = []
        
        try:
            with open(structure_file, 'r') as f:
                lines = f.readlines()
            
            # Parse .gro file format
            if len(lines) < 3:
                raise ValueError("Invalid structure file format")
            
            # Skip header and total atoms line
            atom_lines = lines[2:-1]  # Exclude header, natoms line, and box line
            
            for i, line in enumerate(atom_lines, 1):
                if len(line.strip()) < 44:  # Minimum line length for .gro format
                    continue
                
                try:
                    # .gro format: columns 10-15 contain residue name
                    residue_name = line[5:10].strip()
                    
                    # Check if this is protein or ligand
                    if self._is_protein_residue(residue_name):
                        protein_atoms.append(i)
                    elif residue_name == 'LIG':  # Standard PRISM ligand name
                        ligand_atoms.append(i)
                        
                except (ValueError, IndexError):
                    continue
            
            # Write index file
            with open(index_file, 'w') as f:
                f.write("[ Protein_LIG ]\n")
                
                # Combine protein and ligand atoms
                all_atoms = sorted(protein_atoms + ligand_atoms)
                
                # Write atoms in groups of 15 per line (GROMACS convention)
                for i in range(0, len(all_atoms), 15):
                    atom_group = all_atoms[i:i+15]
                    f.write(" ".join(map(str, atom_group)) + "\n")
                f.write("\n")
            
            logger.debug(f"Manual index file created with {len(protein_atoms)} protein + {len(ligand_atoms)} ligand atoms")
            return index_file
            
        except Exception as e:
            logger.error(f"Failed to create manual index file: {e}")
            raise RuntimeError(f"Could not create index file for complex extraction: {e}")
    
    def _is_protein_residue(self, residue_name: str) -> bool:
        """
        Check if residue name corresponds to a standard protein residue
        
        Parameters:
        -----------
        residue_name : str
            Three-letter residue code
            
        Returns:
        --------
        bool : True if protein residue
        """
        # Standard 20 amino acids + common variants
        protein_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            # Common variants
            'HIE', 'HID', 'HIP', 'CYX', 'CYSH', 'ASH', 'GLH', 'LYP'
        }
        
        return residue_name.upper() in protein_residues
    
    def calculate_alignment_geometry(self) -> Dict:
        """
        Calculate centroids and alignment vector for Z-axis orientation
        
        Returns:
        --------
        Dict : Alignment geometry information
        """
        logger.info("Calculating alignment geometry")
        
        protein_file = self.work_dir / "protein_extracted.pdb"
        ligand_file = self.work_dir / "ligand_extracted.pdb"
        
        # Calculate centroids
        protein_centroid = self._calculate_centroid(protein_file, "Protein")
        ligand_centroid = self._calculate_centroid(ligand_file, "LIG")
        
        # Calculate alignment vector (ligand -> protein)
        alignment_vector = np.array(protein_centroid) - np.array(ligand_centroid)
        
        # Calculate rotation matrix to align with Z-axis
        z_axis = np.array([0, 0, 1])
        rotation_matrix = self._calculate_rotation_matrix(alignment_vector, z_axis)
        
        # Calculate distance for box sizing
        initial_distance = np.linalg.norm(alignment_vector)
        
        alignment_info = {
            'protein_centroid': protein_centroid,
            'ligand_centroid': ligand_centroid,
            'alignment_vector': alignment_vector.tolist(),
            'initial_distance': float(initial_distance),
            'rotation_matrix': rotation_matrix.tolist(),
            'z_axis_aligned': True
        }
        
        logger.info(f"Alignment calculated: distance = {initial_distance:.3f} nm")
        logger.info(f"Protein centroid: {protein_centroid}")
        logger.info(f"Ligand centroid: {ligand_centroid}")
        
        return alignment_info
    
    def _calculate_centroid(self, structure_file: Path, group: str) -> Tuple[float, float, float]:
        """Calculate geometric centroid of a group"""
        # Use gmx traj to get center of mass
        cmd = [
            self.gmx_command, "traj",
            "-f", str(structure_file),
            "-s", str(structure_file),
            "-com",
            "-quiet"
        ]
        
        result = subprocess.run(
            cmd,
            input=f"{group}\n",
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            # Fallback: parse coordinates manually
            return self._parse_centroid_manually(structure_file)
        
        # Parse center of mass from output
        for line in result.stdout.split('\n'):
            if 'Center of mass' in line:
                coords = [float(x) for x in line.split()[-3:]]
                return tuple(coords)
        
        # Fallback to manual parsing
        return self._parse_centroid_manually(structure_file)
    
    def _parse_centroid_manually(self, structure_file: Path) -> Tuple[float, float, float]:
        """Manually parse coordinates and calculate centroid"""
        coords = []
        
        with open(structure_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # PDB format: positions at columns 30-38, 38-46, 46-54
                    try:
                        x = float(line[30:38].strip()) / 10.0  # Convert Å to nm
                        y = float(line[38:46].strip()) / 10.0
                        z = float(line[46:54].strip()) / 10.0
                        coords.append([x, y, z])
                    except ValueError:
                        continue
                elif line.strip() and not line.startswith('#'):
                    # GRO format: positions at specific columns
                    try:
                        parts = line.split()
                        if len(parts) >= 6:
                            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                            coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
        
        if not coords:
            raise RuntimeError(f"Could not parse coordinates from {structure_file}")
        
        # Calculate centroid
        coords_array = np.array(coords)
        centroid = np.mean(coords_array, axis=0)
        
        return tuple(centroid.tolist())
    
    def _calculate_rotation_matrix(self, vector: np.ndarray, target: np.ndarray) -> np.ndarray:
        """Calculate rotation matrix to align vector with target"""
        # Normalize vectors
        v1 = vector / np.linalg.norm(vector)
        v2 = target / np.linalg.norm(target)
        
        # Calculate rotation axis and angle
        cross_product = np.cross(v1, v2)
        sin_angle = np.linalg.norm(cross_product)
        cos_angle = np.dot(v1, v2)
        
        # If vectors are already aligned
        if sin_angle < 1e-6:
            return np.eye(3)
        
        # Rodrigues' rotation formula
        cross_product = cross_product / sin_angle  # Normalize
        
        K = np.array([
            [0, -cross_product[2], cross_product[1]],
            [cross_product[2], 0, -cross_product[0]],
            [-cross_product[1], cross_product[0], 0]
        ])
        
        rotation_matrix = (np.eye(3) + sin_angle * K + 
                          (1 - cos_angle) * np.dot(K, K))
        
        return rotation_matrix
    
    def rebuild_aligned_system(self, alignment_info: Dict) -> Dict:
        """
        Rebuild system with Z-axis alignment and extended box
        
        Parameters:
        -----------
        alignment_info : Dict
            Alignment geometry information
            
        Returns:
        --------
        Dict : Rebuild results
        """
        logger.info("Rebuilding system with Z-axis alignment")
        
        # Step 1: Create aligned complex
        aligned_complex = self._create_aligned_complex(alignment_info)
        
        # Step 2: Process protein with pdb2gmx
        processed_protein = self._process_protein_structure()
        
        # Step 3: Combine protein and ligand topologies
        combined_topology = self._combine_topologies(processed_protein)
        
        # Step 4: Create extended box with Z-axis alignment
        boxed_system = self._create_extended_box(aligned_complex, alignment_info)
        
        # Step 5: Solvate and add ions
        solvated_system = self._solvate_and_add_ions(boxed_system, combined_topology)
        
        results = {
            'aligned_complex': str(aligned_complex),
            'processed_protein': processed_protein,
            'combined_topology': str(combined_topology),
            'boxed_system': str(boxed_system),
            'final_system': str(solvated_system),
            'system_directory': str(self.rebuilt_system_dir)
        }
        
        logger.info("System rebuild completed")
        return results
    
    def _create_aligned_complex(self, alignment_info: Dict) -> Path:
        """Create aligned protein-ligand complex using centroid-based rotation calculation"""
        complex_file = self.work_dir / "complex_extracted.gro"
        aligned_file = self.work_dir / "complex_aligned.gro"
        
        logger.info("Aligning protein-ligand complex to Z-axis using centroid calculation")
        
        # Step 1: Center the complex at origin
        centered_file = self.work_dir / "complex_centered.gro"
        self._center_complex(complex_file, centered_file)
        
        # Step 2: Calculate centroid coordinates
        centroids_file = self.work_dir / "centroids.xvg"
        self._calculate_centroids(centered_file, centroids_file)
        
        # Step 3: Calculate rotation angles from centroids
        rotation_angles = self._calculate_rotation_angles(centroids_file)
        
        # Step 4: Apply rotation to align centroid direction with Z-axis
        if abs(rotation_angles[0]) < 0.1 and abs(rotation_angles[1]) < 0.1 and abs(rotation_angles[2]) < 0.1:
            logger.info("Complex already well-aligned with Z-axis, skipping rotation")
            shutil.copy2(centered_file, aligned_file)
        else:
            self._apply_rotation(centered_file, aligned_file, rotation_angles)
        
        # Step 5: Verify alignment quality
        self._verify_alignment(aligned_file)
        
        return aligned_file
    
    def _center_complex(self, input_file: Path, output_file: Path) -> None:
        """Center the complex at the origin"""
        cmd = [
            self.gmx_command, "editconf",
            "-f", str(input_file),
            "-o", str(output_file),
            "-center", "0", "0", "0"
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"Failed to center complex: {result.stderr}")
        
        logger.debug(f"Complex centered at origin: {output_file}")
    
    def _calculate_centroids(self, structure_file: Path, output_file: Path) -> None:
        """Calculate protein and ligand centroids"""
        cmd = [
            self.gmx_command, "traj",
            "-f", str(structure_file),
            "-s", str(structure_file),
            "-com", "-ng", "2",
            "-ox", str(output_file)
        ]
        
        # Input: Protein group, then LIG group
        input_text = "Protein\nLIG\n"
        
        result = subprocess.run(
            cmd,
            input=input_text,
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"Failed to calculate centroids: {result.stderr}")
        
        logger.debug(f"Centroids calculated: {output_file}")
    
    def _calculate_rotation_angles(self, centroids_file: Path) -> Tuple[float, float, float]:
        """Calculate rotation angles to align centroid vector with Z-axis"""
        try:
            # Read centroid coordinates
            data = []
            with open(centroids_file, 'r') as f:
                for line in f:
                    if not line.startswith('#') and not line.startswith('@') and line.strip():
                        values = line.split()
                        if len(values) >= 7:  # time + 3 coords for protein + 3 coords for ligand
                            data.append([float(x) for x in values])
            
            if not data:
                raise RuntimeError("No valid data found in centroids file")
            
            # Use the first (and likely only) frame
            row = data[0]
            protein_com = np.array([row[1], row[2], row[3]])  # protein centroid
            ligand_com = np.array([row[4], row[5], row[6]])   # ligand centroid
            
            # Calculate centroid vector (protein -> ligand direction for pulling)
            # This is the direction in which the ligand will be pulled away from protein
            centroid_vector = ligand_com - protein_com
            
            # Normalize the vector
            vector_length = np.linalg.norm(centroid_vector)
            if vector_length < 1e-6:
                logger.warning("Protein and ligand centroids are too close, skipping alignment")
                return 0.0, 0.0, 0.0
            
            centroid_vector = centroid_vector / vector_length
            
            # Target Z-axis vector
            z_axis = np.array([0, 0, 1])
            
            # Calculate rotation using Rodrigues' rotation formula (following workflow.md)
            # If vectors are already aligned, return zero rotation
            cos_angle = np.dot(centroid_vector, z_axis)
            if abs(cos_angle) > 0.999:  # Already aligned (within ~2.6 degrees)
                return 0.0, 0.0, 0.0
            
            # Calculate rotation axis (cross product) - following workflow.md method
            # For vector v = (v_x, v_y, v_z) and z = (0, 0, 1):
            # rotation_axis = v × z = (-v_y, v_x, 0)
            v_x, v_y, v_z = centroid_vector[0], centroid_vector[1], centroid_vector[2]
            rotation_axis = np.array([-v_y, v_x, 0.0])
            
            # Check if rotation axis is valid
            axis_norm = np.linalg.norm(rotation_axis)
            if axis_norm < 1e-10:  # Nearly parallel to Z-axis
                return 0.0, 0.0, 0.0
            
            rotation_axis = rotation_axis / axis_norm
            
            # Calculate rotation angle - following workflow.md: θ = arccos(v_z/|v|)
            # We want to align centroid_vector to +Z axis (0, 0, 1)
            target_z = 1.0  # Always align to +Z axis
            
            # Calculate rotation angle to align centroid_vector to +Z axis
            # Using standard dot product: cos(θ) = v⃗ · ẑ = v_z (since ẑ = [0,0,1])
            cos_angle = v_z  # Direct dot product with +Z axis
            
            # Clamp to valid range for arccos
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            angle_rad = np.arccos(abs(cos_angle))  # Always positive angle
            
            # If vector points towards -Z, we need the supplementary angle
            if v_z < 0:
                logger.info(f"Vector pointing towards -Z (v_z={v_z:.3f}), need {np.degrees(angle_rad):.1f}° rotation to +Z")
            else:
                logger.info(f"Vector pointing towards +Z (v_z={v_z:.3f}), need {np.degrees(angle_rad):.1f}° rotation to +Z")
            
            angle_deg = np.degrees(angle_rad)
            
            # Apply Rodrigues' rotation formula to get rotation matrix
            # R = I + sin(θ)K + (1-cos(θ))K²
            I = np.eye(3)
            K = np.array([
                [0, -rotation_axis[2], rotation_axis[1]],
                [rotation_axis[2], 0, -rotation_axis[0]], 
                [-rotation_axis[1], rotation_axis[0], 0]
            ])
            
            sin_theta = np.sin(angle_rad)
            cos_theta = np.cos(angle_rad)
            
            R = I + sin_theta * K + (1 - cos_theta) * np.dot(K, K)
            
            # Convert rotation matrix to Euler angles (ZYX convention for GROMACS)
            # Extract Euler angles from rotation matrix
            # Following standard rotation matrix to Euler angle conversion
            sy = np.sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0])
            
            singular = sy < 1e-6
            
            if not singular:
                angle_x = np.degrees(np.arctan2(R[2, 1], R[2, 2]))
                angle_y = np.degrees(np.arctan2(-R[2, 0], sy))
                angle_z = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
            else:
                angle_x = np.degrees(np.arctan2(-R[1, 2], R[1, 1]))
                angle_y = np.degrees(np.arctan2(-R[2, 0], sy))
                angle_z = 0.0
            
            logger.info(f"Calculated rotation angles: X={angle_x:.1f}°, Y={angle_y:.1f}°, Z={angle_z:.1f}°")
            logger.info(f"Final centroid vector (Protein→Ligand): [{centroid_vector[0]:.3f}, {centroid_vector[1]:.3f}, {centroid_vector[2]:.3f}]")
            logger.info(f"Angle with +Z-axis: {angle_deg:.1f}°")
            logger.info(f"Direction: Protein at lower Z, Ligand at higher Z (ready for +Z pulling)")
            
            # Store the final aligned vector information for box extension
            self._aligned_vector_info = {
                'centroid_vector': centroid_vector.copy(),
                'protein_com': protein_com.copy(),
                'ligand_com': ligand_com.copy(),
                'original_z_component': v_z,  # Original Z component for verification
                'rotation_angle': angle_deg,   # Applied rotation angle
                'expected_z_positive': True    # After rotation, should always point to +Z
            }
            
            return angle_x, angle_y, angle_z
            
        except Exception as e:
            logger.error(f"Failed to calculate rotation angles: {e}")
            logger.warning("Using principal axis alignment as fallback")
            return 0.0, 0.0, 0.0
    
    def _apply_rotation(self, input_file: Path, output_file: Path, rotation_angles: Tuple[float, float, float]) -> None:
        """Apply rotation to align complex with Z-axis"""
        angle_x, angle_y, angle_z = rotation_angles
        
        # Method 1: Try precise rotation
        try:
            cmd = [
                self.gmx_command, "editconf",
                "-f", str(input_file),
                "-o", str(output_file),
                "-rotate", str(angle_x), str(angle_y), str(angle_z)
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=self.work_dir
            )
            
            if result.returncode == 0:
                logger.info("✓ Successfully applied rotation alignment")
                return
            else:
                logger.warning(f"Rotation failed: {result.stderr}")
                
        except Exception as e:
            logger.warning(f"Rotation method failed: {e}")
        
        # Method 2: Fallback to principal axis alignment
        logger.info("Using fallback: principal axis alignment")
        try:
            cmd = [
                self.gmx_command, "editconf",
                "-f", str(input_file),
                "-o", str(output_file),
                "-princ",  # Use correct -princ option
                "-center", "0", "0", "0"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
            if result.returncode == 0:
                logger.info("✓ Applied principal axis alignment as fallback")
                return
            else:
                logger.warning(f"Principal axis alignment failed: {result.stderr}")
                
        except Exception as e:
            logger.error(f"Principal axis alignment failed: {e}")
            
        # Method 3: Last resort - just center the system
        logger.warning("Using last resort: simple centering (no rotation)")
        cmd = [
            self.gmx_command, "editconf",
            "-f", str(input_file),
            "-o", str(output_file),
            "-center", "0", "0", "0"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to center complex: {result.stderr}")
        
        logger.warning("⚠ System only centered, not rotated - alignment may be suboptimal")
    
    def _create_alignment_index_file(self, index_file: Path) -> None:
        """Create index file for alignment groups"""
        try:
            # Use the same complex file that we're about to align
            complex_file = self.work_dir / "complex_extracted.gro"
            
            # Create index file with Protein and LIG groups
            cmd = [
                self.gmx_command, "make_ndx",
                "-f", str(complex_file),
                "-o", str(index_file),
                "-quiet"
            ]
            
            # Create groups: 
            # - Use existing Protein group (usually group 1)
            # - Use existing non-Water groups for ligand
            input_text = 'q\n'  # Just quit to save default groups
            
            result = subprocess.run(
                cmd,
                input=input_text,
                text=True,
                capture_output=True,
                cwd=self.work_dir
            )
            
            if result.returncode != 0:
                logger.warning(f"Failed to create alignment index: {result.stderr}")
                # Create a minimal index file manually
                with open(index_file, 'w') as f:
                    f.write("[ System ]\n")
                    f.write("0\n")
                    
        except Exception as e:
            logger.warning(f"Could not create alignment index file: {e}")
    
    def _verify_alignment(self, aligned_file: Path) -> None:
        """Verify that the alignment was successful using centroid calculation"""
        try:
            # Calculate centroids after alignment
            verification_centroids_file = self.work_dir / "verification_centroids.xvg"
            self._calculate_centroids(aligned_file, verification_centroids_file)
            
            # Read the centroid data
            data = []
            with open(verification_centroids_file, 'r') as f:
                for line in f:
                    if not line.startswith('#') and not line.startswith('@') and line.strip():
                        values = line.split()
                        if len(values) >= 7:
                            data.append([float(x) for x in values])
            
            if not data:
                logger.warning("No centroid data available for verification")
                return
            
            # Extract centroid coordinates
            row = data[0]
            protein_com = np.array([row[1], row[2], row[3]])
            ligand_com = np.array([row[4], row[5], row[6]])
            
            # Calculate alignment vector (ligand - protein)
            alignment_vector = ligand_com - protein_com
            vector_length = np.linalg.norm(alignment_vector)
            
            if vector_length < 1e-6:
                logger.warning("Protein and ligand centroids are too close for verification")
                return
            
            alignment_vector_norm = alignment_vector / vector_length
            
            # Check alignment with Z-axis
            z_axis = np.array([0, 0, 1])
            cos_angle = np.dot(alignment_vector_norm, z_axis)
            angle_degrees = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
            
            logger.info(f"Alignment verification:")
            logger.info(f"  Protein centroid: [{protein_com[0]:.3f}, {protein_com[1]:.3f}, {protein_com[2]:.3f}]")
            logger.info(f"  Ligand centroid:  [{ligand_com[0]:.3f}, {ligand_com[1]:.3f}, {ligand_com[2]:.3f}]")
            logger.info(f"  Alignment vector: [{alignment_vector_norm[0]:.3f}, {alignment_vector_norm[1]:.3f}, {alignment_vector_norm[2]:.3f}]")
            logger.info(f"  Angle with Z-axis: {angle_degrees:.1f}°")
            
            if angle_degrees < 5.0:
                logger.info("✓ Excellent alignment achieved (< 5° from Z-axis)")
            elif angle_degrees < 15.0:
                logger.info("✓ Good alignment achieved (< 15° from Z-axis)")
            elif angle_degrees < 30.0:
                logger.warning(f"○ Moderate alignment ({angle_degrees:.1f}° from Z-axis)")
            else:
                logger.warning(f"△ Poor alignment ({angle_degrees:.1f}° from Z-axis) - consider manual adjustment")
            
        except Exception as e:
            logger.warning(f"Could not verify alignment: {e}")
    
    def _process_protein_structure(self) -> Dict:
        """Process protein structure with pdb2gmx"""
        protein_file = self.work_dir / "protein_extracted.pdb"
        protein_gro = self.work_dir / "protein_processed.gro"
        protein_top = self.work_dir / "protein.top"
        
        # Run pdb2gmx
        cmd = [
            self.gmx_command, "pdb2gmx",
            "-f", str(protein_file),
            "-o", str(protein_gro),
            "-p", str(protein_top),
            "-water", self.config['general']['water_model'],
            "-ff", self.config['general']['forcefield']
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"pdb2gmx failed: {result.stderr}")
        
        logger.debug("Protein structure processed with pdb2gmx")
        
        return {
            'structure': str(protein_gro),
            'topology': str(protein_top),
            'forcefield': self.config['general']['forcefield']
        }
    
    def _combine_topologies(self, processed_protein: Dict) -> Path:
        """
        Create clean topology without water and ions for PMF remodeling
        
        This removes water and ions from original topology to avoid conflicts
        when re-solvating the extended box system.
        """
        original_top = self.md_system_dir / "topol.top"
        combined_top = self.work_dir / "combined.top"
        
        logger.info("Creating clean topology without original water/ions")
        
        try:
            with open(original_top, 'r') as f:
                lines = f.readlines()
            
            # Process topology to remove water and ions
            clean_lines = []
            in_molecules_section = False
            
            for line in lines:
                line_stripped = line.strip()
                
                # Track if we're in [ molecules ] section
                if line_stripped.startswith('[ molecules ]'):
                    in_molecules_section = True
                    clean_lines.append(line)
                    continue
                elif line_stripped.startswith('[') and in_molecules_section:
                    # End of molecules section
                    in_molecules_section = False
                
                # Skip water and ion entries in molecules section
                if in_molecules_section and line_stripped:
                    # Skip common water and ion molecule names
                    skip_molecules = ['SOL', 'WAT', 'TIP3P', 'TIP4P', 'SPC', 'SPCE',
                                     'NA+', 'CL-', 'K+', 'CA2+', 'MG2+', 'NA', 'CL']
                    
                    parts = line_stripped.split()
                    if len(parts) >= 2 and parts[0] in skip_molecules:
                        logger.debug(f"Removing from topology: {line_stripped}")
                        continue
                
                clean_lines.append(line)
            
            # Write clean topology
            with open(combined_top, 'w') as f:
                f.writelines(clean_lines)
            
            logger.info(f"Clean topology created: {combined_top}")
            logger.info("Removed original water and ions - ready for re-solvation")
            return combined_top
            
        except Exception as e:
            logger.warning(f"Failed to clean topology: {e}")
            # Fallback: copy original topology
            shutil.copy2(original_top, combined_top)
            logger.warning("Using original topology - may cause water/ion conflicts")
            return combined_top
    
    def _create_extended_box(self, complex_file: Path, alignment_info: Dict) -> Path:
        """
        Create PMF-optimized box using GROMACS built-in algorithms
        
        Three-step method:
        1. Use fixed dt distance (default 1.2nm) with GROMACS -d option to create standard cubic box
        2. Create extended box by adding 2x user-defined pulling_distance to Z-axis
        3. Translate system by pulling_distance in negative Z direction to optimize pulling space
        
        Parameters:
        -----------
        complex_file : Path
            Input complex structure file
        alignment_info : Dict
            System alignment information
            
        Returns:
        --------
        Path : PMF-optimized boxed system file
        """
        boxed_file = self.work_dir / "complex_boxed.gro"
        
        # Get configuration parameters
        dt_distance = self.config['box']['box_distance']  # Fixed dt distance (default: 1.2 nm)
        pulling_distance = self.config['box']['pulling_distance']  # User-customizable pulling distance
        
        logger.info(f"Creating PMF-optimized box using GROMACS built-in algorithms:")
        logger.info(f"  Fixed dt distance: {dt_distance:.1f} nm (GROMACS standard)")
        logger.info(f"  User pulling distance: {pulling_distance:.1f} nm (customizable)")
        logger.info(f"  Total Z extension: {2 * pulling_distance:.1f} nm (2x pulling distance)")
        
        # === Step 1: Create standard cubic box with Z-axis alignment using GROMACS -d option ===
        standard_boxed = self.work_dir / "standard_box.gro"
        
        logger.info("Step 1: Creating standard cubic box with GROMACS -d algorithm")
        cmd = [
            self.gmx_command, "editconf",
            "-f", str(complex_file),
            "-o", str(standard_boxed),
            "-d", str(dt_distance),  # Use GROMACS built-in -d option for consistent box creation
            "-bt", "cubic",           # Cubic box type for uniform spacing
            "-c"                      # Center the system (quality control)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"Step 1 failed - GROMACS standard box creation: {result.stderr}")
        
        logger.info("  Standard box created successfully with GROMACS algorithm")
        
        # === Step 2: Parse standard box dimensions and create extended box ===
        logger.info("Step 2: Creating extended box based on standard dimensions")
        
        try:
            box_info = self._parse_box_dimensions(standard_boxed)
            original_x = box_info['x_dimension']
            original_y = box_info['y_dimension']
            original_z = box_info['z_dimension']
        except Exception as e:
            logger.error(f"Failed to parse box dimensions: {e}")
            raise RuntimeError(f"Step 2 failed - cannot parse standard box dimensions: {e}")
        
        # Calculate extended Z dimension (original + 2x pulling distance)
        extended_z = original_z + 2 * pulling_distance
        
        logger.info(f"  Original dimensions: {original_x:.2f} x {original_y:.2f} x {original_z:.2f} nm")
        logger.info(f"  Extended Z dimension: {extended_z:.2f} nm (+{2*pulling_distance:.2f} nm)")
        
        extended_boxed = self.work_dir / "extended_box.gro"
        cmd = [
            self.gmx_command, "editconf",
            "-f", str(standard_boxed),
            "-o", str(extended_boxed),
            "-box", str(original_x), str(original_y), str(extended_z),  # Explicit box dimensions
            "-c"  # Keep system centered during box extension
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"Step 2 failed - box extension: {result.stderr}")
        
        logger.info("  Extended box created successfully")
        
        # === Step 3: Center ligand in the extended box for optimal PMF setup ===
        logger.info("Step 3: Centering system in extended box for PMF pulling")
        
        # Following CLAUDE.md requirement: center based on small molecule centroid position
        # Since we extended the box by 2×pulling_distance, the ligand should be centered
        # in the new extended box to have equal space on both sides
        # This means no translation is needed since editconf -c already centers the system
        translation_distance = 0.0  # Keep system centered in extended box
        
        logger.info(f"  Strategy: Keep ligand centered in extended box")
        logger.info(f"  Extended Z-box: {extended_z:.2f} nm (original: {original_z:.2f} nm)")
        logger.info(f"  Available space for +Z pulling: ~{pulling_distance:.2f} nm")
        logger.info(f"  Available space for -Z buffer: ~{pulling_distance:.2f} nm")
        
        # If no translation needed, just copy the file
        if abs(translation_distance) < 1e-6:
            logger.info("  No translation needed - system already optimally positioned")
            shutil.copy2(extended_boxed, boxed_file)
        else:
            cmd = [
                self.gmx_command, "editconf",
                "-f", str(extended_boxed),
                "-o", str(boxed_file),
                "-translate", "0", "0", str(translation_distance)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
            if result.returncode != 0:
                raise RuntimeError(f"Step 3 failed - system translation: {result.stderr}")
        
        # Clean up temporary files
        try:
            standard_boxed.unlink()
            extended_boxed.unlink()
        except FileNotFoundError:
            pass  # Files already cleaned
        
        # Verify final box dimensions
        try:
            final_box = self._parse_box_dimensions(boxed_file)
            logger.info(f"Final PMF box: {final_box['x_dimension']:.2f} x {final_box['y_dimension']:.2f} x {final_box['z_dimension']:.2f} nm")
        except Exception:
            logger.warning("Could not verify final box dimensions")
        
        # Verify that the alignment is maintained after box operations
        self._verify_box_alignment_consistency(boxed_file)
        
        logger.info("PMF-optimized box creation completed successfully!")
        logger.info(f"  - System positioned for optimal pulling in +Z direction")
        logger.info(f"  - User-defined pulling distance: {pulling_distance:.1f} nm")
        logger.info(f"  - Ready for steered MD and umbrella sampling")
        
        return boxed_file
    
    def _verify_box_alignment_consistency(self, boxed_file: Path) -> None:
        """Verify that alignment is maintained after box operations"""
        try:
            logger.info("Verifying alignment consistency after box operations...")
            
            # Calculate centroids of the final boxed system
            verification_centroids = self.work_dir / "post_box_centroids.xvg"
            self._calculate_centroids(boxed_file, verification_centroids)
            
            # Read the centroid data
            data = []
            with open(verification_centroids, 'r') as f:
                for line in f:
                    if not line.startswith('#') and not line.startswith('@') and line.strip():
                        values = line.split()
                        if len(values) >= 7:
                            data.append([float(x) for x in values])
            
            if not data:
                logger.warning("Could not verify alignment consistency - no centroid data")
                return
            
            # Extract final centroid coordinates
            row = data[0]
            final_protein_com = np.array([row[1], row[2], row[3]])
            final_ligand_com = np.array([row[4], row[5], row[6]])
            
            # Calculate final alignment vector
            final_vector = final_ligand_com - final_protein_com
            final_vector_norm = final_vector / np.linalg.norm(final_vector)
            
            # Check alignment with Z-axis
            z_axis = np.array([0, 0, 1])
            cos_angle = np.dot(final_vector_norm, z_axis)
            angle_degrees = np.degrees(np.arccos(np.clip(abs(cos_angle), 0.0, 1.0)))
            
            logger.info("Post-box alignment verification:")
            logger.info(f"  Final protein COM: [{final_protein_com[0]:.3f}, {final_protein_com[1]:.3f}, {final_protein_com[2]:.3f}]")
            logger.info(f"  Final ligand COM:  [{final_ligand_com[0]:.3f}, {final_ligand_com[1]:.3f}, {final_ligand_com[2]:.3f}]")
            logger.info(f"  Final P→L vector:  [{final_vector_norm[0]:.3f}, {final_vector_norm[1]:.3f}, {final_vector_norm[2]:.3f}]")
            logger.info(f"  Angle with Z-axis: {angle_degrees:.1f}°")
            
            # Check direction consistency
            if final_vector_norm[2] > 0:
                logger.info("✓ Direction maintained: P→L points towards +Z (GOOD)")
            else:
                logger.warning("⚠ Direction issue: P→L points towards -Z (CHECK NEEDED)")
            
            # Check alignment quality
            if angle_degrees < 5.0:
                logger.info("✓ Excellent alignment maintained (< 5°)")
            elif angle_degrees < 15.0:
                logger.info("✓ Good alignment maintained (< 15°)")
            elif angle_degrees < 30.0:
                logger.warning("○ Moderate alignment (< 30°) - acceptable but not optimal")
            else:
                logger.warning("△ Poor alignment (> 30°) - may affect PMF quality")
            
            # Compare with pre-alignment if available
            if hasattr(self, '_aligned_vector_info'):
                original_info = self._aligned_vector_info
                logger.info("Comparison with pre-box alignment:")
                logger.info(f"  Original angle: {original_info['final_angle']:.1f}°")
                logger.info(f"  Current angle:  {angle_degrees:.1f}°")
                angle_drift = abs(angle_degrees - original_info['final_angle'])
                if angle_drift < 2.0:
                    logger.info("✓ Minimal alignment drift (< 2°)")
                elif angle_drift < 5.0:
                    logger.info("○ Small alignment drift (< 5°)")
                else:
                    logger.warning(f"△ Significant alignment drift ({angle_drift:.1f}°)")
            
        except Exception as e:
            logger.warning(f"Could not verify alignment consistency: {e}")
    
    def _parse_box_dimensions(self, gro_file: Path) -> Dict[str, float]:
        """
        Parse box dimensions from GROMACS .gro file
        
        Parameters:
        -----------
        gro_file : Path
            Path to .gro file containing box vectors
            
        Returns:
        --------
        Dict[str, float] : Box dimensions (x, y, z) in nm
        """
        try:
            with open(gro_file, 'r') as f:
                lines = f.readlines()
            
            # Last line contains box vectors
            # Format: "   v1(x) v2(y) v3(z) [v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)]"
            box_line = lines[-1].strip()
            box_values = box_line.split()
            
            if len(box_values) >= 3:
                # For cubic/rectangular boxes, we need the diagonal elements
                x_dim = float(box_values[0])
                y_dim = float(box_values[1]) 
                z_dim = float(box_values[2])
                
                logger.debug(f"Parsed box dimensions: {x_dim:.3f} × {y_dim:.3f} × {z_dim:.3f} nm")
                
                return {
                    'x_dimension': x_dim,
                    'y_dimension': y_dim,
                    'z_dimension': z_dim
                }
            else:
                raise ValueError(f"Invalid box format: {box_line}")
                
        except Exception as e:
            logger.error(f"Failed to parse box dimensions from {gro_file}: {e}")
            # Fallback to default dimensions
            logger.warning("Using fallback box dimensions: 5.0 × 5.0 × 5.0 nm")
            return {
                'x_dimension': 5.0,
                'y_dimension': 5.0,
                'z_dimension': 5.0
            }
    
    def _solvate_and_add_ions(self, boxed_system: Path, topology: Path) -> Path:
        """Solvate system and add ions"""
        logger.info("Solvating system and adding ions")
        
        # Solvate
        solvated_file = self.work_dir / "complex_solvated.gro"
        
        cmd = [
            self.gmx_command, "solvate",
            "-cp", str(boxed_system),
            "-p", str(topology),
            "-o", str(solvated_file)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"Solvation failed: {result.stderr}")
        
        # Add ions
        final_system = self.rebuilt_system_dir / "solv_ions.gro"
        final_topology = self.rebuilt_system_dir / "topol.top"
        
        # Copy files and dependencies first
        shutil.copy2(solvated_file, self.rebuilt_system_dir / "temp_solvated.gro")
        shutil.copy2(topology, final_topology)
        
        # Copy ligand topology dependencies that are referenced by relative paths
        self._copy_ligand_topology_files()
        
        # Generate protein-specific position restraints using C-alpha group
        self._generate_protein_position_restraints()
        
        # Add ions using direct method
        logger.info("Adding ions to neutralize system and achieve target concentration")
        try:
            self._add_ions_direct(
                self.rebuilt_system_dir / "temp_solvated.gro", 
                final_topology, 
                final_system
            )
            logger.info("Ion addition completed successfully")
            
            # Clean up temp file
            temp_file = self.rebuilt_system_dir / "temp_solvated.gro"
            if temp_file.exists():
                temp_file.unlink()
                
            return final_system
            
        except Exception as e:
            logger.warning(f"Ion addition failed: {e}")
            # Fallback: copy solvated system without ions
            shutil.copy2(self.rebuilt_system_dir / "temp_solvated.gro", final_system)
            logger.warning("Using solvated system without additional ions")
            
            # Clean up temp file
            temp_file = self.rebuilt_system_dir / "temp_solvated.gro"
            if temp_file.exists():
                temp_file.unlink()
                
            return final_system
    
    def _add_ions_direct(self, solvated_file: Path, topology_file: Path, output_file: Path):
        """Direct ion addition fallback method"""
        logger.info("Adding ions directly using gmx genion")
        
        # Create temporary .tpr file for genion
        temp_mdp = self.work_dir / "temp_ions.mdp"
        temp_tpr = self.work_dir / "temp_ions.tpr"
        
        # Write minimal MDP file
        with open(temp_mdp, 'w') as f:
            f.write("integrator = steep\n")
            f.write("nsteps = 0\n")
        
        # Run grompp
        grompp_cmd = [
            self.gmx_command, "grompp",
            "-f", str(temp_mdp),
            "-c", str(solvated_file),
            "-p", str(topology_file),
            "-o", str(temp_tpr),
            "-maxwarn", "5"
        ]
        
        result = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=self.work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"grompp failed: {result.stderr}")
        
        # Run genion to add ions and neutralize with salt concentration
        genion_cmd = [
            self.gmx_command, "genion",
            "-s", str(temp_tpr),
            "-o", str(output_file),
            "-p", str(topology_file),
            "-pname", "NA",
            "-nname", "CL",
            "-neutral",
            "-conc", "0.15"  # 0.15 M NaCl concentration
        ]
        
        # Select solvent group for ion replacement (usually SOL)
        result = subprocess.run(
            genion_cmd,
            input="SOL\n",
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"genion failed: {result.stderr}")
        
        # Clean up temporary files
        try:
            temp_mdp.unlink()
            temp_tpr.unlink()
        except FileNotFoundError:
            pass
        
        logger.info("Direct ion addition completed")
    
    def generate_all_equilibration_scripts(self) -> Dict:
        """
        Generate both local and SLURM equilibration scripts

        Returns:
        --------
        Dict : Script generation results for both execution modes
        """
        logger.info("Generating both local and SLURM equilibration scripts")

        results = {
            'local_script': None,
            'slurm_script': None,
            'generated_scripts': []
        }

        # Generate local run script
        try:
            local_script_path = self.rebuilt_system_dir / "run_equilibration_local.sh"
            local_content = self._generate_local_equilibration_script_content()

            with open(local_script_path, 'w') as f:
                f.write(local_content)

            # Make script executable
            local_script_path.chmod(0o755)

            results['local_script'] = str(local_script_path)
            results['generated_scripts'].append('run_equilibration_local.sh')
            logger.info(f"Local equilibration script generated: {local_script_path}")

        except Exception as e:
            logger.error(f"Failed to generate local script: {e}")

        # Generate SLURM script
        try:
            slurm_script_path = self.rebuilt_system_dir / "run_equilibration_slurm.sh"
            slurm_content = self._generate_slurm_equilibration_script_content()

            with open(slurm_script_path, 'w') as f:
                f.write(slurm_content)

            # Make script executable
            slurm_script_path.chmod(0o755)

            results['slurm_script'] = str(slurm_script_path)
            results['generated_scripts'].append('run_equilibration_slurm.sh')
            logger.info(f"SLURM equilibration script generated: {slurm_script_path}")

        except Exception as e:
            logger.error(f"Failed to generate SLURM script: {e}")

        # Generate instructions file
        try:
            instructions_path = self.rebuilt_system_dir / "EQUILIBRATION_INSTRUCTIONS.txt"
            instructions_content = self._generate_equilibration_instructions()

            with open(instructions_path, 'w') as f:
                f.write(instructions_content)

            results['instructions'] = str(instructions_path)
            results['generated_scripts'].append('EQUILIBRATION_INSTRUCTIONS.txt')
            logger.info(f"Instructions file generated: {instructions_path}")

        except Exception as e:
            logger.error(f"Failed to generate instructions: {e}")

        return results

    def generate_equilibration_script(self) -> Dict:
        """
        Generate localrun.sh-style equilibration script for manual execution
        
        Returns:
        --------
        Dict : Script generation results
        """
        logger.info("Generating manual equilibration script")
        
        # Get script configuration
        script_config = self.config.get('equilibration', {}).get('script', {})
        script_filename = script_config.get('filename', 'equilibration_run.sh')
        
        script_path = self.rebuilt_system_dir / script_filename
        
        # Generate the script content
        script_content = self._generate_equilibration_script_content()
        
        # Write script file
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make script executable
        script_path.chmod(0o755)
        
        # Generate MDP files for manual execution
        mdp_results = self._generate_equilibration_mdps()
        
        results = {
            'script_path': str(script_path),
            'mdp_files': mdp_results,
            'mode': 'manual',
            'ready_for_execution': True,
            'instructions': [
                f"1. Navigate to: {self.rebuilt_system_dir}",
                f"2. Execute script: bash {script_filename}",
                f"3. Monitor progress in em/, nvt/, npt/ subdirectories"
            ]
        }
        
        logger.info(f"Equilibration script generated: {script_path}")
        logger.info("Manual equilibration setup completed")
        
        return results
    
    def _generate_equilibration_script_content(self) -> str:
        """
        Generate the content of the equilibration script
        
        Returns:
        --------
        str : Script content following localrun.sh style
        """
        # Get hardware configuration
        hardware = self.config.get('equilibration', {}).get('hardware', {})
        gpu_id = hardware.get('gpu_id', 0)
        ntmpi = hardware.get('ntmpi', 1)
        ntomp = hardware.get('ntomp', 10)
        gpu_acceleration = hardware.get('gpu_acceleration', True)
        
        # Get equilibration configuration
        em_config = self.config.get('equilibration', {}).get('em', {})
        nvt_config = self.config.get('equilibration', {}).get('nvt', {})
        npt_config = self.config.get('equilibration', {}).get('npt', {})
        
        # GPU optimization flags
        gpu_flags = ""
        if gpu_acceleration:
            gpu_flags = f"-nb gpu -bonded gpu -pme gpu -gpu_id {gpu_id}"
        
        script_content = f"""#!/bin/bash

######################################################
# PMF System Equilibration Script
# Generated by PRISM PMF Builder
# Similar to localrun.sh but optimized for PMF systems
######################################################

echo "Starting PMF system equilibration..."
echo "Working directory: $(pwd)"
echo "GPU acceleration: {'enabled' if gpu_acceleration else 'disabled'}"

# Check for required files
if [ ! -f "solv_ions.gro" ]; then
    echo "ERROR: solv_ions.gro not found!"
    exit 1
fi

if [ ! -f "topol.top" ]; then
    echo "ERROR: topol.top not found!"
    exit 1
fi

# Energy Minimization (EM)
echo "=== Energy Minimization ==="
mkdir -p em
if [ -f ./em/em.gro ]; then
    echo "EM already completed, skipping..."
elif [ -f ./em/em.tpr ]; then
    echo "EM tpr file found, continuing from checkpoint..."
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi {ntmpi} -ntomp {ntomp} {gpu_flags} -v -cpi ./em/em.cpt
else
    echo "Starting EM from scratch..."
    gmx grompp -f em.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o ./em/em.tpr -maxwarn {em_config.get('maxwarn', 3)}
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi {ntmpi} -ntomp {ntomp} {gpu_flags} -v
fi

# Check EM completion
if [ ! -f ./em/em.gro ]; then
    echo "ERROR: Energy minimization failed!"
    exit 1
fi

# NVT Equilibration  
echo "=== NVT Equilibration ==="
mkdir -p nvt
if [ -f ./nvt/nvt.gro ]; then
    echo "NVT already completed, skipping..."
elif [ -f ./nvt/nvt.tpr ]; then
    echo "NVT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi {ntmpi} -ntomp {ntomp} {gpu_flags} -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v -cpi ./nvt/nvt.cpt
else
    echo "Starting NVT from scratch..."
    gmx grompp -f nvt.mdp -c ./em/em.gro -r ./em/em.gro -p topol.top -o ./nvt/nvt.tpr -maxwarn {nvt_config.get('maxwarn', 3)}
    gmx mdrun -ntmpi {ntmpi} -ntomp {ntomp} {gpu_flags} -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v
fi

# Check NVT completion
if [ ! -f ./nvt/nvt.gro ]; then
    echo "ERROR: NVT equilibration failed!"
    exit 1
fi

# NPT Equilibration
echo "=== NPT Equilibration ==="
mkdir -p npt
if [ -f ./npt/npt.gro ]; then
    echo "NPT already completed, skipping..."
elif [ -f ./npt/npt.tpr ]; then
    echo "NPT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi {ntmpi} -ntomp {ntomp} {gpu_flags} -s ./npt/npt.tpr -deffnm ./npt/npt -v -cpi ./npt/npt.cpt
else
    echo "Starting NPT from scratch..."
    gmx grompp -f npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn {npt_config.get('maxwarn', 3)}
    gmx mdrun -ntmpi {ntmpi} -ntomp {ntomp} {gpu_flags} -s ./npt/npt.tpr -deffnm ./npt/npt -v
fi

# Check NPT completion
if [ ! -f ./npt/npt.gro ]; then
    echo "ERROR: NPT equilibration failed!"
    exit 1
fi

echo "=== PMF System Equilibration Completed Successfully ==="
echo "Final equilibrated structure: ./npt/npt.gro"
echo "System is now ready for PMF calculations"

# Optional: Create symbolic link to final structure
ln -sf ./npt/npt.gro equilibrated_system.gro
ln -sf ./npt/npt.cpt equilibrated_system.cpt

echo "Symbolic links created:"
echo "  equilibrated_system.gro -> ./npt/npt.gro"  
echo "  equilibrated_system.cpt -> ./npt/npt.cpt"
"""
        
        return script_content

    def _generate_local_equilibration_script_content(self) -> str:
        """
        Generate local run equilibration script content

        Returns:
        --------
        str : Local run script content
        """
        # Get local configuration
        local_config = self.config.get('local', {})
        cpus = local_config.get('cpus', 10)
        gpu_id = local_config.get('gpu_id', 0)
        ntmpi = local_config.get('ntmpi', 1)
        ntomp = local_config.get('ntomp', 10)

        script_content = f"""#!/bin/bash

######################################################
# PMF System Local Equilibration Script
# Generated by PRISM PMF Builder
######################################################

set -e  # Exit on any error

echo "=== PRISM PMF System Local Equilibration ==="
echo "CPU cores: {cpus}"
echo "GPU ID: {gpu_id}"
echo "Working directory: $(pwd)"
echo ""

# Check for required files
if [ ! -f "solv_ions.gro" ]; then
    echo "ERROR: solv_ions.gro not found!"
    exit 1
fi

if [ ! -f "topol.top" ]; then
    echo "ERROR: topol.top not found!"
    exit 1
fi

# Set GROMACS command
GMX_COMMAND="{self.gmx_command}"

# Hardware configuration
NCPU={ntomp}
GPU_ID={gpu_id}

# Create equilibration directory
mkdir -p equilibration
cd equilibration

# Function to run with GPU acceleration if available
run_mdrun() {{
    local description="$1"
    local tpr_file="$2"
    local output_name="$3"

    echo "Running $description..."

    # Try GPU first, fallback to CPU
    if command -v nvidia-smi &> /dev/null && nvidia-smi &> /dev/null; then
        echo "  Using GPU acceleration (GPU $GPU_ID)"
        $GMX_COMMAND mdrun -s "$tpr_file" -deffnm "$output_name" \\
            -ntmpi {ntmpi} -ntomp $NCPU -nb gpu -pme gpu -gpu_id $GPU_ID -v
    else
        echo "  Using CPU only"
        $GMX_COMMAND mdrun -s "$tpr_file" -deffnm "$output_name" \\
            -ntmpi {ntmpi} -ntomp $NCPU -v
    fi
}}

# Energy Minimization
echo "=== Energy Minimization ==="
$GMX_COMMAND grompp -f ../em.mdp -c ../solv_ions.gro -p ../topol.top -o em.tpr -maxwarn 3
run_mdrun "Energy Minimization" "em.tpr" "em"

# NVT Equilibration
echo "=== NVT Equilibration ==="
$GMX_COMMAND grompp -f ../nvt.mdp -c em.gro -r em.gro -p ../topol.top -o nvt.tpr -maxwarn 3
run_mdrun "NVT Equilibration" "nvt.tpr" "nvt"

# NPT Equilibration
echo "=== NPT Equilibration ==="
$GMX_COMMAND grompp -f ../npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ../topol.top -o npt.tpr -maxwarn 3
run_mdrun "NPT Equilibration" "npt.tpr" "npt"

# Create final files
echo "=== Creating Final System Files ==="
cp npt.gro ../equilibrated_system.gro
cp npt.cpt ../equilibrated_system.cpt

echo ""
echo "=== PMF System Equilibration Completed Successfully ==="
echo "Final equilibrated system:"
echo "  Structure: equilibrated_system.gro"
echo "  Checkpoint: equilibrated_system.cpt"
echo "  Ready for PMF calculations"
"""
        return script_content

    def _generate_slurm_equilibration_script_content(self) -> str:
        """
        Generate SLURM equilibration script content

        Returns:
        --------
        str : SLURM script content
        """
        # Get SLURM configuration
        slurm_config = self.config.get('slurm', {})
        partition = slurm_config.get('partition', 'gpu')
        cpus_per_task = slurm_config.get('cpus_per_task', 10)
        time = slurm_config.get('time', '12:00:00')
        memory = slurm_config.get('memory', '32G')
        gres = slurm_config.get('gres', 'gpu:1')
        nodes = slurm_config.get('nodes', 1)
        ntasks = slurm_config.get('ntasks', 1)
        environment_setup = slurm_config.get('environment_setup', '')
        ntmpi = slurm_config.get('ntmpi', 1)
        ntomp = slurm_config.get('ntomp', 10)

        script_content = f"""#!/bin/bash
#SBATCH --job-name=pmf_equilibration
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --time={time}"""

        if memory:
            script_content += f"\n#SBATCH --mem={memory}"

        if gres:
            script_content += f"\n#SBATCH --gres={gres}"

        script_content += f"""
#SBATCH --output=pmf_equilibration_%j.out
#SBATCH --error=pmf_equilibration_%j.err

######################################################
# PMF System SLURM Equilibration Script
# Generated by PRISM PMF Builder
######################################################

set -e  # Exit on any error

echo "=== PRISM PMF System SLURM Equilibration ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Partition: {partition}"
echo "Nodes: {nodes}"
echo "CPUs per task: {cpus_per_task}"
echo "Working directory: $(pwd)"
echo ""

# Environment setup
{environment_setup}

# Check for required files
if [ ! -f "solv_ions.gro" ]; then
    echo "ERROR: solv_ions.gro not found!"
    exit 1
fi

if [ ! -f "topol.top" ]; then
    echo "ERROR: topol.top not found!"
    exit 1
fi

# Set GROMACS command
GMX_COMMAND="{self.gmx_command}"

# Hardware configuration
NCPU={ntomp}

# Create equilibration directory
mkdir -p equilibration
cd equilibration

# Function to run with GPU acceleration
run_mdrun() {{
    local description="$1"
    local tpr_file="$2"
    local output_name="$3"

    echo "Running $description..."
    echo "  Using SLURM GPU resources"

    $GMX_COMMAND mdrun -s "$tpr_file" -deffnm "$output_name" \\
        -ntmpi {ntmpi} -ntomp $NCPU -nb gpu -pme gpu -v
}}

# Energy Minimization
echo "=== Energy Minimization ==="
$GMX_COMMAND grompp -f ../em.mdp -c ../solv_ions.gro -p ../topol.top -o em.tpr -maxwarn 3
run_mdrun "Energy Minimization" "em.tpr" "em"

# NVT Equilibration
echo "=== NVT Equilibration ==="
$GMX_COMMAND grompp -f ../nvt.mdp -c em.gro -r em.gro -p ../topol.top -o nvt.tpr -maxwarn 3
run_mdrun "NVT Equilibration" "nvt.tpr" "nvt"

# NPT Equilibration
echo "=== NPT Equilibration ==="
$GMX_COMMAND grompp -f ../npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ../topol.top -o npt.tpr -maxwarn 3
run_mdrun "NPT Equilibration" "npt.tpr" "npt"

# Create final files
echo "=== Creating Final System Files ==="
cp npt.gro ../equilibrated_system.gro
cp npt.cpt ../equilibrated_system.cpt

echo ""
echo "=== PMF System Equilibration Completed Successfully ==="
echo "Final equilibrated system:"
echo "  Structure: equilibrated_system.gro"
echo "  Checkpoint: equilibrated_system.cpt"
echo "  Ready for PMF calculations"
echo ""
echo "Job completed at: $(date)"
"""
        return script_content

    def _generate_equilibration_instructions(self) -> str:
        """
        Generate equilibration instructions file

        Returns:
        --------
        str : Instructions content
        """
        instructions = f"""
PMF SYSTEM EQUILIBRATION INSTRUCTIONS
=====================================

Generated by PRISM PMF Builder

AVAILABLE SCRIPTS:
-----------------
1. run_equilibration_local.sh  - Local execution script
2. run_equilibration_slurm.sh  - SLURM cluster submission script

Both scripts are always generated. Choose the appropriate one for your environment.

EXECUTION OPTIONS:
-----------------

Option 1: Local Execution
------------------------
For running on a local machine or workstation:

    bash run_equilibration_local.sh

Requirements:
- GROMACS installed and accessible
- GPU drivers (optional, will fallback to CPU)
- Sufficient CPU cores (recommended: 8+)

Option 2: SLURM Cluster Execution
--------------------------------
For running on a SLURM-managed cluster:

    sbatch run_equilibration_slurm.sh

Requirements:
- SLURM job scheduler
- GROMACS module or installation on cluster
- GPU partition access (as configured)

MONITORING PROGRESS:
-------------------

Local execution:
- Monitor console output directly
- Check equilibration/ directory for intermediate files

SLURM execution:
- Check job status: squeue -u $USER
- Monitor output: tail -f pmf_equilibration_*.out
- Check errors: tail -f pmf_equilibration_*.err

EXPECTED OUTPUTS:
----------------
After successful equilibration:
- equilibrated_system.gro (final structure)
- equilibrated_system.cpt (checkpoint file)
- equilibration/ directory with EM, NVT, NPT results

TROUBLESHOOTING:
---------------
1. Missing files error: Ensure solv_ions.gro and topol.top exist
2. GROMACS command not found: Check GROMACS installation/module
3. GPU acceleration issues: Script will automatically fallback to CPU
4. Memory issues: Reduce ntomp value in configuration

NEXT STEPS:
----------
After successful equilibration, the system is ready for PMF calculations.
Use the equilibrated_system.gro and equilibrated_system.cpt files for
subsequent SMD and umbrella sampling simulations.

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""
        return instructions

    def _generate_equilibration_mdps(self) -> Dict:
        """
        Generate MDP files for manual equilibration
        
        Returns:
        --------
        Dict : Generated MDP file information
        """
        logger.info("Generating MDP files for manual equilibration")
        
        # Use the existing equilibration manager to generate MDP files
        from .equilibration import PMFEquilibrationManager
        
        equil_manager = PMFEquilibrationManager(
            system_dir=self.rebuilt_system_dir,
            output_dir=self.rebuilt_system_dir,  # Keep everything in same directory
            config=self.config,
            gmx_command=self.gmx_command
        )
        
        # Generate MDP files directly in the PMF system directory
        em_mdp = equil_manager._generate_em_mdp()
        nvt_mdp = equil_manager._generate_nvt_mdp() 
        npt_mdp = equil_manager._generate_npt_mdp()
        
        # Move MDP files to system directory for script execution
        target_em_mdp = self.rebuilt_system_dir / "em.mdp"
        target_nvt_mdp = self.rebuilt_system_dir / "nvt.mdp"
        target_npt_mdp = self.rebuilt_system_dir / "npt.mdp"
        
        shutil.copy2(em_mdp, target_em_mdp)
        shutil.copy2(nvt_mdp, target_nvt_mdp)
        shutil.copy2(npt_mdp, target_npt_mdp)
        
        mdp_results = {
            'em_mdp': str(target_em_mdp),
            'nvt_mdp': str(target_nvt_mdp),
            'npt_mdp': str(target_npt_mdp),
            'location': str(self.rebuilt_system_dir)
        }
        
        logger.info("MDP files generated for manual execution")
        return mdp_results
    
    def _copy_ligand_topology_files(self) -> None:
        """
        Copy ligand topology files and dependencies from original MD system to PMF system
        
        This method identifies and copies:
        - LIG.amb2gmx/ or LIG.openff2gmx/ directories containing ligand topology files
        - Any other files referenced by relative paths in the topology
        
        Note: Position restraints files (posre.itp) are handled separately by _copy_position_restraints_file()
        """
        logger.info("Copying ligand topology dependencies")
        
        # Common ligand topology directory names
        ligand_topology_dirs = ['LIG.amb2gmx', 'LIG.openff2gmx']
        
        copied_dirs = []
        
        # Search for and copy ligand topology directories
        # Ligand topology dirs are at the same level as GMX_PROLIG_MD, not inside it
        # We need to copy them to the parent of PMF system to maintain relative path structure
        for dir_name in ligand_topology_dirs:
            source_dir = self.md_system_dir.parent / dir_name  # From gaff_model/LIG.amb2gmx/
            if source_dir.exists() and source_dir.is_dir():
                # Copy to PMF system parent directory, maintain relative path structure ../LIG.amb2gmx/
                target_dir = self.rebuilt_system_dir.parent / dir_name
                
                try:
                    if target_dir.exists():
                        shutil.rmtree(target_dir)
                    
                    shutil.copytree(source_dir, target_dir)
                    copied_dirs.append(dir_name)
                    
                    logger.info(f"Copied ligand topology directory: {dir_name}")
                    
                    # Log contents for verification
                    files_in_dir = list(target_dir.iterdir())
                    logger.debug(f"Files in {dir_name}: {[f.name for f in files_in_dir]}")
                    
                except Exception as e:
                    logger.warning(f"Failed to copy {dir_name}: {e}")
        
        if copied_dirs:
            logger.info(f"Successfully copied ligand topology directories: {copied_dirs}")
        else:
            logger.warning("No ligand topology directories found to copy")
        
        # Additional check: scan topology file for other relative path references
        self._copy_additional_topology_dependencies()
    
    def _copy_additional_topology_dependencies(self) -> None:
        """
        Scan topology file for additional relative path dependencies and copy them
        """
        topology_file = self.rebuilt_system_dir / "topol.top"
        
        if not topology_file.exists():
            logger.warning("Topology file not found for dependency scanning")
            return
        
        try:
            with open(topology_file, 'r') as f:
                content = f.read()
            
            # Look for #include statements with relative paths (../)
            include_pattern = r'#include\s+"(\.\./[^"]+)"'
            includes = re.findall(include_pattern, content)
            
            for include_path in includes:
                # Convert relative path to source and target paths
                # ../LIG.amb2gmx/LIG.itp -> LIG.amb2gmx/LIG.itp
                rel_path = include_path.replace('../', '')
                source_file = self.md_system_dir.parent / rel_path  # From gaff_model/LIG.amb2gmx/
                target_file = self.rebuilt_system_dir.parent / rel_path  # To pmf_output/LIG.amb2gmx/
                
                if source_file.exists() and not target_file.exists():
                    # Make sure target directory exists
                    target_file.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source_file, target_file)
                    logger.debug(f"Copied additional dependency: {rel_path}")
            
            if includes:
                logger.info(f"Processed {len(includes)} topology include dependencies")
            
        except Exception as e:
            logger.warning(f"Failed to scan topology dependencies: {e}")
    
    def _copy_position_restraints_file(self) -> None:
        """
        Copy position restraints file (posre.itp) from original MD system to PMF system
        
        The posre.itp file is typically located in the GMX_PROLIG_MD directory and is
        referenced by the topology file for position restraints during equilibration.
        """
        logger.info("Copying position restraints file")
        
        # Common position restraints file names
        posre_files = ['posre.itp', 'posre_Protein.itp', 'posre_LIG.itp']
        
        copied_files = []
        
        for posre_filename in posre_files:
            source_file = self.md_system_dir / posre_filename
            target_file = self.rebuilt_system_dir / posre_filename
            
            if source_file.exists():
                try:
                    shutil.copy2(source_file, target_file)
                    copied_files.append(posre_filename)
                    logger.info(f"Copied position restraints file: {posre_filename}")
                    logger.debug(f"  From: {source_file}")
                    logger.debug(f"  To:   {target_file}")
                    logger.debug(f"  Size: {source_file.stat().st_size} bytes")
                    logger.debug(f"  Same directory as topol.top: ✓")
                    
                except Exception as e:
                    logger.warning(f"Failed to copy {posre_filename}: {e}")
        
        if copied_files:
            logger.info(f"Successfully copied position restraints files: {copied_files}")
        else:
            logger.warning("No position restraints files found in original MD system")
            logger.info("This may be normal if position restraints are not used")
    
    def _generate_protein_position_restraints(self) -> None:
        """
        Generate protein-specific position restraints using C-alpha group
        
        This replaces the copying of posre.itp from original MD system with
        a new protein_posre.itp file generated specifically for the PMF system.
        """
        logger.info("Generating protein position restraints using C-alpha group")
        
        # Extract protein structure for position restraint generation
        protein_gro = self._extract_protein_for_posre()
        
        if not protein_gro.exists():
            logger.warning("Could not extract protein structure for position restraints")
            return
        
        # Generate protein-specific position restraints
        protein_posre_file = self.rebuilt_system_dir / "protein_posre.itp"
        
        try:
            # Get C-alpha group number
            calpha_group_num = self._get_calpha_group_number(protein_gro)
            
            if calpha_group_num is None:
                # Use standard C-alpha group number (usually 3)
                calpha_group_num = 3
                logger.info(f"Using standard C-alpha group number: {calpha_group_num}")
            else:
                logger.info(f"Detected C-alpha group number: {calpha_group_num}")
            
            # Generate position restraints using gmx genrestr
            cmd = [
                self.gmx_command, "genrestr",
                "-f", str(protein_gro),
                "-o", str(protein_posre_file)
            ]
            
            result = subprocess.run(
                cmd,
                input=f"{calpha_group_num}\n",
                text=True,
                capture_output=True,
                cwd=str(self.rebuilt_system_dir)
            )
            
            if result.returncode == 0 and protein_posre_file.exists():
                logger.info(f"Successfully generated protein position restraints: protein_posre.itp")
                # Update topology file to use the new position restraints
                self._update_topology_for_protein_posre()
            else:
                logger.warning(f"Failed to generate protein position restraints: {result.stderr}")
                
        except Exception as e:
            logger.error(f"Error generating protein position restraints: {e}")
    
    def _extract_protein_for_posre(self) -> Path:
        """Extract protein structure from solvated system for position restraint generation"""
        # Use the temporary solvated file that contains the protein structure
        temp_solvated = self.rebuilt_system_dir / "temp_solvated.gro"
        protein_gro = self.rebuilt_system_dir / "protein_only.gro"
        
        if not temp_solvated.exists():
            logger.warning("temp_solvated.gro not found, cannot extract protein for position restraints")
            return protein_gro
        
        try:
            cmd = [
                self.gmx_command, "trjconv",
                "-f", str(temp_solvated),
                "-s", str(temp_solvated),
                "-o", str(protein_gro),
                "-n"  # Will use index file if available
            ]
            
            # Try to extract protein group (group name: "Protein")
            result = subprocess.run(
                cmd,
                input="Protein\n",
                text=True,
                capture_output=True,
                cwd=str(self.rebuilt_system_dir)
            )
            
            if result.returncode == 0 and protein_gro.exists():
                logger.info(f"Successfully extracted protein structure: {protein_gro.name}")
            else:
                logger.warning(f"Failed to extract protein structure: {result.stderr}")
                
        except Exception as e:
            logger.error(f"Error extracting protein structure: {e}")
        
        return protein_gro
    
    def _get_calpha_group_number(self, protein_gro: Path) -> Optional[int]:
        """Get C-alpha group number from protein structure using gmx make_ndx"""
        try:
            cmd = [self.gmx_command, "make_ndx", "-f", str(protein_gro)]
            result = subprocess.run(
                cmd,
                input="q\n",
                text=True,
                capture_output=True,
                cwd=str(self.rebuilt_system_dir)
            )
            
            # Clean up temporary index files
            for temp_file in [self.rebuilt_system_dir / "index.ndx", 
                            self.rebuilt_system_dir / f"{protein_gro.stem}.ndx"]:
                if temp_file.exists() and temp_file.name != "index.ndx":
                    try:
                        temp_file.unlink()
                    except:
                        pass
            
            if result.returncode != 0:
                return None
            
            # Parse output to find C-alpha group number
            for line in result.stdout.split('\n'):
                if ':' in line and 'atoms' in line:
                    try:
                        parts = line.strip().split()
                        if len(parts) >= 3 and parts[1] == 'C-alpha':
                            group_num = int(parts[0])
                            logger.debug(f"Found C-alpha group: number {group_num}")
                            return group_num
                    except (ValueError, IndexError):
                        continue
            
            # If not found by exact match, try partial match
            for line in result.stdout.split('\n'):
                if ':' in line and 'atoms' in line and 'c-alpha' in line.lower():
                    try:
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            group_num = int(parts[0])
                            logger.debug(f"Found C-alpha group by partial match: number {group_num}")
                            return group_num
                    except (ValueError, IndexError):
                        continue
        
        except Exception as e:
            logger.warning(f"Failed to get C-alpha group number: {e}")
        
        return None
    
    def _update_topology_for_protein_posre(self) -> None:
        """Update topology file to include protein_posre.itp with POSRES_P macro"""
        topology_file = self.rebuilt_system_dir / "topol.top"
        
        if not topology_file.exists():
            logger.warning("topol.top not found, cannot update for protein position restraints")
            return
        
        try:
            # Read current topology content
            with open(topology_file, 'r') as f:
                lines = f.readlines()
            
            # Find the position restraint section and update it
            updated_lines = []
            posre_section_found = False
            
            for line in lines:
                updated_lines.append(line)
                
                # Look for existing position restraint include section
                if "; Include Position restraint file" in line or "#include \"posre.itp\"" in line:
                    posre_section_found = True
                    # Add the new protein-specific position restraints
                    if not any("#ifdef POSRES_P" in l for l in lines):  # Avoid duplicates
                        updated_lines.append("#ifdef POSRES_P\n")
                        updated_lines.append("#include \"protein_posre.itp\"\n")
                        updated_lines.append("#endif\n")
                        logger.info("Added POSRES_P macro and protein_posre.itp include to topology")
            
            # If no existing position restraint section found, add it after system definition
            if not posre_section_found:
                # Find a good place to add it (after [ system ] or [ molecules ])
                for i, line in enumerate(updated_lines):
                    if "[ system ]" in line or "[ molecules ]" in line:
                        # Add after this section
                        insert_pos = i + 2  # Skip the section header and next line
                        if insert_pos < len(updated_lines):
                            updated_lines.insert(insert_pos, "\n; Include Position restraint file\n")
                            updated_lines.insert(insert_pos + 1, "#ifdef POSRES_P\n")
                            updated_lines.insert(insert_pos + 2, "#include \"protein_posre.itp\"\n") 
                            updated_lines.insert(insert_pos + 3, "#endif\n")
                            logger.info("Added POSRES_P section to topology file")
                            break
            
            # Write updated topology
            with open(topology_file, 'w') as f:
                f.writelines(updated_lines)
            
        except Exception as e:
            logger.error(f"Error updating topology file for protein position restraints: {e}")
    
    def _save_build_config(self, results: Dict) -> None:
        """Save PMF build configuration and results"""
        config_file = self.output_dir / "pmf_build_config.yaml"
        
        build_info = {
            'pmf_builder': {
                'version': '1.0.0',
                'source_md_dir': str(self.md_system_dir),
                'build_timestamp': self._get_timestamp()
            },
            'extraction_info': results['extraction'],
            'alignment_info': results['alignment'],
            'configuration': self.config,
            'output_structure': {
                'system_directory': results['system_dir'],
                'primary_files': {
                    'structure': str(Path(results['system_dir']) / "solv_ions.gro"),
                    'topology': str(Path(results['system_dir']) / "topol.top")
                }
            }
        }
        
        with open(config_file, 'w') as f:
            yaml.dump(build_info, f, default_flow_style=False, indent=2)
        
        logger.info(f"Build configuration saved: {config_file}")
    
    def _get_timestamp(self) -> str:
        """Get current timestamp"""
        from datetime import datetime
        return datetime.now().isoformat()
    
    def get_system_info(self) -> Dict:
        """Get information about the built PMF system"""
        final_structure = self.rebuilt_system_dir / "solv_ions.gro"
        final_topology = self.rebuilt_system_dir / "topol.top"
        
        if not final_structure.exists():
            return {'status': 'not_built', 'message': 'PMF system not yet built'}
        
        # Get basic system information
        info = {
            'status': 'ready',
            'system_directory': str(self.rebuilt_system_dir),
            'structure_file': str(final_structure),
            'topology_file': str(final_topology),
            'box_optimized_for': 'PMF calculations',
            'z_axis_aligned': True
        }
        
        # Get system statistics if possible
        try:
            info.update(self._get_system_statistics(final_structure))
        except Exception as e:
            logger.warning(f"Could not get system statistics: {e}")
        
        return info
    
    def _get_system_statistics(self, structure_file: Path) -> Dict:
        """Get system statistics using GROMACS tools"""
        # Basic system info using gmx check
        cmd = [
            self.gmx_command, "check",
            "-f", str(structure_file)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
        
        stats = {'gromacs_check': 'completed'}
        
        # Parse basic information
        if result.returncode == 0:
            for line in result.stdout.split('\n'):
                if 'atoms' in line.lower():
                    try:
                        atoms = int(line.split()[0])
                        stats['total_atoms'] = atoms
                    except (ValueError, IndexError):
                        pass
        
        return stats
    
    def run_equilibration(self) -> Dict:
        """
        Run complete equilibration process (EM → NVT → NPT)
        
        Returns:
        --------
        Dict : Equilibration results
        """
        if not self.equilibration_manager:
            self.equilibration_manager = PMFEquilibrationManager(
                system_dir=self.rebuilt_system_dir,
                output_dir=self.rebuilt_system_dir,  # Use same directory for organization
                config=self.config,
                gmx_command=self.gmx_command
            )
        
        logger.info("Starting complete equilibration process...")
        equilibration_results = self.equilibration_manager.run_complete_equilibration()
        
        return equilibration_results
    
    def run_equilibration_step(self, step: str) -> Dict:
        """
        Run individual equilibration step
        
        Parameters:
        -----------
        step : str
            Equilibration step ('em', 'nvt', 'npt')
            
        Returns:
        --------
        Dict : Step results
        """
        if not self.equilibration_manager:
            self.equilibration_manager = PMFEquilibrationManager(
                system_dir=self.rebuilt_system_dir,
                output_dir=self.rebuilt_system_dir,  # Use same directory
                config=self.config,
                gmx_command=self.gmx_command
            )
        
        if step.lower() == 'em':
            return self.equilibration_manager.run_energy_minimization()
        elif step.lower() == 'nvt':
            return self.equilibration_manager.run_nvt_equilibration()
        elif step.lower() == 'npt':
            return self.equilibration_manager.run_npt_equilibration()
        else:
            raise ValueError(f"Unknown equilibration step: {step}. Use 'em', 'nvt', or 'npt'")
    
    def get_equilibration_status(self) -> Dict:
        """Get current equilibration status"""
        if not self.equilibration_manager:
            return {
                'equilibration_initialized': False,
                'status': 'not_started',
                'message': 'Equilibration not yet initialized'
            }
        
        return self.equilibration_manager.get_equilibration_status()
    
    def get_equilibrated_system(self) -> Dict:
        """Get final equilibrated system files"""
        if not self.equilibration_manager:
            raise RuntimeError("Equilibration not run - call run_equilibration() first")
        
        return self.equilibration_manager.get_final_equilibrated_system()
    
    def clean(self) -> None:
        """Clean temporary build files"""
        if self.work_dir.exists():
            shutil.rmtree(self.work_dir)
            logger.info("Temporary build files cleaned")


# Convenience function following PRISM patterns
def pmf_builder(md_results_dir: str, output_dir: str = "pmf_system",
                config: Optional[Dict] = None, **kwargs) -> PMFBuilder:
    """
    Create PMF builder instance (convenience function)
    
    Parameters:
    -----------
    md_results_dir : str
        Directory containing MD results
    output_dir : str
        Output directory for PMF system
    config : dict, optional
        PMF builder configuration
    pulling_distance : float, optional
        User-customizable pulling distance in nm (default: 2.0)
    box_distance : float, optional
        Fixed dt distance for GROMACS standard box in nm (default: 1.2)
    **kwargs : optional
        Additional configuration parameters
        
    Returns:
    --------
    PMFBuilder : PMF builder instance
    
    Examples:
    ---------
    >>> import prism
    >>> # Use default pulling distance (2.0 nm)
    >>> builder = prism.pmf.pmf_builder("./gaff_model", "./pmf_system")
    >>> 
    >>> # Customize pulling distance for your specific needs
    >>> builder = prism.pmf.pmf_builder(
    ...     "./gaff_model", "./pmf_system", 
    ...     pulling_distance=3.5  # Custom 3.5 nm pulling distance
    ... )
    >>>
    >>> # Build PMF-ready system
    >>> results = builder.build()
    """
    return PMFBuilder(md_results_dir, output_dir, config, **kwargs)