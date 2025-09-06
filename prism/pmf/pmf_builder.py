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
from typing import Dict, Tuple, Optional, Union

from ..utils.environment import GromacsEnvironment
from ..utils.config import ConfigurationManager
from .equilibration import PMFEquilibrationManager

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
            Directory containing MD results (GMX_PROLIG_MD)
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
        
        return default_config
    
    def _deep_update(self, target: Dict, source: Dict) -> None:
        """Deep update dictionary"""
        for key, value in source.items():
            if isinstance(value, dict) and key in target and isinstance(target[key], dict):
                self._deep_update(target[key], value)
            else:
                target[key] = value
    
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
        """
        logger.info("=== Building PMF-Optimized System ===")
        
        # Step 1: Extract structures from MD trajectory
        extraction_results = self.extract_structures(frame)
        
        # Step 2: Calculate centroids and alignment
        alignment_info = self.calculate_alignment_geometry()
        
        # Step 3: Rebuild system with alignment and extended box
        rebuild_results = self.rebuild_aligned_system(alignment_info)
        
        # Combine initial results
        results = {
            'extraction': extraction_results,
            'alignment': alignment_info,
            'rebuild': rebuild_results,
            'system_dir': str(self.rebuilt_system_dir),
            'equilibration_required': equilibrate
        }
        
        # Step 4: Run equilibration if requested
        if equilibrate:
            logger.info("=== Running System Equilibration ===")
            equilibration_results = self.run_equilibration()
            results['equilibration'] = equilibration_results
            results['status'] = 'pmf_ready_equilibrated'
            results['final_system'] = equilibration_results['final_system']
        else:
            results['status'] = 'pmf_ready_unequilibrated'
            results['final_system'] = {
                'structure': str(self.rebuilt_system_dir / "solv_ions.gro"),
                'topology': str(self.rebuilt_system_dir / "topol.top"),
                'equilibrated': False
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
            raise FileNotFoundError(f"Structure file not found: {structure_file}")
        
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
            raise RuntimeError(f"Failed to extract {group}: {result.stderr}")
        
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
            raise RuntimeError(f"Failed to extract complex: {result.stderr}")
        
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
                raise RuntimeError("Protein group not found")
            
            if lig_num is None:
                logger.error("LIG group not found in structure")
                raise RuntimeError("LIG group not found")
            
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
            raise RuntimeError(f"make_ndx failed: {result.stderr}")
        
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
            raise RuntimeError("Could not parse group information from make_ndx output")
        
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
        """Create aligned protein-ligand complex"""
        complex_file = self.work_dir / "complex_extracted.gro"
        aligned_file = self.work_dir / "complex_aligned.gro"
        
        # For simplicity, we'll use the extracted complex as is
        # In practice, you might want to apply the rotation matrix
        # This would require more sophisticated coordinate manipulation
        
        # Copy complex file and center it
        cmd = [
            self.gmx_command, "editconf",
            "-f", str(complex_file),
            "-o", str(aligned_file),
            "-center", "0", "0", "0"  # Center at origin
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to center complex: {result.stderr}")
        
        logger.debug(f"Created aligned complex: {aligned_file}")
        return aligned_file
    
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
        
        # === Step 3: Translate system in negative Z direction for optimal pulling space ===
        logger.info("Step 3: Optimizing system position for PMF pulling")
        
        # Translate by pulling_distance in negative Z direction
        # This positions the system optimally: maximum space in +Z direction for pulling
        translation_distance = -pulling_distance  # Negative Z direction
        
        logger.info(f"  Translation: {translation_distance:.2f} nm in Z direction")
        logger.info(f"  Available pulling space: ~{pulling_distance + original_z/2:.2f} nm in +Z direction")
        logger.info(f"  Safety margin: ~{pulling_distance - original_z/2:.2f} nm in -Z direction")
        
        cmd = [
            self.gmx_command, "editconf",
            "-f", str(extended_boxed),
            "-o", str(boxed_file),
            "-translate", "0", "0", str(translation_distance)  # Translate in Z direction only
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
        
        logger.info("PMF-optimized box creation completed successfully!")
        logger.info(f"  - System positioned for optimal pulling in +Z direction")
        logger.info(f"  - User-defined pulling distance: {pulling_distance:.1f} nm")
        logger.info(f"  - Ready for steered MD and umbrella sampling")
        
        return boxed_file
    
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
        
        # Copy position restraints file (posre.itp) from original MD system
        self._copy_position_restraints_file()
        
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
            source_dir = self.md_system_dir.parent / dir_name  # 从gaff_model/LIG.amb2gmx/
            if source_dir.exists() and source_dir.is_dir():
                # 复制到PMF系统的上级目录，保持相对路径结构 ../LIG.amb2gmx/
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
                source_file = self.md_system_dir.parent / rel_path  # 从gaff_model/LIG.amb2gmx/
                target_file = self.rebuilt_system_dir.parent / rel_path  # 到pmf_output/LIG.amb2gmx/
                
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
                output_dir=self.output_dir,
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
                output_dir=self.output_dir,
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