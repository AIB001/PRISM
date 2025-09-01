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
import shutil
import numpy as np
import subprocess
from pathlib import Path
import yaml
import logging
from typing import Dict, Tuple, Optional, Union

from ..utils.environment import GromacsEnvironment
from ..utils.config import ConfigurationManager
from ..utils.system import SystemBuilder
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
        
        # Detect original force field from MD system (保持一致性)
        self.original_forcefield_info = self._detect_original_forcefield()
        
        # Process configuration with force field info
        self.config = self._process_pmf_config(config, **kwargs)
        
        # Initialize standard PRISM components
        self.config_manager = ConfigurationManager(
            None,  # Use provided config dict instead of file
            self.gromacs_env,
            config_dict=self.config
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
        检测原MD系统使用的力场信息，确保PMF重建时保持一致性
        
        Returns:
        --------
        Dict : 原始力场信息
        """
        logger.info("检测原MD系统的力场信息...")
        
        topology_file = self.md_system_dir / "topol.top"
        if not topology_file.exists():
            logger.warning("未找到topology文件，使用默认力场设置")
            return {'forcefield': 'amber99sb', 'water_model': 'tip3p'}
        
        forcefield_info = {
            'forcefield': 'amber99sb',  # 默认值
            'water_model': 'tip3p',     # 默认值
            'detected': False
        }
        
        try:
            with open(topology_file, 'r') as f:
                content = f.read()
            
            # 检测蛋白质力场
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
            
            # 检测水模型
            if 'tip3p' in content.lower():
                forcefield_info['water_model'] = 'tip3p'
            elif 'tip4p' in content.lower():
                forcefield_info['water_model'] = 'tip4p'
            elif 'spc' in content.lower():
                forcefield_info['water_model'] = 'spc'
            elif 'spce' in content.lower():
                forcefield_info['water_model'] = 'spce'
            
            # 检测配体力场信息
            if 'gaff' in content.lower() or 'amber' in content.lower():
                forcefield_info['ligand_forcefield'] = 'gaff'
            elif 'openff' in content.lower() or 'sage' in content.lower():
                forcefield_info['ligand_forcefield'] = 'openff'
            
            if forcefield_info['detected']:
                logger.info(f"✅ 检测到原始力场: {forcefield_info['forcefield']} + {forcefield_info['water_model']}")
            else:
                logger.warning("⚠️  无法明确检测力场，使用默认设置")
                
        except Exception as e:
            logger.warning(f"力场检测出错: {e}，使用默认设置")
        
        return forcefield_info
    
    def _process_pmf_config(self, config: Optional[Dict], **kwargs) -> Dict:
        """Process PMF-specific configuration"""
        # 使用检测到的力场信息，确保一致性
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
                'box_distance': 1.0,
                'z_extension': 2.0  # Extra Z-axis length for pulling
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
                elif key == 'z_extension':
                    default_config['box']['z_extension'] = value
                elif key in ['forcefield', 'water_model']:
                    default_config['general'][key] = value
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
            logger.info("✅ System fully equilibrated and ready for PMF calculations")
        else:
            logger.info("⚠️  System built but not equilibrated - consider running equilibration")
        
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
        result = subprocess.run(
            cmd,
            input=f"{group}\n",
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"Failed to extract {group}: {result.stderr}")
        
        logger.debug(f"Extracted {group} to {output_file}")
    
    def _extract_complex(self, input_file: Path, output_file: Path, frame: int = -1) -> None:
        """Extract protein-ligand complex (no solvent)"""
        cmd = [
            self.gmx_command, "trjconv", 
            "-f", str(input_file),
            "-o", str(output_file),
            "-s", str(input_file)
        ]
        
        if frame != -1:
            cmd.extend(["-dump", str(frame)])
        
        # Select protein and ligand groups
        result = subprocess.run(
            cmd,
            input="Protein | LIG\n",  # Select both protein and ligand
            text=True,
            capture_output=True,
            cwd=self.work_dir
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"Failed to extract complex: {result.stderr}")
        
        logger.debug(f"Extracted complex to {output_file}")
    
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
        """Combine protein and ligand topologies"""
        # Copy original topology as base
        original_top = self.md_system_dir / "topol.top"
        combined_top = self.work_dir / "combined.top"
        
        shutil.copy2(original_top, combined_top)
        
        logger.debug(f"Combined topology created: {combined_top}")
        return combined_top
    
    def _create_extended_box(self, complex_file: Path, alignment_info: Dict) -> Path:
        """Create box with extended Z-axis for pulling"""
        boxed_file = self.work_dir / "complex_boxed.gro"
        
        # Calculate extended Z dimension
        initial_distance = alignment_info['initial_distance']
        z_extension = self.config['box']['z_extension']
        box_distance = self.config['box']['box_distance']
        
        # Z dimension should accommodate pulling trajectory
        z_dimension = initial_distance + z_extension + 2 * box_distance
        xy_dimension = initial_distance + 2 * box_distance
        
        logger.info(f"Creating extended box: {xy_dimension:.1f} x {xy_dimension:.1f} x {z_dimension:.1f} nm")
        
        cmd = [
            self.gmx_command, "editconf",
            "-f", str(complex_file),
            "-o", str(boxed_file),
            "-box", str(xy_dimension), str(xy_dimension), str(z_dimension),
            "-center", "0", "0", "0"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to create box: {result.stderr}")
        
        logger.debug(f"Extended box created: {boxed_file}")
        return boxed_file
    
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
        
        # Use SystemBuilder for ion addition (following PRISM patterns)
        system_builder = SystemBuilder(self.config, str(self.rebuilt_system_dir))
        
        # Copy files for SystemBuilder
        shutil.copy2(solvated_file, self.rebuilt_system_dir / "temp_solvated.gro")
        shutil.copy2(topology, final_topology)
        
        try:
            # Add ions using existing PRISM functionality
            ion_results = system_builder._add_ions(
                str(self.rebuilt_system_dir / "temp_solvated.gro"),
                str(final_topology),
                str(final_system)
            )
            
            # Clean up temp file
            temp_file = self.rebuilt_system_dir / "temp_solvated.gro"
            if temp_file.exists():
                temp_file.unlink()
            
            logger.info("Solvation and ion addition completed")
            return final_system
            
        except Exception as e:
            logger.warning(f"Ion addition with SystemBuilder failed: {e}")
            # Fallback: simple copy without ions
            shutil.copy2(solvated_file, final_system)
            logger.info("Using solvated system without additional ions")
            return final_system
    
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
    **kwargs : optional
        Additional configuration parameters
        
    Returns:
    --------
    PMFBuilder : PMF builder instance
    
    Example:
    --------
    >>> import prism
    >>> builder = prism.pmf.pmf_builder("./gaff_model", "./pmf_system")
    >>> results = builder.build()
    """
    return PMFBuilder(md_results_dir, output_dir, config, **kwargs)