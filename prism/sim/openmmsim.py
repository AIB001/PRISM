#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenMM simulator for PRISM MD simulations with enhanced CUDA support and progress bars
"""

import os
import sys
import subprocess
import tempfile
import re
import time
from pathlib import Path
from datetime import datetime, timedelta

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    import numpy as np
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    print("Warning: OpenMM not available. RECOMMENDED INSTALLATION:")
    print("  mamba install -c conda-forge openmm")
    print("  # OR: conda install -c conda-forge openmm")

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    print("Info: tqdm not available for progress bars. Install with: pip install tqdm")


class ProgressReporter:
    """Custom reporter for displaying simulation progress with tqdm"""
    
    def __init__(self, total_steps, update_interval=1000, desc="Simulation"):
        self.total_steps = total_steps
        self.update_interval = update_interval
        self.desc = desc
        self.start_time = None
        self.pbar = None
        
    def __enter__(self):
        if TQDM_AVAILABLE:
            self.pbar = tqdm(total=self.total_steps, desc=self.desc, 
                           unit='steps', unit_scale=True)
        self.start_time = time.time()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.pbar:
            self.pbar.close()
    
    def update(self, current_step):
        if self.pbar:
            increment = current_step - self.pbar.n
            if increment > 0:
                self.pbar.update(increment)
        else:
            # Fallback to simple text progress
            if current_step % (self.total_steps // 20) == 0 or current_step == self.total_steps:
                elapsed = time.time() - self.start_time
                progress = current_step / self.total_steps * 100
                eta = elapsed / (current_step / self.total_steps) - elapsed if current_step > 0 else 0
                print(f"{self.desc}: {progress:.1f}% ({current_step}/{self.total_steps} steps) "
                      f"ETA: {timedelta(seconds=int(eta))}")


class TqdmReporter:
    """OpenMM reporter that updates tqdm progress bar"""
    
    def __init__(self, progress_reporter, reportInterval):
        self.progress_reporter = progress_reporter
        self.reportInterval = reportInterval
        
    def describeNextReport(self, simulation):
        steps = self.reportInterval - simulation.currentStep % self.reportInterval
        return (steps, False, False, False, False, None)
    
    def report(self, simulation, state):
        self.progress_reporter.update(simulation.currentStep)


class OpenMMSimulator:
    """
    OpenMM simulator class for running MD simulations.
    
    This class reads GROMACS topology files and runs equivalent simulations
    using OpenMM with the same parameters as in the MDP files.
    """
    
    def __init__(self, gmx_dir, system_files, ff_dir):
        """
        Initialize OpenMM simulator.
        
        Parameters:
        -----------
        gmx_dir : str
            Path to GMX_PROLIG_MD directory
        system_files : dict
            Dictionary of system file paths
        ff_dir : str
            Path to force field directory
        """
        if not OPENMM_AVAILABLE:
            raise ImportError("OpenMM is not installed. Please install it to use this simulator.")
        
        self.gmx_dir = gmx_dir
        self.system_files = system_files
        self.ff_dir = ff_dir
        
        # Find GROMACS installation
        self.gmx_data_dir = self._find_gromacs_data_dir()
        
        # Parse MDP files to get simulation parameters
        self.mdp_params = self._parse_mdp_files()
        
        # Check CUDA availability at initialization
        self.cuda_info = self._check_cuda_availability()
        
        print("OpenMM simulator initialized")
        print(f"  Platforms available: {self._get_platform_info()}")
        if self.cuda_info['available']:
            print(f"  CUDA: Available (Driver: {self.cuda_info['driver_version']}, "
                  f"Devices: {self.cuda_info['device_count']})")
        else:
            print(f"  CUDA: Not available ({self.cuda_info['reason']})")
        if self.gmx_data_dir:
            print(f"  GROMACS data directory: {self.gmx_data_dir}")
    
    def _check_cuda_availability(self):
        """Check CUDA availability and version"""
        cuda_info = {
            'available': False,
            'driver_version': None,
            'device_count': 0,
            'devices': [],
            'reason': 'Not checked'
        }
        
        try:
            # Check if CUDA platform is available in OpenMM
            for i in range(mm.Platform.getNumPlatforms()):
                platform = mm.Platform.getPlatform(i)
                if platform.getName() == 'CUDA':
                    cuda_info['available'] = True
                    cuda_info['reason'] = 'CUDA platform found'
                    
                    # Try to get CUDA device information
                    try:
                        # Get device count
                        cuda_info['device_count'] = int(platform.getPropertyDefaultValue('DeviceIndex')) + 1
                    except:
                        cuda_info['device_count'] = 1
                    
                    # Try to get driver version using nvidia-smi
                    try:
                        result = subprocess.run(['nvidia-smi', '--query-gpu=driver_version', 
                                               '--format=csv,noheader'], 
                                              capture_output=True, text=True)
                        if result.returncode == 0:
                            cuda_info['driver_version'] = result.stdout.strip().split('\n')[0]
                            
                        # Get device names
                        result = subprocess.run(['nvidia-smi', '--query-gpu=name', 
                                               '--format=csv,noheader'], 
                                              capture_output=True, text=True)
                        if result.returncode == 0:
                            cuda_info['devices'] = result.stdout.strip().split('\n')
                    except:
                        pass
                    
                    break
            
            if not cuda_info['available']:
                cuda_info['reason'] = 'CUDA platform not found in OpenMM'
                
        except Exception as e:
            cuda_info['reason'] = str(e)
        
        return cuda_info
    
    def _diagnose_topology_includes(self, top_file):
        """Diagnose include statements in topology file"""
        print("\nDiagnosing topology file includes:")
        with open(top_file, 'r') as f:
            lines = f.readlines()
        
        includes = []
        for i, line in enumerate(lines, 1):
            if line.strip().startswith('#include'):
                includes.append((i, line.strip()))
        
        if includes:
            print(f"Found {len(includes)} include statements:")
            for line_num, include in includes:
                print(f"  Line {line_num}: {include}")
                
                # Check for syntax errors
                if include.count('"') != 2:
                    print(f"    -> WARNING: Syntax error - unmatched quotes!")
                
                # Extract filename
                match = re.search(r'#include\s+"([^"]+)"', include)
                if match:
                    filename = match.group(1)
                    # Check if it's a force field file
                    if '.ff/' in filename:
                        print(f"    -> Force field include: {filename}")
                        if self.gmx_data_dir:
                            full_path = os.path.join(self.gmx_data_dir, filename)
                            exists = "EXISTS" if os.path.exists(full_path) else "NOT FOUND"
                            print(f"    -> Full path: {full_path} [{exists}]")
                else:
                    # Try to extract filename even with syntax errors
                    match2 = re.search(r'#include\s+"([^"]+)', include)
                    if match2:
                        filename = match2.group(1)
                        print(f"    -> Possibly malformed include: {filename}")
        else:
            print("No include statements found in topology file")
    
    def _find_gromacs_data_dir(self):
        """Find GROMACS data directory containing force field files"""
        try:
            # Find GROMACS binary
            result = subprocess.run(['which', 'gmx'], capture_output=True, text=True)
            if result.returncode == 0:
                gmx_bin = result.stdout.strip()
                print(f"Found GROMACS binary: {gmx_bin}")
                
                # GROMACS data is typically in ../share/gromacs/top relative to bin
                gmx_bin_dir = os.path.dirname(gmx_bin)
                
                # Check if it's a conda installation
                if 'conda' in gmx_bin or 'miniconda' in gmx_bin or 'anaconda' in gmx_bin:
                    print("Detected conda GROMACS installation")
                
                potential_dirs = [
                    os.path.join(os.path.dirname(gmx_bin_dir), 'share', 'gromacs', 'top'),
                    os.path.join(os.path.dirname(gmx_bin_dir), 'share', 'top'),
                    # Conda environments might have it here
                    os.path.join(os.path.dirname(os.path.dirname(gmx_bin)), 'share', 'gromacs', 'top'),
                    # System-wide installations
                    '/usr/share/gromacs/top',
                    '/usr/local/share/gromacs/top',
                    # Alternative conda location
                    os.path.join(os.path.dirname(gmx_bin_dir), 'share', 'gromacs'),
                ]
                
                for dir_path in potential_dirs:
                    if os.path.exists(dir_path):
                        # Verify it contains force field files
                        ff_files = [f for f in os.listdir(dir_path) if f.endswith('.ff')]
                        if ff_files:
                            print(f"Found GROMACS data directory: {dir_path}")
                            print(f"Available force fields: {', '.join(ff_files[:5])}...")
                            return dir_path
        except Exception as e:
            print(f"Warning: Could not find GROMACS installation: {e}")
        
        print("Warning: GROMACS data directory not found")
        return None
    
    def _create_modified_topology(self, original_top):
        """Create a modified topology file with absolute paths for ALL includes"""
        # Read original topology
        with open(original_top, 'r') as f:
            lines = f.readlines()
        
        # Get the directory of the topology file
        top_dir = os.path.dirname(original_top)
        
        print("\nConverting include paths to absolute paths:")
        
        modified_lines = []
        for line in lines:
            # Check if this is an include line
            if line.strip().startswith('#include'):
                # First, fix any syntax errors (missing quotes)
                if line.count('"') == 1:
                    print(f"  Fixing syntax error in: {line.strip()}")
                    line = line.rstrip() + '"\n'
                
                # Now process the include
                match = re.search(r'#include\s+"([^"]+)"', line)
                if match:
                    include_file = match.group(1)
                    new_path = self._resolve_include_path(include_file, top_dir)
                    if new_path != include_file:
                        line = f'#include "{new_path}"\n'
            
            modified_lines.append(line)
        
        # Join lines back together
        modified_content = ''.join(modified_lines)
        
        # Only create a new file if we actually modified something
        if modified_lines != lines:
            # Create temporary file with modified topology
            with tempfile.NamedTemporaryFile(mode='w', suffix='.top', delete=False) as tmp:
                tmp.write(modified_content)
                print(f"\nCreated modified topology: {tmp.name}")
                return tmp.name
        else:
            print("\nNo modifications needed to topology file")
            return original_top
    
    def _resolve_include_path(self, include_file, top_dir):
        """Resolve include file path to absolute path"""
        # Check if it's already an absolute path
        if os.path.isabs(include_file):
            return include_file
        
        # Handle force field includes (contain .ff/)
        if '.ff/' in include_file:
            if self.gmx_data_dir:
                absolute_path = os.path.join(self.gmx_data_dir, include_file)
                if os.path.exists(absolute_path):
                    print(f"  Force field: {include_file} -> {absolute_path}")
                    return absolute_path
                else:
                    # Try alternative locations
                    for base_dir in ['/usr/share/gromacs/top', '/usr/local/share/gromacs/top']:
                        alt_path = os.path.join(base_dir, include_file)
                        if os.path.exists(alt_path):
                            print(f"  Force field: {include_file} -> {alt_path}")
                            return alt_path
        
        # Handle relative paths (../)
        if include_file.startswith('../'):
            # Resolve relative to topology file directory
            absolute_path = os.path.abspath(os.path.join(top_dir, include_file))
            if os.path.exists(absolute_path):
                print(f"  Relative: {include_file} -> {absolute_path}")
                return absolute_path
            else:
                print(f"  ERROR: Could not resolve relative path: {include_file}")
                print(f"    Tried: {absolute_path}")
        
        # Handle files in the same directory as topology
        else:
            # First check in topology directory
            absolute_path = os.path.join(top_dir, include_file)
            if os.path.exists(absolute_path):
                print(f"  Local: {include_file} -> {absolute_path}")
                return absolute_path
            
            # Check in force field directory
            if self.ff_dir:
                ff_path = os.path.join(self.ff_dir, include_file)
                if os.path.exists(ff_path):
                    print(f"  Force field dir: {include_file} -> {ff_path}")
                    return ff_path
            
            # Also check parent directory (for cases like posre.itp)
            parent_dir = os.path.dirname(top_dir)
            parent_path = os.path.join(parent_dir, include_file)
            if os.path.exists(parent_path):
                print(f"  Parent dir: {include_file} -> {parent_path}")
                return parent_path
        
        # If we couldn't resolve it, print a warning but keep original
        print(f"  WARNING: Could not resolve: {include_file}")
        return include_file
    
    def _get_platform_info(self):
        """Get information about available OpenMM platforms"""
        platforms = []
        for i in range(mm.Platform.getNumPlatforms()):
            platform = mm.Platform.getPlatform(i)
            platforms.append(platform.getName())
        return ", ".join(platforms)
    
    def _verify_modified_topology(self, top_file):
        """Verify that all includes in modified topology are absolute paths"""
        with open(top_file, 'r') as f:
            lines = f.readlines()
        
        print("  Checking includes in modified topology:")
        for i, line in enumerate(lines):
            if line.strip().startswith('#include'):
                match = re.search(r'#include\s+"([^"]+)"', line)
                if match:
                    path = match.group(1)
                    if os.path.isabs(path):
                        exists = "EXISTS" if os.path.exists(path) else "NOT FOUND"
                        print(f"    ✓ {path} [{exists}]")
                    else:
                        print(f"    ✗ Still relative: {path}")
    
    def _parse_mdp_files(self):
        """Parse MDP files to extract simulation parameters"""
        params = {}
        
        # Parse each MDP file
        for stage in ['em', 'nvt', 'npt', 'md']:
            mdp_file = self.system_files.get(f'{stage}_mdp')
            if mdp_file and os.path.exists(mdp_file):
                params[stage] = self._parse_mdp(mdp_file)
        
        return params
    
    def _parse_mdp(self, mdp_file):
        """Parse a single MDP file"""
        params = {}
        
        with open(mdp_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith(';'):
                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        # Remove inline comments
                        if ';' in value:
                            value = value.split(';')[0].strip()
                        params[key] = value
        
        return params
    
    def run(self, stages=None, platform='auto', device_index=0, num_threads=10, gmx_data_dir=None, **kwargs):
        """
        Run OpenMM MD simulation.
        
        Parameters:
        -----------
        stages : list, optional
            List of stages to run ['em', 'nvt', 'npt', 'prod']
            If None, runs all stages
        platform : str
            OpenMM platform to use ('auto', 'CUDA', 'OpenCL', 'CPU')
            'auto' will try CUDA first, then fall back to CPU with multiple threads
        device_index : int
            GPU device index (default: 0)
        num_threads : int
            Number of CPU threads to use (default: 10)
        gmx_data_dir : str, optional
            Manual path to GROMACS data directory (e.g., '/usr/share/gromacs/top')
        **kwargs : dict
            Additional parameters
            
        Returns:
        --------
        dict
            Dictionary with paths to output files
        """
        # Override auto-detected GROMACS directory if provided
        if gmx_data_dir:
            print(f"Using manually specified GROMACS data directory: {gmx_data_dir}")
            self.gmx_data_dir = gmx_data_dir
        
        # Default stages
        if stages is None:
            stages = ['em', 'nvt', 'npt', 'prod']
        
        # Setup platform with automatic fallback
        print(f"\nAttempting to use {platform} platform...")
        platform_obj = self._setup_platform(platform, device_index)
        
        # If CUDA failed and we fell back, inform the user
        if platform == 'CUDA' and platform_obj[0].getName() != 'CUDA':
            print(f"\n{'='*60}")
            print(f"Note: Using {platform_obj[0].getName()} platform instead of CUDA")
            print("To fix CUDA compatibility issues:")
            print("1. Check CUDA driver version: nvidia-smi")
            print("2. Reinstall OpenMM with matching CUDA version:")
            print("   mamba install -c conda-forge openmm cudatoolkit=11.7")
            print("   # OR: conda install -c conda-forge openmm cudatoolkit=11.7")
            print("3. Or explicitly use CPU: sim.run(engine='openmm', platform='CPU')")
            print(f"{'='*60}\n")
        
        # Load system from GROMACS files
        print("Loading GROMACS system...")
        gro_file = self.system_files['gro']
        top_file = self.system_files['top']
        
        # Diagnose topology includes
        self._diagnose_topology_includes(top_file)
        
        # Create modified topology with absolute paths
        modified_top = self._create_modified_topology(top_file)
        if modified_top != top_file:
            print("\nTopology file has been modified with absolute paths")
            # Verify the modified topology
            print("\nVerifying modified topology:")
            self._verify_modified_topology(modified_top)
        
        try:
            # Use GromacsGroFile and GromacsTopFile
            gro = app.GromacsGroFile(gro_file)
            
            # Since we've converted all paths to absolute, we don't need includeDir
            # But OpenMM requires it, so we'll use the topology directory
            print(f"\nLoading topology with all absolute paths...")
            
            try:
                top = app.GromacsTopFile(
                    modified_top,
                    periodicBoxVectors=gro.getPeriodicBoxVectors()
                )
                print("Successfully loaded topology file!")
            except Exception as e:
                print("\n" + "="*60)
                print("ERROR: Could not load GROMACS topology file in OpenMM")
                print("="*60)
                print(f"\nError: {str(e)}")
                
                # Try to identify which file is missing
                if "Could not locate #include file:" in str(e):
                    missing_file = str(e).split("Could not locate #include file:")[-1].strip()
                    print(f"\nMissing file: {missing_file}")
                    print("\nThis usually means the path conversion didn't work properly.")
                
                print("\nPossible solutions:")
                print("1. Use GROMACS engine instead (recommended):")
                print("   sim.run(engine='gmx')")
                print("\n2. Check that all required files exist:")
                print(f"   - Topology: {top_file}")
                print(f"   - Ligand files: {self.ff_dir}")
                print(f"   - GROMACS data: {self.gmx_data_dir}")
                print("\n3. Manually copy all required files to one directory")
                print("="*60)
                raise RuntimeError(
                    f"Could not load topology file in OpenMM. "
                    f"Consider using GROMACS engine instead: sim.run(engine='gmx')"
                )
            
            # Create system
            system = top.createSystem(
                nonbondedMethod=app.PME,
                nonbondedCutoff=1.0*unit.nanometer,
                constraints=app.HBonds
            )
            
            # Run stages
            outputs = {}
            positions = gro.positions
            velocities = None
            box_vectors = gro.getPeriodicBoxVectors()
            
            for stage in stages:
                print(f"\n{'='*60}")
                print(f" Running {stage.upper()} stage")
                print(f"{'='*60}")
                
                if stage == 'em':
                    positions, box_vectors = self._run_minimization(
                        system, top.topology, positions, box_vectors, platform_obj
                    )
                    outputs['em'] = {'gro': os.path.join(self.gmx_dir, 'em_openmm.gro')}
                    
                elif stage == 'nvt':
                    positions, velocities, box_vectors = self._run_nvt(
                        system, top.topology, positions, box_vectors, platform_obj
                    )
                    outputs['nvt'] = {'gro': os.path.join(self.gmx_dir, 'nvt_openmm.gro')}
                    
                elif stage == 'npt':
                    positions, velocities, box_vectors = self._run_npt(
                        system, top.topology, positions, velocities, box_vectors, platform_obj
                    )
                    outputs['npt'] = {'gro': os.path.join(self.gmx_dir, 'npt_openmm.gro')}
                    
                elif stage == 'prod':
                    self._run_production(
                        system, top.topology, positions, velocities, box_vectors, platform_obj
                    )
                    outputs['prod'] = {
                        'dcd': os.path.join(self.gmx_dir, 'prod_openmm.dcd'),
                        'log': os.path.join(self.gmx_dir, 'prod_openmm.log'),
                        'gro': os.path.join(self.gmx_dir, 'prod_openmm.gro')
                    }
            
            return outputs
            
        finally:
            # Clean up temporary topology file
            if modified_top != top_file and os.path.exists(modified_top):
                os.unlink(modified_top)
    
    def _setup_platform(self, platform_name, device_index):
        """Setup OpenMM platform with automatic fallback and better diagnostics"""
        try:
            platform = mm.Platform.getPlatformByName(platform_name)
            
            if platform_name in ['CUDA', 'OpenCL']:
                properties = {'DeviceIndex': str(device_index)}
                if platform_name == 'CUDA':
                    properties['Precision'] = 'mixed'
                    # Add more CUDA-specific optimizations
                    properties['UseCpuPme'] = 'false'
                
                # Test if the platform actually works
                try:
                    # Create a minimal test system
                    test_system = mm.System()
                    test_system.addParticle(1.0)
                    test_integrator = mm.VerletIntegrator(1.0)
                    test_context = mm.Context(test_system, test_integrator, platform, properties)
                    
                    # If CUDA, get device info
                    if platform_name == 'CUDA' and self.cuda_info['devices']:
                        device_name = self.cuda_info['devices'][device_index] if device_index < len(self.cuda_info['devices']) else 'Unknown'
                        print(f"Successfully initialized {platform_name} platform on device {device_index}: {device_name}")
                    else:
                        print(f"Successfully initialized {platform_name} platform")
                    
                    del test_context
                    del test_integrator
                    del test_system
                    return platform, properties
                    
                except Exception as e:
                    error_msg = str(e)
                    if "CUDA_ERROR_UNSUPPORTED_PTX_VERSION" in error_msg:
                        print(f"\n⚠️  CUDA version incompatibility detected!")
                        print("Your CUDA driver is too old for this version of OpenMM.")
                        print(f"Current driver version: {self.cuda_info.get('driver_version', 'Unknown')}")
                        print("Falling back to OpenCL platform...")
                        
                        # Try OpenCL
                        if platform_name == 'CUDA':
                            return self._setup_platform('OpenCL', device_index)
                    elif "CUDA_ERROR" in error_msg:
                        print(f"\n⚠️  CUDA error: {error_msg}")
                        print("Falling back to OpenCL platform...")
                        if platform_name == 'CUDA':
                            return self._setup_platform('OpenCL', device_index)
                    else:
                        raise
            else:
                print(f"Using {platform_name} platform (no GPU acceleration)")
                return platform, {}
                
        except Exception as e:
            print(f"⚠️  {platform_name} platform not available: {e}")
            if platform_name != 'CPU':
                print("Falling back to CPU platform...")
                return mm.Platform.getPlatformByName('CPU'), {}
            else:
                raise
    
    def _write_gro_file(self, filename, topology, positions, box_vectors, title="Generated by OpenMM"):
        """Write a GROMACS .gro file"""
        try:
            with open(filename, 'w') as f:
                # Write title
                f.write(f"{title}\n")
                
                # Write number of atoms
                f.write(f"{topology.getNumAtoms():5d}\n")
                
                # Write atom data
                for i, atom in enumerate(topology.atoms()):
                    pos = positions[i].value_in_unit(unit.nanometer)
                    
                    # Format: residue number (5), residue name (5), atom name (5), 
                    # atom number (5), x (8.3f), y (8.3f), z (8.3f)
                    f.write(f"{atom.residue.index + 1:5d}"
                           f"{atom.residue.name[:5]:>5s}"
                           f"{atom.name[:5]:>5s}"
                           f"{i + 1:5d}"
                           f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}\n")
                
                # Write box vectors
                if box_vectors is not None:
                    box = box_vectors.value_in_unit(unit.nanometer)
                    if len(box) == 3:
                        # Rectangular box
                        f.write(f"{box[0][0]:10.5f}{box[1][1]:10.5f}{box[2][2]:10.5f}\n")
                    else:
                        # Triclinic box - write full matrix
                        f.write(f"{box[0][0]:10.5f}{box[1][1]:10.5f}{box[2][2]:10.5f}"
                               f"{box[0][1]:10.5f}{box[0][2]:10.5f}{box[1][0]:10.5f}"
                               f"{box[1][2]:10.5f}{box[2][0]:10.5f}{box[2][1]:10.5f}\n")
        except Exception as e:
            print(f"Warning: Failed to write GRO file: {e}")
            # Fall back to PDB format
            pdb_filename = filename.replace('.gro', '.pdb')
            print(f"Writing PDB file instead: {pdb_filename}")
            with open(pdb_filename, 'w') as f:
                app.PDBFile.writeFile(topology, positions, f)
    
    def _run_minimization(self, system, topology, positions, box_vectors, platform_info):
        """Run energy minimization with progress bar"""
        platform, properties = platform_info
        
        # Get parameters from MDP
        em_params = self.mdp_params.get('em', {})
        max_iterations = int(em_params.get('nsteps', 10000))
        tolerance = float(em_params.get('emtol', 200.0)) * unit.kilojoules_per_mole/unit.nanometer
        
        print(f"Running energy minimization (max {max_iterations} steps)...")
        print(f"Tolerance: {tolerance}")
        
        # Create integrator
        integrator = mm.VerletIntegrator(1.0*unit.femtoseconds)
        
        # Create simulation
        simulation = app.Simulation(topology, system, integrator, platform, properties)
        simulation.context.setPositions(positions)
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        
        # Get initial energy
        state = simulation.context.getState(getEnergy=True)
        initial_energy = state.getPotentialEnergy()
        print(f"Initial potential energy: {initial_energy}")
        
        # Minimize with progress tracking
        with ProgressReporter(max_iterations, desc="Energy Minimization") as progress:
            # We'll do minimization in chunks to show progress
            chunk_size = max(100, max_iterations // 100)
            remaining_iterations = max_iterations
            
            while remaining_iterations > 0:
                current_chunk = min(chunk_size, remaining_iterations)
                simulation.minimizeEnergy(tolerance=tolerance, maxIterations=current_chunk)
                
                # Update progress
                progress.update(max_iterations - remaining_iterations + current_chunk)
                
                # Check if converged
                state = simulation.context.getState(getEnergy=True)
                current_energy = state.getPotentialEnergy()
                
                # Simple convergence check
                if abs((current_energy - initial_energy) / initial_energy) < 1e-6:
                    print(f"\nConverged at step {max_iterations - remaining_iterations + current_chunk}")
                    break
                    
                remaining_iterations -= current_chunk
        
        # Get minimized positions
        state = simulation.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True)
        positions = state.getPositions()
        box_vectors = state.getPeriodicBoxVectors()
        final_energy = state.getPotentialEnergy()
        
        print(f"Final potential energy: {final_energy}")
        print(f"Energy change: {final_energy - initial_energy}")
        
        # Save minimized structure
        output_file = os.path.join(self.gmx_dir, 'em_openmm.gro')
        self._write_gro_file(output_file, topology, positions, box_vectors, 
                            title="Energy minimized structure")
        print(f"Minimized structure saved to: {output_file}")
        
        return positions, box_vectors
    
    def _run_nvt(self, system, topology, positions, box_vectors, platform_info):
        """Run NVT equilibration with progress bar"""
        platform, properties = platform_info
        
        # Get parameters from MDP
        nvt_params = self.mdp_params.get('nvt', {})
        dt = float(nvt_params.get('dt', 0.002)) * unit.picoseconds
        nsteps = int(nvt_params.get('nsteps', 250000))
        temperature = float(nvt_params.get('ref_t', '310').split()[0]) * unit.kelvin
        
        total_time = (nsteps * dt).in_units_of(unit.nanoseconds)
        print(f"Running NVT equilibration")
        print(f"  Steps: {nsteps:,}")
        print(f"  Time: {total_time}")
        print(f"  Temperature: {temperature}")
        
        # Create integrator with thermostat
        integrator = mm.LangevinIntegrator(
            temperature,
            1.0/unit.picosecond,
            dt
        )
        
        # Create simulation
        simulation = app.Simulation(topology, system, integrator, platform, properties)
        simulation.context.setPositions(positions)
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        simulation.context.setVelocitiesToTemperature(temperature)
        
        # Add reporters
        log_file = os.path.join(self.gmx_dir, 'nvt_openmm.log')
        simulation.reporters.append(
            app.StateDataReporter(
                log_file, 5000, step=True, time=True, 
                temperature=True, progress=False,  # Disable built-in progress
                remainingTime=False, speed=True,
                totalSteps=nsteps, separator='\t'
            )
        )
        
        # Add progress bar reporter
        with ProgressReporter(nsteps, desc="NVT Equilibration") as progress:
            simulation.reporters.append(TqdmReporter(progress, 1000))
            
            # Run simulation
            simulation.step(nsteps)
        
        # Get final state
        state = simulation.context.getState(
            getPositions=True, getVelocities=True, enforcePeriodicBox=True
        )
        positions = state.getPositions()
        velocities = state.getVelocities()
        box_vectors = state.getPeriodicBoxVectors()
        
        # Save final structure
        output_file = os.path.join(self.gmx_dir, 'nvt_openmm.gro')
        self._write_gro_file(output_file, topology, positions, box_vectors,
                            title="NVT equilibrated structure")
        print(f"NVT equilibrated structure saved to: {output_file}")
        
        return positions, velocities, box_vectors
    
    def _run_npt(self, system, topology, positions, velocities, box_vectors, platform_info):
        """Run NPT equilibration with progress bar"""
        platform, properties = platform_info
        
        # Get parameters from MDP
        npt_params = self.mdp_params.get('npt', {})
        dt = float(npt_params.get('dt', 0.002)) * unit.picoseconds
        nsteps = int(npt_params.get('nsteps', 250000))
        temperature = float(npt_params.get('ref_t', '310').split()[0]) * unit.kelvin
        pressure = float(npt_params.get('ref_p', 1.0)) * unit.bar
        
        total_time = (nsteps * dt).in_units_of(unit.nanoseconds)
        print(f"Running NPT equilibration")
        print(f"  Steps: {nsteps:,}")
        print(f"  Time: {total_time}")
        print(f"  Temperature: {temperature}")
        print(f"  Pressure: {pressure}")
        
        # Add barostat
        barostat = mm.MonteCarloBarostat(pressure, temperature, 25)
        system.addForce(barostat)
        
        # Create integrator
        integrator = mm.LangevinIntegrator(
            temperature,
            1.0/unit.picosecond,
            dt
        )
        
        # Create simulation
        simulation = app.Simulation(topology, system, integrator, platform, properties)
        simulation.context.setPositions(positions)
        if velocities is not None:
            simulation.context.setVelocities(velocities)
        else:
            simulation.context.setVelocitiesToTemperature(temperature)
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        
        # Add reporters
        log_file = os.path.join(self.gmx_dir, 'npt_openmm.log')
        simulation.reporters.append(
            app.StateDataReporter(
                log_file, 5000, step=True, time=True,
                temperature=True, volume=True, density=True,
                progress=False, remainingTime=False, speed=True,
                totalSteps=nsteps, separator='\t'
            )
        )
        
        # Add progress bar reporter
        with ProgressReporter(nsteps, desc="NPT Equilibration") as progress:
            simulation.reporters.append(TqdmReporter(progress, 1000))
            
            # Run simulation
            simulation.step(nsteps)
        
        # Get final state
        state = simulation.context.getState(
            getPositions=True, getVelocities=True, enforcePeriodicBox=True
        )
        positions = state.getPositions()
        velocities = state.getVelocities()
        box_vectors = state.getPeriodicBoxVectors()
        
        # Save final structure
        output_file = os.path.join(self.gmx_dir, 'npt_openmm.gro')
        self._write_gro_file(output_file, topology, positions, box_vectors,
                            title="NPT equilibrated structure")
        print(f"NPT equilibrated structure saved to: {output_file}")
        
        return positions, velocities, box_vectors
    
    def _run_production(self, system, topology, positions, velocities, box_vectors, platform_info):
        """Run production MD with progress bar"""
        platform, properties = platform_info
        
        # Get parameters from MDP
        md_params = self.mdp_params.get('md', {})
        dt = float(md_params.get('dt', 0.002)) * unit.picoseconds
        nsteps = int(md_params.get('nsteps', 250000000))
        temperature = float(md_params.get('ref_t', '310').split()[0]) * unit.kelvin
        pressure = float(md_params.get('ref_p', 1.0)) * unit.bar
        
        # Output frequency
        nstxtcout = int(md_params.get('nstxtcout', 250000))
        
        total_time = (nsteps * dt).in_units_of(unit.nanoseconds)
        print(f"Running production MD")
        print(f"  Steps: {nsteps:,}")
        print(f"  Time: {total_time}")
        print(f"  Temperature: {temperature}")
        print(f"  Pressure: {pressure}")
        print(f"  Frame output interval: {nstxtcout} steps")
        
        # Add barostat
        barostat = mm.MonteCarloBarostat(pressure, temperature, 25)
        system.addForce(barostat)
        
        # Create integrator
        integrator = mm.LangevinIntegrator(
            temperature,
            1.0/unit.picosecond,
            dt
        )
        
        # Create simulation
        simulation = app.Simulation(topology, system, integrator, platform, properties)
        simulation.context.setPositions(positions)
        if velocities is not None:
            simulation.context.setVelocities(velocities)
        else:
            simulation.context.setVelocitiesToTemperature(temperature)
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        
        # Add reporters
        log_file = os.path.join(self.gmx_dir, 'prod_openmm.log')
        dcd_file = os.path.join(self.gmx_dir, 'prod_openmm.dcd')
        
        simulation.reporters.append(
            app.DCDReporter(dcd_file, nstxtcout)
        )
        simulation.reporters.append(
            app.StateDataReporter(
                log_file, 5000, step=True, time=True,
                temperature=True, volume=True, density=True,
                potentialEnergy=True, kineticEnergy=True,
                progress=False, remainingTime=False, speed=True,
                totalSteps=nsteps, separator='\t'
            )
        )
        
        # Add progress bar reporter with checkpoint saving
        checkpoint_interval = 1000000  # Save checkpoint every 1M steps
        checkpoint_file = os.path.join(self.gmx_dir, 'prod_openmm.chk')
        
        with ProgressReporter(nsteps, desc="Production MD") as progress:
            simulation.reporters.append(TqdmReporter(progress, 1000))
            
            # Run simulation with periodic checkpointing
            steps_completed = 0
            while steps_completed < nsteps:
                steps_to_run = min(checkpoint_interval, nsteps - steps_completed)
                simulation.step(steps_to_run)
                steps_completed += steps_to_run
                
                # Save checkpoint
                if steps_completed % checkpoint_interval == 0:
                    simulation.saveCheckpoint(checkpoint_file)
                    if not TQDM_AVAILABLE:
                        print(f"Checkpoint saved at step {steps_completed}")
        
        # Save final structure
        state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        final_positions = state.getPositions()
        final_box_vectors = state.getPeriodicBoxVectors()
        
        output_gro = os.path.join(self.gmx_dir, 'prod_openmm.gro')
        self._write_gro_file(output_gro, topology, final_positions, final_box_vectors,
                            title="Production MD final structure")
        
        print(f"\n✅ Production MD completed!")
        print(f"  Trajectory: {dcd_file}")
        print(f"  Log file: {log_file}")
        print(f"  Final structure: {output_gro}")
        print(f"  Checkpoint: {checkpoint_file}")