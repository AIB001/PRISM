#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utility functions for PRISM simulation module
"""

import os
from pathlib import Path


def validate_gmx_directory(gmx_dir):
    """
    Validate that a directory contains necessary GROMACS files.
    
    Parameters:
    -----------
    gmx_dir : str
        Path to GMX_PROLIG_MD directory
        
    Returns:
    --------
    bool
        True if directory is valid
        
    Raises:
    -------
    FileNotFoundError
        If directory or required files are missing
    """
    if not os.path.exists(gmx_dir):
        raise FileNotFoundError(f"Directory not found: {gmx_dir}")
    
    # Check for required files
    required_files = ['solv_ions.gro', 'topol.top']
    missing = []
    
    for file in required_files:
        if not os.path.exists(os.path.join(gmx_dir, file)):
            missing.append(file)
    
    if missing:
        raise FileNotFoundError(
            f"Missing required files in {gmx_dir}: {', '.join(missing)}"
        )
    
    return True


def find_forcefield_dir(output_dir):
    """
    Find the force field directory (LIG.amb2gmx or LIG.openff2gmx).
    
    Parameters:
    -----------
    output_dir : str
        Output directory containing GMX_PROLIG_MD, mdps, and force field directory
        
    Returns:
    --------
    str
        Path to force field directory
        
    Raises:
    -------
    FileNotFoundError
        If no force field directory is found
    """
    possible_names = ['LIG.amb2gmx', 'LIG.openff2gmx', 'LIG.gaff2gmx']
    
    # Check in the output directory (same level as GMX_PROLIG_MD)
    for name in possible_names:
        ff_dir = os.path.join(output_dir, name)
        if os.path.exists(ff_dir):
            return ff_dir
    
    # List available directories for better error message
    available_dirs = [d for d in os.listdir(output_dir) 
                     if os.path.isdir(os.path.join(output_dir, d))]
    
    raise FileNotFoundError(
        f"No force field directory found in {output_dir}\n"
        f"Expected one of: {', '.join(possible_names)}\n"
        f"Available directories: {', '.join(available_dirs)}"
    )


def parse_simulation_time(mdp_params):
    """
    Parse simulation time from MDP parameters.
    
    Parameters:
    -----------
    mdp_params : dict
        Dictionary of MDP parameters
        
    Returns:
    --------
    float
        Simulation time in nanoseconds
    """
    dt = float(mdp_params.get('dt', 0.002))  # ps
    nsteps = int(mdp_params.get('nsteps', 0))
    
    total_time_ps = dt * nsteps
    total_time_ns = total_time_ps / 1000.0
    
    return total_time_ns


def check_gpu_available():
    """
    Check if GPU is available for simulation.
    
    Returns:
    --------
    dict
        Dictionary with GPU availability information
    """
    gpu_info = {
        'cuda_available': False,
        'opencl_available': False,
        'devices': []
    }
    
    # Check for CUDA
    try:
        import subprocess
        result = subprocess.run(
            ['nvidia-smi', '--list-gpus'],
            capture_output=True,
            text=True,
            check=False
        )
        if result.returncode == 0:
            gpu_info['cuda_available'] = True
            gpu_info['devices'] = result.stdout.strip().split('\n')
    except:
        pass
    
    # Check for OpenCL
    try:
        import pyopencl as cl
        platforms = cl.get_platforms()
        if platforms:
            gpu_info['opencl_available'] = True
            for platform in platforms:
                devices = platform.get_devices()
                gpu_info['devices'].extend([str(d) for d in devices])
    except:
        pass
    
    return gpu_info


def estimate_simulation_time(mdp_params, system_size):
    """
    Estimate wall-clock time for simulation.
    
    Parameters:
    -----------
    mdp_params : dict
        Dictionary of MDP parameters
    system_size : int
        Number of atoms in the system
        
    Returns:
    --------
    dict
        Estimated times for each stage
    """
    estimates = {}
    
    # Very rough estimates (atoms/microsecond/day)
    gpu_performance = {
        'em': 1000000,    # Very fast for minimization
        'nvt': 100000,    # ~100k atoms/μs/day
        'npt': 90000,     # Slightly slower with barostat
        'prod': 85000     # Production run performance
    }
    
    cpu_performance = {
        'em': 100000,
        'nvt': 5000,
        'npt': 4500,
        'prod': 4000
    }
    
    # Use GPU estimates if available
    gpu_info = check_gpu_available()
    if gpu_info['cuda_available'] or gpu_info['opencl_available']:
        performance = gpu_performance
        hardware = "GPU"
    else:
        performance = cpu_performance
        hardware = "CPU"
    
    for stage, perf in performance.items():
        if stage in mdp_params:
            sim_time_ns = parse_simulation_time(mdp_params[stage])
            
            # Atoms * ns / (atoms/μs/day * 1000)
            wall_time_days = (system_size * sim_time_ns) / (perf * 1000)
            wall_time_hours = wall_time_days * 24
            
            estimates[stage] = {
                'simulation_time_ns': sim_time_ns,
                'estimated_hours': wall_time_hours,
                'hardware': hardware
            }
    
    return estimates


def create_checkpoint_script(gmx_dir, stage='prod'):
    """
    Create a checkpoint restart script for GROMACS.
    
    Parameters:
    -----------
    gmx_dir : str
        Path to GMX_PROLIG_MD directory
    stage : str
        Simulation stage to restart
        
    Returns:
    --------
    str
        Path to created script
    """
    script_content = f"""#!/bin/bash
# Restart {stage} simulation from checkpoint

cd {gmx_dir}

if [ -f ./{stage}/{stage}.cpt ]; then
    echo "Restarting {stage} from checkpoint..."
    gmx mdrun -s ./{stage}/{stage}.tpr \\
              -cpi ./{stage}/{stage}.cpt \\
              -deffnm ./{stage}/{stage} \\
              -ntmpi 1 -ntomp 15 \\
              -nb gpu -bonded gpu -pme gpu \\
              -gpu_id 0 -v
else
    echo "No checkpoint file found at ./{stage}/{stage}.cpt"
    exit 1
fi
"""
    
    script_path = os.path.join(gmx_dir, f'restart_{stage}.sh')
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_path, 0o755)
    return script_path