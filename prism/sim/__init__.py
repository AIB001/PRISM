#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Simulation Module - Run MD simulations using GROMACS or OpenMM
"""

from .base import SimulationModel
from .gmxsim import GMXSimulator
from .openmmsim import OpenMMSimulator

# High-level API function
def model(gmx_dir):
    """
    Create a simulation model from a GROMACS directory.
    
    Parameters:
    -----------
    gmx_dir : str
        Path to the GMX_PROLIG_MD directory
        
    Returns:
    --------
    SimulationModel
        A simulation model object ready to run MD
        
    Examples:
    ---------
    >>> import prism as pm
    >>> sim = pm.model("output/GMX_PROLIG_MD")
    >>> sim.run(engine="gmx")
    """
    return SimulationModel(gmx_dir)

__all__ = [
    "SimulationModel",
    "GMXSimulator", 
    "OpenMMSimulator",
    "model"
]