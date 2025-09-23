#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM - Protein Receptor Interaction Simulation Modeler

A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
"""

__version__ = "1.0.0"
__author__ = "PRISM Development Team"

from .builder import PRISMBuilder
from .core import PRISMSystem
from .sim import model

# High-level API functions
def system(protein_path, ligand_path, config=None, **kwargs):
    """
    Create a protein-ligand system for MD simulation.

    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)
    config : str or dict, optional
        Configuration file path or configuration dictionary
    **kwargs : optional
        Additional parameters (output_dir, ligand_forcefield, forcefield, etc.)

    Returns:
    --------
    PRISMSystem
        A system object with methods to build and run simulations

    Examples:
    ---------
    >>> import prism as pm
    >>> system = pm.system("protein.pdb", "ligand.mol2")
    >>> system.build()
    >>> system.generate_mdp_files()

    >>> # With custom configuration
    >>> system = pm.system("protein.pdb", "ligand.sdf",
    ...                    config="config.yaml",
    ...                    ligand_forcefield="openff")
    >>> system.build()
    """
    return PRISMSystem(protein_path, ligand_path, config=config, **kwargs)

def build_system(protein_path, ligand_path, output_dir="prism_output", **kwargs):
    """
    Build a complete protein-ligand system (one-step function).

    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)
    output_dir : str
        Output directory for generated files
    **kwargs : optional
        Additional parameters (forcefield, water_model, ligand_forcefield, etc.)

    Returns:
    --------
    str
        Path to the output directory

    Examples:
    ---------
    >>> import prism as pm
    >>> output_path = pm.build_system("protein.pdb", "ligand.mol2")

    >>> # With custom force fields
    >>> output_path = pm.build_system("protein.pdb", "ligand.mol2",
    ...                               forcefield="amber99sb-ildn",
    ...                               water_model="tip4p",
    ...                               ligand_forcefield="openff")
    """
    system_obj = PRISMSystem(protein_path, ligand_path, output_dir=output_dir, **kwargs)
    return system_obj.build()

def check_dependencies():
    """
    Check if all required dependencies are available.

    Returns:
    --------
    dict
        Dictionary showing availability of each dependency

    Example:
    --------
    >>> import prism as pm
    >>> deps = pm.check_dependencies()
    >>> print(deps)
    {'gromacs': True, 'pdbfixer': True, 'antechamber': True, 'openff': False}
    """
    from .utils.environment import GromacsEnvironment
    import subprocess

    dependencies = {
        'gromacs': False,
        'pdbfixer': False,
        'antechamber': False,
        'openff': False
    }

    # Check GROMACS
    try:
        env = GromacsEnvironment()
        dependencies['gromacs'] = True
    except:
        pass

    # Check pdbfixer
    try:
        import pdbfixer
        dependencies['pdbfixer'] = True
    except:
        pass

    # Check for AmberTools (antechamber)
    try:
        subprocess.run(['antechamber', '-h'], capture_output=True, check=True)
        dependencies['antechamber'] = True
    except:
        pass

    # Check OpenFF
    try:
        import openff.toolkit
        dependencies['openff'] = True
    except:
        pass

    return dependencies

def list_forcefields():
    """
    List available force fields from GROMACS installation.

    Returns:
    --------
    list
        List of available force field names

    Example:
    --------
    >>> import prism as pm
    >>> ffs = pm.list_forcefields()
    >>> for ff in ffs:
    ...     print(ff)
    amber99sb (index: 1)
    amber99sb-ildn (index: 2)
    amber14sb (index: 3)
    """
    from .utils.environment import GromacsEnvironment

    try:
        env = GromacsEnvironment()
        return env.list_force_fields()
    except Exception as e:
        print(f"Error detecting force fields: {e}")
        return []

def get_version():
    """Get PRISM version."""
    return __version__

# Export main classes and functions
__all__ = [
    "PRISMBuilder",
    "PRISMSystem",
    "system",
    "build_system",
    "model",
    "check_dependencies",
    "list_forcefields",
    "get_version",
    "__version__"
]