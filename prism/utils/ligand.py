#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand Detection Utilities for PRISM

This module provides utilities for automatically detecting and working with ligand
residues in protein-ligand systems for MD analysis and trajectory processing.
"""

import logging
from typing import Optional, List, Union, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

# Common ligand residue names
COMMON_LIGAND_NAMES = ['LIG', 'UNL', 'MOL', 'DRG', 'INH', 'SUB', 'HET']

# Standard biological residues and common solvent/ion names
STANDARD_RESIDUES = {
    # Amino acids
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Nucleotides
    'DA', 'DT', 'DG', 'DC', 'RA', 'RU', 'RG', 'RC',
    # Water and solvents
    'WAT', 'HOH', 'TIP3', 'TIP4', 'SOL', 'DMSO',
    # Common ions
    'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'CU', 'MN'
}

def identify_ligand_residue(traj, manual_ligand: Optional[str] = None):
    """
    Automatically identify ligand residue in a trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory object with topology information
    manual_ligand : str, optional
        Manually specified ligand residue name to look for

    Returns
    -------
    residue or None
        The identified ligand residue object, or None if not found

    Examples
    --------
    >>> import mdtraj as md
    >>> traj = md.load('trajectory.xtc', top='system.pdb')
    >>> ligand = identify_ligand_residue(traj)
    >>> if ligand:
    ...     print(f"Found ligand: {ligand.name}{ligand.resSeq}")
    >>>
    >>> # Manually specify ligand
    >>> ligand = identify_ligand_residue(traj, manual_ligand='LIG')
    """
    if manual_ligand:
        # Look for manually specified ligand
        for residue in traj.topology.residues:
            if residue.name == manual_ligand:
                logger.info(f"Found manually specified ligand: {residue.name}{residue.resSeq}")
                return residue
        logger.warning(f"Manually specified ligand '{manual_ligand}' not found")
        return None

    # Automatic detection
    ligand_residues = []

    # Step 1: Look for common ligand names
    for residue in traj.topology.residues:
        if (residue.name in COMMON_LIGAND_NAMES or
            (residue.name not in STANDARD_RESIDUES and len(list(residue.atoms)) > 5)):
            ligand_residues.append(residue)

    # Step 2: If no obvious ligands found, find largest non-standard residue
    if not ligand_residues:
        logger.info("No obvious ligand residues found, looking for largest non-standard residue")
        largest_residue = None
        max_atoms = 0
        for residue in traj.topology.residues:
            if residue.name not in STANDARD_RESIDUES:
                n_atoms = len(list(residue.atoms))
                if n_atoms > max_atoms:
                    max_atoms = n_atoms
                    largest_residue = residue

        if largest_residue:
            ligand_residues = [largest_residue]
            logger.info(f"Using largest non-standard residue as ligand: {largest_residue.name}{largest_residue.resSeq} ({max_atoms} atoms)")

    if not ligand_residues:
        logger.warning("No ligand residue found automatically")
        return None

    # Return first found ligand
    ligand = ligand_residues[0]
    logger.info(f"Auto-detected ligand: {ligand.name}{ligand.resSeq}")
    return ligand

def get_ligand_selection_string(traj, manual_ligand: Optional[str] = None) -> Optional[str]:
    """
    Get GROMACS selection string for ligand residue.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory object with topology information
    manual_ligand : str, optional
        Manually specified ligand residue name

    Returns
    -------
    str or None
        GROMACS selection string like "r_LIG" or None if no ligand found

    Examples
    --------
    >>> selection = get_ligand_selection_string(traj)
    >>> print(selection)  # "r_LIG"
    """
    ligand = identify_ligand_residue(traj, manual_ligand)
    if ligand:
        return f"r_{ligand.name}"
    return None

def validate_ligand_residue(traj, ligand_name: str) -> bool:
    """
    Check if a ligand residue exists in the trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory object with topology information
    ligand_name : str
        Name of ligand residue to validate

    Returns
    -------
    bool
        True if ligand exists, False otherwise
    """
    for residue in traj.topology.residues:
        if residue.name == ligand_name:
            return True
    return False

def get_ligand_center_atom(traj, manual_ligand: Optional[str] = None) -> Optional[int]:
    """
    Get the best atom index for centering the ligand in PBC processing.

    Tries to find a central heavy atom that's good for PBC centering.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory object with topology information
    manual_ligand : str, optional
        Manually specified ligand residue name

    Returns
    -------
    int or None
        Atom index suitable for PBC centering, or None if no ligand found

    Examples
    --------
    >>> center_atom = get_ligand_center_atom(traj)
    >>> if center_atom is not None:
    ...     print(f"Use atom {center_atom} for PBC centering")
    """
    ligand = identify_ligand_residue(traj, manual_ligand)
    if not ligand:
        return None

    # Get heavy atoms from ligand
    heavy_atoms = []
    for atom in ligand.atoms:
        if atom.element.symbol != 'H':  # Skip hydrogens
            heavy_atoms.append(atom)

    if not heavy_atoms:
        logger.warning(f"No heavy atoms found in ligand {ligand.name}")
        return None

    # For PBC centering, try to find a central atom
    # Simple heuristic: use the middle atom by index
    center_idx = len(heavy_atoms) // 2
    center_atom = heavy_atoms[center_idx]

    logger.info(f"Selected atom {center_atom.index} ({center_atom.name}) for PBC centering")
    return center_atom.index

def list_ligand_candidates(traj) -> List[Tuple[str, int, int]]:
    """
    List all potential ligand residues in the trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory object with topology information

    Returns
    -------
    List[Tuple[str, int, int]]
        List of (residue_name, residue_id, n_atoms) for potential ligands

    Examples
    --------
    >>> candidates = list_ligand_candidates(traj)
    >>> for name, resid, n_atoms in candidates:
    ...     print(f"{name}{resid}: {n_atoms} atoms")
    """
    candidates = []

    for residue in traj.topology.residues:
        if (residue.name in COMMON_LIGAND_NAMES or
            (residue.name not in STANDARD_RESIDUES and len(list(residue.atoms)) > 5)):
            candidates.append((residue.name, residue.resSeq, len(list(residue.atoms))))

    # Sort by atom count (largest first)
    candidates.sort(key=lambda x: x[2], reverse=True)
    return candidates

def get_ligand_info(traj, manual_ligand: Optional[str] = None) -> Optional[dict]:
    """
    Get comprehensive information about the identified ligand.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory object with topology information
    manual_ligand : str, optional
        Manually specified ligand residue name

    Returns
    -------
    dict or None
        Dictionary with ligand information or None if no ligand found

    Examples
    --------
    >>> info = get_ligand_info(traj)
    >>> if info:
    ...     print(f"Ligand: {info['name']}, Atoms: {info['n_atoms']}")
    """
    ligand = identify_ligand_residue(traj, manual_ligand)
    if not ligand:
        return None

    heavy_atoms = [atom for atom in ligand.atoms if atom.element.symbol != 'H']

    return {
        'name': ligand.name,
        'resid': ligand.resSeq,
        'n_atoms': len(list(ligand.atoms)),
        'n_heavy_atoms': len(heavy_atoms),
        'gromacs_selection': f"r_{ligand.name}",
        'center_atom_index': get_ligand_center_atom(traj, manual_ligand),
        'residue_object': ligand
    }