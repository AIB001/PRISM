#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM-FEP Utility Functions

Common utility functions for FEP module.
"""

from typing import List, Tuple
import numpy as np


def calculate_distance(atom1, atom2) -> float:
    """
    Calculate Euclidean distance between two atoms

    Parameters
    ----------
    atom1, atom2 : Atom
        Atom objects with coord attribute

    Returns
    -------
    float
        Distance in Angstroms
    """
    return np.linalg.norm(atom1.coord - atom2.coord)


def get_connected_atoms(atom: "Atom", bonds: List[Tuple]) -> List[str]:
    """
    Get list of connected atom names

    Parameters
    ----------
    atom : Atom
        Central atom
    bonds : List[Tuple]
        List of (atom1_name, atom2_name) tuples

    Returns
    -------
    List[str]
        Names of connected atoms
    """
    connected = []
    for bond in bonds:
        if atom.name in bond:
            other_atom = bond[1] if bond[0] == atom.name else bond[0]
            connected.append(other_atom)
    return connected


def format_charge(charge: float, precision: int = 4) -> str:
    """
    Format charge value for ITP file output

    Parameters
    ----------
    charge : float
        Charge value
    precision : int, optional
        Number of decimal places (default: 4)

    Returns
    -------
    str
        Formatted charge string
    """
    return f"{charge:.{precision}f}"


def validate_atoms(atoms: List) -> bool:
    """
    Validate list of Atom objects

    Parameters
    ----------
    atoms : List
        List of Atom objects

    Returns
    -------
    bool
        True if all atoms are valid, False otherwise

    Raises
    ------
    ValueError
        If atoms are invalid
    """
    if not atoms:
        raise ValueError("Atom list cannot be empty")

    for atom in atoms:
        if not hasattr(atom, "coord") or atom.coord is None:
            raise ValueError(f"Atom {atom} missing coordinate data")
        if not hasattr(atom, "charge"):
            raise ValueError(f"Atom {atom} missing charge data")
        if not hasattr(atom, "atom_type"):
            raise ValueError(f"Atom {atom} missing atom type")

    return True
