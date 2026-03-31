#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Atomic properties for FEP calculations

This module provides atomic property lookups (mass, radius, etc.) for FEP modeling.
"""

from typing import Dict


# Standard atomic masses (in amu)
_ATOMIC_MASSES: Dict[str, float] = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "P": 30.974,
    "S": 32.065,
    "Cl": 35.453,
    "Br": 79.904,
    "I": 126.904,
    # Common metals
    "Mg": 24.305,
    "Ca": 40.078,
    "Zn": 65.38,
    "Fe": 55.845,
    "Cu": 63.546,
    "Mn": 54.938,
}


def get_atomic_mass(element: str) -> float:
    """
    Get standard atomic mass for an element.

    Parameters
    ----------
    element : str
        Element symbol (e.g., 'H', 'C', 'N', 'O', 'Cl', 'Br')

    Returns
    -------
    float
        Atomic mass in atomic mass units (amu)

    Notes
    -----
    Values are from IUPAC atomic weights (most common isotope).
    For elements not in the table, defaults to carbon (12.011 amu).

    Examples
    --------
    >>> get_atomic_mass('H')
    1.008
    >>> get_atomic_mass('Cl')
    35.453
    """
    return _ATOMIC_MASSES.get(element, 12.011)  # Default to carbon


def build_mass_dict_from_atoms(atoms, mass_dict: Dict[str, float] = None) -> Dict[str, float]:
    """
    Build atom_type -> mass mapping dictionary from Atom objects.

    This is a convenience function for building the 'masses' parameter
    dictionary required by HybridTopologyBuilder.

    Parameters
    ----------
    atoms : List[Atom]
        List of Atom objects with element attribute
    mass_dict : Dict[str, float], optional
        Existing mass dictionary to update. If None, creates new dict.

    Returns
    -------
    Dict[str, float]
        Dictionary mapping atom_type -> atomic mass

    Examples
    --------
    >>> from prism.fep.io import read_ligand_from_prism
    >>> atoms = read_ligand_from_prism("LIG.itp", "LIG.gro")
    >>> params = {'masses': build_mass_dict_from_atoms(atoms)}
    """
    if mass_dict is None:
        mass_dict = {}

    for atom in atoms:
        if atom.atom_type not in mass_dict:
            mass_dict[atom.atom_type] = get_atomic_mass(atom.element)

    return mass_dict
