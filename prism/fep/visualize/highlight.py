#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Highlight information generation for FEP atom mapping visualization.

Generates color-coded highlighting for three atom types:
- Common (green): Shared between both ligands
- Transformed (red): Require dummy atom types
- Surrounding (blue): Only need typeB/chargeB
"""

from typing import Dict, List
from rdkit import Chem


def create_highlight_info(
    mol: Chem.Mol,
    common_atoms: List[str],
    transformed_atoms: List[str],
    surrounding_atoms: List[str],
) -> Dict:
    """
    Generate RDKit highlighting information for atom classification.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit Mol object (should have atom name labels)
    common_atoms : List[str]
        Atom names classified as common
    transformed_atoms : List[str]
        Atom names classified as transformed
    surrounding_atoms : List[str]
        Atom names classified as surrounding

    Returns
    -------
    dict
        Dictionary containing:
        - 'mol': RDKit Mol object
        - 'hlist': List of atom indices to highlight
        - 'hmap': Dictionary mapping atom index to RGB color tuple

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles('CCO')
    >>> info = create_highlight_info(mol, ['C1'], [], ['O1'])
    >>> info['hlist']  # [0, 2]
    >>> info['hmap']   # {0: (0.8, 0.9, 0.3), 2: (0.3, 0.6, 1.0)}
    """
    # Build atom name to index mapping
    atom_names = {}
    for atom in mol.GetAtoms():
        name = atom.GetProp("name") if atom.HasProp("name") else ""
        if name:
            atom_names[name] = atom.GetIdx()

    # Define colors for each classification
    colors = {
        "common": (0.8, 0.9, 0.3),  # Green (FEbuilder didn't have this)
        "transformed": (1, 0.3, 0.3),  # Red (FEbuilder didn't have this)
        "surrounding": (0.3, 0.6, 1),  # Blue (FEbuilder didn't have this)
    }

    # Collect highlighted atoms
    hlist = []
    hmap = {}

    for name in common_atoms:
        if name in atom_names:
            idx = atom_names[name]
            hlist.append(idx)
            hmap[idx] = colors["common"]

    for name in transformed_atoms:
        if name in atom_names:
            idx = atom_names[name]
            hlist.append(idx)
            hmap[idx] = colors["transformed"]

    for name in surrounding_atoms:
        if name in atom_names:
            idx = atom_names[name]
            hlist.append(idx)
            hmap[idx] = colors["surrounding"]

    return {
        "mol": mol,
        "hlist": hlist,
        "hmap": hmap,
    }
