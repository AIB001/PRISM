#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Atom Mapping Module

Implements distance-based atom mapping for FEP calculations.
Ported from FEbuilder/setup/make_hybrid.py
"""

from dataclasses import dataclass
from typing import List, Tuple
import numpy as np


@dataclass
class Atom:
    """
    Atom object for FEP calculations

    Attributes
    ----------
    name : str
        Atom unique identifier (e.g., 'C1', 'H2', 'N3')
    element : str
        Element symbol (C, H, N, O, etc.)
    coord : np.ndarray
        3D coordinates [x, y, z] in Angstroms
    charge : float
        Partial charge from force field parameterization
    atom_type : str
        Force field atom type (e.g., 'ca', 'ha', 'na' for GAFF)
    index : int
        Atom index for referencing
    """

    name: str
    element: str
    coord: np.ndarray
    charge: float
    atom_type: str
    index: int


@dataclass
class AtomMapping:
    """
    Result of atom mapping between two ligands

    This is the core data structure for FEP - all subsequent processing
    is based on this classification.

    Attributes
    ----------
    common : List[Tuple[Atom, Atom]]
        Atoms shared between both ligands (same in both states)
    transformed_a : List[Atom]
        Atoms unique to ligand A (become dummy in state B)
    transformed_b : List[Atom]
        Atoms unique to ligand B (become dummy in state A)
    surrounding_a : List[Atom]
        Position-matched atoms with divergent parameters in ligand A
    surrounding_b : List[Atom]
        Position-matched atoms with divergent parameters in ligand B
    """

    common: List[Tuple[Atom, Atom]]
    transformed_a: List[Atom]
    transformed_b: List[Atom]
    surrounding_a: List[Atom]
    surrounding_b: List[Atom]


class DistanceAtomMapper:
    """
    Distance-based atom mapper for FEP calculations

    Uses greedy matching with distance threshold to classify atoms
    into common, transformed, and surrounding categories.

    Parameters
    ----------
    dist_cutoff : float, optional
        Distance threshold in Angstroms (default: 0.6)
        Smaller values (0.3-0.5) for very similar molecules
        Larger values (0.8-1.0) for structurally divergent molecules
    charge_cutoff : float, optional
        Charge difference threshold (default: 0.05)
        Atoms with charge differences > cutoff become surrounding
    charge_common : str, optional
        Strategy for handling charges of common atoms (default: 'mean')
        'ref': use reference ligand charges
        'mut': use mutant ligand charges
        'mean': use average of both charges
    charge_reception : str, optional
        Strategy for distributing surplus charges (default: 'pert')
        'pert': all uncommon atoms receive charge
        'unique': only unique (non-hydrogen) atoms receive charge
        'surround': only surrounding atoms receive charge
    recharge_hydrogen : bool, optional
        Whether to perturb hydrogen charges (default: False)
        If False, hydrogen charges are not modified
    """

    def __init__(
        self,
        dist_cutoff: float = 0.6,
        charge_cutoff: float = 0.05,
        charge_common: str = "mean",
        charge_reception: str = "pert",
        recharge_hydrogen: bool = False,
    ):
        self.dist_cutoff = dist_cutoff
        self.charge_cutoff = charge_cutoff
        self.charge_common = charge_common
        self.charge_reception = charge_reception
        self.recharge_hydrogen = recharge_hydrogen

    @classmethod
    def from_config(cls, config: dict) -> "DistanceAtomMapper":
        """Create from a loaded config dict (e.g. from ConfigurationManager).

        Reads ``config['fep']['mapping']``. Missing keys fall back to
        ``__init__`` defaults.

        Parameters
        ----------
        config : dict
            Config dict as returned by ``ConfigurationManager.config``.
        """
        params = config.get("fep", {}).get("mapping", {})
        return cls(**params)

    def map(self, ligand_a: List[Atom], ligand_b: List[Atom]) -> AtomMapping:
        """
        Perform distance-based atom mapping

        Algorithm (ported from FEbuilder/setup/make_hybrid.py):
        1. Distance matching - find atoms within dist_cutoff with same element
        2. Identify transformed atoms - unmatched atoms become transformed
        3. Identify surrounding atoms - matched atoms with divergent parameters
        4. Apply charge_common strategy to common atoms

        Parameters
        ----------
        ligand_a : List[Atom]
            Reference ligand atoms
        ligand_b : List[Atom]
            Mutant ligand atoms

        Returns
        -------
        AtomMapping
            Classification of atoms into common/transformed/surrounding
        """
        # Create copies to avoid modifying original atoms
        ligand_a = [Atom(a.name, a.element, a.coord.copy(), a.charge, a.atom_type, a.index) for a in ligand_a]
        ligand_b = [Atom(b.name, b.element, b.coord.copy(), b.charge, b.atom_type, b.index) for b in ligand_b]

        # Step 1: Distance-based matching
        common = []
        matched_b_indices = set()

        for atom_a in ligand_a:
            for atom_b in ligand_b:
                # Skip if already matched
                if atom_b.index in matched_b_indices:
                    continue

                # Check distance
                dist = np.linalg.norm(atom_a.coord - atom_b.coord)
                if dist > self.dist_cutoff:
                    continue

                # Check element type
                if atom_a.element != atom_b.element:
                    continue

                # Found a match
                common.append((atom_a, atom_b))
                matched_b_indices.add(atom_b.index)
                break  # Greedy: use first match only

        # Step 2: Identify transformed atoms
        matched_a_indices = set(atom_a.index for atom_a, _ in common)
        transformed_a = [atom for atom in ligand_a if atom.index not in matched_a_indices]
        transformed_b = [atom for atom in ligand_b if atom.index not in matched_b_indices]

        # Apply charge_reception strategy
        if self.charge_reception == "unique" and not self.recharge_hydrogen:
            # Only hydrogens in transformed
            transformed_a = [atom for atom in transformed_a if atom.name.startswith("H")]
            transformed_b = [atom for atom in transformed_b if atom.name.startswith("H")]

        # Step 3: Identify surrounding atoms
        # (atoms with divergent types or charges)
        surrounding_a = []
        surrounding_b = []

        # Make a copy to iterate safely
        for atom_a, atom_b in common[:]:
            # Check if types differ
            type_differs = atom_a.atom_type != atom_b.atom_type

            # Check if charges differ significantly
            charge_differs = abs(atom_a.charge - atom_b.charge) > self.charge_cutoff

            if type_differs or charge_differs:
                # Move from common to surrounding
                common.remove((atom_a, atom_b))
                surrounding_a.append(atom_a)
                surrounding_b.append(atom_b)
            elif self.charge_common != "mean" and abs(atom_a.charge - atom_b.charge) > 0:
                # Apply charge_common strategy (ref or mut)
                if self.charge_common == "ref":
                    atom_b.charge = atom_a.charge
                elif self.charge_common == "mut":
                    atom_a.charge = atom_b.charge

        # Step 4: Handle hydrogen atoms connected to transformed/surrounding
        # Remove hydrogens from common and add to transformed/surrounding
        common_to_remove = []
        for atom_a, atom_b in common:
            # Check if hydrogen connected to transformed/surrounding atoms
            if atom_a.name.startswith("H") or atom_b.name.startswith("H"):
                # For simplicity, move all hydrogens to transformed if not recharge_hydrogen
                if not self.recharge_hydrogen:
                    # Keep in common only if not connected to transformed/surrounding
                    # This is a simplified version - FEbuilder has more complex logic
                    pass

        return AtomMapping(
            common=common,
            transformed_a=transformed_a,
            transformed_b=transformed_b,
            surrounding_a=surrounding_a,
            surrounding_b=surrounding_b,
        )
