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

    def __eq__(self, other):
        """Custom equality comparison for Atom objects"""
        if not isinstance(other, Atom):
            return False
        return (
            self.name == other.name
            and self.element == other.element
            and np.allclose(self.coord, other.coord, atol=1e-6)
            and abs(self.charge - other.charge) < 1e-6
            and self.atom_type == other.atom_type
            and self.index == other.index
        )

    def __hash__(self):
        """Custom hash for Atom objects"""
        # Use tuple of hashable components
        coord_tuple = tuple(self.coord.round(6))  # Round to avoid floating point issues
        return hash((self.name, self.element, coord_tuple, round(self.charge, 6), self.atom_type, self.index))


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

        Note: For OpenFF/GAFF force fields, atom_type check is disabled since
        these force fields use generic types (output_0, output_1, ...) that don't
        encode positional information like CGenFF (CA, CB, CG, etc.).

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

        # Detect if using OpenFF/GAFF (generic atom types)
        # OpenFF types are like: output_0, output_1, ...
        # CGenFF types are like: CA, CB, CG, HA, HB, etc.
        uses_generic_types = self._uses_generic_atom_types(ligand_a, ligand_b)

        if uses_generic_types:
            # For OpenFF/GAFF: skip atom type check, only use distance + element + charge
            pass

        # Step 1: Distance-based matching (greedy, like FEbuilder)
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

                # Check element
                if atom_a.element != atom_b.element:
                    continue

                # Check atom type (skip only for truly generic types like OpenFF's output_X)
                # GAFF: ca, c3, ha, hc, ... (position-specific, should check)
                # CGenFF: CA, CB, CG, HA, HB, ... (position-specific, should check)
                # OpenFF: output_0, output_1, ... (generic, skip check)
                if not uses_generic_types:
                    if atom_a.atom_type != atom_b.atom_type:
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
            # Check if types differ (only for CGenFF, not OpenFF/GAFF)
            if not uses_generic_types:
                type_differs = atom_a.atom_type != atom_b.atom_type
            else:
                type_differs = False  # Skip type check for OpenFF/GAFF

            # Check if charges differ significantly
            charge_differs = abs(atom_a.charge - atom_b.charge) > self.charge_cutoff

            if type_differs or charge_differs:
                # Move from common to surrounding
                common.remove((atom_a, atom_b))
                surrounding_a.append(atom_a)
                surrounding_b.append(atom_b)

        # Note: charge_common strategy is applied by HybridTopologyBuilder, not here
        # DistanceAtomMapper only classifies atoms into common/transformed/surrounding

        # Step 4: Handle hydrogen atoms connected to transformed/surrounding
        # For OpenFF/GAFF: use distance-based proximity to transformed/surrounding
        if uses_generic_types:
            # Remove hydrogens that are close to transformed/surrounding atoms
            common = self._filter_hydrogens_for_generic_types(
                common, transformed_a, transformed_b, surrounding_a, surrounding_b
            )
        else:
            # CGenFF: use bond connectivity (simplified version)
            # TODO: implement full bond-based logic like FEbuilder
            pass

        return AtomMapping(
            common=common,
            transformed_a=transformed_a,
            transformed_b=transformed_b,
            surrounding_a=surrounding_a,
            surrounding_b=surrounding_b,
        )

    def _uses_generic_atom_types(self, ligand_a: List[Atom], ligand_b: List[Atom]) -> bool:
        """
        Detect if force field uses generic/sequential atom types that should skip type checking.

        Generic/Sequential types (skip atom type check):
        - OpenFF: output_0, output_1, ... (generic, not position-specific)
        - OPLS-AA: opls_800, opls_801, ... (sequential numbering, not chemical environment)

        Position-specific types (use atom type check):
        - CGenFF: CA, CB, CG, CD, HA, HB, ... (encode position in molecule)
        - GAFF: ca, c3, ha, hc, ... (encode chemical environment)

        Returns
        -------
        bool
            True if using generic/sequential types (OpenFF/OPLS), False if position-specific (CGenFF/GAFF)
        """
        # Check first few atoms
        for atom in ligand_a[:5]:
            if atom.atom_type.startswith("output_"):
                return True
            if atom.atom_type.startswith("opls_"):
                return True

        for atom in ligand_b[:5]:
            if atom.atom_type.startswith("output_"):
                return True
            if atom.atom_type.startswith("opls_"):
                return True

        return False

    def _filter_hydrogens_for_generic_types(
        self,
        common: List[Tuple[Atom, Atom]],
        transformed_a: List[Atom],
        transformed_b: List[Atom],
        surrounding_a: List[Atom],
        surrounding_b: List[Atom],
    ) -> List[Tuple[Atom, Atom]]:
        """
        For OpenFF/GAFF: remove hydrogens from common if they're close to transformed/surrounding.

        This approximates FEbuilder's bond-based logic using distance.
        """
        # Collect all uncommon atoms
        uncommon_a = transformed_a + surrounding_a
        uncommon_b = transformed_b + surrounding_b

        # Hydrogens to keep in common (not close to uncommon atoms)
        common_to_keep = []

        for atom_a, atom_b in common:
            # Check if this is a hydrogen pair
            if not (atom_a.name.startswith("H") or atom_b.name.startswith("H")):
                # Non-hydrogen: always keep
                common_to_keep.append((atom_a, atom_b))
                continue

            # For hydrogens: check if close to any uncommon atom
            # If close to uncommon, should be in transformed/surrounding, not common
            close_to_uncommon = False

            for uncommon_atom in uncommon_a:
                dist = np.linalg.norm(atom_a.coord - uncommon_atom.coord)
                if dist < 1.5:  # Typical C-H bond length
                    close_to_uncommon = True
                    break

            if not close_to_uncommon:
                # Also check the B side
                for uncommon_atom in uncommon_b:
                    dist = np.linalg.norm(atom_b.coord - uncommon_atom.coord)
                    if dist < 1.5:
                        close_to_uncommon = True
                        break

            if not close_to_uncommon:
                # Hydrogen is not close to any uncommon atom, keep in common
                common_to_keep.append((atom_a, atom_b))

        return common_to_keep
