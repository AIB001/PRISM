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


def normalize_charge_reception(value: str) -> str:
    """Normalize legacy charge reception labels to the current API."""
    if value is None:
        return "surround"
    normalized = str(value).strip().lower()
    if normalized == "pert":
        return "surround"
    return normalized


def _normalize_atom_name(name: str) -> str:
    """Normalize atom names for cross-state matching."""
    name = str(name).strip()
    if name.endswith(("A", "B")) and len(name) > 1:
        return name[:-1]
    return name


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
        'none': keep original charges (no modification)
    charge_reception : str, optional
        Strategy for distributing surplus charges (default: 'surround')
        'unique': only unique (non-hydrogen) atoms receive charge
        'surround': only surrounding atoms receive charge
        'surround_ext': auto-extend to common atoms if charge/atom > 0.02
        'none': no charge redistribution
    recharge_hydrogen : bool, optional
        Whether to perturb hydrogen charges (default: False)
        If False, hydrogen charges are not modified
    """

    def __init__(
        self,
        dist_cutoff: float = 0.6,
        charge_cutoff: float = 0.05,
        charge_common: str = "mean",
        charge_reception: str = "surround",
        recharge_hydrogen: bool = False,
    ):
        if charge_common not in ["ref", "mut", "mean", "none"]:
            raise ValueError(f"Invalid charge_common: {charge_common}")
        self.dist_cutoff = dist_cutoff
        self.charge_cutoff = charge_cutoff
        self.charge_common = charge_common
        self.charge_reception = normalize_charge_reception(charge_reception)
        if self.charge_reception not in ["unique", "surround", "none"]:
            raise ValueError(f"Invalid charge_reception: {charge_reception}")
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
        params = dict(config.get("fep", {}).get("mapping", {}))
        if not params:
            model_cfg = config.get("model", {})
            other_cfg = config.get("other", {})
            params = {
                "dist_cutoff": other_cfg.get("dist_cutoff", 0.6),
                "charge_cutoff": other_cfg.get("charge_cutoff", 0.05),
                "charge_common": model_cfg.get("charge_common", "mean"),
                "charge_reception": model_cfg.get("charge_reception", "surround"),
                "recharge_hydrogen": other_cfg.get("recharge_hydrogen", False),
            }
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
        # Note: We modify the input atom objects in-place to apply charge_common strategy
        # This ensures HTML visualization shows the correct charges

        # Detect if using OpenFF/GAFF (generic atom types)
        # OpenFF types are like: output_0, output_1, ...
        # CGenFF types are like: CA, CB, CG, HA, HB, etc.
        uses_generic_types = self._uses_generic_atom_types(ligand_a, ligand_b)

        if uses_generic_types:
            # For OpenFF/GAFF: skip atom type check, only use distance + element + charge
            pass

        # Step 1: Distance-based matching
        # Build all valid candidates first, then greedily take the globally best
        # pairs in distance order. This avoids order-dependent mismatches where an
        # early acceptable match steals a much better partner from a later atom.
        common = []
        matched_a_indices = set()
        matched_b_indices = set()
        candidates = []

        for atom_a in ligand_a:
            for atom_b in ligand_b:
                dist = np.linalg.norm(atom_a.coord - atom_b.coord)
                if dist > self.dist_cutoff:
                    continue
                if atom_a.element != atom_b.element:
                    continue
                exact_name_match = atom_a.name == atom_b.name
                normalized_name_match = _normalize_atom_name(atom_a.name) == _normalize_atom_name(atom_b.name)
                type_match = atom_a.atom_type == atom_b.atom_type
                # Do not require atom-type identity at the candidate stage.
                # Otherwise atoms like C17 that keep the same position and name
                # but are retyped by CGenFF get forced into transformed_* instead
                # of becoming surrounding_*.
                candidates.append(
                    (
                        0 if exact_name_match else 1,
                        0 if normalized_name_match else 1,
                        0 if type_match else 1,
                        dist,
                        atom_a.index,
                        atom_b.index,
                        atom_a,
                        atom_b,
                    )
                )

        candidates.sort(key=lambda item: item[:6])

        for _, _, _, _, _, _, atom_a, atom_b in candidates:
            if atom_a.index in matched_a_indices or atom_b.index in matched_b_indices:
                continue
            common.append((atom_a, atom_b))
            matched_a_indices.add(atom_a.index)
            matched_b_indices.add(atom_b.index)

        # Step 2: Identify transformed atoms
        matched_a_indices = set(atom_a.index for atom_a, _ in common)
        transformed_a = [atom for atom in ligand_a if atom.index not in matched_a_indices]
        transformed_b = [atom for atom in ligand_b if atom.index not in matched_b_indices]

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

            # Check if charges differ
            if self.charge_common == "none":
                # In 'none' mode: any charge difference → surrounding
                # (no threshold, because we're not adjusting charges anyway)
                charge_differs = abs(atom_a.charge - atom_b.charge) > 1e-6
            else:
                # In ref/mut/mean modes: use charge_cutoff threshold
                # Filters out atoms with charge differences too large to adjust
                charge_differs = abs(atom_a.charge - atom_b.charge) > self.charge_cutoff

            if type_differs or charge_differs:
                # Move from common to surrounding
                common.remove((atom_a, atom_b))
                surrounding_a.append(atom_a)
                surrounding_b.append(atom_b)

        # Save original total charges before applying charge_common
        original_charge_a = sum(a.charge for a in ligand_a)
        original_charge_b = sum(b.charge for b in ligand_b)

        # Step 3.5: Apply charge_common strategy to common atoms
        # This modifies the atom objects in-place so HTML visualization shows correct charges
        if self.charge_common == "ref":
            for atom_a, atom_b in common:
                atom_b.charge = atom_a.charge
        elif self.charge_common == "mut":
            for atom_a, atom_b in common:
                atom_a.charge = atom_b.charge
        elif self.charge_common == "mean":
            for atom_a, atom_b in common:
                ave = (atom_a.charge + atom_b.charge) / 2.0
                atom_a.charge = ave
                atom_b.charge = ave
        elif self.charge_common == "none":
            # Don't modify charges - keep original values
            pass
        else:
            raise ValueError(f"Invalid charge_common: {self.charge_common}")

        # Step 3.6: Charge redistribution to maintain total charge conservation
        # After applying charge_common, the total charge may have changed
        # We need to redistribute the charge difference to maintain neutrality
        # CRITICAL: Only redistribute to non-common atoms to preserve pairing consistency
        if self.charge_common != "none":
            # Calculate current total charges (after charge_common)
            current_charge_a = sum(a.charge for a in ligand_a)
            current_charge_b = sum(b.charge for b in ligand_b)

            # Calculate calibration needed
            calibration_a = original_charge_a - current_charge_a
            calibration_b = original_charge_b - current_charge_b

            # Build modify_list: transformed + surrounding atoms (exclude common atoms)
            # This preserves the charge consistency of paired common atoms
            common_atoms_a = set(a for a, _ in common)
            common_atoms_b = set(b for _, b in common)

            modify_list_a = [a for a in ligand_a if a not in common_atoms_a]
            modify_list_b = [b for b in ligand_b if b not in common_atoms_b]

            # Redistribute charge
            if modify_list_a and abs(calibration_a) > 1e-10:
                ave_a = calibration_a / len(modify_list_a)
                for atom in modify_list_a:
                    atom.charge += ave_a

            if modify_list_b and abs(calibration_b) > 1e-10:
                ave_b = calibration_b / len(modify_list_b)
                for atom in modify_list_b:
                    atom.charge += ave_b

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
        - SwissParam-based (MMFF/MATCH/Both): Use generic atom types from RTF files

        Position-specific types (use atom type check):
        - CGenFF: CA, CB, CG, CD, HA, HB, ... (encode position in molecule)
        - GAFF: ca, c3, ha, hc, ... (encode chemical environment)

        Returns
        -------
        bool
            True if using generic/sequential types (OpenFF/OPLS/SwissParam), False if position-specific (CGenFF/GAFF)
        """
        # Check first few atoms for generic type patterns
        for atom in ligand_a[:5]:
            atom_type = atom.atom_type

            # OpenFF: output_0, output_1, ...
            if atom_type.startswith("output_"):
                return True
            # OPLS-AA: opls_800, opls_801, ...
            if atom_type.startswith("opls_"):
                return True
            # SwissParam (MMFF/MATCH): Check for RTF-style generic types
            # MMFF/MATCH atom types are typically like C_3, H_1, N_3, O_2, etc.
            # which are generic (element_number format), not position-specific
            if "_" in atom_type and len(atom_type) <= 4:
                # Format: C_3, H_1, N_3, O_2, etc. (element_number)
                # These are generic MMFF types, not position-specific like GAFF's ca/c3
                return True

        for atom in ligand_b[:5]:
            atom_type = atom.atom_type

            if atom_type.startswith("output_"):
                return True
            if atom_type.startswith("opls_"):
                return True
            if "_" in atom_type and len(atom_type) <= 4:
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
        Reclassify common hydrogens near uncommon atoms for generic atom types.

        OpenFF/OPLS-style sequential types do not provide reliable bonded-neighbor
        identity for deciding whether a hydrogen should remain in the common core.
        When a hydrogen sits within typical bond distance of a transformed or
        surrounding heavy atom on either side, keep it out of ``common`` and move it
        into ``surrounding`` so downstream visualization and charge handling still
        classify it explicitly instead of dropping it as unknown.
        """
        uncommon_a = transformed_a + surrounding_a
        uncommon_b = transformed_b + surrounding_b
        common_to_keep = []

        surrounding_a_indices = {atom.index for atom in surrounding_a}
        surrounding_b_indices = {atom.index for atom in surrounding_b}

        for atom_a, atom_b in common:
            if atom_a.element != "H" and atom_b.element != "H":
                common_to_keep.append((atom_a, atom_b))
                continue

            close_to_uncommon = False

            for uncommon_atom in uncommon_a:
                dist = np.linalg.norm(atom_a.coord - uncommon_atom.coord)
                if dist < 1.5:
                    close_to_uncommon = True
                    break

            if not close_to_uncommon:
                for uncommon_atom in uncommon_b:
                    dist = np.linalg.norm(atom_b.coord - uncommon_atom.coord)
                    if dist < 1.5:
                        close_to_uncommon = True
                        break

            if close_to_uncommon:
                if atom_a.index not in surrounding_a_indices:
                    surrounding_a.append(atom_a)
                    surrounding_a_indices.add(atom_a.index)
                if atom_b.index not in surrounding_b_indices:
                    surrounding_b.append(atom_b)
                    surrounding_b_indices.add(atom_b.index)
                continue

            common_to_keep.append((atom_a, atom_b))

        return common_to_keep
