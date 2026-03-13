#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Hybrid Topology Construction Module

Builds hybrid topologies with A/B state parameters for GROMACS FEP calculations.
GROMACS uses single-topology approach with typeB/chargeB encoding.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional
from .mapping import AtomMapping


@dataclass
class HybridAtom:
    """
    Hybrid atom containing A/B states for FEP calculations

    Maps directly to a row in GROMACS ITP file.

    Attributes
    ----------
    name : str
        Atom name (original for common, tagged like 'C1A'/'C1B' for transformed)
    index : int
        Atom index (starts from 1, GROMACS convention)
    state_a_type : str
        Force field type for state A
    state_a_charge : float
        Charge for state A
    state_b_type : Optional[str], optional
        Force field type for state B (None for common atoms)
    state_b_charge : Optional[float], optional
        Charge for state B (None for common atoms)
    element : str, optional
        Element symbol (for mass lookup)
    mass : float, optional
        Atomic mass for state A
    mass_b : Optional[float], optional
        Atomic mass for state B (for pmx protein mutation topologies)
    """

    name: str
    index: int
    state_a_type: str
    state_a_charge: float
    state_b_type: Optional[str] = None
    state_b_charge: Optional[float] = None
    element: str = ""
    mass: float = 0.0
    mass_b: Optional[float] = None


class HybridTopologyBuilder:
    """
    Builds hybrid topology for GROMACS FEP calculations

    Creates GROMACS ITP files with typeB/chargeB columns encoding
    both reference and mutant states using single-topology approach.

    Parameters
    ----------
    charge_strategy : str, optional
        Strategy for handling common atom charges (default: 'mean')
        - 'mean': Use average of reference and mutant charges
        - 'ref': Use reference ligand charges
        - 'mut': Use mutant ligand charges

    Notes
    -----
    GROMACS uses single-topology approach: one topology with typeB/chargeB
    columns to encode both states. This is different from NAMD's dual-topology
    approach which requires two separate topologies.
    """

    def __init__(self, charge_strategy: str = "mean", charge_reception: str = "pert", recharge_hydrogen: bool = False):
        if charge_strategy not in ["ref", "mut", "mean"]:
            raise ValueError(f"Invalid charge_strategy: {charge_strategy}")
        if charge_reception not in ["pert", "unique", "surround", "none"]:
            raise ValueError(f"Invalid charge_reception: {charge_reception}")
        self.charge_strategy = charge_strategy
        self.charge_reception = charge_reception
        self.recharge_hydrogen = recharge_hydrogen
        self.hybrid_atoms = []

    def build(
        self, mapping: AtomMapping, params_a: Dict, params_b: Dict, atoms_a: List, atoms_b: List
    ) -> List[HybridAtom]:
        """
        Build hybrid topology with charge redistribution

        Parameters
        ----------
        mapping : AtomMapping
            Atom mapping classification
        params_a : Dict
            Force field parameters for ligand A
            Contains 'masses', 'bonds', 'angles', 'dihedrals', etc.
        params_b : Dict
            Force field parameters for ligand B
        atoms_a : List[Atom]
            Original atoms from ligand A (for charge calculation)
        atoms_b : List[Atom]
            Original atoms from ligand B (for charge calculation)

        Returns
        -------
        List[HybridAtom]
            List of hybrid atoms with typeB/chargeB encoding
        """
        self.hybrid_atoms = []
        index = 1  # GROMACS uses 1-based indexing

        # Get mass dictionaries
        masses_a = params_a.get("masses", {})
        masses_b = params_b.get("masses", {})

        # 1. Common atoms (shared structure, same or similar charges)
        for atom_a, atom_b in mapping.common:
            charge = self._resolve_charge(atom_a.charge, atom_b.charge)

            # Common atoms: typeA = typeB, chargeA = chargeB (after averaging)
            self.hybrid_atoms.append(
                HybridAtom(
                    name=atom_a.name,
                    index=index,
                    state_a_type=atom_a.atom_type,
                    state_a_charge=charge,
                    state_b_type=atom_a.atom_type,  # Same type
                    state_b_charge=charge,  # Same charge
                    element=atom_a.element,
                    mass=masses_a.get(atom_a.atom_type, self._get_default_mass(atom_a.element)),
                )
            )
            index += 1

        # 2. Transformed A atoms (disappear: A -> dummy)
        for atom_a in mapping.transformed_a:
            dummy_type = self._get_dummy_type(atom_a.atom_type)

            self.hybrid_atoms.append(
                HybridAtom(
                    name=f"{atom_a.name}A",  # Add suffix for clarity
                    index=index,
                    state_a_type=atom_a.atom_type,
                    state_a_charge=atom_a.charge,
                    state_b_type=dummy_type,  # Dummy in state B
                    state_b_charge=0.0,
                    element=atom_a.element,
                    mass=masses_a.get(atom_a.atom_type, self._get_default_mass(atom_a.element)),
                )
            )
            index += 1

        # 3. Transformed B atoms (appear: dummy -> B)
        for atom_b in mapping.transformed_b:
            dummy_type = self._get_dummy_type(atom_b.atom_type)

            self.hybrid_atoms.append(
                HybridAtom(
                    name=f"{atom_b.name}B",  # Add suffix for clarity
                    index=index,
                    state_a_type=dummy_type,  # Dummy in state A
                    state_a_charge=0.0,
                    state_b_type=atom_b.atom_type,
                    state_b_charge=atom_b.charge,
                    element=atom_b.element,
                    mass=masses_b.get(atom_b.atom_type, self._get_default_mass(atom_b.element)),
                )
            )
            index += 1

        # 4. Surrounding atoms (same position, different charges/types)
        # CRITICAL: Surrounding atoms MUST have typeB/chargeB because their charges
        # will be different in states A and B after charge redistribution
        for atom_a, atom_b in zip(mapping.surrounding_a, mapping.surrounding_b):
            self.hybrid_atoms.append(
                HybridAtom(
                    name=atom_a.name,
                    index=index,
                    state_a_type=atom_a.atom_type,
                    state_a_charge=atom_a.charge,  # Will be adjusted by redistribution
                    state_b_type=atom_b.atom_type,  # MUST have typeB
                    state_b_charge=atom_b.charge,  # MUST have chargeB (different!)
                    element=atom_a.element,
                    mass=masses_a.get(atom_a.atom_type, self._get_default_mass(atom_a.element)),
                )
            )
            index += 1

        # . Apply charge redistribution to maintain neutrality
        if self.charge_reception != "none":
            self._redistribute_charges(mapping, atoms_a, atoms_b)

        return self.hybrid_atoms

    def _redistribute_charges(self, mapping: AtomMapping, atoms_a: List, atoms_b: List) -> None:
        """
        Redistribute charges to maintain system neutrality

        After applying charge_strategy to common atoms, the total system charge
        changes. This method redistributes the charge difference to selected atoms
        according to charge_reception strategy.

        Parameters
        ----------
        mapping : AtomMapping
            Atom mapping classification
        atoms_a : List[Atom]
            Original atoms from ligand A
        atoms_b : List[Atom]
            Original atoms from ligand B
        """
        # Calculate original total charges
        original_charge_a = sum(atom.charge for atom in atoms_a)
        original_charge_b = sum(atom.charge for atom in atoms_b)

        # Calculate current total charges from hybrid atoms
        current_charge_a = sum(atom.state_a_charge for atom in self.hybrid_atoms)
        current_charge_b = sum(
            atom.state_b_charge if atom.state_b_charge is not None else atom.state_a_charge
            for atom in self.hybrid_atoms
        )

        # Calculate calibration needed
        calibration_a = original_charge_a - current_charge_a
        calibration_b = original_charge_b - current_charge_b

        # Get atoms to modify based on charge_reception strategy
        modify_atoms = self._get_modify_list(mapping)

        # Distribute charge evenly to both states
        if modify_atoms and abs(calibration_a) > 1e-10:
            ave_a = calibration_a / len(modify_atoms)
            for atom in modify_atoms:
                atom.state_a_charge += ave_a

        if modify_atoms and abs(calibration_b) > 1e-10:
            ave_b = calibration_b / len(modify_atoms)
            for atom in modify_atoms:
                atom.state_b_charge += ave_b

    def _get_modify_list(self, mapping: AtomMapping) -> List[HybridAtom]:
        """
        Get list of atoms to modify for charge redistribution

        Following FEbuilder's logic:
        - Common atoms are EXCLUDED from redistribution (they stay the same)
        - Which atoms participate depends on charge_reception:
          * 'pert' (default): transformed + surrounding atoms (all uncommon)
          * 'unique': only non-hydrogen transformed atoms
          * 'surround': only surrounding atoms

        Parameters
        ----------
        mapping : AtomMapping
            Atom mapping classification

        Returns
        -------
        List[HybridAtom]
            List of atoms that should receive redistributed charge
        """
        modify_atoms = []

        # Count atoms in each category
        n_common = len(mapping.common)
        n_transformed_a = len(mapping.transformed_a)
        n_transformed_b = len(mapping.transformed_b)
        n_surrounding = len(mapping.surrounding_a)

        if self.charge_reception == "surround":
            # Only surrounding atoms participate
            start_idx = n_common + n_transformed_a + n_transformed_b
            end_idx = start_idx + n_surrounding
            atom_indices = range(start_idx, end_idx)

        elif self.charge_reception == "unique":
            # Only non-hydrogen transformed atoms participate
            # transformed atoms are at indices [n_common, n_common + n_transformed_a + n_transformed_b)
            start_idx = n_common
            end_idx = n_common + n_transformed_a + n_transformed_b
            for idx in range(start_idx, end_idx):
                atom = self.hybrid_atoms[idx]
                if atom.element != "H":
                    modify_atoms.append(atom)
            return modify_atoms  # Early return

        else:  # 'pert' (default)
            # All uncommon atoms participate (transformed + surrounding)
            start_idx = n_common
            end_idx = n_common + n_transformed_a + n_transformed_b + n_surrounding
            atom_indices = range(start_idx, end_idx)

        for atom in [self.hybrid_atoms[i] for i in atom_indices]:
            # Apply recharge_hydrogen filter
            if atom.element != "H" or self.recharge_hydrogen:
                modify_atoms.append(atom)

        return modify_atoms

    def _resolve_charge(self, charge_a: float, charge_b: float) -> float:
        """
        Resolve charge for common atoms

        Parameters
        ----------
        charge_a : float
            Charge from ligand A
        charge_b : float
            Charge from ligand B

        Returns
        -------
        float
            Resolved charge based on strategy
        """
        if self.charge_strategy == "ref":
            return charge_a
        elif self.charge_strategy == "mut":
            return charge_b
        else:  # mean
            return (charge_a + charge_b) / 2.0

    def _get_dummy_type(self, original_type: str) -> str:
        """
        Get dummy atom type for transformed atoms

        For GROMACS FEP, dummy atoms need special atom types.
        Convention: DUM_<original_type>

        Parameters
        ----------
        original_type : str
            Original atom type

        Returns
        -------
        str
            Dummy atom type
        """
        # Simple convention: DUM_<type>
        # In practice, this should match the force field's dummy type definitions
        return f"DUM_{original_type}"

    def _get_default_mass(self, element: str) -> float:
        """
        Get default atomic mass for element

        Parameters
        ----------
        element : str
            Element symbol (H, C, N, O, etc.)

        Returns
        -------
        float
            Atomic mass in amu
        """
        # Standard atomic masses
        masses = {
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
        }
        return masses.get(element, 12.011)  # Default to carbon if unknown
