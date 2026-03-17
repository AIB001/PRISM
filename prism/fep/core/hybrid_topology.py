#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Hybrid Topology Construction Module

Builds hybrid topologies with A/B state parameters for GROMACS FEP calculations.
GROMACS uses single-topology approach with typeB/chargeB encoding.
"""

from dataclasses import dataclass
from typing import List, Dict, Optional
from .mapping import AtomMapping, normalize_charge_reception


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
        - 'none': Keep original charges (no modification)

    charge_reception : str, optional
        Strategy for charge redistribution (default: 'surround')
        - 'surround': Only surrounding atoms receive redistributed charge
        - 'unique': Surrounding + transformed atoms (state-specific)
        - 'none': No charge redistribution

    recharge_hydrogen : bool, optional
        Whether to include hydrogen atoms in charge redistribution (default: False)

    Notes
    -----
    GROMACS uses single-topology approach: one topology with typeB/chargeB
    columns to encode both states. This is different from NAMD's dual-topology
    approach which requires two separate topologies.

    Charge redistribution is necessary after applying charge_strategy to common atoms
    to maintain system neutrality. The redistribution is state-specific:
    - State A: Only affects atoms that exist in state A (excludes transformed B)
    - State B: Only affects atoms that exist in state B (excludes transformed A)
    """

    def __init__(
        self, charge_strategy: str = "mean", charge_reception: str = "surround", recharge_hydrogen: bool = False
    ):
        if charge_strategy not in ["ref", "mut", "mean", "none"]:
            raise ValueError(f"Invalid charge_strategy: {charge_strategy}")
        charge_reception = normalize_charge_reception(charge_reception)
        if charge_reception not in ["unique", "surround", "surround_ext", "none"]:
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
            if self.charge_strategy == "none":
                # Keep original charges (different for A/B)
                # MUST have typeB/chargeB because charges are different
                self.hybrid_atoms.append(
                    HybridAtom(
                        name=atom_a.name,
                        index=index,
                        state_a_type=atom_a.atom_type,
                        state_a_charge=atom_a.charge,  # Original A charge
                        state_b_type=atom_b.atom_type,  # MUST have typeB
                        state_b_charge=atom_b.charge,  # Original B charge (different!)
                        element=atom_a.element,
                        mass=masses_a.get(atom_a.atom_type, self._get_default_mass(atom_a.element)),
                    )
                )
            else:
                # Use resolved charge (same for A/B)
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

        # Apply charge redistribution to maintain neutrality
        # Skip for 'none' mode
        if self.charge_strategy != "none" and self.charge_reception != "none":
            self._redistribute_charges(mapping, atoms_a, atoms_b)

        return self.hybrid_atoms

    def _redistribute_charges(self, mapping: AtomMapping, atoms_a: List, atoms_b: List) -> None:
        """
        Redistribute charges to maintain system neutrality

        After applying charge_strategy to common atoms, the total system charge
        changes. This method redistributes the charge difference to selected atoms
        according to charge_reception strategy.

        CRITICAL: For GROMACS single-topology FEP:
        - State A (lambda=0): Uses state_a_charge for all atoms (transformed B atoms are dummy, charge=0)
        - State B (lambda=1): Uses state_b_charge for all atoms (transformed A atoms are dummy, charge=0)
        - So state A should match ligand A's total charge
        - And state B should match ligand B's total charge

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
        # State A: sum of all state_a_charge (transformed B atoms contribute 0)
        current_charge_a = sum(atom.state_a_charge for atom in self.hybrid_atoms)

        # State B: sum of state_b_charge (or state_a_charge if state_b is None)
        # For transformed A atoms, state_b_charge should be 0 (dummy)
        current_charge_b = sum(
            atom.state_b_charge if atom.state_b_charge is not None else atom.state_a_charge
            for atom in self.hybrid_atoms
        )

        # Calculate calibration needed
        calibration_a = original_charge_a - current_charge_a
        calibration_b = original_charge_b - current_charge_b

        # Get separate lists for state A and state B redistribution
        # State A: transformed A atoms can participate (they exist in state A)
        # State B: transformed B atoms can participate (they exist in state B)
        modify_atoms_a = self._get_modify_list_for_state(mapping, "A")
        modify_atoms_b = self._get_modify_list_for_state(mapping, "B")

        # Distribute charge to state A
        if modify_atoms_a and abs(calibration_a) > 1e-10:
            ave_a = calibration_a / len(modify_atoms_a)

            # surround_ext: if charge per atom > 0.02, extend to common atoms
            if self.charge_reception == "surround_ext" and abs(ave_a) >= 0.02:
                # Add common atoms to modify list
                common_atoms = self.hybrid_atoms[: len(mapping.common)]
                modify_atoms_a.extend(common_atoms)
                ave_a = calibration_a / len(modify_atoms_a)

            for atom in modify_atoms_a:
                atom.state_a_charge += ave_a

        # Distribute charge to state B
        if modify_atoms_b and abs(calibration_b) > 1e-10:
            ave_b = calibration_b / len(modify_atoms_b)

            # surround_ext: if charge per atom > 0.02, extend to common atoms
            if self.charge_reception == "surround_ext" and abs(ave_b) >= 0.02:
                # Add common atoms to modify list
                common_atoms = self.hybrid_atoms[: len(mapping.common)]
                modify_atoms_b.extend(common_atoms)
                ave_b = calibration_b / len(modify_atoms_b)

            for atom in modify_atoms_b:
                if atom.state_b_charge is not None:
                    atom.state_b_charge += ave_b

    def _get_modify_list_for_state(self, mapping: AtomMapping, state: str) -> List[HybridAtom]:
        """
        Get list of atoms to modify for charge redistribution (state-specific)

        For GROMACS single-topology FEP, we need to handle state A and state B separately:
        - State A: transformed A atoms exist, transformed B atoms are dummy (exclude)
        - State B: transformed B atoms exist, transformed A atoms are dummy (exclude)
        - Common atoms: EXCLUDED from redistribution (already handled by charge_strategy)
        - Surrounding atoms: always participate when available

        CRITICAL: If no surrounding atoms exist, fallback to transformed atoms (non-dummy state only).
        If still empty after applying recharge_hydrogen filter, ignore the filter and use all available atoms.

        Parameters
        ----------
        mapping : AtomMapping
            Atom mapping classification
        state : str
            Either 'A' or 'B', specifies which state to redistribute

        Returns
        -------
        List[HybridAtom]
            List of atoms that should receive redistributed charge for this state
        """
        modify_atoms = []

        # Count atoms in each category
        n_common = len(mapping.common)
        n_transformed_a = len(mapping.transformed_a)
        n_transformed_b = len(mapping.transformed_b)
        n_surrounding = len(mapping.surrounding_a)

        # Indices for different atom categories in hybrid_atoms list
        # Order: common, transformed_a, transformed_b, surrounding
        common_end = n_common
        transformed_a_end = n_common + n_transformed_a
        transformed_b_end = transformed_a_end + n_transformed_b
        surrounding_end = transformed_b_end + n_surrounding

        if self.charge_reception == "none":
            return []

        # Build atom indices list based on charge_reception mode
        atom_indices = []

        if self.charge_reception in ["surround", "surround_ext"]:
            # Only surrounding atoms participate (initially)
            # surround_ext will extend to common atoms if needed (in _redistribute_charges)
            # If no surrounding atoms, fallback to transformed atoms (non-dummy state)
            if n_surrounding > 0:
                atom_indices = range(transformed_b_end, surrounding_end)
            else:
                # Fallback: use transformed atoms (state-specific)
                if state == "A":
                    atom_indices = range(common_end, transformed_a_end)  # transformed A
                else:  # state == 'B'
                    atom_indices = range(transformed_a_end, transformed_b_end)  # transformed B

        elif self.charge_reception == "unique":
            # Surrounding + transformed atoms (state-specific)
            if state == "A":
                # Surrounding + transformed A (exclude transformed B because they're dummy in state A)
                atom_indices = list(range(common_end, transformed_a_end))
                atom_indices.extend(range(transformed_b_end, surrounding_end))
            else:  # state == 'B'
                # Surrounding + transformed B (exclude transformed A because they're dummy in state B)
                atom_indices = list(range(transformed_a_end, transformed_b_end))
                atom_indices.extend(range(transformed_b_end, surrounding_end))

        else:
            raise ValueError(f"Invalid charge_reception: {self.charge_reception}")

        for atom in [self.hybrid_atoms[i] for i in atom_indices]:
            # Apply recharge_hydrogen filter
            if atom.element != "H" or self.recharge_hydrogen:
                modify_atoms.append(atom)

        # CRITICAL: If modify list is empty (e.g., all are H and recharge_hydrogen=False),
        # fallback to include all atoms from atom_indices (ignore recharge_hydrogen)
        if not modify_atoms and atom_indices:
            for atom in [self.hybrid_atoms[i] for i in atom_indices]:
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
