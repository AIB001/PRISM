#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dual Topology Construction Module

Builds hybrid topologies with A/B state parameters for GROMACS FEP calculations.
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


class DualTopologyBuilder:
    """
    Builds dual-topology hybrid systems for FEP calculations

    Creates GROMACS ITP files with typeB/chargeB columns encoding
    both reference and mutant states.

    Parameters
    ----------
    charge_strategy : str, optional
        Strategy for handling common atom charges (default: 'mean')
        - 'mean': Use average of reference and mutant charges
        - 'ref': Use reference ligand charges
        - 'mut': Use mutant ligand charges

    Notes
    -----
    Key design decision: GROMACS does NOT need parameter value merging
    like NAMD/CHARMM. We only need to map atom indices - GROMACS
    automatically selects A or B state parameters based on typeB/chargeB.
    """

    def __init__(self, charge_strategy: str = "mean"):
        if charge_strategy not in ["ref", "mut", "mean"]:
            raise ValueError(f"Invalid charge_strategy: {charge_strategy}")
        self.charge_strategy = charge_strategy

    def build(self, mapping: AtomMapping, params_a: Dict, params_b: Dict) -> List[HybridAtom]:
        """
        Build dual-topology hybrid atoms

        Parameters
        ----------
        mapping : AtomMapping
            Atom mapping classification
        params_a : Dict
            Force field parameters for ligand A
            Contains 'masses', 'bonds', 'angles', 'dihedrals', etc.
        params_b : Dict
            Force field parameters for ligand B

        Returns
        -------
        List[HybridAtom]
            List of hybrid atoms with typeB/chargeB encoding

        Raises
        ------
        NotImplementedError
            To be implemented in next phase (FEbuilder port)
        """
        raise NotImplementedError("Dual topology construction to be ported from FEbuilder/setup/make_hybrid.py")

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
