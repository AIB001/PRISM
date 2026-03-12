#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GROMACS ITP File Builder

Generates GROMACS ITP files with typeB/chargeB encoding for FEP calculations.
"""

from typing import List, Dict
from ..core.dual_topology import HybridAtom


class ITPBuilder:
    """
    GROMACS ITP file generator for FEP calculations

    Generates hybrid topology files with A/B state parameters encoded
    in typeB/chargeB columns following GROMACS FEP format.

    Parameters
    ----------
    hybrid_atoms : List[HybridAtom]
        List of hybrid atoms with A/B state parameters
    hybrid_params : Dict[str, List]
        Hybrid parameters dictionary containing:
        - 'bonds': List of bond parameters
        - 'angles': List of angle parameters
        - 'dihedrals': List of dihedral parameters
        - 'impropers': List of improper parameters

    Notes
    -----
    GROMACS FEP format requirements:
    - [ atoms ] section needs typeB/chargeB columns (massB for pmx)
    - Common atoms: type = typeB, charge = chargeB (write once)
    - Transformed atoms: typeA ≠ typeB (write both states)
    - Virtual atoms: use dummy types (DUM_*) from forcefield
    """

    def __init__(self, hybrid_atoms: List[HybridAtom], hybrid_params: Dict[str, List]):
        self.hybrid_atoms = hybrid_atoms
        self.hybrid_params = hybrid_params

    def write_itp(self, output_path: str, molecule_name: str = "HYB"):
        """
        Write hybrid topology to GROMACS ITP file

        File structure:
        1. Header comments
        2. [ moleculetype ]
        3. [ atoms ] (with typeB/chargeB columns)
        4. [ bonds ]
        5. [ pairs ]
        6. [ angles ]
        7. [ dihedrals ]
        8. [ impropers ]

        Parameters
        ----------
        output_path : str
            Path to output ITP file
        molecule_name : str, optional
            Molecule name (default: "HYB")

        Raises
        ------
        NotImplementedError
            To be implemented in next phase
        """
        raise NotImplementedError("ITP file generation to be implemented")

    def _get_atoms_section(self) -> str:
        """
        Generate [ atoms ] section with typeB/chargeB columns

        This is the most critical section for GROMACS FEP format.

        Returns
        -------
        str
            Formatted [ atoms ] section

        Notes
        -----
        GROMACS format (important):
        ;   nr       type  resnr residue  atom   cgnr    charge   mass    typeB   chargeB   massB
        Common atoms: only write first 8 columns
        Transformed atoms: write all 10 columns
        """
        raise NotImplementedError("_get_atoms_section to be implemented")

    def _get_bonds_section(self) -> str:
        """Generate [ bonds ] section"""
        raise NotImplementedError("_get_bonds_section to be implemented")

    def _get_angles_section(self) -> str:
        """Generate [ angles ] section"""
        raise NotImplementedError("_get_angles_section to be implemented")

    def _get_dihedrals_section(self) -> str:
        """Generate [ dihedrals ] section"""
        raise NotImplementedError("_get_dihedrals_section to be implemented")

    def _get_impropers_section(self) -> str:
        """Generate [ impropers ] section"""
        raise NotImplementedError("_get_impropers_section to be implemented")
