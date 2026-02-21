"""Tests for prism.utils.ligand detection constants and logic."""

import pytest
from prism.utils.ligand import (
    COMMON_LIGAND_NAMES,
    STANDARD_RESIDUES,
)


class TestCommonLigandNames:
    """Test ligand name constants."""

    def test_lig_is_common(self):
        """LIG should be a common ligand name (PRISM default)."""
        assert 'LIG' in COMMON_LIGAND_NAMES

    def test_unl_is_common(self):
        """UNL (unknown ligand) should be recognized."""
        assert 'UNL' in COMMON_LIGAND_NAMES

    def test_mol_is_common(self):
        assert 'MOL' in COMMON_LIGAND_NAMES


class TestStandardResidueExclusion:
    """Test that standard residues are properly defined for exclusion."""

    def test_all_20_amino_acids_present(self):
        """All 20 standard amino acids should be in STANDARD_RESIDUES."""
        amino_acids = [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        ]
        for aa in amino_acids:
            assert aa in STANDARD_RESIDUES, f"{aa} missing from STANDARD_RESIDUES"

    def test_water_residues(self):
        """Water residue names should be standard."""
        for water in ['WAT', 'HOH', 'SOL']:
            assert water in STANDARD_RESIDUES, f"{water} missing from STANDARD_RESIDUES"

    def test_common_ions(self):
        """Common ion names should be standard."""
        for ion in ['NA', 'CL', 'K', 'MG', 'CA']:
            assert ion in STANDARD_RESIDUES, f"{ion} missing from STANDARD_RESIDUES"

    def test_ligand_not_standard(self):
        """Ligand names should NOT be in STANDARD_RESIDUES."""
        assert 'LIG' not in STANDARD_RESIDUES
        assert 'UNL' not in STANDARD_RESIDUES
        assert 'MOL' not in STANDARD_RESIDUES

    def test_nucleotides_present(self):
        """Nucleotide residues should be standard."""
        for nt in ['DA', 'DT', 'DG', 'DC']:
            assert nt in STANDARD_RESIDUES, f"{nt} missing from STANDARD_RESIDUES"
