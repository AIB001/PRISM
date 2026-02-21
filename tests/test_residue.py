"""Tests for prism.utils.residue amino acid code conversion."""

import pytest
from prism.utils.residue import (
    convert_residue_to_1letter,
    convert_residue_to_3letter,
    AA_3TO1,
    AA_1TO3,
)


class TestThreeToOneLetter:
    """Test 3-letter to 1-letter residue conversion."""

    def test_asp(self):
        assert convert_residue_to_1letter("ASP618") == "D618"

    def test_ser(self):
        assert convert_residue_to_1letter("SER759") == "S759"

    def test_glu(self):
        assert convert_residue_to_1letter("GLU123") == "E123"

    def test_ala(self):
        assert convert_residue_to_1letter("ALA1") == "A1"

    def test_already_one_letter(self):
        """Already 1-letter format should pass through unchanged."""
        assert convert_residue_to_1letter("D618") == "D618"

    def test_unknown_residue(self):
        """Unknown 3-letter codes should pass through unchanged."""
        assert convert_residue_to_1letter("XYZ123") == "XYZ123"

    def test_no_number(self):
        """3-letter code without number should pass through unchanged."""
        assert convert_residue_to_1letter("ALA") == "ALA"

    def test_all_standard_amino_acids(self):
        """All 20 standard amino acids should convert correctly."""
        for three, one in AA_3TO1.items():
            assert convert_residue_to_1letter(f"{three}42") == f"{one}42"


class TestOneToThreeLetter:
    """Test 1-letter to 3-letter residue conversion."""

    def test_d_to_asp(self):
        assert convert_residue_to_3letter("D618") == "ASP618"

    def test_s_to_ser(self):
        assert convert_residue_to_3letter("S759") == "SER759"

    def test_e_to_glu(self):
        assert convert_residue_to_3letter("E123") == "GLU123"

    def test_already_three_letter(self):
        """Already 3-letter format should pass through unchanged."""
        assert convert_residue_to_3letter("ASP618") == "ASP618"

    def test_unknown_letter(self):
        """Unknown 1-letter code should pass through unchanged."""
        assert convert_residue_to_3letter("X123") == "X123"

    def test_no_number(self):
        """Single letter without number should pass through unchanged."""
        assert convert_residue_to_3letter("D") == "D"

    def test_all_standard_amino_acids(self):
        """All 20 standard amino acids should convert correctly."""
        for one, three in AA_1TO3.items():
            assert convert_residue_to_3letter(f"{one}42") == f"{three}42"


class TestRoundtrip:
    """Test roundtrip conversion."""

    def test_3_to_1_to_3(self):
        """3-letter -> 1-letter -> 3-letter should recover original."""
        original = "ASP618"
        one_letter = convert_residue_to_1letter(original)
        recovered = convert_residue_to_3letter(one_letter)
        assert recovered == original

    def test_1_to_3_to_1(self):
        """1-letter -> 3-letter -> 1-letter should recover original."""
        original = "D618"
        three_letter = convert_residue_to_3letter(original)
        recovered = convert_residue_to_1letter(three_letter)
        assert recovered == original
