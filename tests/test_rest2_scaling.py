"""Tests for prism.rest2.partial_tempering scaling functions."""

import pytest
import math
from prism.rest2.partial_tempering import (
    _is_solute_type,
    _classify_interaction,
    _get_scale_factor,
    _strip_comment,
    _looks_numeric,
)


class TestIsSoluteType:
    """Test atom type solute classification."""

    def test_solute_type_with_underscore(self):
        """Atom types ending in '_' are solute."""
        assert _is_solute_type("CT_") is True

    def test_solute_type_double_underscore(self):
        assert _is_solute_type("CA__") is True

    def test_solvent_type_no_underscore(self):
        """Atom types without trailing '_' are solvent."""
        assert _is_solute_type("CT") is False

    def test_solvent_type_middle_underscore(self):
        """Underscore in middle doesn't make it solute."""
        assert _is_solute_type("C_T") is False

    def test_empty_string(self):
        assert _is_solute_type("") is False

    def test_single_underscore(self):
        assert _is_solute_type("_") is True


class TestClassifyInteraction:
    """Test bonded interaction classification."""

    def test_all_solute(self):
        """All atoms are solute -> all_solute."""
        atom_is_solute = {1: True, 2: True, 3: True}
        result = _classify_interaction([1, 2, 3], atom_is_solute)
        assert result == "all_solute"

    def test_all_solvent(self):
        """No atoms are solute -> all_solvent."""
        atom_is_solute = {1: False, 2: False}
        result = _classify_interaction([1, 2], atom_is_solute)
        assert result == "all_solvent"

    def test_mixed(self):
        """Some solute, some solvent -> mixed."""
        atom_is_solute = {1: True, 2: False}
        result = _classify_interaction([1, 2], atom_is_solute)
        assert result == "mixed"

    def test_single_solute_atom(self):
        """Single solute atom -> all_solute."""
        atom_is_solute = {5: True}
        result = _classify_interaction([5], atom_is_solute)
        assert result == "all_solute"

    def test_unknown_atom_treated_as_solvent(self):
        """Atoms not in dict default to False (solvent)."""
        atom_is_solute = {1: True}
        result = _classify_interaction([1, 99], atom_is_solute)
        assert result == "mixed"


class TestGetScaleFactor:
    """Test REST2 scale factor selection."""

    def test_all_solute_gets_lambda(self):
        """All-solute interactions scale by lambda."""
        lam = 0.7
        result = _get_scale_factor("all_solute", lam, math.sqrt(lam))
        assert result == pytest.approx(lam)

    def test_mixed_gets_sqrt_lambda(self):
        """Mixed interactions scale by sqrt(lambda)."""
        lam = 0.7
        sqrt_lam = math.sqrt(lam)
        result = _get_scale_factor("mixed", lam, sqrt_lam)
        assert result == pytest.approx(sqrt_lam)

    def test_all_solvent_gets_one(self):
        """All-solvent interactions are unscaled (factor = 1.0)."""
        lam = 0.7
        result = _get_scale_factor("all_solvent", lam, math.sqrt(lam))
        assert result == pytest.approx(1.0)

    def test_lambda_one_gives_no_scaling(self):
        """When lambda=1.0, all scale factors should be 1.0."""
        lam = 1.0
        assert _get_scale_factor("all_solute", lam, math.sqrt(lam)) == pytest.approx(1.0)
        assert _get_scale_factor("mixed", lam, math.sqrt(lam)) == pytest.approx(1.0)
        assert _get_scale_factor("all_solvent", lam, math.sqrt(lam)) == pytest.approx(1.0)

    def test_charge_scaling(self):
        """Charges scale by sqrt(lambda): charge * sqrt(lambda)."""
        lam = 0.5
        original_charge = 0.42
        scaled_charge = original_charge * math.sqrt(lam)
        assert scaled_charge == pytest.approx(original_charge * math.sqrt(0.5))


class TestStripComment:
    """Test line comment stripping."""

    def test_line_with_semicolon_comment(self):
        data, comment = _strip_comment("  1  CT  0.1234  ; some comment")
        assert "some comment" in comment
        assert "1" in data

    def test_line_without_comment(self):
        data, comment = _strip_comment("  1  CT  0.1234")
        assert comment == ""
        assert "CT" in data

    def test_empty_line(self):
        data, comment = _strip_comment("")
        assert data == ""


class TestLooksNumeric:
    """Test numeric detection."""

    def test_integer(self):
        assert _looks_numeric("42") is True

    def test_float(self):
        assert _looks_numeric("3.14") is True

    def test_scientific(self):
        assert _looks_numeric("1.23e-05") is True

    def test_negative(self):
        assert _looks_numeric("-2.5") is True

    def test_word(self):
        assert _looks_numeric("hello") is False

    def test_atom_type(self):
        assert _looks_numeric("CT") is False
