"""Tests for prism.utils.cleaner.ProteinCleaner."""

import pytest
import numpy as np
from prism.utils.cleaner import ProteinCleaner


@pytest.fixture
def cleaner():
    """Create a ProteinCleaner with default settings."""
    return ProteinCleaner(verbose=False)


@pytest.fixture
def cleaner_keep_all():
    """Create a ProteinCleaner in keep_all mode."""
    return ProteinCleaner(ion_mode="keep_all", verbose=False)


@pytest.fixture
def cleaner_remove_all():
    """Create a ProteinCleaner in remove_all mode."""
    return ProteinCleaner(ion_mode="remove_all", verbose=False)


class TestIsMetalOrIon:
    """Test metal/ion identification."""

    def test_zinc_detected(self, cleaner):
        assert cleaner._is_metal_or_ion("ZN", "ZN") is True

    def test_magnesium_detected(self, cleaner):
        assert cleaner._is_metal_or_ion("MG", "MG") is True

    def test_calcium_detected(self, cleaner):
        assert cleaner._is_metal_or_ion("CA", "CA") is True

    def test_iron_detected(self, cleaner):
        assert cleaner._is_metal_or_ion("FE", "FE") is True

    def test_sodium_detected(self, cleaner):
        assert cleaner._is_metal_or_ion("NA", "NA") is True

    def test_chloride_detected(self, cleaner):
        assert cleaner._is_metal_or_ion("CL", "CL") is True

    def test_regular_residue_not_metal(self, cleaner):
        """Normal amino acid atoms (not CA/ambiguous) should not be metal."""
        assert cleaner._is_metal_or_ion("ALA", "N") is False
        assert cleaner._is_metal_or_ion("ALA", "C") is False
        assert cleaner._is_metal_or_ion("ALA", "O") is False

    def test_water_not_metal(self, cleaner):
        assert cleaner._is_metal_or_ion("HOH", "O") is False


class TestShouldKeepMetalSmart:
    """Test smart metal retention logic."""

    def test_keep_structural_zinc(self, cleaner):
        """Zinc is structural and should be kept."""
        assert cleaner._should_keep_metal_smart("ZN", "ZN") is True

    def test_keep_structural_magnesium(self, cleaner):
        """Magnesium is structural and should be kept."""
        assert cleaner._should_keep_metal_smart("MG", "MG") is True

    def test_keep_structural_calcium(self, cleaner):
        """Calcium is structural and should be kept."""
        assert cleaner._should_keep_metal_smart("CA", "CA") is True

    def test_remove_sodium(self, cleaner):
        """Sodium is non-structural and should be removed."""
        assert cleaner._should_keep_metal_smart("NA", "NA") is False

    def test_remove_chloride(self, cleaner):
        """Chloride is non-structural and should be removed."""
        assert cleaner._should_keep_metal_smart("CL", "CL") is False

    def test_remove_potassium(self, cleaner):
        """Potassium is non-structural and should be removed."""
        assert cleaner._should_keep_metal_smart("K", "K") is False


class TestCrystallizationArtifacts:
    """Test crystallization artifact identification."""

    def test_glycerol_is_artifact(self):
        assert "GOL" in ProteinCleaner.CRYSTALLIZATION_ARTIFACTS

    def test_ethylene_glycol_is_artifact(self):
        assert "EDO" in ProteinCleaner.CRYSTALLIZATION_ARTIFACTS

    def test_peg_is_artifact(self):
        assert "PEG" in ProteinCleaner.CRYSTALLIZATION_ARTIFACTS

    def test_nag_is_artifact(self):
        assert "NAG" in ProteinCleaner.CRYSTALLIZATION_ARTIFACTS

    def test_standard_residue_not_artifact(self):
        assert "ALA" not in ProteinCleaner.CRYSTALLIZATION_ARTIFACTS
        assert "GLY" not in ProteinCleaner.CRYSTALLIZATION_ARTIFACTS


class TestWaterIdentification:
    """Test water residue identification."""

    def test_hoh_is_water(self):
        assert "HOH" in ProteinCleaner.WATER_NAMES

    def test_wat_is_water(self):
        assert "WAT" in ProteinCleaner.WATER_NAMES

    def test_sol_is_water(self):
        assert "SOL" in ProteinCleaner.WATER_NAMES


class TestMinDistance:
    """Test minimum distance calculation."""

    def test_zero_distance(self, cleaner):
        """Point at origin, coords include origin."""
        point = np.array([0.0, 0.0, 0.0])
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        assert cleaner._calculate_min_distance(point, coords) == pytest.approx(0.0)

    def test_known_distance(self, cleaner):
        """Known distance between points."""
        point = np.array([3.0, 4.0, 0.0])
        coords = np.array([[0.0, 0.0, 0.0]])
        assert cleaner._calculate_min_distance(point, coords) == pytest.approx(5.0)

    def test_multiple_coords_returns_min(self, cleaner):
        """Should return minimum distance among multiple coordinates."""
        point = np.array([1.0, 0.0, 0.0])
        coords = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
        assert cleaner._calculate_min_distance(point, coords) == pytest.approx(0.5)
