"""Tests for prism.pmf.alignment helper functions and PMFAligner geometry."""

import pytest
import numpy as np
from prism.pmf.alignment import (
    _is_heavy,
    _spherical_to_cartesian,
    _cartesian_to_spherical,
    _compute_energy,
    PMFAligner,
)


class TestIsHeavy:
    """Test heavy atom detection."""

    def test_carbon_is_heavy(self):
        assert _is_heavy('C') is True

    def test_nitrogen_is_heavy(self):
        assert _is_heavy('N') is True

    def test_oxygen_is_heavy(self):
        assert _is_heavy('O') is True

    def test_hydrogen_not_heavy(self):
        assert _is_heavy('H') is False

    def test_deuterium_not_heavy(self):
        assert _is_heavy('D') is False

    def test_lowercase_hydrogen(self):
        assert _is_heavy('h') is False

    def test_sulfur_is_heavy(self):
        assert _is_heavy('S') is True

    def test_zinc_is_heavy(self):
        assert _is_heavy('Zn') is True


class TestSphericalCartesianConversion:
    """Test spherical <-> Cartesian coordinate conversion roundtrip."""

    def test_z_axis(self):
        """theta=0 should give [0, 0, 1] (Z-axis)."""
        v = _spherical_to_cartesian(0.0, 0.0)
        np.testing.assert_allclose(v, [0, 0, 1], atol=1e-10)

    def test_x_axis(self):
        """theta=pi/2, phi=0 should give [1, 0, 0] (X-axis)."""
        v = _spherical_to_cartesian(np.pi / 2, 0.0)
        np.testing.assert_allclose(v, [1, 0, 0], atol=1e-10)

    def test_y_axis(self):
        """theta=pi/2, phi=pi/2 should give [0, 1, 0] (Y-axis)."""
        v = _spherical_to_cartesian(np.pi / 2, np.pi / 2)
        np.testing.assert_allclose(v, [0, 1, 0], atol=1e-10)

    def test_negative_z_axis(self):
        """theta=pi should give [0, 0, -1] (-Z-axis)."""
        v = _spherical_to_cartesian(np.pi, 0.0)
        np.testing.assert_allclose(v, [0, 0, -1], atol=1e-10)

    def test_roundtrip_random_directions(self):
        """Converting to spherical and back should recover the original vector."""
        rng = np.random.RandomState(42)
        for _ in range(20):
            # Generate random unit vector
            v = rng.randn(3)
            v = v / np.linalg.norm(v)

            theta, phi = _cartesian_to_spherical(v)
            v_recovered = _spherical_to_cartesian(theta, phi)
            np.testing.assert_allclose(v_recovered, v, atol=1e-10)

    def test_unit_vector_output(self):
        """Output of _spherical_to_cartesian should always be unit length."""
        for theta in np.linspace(0, np.pi, 10):
            for phi in np.linspace(0, 2 * np.pi, 10):
                v = _spherical_to_cartesian(theta, phi)
                assert np.linalg.norm(v) == pytest.approx(1.0, abs=1e-10)


class TestComputeEnergy:
    """Test PMF pulling direction energy computation."""

    def test_energy_is_negative(self):
        """Energy should be negative (it's -sum of distances)."""
        lig_coords = np.array([[0.0, 0.0, 0.0]])
        pocket_coords = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        direction = np.array([0.0, 0.0, 1.0])
        energy = _compute_energy(lig_coords, pocket_coords, direction)
        assert energy < 0

    def test_energy_varies_with_direction(self):
        """Energy should change with different pulling directions."""
        lig_coords = np.array([[0.0, 0.0, 0.0]])
        pocket_coords = np.array([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0]])

        e_z = _compute_energy(lig_coords, pocket_coords, np.array([0.0, 0.0, 1.0]))
        e_x = _compute_energy(lig_coords, pocket_coords, np.array([1.0, 0.0, 0.0]))
        # Pulling along Z should have different clearance than along X
        assert e_z != e_x


class TestRotationMatrix:
    """Test PMFAligner._calculate_rotation_matrix."""

    @pytest.fixture
    def aligner(self):
        return PMFAligner(pocket_cutoff=4.0, verbose=False)

    def test_identity_when_aligned(self, aligner):
        """Same vector should give identity matrix."""
        v = np.array([0.0, 0.0, 1.0])
        R = aligner._calculate_rotation_matrix(v, v)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-10)

    def test_90_degree_rotation(self, aligner):
        """X-axis to Z-axis should be a 90-degree rotation."""
        v_from = np.array([1.0, 0.0, 0.0])
        v_to = np.array([0.0, 0.0, 1.0])
        R = aligner._calculate_rotation_matrix(v_from, v_to)

        # Apply rotation to v_from
        rotated = R @ v_from
        np.testing.assert_allclose(rotated, v_to, atol=1e-10)

    def test_180_degree_rotation(self, aligner):
        """Opposite vectors should produce a valid 180-degree rotation."""
        v_from = np.array([0.0, 0.0, 1.0])
        v_to = np.array([0.0, 0.0, -1.0])
        R = aligner._calculate_rotation_matrix(v_from, v_to)

        rotated = R @ v_from
        np.testing.assert_allclose(rotated, v_to, atol=1e-10)

    def test_rotation_is_orthogonal(self, aligner):
        """Rotation matrix should be orthogonal (R^T R = I)."""
        v_from = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
        v_to = np.array([0.0, 0.0, 1.0])
        R = aligner._calculate_rotation_matrix(v_from, v_to)

        np.testing.assert_allclose(R.T @ R, np.eye(3), atol=1e-10)
        np.testing.assert_allclose(np.linalg.det(R), 1.0, atol=1e-10)

    def test_arbitrary_rotation(self, aligner):
        """Arbitrary vectors should produce correct rotation."""
        v_from = np.array([1.0, 2.0, 3.0])
        v_from = v_from / np.linalg.norm(v_from)
        v_to = np.array([0.0, 0.0, 1.0])

        R = aligner._calculate_rotation_matrix(v_from, v_to)
        rotated = R @ v_from
        np.testing.assert_allclose(rotated, v_to, atol=1e-10)
