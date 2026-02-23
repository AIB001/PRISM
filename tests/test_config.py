"""Tests for prism.utils.config.ConfigurationManager."""

import pytest
from prism.utils.config import ConfigurationManager


class TestConfigDefaults:
    """Test default configuration generation."""

    @pytest.fixture
    def cm(self):
        """Create a ConfigurationManager without running __init__."""
        obj = ConfigurationManager.__new__(ConfigurationManager)
        obj.gromacs_env = None  # _get_default_config checks this
        obj.forcefield_name = "amber99sb"
        obj.water_model_name = "tip3p"
        return obj

    def test_default_config_has_required_sections(self, cm):
        """Default config should contain all required sections."""
        defaults = cm._get_default_config()

        required_sections = [
            "general",
            "box",
            "simulation",
            "ions",
            "constraints",
            "energy_minimization",
            "output",
            "electrostatics",
            "vdw",
            "temperature_coupling",
            "pressure_coupling",
        ]
        for section in required_sections:
            assert section in defaults, f"Missing section: {section}"

    def test_default_temperature(self, cm):
        """Default temperature should be 310 K."""
        defaults = cm._get_default_config()
        assert defaults["simulation"]["temperature"] == 310

    def test_default_production_time(self, cm):
        """Default production time should be 500 ns."""
        defaults = cm._get_default_config()
        assert defaults["simulation"]["production_time_ns"] == 500

    def test_default_box_distance(self, cm):
        """Default box distance should be 1.5 nm."""
        defaults = cm._get_default_config()
        assert defaults["box"]["distance"] == 1.5

    def test_default_ion_concentration(self, cm):
        """Default ion concentration should be 0.15 M."""
        defaults = cm._get_default_config()
        assert defaults["ions"]["concentration"] == 0.15


class TestConfigMerge:
    """Test configuration merging logic."""

    def test_merge_preserves_defaults(self):
        """User config that doesn't override should preserve defaults."""
        cm = ConfigurationManager.__new__(ConfigurationManager)
        default = {"box": {"distance": 1.5, "shape": "cubic"}, "simulation": {"temperature": 310}}
        user = {"simulation": {"temperature": 300}}

        merged = cm._merge_configs(default, user)
        assert merged["box"]["distance"] == 1.5
        assert merged["box"]["shape"] == "cubic"
        assert merged["simulation"]["temperature"] == 300

    def test_merge_user_overrides(self):
        """User values should override defaults."""
        cm = ConfigurationManager.__new__(ConfigurationManager)
        default = {"simulation": {"temperature": 310, "dt": 0.002}}
        user = {"simulation": {"temperature": 300}}

        merged = cm._merge_configs(default, user)
        assert merged["simulation"]["temperature"] == 300
        assert merged["simulation"]["dt"] == 0.002

    def test_merge_nested_dicts(self):
        """Nested dicts should be recursively merged."""
        cm = ConfigurationManager.__new__(ConfigurationManager)
        default = {"a": {"b": {"c": 1, "d": 2}}}
        user = {"a": {"b": {"c": 10}}}

        merged = cm._merge_configs(default, user)
        assert merged["a"]["b"]["c"] == 10
        assert merged["a"]["b"]["d"] == 2

    def test_merge_new_keys_from_user(self):
        """User config can introduce new keys."""
        cm = ConfigurationManager.__new__(ConfigurationManager)
        default = {"a": 1}
        user = {"a": 1, "b": 2}

        merged = cm._merge_configs(default, user)
        assert merged["b"] == 2

    def test_merge_empty_user(self):
        """Empty user config should return defaults."""
        cm = ConfigurationManager.__new__(ConfigurationManager)
        default = {"a": 1, "b": {"c": 2}}
        user = {}

        merged = cm._merge_configs(default, user)
        assert merged == default
