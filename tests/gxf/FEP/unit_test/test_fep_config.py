#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unit tests for FEPConfig class
"""

import pytest
import tempfile
from pathlib import Path

from prism.fep.config import FEPConfig


class TestFEPConfig:
    """Test FEPConfig parameter loading and defaults"""

    def test_simulation_params_defaults(self):
        """Test simulation parameters use correct defaults"""
        # Create minimal fep.yaml
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("mapping:\n  charge_common: mean\n")

            config = FEPConfig(tmpdir)
            params = config.get_simulation_params()

            assert params["equilibration_nvt_time_ps"] == 500
            assert params["equilibration_npt_time_ps"] == 500
            assert params["production_time_ns"] == 5.0
            assert params["dt"] == 0.002
            assert params["temperature"] == 310
            assert params["pressure"] == 1.0

    def test_soft_core_params_defaults(self):
        """Test soft-core parameters use correct defaults"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("mapping:\n  charge_common: mean\n")

            config = FEPConfig(tmpdir)
            params = config.get_soft_core_params()

            assert params["alpha"] == 0.5
            assert params["sigma"] == 0.3

    def test_lambda_params_defaults(self):
        """Test lambda parameters use correct defaults"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("mapping:\n  charge_common: mean\n")

            config = FEPConfig(tmpdir)
            params = config.get_lambda_params()

            assert params["windows"] == 32

    def test_electrostatics_params_defaults(self):
        """Test electrostatics parameters use correct defaults"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("mapping:\n  charge_common: mean\n")

            config = FEPConfig(tmpdir)
            params = config.get_electrostatics_params()

            assert params["rcoulomb"] == 1.0

    def test_vdw_params_defaults(self):
        """Test VDW parameters use correct defaults"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("mapping:\n  charge_common: mean\n")

            config = FEPConfig(tmpdir)
            params = config.get_vdw_params()

            assert params["rvdw"] == 1.0

    def test_output_params_defaults(self):
        """Test output parameters use correct defaults"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("mapping:\n  charge_common: mean\n")

            config = FEPConfig(tmpdir)
            params = config.get_output_params()

            assert params["trajectory_interval_ps"] == 500
            assert params["energy_interval_ps"] == 10
            assert params["log_interval_ps"] == 10

    def test_get_all_mdp_params(self):
        """Test get_all_mdp_params returns correct structure"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("mapping:\n  charge_common: mean\n")

            config = FEPConfig(tmpdir)
            params = config.get_all_mdp_params()

            # Check structure
            assert "simulation" in params
            assert "electrostatics" in params
            assert "vdw" in params
            assert "output" in params

            # Check nested structure
            assert "temperature" in params["simulation"]
            assert "rcoulomb" in params["electrostatics"]
            assert "rvdw" in params["vdw"]
            assert "trajectory_interval_ps" in params["output"]

    def test_custom_simulation_params(self):
        """Test custom simulation parameters override defaults"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("""
mapping:
  charge_common: mean
simulation:
  temperature: 320
  production_time_ns: 10.0
""")

            config = FEPConfig(tmpdir)
            params = config.get_simulation_params()

            assert params["temperature"] == 320
            assert params["production_time_ns"] == 10.0
            # Other params should still use defaults
            assert params["pressure"] == 1.0
            assert params["dt"] == 0.002

    def test_custom_lambda_windows(self):
        """Test custom lambda windows parameter"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("""
mapping:
  charge_common: mean
lambda:
  windows: 21
""")

            config = FEPConfig(tmpdir)
            params = config.get_lambda_params()

            assert params["windows"] == 21

    def test_real_fep_yaml_files(self):
        """Test loading from real fep.yaml files"""
        # Test oMeEtPh-EtPh configuration
        config = FEPConfig("tests/gxf/FEP/unit_test/oMeEtPh-EtPh")

        # Should load all parameter sections
        sim = config.get_simulation_params()
        sc = config.get_soft_core_params()
        lam = config.get_lambda_params()

        # Check values from the actual file
        assert sim["temperature"] == 310
        assert sc["alpha"] == 0.5
        assert lam["windows"] == 32

    def test_backward_compatibility(self):
        """Test that configs without new sections still work"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create an old-style fep.yaml with only mapping section
            fep_file = Path(tmpdir) / "fep.yaml"
            fep_file.write_text("""
mapping:
  dist_cutoff: 0.6
  charge_cutoff: 0.05
  charge_common: mean
  charge_reception: surround
""")

            config = FEPConfig(tmpdir)

            # All getters should work with defaults
            config.get_simulation_params()
            config.get_soft_core_params()
            config.get_lambda_params()
            config.get_electrostatics_params()
            config.get_vdw_params()
            config.get_output_params()
            config.get_all_mdp_params()

            # Should not raise any exceptions


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
