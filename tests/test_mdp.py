"""Tests for prism.utils.mdp.MDPGenerator parameter calculations."""

import os
import pytest
import tempfile
from prism.utils.mdp import MDPGenerator


@pytest.fixture
def mdp_gen(sample_config, tmp_path):
    """Create an MDPGenerator with sample config and temp directory."""
    return MDPGenerator(sample_config, str(tmp_path))


class TestMDPGeneration:
    """Test MDP file generation and parameter calculations."""

    def test_mdp_dir_created(self, mdp_gen, tmp_path):
        """MDP directory should be created during generate_all."""
        mdp_gen.generate_all()
        mdp_dir = os.path.join(str(tmp_path), 'mdps')
        assert os.path.isdir(mdp_dir)

    def test_em_mdp_created(self, mdp_gen, tmp_path):
        """Energy minimization MDP should be generated."""
        mdp_gen.generate_all()
        em_mdp = os.path.join(str(tmp_path), 'mdps', 'em.mdp')
        assert os.path.isfile(em_mdp)

    def test_nvt_mdp_created(self, mdp_gen, tmp_path):
        """NVT equilibration MDP should be generated."""
        mdp_gen.generate_all()
        nvt_mdp = os.path.join(str(tmp_path), 'mdps', 'nvt.mdp')
        assert os.path.isfile(nvt_mdp)

    def test_npt_mdp_created(self, mdp_gen, tmp_path):
        """NPT equilibration MDP should be generated."""
        mdp_gen.generate_all()
        npt_mdp = os.path.join(str(tmp_path), 'mdps', 'npt.mdp')
        assert os.path.isfile(npt_mdp)

    def test_production_mdp_created(self, mdp_gen, tmp_path):
        """Production MD MDP should be generated."""
        mdp_gen.generate_all()
        md_mdp = os.path.join(str(tmp_path), 'mdps', 'md.mdp')
        assert os.path.isfile(md_mdp)


class TestMDPParameterCalculations:
    """Test step count and interval calculations in generated MDP files."""

    def test_production_nsteps(self, mdp_gen, tmp_path):
        """Production MD nsteps = production_time_ns * 1000 / dt."""
        mdp_gen.generate_all()
        md_mdp = os.path.join(str(tmp_path), 'mdps', 'md.mdp')

        with open(md_mdp, 'r') as f:
            content = f.read()

        # 500 ns / 0.002 ps = 250,000,000 steps
        expected_steps = int(500 * 1000 / 0.002)
        assert f"nsteps" in content
        # Find the nsteps value
        for line in content.split('\n'):
            if line.strip().startswith('nsteps'):
                parts = line.split('=')
                if len(parts) == 2:
                    val = int(parts[1].strip().split()[0])
                    assert val == expected_steps
                    break

    def test_nvt_temperature(self, mdp_gen, tmp_path):
        """NVT MDP should contain correct temperature."""
        mdp_gen.generate_all()
        nvt_mdp = os.path.join(str(tmp_path), 'mdps', 'nvt.mdp')

        with open(nvt_mdp, 'r') as f:
            content = f.read()

        assert '310' in content  # Temperature from config

    def test_em_integrator(self, mdp_gen, tmp_path):
        """EM MDP should use steep integrator."""
        mdp_gen.generate_all()
        em_mdp = os.path.join(str(tmp_path), 'mdps', 'em.mdp')

        with open(em_mdp, 'r') as f:
            content = f.read()

        assert 'steep' in content

    def test_dt_value(self, mdp_gen, tmp_path):
        """MDP files should contain correct timestep."""
        mdp_gen.generate_all()
        nvt_mdp = os.path.join(str(tmp_path), 'mdps', 'nvt.mdp')

        with open(nvt_mdp, 'r') as f:
            content = f.read()

        assert '0.002' in content
