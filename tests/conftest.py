"""Shared fixtures and test data paths for PRISM unit tests."""

import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    """Path to test data directory."""
    return Path(__file__).parent.parent / "test" / "4xb4"


@pytest.fixture
def sample_config():
    """Minimal PRISM configuration dict for testing."""
    return {
        'general': {'overwrite': False},
        'forcefield': {
            'index': 1,
            'custom_forcefields': {1: {'name': 'amber99sb', 'dir': 'amber99sb.ff', 'path': None}},
        },
        'water_model': {
            'index': 1,
            'custom_water_models': {1: {'name': 'tip3p'}},
        },
        'simulation': {
            'temperature': 310,
            'pressure': 1.0,
            'pH': 7.0,
            'ligand_charge': 0,
            'production_time_ns': 500,
            'dt': 0.002,
            'equilibration_nvt_time_ps': 500,
            'equilibration_npt_time_ps': 500,
        },
        'box': {
            'distance': 1.5,
            'shape': 'cubic',
            'center': True,
        },
        'ions': {
            'neutral': True,
            'concentration': 0.15,
            'positive_ion': 'NA',
            'negative_ion': 'CL',
        },
        'constraints': {
            'algorithm': 'lincs',
            'type': 'h-bonds',
            'lincs_iter': 1,
            'lincs_order': 4,
        },
        'energy_minimization': {
            'integrator': 'steep',
            'emtol': 200.0,
            'emstep': 0.01,
            'nsteps': 10000,
        },
        'output': {
            'trajectory_interval_ps': 500,
            'energy_interval_ps': 10,
            'log_interval_ps': 10,
            'compressed_trajectory': True,
        },
        'electrostatics': {
            'coulombtype': 'PME',
            'rcoulomb': 1.0,
            'pme_order': 4,
            'fourierspacing': 0.16,
        },
        'vdw': {
            'rvdw': 1.0,
            'dispcorr': 'EnerPres',
        },
        'temperature_coupling': {
            'tcoupl': 'V-rescale',
            'tc_grps': ['Protein', 'Non-Protein'],
            'tau_t': [0.1, 0.1],
        },
        'pressure_coupling': {
            'pcoupl': 'C-rescale',
            'pcoupltype': 'isotropic',
            'tau_p': 1.0,
            'compressibility': 4.5e-05,
        },
    }
