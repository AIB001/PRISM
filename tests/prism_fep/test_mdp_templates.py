from pathlib import Path

from prism.fep.gromacs.mdp_templates import write_fep_mdps


def test_write_fep_mdps_uses_conservative_equilibration_dt_for_charmm(tmp_path: Path):
    config = {
        "forcefield": {"name": "charmm36m-mut"},
        "simulation": {
            "dt": 0.002,
            "equilibration_nvt_time_ps": 500,
            "equilibration_npt_time_ps": 500,
            "per_window_npt_time_ps": 100,
            "production_time_ns": 2.0,
            "temperature": 310,
            "pressure": 1.0,
        },
        "fep": {"lambda_windows": 11},
    }

    write_fep_mdps(str(tmp_path), config=config, leg_name="bound")

    nvt = (tmp_path / "nvt.mdp").read_text()
    npt = (tmp_path / "npt.mdp").read_text()

    assert "dt                  = 0.001" in nvt
    assert "nsteps              = 500000" in nvt
    assert "dt                  = 0.001" in npt
    assert "nsteps              = 500000" in npt


def test_write_fep_mdps_preserves_equilibration_dt_for_non_charmm(tmp_path: Path):
    config = {
        "forcefield": {"name": "amber99sb"},
        "simulation": {
            "dt": 0.002,
            "equilibration_nvt_time_ps": 500,
            "equilibration_npt_time_ps": 500,
            "per_window_npt_time_ps": 100,
            "production_time_ns": 2.0,
            "temperature": 310,
            "pressure": 1.0,
        },
        "fep": {"lambda_windows": 11},
    }

    write_fep_mdps(str(tmp_path), config=config, leg_name="bound")

    nvt = (tmp_path / "nvt.mdp").read_text()
    npt = (tmp_path / "npt.mdp").read_text()

    assert "dt                  = 0.002" in nvt
    assert "nsteps              = 250000" in nvt
    assert "dt                  = 0.002" in npt
    assert "nsteps              = 250000" in npt
