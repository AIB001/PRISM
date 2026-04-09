#!/usr/bin/env python3

from pathlib import Path

import pytest

from prism.fep.gromacs.mdp_templates import _detect_hybrid_runtime_guard, _generate_distribution, write_fep_mdps


def _write_hybrid_itp(path: Path, zeroized_h_bond: bool) -> None:
    bond_line = (
        "    1    2    1 0.000000 0.000000 0.109000 284512.000000\n"
        if zeroized_h_bond
        else "    1    2    1 0.109000 284512.000000 0.109000 284512.000000\n"
    )
    path.write_text(
        "[ moleculetype ]\nHYB 3\n\n"
        "[ atoms ]\n"
        "    1   C261   1 LIG   C1   1  0.0 12.011\n"
        "    2   HGA2   1 LIG   H1   1  0.0  1.008\n\n"
        "[ bonds ]\n" + bond_line,
    )


def test_detect_hybrid_runtime_guard_and_write_mdps(tmp_path):
    mdps_dir = tmp_path / "bound" / "repeat1" / "mdps"
    hybrid_dir = tmp_path / "bound" / "common" / "hybrid"
    mdps_dir.mkdir(parents=True)
    hybrid_dir.mkdir(parents=True)
    _write_hybrid_itp(hybrid_dir / "hybrid.itp", zeroized_h_bond=True)

    needs_guard, zeroized = _detect_hybrid_runtime_guard(mdps_dir)
    assert needs_guard is True
    assert zeroized == 1

    write_fep_mdps(str(mdps_dir), config={"forcefield": {"name": "charmm36-jul2022"}}, leg_name="bound")
    nvt = (mdps_dir / "nvt.mdp").read_text()
    npt = (mdps_dir / "npt.mdp").read_text()
    assert "dt                  = 0.0005" in nvt
    assert "constraints             = none" in nvt
    assert "zeroized state" in nvt
    assert "dt                  = 0.0005" in npt
    assert "constraints             = none" in npt


def test_write_mdps_preserves_hbonds_when_no_zeroized_h_bonds(tmp_path):
    mdps_dir = tmp_path / "bound" / "repeat1" / "mdps"
    hybrid_dir = tmp_path / "bound" / "common" / "hybrid"
    mdps_dir.mkdir(parents=True)
    hybrid_dir.mkdir(parents=True)
    _write_hybrid_itp(hybrid_dir / "hybrid.itp", zeroized_h_bond=False)

    needs_guard, zeroized = _detect_hybrid_runtime_guard(mdps_dir)
    assert needs_guard is False
    assert zeroized == 0

    write_fep_mdps(str(mdps_dir), config={"forcefield": {"name": "amber14sb_OL15"}}, leg_name="bound")
    nvt = (mdps_dir / "nvt.mdp").read_text()
    assert "dt                  = 0.002" in nvt
    assert "constraints             = h-bonds" in nvt


def test_quadratic_distribution_exponent_controls_endpoint_density():
    default = _generate_distribution(11, "quadratic", quadratic_exponent=2.0)
    stronger = _generate_distribution(11, "quadratic", quadratic_exponent=4.0)

    # Higher exponent should pack points more tightly near the endpoints and
    # therefore create larger gaps around the center.
    assert stronger[1] < default[1]
    assert stronger[-2] > default[-2]
    assert (stronger[6] - stronger[5]) > (default[6] - default[5])
    assert default[5] == pytest.approx(0.5)
    assert stronger[5] == pytest.approx(0.5)


def test_write_mdps_honors_lambda_quadratic_exponent(tmp_path):
    mdps_dir = tmp_path / "bound" / "repeat1" / "mdps"
    hybrid_dir = tmp_path / "bound" / "common" / "hybrid"
    mdps_dir.mkdir(parents=True)
    hybrid_dir.mkdir(parents=True)
    _write_hybrid_itp(hybrid_dir / "hybrid.itp", zeroized_h_bond=False)

    write_fep_mdps(
        str(mdps_dir),
        lambda_distribution="quadratic",
        lambda_windows=11,
        config={"lambda": {"quadratic_exponent": 4.0}, "forcefield": {"name": "amber14sb_OL15"}},
        leg_name="bound",
    )
    schedule = (mdps_dir / "lambda_schedule.json").read_text()
    assert '"quadratic_exponent": 4.0' in schedule


def test_write_fep_mdps_propagates_output_nstdhdl(tmp_path: Path):
    """Custom output.nstdhdl should propagate into all FEP-specific MDPs."""
    mdp_dir = tmp_path / "bound" / "repeat1" / "mdps"
    write_fep_mdps(
        output_dir=str(mdp_dir),
        config={
            "output": {
                "trajectory_interval_ps": 500,
                "energy_interval_ps": 10,
                "log_interval_ps": 10,
                "nstdhdl": 25,
            }
        },
        leg_name="bound",
    )

    prod_text = (mdp_dir / "prod_00.mdp").read_text()
    em_short_text = (mdp_dir / "em_short_00.mdp").read_text()
    npt_short_text = (mdp_dir / "npt_short_00.mdp").read_text()

    assert "nstdhdl                  = 25" in prod_text
    assert "nstdhdl                  = 25" in em_short_text
    assert "nstdhdl                  = 25" in npt_short_text
    assert "separate-dhdl-file       = yes" in prod_text


def test_write_fep_mdps_rejects_negative_nstdhdl(tmp_path: Path):
    """Negative dhdl output cadence should be rejected early."""
    mdp_dir = tmp_path / "unbound" / "repeat1" / "mdps"

    with pytest.raises(ValueError, match="output.nstdhdl"):
        write_fep_mdps(
            output_dir=str(mdp_dir),
            config={"output": {"nstdhdl": -1}},
            leg_name="unbound",
        )
