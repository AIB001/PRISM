"""Regression tests for the membrane builder and its public integrations."""

from __future__ import annotations

import importlib
import json
import struct
import sys
import urllib.request
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from prism.builder import cli as builder_cli
from prism.membrane.builder import MembraneBuilder
from prism.membrane.config import (
    DEFAULT_PRODUCTION_NS,
    MEMBRANE_TIMESTEP_PS,
    MembraneConfig,
    parse_membrane_cli,
    parse_membrane_yaml,
)
from prism.membrane.membrane_mdp import write_membrane_mdps
from prism.membrane.orient import fetch_opm, insert_ter_at_chain_breaks, orient_protein
from prism.membrane import packmol_memgen


# Load only the MCP build module. The production prism.mcp package eagerly
# imports the optional FastMCP dependency and every MCP tool, neither of which is
# needed by these focused integration tests.
_MCP_DIR = Path(__file__).resolve().parents[1] / "prism" / "mcp"
_MCP_PACKAGE_NAME = "_prism_mcp_under_test"
_mcp_package = ModuleType(_MCP_PACKAGE_NAME)
_mcp_package.__path__ = [str(_MCP_DIR)]
sys.modules[_MCP_PACKAGE_NAME] = _mcp_package
mcp_build = importlib.import_module(f"{_MCP_PACKAGE_NAME}.build")


def _option_value(command, option):
    return command[command.index(option) + 1]


def _write_minimal_amber_topology(path, atom_count=2):
    path.write_text(
        "%VERSION VERSION_STAMP = V0001.000\n"
        "%FLAG POINTERS\n"
        "%FORMAT(10I8)\n"
        f"{atom_count:8d}\n"
    )


def _write_minimal_amber_coordinates(path, atom_count=2):
    path.write_text(f"test coordinates\n{atom_count:6d}\n")


def _write_minimal_amber_netcdf_restart(path, atom_count=2):
    # Minimal CDF-2 header through the AMBER ``atom`` dimension. The production
    # selector needs no coordinate payload to validate NATOM.
    path.write_bytes(
        b"CDF\x02"
        + struct.pack(">I", 0)
        + struct.pack(">II", 10, 1)
        + struct.pack(">I", 4)
        + b"atom"
        + struct.pack(">I", atom_count)
    )


def test_membrane_config_preserves_runtime_settings_and_rejects_bad_ratios():
    config = parse_membrane_yaml(
        {
            "lipids": ["popc", "chl1"],
            "ratio": [3, 1],
            "orient": "OPM",
            "pdb_id": "4XB4",
            "positive_ion": "K+",
            "negative_ion": "Cl-",
            "production_ns": 25,
        },
        temperature=310,
        production_ns=500,
    )

    assert config.lipids == ["POPC", "CHL1"]
    assert config.orient == "opm"
    assert config.positive_ion == "K+"
    assert config.production_ns == 25
    assert config.resolved_ratio() == [3, 1]

    with pytest.raises(ValueError, match="must match lipid count"):
        MembraneConfig(lipids=["POPC", "CHL1"], ratio=[1]).resolved_ratio()

    # Presence of an empty YAML section means "enabled with defaults", not
    # the same thing as an absent/null section.
    assert parse_membrane_yaml({}).enabled is True


def test_membrane_public_api_uses_the_canonical_production_default(tmp_path):
    assert DEFAULT_PRODUCTION_NS == 500.0
    assert MEMBRANE_TIMESTEP_PS == 0.002
    assert MembraneConfig().production_ns == DEFAULT_PRODUCTION_NS
    assert parse_membrane_yaml({}).production_ns == DEFAULT_PRODUCTION_NS
    config = parse_membrane_cli(
        membrane=True,
        lipid=None,
        lipid_ratio=None,
        orient=None,
        pdb_id=None,
        lipid_ff=None,
    )
    assert config.production_ns == DEFAULT_PRODUCTION_NS
    assert config.lipids == ["POPC"]
    assert config.packmol_forcefield_pair() == ("ff14SB", "tip3p")

    mdps = write_membrane_mdps(str(tmp_path))
    assert "nsteps                  = 250000000" in Path(mdps["md"]).read_text()
    assert "pcoupltype              = semiisotropic" in Path(mdps["npt"]).read_text()


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ({"enabled": "false"}, "enabled must be a boolean"),
        ({"minimize": "false"}, "minimize must be a boolean"),
        ({"ratio": [1.5]}, "floats and strings are not allowed"),
        ({"ratio": ["2"]}, "floats and strings are not allowed"),
        ({"orient": "sideways"}, "Unknown membrane orientation method"),
        ({"lipids": []}, "At least one membrane lipid"),
    ],
)
def test_membrane_yaml_rejects_lossy_or_ambiguous_values(payload, message):
    with pytest.raises(ValueError, match=message):
        parse_membrane_yaml(payload)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("salt_concentration", float("nan")),
        ("water_thickness", float("inf")),
        ("xy_distance", float("-inf")),
        ("temperature", float("nan")),
        ("production_ns", float("inf")),
    ],
)
def test_membrane_yaml_rejects_nonfinite_numbers(field, value):
    with pytest.raises(ValueError, match="finite"):
        parse_membrane_yaml({field: value})


def test_invalid_direct_config_never_launches_packmol_or_creates_output(monkeypatch, tmp_path):
    builder = MembraneBuilder(
        MembraneConfig(enabled=True, salt_concentration=float("nan")),
        verbose=False,
    )
    monkeypatch.setattr(
        packmol_memgen,
        "run",
        lambda *_args, **_kwargs: pytest.fail("invalid config reached PACKMOL-Memgen"),
    )
    output = tmp_path / "invalid-output"

    with pytest.raises(ValueError, match="NaN/Inf"):
        builder.build(str(tmp_path / "protein.pdb"), str(output))

    assert not output.exists()


@pytest.mark.parametrize(
    (
        "membrane_yaml",
        "membrane_section",
        "expected_forcefield",
        "expected_water",
        "expected_runtime",
    ),
    [
        (
            "membrane:\n  enabled: true\n",
            {"enabled": True},
            "amber14sb",
            "tip3p",
            {
                "temperature": 310,
                "production_ns": 500,
                "salt_concentration": 0.15,
                "positive_ion": "Na+",
                "negative_ion": "Cl-",
            },
        ),
        (
            "membrane:\n  enabled: true\n  protein_forcefield: amber19sb\n  water_model: opc\n",
            {
                "enabled": True,
                "protein_forcefield": "amber19sb",
                "water_model": "opc",
            },
            "amber19sb",
            "opc",
            {
                "temperature": 310,
                "production_ns": 500,
                "salt_concentration": 0.15,
                "positive_ion": "Na+",
                "negative_ion": "Cl-",
            },
        ),
        (
            "simulation:\n  temperature: 315\n  production_ns: 25\n"
            "ions:\n  concentration: 0.2\n  positive_ion: K+\n  negative_ion: Br-\n"
            "membrane:\n  enabled: true\n",
            {"enabled": True},
            "amber14sb",
            "tip3p",
            {
                "temperature": 315,
                "production_ns": 25,
                "salt_concentration": 0.2,
                "positive_ion": "K+",
                "negative_ion": "Br-",
            },
        ),
    ],
)
def test_prismbuilder_python_api_detects_yaml_membrane_before_forcefield_defaults(
    monkeypatch,
    tmp_path,
    membrane_yaml,
    membrane_section,
    expected_forcefield,
    expected_water,
    expected_runtime,
):
    from prism.builder import core as builder_core

    config_path = tmp_path / "membrane.yaml"
    config_path.write_text(membrane_yaml, encoding="utf-8")
    captured = {}

    class FakeEnvironment:
        pass

    class FakeConfigurationManager:
        def __init__(
            self,
            _config_path,
            _gromacs_env,
            forcefield_name="amber99sb",
            water_model_name="tip3p",
        ):
            captured["forcefield_name"] = forcefield_name
            captured["water_model_name"] = water_model_name
            self.config = {
                "general": {"overwrite": False},
                "forcefield": {
                    "index": 1,
                    "custom_forcefields": {1: {"name": forcefield_name}},
                },
                "water_model": {
                    "index": 1,
                    "custom_water_models": {1: {"name": water_model_name}},
                },
                "simulation": {"temperature": 310, "production_time_ns": 500},
                "membrane": dict(membrane_section),
                "box": {"distance": 1.5},
                "ions": {},
                "energy_minimization": {},
            }

    monkeypatch.setattr(builder_core, "GromacsEnvironment", FakeEnvironment)
    monkeypatch.setattr(builder_core, "ConfigurationManager", FakeConfigurationManager)
    monkeypatch.setattr(builder_core, "MDPGenerator", lambda *_args, **_kwargs: SimpleNamespace())
    monkeypatch.setattr(builder_core, "SystemBuilder", lambda *_args, **_kwargs: SimpleNamespace())
    monkeypatch.setattr(builder_core.PRISMBuilder, "_print_initialization_info", lambda _self: None)

    builder = builder_core.PRISMBuilder(
        protein_path=str(tmp_path / "protein.pdb"),
        ligand_paths=[],
        output_dir=str(tmp_path / "output"),
        config_path=str(config_path),
        forcefield=None,
        water_model=None,
    )

    assert captured["forcefield_name"] == expected_forcefield
    assert captured["water_model_name"] == expected_water
    assert builder.membrane_mode is True
    assert builder.membrane_config.protein_forcefield == expected_forcefield
    assert builder.membrane_config.water_model == expected_water
    for field_name, expected_value in expected_runtime.items():
        assert getattr(builder.membrane_config, field_name) == expected_value


def test_prismbuilder_membrane_does_not_require_packmol_forcefield_in_gromacs(
    monkeypatch, tmp_path
):
    """PACKMOL-Memgen FF pairs must not be rejected by GROMACS discovery.

    The membrane backend produces a self-contained AMBER topology and converts
    it with ParmEd.  Requiring the same protein/water pair to exist as a local
    GROMACS ``.ff`` directory prevents valid ff19SB/OPC builds before the
    membrane backend is reached.
    """
    from prism.builder import core as builder_core

    class GromacsWithoutPackmolForcefields:
        gmx_command = "gmx"

        def get_force_field_config(self):
            pytest.fail("membrane initialization queried GROMACS force fields")

    monkeypatch.setattr(
        builder_core, "GromacsEnvironment", GromacsWithoutPackmolForcefields
    )
    monkeypatch.setattr(
        builder_core, "MDPGenerator", lambda *_args, **_kwargs: SimpleNamespace()
    )
    monkeypatch.setattr(
        builder_core, "SystemBuilder", lambda *_args, **_kwargs: SimpleNamespace()
    )
    monkeypatch.setattr(
        builder_core.PRISMBuilder, "_print_initialization_info", lambda _self: None
    )

    membrane_config = MembraneConfig(
        enabled=True,
        protein_forcefield="amber19sb",
        water_model="opc",
    )
    builder = builder_core.PRISMBuilder(
        protein_path=str(tmp_path / "protein.pdb"),
        ligand_paths=[],
        output_dir=str(tmp_path / "output"),
        membrane_config=membrane_config,
    )

    assert builder.forcefield == {
        "name": "amber19sb",
        "dir": None,
        "source": "PACKMOL-Memgen",
    }
    assert builder.water_model == {"name": "opc", "source": "PACKMOL-Memgen"}
    assert builder.membrane_config.packmol_forcefield_pair() == ("ff19SB", "opc")


@pytest.mark.parametrize(
    ("invalid_config", "message"),
    [
        (
            MembraneConfig(enabled=True, salt_concentration=float("nan")),
            "NaN/Inf",
        ),
        (MembraneConfig(enabled=""), "enabled must be a boolean"),
    ],
)
def test_prismbuilder_rejects_invalid_direct_membrane_config_before_output(
    monkeypatch,
    tmp_path,
    invalid_config,
    message,
):
    from prism.builder import core as builder_core

    class FakeEnvironment:
        pass

    class FakeConfigurationManager:
        def __init__(
            self,
            _config_path,
            _gromacs_env,
            forcefield_name="amber99sb",
            water_model_name="tip3p",
        ):
            self.config = {
                "general": {"overwrite": False},
                "forcefield": {
                    "index": 1,
                    "custom_forcefields": {1: {"name": forcefield_name}},
                },
                "water_model": {
                    "index": 1,
                    "custom_water_models": {1: {"name": water_model_name}},
                },
                "simulation": {"temperature": 310, "production_time_ns": 500},
                "box": {"distance": 1.5},
                "ions": {},
                "energy_minimization": {},
            }

    monkeypatch.setattr(builder_core, "GromacsEnvironment", FakeEnvironment)
    monkeypatch.setattr(builder_core, "ConfigurationManager", FakeConfigurationManager)
    monkeypatch.setattr(builder_core, "MDPGenerator", lambda *_args, **_kwargs: SimpleNamespace())
    monkeypatch.setattr(builder_core, "SystemBuilder", lambda *_args, **_kwargs: SimpleNamespace())
    monkeypatch.setattr(builder_core.PRISMBuilder, "_print_initialization_info", lambda _self: None)

    output = tmp_path / "invalid-python-api-output"
    with pytest.raises(ValueError, match=message):
        builder_core.PRISMBuilder(
            protein_path=str(tmp_path / "protein.pdb"),
            ligand_paths=[],
            output_dir=str(output),
            membrane_config=invalid_config,
        )

    assert not output.exists()


@pytest.mark.parametrize(
    ("control", "message"),
    [
        ("pressure: 2.0", r"simulation\.pressure"),
        ("nvt_ps: 250", r"simulation\.nvt_ps"),
    ],
)
def test_prismbuilder_python_api_rejects_unconsumed_membrane_yaml_controls(
    tmp_path,
    control,
    message,
):
    from prism.builder import core as builder_core

    config_path = tmp_path / "unsupported-membrane.yaml"
    config_path.write_text(
        f"simulation:\n  {control}\nmembrane:\n  enabled: true\n",
        encoding="utf-8",
    )
    output = tmp_path / "unsupported-python-api-output"

    with pytest.raises(ValueError, match=message):
        builder_core.PRISMBuilder(
            protein_path=str(tmp_path / "protein.pdb"),
            ligand_paths=[],
            output_dir=str(output),
            config_path=str(config_path),
        )

    assert not output.exists()


@pytest.mark.parametrize(
    ("prefix", "message"),
    [
        ("forcefield: amber19sb\n", r"forcefield.*mapping"),
        ("forcefield:\n  name: 123\n", r"protein force field.*non-empty string"),
        ("simulation: []\n", r"simulation.*mapping"),
        ("ions: 0\n", r"ions.*mapping"),
    ],
)
def test_prismbuilder_python_api_rejects_malformed_membrane_yaml_dependencies(
    tmp_path,
    prefix,
    message,
):
    from prism.builder import core as builder_core

    config_path = tmp_path / "malformed-membrane.yaml"
    config_path.write_text(
        prefix + "membrane:\n  enabled: true\n",
        encoding="utf-8",
    )
    output = tmp_path / "malformed-python-api-output"

    with pytest.raises(ValueError, match=message):
        builder_core.PRISMBuilder(
            protein_path=str(tmp_path / "protein.pdb"),
            ligand_paths=[],
            output_dir=str(output),
            config_path=str(config_path),
        )

    assert not output.exists()


def test_packmol_command_selects_lipid21_and_normalizes_cli_ions():
    config = parse_membrane_cli(
        membrane=True,
        lipid=["POPC", "CHL1"],
        lipid_ratio=[3, 1],
        orient="OPM",
        pdb_id="4xb4",
        lipid_ff="lipid21",
        positive_ion="NA",
        negative_ion="CL",
    )

    command = packmol_memgen.build_command("protein.pdb", config)

    assert _option_value(command, "--lipids") == "POPC:CHL1"
    assert _option_value(command, "--ratio") == "3:1"
    assert _option_value(command, "--fflip") == "lipid21"
    assert _option_value(command, "--salt_c") == "Na+"
    assert _option_value(command, "--salt_a") == "Cl-"
    assert "--preoriented" in command
    assert "--fprot" not in command
    assert "--ffprot" not in command


def test_packmol_command_propagates_nondefault_protein_and_water_forcefields():
    config = MembraneConfig(
        enabled=True,
        lipids=["POPC", "CHL1"],
        ratio=[3, 1],
        protein_forcefield="amber19sb",
        water_model="opc",
    )

    command = packmol_memgen.build_command(
        "protein.pdb",
        config,
        protein_ff_option="--ffprot",
    )

    assert _option_value(command, "--ffprot") == "ff19SB"
    assert _option_value(command, "--ffwat") == "opc"
    assert _option_value(command, "--lipids") == "POPC:CHL1"
    assert _option_value(command, "--ratio") == "3:1"

    config.water_model = "tip3p"
    with pytest.raises(ValueError, match="pairing"):
        packmol_memgen.build_command("protein.pdb", config)


def test_packmol_detects_forcefield_flag_spelling(monkeypatch):
    monkeypatch.setattr(
        packmol_memgen.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(
            returncode=0,
            stdout="options: --ffprot PROTEIN --ffwat WATER",
            stderr="",
        ),
    )

    assert packmol_memgen._detect_forcefield_options() == ("--ffprot", "--ffwat")


def test_packmol_command_refuses_unimplemented_parametrization_backend():
    config = MembraneConfig(enabled=True, lipid_ff="charmm36")

    with pytest.raises(ValueError, match="does not support"):
        packmol_memgen.build_command("protein.pdb", config)


def test_packmol_result_never_combines_mismatched_amber_files(monkeypatch, tmp_path):
    (tmp_path / "old_a.top").write_text("topology")
    (tmp_path / "old_b.crd").write_text("coordinates")
    (tmp_path / "topol.top").write_text("PRISM GROMACS topology")

    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)
    monkeypatch.setattr(
        packmol_memgen.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=0, stdout="", stderr=""),
    )

    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=False),
        str(tmp_path),
    )

    assert result.success is False
    assert result.prmtop is None
    assert result.inpcrd is None


def test_packmol_result_does_not_reuse_stale_matching_pair(monkeypatch, tmp_path):
    (tmp_path / "old_system.top").write_text("topology")
    (tmp_path / "old_system.crd").write_text("coordinates")

    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)
    monkeypatch.setattr(
        packmol_memgen.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=0, stdout="", stderr=""),
    )

    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=False),
        str(tmp_path),
    )

    assert result.success is False
    assert result.prmtop is None
    assert result.inpcrd is None


def test_packmol_result_uses_matching_amber_pair(monkeypatch, tmp_path):
    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)

    def fake_run(*_args, cwd, **_kwargs):
        _write_minimal_amber_topology(Path(cwd, "bilayer_system.top"))
        _write_minimal_amber_coordinates(Path(cwd, "bilayer_system.crd"))
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(packmol_memgen.subprocess, "run", fake_run)
    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=False),
        str(tmp_path),
    )

    assert result.success is True
    assert Path(result.prmtop).name == "bilayer_system.top"
    assert Path(result.inpcrd).name == "bilayer_system.crd"


def test_packmol_run_absolutizes_relative_input_before_changing_cwd(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    protein = Path("protein with space.pdb")
    protein.write_text("ATOM\n")
    captured = {}
    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)

    def fake_run(command, cwd, **_kwargs):
        captured["command"] = command
        captured["cwd"] = cwd
        _write_minimal_amber_topology(Path(cwd, "bilayer_system.top"))
        _write_minimal_amber_coordinates(Path(cwd, "bilayer_system.crd"))
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(packmol_memgen.subprocess, "run", fake_run)
    result = packmol_memgen.run(
        str(protein),
        MembraneConfig(enabled=True, minimize=False),
        "work dir",
    )

    command_pdb = _option_value(captured["command"], "--pdb")
    assert Path(command_pdb).is_absolute()
    assert Path(command_pdb) == protein.resolve()
    assert Path(captured["cwd"]).is_absolute()
    assert result.success is True


def test_packmol_result_does_not_report_stale_packed_pdb(monkeypatch, tmp_path):
    stale_pdb = tmp_path / "bilayer_old.pdb"
    stale_pdb.write_text("stale")
    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)

    def fake_run(*_args, cwd, **_kwargs):
        _write_minimal_amber_topology(Path(cwd, "bilayer_system.top"))
        _write_minimal_amber_coordinates(Path(cwd, "bilayer_system.crd"))
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(packmol_memgen.subprocess, "run", fake_run)
    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=False),
        str(tmp_path),
    )

    assert result.success is True
    assert result.packed_pdb is None


def test_packmol_nonzero_exit_never_exposes_new_amber_pair(monkeypatch, tmp_path):
    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)

    def fake_run(*_args, cwd, **_kwargs):
        _write_minimal_amber_topology(Path(cwd, "bilayer_system.top"))
        _write_minimal_amber_coordinates(Path(cwd, "bilayer_system.crd"))
        return SimpleNamespace(returncode=2, stdout="", stderr="failed")

    monkeypatch.setattr(packmol_memgen.subprocess, "run", fake_run)
    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=False),
        str(tmp_path),
    )

    assert result.success is False
    assert result.prmtop is None
    assert result.inpcrd is None


def test_packmol_result_prefers_final_unrestrained_minimized_restart(monkeypatch, tmp_path):
    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)

    def fake_run(*_args, cwd, **_kwargs):
        _write_minimal_amber_topology(Path(cwd, "bilayer_protein_oriented_lipid.top"), 5)
        _write_minimal_amber_coordinates(Path(cwd, "bilayer_protein_oriented_lipid.crd"), 5)
        _write_minimal_amber_netcdf_restart(Path(cwd, "min.restrt"), 5)
        _write_minimal_amber_netcdf_restart(Path(cwd, "bilayer_protein_oriented_min.restrt"), 5)
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(packmol_memgen.subprocess, "run", fake_run)
    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=True),
        str(tmp_path),
    )

    assert result.success is True
    assert Path(result.prmtop).name == "bilayer_protein_oriented_lipid.top"
    assert Path(result.inpcrd).name == "bilayer_protein_oriented_min.restrt"


def test_packmol_result_rejects_minimized_restart_with_wrong_atom_count(monkeypatch, tmp_path):
    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)

    def fake_run(*_args, cwd, **_kwargs):
        _write_minimal_amber_topology(Path(cwd, "bilayer_protein_oriented_lipid.top"), 5)
        _write_minimal_amber_coordinates(Path(cwd, "bilayer_protein_oriented_lipid.crd"), 5)
        _write_minimal_amber_netcdf_restart(Path(cwd, "bilayer_protein_oriented_min.restrt"), 4)
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(packmol_memgen.subprocess, "run", fake_run)
    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=True),
        str(tmp_path),
    )

    assert result.success is False
    assert result.prmtop is None
    assert result.inpcrd is None


def test_packmol_minimize_does_not_fall_back_when_final_restart_is_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)

    def fake_run(*_args, cwd, **_kwargs):
        _write_minimal_amber_topology(Path(cwd, "bilayer_protein_oriented_lipid.top"), 5)
        _write_minimal_amber_coordinates(Path(cwd, "bilayer_protein_oriented_lipid.crd"), 5)
        _write_minimal_amber_netcdf_restart(Path(cwd, "min.restrt"), 5)
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(packmol_memgen.subprocess, "run", fake_run)
    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=True),
        str(tmp_path),
    )

    assert result.success is False
    assert result.prmtop is None
    assert result.inpcrd is None


def test_packmol_minimize_does_not_reuse_stale_final_restart(monkeypatch, tmp_path):
    _write_minimal_amber_netcdf_restart(
        tmp_path / "bilayer_protein_oriented_min.restrt",
        5,
    )
    monkeypatch.setattr(packmol_memgen, "is_available", lambda: True)

    def fake_run(*_args, cwd, **_kwargs):
        _write_minimal_amber_topology(Path(cwd, "bilayer_protein_oriented_lipid.top"), 5)
        _write_minimal_amber_coordinates(Path(cwd, "bilayer_protein_oriented_lipid.crd"), 5)
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(packmol_memgen.subprocess, "run", fake_run)
    result = packmol_memgen.run(
        "protein.pdb",
        MembraneConfig(enabled=True, minimize=True),
        str(tmp_path),
    )

    assert result.success is False
    assert result.prmtop is None
    assert result.inpcrd is None


def test_membrane_mdp_honors_production_duration_and_normalizes_family(tmp_path):
    files = write_membrane_mdps(
        str(tmp_path),
        temperature=310,
        ff_family="CHARMM",
        production_ns=2.0,
    )

    production = Path(files["md"]).read_text()
    assert "nsteps                  = 1000000" in production
    assert "vdw-modifier            = Force-switch" in production
    assert "pcoupltype              = semiisotropic" in production

    with pytest.raises(ValueError, match="positive"):
        write_membrane_mdps(str(tmp_path / "invalid"), production_ns=0)


def test_substep_membrane_duration_fails_before_creating_mdp_directory(tmp_path):
    output = tmp_path / "substep"

    with pytest.raises(ValueError, match="at least one.*integration step"):
        write_membrane_mdps(str(output), production_ns=1e-6)

    assert not output.exists()

    with pytest.raises(ValueError, match="at least one.*integration step"):
        parse_membrane_yaml({"production_ns": 1e-6})


def test_fetch_opm_accepts_a_valid_pdb_after_a_long_header(monkeypatch, tmp_path):
    payload = b"REMARK metadata\n" * 500 + b"ATOM      1  CA  ALA A   1\n"

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return payload

    monkeypatch.setattr(urllib.request, "urlopen", lambda *_args, **_kwargs: Response())
    output = tmp_path / "opm.pdb"

    assert fetch_opm("4xb4", str(output)) == str(output)
    assert output.read_bytes() == payload

    payload = b"<html><body>OPM error: ATOM record unavailable</body></html>"
    rejected = tmp_path / "opm_error.pdb"
    assert fetch_opm("missing", str(rejected)) is None
    assert not rejected.exists()


def _backbone_residue(serial, resname, chain, resid, n, ca, c):
    lines = []
    for atom, coordinate, element in (("N", n, "N"), ("CA", ca, "C"), ("C", c, "C")):
        x, y, z = coordinate
        lines.append(
            f"ATOM  {serial:5d} {atom:^4s} {resname:>3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        )
        serial += 1
    return "".join(lines), serial


def test_chain_break_preparation_inserts_ter_for_numbering_and_geometry_gaps(tmp_path):
    source = tmp_path / "broken.pdb"
    output = tmp_path / "prepared.pdb"
    blocks = []
    serial = 1
    block, serial = _backbone_residue(
        serial, "ALA", "A", 1, (0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)
    )
    blocks.append(block)
    # Consecutive and peptide-bond close: no TER here.
    block, serial = _backbone_residue(
        serial, "GLY", "A", 2, (2.3, 0.0, 0.0), (2.8, 0.0, 0.0), (3.3, 0.0, 0.0)
    )
    blocks.append(block)
    # Residues 3-10 are absent, matching the real 2K0L A165 -> A177 defect.
    block, serial = _backbone_residue(
        serial, "THR", "A", 11, (20.0, 0.0, 0.0), (20.5, 0.0, 0.0), (21.0, 0.0, 0.0)
    )
    blocks.append(block)
    # Consecutive numbering but impossible C-N distance also requires a break.
    block, serial = _backbone_residue(
        serial, "VAL", "A", 12, (30.0, 0.0, 0.0), (30.5, 0.0, 0.0), (31.0, 0.0, 0.0)
    )
    blocks.append(block)
    source.write_text("".join(blocks) + "END\n", encoding="ascii")

    assert insert_ter_at_chain_breaks(str(source), str(output)) == 2
    lines = output.read_text(encoding="ascii").splitlines()
    assert sum(line.startswith("TER") for line in lines) == 2
    assert lines.index("TER") == 6  # after two complete, genuinely connected residues


def test_preoriented_path_inserts_chain_break_ter_in_place_and_reports_it(tmp_path):
    source = tmp_path / "broken.pdb"
    first, serial = _backbone_residue(
        1, "TRP", "A", 165, (0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)
    )
    second, _ = _backbone_residue(
        serial, "THR", "A", 177, (15.0, 0.0, 0.0), (15.5, 0.0, 0.0), (16.0, 0.0, 0.0)
    )
    source.write_text(first + second + "END\n", encoding="ascii")

    result = orient_protein(str(source), str(source), "preoriented")

    assert result["oriented_pdb"] == str(source)
    assert "Inserted 1 TER" in result["note"]
    assert source.read_text(encoding="ascii").splitlines().count("TER") == 1


def test_chain_break_preparation_keeps_close_same_chain_residues_despite_numbering_gap(tmp_path):
    source = tmp_path / "renumbered.pdb"
    output = tmp_path / "prepared.pdb"
    first, serial = _backbone_residue(
        1, "ALA", "A", 1, (0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)
    )
    second, _ = _backbone_residue(
        serial, "GLY", "A", 11, (2.3, 0.0, 0.0), (2.8, 0.0, 0.0), (3.3, 0.0, 0.0)
    )
    source.write_text(first + second + "END\n", encoding="ascii")

    assert insert_ter_at_chain_breaks(str(source), str(output)) == 0
    assert "TER" not in output.read_text(encoding="ascii").splitlines()


def test_chain_break_preparation_uses_numbering_gap_when_geometry_is_missing(tmp_path):
    source = tmp_path / "numbering_fallback.pdb"
    output = tmp_path / "prepared.pdb"
    first, serial = _backbone_residue(
        1, "ALA", "A", 1, (0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)
    )
    second, _ = _backbone_residue(
        serial, "GLY", "A", 11, (2.3, 0.0, 0.0), (2.8, 0.0, 0.0), (3.3, 0.0, 0.0)
    )
    first_without_c = "".join(
        line
        for line in first.splitlines(keepends=True)
        if line[12:16].strip() != "C"
    )
    source.write_text(first_without_c + second + "END\n", encoding="ascii")

    assert insert_ter_at_chain_breaks(str(source), str(output)) == 1
    assert output.read_text(encoding="ascii").splitlines().count("TER") == 1


def test_chain_break_preparation_splits_changed_chain_with_incomplete_backbone(tmp_path):
    source = tmp_path / "changed_chain.pdb"
    output = tmp_path / "prepared.pdb"
    first, serial = _backbone_residue(
        1, "ALA", "A", 1, (0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)
    )
    second, _ = _backbone_residue(
        serial, "GLY", "B", 1, (2.3, 0.0, 0.0), (2.8, 0.0, 0.0), (3.3, 0.0, 0.0)
    )

    def without_ca(block):
        return "".join(
            line
            for line in block.splitlines(keepends=True)
            if line[12:16].strip() != "CA"
        )

    source.write_text(without_ca(first) + without_ca(second) + "END\n", encoding="ascii")

    assert insert_ter_at_chain_breaks(str(source), str(output)) == 1
    assert output.read_text(encoding="ascii").splitlines().count("TER") == 1


def test_chain_break_preparation_preserves_existing_ter_idempotently(tmp_path):
    source = tmp_path / "terminated.pdb"
    first, serial = _backbone_residue(
        1, "ALA", "A", 1, (0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)
    )
    second, _ = _backbone_residue(
        serial, "GLY", "B", 1, (2.3, 0.0, 0.0), (2.8, 0.0, 0.0), (3.3, 0.0, 0.0)
    )
    source.write_text(first + "TER\n" + second + "END\n", encoding="ascii")

    assert insert_ter_at_chain_breaks(str(source), str(source)) == 0
    assert insert_ter_at_chain_breaks(str(source), str(source)) == 0
    assert source.read_text(encoding="ascii").splitlines().count("TER") == 1


def test_orientation_rejects_missing_opm_id_and_unknown_method(tmp_path):
    input_pdb = tmp_path / "protein.pdb"
    input_pdb.write_text("ATOM\n")

    with pytest.raises(ValueError, match="requires a PDB id"):
        orient_protein(str(input_pdb), str(tmp_path / "out.pdb"), "opm")
    with pytest.raises(ValueError, match="Unsupported"):
        orient_protein(str(input_pdb), str(tmp_path / "out.pdb"), "typo")


def test_orientation_never_treats_failed_opm_or_unimplemented_ppm_as_preoriented(monkeypatch, tmp_path):
    input_pdb = tmp_path / "protein.pdb"
    input_pdb.write_text("ATOM\n")
    output_pdb = tmp_path / "oriented.pdb"
    monkeypatch.setattr("prism.membrane.orient.fetch_opm", lambda *_args, **_kwargs: None)

    with pytest.raises(RuntimeError, match="OPM fetch failed"):
        orient_protein(str(input_pdb), str(output_pdb), "opm", "4xb4")
    with pytest.raises(NotImplementedError, match="PPM"):
        orient_protein(str(input_pdb), str(output_pdb), "ppm")

    assert not output_pdb.exists()


def test_membrane_index_requires_all_thermostat_groups(tmp_path):
    builder = MembraneBuilder(MembraneConfig(enabled=True), verbose=False)
    atoms = [
        SimpleNamespace(residue=SimpleNamespace(name="ALA")),
        SimpleNamespace(residue=SimpleNamespace(name="POPC")),
    ]
    index_path = tmp_path / "index.ndx"

    with pytest.raises(ValueError, match="SOLV"):
        builder._write_membrane_index(SimpleNamespace(atoms=atoms), str(index_path))
    assert not index_path.exists()


def test_membrane_index_never_classifies_unknown_residue_as_lipid(tmp_path):
    builder = MembraneBuilder(MembraneConfig(enabled=True), verbose=False)
    atoms = [
        SimpleNamespace(residue=SimpleNamespace(name="ALA")),
        SimpleNamespace(residue=SimpleNamespace(name="POPC")),
        SimpleNamespace(residue=SimpleNamespace(name="WAT")),
        SimpleNamespace(residue=SimpleNamespace(name="MSE")),
    ]
    index_path = tmp_path / "index.ndx"

    with pytest.raises(ValueError, match="MSE"):
        builder._write_membrane_index(SimpleNamespace(atoms=atoms), str(index_path))
    assert not index_path.exists()


@pytest.mark.parametrize(
    ("lipids", "amber_residue_names"),
    [
        (["POPC", "CHL1"], ["PA", "PC", "OL", "CHL"]),
        (["POPE", "DOPC"], ["PA", "PE", "OL", "PC"]),
    ],
)
def test_membrane_index_classifies_amber_lipid_fragments(
    tmp_path, lipids, amber_residue_names
):
    builder = MembraneBuilder(
        MembraneConfig(enabled=True, lipids=lipids), verbose=False
    )
    residue_names = ["ALA", *amber_residue_names, "WAT"]
    atoms = [
        SimpleNamespace(residue=SimpleNamespace(name=name))
        for name in residue_names
    ]
    index_path = tmp_path / "index.ndx"

    builder._write_membrane_index(SimpleNamespace(atoms=atoms), str(index_path))

    groups = {}
    current = None
    for line in index_path.read_text().splitlines():
        if line.startswith("["):
            current = line.strip("[] ")
            groups[current] = []
        elif line.strip():
            groups[current].extend(int(value) for value in line.split())

    assert groups["SOLU"] == [1]
    assert groups["MEMB"] == list(range(2, 2 + len(amber_residue_names)))
    assert groups["SOLV"] == [len(residue_names)]


def test_parmed_conversion_does_not_accept_unchanged_stale_outputs(monkeypatch, tmp_path):
    top = tmp_path / "topol.top"
    gro = tmp_path / "system.gro"
    top.write_text("stale topology")
    gro.write_text("stale coordinates")

    class NoOpParm:
        atoms = []

        def save(self, *_args, **_kwargs):
            return None

    monkeypatch.setitem(
        sys.modules,
        "parmed",
        SimpleNamespace(load_file=lambda *_args, **_kwargs: NoOpParm()),
    )
    result = SimpleNamespace(warnings=[], ndx=None)
    builder = MembraneBuilder(MembraneConfig(enabled=True), verbose=False)

    assert builder._amber_to_gromacs(
        "new.prmtop",
        "new.inpcrd",
        str(top),
        str(gro),
        result,
    ) is False
    assert any("stale files" in warning for warning in result.warnings)


def test_membrane_builder_scaffold_uses_configured_production_time(monkeypatch, tmp_path):
    input_pdb = tmp_path / "protein.pdb"
    input_pdb.write_text("ATOM\n")
    config = MembraneConfig(enabled=True, production_ns=3.0)
    builder = MembraneBuilder(config, verbose=False)
    monkeypatch.setattr(builder, "capability_report", lambda: ["AmberTools unavailable"])

    result = builder.build(str(input_pdb), str(tmp_path / "output"))

    assert result.success is False
    assert result.oriented_pdb
    assert "nsteps                  = 1500000" in Path(result.mdps["md"]).read_text()


class _ToolRegistry:
    def __init__(self):
        self.tools = {}

    def tool(self):
        def decorator(function):
            self.tools[function.__name__] = function
            return function

        return decorator


def test_membrane_entrypoints_reject_unapplied_protocol_controls(monkeypatch, tmp_path):
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM\n")
    output = tmp_path / "output"

    class UnexpectedBuilder:
        def __init__(self, *_args, **_kwargs):
            raise AssertionError("builder must not run for an unsupported membrane control")

    monkeypatch.setattr(builder_cli, "PRISMBuilder", UnexpectedBuilder)
    monkeypatch.setattr(
        sys,
        "argv",
        ["prism", str(protein), "--membrane", "--pressure", "2.0", "-o", str(output)],
    )
    with pytest.raises(SystemExit):
        builder_cli.main()
    assert not output.exists()

    config = tmp_path / "unsupported.yaml"
    config.write_text(
        "simulation:\n  pressure: 2.0\nmembrane:\n  enabled: true\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["prism", str(protein), "--config", str(config), "-o", str(output)],
    )
    with pytest.raises(SystemExit):
        builder_cli.main()
    assert not output.exists()

    registry = _ToolRegistry()
    mcp_build.register(registry)
    response = json.loads(
        registry.tools["build_system"](
            protein_path=str(protein.resolve()),
            ligand_paths="",
            output_dir=str(output),
            membrane=True,
            pressure=2.0,
        )
    )
    assert response["success"] is False
    assert any("pressure" in error for error in response["errors"])
    assert not output.exists()


def test_mcp_membrane_build_allows_no_ligand_and_reports_partial(monkeypatch, tmp_path):
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM\n")
    captured = {}

    class FakeBuilder:
        def __init__(self, **kwargs):
            captured.update(kwargs)
            self.output_dir = kwargs["output_dir"]
            self.config = {
                "simulation": {},
                "ligand_forcefield": {},
                "box": {},
                "ions": {},
                "energy_minimization": {},
            }

        def run(self):
            output = Path(self.output_dir)
            system_dir = output / "GMX_PROLIG_MEMB"
            system_dir.mkdir(parents=True)
            # Simulate complete files left by an older run. A partial response
            # must not present them as artifacts of the current invocation.
            (system_dir / "system.gro").write_text("stale")
            (system_dir / "topol.top").write_text("stale")
            (system_dir / "index.ndx").write_text("stale")
            mdps = output / "mdps"
            mdps.mkdir()
            for name in ("em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"):
                (mdps / name).write_text("scaffold")
            self.membrane_result = SimpleNamespace(
                success=False,
                gro=None,
                top=None,
                ndx=None,
                mdps={stage: str(mdps / f"{stage}.mdp") for stage in ("em", "nvt", "npt", "md")},
                note="Membrane build prerequisites missing.",
                warnings=["packmol-memgen not found"],
            )
            return self.output_dir

    import prism.builder

    monkeypatch.setattr(prism.builder, "PRISMBuilder", FakeBuilder)
    monkeypatch.setattr(mcp_build, "_ensure_prism_importable", lambda: None)
    registry = _ToolRegistry()
    mcp_build.register(registry)

    response = json.loads(
        registry.tools["build_system"](
            protein_path=str(protein.resolve()),
            ligand_paths="",
            output_dir=str(tmp_path / "output"),
            membrane=True,
            production_ns=42,
            lipid_ff="lipid17",
        )
    )

    assert captured["ligand_paths"] == []
    assert captured["membrane_config"].production_ns == 42
    assert captured["membrane_config"].lipid_ff == "lipid17"
    assert captured["membrane_config"].protein_forcefield == "amber14sb"
    assert captured["membrane_config"].water_model == "tip3p"
    assert response["success"] is False
    assert response["partial"] is True
    assert response["gmx_dir"].endswith("GMX_PROLIG_MEMB")
    assert "localrun.sh" not in response["message"]
    assert response["warnings"] == ["packmol-memgen not found"]
    assert response["parameters"]["lipid_ff"] == "lipid17"
    assert not {"system_coordinates", "topology", "index"} & response["files"].keys()
    assert set(response["files"]) == {"em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"}



def test_mcp_membrane_build_rejects_silently_ignored_ligand(tmp_path):
    protein = tmp_path / "protein.pdb"
    ligand = tmp_path / "ligand.sdf"
    protein.write_text("ATOM\n")
    ligand.write_text("ligand\n")
    registry = _ToolRegistry()
    mcp_build.register(registry)

    response = json.loads(
        registry.tools["build_system"](
            protein_path=str(protein.resolve()),
            ligand_paths=str(ligand.resolve()),
            membrane=True,
        )
    )

    assert response["success"] is False
    assert any("not yet incorporated" in error for error in response["errors"])


def test_cli_membrane_build_allows_protein_only(monkeypatch, tmp_path):
    captured = {}

    class FakeBuilder:
        def __init__(self, protein_path, ligand_paths, output_dir, **kwargs):
            captured.update(
                protein_path=protein_path,
                ligand_paths=ligand_paths,
                output_dir=output_dir,
                **kwargs,
            )

        def run(self):
            return captured["output_dir"]

    monkeypatch.setattr(builder_cli, "PRISMBuilder", FakeBuilder)
    monkeypatch.setattr(
        sys,
        "argv",
        ["prism", str(tmp_path / "protein.pdb"), "--membrane", "-o", str(tmp_path / "output")],
    )

    builder_cli.main()

    assert captured["ligand_paths"] == []
    assert captured["forcefield"] == "amber14sb"
    assert captured["water_model"] == "tip3p"
    assert captured["membrane_config"].enabled is True
    assert captured["membrane_config"].production_ns == 500


def test_cli_allows_membrane_enabled_only_in_yaml_without_a_ligand(monkeypatch, tmp_path):
    captured = {}
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM\n")
    config = tmp_path / "prism.yaml"
    config.write_text("membrane:\n  enabled: true\n", encoding="utf-8")

    class FakeBuilder:
        def __init__(self, protein_path, ligand_paths, output_dir, **kwargs):
            captured.update(
                protein_path=protein_path,
                ligand_paths=ligand_paths,
                output_dir=output_dir,
                **kwargs,
            )

        def run(self):
            return captured["output_dir"]

    monkeypatch.setattr(builder_cli, "PRISMBuilder", FakeBuilder)
    monkeypatch.setattr(
        sys,
        "argv",
        ["prism", str(protein), "--config", str(config), "-o", str(tmp_path / "output")],
    )

    builder_cli.main()

    assert captured["ligand_paths"] == []
    assert captured["config_path"] == str(config)
    assert captured["membrane_config"].enabled is True
    assert captured["membrane_config"].protein_forcefield == "amber14sb"
    assert captured["membrane_config"].water_model == "tip3p"


def test_cli_membrane_override_merges_into_an_empty_yaml_file(monkeypatch, tmp_path):
    captured = {}
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM\n")
    config = tmp_path / "empty.yaml"
    config.write_text("", encoding="utf-8")

    class FakeBuilder:
        def __init__(self, protein_path, ligand_paths, output_dir, **kwargs):
            import yaml

            captured.update(
                protein_path=protein_path,
                ligand_paths=ligand_paths,
                output_dir=output_dir,
                **kwargs,
            )
            with open(kwargs["config_path"], "r", encoding="utf-8") as config_handle:
                captured["effective_yaml"] = yaml.safe_load(config_handle)

        def run(self):
            return captured["output_dir"]

    monkeypatch.setattr(builder_cli, "PRISMBuilder", FakeBuilder)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prism",
            str(protein),
            "--config",
            str(config),
            "--membrane",
            "--temperature",
            "310",
            "-o",
            str(tmp_path / "output"),
        ],
    )

    builder_cli.main()

    assert captured["membrane_config"].temperature == 310
    assert captured["effective_yaml"]["simulation"]["temperature"] == 310


def test_cli_yaml_membrane_values_survive_argparse_defaults(monkeypatch, tmp_path):
    captured = {}
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM\n")
    config = tmp_path / "prism.yaml"
    config.write_text(
        """forcefield:
  name: amber19sb
water_model:
  name: opc
simulation:
  temperature: 315
  production_time_ns: 25
ions:
  concentration: 0.2
  positive_ion: K+
  negative_ion: Cl-
membrane:
  enabled: true
  lipids: [POPE]
  ratio: [2]
  orient: opm
  pdb_id: 4xb4
  lipid_ff: lipid17
  water_thickness: 22
  minimize: false
""",
        encoding="utf-8",
    )

    class FakeBuilder:
        def __init__(self, protein_path, ligand_paths, output_dir, **kwargs):
            captured.update(
                protein_path=protein_path,
                ligand_paths=ligand_paths,
                output_dir=output_dir,
                **kwargs,
            )

        def run(self):
            return captured["output_dir"]

    monkeypatch.setattr(builder_cli, "PRISMBuilder", FakeBuilder)
    monkeypatch.setattr(
        sys,
        "argv",
        ["prism", str(protein), "--config", str(config), "-o", str(tmp_path / "output")],
    )

    builder_cli.main()

    membrane = captured["membrane_config"]
    assert captured["forcefield"] == "amber19sb"
    assert captured["water_model"] == "opc"
    assert membrane.lipids == ["POPE"]
    assert membrane.ratio == [2]
    assert membrane.orient == "opm"
    assert membrane.pdb_id == "4xb4"
    assert membrane.lipid_ff == "lipid17"
    assert membrane.temperature == 315
    assert membrane.production_ns == 25
    assert membrane.salt_concentration == 0.2
    assert membrane.positive_ion == "K+"
    assert membrane.water_thickness == 22
    assert membrane.minimize is False


def test_cli_explicit_membrane_options_override_yaml_even_at_parser_defaults(monkeypatch, tmp_path):
    captured = {}
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM\n")
    config = tmp_path / "prism.yaml"
    config.write_text(
        """forcefield:
  name: amber19sb
water_model:
  name: opc
simulation:
  temperature: 315
  production_time_ns: 25
membrane:
  enabled: true
  lipids: [POPE]
  ratio: [2]
  orient: opm
  pdb_id: 4xb4
  lipid_ff: lipid17
  water_thickness: 22
  minimize: false
""",
        encoding="utf-8",
    )

    class FakeBuilder:
        def __init__(self, protein_path, ligand_paths, output_dir, **kwargs):
            import yaml

            captured.update(
                protein_path=protein_path,
                ligand_paths=ligand_paths,
                output_dir=output_dir,
                **kwargs,
            )
            with open(kwargs["config_path"], "r", encoding="utf-8") as config_handle:
                captured["effective_yaml"] = yaml.safe_load(config_handle)

        def run(self):
            return captured["output_dir"]

    monkeypatch.setattr(builder_cli, "PRISMBuilder", FakeBuilder)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prism",
            str(protein),
            "--config",
            str(config),
            "--membrane",
            "--forcefield",
            "amber14sb",
            "--water",
            "tip3p",
            "--temperature",
            "310",
            "--production-ns",
            "500",
            "--membrane-orient",
            "memembed",
            "--lipid-ff",
            "lipid21",
            "--lipid",
            "POPC",
            "--lipid-ratio",
            "3",
            "-o",
            str(tmp_path / "output"),
        ],
    )

    builder_cli.main()

    membrane = captured["membrane_config"]
    assert captured["forcefield"] == "amber14sb"
    assert captured["water_model"] == "tip3p"
    assert membrane.lipids == ["POPC"]
    assert membrane.ratio == [3]
    assert membrane.orient == "memembed"
    assert membrane.lipid_ff == "lipid21"
    assert membrane.temperature == 310
    assert membrane.production_ns == 500
    assert captured["effective_yaml"]["simulation"]["temperature"] == 310
    assert captured["effective_yaml"]["simulation"]["production_time_ns"] == 500
    assert captured["effective_yaml"]["forcefield"]["name"] == "amber14sb"
    assert captured["effective_yaml"]["water_model"]["name"] == "tip3p"
    # Unoverridden YAML-only fields are preserved.
    assert membrane.water_thickness == 22
    assert membrane.minimize is False




def test_cli_membrane_build_rejects_silently_ignored_ligand(monkeypatch, tmp_path):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prism",
            str(tmp_path / "protein.pdb"),
            str(tmp_path / "ligand.sdf"),
            "--membrane",
        ],
    )

    with pytest.raises(SystemExit):
        builder_cli.main()
