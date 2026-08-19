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
from prism.membrane.builder import (
    MEMBRANE_SYSTEM_DIRNAME,
    MembraneBuilder,
    MembraneBuildResult,
)
from prism.membrane.config import (
    DEFAULT_PRODUCTION_NS,
    MEMBRANE_EQUIL_TIMESTEP_PS,
    MEMBRANE_TIMESTEP_PS,
    MembraneConfig,
    gaff2_selected,
    parse_membrane_cli,
    parse_membrane_yaml,
)
from prism.membrane.membrane_mdp import MEMBRANE_MDP_DIRNAME, write_membrane_mdps
from prism.membrane.orient import (
    fetch_opm,
    insert_ter_at_chain_breaks,
    membrane_alignment_report,
    orient_protein,
    validate_membrane_alignment,
)
from prism.membrane import builder as membrane_builder
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


#: The staged membrane protocol, in execution order. Equilibration releases
#: position restraints across nvt -> npt -> npt2 before md runs unrestrained.
_MEMBRANE_STAGES = ("em", "nvt", "npt", "npt2", "md")

#: The equilibration stages that carry a ``define`` line, in release order.
_RESTRAINED_STAGES = ("nvt", "npt", "npt2")


def _option_value(command, option):
    return command[command.index(option) + 1]


def _mdp_option(text, key):
    """Return the value of ``key`` in an .mdp file's text, or None."""
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if "=" not in line:
            continue
        name, value = line.split("=", 1)
        if name.strip().lower().replace("_", "-") == key:
            return value.strip()
    return None


def _define_macros(text):
    """Parse an .mdp ``define`` line into ``{macro: value}``.

    A bare guard such as ``-DPOSRES`` maps to ``None``; ``-DPOSRES_FC_BB=1000.0``
    maps to the string ``"1000.0"``. Returns ``{}`` when there is no define line,
    which is how an unrestrained stage is distinguished from a restrained one.
    """
    define = _mdp_option(text, "define")
    if define is None:
        return {}
    macros = {}
    for token in define.split():
        assert token.startswith("-D"), f"unexpected define token {token!r}"
        macro, separator, value = token[2:].partition("=")
        macros[macro] = value if separator else None
    return macros


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


@pytest.mark.parametrize("supplied", ["ligand_frcmod", "ligand_lib"])
def test_membrane_config_rejects_a_half_supplied_ligand_parameter_pair(supplied):
    """PACKMOL-Memgen's ``--ligand_param`` takes both files or neither.

    Half a pair is the dangerous case: ``has_ligand()`` goes False, the ligand
    flags are never emitted, and the bilayer is built around a protein whose
    ligand has quietly vanished. Fail loudly instead.
    """
    half = MembraneConfig(enabled=True, **{supplied: f"/params/lig.{supplied}"})

    assert half.has_ligand() is False
    with pytest.raises(ValueError, match=f"got only {supplied}"):
        half.raise_for_errors()

    matched = MembraneConfig(
        enabled=True,
        ligand_frcmod="/params/lig.frcmod",
        ligand_lib="/params/lig.lib",
    )
    matched.raise_for_errors()
    assert matched.has_ligand() is True


def test_build_command_emits_ligand_options_only_when_a_ligand_is_configured(tmp_path):
    """Ligand flags appear as a set, and a ligand-free build is untouched."""
    frcmod = tmp_path / "lig.frcmod"
    lib = tmp_path / "lig.lib"

    baseline = packmol_memgen.build_command("protein.pdb", MembraneConfig(enabled=True))

    assert "--keepligs" not in baseline
    assert "--ligand_param" not in baseline
    assert "--gaff2" not in baseline

    # --gaff2 selects the leaprc PACKMOL-Memgen sources for the ligand, so it
    # tracks the GAFF generation that actually typed frcmod/lib rather than
    # riding along with every ligand.  A mismatch is silent (GAFF's atom types
    # are a subset of GAFF2's, same masses) and only corrupts force constants,
    # so both generations are pinned here.
    for gaff2, expected in ((False, False), (True, True)):
        command = packmol_memgen.build_command(
            "protein.pdb",
            MembraneConfig(
                enabled=True,
                ligand_frcmod=str(frcmod),
                ligand_lib=str(lib),
                gaff2=gaff2,
            ),
        )

        assert "--keepligs" in command
        assert _option_value(command, "--ligand_param") == f"{frcmod}:{lib}"
        assert ("--gaff2" in command) is expected

        # Configuring a ligand adds exactly the ligand options and perturbs
        # nothing else about the packing/parametrization request.
        ligand_tokens = {
            "--keepligs",
            "--gaff2",
            "--ligand_param",
            _option_value(command, "--ligand_param"),
        }
        assert [token for token in command if token not in ligand_tokens] == baseline

    # The standalone default follows PRISM's own ``-lff`` default of plain GAFF,
    # which is also PACKMOL-Memgen's documented leaprc default.
    assert MembraneConfig(enabled=True).gaff2 is False
    assert gaff2_selected("gaff2") is True
    assert gaff2_selected("gaff") is False
    assert gaff2_selected("openff") is False
    assert gaff2_selected(None) is False


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


def test_membrane_mdps_emit_every_stage_of_the_staged_protocol(tmp_path):
    files = write_membrane_mdps(str(tmp_path))

    assert set(files) == set(_MEMBRANE_STAGES)
    for stage in _MEMBRANE_STAGES:
        assert Path(files[stage]).name == f"{stage}.mdp"
        assert Path(files[stage]).is_file()


def test_membrane_restraints_decay_monotonically_and_production_is_unrestrained(tmp_path):
    """Restraints are released gradually; production carries none.

    Dropping restraints in one step re-introduces exactly the shock the staged
    protocol exists to avoid, and any restraint surviving into production would
    bias the run everything else exists to set up.
    """
    files = write_membrane_mdps(str(tmp_path))
    macros = {
        stage: _define_macros(Path(files[stage]).read_text())
        for stage in _MEMBRANE_STAGES
    }

    # Minimisation is deliberately unrestrained: restraining the very atoms that
    # are clashing turns a fixable packing clash into a LINCS failure later.
    assert macros["em"] == {}
    assert macros["md"] == {}

    for name in ("POSRES_FC_BB", "POSRES_FC_SC"):
        values = [float(macros[stage][name]) for stage in _RESTRAINED_STAGES]
        assert values[0] > values[1] > values[2], f"{name} does not decay: {values}"
        assert values[0] > 0

    # Lipids are freed one stage before the protein backbone, so their constant
    # decays and then the whole restraint block is switched off.
    lipid = [macros[stage].get("POSRES_FC_LIPID") for stage in _RESTRAINED_STAGES]
    assert float(lipid[0]) > float(lipid[1]) > 0
    assert lipid[2] is None
    assert "POSRES_LIPID" not in macros["npt2"]

    # The final equilibration stage is backbone-only: side chains are released
    # to zero while the fold itself is still held.
    assert float(macros["npt2"]["POSRES_FC_SC"]) == 0.0
    assert float(macros["npt2"]["POSRES_FC_BB"]) > 0.0


def test_membrane_define_line_pairs_every_force_constant_with_its_guard(tmp_path):
    """An unpaired macro is a grompp failure, not a cosmetic inconsistency.

    The restraint .itp files reference ``POSRES_FC_BB``/``POSRES_FC_SC`` inside
    ``#ifdef POSRES`` and ``POSRES_FC_LIPID`` inside ``#ifdef POSRES_LIPID``.
    Emitting a guard without its force constant leaves an unsubstituted token in
    the topology and grompp dies; emitting a force constant without its guard
    silently restrains nothing. So the two must always travel together.
    """
    files = write_membrane_mdps(str(tmp_path))
    guards = {
        "POSRES": ("POSRES_FC_BB", "POSRES_FC_SC"),
        "POSRES_LIPID": ("POSRES_FC_LIPID",),
    }
    constants = {name for group in guards.values() for name in group}

    for stage in _MEMBRANE_STAGES:
        macros = _define_macros(Path(files[stage]).read_text())
        for guard, group in guards.items():
            for constant in group:
                assert (guard in macros) == (constant in macros), (
                    f"{stage}.mdp pairs {guard} with {constant} inconsistently: "
                    f"{sorted(macros)}"
                )
        # A typo'd macro name would restrain nothing and never be noticed, so no
        # macro outside the known contract may appear.
        assert set(macros) <= set(guards) | constants, f"{stage}.mdp: {sorted(macros)}"
        for macro, value in macros.items():
            if macro in guards:
                assert value is None, f"guard {macro} must be bare, got {value!r}"
            else:
                assert value is not None and float(value) >= 0


def test_membrane_equilibration_starts_at_1fs_and_production_runs_at_2fs(tmp_path):
    """A freshly packed bilayer handed straight to 2 fs dynamics is a classic
    LINCS failure, so the restrained stages integrate at 1 fs and only step up
    once the restraints that make 2 fs safe have done their job.
    """
    assert MEMBRANE_EQUIL_TIMESTEP_PS == 0.001
    assert MEMBRANE_TIMESTEP_PS == 0.002

    files = write_membrane_mdps(str(tmp_path), production_ns=1.0)

    def option(stage, key):
        return _mdp_option(Path(files[stage]).read_text(), key)

    assert float(option("nvt", "dt")) == MEMBRANE_EQUIL_TIMESTEP_PS
    assert float(option("npt", "dt")) == MEMBRANE_EQUIL_TIMESTEP_PS
    assert float(option("npt2", "dt")) == MEMBRANE_TIMESTEP_PS
    assert float(option("md", "dt")) == MEMBRANE_TIMESTEP_PS

    # Minimisation does not integrate at all, so it carries no timestep.
    assert option("em", "integrator") == "steep"
    assert option("em", "dt") is None

    # Steps are derived from a physical duration, so the smaller equilibration
    # timestep preserves each stage's wall-clock length instead of halving it.
    def duration_ps(stage):
        return int(option(stage, "nsteps")) * float(option(stage, "dt"))

    assert duration_ps("nvt") == pytest.approx(250.0)
    assert duration_ps("npt") == pytest.approx(500.0)
    # npt2 is the stage that lets area per lipid converge before Parrinello-Rahman
    # takes over, so it is deliberately the longest equilibration stage.
    assert duration_ps("npt2") == pytest.approx(5000.0)
    assert duration_ps("md") == pytest.approx(1000.0)  # the requested 1 ns

    # The barostat runs for substantially longer than the volume relaxation time,
    # because area per lipid is the slow variable and production starts on a
    # barostat with no dissipation.
    assert duration_ps("npt") + duration_ps("npt2") >= 5000.0


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


def _write_z_stack_pdb(path, z_values, resname="ALA"):
    """Write a PDB whose protein atoms sit at the given z coordinates."""
    lines = []
    for serial, z in enumerate(z_values, start=1):
        lines.append(
            f"ATOM  {serial:5d} {'CA':^4s} {resname:>3s} A{serial:4d}    "
            f"{0.0:8.3f}{0.0:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {'C':>2s}\n"
        )
    path.write_text("".join(lines) + "END\n", encoding="ascii")
    return path


#: A transmembrane segment in the frame PACKMOL-Memgen builds for: it crosses
#: z = 0 and spans far enough to be judged by the two-sided leaflet test.
_TRANSMEMBRANE_Z = [float(value) for value in range(-20, 21)]


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


def test_validate_membrane_alignment_rejects_a_structure_off_the_midplane(tmp_path):
    """PACKMOL-Memgen builds the bilayer at z = 0 and leaves the solute where it
    is, so an off-midplane input packs into bulk water and "succeeds" as a
    scientifically worthless system. It must be refused instead.
    """
    offset = _write_z_stack_pdb(
        tmp_path / "offset.pdb", [60.0 + value for value in _TRANSMEMBRANE_Z]
    )

    # The report measures and describes; only the validator enforces.
    report = membrane_alignment_report(str(offset))
    assert report["measurable"] is True
    assert report["aligned"] is False

    with pytest.raises(ValueError, match="not aligned to the membrane frame"):
        validate_membrane_alignment(str(offset))


def test_validate_membrane_alignment_rejects_a_protein_resting_on_one_leaflet(tmp_path):
    """Spans far enough to be judged by the two-sided test, but every atom near
    the midplane is above it: this protein sits ON the bilayer, not through it.
    """
    grazing = _write_z_stack_pdb(
        tmp_path / "grazing.pdb", [25.0 + value for value in _TRANSMEMBRANE_Z]
    )

    with pytest.raises(ValueError, match="one-sided"):
        validate_membrane_alignment(str(grazing))


def test_validate_membrane_alignment_accepts_a_centred_membrane_protein(tmp_path):
    centred = _write_z_stack_pdb(tmp_path / "centred.pdb", _TRANSMEMBRANE_Z)

    report = validate_membrane_alignment(str(centred))

    assert report["aligned"] is True
    assert report["measurable"] is True
    assert report["straddles_midplane"] is True


def test_validate_membrane_alignment_accepts_a_large_one_sided_extramembrane_domain(tmp_path):
    """The false positive that actually matters.

    Real membrane proteins carry big soluble domains on one side (a GPCR
    ectodomain can be 100 A tall). A whole-molecule centre test would drag this
    structure's centre far above the midplane and reject it, even though its
    membrane-spanning region is exactly where the bilayer will be built. The
    verdict must therefore be taken over the membrane slab, not the molecule.
    """
    ectodomain = [30.0 + 0.5 * index for index in range(180)]
    structure = _write_z_stack_pdb(
        tmp_path / "ectodomain.pdb", _TRANSMEMBRANE_Z + ectodomain
    )

    report = validate_membrane_alignment(str(structure))

    assert report["aligned"] is True
    assert report["straddles_midplane"] is True
    assert report["slab_median_z"] == pytest.approx(0.0)
    # The naive whole-molecule measure is far outside tolerance, which is
    # precisely why the verdict must not be based on it.
    assert abs(report["median_z"]) > report["tolerance"]


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
        SimpleNamespace(residue=SimpleNamespace(name="ZZZ")),
    ]
    index_path = tmp_path / "index.ndx"

    # ZZZ, not MSE: selenomethionine is a recognised protein residue and must
    # classify -- it appears in ordinary crystal structures, and refusing it
    # used to kill a build at index generation, after packing and solvation had
    # already run. The property under test is the fail-closed one, so the sample
    # has to be a name nothing recognises.
    with pytest.raises(ValueError, match="ZZZ"):
        builder._write_membrane_index(SimpleNamespace(atoms=atoms), str(index_path))
    assert not index_path.exists()

    # The residue that motivated the vocabulary fix now classifies instead.
    atoms[3] = SimpleNamespace(residue=SimpleNamespace(name="MSE"))
    builder._write_membrane_index(SimpleNamespace(atoms=atoms), str(index_path))
    assert index_path.exists()


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

    # A writer that silently produces nothing is the failure this guards against,
    # so stub the WRITERS rather than ParmEd itself. Faking the parmed module
    # instead would exercise the topology machinery against an object that is not
    # a Structure, and fail for that reason long before the staleness check --
    # testing the stub, not the guarantee.
    class NoOpParm:
        atoms = []
        dihedrals = []

        def save(self, *_args, **_kwargs):
            return None

    monkeypatch.setitem(
        sys.modules,
        "parmed",
        SimpleNamespace(load_file=lambda *_args, **_kwargs: NoOpParm()),
    )
    monkeypatch.setattr(
        membrane_builder,
        "save_gromacs_topology",
        lambda *_args, **_kwargs: SimpleNamespace(describe=lambda: "stubbed"),
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


def test_membrane_system_directory_is_named_for_a_full_md_system(monkeypatch, tmp_path):
    """The membrane system directory mirrors the soluble ``GMX_PROLIG_MD``.

    Both workflows now hand back the same shape of deliverable (topology,
    coordinates, index, restraints, run script), so the name says "a runnable MD
    system" rather than naming a membrane scaffold.
    """
    assert MEMBRANE_SYSTEM_DIRNAME == "GMX_PROLIG_MEMB"

    input_pdb = tmp_path / "protein.pdb"
    input_pdb.write_text("ATOM\n")
    builder = MembraneBuilder(MembraneConfig(enabled=True), verbose=False)
    monkeypatch.setattr(builder, "missing_capabilities", lambda: ["AmberTools unavailable"])
    output = tmp_path / "output"

    result = builder.build(str(input_pdb), str(output))

    assert Path(result.system_dir) == output / MEMBRANE_SYSTEM_DIRNAME
    assert Path(result.system_dir).is_dir()
    # Exactly one system directory, under the current name: a rename that left
    # the old directory behind would show up here as a stray sibling.  The MDP
    # directory is deliberately *not* the soluble path's ``mdps``: these MDPs
    # name membrane thermostat groups and are unusable without ``-n index.ndx``
    # (see MEMBRANE_MDP_DIRNAME), so a generic consumer must not find them there.
    assert {entry.name for entry in output.iterdir() if entry.is_dir()} == {
        MEMBRANE_SYSTEM_DIRNAME,
        MEMBRANE_MDP_DIRNAME,
    }
    assert MEMBRANE_MDP_DIRNAME != "mdps"
    # A capability-gated scaffold is not a runnable system, so it carries
    # neither restraints nor a driver script.
    assert result.posres == {}
    assert result.runscript is None
    assert result.success is False


def test_membrane_builder_scaffold_uses_configured_production_time(monkeypatch, tmp_path):
    input_pdb = tmp_path / "protein.pdb"
    input_pdb.write_text("ATOM\n")
    config = MembraneConfig(enabled=True, production_ns=3.0)
    builder = MembraneBuilder(config, verbose=False)
    monkeypatch.setattr(builder, "missing_capabilities", lambda: ["AmberTools unavailable"])

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
            system_dir = output / MEMBRANE_SYSTEM_DIRNAME
            system_dir.mkdir(parents=True)
            # Simulate complete files left by an older run. A partial response
            # must not present them as artifacts of the current invocation --
            # including the run script and restraints, which are what a caller
            # would otherwise pick up and try to execute.
            (system_dir / "system.gro").write_text("stale")
            (system_dir / "topol.top").write_text("stale")
            (system_dir / "index.ndx").write_text("stale")
            (system_dir / "localrun.sh").write_text("stale")
            (system_dir / "posre_solu.itp").write_text("stale")
            mdps = output / "mdps"
            mdps.mkdir()
            for name in _MEMBRANE_STAGES:
                (mdps / f"{name}.mdp").write_text("scaffold")
            self.membrane_result = SimpleNamespace(
                success=False,
                gro=None,
                top=None,
                ndx=None,
                posres={},
                runscript=None,
                mdps={stage: str(mdps / f"{stage}.mdp") for stage in _MEMBRANE_STAGES},
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
    assert response["gmx_dir"].endswith(MEMBRANE_SYSTEM_DIRNAME)
    assert "localrun.sh" not in response["message"]
    assert response["warnings"] == ["packmol-memgen not found"]
    assert response["parameters"]["lipid_ff"] == "lipid17"
    assert not {
        "system_coordinates", "topology", "index", "run_script"
    } & response["files"].keys()
    assert set(response["files"]) == {f"{stage}.mdp" for stage in _MEMBRANE_STAGES}



def test_mcp_membrane_build_reports_the_run_script_and_restraints(monkeypatch, tmp_path):
    """A complete membrane build reports the same shape as the soluble path.

    The membrane system is runnable now, so the caller is told how to run it
    (``localrun.sh``) and handed the restraint files grompp needs, instead of a
    bare list of MDPs and a hand-assembled grompp invocation.
    """
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM\n")
    restraints = ("posre_solu.itp", "posre_memb.itp")

    class FakeBuilder:
        def __init__(self, **kwargs):
            self.output_dir = kwargs["output_dir"]
            self.config = {
                "simulation": {},
                "ligand_forcefield": {},
                "box": {},
                "ions": {},
                "energy_minimization": {},
            }

        def run(self):
            root = Path(self.output_dir)
            system_dir = root / MEMBRANE_SYSTEM_DIRNAME
            system_dir.mkdir(parents=True)
            mdp_dir = root / "mdps"
            mdp_dir.mkdir()
            for stage in _MEMBRANE_STAGES:
                (mdp_dir / f"{stage}.mdp").write_text("stage")
            for name in ("system.gro", "topol.top", "index.ndx", "localrun.sh", *restraints):
                (system_dir / name).write_text("built")
            self.membrane_result = SimpleNamespace(
                success=True,
                gro=str(system_dir / "system.gro"),
                top=str(system_dir / "topol.top"),
                ndx=str(system_dir / "index.ndx"),
                oriented_pdb=None,
                posres={name: str(system_dir / name) for name in restraints},
                runscript=str(system_dir / "localrun.sh"),
                mdps={stage: str(mdp_dir / f"{stage}.mdp") for stage in _MEMBRANE_STAGES},
                note="Membrane system built.",
                warnings=[],
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
        )
    )

    assert response["success"] is True
    assert response["partial"] is False
    assert response["gmx_dir"].endswith(MEMBRANE_SYSTEM_DIRNAME)
    assert Path(response["files"]["run_script"]).name == "localrun.sh"
    assert set(restraints) <= response["files"].keys()
    assert {f"{stage}.mdp" for stage in _MEMBRANE_STAGES} <= response["files"].keys()
    assert {"system_coordinates", "topology", "index"} <= response["files"].keys()
    # The soluble path tells the user to run localrun.sh; so must this one.
    assert "bash localrun.sh" in response["message"]


def test_mcp_membrane_build_accepts_a_ligand(monkeypatch, tmp_path):
    """Membrane mode now embeds a ligand, so the tool must not pre-reject one.

    The rejection used to live here *and* in ``prism.builder.core``; the
    duplicate was the one that ran first, so a membrane+ligand request never
    reached the builder at all. This asserts the request is forwarded instead.
    """
    protein = tmp_path / "protein.pdb"
    ligand = tmp_path / "ligand.sdf"
    protein.write_text("ATOM\n")
    ligand.write_text("ligand\n")
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
            self.membrane_result = SimpleNamespace(
                success=False,
                gro=None,
                top=None,
                ndx=None,
                posres={},
                runscript=None,
                mdps={},
                note="Membrane build prerequisites missing.",
                warnings=[],
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
            ligand_paths=str(ligand.resolve()),
            output_dir=str(tmp_path / "output"),
            membrane=True,
        )
    )

    # The ligand reached the builder rather than being refused up front.
    assert captured["ligand_paths"] == str(ligand.resolve())
    assert captured["membrane_config"].enabled is True
    # No validation error at all, and specifically not the retired message.
    assert "errors" not in response
    assert "not yet incorporated" not in json.dumps(response)


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




def test_cli_membrane_build_rejects_multiple_ligands(monkeypatch, tmp_path, capsys):
    """PACKMOL-Memgen's ``--ligand_param`` takes one FRCMOD:LIB pair.

    Accepting several ligands could only mean silently dropping all but one, so
    this fails at argument-parse time rather than minutes into packing.
    """
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prism",
            "-pf", str(tmp_path / "protein.pdb"),
            "-lf", str(tmp_path / "ligand1.mol2"),
            "-lf", str(tmp_path / "ligand2.mol2"),
            "--membrane",
        ],
    )

    with pytest.raises(SystemExit):
        builder_cli.main()

    message = capsys.readouterr().err
    assert "single ligand" in message
    assert "2 were given" in message
    assert "not yet incorporated" not in message


def test_cli_membrane_build_rejects_non_amber_ligand_forcefield(monkeypatch, tmp_path, capsys):
    """The bilayer is parameterized in AMBER, so the ligand must be too.

    OpenFF/OPLS/CGenFF generators emit GROMACS-side parameters only; there is no
    frcmod/lib to hand PACKMOL-Memgen, and tleap would happily build a
    ligand-free system without complaining.
    """
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prism",
            "-pf", str(tmp_path / "protein.pdb"),
            "-lf", str(tmp_path / "ligand.mol2"),
            "--membrane",
            "--ligand-forcefield", "openff",
        ],
    )

    with pytest.raises(SystemExit):
        builder_cli.main()

    message = capsys.readouterr().err
    assert "openff" in message
    assert "--ligand-forcefield gaff" in message


# --------------------------------------------------------------------------- #
# Regression guards for silent-failure modes found in adversarial review.
# Each of these was a real bug: grompp accepts every one of the broken states
# below without a diagnostic, so only PRISM can catch them.
# --------------------------------------------------------------------------- #

_TWO_BLOCK_TOPOLOGY = """\
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
; Name            nrexcl
system1          3

[ atoms ]
;   nr  type resnr residue atom cgnr charge mass
     1     N     1     MET     N    1  0.1  14.010000
     2    CT     1     MET    CA    2  0.0  12.010000
     3     C     1     MET     C    3  0.5  12.010000
     4     O     1     MET     O    4 -0.5  16.000000
     5    CT     1     MET    CB    5  0.0  12.010000
     6    HC     1     MET   HB1    6  0.0   1.008000

[ moleculetype ]
; Name            nrexcl
system2          3

[ atoms ]
;   nr  type resnr residue atom cgnr charge mass
     1    cD     1      PA  C116    1  0.0  12.010000
     2    pA     2      PC   P31    2  1.2  30.970000

[ system ]
membrane

[ molecules ]
system1 1
system2 4
"""


def _membrane_builder_with_topology(tmp_path):
    """A builder plus a two-moleculetype topology (one protein, one POPC)."""
    top = tmp_path / "topol.top"
    top.write_text(_TWO_BLOCK_TOPOLOGY)
    builder = MembraneBuilder(
        MembraneConfig(enabled=True, lipids=["POPC"]), verbose=False
    )
    return builder, top


def test_position_restraints_classify_backbone_sidechain_and_phosphorus(tmp_path):
    builder, top = _membrane_builder_with_topology(tmp_path)
    result = MembraneBuildResult(system_dir=str(tmp_path))

    builder._write_position_restraints(str(top), str(tmp_path), result)

    assert result.warnings == []

    def _rows(path):
        body = Path(path).read_text().split("[ position_restraints ]")[1]
        return [
            line.split()
            for line in body.strip().splitlines()
            if line.strip() and not line.lstrip().startswith(";")
        ]

    solu = _rows(tmp_path / "posre_solu.itp")
    # N/CA/C/O are backbone; CB is sidechain; HB1 is hydrogen and is skipped
    # because LINCS already ties it to its heavy atom.
    assert [row[0] for row in solu] == ["1", "2", "3", "4", "5"]
    assert [row[2] for row in solu] == ["POSRES_FC_BB"] * 4 + ["POSRES_FC_SC"]
    # every row restrains all three dimensions with the same constant
    assert all(row[2] == row[3] == row[4] for row in solu)

    rows = _rows(tmp_path / "posre_memb.itp")
    assert len(rows) == 1, "only the headgroup phosphorus should be restrained"
    assert rows[0][0] == "2"
    # z-only: restraining x/y would freeze lateral diffusion and stop the
    # bilayer relaxing its area per lipid.
    assert rows[0][2:] == ["0", "0", "POSRES_FC_LIPID"]


def test_position_restraints_refuse_second_injection(tmp_path):
    """Re-injecting would ADD a section, multiplying every force constant.

    grompp accepts a doubled [ position_restraints ] block silently, so the
    only place this can be caught is here.
    """
    builder, top = _membrane_builder_with_topology(tmp_path)

    first = MembraneBuildResult(system_dir=str(tmp_path))
    builder._write_position_restraints(str(top), str(tmp_path), first)
    assert first.warnings == []

    second = MembraneBuildResult(system_dir=str(tmp_path))
    builder._write_position_restraints(str(top), str(tmp_path), second)

    assert second.warnings, "a second injection must be refused"
    assert "already includes" in second.warnings[0]
    assert top.read_text().count('#include "posre_solu.itp"') == 1


def test_position_restraints_reject_unparsable_atom_columns(tmp_path):
    """A protein block yielding no restrainable atom means we misparsed it.

    Failing closed here is what stops an unrestrained run being reported as
    restrained.
    """
    builder, _ = _membrane_builder_with_topology(tmp_path)
    top = tmp_path / "short.top"
    # The optional mass column is dropped from every protein atom row. That is
    # legal GROMACS, but it leaves the classifier with nothing to weigh, and
    # silently skipping the block would produce an unrestrained protein.
    top.write_text(
        "[ moleculetype ]\n"
        "; Name            nrexcl\n"
        "system1          3\n"
        "\n"
        "[ atoms ]\n"
        ";   nr  type resnr residue atom cgnr charge\n"
        "     1     N     1     MET     N    1  0.1\n"
        "     2    CT     1     MET    CA    2  0.0\n"
        "\n"
        "[ system ]\n"
        "membrane\n"
    )
    result = MembraneBuildResult(system_dir=str(tmp_path))

    builder._write_position_restraints(str(top), str(tmp_path), result)

    assert result.warnings
    assert not (tmp_path / "posre_solu.itp").exists()


def test_packing_skip_is_detected_as_failure():
    """PACKMOL-Memgen exits 0 when it reuses an existing bilayer."""
    assert packmol_memgen._packing_was_skipped(
        "Packed PDB bilayer_protein_oriented.pdb found. Skipping PACKMOL"
    )
    assert not packmol_memgen._packing_was_skipped("Building system...\nDone.")


def test_ligand_absence_in_packed_system_is_detected(tmp_path):
    """A ligand dropped during packing yields a valid but wrong bilayer."""
    packed = tmp_path / "bilayer.pdb"
    packed.write_text(
        "ATOM      1  N   MET A   1       0.000   0.000   0.000\n"
        "ATOM      2  P31 PC      2       0.000   0.000  18.000\n"
    )
    assert packmol_memgen._ligand_absent_from(str(packed), "LIG")

    packed.write_text(
        "ATOM      1  N   MET A   1       0.000   0.000   0.000\n"
        "HETATM    2  C1  LIG     2       1.000   0.000   0.000\n"
    )
    assert not packmol_memgen._ligand_absent_from(str(packed), "LIG")


def test_sterols_are_restrained_by_their_hydroxyl_oxygen(tmp_path):
    """Cholesterol has no phosphate, so a P-only rule leaves it free.

    In a POPC:CHL1 bilayer that would hold the phospholipids while the sterol
    fraction moves unrestrained -- the bilayer deforms around it during
    equilibration, and nothing in GROMACS reports it.
    """
    top = tmp_path / "topol.top"
    top.write_text(
        _TWO_BLOCK_TOPOLOGY.replace(
            "[ system ]",
            "[ moleculetype ]\n"
            "; Name            nrexcl\n"
            "system3          3\n"
            "\n"
            "[ atoms ]\n"
            "     1    cA     1     CHL    C1    1  0.0  12.010000\n"
            "     2    oH     1     CHL    O1    2 -0.6  16.000000\n"
            "     3    hO     1     CHL    H1    3  0.4   1.008000\n"
            "\n"
            "[ system ]",
        )
    )
    builder = MembraneBuilder(
        MembraneConfig(enabled=True, lipids=["POPC", "CHL1"], ratio=[9, 1]), verbose=False
    )
    result = MembraneBuildResult(system_dir=str(tmp_path))

    builder._write_position_restraints(str(top), str(tmp_path), result)

    assert result.warnings == []
    # One restraint file per lipid moleculetype: the phospholipid and the sterol.
    assert "posre_memb.itp" in result.posres
    assert "posre_memb2.itp" in result.posres

    sterol = Path(result.posres["posre_memb2.itp"]).read_text()
    rows = [
        line.split()
        for line in sterol.split("[ position_restraints ]")[1].strip().splitlines()
        if line.strip() and not line.lstrip().startswith(";")
    ]
    assert len(rows) == 1, "the single hydroxyl oxygen anchors the sterol"
    assert rows[0][0] == "2"  # O1, not the carbon and not the hydrogen
    assert rows[0][2:] == ["0", "0", "POSRES_FC_LIPID"]  # z only, as for phosphate


# --------------------------------------------------------------------------- #
# Fast (GROMACS-solvation) path: the solvent parameter blocks
#
# A dry pack has no water and no ions in its prmtop, so the converted topology
# carries neither the [ atomtypes ] entries nor the [ moleculetype ] blocks that
# gmx solvate and gmx genion are about to reference. These tests pin the
# contract between prism.membrane.solvent_params and its two consumers
# (builder._probe_solvent_parameters / _fast_solvate and solvate.MembraneSolvator).
# --------------------------------------------------------------------------- #

from prism.membrane import solvent_params as _solvent_params
from prism.membrane.solvate import WATER_MOLECULETYPE, read_topology

#: A dry topology in exactly the shape ParmEd emits for a packed bilayer.
_DRY_TOP = """[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             1            0.83333333

[ atomtypes ]
; name    at.num    mass    charge ptype  sigma      epsilon
N3             7  14.010000  0.00000000  A     0.32499985        0.71128
CX             6  12.010000  0.00000000  A     0.33996695      0.4577296

[ moleculetype ]
; Name            nrexcl
system1          3

[ atoms ]
    1         N3      1    ALA      N      1 -0.40000000  14.010000
    2         CX      1    ALA     CA      2  0.40000000  12.010000

[ system ]
; Name
Generic title

[ molecules ]
; Compound       #mols
system1              1
"""


def _fake_blocks(**overrides):
    """SolventBlocks equivalent to a tip3p/Na+/Cl- tleap probe."""
    defaults = dict(
        water_moleculetype="SOL",
        water_sites=3,
        cation_moleculetype="Na+",
        anion_moleculetype="Cl-",
        atomtypes=(
            ("Na+", "Na+           11  22.990000  0.00000000  A     0.24392807     0.36584603\n"),
            ("Cl-", "Cl-           17  35.450000  0.00000000  A      0.4477657     0.14891274\n"),
            ("OW", "OW             8  16.000000  0.00000000  A     0.31507524       0.635968\n"),
            ("HW", "HW             1   1.008000  0.00000000  A              0              0\n"),
        ),
        moleculetypes=(
            (
                "Na+",
                "[ moleculetype ]\nNa+          3\n\n[ atoms ]\n"
                "    1        Na+      1    Na+    Na+      1 1.00000000  22.990000\n",
            ),
            (
                "Cl-",
                "[ moleculetype ]\nCl-          3\n\n[ atoms ]\n"
                "    1        Cl-      1    Cl-    Cl-      1 -1.00000000  35.450000\n",
            ),
            (
                "SOL",
                "[ moleculetype ]\nSOL          3\n\n[ atoms ]\n"
                "    1         OW      1    SOL      O      1 -0.83400000  16.000000\n"
                "    2         HW      1    SOL     H1      2 0.41700000   1.008000\n"
                "    3         HW      1    SOL     H2      3 0.41700000   1.008000\n"
                "\n[ settles ]\n1     1   0.09572000   0.15136000\n"
                "\n[ exclusions ]\n1  2  3\n2  1  3\n3  1  2\n",
            ),
        ),
        comb_rule=2,
        water_model="tip3p",
        charges={"Na+": 1.0, "Cl-": -1.0, "SOL": 0.0},
        source="tleap leaprc.water.tip3p (Na+/Cl-)",
    )
    defaults.update(overrides)
    return _solvent_params.SolventBlocks(**defaults)


def test_solvent_blocks_merge_is_visible_to_the_topology_parser(tmp_path):
    """The blocks must be inlined, not #included.

    solvate.read_topology skips preprocessor lines, so an #included
    moleculetype would be invisible to exactly the guards that exist to prove
    gmx genion has something to resolve -pname against.
    """
    top = tmp_path / "system.top"
    top.write_text(_DRY_TOP)

    _solvent_params.merge_solvent_blocks(str(top), _fake_blocks())

    view = read_topology(str(top))
    assert WATER_MOLECULETYPE in view.moleculetypes
    assert view.moleculetypes[WATER_MOLECULETYPE].natoms == 3
    assert {"Na+", "Cl-"} <= set(view.moleculetypes)
    assert view.moleculetypes["Na+"].charge == pytest.approx(1.0)
    assert view.moleculetypes["Cl-"].charge == pytest.approx(-1.0)
    # gmx solvate and gmx genion own [ molecules ]; a pre-seeded entry would be
    # counted twice.
    assert view.molecules == (("system1", 1),)


def test_solvent_blocks_merge_is_idempotent(tmp_path):
    """A retry must not duplicate a moleculetype; read_topology rejects those."""
    top = tmp_path / "system.top"
    top.write_text(_DRY_TOP)

    _solvent_params.merge_solvent_blocks(str(top), _fake_blocks())
    once = top.read_text()
    _solvent_params.merge_solvent_blocks(str(top), _fake_blocks())

    assert top.read_text() == once
    read_topology(str(top))  # raises on a duplicate [ moleculetype ]


def test_solvent_blocks_refuse_a_clashing_atomtype(tmp_path):
    """Reusing a same-named type with different LJ parameters is silent and wrong."""
    top = tmp_path / "system.top"
    top.write_text(
        _DRY_TOP.replace(
            "CX             6  12.010000  0.00000000  A     0.33996695      0.4577296",
            "CX             6  12.010000  0.00000000  A     0.33996695      0.4577296\n"
            "OW             8  16.000000  0.00000000  A     0.99999999       0.111111",
        )
    )

    with pytest.raises(_solvent_params.SolventParameterError, match="different"):
        _solvent_params.merge_solvent_blocks(str(top), _fake_blocks())


def test_solvent_blocks_refuse_a_mismatched_combination_rule(tmp_path):
    """sigma/epsilon and C6/C12 are not interchangeable."""
    top = tmp_path / "system.top"
    top.write_text(_DRY_TOP.replace("1               2               yes", "1               1               yes"))

    with pytest.raises(_solvent_params.SolventParameterError, match="combination rule"):
        _solvent_params.merge_solvent_blocks(str(top), _fake_blocks())


@pytest.mark.parametrize(
    "overrides, match",
    [
        ({"water_moleculetype": "WAT"}, "SOL"),
        ({"charges": {"Na+": -1.0, "Cl-": -1.0, "SOL": 0.0}}, "wrong sign"),
        ({"charges": {"Na+": 1.0, "Cl-": -1.0, "SOL": 0.6}}, "neutral"),
        ({"comb_rule": 1}, "combination rule"),
    ],
)
def test_assert_solvent_blocks_rejects_unusable_probes(overrides, match):
    """A mis-mapped ion unit or a WAT-named water only fails much later, if ever."""
    with pytest.raises(_solvent_params.SolventParameterError, match=match):
        _solvent_params.assert_solvent_blocks(_fake_blocks(**overrides))


def test_builder_fast_path_modules_are_importable():
    """builder._fast_path_modules silently disables the fast path on ImportError.

    A missing solvent_params module therefore does not fail a build -- it makes
    every build take the slow PACKMOL route, which is the regression this pins.
    """
    solvate_module, solvent_params_module = MembraneBuilder._fast_path_modules()
    assert solvate_module.WATER_MOLECULETYPE == "SOL"
    for required in ("solvent_parameter_blocks", "assert_solvent_blocks", "merge_solvent_blocks"):
        assert hasattr(solvent_params_module, required), required


def test_membrane_builder_reports_no_fast_path_blockers_when_tools_exist(monkeypatch):
    """The blocker list must not name the fast-path modules once they exist."""
    builder = MembraneBuilder(MembraneConfig(enabled=True, solvate="gromacs"), verbose=False)
    monkeypatch.setattr(MembraneBuilder, "_gmx_command", staticmethod(lambda: "gmx"))
    blockers = builder._fast_path_blockers()
    assert not any("solvation modules are unavailable" in b for b in blockers), blockers


# --------------------------------------------------------------------------- #
# The default bilayer route: prism.membrane.patch
#
# bilayer_source defaults to "auto", and on a POPC system "auto" resolves to
# this module: tile a bundled pre-equilibrated Lipid21 patch, drop the oriented
# solute in unmoved, delete the lipids that overlap it. Nothing downstream
# re-derives any of that. A swapped asset, a mis-tiled lattice or an off-by-one
# in the clash deletion all yield a system that still passes grompp and still
# runs -- the damage only surfaces as wrong physics, long after the build. The
# tests below therefore target the failures that are silent, not the ones that
# crash.
# --------------------------------------------------------------------------- #

import hashlib
import os
import shutil
import subprocess

import numpy as np

from prism.membrane import patch as _patch
from prism.membrane import runscript as _runscript

#: Resolved from the package rather than through ``patch_data_dir()`` so that a
#: PRISM_MEMBRANE_PATCH_DIR left set in the developer's environment cannot
#: quietly redirect collection to somebody else's assets.
_BUNDLED_PATCH_DIR = (
    Path(_patch.__file__).resolve().parent.parent / "data" / "membrane_patches"
)
_BUNDLED_MANIFEST = json.loads(
    (_BUNDLED_PATCH_DIR / "manifest.json").read_text(encoding="utf-8")
)
_BUNDLED_PATCH_KEYS = sorted(_BUNDLED_MANIFEST["patches"])

#: GRO stores nm to three decimals, i.e. 0.01 A. Any geometry the manifest
#: recorded from full-precision arrays can therefore come back from the file up
#: to half a quantum different per coordinate. 0.05 A is comfortably above that
#: and still ~50x smaller than the gap between any two bundled patches.
_GRO_QUANTISATION_A = 0.05

#: A patch carrying more than one lipid species, so the tests that care about
#: the species ordering contract have an ordering available to get wrong.
_MIXED_PATCH = "POPC_CHOL15"


@pytest.fixture
def bundled_patches(monkeypatch):
    """Force the bundled assets, ignoring any patch dir set in the environment."""
    monkeypatch.delenv(_patch.PATCH_DIR_ENV, raising=False)
    return _BUNDLED_PATCH_DIR


def _write_solute_pdb(path, coords, resname="LIG", element="C"):
    """Write an all-carbon solute PDB at the given coordinates.

    Column layout matches ``_write_z_stack_pdb`` above; the serial wraps modulo
    100000 so a large synthetic solute cannot push the coordinate fields out of
    their columns.
    """
    lines = []
    for serial, (x, y, z) in enumerate(coords, start=1):
        lines.append(
            f"ATOM  {serial % 100000:5d} {'C':^4s} {resname:>3s} A{1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        )
    path.write_text("".join(lines) + "TER\nEND\n", encoding="ascii")
    return path


def _molecules_within(assembly, points, cutoff):
    """Molecule indices owning a heavy atom within ``cutoff`` of any point.

    A deliberately naive minimum-image brute force, independent of the cell
    list / KD-tree the production search picks at runtime. Affordable because
    the tests that use it place only a handful of solute atoms.
    """
    lx, ly = float(assembly.box[0]), float(assembly.box[1])
    heavy = assembly.heavy_mask
    xyz = assembly.xyz[heavy]
    molecule = assembly.atom_molecule[heavy]
    hit = np.zeros(xyz.shape[0], dtype=bool)
    for point in np.atleast_2d(np.asarray(points, dtype=float)):
        delta = xyz - point
        delta[:, 0] -= lx * np.round(delta[:, 0] / lx)
        delta[:, 1] -= ly * np.round(delta[:, 1] / ly)
        hit |= (delta * delta).sum(axis=1) < cutoff * cutoff
    return set(int(m) for m in np.unique(molecule[hit]))


def _closest_intermolecular_gap(assembly, plane_x, halfwidth=3.5):
    """Closest heavy-atom approach between *different* molecules near a plane.

    Restricted to a slab around ``plane_x`` so the O(n^2) sweep stays cheap,
    and periodic in x/y so a contact across the box face counts. Returns
    (min_distance, n_pairs_below_the_default_clash_cutoff).
    """
    lx, ly = float(assembly.box[0]), float(assembly.box[1])
    heavy = assembly.heavy_mask
    xyz = assembly.xyz[heavy]
    molecule = assembly.atom_molecule[heavy]

    offset = xyz[:, 0] - plane_x
    offset -= lx * np.round(offset / lx)
    near = np.abs(offset) < halfwidth
    xyz, molecule = xyz[near], molecule[near]

    delta = xyz[:, None, :] - xyz[None, :, :]
    delta[:, :, 0] -= lx * np.round(delta[:, :, 0] / lx)
    delta[:, :, 1] -= ly * np.round(delta[:, :, 1] / ly)
    distance = np.sqrt((delta * delta).sum(axis=2))
    distance[molecule[:, None] == molecule[None, :]] = np.inf
    return float(distance.min()), int((distance < _patch.DEFAULT_CLASH_CUTOFF).sum() // 2)


# -- the bundled assets ----------------------------------------------------- #


@pytest.mark.parametrize("key", _BUNDLED_PATCH_KEYS)
def test_every_bundled_patch_matches_the_geometry_its_manifest_declares(
    key, bundled_patches
):
    """The manifest is the only description of an asset anyone ever reads.

    ``required_tiling`` sizes the box from the manifest's box, the builder sizes
    the water slab from its thickness, and the topology's ``[molecules]`` counts
    come from its composition. If the file and the manifest ever disagree -- a
    regenerated asset, a half-applied edit, two patches swapped -- every one of
    those is computed from the wrong bilayer and nothing checks back. So the
    declared numbers are re-derived here from the file itself.
    """
    entry = _BUNDLED_MANIFEST["patches"][key]
    path = bundled_patches / entry["filename"]
    assert path.is_file(), f"{key} is in the manifest but its file is missing"

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == entry["sha256"], f"{key}: bundled file does not match its checksum"

    assembly = _patch.load_patch(key)
    assert assembly.n_atoms == entry["n_atoms"]
    assert assembly.n_molecules == entry["n_molecules"]
    assert assembly.species_counts() == {n: c for n, c in entry["composition"]}
    assert assembly.box.tolist() == pytest.approx(entry["box"], abs=1e-4)

    geometry = _patch.patch_geometry(assembly)
    measured = entry["measured"]
    assert geometry.n_lower == measured["n_lower"]
    assert geometry.n_upper == measured["n_upper"]
    assert geometry.area_per_lipid == pytest.approx(
        measured["area_per_lipid"], abs=_GRO_QUANTISATION_A
    )
    assert geometry.phosphate_thickness == pytest.approx(
        measured["phosphate_thickness"], abs=_GRO_QUANTISATION_A
    )
    assert geometry.core_half_width == pytest.approx(
        measured["core_half_width"], abs=_GRO_QUANTISATION_A
    )

    # A patch is only usable if both leaflets carry the same number of lipids;
    # an imbalance builds a permanent spontaneous curvature into the system.
    assert geometry.n_lower == geometry.n_upper, "leaflets are not balanced"

    # Every alias has to reach the same asset, or a user asking for "PC" gets a
    # different bilayer from the one asking for "POPC".
    for alias in entry.get("aliases", ()):
        assert _patch.load_patch(alias, verify=False).key == key


def test_a_tampered_patch_is_refused_rather_than_tiled(tmp_path, monkeypatch):
    """Corruption must stop the build, because its symptoms are not local.

    A patch is replicated across the whole box, so a damaged one does not fail
    in one place -- it produces a bilayer that is uniformly wrong and entirely
    plausible-looking. The checksum is the first guard and the structural
    checks are the second, and the second matters on its own: an asset swapped
    *together with* its checksum (a bad regeneration, a mis-merged manifest)
    passes the first one.
    """
    monkeypatch.setenv(_patch.PATCH_DIR_ENV, str(tmp_path))
    shutil.copy(_BUNDLED_PATCH_DIR / "manifest.json", tmp_path / "manifest.json")
    shutil.copy(_BUNDLED_PATCH_DIR / "POPC.gro.gz", tmp_path / "POPC.gro.gz")
    assert _patch.load_patch("POPC").n_atoms == 17152  # the copy is good

    # 1. A flipped byte: caught by the checksum.
    payload = bytearray((tmp_path / "POPC.gro.gz").read_bytes())
    payload[-1] ^= 0xFF
    (tmp_path / "POPC.gro.gz").write_bytes(bytes(payload))
    with pytest.raises(_patch.PatchIntegrityError, match="[Cc]hecksum"):
        _patch.load_patch("POPC")

    # 2. A different lipid published under POPC's name *and* POPC's checksum
    #    updated to match -- the checksum now agrees and the asset is still the
    #    wrong bilayer. Only the structural cross-check can see this.
    shutil.copy(_BUNDLED_PATCH_DIR / "POPE.gro.gz", tmp_path / "POPC.gro.gz")
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    manifest["patches"]["POPC"]["sha256"] = hashlib.sha256(
        (tmp_path / "POPC.gro.gz").read_bytes()
    ).hexdigest()
    (tmp_path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(_patch.PatchIntegrityError):
        _patch.load_patch("POPC")


def test_an_unknown_lipid_names_the_patches_that_do_exist(bundled_patches):
    """The failure has to be actionable: the alternative is the 30-minute route."""
    with pytest.raises(_patch.PatchNotAvailableError) as excinfo:
        _patch.load_patch("DOPC")
    message = str(excinfo.value)
    assert "POPC" in message and "PACKMOL" in message


# -- tiling ----------------------------------------------------------------- #


def test_tiling_scales_the_cell_without_disturbing_the_equilibrated_structure(
    bundled_patches,
):
    """Integer replication has to be exactly that: same material, bigger box.

    Everything downstream trusts this. ``[molecules]`` counts are the tiled
    counts, the water slab is sized from the tiled box, and the whole reason
    for tiling rather than packing is that the lipids stay equilibrated. A
    tiler that stretched, rotated or dropped molecules would still emit a
    loadable bilayer at the wrong density, which no later step re-measures.
    """
    patch = _patch.load_patch("POPC")
    nx, ny = 2, 3
    base = _patch.tile_patch(patch, nx=1, ny=1, center_xy=(0.0, 0.0))
    tiled = _patch.tile_patch(patch, nx=nx, ny=ny, center_xy=(0.0, 0.0))
    replicas = nx * ny

    # The box grows in xy only; z is the bilayer normal and must not change.
    assert tiled.box[0] == pytest.approx(patch.box[0] * nx)
    assert tiled.box[1] == pytest.approx(patch.box[1] * ny)
    assert tiled.box[2] == pytest.approx(patch.box[2])

    assert tiled.n_molecules == patch.n_molecules * replicas
    assert tiled.n_atoms == patch.n_atoms * replicas
    assert tiled.n_tiles == (nx, ny)

    # Per species *and* per leaflet, so a tiler that replicated one leaflet
    # twice and the other not at all cannot pass on totals alone.
    for species, (lower, upper) in patch.leaflet_counts().items():
        assert tiled.leaflet_counts()[species] == (lower * replicas, upper * replicas)
    assert sum(v for v, _ in tiled.leaflet_counts().values()) == sum(
        u for _, u in tiled.leaflet_counts().values()
    )

    # Area per lipid is the density check: identical before and after.
    assert _patch.patch_geometry(tiled).area_per_lipid == pytest.approx(
        _patch.patch_geometry(patch).area_per_lipid, abs=1e-6
    )

    # Every molecule's internal geometry survives untouched. Comparing the
    # sorted distribution of atom-minus-own-centroid vectors needs no molecule
    # correspondence, and any rotation, rescaling or distortion moves it.
    #
    # The reference is the freshly loaded patch, never a 1x1 tiling: a
    # distortion introduced inside tile_patch would be applied to a tiled
    # reference as well, and the two would then agree while both were wrong.
    def internal(assembly):
        centroids = assembly.molecule_centroids()
        counts = np.diff(assembly.mol_offsets)
        return assembly.xyz - np.repeat(centroids, counts, axis=0)

    assert np.allclose(
        np.sort(internal(tiled), axis=0),
        np.sort(np.tile(internal(patch), (replicas, 1)), axis=0),
        atol=1e-9,
    )
    # A 1x1 tiling is the identity on everything except placement.
    assert base.n_molecules == patch.n_molecules
    assert base.box.tolist() == pytest.approx(patch.box.tolist())


def test_tiling_leaves_the_seams_indistinguishable_from_the_bulk(bundled_patches):
    """The seam is where a wrong lattice vector would show, and only there.

    Replicating molecules that were never wrapped into the cell, or trimming a
    tiled patch back to an exact target box, both leave the interior perfect
    and wreck the joins -- hard overlaps that minimisation cannot undo, or a
    low-density crack the lipids never close. Both survive grompp. So the seam
    plane is measured against an interior plane of the same cell: the two must
    look the same, and neither may contain a clash.
    """
    patch = _patch.load_patch("POPC")
    lx = float(patch.box[0])
    tiled = _patch.tile_patch(patch, nx=2, ny=2, center_xy=(0.0, 0.0))
    origin_x = 0.0 - 0.5 * float(tiled.box[0])

    seam_gap, seam_clashes = _closest_intermolecular_gap(tiled, origin_x + lx)
    bulk_gap, bulk_clashes = _closest_intermolecular_gap(tiled, origin_x + 0.5 * lx)

    assert seam_clashes == 0, f"{seam_clashes} lipid-lipid clashes at the 2x2 seam"
    assert bulk_clashes == 0, "the untiled interior already contains clashes"
    assert seam_gap >= _patch.DEFAULT_CLASH_CUTOFF
    # Self-normalising: no absolute number to drift, just "the join is as good
    # as the material it joins".
    assert seam_gap > 0.9 * bulk_gap, (
        f"seam closest approach {seam_gap:.3f} A is worse than the bulk "
        f"{bulk_gap:.3f} A -- the tiles do not line up"
    )

    # The same patch handed in one cell over has to tile to the same bilayer.
    #
    # The bundled assets already sit inside the primary cell, so this cannot be
    # provoked with them -- but PRISM_MEMBRANE_PATCH_DIR exists precisely so a
    # user can supply their own, and a deposited frame is typically not wrapped
    # at all (atoms several box lengths from the origin). Replicating one of
    # those verbatim stacks tiles on top of each other and opens voids
    # elsewhere: a bilayer at locally double density that still loads, still
    # grompps, and is not recoverable by minimisation.
    counts = np.diff(patch.mol_offsets)
    offsets = np.zeros((patch.n_molecules, 3))
    offsets[::2, 0] = float(patch.box[0])          # half of them a cell to the right
    offsets[1::3, 1] = -2.0 * float(patch.box[1])  # and some two cells down
    displaced = _patch.LipidAssembly(
        atom_names=patch.atom_names,
        res_names=patch.res_names,
        res_local=patch.res_local,
        elements=patch.elements,
        xyz=patch.xyz + np.repeat(offsets, counts, axis=0),
        mol_offsets=patch.mol_offsets,
        mol_species=patch.mol_species,
        mol_leaflet=patch.mol_leaflet,
        box=patch.box.copy(),
        key=patch.key,
        core_half_width=patch.core_half_width,
    )
    from_displaced = _patch.tile_patch(displaced, nx=2, ny=2, center_xy=(0.0, 0.0))
    assert np.allclose(from_displaced.xyz, tiled.xyz, atol=1e-6), (
        "a patch translated by whole cell vectors tiled to a different bilayer"
    )


def test_required_tiling_never_returns_a_box_smaller_than_asked_for(bundled_patches):
    """Under-tiling by one is invisible: the build succeeds with the protein
    too close to its own periodic image."""
    patch = _patch.load_patch("POPC")
    lx, ly = float(patch.box[0]), float(patch.box[1])
    for target_x, target_y in [(1.0, 1.0), (lx, ly), (lx + 0.01, ly + 0.01), (250.0, 130.0)]:
        nx, ny = _patch.required_tiling(patch, target_x, target_y)
        assert nx * lx >= target_x - 1e-6 and ny * ly >= target_y - 1e-6
        # and minimal: one tile fewer would not have covered it
        assert (nx - 1) * lx < target_x and (ny - 1) * ly < target_y
    with pytest.raises(ValueError):
        _patch.required_tiling(patch, 0.0, 10.0)


# -- insertion and clash deletion ------------------------------------------- #


def test_insert_solute_refuses_a_solute_wider_than_the_tiled_bilayer(
    tmp_path, bundled_patches
):
    """This guard was dead once, and a dead guard here is not a crash.

    The clearance was computed with lo/hi transposed and against a box origin
    taken from ``lipids.xyz.min()``, so it returned roughly a box length and
    never fired. What gets through is a protein sticking out of its own
    bilayer: it solvates, it grompps, it runs, and the protein interacts with
    its periodic image through water for the whole trajectory. The two cases
    below are the ones the transposed form could not tell apart.
    """
    patch = _patch.load_patch("POPC")
    tiled = _patch.tile_patch(patch, nx=1, ny=1, center_xy=(0.0, 0.0))
    lx = float(tiled.box[0])

    too_wide = _patch.read_solute(
        _write_solute_pdb(tmp_path / "wide.pdb", [(-lx, 0.0, 0.0), (lx, 0.0, 0.0)])
    )
    with pytest.raises(_patch.MembranePatchError, match="does not fit|sticks out"):
        _patch.insert_solute(tiled, too_wide)

    # Comfortably inside: must still be accepted, or the guard is merely noisy.
    fits = _patch.read_solute(
        _write_solute_pdb(
            tmp_path / "fits.pdb", [(-0.25 * lx, 0.0, 0.0), (0.25 * lx, 0.0, 0.0)]
        )
    )
    assert _patch.insert_solute(tiled, fits).n_atoms == tiled.n_atoms + 2

    # The margin is the protein-to-image separation, so it has to bite before
    # the solute actually touches the face.
    with pytest.raises(_patch.MembranePatchError, match="does not fit"):
        _patch.insert_solute(tiled, fits, margin=0.3 * lx)

    # The check must depend on the box and the solute only. Offsetting the
    # slab's centre is exactly what the old ``lipids.xyz.min()`` origin was
    # sensitive to, and it must make no difference here.
    offset = _patch.tile_patch(patch, nx=1, ny=1, center_xy=(500.0, -500.0))
    with pytest.raises(_patch.MembranePatchError):
        _patch.insert_solute(offset, too_wide)


def test_clash_deletion_removes_exactly_the_lipids_that_overlap(
    tmp_path, bundled_patches
):
    """Deleting one lipid too few leaves an overlap no minimiser can undo;
    one too many leaves a vacuum bubble that closes over the first ns and
    changes the local area per lipid. Neither is visible in any log, so the
    surviving set is compared atom-for-atom against an independent search.
    """
    patch = _patch.load_patch("POPC")
    tiled = _patch.tile_patch(patch, nx=1, ny=1, center_xy=(0.0, 0.0))
    probes = [(0.0, 0.0, 0.0), (10.0, -8.0, 6.0)]
    solute = _patch.read_solute(_write_solute_pdb(tmp_path / "probe.pdb", probes))

    system = _patch.insert_solute(tiled, solute)
    cutoff = 6.0
    result = _patch.delete_clashing_lipids(system, cutoff)

    doomed = _molecules_within(tiled, probes, cutoff)
    assert doomed, "the probe must overlap something for this test to mean anything"
    keep = np.ones(tiled.n_molecules, dtype=bool)
    keep[sorted(doomed)] = False
    expected = tiled.select_molecules(keep)

    assert result.deletion.n_deleted == len(doomed)
    assert result.lipids.n_molecules == expected.n_molecules
    assert np.array_equal(result.lipids.mol_species, expected.mol_species)
    assert np.array_equal(result.lipids.mol_leaflet, expected.mol_leaflet)
    assert np.array_equal(result.lipids.xyz, expected.xyz)
    # The report is what the builder logs and what the topology is written
    # from, so it has to agree with the coordinates it describes.
    assert result.deletion.surviving_by_species == expected.species_counts()
    assert result.deletion.surviving_by_leaflet == (
        int(np.count_nonzero(expected.mol_leaflet < 0)),
        int(np.count_nonzero(expected.mol_leaflet > 0)),
    )
    assert result.deletion.residual_contacts == 0
    assert result.deletion.min_heavy_distance >= cutoff


def test_clash_deletion_of_a_solute_that_overlaps_nothing_is_a_no_op(
    tmp_path, bundled_patches
):
    """A soluble domain far above the bilayer deletes nothing, and that is a
    result rather than an error. Raising here (or silently deleting a lipid
    anyway) would break every peripheral-protein system."""
    patch = _patch.load_patch("POPC")
    tiled = _patch.tile_patch(patch, nx=1, ny=1, center_xy=(0.0, 0.0))
    solute = _patch.read_solute(
        _write_solute_pdb(tmp_path / "far.pdb", [(0.0, 0.0, 400.0)])
    )

    result = _patch.delete_clashing_lipids(_patch.insert_solute(tiled, solute))

    assert result.deletion.n_deleted == 0
    assert result.deletion.deleted_by_species == {}
    assert result.deletion.deleted_by_leaflet == (0, 0)
    assert result.lipids.n_molecules == tiled.n_molecules
    assert np.array_equal(result.lipids.xyz, tiled.xyz)
    assert result.deletion.residual_contacts == 0


def test_clash_deletion_can_empty_a_leaflet_and_says_so(tmp_path, bundled_patches):
    """The pathological input still has to produce an honest report.

    A solute lying flat across one leaflet -- a badly oriented protein, an OPM
    fetch that returned the wrong frame -- wipes that leaflet. The build must
    not pretend otherwise, because a monolayer is the one outcome a caller has
    to be able to detect and reject.
    """
    patch = _patch.load_patch("POPC")
    tiled = _patch.tile_patch(patch, nx=1, ny=1, center_xy=(0.0, 0.0))
    half_x = 0.5 * float(tiled.box[0]) - 1.0
    half_y = 0.5 * float(tiled.box[1]) - 1.0
    grid = [
        (x, y, z)
        for x in np.arange(-half_x, half_x, 2.0)
        for y in np.arange(-half_y, half_y, 2.0)
        for z in (8.0, 14.0, 20.0)
    ]
    solute = _patch.read_solute(_write_solute_pdb(tmp_path / "sheet.pdb", grid))

    result = _patch.delete_clashing_lipids(_patch.insert_solute(tiled, solute))

    lower_before, upper_before = patch.leaflet_counts()["POPC"]
    assert result.deletion.deleted_by_leaflet[1] == upper_before
    assert result.deletion.surviving_by_leaflet[1] == 0, "upper leaflet must be gone"
    assert result.deletion.surviving_by_leaflet[0] > 0, "lower leaflet must survive"
    assert result.deletion.n_deleted == (
        lower_before + upper_before - result.lipids.n_molecules
    )
    assert result.deletion.residual_contacts == 0


# -- the ordering contract and determinism ---------------------------------- #


def test_emitted_molecule_order_matches_the_declared_composition(
    tmp_path, bundled_patches
):
    """``MembraneComposition.species`` is what a ``[molecules]`` section copies.

    GROMACS reads ``[molecules]`` as a run-length encoding of the coordinate
    file: it never checks the names against the coordinates. If the emitted
    order and the declared order drift apart, every atom past the drift is
    assigned the wrong parameters, and the system still builds, still grompps
    and still runs. It is the worst failure in this module and the quietest.
    """
    solute = _patch.read_solute(
        _write_solute_pdb(tmp_path / "core.pdb", [(0.0, 0.0, 0.0)])
    )
    system = _patch.build_dry_membrane_system(solute, lipid=_MIXED_PATCH, xy_margin=4.0)
    composition = system.composition()
    lipids = system.lipids

    assert composition.solute_atoms == solute.n_atoms
    assert composition.n_atoms == solute.n_atoms + lipids.n_atoms

    # 1. Species appear as contiguous runs, in the declared order and count.
    runs = []
    for name in lipids.mol_species.tolist():
        if not runs or runs[-1][0] != name:
            runs.append([name, 0])
        runs[-1][1] += 1
    assert tuple((n, c) for n, c in runs) == composition.species

    # 2. Within each species, the lower leaflet precedes the upper.
    start = 0
    for name, count in composition.species:
        leaflets = lipids.mol_leaflet[start:start + count]
        boundary = int(np.count_nonzero(leaflets < 0))
        assert (leaflets[:boundary] < 0).all() and (leaflets[boundary:] > 0).all(), (
            f"{name}: leaflets are interleaved, not lower-then-upper"
        )
        assert composition.per_leaflet[name] == (boundary, count - boundary)
        start += count

    # 3. Molecules are contiguous, which is what GROMACS actually requires.
    assert np.all(np.diff(lipids.mol_offsets) > 0)
    assert lipids.mol_offsets[-1] == lipids.n_atoms

    # 4. The charge the caller has to neutralise is derived from those counts.
    assert composition.formal_charge == sum(
        _patch.LIPID21_FORMAL_CHARGE[name] * count for name, count in composition.species
    )

    # 5. The written coordinates carry the same species runs in the same order.
    written = system.write_gro(tmp_path / "dry.gro")
    body = written.read_text(encoding="utf-8").splitlines()[2:-1][solute.n_atoms:]
    residues = [line[5:10].strip() for line in body]
    seen = []
    for name in residues:
        if not seen or seen[-1] != name:
            seen.append(name)
    for name, _ in composition.species:
        assert set(_patch.LIPID21_RESIDUE_PATTERNS[name]) <= set(seen)


def test_the_solute_is_emitted_exactly_where_it_was_read(tmp_path, bundled_patches):
    """The bilayer is positioned around the solute, never the other way round.

    The solute arrives already oriented in the membrane frame -- from OPM, or
    from PRISM's own alignment -- and that orientation is the only thing making
    the system a membrane protein rather than a protein in a lipid box. A
    translation applied here would be undetectable downstream: the system still
    looks like a protein in a bilayer, just not the one that was validated.
    """
    # Chosen to be exact in %8.3f, so "unmoved" is asserted as equality rather
    # than as a tolerance a small systematic drift could hide inside.
    coords = [(3.5, -7.25, 1.125), (-12.0, 4.5, -9.75)]
    source = _write_solute_pdb(tmp_path / "solute.pdb", coords)
    solute = _patch.read_solute(source)
    assert solute.xyz.tolist() == [list(c) for c in coords]

    system = _patch.build_dry_membrane_system(solute, lipid="POPC", xy_margin=10.0)
    assert system.solute.xyz.tolist() == [list(c) for c in coords]

    written = system.write_pdb(tmp_path / "dry.pdb").read_text(encoding="utf-8")
    emitted = [
        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        for line in written.splitlines()
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "LIG"
    ]
    assert emitted == coords


def test_the_dry_build_is_byte_identical_across_processes(tmp_path, bundled_patches):
    """No optimiser, no randomness -- so the same inputs owe the same bytes.

    The entire case for this route over PACKMOL is that it is a deterministic
    geometric construction rather than a stochastic search. If that quietly
    stopped being true, the route would keep working and merely stop being
    reproducible, which nothing would ever report. Separate interpreters with
    different hash seeds also catch the usual way determinism is lost here: a
    set or dict iteration order leaking into the molecule ordering. The
    two-species patch is deliberate: with a single lipid there is no species
    ordering for a hash seed to permute, so a pure patch could not detect it.
    """
    source = _write_solute_pdb(
        tmp_path / "solute.pdb", [(0.0, 0.0, 0.0), (6.0, 4.0, -3.0)]
    )

    def digests(target):
        system = _patch.build_dry_membrane_system(
            source, lipid=_MIXED_PATCH, xy_margin=8.0
        )
        return [
            hashlib.sha256(system.write_gro(target / "dry.gro").read_bytes()).hexdigest(),
            hashlib.sha256(system.write_pdb(target / "dry.pdb").read_bytes()).hexdigest(),
        ]

    first = digests(tmp_path / "a")
    assert digests(tmp_path / "b") == first, "two builds in one process disagree"

    # patch.py imports nothing from prism itself, so it loads standalone in a
    # fraction of the time the whole package takes.
    driver = tmp_path / "driver.py"
    driver.write_text(
        "import hashlib, importlib.util, sys\n"
        f"spec = importlib.util.spec_from_file_location('p', {str(Path(_patch.__file__))!r})\n"
        "m = importlib.util.module_from_spec(spec)\n"
        "sys.modules['p'] = m\n"
        "spec.loader.exec_module(m)\n"
        f"s = m.build_dry_membrane_system({str(source)!r}, lipid={_MIXED_PATCH!r}, "
        "xy_margin=8.0)\n"
        "for name, write in (('dry.gro', s.write_gro), ('dry.pdb', s.write_pdb)):\n"
        "    print(hashlib.sha256(write(sys.argv[1] + '/' + name).read_bytes()).hexdigest())\n",
        encoding="ascii",
    )
    for seed in ("0", "1234567"):
        out = tmp_path / f"seed{seed}"
        out.mkdir()
        environ = dict(os.environ, PYTHONHASHSEED=seed, PYTHONNOUSERSITE="1")
        environ.pop(_patch.PATCH_DIR_ENV, None)
        completed = subprocess.run(
            [sys.executable, str(driver), str(out)],
            capture_output=True, text=True, env=environ, check=True,
        )
        assert completed.stdout.split() == first, (
            f"PYTHONHASHSEED={seed} produced different bytes"
        )


# --------------------------------------------------------------------------- #
# The membrane driver script: prism.membrane.runscript
#
# This module exists for one reason. The membrane MDPs thermostat on
# ``tc-grps = SOLU MEMB SOLV``, and those three are PRISM-defined index groups
# rather than GROMACS built-ins, so every single grompp call needs
# ``-n index.ndx`` or it aborts. The soluble path's generator passes none. A
# regression here does not corrupt a system -- it produces a run script that
# cannot start at all, or worse, one that starts and thermostats the wrong
# atoms.
# --------------------------------------------------------------------------- #

#: A generated script is only meaningful if a shell can parse it.
_HAVE_BASH = shutil.which("bash") is not None


def _grompp_calls(script_text):
    """Every ``gmx grompp`` command line in a generated script."""
    return [
        line.strip()
        for line in script_text.splitlines()
        if "grompp" in line and '"$GMX"' in line
    ]


def _command_options(line):
    """``-flag value`` pairs of a shell command line, ignoring bare words."""
    tokens = line.split()
    options = {}
    for index, token in enumerate(tokens):
        if token.startswith("-") and index + 1 < len(tokens):
            following = tokens[index + 1]
            if not following.startswith("-"):
                options[token] = following
    return options


def _write_stage_mdps(directory, stages, define=None):
    """Minimal .mdp files for stage-plumbing tests that do not need real physics."""
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    for stage in stages:
        path = directory / f"{stage}.mdp"
        text = "integrator = md\nnsteps = 10\n"
        if define and stage in define:
            text += f"define = {define[stage]}\n"
        path.write_text(text, encoding="ascii")
        paths[stage] = str(path)
    return paths


def test_every_grompp_call_passes_the_index_file(tmp_path):
    """The single reason this module exists, over the soluble generator.

    ``tc-grps = SOLU MEMB SOLV`` names PRISM's own index groups, so a grompp
    without ``-n`` stops with "Group SOLU referenced in the .mdp file was not
    found in the index file" and the overnight job never starts. One missed
    call is enough, and it is the *last* stage that is most likely to be
    missed and least likely to be noticed before production.
    """
    mdps = write_membrane_mdps(str(tmp_path))
    system_dir = tmp_path / MEMBRANE_SYSTEM_DIRNAME
    script = Path(_runscript.write_membrane_runscript(str(system_dir), mdps))
    text = script.read_text(encoding="utf-8")

    calls = _grompp_calls(text)
    assert len(calls) == len(mdps), "every stage must grompp exactly once"
    for call in calls:
        options = _command_options(call)
        assert options.get("-n") == _runscript.INDEX_NAME, call
        assert options.get("-p") == _runscript.TOPOL_NAME, call

    # The names are parameters, not literals: the builder is free to rename the
    # index file, and a hardcoded "index.ndx" would then point at nothing.
    custom = Path(
        _runscript.write_membrane_runscript(
            str(tmp_path / "custom"), mdps,
            gro="built.gro", top="membrane.top", ndx="groups.ndx",
        )
    ).read_text(encoding="utf-8")
    for call in _grompp_calls(custom):
        options = _command_options(call)
        assert options.get("-n") == "groups.ndx", call
        assert options.get("-p") == "membrane.top", call
    assert _command_options(_grompp_calls(custom)[0])["-c"] == "built.gro"


def test_every_grompp_call_passes_a_restraint_reference(tmp_path):
    """GROMACS >= 2018 refuses to build a .tpr when restraints are live and no
    ``-r`` was given, and the staged protocol keeps ``-DPOSRES`` on through
    every equilibration stage. ``-r`` also has to name the *previous* stage's
    output: referencing the original build instead would drag every atom back
    to its packed position at each stage, undoing the equilibration it just
    paid for while reporting nothing.
    """
    mdps = write_membrane_mdps(str(tmp_path))
    system_dir = tmp_path / MEMBRANE_SYSTEM_DIRNAME
    text = Path(
        _runscript.write_membrane_runscript(str(system_dir), mdps)
    ).read_text(encoding="utf-8")

    calls = _grompp_calls(text)
    for call in calls:
        options = _command_options(call)
        assert "-r" in options, f"no restraint reference: {call}"
        assert options["-r"] == options["-c"], (
            f"-r must restrain to the coordinates the stage starts from: {call}"
        )

    # The restrained stages specifically -- the ones grompp would reject.
    restrained = [
        stage for stage in mdps
        if "POSRES" in (_mdp_option(Path(mdps[stage]).read_text(), "define") or "").upper()
    ]
    assert set(restrained) == set(_RESTRAINED_STAGES), restrained
    for stage in restrained:
        matching = [c for c in calls if f"-f " in c and f"{stage}.mdp" in c]
        assert matching, f"no grompp call for restrained stage {stage}"
        assert all("-r " in c for c in matching)

    # The chain: stage N starts from stage N-1's output, and only the first
    # stage reads the freshly built system.
    assert _command_options(calls[0])["-c"] == _runscript.SYSTEM_GRO_NAME
    for previous, call in zip(calls, calls[1:]):
        produced = _command_options(previous)["-o"]
        assert _command_options(call)["-c"] == produced.replace(".tpr", ".gro"), call


@pytest.mark.parametrize(
    "declared, expected",
    [
        # Declaration order must not survive into execution order.
        (["md", "npt2", "em", "npt", "nvt"], ["em", "nvt", "npt", "npt2", "md"]),
        (["md", "npt2", "npt", "nvt", "em"], ["em", "nvt", "npt", "npt2", "md"]),
        # Nor may an alphabetic one: this is the npt10 trap. Sorted lexically,
        # "npt10" lands between "npt" and "npt2" and a tenth NPT stage would
        # silently run second.
        (["npt10", "npt9", "npt2", "npt"], ["npt", "npt2", "npt9", "npt10"]),
        (["npt", "npt2", "npt10"], ["npt", "npt2", "npt10"]),
        # Phases, not a hardcoded stage list: the protocol is being reworked.
        (["production", "minimization", "heat", "equil", "nvt"],
         ["minimization", "heat", "nvt", "equil", "production"]),
        # Case is a spelling, not a phase.
        (["MD", "NPT2", "EM", "NVT"], ["EM", "NVT", "NPT2", "MD"]),
        # An unrecognised stage is almost certainly extra equilibration, so it
        # goes after what we do recognise but never displaces either anchor.
        (["md", "weird", "em", "nvt"], ["em", "nvt", "weird", "md"]),
    ],
)
def test_stage_order_follows_physical_phase_not_declaration_or_alphabet(
    declared, expected
):
    """Order is the protocol. Running npt before nvt, or production before the
    restraints have been released, produces a trajectory that looks entirely
    normal and started from an unequilibrated bilayer."""
    assert _runscript.order_membrane_stages(declared) == expected


def test_minimisation_stays_first_and_production_last_for_any_input_order(tmp_path):
    """The two anchors are structural, not conventional.

    Minimisation is the only stage that can consume the freshly packed
    coordinates -- it is what relieves the insertion clashes -- and production
    is the unrestrained run everything else exists to set up. Either one out of
    place still runs to completion.
    """
    stages = ["md", "npt2", "em", "npt", "nvt"]
    for rotation in range(len(stages)):
        declared = stages[rotation:] + stages[:rotation]
        ordered = _runscript.order_membrane_stages(declared)
        assert ordered[0] == "em"
        assert ordered[-1] == "md"

    # And the emitted script agrees with the ordering function.
    mdps = _write_stage_mdps(tmp_path / "mdps", ["md", "npt2", "em", "npt", "nvt"])
    text = Path(
        _runscript.write_membrane_runscript(str(tmp_path / "sys"), mdps)
    ).read_text(encoding="utf-8")
    used = [_command_options(c)["-f"].rsplit("/", 1)[-1] for c in _grompp_calls(text)]
    assert used == ["em.mdp", "nvt.mdp", "npt.mdp", "npt2.mdp", "md.mdp"]


def test_production_is_written_where_every_downstream_tool_looks(tmp_path):
    """prism.analysis, prism.sim and prism.mmpbsa all discover a run by looking
    for ``<system dir>/prod/md.xtc``. Production written anywhere else is not a
    broken run -- it is a finished run that is invisible."""
    mdps = _write_stage_mdps(tmp_path / "mdps", ["em", "nvt", "npt", "npt2", "md"])
    text = Path(
        _runscript.write_membrane_runscript(str(tmp_path / "sys"), mdps)
    ).read_text(encoding="utf-8")

    final = _grompp_calls(text)[-1]
    assert _command_options(final)["-o"] == "./prod/md.tpr"
    assert "-deffnm ./prod/md " in text

    # Equilibration stays in its own per-stage directory, so a rerun of one
    # stage cannot overwrite another's outputs.
    for stage in ("em", "nvt", "npt", "npt2"):
        assert f"-deffnm ./{stage}/{stage} " in text


@pytest.mark.skipif(not _HAVE_BASH, reason="needs a bash to parse the generated script")
def test_the_generated_script_is_valid_shell(tmp_path):
    """A syntax error is a whole queue slot wasted before anything runs, and
    the script is assembled from string fragments, quoted echo lines and a
    trap, so it is exactly the sort of thing that breaks silently on an edit."""
    mdps = write_membrane_mdps(str(tmp_path))
    system_dir = tmp_path / MEMBRANE_SYSTEM_DIRNAME
    script = Path(_runscript.write_membrane_runscript(str(system_dir), mdps))

    for candidate in (script, script.parent / _runscript.RUNSCRIPT_BACKUP_NAME):
        completed = subprocess.run(
            ["bash", "-n", str(candidate)], capture_output=True, text=True
        )
        assert completed.returncode == 0, f"{candidate.name}: {completed.stderr}"

    # The backup exists because prism.sim overwrites localrun.sh with the
    # soluble script, which has no -n and cannot grompp this system. It has to
    # be a real copy, and the sentinel is how a generator recognises a membrane
    # script it did not author.
    backup = script.parent / _runscript.RUNSCRIPT_BACKUP_NAME
    assert backup.read_text(encoding="utf-8") == script.read_text(encoding="utf-8")
    assert script.read_text(encoding="utf-8").splitlines()[1] == _runscript.RUNSCRIPT_SENTINEL
    assert os.access(script, os.X_OK)


@pytest.mark.parametrize(
    "hostile",
    [
        "; rm -rf /",       # command separator
        "$(id)",            # command substitution
        "`id`",             # legacy command substitution
        "../../etc",        # path traversal: stage names become directories
        "em nvt",           # word splitting
        "-rf",              # reads as an option wherever it is interpolated
        "",                 # empty
        "a;b",
        "*",                # glob
    ],
)
def test_hostile_stage_names_are_refused_rather_than_quoted(hostile, tmp_path):
    """Stage names become directory names, file stems and bare shell words --
    including the argument to an ``rm -rf``. They are rejected outright instead
    of being defensively quoted at each of the dozen sites that interpolate
    them, because one missed site is a deleted working directory."""
    mdps = _write_stage_mdps(tmp_path / "mdps", ["em"])
    with pytest.raises(ValueError, match="stage name"):
        _runscript.write_membrane_runscript(
            str(tmp_path / "sys"), {hostile: mdps["em"]}
        )
    assert not (tmp_path / "sys" / _runscript.RUNSCRIPT_NAME).exists()


def test_an_empty_protocol_is_refused(tmp_path):
    """Writing a script with no stages would produce a driver that exits 0
    having run nothing -- indistinguishable from a completed build."""
    with pytest.raises(ValueError, match="empty"):
        _runscript.write_membrane_runscript(str(tmp_path / "sys"), {})


def test_the_numpy_fallback_finds_the_same_clashes_as_the_kd_tree(
    tmp_path, bundled_patches, monkeypatch
):
    """Two neighbour searches, one answer -- and nothing checks that.

    ``setup.py`` requires numpy and lists scipy only as an extra, so the
    dependency-free cell list is not a rare fallback: it is what a base PRISM
    install actually runs. The two implementations decide which lipids get
    deleted, so a divergence between them means two users build measurably
    different bilayers from identical inputs, with nothing in either log to
    say so. Cross-checking them against each other is the only place that can
    be caught.
    """
    pytest.importorskip("scipy.spatial", reason="nothing to cross-check without scipy")

    patch = _patch.load_patch("POPC")
    tiled = _patch.tile_patch(patch, nx=1, ny=1, center_xy=(0.0, 0.0))
    # Deliberately near a box face, so the periodic wrap-around that the two
    # searches implement differently is actually exercised.
    probes = [(0.0, 0.0, 0.0), (28.0, -29.0, 5.0), (-30.0, 30.0, -4.0)]
    solute = _patch.read_solute(_write_solute_pdb(tmp_path / "probe.pdb", probes))
    system = _patch.insert_solute(tiled, solute)
    cutoff = 7.0

    with_scipy = _patch.delete_clashing_lipids(system, cutoff)

    # A None entry in sys.modules makes "from scipy.spatial import cKDTree"
    # raise ImportError, which is exactly the state of a base install.
    monkeypatch.setitem(sys.modules, "scipy.spatial", None)
    without_scipy = _patch.delete_clashing_lipids(system, cutoff)

    assert without_scipy.deletion.n_deleted == with_scipy.deletion.n_deleted
    assert without_scipy.deletion.deleted_by_species == with_scipy.deletion.deleted_by_species
    assert without_scipy.deletion.deleted_by_leaflet == with_scipy.deletion.deleted_by_leaflet
    assert np.array_equal(without_scipy.lipids.xyz, with_scipy.lipids.xyz)
    assert without_scipy.deletion.min_heavy_distance == pytest.approx(
        with_scipy.deletion.min_heavy_distance, abs=1e-9
    )
    # Both must actually have deleted something, or the agreement is vacuous.
    assert with_scipy.deletion.n_deleted > 0
