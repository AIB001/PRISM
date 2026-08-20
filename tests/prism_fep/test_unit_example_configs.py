"""Regression checks for publication-ready PRISM-FEP example inputs."""

from pathlib import Path

import yaml


EXAMPLE_ROOT = Path(__file__).parent / "unit_test"


def _load_yaml(path: Path) -> dict:
    data = yaml.safe_load(path.read_text())
    assert isinstance(data, dict), f"{path} must contain a YAML mapping"
    return data


def test_example_yaml_files_are_valid_mappings():
    configs = sorted(EXAMPLE_ROOT.glob("*/configs/*.yaml"))
    assert configs
    for config in configs:
        _load_yaml(config)


def test_example_default_profiles_and_input_references():
    case_42 = EXAMPLE_ROOT / "42-38"
    case_p38 = EXAMPLE_ROOT / "p38-19-24"
    case_ome = EXAMPLE_ROOT / "oMeEtPh-EtPh"

    for input_path in (
        case_42 / "input" / "receptor.pdb",
        case_42 / "input" / "42.mol2",
        case_42 / "input" / "38.mol2",
        case_p38 / "input" / "protein.pdb",
        case_p38 / "input" / "19.mol2",
        case_p38 / "input" / "24.mol2",
        case_ome / "input" / "receptor.pdb",
        case_ome / "input" / "oMeEtPh.mol2",
        case_ome / "input" / "EtPh.mol2",
    ):
        assert input_path.exists(), input_path

    assert _load_yaml(case_42 / "configs" / "config.yaml")["forcefield"]["type"] == "gaff2"
    assert _load_yaml(case_42 / "configs" / "fep.yaml")["forcefield"]["type"] == "gaff2"
    assert _load_yaml(case_42 / "configs" / "config_charmm.yaml")["forcefield"]["type"] == "cgenff"
    assert _load_yaml(case_42 / "configs" / "fep_charmm.yaml")["forcefield"]["type"] == "cgenff"
