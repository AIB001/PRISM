#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Configuration file parser for PRISM-FEP

Supports both:
- Legacy FEbuilder-compatible config.conf files
- New YAML-based config.yaml + fep.yaml architecture
"""

import configparser
from pathlib import Path
from typing import Dict, Any
import yaml


class FEPConfig:
    """
    Unified FEP configuration manager

    Loads configuration from:
    - config.yaml: General PRISM parameters (force field, temperature, etc.)
    - fep.yaml: FEP-specific parameters (mapping cutoffs, charge strategies)
    """

    def __init__(self, work_dir: str):
        """
        Initialize FEP configuration from work directory

        Parameters
        ----------
        work_dir : str
            Working directory containing config.yaml and fep.yaml
        """
        self.work_dir = Path(work_dir)
        self.config_yaml = self.work_dir / "config.yaml"
        self.fep_yaml = self.work_dir / "fep.yaml"

        # Load configurations
        self.general_config = self._load_yaml(self.config_yaml)
        self.fep_config = self._load_yaml(self.fep_yaml)

    def _load_yaml(self, yaml_file: Path) -> Dict[str, Any]:
        """Load YAML file"""
        if yaml_file.exists():
            with open(yaml_file) as f:
                return yaml.safe_load(f) or {}
        return {}

    def get_forcefield_type(self) -> str:
        """Get force field type from config.yaml"""
        return self.general_config.get("forcefield", {}).get("type", "gaff")

    def get_forcefield_params(self) -> Dict[str, Any]:
        """Get force field parameters from config.yaml"""
        return self.general_config.get("forcefield", {}).get("params", {})

    def get_mapping_params(self) -> Dict[str, Any]:
        """
        Get mapping parameters from fep.yaml

        Returns dict with keys: dist_cutoff, charge_cutoff, charge_common, charge_reception
        """
        defaults = {"dist_cutoff": 0.6, "charge_cutoff": 0.05, "charge_common": "mean", "charge_reception": "pert"}
        mapping_config = self.fep_config.get("mapping", {})
        # Merge with defaults
        return {**defaults, **mapping_config}

    def get_html_config(self) -> Dict[str, Any]:
        """Get configuration for HTML visualization"""
        return {
            "forcefield": {"type": self.get_forcefield_type(), "params": self.get_forcefield_params()},
            "fep": {"mapping": self.get_mapping_params()},
        }

    def __repr__(self):
        return f"FEPConfig(work_dir={self.work_dir})"


# Legacy functions for backward compatibility


def read_fep_config(config_file: str) -> Dict[str, Any]:
    """
    Read FEP configuration file

    Supports FEbuilder-compatible config.conf format:

    [Model]
    charge_common = ref|mut|mean
    charge_reception = pert|unique|surround
    distance = 12
    recharge_hydrogen = False

    [Other]
    dist_cutoff = 0.6
    charge_cutoff = 0.05

    Parameters
    ----------
    config_file : str
        Path to configuration file

    Returns
    -------
    dict
        Configuration parameters

    Examples
    --------
    >>> config = read_fep_config('tests/gxf/FEP/unit_test/25-36/config.conf')
    >>> mapper = DistanceAtomMapper(
    ...     dist_cutoff=config['dist_cutoff'],
    ...     charge_cutoff=config['charge_cutoff'],
    ...     charge_common=config['charge_common'],
    ...     charge_reception=config['charge_reception']
    ... )
    """
    if not Path(config_file).exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")

    config = configparser.ConfigParser()
    config.read(config_file)

    # Default values
    params = {
        "dist_cutoff": 0.6,
        "charge_cutoff": 0.05,
        "charge_common": "mean",
        "charge_reception": "pert",
        "recharge_hydrogen": False,
        "distance": 10,  # Solvation distance
    }

    # Read [Model] section
    if "Model" in config:
        model = config["Model"]
        if "charge_common" in model:
            params["charge_common"] = model["charge_common"]
        if "charge_reception" in model:
            params["charge_reception"] = model["charge_reception"]
        if "distance" in model:
            params["distance"] = float(model["distance"])

    # Read [Other] section
    if "Other" in config:
        other = config["Other"]
        if "dist_cutoff" in other:
            params["dist_cutoff"] = float(other["dist_cutoff"])
        if "charge_cutoff" in other:
            params["charge_cutoff"] = float(other["charge_cutoff"])
        if "recharge_hydrogen" in other:
            params["recharge_hydrogen"] = other["recharge_hydrogen"].lower() in ["true", "yes", "1"]

    return params


def write_fep_config(config_file: str, params: Dict[str, Any]):
    """
    Write FEP configuration file

    Parameters
    ----------
    config_file : str
        Path to output configuration file
    params : dict
        Configuration parameters
    """
    config = configparser.ConfigParser()

    # [Model] section
    config.add_section("Model")
    if "charge_common" in params:
        config.set("Model", "charge_common", params["charge_common"])
    if "charge_reception" in params:
        config.set("Model", "charge_reception", params["charge_reception"])
    if "distance" in params:
        config.set("Model", "distance", str(params["distance"]))

    # [Other] section
    config.add_section("Other")
    if "dist_cutoff" in params:
        config.set("Other", "dist_cutoff", str(params["dist_cutoff"]))
    if "charge_cutoff" in params:
        config.set("Other", "charge_cutoff", str(params["charge_cutoff"]))
    if "recharge_hydrogen" in params:
        config.set("Other", "recharge_hydrogen", str(params["recharge_hydrogen"]))

    with open(config_file, "w") as f:
        config.write(f)
