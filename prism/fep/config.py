#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Configuration file parser for PRISM-FEP

Supports FEbuilder-compatible config.conf files.
"""

import configparser
from pathlib import Path
from typing import Dict, Any


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
