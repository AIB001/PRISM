#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Configuration file parser for PRISM-FEP

Supports both:
- Legacy FEbuilder-compatible config.conf files
- New YAML-based config.yaml + fep.yaml architecture
"""

import yaml
import configparser
from pathlib import Path
from typing import Dict, Any


DEFAULT_MAPPING_PARAMS = {
    "dist_cutoff": 0.6,
    "charge_cutoff": 0.05,
    "charge_common": "mean",
    "charge_reception": "surround",
}

DEFAULT_SIMULATION_PARAMS = {
    "equilibration_nvt_time_ps": 500,
    "equilibration_npt_time_ps": 500,
    "production_time_ns": 5.0,
    "dt": 0.002,
    "temperature": 310,
    "pressure": 1.0,
}

DEFAULT_SOFT_CORE_PARAMS = {"alpha": 0.5, "sigma": 0.3}

DEFAULT_LAMBDA_PARAMS = {
    "strategy": "decoupled",
    "distribution": "nonlinear",
    "windows": 32,
    "coul_windows": 12,
    "vdw_windows": 20,
}

DEFAULT_ELECTROSTATICS_PARAMS = {"rcoulomb": 1.0}
DEFAULT_VDW_PARAMS = {"rvdw": 1.0}
DEFAULT_OUTPUT_PARAMS = {
    "trajectory_interval_ps": 500,
    "energy_interval_ps": 10,
    "log_interval_ps": 10,
}
DEFAULT_EXECUTION_PARAMS = {
    "mode": "standard",
    "num_gpus": None,
    "parallel_windows": None,
    "total_cpus": None,
    "omp_threads": None,
    "use_gpu_pme": True,
}

DEFAULT_LEGACY_PARAMS = {
    "dist_cutoff": 0.6,
    "charge_cutoff": 0.05,
    "charge_common": "mean",
    "charge_reception": "surround",
    "recharge_hydrogen": False,
    "distance": 10,
}


def normalize_charge_reception(value: Any) -> str:
    """Normalize legacy charge_reception labels to the current API."""
    if value is None:
        return "surround"
    normalized = str(value).strip().lower()
    if normalized == "pert":
        return "surround"
    return normalized


class FEPConfig:
    """
    Unified FEP configuration manager

    Loads configuration from:
    - config.yaml: General PRISM parameters (force field, temperature, etc.)
    - fep.yaml: FEP-specific parameters (mapping cutoffs, charge strategies)
    """

    def __init__(self, work_dir: str, config_file: str = "config.yaml", fep_file: str = "fep.yaml"):
        """
        Initialize FEP configuration from work directory

        Parameters
        ----------
        work_dir : str
            Working directory containing the config YAML files.
        config_file : str, optional
            Filename of the general config YAML (default: ``config.yaml``).
        fep_file : str, optional
            Filename of the FEP-specific YAML (default: ``fep.yaml``).
        """
        self.work_dir = Path(work_dir)
        self.config_yaml = self.work_dir / config_file
        self.fep_yaml = self.work_dir / fep_file

        # Load configurations
        self.general_config = self._load_yaml(self.config_yaml)
        self.fep_config = self._load_yaml(self.fep_yaml)

    def _load_yaml(self, yaml_file: Path) -> Dict[str, Any]:
        """Load YAML file."""
        if yaml_file.exists():
            with open(yaml_file) as f:
                return yaml.safe_load(f) or {}
        return {}

    def _merged_fep_section(self, section: str, defaults: Dict[str, Any]) -> Dict[str, Any]:
        """Merge one FEP YAML section with its canonical defaults."""
        return {**defaults, **self.fep_config.get(section, {})}

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
        params = self._merged_fep_section("mapping", DEFAULT_MAPPING_PARAMS)
        params["charge_reception"] = normalize_charge_reception(params.get("charge_reception"))
        return params

    def get_html_config(self) -> Dict[str, Any]:
        """Get configuration for HTML visualization"""
        return {
            "forcefield": {"type": self.get_forcefield_type(), "params": self.get_forcefield_params()},
            "fep": {"mapping": self.get_mapping_params()},
        }

    def get_simulation_params(self) -> Dict[str, Any]:
        """
        Get simulation parameters from fep.yaml

        Returns dict with keys:
        - equilibration_nvt_time_ps: NVT equilibration time in ps
        - equilibration_npt_time_ps: NPT equilibration time in ps
        - production_time_ns: Production time per window in ns
        - dt: Time step in ps
        - temperature: Temperature in K
        - pressure: Pressure in bar
        """
        return self._merged_fep_section("simulation", DEFAULT_SIMULATION_PARAMS)

    def get_soft_core_params(self) -> Dict[str, float]:
        """
        Get soft-core parameters from fep.yaml

        Returns dict with keys:
        - alpha: Soft-core alpha parameter (0.3-0.5 recommended)
        - sigma: Soft-core sigma parameter in nm
        """
        return self._merged_fep_section("soft_core", DEFAULT_SOFT_CORE_PARAMS)

    def get_lambda_params(self) -> Dict[str, Any]:
        """
        Get lambda window configuration from fep.yaml

        Returns dict with keys:
        - strategy: Lambda strategy ('coupled', 'decoupled', or 'custom')
        - distribution: Distribution type ('linear', 'nonlinear', or 'quadratic')
        - windows: Total number of lambda windows (for coupled mode)
        - coul_windows: Number of windows for electrostatics stage (for decoupled mode)
        - vdw_windows: Number of windows for VDW stage (for decoupled mode)
        - custom_coul_lambdas: Optional custom coul lambda array
        - custom_vdw_lambdas: Optional custom vdw lambda array
        """
        lambda_config = self.fep_config.get("lambda", {})
        params = self._merged_fep_section("lambda", DEFAULT_LAMBDA_PARAMS)

        # Handle custom lambda arrays
        if "custom_coul_lambdas" in lambda_config:
            params["custom_coul_lambdas"] = lambda_config["custom_coul_lambdas"]
        if "custom_vdw_lambdas" in lambda_config:
            params["custom_vdw_lambdas"] = lambda_config["custom_vdw_lambdas"]

        return params

    def get_electrostatics_params(self) -> Dict[str, Any]:
        """
        Get electrostatics parameters from fep.yaml

        Returns dict with keys:
        - rcoulomb: Coulomb cutoff in nm
        """
        return self._merged_fep_section("electrostatics", DEFAULT_ELECTROSTATICS_PARAMS)

    def get_vdw_params(self) -> Dict[str, Any]:
        """
        Get van der Waals parameters from fep.yaml

        Returns dict with keys:
        - rvdw: VDW cutoff in nm
        """
        return self._merged_fep_section("vdw", DEFAULT_VDW_PARAMS)

    def get_output_params(self) -> Dict[str, Any]:
        """
        Get output control parameters from fep.yaml

        Returns dict with keys:
        - trajectory_interval_ps: Trajectory output interval in ps
        - energy_interval_ps: Energy output interval in ps
        - log_interval_ps: Log output interval in ps
        """
        return self._merged_fep_section("output", DEFAULT_OUTPUT_PARAMS)

    def get_execution_params(self) -> Dict[str, Any]:
        """
        Get execution-layer parameters from fep.yaml.

        These settings only affect generated run scripts. They do not change
        topology generation, MDP generation, or analysis.

        Returns dict with keys:
        - mode: Execution mode ('standard' or 'repex')
        - num_gpus: Total GPUs available to the generated scripts
        - parallel_windows: Number of concurrent windows in standard mode
        - omp_threads: OpenMP threads per worker/rank
        - use_gpu_pme: Whether to request GPU PME in mdrun commands
        """
        params = self._merged_fep_section("execution", DEFAULT_EXECUTION_PARAMS)
        params["mode"] = str(params.get("mode", "standard")).strip().lower()
        if params["mode"] not in {"standard", "repex"}:
            raise ValueError(f"Invalid execution.mode: {params['mode']}")
        return params

    def get_all_mdp_params(self) -> Dict[str, Any]:
        """
        Get all parameters needed for write_fep_mdps()

        Returns a dict matching the config parameter format expected by mdp_templates.py
        """
        return {
            "simulation": self.get_simulation_params(),
            "electrostatics": self.get_electrostatics_params(),
            "vdw": self.get_vdw_params(),
            "output": self.get_output_params(),
        }

    def to_dict(self) -> Dict[str, Any]:
        """Return the canonical merged FEP configuration view."""
        return {
            "forcefield": {"type": self.get_forcefield_type(), "params": self.get_forcefield_params()},
            "mapping": self.get_mapping_params(),
            "simulation": self.get_simulation_params(),
            "soft_core": self.get_soft_core_params(),
            "lambda": self.get_lambda_params(),
            "electrostatics": self.get_electrostatics_params(),
            "vdw": self.get_vdw_params(),
            "output": self.get_output_params(),
            "execution": self.get_execution_params(),
        }

    def __repr__(self):
        return f"FEPConfig(work_dir={self.work_dir})"


# Legacy functions for backward compatibility


def read_fep_config(config_file: str) -> Dict[str, Any]:
    """
    Read FEP configuration file

    Supports FEbuilder-compatible config.conf format:

    [Model]
    charge_common = ref|mut|mean|none
    charge_reception = unique|surround|surround_ext|none
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

    params = dict(DEFAULT_LEGACY_PARAMS)

    # Read [Model] section
    if "Model" in config:
        model = config["Model"]
        if "charge_common" in model:
            params["charge_common"] = model["charge_common"]
        if "charge_reception" in model:
            params["charge_reception"] = normalize_charge_reception(model["charge_reception"])
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

    params["charge_reception"] = normalize_charge_reception(params.get("charge_reception"))
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
