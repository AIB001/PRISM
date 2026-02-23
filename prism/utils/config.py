#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Configuration - Configuration loading and management
"""

import os
import yaml


class ConfigurationManager:
    """Handle configuration loading and management"""

    def __init__(self, config_path=None, gromacs_env=None, forcefield_name="amber99sb", water_model_name="tip3p"):
        self.config_path = config_path
        self.gromacs_env = gromacs_env
        self.forcefield_name = forcefield_name
        self.water_model_name = water_model_name
        self.config = self._load_config()
        self._process_forcefield_names()

    def _load_config(self):
        """Load configuration from YAML file or use defaults"""
        default_config = self._get_default_config()

        if self.config_path and os.path.exists(self.config_path):
            print(f"Loading configuration from: {self.config_path}")
            with open(self.config_path, "r") as f:
                user_config = yaml.safe_load(f)

            # Merge user config with defaults
            config = self._merge_configs(default_config, user_config)
        else:
            if self.config_path:
                print(f"Config file not found: {self.config_path}. Using default configuration.")
            else:
                print("No config file specified. Using default configuration.")
            config = default_config

        return config

    def _get_default_config(self):
        """Get default configuration with detected force fields"""
        # Get force fields from GROMACS environment if available
        if self.gromacs_env:
            force_fields = self.gromacs_env.get_force_field_config()
            gmx_command = self.gromacs_env.gmx_command
        else:
            # Fallback to hardcoded defaults
            force_fields = {
                1: {"name": "amber99sb", "dir": "amber99sb.ff"},
                2: {"name": "amber99sb-ildn", "dir": "amber99sb-ildn.ff"},
                3: {"name": "amber03", "dir": "amber03.ff"},
                4: {"name": "amber14sb", "dir": "amber14sb.ff"},
                5: {"name": "charmm27", "dir": "charmm27.ff"},
                6: {"name": "oplsaa", "dir": "oplsaa.ff"},
                7: {"name": "gromos54a7", "dir": "gromos54a7.ff"},
            }
            gmx_command = "gmx"

        return {
            "general": {"overwrite": False, "gmx_command": gmx_command},
            "forcefield": {
                "name": self.forcefield_name,  # Use name instead of index
                "index": 1,  # Will be determined from name
                "custom_forcefields": force_fields,
            },
            "water_model": {
                "name": self.water_model_name,  # Use name instead of index
                "index": 1,  # Will be determined from name
                "custom_water_models": {
                    1: {"name": "tip3p"},
                    2: {"name": "tip4p"},
                    3: {"name": "spc"},
                    4: {"name": "spce"},
                    5: {"name": "none"},
                },
            },
            "box": {
                "distance": 1.5,  # nm
                "shape": "cubic",  # cubic, dodecahedron, octahedron
                "center": True,
            },
            "simulation": {
                "temperature": 310,  # K
                "pressure": 1.0,  # bar
                "pH": 7.0,
                "ligand_charge": 0,  # Default ligand charge
                "production_time_ns": 500,  # ns
                "dt": 0.002,  # ps
                "equilibration_nvt_time_ps": 500,  # ps
                "equilibration_npt_time_ps": 500,  # ps
            },
            "ions": {
                "neutral": True,
                "concentration": 0.15,  # M
                "positive_ion": "NA",
                "negative_ion": "CL",
            },
            "constraints": {
                "algorithm": "lincs",
                "type": "h-bonds",  # none, h-bonds, all-bonds
                "lincs_iter": 1,
                "lincs_order": 4,
            },
            "energy_minimization": {
                "integrator": "steep",
                "emtol": 200.0,  # kJ/mol/nm
                "emstep": 0.01,
                "nsteps": 10000,
            },
            "output": {
                "trajectory_interval_ps": 500,  # ps
                "energy_interval_ps": 10,  # ps
                "log_interval_ps": 10,  # ps
                "compressed_trajectory": True,
            },
            "electrostatics": {
                "coulombtype": "PME",
                "rcoulomb": 1.0,  # nm
                "pme_order": 4,
                "fourierspacing": 0.16,  # nm
            },
            "vdw": {
                "rvdw": 1.0,  # nm
                "dispcorr": "EnerPres",
            },
            "temperature_coupling": {
                "tcoupl": "V-rescale",
                "tc_grps": ["Protein", "Non-Protein"],
                "tau_t": [0.1, 0.1],  # ps
            },
            "pressure_coupling": {
                "pcoupl": "C-rescale",
                "pcoupltype": "isotropic",
                "tau_p": 1.0,  # ps
                "compressibility": 4.5e-5,  # bar^-1
            },
        }

    def _merge_configs(self, default, user):
        """Recursively merge user config with default config"""
        if isinstance(default, dict):
            if not isinstance(user, dict):
                return default
            merged = default.copy()
            for key, value in user.items():
                if key in merged:
                    merged[key] = self._merge_configs(merged[key], value)
                else:
                    merged[key] = value
            return merged
        else:
            return user

    def _process_forcefield_names(self):
        """Convert force field names to indices"""
        if not self.gromacs_env:
            return

        # Handle force field name
        ff_name = self.config["forcefield"].get("name")
        if ff_name and isinstance(ff_name, str):
            ff_index = self.gromacs_env.get_force_field_index(ff_name)
            if ff_index:
                self.config["forcefield"]["index"] = ff_index
                print(f"Force field '{ff_name}' mapped to index {ff_index}")
            else:
                print(f"Warning: Force field '{ff_name}' not found. Available force fields:")
                for ff in self.gromacs_env.list_force_fields():
                    print(f"  - {ff}")
                raise ValueError(f"Force field '{ff_name}' not found in GROMACS installation")

        # Update water models for the selected force field
        self.update_water_models(self.config["forcefield"]["index"])

        # Handle water model name
        water_name = self.config["water_model"].get("name")
        if water_name and isinstance(water_name, str):
            ff_index = self.config["forcefield"]["index"]
            wm_index = self.gromacs_env.get_water_model_index(ff_index, water_name)
            if wm_index:
                self.config["water_model"]["index"] = wm_index
                print(f"Water model '{water_name}' mapped to index {wm_index}")
            else:
                print(f"Warning: Water model '{water_name}' not found for force field. Available water models:")
                for wm in self.gromacs_env.list_water_models(ff_index):
                    print(f"  - {wm}")
                raise ValueError(f"Water model '{water_name}' not found for the selected force field")

    def update_water_models(self, ff_index):
        """Update water models based on selected force field"""
        if self.gromacs_env:
            water_models = self.gromacs_env.get_water_models_for_forcefield(ff_index)
            if water_models:
                self.config["water_model"]["custom_water_models"] = water_models

    def save_config(self, output_path):
        """Save configuration to file"""
        with open(output_path, "w") as f:
            yaml.dump(self.config, f, default_flow_style=False, sort_keys=False)

    def update_forcefield_by_name(self, ff_name):
        """Update force field by name"""
        if not self.gromacs_env:
            print("Warning: Cannot update force field without GROMACS environment")
            return

        ff_index = self.gromacs_env.get_force_field_index(ff_name)
        if ff_index:
            self.config["forcefield"]["name"] = ff_name
            self.config["forcefield"]["index"] = ff_index
            self.update_water_models(ff_index)
            print(f"Updated force field to '{ff_name}' (index {ff_index})")
        else:
            raise ValueError(f"Force field '{ff_name}' not found")

    def update_water_model_by_name(self, water_name):
        """Update water model by name"""
        if not self.gromacs_env:
            print("Warning: Cannot update water model without GROMACS environment")
            return

        ff_index = self.config["forcefield"]["index"]
        wm_index = self.gromacs_env.get_water_model_index(ff_index, water_name)
        if wm_index:
            self.config["water_model"]["name"] = water_name
            self.config["water_model"]["index"] = wm_index
            print(f"Updated water model to '{water_name}' (index {wm_index})")
        else:
            raise ValueError(f"Water model '{water_name}' not found for current force field")
