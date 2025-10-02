#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Configuration - Enhanced configuration loading, validation, and management

Provides comprehensive configuration management with validation, templates,
and migration capabilities for PRISM molecular dynamics simulations.
"""

import os
import yaml
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List, Union
from .config_validator import ConfigurationValidator, ValidationLevel, validate_config_file
from .config_templates import ConfigTemplateGenerator, ConfigTemplate


logger = logging.getLogger(__name__)


class ConfigurationManager:
    """Enhanced configuration loading, validation, and management"""

    def __init__(self, config_path=None, gromacs_env=None,
                 forcefield_name="amber99sb", water_model_name="tip3p",
                 validate=True, validation_level=ValidationLevel.STANDARD):
        """
        Initialize configuration manager with validation
        
        Parameters:
        -----------
        config_path : str, optional
            Path to configuration file
        gromacs_env : GromacsEnvironment, optional
            GROMACS environment for force field detection
        forcefield_name : str
            Default force field name
        water_model_name : str
            Default water model name
        validate : bool
            Whether to validate configuration
        validation_level : ValidationLevel
            Level of validation to perform
        """
        self.config_path = config_path
        self.gromacs_env = gromacs_env
        self.forcefield_name = forcefield_name
        self.water_model_name = water_model_name
        self.validate_config = validate
        self.validation_level = validation_level
        self.validator = ConfigurationValidator() if validate else None
        self.template_generator = ConfigTemplateGenerator()
        
        # Load and validate configuration
        self.config = self._load_config()
        self._process_forcefield_names()
        
        # Validate if requested
        if self.validate_config:
            self._validate_configuration()

    def _load_config(self):
        """Load configuration from YAML file or use defaults"""
        default_config = self._get_default_config()

        if self.config_path and os.path.exists(self.config_path):
            print(f"Loading configuration from: {self.config_path}")
            with open(self.config_path, 'r') as f:
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
                7: {"name": "gromos54a7", "dir": "gromos54a7.ff"}
            }
            gmx_command = 'gmx'

        return {
            'general': {
                'overwrite': False,
                'gmx_command': gmx_command
            },
            'forcefield': {
                'name': self.forcefield_name,  # Use name instead of index
                'index': 1,  # Will be determined from name
                'custom_forcefields': force_fields
            },
            'water_model': {
                'name': self.water_model_name,  # Use name instead of index
                'index': 1,  # Will be determined from name
                'custom_water_models': {
                    1: {"name": "tip3p"},
                    2: {"name": "tip4p"},
                    3: {"name": "spc"},
                    4: {"name": "spce"},
                    5: {"name": "none"}
                }
            },
            'box': {
                'distance': 1.5,  # nm
                'shape': 'cubic',  # cubic, dodecahedron, octahedron
                'center': True
            },
            'simulation': {
                'temperature': 310,  # K
                'pressure': 1.0,  # bar
                'pH': 7.0,
                'ligand_charge': 0,  # Default ligand charge
                'production_time_ns': 500,  # ns
                'dt': 0.002,  # ps
                'equilibration_nvt_time_ps': 500,  # ps
                'equilibration_npt_time_ps': 500  # ps
            },
            'ions': {
                'neutral': True,
                'concentration': 0.15,  # M
                'positive_ion': 'NA',
                'negative_ion': 'CL'
            },
            'constraints': {
                'algorithm': 'lincs',
                'type': 'h-bonds',  # none, h-bonds, all-bonds
                'lincs_iter': 1,
                'lincs_order': 4
            },
            'energy_minimization': {
                'integrator': 'steep',
                'emtol': 200.0,  # kJ/mol/nm
                'emstep': 0.01,
                'nsteps': 10000
            },
            'output': {
                'trajectory_interval_ps': 500,  # ps
                'energy_interval_ps': 10,  # ps
                'log_interval_ps': 10,  # ps
                'compressed_trajectory': True
            },
            'electrostatics': {
                'coulombtype': 'PME',
                'rcoulomb': 1.0,  # nm
                'pme_order': 4,
                'fourierspacing': 0.16  # nm
            },
            'vdw': {
                'rvdw': 1.0,  # nm
                'dispcorr': 'EnerPres'
            },
            'temperature_coupling': {
                'tcoupl': 'V-rescale',
                'tc_grps': ['Protein', 'Non-Protein'],
                'tau_t': [0.1, 0.1],  # ps
            },
            'pressure_coupling': {
                'pcoupl': 'C-rescale',
                'pcoupltype': 'isotropic',
                'tau_p': 1.0,  # ps
                'compressibility': 4.5e-5  # bar^-1
            }
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
        ff_name = self.config['forcefield'].get('name')
        if ff_name and isinstance(ff_name, str):
            ff_index = self.gromacs_env.get_force_field_index(ff_name)
            if ff_index:
                self.config['forcefield']['index'] = ff_index
                print(f"Force field '{ff_name}' mapped to index {ff_index}")
            else:
                print(f"Warning: Force field '{ff_name}' not found. Available force fields:")
                for ff in self.gromacs_env.list_force_fields():
                    print(f"  - {ff}")
                raise ValueError(f"Force field '{ff_name}' not found in GROMACS installation")

        # Update water models for the selected force field
        self.update_water_models(self.config['forcefield']['index'])

        # Handle water model name
        water_name = self.config['water_model'].get('name')
        if water_name and isinstance(water_name, str):
            ff_index = self.config['forcefield']['index']
            wm_index = self.gromacs_env.get_water_model_index(ff_index, water_name)
            if wm_index:
                self.config['water_model']['index'] = wm_index
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
                self.config['water_model']['custom_water_models'] = water_models

    def save_config(self, output_path):
        """Save configuration to file"""
        with open(output_path, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False, sort_keys=False)

    def update_forcefield_by_name(self, ff_name):
        """Update force field by name"""
        if not self.gromacs_env:
            print("Warning: Cannot update force field without GROMACS environment")
            return

        ff_index = self.gromacs_env.get_force_field_index(ff_name)
        if ff_index:
            self.config['forcefield']['name'] = ff_name
            self.config['forcefield']['index'] = ff_index
            self.update_water_models(ff_index)
            print(f"Updated force field to '{ff_name}' (index {ff_index})")
        else:
            raise ValueError(f"Force field '{ff_name}' not found")

    def update_water_model_by_name(self, water_name):
        """Update water model by name"""
        if not self.gromacs_env:
            print("Warning: Cannot update water model without GROMACS environment")
            return

        ff_index = self.config['forcefield']['index']
        wm_index = self.gromacs_env.get_water_model_index(ff_index, water_name)
        if wm_index:
            self.config['water_model']['name'] = water_name
            self.config['water_model']['index'] = wm_index
            print(f"Updated water model to '{water_name}' (index {wm_index})")
        else:
            raise ValueError(f"Water model '{water_name}' not found for current force field")
    
    def _validate_configuration(self):
        """Validate current configuration"""
        if not self.validator:
            return
        
        logger.debug(f"Validating configuration with {self.validation_level.value} level")
        is_valid, report = self.validator.validate_config(self.config, self.validation_level)
        
        if not is_valid:
            error_summary = f"Configuration validation failed with {len(report['errors'])} errors"
            logger.error(error_summary)
            
            # Log first few errors
            for i, error in enumerate(report['errors'][:3], 1):
                logger.error(f"  {i}. {error.get('message', 'Unknown error')}")
            
            if len(report['errors']) > 3:
                logger.error(f"  ... and {len(report['errors']) - 3} more errors")
            
            # Optionally raise exception for invalid config
            if self.validation_level == ValidationLevel.STRICT:
                raise ValueError(f"Configuration validation failed: {error_summary}")
        
        # Log warnings
        if report['warnings']:
            logger.warning(f"Configuration validation completed with {len(report['warnings'])} warnings")
            for warning in report['warnings'][:2]:  # Show first 2 warnings
                logger.warning(f"  - {warning.get('message', 'Unknown warning')}")
        
        logger.info("Configuration validation completed successfully")
    
    def get_validation_report(self) -> str:
        """Get detailed validation report"""
        if not self.validator:
            return "Validation is disabled for this configuration manager"
        
        return self.validator.generate_validation_report(self.config, self.validation_level)
    
    def validate_and_fix(self) -> Dict[str, Any]:
        """Validate configuration and suggest fixes"""
        if not self.validator:
            return {"status": "validation_disabled"}
        
        is_valid, report = self.validator.validate_config(self.config, self.validation_level)
        
        fixes_applied = []
        warnings_noted = []
        
        # Apply automatic fixes for common issues
        for error in report.get('errors', []):
            error_code = error.get('error_code')
            field_path = error.get('field_path', '')
            
            if error_code == 'TEMPERATURE_INCONSISTENCY':
                # Fix temperature inconsistencies
                target_temp = 310.0  # Default body temperature
                self._fix_temperature_consistency(target_temp)
                fixes_applied.append(f"Set all temperatures to {target_temp} K")
            
            elif error_code == 'TIMESTEP_INCONSISTENCY':
                # Fix timestep inconsistencies
                target_dt = 0.002
                self._fix_timestep_consistency(target_dt)
                fixes_applied.append(f"Set all timesteps to {target_dt} ps")
        
        # Note warnings for user attention
        for warning in report.get('warnings', []):
            warnings_noted.append(warning.get('message', 'Unknown warning'))
        
        return {
            "status": "completed",
            "original_valid": is_valid,
            "fixes_applied": fixes_applied,
            "warnings": warnings_noted,
            "validation_report": self.get_validation_report()
        }
    
    def _fix_temperature_consistency(self, target_temp: float):
        """Fix temperature consistency across sections"""
        sections_to_fix = ['simulation', 'smd', 'umbrella', 'analysis']
        
        for section in sections_to_fix:
            if section in self.config and isinstance(self.config[section], dict):
                if 'temperature' in self.config[section]:
                    self.config[section]['temperature'] = target_temp
        
        # Fix equilibration subsections
        if 'equilibration' in self.config:
            for subsection in ['nvt', 'npt']:
                if (subsection in self.config['equilibration'] and 
                    isinstance(self.config['equilibration'][subsection], dict)):
                    self.config['equilibration'][subsection]['temperature'] = target_temp
    
    def _fix_timestep_consistency(self, target_dt: float):
        """Fix timestep consistency across sections"""
        sections_to_fix = ['smd', 'umbrella']
        
        for section in sections_to_fix:
            if section in self.config and isinstance(self.config[section], dict):
                if 'dt' in self.config[section]:
                    self.config[section]['dt'] = target_dt
        
        # Fix equilibration subsections
        if 'equilibration' in self.config:
            for subsection in ['nvt', 'npt']:
                if (subsection in self.config['equilibration'] and 
                    isinstance(self.config['equilibration'][subsection], dict)):
                    self.config['equilibration'][subsection]['dt'] = target_dt
    
    def create_template(self, template_type: str, custom_params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Create configuration from template
        
        Parameters:
        -----------
        template_type : str
            Template type ('fast', 'standard', 'accurate', 'production')
        custom_params : dict, optional
            Custom parameters to override
            
        Returns:
        --------
        Dict : Generated configuration
        """
        template_map = {
            'fast': ConfigTemplate.PMF_FAST,
            'standard': ConfigTemplate.PMF_STANDARD,
            'accurate': ConfigTemplate.PMF_ACCURATE,
            'production': ConfigTemplate.PMF_PRODUCTION,
            'drug_discovery': ConfigTemplate.DRUG_DISCOVERY,
            'membrane': ConfigTemplate.MEMBRANE_PROTEIN
        }
        
        if template_type not in template_map:
            raise ValueError(f"Unknown template: {template_type}. Available: {list(template_map.keys())}")
        
        template_enum = template_map[template_type]
        return self.template_generator.generate_template(template_enum, custom_params)
    
    def save_template(self, template_type: str, output_path: str, 
                     custom_params: Optional[Dict[str, Any]] = None):
        """Save configuration template to file"""
        template_map = {
            'fast': ConfigTemplate.PMF_FAST,
            'standard': ConfigTemplate.PMF_STANDARD,
            'accurate': ConfigTemplate.PMF_ACCURATE,
            'production': ConfigTemplate.PMF_PRODUCTION,
            'drug_discovery': ConfigTemplate.DRUG_DISCOVERY,
            'membrane': ConfigTemplate.MEMBRANE_PROTEIN
        }
        
        if template_type not in template_map:
            raise ValueError(f"Unknown template: {template_type}")
        
        template_enum = template_map[template_type]
        self.template_generator.save_template(template_enum, output_path, custom_params)
        logger.info(f"Template '{template_type}' saved to {output_path}")
    
    def migrate_config(self, target_version: str = "latest") -> Dict[str, Any]:
        """
        Migrate configuration to newer version format
        
        Parameters:
        -----------
        target_version : str
            Target configuration version
            
        Returns:
        --------
        Dict : Migration report
        """
        logger.info(f"Migrating configuration to version {target_version}")
        
        migration_report = {
            "original_config": self.config.copy(),
            "migrations_applied": [],
            "warnings": [],
            "success": True
        }
        
        # Example migrations - add more as needed
        if target_version == "latest" or target_version >= "2.0":
            # Migrate old parameter names
            migrations = [
                ("smd.pullf_rate", "smd.pull_rate"),
                ("umbrella.force_const", "umbrella.force_constant"),
                ("analysis.energy_units", "analysis.energy_unit")
            ]
            
            for old_path, new_path in migrations:
                if self._migrate_parameter(old_path, new_path):
                    migration_report["migrations_applied"].append(f"{old_path} -> {new_path}")
            
            # Add missing required sections
            if "output" not in self.config:
                self.config["output"] = {
                    "log_level": "INFO",
                    "save_configs": True,
                    "validate_inputs": True
                }
                migration_report["migrations_applied"].append("Added 'output' section")
        
        # Validate migrated configuration
        if self.validator:
            is_valid, _ = self.validator.validate_config(self.config, ValidationLevel.BASIC)
            if not is_valid:
                migration_report["success"] = False
                migration_report["warnings"].append("Migrated configuration failed validation")
        
        logger.info(f"Migration completed. Applied {len(migration_report['migrations_applied'])} changes")
        return migration_report
    
    def _migrate_parameter(self, old_path: str, new_path: str) -> bool:
        """Migrate a single parameter from old path to new path"""
        old_parts = old_path.split('.')
        new_parts = new_path.split('.')
        
        # Check if old parameter exists
        current = self.config
        for part in old_parts[:-1]:
            if part not in current or not isinstance(current[part], dict):
                return False
            current = current[part]
        
        if old_parts[-1] not in current:
            return False
        
        old_value = current[old_parts[-1]]
        
        # Create new parameter location
        current = self.config
        for part in new_parts[:-1]:
            if part not in current:
                current[part] = {}
            current = current[part]
        
        # Move value and remove old
        current[new_parts[-1]] = old_value
        del self.config[old_parts[-2]][old_parts[-1]]
        
        return True
    
    def get_summary(self) -> str:
        """Get configuration summary"""
        lines = []
        lines.append("PRISM Configuration Summary")
        lines.append("=" * 40)
        
        if self.config_path:
            lines.append(f"Source: {self.config_path}")
        else:
            lines.append("Source: Default configuration")
        
        # Basic settings
        if 'forcefield' in self.config:
            ff_name = self.config['forcefield'].get('name', 'Unknown')
            lines.append(f"Force Field: {ff_name}")
        
        if 'water_model' in self.config:
            wm_name = self.config['water_model'].get('name', 'Unknown')
            lines.append(f"Water Model: {wm_name}")
        
        if 'simulation' in self.config:
            sim = self.config['simulation']
            temp = sim.get('temperature', 'Unknown')
            lines.append(f"Temperature: {temp} K")
        
        # PMF specific
        if 'smd' in self.config:
            smd = self.config['smd']
            pull_rate = smd.get('pull_rate', 'Unknown')
            lines.append(f"SMD Pull Rate: {pull_rate} nm/ps")
        
        if 'umbrella' in self.config:
            umb = self.config['umbrella']
            prod_time = umb.get('production_time_ps', 'Unknown')
            lines.append(f"Umbrella Production Time: {prod_time} ps")
        
        # Validation status
        if self.validate_config and self.validator:
            is_valid, report = self.validator.validate_config(self.config, ValidationLevel.BASIC)
            status = "✅ Valid" if is_valid else "❌ Invalid"
            lines.append(f"Validation Status: {status}")
        
        lines.append("=" * 40)
        return "\n".join(lines)