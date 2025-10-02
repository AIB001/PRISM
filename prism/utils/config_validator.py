#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Configuration Validator

Provides comprehensive validation for PRISM configuration files with schema validation,
type checking, value range validation, and dependency verification.
"""

import os
import re
import yaml
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class ValidationLevel(Enum):
    """Configuration validation levels"""
    BASIC = "basic"          # Basic type and required field checking
    STANDARD = "standard"    # Standard validation with range checks
    STRICT = "strict"        # Strict validation with full dependency checks


class ConfigValidationError(Exception):
    """Configuration validation error"""
    
    def __init__(self, message: str, field_path: str = None, 
                 suggestions: List[str] = None, error_code: str = None):
        self.message = message
        self.field_path = field_path
        self.suggestions = suggestions or []
        self.error_code = error_code
        super().__init__(message)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'message': self.message,
            'field_path': self.field_path,
            'suggestions': self.suggestions,
            'error_code': self.error_code
        }


class ConfigurationValidator:
    """Comprehensive configuration validator for PRISM"""
    
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.schema = self._load_validation_schema()
    
    def _load_validation_schema(self) -> Dict[str, Any]:
        """Load validation schema for PRISM configurations"""
        return {
            # General settings schema
            'general': {
                'type': 'object',
                'properties': {
                    'overwrite': {'type': 'boolean'},
                    'gmx_command': {'type': 'string', 'pattern': r'^gmx|gromacs'},
                    'project_name': {'type': 'string', 'minLength': 1}
                }
            },
            
            # Force field schema
            'forcefield': {
                'type': 'object',
                'required': ['name'],
                'properties': {
                    'name': {
                        'type': 'string', 
                        'enum': ['amber99sb', 'amber99sb-ildn', 'amber03', 'amber14sb', 
                                'charmm27', 'oplsaa', 'gromos54a7']
                    },
                    'index': {'type': 'integer', 'minimum': 1, 'maximum': 20}
                }
            },
            
            # Water model schema
            'water_model': {
                'type': 'object',
                'required': ['name'],
                'properties': {
                    'name': {
                        'type': 'string',
                        'enum': ['tip3p', 'tip4p', 'spc', 'spce', 'none']
                    },
                    'index': {'type': 'integer', 'minimum': 1, 'maximum': 10}
                }
            },
            
            # Simulation parameters schema
            'simulation': {
                'type': 'object',
                'properties': {
                    'temperature': {
                        'type': 'number', 
                        'minimum': 250, 
                        'maximum': 450,
                        'description': 'Temperature in Kelvin'
                    },
                    'pressure': {
                        'type': 'number',
                        'minimum': 0.5,
                        'maximum': 10.0,
                        'description': 'Pressure in bar'
                    },
                    'pH': {'type': 'number', 'minimum': 0, 'maximum': 14},
                    'production_time_ns': {
                        'type': 'number', 
                        'minimum': 1, 
                        'maximum': 10000
                    },
                    'dt': {
                        'type': 'number',
                        'minimum': 0.001,
                        'maximum': 0.005,
                        'description': 'Time step in ps'
                    }
                }
            },
            
            # PMF-specific schemas
            'smd': {
                'type': 'object',
                'properties': {
                    'nsteps': {'type': 'integer', 'minimum': 100000, 'maximum': 50000000},
                    'dt': {'type': 'number', 'minimum': 0.001, 'maximum': 0.005},
                    'pull_rate': {
                        'type': 'number', 
                        'minimum': 0.001, 
                        'maximum': 0.1,
                        'description': 'Pulling rate in nm/ps'
                    },
                    'pull_k': {
                        'type': 'number',
                        'minimum': 100,
                        'maximum': 10000,
                        'description': 'Force constant in kJ/mol/nm²'
                    },
                    'temperature': {'type': 'number', 'minimum': 250, 'maximum': 450}
                }
            },
            
            'umbrella': {
                'type': 'object',
                'properties': {
                    'sample_interval_near': {
                        'type': 'number',
                        'minimum': 0.05,
                        'maximum': 0.5,
                        'description': 'Sampling interval for near distances in nm'
                    },
                    'sample_interval_far': {
                        'type': 'number',
                        'minimum': 0.1,
                        'maximum': 1.0,
                        'description': 'Sampling interval for far distances in nm'
                    },
                    'force_constant': {
                        'type': 'number',
                        'minimum': 100,
                        'maximum': 10000,
                        'description': 'Umbrella force constant in kJ/mol/nm²'
                    },
                    'production_time_ps': {
                        'type': 'number',
                        'minimum': 5000,
                        'maximum': 100000,
                        'description': 'Production time per window in ps'
                    }
                }
            },
            
            'analysis': {
                'type': 'object',
                'properties': {
                    'begin_time_ps': {
                        'type': 'number',
                        'minimum': 0,
                        'maximum': 50000,
                        'description': 'Skip time for analysis in ps'
                    },
                    'bootstrap_iterations': {
                        'type': 'integer',
                        'minimum': 10,
                        'maximum': 1000
                    },
                    'energy_unit': {
                        'type': 'string',
                        'enum': ['kJ', 'kcal', 'kCal', 'kT']
                    },
                    'temperature': {'type': 'number', 'minimum': 250, 'maximum': 450}
                }
            }
        }
    
    def validate_config(self, config: Dict[str, Any], 
                       level: ValidationLevel = ValidationLevel.STANDARD) -> Tuple[bool, Dict[str, Any]]:
        """
        Validate configuration with specified validation level
        
        Parameters:
        -----------
        config : dict
            Configuration dictionary to validate
        level : ValidationLevel
            Validation strictness level
            
        Returns:
        --------
        Tuple[bool, Dict]
            (is_valid, validation_report)
        """
        self.errors = []
        self.warnings = []
        
        logger.info(f"Validating configuration with {level.value} level")
        
        # Basic validation
        if level.value in ['basic', 'standard', 'strict']:
            self._validate_basic_structure(config)
            self._validate_data_types(config)
        
        # Standard validation
        if level.value in ['standard', 'strict']:
            self._validate_value_ranges(config)
            self._validate_parameter_consistency(config)
        
        # Strict validation
        if level.value == 'strict':
            self._validate_dependencies(config)
            self._validate_scientific_constraints(config)
        
        is_valid = len(self.errors) == 0
        
        validation_report = {
            'is_valid': is_valid,
            'level': level.value,
            'errors': [error.to_dict() if hasattr(error, 'to_dict') else str(error) for error in self.errors],
            'warnings': [warning.to_dict() if hasattr(warning, 'to_dict') else str(warning) for warning in self.warnings],
            'summary': {
                'total_errors': len(self.errors),
                'total_warnings': len(self.warnings),
                'validated_sections': list(config.keys())
            }
        }
        
        return is_valid, validation_report
    
    def _validate_basic_structure(self, config: Dict[str, Any]):
        """Validate basic configuration structure"""
        logger.debug("Validating basic configuration structure")
        
        # Check for empty configuration
        if not config:
            self.errors.append(ConfigValidationError(
                "Configuration is empty",
                error_code="CONFIG_EMPTY"
            ))
            return
        
        # Check for required top-level sections in PMF configs
        if 'smd' in config or 'umbrella' in config or 'analysis' in config:
            # This appears to be a PMF configuration
            pmf_required = ['simulation']
            for section in pmf_required:
                if section not in config:
                    self.errors.append(ConfigValidationError(
                        f"Required PMF section '{section}' is missing",
                        field_path=section,
                        suggestions=[f"Add '{section}' section to configuration"],
                        error_code="MISSING_PMF_SECTION"
                    ))
    
    def _validate_data_types(self, config: Dict[str, Any]):
        """Validate data types according to schema"""
        logger.debug("Validating data types")
        
        for section_name, section_config in config.items():
            if section_name in self.schema:
                section_schema = self.schema[section_name]
                self._validate_section_types(section_config, section_schema, section_name)
    
    def _validate_section_types(self, section: Any, schema: Dict[str, Any], path: str):
        """Validate types for a configuration section"""
        if schema.get('type') == 'object' and isinstance(section, dict):
            properties = schema.get('properties', {})
            
            for field_name, field_value in section.items():
                field_path = f"{path}.{field_name}"
                
                if field_name in properties:
                    field_schema = properties[field_name]
                    self._validate_field_type(field_value, field_schema, field_path)
        elif schema.get('type') == 'object' and not isinstance(section, dict):
            self.errors.append(ConfigValidationError(
                f"Expected object for section '{path}', got {type(section).__name__}",
                field_path=path,
                error_code="TYPE_MISMATCH"
            ))
    
    def _validate_field_type(self, value: Any, schema: Dict[str, Any], path: str):
        """Validate individual field type"""
        expected_type = schema.get('type')
        
        if expected_type == 'string' and not isinstance(value, str):
            self.errors.append(ConfigValidationError(
                f"Expected string for '{path}', got {type(value).__name__}",
                field_path=path,
                error_code="TYPE_MISMATCH"
            ))
        elif expected_type == 'number' and not isinstance(value, (int, float)):
            self.errors.append(ConfigValidationError(
                f"Expected number for '{path}', got {type(value).__name__}",
                field_path=path,
                error_code="TYPE_MISMATCH"
            ))
        elif expected_type == 'integer' and not isinstance(value, int):
            self.errors.append(ConfigValidationError(
                f"Expected integer for '{path}', got {type(value).__name__}",
                field_path=path,
                error_code="TYPE_MISMATCH"
            ))
        elif expected_type == 'boolean' and not isinstance(value, bool):
            self.errors.append(ConfigValidationError(
                f"Expected boolean for '{path}', got {type(value).__name__}",
                field_path=path,
                error_code="TYPE_MISMATCH"
            ))
    
    def _validate_value_ranges(self, config: Dict[str, Any]):
        """Validate value ranges and constraints"""
        logger.debug("Validating value ranges")
        
        for section_name, section_config in config.items():
            if section_name in self.schema and isinstance(section_config, dict):
                section_schema = self.schema[section_name]
                properties = section_schema.get('properties', {})
                
                for field_name, field_value in section_config.items():
                    if field_name in properties:
                        field_schema = properties[field_name]
                        field_path = f"{section_name}.{field_name}"
                        self._validate_field_range(field_value, field_schema, field_path)
    
    def _validate_field_range(self, value: Any, schema: Dict[str, Any], path: str):
        """Validate field value ranges"""
        # Enum validation
        if 'enum' in schema and value not in schema['enum']:
            self.errors.append(ConfigValidationError(
                f"Invalid value '{value}' for '{path}'. Valid values: {schema['enum']}",
                field_path=path,
                suggestions=[f"Use one of: {', '.join(map(str, schema['enum']))}"],
                error_code="INVALID_ENUM_VALUE"
            ))
        
        # Numeric range validation
        if isinstance(value, (int, float)):
            if 'minimum' in schema and value < schema['minimum']:
                self.errors.append(ConfigValidationError(
                    f"Value {value} for '{path}' is below minimum {schema['minimum']}",
                    field_path=path,
                    suggestions=[f"Use value >= {schema['minimum']}"],
                    error_code="VALUE_BELOW_MINIMUM"
                ))
            
            if 'maximum' in schema and value > schema['maximum']:
                self.errors.append(ConfigValidationError(
                    f"Value {value} for '{path}' exceeds maximum {schema['maximum']}",
                    field_path=path,
                    suggestions=[f"Use value <= {schema['maximum']}"],
                    error_code="VALUE_ABOVE_MAXIMUM"
                ))
        
        # String pattern validation
        if isinstance(value, str) and 'pattern' in schema:
            pattern = schema['pattern']
            if not re.match(pattern, value):
                self.errors.append(ConfigValidationError(
                    f"Value '{value}' for '{path}' does not match required pattern",
                    field_path=path,
                    suggestions=[f"Value must match pattern: {pattern}"],
                    error_code="PATTERN_MISMATCH"
                ))
    
    def _validate_parameter_consistency(self, config: Dict[str, Any]):
        """Validate consistency between related parameters"""
        logger.debug("Validating parameter consistency")
        
        # Temperature consistency across sections
        temperatures = {}
        for section in ['simulation', 'smd', 'umbrella', 'equilibration', 'analysis']:
            if section in config:
                section_config = config[section]
                if isinstance(section_config, dict):
                    if 'temperature' in section_config:
                        temperatures[section] = section_config['temperature']
                    # Check nested equilibration sections
                    if section == 'equilibration':
                        for sub_section in ['nvt', 'npt']:
                            if sub_section in section_config and isinstance(section_config[sub_section], dict):
                                if 'temperature' in section_config[sub_section]:
                                    temperatures[f"{section}.{sub_section}"] = section_config[sub_section]['temperature']
        
        if len(set(temperatures.values())) > 1:
            self.warnings.append(ConfigValidationError(
                f"Temperature inconsistency detected: {temperatures}",
                suggestions=["Ensure all temperature values are consistent across sections"],
                error_code="TEMPERATURE_INCONSISTENCY"
            ))
        
        # Time step consistency
        time_steps = {}
        for section in ['smd', 'umbrella', 'equilibration']:
            if section in config and isinstance(config[section], dict):
                if 'dt' in config[section]:
                    time_steps[section] = config[section]['dt']
                # Check nested sections
                if section == 'equilibration':
                    for sub_section in ['nvt', 'npt']:
                        if (sub_section in config[section] and 
                            isinstance(config[section][sub_section], dict) and
                            'dt' in config[section][sub_section]):
                            time_steps[f"{section}.{sub_section}"] = config[section][sub_section]['dt']
        
        if len(set(time_steps.values())) > 1:
            self.warnings.append(ConfigValidationError(
                f"Time step inconsistency detected: {time_steps}",
                suggestions=["Consider using consistent time steps across simulation sections"],
                error_code="TIMESTEP_INCONSISTENCY"
            ))
    
    def _validate_dependencies(self, config: Dict[str, Any]):
        """Validate dependencies between configuration sections"""
        logger.debug("Validating configuration dependencies")
        
        # PMF workflow dependencies
        if 'umbrella' in config and 'smd' not in config:
            self.warnings.append(ConfigValidationError(
                "Umbrella sampling configuration found without SMD configuration",
                suggestions=["Consider adding SMD configuration for complete PMF workflow"],
                error_code="MISSING_SMD_FOR_UMBRELLA"
            ))
        
        if 'analysis' in config and 'umbrella' not in config:
            self.errors.append(ConfigValidationError(
                "Analysis configuration requires umbrella sampling configuration",
                suggestions=["Add umbrella sampling configuration before analysis"],
                error_code="MISSING_UMBRELLA_FOR_ANALYSIS"
            ))
        
        # Force field and water model compatibility
        if 'forcefield' in config and 'water_model' in config:
            ff_name = config['forcefield'].get('name')
            wm_name = config['water_model'].get('name')
            
            # Check known incompatibilities
            incompatible_combinations = [
                ('charmm27', 'tip4p'),  # Example incompatibility
            ]
            
            if (ff_name, wm_name) in incompatible_combinations:
                self.warnings.append(ConfigValidationError(
                    f"Potential incompatibility between force field '{ff_name}' and water model '{wm_name}'",
                    suggestions=["Verify force field and water model compatibility"],
                    error_code="FF_WATER_INCOMPATIBILITY"
                ))
    
    def _validate_scientific_constraints(self, config: Dict[str, Any]):
        """Validate scientific and physical constraints"""
        logger.debug("Validating scientific constraints")
        
        # Check SMD pulling rate - too fast can cause artifacts
        if 'smd' in config and isinstance(config['smd'], dict):
            pull_rate = config['smd'].get('pull_rate')
            if pull_rate and pull_rate > 0.01:
                self.warnings.append(ConfigValidationError(
                    f"SMD pulling rate {pull_rate} nm/ps may be too fast and cause artifacts",
                    field_path="smd.pull_rate",
                    suggestions=["Consider using pulling rate <= 0.01 nm/ps for more accurate results"],
                    error_code="SMD_RATE_TOO_FAST"
                ))
        
        # Check umbrella sampling time - insufficient sampling can lead to poor convergence
        if 'umbrella' in config and isinstance(config['umbrella'], dict):
            prod_time = config['umbrella'].get('production_time_ps')
            if prod_time and prod_time < 10000:  # Less than 10 ns
                self.warnings.append(ConfigValidationError(
                    f"Umbrella production time {prod_time/1000:.1f} ns may be insufficient for convergence",
                    field_path="umbrella.production_time_ps",
                    suggestions=["Consider using at least 10-20 ns production time per window"],
                    error_code="UMBRELLA_TIME_TOO_SHORT"
                ))
        
        # Check analysis begin time vs umbrella equilibration
        if ('analysis' in config and 'umbrella' in config and 
            isinstance(config['analysis'], dict) and isinstance(config['umbrella'], dict)):
            
            begin_time = config['analysis'].get('begin_time_ps', 0)
            equil_time = config['umbrella'].get('equilibration_time_ps', 0)
            
            if begin_time < equil_time:
                self.warnings.append(ConfigValidationError(
                    f"Analysis begin time ({begin_time} ps) is less than umbrella equilibration time ({equil_time} ps)",
                    suggestions=["Set analysis begin_time >= umbrella equilibration_time"],
                    error_code="ANALYSIS_START_TOO_EARLY"
                ))
    
    def generate_validation_report(self, config: Dict[str, Any], 
                                 level: ValidationLevel = ValidationLevel.STANDARD) -> str:
        """Generate detailed validation report"""
        is_valid, report = self.validate_config(config, level)
        
        lines = []
        lines.append("="*60)
        lines.append("PRISM CONFIGURATION VALIDATION REPORT")
        lines.append("="*60)
        lines.append(f"Validation Level: {level.value.upper()}")
        lines.append(f"Overall Status: {'✅ VALID' if is_valid else '❌ INVALID'}")
        lines.append(f"Sections Validated: {', '.join(report['summary']['validated_sections'])}")
        lines.append("")
        
        if report['errors']:
            lines.append(f"❌ ERRORS ({len(report['errors'])})")
            lines.append("-" * 40)
            for i, error in enumerate(report['errors'], 1):
                lines.append(f"{i}. {error.get('message', 'Unknown error')}")
                if error.get('field_path'):
                    lines.append(f"   Field: {error['field_path']}")
                if error.get('suggestions'):
                    lines.append(f"   Suggestions: {'; '.join(error['suggestions'])}")
                lines.append("")
        
        if report['warnings']:
            lines.append(f"⚠️  WARNINGS ({len(report['warnings'])})")
            lines.append("-" * 40)
            for i, warning in enumerate(report['warnings'], 1):
                lines.append(f"{i}. {warning.get('message', 'Unknown warning')}")
                if warning.get('suggestions'):
                    lines.append(f"   Suggestions: {'; '.join(warning['suggestions'])}")
                lines.append("")
        
        if is_valid:
            lines.append("✅ Configuration is valid and ready for use!")
        else:
            lines.append("❌ Please fix the errors before proceeding.")
        
        lines.append("="*60)
        
        return "\n".join(lines)


def validate_config_file(config_path: str, 
                        level: ValidationLevel = ValidationLevel.STANDARD) -> Tuple[bool, Dict[str, Any]]:
    """
    Validate configuration file
    
    Parameters:
    -----------
    config_path : str
        Path to configuration file
    level : ValidationLevel
        Validation level
        
    Returns:
    --------
    Tuple[bool, Dict]
        (is_valid, validation_report)
    """
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        validator = ConfigurationValidator()
        return validator.validate_config(config, level)
        
    except FileNotFoundError:
        return False, {
            'is_valid': False,
            'errors': [{'message': f"Configuration file not found: {config_path}"}],
            'warnings': [],
            'summary': {'total_errors': 1, 'total_warnings': 0}
        }
    except yaml.YAMLError as e:
        return False, {
            'is_valid': False,
            'errors': [{'message': f"YAML parsing error: {e}"}],
            'warnings': [],
            'summary': {'total_errors': 1, 'total_warnings': 0}
        }


# Convenience functions
def quick_validate(config: Dict[str, Any]) -> bool:
    """Quick validation with basic level"""
    validator = ConfigurationValidator()
    is_valid, _ = validator.validate_config(config, ValidationLevel.BASIC)
    return is_valid


def validate_pmf_config(config: Dict[str, Any]) -> Tuple[bool, str]:
    """Validate PMF configuration and return report"""
    validator = ConfigurationValidator()
    is_valid, _ = validator.validate_config(config, ValidationLevel.STANDARD)
    report = validator.generate_validation_report(config, ValidationLevel.STANDARD)
    return is_valid, report