#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Configuration Templates Generator

Provides template generation for different types of PRISM calculations
with appropriate default values and documentation.
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class ConfigTemplate(Enum):
    """Available configuration templates"""
    BASIC_MD = "basic_md"                # Basic MD simulation
    PMF_FAST = "pmf_fast"                # Fast PMF calculation (shorter times)
    PMF_STANDARD = "pmf_standard"        # Standard PMF calculation
    PMF_ACCURATE = "pmf_accurate"        # High-accuracy PMF calculation
    PMF_PRODUCTION = "pmf_production"    # Production-ready PMF calculation
    DRUG_DISCOVERY = "drug_discovery"    # Drug discovery optimized
    MEMBRANE_PROTEIN = "membrane_protein" # Membrane protein systems


class ConfigTemplateGenerator:
    """Generate configuration templates for different use cases"""
    
    def __init__(self):
        self.base_config = self._get_base_configuration()
    
    def _get_base_configuration(self) -> Dict[str, Any]:
        """Base configuration with common settings"""
        return {
            'general': {
                'overwrite': False,
                'gmx_command': 'gmx',
                'project_name': 'prism_simulation'
            },
            'forcefield': {
                'name': 'amber99sb-ildn',
                'protein_forcefield': 'amber99sb-ildn',
                'ligand_forcefield': 'gaff'
            },
            'water_model': {
                'name': 'tip3p'
            },
            'box': {
                'distance': 1.2,
                'shape': 'cubic',
                'center': True
            },
            'simulation': {
                'temperature': 310.0,
                'pressure': 1.0,
                'pH': 7.0
            },
            'ions': {
                'neutral': True,
                'concentration': 0.15,
                'positive_ion': 'NA',
                'negative_ion': 'CL'
            },
            'constraints': {
                'algorithm': 'lincs',
                'type': 'h-bonds',
                'lincs_iter': 1,
                'lincs_order': 4
            },
            'output': {
                'log_level': 'INFO',
                'save_configs': True,
                'validate_inputs': True
            }
        }
    
    def generate_template(self, template_type: ConfigTemplate, 
                         custom_params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Generate configuration template
        
        Parameters:
        -----------
        template_type : ConfigTemplate
            Type of template to generate
        custom_params : dict, optional
            Custom parameters to override defaults
            
        Returns:
        --------
        Dict : Generated configuration
        """
        config = self.base_config.copy()
        
        # Apply template-specific settings
        if template_type == ConfigTemplate.BASIC_MD:
            config.update(self._get_basic_md_config())
        elif template_type == ConfigTemplate.PMF_FAST:
            config.update(self._get_pmf_fast_config())
        elif template_type == ConfigTemplate.PMF_STANDARD:
            config.update(self._get_pmf_standard_config())
        elif template_type == ConfigTemplate.PMF_ACCURATE:
            config.update(self._get_pmf_accurate_config())
        elif template_type == ConfigTemplate.PMF_PRODUCTION:
            config.update(self._get_pmf_production_config())
        elif template_type == ConfigTemplate.DRUG_DISCOVERY:
            config.update(self._get_drug_discovery_config())
        elif template_type == ConfigTemplate.MEMBRANE_PROTEIN:
            config.update(self._get_membrane_protein_config())
        else:
            raise ValueError(f"Unknown template type: {template_type}")
        
        # Apply custom parameters
        if custom_params:
            config = self._deep_merge(config, custom_params)
        
        # Add metadata
        config['_template_info'] = {
            'template_type': template_type.value,
            'generated_by': 'PRISM ConfigTemplateGenerator',
            'description': self._get_template_description(template_type)
        }
        
        return config
    
    def _get_basic_md_config(self) -> Dict[str, Any]:
        """Basic MD simulation configuration"""
        return {
            'simulation': {
                'production_time_ns': 100,
                'dt': 0.002,
                'equilibration_nvt_time_ps': 1000,
                'equilibration_npt_time_ps': 1000
            },
            'energy_minimization': {
                'integrator': 'steep',
                'emtol': 1000.0,
                'emstep': 0.01,
                'nsteps': 50000
            },
            'output': {
                'trajectory_interval_ps': 1000,
                'energy_interval_ps': 100,
                'log_interval_ps': 100
            }
        }
    
    def _get_pmf_fast_config(self) -> Dict[str, Any]:
        """Fast PMF calculation (for testing/development)"""
        return {
            'builder': {
                'rebuild_system': True,
                'frame_number': -1,
                'z_extension': 2.0,
                'box_type': 'cubic'
            },
            'equilibration': {
                'em': {
                    'integrator': 'steep',
                    'emtol': 100.0,
                    'nsteps': 25000,
                    'posres_protein': 1000.0,
                    'posres_ligand': 500.0
                },
                'nvt': {
                    'nsteps': 25000,  # 50 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'posres_protein': 1000.0,
                    'posres_ligand': 500.0
                },
                'npt': {
                    'nsteps': 50000,  # 100 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'pressure': 1.0,
                    'posres_protein': 400.0,
                    'posres_ligand': 200.0
                }
            },
            'smd': {
                'nsteps': 1000000,  # 2 ns
                'dt': 0.002,
                'pull_rate': 0.01,  # Faster for testing
                'pull_k': 1000.0,
                'temperature': 310.0,
                'nstpull': 50
            },
            'umbrella': {
                'sample_interval_near': 0.2,
                'sample_interval_far': 0.3,
                'force_constant': 1000.0,
                'production_time_ps': 5000,  # 5 ns per window
                'equilibration_time_ps': 1000  # 1 ns equilibration
            },
            'analysis': {
                'begin_time_ps': 1000,
                'bootstrap_iterations': 25,
                'energy_unit': 'kcal'
            }
        }
    
    def _get_pmf_standard_config(self) -> Dict[str, Any]:
        """Standard PMF calculation (balanced accuracy/speed)"""
        return {
            'builder': {
                'rebuild_system': True,
                'frame_number': -1,
                'z_extension': 2.5,
                'box_type': 'cubic',
                'align_axis': 'z',
                'center_system': True
            },
            'equilibration': {
                'em': {
                    'integrator': 'steep',
                    'emtol': 10.0,
                    'nsteps': 50000,
                    'posres_protein': 1000.0,
                    'posres_ligand': 500.0
                },
                'nvt': {
                    'nsteps': 50000,  # 100 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'posres_protein': 1000.0,
                    'posres_ligand': 500.0
                },
                'npt': {
                    'nsteps': 100000,  # 200 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'pressure': 1.0,
                    'posres_protein': 400.0,
                    'posres_ligand': 200.0
                }
            },
            'smd': {
                'nsteps': 2500000,  # 5 ns
                'dt': 0.002,
                'pull_rate': 0.005,
                'pull_k': 1000.0,
                'temperature': 310.0,
                'nstpull': 50
            },
            'umbrella': {
                'sample_interval_near': 0.1,
                'sample_interval_far': 0.2,
                'force_constant': 1000.0,
                'production_time_ps': 15000,  # 15 ns per window
                'equilibration_time_ps': 2000  # 2 ns equilibration
            },
            'analysis': {
                'begin_time_ps': 2000,
                'bootstrap_iterations': 50,
                'energy_unit': 'kcal',
                'generate_plots': True
            }
        }
    
    def _get_pmf_accurate_config(self) -> Dict[str, Any]:
        """High-accuracy PMF calculation (longer simulations)"""
        return {
            'builder': {
                'rebuild_system': True,
                'frame_number': -1,
                'z_extension': 3.0,
                'box_type': 'cubic',
                'align_axis': 'z',
                'center_system': True
            },
            'equilibration': {
                'em': {
                    'integrator': 'steep',
                    'emtol': 1.0,
                    'nsteps': 100000,
                    'posres_protein': 1000.0,
                    'posres_ligand': 500.0
                },
                'nvt': {
                    'nsteps': 100000,  # 200 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'posres_protein': 1000.0,
                    'posres_ligand': 500.0
                },
                'npt': {
                    'nsteps': 250000,  # 500 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'pressure': 1.0,
                    'posres_protein': 400.0,
                    'posres_ligand': 200.0
                }
            },
            'smd': {
                'nsteps': 5000000,  # 10 ns
                'dt': 0.002,
                'pull_rate': 0.002,  # Slower for accuracy
                'pull_k': 1000.0,
                'temperature': 310.0,
                'nstpull': 25  # More frequent output
            },
            'umbrella': {
                'sample_interval_near': 0.08,
                'sample_interval_far': 0.15,
                'force_constant': 1000.0,
                'production_time_ps': 30000,  # 30 ns per window
                'equilibration_time_ps': 5000  # 5 ns equilibration
            },
            'analysis': {
                'begin_time_ps': 5000,
                'bootstrap_iterations': 100,
                'energy_unit': 'kcal',
                'generate_plots': True,
                'tolerance': 1e-8,
                'max_iterations': 2000
            }
        }
    
    def _get_pmf_production_config(self) -> Dict[str, Any]:
        """Production-ready PMF calculation (publication quality)"""
        return {
            'builder': {
                'rebuild_system': True,
                'frame_number': -1,
                'z_extension': 3.0,
                'box_type': 'cubic',
                'align_axis': 'z',
                'center_system': True
            },
            'equilibration': {
                'em': {
                    'integrator': 'steep',
                    'emtol': 1.0,
                    'nsteps': 100000,
                    'posres_protein': 1000.0,
                    'posres_ligand': 500.0
                },
                'nvt': {
                    'nsteps': 250000,  # 500 ps
                    'dt': 0.002,
                    'temperature': 310.0,
                    'posres_protein': 500.0,
                    'posres_ligand': 250.0
                },
                'npt': {
                    'nsteps': 500000,  # 1 ns
                    'dt': 0.002,
                    'temperature': 310.0,
                    'pressure': 1.0,
                    'posres_protein': 200.0,
                    'posres_ligand': 100.0
                }
            },
            'smd': {
                'nsteps': 10000000,  # 20 ns
                'dt': 0.002,
                'pull_rate': 0.001,  # Very slow for accuracy
                'pull_k': 1000.0,
                'temperature': 310.0,
                'nstpull': 25
            },
            'umbrella': {
                'sample_interval_near': 0.05,
                'sample_interval_far': 0.1,
                'force_constant': 1000.0,
                'production_time_ps': 50000,  # 50 ns per window
                'equilibration_time_ps': 10000  # 10 ns equilibration
            },
            'analysis': {
                'begin_time_ps': 10000,
                'bootstrap_iterations': 200,
                'energy_unit': 'kcal',
                'generate_plots': True,
                'tolerance': 1e-8,
                'max_iterations': 5000
            },
            'advanced': {
                'use_gpu': True,
                'max_retries': 5,
                'continue_on_error': False
            }
        }
    
    def _get_drug_discovery_config(self) -> Dict[str, Any]:
        """Drug discovery optimized configuration"""
        config = self._get_pmf_standard_config()
        config.update({
            'simulation': {
                'temperature': 310.0,  # Body temperature
                'pressure': 1.0,
                'pH': 7.4  # Physiological pH
            },
            'ions': {
                'neutral': True,
                'concentration': 0.15,  # Physiological salt concentration
                'positive_ion': 'NA',
                'negative_ion': 'CL'
            },
            'smd': {
                'nsteps': 2500000,  # 5 ns - balanced for throughput
                'pull_rate': 0.005,  # Standard rate
                'temperature': 310.0
            },
            'umbrella': {
                'production_time_ps': 20000,  # 20 ns - good for drug discovery
                'equilibration_time_ps': 2000
            },
            'analysis': {
                'energy_unit': 'kcal',  # Preferred in drug discovery
                'generate_plots': True,
                'bootstrap_iterations': 50
            }
        })
        return config
    
    def _get_membrane_protein_config(self) -> Dict[str, Any]:
        """Membrane protein specific configuration"""
        config = self._get_pmf_standard_config()
        config.update({
            'box': {
                'distance': 2.0,  # Larger box for membrane
                'shape': 'cubic'
            },
            'simulation': {
                'temperature': 310.0,
                'pressure': 1.0
            },
            'equilibration': {
                'em': {
                    'emtol': 10.0,
                    'nsteps': 100000,
                    'posres_protein': 2000.0,  # Stronger restraints for membrane
                    'posres_ligand': 1000.0
                },
                'nvt': {
                    'nsteps': 100000,  # Longer equilibration
                    'posres_protein': 1000.0,
                    'posres_ligand': 500.0
                },
                'npt': {
                    'nsteps': 250000,  # Extended NPT
                    'pcoupl': 'semiisotropic',  # Better for membranes
                    'posres_protein': 500.0,
                    'posres_ligand': 250.0
                }
            },
            'smd': {
                'pull_rate': 0.002,  # Slower for membrane systems
                'nsteps': 5000000  # Longer SMD
            },
            'umbrella': {
                'production_time_ps': 25000,  # Longer sampling
                'equilibration_time_ps': 5000
            }
        })
        return config
    
    def _get_template_description(self, template_type: ConfigTemplate) -> str:
        """Get description for template type"""
        descriptions = {
            ConfigTemplate.BASIC_MD: "Basic molecular dynamics simulation configuration",
            ConfigTemplate.PMF_FAST: "Fast PMF calculation for testing and development",
            ConfigTemplate.PMF_STANDARD: "Standard PMF calculation with balanced accuracy and speed",
            ConfigTemplate.PMF_ACCURATE: "High-accuracy PMF calculation with extended sampling",
            ConfigTemplate.PMF_PRODUCTION: "Production-ready PMF calculation for publication",
            ConfigTemplate.DRUG_DISCOVERY: "Drug discovery optimized PMF calculation",
            ConfigTemplate.MEMBRANE_PROTEIN: "Membrane protein specific PMF configuration"
        }
        return descriptions.get(template_type, "Unknown template type")
    
    def _deep_merge(self, base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
        """Deep merge two dictionaries"""
        result = base.copy()
        
        for key, value in override.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = self._deep_merge(result[key], value)
            else:
                result[key] = value
        
        return result
    
    def save_template(self, template_type: ConfigTemplate, 
                     output_path: str, custom_params: Optional[Dict[str, Any]] = None):
        """
        Generate and save configuration template to file
        
        Parameters:
        -----------
        template_type : ConfigTemplate
            Type of template to generate
        output_path : str
            Output file path
        custom_params : dict, optional
            Custom parameters to override
        """
        config = self.generate_template(template_type, custom_params)
        
        # Add header comment
        header = self._generate_header_comment(template_type)
        
        # Save to file
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            f.write(header)
            yaml.dump(config, f, default_flow_style=False, sort_keys=False, indent=2)
        
        logger.info(f"Configuration template saved: {output_file}")
    
    def _generate_header_comment(self, template_type: ConfigTemplate) -> str:
        """Generate header comment for configuration file"""
        description = self._get_template_description(template_type)
        
        header = f"""# PRISM Configuration File
# Template Type: {template_type.value.upper()}
# Description: {description}
# Generated by PRISM ConfigTemplateGenerator
#
# This configuration file contains optimized parameters for {template_type.value} calculations.
# You can modify the values below to suit your specific system and requirements.
#
# For more information about PRISM configuration, see the documentation.

"""
        return header
    
    def list_templates(self) -> Dict[str, str]:
        """List available templates with descriptions"""
        return {
            template.value: self._get_template_description(template)
            for template in ConfigTemplate
        }


# Convenience functions
def create_pmf_config(config_type: str = "standard", 
                     output_file: Optional[str] = None,
                     custom_params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Create PMF configuration
    
    Parameters:
    -----------
    config_type : str
        Configuration type: 'fast', 'standard', 'accurate', 'production'
    output_file : str, optional
        Output file path to save configuration
    custom_params : dict, optional
        Custom parameters to override
        
    Returns:
    --------
    Dict : Generated configuration
    """
    generator = ConfigTemplateGenerator()
    
    template_map = {
        'fast': ConfigTemplate.PMF_FAST,
        'standard': ConfigTemplate.PMF_STANDARD,
        'accurate': ConfigTemplate.PMF_ACCURATE,
        'production': ConfigTemplate.PMF_PRODUCTION
    }
    
    if config_type not in template_map:
        raise ValueError(f"Unknown config type: {config_type}. Available: {list(template_map.keys())}")
    
    template_type = template_map[config_type]
    config = generator.generate_template(template_type, custom_params)
    
    if output_file:
        generator.save_template(template_type, output_file, custom_params)
    
    return config


def create_drug_discovery_config(output_file: Optional[str] = None,
                                ligand_charge: int = 0,
                                custom_params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Create drug discovery optimized configuration
    
    Parameters:
    -----------
    output_file : str, optional
        Output file path
    ligand_charge : int
        Ligand charge
    custom_params : dict, optional
        Custom parameters
        
    Returns:
    --------
    Dict : Generated configuration
    """
    generator = ConfigTemplateGenerator()
    
    # Add drug discovery specific parameters
    dd_params = {
        'simulation': {
            'ligand_charge': ligand_charge,
            'temperature': 310.0,  # Body temperature
            'pH': 7.4  # Physiological pH
        }
    }
    
    if custom_params:
        dd_params = generator._deep_merge(dd_params, custom_params)
    
    config = generator.generate_template(ConfigTemplate.DRUG_DISCOVERY, dd_params)
    
    if output_file:
        generator.save_template(ConfigTemplate.DRUG_DISCOVERY, output_file, dd_params)
    
    return config