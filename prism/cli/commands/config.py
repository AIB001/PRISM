#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM CLI Configuration Commands

Command handlers for configuration management, validation, and interactive setup.
"""

import yaml
import json
from pathlib import Path
from typing import Dict, Any, List, Optional
from argparse import ArgumentParser, Namespace

from ...utils.logging_system import PrismLogger
from ..utils import CLIFormatter, TableFormatter, confirm_action


class ConfigCommands:
    """Handler for configuration-related CLI commands"""
    
    def __init__(self, context):
        self.context = context
        self.logger = PrismLogger("config_cli")
        self.formatter = CLIFormatter()
    
    def setup_parser(self, parser: ArgumentParser):
        """Setup config command parser"""
        subparsers = parser.add_subparsers(dest='config_action', help='Configuration actions')
        
        # Show configuration
        show_parser = subparsers.add_parser('show', help='Show current configuration')
        show_parser.add_argument('--section', help='Show specific configuration section')
        show_parser.add_argument('--key', help='Show specific configuration key')
        
        # Validate configuration
        validate_parser = subparsers.add_parser('validate', help='Validate configuration file')
        validate_parser.add_argument('config_file', nargs='?', help='Configuration file to validate')
        validate_parser.add_argument('--strict', action='store_true', help='Use strict validation')
        
        # Set configuration values
        set_parser = subparsers.add_parser('set', help='Set configuration value')
        set_parser.add_argument('key', help='Configuration key (e.g., workflow.default_priority)')
        set_parser.add_argument('value', help='Configuration value')
        set_parser.add_argument('--global', action='store_true', help='Set global configuration')
        
        # Get configuration values
        get_parser = subparsers.add_parser('get', help='Get configuration value')
        get_parser.add_argument('key', help='Configuration key')
        
        # Reset configuration
        reset_parser = subparsers.add_parser('reset', help='Reset configuration to defaults')
        reset_parser.add_argument('--section', help='Reset specific section')
        reset_parser.add_argument('--force', action='store_true', help='Force reset without confirmation')
        
        # Initialize configuration
        init_parser = subparsers.add_parser('init', help='Initialize configuration interactively')
        init_parser.add_argument('--template', choices=['basic', 'advanced', 'cluster'],
                               default='basic', help='Configuration template')
        
        # Import/Export configuration
        export_parser = subparsers.add_parser('export', help='Export configuration')
        export_parser.add_argument('output_file', help='Output file path')
        export_parser.add_argument('--format', choices=['yaml', 'json'], default='yaml', 
                                 help='Export format')
        
        import_parser = subparsers.add_parser('import', help='Import configuration')
        import_parser.add_argument('config_file', help='Configuration file to import')
        import_parser.add_argument('--merge', action='store_true', 
                                 help='Merge with existing configuration')
        
        # Template management
        template_parser = subparsers.add_parser('template', help='Configuration template management')
        template_subparsers = template_parser.add_subparsers(dest='template_action')
        
        template_list = template_subparsers.add_parser('list', help='List configuration templates')
        template_show = template_subparsers.add_parser('show', help='Show template details')
        template_show.add_argument('template_name', help='Template name')
    
    def handle(self, args: Namespace) -> int:
        """Handle configuration commands"""
        try:
            if args.config_action == 'show':
                return self._show_config(args)
            elif args.config_action == 'validate':
                return self._validate_config(args)
            elif args.config_action == 'set':
                return self._set_config(args)
            elif args.config_action == 'get':
                return self._get_config(args)
            elif args.config_action == 'reset':
                return self._reset_config(args)
            elif args.config_action == 'init':
                return self._init_config(args)
            elif args.config_action == 'export':
                return self._export_config(args)
            elif args.config_action == 'import':
                return self._import_config(args)
            elif args.config_action == 'template':
                return self._handle_templates(args)
            else:
                self.formatter.error("No configuration action specified")
                return 1
                
        except Exception as e:
            self.formatter.error(f"Configuration command failed: {e}")
            if self.context.verbose:
                self.logger.exception("Configuration command error")
            return 1
    
    def _show_config(self, args: Namespace) -> int:
        """Show current configuration"""
        config = self._load_current_config()
        
        if args.key:
            # Show specific key
            value = self._get_nested_value(config, args.key)
            if value is not None:
                if isinstance(value, (dict, list)):
                    print(yaml.dump({args.key: value}, default_flow_style=False))
                else:
                    print(f"{args.key}: {value}")
            else:
                self.formatter.error(f"Configuration key not found: {args.key}")
                return 1
        elif args.section:
            # Show specific section
            section = config.get(args.section)
            if section is not None:
                print(yaml.dump({args.section: section}, default_flow_style=False))
            else:
                self.formatter.error(f"Configuration section not found: {args.section}")
                return 1
        else:
            # Show all configuration
            self.formatter.header("Current Configuration")
            
            if self.context.output_format == 'json':
                print(json.dumps(config, indent=2))
            else:
                print(yaml.dump(config, default_flow_style=False))
        
        return 0
    
    def _validate_config(self, args: Namespace) -> int:
        """Validate configuration file"""
        if args.config_file:
            config_path = Path(args.config_file)
        else:
            config_path = self._get_default_config_path()
        
        if not config_path.exists():
            self.formatter.error(f"Configuration file not found: {config_path}")
            return 1
        
        self.formatter.info(f"Validating configuration: {config_path}")
        
        try:
            # Load and parse configuration
            with open(config_path, 'r') as f:
                if config_path.suffix in ['.yaml', '.yml']:
                    config = yaml.safe_load(f)
                else:
                    config = json.load(f)
            
            # Validate configuration structure
            errors = self._validate_config_structure(config, strict=args.strict)
            
            if errors:
                self.formatter.error(f"Configuration validation failed ({len(errors)} errors):")
                for error in errors:
                    print(f"  • {error}")
                return 1
            else:
                self.formatter.success("Configuration validation passed")
                return 0
                
        except Exception as e:
            self.formatter.error(f"Failed to validate configuration: {e}")
            return 1
    
    def _set_config(self, args: Namespace) -> int:
        """Set configuration value"""
        config = self._load_current_config()
        
        # Parse value
        value = self._parse_config_value(args.value)
        
        # Set nested value
        self._set_nested_value(config, args.key, value)
        
        # Save configuration
        config_path = self._get_config_path(args)
        self._save_config(config, config_path)
        
        self.formatter.success(f"Set {args.key} = {value}")
        return 0
    
    def _get_config(self, args: Namespace) -> int:
        """Get configuration value"""
        config = self._load_current_config()
        
        value = self._get_nested_value(config, args.key)
        
        if value is not None:
            if isinstance(value, (dict, list)):
                if self.context.output_format == 'json':
                    print(json.dumps(value, indent=2))
                else:
                    print(yaml.dump(value, default_flow_style=False))
            else:
                print(value)
        else:
            self.formatter.error(f"Configuration key not found: {args.key}")
            return 1
        
        return 0
    
    def _reset_config(self, args: Namespace) -> int:
        """Reset configuration to defaults"""
        if not args.force:
            if args.section:
                message = f"Reset section '{args.section}' to defaults?"
            else:
                message = "Reset entire configuration to defaults?"
            
            if not confirm_action(message):
                self.formatter.info("Reset cancelled")
                return 0
        
        # Load default configuration
        default_config = self._get_default_config()
        
        if args.section:
            # Reset specific section
            current_config = self._load_current_config()
            if args.section in default_config:
                current_config[args.section] = default_config[args.section]
                config_to_save = current_config
                self.formatter.success(f"Reset section '{args.section}' to defaults")
            else:
                self.formatter.error(f"Unknown configuration section: {args.section}")
                return 1
        else:
            # Reset entire configuration
            config_to_save = default_config
            self.formatter.success("Reset configuration to defaults")
        
        # Save configuration
        config_path = self._get_default_config_path()
        self._save_config(config_to_save, config_path)
        
        return 0
    
    def _init_config(self, args: Namespace) -> int:
        """Initialize configuration interactively"""
        self.formatter.header(f"Interactive Configuration Setup ({args.template} template)")
        
        config = self._get_template_config(args.template)
        
        # Interactive configuration
        if args.template == 'basic':
            config = self._interactive_basic_config(config)
        elif args.template == 'advanced':
            config = self._interactive_advanced_config(config)
        elif args.template == 'cluster':
            config = self._interactive_cluster_config(config)
        
        # Save configuration
        config_path = self._get_default_config_path()
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        self._save_config(config, config_path)
        
        self.formatter.success(f"Configuration saved to: {config_path}")
        return 0
    
    def _export_config(self, args: Namespace) -> int:
        """Export configuration"""
        config = self._load_current_config()
        output_path = Path(args.output_file)
        
        try:
            with open(output_path, 'w') as f:
                if args.format == 'json':
                    json.dump(config, f, indent=2)
                else:
                    yaml.dump(config, f, default_flow_style=False)
            
            self.formatter.success(f"Configuration exported to: {output_path}")
            return 0
            
        except Exception as e:
            self.formatter.error(f"Failed to export configuration: {e}")
            return 1
    
    def _import_config(self, args: Namespace) -> int:
        """Import configuration"""
        import_path = Path(args.config_file)
        
        if not import_path.exists():
            self.formatter.error(f"Configuration file not found: {import_path}")
            return 1
        
        try:
            # Load configuration to import
            with open(import_path, 'r') as f:
                if import_path.suffix in ['.yaml', '.yml']:
                    new_config = yaml.safe_load(f)
                else:
                    new_config = json.load(f)
            
            if args.merge:
                # Merge with existing configuration
                current_config = self._load_current_config()
                merged_config = self._merge_configs(current_config, new_config)
                config_to_save = merged_config
                action = "merged"
            else:
                # Replace configuration
                config_to_save = new_config
                action = "imported"
            
            # Save configuration
            config_path = self._get_default_config_path()
            self._save_config(config_to_save, config_path)
            
            self.formatter.success(f"Configuration {action} from: {import_path}")
            return 0
            
        except Exception as e:
            self.formatter.error(f"Failed to import configuration: {e}")
            return 1
    
    def _handle_templates(self, args: Namespace) -> int:
        """Handle template management"""
        if args.template_action == 'list':
            return self._list_templates()
        elif args.template_action == 'show':
            return self._show_template(args.template_name)
        else:
            self.formatter.error("No template action specified")
            return 1
    
    def _list_templates(self) -> int:
        """List available configuration templates"""
        templates = [
            {
                'name': 'basic',
                'description': 'Basic single-node configuration',
                'use_case': 'Development and small-scale calculations'
            },
            {
                'name': 'advanced',
                'description': 'Advanced configuration with optimization',
                'use_case': 'Production single-node systems'
            },
            {
                'name': 'cluster',
                'description': 'Multi-node cluster configuration',
                'use_case': 'High-performance distributed computing'
            }
        ]
        
        table = TableFormatter(['Template', 'Description', 'Use Case'])
        
        for template in templates:
            table.add_row([
                template['name'],
                template['description'],
                template['use_case']
            ])
        
        self.formatter.header("Available Configuration Templates")
        table.print(self.context.output_format)
        
        return 0
    
    def _show_template(self, template_name: str) -> int:
        """Show template details"""
        template_config = self._get_template_config(template_name)
        
        if not template_config:
            self.formatter.error(f"Template not found: {template_name}")
            return 1
        
        self.formatter.header(f"Configuration Template: {template_name}")
        
        if self.context.output_format == 'json':
            print(json.dumps(template_config, indent=2))
        else:
            print(yaml.dump(template_config, default_flow_style=False))
        
        return 0
    
    # Helper methods
    def _load_current_config(self) -> Dict[str, Any]:
        """Load current configuration"""
        config_path = self.context.config_path or self._get_default_config_path()
        
        if config_path.exists():
            try:
                with open(config_path, 'r') as f:
                    if config_path.suffix in ['.yaml', '.yml']:
                        return yaml.safe_load(f) or {}
                    else:
                        return json.load(f)
            except Exception as e:
                self.logger.warning(f"Failed to load config from {config_path}: {e}")
        
        return self._get_default_config()
    
    def _get_default_config_path(self) -> Path:
        """Get default configuration file path"""
        return Path.home() / '.prism' / 'config.yaml'
    
    def _get_config_path(self, args: Namespace) -> Path:
        """Get configuration file path"""
        if hasattr(args, 'global') and getattr(args, 'global'):
            return Path('/etc/prism/config.yaml')
        else:
            return self._get_default_config_path()
    
    def _get_default_config(self) -> Dict[str, Any]:
        """Get default configuration"""
        return {
            'system': {
                'max_concurrent_workflows': 10,
                'default_priority': 'normal',
                'log_level': 'INFO',
                'working_directory': str(Path.cwd())
            },
            'workflow': {
                'default_template': 'pmf',
                'auto_submit': False,
                'retry_failed_tasks': True,
                'max_retries': 3
            },
            'resources': {
                'default_cpu_cores': 4,
                'default_memory_gb': 8.0,
                'max_memory_utilization': 0.9,
                'resource_check_interval': 30
            },
            'monitoring': {
                'enabled': True,
                'dashboard_update_interval': 5,
                'alert_email': None,
                'log_retention_days': 30
            },
            'storage': {
                'data_directory': str(Path.cwd() / 'prism_data'),
                'backup_enabled': True,
                'compression_enabled': True,
                'cleanup_old_files': True
            }
        }
    
    def _get_template_config(self, template_name: str) -> Optional[Dict[str, Any]]:
        """Get template configuration"""
        templates = {
            'basic': {
                'system': {
                    'max_concurrent_workflows': 5,
                    'default_priority': 'normal',
                    'log_level': 'INFO'
                },
                'workflow': {
                    'default_template': 'pmf',
                    'auto_submit': False,
                    'retry_failed_tasks': True,
                    'max_retries': 2
                },
                'resources': {
                    'default_cpu_cores': 2,
                    'default_memory_gb': 4.0
                },
                'monitoring': {
                    'enabled': True,
                    'dashboard_update_interval': 10
                }
            },
            'advanced': {
                'system': {
                    'max_concurrent_workflows': 20,
                    'default_priority': 'normal',
                    'log_level': 'INFO',
                    'performance_optimization': True
                },
                'workflow': {
                    'default_template': 'pmf',
                    'auto_submit': True,
                    'retry_failed_tasks': True,
                    'max_retries': 5,
                    'intelligent_scheduling': True
                },
                'resources': {
                    'default_cpu_cores': 8,
                    'default_memory_gb': 16.0,
                    'adaptive_resource_allocation': True
                },
                'monitoring': {
                    'enabled': True,
                    'dashboard_update_interval': 5,
                    'advanced_metrics': True,
                    'predictive_alerts': True
                }
            },
            'cluster': {
                'system': {
                    'max_concurrent_workflows': 100,
                    'default_priority': 'normal',
                    'log_level': 'INFO',
                    'distributed_mode': True
                },
                'workflow': {
                    'default_template': 'pmf',
                    'auto_submit': True,
                    'retry_failed_tasks': True,
                    'max_retries': 5,
                    'load_balancing': True
                },
                'resources': {
                    'default_cpu_cores': 16,
                    'default_memory_gb': 32.0,
                    'cluster_resource_management': True
                },
                'cluster': {
                    'node_discovery': 'auto',
                    'load_balancing_strategy': 'resource_aware',
                    'fault_tolerance': True,
                    'auto_scaling': True
                },
                'monitoring': {
                    'enabled': True,
                    'dashboard_update_interval': 2,
                    'cluster_monitoring': True,
                    'distributed_logging': True
                }
            }
        }
        
        return templates.get(template_name)
    
    def _validate_config_structure(self, config: Dict[str, Any], strict: bool = False) -> List[str]:
        """Validate configuration structure"""
        errors = []
        
        # Required sections
        required_sections = ['system', 'workflow', 'resources']
        for section in required_sections:
            if section not in config:
                errors.append(f"Missing required section: {section}")
        
        # Validate system section
        if 'system' in config:
            system = config['system']
            
            if 'max_concurrent_workflows' in system:
                if not isinstance(system['max_concurrent_workflows'], int) or system['max_concurrent_workflows'] <= 0:
                    errors.append("system.max_concurrent_workflows must be a positive integer")
            
            if 'default_priority' in system:
                valid_priorities = ['low', 'normal', 'high', 'critical']
                if system['default_priority'] not in valid_priorities:
                    errors.append(f"system.default_priority must be one of: {', '.join(valid_priorities)}")
        
        # Validate resources section
        if 'resources' in config:
            resources = config['resources']
            
            if 'default_cpu_cores' in resources:
                if not isinstance(resources['default_cpu_cores'], int) or resources['default_cpu_cores'] <= 0:
                    errors.append("resources.default_cpu_cores must be a positive integer")
            
            if 'default_memory_gb' in resources:
                if not isinstance(resources['default_memory_gb'], (int, float)) or resources['default_memory_gb'] <= 0:
                    errors.append("resources.default_memory_gb must be a positive number")
        
        return errors
    
    def _get_nested_value(self, config: Dict[str, Any], key: str) -> Any:
        """Get nested configuration value using dot notation"""
        keys = key.split('.')
        value = config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return None
        
        return value
    
    def _set_nested_value(self, config: Dict[str, Any], key: str, value: Any):
        """Set nested configuration value using dot notation"""
        keys = key.split('.')
        current = config
        
        for k in keys[:-1]:
            if k not in current:
                current[k] = {}
            current = current[k]
        
        current[keys[-1]] = value
    
    def _parse_config_value(self, value_str: str) -> Any:
        """Parse configuration value from string"""
        # Try to parse as JSON first (handles booleans, numbers, lists, etc.)
        try:
            return json.loads(value_str)
        except json.JSONDecodeError:
            # Return as string if JSON parsing fails
            return value_str
    
    def _save_config(self, config: Dict[str, Any], config_path: Path):
        """Save configuration to file"""
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, indent=2)
    
    def _merge_configs(self, base: Dict[str, Any], overlay: Dict[str, Any]) -> Dict[str, Any]:
        """Merge two configuration dictionaries"""
        result = base.copy()
        
        for key, value in overlay.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = self._merge_configs(result[key], value)
            else:
                result[key] = value
        
        return result
    
    def _interactive_basic_config(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Interactive basic configuration setup"""
        print("\nBasic Configuration Setup:")
        print("Press Enter to use default values shown in brackets")
        
        # System settings
        print(f"\n{Colors.BOLD}System Settings:{Colors.RESET}")
        
        max_workflows = input(f"Max concurrent workflows [{config['system']['max_concurrent_workflows']}]: ").strip()
        if max_workflows:
            config['system']['max_concurrent_workflows'] = int(max_workflows)
        
        # Resource settings
        print(f"\n{Colors.BOLD}Resource Settings:{Colors.RESET}")
        
        cpu_cores = input(f"Default CPU cores [{config['resources']['default_cpu_cores']}]: ").strip()
        if cpu_cores:
            config['resources']['default_cpu_cores'] = int(cpu_cores)
        
        memory_gb = input(f"Default memory (GB) [{config['resources']['default_memory_gb']}]: ").strip()
        if memory_gb:
            config['resources']['default_memory_gb'] = float(memory_gb)
        
        return config
    
    def _interactive_advanced_config(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Interactive advanced configuration setup"""
        # Start with basic config
        config = self._interactive_basic_config(config)
        
        # Advanced settings
        print(f"\n{Colors.BOLD}Advanced Settings:{Colors.RESET}")
        
        auto_submit = input("Auto-submit workflows? [y/N]: ").strip().lower()
        config['workflow']['auto_submit'] = auto_submit in ('y', 'yes')
        
        alert_email = input("Alert email address [none]: ").strip()
        if alert_email:
            config['monitoring']['alert_email'] = alert_email
        
        return config
    
    def _interactive_cluster_config(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Interactive cluster configuration setup"""
        # Start with advanced config
        config = self._interactive_advanced_config(config)
        
        # Cluster settings
        print(f"\n{Colors.BOLD}Cluster Settings:{Colors.RESET}")
        
        load_strategy = input("Load balancing strategy [resource_aware]: ").strip()
        if load_strategy:
            config['cluster']['load_balancing_strategy'] = load_strategy
        
        auto_scaling = input("Enable auto-scaling? [Y/n]: ").strip().lower()
        config['cluster']['auto_scaling'] = auto_scaling not in ('n', 'no')
        
        return config


# Import Colors for use in this module
from ..utils import Colors