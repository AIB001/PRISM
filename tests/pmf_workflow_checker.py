#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF Workflow Comprehensive Checker

This script provides a complete validation and testing framework for the PRISM PMF module.
It checks all components, APIs, configurations, and provides debugging information.

Usage:
    python pmf_workflow_checker.py
    python pmf_workflow_checker.py --test-system ./gaff_model
    python pmf_workflow_checker.py --verbose --check-imports
"""

import sys
import os
import traceback
from pathlib import Path
import argparse

# Add PRISM to path
sys.path.insert(0, str(Path(__file__).parents[3]))

class PMFWorkflowChecker:
    """Comprehensive PMF workflow validation and testing"""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.results = {
            'import_check': False,
            'api_check': False,
            'config_check': False,
            'examples_check': False,
            'documentation_check': False,
            'system_check': False,
            'errors': [],
            'warnings': []
        }

    def log(self, message, level='INFO'):
        """Log message with level"""
        if self.verbose or level in ['ERROR', 'WARNING']:
            print(f"[{level}] {message}")

    def check_imports(self):
        """Check all PMF module imports"""
        self.log("="*60)
        self.log("1. CHECKING PMF MODULE IMPORTS")
        self.log("="*60)

        try:
            # Test basic import
            import prism.pmf as pmf
            self.log("✓ Basic prism.pmf import successful")

            # Test core components
            core_imports = [
                ('PMFSystem', pmf.PMFSystem),
                ('pmf_system', pmf.pmf_system),
                ('PMFRunner', pmf.PMFRunner),
                ('run_pmf_workflow', pmf.run_pmf_workflow),
                ('create_pmf_config', pmf.create_pmf_config),
                ('PMFBuilder', pmf.PMFBuilder),
                ('SMDManager', pmf.SMDManager),
                ('UmbrellaManager', pmf.UmbrellaManager),
                ('PMFAnalyzer', pmf.PMFAnalyzer),
            ]

            for name, obj in core_imports:
                if obj is not None:
                    self.log(f"✓ {name} import successful")
                else:
                    self.log(f"✗ {name} import failed", 'ERROR')
                    self.results['errors'].append(f"Failed to import {name}")

            # Test module info
            info = pmf.get_pmf_info()
            self.log(f"✓ PMF module version: {info['version']}")
            self.log(f"✓ Architecture: {info['architecture']}")

            # Test optional components
            try:
                optimizer = pmf.PMFWorkflowOptimizer
                self.log("✓ Workflow optimizer available")
            except AttributeError:
                self.log("⚠ Workflow optimizer not available", 'WARNING')
                self.results['warnings'].append("Workflow optimizer module missing")

            self.results['import_check'] = True
            self.log("✓ All core imports successful")

        except Exception as e:
            self.log(f"✗ Import check failed: {e}", 'ERROR')
            self.results['errors'].append(f"Import check failed: {e}")
            if self.verbose:
                traceback.print_exc()

    def check_api_interfaces(self):
        """Check API interfaces and basic functionality"""
        self.log("\\n" + "="*60)
        self.log("2. CHECKING API INTERFACES")
        self.log("="*60)

        try:
            import prism.pmf as pmf

            # Test configuration creation
            self.log("Testing configuration creation...")
            test_config_path = "/tmp/test_pmf_config.yaml"

            templates = ['default', 'fast', 'accurate']
            for template in templates:
                try:
                    pmf.create_pmf_config(test_config_path, template=template)
                    if Path(test_config_path).exists():
                        self.log(f"✓ {template} configuration template created")
                        os.remove(test_config_path)
                    else:
                        self.log(f"✗ {template} configuration not created", 'ERROR')
                except Exception as e:
                    self.log(f"✗ {template} config creation failed: {e}", 'ERROR')

            # Test PMFSystem initialization
            self.log("Testing PMFSystem initialization...")
            try:
                system = pmf.pmf_system(
                    system_dir="./test_system",
                    output_dir="./test_output",
                    pulling_distance=2.0
                )
                self.log("✓ PMFSystem initialization successful")

                # Test status check (should work even without actual system)
                try:
                    status = system.get_status()
                    self.log(f"✓ Status check successful: {status.get('current_stage', 'unknown')}")
                except Exception as e:
                    self.log(f"⚠ Status check warning: {e}", 'WARNING')

            except Exception as e:
                self.log(f"✗ PMFSystem initialization failed: {e}", 'ERROR')
                self.results['errors'].append(f"PMFSystem init failed: {e}")

            # Test PMFRunner initialization
            self.log("Testing PMFRunner initialization...")
            try:
                runner = pmf.PMFRunner()
                self.log("✓ PMFRunner initialization successful")
            except Exception as e:
                self.log(f"✗ PMFRunner initialization failed: {e}", 'ERROR')

            # Test individual managers
            managers = [
                ('SMDManager', pmf.SMDManager),
                ('UmbrellaManager', pmf.UmbrellaManager),
                ('PMFAnalyzer', pmf.PMFAnalyzer),
                ('PMFBuilder', pmf.PMFBuilder)
            ]

            for name, manager_class in managers:
                try:
                    manager = manager_class()
                    self.log(f"✓ {name} initialization successful")
                except Exception as e:
                    self.log(f"✗ {name} initialization failed: {e}", 'ERROR')

            self.results['api_check'] = True
            self.log("✓ API interface check completed")

        except Exception as e:
            self.log(f"✗ API check failed: {e}", 'ERROR')
            self.results['errors'].append(f"API check failed: {e}")

    def check_configuration_system(self):
        """Check configuration file system"""
        self.log("\\n" + "="*60)
        self.log("3. CHECKING CONFIGURATION SYSTEM")
        self.log("="*60)

        try:
            import prism.pmf as pmf
            import yaml

            # Test config creation and loading
            test_dir = Path("/tmp/pmf_config_test")
            test_dir.mkdir(exist_ok=True)

            templates = ['default', 'fast', 'accurate']
            for template in templates:
                config_file = test_dir / f"{template}_config.yaml"

                try:
                    pmf.create_pmf_config(str(config_file), template=template)

                    if config_file.exists():
                        # Try to load and validate config
                        with open(config_file, 'r') as f:
                            config_data = yaml.safe_load(f)

                        # Check required sections
                        required_sections = ['builder', 'smd', 'umbrella', 'analysis']
                        missing_sections = []

                        for section in required_sections:
                            if section not in config_data:
                                missing_sections.append(section)

                        if missing_sections:
                            self.log(f"⚠ {template} config missing sections: {missing_sections}", 'WARNING')
                        else:
                            self.log(f"✓ {template} configuration valid")

                        # Check key parameters
                        if 'smd' in config_data:
                            smd_params = ['pull_rate', 'pull_k', 'nsteps']
                            for param in smd_params:
                                if param in config_data['smd']:
                                    self.log(f"  ✓ SMD parameter {param}: {config_data['smd'][param]}")

                        config_file.unlink()  # Clean up

                except Exception as e:
                    self.log(f"✗ {template} config test failed: {e}", 'ERROR')

            # Clean up test directory
            test_dir.rmdir()

            self.results['config_check'] = True
            self.log("✓ Configuration system check completed")

        except Exception as e:
            self.log(f"✗ Configuration check failed: {e}", 'ERROR')
            self.results['errors'].append(f"Configuration check failed: {e}")

    def check_examples(self):
        """Check example files"""
        self.log("\\n" + "="*60)
        self.log("4. CHECKING EXAMPLE FILES")
        self.log("="*60)

        try:
            examples_dir = Path(__file__).parent

            example_files = [
                'PMF_tutorial.ipynb',
                'smd_user_custom_example.py',
                '../api_examples.py'
            ]

            for example_file in example_files:
                example_path = examples_dir / example_file

                if example_path.exists():
                    self.log(f"✓ Example found: {example_file}")

                    # Check if Python file can be imported (syntax check)
                    if example_file.endswith('.py'):
                        try:
                            with open(example_path, 'r') as f:
                                content = f.read()

                            # Basic syntax check
                            compile(content, str(example_path), 'exec')
                            self.log(f"  ✓ Syntax valid")
                        except SyntaxError as e:
                            self.log(f"  ✗ Syntax error: {e}", 'ERROR')
                        except Exception as e:
                            self.log(f"  ⚠ Could not validate: {e}", 'WARNING')

                    elif example_file.endswith('.ipynb'):
                        # Check notebook structure
                        try:
                            import json
                            with open(example_path, 'r') as f:
                                notebook = json.load(f)

                            if 'cells' in notebook:
                                cell_count = len(notebook['cells'])
                                self.log(f"  ✓ Notebook with {cell_count} cells")
                            else:
                                self.log(f"  ⚠ Invalid notebook structure", 'WARNING')
                        except Exception as e:
                            self.log(f"  ⚠ Could not validate notebook: {e}", 'WARNING')

                else:
                    self.log(f"✗ Example missing: {example_file}", 'ERROR')
                    self.results['errors'].append(f"Missing example: {example_file}")

            self.results['examples_check'] = True
            self.log("✓ Example files check completed")

        except Exception as e:
            self.log(f"✗ Examples check failed: {e}", 'ERROR')
            self.results['errors'].append(f"Examples check failed: {e}")

    def check_documentation(self):
        """Check documentation completeness"""
        self.log("\\n" + "="*60)
        self.log("5. CHECKING DOCUMENTATION")
        self.log("="*60)

        try:
            pmf_dir = Path(__file__).parent.parent

            # Check README
            readme_path = pmf_dir / "README.md"
            if readme_path.exists():
                with open(readme_path, 'r') as f:
                    readme_content = f.read()

                # Check for key sections
                required_sections = [
                    'Quick Start',
                    'Configuration',
                    'API Reference',
                    'Troubleshooting'
                ]

                for section in required_sections:
                    if section.lower() in readme_content.lower():
                        self.log(f"✓ README section found: {section}")
                    else:
                        self.log(f"⚠ README section missing: {section}", 'WARNING')

                self.log(f"✓ README.md exists ({len(readme_content)} characters)")
            else:
                self.log("✗ README.md missing", 'ERROR')
                self.results['errors'].append("README.md missing")

            # Check docstrings in key modules
            key_modules = [
                '__init__.py',
                'runner.py',
                'core/system.py',
                'methods/smd.py'
            ]

            for module in key_modules:
                module_path = pmf_dir / module
                if module_path.exists():
                    with open(module_path, 'r') as f:
                        content = f.read()

                    if '"""' in content:
                        self.log(f"✓ {module} has docstrings")
                    else:
                        self.log(f"⚠ {module} missing docstrings", 'WARNING')
                else:
                    self.log(f"✗ {module} not found", 'ERROR')

            self.results['documentation_check'] = True
            self.log("✓ Documentation check completed")

        except Exception as e:
            self.log(f"✗ Documentation check failed: {e}", 'ERROR')
            self.results['errors'].append(f"Documentation check failed: {e}")

    def check_test_system(self, system_path=None):
        """Check with a real test system if provided"""
        self.log("\\n" + "="*60)
        self.log("6. CHECKING TEST SYSTEM INTEGRATION")
        self.log("="*60)

        if system_path is None:
            self.log("⚠ No test system provided, skipping integration test", 'WARNING')
            self.results['system_check'] = True
            return

        try:
            import prism.pmf as pmf

            system_path = Path(system_path)
            if not system_path.exists():
                self.log(f"✗ Test system not found: {system_path}", 'ERROR')
                self.results['errors'].append(f"Test system not found: {system_path}")
                return

            self.log(f"Testing with system: {system_path}")

            # Check required files
            required_files = [
                'GMX_PROLIG_MD/solv_ions.gro',
                'GMX_PROLIG_MD/topol.top'
            ]

            for file_path in required_files:
                full_path = system_path / file_path
                if full_path.exists():
                    self.log(f"✓ Required file found: {file_path}")
                else:
                    self.log(f"✗ Required file missing: {file_path}", 'ERROR')
                    self.results['errors'].append(f"Missing file: {file_path}")

            # Test PMF system creation
            try:
                test_output = "/tmp/pmf_test_output"
                system = pmf.pmf_system(
                    system_dir=str(system_path),
                    output_dir=test_output,
                    pulling_distance=2.0
                )

                self.log("✓ PMF system created successfully")

                # Test status check
                status = system.get_status()
                self.log(f"✓ System status: {status.get('current_stage', 'unknown')}")

                # Test validation
                try:
                    from prism.pmf.analysis import validate_md_results
                    validation_result = validate_md_results(str(system_path))
                    self.log(f"✓ System validation: {validation_result.get('status', 'unknown')}")
                except Exception as e:
                    self.log(f"⚠ Validation warning: {e}", 'WARNING')

                # Clean up
                import shutil
                if Path(test_output).exists():
                    shutil.rmtree(test_output)

            except Exception as e:
                self.log(f"✗ PMF system test failed: {e}", 'ERROR')
                self.results['errors'].append(f"System test failed: {e}")

            self.results['system_check'] = True
            self.log("✓ Test system integration completed")

        except Exception as e:
            self.log(f"✗ System check failed: {e}", 'ERROR')
            self.results['errors'].append(f"System check failed: {e}")

    def generate_report(self):
        """Generate comprehensive test report"""
        self.log("\\n" + "="*60)
        self.log("PMF WORKFLOW CHECK SUMMARY")
        self.log("="*60)

        total_checks = len([k for k in self.results.keys() if k.endswith('_check')])
        passed_checks = sum([v for k, v in self.results.items() if k.endswith('_check')])

        self.log(f"Total checks: {total_checks}")
        self.log(f"Passed checks: {passed_checks}")
        self.log(f"Success rate: {passed_checks/total_checks*100:.1f}%")

        check_status = {
            'import_check': 'Module Imports',
            'api_check': 'API Interfaces',
            'config_check': 'Configuration System',
            'examples_check': 'Example Files',
            'documentation_check': 'Documentation',
            'system_check': 'System Integration'
        }

        self.log("\\nDetailed Results:")
        for check, description in check_status.items():
            status = "✓ PASS" if self.results[check] else "✗ FAIL"
            self.log(f"  {description}: {status}")

        if self.results['errors']:
            self.log(f"\\n❌ ERRORS ({len(self.results['errors'])}):")
            for error in self.results['errors']:
                self.log(f"  - {error}")

        if self.results['warnings']:
            self.log(f"\\n⚠️  WARNINGS ({len(self.results['warnings'])}):")
            for warning in self.results['warnings']:
                self.log(f"  - {warning}")

        if passed_checks == total_checks:
            self.log("\\n🎉 ALL CHECKS PASSED! PMF workflow is ready for use.")
        elif passed_checks >= total_checks * 0.8:
            self.log("\\n✅ MOSTLY WORKING! Some issues found but core functionality intact.")
        else:
            self.log("\\n❌ SIGNIFICANT ISSUES! Please address errors before using PMF workflow.")

        return self.results

    def run_complete_check(self, system_path=None):
        """Run all checks"""
        self.log("Starting comprehensive PMF workflow check...")
        self.log(f"Working directory: {Path.cwd()}")
        self.log(f"Script location: {Path(__file__).parent}")

        self.check_imports()
        self.check_api_interfaces()
        self.check_configuration_system()
        self.check_examples()
        self.check_documentation()
        self.check_test_system(system_path)

        return self.generate_report()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='PMF Workflow Comprehensive Checker')
    parser.add_argument('--test-system', type=str, help='Path to test MD system')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--check-imports', action='store_true', help='Only check imports')

    args = parser.parse_args()

    checker = PMFWorkflowChecker(verbose=args.verbose)

    if args.check_imports:
        checker.check_imports()
        checker.generate_report()
    else:
        results = checker.run_complete_check(args.test_system)

        # Exit with appropriate code
        if results['errors']:
            sys.exit(1)
        elif results['warnings']:
            sys.exit(2)
        else:
            sys.exit(0)


if __name__ == "__main__":
    main()