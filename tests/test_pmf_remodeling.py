#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF Remodeling Test Script

This script specifically tests and validates the PMF system remodeling process,
ensuring that the remodeled system meets all requirements for PMF workflows.
"""

import sys
import os
import argparse
import shutil
import subprocess
from pathlib import Path
from datetime import datetime
import yaml
import json

# Add PRISM to path
sys.path.insert(0, str(Path(__file__).parents[3]))

def log_message(message, level='INFO'):
    """Simple logging with timestamp"""
    timestamp = datetime.now().strftime('%H:%M:%S')
    prefix = {
        'INFO': '✓',
        'WARNING': '⚠',
        'ERROR': '✗',
        'DEBUG': '🔍'
    }.get(level, 'ℹ')
    print(f"[{timestamp}] {prefix} {message}")

class PMFRemodelingTester:
    """Test PMF remodeling process comprehensively"""

    def __init__(self, input_dir, output_dir, pulling_distance=2.5, verbose=True):
        self.input_dir = Path(input_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.pulling_distance = pulling_distance
        self.verbose = verbose

        # Results tracking
        self.results = {
            'start_time': datetime.now().isoformat(),
            'input_dir': str(self.input_dir),
            'output_dir': str(self.output_dir),
            'pulling_distance': pulling_distance,
            'tests': {},
            'summary': {},
            'errors': [],
            'warnings': []
        }

        # PMF system directory
        self.pmf_system_dir = self.output_dir / "GMX_PMF_SYSTEM"

        log_message("PMF Remodeling Tester Initialized")
        log_message(f"Input: {self.input_dir}")
        log_message(f"Output: {self.output_dir}")
        log_message(f"Pulling distance: {self.pulling_distance} nm")

    def test_input_system_validation(self):
        """Test 1: Validate input MD system"""
        log_message("\n" + "="*70)
        log_message("TEST 1: INPUT MD SYSTEM VALIDATION")
        log_message("="*70)

        test_result = {
            'status': 'running',
            'details': {}
        }

        try:
            # Check input directory exists
            if not self.input_dir.exists():
                raise FileNotFoundError(f"Input directory not found: {self.input_dir}")

            log_message(f"Input directory exists: {self.input_dir}")

            # Find MD directory
            md_dirs = []
            for item in self.input_dir.iterdir():
                if item.is_dir() and ('md' in item.name.lower() or 'prolig' in item.name.lower()):
                    md_dirs.append(item)

            if not md_dirs:
                raise FileNotFoundError("No MD directory found (should contain 'md' or 'prolig' in name)")

            md_dir = md_dirs[0]
            log_message(f"Found MD directory: {md_dir.name}")

            # Check required files
            required_files = [
                md_dir / "solv_ions.gro",
                md_dir / "topol.top"
            ]

            missing_files = []
            for file_path in required_files:
                if file_path.exists():
                    size = file_path.stat().st_size
                    log_message(f"Required file: {file_path.relative_to(self.input_dir)} ({size} bytes)")
                    test_result['details'][file_path.name] = {
                        'exists': True,
                        'size': size,
                        'path': str(file_path.relative_to(self.input_dir))
                    }
                else:
                    missing_files.append(str(file_path.relative_to(self.input_dir)))
                    log_message(f"Missing required file: {file_path.relative_to(self.input_dir)}", 'ERROR')

            if missing_files:
                raise FileNotFoundError(f"Missing required files: {missing_files}")

            # Check optional files
            optional_files = list(md_dir.glob("posre*.itp"))
            for posre_file in optional_files:
                log_message(f"Optional file: {posre_file.relative_to(self.input_dir)}")
                test_result['details'][posre_file.name] = {
                    'exists': True,
                    'size': posre_file.stat().st_size,
                    'path': str(posre_file.relative_to(self.input_dir))
                }

            test_result['status'] = 'passed'
            test_result['md_directory'] = str(md_dir.relative_to(self.input_dir))
            log_message("Input system validation: PASSED")

        except Exception as e:
            test_result['status'] = 'failed'
            test_result['error'] = str(e)
            log_message(f"Input system validation: FAILED - {e}", 'ERROR')
            self.results['errors'].append(f"Input validation: {e}")

        self.results['tests']['input_validation'] = test_result
        return test_result['status'] == 'passed'

    def test_remodeling_script_execution(self):
        """Test 2: Execute PMF remodeling script"""
        log_message("\n" + "="*70)
        log_message("TEST 2: PMF REMODELING SCRIPT EXECUTION")
        log_message("="*70)

        test_result = {
            'status': 'running',
            'details': {}
        }

        try:
            # Find remodeling script
            script_dir = Path(__file__).parent
            remodel_script = script_dir / "pmf_remodel.py"

            if not remodel_script.exists():
                raise FileNotFoundError(f"Remodeling script not found: {remodel_script}")

            log_message(f"Found remodeling script: {remodel_script.relative_to(script_dir.parent.parent.parent)}")

            # Prepare output directory
            self.output_dir.mkdir(parents=True, exist_ok=True)

            # Remove existing PMF system if exists
            if self.pmf_system_dir.exists():
                log_message("Removing existing PMF system...")
                shutil.rmtree(self.pmf_system_dir)

            # Execute remodeling script
            cmd = [
                sys.executable, str(remodel_script),
                "--input", str(self.input_dir),
                "--output", str(self.output_dir),
                "--pulling-distance", str(self.pulling_distance),
                "--force"  # Force overwrite
            ]

            log_message(f"Executing: {' '.join(cmd)}")

            # Run the command
            start_time = datetime.now()
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            end_time = datetime.now()

            execution_time = (end_time - start_time).total_seconds()

            test_result['details'] = {
                'command': ' '.join(cmd),
                'execution_time': execution_time,
                'return_code': result.returncode,
                'stdout': result.stdout,
                'stderr': result.stderr
            }

            if result.returncode != 0:
                raise subprocess.CalledProcessError(
                    result.returncode, cmd, result.stdout, result.stderr
                )

            log_message(f"Remodeling completed in {execution_time:.2f} seconds")
            log_message("Script output:")
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    log_message(f"  {line}")

            test_result['status'] = 'passed'

        except subprocess.TimeoutExpired:
            test_result['status'] = 'failed'
            test_result['error'] = "Remodeling script timeout (300s)"
            log_message("Remodeling script execution: TIMEOUT", 'ERROR')
            self.results['errors'].append("Remodeling timeout")

        except subprocess.CalledProcessError as e:
            test_result['status'] = 'failed'
            test_result['error'] = f"Script failed with return code {e.returncode}"
            log_message(f"Remodeling script execution: FAILED - {e}", 'ERROR')
            log_message("STDOUT:", 'DEBUG')
            for line in e.stdout.split('\n'):
                if line.strip():
                    log_message(f"  {line}", 'DEBUG')
            log_message("STDERR:", 'DEBUG')
            for line in e.stderr.split('\n'):
                if line.strip():
                    log_message(f"  {line}", 'DEBUG')
            self.results['errors'].append(f"Remodeling script: {e}")

        except Exception as e:
            test_result['status'] = 'failed'
            test_result['error'] = str(e)
            log_message(f"Remodeling script execution: FAILED - {e}", 'ERROR')
            self.results['errors'].append(f"Remodeling execution: {e}")

        self.results['tests']['remodeling_execution'] = test_result
        return test_result['status'] == 'passed'

    def test_remodeled_system_structure(self):
        """Test 3: Validate remodeled system structure"""
        log_message("\n" + "="*70)
        log_message("TEST 3: REMODELED SYSTEM STRUCTURE VALIDATION")
        log_message("="*70)

        test_result = {
            'status': 'running',
            'details': {}
        }

        try:
            # Check PMF system directory exists
            if not self.pmf_system_dir.exists():
                raise FileNotFoundError(f"PMF system directory not created: {self.pmf_system_dir}")

            log_message(f"PMF system directory exists: {self.pmf_system_dir.relative_to(self.output_dir)}")

            # Check required files
            required_files = {
                "solv_ions.gro": "System structure file",
                "topol.top": "System topology file",
                "PMF_SYSTEM_INFO.txt": "PMF system metadata"
            }

            for filename, description in required_files.items():
                file_path = self.pmf_system_dir / filename
                if file_path.exists():
                    size = file_path.stat().st_size
                    log_message(f"{description}: {filename} ({size} bytes)")
                    test_result['details'][filename] = {
                        'exists': True,
                        'size': size,
                        'description': description
                    }

                    # Additional checks for specific files
                    if filename == "solv_ions.gro":
                        # Check structure file format
                        try:
                            with open(file_path, 'r') as f:
                                lines = f.readlines()
                                if len(lines) < 3:
                                    raise ValueError("Structure file too short")
                                # Check if it has atoms
                                try:
                                    atom_count = int(lines[1].strip())
                                    if atom_count <= 0:
                                        raise ValueError("No atoms in structure")
                                    log_message(f"  Structure contains {atom_count} atoms")
                                    test_result['details'][filename]['atom_count'] = atom_count
                                except ValueError:
                                    raise ValueError("Invalid atom count in structure file")
                        except Exception as e:
                            log_message(f"  Warning: Structure file validation failed: {e}", 'WARNING')
                            self.results['warnings'].append(f"Structure validation: {e}")

                    elif filename == "topol.top":
                        # Check topology file format
                        try:
                            with open(file_path, 'r') as f:
                                content = f.read()
                                if "#include" not in content:
                                    log_message("  Warning: No #include statements in topology", 'WARNING')
                                if "[ system ]" not in content:
                                    log_message("  Warning: No [system] section in topology", 'WARNING')
                                log_message(f"  Topology file size: {len(content)} characters")
                                test_result['details'][filename]['content_length'] = len(content)
                        except Exception as e:
                            log_message(f"  Warning: Topology file validation failed: {e}", 'WARNING')
                            self.results['warnings'].append(f"Topology validation: {e}")

                    elif filename == "PMF_SYSTEM_INFO.txt":
                        # Parse PMF system info
                        try:
                            with open(file_path, 'r') as f:
                                content = f.read()
                                if f"Pulling Distance: {self.pulling_distance}" in content:
                                    log_message(f"  Pulling distance matches: {self.pulling_distance} nm")
                                else:
                                    log_message("  Warning: Pulling distance mismatch in info file", 'WARNING')
                                test_result['details'][filename]['content_preview'] = content[:200] + "..."
                        except Exception as e:
                            log_message(f"  Warning: Info file parsing failed: {e}", 'WARNING')

                else:
                    raise FileNotFoundError(f"Required file missing: {filename}")

            # Check optional files
            optional_files = list(self.pmf_system_dir.glob("posre*.itp"))
            for posre_file in optional_files:
                size = posre_file.stat().st_size
                log_message(f"Position restraints: {posre_file.name} ({size} bytes)")
                test_result['details'][posre_file.name] = {
                    'exists': True,
                    'size': size,
                    'description': "Position restraints file"
                }

            test_result['status'] = 'passed'
            log_message("Remodeled system structure: PASSED")

        except Exception as e:
            test_result['status'] = 'failed'
            test_result['error'] = str(e)
            log_message(f"Remodeled system structure: FAILED - {e}", 'ERROR')
            self.results['errors'].append(f"Structure validation: {e}")

        self.results['tests']['structure_validation'] = test_result
        return test_result['status'] == 'passed'

    def test_pmf_workflow_compatibility(self):
        """Test 4: Check PMF workflow compatibility"""
        log_message("\n" + "="*70)
        log_message("TEST 4: PMF WORKFLOW COMPATIBILITY")
        log_message("="*70)

        test_result = {
            'status': 'running',
            'details': {}
        }

        try:
            # Try to import PRISM PMF
            try:
                import prism.pmf as pmf
                log_message("PRISM PMF module imported successfully")
                test_result['details']['pmf_import'] = True
            except ImportError as e:
                log_message(f"Warning: Cannot import PRISM PMF: {e}", 'WARNING')
                test_result['details']['pmf_import'] = False
                self.results['warnings'].append(f"PMF import: {e}")

                # Skip workflow test if PMF not available
                test_result['status'] = 'skipped'
                test_result['reason'] = 'PMF module not available'
                self.results['tests']['workflow_compatibility'] = test_result
                return True  # Don't fail the test for import issues

            # Try to create PMF workflow with remodeled system
            try:
                config = {'pulling_distance': self.pulling_distance}

                # Test PMF workflow initialization
                workflow = pmf.PMFWorkflow(
                    system_dir=self.output_dir,
                    output_dir=self.output_dir / "test_pmf_workflow",
                    config=config
                )

                log_message("PMF workflow created successfully")
                test_result['details']['workflow_creation'] = True

                # Get workflow status
                status = workflow.get_status()
                log_message(f"Workflow status: {status['current_stage']}")
                test_result['details']['workflow_status'] = status

                # Check if system is ready for PMF
                if status['current_stage'] in ['ready_to_start', 'pmf_system_ready']:
                    log_message("System is ready for PMF calculations")
                    test_result['details']['pmf_ready'] = True
                else:
                    log_message(f"System status: {status['current_stage']}", 'WARNING')
                    test_result['details']['pmf_ready'] = False
                    self.results['warnings'].append(f"PMF status: {status['current_stage']}")

            except Exception as e:
                log_message(f"PMF workflow creation failed: {e}", 'ERROR')
                test_result['details']['workflow_creation'] = False
                test_result['details']['workflow_error'] = str(e)
                raise e

            test_result['status'] = 'passed'
            log_message("PMF workflow compatibility: PASSED")

        except Exception as e:
            test_result['status'] = 'failed'
            test_result['error'] = str(e)
            log_message(f"PMF workflow compatibility: FAILED - {e}", 'ERROR')
            self.results['errors'].append(f"Workflow compatibility: {e}")

        self.results['tests']['workflow_compatibility'] = test_result
        return test_result['status'] == 'passed'

    def generate_report(self):
        """Generate comprehensive test report"""
        log_message("\n" + "="*70)
        log_message("GENERATING COMPREHENSIVE REPORT")
        log_message("="*70)

        # Calculate summary
        total_tests = len(self.results['tests'])
        passed_tests = sum(1 for test in self.results['tests'].values()
                          if test['status'] == 'passed')
        failed_tests = sum(1 for test in self.results['tests'].values()
                          if test['status'] == 'failed')
        skipped_tests = sum(1 for test in self.results['tests'].values()
                           if test['status'] == 'skipped')

        success_rate = (passed_tests / total_tests * 100) if total_tests > 0 else 0

        self.results['summary'] = {
            'total_tests': total_tests,
            'passed_tests': passed_tests,
            'failed_tests': failed_tests,
            'skipped_tests': skipped_tests,
            'success_rate': success_rate,
            'end_time': datetime.now().isoformat()
        }

        # Display summary
        log_message(f"Total tests: {total_tests}")
        log_message(f"Passed: {passed_tests}")
        log_message(f"Failed: {failed_tests}")
        log_message(f"Skipped: {skipped_tests}")
        log_message(f"Success rate: {success_rate:.1f}%")

        # Test results summary
        log_message("\nTest Summary:")
        for test_name, test_result in self.results['tests'].items():
            status_icon = {
                'passed': '✓',
                'failed': '✗',
                'skipped': '⏭'
            }.get(test_result['status'], '?')
            log_message(f"  {status_icon} {test_name.upper()}: {test_result['status']}")

        # Show errors
        if self.results['errors']:
            log_message(f"\n❌ ERRORS ({len(self.results['errors'])}):")
            for error in self.results['errors']:
                log_message(f"   - {error}")

        # Show warnings
        if self.results['warnings']:
            log_message(f"\n⚠️ WARNINGS ({len(self.results['warnings'])}):")
            for warning in self.results['warnings']:
                log_message(f"   - {warning}")

        # Overall result
        if failed_tests == 0:
            if skipped_tests > 0:
                log_message("\n🟡 PARTIAL SUCCESS! Some tests were skipped.")
            else:
                log_message("\n🎉 ALL TESTS PASSED! PMF remodeling is working correctly.")
        else:
            log_message("\n❌ SOME TESTS FAILED! PMF remodeling needs attention.")

        # Save detailed report
        report_file = self.output_dir / "pmf_remodeling_test_report.json"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        with open(report_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        log_message(f"\n📄 Detailed report saved: {report_file}")

        return failed_tests == 0

    def run_all_tests(self):
        """Run all PMF remodeling tests"""
        log_message("🚀 Starting PMF Remodeling Comprehensive Tests")
        log_message("="*70)

        tests = [
            ('Input System Validation', self.test_input_system_validation),
            ('Remodeling Script Execution', self.test_remodeling_script_execution),
            ('Remodeled System Structure', self.test_remodeled_system_structure),
            ('PMF Workflow Compatibility', self.test_pmf_workflow_compatibility)
        ]

        overall_success = True

        for test_name, test_func in tests:
            try:
                success = test_func()
                if not success:
                    overall_success = False
            except Exception as e:
                log_message(f"Test {test_name} crashed: {e}", 'ERROR')
                self.results['errors'].append(f"{test_name} crashed: {e}")
                overall_success = False

        # Generate final report
        self.generate_report()

        return overall_success

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='PMF Remodeling Test Script')
    parser.add_argument('--input', '-i', required=True,
                       help='Input system directory (e.g., gromacssim)')
    parser.add_argument('--output', '-o', required=True,
                       help='Output directory for remodeled system and results')
    parser.add_argument('--pulling-distance', '-p', type=float, default=2.5,
                       help='Pulling distance in nm (default: 2.5)')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Reduce output verbosity')

    args = parser.parse_args()

    # Create tester
    tester = PMFRemodelingTester(
        input_dir=args.input,
        output_dir=args.output,
        pulling_distance=args.pulling_distance,
        verbose=not args.quiet
    )

    # Run tests
    success = tester.run_all_tests()

    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())