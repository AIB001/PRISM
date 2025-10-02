#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Core Testing Framework

Provides a comprehensive testing framework specifically designed for PRISM
molecular dynamics and PMF calculations with scientific computing considerations.
"""

import os
import sys
import time
import unittest
import traceback
from pathlib import Path
from typing import Dict, Any, List, Optional, Callable, Union
from dataclasses import dataclass, asdict
from enum import Enum
from contextlib import contextmanager
import json
import tempfile
import shutil

from ..utils.logging_system import PrismLogger, LogLevel, EventType


class TestType(Enum):
    """Types of tests in PRISM testing framework"""
    UNIT = "unit"
    INTEGRATION = "integration"
    PERFORMANCE = "performance"
    REGRESSION = "regression"
    SYSTEM = "system"
    ACCEPTANCE = "acceptance"


class TestStatus(Enum):
    """Test execution status"""
    PENDING = "pending"
    RUNNING = "running"
    PASSED = "passed"
    FAILED = "failed"
    SKIPPED = "skipped"
    ERROR = "error"


@dataclass
class TestResult:
    """Comprehensive test result structure"""
    test_name: str
    test_type: TestType
    status: TestStatus
    execution_time: float
    memory_usage: Optional[float] = None
    assertions_count: int = 0
    error_message: Optional[str] = None
    error_traceback: Optional[str] = None
    output_data: Dict[str, Any] = None
    performance_metrics: Dict[str, float] = None
    
    def __post_init__(self):
        if self.output_data is None:
            self.output_data = {}
        if self.performance_metrics is None:
            self.performance_metrics = {}
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        data = asdict(self)
        data['test_type'] = self.test_type.value
        data['status'] = self.status.value
        return data
    
    def is_successful(self) -> bool:
        """Check if test was successful"""
        return self.status in [TestStatus.PASSED, TestStatus.SKIPPED]


class PrismTestCase(unittest.TestCase):
    """Enhanced test case for PRISM with scientific computing features"""
    
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)
        self.test_type = TestType.UNIT
        self.logger = PrismLogger(f"test.{self.__class__.__name__}")
        self.test_start_time = None
        self.test_temp_dir = None
        self.performance_metrics = {}
        self.assertions_count = 0
        
    def setUp(self):
        """Enhanced setup with logging and temp directory"""
        self.test_start_time = time.time()
        self.test_temp_dir = Path(tempfile.mkdtemp(prefix=f"prism_test_{self._testMethodName}_"))
        self.logger.info(f"Starting test: {self._testMethodName}")
        
    def tearDown(self):
        """Enhanced teardown with cleanup"""
        execution_time = time.time() - self.test_start_time if self.test_start_time else 0
        
        # Cleanup temp directory
        if self.test_temp_dir and self.test_temp_dir.exists():
            try:
                shutil.rmtree(self.test_temp_dir)
            except Exception as e:
                self.logger.warning(f"Failed to cleanup temp dir: {e}")
        
        self.logger.info(f"Completed test: {self._testMethodName} in {execution_time:.3f}s")
    
    def assertFloatEqual(self, first: float, second: float, tolerance: float = 1e-6,
                        msg: Optional[str] = None):
        """Assert floating point equality with tolerance"""
        self.assertions_count += 1
        if abs(first - second) > tolerance:
            standard_msg = f"{first} != {second} within tolerance {tolerance}"
            msg = self._formatMessage(msg, standard_msg)
            raise self.failureException(msg)
    
    def assertEnergyEqual(self, first: float, second: float, tolerance: float = 0.1,
                         msg: Optional[str] = None):
        """Assert energy values equality (kcal/mol typically)"""
        self.assertions_count += 1
        self.assertFloatEqual(first, second, tolerance, msg)
    
    def assertDistanceEqual(self, first: float, second: float, tolerance: float = 0.01,
                           msg: Optional[str] = None):
        """Assert distance values equality (nm typically)"""
        self.assertions_count += 1
        self.assertFloatEqual(first, second, tolerance, msg)
    
    def assertFileExists(self, file_path: Union[str, Path], msg: Optional[str] = None):
        """Assert that a file exists"""
        self.assertions_count += 1
        file_path = Path(file_path)
        if not file_path.exists():
            standard_msg = f"File does not exist: {file_path}"
            msg = self._formatMessage(msg, standard_msg)
            raise self.failureException(msg)
    
    def assertFileNotEmpty(self, file_path: Union[str, Path], msg: Optional[str] = None):
        """Assert that a file exists and is not empty"""
        self.assertFileExists(file_path, msg)
        self.assertions_count += 1
        file_path = Path(file_path)
        if file_path.stat().st_size == 0:
            standard_msg = f"File is empty: {file_path}"
            msg = self._formatMessage(msg, standard_msg)
            raise self.failureException(msg)
    
    def assertValidTrajectory(self, trajectory_path: Union[str, Path], 
                             topology_path: Union[str, Path],
                             min_frames: int = 1, msg: Optional[str] = None):
        """Assert that trajectory file is valid"""
        self.assertFileExists(trajectory_path, msg)
        self.assertFileExists(topology_path, msg)
        
        # Additional trajectory validation would go here
        # For now, just check file existence and size
        self.assertFileNotEmpty(trajectory_path, msg)
        self.assertFileNotEmpty(topology_path, msg)
    
    def assertValidPMFResult(self, pmf_result: Dict[str, Any], msg: Optional[str] = None):
        """Assert that PMF calculation result is valid"""
        self.assertions_count += 1
        required_keys = ['binding_energy', 'pmf_profile', 'error_estimate']
        
        for key in required_keys:
            if key not in pmf_result:
                standard_msg = f"PMF result missing required key: {key}"
                msg = self._formatMessage(msg, standard_msg)
                raise self.failureException(msg)
        
        # Validate binding energy is reasonable
        binding_energy = pmf_result['binding_energy']
        if not isinstance(binding_energy, (int, float)):
            standard_msg = f"Binding energy must be numeric, got {type(binding_energy)}"
            msg = self._formatMessage(msg, standard_msg)
            raise self.failureException(msg)
        
        # Reasonable range for binding energies (-50 to +10 kcal/mol)
        if not (-50 <= binding_energy <= 10):
            standard_msg = f"Binding energy {binding_energy} outside reasonable range (-50 to +10 kcal/mol)"
            msg = self._formatMessage(msg, standard_msg)
            raise self.failureException(msg)
    
    @contextmanager
    def assertExecutionTime(self, max_seconds: float):
        """Context manager to assert execution time"""
        start_time = time.time()
        yield
        execution_time = time.time() - start_time
        self.performance_metrics['execution_time'] = execution_time
        
        if execution_time > max_seconds:
            raise self.failureException(
                f"Execution took {execution_time:.3f}s, expected <= {max_seconds}s"
            )
    
    def create_test_file(self, filename: str, content: str = "") -> Path:
        """Create a test file in the temp directory"""
        if not self.test_temp_dir:
            raise RuntimeError("Test temp directory not available")
        
        file_path = self.test_temp_dir / filename
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(file_path, 'w') as f:
            f.write(content)
        
        return file_path
    
    def create_mock_system(self, system_name: str = "test_system") -> Path:
        """Create a mock molecular system for testing"""
        system_dir = self.test_temp_dir / system_name
        system_dir.mkdir(exist_ok=True)
        
        # Create mock files
        files_to_create = {
            'system.gro': self._mock_gro_content(),
            'system.top': self._mock_top_content(),
            'index.ndx': self._mock_ndx_content(),
            'production.mdp': self._mock_mdp_content()
        }
        
        for filename, content in files_to_create.items():
            with open(system_dir / filename, 'w') as f:
                f.write(content)
        
        return system_dir
    
    def _mock_gro_content(self) -> str:
        """Generate mock GRO file content"""
        return """Test system
    2
    1PRO     CA    1   1.000   1.000   1.000
    2LIG     C1    2   2.000   2.000   2.000
   3.00000   3.00000   3.00000
"""
    
    def _mock_top_content(self) -> str:
        """Generate mock TOP file content"""
        return """; Mock topology file
#include "forcefield.itp"

[ system ]
Test System

[ molecules ]
Protein 1
LIG     1
"""
    
    def _mock_ndx_content(self) -> str:
        """Generate mock index file content"""
        return """[ System ]
1 2

[ Protein ]
1

[ LIG ]
2
"""
    
    def _mock_mdp_content(self) -> str:
        """Generate mock MDP file content"""
        return """; Mock MDP file
integrator = md
dt = 0.002
nsteps = 1000
"""


class PrismTestSuite:
    """Test suite manager for PRISM tests"""
    
    def __init__(self, name: str):
        self.name = name
        self.test_cases: List[PrismTestCase] = []
        self.test_results: List[TestResult] = []
        self.logger = PrismLogger(f"test_suite.{name}")
        
    def add_test(self, test_case: Union[PrismTestCase, unittest.TestCase]):
        """Add a test case to the suite"""
        self.test_cases.append(test_case)
    
    def add_tests_from_module(self, module):
        """Add all test cases from a module"""
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromModule(module)
        
        for test_group in suite:
            if hasattr(test_group, '_tests'):
                for test in test_group._tests:
                    self.add_test(test)
    
    def run_tests(self, runner: Optional['PrismTestRunner'] = None) -> List[TestResult]:
        """Run all tests in the suite"""
        if runner is None:
            runner = PrismTestRunner()
        
        self.logger.info(f"Running test suite: {self.name} ({len(self.test_cases)} tests)")
        results = []
        
        for test_case in self.test_cases:
            result = runner.run_single_test(test_case)
            results.append(result)
            self.test_results.append(result)
        
        return results
    
    def get_summary(self) -> Dict[str, Any]:
        """Get test suite summary"""
        if not self.test_results:
            return {"message": "No tests run yet"}
        
        total = len(self.test_results)
        passed = sum(1 for r in self.test_results if r.status == TestStatus.PASSED)
        failed = sum(1 for r in self.test_results if r.status == TestStatus.FAILED)
        skipped = sum(1 for r in self.test_results if r.status == TestStatus.SKIPPED)
        errors = sum(1 for r in self.test_results if r.status == TestStatus.ERROR)
        
        total_time = sum(r.execution_time for r in self.test_results)
        
        return {
            'suite_name': self.name,
            'total_tests': total,
            'passed': passed,
            'failed': failed,
            'skipped': skipped,
            'errors': errors,
            'success_rate': (passed / total) * 100 if total > 0 else 0,
            'total_execution_time': total_time,
            'average_execution_time': total_time / total if total > 0 else 0
        }


class PrismTestRunner:
    """Enhanced test runner for PRISM tests"""
    
    def __init__(self, output_dir: Optional[Path] = None):
        self.output_dir = output_dir or Path("test_results")
        self.output_dir.mkdir(exist_ok=True)
        self.logger = PrismLogger("test_runner")
        
    def run_single_test(self, test_case: unittest.TestCase) -> TestResult:
        """Run a single test case"""
        test_name = f"{test_case.__class__.__name__}.{test_case._testMethodName}"
        start_time = time.time()
        
        self.logger.info(f"Running test: {test_name}")
        
        # Create test result
        result = TestResult(
            test_name=test_name,
            test_type=getattr(test_case, 'test_type', TestType.UNIT),
            status=TestStatus.RUNNING,
            execution_time=0.0
        )
        
        try:
            # Run the test
            test_case.setUp()
            test_method = getattr(test_case, test_case._testMethodName)
            test_method()
            test_case.tearDown()
            
            # Test passed
            result.status = TestStatus.PASSED
            result.assertions_count = getattr(test_case, 'assertions_count', 0)
            result.performance_metrics = getattr(test_case, 'performance_metrics', {})
            
        except unittest.SkipTest as e:
            result.status = TestStatus.SKIPPED
            result.error_message = str(e)
            
        except AssertionError as e:
            result.status = TestStatus.FAILED
            result.error_message = str(e)
            result.error_traceback = traceback.format_exc()
            
        except Exception as e:
            result.status = TestStatus.ERROR
            result.error_message = str(e)
            result.error_traceback = traceback.format_exc()
        
        finally:
            result.execution_time = time.time() - start_time
            
            # Try to get memory usage
            try:
                import psutil
                process = psutil.Process()
                result.memory_usage = process.memory_info().rss / 1024 / 1024  # MB
            except ImportError:
                pass
        
        # Log result
        status_emoji = {
            TestStatus.PASSED: "✅",
            TestStatus.FAILED: "❌", 
            TestStatus.SKIPPED: "⏭️",
            TestStatus.ERROR: "💥"
        }
        
        emoji = status_emoji.get(result.status, "❓")
        self.logger.info(f"{emoji} {test_name}: {result.status.value} ({result.execution_time:.3f}s)")
        
        if result.error_message:
            self.logger.error(f"   Error: {result.error_message}")
        
        return result
    
    def run_suite(self, test_suite: PrismTestSuite) -> List[TestResult]:
        """Run a complete test suite"""
        self.logger.info(f"Starting test suite: {test_suite.name}")
        results = test_suite.run_tests(self)
        
        # Generate report
        self._generate_test_report(test_suite, results)
        
        return results
    
    def _generate_test_report(self, test_suite: PrismTestSuite, results: List[TestResult]):
        """Generate test report"""
        summary = test_suite.get_summary()
        
        # JSON report
        json_report = {
            'summary': summary,
            'results': [result.to_dict() for result in results],
            'generated_at': time.time(),
            'generated_by': 'PrismTestRunner'
        }
        
        json_file = self.output_dir / f"{test_suite.name}_report.json"
        with open(json_file, 'w') as f:
            json.dump(json_report, f, indent=2)
        
        # HTML report
        html_report = self._generate_html_report(summary, results)
        html_file = self.output_dir / f"{test_suite.name}_report.html"
        with open(html_file, 'w') as f:
            f.write(html_report)
        
        self.logger.info(f"Test reports generated: {json_file}, {html_file}")
    
    def _generate_html_report(self, summary: Dict[str, Any], results: List[TestResult]) -> str:
        """Generate HTML test report"""
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>PRISM Test Report - {summary['suite_name']}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 10px; border-radius: 5px; }}
                .summary {{ margin: 20px 0; }}
                .test-result {{ margin: 10px 0; padding: 10px; border-left: 4px solid #ccc; }}
                .passed {{ border-color: #4CAF50; background-color: #f1f8e9; }}
                .failed {{ border-color: #F44336; background-color: #ffebee; }}
                .skipped {{ border-color: #FF9800; background-color: #fff3e0; }}
                .error {{ border-color: #9C27B0; background-color: #f3e5f5; }}
                .metrics {{ font-size: 0.9em; color: #666; }}
                pre {{ background-color: #f5f5f5; padding: 10px; border-radius: 3px; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>PRISM Test Report</h1>
                <h2>{summary['suite_name']}</h2>
            </div>
            
            <div class="summary">
                <h3>Summary</h3>
                <p><strong>Total Tests:</strong> {summary['total_tests']}</p>
                <p><strong>Passed:</strong> {summary['passed']} ({summary['success_rate']:.1f}%)</p>
                <p><strong>Failed:</strong> {summary['failed']}</p>
                <p><strong>Skipped:</strong> {summary['skipped']}</p>
                <p><strong>Errors:</strong> {summary['errors']}</p>
                <p><strong>Total Time:</strong> {summary['total_execution_time']:.3f}s</p>
            </div>
            
            <div class="results">
                <h3>Test Results</h3>
        """
        
        for result in results:
            status_class = result.status.value
            html += f"""
                <div class="test-result {status_class}">
                    <h4>{result.test_name}</h4>
                    <p><strong>Status:</strong> {result.status.value.upper()}</p>
                    <p><strong>Execution Time:</strong> {result.execution_time:.3f}s</p>
                    <p><strong>Assertions:</strong> {result.assertions_count}</p>
            """
            
            if result.error_message:
                html += f"<p><strong>Error:</strong> {result.error_message}</p>"
            
            if result.error_traceback:
                html += f"<pre>{result.error_traceback}</pre>"
            
            html += "</div>"
        
        html += """
            </div>
        </body>
        </html>
        """
        
        return html


# Convenience functions
def create_test_suite(name: str) -> PrismTestSuite:
    """Create a new test suite"""
    return PrismTestSuite(name)


def run_tests(test_modules: List, output_dir: Optional[Path] = None) -> List[TestResult]:
    """Run tests from multiple modules"""
    suite = PrismTestSuite("comprehensive_tests")
    
    for module in test_modules:
        suite.add_tests_from_module(module)
    
    runner = PrismTestRunner(output_dir)
    return runner.run_suite(suite)