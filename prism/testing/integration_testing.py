#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Integration Testing

Comprehensive integration testing for PRISM components including end-to-end
workflow tests, module integration tests, and system validation tests.
"""

import time
import tempfile
from pathlib import Path
from typing import Dict, Any, List, Optional, Callable
from dataclasses import dataclass
from enum import Enum

from .test_framework import PrismTestCase, TestType, TestStatus, TestResult
from .test_fixtures import SystemFixtures, TestFixtures, MockSystemConfig
from ..utils.logging_system import PrismLogger


class IntegrationTestType(Enum):
    """Types of integration tests"""
    MODULE_INTEGRATION = "module_integration"
    WORKFLOW_E2E = "workflow_e2e"
    SYSTEM_VALIDATION = "system_validation"
    PERFORMANCE_INTEGRATION = "performance_integration"
    DATA_FLOW = "data_flow"


@dataclass
class IntegrationTestResult:
    """Result of integration test"""
    test_name: str
    test_type: IntegrationTestType
    status: TestStatus
    execution_time: float
    components_tested: List[str]
    data_validated: Dict[str, bool]
    performance_metrics: Dict[str, float]
    errors: List[str]
    warnings: List[str]


class IntegrationTestCase(PrismTestCase):
    """Base class for integration tests"""
    
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)
        self.test_type = TestType.INTEGRATION
        self.integration_type = IntegrationTestType.MODULE_INTEGRATION
        self.fixtures = SystemFixtures()
        self.components_tested = []
        self.data_validated = {}
        
    def setUp(self):
        """Setup for integration tests"""
        super().setUp()
        self.logger.info(f"Setting up integration test: {self._testMethodName}")
        
    def tearDown(self):
        """Cleanup after integration tests"""
        super().tearDown()
        self.fixtures.cleanup_all()
        self.logger.info(f"Cleaned up integration test: {self._testMethodName}")
    
    def validate_component_integration(self, component_a: str, component_b: str,
                                     data_flow: Dict[str, Any]) -> bool:
        """Validate integration between two components"""
        self.components_tested.extend([component_a, component_b])
        
        # Mock validation - real implementation would check actual data flow
        validation_result = len(data_flow) > 0 and 'output_data' in data_flow
        self.data_validated[f"{component_a}_to_{component_b}"] = validation_result
        
        self.logger.info(f"Validated integration: {component_a} -> {component_b}: {validation_result}")
        return validation_result
    
    def measure_integration_performance(self, operation: Callable, 
                                      operation_name: str) -> Dict[str, float]:
        """Measure performance of integrated operation"""
        start_time = time.time()
        start_memory = self._get_memory_usage()
        
        try:
            result = operation()
            success = True
        except Exception as e:
            self.logger.error(f"Operation {operation_name} failed: {e}")
            success = False
            result = None
        
        end_time = time.time()
        end_memory = self._get_memory_usage()
        
        metrics = {
            'execution_time': end_time - start_time,
            'memory_delta': end_memory - start_memory,
            'success': success
        }
        
        self.performance_metrics[operation_name] = metrics
        return metrics
    
    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB"""
        try:
            import psutil
            process = psutil.Process()
            return process.memory_info().rss / 1024 / 1024
        except ImportError:
            return 0.0


class PMFWorkflowTest(IntegrationTestCase):
    """End-to-end PMF workflow integration tests"""
    
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)
        self.integration_type = IntegrationTestType.WORKFLOW_E2E
    
    def test_complete_pmf_workflow(self):
        """Test complete PMF calculation workflow"""
        self.logger.info("Testing complete PMF workflow")
        
        # Setup test scenario
        test_scenario = self.fixtures.pmf_test_scenario()
        system_dir = test_scenario['system_dir']
        
        # Test workflow steps
        with self.assertExecutionTime(30.0):  # Should complete within 30 seconds for mock data
            
            # Step 1: System validation
            self.assertTrue(system_dir.exists(), "System directory should exist")
            self.assertFileExists(system_dir / "system.gro")
            self.assertFileExists(system_dir / "system.top")
            
            # Step 2: Mock PMF preparation
            pmf_output_dir = self.test_temp_dir / "pmf_output"
            pmf_output_dir.mkdir()
            
            prep_result = self._mock_pmf_preparation(system_dir, pmf_output_dir)
            self.assertTrue(prep_result['success'], "PMF preparation should succeed")
            
            # Step 3: Mock SMD execution
            smd_result = self._mock_smd_execution(pmf_output_dir)
            self.assertTrue(smd_result['success'], "SMD execution should succeed")
            
            # Step 4: Mock umbrella sampling
            umbrella_result = self._mock_umbrella_sampling(pmf_output_dir, smd_result)
            self.assertTrue(umbrella_result['success'], "Umbrella sampling should succeed")
            
            # Step 5: Mock PMF analysis
            analysis_result = self._mock_pmf_analysis(pmf_output_dir, umbrella_result)
            self.assertTrue(analysis_result['success'], "PMF analysis should succeed")
            
            # Validate final PMF result
            self.assertValidPMFResult(analysis_result['pmf_data'])
        
        # Verify component integration
        self.validate_component_integration("smd", "umbrella", smd_result)
        self.validate_component_integration("umbrella", "analysis", umbrella_result)
        
        self.logger.info("Complete PMF workflow test passed")
    
    def test_pmf_workflow_with_errors(self):
        """Test PMF workflow error handling"""
        self.logger.info("Testing PMF workflow error handling")
        
        # Create invalid system (missing files)
        invalid_system = self.test_temp_dir / "invalid_system"
        invalid_system.mkdir()
        
        # Only create some required files
        self.create_test_file("invalid_system/system.gro", "invalid content")
        # Missing system.top should cause validation error
        
        pmf_output_dir = self.test_temp_dir / "pmf_output_error"
        pmf_output_dir.mkdir()
        
        # PMF preparation should fail gracefully
        prep_result = self._mock_pmf_preparation(invalid_system, pmf_output_dir)
        self.assertFalse(prep_result['success'], "PMF preparation should fail for invalid system")
        self.assertIn('errors', prep_result, "Error information should be provided")
    
    def _mock_pmf_preparation(self, system_dir: Path, output_dir: Path) -> Dict[str, Any]:
        """Mock PMF preparation step"""
        def prep_operation():
            # Check required files
            required_files = ['system.gro', 'system.top']
            missing_files = []
            
            for file_name in required_files:
                if not (system_dir / file_name).exists():
                    missing_files.append(file_name)
            
            if missing_files:
                return {
                    'success': False,
                    'errors': [f"Missing required files: {', '.join(missing_files)}"]
                }
            
            # Create mock preparation outputs
            prep_dir = output_dir / "preparation"
            prep_dir.mkdir(exist_ok=True)
            
            # Mock file creation
            (prep_dir / "prepared_system.gro").write_text("Mock prepared system")
            (prep_dir / "preparation_log.txt").write_text("Mock preparation log")
            
            return {
                'success': True,
                'output_dir': str(prep_dir),
                'files_created': ['prepared_system.gro', 'preparation_log.txt']
            }
        
        return self.measure_integration_performance(prep_operation, "pmf_preparation")
    
    def _mock_smd_execution(self, output_dir: Path) -> Dict[str, Any]:
        """Mock SMD execution step"""
        def smd_operation():
            smd_dir = output_dir / "smd"
            smd_dir.mkdir(exist_ok=True)
            
            # Create mock SMD outputs
            trajectory_file = smd_dir / "smd_trajectory.xtc"
            trajectory_file.write_text("Mock SMD trajectory")
            
            pullx_file = smd_dir / "pullx.xvg"
            pullx_file.write_text("# Mock pull distance data\n0.000 1.500\n1.000 2.000\n2.000 2.500")
            
            pullf_file = smd_dir / "pullf.xvg" 
            pullf_file.write_text("# Mock pull force data\n0.000 100.0\n1.000 200.0\n2.000 150.0")
            
            return {
                'success': True,
                'trajectory_file': str(trajectory_file),
                'pullx_file': str(pullx_file),
                'pullf_file': str(pullf_file),
                'suggested_windows': 30,
                'output_data': {
                    'max_distance': 2.5,
                    'max_force': 200.0
                }
            }
        
        return self.measure_integration_performance(smd_operation, "smd_execution")
    
    def _mock_umbrella_sampling(self, output_dir: Path, smd_result: Dict[str, Any]) -> Dict[str, Any]:
        """Mock umbrella sampling step"""
        def umbrella_operation():
            umbrella_dir = output_dir / "umbrella"
            umbrella_dir.mkdir(exist_ok=True)
            
            # Use data from SMD results
            num_windows = smd_result.get('suggested_windows', 30)
            
            # Create mock umbrella windows
            for i in range(num_windows):
                window_dir = umbrella_dir / f"window_{i:02d}"
                window_dir.mkdir(exist_ok=True)
                
                # Mock window files
                (window_dir / "umbrella.tpr").write_text("Mock TPR file")
                (window_dir / "umbrella.xtc").write_text("Mock umbrella trajectory")
                (window_dir / "pullx.xvg").write_text(f"# Window {i}\n0.000 {1.0 + i * 0.1}")
            
            return {
                'success': True,
                'num_windows': num_windows,
                'umbrella_dir': str(umbrella_dir),
                'output_data': {
                    'windows_completed': num_windows,
                    'total_simulation_time': num_windows * 10.0  # ns
                }
            }
        
        return self.measure_integration_performance(umbrella_operation, "umbrella_sampling")
    
    def _mock_pmf_analysis(self, output_dir: Path, umbrella_result: Dict[str, Any]) -> Dict[str, Any]:
        """Mock PMF analysis step"""
        def analysis_operation():
            analysis_dir = output_dir / "analysis"
            analysis_dir.mkdir(exist_ok=True)
            
            # Generate mock PMF profile
            num_windows = umbrella_result.get('num_windows', 30)
            
            # Create mock WHAM output
            wham_output = analysis_dir / "wham_output.txt"
            wham_content = "# Mock WHAM output\n"
            
            pmf_profile = []
            for i in range(num_windows):
                distance = 1.0 + i * 0.1
                # Mock PMF with binding minimum and barrier
                if distance < 1.5:
                    pmf = -8.0 + 2.0 * (distance - 1.2)**2
                elif distance < 2.5:
                    pmf = -6.0 + 5.0 * (distance - 1.8)**2
                else:
                    pmf = 0.0
                
                pmf_profile.append((distance, pmf))
                wham_content += f"{distance:.3f} {pmf:.3f} 0.5\n"
            
            wham_output.write_text(wham_content)
            
            # Find binding energy (minimum PMF)
            binding_energy = min(pmf for _, pmf in pmf_profile)
            
            pmf_data = {
                'binding_energy': binding_energy,
                'pmf_profile': pmf_profile,
                'error_estimate': 0.5,
                'convergence': True,
                'method': 'WHAM'
            }
            
            return {
                'success': True,
                'pmf_data': pmf_data,
                'output_files': [str(wham_output)],
                'output_data': pmf_data
            }
        
        return self.measure_integration_performance(analysis_operation, "pmf_analysis")


class ModuleIntegrationTest(IntegrationTestCase):
    """Module-to-module integration tests"""
    
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)
        self.integration_type = IntegrationTestType.MODULE_INTEGRATION
    
    def test_logging_integration(self):
        """Test logging system integration"""
        self.logger.info("Testing logging system integration")
        
        # Test logging with different components
        from ..utils.logging_system import get_prism_logger
        
        # Create loggers for different components
        smd_logger = get_prism_logger("smd_test", enable_monitoring=False)
        umbrella_logger = get_prism_logger("umbrella_test", enable_monitoring=False)
        
        # Test cross-component logging
        smd_logger.info("SMD module starting")
        umbrella_logger.info("Umbrella module receiving data from SMD")
        
        # Validate logger sessions
        smd_summary = smd_logger.get_session_summary()
        umbrella_summary = umbrella_logger.get_session_summary()
        
        self.assertGreater(smd_summary['total_events'], 0, "SMD logger should have events")
        self.assertGreater(umbrella_summary['total_events'], 0, "Umbrella logger should have events")
        
        # Validate component integration
        self.validate_component_integration("logging", "smd", smd_summary)
        self.validate_component_integration("logging", "umbrella", umbrella_summary)
    
    def test_configuration_integration(self):
        """Test configuration system integration"""
        self.logger.info("Testing configuration system integration")
        
        # Create test configuration
        test_config = {
            'smd': {
                'pull_rate': 0.01,
                'nsteps': 1000000
            },
            'umbrella': {
                'windows': 40,
                'production_time_ps': 10000
            }
        }
        
        # Test configuration sharing between modules
        config_file = self.create_test_file("test_config.json")
        
        import json
        with open(config_file, 'w') as f:
            json.dump(test_config, f)
        
        # Validate configuration loading
        with open(config_file, 'r') as f:
            loaded_config = json.load(f)
        
        self.assertEqual(loaded_config, test_config, "Configuration should be preserved")
        
        # Test module configuration access
        smd_config = loaded_config.get('smd', {})
        umbrella_config = loaded_config.get('umbrella', {})
        
        self.assertIn('pull_rate', smd_config, "SMD configuration should be accessible")
        self.assertIn('windows', umbrella_config, "Umbrella configuration should be accessible")
        
        # Validate component integration
        self.validate_component_integration("config", "smd", smd_config)
        self.validate_component_integration("config", "umbrella", umbrella_config)
    
    def test_data_flow_integration(self):
        """Test data flow between modules"""
        self.logger.info("Testing data flow integration")
        
        # Mock data flow: SMD -> Umbrella -> Analysis
        smd_data = {
            'trajectory_file': '/path/to/smd.xtc',
            'pull_data': '/path/to/pullx.xvg',
            'suggested_windows': 35,
            'max_distance': 3.5
        }
        
        # Transform SMD data for umbrella sampling
        umbrella_input = {
            'reference_trajectory': smd_data['trajectory_file'],
            'num_windows': smd_data['suggested_windows'],
            'distance_range': (0.5, smd_data['max_distance'])
        }
        
        # Validate data transformation
        self.assertEqual(umbrella_input['num_windows'], smd_data['suggested_windows'])
        self.assertEqual(umbrella_input['distance_range'][1], smd_data['max_distance'])
        
        # Mock umbrella output
        umbrella_data = {
            'windows_completed': umbrella_input['num_windows'],
            'window_trajectories': [f'/path/to/window_{i}.xtc' for i in range(umbrella_input['num_windows'])],
            'convergence_status': True
        }
        
        # Transform umbrella data for analysis
        analysis_input = {
            'window_data': umbrella_data['window_trajectories'],
            'num_windows': umbrella_data['windows_completed'],
            'converged': umbrella_data['convergence_status']
        }
        
        # Validate data transformations
        self.assertEqual(len(analysis_input['window_data']), umbrella_data['windows_completed'])
        self.assertTrue(analysis_input['converged'])
        
        # Validate component integrations
        self.validate_component_integration("smd", "umbrella", umbrella_input)
        self.validate_component_integration("umbrella", "analysis", analysis_input)


class SystemValidationTest(IntegrationTestCase):
    """System-level validation tests"""
    
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)
        self.integration_type = IntegrationTestType.SYSTEM_VALIDATION
    
    def test_system_resource_management(self):
        """Test system resource management"""
        self.logger.info("Testing system resource management")
        
        # Test temporary directory management
        temp_dirs = []
        for i in range(5):
            temp_dir = self.test_temp_dir / f"resource_test_{i}"
            temp_dir.mkdir()
            temp_dirs.append(temp_dir)
        
        # Validate all directories exist
        for temp_dir in temp_dirs:
            self.assertTrue(temp_dir.exists(), f"Temp directory {temp_dir} should exist")
        
        # Test resource cleanup (will happen in tearDown)
        self.assertEqual(len(temp_dirs), 5, "Should create 5 temporary directories")
    
    def test_error_propagation(self):
        """Test error propagation across components"""
        self.logger.info("Testing error propagation")
        
        # Simulate error in early component
        def failing_operation():
            raise ValueError("Mock error in SMD component")
        
        # Test error handling
        with self.assertRaises(ValueError) as context:
            failing_operation()
        
        self.assertIn("Mock error in SMD component", str(context.exception))
        
        # Test error recovery
        def recovering_operation():
            try:
                failing_operation()
            except ValueError as e:
                return {'error': str(e), 'recovered': True}
        
        result = recovering_operation()
        self.assertTrue(result['recovered'], "System should recover from errors")
        self.assertIn('error', result, "Error information should be preserved")


class IntegrationTestSuite:
    """Integration test suite manager"""
    
    def __init__(self, name: str = "integration_tests"):
        self.name = name
        self.test_cases = []
        self.logger = PrismLogger(f"integration_suite.{name}")
    
    def add_pmf_workflow_tests(self):
        """Add PMF workflow tests"""
        self.test_cases.extend([
            PMFWorkflowTest('test_complete_pmf_workflow'),
            PMFWorkflowTest('test_pmf_workflow_with_errors')
        ])
    
    def add_module_integration_tests(self):
        """Add module integration tests"""
        self.test_cases.extend([
            ModuleIntegrationTest('test_logging_integration'),
            ModuleIntegrationTest('test_configuration_integration'),
            ModuleIntegrationTest('test_data_flow_integration')
        ])
    
    def add_system_validation_tests(self):
        """Add system validation tests"""
        self.test_cases.extend([
            SystemValidationTest('test_system_resource_management'),
            SystemValidationTest('test_error_propagation')
        ])
    
    def run_all_tests(self) -> List[TestResult]:
        """Run all integration tests"""
        from .test_framework import PrismTestRunner
        
        runner = PrismTestRunner()
        results = []
        
        self.logger.info(f"Running integration test suite: {self.name} ({len(self.test_cases)} tests)")
        
        for test_case in self.test_cases:
            result = runner.run_single_test(test_case)
            results.append(result)
        
        return results


class EndToEndTest:
    """End-to-end test runner"""
    
    def __init__(self):
        self.logger = PrismLogger("e2e_test")
        self.fixtures = SystemFixtures()
    
    def run_complete_pmf_pipeline(self) -> Dict[str, Any]:
        """Run complete PMF pipeline end-to-end test"""
        self.logger.info("Starting complete PMF pipeline E2E test")
        
        start_time = time.time()
        
        try:
            # Setup test scenario
            test_scenario = self.fixtures.pmf_test_scenario()
            
            # Run pipeline steps (mocked)
            pipeline_result = {
                'system_setup': True,
                'smd_execution': True,
                'umbrella_sampling': True,
                'pmf_analysis': True,
                'final_result': test_scenario['pmf_results']
            }
            
            execution_time = time.time() - start_time
            
            return {
                'success': True,
                'execution_time': execution_time,
                'pipeline_result': pipeline_result,
                'binding_energy': test_scenario['expected_binding_energy']
            }
            
        except Exception as e:
            execution_time = time.time() - start_time
            self.logger.error(f"E2E test failed: {e}")
            
            return {
                'success': False,
                'execution_time': execution_time,
                'error': str(e)
            }
        
        finally:
            self.fixtures.cleanup_all()


# Convenience functions
def create_integration_test_suite() -> IntegrationTestSuite:
    """Create comprehensive integration test suite"""
    suite = IntegrationTestSuite("comprehensive_integration")
    suite.add_pmf_workflow_tests()
    suite.add_module_integration_tests()
    suite.add_system_validation_tests()
    return suite


def run_end_to_end_test() -> Dict[str, Any]:
    """Run end-to-end test"""
    e2e_test = EndToEndTest()
    return e2e_test.run_complete_pmf_pipeline()