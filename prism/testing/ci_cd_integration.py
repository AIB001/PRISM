#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM CI/CD Integration and Automation

Provides automation tools for continuous integration, continuous deployment,
and automated testing workflows for PRISM development.
"""

import os
import json
import yaml
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, asdict
from enum import Enum
import time

from ..utils.logging_system import PrismLogger


class CIStage(Enum):
    """CI/CD pipeline stages"""
    LINT = "lint"
    UNIT_TESTS = "unit_tests"
    INTEGRATION_TESTS = "integration_tests"
    PERFORMANCE_TESTS = "performance_tests"
    QUALITY_CHECKS = "quality_checks"
    SECURITY_SCAN = "security_scan"
    BUILD = "build"
    DEPLOY = "deploy"


class PipelineStatus(Enum):
    """Pipeline execution status"""
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class PipelineStage:
    """CI/CD pipeline stage configuration"""
    name: str
    stage_type: CIStage
    commands: List[str]
    dependencies: List[str]
    timeout_minutes: int = 30
    allow_failure: bool = False
    artifacts: List[str] = None
    environment: Dict[str, str] = None
    
    def __post_init__(self):
        if self.artifacts is None:
            self.artifacts = []
        if self.environment is None:
            self.environment = {}


@dataclass
class PipelineResult:
    """Result of pipeline stage execution"""
    stage_name: str
    status: PipelineStatus
    execution_time: float
    output: str
    error_output: str
    exit_code: int
    artifacts_created: List[str]


class GitHubActionsGenerator:
    """Generates GitHub Actions workflow files"""
    
    def __init__(self):
        self.logger = PrismLogger("github_actions_generator")
    
    def generate_test_workflow(self, output_file: Path) -> None:
        """Generate GitHub Actions workflow for testing"""
        
        workflow = {
            'name': 'PRISM Tests',
            'on': {
                'push': {
                    'branches': ['main', 'develop']
                },
                'pull_request': {
                    'branches': ['main', 'develop']
                }
            },
            'jobs': {
                'test': {
                    'runs-on': 'ubuntu-latest',
                    'strategy': {
                        'matrix': {
                            'python-version': ['3.8', '3.9', '3.10', '3.11']
                        }
                    },
                    'steps': [
                        {
                            'uses': 'actions/checkout@v3'
                        },
                        {
                            'name': 'Set up Python ${{ matrix.python-version }}',
                            'uses': 'actions/setup-python@v4',
                            'with': {
                                'python-version': '${{ matrix.python-version }}'
                            }
                        },
                        {
                            'name': 'Install dependencies',
                            'run': '\n'.join([
                                'python -m pip install --upgrade pip',
                                'pip install -r requirements.txt',
                                'pip install -r requirements-dev.txt'
                            ])
                        },
                        {
                            'name': 'Lint with flake8',
                            'run': '\n'.join([
                                'flake8 prism --count --select=E9,F63,F7,F82 --show-source --statistics',
                                'flake8 prism --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics'
                            ])
                        },
                        {
                            'name': 'Type check with mypy',
                            'run': 'mypy prism --ignore-missing-imports'
                        },
                        {
                            'name': 'Test with pytest',
                            'run': '\n'.join([
                                'pytest tests/ --cov=prism --cov-report=xml --cov-report=html',
                                'python -m prism.testing.integration_testing'
                            ])
                        },
                        {
                            'name': 'Upload coverage to Codecov',
                            'uses': 'codecov/codecov-action@v3',
                            'with': {
                                'file': './coverage.xml',
                                'flags': 'unittests',
                                'name': 'codecov-umbrella'
                            }
                        },
                        {
                            'name': 'Archive test results',
                            'uses': 'actions/upload-artifact@v3',
                            'with': {
                                'name': 'test-results-${{ matrix.python-version }}',
                                'path': 'test_results/'
                            }
                        }
                    ]
                },
                'quality-check': {
                    'runs-on': 'ubuntu-latest',
                    'steps': [
                        {
                            'uses': 'actions/checkout@v3'
                        },
                        {
                            'name': 'Set up Python',
                            'uses': 'actions/setup-python@v4',
                            'with': {
                                'python-version': '3.10'
                            }
                        },
                        {
                            'name': 'Install dependencies',
                            'run': '\n'.join([
                                'python -m pip install --upgrade pip',
                                'pip install -r requirements.txt',
                                'pip install -r requirements-dev.txt'
                            ])
                        },
                        {
                            'name': 'Run quality checks',
                            'run': 'python -m prism.testing.quality_assurance'
                        },
                        {
                            'name': 'Upload quality report',
                            'uses': 'actions/upload-artifact@v3',
                            'with': {
                                'name': 'quality-report',
                                'path': 'quality_report.html'
                            }
                        }
                    ]
                }
            }
        }
        
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            yaml.dump(workflow, f, default_flow_style=False, sort_keys=False)
        
        self.logger.info(f"Generated GitHub Actions workflow: {output_file}")
    
    def generate_performance_workflow(self, output_file: Path) -> None:
        """Generate performance testing workflow"""
        
        workflow = {
            'name': 'PRISM Performance Tests',
            'on': {
                'schedule': [
                    {'cron': '0 2 * * *'}  # Run daily at 2 AM
                ],
                'workflow_dispatch': {}  # Allow manual triggering
            },
            'jobs': {
                'performance': {
                    'runs-on': 'ubuntu-latest',
                    'steps': [
                        {
                            'uses': 'actions/checkout@v3'
                        },
                        {
                            'name': 'Set up Python',
                            'uses': 'actions/setup-python@v4',
                            'with': {
                                'python-version': '3.10'
                            }
                        },
                        {
                            'name': 'Install dependencies',
                            'run': '\n'.join([
                                'python -m pip install --upgrade pip',
                                'pip install -r requirements.txt',
                                'pip install -r requirements-dev.txt'
                            ])
                        },
                        {
                            'name': 'Run performance benchmarks',
                            'run': 'python -m prism.testing.performance_testing'
                        },
                        {
                            'name': 'Compare with baseline',
                            'run': '\n'.join([
                                'python scripts/compare_performance.py',
                                'python scripts/update_baseline.py'
                            ])
                        },
                        {
                            'name': 'Upload benchmark results',
                            'uses': 'actions/upload-artifact@v3',
                            'with': {
                                'name': 'benchmark-results',
                                'path': 'benchmark_results/'
                            }
                        }
                    ]
                }
            }
        }
        
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            yaml.dump(workflow, f, default_flow_style=False, sort_keys=False)
        
        self.logger.info(f"Generated performance workflow: {output_file}")


class JenkinsGenerator:
    """Generates Jenkins pipeline files"""
    
    def __init__(self):
        self.logger = PrismLogger("jenkins_generator")
    
    def generate_jenkinsfile(self, output_file: Path) -> None:
        """Generate Jenkinsfile for CI/CD pipeline"""
        
        jenkinsfile_content = '''
pipeline {
    agent any
    
    environment {
        PYTHON_VERSION = '3.10'
        PRISM_ENV = 'test'
    }
    
    stages {
        stage('Checkout') {
            steps {
                checkout scm
            }
        }
        
        stage('Setup Environment') {
            steps {
                sh '''
                    python -m venv venv
                    . venv/bin/activate
                    pip install --upgrade pip
                    pip install -r requirements.txt
                    pip install -r requirements-dev.txt
                '''
            }
        }
        
        stage('Lint') {
            steps {
                sh '''
                    . venv/bin/activate
                    flake8 prism --count --select=E9,F63,F7,F82 --show-source --statistics
                    flake8 prism --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics
                '''
            }
        }
        
        stage('Unit Tests') {
            steps {
                sh '''
                    . venv/bin/activate
                    pytest tests/unit/ --junitxml=test-results/unit-tests.xml --cov=prism --cov-report=xml
                '''
            }
            post {
                always {
                    junit 'test-results/unit-tests.xml'
                    publishCoverage adapters: [coberturaAdapter('coverage.xml')], sourceFileResolver: sourceFiles('STORE_LAST_BUILD')
                }
            }
        }
        
        stage('Integration Tests') {
            steps {
                sh '''
                    . venv/bin/activate
                    pytest tests/integration/ --junitxml=test-results/integration-tests.xml
                '''
            }
            post {
                always {
                    junit 'test-results/integration-tests.xml'
                }
            }
        }
        
        stage('Quality Analysis') {
            parallel {
                stage('Code Quality') {
                    steps {
                        sh '''
                            . venv/bin/activate
                            python -m prism.testing.quality_assurance
                        '''
                    }
                    post {
                        always {
                            publishHTML([
                                allowMissing: false,
                                alwaysLinkToLastBuild: true,
                                keepAll: true,
                                reportDir: 'quality_reports',
                                reportFiles: 'quality_report.html',
                                reportName: 'Code Quality Report'
                            ])
                        }
                    }
                }
                
                stage('Performance Tests') {
                    steps {
                        sh '''
                            . venv/bin/activate
                            python -m prism.testing.performance_testing
                        '''
                    }
                    post {
                        always {
                            publishHTML([
                                allowMissing: false,
                                alwaysLinkToLastBuild: true,
                                keepAll: true,
                                reportDir: 'benchmark_results',
                                reportFiles: 'benchmark_report.html',
                                reportName: 'Performance Report'
                            ])
                        }
                    }
                }
            }
        }
        
        stage('Build Package') {
            steps {
                sh '''
                    . venv/bin/activate
                    python setup.py sdist bdist_wheel
                '''
            }
            post {
                success {
                    archiveArtifacts artifacts: 'dist/*', fingerprint: true
                }
            }
        }
        
        stage('Deploy to Test') {
            when {
                branch 'develop'
            }
            steps {
                sh '''
                    . venv/bin/activate
                    pip install dist/*.whl
                    python -c "import prism; print(prism.__version__)"
                '''
            }
        }
    }
    
    post {
        always {
            cleanWs()
        }
        success {
            emailext (
                subject: "SUCCESS: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]'",
                body: "Build succeeded: ${env.BUILD_URL}",
                to: "${env.CHANGE_AUTHOR_EMAIL}"
            )
        }
        failure {
            emailext (
                subject: "FAILURE: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]'",
                body: "Build failed: ${env.BUILD_URL}",
                to: "${env.CHANGE_AUTHOR_EMAIL}"
            )
        }
    }
}
'''
        
        with open(output_file, 'w') as f:
            f.write(jenkinsfile_content.strip())
        
        self.logger.info(f"Generated Jenkinsfile: {output_file}")


class PipelineExecutor:
    """Local pipeline executor for testing CI/CD configurations"""
    
    def __init__(self, working_dir: Path):
        self.working_dir = working_dir
        self.logger = PrismLogger("pipeline_executor")
        self.stages: List[PipelineStage] = []
        self.results: List[PipelineResult] = []
    
    def add_stage(self, stage: PipelineStage) -> None:
        """Add a pipeline stage"""
        self.stages.append(stage)
        self.logger.info(f"Added pipeline stage: {stage.name}")
    
    def execute_pipeline(self) -> Dict[str, Any]:
        """Execute the complete pipeline"""
        self.logger.info(f"Starting pipeline execution with {len(self.stages)} stages")
        
        start_time = time.time()
        overall_status = PipelineStatus.SUCCESS
        
        for stage in self.stages:
            # Check dependencies
            if not self._check_dependencies(stage):
                self.logger.error(f"Dependencies not met for stage: {stage.name}")
                overall_status = PipelineStatus.FAILED
                break
            
            # Execute stage
            result = self._execute_stage(stage)
            self.results.append(result)
            
            if result.status == PipelineStatus.FAILED and not stage.allow_failure:
                overall_status = PipelineStatus.FAILED
                self.logger.error(f"Pipeline failed at stage: {stage.name}")
                break
        
        total_time = time.time() - start_time
        
        pipeline_summary = {
            'status': overall_status.value,
            'total_time': total_time,
            'stages_executed': len(self.results),
            'stages_passed': len([r for r in self.results if r.status == PipelineStatus.SUCCESS]),
            'stages_failed': len([r for r in self.results if r.status == PipelineStatus.FAILED]),
            'results': [self._result_to_dict(r) for r in self.results]
        }
        
        self.logger.info(f"Pipeline completed: {overall_status.value} in {total_time:.2f}s")
        return pipeline_summary
    
    def _execute_stage(self, stage: PipelineStage) -> PipelineResult:
        """Execute a single pipeline stage"""
        self.logger.info(f"Executing stage: {stage.name}")
        
        start_time = time.time()
        artifacts_created = []
        
        try:
            # Set environment variables
            env = os.environ.copy()
            env.update(stage.environment)
            
            # Execute commands
            output_lines = []
            error_lines = []
            final_exit_code = 0
            
            for command in stage.commands:
                self.logger.info(f"Running command: {command}")
                
                result = subprocess.run(
                    command,
                    shell=True,
                    cwd=self.working_dir,
                    env=env,
                    capture_output=True,
                    text=True,
                    timeout=stage.timeout_minutes * 60
                )
                
                output_lines.append(f"=== Command: {command} ===")
                output_lines.append(result.stdout)
                
                if result.stderr:
                    error_lines.append(f"=== Command: {command} ===")
                    error_lines.append(result.stderr)
                
                if result.returncode != 0:
                    final_exit_code = result.returncode
                    if not stage.allow_failure:
                        break
            
            # Check for artifacts
            for artifact_pattern in stage.artifacts:
                artifact_path = self.working_dir / artifact_pattern
                if artifact_path.exists():
                    artifacts_created.append(str(artifact_path))
            
            execution_time = time.time() - start_time
            status = PipelineStatus.SUCCESS if final_exit_code == 0 else PipelineStatus.FAILED
            
            return PipelineResult(
                stage_name=stage.name,
                status=status,
                execution_time=execution_time,
                output='\n'.join(output_lines),
                error_output='\n'.join(error_lines),
                exit_code=final_exit_code,
                artifacts_created=artifacts_created
            )
            
        except subprocess.TimeoutExpired:
            execution_time = time.time() - start_time
            return PipelineResult(
                stage_name=stage.name,
                status=PipelineStatus.FAILED,
                execution_time=execution_time,
                output="",
                error_output=f"Stage timed out after {stage.timeout_minutes} minutes",
                exit_code=-1,
                artifacts_created=artifacts_created
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            return PipelineResult(
                stage_name=stage.name,
                status=PipelineStatus.FAILED,
                execution_time=execution_time,
                output="",
                error_output=f"Stage execution failed: {str(e)}",
                exit_code=-1,
                artifacts_created=artifacts_created
            )
    
    def _check_dependencies(self, stage: PipelineStage) -> bool:
        """Check if stage dependencies are met"""
        for dependency in stage.dependencies:
            # Check if dependency stage passed
            dependency_result = next((r for r in self.results if r.stage_name == dependency), None)
            if not dependency_result or dependency_result.status != PipelineStatus.SUCCESS:
                return False
        return True
    
    def _result_to_dict(self, result: PipelineResult) -> Dict[str, Any]:
        """Convert result to dictionary"""
        return {
            'stage_name': result.stage_name,
            'status': result.status.value,
            'execution_time': result.execution_time,
            'exit_code': result.exit_code,
            'artifacts_created': result.artifacts_created,
            'has_output': len(result.output) > 0,
            'has_errors': len(result.error_output) > 0
        }


class CICDManager:
    """Main CI/CD configuration manager"""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.logger = PrismLogger("cicd_manager")
        self.github_generator = GitHubActionsGenerator()
        self.jenkins_generator = JenkinsGenerator()
    
    def setup_github_actions(self) -> None:
        """Setup GitHub Actions workflows"""
        workflows_dir = self.project_root / ".github" / "workflows"
        
        # Generate test workflow
        self.github_generator.generate_test_workflow(workflows_dir / "tests.yml")
        
        # Generate performance workflow
        self.github_generator.generate_performance_workflow(workflows_dir / "performance.yml")
        
        self.logger.info("GitHub Actions workflows configured")
    
    def setup_jenkins(self) -> None:
        """Setup Jenkins pipeline"""
        jenkinsfile = self.project_root / "Jenkinsfile"
        self.jenkins_generator.generate_jenkinsfile(jenkinsfile)
        
        self.logger.info("Jenkins pipeline configured")
    
    def create_test_pipeline(self) -> PipelineExecutor:
        """Create a test pipeline for local execution"""
        executor = PipelineExecutor(self.project_root)
        
        # Lint stage
        lint_stage = PipelineStage(
            name="lint",
            stage_type=CIStage.LINT,
            commands=[
                "python -m flake8 prism --count --statistics",
                "python -m mypy prism --ignore-missing-imports"
            ],
            dependencies=[],
            timeout_minutes=5
        )
        executor.add_stage(lint_stage)
        
        # Unit tests stage
        unit_test_stage = PipelineStage(
            name="unit_tests",
            stage_type=CIStage.UNIT_TESTS,
            commands=[
                "python -m pytest tests/unit/ -v --tb=short"
            ],
            dependencies=["lint"],
            timeout_minutes=10,
            artifacts=["test_results/unit_tests.xml"]
        )
        executor.add_stage(unit_test_stage)
        
        # Integration tests stage
        integration_test_stage = PipelineStage(
            name="integration_tests",
            stage_type=CIStage.INTEGRATION_TESTS,
            commands=[
                "python -m prism.testing.integration_testing"
            ],
            dependencies=["unit_tests"],
            timeout_minutes=15,
            artifacts=["test_results/integration_tests.xml"]
        )
        executor.add_stage(integration_test_stage)
        
        # Quality checks stage
        quality_stage = PipelineStage(
            name="quality_checks",
            stage_type=CIStage.QUALITY_CHECKS,
            commands=[
                "python -m prism.testing.quality_assurance"
            ],
            dependencies=["unit_tests"],
            timeout_minutes=10,
            artifacts=["quality_report.html"],
            allow_failure=True
        )
        executor.add_stage(quality_stage)
        
        return executor
    
    def generate_requirements_files(self) -> None:
        """Generate requirements files for CI/CD"""
        
        # Main requirements
        requirements = [
            "numpy>=1.20.0",
            "scipy>=1.7.0",
            "matplotlib>=3.4.0",
            "pandas>=1.3.0",
            "pyyaml>=5.4.0",
            "click>=8.0.0",
            "tqdm>=4.60.0"
        ]
        
        requirements_file = self.project_root / "requirements.txt"
        with open(requirements_file, 'w') as f:
            f.write('\n'.join(requirements))
        
        # Development requirements
        dev_requirements = [
            "pytest>=6.2.0",
            "pytest-cov>=2.12.0",
            "flake8>=3.9.0",
            "mypy>=0.900",
            "black>=21.0.0",
            "isort>=5.9.0",
            "pylint>=2.9.0",
            "coverage>=5.5.0",
            "psutil>=5.8.0"
        ]
        
        dev_requirements_file = self.project_root / "requirements-dev.txt"
        with open(dev_requirements_file, 'w') as f:
            f.write('\n'.join(dev_requirements))
        
        self.logger.info("Requirements files generated")


# Convenience functions
def setup_ci_cd(project_root: Path) -> CICDManager:
    """Setup CI/CD for PRISM project"""
    manager = CICDManager(project_root)
    manager.generate_requirements_files()
    manager.setup_github_actions()
    manager.setup_jenkins()
    return manager


def test_local_pipeline(project_root: Path) -> Dict[str, Any]:
    """Test pipeline locally"""
    manager = CICDManager(project_root)
    pipeline = manager.create_test_pipeline()
    return pipeline.execute_pipeline()