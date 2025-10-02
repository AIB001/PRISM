#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Testing Framework

Comprehensive testing system for PRISM components including unit tests,
integration tests, performance benchmarks, and quality assurance tools.
"""

from .test_framework import (
    PrismTestCase, PrismTestSuite, PrismTestRunner,
    TestType, TestStatus, TestResult
)
from .test_fixtures import (
    TestFixtures, MockDataGenerator, SystemFixtures
)
from .performance_testing import (
    PerformanceBenchmark, BenchmarkSuite, BenchmarkRunner,
    PerformanceMetrics, BenchmarkResult
)
from .quality_assurance import (
    CodeQualityChecker, StaticAnalyzer, CoverageAnalyzer,
    QualityReport, QualityMetrics
)
from .integration_testing import (
    IntegrationTestSuite, EndToEndTest, WorkflowTest
)

__all__ = [
    # Core testing framework
    "PrismTestCase",
    "PrismTestSuite", 
    "PrismTestRunner",
    "TestType",
    "TestStatus",
    "TestResult",
    
    # Test fixtures and data
    "TestFixtures",
    "MockDataGenerator",
    "SystemFixtures",
    
    # Performance testing
    "PerformanceBenchmark",
    "BenchmarkSuite",
    "BenchmarkRunner", 
    "PerformanceMetrics",
    "BenchmarkResult",
    
    # Quality assurance
    "CodeQualityChecker",
    "StaticAnalyzer",
    "CoverageAnalyzer",
    "QualityReport",
    "QualityMetrics",
    
    # Integration testing
    "IntegrationTestSuite",
    "EndToEndTest",
    "WorkflowTest",
]