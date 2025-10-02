#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Module

This module provides comprehensive PMF calculation capabilities following
PRISM architectural patterns, with support for both automated and manual
execution modes.

Features:
- Integrated PMF workflow management
- SMD and umbrella sampling preparation
- WHAM analysis with error estimation
- Advanced visualization and reporting
- Flexible configuration system
"""

# Core classes following PRISM patterns
from .core import PMFSystem, pmf_system, PMFWorkflow
from .builders import PMFBuilder, pmf_builder, PMFEquilibrationManager
from .methods import SMDManager, SMDBuilder, UmbrellaManager, WindowSelector
from .analysis import PMFAnalyzer, WhamRunner, MDResultsValidator, validate_md_results, ValidationLevel
from .runner import PMFRunner, run_pmf_workflow, run_pmf_step, create_pmf_config

# Modular interface system
from .interfaces import (
    ModulePhase, ModuleStatus, ModuleResult, ModuleInterface, PMFModuleInterface,
    ModuleRegistry, WorkflowOrchestrator, create_pmf_module, setup_pmf_workflow
)

# Workflow optimizer removed - was incomplete and unused

__all__ = [
    # High-level interface (main entry points)
    "PMFSystem",
    "pmf_system",

    # Unified workflow runner (recommended API)
    "PMFRunner",
    "run_pmf_workflow",
    "run_pmf_step",
    "create_pmf_config",

    # Workflow management
    "PMFWorkflow",

    # System building
    "PMFBuilder",
    "pmf_builder",
    "PMFEquilibrationManager",

    # Component managers
    "SMDManager",
    "SMDBuilder",
    "UmbrellaManager",
    "WindowSelector",
    "PMFAnalyzer",
    "WhamRunner",

    # Validation and utilities
    "MDResultsValidator",
    "validate_md_results",
    "ValidationLevel",
    "get_pmf_info",

    # Modular interface system
    "ModulePhase",
    "ModuleStatus",
    "ModuleResult",
    "ModuleInterface",
    "PMFModuleInterface",
    "ModuleRegistry",
    "WorkflowOrchestrator",
    "create_pmf_module",
    "setup_pmf_workflow",
]


def get_pmf_info():
    """Get information about PRISM PMF module"""
    info = {
        'version': '1.1.0',
        'architecture': 'PRISM-compatible with enhanced box remodeling',
        'execution_modes': ['automated', 'step_by_step', 'pmf_optimized'],
        'supported_forcefields': ['GAFF', 'OpenFF'],
        'pmf_enhancements': [
            'Z-axis alignment for optimal pulling',
            'Extended box design for PMF calculations',
            'Protein+LIG extraction and remodeling',
            'Equilibration scripts (localrun.sh-style)',
            'PMF-optimized SMD protocols'
        ],
        'analysis_features': [
            'SMD preparation and execution',
            'Adaptive umbrella window selection',
            'WHAM analysis with bootstrap errors',
            'PMF visualization and reporting',
            'Binding energy calculation'
        ],
        'key_features': [
            'Automated box extension based on pulling distance',
            'Quality control through force field detection',
            'GROMACS built-in algorithms for box creation',
            'Support for equilibrated and unequilibrated systems'
        ],
        'integration': 'Seamless PRISM workflow integration'
    }
    return info