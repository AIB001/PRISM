#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Module

Two-API Design for PMF Calculations:
1. Remodeling API - System preparation and rebuilding
2. Runner API - Complete PMF workflow execution

Usage:
    # Step 1: Remodel MD system for PMF
    from prism.pmf import PMFBuilder
    builder = PMFBuilder(md_results_dir, output_dir)
    builder.build(equilibrate=True)

    # Step 2: Run PMF workflow
    from prism.pmf import PMFRunner
    runner = PMFRunner(config)
    results = runner.run_complete_workflow(system_dir, output_dir)
"""

# ============================================================================
# PRIMARY APIs - Only use these two entry points
# ============================================================================

# 1. Remodeling API - System preparation
from .builders import PMFBuilder

# 2. Runner API - PMF workflow execution
from .runner import PMFRunner, run_pmf_workflow, create_pmf_config

# ============================================================================
# Internal components - Do not use directly (for advanced users only)
# ============================================================================
from .core import PMFSystem, PMFWorkflow
from .builders import PMFEquilibrationManager
from .methods import SMDManager, UmbrellaManager
from .analysis import PMFAnalyzer, MDResultsValidator, validate_md_results, ValidationLevel

__all__ = [
    # ===== PRIMARY APIs (RECOMMENDED) =====
    # Remodeling API
    "PMFBuilder",

    # Runner API
    "PMFRunner",
    "run_pmf_workflow",
    "create_pmf_config",

    # ===== INTERNAL COMPONENTS (Advanced use only) =====
    "PMFSystem",
    "PMFWorkflow",
    "PMFEquilibrationManager",
    "SMDManager",
    "UmbrellaManager",
    "PMFAnalyzer",
    "MDResultsValidator",
    "validate_md_results",
    "ValidationLevel",
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