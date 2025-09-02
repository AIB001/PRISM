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
from .core import PMFSystem, pmf_system
from .workflow import PMFWorkflow  
from .pmf_builder import PMFBuilder, pmf_builder
from .smd import SMDManager, SMDBuilder
from .umbrella import UmbrellaManager, WindowSelector
from .analyzer import PMFAnalyzer, WhamRunner
from .equilibration import PMFEquilibrationManager
from .runner import PMFRunner, run_pmf_workflow, run_pmf_step, create_pmf_config

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
    
    # Component managers
    "SMDManager", 
    "SMDBuilder",
    "UmbrellaManager",
    "WindowSelector", 
    "PMFAnalyzer",
    "WhamRunner",
    "PMFEquilibrationManager",
    
    # Utility functions
    "get_pmf_info",
]


def get_pmf_info():
    """Get information about PRISM PMF module"""
    info = {
        'version': '1.0.0',
        'architecture': 'PRISM-compatible',
        'execution_modes': ['automated', 'step_by_step'],
        'supported_forcefields': ['GAFF', 'OpenFF'],
        'analysis_features': [
            'SMD preparation and execution',
            'Adaptive umbrella window selection',
            'WHAM analysis with bootstrap errors',
            'PMF visualization and reporting',
            'Binding energy calculation'
        ],
        'integration': 'Seamless PRISM workflow integration'
    }
    return info