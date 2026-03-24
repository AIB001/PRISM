#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM-FEP Analysis Module

Analysis tools for FEP results including XVG parsing, free energy estimators,
and comprehensive HTML report generation.

Main Components:
--------------
- FEPAnalyzer: Main analyzer class for complete FEP workflow
- FEPMultiEstimatorAnalyzer: Multi-estimator analysis with comparison
- HTMLReportGenerator: Generate publication-ready HTML reports
- MultiEstimatorReportGenerator: Tabbed reports for multiple estimators

Usage:
------
    from prism.fep.analysis import FEPAnalyzer

    analyzer = FEPAnalyzer(
        bound_dir='fep_project/bound',
        unbound_dir='fep_project/unbound',
        temperature=310.0,
        estimator='MBAR'
    )

    results = analyzer.analyze()
    report_path = analyzer.generate_html_report('report.html')

Command Line:
-------------
    prism --fep-analyze \\
        --bound-dir fep_project/bound \\
        --unbound-dir fep_project/unbound \\
        --output report.html \\
        --estimator MBAR

Dependencies:
------------
- alchemlyb: Free energy analysis library
- pandas, numpy: Data manipulation
- plotly: Interactive plotting (for HTML reports)

Install:
    pip install alchemlyb pandas numpy plotly

Author: PRISM Team
"""

# Core data models
from .core import FEResults, MultiEstimatorResults

# Analyzers
from .analyzers import FEPAnalyzer, FEPMultiEstimatorAnalyzer

# Report generators
from .reports import (
    HTMLReportGenerator,
    MultiEstimatorReportGenerator,
    MultiBackendReportGenerator,
)

# Other components
from . import cli
from . import plots
from .core import convergence, estimators, profiles, xvg_parser

__all__ = [
    # Core models
    "FEResults",
    "MultiEstimatorResults",
    # Analyzers
    "FEPAnalyzer",
    "FEPMultiEstimatorAnalyzer",
    # Reports
    "HTMLReportGenerator",
    "MultiEstimatorReportGenerator",
    "MultiBackendReportGenerator",
    # Other
    "cli",
    "convergence",
    "estimators",
    "plots",
    "profiles",
    "xvg_parser",
]
