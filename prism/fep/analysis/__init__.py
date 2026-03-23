#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM-FEP Analysis Module

Analysis tools for FEP results including XVG parsing, free energy estimators,
and comprehensive HTML report generation.

Main Components:
--------------
- FEPAnalyzer: Main analyzer class for complete FEP workflow
- GROMACSParser: Parse GROMACS output files (xvg, edr)
- FreeEnergyEstimator: Compute ΔG using BAR, MBAR, or TI
- HTMLReportGenerator: Generate publication-ready HTML reports

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

from .analyzer import FEPAnalyzer, FEResults, MultiEstimatorResults, FEPMultiEstimatorAnalyzer
from .xvg_parser import XVGParser
from .estimators import FEstimator
from .profiles import build_lambda_profiles, extract_lambda_data
from .report import HTMLReportGenerator, MultiEstimatorReportGenerator
from .cli import main

__all__ = [
    "FEPAnalyzer",
    "FEResults",
    "MultiEstimatorResults",
    "FEPMultiEstimatorAnalyzer",
    "XVGParser",
    "FEstimator",
    "build_lambda_profiles",
    "extract_lambda_data",
    "HTMLReportGenerator",
    "MultiEstimatorReportGenerator",
    "main",
]

__version__ = "0.1.0"
