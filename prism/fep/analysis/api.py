#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Stable public API for PRISM-FEP analysis.

External callers should import analysis entrypoints from here (or from
``prism.fep.analysis``) rather than reaching into ``core``/``reports``/
``analyzers`` subpackages directly.
"""

from .analyzers import FEPAnalyzer, FEPMultiEstimatorAnalyzer
from .reports import HTMLReportGenerator, MultiBackendReportGenerator, MultiEstimatorReportGenerator
from .estimators import FEstimator
from .xvg_parser import XVGParser, find_xvg_file, parse_leg_data

__all__ = [
    "FEPAnalyzer",
    "FEPMultiEstimatorAnalyzer",
    "HTMLReportGenerator",
    "MultiBackendReportGenerator",
    "MultiEstimatorReportGenerator",
    "FEstimator",
    "XVGParser",
    "find_xvg_file",
    "parse_leg_data",
]
