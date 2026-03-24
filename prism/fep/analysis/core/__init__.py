#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP Analysis Core

This module provides core data models and utility functions for FEP analysis.
"""

# Data models
from .models import FEResults, MultiEstimatorResults

# Utility modules
from . import estimators
from . import convergence
from . import profiles
from . import xvg_parser

__all__ = [
    # Data models
    "FEResults",
    "MultiEstimatorResults",
    # Utility modules
    "estimators",
    "convergence",
    "profiles",
    "xvg_parser",
]
