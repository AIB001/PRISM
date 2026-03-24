#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP Analyzers - Single and Multi-Estimator Analysis

This module provides analyzer classes for FEP calculations.
"""

from .single import FEPAnalyzer
from .multi import FEPMultiEstimatorAnalyzer

__all__ = ["FEPAnalyzer", "FEPMultiEstimatorAnalyzer"]
