#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Backward-compatible exports for FEP estimator helpers.

This wrapper preserves the older import path:
`prism.fep.analysis.estimators`.
"""

from .core.estimators import (
    FEstimator,
    compute_free_energy_alchemlyb,
    compute_free_energy_gmx_bar,
    summarize_repeat_results,
)

__all__ = [
    "FEstimator",
    "compute_free_energy_alchemlyb",
    "compute_free_energy_gmx_bar",
    "summarize_repeat_results",
]
