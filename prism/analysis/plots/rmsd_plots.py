#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD/RMSF plotting utilities (Backward Compatibility Layer).

This module now serves as a facade, importing from `prism.analysis.plots.rmsd/`.

For new code, prefer importing directly from submodules:
    from prism.analysis.plots.rmsd import plot_rmsd_simple_timeseries

This file maintains backward compatibility for existing code.
"""

# Re-export all functions from rmsd submodule
from .rmsd import *

__all__ = [
    'plot_rmsd_simple_timeseries',
    'plot_rmsd_analysis',
    'plot_rmsf_analysis',
    'plot_rmsf_with_auto_chains',
    'plot_multi_trajectory_rmsf_comprehensive',
]
