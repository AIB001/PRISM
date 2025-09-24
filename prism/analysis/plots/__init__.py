#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRISM Analysis Plotting Module

Centralized plotting utilities for molecular dynamics trajectory analysis.
"""

from .basic import BasicPlotter, Visualizer
from .structural_plots import (
    plot_violin_comparison, plot_ramachandran, plot_dihedral_time_series,
    plot_sasa_comparison, plot_property_distribution, plot_rmsd_time_series,
    plot_rmsf_per_residue, plot_rmsd_rmsf_combined, plot_multi_chain_rmsf,
    plot_separate_rmsd, plot_multi_repeat_ligand_rmsd
)
from .comparison_plots import (
    plot_multi_system_comparison, plot_correlation_matrix,
    plot_statistical_comparison, plot_difference_analysis
)

__all__ = [
    # Basic plotting
    'BasicPlotter',
    'Visualizer',  # Legacy compatibility

    # Structural plotting
    'plot_violin_comparison',
    'plot_ramachandran',
    'plot_dihedral_time_series',
    'plot_sasa_comparison',
    'plot_property_distribution',
    'plot_rmsd_time_series',
    'plot_rmsf_per_residue',
    'plot_rmsd_rmsf_combined',
    'plot_multi_chain_rmsf',
    'plot_separate_rmsd',
    'plot_multi_repeat_ligand_rmsd',

    # Comparison plotting
    'plot_multi_system_comparison',
    'plot_correlation_matrix',
    'plot_statistical_comparison',
    'plot_difference_analysis'
]