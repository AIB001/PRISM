#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact analysis plotting utilities (Backward Compatibility Layer).

This module now serves as a facade, importing from modularized submodules
under `prism.analysis.plots.contact/` for better organization.

Original file was 1447 lines. Now refactored into:
- contact_violin_plots.py: Contact violin plots
- contact_heatmap_plots.py: Contact heatmaps
- contact_timeseries_plots.py: Contact timeseries
- contact_analysis_plots.py: Contact analysis plots

For new code, prefer importing directly from submodules:
    from prism.analysis.plots.contact import plot_contact_analysis

This file maintains backward compatibility for existing code.
"""

# Re-export all functions from submodules for backward compatibility
from .contact import *

__all__ = [
    # Violin plots
    "plot_comparison_contact_distances_violin",
    "plot_key_residue_contact_violin",
    "plot_residue_contact_numbers_violin",
    "plot_contact_distances_violin",
    # Heatmap plots
    "plot_contact_heatmap_annotated",
    # Timeseries plots
    "plot_contact_analysis",
    "plot_contact_numbers_timeseries",
    "plot_residue_distance_timeseries",
    # Analysis plots
    "plot_grouped_contact_bars",
    "plot_key_residue_contacts",
    "plot_key_residue_contact_distribution",
    "plot_contact_distances_raincloud",
]
