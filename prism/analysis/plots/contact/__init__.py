#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact analysis plotting subpackage.

Hierarchical organization by visualization type.
"""

# Import from core
from .core.analysis import (
    plot_contact_analysis,
    plot_grouped_contact_bars,
    plot_key_residue_contacts,
    plot_key_residue_contact_distribution,
    plot_contact_distances_raincloud,
)

# Import from violin
from .violin.violin import (
    plot_comparison_contact_distances_violin,
    plot_key_residue_contact_violin,
    plot_residue_contact_numbers_violin,
    plot_contact_distances_violin,
)

# Import from heatmap
from .heatmap.probability import plot_contact_heatmap_annotated

# Import from timeseries
from .timeseries.contacts import plot_contact_numbers_timeseries, plot_residue_distance_timeseries

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
