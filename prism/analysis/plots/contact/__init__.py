#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact analysis plotting subpackage.

Provides modular plotting functions for contact analysis,
organized by visualization type.
"""

# Import all functions from submodules
from .contact_violin_plots import (
    plot_comparison_contact_distances_violin,
    plot_key_residue_contact_violin,
    plot_residue_contact_numbers_violin,
    plot_contact_distances_violin
)

from .contact_heatmap_plots import (
    plot_contact_heatmap_annotated
)

from .contact_analysis_plots import (
    plot_contact_analysis
)
from .contact_timeseries_plots import (
    plot_contact_numbers_timeseries,
    plot_residue_distance_timeseries
)

from .contact_analysis_plots import (
    plot_grouped_contact_bars,
    plot_key_residue_contacts,
    plot_key_residue_contact_distribution,
    plot_contact_distances_raincloud
)

__all__ = [
    # Violin plots
    'plot_comparison_contact_distances_violin',
    'plot_key_residue_contact_violin',
    'plot_residue_contact_numbers_violin',
    'plot_contact_distances_violin',

    # Heatmap plots
    'plot_contact_heatmap_annotated',

    # Timeseries plots
    'plot_contact_analysis',
    'plot_contact_numbers_timeseries',
    'plot_residue_distance_timeseries',

    # Analysis plots
    'plot_grouped_contact_bars',
    'plot_key_residue_contacts',
    'plot_key_residue_contact_distribution',
    'plot_contact_distances_raincloud',
]
