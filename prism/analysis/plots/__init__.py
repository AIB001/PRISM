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
    plot_separate_rmsd, plot_multi_repeat_ligand_rmsd,
    plot_contact_probability_barplot, plot_contact_probability_heatmap,
    plot_contact_distance_distribution, plot_hydrogen_bond_analysis
)
from .rmsd_plots import (
    plot_rmsd_simple_timeseries, plot_rmsd_analysis, plot_rmsf_analysis,
    plot_rmsf_with_auto_chains
)
from .comparison_plots import (
    plot_multi_system_comparison, plot_correlation_matrix,
    plot_statistical_comparison, plot_difference_analysis
)
from .publication_utils import (
    get_publication_style, get_color_palette, setup_publication_figure,
    fix_rotated_labels, add_statistical_annotations, style_axes_for_publication,
    save_publication_figure, create_publication_colormap,
    PUBLICATION_FONTS, PUBLICATION_COLORS
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

    # RMSD plotting
    'plot_rmsd_simple_timeseries',
    'plot_rmsd_analysis',
    'plot_rmsf_analysis',
    'plot_rmsf_with_auto_chains',

    # Contact analysis plotting
    'plot_contact_probability_barplot',
    'plot_contact_probability_heatmap',
    'plot_contact_distance_distribution',
    'plot_hydrogen_bond_analysis',

    # Comparison plotting
    'plot_multi_system_comparison',
    'plot_correlation_matrix',
    'plot_statistical_comparison',
    'plot_difference_analysis',

    # Publication utilities
    'get_publication_style',
    'get_color_palette',
    'setup_publication_figure',
    'fix_rotated_labels',
    'add_statistical_annotations',
    'style_axes_for_publication',
    'save_publication_figure',
    'create_publication_colormap',
    'PUBLICATION_FONTS',
    'PUBLICATION_COLORS'
]