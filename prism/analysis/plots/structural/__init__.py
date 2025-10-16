#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structural visualization subpackage.

This module provides modular plotting functions for structural analysis,
organized by functionality.
"""

# Import all functions from submodules
from .rmsd_rmsf_plots import (
    plot_rmsd_time_series,
    plot_rmsf_per_residue,
    plot_rmsd_rmsf_combined,
    plot_multi_chain_rmsf,
    plot_separate_rmsd,
    plot_multi_repeat_ligand_rmsd,
    plot_multi_chain_rmsf_example_style
)

from .contact_hbond_plots import (
    generate_publication_contact_plots,
    plot_contact_probability_barplot,
    plot_contact_probability_heatmap,
    plot_contact_distance_distribution,
    plot_hydrogen_bond_analysis,
    plot_key_residue_distances,
    plot_hydrogen_bond_stability
)

from .dihedral_sasa_plots import (
    plot_ramachandran,
    plot_dihedral_time_series,
    plot_sasa_comparison,
    plot_property_distribution
)

from .specialized_plots import (
    plot_violin_comparison,
    plot_distance_time_series,
    plot_magnesium_coordination
)

__all__ = [
    # RMSD/RMSF plots
    'plot_rmsd_time_series',
    'plot_rmsf_per_residue',
    'plot_rmsd_rmsf_combined',
    'plot_multi_chain_rmsf',
    'plot_separate_rmsd',
    'plot_multi_repeat_ligand_rmsd',
    'plot_multi_chain_rmsf_example_style',

    # Contact/H-bond plots
    'generate_publication_contact_plots',
    'plot_contact_probability_barplot',
    'plot_contact_probability_heatmap',
    'plot_contact_distance_distribution',
    'plot_hydrogen_bond_analysis',
    'plot_key_residue_distances',
    'plot_hydrogen_bond_stability',

    # Dihedral/SASA plots
    'plot_ramachandran',
    'plot_dihedral_time_series',
    'plot_sasa_comparison',
    'plot_property_distribution',

    # Specialized plots
    'plot_violin_comparison',
    'plot_distance_time_series',
    'plot_magnesium_coordination',
]
