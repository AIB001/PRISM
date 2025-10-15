#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structural analysis plotting utilities (Backward Compatibility Layer).

This module now serves as a facade, importing from modularized submodules
under `prism.analysis.plots.structural/` for better organization.

Original file was 2340 lines. Now refactored into:
- rmsd_rmsf_plots.py (621 lines)
- contact_hbond_plots.py (1103 lines)
- dihedral_sasa_plots.py (299 lines)
- specialized_plots.py (286 lines)

For new code, prefer importing directly from submodules:
    from prism.analysis.plots.structural import plot_rmsd_time_series

This file maintains backward compatibility for existing code.
"""

# Re-export all functions from submodules for backward compatibility
from .structural import *

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
