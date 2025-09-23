#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein-Ligand Analysis Module - Trajectory analysis for protein-ligand systems
"""

from .config import AnalysisConfig, convert_numpy_types
from .trajectory import TrajectoryManager
from .multisys import MultiSystemAnalyzer
from .export import DataExporter
from .analyzer import IntegratedProteinLigandAnalyzer, analyze_and_visualize

# Import calculation submodule
from . import calc
from .calc import (
    ContactAnalyzer, DistanceCalculator, HBondAnalyzer, DistanceAnalyzer,
    RMSDAnalyzer, ClusteringAnalyzer, StructuralAnalyzer, SASAAnalyzer, DihedralAnalyzer
)

# Import plotting submodule
from . import plots
from .plots import (
    BasicPlotter, Visualizer,
    plot_violin_comparison, plot_ramachandran, plot_dihedral_time_series,
    plot_sasa_comparison, plot_property_distribution,
    plot_rmsd_time_series, plot_rmsf_per_residue, plot_rmsd_rmsf_combined,
    plot_multi_system_comparison, plot_correlation_matrix,
    plot_statistical_comparison, plot_difference_analysis
)

# Import contact visualization module (renamed from visualization)
from . import contact
from .contact import generate_html, HTMLGenerator

__all__ = [
    # Main interface
    'IntegratedProteinLigandAnalyzer',
    'analyze_and_visualize',

    # Configuration
    'AnalysisConfig',
    'convert_numpy_types',

    # Core analyzers
    'ContactAnalyzer',
    'HBondAnalyzer',
    'DistanceAnalyzer',
    'MultiSystemAnalyzer',

    # New analysis modules
    'RMSDAnalyzer',
    'ClusteringAnalyzer',
    'StructuralAnalyzer',
    'SASAAnalyzer',
    'DihedralAnalyzer',

    # Utilities
    'TrajectoryManager',
    'DistanceCalculator',
    'DataExporter',

    # Calculation submodule
    'calc',

    # Plotting submodule and functions
    'plots',
    'BasicPlotter',
    'Visualizer',
    'plot_violin_comparison',
    'plot_ramachandran',
    'plot_dihedral_time_series',
    'plot_sasa_comparison',
    'plot_property_distribution',
    'plot_rmsd_time_series',
    'plot_rmsf_per_residue',
    'plot_rmsd_rmsf_combined',
    'plot_multi_system_comparison',
    'plot_correlation_matrix',
    'plot_statistical_comparison',
    'plot_difference_analysis',

    # Contact visualization module (renamed from visualization)
    'contact',
    'generate_html',
    'HTMLGenerator'
]