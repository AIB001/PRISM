#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein-Ligand Analysis Module - Trajectory analysis for protein-ligand systems
"""

import os
import multiprocessing

# Enable automatic OpenMP parallelization for MDTraj calculations
# Set OMP_NUM_THREADS to use all available CPU cores by default
if 'OMP_NUM_THREADS' not in os.environ:
    n_cores = multiprocessing.cpu_count()
    os.environ['OMP_NUM_THREADS'] = str(n_cores)
    print(f"PRISM: Enabled OpenMP parallelization with {n_cores} CPU cores")

from .config import AnalysisConfig, convert_numpy_types
from .trajectory import TrajectoryManager
from .trajectory_processor import TrajectoryProcessor, process_trajectory_simple, batch_process_trajectories
# from .multisys import MultiSystemAnalyzer  # TODO: Fix MDTraj conversion
from .export import DataExporter
from .analyzer import IntegratedProteinLigandAnalyzer, analyze_and_visualize

# Import calculation submodule
from . import calc
from .calc import (
    ContactAnalyzer, HBondAnalyzer, DistanceAnalyzer,
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
    # 'MultiSystemAnalyzer',  # TODO: Fix MDTraj conversion

    # New analysis modules
    'RMSDAnalyzer',
    'ClusteringAnalyzer',
    'StructuralAnalyzer',
    'SASAAnalyzer',
    'DihedralAnalyzer',

    # Utilities
    'TrajectoryManager',
    'TrajectoryProcessor',
    'process_trajectory_simple',
    'batch_process_trajectories',
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