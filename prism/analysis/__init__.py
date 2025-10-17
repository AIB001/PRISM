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

# Core infrastructure components (maintain backward compatibility)
from .core.config import AnalysisConfig, convert_numpy_types
from .core.trajectory import TrajectoryManager
from .core.trajectory_processor import TrajectoryProcessor, process_trajectory_simple, batch_process_trajectories
# from .core.multisys import MultiSystemAnalyzer  # TODO: Fix MDTraj conversion
from .core.export import DataExporter
from .core.analyzer import IntegratedProteinLigandAnalyzer, analyze_and_visualize

# Calculation modules (maintain backward compatibility)
from .calc.contacts import ContactAnalyzer
from .calc.hbonds import HBondAnalyzer
from .calc.distance import DistanceAnalyzer
from .calc.rmsd import RMSDAnalyzer
from .calc.clustering import ClusteringAnalyzer
from .calc.structural import StructuralAnalyzer
from .calc.sasa import SASAAnalyzer
from .calc.dihedral import DihedralAnalyzer

# Re-export core components at top level for backward compatibility

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