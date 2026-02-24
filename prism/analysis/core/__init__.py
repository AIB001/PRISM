#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRISM Analysis Core Module

Core infrastructure components for molecular dynamics trajectory analysis.
"""

from .analyzer import IntegratedProteinLigandAnalyzer, analyze_and_visualize
from .config import AnalysisConfig, convert_numpy_types
from .export import DataExporter

# from .multisys import MultiSystemAnalyzer  # TODO: Fix MDTraj conversion
from .parallel import ParallelProcessor
from .trajectory import TrajectoryManager
from .trajectory_processor import TrajectoryProcessor, process_trajectory_simple, batch_process_trajectories
from .visualize import Visualizer

# Create default processor instance
default_processor = ParallelProcessor()

__all__ = [
    # Main analyzer
    "IntegratedProteinLigandAnalyzer",
    "analyze_and_visualize",
    # Configuration
    "AnalysisConfig",
    "convert_numpy_types",
    # Export utilities
    "DataExporter",
    # Multi-system analysis
    # 'MultiSystemAnalyzer',  # TODO: Fix MDTraj conversion
    # Parallel processing
    "ParallelProcessor",
    "default_processor",
    # Trajectory handling
    "TrajectoryManager",
    "TrajectoryProcessor",
    "process_trajectory_simple",
    "batch_process_trajectories",
    # Visualization utilities
    "Visualizer",
]
