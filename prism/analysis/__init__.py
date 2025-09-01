# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# """
# PRISM Analysis Module - Trajectory analysis for protein-ligand systems
# """

# from .traj_analysis import TrajAnalysis
# from .contact import ContactAnalyzer
# from .hbond import HydrogenBondAnalyzer
# from .distance import DistanceAnalyzer
# from .io import DataExporter, ReportGenerator
# from .config import AnalysisConfig

# from . import visualization
# from .visualization import generate_html, HTMLGenerator

# __all__ = [
#     'TrajAnalysis',
#     'ContactAnalyzer',
#     'HydrogenBondAnalyzer',
#     'DistanceAnalyzer',
#     'DataExporter',
#     'ReportGenerator',
#     'AnalysisConfig',
#     'visualization' 
# ]

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein-Ligand Analysis Module - Trajectory analysis for protein-ligand systems
"""

from .config import AnalysisConfig, convert_numpy_types
from .trajectory import TrajectoryManager
from .contacts import ContactAnalyzer, DistanceCalculator
from .hbonds import HBondAnalyzer
from .distance import DistanceAnalyzer
from .multisys import MultiSystemAnalyzer
from .export import DataExporter
from .visualize import Visualizer
from .analyzer import IntegratedProteinLigandAnalyzer, analyze_and_visualize

# Import visualization submodule
from . import visualization
from .visualization import generate_html, HTMLGenerator

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
    
    # Utilities
    'TrajectoryManager',
    'DistanceCalculator',
    'DataExporter',
    'Visualizer',
    
    # HTML Visualization module
    'visualization',
    'generate_html',
    'HTMLGenerator'
]
