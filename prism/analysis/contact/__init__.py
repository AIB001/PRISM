#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Visualization Module - Enhanced HTML visualization for trajectory analysis
"""

from .htmlgen import generate_html, HTMLGenerator
from .contact_analyzer import FastContactAnalyzer
from .data_processor import DataProcessor
from .visualization_generator import VisualizationGenerator
from .html_builder import HTMLBuilder
from .interactions import (
    InteractionTyper,
    INTERACTION_TYPES,
    INTERACTION_STYLE,
    RESIDUE_CLASS_COLORS,
    residue_class,
)

__all__ = [
    "generate_html",
    "HTMLGenerator",
    "FastContactAnalyzer",
    "DataProcessor",
    "VisualizationGenerator",
    "HTMLBuilder",
    "InteractionTyper",
    "INTERACTION_TYPES",
    "INTERACTION_STYLE",
    "RESIDUE_CLASS_COLORS",
    "residue_class",
]
