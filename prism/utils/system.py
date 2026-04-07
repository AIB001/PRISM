#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM System Builder - Handles GROMACS system setup, solvation, and ion addition.

This module has been refactored into a modular package structure.
For backward compatibility, SystemBuilder is re-exported from the system package.

New structure:
- prism.utils.system.base: Base class with common utilities
- prism.utils.system.protein: Protein processing
- prism.utils.system.topology: Topology generation and fixing
- prism.utils.system.coordinates: Coordinate file processing
- prism.utils.system.solvation: Solvation and ion addition
- prism.utils.system.metals: Metal ion processing
- prism.utils.system.builder: Main SystemBuilder class
"""

# Re-export SystemBuilder for backward compatibility
from .system import SystemBuilder

__all__ = ["SystemBuilder"]
