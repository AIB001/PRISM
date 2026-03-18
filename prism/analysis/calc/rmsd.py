#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD calculation for MD trajectories.

Based on SciDraft-Studio implementation with PRISM integration.

This module has been reorganized into a package structure.
All classes are re-exported for backward compatibility.
"""

# Re-export RMSDAnalyzer for backward compatibility
from .rmsd import RMSDAnalyzer

__all__ = ["RMSDAnalyzer"]
