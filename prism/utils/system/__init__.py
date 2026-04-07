#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM System Builder - Modular system building for GROMACS simulations.

This package provides a modular architecture for building protein-ligand systems:
- base: Base class with common utilities
- protein: Protein processing (pdbfixer, PROPKA, histidine handling)
- topology: Topology generation and fixing (pdb2gmx, multi-ligand support)
- coordinates: Coordinate file processing (combining, box creation)
- solvation: Solvation and ion addition
- metals: Metal ion processing
- builder: Main SystemBuilder class orchestrating the workflow
"""

from .builder import SystemBuilder

__all__ = ["SystemBuilder"]
