#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Protein Cleaner - Backward compatibility wrapper.

This module has been reorganized into a package structure.
All classes and functions are re-exported for backward compatibility.
"""

from .cleaner import ProteinCleaner, clean_protein_pdb, fix_terminal_atoms

__all__ = ["ProteinCleaner", "clean_protein_pdb", "fix_terminal_atoms"]
