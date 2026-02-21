#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Module

Provides structure alignment and analysis for PMF (Potential of Mean Force)
calculations.  The PMF workflow (SMD, umbrella sampling, WHAM) is driven
by prism.builder and prism.utils.mdp.

Usage:
    # Via CLI (recommended):
    prism --pmf protein.pdb ligand.mol2 -o output

    # Direct alignment access (advanced):
    from prism.pmf.alignment import PMFAligner

    # Analyze PMF results:
    from prism.pmf.analysis import PMFAnalyzer
    analyzer = PMFAnalyzer("output/GMX_PROLIG_PMF")
    report = analyzer.generate_report()
"""

from .alignment import PMFAligner
from .analysis import PMFAnalyzer

__all__ = [
    "PMFAligner",
    "PMFAnalyzer",
]
