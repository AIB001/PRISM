#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM-FEP Core Module

Core algorithms for FEP calculations including atom mapping and dual topology construction.
"""

from .mapping import Atom, AtomMapping, DistanceAtomMapper
from .dual_topology import HybridAtom, DualTopologyBuilder

__all__ = [
    "Atom",
    "AtomMapping",
    "DistanceAtomMapper",
    "HybridAtom",
    "DualTopologyBuilder",
]
