#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM-FEP Core Module

Core algorithms for FEP calculations including atom mapping and dual topology construction.
"""

from .mapping import Atom, AtomMapping, DistanceAtomMapper
from .hybrid_topology import HybridAtom, HybridTopologyBuilder

__all__ = [
    "Atom",
    "AtomMapping",
    "DistanceAtomMapper",
    "HybridAtom",
    "HybridTopologyBuilder",
]
