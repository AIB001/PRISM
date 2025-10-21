#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM MM/PBSA Module

This module provides MM/PBSA (Molecular Mechanics/Poisson-Boltzmann Surface Area)
calculations for protein-ligand binding free energy estimation.

Two modes are supported:
1. Single-frame: Calculate MM/PBSA for a single structure (e.g., docking pose)
2. Trajectory: Calculate MM/PBSA along MD trajectory (standard workflow)
"""

from .single_frame import SingleFrameMMPBSA

__all__ = ['SingleFrameMMPBSA']
