#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Utils - Utility modules for PRISM Builder
"""

from .environment import GromacsEnvironment
from .config import ConfigurationManager
from .mdp import MDPGenerator
from .system import SystemBuilder
from .ligand import (
    identify_ligand_residue, get_ligand_selection_string,
    validate_ligand_residue, get_ligand_center_atom,
    list_ligand_candidates, get_ligand_info
)
from .topology import (
    detect_protein_chains, detect_nucleic_chains, detect_all_chains,
    get_primary_protein_chain
)

__all__ = [
    'GromacsEnvironment',
    'ConfigurationManager',
    'MDPGenerator',
    'SystemBuilder',
    # Ligand detection utilities
    'identify_ligand_residue',
    'get_ligand_selection_string',
    'validate_ligand_residue',
    'get_ligand_center_atom',
    'list_ligand_candidates',
    'get_ligand_info',
    # Topology analysis utilities
    'detect_protein_chains',
    'detect_nucleic_chains',
    'detect_all_chains',
    'get_primary_protein_chain',
]