#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM-FEP Module

Free Energy Perturbation (FEP) module for relative binding free energy calculations
between structurally similar ligands using GROMACS.

This module implements:
- Distance-based atom mapping for hybrid topology construction
- Single-topology GROMACS format output (typeB/chargeB encoding)
- FEP simulation setup and analysis tools
"""

from prism.fep.core.mapping import Atom, AtomMapping, DistanceAtomMapper
from prism.fep.core.hybrid_topology import HybridAtom, HybridTopologyBuilder, LigandTopologyInput
from prism.fep.gromacs.itp_builder import ITPBuilder
from prism.fep.gromacs.mdp_templates import FEP_PROD_MDP, write_fep_mdps
from prism.fep.modeling import FEPScaffoldBuilder
from prism.fep.analysis import XVGParser, FEstimator
from prism.fep.io import read_ligand_from_prism, read_mol2_atoms, read_rtf_for_fep
from prism.fep.config import FEPConfig, read_fep_config, write_fep_config

__version__ = "0.1.0"
__all__ = [
    # Core data structures
    "Atom",
    "AtomMapping",
    "DistanceAtomMapper",
    "HybridAtom",
    "HybridTopologyBuilder",
    "LigandTopologyInput",
    # File I/O
    "read_ligand_from_prism",
    "read_mol2_atoms",
    "read_rtf_for_fep",
    # Configuration
    "FEPConfig",
    "read_fep_config",
    "write_fep_config",
    # GROMACS output
    "ITPBuilder",
    "FEP_PROD_MDP",
    "write_fep_mdps",
    "FEPScaffoldBuilder",
    # Analysis
    "XVGParser",
    "FEstimator",
]
