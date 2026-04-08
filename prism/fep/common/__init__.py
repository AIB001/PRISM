#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Shared infrastructure for the FEP package.

This package hosts cross-cutting utilities that are reused across the FEP
domain, modeling, workflow, and visualization layers.
"""

from prism.fep.common.coordinates import (
    ANGSTROM_PER_NM,
    angstrom_to_nm,
    angstrom_xyz_to_nm,
    nm_to_angstrom,
    nm_xyz_to_angstrom,
)
from prism.fep.common.config import FEPConfig, read_fep_config, write_fep_config
from prism.fep.common.io import (
    read_hybrid_topology_dual,
    read_ligand_from_prism,
    read_mol2_atoms,
    read_rtf_for_fep,
    write_ligand_to_pdb,
)
from prism.fep.common.naming import (
    FEPSystemNamer,
    generate_fep_system_name,
    parse_fep_system_name,
    validate_fep_system_name,
)
from prism.fep.common.utils import calculate_distance, format_charge, get_connected_atoms, validate_atoms

__all__ = [
    "ANGSTROM_PER_NM",
    "angstrom_to_nm",
    "angstrom_xyz_to_nm",
    "nm_to_angstrom",
    "nm_xyz_to_angstrom",
    "FEPConfig",
    "read_fep_config",
    "write_fep_config",
    "read_hybrid_topology_dual",
    "read_ligand_from_prism",
    "read_mol2_atoms",
    "read_rtf_for_fep",
    "write_ligand_to_pdb",
    "FEPSystemNamer",
    "generate_fep_system_name",
    "parse_fep_system_name",
    "validate_fep_system_name",
    "calculate_distance",
    "format_charge",
    "get_connected_atoms",
    "validate_atoms",
]
