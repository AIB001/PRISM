"""
Common utilities for PRISM forcefield module.

This module contains shared functionality used across multiple forcefield
generators to avoid code duplication.

Modules:
    io: Unified coordinate file I/O (GRO, MOL2, PDB)
    units: Unit conversion functions (Å ↔ nm)
    itp_parser: Unified ITP file parser
"""

from .io import (
    read_gro_coordinates,
    read_mol2_coordinates,
    read_pdb_coordinates,
    write_gro_coordinates,
    detect_coordinate_file_format,
)
from .units import (
    angstrom_to_nm,
    nm_to_angstrom,
    angstrom_xyz_to_nm,
    nm_xyz_to_angstrom,
)
from .itp_parser import ITPParser

__all__ = [
    # Coordinate I/O
    "read_gro_coordinates",
    "read_mol2_coordinates",
    "read_pdb_coordinates",
    "write_gro_coordinates",
    "detect_coordinate_file_format",
    # Unit conversion
    "angstrom_to_nm",
    "nm_to_angstrom",
    "angstrom_xyz_to_nm",
    "nm_xyz_to_angstrom",
    # ITP parsing
    "ITPParser",
    "extract_section_from_content",
]

# Export static method directly for convenience
extract_section_from_content = ITPParser.extract_section_from_content
