"""
REST2 (Replica Exchange with Solute Tempering v2) module for PRISM.

Generates REST2 replica exchange setups from standard GROMACS MD directories
produced by the PRISM builder workflow.
"""

from .rest2_workflow import REST2Workflow
from .topology import (
    parse_gro,
    find_pocket_residues,
    merge_topology,
    parse_molecules_section,
    count_atoms_per_moltype,
    mark_solute_atoms,
)
from .partial_tempering import partial_tempering
from .builder import build_rest2, calculate_temperature_ladder

__all__ = [
    "REST2Workflow",
    "parse_gro",
    "find_pocket_residues",
    "merge_topology",
    "parse_molecules_section",
    "count_atoms_per_moltype",
    "mark_solute_atoms",
    "partial_tempering",
    "build_rest2",
    "calculate_temperature_ladder",
]
