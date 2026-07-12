#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM membrane module — all-atom protein-in-bilayer system building for GROMACS.

Pipeline (Tier-A, AMBER-consistent): orient (OPM/PPM/preoriented) ->
PACKMOL-Memgen bilayer (Lipid21) -> ParmEd AMBER->GROMACS -> semiisotropic
membrane MDP protocol. All external-tool steps are runtime-gated with an
actionable capability report.

Typical use::

    from prism.membrane import MembraneConfig, MembraneBuilder, parse_membrane_cli
    cfg = parse_membrane_cli(membrane=True, lipid=["POPC"], lipid_ratio=None,
                             orient="opm", pdb_id="4xb4", lipid_ff="lipid21")
    builder = MembraneBuilder(cfg)
    print(builder.capability_report())          # check prerequisites
    result = builder.build("protein.pdb", "out") # build (gated)
"""

from .config import MembraneConfig, parse_membrane_cli, parse_membrane_yaml, COMMON_LIPIDS, LIPID_FF_FAMILY
from .membrane_mdp import write_membrane_mdps
from .orient import orient_protein, fetch_opm, strip_opm_dummy_atoms
from .packmol_memgen import build_command, is_available as packmol_memgen_available
from .builder import MembraneBuilder, MembraneBuildResult

__all__ = [
    "MembraneConfig",
    "parse_membrane_cli",
    "parse_membrane_yaml",
    "COMMON_LIPIDS",
    "LIPID_FF_FAMILY",
    "write_membrane_mdps",
    "orient_protein",
    "fetch_opm",
    "strip_opm_dummy_atoms",
    "build_command",
    "packmol_memgen_available",
    "MembraneBuilder",
    "MembraneBuildResult",
]
