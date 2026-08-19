#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM membrane module — all-atom protein-in-bilayer system building for GROMACS.

Pipeline (Tier-A, AMBER-consistent): orient (OPM/PPM/preoriented) ->
PACKMOL-Memgen bilayer (Lipid21) -> ParmEd AMBER->GROMACS -> semiisotropic
membrane MDP protocol -> ``localrun.sh`` driver script. All external-tool steps
are runtime-gated with an actionable capability report.

Typical use::

    from prism.membrane import MembraneConfig, MembraneBuilder, parse_membrane_cli
    cfg = parse_membrane_cli(membrane=True, lipid=["POPC"], lipid_ratio=None,
                             orient="opm", pdb_id="4xb4", lipid_ff="lipid21")
    builder = MembraneBuilder(cfg)
    print(builder.capability_report())          # check prerequisites
    result = builder.build("protein.pdb", "out") # build (gated)

The names re-exported below are the supported surface: everything a consumer
outside this package needs, including ``MEMBRANE_SYSTEM_DIRNAME`` so nobody has
to hardcode the output directory or deep-import ``prism.membrane.builder`` for
it. The PACKMOL-Memgen backend itself (command construction, invocation, log
parsing) is deliberately *not* re-exported piecemeal — reach it through
``prism.membrane.packmol_memgen``, whose names only make sense qualified by that
module.

Every module-level import in this package is stdlib; the one third-party
dependency (ParmEd) is imported inside the function that needs it. That is what
makes ``import prism.membrane`` cheap and side-effect free, and therefore safe
to do eagerly from ``prism/__init__.py``. Keep it that way.
"""

from .config import (
    COMMON_LIPIDS,
    DEFAULT_PRODUCTION_NS,
    LIPID_FF_FAMILY,
    MEMBRANE_TIMESTEP_PS,
    MembraneConfig,
    parse_membrane_cli,
    parse_membrane_yaml,
)
from .membrane_mdp import write_membrane_mdps
from .orient import (
    fetch_opm,
    membrane_alignment_report,
    orient_protein,
    strip_opm_dummy_atoms,
    validate_membrane_alignment,
)
from .packmol_memgen import is_available as packmol_memgen_available
from .runscript import RUNSCRIPT_NAME, write_membrane_runscript
from .builder import MEMBRANE_SYSTEM_DIRNAME, MembraneBuilder, MembraneBuildResult

__all__ = [
    "MembraneConfig",
    "parse_membrane_cli",
    "parse_membrane_yaml",
    "COMMON_LIPIDS",
    "DEFAULT_PRODUCTION_NS",
    "MEMBRANE_TIMESTEP_PS",
    "LIPID_FF_FAMILY",
    "write_membrane_mdps",
    "orient_protein",
    "fetch_opm",
    "strip_opm_dummy_atoms",
    "validate_membrane_alignment",
    "membrane_alignment_report",
    "packmol_memgen_available",
    "RUNSCRIPT_NAME",
    "write_membrane_runscript",
    "MEMBRANE_SYSTEM_DIRNAME",
    "MembraneBuilder",
    "MembraneBuildResult",
]
