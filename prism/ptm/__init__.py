#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PTM module — post-translational modification staging for GROMACS.

This package adds modified amino acids (phosphorylation, methylation,
acetylation, ...) and disulfide-bond handling to PRISM's ``pdb2gmx``-based
build pipeline. It supports both force-field families PRISM ships:

* **CHARMM36** — the bundled ``charmm36-jul2022.ff`` already contains the PTM
  residues (SEP/TPO/PTR/SP1/SP2/TP1/TP2/MLZ/MLY/M3L/ALY/2MR); staging registers
  them as ``Protein`` and normalises residue/atom names so ``pdb2gmx`` builds
  them directly — fully local, no external tools.
* **AMBER** — phosphorylation via ``tleap`` + ``phosaa14SB``/``phosaa19SB`` →
  ParmEd/ACPYPE (reuses PRISM's ``amb2gmx`` conversion path).

Disulfides are native to ``pdb2gmx`` and require no parameters.

Typical use::

    from prism.ptm import PTMConfig, PTMStager, parse_ptm_cli
    cfg = parse_ptm_cli(["A:145:SEP", "A:32:TPO:-1"], ssbond="auto", phospho_ff="phosaa19SB")
    stager = PTMStager(cfg, forcefield_name="charmm36-jul2022", forcefield_dir=ff_dir)
    result = stager.apply("protein_clean.pdb", "protein_ptm.pdb", work_dir=model_dir)
"""

from .catalog import PTMDef, get_ptm, is_known_ptm, iter_ptm_defs, supported_codes
from .disulfides import Disulfide, detect_disulfides, bonded_cys_resname
from .spec import PTMConfig, PTMRequest, parse_ptm_cli, parse_ptm_yaml
from .stager import PTMStager, PTMResult, promote_requested_ptm_records

__all__ = [
    "PTMDef",
    "get_ptm",
    "is_known_ptm",
    "iter_ptm_defs",
    "supported_codes",
    "Disulfide",
    "detect_disulfides",
    "bonded_cys_resname",
    "PTMConfig",
    "PTMRequest",
    "parse_ptm_cli",
    "parse_ptm_yaml",
    "PTMStager",
    "PTMResult",
    "promote_requested_ptm_records",
]
