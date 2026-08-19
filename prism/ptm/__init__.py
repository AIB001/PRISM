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
  them directly — fully local, no external tools. This is the fully automated
  path for every catalogued modification.
* **AMBER** — methylation/acetylation stage exactly as above. Phosphorylation
  is the documented ``tleap`` + ``phosaa14SB``/``phosaa19SB`` -> ParmEd/ACPYPE
  route (PRISM's existing ``amb2gmx`` conversion path), but it is *not* wired
  into the default ``pdb2gmx`` build: :class:`PTMStager` gates on it and raises
  with an actionable message pointing at ``charmm36-jul2022``.

Disulfides are native to ``pdb2gmx`` and require no parameters.

Typical use::

    from prism.ptm import PTMConfig, PTMStager, parse_ptm_cli
    cfg = parse_ptm_cli(["A:145:SEP", "A:32:TPO:-1"], ssbond="auto", phospho_ff="phosaa19SB")
    stager = PTMStager(cfg, forcefield_name="charmm36-jul2022", forcefield_dir=ff_dir)
    result = stager.apply("protein_clean.pdb", "protein_ptm.pdb", work_dir=model_dir)

The names re-exported below are the supported surface: everything a consumer
outside this package needs in order to describe, validate and stage a
modification, plus the result types they get back (``PTMResult`` carries
``Disulfide`` objects, which in turn carry ``CysSG`` atoms, so all three are
exported rather than left as deep imports). The staging internals — force-field
family classification, PDB column surgery, ``residuetypes.dat`` discovery — are
deliberately *not* re-exported piecemeal; reach them through
``prism.ptm.stager``, whose names only make sense qualified by that module.

Every import in this package is stdlib; there is no third-party dependency at
all, not even a lazily imported one. That is what makes ``import prism.ptm``
cheap and side-effect free, and therefore safe to do eagerly from a CLI
argument handler or from the middle of a build. Keep it that way: the external
tools this module targets (``pdb2gmx``, ``tleap``) are invoked as subprocesses
and probed with ``shutil.which``, never imported.
"""

from .catalog import PTMDef, get_ptm, is_known_ptm, iter_ptm_defs, supported_codes
from .disulfides import (
    CysSG,
    Disulfide,
    apply_disulfide_renaming,
    bonded_cys_resname,
    detect_disulfides,
)
from .spec import PTMConfig, PTMRequest, parse_ptm_cli, parse_ptm_yaml
from .stager import PTMResult, PTMStager, promote_requested_ptm_records

__all__ = [
    "PTMDef",
    "get_ptm",
    "is_known_ptm",
    "iter_ptm_defs",
    "supported_codes",
    "CysSG",
    "Disulfide",
    "apply_disulfide_renaming",
    "bonded_cys_resname",
    "detect_disulfides",
    "PTMConfig",
    "PTMRequest",
    "parse_ptm_cli",
    "parse_ptm_yaml",
    "PTMResult",
    "PTMStager",
    "promote_requested_ptm_records",
]
