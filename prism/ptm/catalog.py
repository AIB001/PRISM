#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Catalog of supported post-translational modifications (PTMs).

Each PTM is keyed by its wwPDB Chemical Component Dictionary (CCD) 3-letter
code (the code that normally appears in the ``HETATM`` records of a deposited
structure). For every modification we record how to realise it in an all-atom
GROMACS ``pdb2gmx`` pipeline for both supported force-field families:

* **CHARMM36** — the bundled ``charmm36-jul2022.ff`` already contains the PTM
  residues in ``aminoacids.rtp`` (verified: SEP, TPO, PTR, SP1/SP2, TP1/TP2,
  MLZ, MLY, M3L, ALY, 2MR). The only requirements are that the residue is
  registered as ``Protein`` in a ``residuetypes.dat`` visible to ``pdb2gmx`` and
  that the heavy-atom names match the ``.rtp`` entry.
* **AMBER** — phosphorylation uses the ``phosaa14SB`` / ``phosaa19SB`` parameter
  sets shipped with AmberTools (``leaprc.phosaa19SB``); the modified residue is
  built with ``tleap`` and converted to GROMACS with ParmEd/ACPYPE through
  PRISM's existing ``amb2gmx`` path.

References (see PRISM_Capability_Expansion_Roadmap.md for full citations):
    - Raguette et al. JCTC 2024 (phosaa14SB/phosaa19SB)
    - Croitoru et al. JCTC 2021 (CHARMM36 non-standard amino acids)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class PTMDef:
    """Definition of a single post-translational modification."""

    code: str                       # CCD code, e.g. "SEP"
    name: str                       # human-readable name
    category: str                   # phosphorylation | methylation | acetylation | ...
    parent: str                     # parent standard residue (3-letter)
    default_charge: int             # net side-chain charge in the default state

    #: CHARMM36 .rtp residue name (None => not available in CHARMM36 path).
    charmm_resname: Optional[str] = None
    #: Alternative CHARMM rtp names for other protonation states, by charge.
    charmm_charge_variants: Dict[int, str] = field(default_factory=dict)

    #: AMBER tleap residue name (None => not available in the AMBER phosaa path).
    amber_resname: Optional[str] = None
    #: AMBER leaprc that must be sourced to parameterize this residue.
    amber_leaprc: Optional[str] = None
    #: Alternative AMBER residue names by charge.
    amber_charge_variants: Dict[int, str] = field(default_factory=dict)

    #: Optional atom-name normalization {input_name: rtp_name}. pdb2gmx matches
    #: residues by atom name, so deposited names that differ from the .rtp entry
    #: must be renamed first.
    atom_aliases: Dict[str, str] = field(default_factory=dict)

    #: Whether validated parameters exist. If False, the modification is only a
    #: structural placeholder and requires user-supplied parameters.
    validated: bool = True
    note: str = ""

    def charmm_name_for_charge(self, charge: Optional[int]) -> Optional[str]:
        if charge is not None and charge in self.charmm_charge_variants:
            return self.charmm_charge_variants[charge]
        return self.charmm_resname

    def amber_name_for_charge(self, charge: Optional[int]) -> Optional[str]:
        if charge is not None and charge in self.amber_charge_variants:
            return self.amber_charge_variants[charge]
        return self.amber_resname


# --------------------------------------------------------------------------- #
# The catalog                                                                  #
# --------------------------------------------------------------------------- #

_CATALOG: Tuple[PTMDef, ...] = (
    # ---- Phosphorylation --------------------------------------------------
    PTMDef(
        code="SEP", name="Phosphoserine", category="phosphorylation", parent="SER",
        default_charge=-2,
        charmm_resname="SEP", charmm_charge_variants={-2: "SEP", -1: "SP1"},
        amber_resname="SEP", amber_leaprc="leaprc.phosaa19SB",
        amber_charge_variants={-2: "SEP", -1: "S1P"},
    ),
    PTMDef(
        code="TPO", name="Phosphothreonine", category="phosphorylation", parent="THR",
        default_charge=-2,
        charmm_resname="TPO", charmm_charge_variants={-2: "TPO", -1: "TP1"},
        amber_resname="TPO", amber_leaprc="leaprc.phosaa19SB",
        amber_charge_variants={-2: "TPO", -1: "T1P"},
    ),
    PTMDef(
        code="PTR", name="Phosphotyrosine", category="phosphorylation", parent="TYR",
        default_charge=-2,
        charmm_resname="PTR", charmm_charge_variants={-2: "PTR"},
        amber_resname="PTR", amber_leaprc="leaprc.phosaa19SB",
        amber_charge_variants={-2: "PTR", -1: "Y1P"},
    ),
    PTMDef(
        code="PTM", name="Phosphohistidine (3-pHis)", category="phosphorylation", parent="HIS",
        default_charge=-1,
        charmm_resname=None,
        amber_resname=None, amber_leaprc="leaprc.phosaa19SB",
        validated=False,
        note="Phospho-His naming is version-dependent; verify the leaprc/rtp residue name before use.",
    ),
    # ---- Lysine methylation ----------------------------------------------
    PTMDef(
        code="MLZ", name="N6-methyl-lysine", category="methylation", parent="LYS",
        default_charge=1, charmm_resname="MLZ",
        note="Monomethyl-lysine (CHARMM36 / Croitoru 2021).",
    ),
    PTMDef(
        code="MLY", name="N6,N6-dimethyl-lysine", category="methylation", parent="LYS",
        default_charge=1, charmm_resname="MLY",
    ),
    PTMDef(
        code="M3L", name="N6,N6,N6-trimethyl-lysine", category="methylation", parent="LYS",
        default_charge=1, charmm_resname="M3L",
    ),
    # ---- Lysine acetylation ----------------------------------------------
    PTMDef(
        code="ALY", name="N6-acetyl-lysine", category="acetylation", parent="LYS",
        default_charge=0, charmm_resname="ALY",
    ),
    # ---- Arginine methylation --------------------------------------------
    PTMDef(
        code="2MR", name="Asymmetric dimethyl-arginine", category="methylation", parent="ARG",
        default_charge=1, charmm_resname="2MR",
    ),
    # ---- Modifications without canonical validated parameters ------------
    PTMDef(
        code="HYP", name="4-hydroxyproline", category="hydroxylation", parent="PRO",
        default_charge=0, charmm_resname=None, validated=False,
        note="No canonical AMBER/CHARMM36 parameters in the bundled FFs; "
        "use Vienna-PTM (GROMOS) or custom CGenFF/ffTK parameters.",
    ),
    PTMDef(
        code="CIR", name="Citrulline", category="deimination", parent="ARG",
        default_charge=0, charmm_resname=None, validated=False,
        note="No canonical validated parameter set; requires de novo CGenFF+ffTK "
        "(CHARMM) or GAFF+QM (AMBER) parameterization.",
    ),
)

_BY_CODE: Dict[str, PTMDef] = {d.code: d for d in _CATALOG}


def iter_ptm_defs() -> Tuple[PTMDef, ...]:
    """Return all catalogued PTM definitions."""
    return _CATALOG


def get_ptm(code: str) -> Optional[PTMDef]:
    """Look up a PTM definition by CCD code (case-insensitive)."""
    return _BY_CODE.get(str(code).strip().upper())


def is_known_ptm(code: str) -> bool:
    return str(code).strip().upper() in _BY_CODE


def supported_codes(ff_family: Optional[str] = None) -> List[str]:
    """List PTM codes supported (optionally for a given force-field family).

    ``ff_family`` is matched loosely against ``"amber"`` / ``"charmm"``.
    """
    fam = (ff_family or "").lower()
    out = []
    for d in _CATALOG:
        if not d.validated:
            continue
        if "charmm" in fam and not d.charmm_resname:
            continue
        if "amber" in fam and not d.amber_resname:
            continue
        out.append(d.code)
    return out
