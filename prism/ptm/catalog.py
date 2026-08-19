#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Catalog of supported post-translational modifications (PTMs).

Each PTM is keyed by its wwPDB Chemical Component Dictionary (CCD) 3-letter
code (the code that normally appears in the ``HETATM`` records of a deposited
structure). For every modification we record how to realise it in an all-atom
GROMACS ``pdb2gmx`` pipeline for both supported force-field families:

* **CHARMM36** — the bundled ``charmm36-jul2022.ff`` already contains the PTM
  residues in ``aminoacids.rtp`` (verified: SEP/SP1/SP2, TPO/THP1/THP2,
  PTR/TP1/TP2, MLZ, MLY, M3L, ALY, 2MR). SEP/TPO/PTR and SP1/THP1/TP1 are
  duplicate monoanionic entries that differ only in phosphate atom naming, so
  the catalog emits SEP/TPO/PTR for charge -1 and records SP1/THP1/TP1 as
  recognised equivalents. The only requirements are that the residue is
  registered as ``Protein`` in a ``residuetypes.dat`` visible to ``pdb2gmx`` and
  that the heavy-atom names match the ``.rtp`` entry.
* **AMBER** — phosphorylation uses the ``phosaa14SB`` / ``phosaa19SB`` parameter
  sets shipped with AmberTools (``leaprc.phosaa19SB``); the modified residue is
  built with ``tleap`` and converted to GROMACS with ParmEd/ACPYPE through
  PRISM's existing ``amb2gmx`` path.

This module is the single source of truth for PTM residue names. Any other
module that needs to recognise a modified amino acid (thermostat-group
classification, metal-vs-protein discrimination, position-restraint selections)
must derive its table from :func:`ptm_resnames` rather than re-listing the
codes, exactly as ligand-backend metadata is derived from
``prism.forcefield.registry`` rather than re-listed per caller.

References (see PRISM_Capability_Expansion_Roadmap.md for full citations):
    - Raguette et al. JCTC 2024 (phosaa14SB/phosaa19SB)
    - Croitoru et al. JCTC 2021 (CHARMM36 non-standard amino acids)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Set, Tuple


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
    #: CHARMM .rtp entries that describe the same chemical state as one of the
    #: names above under a different atom-naming scheme. They are never emitted
    #: by PRISM, but a structure prepared elsewhere may already carry them, so
    #: they belong in the residue-name table exported by :func:`ptm_resnames`.
    charmm_equivalent_resnames: Tuple[str, ...] = ()

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
    #: CHARMM-specific aliases keyed by the selected .rtp residue name.
    #: Charge variants can use atom names that differ from the CCD entry.
    charmm_atom_aliases: Dict[str, Dict[str, str]] = field(default_factory=dict)
    #: Modification-specific heavy atoms required by each automated target
    #: residue. These signatures prevent a plain parent residue from being
    #: relabelled as a PTM without the atoms that make up the modification.
    required_heavy_atoms: Dict[str, Tuple[str, ...]] = field(default_factory=dict)

    #: Whether validated parameters exist. If False, the modification is only a
    #: structural placeholder and requires user-supplied parameters.
    validated: bool = True
    note: str = ""

    def charmm_name_for_charge(self, charge: Optional[int]) -> Optional[str]:
        if charge is None:
            return self.charmm_resname
        if charge in self.charmm_charge_variants:
            return self.charmm_charge_variants[charge]
        if charge == self.default_charge:
            return self.charmm_resname
        return None

    def amber_name_for_charge(self, charge: Optional[int]) -> Optional[str]:
        if charge is None:
            return self.amber_resname
        if charge in self.amber_charge_variants:
            return self.amber_charge_variants[charge]
        if charge == self.default_charge:
            return self.amber_resname
        return None

    def atom_aliases_for(self, family: str, resname: str) -> Dict[str, str]:
        """Return atom aliases for the selected force-field residue."""
        aliases = dict(self.atom_aliases)
        if family == "charmm":
            aliases.update(self.charmm_atom_aliases.get(resname, {}))
        return aliases

    def required_heavy_atoms_for(self, resname: str) -> Tuple[str, ...]:
        """Return the modification-specific heavy atoms for ``resname``."""
        return self.required_heavy_atoms.get(resname, ())

    def resnames(self, ff_family: Optional[str] = None) -> FrozenSet[str]:
        """Residue names this modification can carry in a built system.

        Includes the CCD code, every force-field residue name and every
        protonation-state variant, plus the CHARMM ``.rtp`` entries that are
        equivalent to one of them. ``ff_family`` is matched loosely against
        ``"amber"`` / ``"charmm"``; ``None`` returns the union. The CCD code is
        always included: a deposited structure carries it whichever force field
        is selected afterwards. Use :func:`supported_codes` instead when the
        question is which PTMs a family can actually build.
        """
        fam = (ff_family or "").lower()
        want_charmm = not fam or "charmm" in fam
        # Non-CHARMM families (amber, opls) all resolve through the AMBER
        # names, matching how PTMStager classifies a force field.
        want_amber = "charmm" not in fam

        names = {self.code}
        if want_charmm:
            names.add(self.charmm_resname)
            names.update(self.charmm_charge_variants.values())
            names.update(self.charmm_equivalent_resnames)
        if want_amber:
            names.add(self.amber_resname)
            names.update(self.amber_charge_variants.values())
        return frozenset(name.strip().upper() for name in names if name)

    def known_source_resnames(self) -> FrozenSet[str]:
        """Residue names that may legitimately denote this PTM in an input PDB.

        This is :meth:`resnames` plus the unmodified parent residue, which is
        what a deposited structure carries when the modification is described
        only by the requested PTM code.
        """
        return self.resnames() | {self.parent.strip().upper()}


# --------------------------------------------------------------------------- #
# The catalog                                                                  #
# --------------------------------------------------------------------------- #

_CATALOG: Tuple[PTMDef, ...] = (
    # ---- Phosphorylation --------------------------------------------------
    PTMDef(
        code="SEP", name="Phosphoserine", category="phosphorylation", parent="SER",
        default_charge=-2,
        # SEP is monoanionic in the bundled .rtp; SP2 is the dianion.
        # SP1 is the same monoanion written with OT/HT instead of O3P/H3T.
        charmm_resname="SP2", charmm_charge_variants={-2: "SP2", -1: "SEP"},
        charmm_equivalent_resnames=("SP1",),
        charmm_atom_aliases={
            "SP2": {"O3P": "OT"},
            "SEP": {"OT": "O3P", "HT": "H3T"},
        },
        required_heavy_atoms={
            "SP2": ("P", "O1P", "O2P", "OT"),
            "SEP": ("P", "O1P", "O2P", "O3P"),
        },
        amber_resname="SEP", amber_leaprc="leaprc.phosaa19SB",
        amber_charge_variants={-2: "SEP", -1: "S1P"},
    ),
    PTMDef(
        code="TPO", name="Phosphothreonine", category="phosphorylation", parent="THR",
        default_charge=-2,
        # TPO is monoanionic. THP1/THP2 are phosphothreonine variants;
        # TP1/TP2 in this force field are phosphotyrosine variants.
        # THP1 is the same monoanion as TPO written with OT/HT.
        charmm_resname="THP2", charmm_charge_variants={-2: "THP2", -1: "TPO"},
        charmm_equivalent_resnames=("THP1",),
        charmm_atom_aliases={
            "THP2": {"O3P": "OT"},
            "TPO": {"OT": "O3P", "HT": "H3T"},
        },
        required_heavy_atoms={
            "THP2": ("P", "O1P", "O2P", "OT"),
            "TPO": ("P", "O1P", "O2P", "O3P"),
        },
        amber_resname="TPO", amber_leaprc="leaprc.phosaa19SB",
        amber_charge_variants={-2: "TPO", -1: "T1P"},
    ),
    PTMDef(
        code="PTR", name="Phosphotyrosine", category="phosphorylation", parent="TYR",
        default_charge=-2,
        # PTR is monoanionic in the bundled .rtp; TP2 is the dianion.
        # TP1 is the same monoanion written with P1/O2/H2 instead of P/O3P/H3T.
        charmm_resname="TP2", charmm_charge_variants={-2: "TP2", -1: "PTR"},
        charmm_equivalent_resnames=("TP1",),
        charmm_atom_aliases={
            "TP2": {"P": "P1", "O3P": "O2", "O1P": "O3", "O2P": "O4"},
            "PTR": {"P1": "P", "O2": "O3P", "H2": "H3T", "O3": "O1P", "O4": "O2P"},
        },
        required_heavy_atoms={
            "TP2": ("P1", "O2", "O3", "O4"),
            "PTR": ("P", "O1P", "O2P", "O3P"),
        },
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
        required_heavy_atoms={"MLZ": ("CM",)},
        note="Monomethyl-lysine (CHARMM36 / Croitoru 2021).",
    ),
    PTMDef(
        code="MLY", name="N6,N6-dimethyl-lysine", category="methylation", parent="LYS",
        default_charge=1, charmm_resname="MLY",
        required_heavy_atoms={"MLY": ("CH1", "CH2")},
    ),
    PTMDef(
        code="M3L", name="N6,N6,N6-trimethyl-lysine", category="methylation", parent="LYS",
        default_charge=1, charmm_resname="M3L",
        required_heavy_atoms={"M3L": ("CM1", "CM2", "CM3")},
    ),
    # ---- Lysine acetylation ----------------------------------------------
    PTMDef(
        code="ALY", name="N6-acetyl-lysine", category="acetylation", parent="LYS",
        default_charge=0, charmm_resname="ALY",
        required_heavy_atoms={"ALY": ("CH", "OH", "CH3")},
    ),
    # ---- Arginine methylation --------------------------------------------
    PTMDef(
        code="2MR", name="Symmetric dimethyl-arginine", category="methylation", parent="ARG",
        default_charge=1, charmm_resname="2MR",
        required_heavy_atoms={"2MR": ("CQ1", "CQ2")},
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


def ptm_resnames(
    ff_family: Optional[str] = None,
    validated_only: bool = True,
    include_parents: bool = False,
) -> FrozenSet[str]:
    """Every residue name a catalogued PTM can appear under.

    This is the table other PRISM modules must use instead of hand-maintaining
    their own list of modified-residue names (thermostat-group classification,
    metal-vs-protein discrimination, position-restraint selections). Keeping it
    here means adding a PTM to :data:`_CATALOG` updates every consumer.

    ``ff_family`` is matched loosely against ``"amber"`` / ``"charmm"``;
    ``None`` returns the union over both families. Unvalidated placeholders are
    excluded by default because staging refuses to build them.
    ``include_parents`` additionally returns the unmodified parent residues,
    which is what a deposited structure carries before staging.
    """
    names: Set[str] = set()
    for definition in _CATALOG:
        if validated_only and not definition.validated:
            continue
        names |= definition.resnames(ff_family)
        if include_parents:
            names.add(definition.parent.strip().upper())
    return frozenset(names)


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
