#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Membrane-build configuration.

Holds the lipid composition, orientation source, and bilayer/solvent box
settings for an all-atom protein-in-membrane build. Parsed from CLI flags or a
YAML ``membrane:`` section.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import MISSING, dataclass, field, fields
from numbers import Integral, Real
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


# Common lipids understood by PACKMOL-Memgen (AMBER Lipid21 names). This is a
# convenience whitelist for validation/help; PACKMOL-Memgen supports many more.
COMMON_LIPIDS = {
    "POPC", "POPE", "POPS", "POPG", "POPA", "POPI",
    "DOPC", "DOPE", "DOPS", "DPPC", "DMPC", "DLPC",
    "CHL1", "PSM",  # cholesterol / sphingomyelin
}

# PACKMOL-Memgen keys its lipid table by exactly one name per molecule, so a
# common synonym is not a "maybe valid" name that deserves a warning -- it is a
# hard failure inside the tool (memgen.parm has no ``CHOL`` entry; cholesterol
# is ``CHL1``). Translate the synonyms PRISM already recognises elsewhere.
#
# ``CHL`` is here for a second reason.  Cholesterol has two spellings in this
# codebase and they are both correct in their own vocabulary: ``CHL1`` is the
# *molecule* name PACKMOL-Memgen keys on, and ``CHL`` is the Lipid21 *residue*
# name, which is what the bundled patch manifest and
# :data:`prism.membrane.patch.LIPID21_RESIDUE_PATTERNS` use.  Config speaks the
# PACKMOL-Memgen vocabulary, so the residue spelling is folded into it on the
# way in -- otherwise the same physical bilayer has two composition keys and
# only one of them finds the patch.
LIPID_NAME_ALIASES = {"CHOL": "CHL1", "CHL": "CHL1"}

#: Shortest PDB coordinate record this package can read: columns 31-54 hold x, y
#: and z, so anything shorter has no coordinates to place and no chain/residue
#: columns to renumber against.  Shared so that the modules which *count*
#: residues and the module which *detects* disulfides in them agree on exactly
#: which records exist -- disagreeing by one record shifts every residue ordinal
#: after it, and a shifted ordinal declares a bond against the wrong residue.
PDB_MIN_COORDINATE_WIDTH = 54


# --------------------------------------------------------------------------- #
# Canonical residue-name vocabulary
# --------------------------------------------------------------------------- #
#
# Every module in this package answers some question of the form "what kind of
# thing is this residue?", and until these tables existed each one carried its
# own hand-maintained answer.  The copies drifted, and the drift was only ever
# discovered by a build that died hours in:
#
#   * ``orient`` knew ``MSE`` (selenomethionine -- routine in crystal
#     structures) while the builder's index writer did not, so an MSE structure
#     oriented, packed, parametrized and solvated and *then* failed thermostat
#     classification, which has no fallback and no partial credit.
#   * ``orient`` did not know ``SEP``/``TPO``/``PTR``, which is exactly what
#     PRISM's own :mod:`prism.ptm` produces, so alignment was measured on an
#     atom set that silently omitted the phosphorylated residues.
#   * no two of the four water-name sets in this package agreed with each other.
#
# This module is their home because it has no intra-package imports, so every
# other module in ``prism.membrane`` can consume it without an import cycle.
#
# These are deliberately *building blocks* rather than one flat "known residues"
# set, because the call sites ask genuinely different questions:
#
#   orient        -- "is this atom part of the solute whose position I should
#                     measure?"  Wants the protein and anything covalently part
#                     of it, and must exclude the surrounding bath.
#   index writer  -- "which thermostat group does this atom belong to?"  Must
#                     partition *every* atom, and fails closed on anything it
#                     cannot place rather than guessing.
#   solvate       -- "is this a solvent molecule I may delete or refuse?"
#                     Wants water and free ions only, never the solute.
#
# A single flat set satisfies none of the three, so each call site composes the
# union it needs from these pieces and documents which question it is asking.

#: The 20 residues every force field names identically.
STANDARD_AMINO_ACID_RESIDUES = frozenset(
    """
    ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
    """.split()
)

#: Protonation / tautomer spellings emitted by AMBER, CHARMM and GROMACS for the
#: residues above.  Same chemistry, different name depending on which pipeline
#: wrote the file -- and membrane inputs reach PRISM from all three.
PROTONATION_VARIANT_RESIDUES = frozenset(
    """
    ASH ASX CYM CYX GLH GLX HID HIE HIP HSD HSE HSP LYN
    """.split()
)

#: Modified residues that are still part of the protein chain.  ``MSE`` is the
#: one that motivated this table (see above).  The phospho- and methyl-lysine
#: names are the ones :mod:`prism.ptm` writes, so PRISM must recognise its own
#: output; ``SEC``/``PYL`` are the two genetically encoded extras.
MODIFIED_AMINO_ACID_RESIDUES = frozenset(
    """
    MSE SEC PYL
    SEP TPO PTR S1P T1P Y1P SP1 SP2 THP1 THP2 TP1 TP2
    MLZ MLY M3L ALY 2MR
    """.split()
)

#: Terminal caps.  tleap writes these for any truncated construct and for
#: ``--nter``/``--cter``, so leaving them out made every capped peptide fail
#: index generation -- after the packing had already run.
TERMINAL_CAP_RESIDUES = frozenset("ACE NME NHE NH2 FOR".split())

#: Everything that is part of a protein chain.  This is the union the index
#: writer uses for its SOLU group and the base of what ``orient`` measures.
PROTEIN_RESIDUE_NAMES = frozenset(
    STANDARD_AMINO_ACID_RESIDUES
    | PROTONATION_VARIANT_RESIDUES
    | MODIFIED_AMINO_ACID_RESIDUES
    | TERMINAL_CAP_RESIDUES
)

#: A polymer residue whose chemical identity the depositor did not know.
#: Deliberately NOT in :data:`PROTEIN_RESIDUE_NAMES`: ``orient`` should still
#: measure it (it is part of the solute, and mismeasuring the solute is the
#: worse error), but the index writer must not silently assign an unidentified
#: molecule to the protein thermostat group.  In practice tleap refuses ``UNK``
#: long before the index is written, so this asymmetry costs nothing.
UNKNOWN_POLYMER_RESIDUE_NAMES = frozenset({"UNK"})

#: Residue names ``orient`` measures as solute.  Its question is "where is the
#: thing I am aligning to the bilayer", so an unidentified polymer residue
#: counts and a free ion does not.
SOLUTE_MEASUREMENT_RESIDUE_NAMES = frozenset(
    PROTEIN_RESIDUE_NAMES | UNKNOWN_POLYMER_RESIDUE_NAMES
)

#: Water, in every spelling ``gmx solvate``, ParmEd, tleap and the PDB use.
WATER_RESIDUE_NAMES = frozenset(
    """
    SOL WAT HOH DOD H2O
    TIP TIP3 TIP3P T3P TIP4 TIP4P T4P TIP5P T5P
    SPC SPCE OPC OPC3
    """.split()
)

#: Monatomic ions that belong to the *bath* -- the ones ``gmx genion`` and
#: PACKMOL-Memgen place to neutralise and salt the box.  Every spelling that can
#: reach a coordinate file is listed, because the configured ion name propagates
#: verbatim into the residue name (``positive_ion: "Li+"`` produces a residue
#: literally named ``LI+``), and the index writer has to classify whatever
#: arrives rather than raise on it.
COUNTER_ION_RESIDUE_NAMES = frozenset(
    """
    NA NA+ SOD SODIUM
    K K+ POT POTASSIUM
    LI LI+ RB RB+ CS CS+
    MG MG2 MG2+ CA CA2 CA2+ CAL ZN ZN2 ZN2+
    CL CL- CLA CHLORIDE
    BR BR- I I- IOD F F- ION
    """.split()
)

#: Transition metals that are normally *structural* -- part of the solute, not
#: the bath.  Kept out of :data:`SOLVENT_RESIDUE_NAMES` on purpose: a buried
#: catalytic iron coupled to the solvent thermostat is a physics error, so the
#: index writer should fail closed and make the user say what it is.
#:
#: The divalent ions above (MG/CA/ZN) are the genuinely ambiguous case -- they
#: are both common counter-ions and common structural cofactors.  They are
#: treated as bath ions here because that is what PRISM has always done and
#: because they are selectable as ``positive_ion``; a structural Zn in a zinc
#: finger is therefore still misgrouped.  That is a pre-existing limitation this
#: table makes visible rather than one it introduces.
STRUCTURAL_METAL_RESIDUE_NAMES = frozenset("FE FE2 FE3 MN CU CU1 CD CO NI HG".split())

#: Water + bath ions: "is this a solvent molecule?"  This is the set the index
#: writer puts in SOLV and the set :mod:`prism.membrane.solvate` uses to refuse
#: to solvate a system that is already wet.
SOLVENT_RESIDUE_NAMES = frozenset(WATER_RESIDUE_NAMES | COUNTER_ION_RESIDUE_NAMES)

#: Lipid molecule names plus the Lipid21/CHARMM *fragment* residues a packed
#: bilayer is actually written with.
LIPID_RESIDUE_NAMES = frozenset(
    """
    POPC POPE POPS POPG POPA POPI DOPC DOPE DOPS DPPC DPPE DMPC DLPC DSPC SDPC
    CHL CHL1 CHOL PSM SSM DHA
    OL PA PC PE PS PGR PI SA ST MY LAL SPM PH-
    """.split()
)

#: Everything that is demonstrably *not* solute.  ``orient`` uses this as the
#: negative test when a structure uses none of the protein names it knows: an
#: unrecognised protein must still be measured rather than waved through, but a
#: water or a lipid must never dilute the measurement.
NON_SOLUTE_RESIDUE_NAMES = frozenset(
    SOLVENT_RESIDUE_NAMES | STRUCTURAL_METAL_RESIDUE_NAMES | LIPID_RESIDUE_NAMES
)

#: Which Lipid21 *fragment* residues each lipid molecule is built from.
#:
#: AMBER Lipid17/21 does not keep the input molecule name: one POPC is written
#: as PA (palmitoyl) + PC (headgroup) + OL (oleoyl).  This is the *unordered
#: membership* form -- "is this residue part of that lipid?" -- which is the
#: only question the index writer and the lipid tally ask.
#:
#: :data:`prism.membrane.patch.LIPID21_RESIDUE_PATTERNS` holds the same
#: chemistry in an *ordered tuple* form because it answers a different question:
#: it has to emit the fragments of one molecule in file order, so it records
#: stoichiometry and position (``DMPC`` is ``("MY", "PC", "MY")`` there -- two
#: myristoyl tails around one headgroup -- and ``{"MY", "PC"}`` here, because a
#: set cannot and should not carry the multiplicity).  The relationship between
#: the two is exactly
#:
#:     set(LIPID21_RESIDUE_PATTERNS[x]) == LIPID_FRAGMENT_RESIDUES[LIPID_NAME_ALIASES.get(x, x)]
#:
#: for every lipid both tables know -- verified for all six of ``patch``'s keys.
#: The alias lookup is load-bearing rather than decorative: ``patch`` keys
#: cholesterol on the Lipid21 *residue* name ``CHL`` and this module keys it on
#: the PACKMOL-Memgen *molecule* name ``CHL1``, which is the same dual spelling
#: :data:`LIPID_NAME_ALIASES` above exists to reconcile.  They are deliberately
#: NOT merged, because
#: collapsing the ordered form would lose the stoichiometry the patch route
#: needs and inflating the set form would imply an order the index writer must
#: not depend on.  ``patch`` covers only the lipids it ships a patch for, so
#: this table is the wider of the two.
LIPID_FRAGMENT_RESIDUES = {
    "POPC": frozenset({"PA", "PC", "OL"}),
    "POPE": frozenset({"PA", "PE", "OL"}),
    "POPS": frozenset({"PA", "PS", "OL"}),
    "POPG": frozenset({"PA", "PGR", "OL"}),
    "POPA": frozenset({"PA", "PH-", "OL"}),
    "DOPC": frozenset({"PC", "OL"}),
    "DOPE": frozenset({"PE", "OL"}),
    "DOPS": frozenset({"PS", "OL"}),
    "DPPC": frozenset({"PA", "PC"}),
    "DMPC": frozenset({"MY", "PC"}),
    "DLPC": frozenset({"LAL", "PC"}),
    "CHL1": frozenset({"CHL"}),
    "CHOL": frozenset({"CHL"}),
    "PSM": frozenset({"PA", "SA", "SPM"}),
}

#: Counter-ion names PRISM can carry all the way to a finished system.
#:
#: This is narrower than the set of names that *parse*, and the gap is measured
#: rather than assumed.  ``prism.membrane.solvent_params`` asks tleap for an ion
#: unit by name; the names it asks for are only real units in AmberTools'
#: ``atomic_ions.lib`` for sodium, potassium and chloride.  For every other ion
#: tleap accepts the script, reports ``Errors = 0``, and silently returns a
#: topology with the ion missing (measured against AmberTools' teLeap: a probe
#: for ``Mg+`` yields a topology whose only residue is ``WAT``).  That is caught
#: downstream by ``assert_solvent_blocks`` -- but only after the bilayer has been
#: packed and parametrized, which is the expensive part of the build.
#:
#: Kept as an advisory rather than a hard config error on purpose: PACKMOL-Memgen
#: has its own ion table and remaps some of these itself, so a wet PACKMOL build
#: may well succeed with an ion the GROMACS solvation route cannot carry.
#: Refusing at config time would block a route that works.
GROMACS_ROUTE_SUPPORTED_IONS = frozenset(
    """
    NA NA+ SOD SODIUM
    K K+ POT POTASSIUM
    CL CL- CLA CHLORIDE
    """.split()
)

#: Ion spellings PRISM recognises at all (mirrors ``solvent_params._ION_UNITS``
#: keys).  A name outside this set is not merely unsupported downstream, it is
#: unrecognised, and that IS worth a hard config error -- it is almost always a
#: typo, and the failure would otherwise be an opaque tleap message.
RECOGNISED_ION_NAMES = frozenset(
    """
    NA NA+ SOD SODIUM K K+ POT POTASSIUM LI LI+ RB RB+ CS CS+
    MG MG2+ CA CA2+ CAL ZN ZN2+
    CL CL- CLA CHLORIDE BR BR- I I- IOD F F-
    """.split()
)


# Lipid force fields and the protein-FF families they pair with.
LIPID_FF_FAMILY = {
    "lipid17": "amber",     # AMBER Lipid17 (legacy PACKMOL-Memgen option)
    "lipid21": "amber",     # AMBER Lipid21 (PRISM default; selected explicitly)
    "slipids": "amber",     # Slipids-2020 (AMBER-compatible)
    "charmm36": "charmm",   # CHARMM36 lipids
}

# ``--parametrize`` in PACKMOL-Memgen produces an AMBER topology.  PRISM can
# therefore select the AMBER lipid force fields directly, but the advertised
# Slipids/CHARMM entries need a separate parameterisation backend and must not
# be silently passed off as Lipid21.
PACKMOL_PARAMETRIZED_LIPID_FFS = {"lipid17", "lipid21"}

# PACKMOL-Memgen does not consume GROMACS force-field directory names.  Keep
# the translation in one place so a membrane request cannot silently report a
# GROMACS force field while PACKMOL-Memgen parameterizes with its defaults.
PACKMOL_PROTEIN_FF_ALIASES = {
    "amber14sb": "ff14SB",
    "ff14sb": "ff14SB",
    "amber19sb": "ff19SB",
    "ff19sb": "ff19SB",
}
PACKMOL_WATER_MODEL_ALIASES = {
    "tip3p": "tip3p",
    "opc": "opc",
}
PACKMOL_FORCEFIELD_PAIRS = {
    ("ff14SB", "tip3p"),
    ("ff19SB", "opc"),
}

#: PACKMOL-Memgen writes an intermediate structure every N GENCAN loops and
#: refuses any iteration budget that is not strictly larger than N.
PACKMOL_WRITEOUT_FREQUENCY = 10

#: Backends for placing the aqueous phase (water + ions).
#:
#: ``packmol``  -- PACKMOL-Memgen packs water and ions together with the lipids.
#:                 The original path; still correct, and still the only one that
#:                 can place water on a bilayer PACKMOL-Memgen also packed.
#: ``gromacs``  -- PACKMOL packs protein+lipids only ("dry"); ``gmx solvate`` and
#:                 ``gmx genion`` add the aqueous phase afterwards.  PACKMOL then
#:                 places ~10^2 molecules instead of ~10^5, which is the whole
#:                 point: on a class-C GPCR homodimer the wet pack never finished.
#: ``auto``     -- the default: pick ``gromacs`` only when the projected water
#:                 count makes the wet pack impractical (see
#:                 ``solvate_auto_threshold``), and ``packmol`` otherwise.
SOLVATE_MODES = {"packmol", "gromacs", "auto"}

#: ``gmx solvate`` reads ``vdwradii.dat`` from its working directory and inserts
#: water wherever the radii allow.  Stock GROMACS gives carbon 0.17 nm, which is
#: far too small for the gaps between freshly packed acyl tails: water lands in
#: the hydrophobic core.  0.375 nm is the value the GROMACS membrane tutorials
#: use for exactly this reason.
#:
#: Inflating carbon is a *mitigation*, never a fix -- measured on a PACKMOL-packed
#: POPC bilayer it is asymptotic (0.17 nm leaves ~1400-2400 core waters, 0.375 nm
#: still leaves ~600-1100, and pushing past 0.5 nm starts starving the bulk before
#: it reaches zero).  Only the geometric core purge reaches PACKMOL's exact 0,
#: which is why ``core_water_purge`` defaults to True and disabling it is warned
#: about rather than merely documented.
DEFAULT_VDWRADII_CARBON_NM = 0.375
#: Stock GROMACS ``vdwradii.dat`` carbon radius.  Anything at or below this does
#: no inflation at all and hands the entire core-water problem to the purge.
GROMACS_DEFAULT_CARBON_RADIUS_NM = 0.17
#: Above this the dead shell carved around every solute carbon starts to
#: under-fill the bulk water noticeably (advisory).
VDWRADII_CARBON_BULK_WARNING_NM = 0.5
#: Above this the box would be so starved of water that the build is wrong
#: rather than merely suboptimal (fatal).
VDWRADII_CARBON_MAX_NM = 1.0

#: Water count above which ``solvate="auto"`` switches to the GROMACS path.
#: Below it the wet pack completes in a tolerable time and there is no reason to
#: leave the route PACKMOL-Memgen does in one run.
#:
#: The number is a real discriminator rather than decoration, and it is worth
#: recording which side each measured build falls on.  mGluR7 (PDB 9OMO, 21 348
#: protein atoms) projects ~158 800 waters -- 5x over -- and takes the GROMACS
#: path; a minimal single-patch box (128 lipids, ~64 x 64 A, 17.5 A of water
#: either side) projects a few thousand and stays on PACKMOL.  So both branches
#: are reachable from realistic inputs, which is the only property that makes a
#: threshold worth having.
#:
#: It is consulted less often than it looks, though: when the bilayer itself
#: comes from a patch, :meth:`MembraneConfig.resolved_solvate_mode` returns
#: ``gromacs`` without measuring anything, because a tiled patch is dry and
#: there is no PACKMOL run to place water in.  This threshold therefore only
#: decides builds whose bilayer PACKMOL-Memgen packed -- i.e. compositions with
#: no bundled patch, which is exactly where the wet pack is still an option.
DEFAULT_SOLVATE_AUTO_THRESHOLD = 30000

#: Backends for building the **bilayer itself**.  This is a different axis from
#: :data:`SOLVATE_MODES`, which chooses who places the *aqueous phase*: this
#: setting is about the lipids, that one is about the water.  Both can vary
#: independently except for the one combination ruled out below.
#:
#: ``packmol`` -- PACKMOL-Memgen packs the lipids around the solute from single
#:                lipid conformers.  The original path, and the only option for
#:                a composition with no bundled patch.
#: ``patch``   -- tile a bundled, pre-equilibrated Lipid21 bilayer patch in xy,
#:                insert the oriented solute at z = 0, and delete the lipids it
#:                overlaps.  Seconds instead of tens of minutes, and it does not
#:                depend on PACKMOL converging: at the physical area per lipid,
#:                PACKMOL's rigid-body packing cannot satisfy its own 2.0 A
#:                tolerance and simply never converges.
#: ``auto``    -- the default: ``patch`` when this composition has a bundled
#:                patch and nothing else rules it out, otherwise ``packmol``.
#:                Never fails for a composition ``packmol`` would have accepted,
#:                because every reason to decline is enumerated up front by
#:                :meth:`MembraneConfig.bilayer_refusals` and the builder adds a
#:                second, installation-level gate on top of it.
BILAYER_SOURCES = {"packmol", "patch", "auto"}

#: Provenance of the bundled bilayer patches, reproduced in errors and reports
#: because a build that silently substitutes someone else's equilibrated lipids
#: for your own packing has to say whose they are.
BILAYER_PATCH_PROVENANCE = "Zenodo doi:10.5281/zenodo.14776136 (CC-BY-4.0)"

#: Environment variable that points :mod:`prism.membrane.patch` at a directory
#: of extra or overriding patches.  Mirrored here rather than imported because
#: that module pulls in :mod:`numpy`, and this one is imported by the CLI on
#: every run; the two must name the same variable or config would answer "is the
#: patch route available" from a different set of assets than the one that
#: actually gets tiled.
BUNDLED_PATCH_DIR_ENV = "PRISM_MEMBRANE_PATCH_DIR"
_PATCH_MANIFEST_NAME = "manifest.json"
_PATCH_SUBDIR = ("data", "membrane_patches")

#: Lipid force field the bundled patches were equilibrated under, used when the
#: manifest does not say.  A patch is a box of coordinates equilibrated under a
#: *particular* force field, so this is part of the lookup key and not a detail:
#: tiling Lipid21 coordinates and parameterising them as CHARMM36 produces a
#: system grompp accepts and that is not the force field either party asked for.
BUNDLED_PATCH_DEFAULT_LIPID_FF = "lipid21"


def _canonical_composition_key(
    lipid_ff: str,
    names: Sequence[str],
    counts: Sequence[int],
) -> Optional[Tuple[str, Tuple[str, ...], Tuple[int, ...]]]:
    """Return the canonical ``(lipid_ff, lipids, ratio)`` key for a composition.

    One canonicaliser, used by both sides of the patch lookup: the table is
    derived through it from the shipped manifest, and a user's request is
    reduced through it by :func:`bundled_patch_id`.  That is deliberate --  a
    lookup whose two sides normalise separately is a lookup that silently stops
    matching the day one side learns a new spelling.

    Duplicates are merged, synonyms folded through :data:`LIPID_NAME_ALIASES`,
    names sorted, and the ratio divided down to lowest terms, so
    ``POPC:CHL1 = 108:20``, ``CHL1:POPC = 5:27`` and the manifest's own
    ``[["CHL", 20], ["POPC", 108]]`` all reduce to one key.
    """
    if not isinstance(lipid_ff, str) or not lipid_ff.strip():
        return None
    if len(names) != len(counts) or not names:
        return None

    merged: Dict[str, int] = {}
    for name, count in zip(names, counts):
        if not isinstance(name, str) or not name.strip():
            return None
        if isinstance(count, bool) or not isinstance(count, Integral) or int(count) <= 0:
            return None
        lipid = name.strip().upper()
        lipid = LIPID_NAME_ALIASES.get(lipid, lipid)
        merged[lipid] = merged.get(lipid, 0) + int(count)

    ordered = sorted(merged.items())
    divisor = 0
    for _name, count in ordered:
        divisor = math.gcd(divisor, count)
    divisor = divisor or 1
    return (
        lipid_ff.strip().lower(),
        tuple(name for name, _count in ordered),
        tuple(count // divisor for _name, count in ordered),
    )


def _patch_manifest_path() -> Path:
    """Return the manifest path, honouring :data:`BUNDLED_PATCH_DIR_ENV`."""
    override = os.environ.get(BUNDLED_PATCH_DIR_ENV)
    if override:
        return Path(override).expanduser() / _PATCH_MANIFEST_NAME
    return Path(__file__).resolve().parent.parent.joinpath(
        *_PATCH_SUBDIR, _PATCH_MANIFEST_NAME
    )


def _load_bundled_patches() -> Tuple[Dict[Tuple[str, Tuple[str, ...], Tuple[int, ...]], str], str]:
    """Derive the composition -> patch-id table from the shipped manifest.

    Returns ``(table, problem)``; ``problem`` is empty when the manifest was
    read cleanly and otherwise says what went wrong, so a refusal can explain
    that the patch route is unavailable because the *install* is incomplete
    rather than because the composition is unsupported.

    Cheap and offline by construction: this reads one small JSON file and never
    opens a coordinate file.  It also never raises -- a PRISM installed as a
    wheel without package data must still import and still build the PACKMOL
    route, which needs no patches at all.
    """
    path = _patch_manifest_path()
    try:
        with open(path, "r", encoding="utf-8") as handle:
            raw = json.load(handle)
    except FileNotFoundError:
        return {}, f"no bundled patch manifest was found at {path}"
    except (OSError, ValueError, UnicodeDecodeError) as exc:
        return {}, f"the bundled patch manifest at {path} could not be read ({exc})"

    if not isinstance(raw, dict) or not isinstance(raw.get("patches"), dict):
        return {}, f"the bundled patch manifest at {path} lists no 'patches' mapping"

    default_ff = raw.get("lipid_ff", BUNDLED_PATCH_DEFAULT_LIPID_FF)
    table: Dict[Tuple[str, Tuple[str, ...], Tuple[int, ...]], str] = {}
    unreadable: List[str] = []
    collisions: List[str] = []
    for key, entry in sorted(raw["patches"].items()):
        composition = entry.get("composition") if isinstance(entry, dict) else None
        if not isinstance(composition, (list, tuple)) or not composition:
            unreadable.append(str(key))
            continue
        try:
            names = [str(pair[0]) for pair in composition]
            counts = [int(pair[1]) for pair in composition]
        except (TypeError, ValueError, IndexError, KeyError):
            unreadable.append(str(key))
            continue
        composition_key = _canonical_composition_key(
            str(entry.get("lipid_ff", default_ff)), names, counts
        )
        if composition_key is None:
            unreadable.append(str(key))
            continue
        if composition_key in table:
            # Two patches for one composition: the request cannot say which is
            # meant, so picking one silently would hand the user someone else's
            # bilayer under the name of the one they asked for.  Iteration is
            # sorted, so the winner is at least deterministic, and the clash is
            # reported rather than absorbed.
            collisions.append(f"{table[composition_key]}/{key}")
            continue
        table[composition_key] = str(key)

    problems = []
    if unreadable:
        problems.append(
            f"{len(unreadable)} entr{'y' if len(unreadable) == 1 else 'ies'} in the bundled "
            f"patch manifest at {path} could not be interpreted "
            f"({', '.join(sorted(unreadable))})"
        )
    if collisions:
        problems.append(
            f"the bundled patch manifest at {path} maps more than one patch to the same "
            f"composition ({', '.join(sorted(collisions))}); the first of each pair is used"
        )
    return table, "; ".join(problems)


#: Compositions that ship with a pre-equilibrated patch, keyed by
#: ``(lipid force field, lipid names, ratio)`` -- canonicalised by
#: :func:`_canonical_composition_key`, so the ratio is in lowest terms and the
#: lipids are sorted.  The value is the patch identifier
#: :mod:`prism.membrane.patch` resolves to a bundled file, which is the
#: manifest's own key for it.
#:
#: **Derived, never restated.**  ``prism/data/membrane_patches/manifest.json`` is
#: the authority on which patches exist -- it is what ships the files, their
#: checksums and their compositions -- so this table is read from it at import.
#: The previous hand-maintained copy listed one composition while the package
#: shipped seven, and every request for the other six was refused with a hint
#: generated from the stale copy: a duplicate of a shipped data file drifts, and
#: this one did.
#:
#: It deliberately carries no geometry (periodicity, lipid count, area per
#: lipid): those are properties of the bundled file and must be read from it.
#: Answering "is the patch route available" from the *request alone* -- at config
#: time, before orientation and parameterisation have been paid for -- needs only
#: the composition, which is why this much is loaded eagerly and no more.
BUNDLED_BILAYER_PATCHES, _BUNDLED_PATCH_MANIFEST_PROBLEM = _load_bundled_patches()


def _describe_bundled_patches() -> str:
    """Render the bundled compositions as the requests that select them."""
    if not BUNDLED_BILAYER_PATCHES:
        # Say *why* there are none.  "no patch is bundled for POPC (available:
        # none)" sends the user looking for a composition problem when the real
        # answer is that the install shipped without its package data.
        if _BUNDLED_PATCH_MANIFEST_PROBLEM:
            return f"none -- {_BUNDLED_PATCH_MANIFEST_PROBLEM}"
        return "none"
    described = []
    for (ff, names, counts), key in BUNDLED_BILAYER_PATCHES.items():
        composition = f"{ff} {'/'.join(names)}"
        if len(names) > 1:
            composition += f" ratio {':'.join(str(count) for count in counts)}"
        described.append(f"{composition} ({key})")
    return ", ".join(sorted(described))

#: Protein-to-lipid heavy-atom distance (Angstrom) below which a lipid is deleted
#: to make room for the solute.  Whole lipids go, not individual atoms: half a
#: lipid is not a molecule.  3.2 A is the prototype's measured value -- it is
#: roughly a heavy-atom contact distance, so it clears the lipids the solute
#: actually occupies without carving a wider hole than the solute fills.
DEFAULT_PATCH_CLASH_A = 3.2
#: Lipid-to-lipid heavy-atom distance (Angstrom) used to repair the seams where
#: two tiles meet.  Tighter than the protein-lipid criterion on purpose: the
#: patch is already equilibrated, so genuine lipid-lipid contacts at ~3 A are
#: correct structure and must not be deleted; only true overlaps are.
DEFAULT_PATCH_SEAM_CLASH_A = 2.55
#: Below this the deleted shell is thinner than a heavy-atom contact, so lipids
#: are left overlapping the solute and the first minimisation step has to undo it
#: (advisory).
PATCH_CLASH_MIN_WARNING_A = 2.0
#: Above this the hole cut for the solute is far wider than the solute, leaving a
#: vacuum gap at the protein-lipid interface that no equilibration closes on a
#: useful timescale (fatal).
PATCH_CLASH_MAX_A = 10.0

#: The flat shrink factor ``openmm.app.Modeller.addMembrane`` hardcodes, and the
#: one Wolf et al. (2010) published.  Named here so a request to "do what OpenMM
#: does" has a spelling, and so the reason *not* to copy it has somewhere to
#: live -- see :attr:`MembraneConfig.patch_shrink_factor`, which records what a
#: flat 0.5 measured on a real protein.
#:
#: Mirrored from :data:`prism.membrane.patch.OPENMM_SHRINK_FACTOR` rather than
#: imported, for the same reason :data:`BUNDLED_PATCH_DIR_ENV` is: that module
#: pulls in :mod:`numpy` and this one is imported by the CLI on every run.  The
#: two must hold the same number.
OPENMM_SHRINK_FACTOR = 0.5

#: Range :func:`prism.membrane.patch.area_budget_shrink_factor` clamps its
#: solved factor to, mirrored from ``patch.SHRINK_FACTOR_BOUNDS`` (see above for
#: why it is copied rather than imported).
#:
#: Advisory here, not fatal.  ``patch.py`` clamps its own *solved* factor to this
#: interval, but it accepts any explicit float in (0, 1), and a methods study
#: comparing against OpenMM's 0.5 or sweeping the factor is a legitimate thing to
#: want.  So a value outside these bounds is warned about with the reason -- below
#: the floor the survivors end up so deep inside the real solute that the
#: grow-back has to drag them through it, above the ceiling the shrink does
#: nothing plain deletion did not already do -- and still honoured.
SHRINK_FACTOR_WARNING_BOUNDS = (0.45, 0.95)

#: Spellings that mean "no shrink", normalised to ``None`` on the way in.
#:
#: ``none`` is the one that earns this table.  Unquoted ``none`` in YAML is the
#: *string* ``"none"``, not a null -- YAML spells null as ``null``, ``~`` or an
#: empty value -- so a Python-minded user writing ``patch_shrink_factor: none``
#: hands this module a string.  Rejecting it would be technically defensible and
#: practically obtuse: they asked for the default, in the most natural spelling
#: available, and a config error is a poor reward.  ``prism.membrane.builder``
#: accepts the same set, so normalising here means the stored config is canonical
#: (``None``) and the builder never has to re-interpret a string PRISM already
#: understood.
_SHRINK_FACTOR_OFF_WORDS = frozenset({"", "none", "off", "false", "no"})

#: Lipids that carry no phosphorus.  The GROMACS path locates the two leaflets
#: from lipid phosphorus atoms, so a composition built only from these has no
#: leaflet planes to find, no way to define the hydrophobic slab, and therefore
#: no way to purge core water.  Such a request is pushed back to the PACKMOL
#: path instead of being solvated blind.  Taken from PACKMOL-Memgen's
#: ``memgen.parm``, where the sterols are the only phosphorus-free entries.
#:
#: The patch route needs the same set for its own reason: after tiling it has to
#: say which leaflet each surviving lipid belongs to, and it identifies leaflets
#: from the phosphorus positions too.  A pure-sterol composition gives it nothing
#: to key on, so it is refused there as well -- see :meth:`bilayer_refusals`.
#:
#: Cholesterol is listed under both of its spellings on purpose.  ``CHL1`` is the
#: PACKMOL-Memgen molecule name and ``CHL`` the Lipid21 residue name used by the
#: bundled patches and by :mod:`prism.membrane.patch`; :data:`LIPID_NAME_ALIASES`
#: folds the latter into the former on the way in, but this guard decides whether
#: a bilayer gets solvated blind, so it holds for whichever vocabulary reaches it
#: rather than relying on a normalisation step upstream having happened.
PHOSPHORUS_FREE_LIPIDS = {"CHL1", "CHL", "ERG", "CAM", "SIT", "STI"}

#: How the solute is put into the bilayer frame.
#:
#: ``auto``        -- the default: measure the input against the membrane frame
#:                    and run MEMEMBED only if it is not already framed.  A
#:                    structure that is already framed is passed through with its
#:                    coordinates untouched, so a hand-prepared input still gets
#:                    exactly what was prepared; anything else -- notably a plain
#:                    RCSB download, which is in the crystal frame -- is oriented
#:                    instead of being refused.  The measurement uses
#:                    ``center_tolerance``.
#: ``preoriented`` -- assert that the input is already membrane-framed.  Still
#:                    the right answer when you want the assertion itself: it
#:                    never moves an atom, and with ``validate_orientation`` it
#:                    fails rather than silently orienting a structure you
#:                    believe is already correct.
#: ``memembed``    -- always run MEMEMBED, even on an input that is already framed.
#: ``opm``         -- fetch a pre-oriented structure from OPM by ``pdb_id``.
#: ``ppm``         -- recognised but not automated; see
#:                    :data:`UNAUTOMATED_ORIENTATION_METHODS`.
#:
#: **This set is duplicated as string literals elsewhere and the copies have to
#: agree.**  ``prism/builder/cli.py`` imports it, so its ``--membrane-orient``
#: choices follow automatically; but ``prism/membrane/packmol_memgen.py`` and
#: ``prism/mcp/build.py`` each hardcode ``{"preoriented", "opm", "ppm",
#: "memembed"}`` and reject anything else.  Adding a method here therefore does
#: not add it there: until those two literals learn ``"auto"``, a build that
#: reaches PACKMOL-Memgen with ``orient="auto"`` dies on
#: ``Unsupported membrane orientation method``.  For both of them the correct
#: handling is to treat ``"auto"`` exactly as ``"preoriented"``, because PRISM
#: resolves the route and writes the oriented structure itself before either is
#: reached -- whatever they receive has already been put in the membrane frame.
ORIENTATION_METHODS = {"opm", "ppm", "preoriented", "memembed", "auto"}
#: Orientation methods PRISM recognises by name but cannot run itself.  They
#: stay in ``ORIENTATION_METHODS`` so the request is rejected with the actual
#: workaround instead of being reported as a typo, and so the rejection happens
#: during validation rather than hours later inside the build.
UNAUTOMATED_ORIENTATION_METHODS = {
    "ppm": (
        "orient=ppm is not automated by PRISM; run the PPM web server or the standalone "
        "'immers' binary, then pass its output with orient=preoriented."
    ),
}
DEFAULT_PRODUCTION_NS = 500.0
#: Bilayers are equilibrated well above the gel-transition temperature of the
#: lipid, so the membrane default is warmer than PRISM's soluble-system default.
#: Declared once here and read by both parsers via :func:`membrane_defaults`.
DEFAULT_MEMBRANE_TEMPERATURE_K = 303.15
MEMBRANE_TIMESTEP_PS = 0.002
#: Timestep for the first (restrained) equilibration stages.  A freshly packed
#: bilayer has poor local geometry; starting straight at 2 fs is a common cause
#: of LINCS failures, so the early stages integrate at 1 fs.
MEMBRANE_EQUIL_TIMESTEP_PS = 0.001
#: A ``preoriented`` input must be centred on the bilayer midplane (z = 0),
#: because PACKMOL-Memgen always builds the bilayer there and leaves the solute
#: at its input coordinates.  Structures further off than this are rejected
#: rather than silently packed outside the membrane.  (The patch route places
#: its bilayer at z = 0 too, so the requirement is the same on both.)
#:
#: Under the ``auto`` default this same tolerance decides a *route* and not just
#: a verdict: a solute that measures within it is passed through untouched, and
#: one that does not is handed to MEMEMBED.  Widening it therefore does not only
#: relax a check -- it makes ``auto`` accept a frame it would otherwise have
#: corrected.
DEFAULT_CENTER_TOLERANCE_A = 15.0


def _require_bool(value, field_name: str) -> bool:
    """Return a real boolean without accepting truthy strings or integers."""
    if type(value) is not bool:
        raise ValueError(f"Membrane {field_name} must be a boolean (true or false).")
    return value


def _require_finite_number(
    value,
    field_name: str,
    *,
    allow_zero: bool = False,
) -> float:
    """Parse a YAML/CLI number without silently accepting strings, NaN, or Inf."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"Membrane {field_name} must be a finite number.")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Membrane {field_name} must be finite (NaN/Inf are not allowed).")
    if number < 0 or (number == 0 and not allow_zero):
        qualifier = "non-negative" if allow_zero else "positive"
        raise ValueError(f"Membrane {field_name} must be {qualifier}.")
    return number


def _finite_or_none(value) -> Optional[float]:
    """Return ``value`` as a float, or ``None`` when it is not a finite number.

    Range checks and advisory warnings run alongside the type checks rather than
    after them, so they need a way to say "this value is not mine to complain
    about" instead of raising a second, redundant error.
    """
    if isinstance(value, bool) or not isinstance(value, Real):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _require_positive_int(value, field_name: str) -> int:
    """Parse a count without accepting bools, floats, or numeric strings.

    Mirrors :func:`_require_ratio`'s strictness: a threshold silently truncated
    from ``3e4`` or read from the string ``"30000"`` would compare against water
    counts and pick a solvation backend, so a lossy cast is not acceptable here.
    """
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"Membrane {field_name} must be a positive integer.")
    number = int(value)
    if number < 1:
        raise ValueError(f"Membrane {field_name} must be a positive integer.")
    return number


def membrane_production_steps(value) -> int:
    """Return production steps, rejecting durations shorter than one step."""
    duration_ns = _require_finite_number(value, "production duration")
    steps = int(duration_ns * 1000 / MEMBRANE_TIMESTEP_PS)
    if steps < 1:
        raise ValueError(
            "Membrane production duration must cover at least one "
            f"{MEMBRANE_TIMESTEP_PS:g} ps integration step."
        )
    return steps


def _require_ratio(value) -> List[int]:
    """Parse an integer lipid-ratio sequence without lossy ``int(...)`` casts."""
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise ValueError("Membrane lipid ratio must be a sequence of positive integers.")
    ratio = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, Integral):
            raise ValueError("Membrane lipid ratios must be positive integers; floats and strings are not allowed.")
        item = int(item)
        if item <= 0:
            raise ValueError("Membrane lipid ratios must be positive integers.")
        ratio.append(item)
    return ratio


def _require_lipids(value) -> List[str]:
    if isinstance(value, str):
        value = [value]
    if isinstance(value, (bytes, dict)) or not isinstance(value, (list, tuple)):
        raise ValueError("Membrane lipids must be a lipid name or a sequence of lipid names.")
    lipids = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise ValueError("Membrane lipid names must be non-empty strings.")
        lipid = item.strip().upper()
        if "/" in lipid:
            # PACKMOL-Memgen's ``A//B`` leaflet syntax needs a matching
            # per-leaflet ``--ratio``, which MembraneConfig cannot express: the
            # request would be packed with a single-leaflet ratio and the
            # per-leaflet names would then never match any Amber residue.
            raise ValueError(
                f"Asymmetric bilayer syntax in membrane lipid name {item!r} is not supported; "
                "specify one symmetric lipid composition."
            )
        if ":" in lipid or any(character.isspace() for character in lipid):
            raise ValueError(f"Invalid membrane lipid name {item!r}.")
        lipids.append(LIPID_NAME_ALIASES.get(lipid, lipid))
    if not lipids:
        raise ValueError("At least one membrane lipid must be specified.")
    return lipids


def _require_choice(value, field_name: str, choices) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Membrane {field_name} must be one of: {', '.join(sorted(choices))}.")
    normalized = value.strip().lower()
    if normalized not in choices:
        raise ValueError(
            f"Unknown membrane {field_name} {value!r}; choose one of: {', '.join(sorted(choices))}."
        )
    return normalized


def _require_nonempty_string(value, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Membrane {field_name} must be a non-empty string.")
    return value.strip()


#: What the number means, quoted into every rejection.  A user who typed a bad
#: value needs to know what a good one would have done, and "must be in (0, 1)"
#: does not tell them that.
_SHRINK_FACTOR_MEANING = (
    "patch_shrink_factor is the lateral scale the solute is temporarily shrunk to "
    "before the overlapping lipids are deleted (0.7 = 70% of its real width in the "
    "membrane plane); the solute is then grown back in steps, pushing the survivors "
    "out into contact"
)


def _require_shrink_factor(value):
    """Normalise ``patch_shrink_factor`` to ``None``, ``"auto"`` or a float.

    The three forms :func:`prism.membrane.patch.delete_clashing_lipids` accepts,
    and nothing else.  Vocabulary is normalised (case, surrounding whitespace,
    the off-spellings in :data:`_SHRINK_FACTOR_OFF_WORDS`) exactly as
    :func:`_require_choice` and :func:`_require_lipids` normalise theirs, so the
    stored config carries one spelling per meaning.

    A *quoted number* is refused rather than coerced.  That is the same line this
    module already draws for ``solvate_auto_threshold`` (``"30000"`` is rejected,
    not parsed): normalising a vocabulary is safe, casting a string to a number
    hides a type confusion.  It matters more here than usual because this field's
    string domain is already occupied -- ``"auto"`` is a real value -- so
    accepting ``"0.7"`` too would mean a mistyped ``"0..7"`` or ``"atuo"`` has to
    be told apart from a valid string by a float parse that happens somewhere
    else.  Refusing here puts the message next to the typo.
    """
    if value is None:
        return None

    low, high = SHRINK_FACTOR_WARNING_BOUNDS
    advice = (
        f"Use 'auto' (recommended -- PRISM solves the factor from the area budget), "
        f"a number strictly between 0 and 1 (the useful range is about "
        f"{low:g}-{high:g}; OpenMM's own flat value is {OPENMM_SHRINK_FACTOR:g}), "
        f"or null to leave the shrink off and delete lipids against the full-size "
        f"solute."
    )

    if isinstance(value, str):
        text = value.strip().lower()
        if text in _SHRINK_FACTOR_OFF_WORDS:
            return None
        if text == "auto":
            return "auto"
        raise ValueError(
            f"Membrane patch_shrink_factor {value!r} is not a recognised setting "
            f"(a quoted number such as '0.7' is not accepted either -- write it "
            f"unquoted, as a YAML float). {_SHRINK_FACTOR_MEANING}. {advice}"
        )

    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(
            f"Membrane patch_shrink_factor must be a number, 'auto', or null; got "
            f"{type(value).__name__}. {_SHRINK_FACTOR_MEANING}. {advice}"
        )

    number = float(value)
    if not math.isfinite(number):
        raise ValueError(
            "Membrane patch_shrink_factor must be finite (NaN/Inf are not allowed). "
            f"{_SHRINK_FACTOR_MEANING}. {advice}"
        )

    # The interval is open at BOTH ends, and each end is closed for its own
    # reason.  At 1.0 ``patch.py`` skips the scaling entirely and builds no
    # growth plan, so the request is honoured by doing nothing -- a no-op that
    # reports success is the worst of the three outcomes here, because the user
    # believes the algorithm ran.  At or below 0 there is no solute left to
    # delete lipids against.
    if number >= 1.0:
        what = (
            "expands it instead"
            if number > 1.0
            else "leaves it at full size, so the shrink-and-grow algorithm would run "
            "its minimisations and report success having produced exactly what plain "
            "deletion produces"
        )
        raise ValueError(
            f"Membrane patch_shrink_factor must be less than 1; {number:g} does not "
            f"shrink the solute -- it {what}. Either way the lipids are deleted against "
            f"a solute at least its full size, so there is no vacuum annulus left for "
            f"the grow-back to close. {_SHRINK_FACTOR_MEANING}. {advice}"
        )
    if number <= 0.0:
        raise ValueError(
            f"Membrane patch_shrink_factor must be greater than 0; {number:g} does not "
            f"shrink the solute but collapses it to a point (or, below 0, mirrors it), "
            f"leaving no footprint to delete lipids against. {_SHRINK_FACTOR_MEANING}. "
            f"{advice}"
        )
    return number


def bundled_patch_id(lipid_ff, lipids, ratio=None) -> Optional[str]:
    """Return the bundled patch for this composition, or ``None`` if there is none.

    The lookup is canonical, so the same physical bilayer is recognised however
    it was spelled: lipid names are normalised through
    :data:`LIPID_NAME_ALIASES`, sorted, and the ratio is carried along with them
    and reduced to lowest terms.  ``POPC:CHL1 = 3:1`` and ``CHL1:POPC = 2:6``
    therefore hit the same entry, and a plain single-lipid request matches
    whatever ratio it was given (with one lipid, the ratio says nothing).

    Returns ``None`` -- rather than raising -- for a malformed composition, so
    the "no patch for this" diagnostic never masks the "this is not a valid lipid
    list" one.
    """
    try:
        names = _require_lipids(lipids)
        counts = _require_ratio(ratio)
    except ValueError:
        return None
    if not counts:
        counts = [1] * len(names)
    if len(counts) != len(names):
        return None

    # Canonicalised by the same function that derived the table from the
    # manifest, so a spelling the table understands is a spelling a request can
    # use.  (Duplicates merge on the way through: ["POPC", "POPC"] with ratio
    # [1, 1] is a pure POPC bilayer and matches the pure POPC patch.)
    key = _canonical_composition_key(lipid_ff, names, counts)
    if key is None:
        return None
    return BUNDLED_BILAYER_PATCHES.get(key)


def bilayer_patch_citation(patch_id: Optional[str] = None) -> str:
    """Return the citation for a build whose bilayer *was* tiled from a patch.

    Past tense, and deliberately not emitted by
    :meth:`MembraneConfig.validation_warnings`: a citation is a claim about what
    ran, and a configuration only says what was requested.  The patch route can
    still decline at build time -- a missing or corrupt asset, a capability gate,
    an unorientable solute -- and a build that fell back to PACKMOL packing must
    not carry an instruction to cite equilibrated coordinates it never opened.

    Call this from whatever observed the tiling actually happen, and pass the
    patch identifier it used so the record names the specific asset rather than
    the collection.  :mod:`prism.membrane.patch` independently writes the full
    licence attribution into the structures it emits (PDB ``REMARK`` records and
    ``MembraneComposition.citation``), so the CC-BY obligation travels with the
    artefacts whether or not this sentence is surfaced in a report.
    """
    asset = ""
    if isinstance(patch_id, str) and patch_id.strip():
        asset = f" '{patch_id.strip()}'"
    return (
        f"This system's bilayer was tiled from the bundled pre-equilibrated patch{asset} "
        f"({BILAYER_PATCH_PROVENANCE}) rather than packed for this system; cite it "
        "alongside the lipid force field."
    )


@dataclass
class MembraneConfig:
    """Configuration for building a protein-membrane system."""

    enabled: bool = False
    #: Lipid composition (one or more lipid names).
    lipids: List[str] = field(default_factory=lambda: ["POPC"])
    #: Molar ratio per lipid (parallel to ``lipids``); defaults to equal parts.
    ratio: List[int] = field(default_factory=list)
    #: Orientation source -- see :data:`ORIENTATION_METHODS`.  Defaults to
    #: ``auto``: it keeps an already-framed input byte-for-byte, and orients one
    #: that is not, so ``prism -pf protein.pdb --membrane`` works on a structure
    #: straight from RCSB.  ``preoriented`` -- the previous default -- refused
    #: every such structure, which made the no-YAML command unusable for the
    #: overwhelmingly common case of a PDB in the crystal frame.
    orient: str = "auto"
    #: PDB id to fetch a pre-oriented structure from the OPM database.
    pdb_id: Optional[str] = None
    #: Lipid force field (selects the protein-FF family pairing).
    lipid_ff: str = "lipid21"
    #: Protein and water force fields used by PACKMOL-Memgen parameterization.
    #: These are synchronized with PRISMBuilder's resolved global settings, so
    #: they carry PRISM's vocabulary ("amber14sb"), not PACKMOL-Memgen's
    #: ("ff14SB"); ``packmol_forcefield_pair`` does the translation.
    protein_forcefield: str = "amber14sb"
    water_model: str = "tip3p"
    #: Salt concentration (mol/L) for the aqueous phase.
    salt_concentration: float = 0.15
    positive_ion: str = "Na+"
    negative_ion: str = "Cl-"
    #: Water-layer thickness above/below the bilayer (Angstrom; --dist_wat).
    water_thickness: float = 17.5
    #: Minimum protein-to-box distance in the membrane plane (Angstrom; --dist).
    #: PACKMOL-Memgen's flag, and PACKMOL-Memgen's only: the patch route sizes
    #: its tiling from the bundled patch's own periodicity and honours
    #: ``patch_xy_margin`` instead, because a tile is only clash-free at integer
    #: multiples of that period.  Since ``bilayer_source`` now defaults to
    #: ``auto`` this knob is inert on the route a default build usually takes,
    #: so a non-default value that will not be applied is warned about rather
    #: than quietly dropped -- see :meth:`validation_warnings`.
    xy_distance: float = 17.0
    #: Run a short AMBER (sander/pmemd) minimization on the packed system
    #: before conversion (works around GROMACS-crash regressions).
    minimize: bool = True
    temperature: float = DEFAULT_MEMBRANE_TEMPERATURE_K
    production_ns: float = DEFAULT_PRODUCTION_NS
    #: Verify that a ``preoriented`` input is actually centred on the bilayer
    #: midplane before packing.  Disable only when the caller has already
    #: validated the frame of reference.
    validate_orientation: bool = True
    #: Allowed deviation (Angstrom) of the membrane-spanning region from z = 0.
    center_tolerance: float = DEFAULT_CENTER_TOLERANCE_A
    #: Position-restraint force constants (kJ/mol/nm^2) used by the staged
    #: equilibration protocol.  Restraints are released over several stages.
    posres_fc_backbone: float = 1000.0
    posres_fc_sidechain: float = 500.0
    posres_fc_lipid: float = 1000.0
    #: AMBER parameters for a ligand embedded together with the protein.
    #: ``ligand_frcmod``/``ligand_lib`` feed PACKMOL-Memgen's ``--ligand_param``.
    ligand_frcmod: Optional[str] = None
    ligand_lib: Optional[str] = None
    ligand_resname: str = "LIG"
    #: Use GAFF2 (rather than GAFF) when ligand parameters are supplied. This
    #: selects the leaprc PACKMOL-Memgen sources for the ligand, so it must name
    #: the generation that actually typed ``ligand_frcmod``/``ligand_lib``; a
    #: mismatch is silent (GAFF's atom types are a subset of GAFF2's, with the
    #: same masses) and only shifts force constants. It is therefore derived
    #: from the ligand force field -- see :func:`gaff2_selected` -- and the
    #: standalone default follows what both PACKMOL-Memgen and PRISM's ``-lff``
    #: default to on their own: plain GAFF.
    gaff2: bool = False
    #: Wall-clock budget (seconds) for the PACKMOL-Memgen pack+parametrize+
    #: minimize step.  Packing cost grows with box volume and with how awkward
    #: the solute's surface is, so a large complex -- or a slow filesystem --
    #: can legitimately exceed the default.  A timeout discards hours of work,
    #: so this is configurable rather than hardcoded.
    packmol_timeout: float = 7200.0
    #: GENCAN iteration budgets for PACKMOL.  ``packing_nloop`` covers the
    #: per-species stage (PACKMOL-Memgen default 20) and ``packing_nloop_all``
    #: the final all-together stage (default 100), which dominates runtime
    #: because every iteration evaluates overlaps across the whole system.
    #:
    #: PRISM overrides the all-together default because the extra iterations
    #: measurably buy nothing.  On a 56 064-atom POPC bilayer around T4 lysozyme
    #: plus a GAFF2 ligand, ``nloop_all`` of 20 and of 100 produced a
    #: BYTE-IDENTICAL packed structure (same md5), in 35 min versus 2 h 54 min:
    #: PACKMOL plateaus around loop 21, reports ENDED WITHOUT PERFECT PACKING in
    #: both cases, and keeps the best structure found, so loops 22-101 only spin.
    #: Dropping to 11 is NOT safe -- worst inter-residue contact degrades from
    #: 0.42 A to 0.09 A (essentially superimposed atoms, which is what makes a
    #: first-step LINCS failure likely) and severe clashes triple.
    #:
    #: This is single-system evidence, so the knob stays: raise it for an
    #: unusually awkward solute, and note that PACKMOL-Memgen rejects any value
    #: not greater than its write-out frequency.  ``None`` restores
    #: PACKMOL-Memgen's own default for the per-species stage.
    packing_nloop: Optional[int] = None
    packing_nloop_all: Optional[int] = 20
    #: Which backend builds the **bilayer** -- see :data:`BILAYER_SOURCES`.  Not
    #: to be confused with ``solvate`` directly below, which chooses who places
    #: the **aqueous phase**; this one is about the lipids.
    #:
    #: Defaults to ``auto``, which takes the patch route when this composition
    #: has a bundled patch and falls back to PACKMOL packing when it does not.
    #: The promotion is measured, not assumed: on mGluR7 (PDB 9OMO, 587 528
    #: atoms) the patch/GROMACS pair finished in 75 s with five clean grompp
    #: stages, while the PACKMOL pair was still running at 600 s and, on a
    #: class-C GPCR, its wet pack never converged at all.  For a composition
    #: with no patch ``auto`` is byte-for-byte the old default.
    bilayer_source: str = "auto"
    #: Protein-lipid heavy-atom distance (Angstrom) below which the patch route
    #: deletes a lipid to make room for the solute.
    patch_clash_distance: float = DEFAULT_PATCH_CLASH_A
    #: Lipid-lipid heavy-atom distance (Angstrom) at which the patch route treats
    #: two tiled lipids as overlapping and repairs the seam between tiles.
    patch_seam_clash_distance: float = DEFAULT_PATCH_SEAM_CLASH_A
    #: Lipid margin (Angstrom) kept around the solute in the membrane plane when
    #: sizing the tiled patch.  ``None`` leaves the choice to
    #: :mod:`prism.membrane.patch`, which sizes it from the bundled patch's own
    #: periodicity -- the tiling is only clash-free at integer multiples of that
    #: period, so a margin this module invented could not be honoured exactly
    #: anyway.  Distinct from ``xy_distance``, which is PACKMOL-Memgen's
    #: solute-to-box distance.
    patch_xy_margin: Optional[float] = None
    #: Shrink-and-grow lipid packing on the patch route -- the algorithm of Wolf
    #: et al., *J. Comp. Chem.* **31**, 2169-2174 (2010), the one
    #: ``openmm.app.Modeller.addMembrane`` uses.  ``None`` (off), ``"auto"``, or
    #: a float strictly between 0 and 1.
    #:
    #: **What it fixes.**  Tiling a pre-equilibrated patch and deleting the
    #: lipids the solute overlaps gets the *count* right but leaves every
    #: survivor on the lattice position it held in a protein-free bilayer, so the
    #: bilayer's inner edge stands off the protein: a vacuum annulus exactly where
    #: the first shell of annular lipid belongs.  The algorithm shrinks the solute
    #: laterally, deletes against the *shrunken* copy (so fewer lipids go, and the
    #: ones that stay sit inside the real footprint), then grows the solute back
    #: in increments with a minimisation between each, pushing the survivors out
    #: into real van der Waals contact.  Measured on mGluR7 (PDB 9OMO): rim gap
    #: 7.56 -> 4.29 A, lipids within 4.5 A of the protein 13 -> 49, area per lipid
    #: 64.02 -> 60.62 A^2 against the patch's own 60.43 -- i.e. *toward* the
    #: equilibrated density, not away from it -- with zero sub-2.0 A overlaps and
    #: a converging EM.
    #:
    #: **Why it is on by default.**  It was opt-in while three things were true;
    #: all three were settled by measurement on 6KQI, with the xy box fix in
    #: place so the numbers are not confounded by the global vacuum stripe that
    #: bug produced:
    #:
    #: 1. *The atom-order blocker is closed.*  The solute is rewritten into
    #:    ``[ molecules ]`` order before the bilayer is built around it, so the
    #:    coordinates the grow-back minimises and the topology it minimises
    #:    against agree by construction rather than by luck.
    #: 2. *The rim annulus is real and plain deletion does not fix it.*  With the
    #:    box corrected, plain deletion still leaves the rim at 6.30 A -- van der
    #:    Waals contact is ~3 A -- with 2 lipids within 4.5 A of the solute.  The
    #:    grow-back gives 3.26 A and 42 lipids, and area per lipid moves 65.8 ->
    #:    61.7 A^2, *toward* the patch's equilibrated 60.43 rather than away.
    #:    Nothing is traded for the improvement.
    #: 3. *The cost is not what it looked like.*  A protein-only 6KQI build with
    #:    the grow-back takes 244 s end to end -- faster than OpenMM's 716 s for
    #:    the same job, because minimising to convergence per increment allows a
    #:    dozen large steps where OpenMM's 10 K dynamics needs ~60 small ones.
    #:    The earlier 821 s figure was the grow-back measured alone on the larger
    #:    mGluR7 dimer, not a build's total.
    #:
    #: Set it to ``None`` for plain clash deletion -- appreciably faster, and
    #: sound as a starting structure, but its rim closes during equilibration
    #: instead of at build time.
    #:
    #: **``"auto"`` is not a vague convenience, and 0.5 is not its safe
    #: equivalent.**  OpenMM shrinks by a flat 0.5 whatever the protein, which is
    #: fine only while the solute is small next to the patch: the lipids deleted
    #: scale with the *scaled* footprint, about a quarter of the true one at 0.5.
    #: On the mGluR7 dimer (transmembrane cross-section 2405 A^2 in a 23205 A^2
    #: box) a flat 0.5 deletes 67 lipids where the area budget calls for 88, and
    #: the free area per lipid lands at 57.9 A^2 against the patch's 60.4 -- a 4%
    #: compression of a bilayer that arrived already equilibrated.  ``"auto"``
    #: solves the factor from the area actually left for the survivors
    #: ((box area - cross-section) / (survivors / 2) == the patch's own APL) and
    #: returns 0.700 on this protein; it reduces to roughly 0.5 for a small one.
    #: So prefer ``"auto"``; reach for :data:`OPENMM_SHRINK_FACTOR` only to
    #: reproduce OpenMM or Wolf et al. deliberately, not as a default-by-analogy.
    patch_shrink_factor: Optional[float | str] = "auto"
    # ``allow_mixed_14_scaling`` was removed here.  It converted a system whose
    # AMBER 1-4 scaling is not uniform by applying the most common factor to
    # every 1-4 pair anyway -- knowingly wrong physics, offered because the
    # alternative was a hard stop.  The converter now writes explicit per-pair
    # ``[ pairs ]`` V/W with ``gen_pairs = no`` (see
    # :mod:`prism.membrane.solvent_params`), so Lipid21's ``SCNB=6.0``
    # ``cD-cD-cD-cD`` acyl-tail torsion is carried exactly alongside the
    # ``SCNB=2.0`` of everything else and no global ``fudgeLJ`` has to be wrong.
    # The approximation this enabled therefore no longer exists to be requested,
    # and its warning ("...will be scaled by a single global fudgeLJ...the
    # acyl-tail 1-4 pairs by 3x") had become a false claim about the emitted
    # force field.  :func:`parse_membrane_yaml` rejects an explicit ``true`` so
    # a stale config is told, rather than silently getting something better than
    # it asked for.
    #: Which backend places water and ions -- see :data:`SOLVATE_MODES`.
    #:
    #: Defaults to ``auto``.  This knob used to say the GROMACS path was "opt-in
    #: until it has regression-passed on more than one system"; that condition
    #: has since been met, so the sentence is gone rather than left standing
    #: next to a default it no longer describes.  The systems it was waiting
    #: for: a 1371-atom test protein (163 836 atoms, 30.5 s) and mGluR7 (587 528
    #: atoms, 75 s), the latter verified end to end -- ``gen-pairs = no`` with
    #: explicit per-pair 1-4 LJ, SOLU/MEMB/SOLV index groups, four restraint
    #: files, and all five grompp stages passing.
    #:
    #: Two systems is a promotion argument, not a proof, which is why ``auto``
    #: rather than ``gromacs`` is the default: it leaves every composition the
    #: GROMACS path cannot serve on PACKMOL automatically, and an explicit
    #: ``solvate: packmol`` still gets PACKMOL, unchanged, on every input.
    solvate: str = "auto"
    #: Water count above which ``solvate="auto"`` selects the GROMACS path.
    solvate_auto_threshold: int = DEFAULT_SOLVATE_AUTO_THRESHOLD
    #: Delete water that ``gmx solvate`` placed inside the hydrophobic slab.
    #: Mandatory for a correct bilayer (an inflated carbon radius alone never
    #: reaches zero core waters); the knob exists only so that a methods study
    #: can quantify what the purge removes.
    core_water_purge: bool = True
    #: Carbon radius (nm) written into the ``vdwradii.dat`` used by
    #: ``gmx solvate`` on the GROMACS path.
    vdwradii_carbon: float = DEFAULT_VDWRADII_CARBON_NM
    #: Run a zero-step ``grompp``/``mdrun`` energy gate on the finished system.
    #: A scrambled topology/coordinate pairing still produces a valid ``.tpr``
    #: behind a single warning, so the cheapest reliable detector of a broken
    #: build is its potential energy -- worth ~10-30 s at the end of a build
    #: that took tens of minutes.
    preflight_energy: bool = True
    #: Extra disulfide bonds to declare to tleap, as pairs of residue numbers
    #: *in tleap's own numbering of the solute*.
    #:
    #: Both bilayer routes auto-detect disulfides from SG-SG distances (see
    #: :func:`prism.ptm.disulfides.detect_disulfides`), so this is an override
    #: for the cases geometry cannot settle: a bond the deposited coordinates
    #: place just outside the distance cutoff, or one that has to be declared
    #: for a model built without its partner in place.  Declared pairs are
    #: merged with the detected ones and de-duplicated.
    #:
    #: This field exists because both routes already read it -- previously via a
    #: ``getattr(config, "disulfide_pairs", None)`` that could never find
    #: anything, so the PACKMOL route silently declared no disulfides at all.
    disulfide_pairs: Optional[List[Tuple[int, int]]] = None

    def has_ligand(self) -> bool:
        """True when AMBER ligand parameters are configured for this build."""
        return bool(self.ligand_frcmod and self.ligand_lib)

    def family(self) -> str:
        """Protein force-field family implied by the lipid force field."""
        lipid_ff = self.lipid_ff.lower() if isinstance(self.lipid_ff, str) else ""
        try:
            return LIPID_FF_FAMILY[lipid_ff]
        except KeyError as exc:
            raise ValueError(f"Unknown lipid force field {self.lipid_ff!r}.") from exc

    def resolved_lipids(self) -> List[str]:
        """Validated lipid names in the spelling PACKMOL-Memgen understands.

        Directly constructed configs never pass through a parser, so this is
        where a synonym such as ``CHOL`` becomes ``CHL1`` for them.
        """
        return _require_lipids(self.lipids)

    def resolved_ratio(self) -> List[int]:
        if not isinstance(self.lipids, (list, tuple)) or not self.lipids:
            raise ValueError("At least one membrane lipid must be specified.")
        ratio = _require_ratio(self.ratio)
        if not ratio:
            return [1] * len(self.lipids)
        if len(ratio) != len(self.lipids):
            raise ValueError(
                f"Lipid ratio count ({len(ratio)}) must match lipid count ({len(self.lipids)})."
            )
        return ratio

    def packmol_forcefield_pair(self):
        """Return validated PACKMOL-Memgen protein/water force-field names."""
        protein_key = str(self.protein_forcefield or "").strip().lower()
        water_key = str(self.water_model or "").strip().lower()
        try:
            protein_ff = PACKMOL_PROTEIN_FF_ALIASES[protein_key]
        except KeyError as exc:
            supported = ", ".join(sorted({"amber14sb", "amber19sb"}))
            raise ValueError(
                f"Unsupported membrane protein force field {self.protein_forcefield!r}; "
                f"supported PRISM names: {supported}."
            ) from exc
        try:
            water_ff = PACKMOL_WATER_MODEL_ALIASES[water_key]
        except KeyError as exc:
            supported = ", ".join(sorted(PACKMOL_WATER_MODEL_ALIASES))
            raise ValueError(
                f"Unsupported membrane water model {self.water_model!r}; supported models: {supported}."
            ) from exc
        if (protein_ff, water_ff) not in PACKMOL_FORCEFIELD_PAIRS:
            raise ValueError(
                "Unsupported membrane protein/water force-field pairing "
                f"{self.protein_forcefield!r}/{self.water_model!r}; use "
                "amber14sb/tip3p or amber19sb/opc."
            )
        return protein_ff, water_ff

    def bundled_patch(self) -> Optional[str]:
        """Return the bundled bilayer patch for this composition, or ``None``."""
        return bundled_patch_id(self.lipid_ff, self.lipids, self.ratio)

    def bilayer_refusals(self) -> List[str]:
        """Return reasons this request cannot build its bilayer from a patch.

        Deliberately shaped like :meth:`solvate_refusals`: a list of causes in
        prose, cheap to compute from the request alone.  How they are treated
        differs, though, and the difference is the point -- for
        ``bilayer_source="auto"`` a refusal is a silent, reported downgrade to
        PACKMOL, while for an explicit ``bilayer_source="patch"`` it is a fatal
        configuration error.  Asking for the patch route and being given the
        30-minute one without comment is not an acceptable answer to an explicit
        request.
        """
        refusals = []

        try:
            lipids = _require_lipids(self.lipids)
        except ValueError:
            # Already a validation error in its own right; not repeated here.
            lipids = []

        lipid_ff = self.lipid_ff.strip().lower() if isinstance(self.lipid_ff, str) else ""

        if lipids and self.bundled_patch() is None:
            composition = "/".join(lipids)
            try:
                requested_ratio = _require_ratio(self.ratio)
            except ValueError:
                requested_ratio = []
            if len(lipids) > 1 and len(requested_ratio) == len(lipids):
                composition += f" ratio {':'.join(str(count) for count in requested_ratio)}"
            refusals.append(
                f"no pre-equilibrated bilayer patch is bundled for {composition} under "
                f"lipid_ff={lipid_ff or self.lipid_ff!r} (available: "
                f"{_describe_bundled_patches()})"
            )

        # Sterols carry no phosphorus, and the patch route identifies which
        # leaflet each surviving lipid belongs to from phosphorus positions.
        if lipids and all(lipid in PHOSPHORUS_FREE_LIPIDS for lipid in lipids):
            refusals.append(
                "the lipid composition contains no phosphorus, so the leaflets of the tiled "
                "patch cannot be identified"
            )

        # The other axis, and the reason this refusal lives here rather than
        # only in ``validation_errors``.  PACKMOL-Memgen places lipids and water
        # in a single run; the patch route never invokes it, so it cannot honour
        # an explicit ``solvate='packmol'``.  Treating that as a reason the patch
        # route is unavailable is what lets ``bilayer_source='auto'`` keep its
        # promise -- it downgrades to packing and the user gets the PACKMOL water
        # they asked for, instead of being moved onto the GROMACS solvator in
        # silence.  For an explicit ``bilayer_source='patch'`` the same refusal
        # is fatal, which is correct: that request contradicts itself.
        if isinstance(self.solvate, str) and self.solvate.strip().lower() == "packmol":
            refusals.append(
                "solvate='packmol' asks PACKMOL-Memgen to place the water and ions, and the "
                "patch route never runs it (use solvate='gromacs', or 'auto')"
            )

        return refusals

    def resolved_bilayer_source(self) -> str:
        """Return the bilayer backend to use: ``"packmol"`` or ``"patch"``.

        ``auto`` resolves to ``patch`` only when nothing in
        :meth:`bilayer_refusals` rules it out, so it can never turn a request
        that PACKMOL would have built into one that fails.
        """
        source = _require_choice(self.bilayer_source, "bilayer source", BILAYER_SOURCES)
        if source == "auto":
            return "packmol" if self.bilayer_refusals() else "patch"
        return source

    def solvate_refusals(self) -> List[str]:
        """Return reasons this composition cannot use the GROMACS solvation path.

        These are *not* validation errors: an unsupported composition still
        builds, it just builds the slow way.  They are collected here rather
        than discovered inside the build because both are cheap to check from
        the request alone, and because a refusal has to be reported to the user
        -- silently solvating a membrane whose leaflets cannot be located would
        put water in the hydrophobic core, which is precisely the failure the
        fast path is designed to avoid.
        """
        refusals = []

        try:
            lipids = _require_lipids(self.lipids)
        except ValueError:
            # Malformed composition is already a validation error; do not
            # duplicate it as a refusal.
            lipids = []
        if lipids and all(lipid in PHOSPHORUS_FREE_LIPIDS for lipid in lipids):
            refusals.append(
                "the lipid composition contains no phosphorus, so the leaflet planes the "
                "core-water purge needs cannot be located"
            )

        lipid_ff = self.lipid_ff.strip().lower() if isinstance(self.lipid_ff, str) else ""
        if lipid_ff and lipid_ff not in PACKMOL_PARAMETRIZED_LIPID_FFS:
            # The dry pack still ends in ``--parametrize``; a lipid force field
            # that PACKMOL-Memgen cannot parameterize has no AMBER topology to
            # convert, let alone to append solvent parameters to.
            refusals.append(
                f"lipid_ff={self.lipid_ff!r} is not parameterized by PACKMOL-Memgen"
            )

        return refusals

    def resolved_solvate_mode(self, estimated_waters: Optional[int] = None) -> str:
        """Return the solvation backend to use: ``"packmol"`` or ``"gromacs"``.

        ``estimated_waters`` is the projected water count for this system; it is
        only consulted for ``solvate="auto"``, and then only when the bilayer
        did not already force the answer (see below).  When it is unknown,
        ``auto`` resolves to ``packmol``: an unmeasured system is not evidence
        that the fast path is needed, so the answer falls back to the route that
        does not depend on having measured anything.  That is why ``auto`` is a
        safe default even though it is a decision -- the branch it takes without
        evidence is the one that was the default before it.

        On the patch route the answer is forced: PACKMOL-Memgen packs the lipids
        and the water in a single run, so a build that never invokes it has no
        PACKMOL solvation step to fall back to.  That makes ``solvate="packmol"``
        and a resolved patch bilayer contradictory, and the contradiction is
        never resolved by quietly substituting one for the other:

        * ``bilayer_source="auto"`` counts an explicit ``solvate="packmol"`` as a
          reason the patch route is unavailable (see :meth:`bilayer_refusals`),
          so it downgrades to PACKMOL packing and the requested PACKMOL
          solvation is honoured exactly as asked;
        * ``bilayer_source="patch"`` with ``solvate="packmol"`` is a fatal
          configuration error in :meth:`validation_errors`, and raises here for
          a config that was built directly and never validated.

        So a returned ``"gromacs"`` on the patch route only ever means the caller
        asked for GROMACS solvation or left the choice open -- this method
        enforces that rather than assuming it, because the previous version
        asserted the same invariant in prose while ``auto`` walked straight past
        it and silently moved packmol-solvation requests onto ``gmx solvate``.
        """
        mode = _require_choice(self.solvate, "solvation mode", SOLVATE_MODES)

        if self.resolved_bilayer_source() == "patch":
            if mode == "packmol":
                raise ValueError(
                    "bilayer_source='patch' is incompatible with solvate='packmol': the patch "
                    "route does not run PACKMOL-Memgen, so there is nothing to place the water "
                    "and ions. Use solvate='gromacs' (or 'auto')."
                )
            return "gromacs"

        if mode == "auto":
            if estimated_waters is None:
                mode = "packmol"
            else:
                if isinstance(estimated_waters, bool) or not isinstance(estimated_waters, Real):
                    raise ValueError("Estimated water count must be a number.")
                waters = float(estimated_waters)
                if not math.isfinite(waters) or waters < 0:
                    raise ValueError("Estimated water count must be finite and non-negative.")
                threshold = _require_positive_int(
                    self.solvate_auto_threshold, "auto-solvation threshold"
                )
                mode = "gromacs" if waters > threshold else "packmol"

        # A refusal outranks the request.  Downgrading (rather than raising) is
        # deliberate: the user asked for a membrane, and the PACKMOL path still
        # produces the correct answer.  The caller is expected to surface
        # ``solvate_refusals()`` so the downgrade is never invisible.
        if mode == "gromacs" and self.solvate_refusals():
            return "packmol"
        return mode

    def validation_errors(self) -> List[str]:
        """Return fatal configuration errors that must block external tools."""
        problems = []

        if type(self.enabled) is not bool:
            problems.append("Membrane enabled must be a boolean (true or false).")
        if type(self.minimize) is not bool:
            problems.append("Membrane minimize must be a boolean (true or false).")

        try:
            lipids = _require_lipids(self.lipids)
        except ValueError as exc:
            problems.append(str(exc))
            lipids = []

        try:
            orient = _require_choice(self.orient, "orientation method", ORIENTATION_METHODS)
        except ValueError as exc:
            problems.append(str(exc))
            orient = ""

        try:
            # Checked for its own sake.  How it interacts with the bilayer route
            # is decided in one place, ``bilayer_refusals()``, and reported below.
            _require_choice(self.solvate, "solvation mode", SOLVATE_MODES)
        except ValueError as exc:
            problems.append(str(exc))
        try:
            bilayer_source = _require_choice(
                self.bilayer_source, "bilayer source", BILAYER_SOURCES
            )
        except ValueError as exc:
            problems.append(str(exc))
            bilayer_source = ""
        try:
            _require_positive_int(self.solvate_auto_threshold, "auto-solvation threshold")
        except ValueError as exc:
            problems.append(str(exc))

        # An explicit patch request for a composition that has no patch is fatal
        # here, which is the whole reason this check exists at config time: the
        # alternative is discovering it after orientation, ligand
        # parameterisation and a tleap run have already been paid for.
        if bilayer_source == "patch":
            for refusal in self.bilayer_refusals():
                problems.append(
                    f"bilayer_source='patch' cannot be used because {refusal}; "
                    "use bilayer_source='packmol' (or 'auto', which falls back to it)."
                )
            # PACKMOL-Memgen places lipids and water in one run, so asking it for
            # the water on a route that never runs it is asking a tool that is
            # not in the pipeline.  That is one of ``bilayer_refusals()`` above --
            # stated once, so the two axes cannot be checked here and resolved
            # differently in ``resolved_bilayer_source``.
            #
            # ...and if GROMACS solvation is refused for this composition too,
            # then no backend can place the aqueous phase at all.
            for refusal in self.solvate_refusals():
                problems.append(
                    "bilayer_source='patch' requires GROMACS solvation, but that is unavailable "
                    f"because {refusal}. Use bilayer_source='packmol'."
                )
        if orient in UNAUTOMATED_ORIENTATION_METHODS:
            problems.append(UNAUTOMATED_ORIENTATION_METHODS[orient])

        lipid_ff = self.lipid_ff.strip().lower() if isinstance(self.lipid_ff, str) else ""
        if lipid_ff not in LIPID_FF_FAMILY:
            problems.append(f"Unknown lipid force field {self.lipid_ff!r}.")
        elif lipid_ff not in PACKMOL_PARAMETRIZED_LIPID_FFS:
            problems.append(
                "The automated membrane backend does not support "
                f"lipid_ff={self.lipid_ff!r}; use lipid17 or lipid21."
            )

        if orient == "opm" and (not isinstance(self.pdb_id, str) or not self.pdb_id.strip()):
            problems.append("orient=opm requires a non-empty membrane PDB id.")
        elif self.pdb_id is not None and (not isinstance(self.pdb_id, str) or not self.pdb_id.strip()):
            problems.append("Membrane PDB id must be a non-empty string when provided.")

        try:
            ratio = _require_ratio(self.ratio)
            if ratio and lipids and len(ratio) != len(lipids):
                problems.append(
                    f"Lipid ratio count ({len(ratio)}) must match lipid count ({len(lipids)})."
                )
        except ValueError as exc:
            problems.append(str(exc))

        for value, field_name, allow_zero in (
            (self.salt_concentration, "salt concentration", True),
            (self.water_thickness, "water thickness", False),
            (self.xy_distance, "protein-to-box distance", False),
            (self.temperature, "temperature", False),
            (self.center_tolerance, "centering tolerance", False),
            (self.packmol_timeout, "packmol timeout", False),
            (self.posres_fc_backbone, "backbone restraint force constant", True),
            (self.posres_fc_sidechain, "sidechain restraint force constant", True),
            (self.posres_fc_lipid, "lipid restraint force constant", True),
            (self.vdwradii_carbon, "carbon van der Waals radius", False),
            (self.patch_clash_distance, "patch clash distance", False),
            (self.patch_seam_clash_distance, "patch seam clash distance", False),
        ):
            try:
                _require_finite_number(value, field_name, allow_zero=allow_zero)
            except ValueError as exc:
                problems.append(str(exc))

        # A clash radius this large deletes lipids far beyond the solute's own
        # footprint, leaving a vacuum annulus at the protein-lipid interface that
        # equilibration will not close.
        clash = _finite_or_none(self.patch_clash_distance)
        if clash is not None and clash > PATCH_CLASH_MAX_A:
            problems.append(
                "Membrane patch clash distance must not exceed "
                f"{PATCH_CLASH_MAX_A:g} A; a wider cut leaves a vacuum gap around the solute "
                "instead of a lipid-packed interface."
            )
        seam = _finite_or_none(self.patch_seam_clash_distance)
        if seam is not None and seam > PATCH_CLASH_MAX_A:
            problems.append(
                "Membrane patch seam clash distance must not exceed "
                f"{PATCH_CLASH_MAX_A:g} A; at that separation every lipid in the patch counts as "
                "clashing with its neighbours and the bilayer would be deleted."
            )

        # A carbon radius this large carves a >1.1 nm dead shell around every
        # solute carbon, so ``gmx solvate`` would leave the aqueous phase
        # grossly under-filled -- a wrong system, not a conservative one.
        carbon_radius = _finite_or_none(self.vdwradii_carbon)
        if carbon_radius is not None and carbon_radius > VDWRADII_CARBON_MAX_NM:
            problems.append(
                "Membrane carbon van der Waals radius must not exceed "
                f"{VDWRADII_CARBON_MAX_NM:g} nm; larger radii starve the water box."
            )

        for value, field_name in (
            (self.validate_orientation, "orientation validation"),
            (self.gaff2, "gaff2"),
            (self.core_water_purge, "core water purge"),
            (self.preflight_energy, "preflight energy check"),
        ):
            if type(value) is not bool:
                problems.append(f"Membrane {field_name} must be a boolean (true or false).")

        # ``None`` means "let the patch module size it from the bundled patch's
        # periodicity"; a supplied value still has to be a real distance.
        if self.patch_xy_margin is not None:
            try:
                _require_finite_number(self.patch_xy_margin, "patch xy margin", allow_zero=True)
            except ValueError as exc:
                problems.append(str(exc))

        # Both parsers already normalise this, so reaching a bad value here means
        # the config was constructed directly -- which is exactly the path that
        # would otherwise carry ``patch_shrink_factor=1.0`` into a build that
        # reports success having done nothing.
        try:
            _require_shrink_factor(self.patch_shrink_factor)
        except ValueError as exc:
            problems.append(str(exc))

        # Ligand parameters are only meaningful as a matched frcmod/lib pair:
        # PACKMOL-Memgen's --ligand_param takes both or neither.
        supplied = [name for name, value in (
            ("ligand_frcmod", self.ligand_frcmod), ("ligand_lib", self.ligand_lib)
        ) if value is not None]
        if len(supplied) == 1:
            problems.append(
                f"Membrane ligand parameters need both ligand_frcmod and ligand_lib; got only {supplied[0]}."
            )
        for value, field_name in (
            (self.ligand_frcmod, "ligand frcmod"),
            (self.ligand_lib, "ligand lib"),
        ):
            if value is not None and (not isinstance(value, str) or not value.strip()):
                problems.append(f"Membrane {field_name} must be a non-empty path.")
        for value, field_name in (
            (self.packing_nloop, "packing nloop"),
            (self.packing_nloop_all, "all-together packing nloop"),
        ):
            if value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
                problems.append(f"Membrane {field_name} must be a positive integer.")
            elif int(value) <= PACKMOL_WRITEOUT_FREQUENCY:
                # PACKMOL-Memgen rejects a budget that is not larger than its
                # write-out frequency, but only after orientation and ligand
                # parameterization have already run. Catching it here saves
                # that work.
                problems.append(
                    f"Membrane {field_name} must exceed PACKMOL-Memgen's write-out frequency "
                    f"({PACKMOL_WRITEOUT_FREQUENCY}); use at least {PACKMOL_WRITEOUT_FREQUENCY + 1}."
                )

        if self.ligand_frcmod or self.ligand_lib:
            try:
                _require_nonempty_string(self.ligand_resname, "ligand residue name")
            except ValueError as exc:
                problems.append(str(exc))

        try:
            membrane_production_steps(self.production_ns)
        except ValueError as exc:
            problems.append(str(exc))

        for value, field_name in (
            (self.positive_ion, "positive ion"),
            (self.negative_ion, "negative ion"),
        ):
            try:
                name = _require_nonempty_string(value, field_name)
            except ValueError as exc:
                problems.append(str(exc))
                continue
            # An unrecognised spelling is almost always a typo, and left alone it
            # surfaces as an opaque tleap message after the pack.  A name PRISM
            # *recognises* but cannot carry end to end is a different matter and
            # is only warned about -- see validation_warnings.
            if name.upper() not in RECOGNISED_ION_NAMES:
                problems.append(
                    f"Membrane {field_name} {name!r} is not a recognised ion name. "
                    "Known names: " + ", ".join(sorted(RECOGNISED_ION_NAMES)) + "."
                )

        problems.extend(self._disulfide_pair_errors())

        if lipid_ff in PACKMOL_PARAMETRIZED_LIPID_FFS:
            try:
                self.packmol_forcefield_pair()
            except ValueError as exc:
                problems.append(str(exc))
        return problems

    def _disulfide_pair_errors(self) -> List[str]:
        """Type/range problems in :attr:`disulfide_pairs`.

        Checked here rather than only where the ``bond`` lines are rendered so a
        malformed override is reported alongside every other configuration
        problem, before any packing happens.  The rendering side validates again
        -- it also sees the auto-detected pairs, which never come through here.
        """
        pairs = self.disulfide_pairs
        if pairs is None:
            return []
        if isinstance(pairs, (str, bytes)) or not isinstance(pairs, (list, tuple)):
            return [
                "Membrane disulfide_pairs must be a list of (residue, residue) "
                f"number pairs; got {type(pairs).__name__}."
            ]
        problems = []
        for pair in pairs:
            if isinstance(pair, (str, bytes)) or not isinstance(pair, (list, tuple)):
                problems.append(
                    f"Membrane disulfide_pairs entry must be a pair of residue numbers: {pair!r}"
                )
                continue
            if len(pair) != 2:
                problems.append(
                    f"Membrane disulfide_pairs entry must hold exactly two residue "
                    f"numbers: {pair!r}"
                )
                continue
            if any(isinstance(v, bool) or not isinstance(v, Integral) or int(v) < 1 for v in pair):
                problems.append(
                    f"Membrane disulfide_pairs residue numbers must be positive integers: {pair!r}"
                )
                continue
            if int(pair[0]) == int(pair[1]):
                problems.append(
                    f"Membrane disulfide_pairs entry bonds a residue to itself: {pair!r}"
                )
        return problems

    def validation_warnings(self) -> List[str]:
        """Return advisory messages which do not make a build scientifically invalid."""
        warnings = []

        # Surfaced here, at config time, rather than left to fail after the
        # bilayer has been packed and parametrized.  See
        # GROMACS_ROUTE_SUPPORTED_IONS for the measurement behind this list.
        for value, field_name in (
            (self.positive_ion, "positive ion"),
            (self.negative_ion, "negative ion"),
        ):
            if not isinstance(value, str) or not value.strip():
                continue
            name = value.strip()
            if name.upper() in RECOGNISED_ION_NAMES and (
                name.upper() not in GROMACS_ROUTE_SUPPORTED_IONS
            ):
                warnings.append(
                    f"Membrane {field_name} '{name}' cannot be carried through the GROMACS "
                    "solvation route: the ion unit PRISM asks tleap for does not exist in "
                    "AmberTools' atomic_ions.lib, and tleap answers with a topology that "
                    "silently omits the ion (this is caught later, but only after the "
                    "bilayer has been packed and parametrized). Ions verified end to end "
                    "are sodium, potassium and chloride "
                    f"({', '.join(sorted(GROMACS_ROUTE_SUPPORTED_IONS))}). "
                    "A wet PACKMOL-Memgen build may still succeed with this ion."
                )

        if isinstance(self.lipids, (list, tuple)):
            for lipid in self.lipids:
                if isinstance(lipid, str) and lipid.strip().upper() not in COMMON_LIPIDS:
                    warnings.append(
                        f"Lipid '{lipid}' is not in the common-lipid whitelist "
                        "(it may still be valid for PACKMOL-Memgen)."
                    )

        requested_orient = self.orient.strip().lower() if isinstance(self.orient, str) else ""
        if requested_orient == "auto" and self.validate_orientation is False:
            # The one combination in which the ``auto`` default is not a
            # superset of the ``preoriented`` one it replaced.
            #
            # ``preoriented`` + ``validate_orientation=False`` was a promise in
            # both directions: PRISM neither checked the frame nor touched it.
            # ``auto`` still measures -- it has to, because the measurement is
            # how it picks a route -- so an input this module judges un-framed
            # is now handed to MEMEMBED, and a deliberately prepared frame that
            # PRISM measures differently (a monotopic or peripheral protein is
            # the realistic case, since the membrane-spanning region it looks
            # for is not there to find) gets reoriented instead of preserved.
            #
            # Not an error: the reorientation is the right answer far more often
            # than it is wrong, and the build still succeeds.  But it is the one
            # way the new default can hand back coordinates the old one would
            # have left alone, so it is said out loud rather than discovered by
            # diffing the output.
            warnings.append(
                "orient='auto' with validate_orientation=false: 'auto' measures the input "
                "against the membrane frame in order to choose its route, so switching "
                "validation off does not stop it. If PRISM judges this structure un-framed it "
                "will be reoriented with MEMEMBED rather than passed through. To keep the "
                "coordinates exactly as prepared, use orient='preoriented'."
            )

        requested_mode = self.solvate.strip().lower() if isinstance(self.solvate, str) else ""

        if self.core_water_purge is False:
            # The purge is what makes GROMACS solvation of a bilayer defensible.
            # Without it the system ships with water in the hydrophobic core --
            # a result, not a nuisance -- so say so in those words.
            warnings.append(
                "core_water_purge is disabled: with solvate='gromacs' the finished system "
                "will contain water inside the hydrophobic core of the bilayer, which is "
                "physically wrong. Leave it enabled unless you are deliberately measuring "
                "what the purge removes."
            )

        carbon_radius = _finite_or_none(self.vdwradii_carbon)
        if carbon_radius is not None and carbon_radius <= GROMACS_DEFAULT_CARBON_RADIUS_NM:
            warnings.append(
                f"vdwradii_carbon={carbon_radius:g} nm is at or below stock GROMACS "
                f"({GROMACS_DEFAULT_CARBON_RADIUS_NM:g} nm), so 'gmx solvate' will flood the "
                f"acyl-tail region and the core-water purge has to remove all of it; "
                f"{DEFAULT_VDWRADII_CARBON_NM:g} nm is the tutorial value."
            )
        elif carbon_radius is not None and carbon_radius > VDWRADII_CARBON_BULK_WARNING_NM:
            warnings.append(
                f"vdwradii_carbon={carbon_radius:g} nm is well above the tutorial value "
                f"({DEFAULT_VDWRADII_CARBON_NM:g} nm); such a radius starts to under-fill the "
                "bulk water without removing the last of the core water."
            )

        if requested_mode in ("gromacs", "auto"):
            for refusal in self.solvate_refusals():
                warnings.append(
                    f"solvate={requested_mode!r} falls back to PACKMOL solvation because "
                    f"{refusal}."
                )

        requested_source = (
            self.bilayer_source.strip().lower() if isinstance(self.bilayer_source, str) else ""
        )
        if requested_source == "auto":
            # 'patch' with a refusal is fatal (see validation_errors); 'auto' with
            # a refusal is the downgrade this warning exists to make visible.
            for refusal in self.bilayer_refusals():
                warnings.append(
                    f"bilayer_source='auto' falls back to PACKMOL packing because {refusal}."
                )
        if requested_source in ("patch", "auto") and not self.bilayer_refusals():
            # A *plan*, in the future tense, carrying no citation.
            #
            # This used to assert the provenance outright -- "the bilayer will be
            # tiled from a bundled patch (Zenodo doi:...); cite it" -- from the
            # configuration alone, before anything knew which route would run.
            # The builder copies these warnings into its result verbatim, so a
            # finished PACKMOL-route build (the patch route can still decline at
            # build time: missing assets, a failed checksum, a capability gate)
            # carried a false statement about where its lipids came from and an
            # instruction to cite a dataset it had never opened.  A citation
            # belongs to what ran, so it now lives in
            # :func:`bilayer_patch_citation`, which is past tense and is the
            # caller's to emit once the route has actually produced a bilayer.
            #
            # The phrase "will be" is load-bearing, not style.  The builder
            # suppresses a stale plan by matching this message on
            # "patch" + "bilayer" + "will be" (MembraneBuilder._claims_patch_route)
            # before copying it into a result whose route turned out to be
            # PACKMOL.  An earlier rewording dropped the phrase, so the filter
            # stopped matching and the plan survived into builds that fell back.
            # Keep all three tokens in any future rewording of this sentence.
            #
            # Note that under the ``auto`` default this fires on every build of a
            # composition that has a bundled patch -- which is the point (the
            # lipids really do come from someone else's equilibration, and that
            # has to be disclosed), but it does mean a non-empty ``warnings`` is
            # no longer a useful health signal on the default path.  The right
            # home for it is the builder's ``info`` channel plus a post-build
            # :func:`bilayer_patch_citation`; that is the caller's call to make.
            warnings.append(
                f"bilayer_source={requested_source!r} selects the bundled pre-equilibrated "
                "patch route for this composition: the bilayer will be tiled from third-party "
                "equilibrated coordinates rather than packed for this system, and GROMACS will "
                "place the water and ions. Those coordinates are under a citation-required "
                "licence. Cite what the finished build reports rather than this plan -- the "
                "route can still fall back to packing at build time."
            )

        # ``xy_distance`` is PACKMOL-Memgen's ``--dist`` and nothing else reads
        # it.  Before ``bilayer_source`` defaulted to ``auto`` that was fine --
        # the default route was the one that used it.  Now a plain build of a
        # composition with a bundled patch tiles instead, and a user who widened
        # the in-plane margin would get the patch's periodicity regardless, with
        # nothing anywhere saying their number had been ignored.  Only a
        # *non-default* value is reported, so the warning marks a request that
        # was actually made rather than firing on every build.
        xy_distance = _finite_or_none(self.xy_distance)
        default_xy_distance = membrane_defaults()["xy_distance"]
        if xy_distance is not None and xy_distance != default_xy_distance:
            try:
                tiles = self.resolved_bilayer_source() == "patch"
            except ValueError:
                tiles = False
            if tiles:
                warnings.append(
                    f"xy_distance={xy_distance:g} A is a PACKMOL-Memgen setting (--dist) and this "
                    "build tiles a bundled patch instead, so it will not be applied: the tiled "
                    "bilayer is sized from the patch's own periodicity. Use patch_xy_margin to "
                    "widen the in-plane margin on the patch route, or bilayer_source='packmol' "
                    "to get the packing route this value belongs to."
                )

        clash = _finite_or_none(self.patch_clash_distance)
        seam = _finite_or_none(self.patch_seam_clash_distance)
        if clash is not None and clash < PATCH_CLASH_MIN_WARNING_A:
            warnings.append(
                f"patch_clash_distance={clash:g} A is below a heavy-atom contact distance "
                f"({PATCH_CLASH_MIN_WARNING_A:g} A), so lipids overlapping the solute will "
                f"survive the cut; {DEFAULT_PATCH_CLASH_A:g} A is the validated value."
            )
        if clash is not None and seam is not None and seam > clash:
            # The patch is already equilibrated, so its own lipid-lipid contacts
            # are correct structure.  A seam criterion looser than the
            # protein-lipid one deletes that structure at every tile boundary.
            warnings.append(
                f"patch_seam_clash_distance={seam:g} A exceeds patch_clash_distance={clash:g} A: "
                "the patch's own equilibrated lipid-lipid contacts would be treated as clashes "
                "and deleted at every tile seam."
            )

        try:
            shrink = _require_shrink_factor(self.patch_shrink_factor)
        except ValueError:
            # Already fatal in validation_errors, with the full explanation.
            shrink = None
        if isinstance(shrink, float):
            low, high = SHRINK_FACTOR_WARNING_BOUNDS
            if shrink < low:
                warnings.append(
                    f"patch_shrink_factor={shrink:g} is below the range PRISM solves within "
                    f"({low:g}-{high:g}): the lipids are deleted against so small a copy of "
                    "the solute that the survivors end up deep inside its real footprint, and "
                    "the grow-back has to drag them out through it. Use 'auto' to size the "
                    "shrink from the area budget instead."
                )
            elif shrink > high:
                warnings.append(
                    f"patch_shrink_factor={shrink:g} is above the range PRISM solves within "
                    f"({low:g}-{high:g}): a shrink this slight deletes very nearly the lipids "
                    "plain deletion would, so the grow-back's minimisations are paid for and "
                    "buy almost nothing."
                )
        if shrink is not None:
            # Same shape as the xy_distance warning above, and the same reason:
            # a knob that belongs to one route should say so when the build is
            # taking the other, rather than being dropped in silence.
            try:
                packs = self.resolved_bilayer_source() == "packmol"
            except ValueError:
                packs = False
            if packs:
                warnings.append(
                    f"patch_shrink_factor={shrink!r} applies to the patch route and this build "
                    "packs its bilayer with PACKMOL-Memgen instead, so it will not be applied. "
                    "PACKMOL packs lipids around the solute directly and has no tiled lattice "
                    "to close in, which is the gap the shrink-and-grow exists to close."
                )

        return warnings

    def raise_for_errors(self) -> None:
        errors = self.validation_errors()
        if errors:
            raise ValueError("Invalid membrane configuration: " + "; ".join(errors))


def membrane_defaults() -> dict:
    """Return the default value of every :class:`MembraneConfig` field.

    Both parsers read their fallbacks from here (and a generated YAML template
    can too), so each default is written once -- in the dataclass -- instead of
    being hand-copied into every parser, where the copies drift.
    """
    defaults = {}
    for spec in fields(MembraneConfig):
        # Containers are rebuilt per call: the caller may hand one to a config.
        if spec.default_factory is not MISSING:
            defaults[spec.name] = spec.default_factory()
        else:
            defaults[spec.name] = spec.default
    return defaults


def gaff2_selected(ligand_forcefield: Optional[str]) -> bool:
    """True when the ligand force field is GAFF2 rather than GAFF.

    PACKMOL-Memgen's ``--gaff2`` picks a leaprc, so it has to follow the force
    field that actually typed the ligand; anything else (including "unknown")
    leaves PACKMOL-Memgen on its documented ``leaprc.gaff`` default.
    """
    return isinstance(ligand_forcefield, str) and ligand_forcefield.strip().lower() == "gaff2"


def parse_membrane_cli(
    membrane: bool,
    lipid: Optional[List[str]],
    lipid_ratio: Optional[List[int]],
    orient: Optional[str],
    pdb_id: Optional[str],
    lipid_ff: Optional[str],
    salt_concentration: Optional[float] = None,
    temperature: Optional[float] = None,
    production_ns: Optional[float] = None,
    positive_ion: Optional[str] = None,
    negative_ion: Optional[str] = None,
    protein_forcefield: Optional[str] = None,
    water_model: Optional[str] = None,
    ligand_forcefield: Optional[str] = None,
    solvate: Optional[str] = None,
    core_water_purge: Optional[bool] = None,
    vdwradii_carbon: Optional[float] = None,
    bilayer_source: Optional[str] = None,
    patch_clash_distance: Optional[float] = None,
    patch_seam_clash_distance: Optional[float] = None,
    patch_xy_margin: Optional[float] = None,
    patch_shrink_factor: Optional[float | str] = None,
) -> MembraneConfig:
    """Build a :class:`MembraneConfig` from CLI arguments.

    ``ligand_forcefield`` is the caller's ``-lff`` selection. It must be passed
    whenever the build can carry a ligand, because it decides which GAFF
    generation PACKMOL-Memgen parameterizes that ligand with.

    ``solvate``/``core_water_purge``/``vdwradii_carbon`` back
    ``--membrane-solvate``, ``--membrane-no-core-purge`` and
    ``--membrane-vdw-carbon``; ``bilayer_source``/``patch_clash_distance``/
    ``patch_seam_clash_distance`` back ``--membrane-bilayer``,
    ``--membrane-patch-clash`` and ``--membrane-patch-seam-clash``.  Every one of
    them is ``None`` when the flag was not given, and ``None`` means "take the
    dataclass default" -- it does not mean "keep the old behaviour".

    That distinction used to be moot and is not any more.  ``orient``,
    ``bilayer_source`` and ``solvate`` all default to ``auto`` now, so a caller
    that passes ``None`` for them is asking for the measured route rather than
    for PACKMOL packing and a ``preoriented`` assertion.  A caller that wants
    the pre-``auto`` behaviour has to name it: ``orient="preoriented"``,
    ``bilayer_source="packmol"``, ``solvate="packmol"``.
    """
    defaults = membrane_defaults()

    def _or_default(value, field_name):
        return defaults[field_name] if value is None else value

    config = MembraneConfig(
        enabled=_require_bool(membrane, "enabled"),
        lipids=_require_lipids(lipid or defaults["lipids"]),
        ratio=_require_ratio(lipid_ratio),
        orient=_require_choice(
            orient or defaults["orient"], "orientation method", ORIENTATION_METHODS
        ),
        pdb_id=pdb_id.strip() if isinstance(pdb_id, str) else pdb_id,
        lipid_ff=_require_choice(
            lipid_ff or defaults["lipid_ff"], "lipid force field", LIPID_FF_FAMILY
        ),
        salt_concentration=_require_finite_number(
            _or_default(salt_concentration, "salt_concentration"),
            "salt concentration",
            allow_zero=True,
        ),
        positive_ion=_require_nonempty_string(
            _or_default(positive_ion, "positive_ion"), "positive ion"
        ),
        negative_ion=_require_nonempty_string(
            _or_default(negative_ion, "negative_ion"), "negative ion"
        ),
        protein_forcefield=_require_nonempty_string(
            _or_default(protein_forcefield, "protein_forcefield"), "protein force field"
        ),
        water_model=_require_nonempty_string(_or_default(water_model, "water_model"), "water model"),
        temperature=_require_finite_number(_or_default(temperature, "temperature"), "temperature"),
        production_ns=_require_finite_number(
            _or_default(production_ns, "production_ns"), "production duration"
        ),
        gaff2=gaff2_selected(ligand_forcefield),
        solvate=_require_choice(
            _or_default(solvate, "solvate"), "solvation mode", SOLVATE_MODES
        ),
        core_water_purge=_require_bool(
            _or_default(core_water_purge, "core_water_purge"), "core water purge"
        ),
        vdwradii_carbon=_require_finite_number(
            _or_default(vdwradii_carbon, "vdwradii_carbon"), "carbon van der Waals radius"
        ),
        bilayer_source=_require_choice(
            _or_default(bilayer_source, "bilayer_source"), "bilayer source", BILAYER_SOURCES
        ),
        patch_clash_distance=_require_finite_number(
            _or_default(patch_clash_distance, "patch_clash_distance"), "patch clash distance"
        ),
        patch_seam_clash_distance=_require_finite_number(
            _or_default(patch_seam_clash_distance, "patch_seam_clash_distance"),
            "patch seam clash distance",
        ),
        # Passed through unvalidated: ``validation_errors`` owns the "None means
        # the patch module decides" case, exactly as it does for packing_nloop.
        patch_xy_margin=_or_default(patch_xy_margin, "patch_xy_margin"),
        # Normalised here (not merely passed through) because this field has a
        # *vocabulary* -- 'auto', and the off-spellings -- and normalising it at
        # the boundary is what ``solvate`` and ``bilayer_source`` do too.  The
        # ``None`` is argparse's "flag not given", which is why it goes through
        # ``_or_default`` rather than straight in: the dataclass default is
        # ``"auto"``, so passing the sentinel on would turn the shrink-and-grow
        # off for every CLI build while the YAML path kept it on.
        patch_shrink_factor=_require_shrink_factor(
            _or_default(patch_shrink_factor, "patch_shrink_factor")
        ),
    )
    config.raise_for_errors()
    return config


def parse_membrane_yaml(
    cfg: Optional[dict],
    temperature: Optional[float] = None,
    production_ns: Optional[float] = None,
) -> MembraneConfig:
    """Build a :class:`MembraneConfig` from a parsed YAML ``membrane:`` mapping."""
    if cfg is None:
        return MembraneConfig(enabled=False)
    if not isinstance(cfg, dict):
        raise ValueError("The YAML 'membrane' section must be a mapping.")

    defaults = membrane_defaults()
    if temperature is None:
        temperature = defaults["temperature"]
    if production_ns is None:
        production_ns = defaults["production_ns"]

    # Retired knob.  ``false`` asked for nothing and is accepted in silence, so
    # an old config keeps working; ``true`` asked for a specific approximation
    # that no longer exists, and being handed the exact conversion instead
    # without being told is how a config comes to describe a build it did not
    # produce.
    if "allow_mixed_14_scaling" in cfg and _require_bool(
        cfg["allow_mixed_14_scaling"], "allow mixed 1-4 scaling"
    ):
        raise ValueError(
            "Membrane allow_mixed_14_scaling has been retired. It applied one global "
            "fudgeLJ to every 1-4 pair of a system whose AMBER SCNB is not uniform; the "
            "converter now emits explicit per-pair [ pairs ] parameters with gen_pairs=no, "
            "so Lipid21's SCNB=6.0 acyl-tail torsion is converted exactly and the "
            "approximation is neither available nor needed. Remove the key."
        )

    pdb_id = cfg.get("pdb_id")
    if isinstance(pdb_id, str):
        pdb_id = pdb_id.strip()

    if "lipids" in cfg:
        lipid_value = cfg["lipids"]
    elif "lipid" in cfg:
        lipid_value = cfg["lipid"]
    else:
        lipid_value = defaults["lipids"]

    config = MembraneConfig(
        enabled=_require_bool(cfg.get("enabled", True), "enabled"),
        lipids=_require_lipids(lipid_value),
        ratio=_require_ratio(cfg.get("ratio", defaults["ratio"])),
        orient=_require_choice(
            cfg.get("orient", defaults["orient"]), "orientation method", ORIENTATION_METHODS
        ),
        pdb_id=pdb_id,
        lipid_ff=_require_choice(
            cfg.get("lipid_ff", defaults["lipid_ff"]), "lipid force field", LIPID_FF_FAMILY
        ),
        protein_forcefield=_require_nonempty_string(
            cfg.get("protein_forcefield", cfg.get("forcefield", defaults["protein_forcefield"])),
            "protein force field",
        ),
        water_model=_require_nonempty_string(
            cfg.get("water_model", cfg.get("water", defaults["water_model"])), "water model"
        ),
        salt_concentration=_require_finite_number(
            cfg.get("salt_concentration", defaults["salt_concentration"]),
            "salt concentration",
            allow_zero=True,
        ),
        positive_ion=_require_nonempty_string(
            cfg.get("positive_ion", defaults["positive_ion"]), "positive ion"
        ),
        negative_ion=_require_nonempty_string(
            cfg.get("negative_ion", defaults["negative_ion"]), "negative ion"
        ),
        water_thickness=_require_finite_number(
            cfg.get("water_thickness", defaults["water_thickness"]), "water thickness"
        ),
        xy_distance=_require_finite_number(
            cfg.get("xy_distance", defaults["xy_distance"]), "protein-to-box distance"
        ),
        minimize=_require_bool(cfg.get("minimize", defaults["minimize"]), "minimize"),
        validate_orientation=_require_bool(
            cfg.get("validate_orientation", defaults["validate_orientation"]),
            "orientation validation",
        ),
        center_tolerance=_require_finite_number(
            cfg.get("center_tolerance", defaults["center_tolerance"]), "centering tolerance"
        ),
        posres_fc_backbone=_require_finite_number(
            cfg.get("posres_fc_backbone", defaults["posres_fc_backbone"]),
            "backbone restraint force constant",
            allow_zero=True,
        ),
        posres_fc_sidechain=_require_finite_number(
            cfg.get("posres_fc_sidechain", defaults["posres_fc_sidechain"]),
            "sidechain restraint force constant",
            allow_zero=True,
        ),
        posres_fc_lipid=_require_finite_number(
            cfg.get("posres_fc_lipid", defaults["posres_fc_lipid"]),
            "lipid restraint force constant",
            allow_zero=True,
        ),
        ligand_frcmod=cfg.get("ligand_frcmod"),
        ligand_lib=cfg.get("ligand_lib"),
        ligand_resname=cfg.get("ligand_resname", defaults["ligand_resname"]),
        gaff2=_require_bool(cfg.get("gaff2", defaults["gaff2"]), "gaff2"),
        packmol_timeout=_require_finite_number(
            cfg.get("packmol_timeout", defaults["packmol_timeout"]), "packmol timeout"
        ),
        packing_nloop=cfg.get("packing_nloop", defaults["packing_nloop"]),
        packing_nloop_all=cfg.get("packing_nloop_all", defaults["packing_nloop_all"]),
        solvate=_require_choice(
            cfg.get("solvate", defaults["solvate"]), "solvation mode", SOLVATE_MODES
        ),
        bilayer_source=_require_choice(
            cfg.get("bilayer_source", defaults["bilayer_source"]),
            "bilayer source",
            BILAYER_SOURCES,
        ),
        patch_clash_distance=_require_finite_number(
            cfg.get("patch_clash_distance", defaults["patch_clash_distance"]),
            "patch clash distance",
        ),
        patch_seam_clash_distance=_require_finite_number(
            cfg.get("patch_seam_clash_distance", defaults["patch_seam_clash_distance"]),
            "patch seam clash distance",
        ),
        # Left unvalidated here on purpose: ``None`` is a meaningful value for
        # this field, and ``validation_errors`` checks it only when it is not.
        patch_xy_margin=cfg.get("patch_xy_margin", defaults["patch_xy_margin"]),
        # Unlike patch_xy_margin, this one IS normalised at the boundary: it has
        # a vocabulary ('auto', plus the off-spellings), and unquoted ``none`` in
        # YAML arrives as the string "none" rather than as a null.
        patch_shrink_factor=_require_shrink_factor(
            cfg.get("patch_shrink_factor", defaults["patch_shrink_factor"])
        ),
        # Handed over unvalidated on purpose: ``validation_errors`` owns the
        # integer check, exactly as it does for the packing nloops.
        solvate_auto_threshold=cfg.get(
            "solvate_auto_threshold", defaults["solvate_auto_threshold"]
        ),
        core_water_purge=_require_bool(
            cfg.get("core_water_purge", defaults["core_water_purge"]), "core water purge"
        ),
        vdwradii_carbon=_require_finite_number(
            cfg.get("vdwradii_carbon", defaults["vdwradii_carbon"]),
            "carbon van der Waals radius",
        ),
        preflight_energy=_require_bool(
            cfg.get("preflight_energy", defaults["preflight_energy"]), "preflight energy check"
        ),
        temperature=_require_finite_number(cfg.get("temperature", temperature), "temperature"),
        production_ns=_require_finite_number(
            cfg.get("production_ns", cfg.get("production_time_ns", production_ns)),
            "production duration",
        ),
    )
    config.raise_for_errors()
    return config
