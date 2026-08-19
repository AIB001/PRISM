"""Water/ion force-field blocks for the GROMACS membrane solvation path.

Why this module exists
----------------------
On the PACKMOL route tleap parametrises the *finished* system, so the water
and the ions arrive in the topology as a side effect of the same run that
parametrises the protein and the lipids.  The fast route deliberately packs a
**dry** system, so its prmtop contains no water and no ions -- and therefore
neither the ``[ atomtypes ]`` entries nor the ``[ moleculetype ]`` blocks that
``gmx solvate`` and ``gmx genion`` are about to reference.  ``gmx solvate``
appends ``SOL <n>`` to ``[ molecules ]`` without checking that ``SOL`` is
defined; the failure surfaces two stages later as an unreadable grompp error.

The blocks are *probed from tleap* rather than tabulated in Python.  A table
would answer an ``opc`` request with TIP3P parameters and nothing downstream
would object: the atom count matches, grompp is happy, and the run is silently
wrong.  Asking tleap for exactly the requested water model and ion pair costs
about a second and makes that class of mistake impossible -- the parameters
that reach the topology are, byte for byte, the ones the validated PACKMOL
route would have produced (verified against a full tleap-parametrised system:
the ``OW``/``HW``/``Na+``/``Cl-`` atomtype lines are identical).

What is renamed, and why
------------------------
tleap/ParmEd name the water moleculetype after its residue, ``WAT``.
``gmx solvate`` stamps down the residue name baked into the ``-cs`` box, which
is ``SOL``.  The fast path's purge, its index generation and its
``[ molecules ]`` bookkeeping all key on ``SOL``, so the water block is
renamed here -- once, at the source -- instead of being special-cased in three
later places.  Ion blocks are renamed to the exact strings the caller passes as
``-pname``/``-nname``, because ``genion`` resolves those against
``[ moleculetype ]`` names literally.

Only ``[ atomtypes ]`` and the three ``[ moleculetype ]`` blocks are merged.
``[ molecules ]`` is left alone: ``gmx solvate`` and ``gmx genion`` write their
own counts there, and a pre-seeded entry would be double-counted.

The blocks are inserted **inline** rather than via ``#include``.  That is not a
style choice: the fast path's guards parse the topology with
:func:`prism.membrane.solvate.read_topology`, which skips preprocessor lines,
so an ``#include``d moleculetype would be invisible to exactly the checks that
exist to prove it is present.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = [
    "SolventParameterError",
    "SolventBlocks",
    "solvent_parameter_blocks",
    "assert_solvent_blocks",
    "merge_solvent_blocks",
    "WATER_MOLECULETYPE",
    # AMBER 1-4 scaling -> GROMACS [ pairs ].  Public and reusable on purpose:
    # every ParmEd ``save(format="gromacs")`` in PRISM needs it, and a second
    # copy of this logic would be a second chance to get it wrong.
    "AmberFudgeError",
    "Amber14ConflictError",
    "AmberFudge",
    "derive_amber_fudge",
    "attach_explicit_14_pairs",
    "apply_amber_fudge_to_topology",
    "save_gromacs_topology",
    "topology_defaults",
    "topology_pair_census",
]

#: The moleculetype name every downstream stage expects for water.  Kept in
#: step with :data:`prism.membrane.solvate.WATER_MOLECULETYPE`; imported
#: defensively so this module still works if solvate.py cannot be imported.
try:  # pragma: no cover - trivial fallback
    from .solvate import WATER_MOLECULETYPE  # type: ignore
except Exception:  # pragma: no cover
    WATER_MOLECULETYPE = "SOL"

#: ParmEd's "one torsion, several periodicities" container.  Resolved once here
#: rather than inside the dihedral walk: that walk runs ~10^5 times per bilayer,
#: and re-executing even a cached ``from ... import`` per iteration cost more
#: than everything else in the scan put together.  ParmEd stays an optional
#: dependency -- the fallback makes ``isinstance`` simply never match, which is
#: the right answer for a structure that cannot have come from ParmEd anyway.
try:  # pragma: no cover - environment dependent
    from parmed.topologyobjects import DihedralTypeList as _DihedralTypeList
except Exception:  # pragma: no cover
    _DihedralTypeList = ()  # type: ignore[assignment,misc]

#: Water model -> the ``leaprc`` that defines its ``WAT`` unit and the matching
#: ion parameter set.  Sourcing the water-model-specific leaprc is what makes
#: the ion parameters consistent with the water: Joung-Cheatham ions are
#: parametrised *per water model*, so ``frcmod.ionsjc_tip3p`` and
#: ``frcmod.ionsjc_opc`` are different numbers for the same ion.
_WATER_LEAPRC = {
    "tip3p": "leaprc.water.tip3p",
    "tip3p_original": "leaprc.water.tip3p",
    "spce": "leaprc.water.spce",
    "spc": "leaprc.water.spce",
    "opc": "leaprc.water.opc",
    "opc3": "leaprc.water.opc3",
    "tip4pew": "leaprc.water.tip4pew",
    "fb3": "leaprc.water.fb3",
    "fb4": "leaprc.water.fb4",
}

#: Caller-supplied ion name -> the unit name tleap actually defines in
#: ``atomic_ions.lib``.  PRISM's config spells them ``Na+``/``Cl-`` (which are
#: already tleap units), but CHARMM-style (``SOD``/``CLA``) and bare
#: (``NA``/``CL``) spellings reach here from user configs and YAML.
_ION_UNITS = {
    "NA": "Na+", "NA+": "Na+", "SOD": "Na+", "SODIUM": "Na+",
    "K": "K+", "K+": "K+", "POT": "K+", "POTASSIUM": "K+",
    "LI": "Li+", "LI+": "Li+",
    "RB": "Rb+", "RB+": "Rb+",
    "CS": "Cs+", "CS+": "Cs+",
    "MG": "Mg+", "MG2+": "Mg+",
    "CA": "Ca2+", "CA2+": "Ca2+", "CAL": "Ca2+",
    "ZN": "Zn2+", "ZN2+": "Zn2+",
    "CL": "Cl-", "CL-": "Cl-", "CLA": "Cl-", "CHLORIDE": "Cl-",
    "BR": "Br-", "BR-": "Br-",
    "I": "I-", "I-": "I-", "IOD": "I-",
    "F": "F-", "F-": "F-",
}

_DIRECTIVE_RE = re.compile(r"^\s*\[\s*([A-Za-z_0-9]+)\s*\]")

#: Relative tolerance when deciding whether an atomtype already present in the
#: host topology is *the same* type as the probed one.  ParmEd writes both, so
#: they agree to the last digit in practice; the tolerance only absorbs a
#: formatting difference between AmberTools releases.
_PARAM_RTOL = 1e-6

_PACKMOL_HINT = (
    "Re-run with --membrane-solvate packmol to use the validated "
    "PACKMOL-Memgen solvation path instead."
)


class SolventParameterError(RuntimeError):
    """Raised when the solvent parameters cannot be probed or merged."""


@dataclass(frozen=True)
class SolventBlocks:
    """Force-field text for water + the two ions, ready to merge.

    Attributes are the names :mod:`prism.membrane.builder` reports and
    :mod:`prism.membrane.solvate` checks against; the two text tuples are the
    payload :func:`merge_solvent_blocks` splices into the topology.
    """

    water_moleculetype: str
    water_sites: int
    cation_moleculetype: str
    anion_moleculetype: str
    #: (type name, verbatim ``[ atomtypes ]`` line)
    atomtypes: Tuple[Tuple[str, str], ...]
    #: (moleculetype name, full block text including the ``[ moleculetype ]``
    #: directive and everything up to the next one)
    moleculetypes: Tuple[Tuple[str, str], ...]
    #: Combination rule of the probe topology; must match the host's.
    comb_rule: int
    water_model: str
    #: Net charge of each moleculetype, used by :func:`assert_solvent_blocks`.
    charges: Dict[str, float] = field(default_factory=dict)
    source: str = ""


# --------------------------------------------------------------------------- #
# AMBER 1-4 scaling -> GROMACS [ pairs ]
# --------------------------------------------------------------------------- #
#
# AMBER scales 1-4 non-bonded interactions down by SCNB (Lennard-Jones) and
# SCEE (Coulomb), both stored *per dihedral type* in the prmtop.  GROMACS'
# ``gen-pairs = yes`` has no per-dihedral equivalent: it generates every 1-4 LJ
# interaction from the atom types and scales all of them by the single global
# ``fudgeLJ`` in ``[ defaults ]``.
#
# ParmEd 4.3.1 fills that field asymmetrically.  ``GromacsTopologyFile.
# from_structure`` derives ``fudgeQQ`` from SCEE and ``fudgeLJ`` from SCNB, but
# only when the value is unique; when the set of SCNB values is empty *or has
# more than one member* it silently falls back to ``fudgeLJ = 1.0`` -- while the
# same non-uniformity in SCEE raises ``GromacsError``.  ``1.0`` means "do not
# scale", so every 1-4 LJ interaction comes out SCNB times too strong.  grompp
# accepts it without a murmur: 1.0 is a legal fudge factor, merely the wrong one.
#
# **AMBER Lipid21 is non-uniform by design.**  Its saturated-tail torsion
# ``cD-cD-cD-cD`` carries ``SCNB=6.0`` (``dat/leap/parm/lipid21.dat`` lines
# 118-122) while every other torsion uses the AMBER-standard ``SCNB=2.0``.  So
# *no* global ``fudgeLJ`` is correct for any Lipid21 bilayer: ``1.0`` makes the
# 2.0 pairs twice and the 6.0 pairs six times too strong, and ``0.5`` fixes the
# majority while leaving the tails three times too strong.
#
# The fix is therefore not a better global constant -- it is to stop using one.
# :func:`attach_explicit_14_pairs` writes the *combined* 1-4 LJ parameters onto
# every ``[ pairs ]`` row and sets ``gen-pairs = no``, which is the only GROMACS
# construct that can express per-pair scaling:
#
#     [ pairs ]
#     ;   ai     aj funct       sigma1-4      epsilon1-4
#          44     50     1    0.300054688     0.069194303
#
# ``epsilon1-4`` already carries the division by that pair's own SCNB, so the
# 2.0 and the 6.0 pairs coexist in one topology with no global compromise.  The
# file stays small because GROMACS de-duplicates per ``[ moleculetype ]``: one
# POPC's worth of rows is written once and reused by ``[ molecules ]``.
#
# ``fudgeQQ`` cannot be handled the same way and does not need to be.  A funct-1
# ``[ pairs ]`` row carries sigma and epsilon only -- there is no per-pair charge
# scale column -- so the Coulomb 1-4 stays global.  Every AMBER force field in
# use here (protein, lipid, GAFF, water) sets SCEE=1.2 throughout, so a global
# ``fudgeQQ = 1/1.2`` is exact; non-uniform SCEE is still refused, loudly,
# because it genuinely cannot be represented.
#
#
# THE SCNB CONFLICT, AND WHY IT IS NOT A JUDGEMENT CALL
# ----------------------------------------------------
# One 1-4 atom pair can be reached by several dihedrals -- through a ring, or
# through the multiple Fourier terms AMBER stores as separate dihedral entries.
# Emitting one ``[ pairs ]`` row per dihedral would double-count the interaction,
# so exactly one dihedral per pair must own it.  If two candidate dihedrals
# disagreed about SCNB, whichever one happened to be listed first in the prmtop
# would silently decide the physics.
#
# PRISM does not invent an arbitration rule for that, because AMBER already has
# one and it is recorded in the file: the prmtop's per-dihedral "end-group
# interactions ignored" flag, which ParmEd exposes as ``Dihedral.ignore_end`` and
# recomputes in ``Structure.update_dihedral_exclusions()``.  That method walks
# the dihedrals in order and sets ``ignore_end = True`` on every one whose
# ``(atom1, atom4)`` pair has already been claimed.  The dihedral left with
# ``ignore_end = False`` is precisely the one whose SCNB *sander and pmemd*
# themselves use for that pair.  Deriving the pair list from those survivors is
# therefore not "first one wins" -- it is "the same one AMBER wins with", and it
# makes the GROMACS pair set identical, row for row, to what ParmEd's own writer
# would have emitted.
#
# That leaves one residual question worth asking rather than assuming: do the
# *discarded* dihedrals ever disagree with the survivor about SCNB?  If they did,
# AMBER's own answer would depend on prmtop ordering and the value would deserve
# a warning even though PRISM had faithfully reproduced it.  Measured on a real
# 57 271-atom Lipid21 POPC bilayer (183 POPC + protein + ions + water, 109 225
# dihedrals):
#
#     71 920 1-4 pairs emitted   (68 077 at SCNB=2.0, 3 843 at SCNB=6.0)
#     12 942 of those 71 920 pairs (18 %) are reached by 2-5 distinct dihedrals
#          0 of those 12 942 disagree about SCNB
#          0 disagree about SCEE
#
# So the conflict does not occur -- not vacuously, but across ~13 000 genuine
# opportunities.  This is structural rather than lucky: the SCNB=6.0 torsions are
# the saturated acyl tails, which are unbranched chains, so their 1-4 pairs have
# no second path; multiple coverage comes from rings and Fourier series, and
# within either of those all the terms come from the same force field with the
# same SCNB.  Because it is measured rather than guaranteed, the policy is
# **detect and refuse, never average and never pick**: :func:`_scan_14_pairs`
# checks every covering dihedral, not just the surviving one, and raises
# :class:`Amber14ConflictError` if any pair is ever ambiguous.  A system that
# trips it is not silently mis-parametrised, it fails to build.
#
# Nothing here hardcodes 0.5, or 6.0, or 1.2.  Every number is read out of the
# parameters that were actually loaded, written into the file, and then read back
# out of the file that was actually written.

#: How close two SCNB (or SCEE) values must be to count as the same value.
#: Deliberately the same ``'%.5f'`` bucketing ParmEd itself uses to decide
#: uniqueness, so PRISM and ParmEd can never disagree about whether a structure
#: is uniform.
_FUDGE_VALUE_FORMAT = "%.5f"

#: Tolerance for "the number in the file is the number we asked for".  The
#: comparison is against the *rendered* text (see :func:`_format_fudge`), so this
#: only has to absorb a re-parse; it stays orders of magnitude tighter than any
#: difference that could matter physically -- the failure being guarded against
#: is 1.0 versus 0.5, which is 100 % off, not a rounding artefact.
_FUDGE_VERIFY_ATOL = 1e-9

#: How far the rendered decimal is allowed to sit from the exact ``1/SCEE``.
#: ``1/1.2`` is not representable, so ``fudgeQQ`` is always a truncation; 1e-6
#: relative is far below the precision GROMACS' own force field files carry.
_FUDGE_RENDER_RTOL = 1e-6


class AmberFudgeError(SolventParameterError):
    """Raised when AMBER 1-4 scaling cannot be expressed as a global fudge factor.

    A subclass of :class:`SolventParameterError` so existing solvent-path error
    handling still catches it, but distinct so a caller that must not downgrade
    it into a warning can single it out.  :mod:`prism.membrane.builder` in
    particular converts ParmEd failures into warnings and returns ``False``;
    this one is a *wrong-physics* condition, not a tool fault, and must be
    re-raised as a build error instead.
    """


class Amber14ConflictError(AmberFudgeError):
    """Raised when one 1-4 atom pair is claimed by dihedrals with different SCNB/SCEE.

    Explicit per-pair parameters can express any *single* scale factor per pair,
    but not two at once.  When the prmtop offers two, the value AMBER itself uses
    depends on which dihedral it happens to list first, so there is no answer
    PRISM could pick that would be reproducibly right.  Measured frequency on a
    real Lipid21 bilayer: zero, over ~13 000 multiply-covered pairs (see the
    module comment).  Kept as a hard failure so that "zero" stays a fact about
    every system PRISM builds rather than a fact about the one it was measured on.
    """


@dataclass(frozen=True)
class AmberFudge:
    """The 1-4 scaling of an AMBER structure, in GROMACS' vocabulary.

    ``fudge_lj``/``fudge_qq`` are ``None`` when the structure pins down no value
    -- which happens only when it has no 1-4 pairs at all, and is therefore
    inert rather than unknown (:func:`derive_amber_fudge` proves that rather
    than assuming it).
    """

    #: ``1 / SCNB``, or ``None`` when the structure has no 1-4 LJ pairs.
    fudge_lj: Optional[float]
    #: ``1 / SCEE``, or ``None`` when the structure has no 1-4 pairs.
    fudge_qq: Optional[float]
    #: Every distinct SCNB value found, ascending.  More than one entry means
    #: :func:`derive_amber_fudge` was called with ``strict=False``.
    scnb_values: Tuple[float, ...] = ()
    scee_values: Tuple[float, ...] = ()
    #: Number of 1-4 pairs GROMACS will generate -- counted with exactly the
    #: filter ParmEd's writer uses, so it is the count of ``[ pairs ]`` lines.
    pair_count: int = 0
    #: ``SCNB value -> (number of 1-4 pairs, representative atom-type quartet)``.
    #: Makes the non-uniform message actionable, and is what lets a build log say
    #: "3 843 of 71 920 pairs are lipid tails" instead of just "mixed".
    scnb_detail: Dict[float, Tuple[int, str]] = field(default_factory=dict)
    #: The ``gen-pairs`` setting this scaling requires.  ``'no'`` means the 1-4
    #: LJ parameters are written per pair and ``fudge_lj`` is inert (see below);
    #: ``'yes'`` means GROMACS generates them and ``fudge_lj`` is load-bearing.
    gen_pairs: str = "yes"
    #: Number of ``[ pairs ]`` rows carrying explicit sigma/epsilon columns.
    #: Zero unless :func:`attach_explicit_14_pairs` ran.
    explicit_pairs: int = 0

    @property
    def uniform(self) -> bool:
        """True when a single global ``fudgeLJ``/``fudgeQQ`` is faithful."""
        return len(self.scnb_values) <= 1 and len(self.scee_values) <= 1

    @property
    def explicit(self) -> bool:
        """True when 1-4 LJ is carried per pair rather than by ``fudgeLJ``."""
        return self.gen_pairs == "no" and self.explicit_pairs > 0

    def describe(self) -> str:
        """One-line human summary for build logs."""
        scnb = ",".join(f"{v:g}" for v in self.scnb_values) or "n/a"
        scee = ",".join(f"{v:g}" for v in self.scee_values) or "n/a"
        if self.explicit:
            # The per-SCNB split is the whole reason this path exists, so it is
            # the thing the log line leads with.
            split = ", ".join(
                f"{count} at SCNB={value:g}"
                for value, (count, _quartet) in sorted(self.scnb_detail.items())
            )
            return (
                f"explicit 1-4 LJ on {self.explicit_pairs} pairs (gen-pairs=no; {split}) "
                f"fudgeQQ={_format_fudge(self.fudge_qq)} (SCEE={scee}); "
                f"fudgeLJ={_format_fudge(self.fudge_lj)} written but unused"
            )
        return (
            f"fudgeLJ={_format_fudge(self.fudge_lj)} (SCNB={scnb}) "
            f"fudgeQQ={_format_fudge(self.fudge_qq)} (SCEE={scee}) "
            f"over {self.pair_count} 1-4 pairs"
        )


@dataclass
class _Pairs14Scan:
    """The 1-4 pair set of an AMBER structure, and the SCNB/SCEE behind each pair.

    ``records`` is ``(atom1, atom2, scnb, scee)`` for exactly the pairs GROMACS
    will apply, in the order ParmEd's writer emits them.  Everything that needs
    to reason about 1-4 scaling goes through this one walk, so the count used to
    *report* the scaling and the count used to *write* it cannot drift apart.
    """

    records: List[Tuple[object, object, float, float]] = field(default_factory=list)
    scnb_stats: Dict[str, List] = field(default_factory=dict)
    scee_stats: Dict[str, List] = field(default_factory=dict)
    #: 1-4 pairs whose dihedral type carries no usable SCNB at all.
    unscaled_pairs: int = 0
    #: Number of *distinct* 1-4 pairs reached by more than one dihedral -- the
    #: population in which an SCNB conflict could have shown up.  Reported so
    #: that "no conflicts" is visibly a measurement and not an empty set.
    multiply_covered: int = 0

    @property
    def pair_count(self) -> int:
        return len(self.records)


def _quartet(dihedral) -> str:
    """The four atom types of *dihedral*, for error messages and logs.

    Knowing the offending torsion is ``cD-cD-cD-cD`` is what turns "mixed SCNB"
    into "the lipid tails".  Built lazily -- only for the first dihedral in each
    distinct SCNB bucket -- because doing it for all ~10^5 dihedrals of a bilayer
    costs more than the rest of the walk put together.
    """
    return "-".join(
        str(getattr(atom, "type", "?"))
        for atom in (dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4)
    )


def _bump(stats: Dict[str, List], value: float, dihedral) -> None:
    bucket = _FUDGE_VALUE_FORMAT % value
    row = stats.get(bucket)
    if row is None:
        stats[bucket] = [value, 1, _quartet(dihedral)]
    else:
        row[1] += 1


def _dihedral_scale_factors(dihedral) -> Tuple[Optional[float], Optional[float]]:
    """Return the single (SCNB, SCEE) *this* dihedral specifies, or ``None`` each.

    A ``DihedralTypeList`` holds one entry per periodicity of the same torsion.
    AMBER stores one scale factor per dihedral, so those entries agree in every
    real force field -- but they are separate objects and nothing enforces it, so
    disagreement is detected here rather than resolved by taking entry zero.
    """
    dtype = dihedral.type
    terms = dtype if isinstance(dtype, _DihedralTypeList) else [dtype]
    scnb: Optional[float] = None
    scee: Optional[float] = None
    for term in terms:
        # ``if value`` (not ``is not None``) mirrors ParmEd: AMBER writes a 0.0
        # scale factor for dihedral types whose 1-4 term is switched off, and
        # 1/0 is not a scale factor.
        raw_scnb = getattr(term, "scnb", None)
        if raw_scnb:
            value = float(raw_scnb)
            if scnb is not None and not _same_scale(scnb, value):
                raise Amber14ConflictError(
                    f"The dihedral {_quartet(dihedral)} carries more than one SCNB across its "
                    f"periodicity terms ({scnb:g} and {value:g}). A single 1-4 pair cannot be "
                    "given two Lennard-Jones scale factors, and choosing one of them would make "
                    "the result depend on the order the terms happen to be stored in."
                )
            scnb = value
        raw_scee = getattr(term, "scee", None)
        if raw_scee:
            value = float(raw_scee)
            if scee is not None and not _same_scale(scee, value):
                raise Amber14ConflictError(
                    f"The dihedral {_quartet(dihedral)} carries more than one SCEE across its "
                    f"periodicity terms ({scee:g} and {value:g})."
                )
            scee = value
    return scnb, scee


def _same_scale(a: float, b: float) -> bool:
    """Whether two scale factors are the same value, at ParmEd's own resolution."""
    return (_FUDGE_VALUE_FORMAT % a) == (_FUDGE_VALUE_FORMAT % b)


def _fully_excluded(atom, cache: Dict[int, frozenset]) -> frozenset:
    """Identities of the atoms 1-2 or 1-3 bonded to *atom*, memoised.

    ``Atom.bond_partners`` and ``Atom.angle_partners`` are ParmEd *properties*
    that rebuild a set and then ``sorted()`` it on every single access.  The 1-4
    walk asks about both for every dihedral -- ~2 x 10^5 times on a bilayer, each
    one dragging in a sort whose comparison key is itself the ``Atom.idx``
    property -- which profiled as the single largest cost in this module by a
    wide margin.  Only membership is needed, never order, so each atom's
    exclusions are built once and reused.

    The private ``_bond_partners``/``_angle_partners`` lists are preferred purely
    to skip that sort; the public properties are used when they are absent, so a
    ParmEd that renames them loses the speed and keeps the answer.
    """
    key = id(atom)
    found = cache.get(key)
    if found is None:
        bonded = getattr(atom, "_bond_partners", None)
        angled = getattr(atom, "_angle_partners", None)
        if bonded is None or angled is None:  # pragma: no cover - ParmEd fallback
            bonded, angled = atom.bond_partners, atom.angle_partners
        found = frozenset([id(a) for a in bonded] + [id(a) for a in angled])
        cache[key] = found
    return found


def _scan_14_pairs(parm) -> _Pairs14Scan:
    """Walk *parm*'s dihedrals and return the 1-4 pairs GROMACS will apply.

    The filter is ParmEd's writer's, exactly: skip ``improper`` and typeless
    dihedrals, skip those whose 1-4 partners are already bonded or angle
    partners, and honour ``ignore_end`` -- which
    ``Structure.update_dihedral_exclusions()`` is called first to bring up to
    date, so the flags filtered on here are the flags the writer would filter on.

    Every dihedral that *covers* a pair is examined, not just the one that owns
    it, so that a disagreement between the owner and a discarded duplicate is
    caught instead of being decided by prmtop ordering.  See the module comment
    for why that arbitration belongs to AMBER rather than to PRISM, and for the
    measured frequency of the conflict.
    """
    dihedrals = getattr(parm, "dihedrals", None)
    if dihedrals is None:
        raise AmberFudgeError(
            f"Cannot derive 1-4 scaling from {type(parm).__name__}: it has no dihedral list."
        )

    # ParmEd's writer calls this immediately before emitting [ pairs ].
    updater = getattr(parm, "update_dihedral_exclusions", None)
    if callable(updater):
        try:
            updater()
        except Exception:  # pragma: no cover - defensive
            pass

    scan = _Pairs14Scan()
    # key -> (scnb, scee, dihedral) of the first dihedral seen covering that pair,
    # whether or not that dihedral is the one that owns the 1-4 term.  The
    # dihedral is kept rather than its rendered atom-type quartet because there
    # is one entry per 1-4 pair -- tens of thousands per bilayer -- and the
    # quartet is only ever read on the error path.
    coverage: Dict[Tuple[int, int], Tuple[Optional[float], Optional[float], object]] = {}
    owned: set = set()
    multiply_covered: set = set()

    excluded_cache: Dict[int, frozenset] = {}
    for dihedral in dihedrals:
        if dihedral.type is None or getattr(dihedral, "improper", False):
            continue
        first, last = dihedral.atom1, dihedral.atom4
        # The writer drops 1-4 "pairs" that are also 1-2 or 1-3 pairs: those are
        # fully excluded, so they are not 1-4 interactions at all.
        if id(first) in _fully_excluded(last, excluded_cache):
            continue

        scnb, scee = _dihedral_scale_factors(dihedral)
        first_idx, last_idx = first.idx, last.idx
        key = (first_idx, last_idx) if first_idx < last_idx else (last_idx, first_idx)
        seen = coverage.get(key)
        if seen is None:
            coverage[key] = (scnb, scee, dihedral)
        else:
            multiply_covered.add(key)
            _assert_no_conflict(key, seen, (scnb, scee), dihedral, parm)

        if getattr(dihedral, "ignore_end", False):
            continue
        if key in owned:
            # update_dihedral_exclusions() is supposed to make this unreachable;
            # if it ever fires, every affected 1-4 interaction would be counted
            # twice, so it is a hard error rather than a de-duplication.
            raise Amber14ConflictError(
                f"Two dihedrals both claim the 1-4 pair (atoms {key[0] + 1}, {key[1] + 1}) with "
                "ignore_end unset, so its Lennard-Jones and Coulomb interactions would be "
                "counted twice. The prmtop's end-group exclusion flags are inconsistent."
            )
        owned.add(key)

        if scnb:
            _bump(scan.scnb_stats, scnb, dihedral)
        else:
            scan.unscaled_pairs += 1
        if scee:
            _bump(scan.scee_stats, scee, dihedral)
        scan.records.append((first, last, scnb or 0.0, scee or 0.0))

    scan.multiply_covered = len(multiply_covered)
    return scan


def _assert_no_conflict(key, seen, incoming, dihedral, parm) -> None:
    """Refuse a 1-4 pair whose covering dihedrals disagree about SCNB or SCEE."""
    for label, old, new in (("SCNB", seen[0], incoming[0]), ("SCEE", seen[1], incoming[1])):
        if old is None or new is None or _same_scale(old, new):
            continue
        try:
            residues = "/".join(
                sorted({parm.atoms[key[0]].residue.name, parm.atoms[key[1]].residue.name})
            )
        except Exception:  # pragma: no cover - defensive
            residues = "?"
        raise Amber14ConflictError(
            f"The 1-4 pair (atoms {key[0] + 1}, {key[1] + 1}, residue(s) {residues}) is reached by "
            f"two dihedrals that disagree about {label}: {_quartet(seen[2])} specifies {old:g} while "
            f"{_quartet(dihedral)} specifies {new:g}. AMBER resolves this by whichever torsion the "
            "prmtop happens to list first, so the 1-4 energy of this pair is not reproducibly "
            "defined and PRISM will not guess at it. This has never been observed in an AMBER "
            "force field (measured: 0 conflicts over ~13 000 multiply-covered pairs in a Lipid21 "
            "bilayer); if it reaches you, the prmtop or its parameter set is the thing to inspect."
        )


def attach_explicit_14_pairs(parm, *, context: str = "") -> AmberFudge:
    """Give *parm* explicit per-pair 1-4 LJ parameters, and describe what it got.

    Builds one :class:`~parmed.topologyobjects.NonbondedException` per 1-4 pair,
    carrying the *combined* Lennard-Jones parameters with that pair's own SCNB
    already divided into epsilon, and assigns them to ``parm.adjusts``.  ParmEd's
    GROMACS writer emits ``adjusts`` as ``[ pairs ]`` rows with explicit sigma and
    epsilon columns, which is what lets Lipid21's SCNB=2.0 and SCNB=6.0 torsions
    coexist in one topology instead of being averaged into a wrong global
    ``fudgeLJ``.

    The structure is *mutated*.  :func:`save_gromacs_topology` restores the
    previous ``adjusts`` afterwards so that writing a topology is not also a
    silent edit of the caller's object; a caller using this directly owns that.

    Returns the :class:`AmberFudge` describing the result, with ``gen_pairs ==
    'no'`` and ``explicit_pairs`` set.  Raises :class:`Amber14ConflictError` if
    any pair's scaling is ambiguous and :class:`AmberFudgeError` if the Coulomb
    scaling cannot be represented (see below).
    """
    scan = _scan_14_pairs(parm)
    scnb_values = _sorted_values(scan.scnb_stats)
    scee_values = _sorted_values(scan.scee_stats)
    scnb_detail = {value: (count, quartet) for value, count, quartet in _stat_rows(scan.scnb_stats)}

    if not scan.records:
        # No dihedrals means no 1-4 pairs, so there is nothing to scale and no
        # question to answer.  The solvent probe (water + two ions) lands here.
        return AmberFudge(
            fudge_lj=None,
            fudge_qq=None,
            scnb_values=scnb_values,
            scee_values=scee_values,
            pair_count=0,
            scnb_detail=scnb_detail,
            gen_pairs="yes",
            explicit_pairs=0,
        )

    # A funct-1 [ pairs ] row has columns for sigma and epsilon and nothing else:
    # there is no per-pair charge-scale column, so unlike SCNB the Coulomb factor
    # really must be global.  Every AMBER force field PRISM uses sets SCEE=1.2
    # throughout, so this is a guard rather than a limitation.
    if len(scee_values) > 1:
        raise AmberFudgeError(_nonuniform_message("SCEE", "fudgeQQ", scan.scee_stats, scan.pair_count))
    if scan.unscaled_pairs:
        raise AmberFudgeError(
            f"{scan.unscaled_pairs} of {scan.pair_count} 1-4 pairs come from dihedral types that "
            "specify no SCNB scale factor, so the Lennard-Jones strength AMBER intends for them is "
            f"undefined{_suffix(context)}."
        )
    if not scee_values:
        raise AmberFudgeError(
            f"This structure generates {scan.pair_count} 1-4 pairs but pins down no SCEE scale "
            f"factor, so the 1-4 Coulomb scaling GROMACS should apply is unknown{_suffix(context)}."
        )

    # Imported here rather than at the top of the function so that a structure
    # with no dihedrals is rejected on its own merits: the "has no dihedral list"
    # error is about the structure, and it must not be pre-empted by an
    # ImportError from a ParmEd that was never going to be needed.
    try:
        from parmed import TrackedList
        from parmed.topologyobjects import NonbondedException, NonbondedExceptionType
    except Exception as exc:  # pragma: no cover - environment dependent
        raise SolventParameterError(
            f"ParmEd's nonbonded-exception classes are unavailable, so explicit 1-4 parameters "
            f"cannot be built{_suffix(context)}: {exc}"
        ) from exc

    scee = scee_values[0]
    chgscale = 1.0 / scee
    # OPLS-style structures combine sigma geometrically; AMBER uses
    # Lorentz-Berthelot.  ParmEd stores ``atom.rmin`` as Rmin/2, so the
    # Lorentz-Berthelot combination is a plain sum.  The ``_14`` variants are
    # used because CHAMBER-derived structures carry a separate 1-4 Lennard-Jones
    # set; for a plain AMBER prmtop they are the same numbers.
    geometric = str(getattr(parm, "combining_rule", "lorentz")).lower() == "geometric"

    adjusts = TrackedList()
    adjust_types = TrackedList()
    interned: Dict[Tuple[int, int], NonbondedExceptionType] = {}

    for atom1, atom2, scnb, _scee in scan.records:
        r1 = float(getattr(atom1, "rmin_14", None) or atom1.rmin)
        r2 = float(getattr(atom2, "rmin_14", None) or atom2.rmin)
        e1 = float(getattr(atom1, "epsilon_14", None) or atom1.epsilon)
        e2 = float(getattr(atom2, "epsilon_14", None) or atom2.epsilon)
        rmin = 2.0 * math.sqrt(r1 * r2) if geometric else r1 + r2
        epsilon = math.sqrt(e1 * e2) / scnb
        # Distinct atom-type combinations number in the hundreds while pairs
        # number in the tens of thousands, so the types are interned: it keeps
        # ParmEd's parameter-set construction and its [ pairtypes ] section
        # proportional to the force field rather than to the system.  The key is
        # the rendered value, so two types are shared only when they would have
        # written identical columns.
        key = ("%.9f" % rmin, "%.9f" % epsilon)
        pair_type = interned.get(key)
        if pair_type is None:
            pair_type = NonbondedExceptionType(rmin, epsilon, chgscale, list=adjust_types)
            adjust_types.append(pair_type)
            interned[key] = pair_type
        adjusts.append(NonbondedException(atom1, atom2, pair_type))

    parm.adjusts = adjusts
    parm.adjust_types = adjust_types

    return AmberFudge(
        # Inert but not arbitrary: with ``gen-pairs = no`` GROMACS never reads
        # fudgeLJ, so this records the value the *majority* of pairs imply.  It is
        # the least-wrong number to leave behind if the explicit columns are ever
        # stripped, and it is never 1.0 -- the "do not scale at all" value this
        # module exists to keep out of PRISM's topologies.
        fudge_lj=1.0 / max(_stat_rows(scan.scnb_stats), key=lambda row: row[1])[0],
        fudge_qq=chgscale,
        scnb_values=scnb_values,
        scee_values=scee_values,
        pair_count=scan.pair_count,
        scnb_detail=scnb_detail,
        gen_pairs="no",
        explicit_pairs=len(adjusts),
    )


def _suffix(context: str) -> str:
    return f" ({context})" if context else ""


def derive_amber_fudge(parm, *, strict: bool = True) -> AmberFudge:
    """Read GROMACS ``fudgeLJ``/``fudgeQQ`` off an AMBER structure's own SCNB/SCEE.

    *parm* is a ParmEd structure (typically an ``AmberParm`` from
    ``pmd.load_file(prmtop)``).  The SCNB/SCEE values are collected over exactly
    the dihedrals ParmEd's GROMACS writer turns into ``[ pairs ]`` lines --
    skipping ``ignore_end`` and ``improper`` dihedrals and 1-4 "pairs" that are
    really bonded or angle partners -- so the derived factors apply to precisely
    the pairs that end up in the file, no more and no fewer.

    With ``strict=True`` (the default) a non-uniform SCNB or SCEE raises
    :class:`AmberFudgeError`, because no single global factor can represent it
    and silently choosing one member of the set is how 1-4 interactions end up
    several times too strong.  ``strict=False`` returns the full value sets
    instead of raising, for callers that want to *report* the situation.

    This answers the narrower question "is a single global ``fudgeLJ`` faithful
    to this structure?", which is still worth asking -- of a topology written
    with ``gen-pairs = yes``, or of a force field one is about to trust.  It is
    no longer how PRISM *writes* topologies: :func:`save_gromacs_topology` emits
    explicit per-pair parameters, which represent mixed SCNB exactly and so have
    no uniformity requirement to enforce.
    """
    scan = _scan_14_pairs(parm)
    scnb_stats = scan.scnb_stats
    scee_stats = scan.scee_stats
    pair_count = scan.pair_count
    unscaled_pairs = scan.unscaled_pairs
    scnb_values = _sorted_values(scnb_stats)
    scee_values = _sorted_values(scee_stats)

    if strict:
        if len(scnb_values) > 1:
            raise AmberFudgeError(_nonuniform_message("SCNB", "fudgeLJ", scnb_stats, pair_count))
        if len(scee_values) > 1:
            raise AmberFudgeError(_nonuniform_message("SCEE", "fudgeQQ", scee_stats, pair_count))
        # 1-4 pairs whose dihedral type carries no SCNB at all would be written
        # into [ pairs ] and then scaled by whatever the other dihedrals voted
        # for -- a different wrong answer from the one this module exists to fix.
        if unscaled_pairs and scnb_values:
            raise AmberFudgeError(
                f"{unscaled_pairs} of {pair_count} 1-4 pairs come from dihedral types with no "
                f"SCNB scale factor, but the rest specify SCNB={scnb_values[0]:g}. GROMACS would "
                "scale all of them alike. Explicit per-pair parameters in [ pairs ] are needed to "
                "represent this system."
            )

    fudge_lj = 1.0 / scnb_values[0] if len(scnb_values) == 1 else None
    fudge_qq = 1.0 / scee_values[0] if len(scee_values) == 1 else None

    if strict and pair_count and fudge_lj is None:
        raise AmberFudgeError(
            f"This structure generates {pair_count} 1-4 pairs but pins down no SCNB scale factor, "
            "so the 1-4 Lennard-Jones scaling GROMACS should apply cannot be determined."
        )

    return AmberFudge(
        fudge_lj=fudge_lj,
        fudge_qq=fudge_qq,
        scnb_values=scnb_values,
        scee_values=scee_values,
        pair_count=pair_count,
        scnb_detail={
            value: (count, quartet) for value, count, quartet in _stat_rows(scnb_stats)
        },
    )


def _sorted_values(stats: Dict[str, List]) -> Tuple[float, ...]:
    return tuple(sorted(row[0] for row in stats.values()))


def _stat_rows(stats: Dict[str, List]) -> List[Tuple[float, int, str]]:
    return sorted((row[0], row[1], row[2]) for row in stats.values())


def _nonuniform_message(
    amber_name: str, gromacs_name: str, stats: Dict[str, List], pair_count: int
) -> str:
    """Explain a non-uniform 1-4 scale factor and what representing it requires."""
    rows = _stat_rows(stats)
    breakdown = "; ".join(
        f"{amber_name}={value:g} on {count} 1-4 term(s), e.g. {quartet}"
        for value, count, quartet in rows
    )
    biggest = max(rows, key=lambda row: row[1])
    common = (
        f"This structure's AMBER {amber_name} is not uniform, so it cannot be written as a single "
        f"GROMACS {gromacs_name}: {len(rows)} distinct {amber_name} values are in use "
        f"({breakdown}) across {pair_count} generated 1-4 pairs. GROMACS applies one global "
        f"{gromacs_name} to every pair, so picking {amber_name}={biggest[0]:g} (the most common) "
        f"would rescale every other pair by the wrong factor -- silently, because any value in "
        f"[ defaults ] is legal."
    )
    if amber_name == "SCNB":
        # Reachable from derive_amber_fudge(strict=True) only.  The writer no
        # longer has this problem, so the message says where the answer is rather
        # than presenting it as a dead end.
        return common + (
            " Use save_gromacs_topology(), which writes explicit per-pair parameters "
            "('gen-pairs = no' plus sigma/epsilon columns on every [ pairs ] row) and represents "
            "mixed SCNB exactly. Note that AMBER Lipid21 is non-uniform by design -- its "
            "cD-cD-cD-cD tail torsion is fit with SCNB=6.0 while the rest of the lipid uses "
            "SCNB=2.0 -- so this is expected of any Lipid17/Lipid21 bilayer and is not a corrupt "
            "prmtop."
        )
    return common + (
        " Unlike SCNB this cannot be moved onto the individual [ pairs ] rows: a funct-1 pair "
        "carries sigma and epsilon and no charge-scale column, so the 1-4 Coulomb factor is "
        "necessarily global. Every AMBER force field in normal use sets SCEE=1.2 throughout, so a "
        "structure that reaches this point has had incompatible parameter sets combined."
    )


def save_gromacs_topology(
    parm,
    out_top: str,
    *,
    overwrite: bool = True,
    context: str = "",
) -> AmberFudge:
    """Convert *parm* to a GROMACS topology at *out_top* with correct 1-4 scaling.

    Drop-in replacement for ``parm.save(out_top, format="gromacs",
    overwrite=True)``, with the same signature it has always had.  What is added
    is that each 1-4 Lennard-Jones interaction is written with **its own** scale
    factor applied, instead of every one of them being scaled by a single global
    ``fudgeLJ`` that cannot be right for a Lipid21 bilayer whatever value it
    takes.  Concretely the topology comes out with ``gen-pairs = no`` and sigma /
    epsilon columns on every ``[ pairs ]`` row.

    All of it is **verified by reading the written file back**: that the
    ``[ defaults ]`` line says what it was told to say, that no ``[ pairs ]`` row
    was left without parameters, and that the number of 1-4 pairs the topology
    actually contains -- expanded through ``[ molecules ]``, because GROMACS
    stores each ``[ moleculetype ]`` once -- equals the number derived from
    *parm*.  A topology can be wrong while every API call succeeded, so success
    is defined as what is in the file.

    Returns the :class:`AmberFudge` describing what was applied, so the caller
    can log it; ``AmberFudge.describe()`` is a ready-made one-liner.  Raises
    :class:`Amber14ConflictError` if any pair's scaling is ambiguous, and
    :class:`AmberFudgeError` if the 1-4 *Coulomb* scaling is non-uniform (which,
    unlike SCNB, GROMACS genuinely cannot represent per pair).

    *context* is a short phrase naming the topology ("the dry membrane
    topology") used to make errors locatable.
    """
    if os.path.exists(out_top) and not overwrite:
        raise SolventParameterError(
            f"{_where(out_top, context)} exists and overwrite=False."
        )

    previous_adjusts = getattr(parm, "adjusts", None)
    previous_adjust_types = getattr(parm, "adjust_types", None)
    if previous_adjusts:
        # The structure already carries explicit nonbonded exceptions (a CHARMM
        # or already-converted structure).  Those are authoritative -- rebuilding
        # them from dihedrals would discard whatever the source encoded -- so
        # they are written as they stand and only surveyed for the log.
        fudge = derive_amber_fudge(parm, strict=False)
        fudge = AmberFudge(
            fudge_lj=fudge.fudge_lj,
            fudge_qq=fudge.fudge_qq,
            scnb_values=fudge.scnb_values,
            scee_values=fudge.scee_values,
            pair_count=fudge.pair_count,
            scnb_detail=fudge.scnb_detail,
            gen_pairs="no",
            explicit_pairs=len(previous_adjusts),
        )
    else:
        fudge = attach_explicit_14_pairs(parm, context=context)

    try:
        _write_gromacs_topology(parm, out_top, fudge, context=context)
    finally:
        # Writing a topology should not also be a silent edit of the caller's
        # structure: a caller that goes on to ``parm.save(gro)`` or re-saves the
        # prmtop gets back exactly the object it handed over.
        if previous_adjusts is not None:
            parm.adjusts = previous_adjusts
        if previous_adjust_types is not None:
            parm.adjust_types = previous_adjust_types

    if not os.path.isfile(out_top):
        raise SolventParameterError(
            f"ParmEd returned without writing {_where(out_top, context)}."
        )
    apply_amber_fudge_to_topology(out_top, fudge, context=context)
    return fudge


def _write_gromacs_topology(parm, out_top: str, fudge: AmberFudge, *, context: str) -> None:
    """Run ParmEd's GROMACS writer with ``[ defaults ]`` pinned to *fudge*.

    Not ``parm.save(format="gromacs")``: that routes through
    ``GromacsTopologyFile.from_structure``, whose ``[ defaults ]`` handling is
    the thing being corrected.  ``from_structure`` sets ``gen-pairs = no`` by
    itself when the structure has ``adjusts`` -- but only on the branch where the
    structure has no ``defaults`` attribute of its own, and a caller may well
    have set one.  Pinning the three fields afterwards makes the outcome
    independent of which branch ran.
    """
    try:
        from parmed.gromacs import GromacsTopologyFile
    except Exception as exc:  # pragma: no cover - environment dependent
        raise SolventParameterError(
            f"ParmEd's GROMACS writer is unavailable, so {_where(out_top, context)} cannot be "
            f"written: {exc}"
        ) from exc

    try:
        gmx = GromacsTopologyFile.from_structure(parm)
    except Exception as exc:
        raise SolventParameterError(
            f"ParmEd could not convert {_where(out_top, context)} to GROMACS format: {exc}"
        ) from exc

    gmx.defaults.gen_pairs = fudge.gen_pairs
    if fudge.fudge_lj is not None:
        gmx.defaults.fudgeLJ = float(fudge.fudge_lj)
    if fudge.fudge_qq is not None:
        gmx.defaults.fudgeQQ = float(fudge.fudge_qq)
    # 'inline' is ParmEd's default and is what puts the sigma/epsilon columns on
    # the [ pairs ] rows themselves; with parameters written to a separate file,
    # rows whose parameters match the shared [ pairtypes ] entry are emitted bare
    # and would silently pick up the wrong pair's epsilon under mixed SCNB.
    gmx.write(out_top, parameters="inline")


def apply_amber_fudge_to_topology(
    top_path: str, fudge: AmberFudge, *, context: str = ""
) -> bool:
    """Force *top_path*'s ``[ defaults ]`` to carry *fudge*, then verify it landed.

    Returns True when the line was changed.  Separate from
    :func:`save_gromacs_topology` so an already-written topology (one converted
    by an earlier PRISM release, or by ``parmed -i``) can be repaired and
    re-verified without redoing the conversion.

    ``nbfunc`` and ``comb-rule`` are preserved exactly: this function's business
    is the 1-4 scaling, and silently changing the combination rule would rewrite
    every Lennard-Jones interaction in the file.  ``gen-pairs`` is set from
    *fudge* when *fudge* pins it down (explicit per-pair parameters require
    ``no``) and preserved otherwise.
    """
    top_path = os.path.abspath(top_path)
    if not os.path.isfile(top_path):
        raise SolventParameterError(f"Topology not found: {top_path}")
    with open(top_path) as handle:
        lines = handle.readlines()

    index, fields = _locate_defaults(lines)
    if index is None:
        raise SolventParameterError(
            f"{_where(top_path, context)} has no [ defaults ] section, so its 1-4 scaling cannot "
            "be set."
        )

    bare_pairs, parameterised_pairs = _count_pair_lines(lines)
    if fudge.fudge_lj is None:
        # Undetermined *and* inert is fine -- undetermined while GROMACS is about
        # to generate 1-4 LJ interactions is not.
        if bare_pairs:
            raise AmberFudgeError(
                f"{_where(top_path, context)} contains {bare_pairs} [ pairs ] entries with no "
                "explicit parameters, so GROMACS will generate their 1-4 Lennard-Jones "
                "interactions and scale them by fudgeLJ -- but no SCNB scale factor could be "
                "derived from the AMBER parameters, so the correct fudgeLJ is unknown."
            )
        return False

    nbfunc = fields[0] if len(fields) > 0 else "1"
    comb_rule = fields[1] if len(fields) > 1 else "2"
    gen_pairs = fudge.gen_pairs or (fields[2] if len(fields) > 2 else "yes")
    fudge_qq = fudge.fudge_qq
    if fudge_qq is None:
        # Keep whatever ParmEd derived rather than inventing one; fudgeQQ is only
        # None when there is nothing for it to scale.
        fudge_qq = _as_float(fields[4]) if len(fields) > 4 else None

    replacement = "%-15s %-15s %-15s %-12s %s\n" % (
        nbfunc,
        comb_rule,
        gen_pairs,
        _format_fudge(fudge.fudge_lj),
        _format_fudge(fudge_qq),
    )
    changed = lines[index] != replacement
    if changed:
        lines[index] = replacement
        with open(top_path, "w") as handle:
            handle.writelines(lines)

    # Read the file back: the point of this module is that a topology can be
    # wrong while every API call succeeded, so success is defined as "the number
    # is in the file", not "the write did not raise".
    written = topology_defaults(top_path)
    if written is None:
        raise SolventParameterError(
            f"{_where(top_path, context)} lost its [ defaults ] section while its 1-4 scaling "
            "was being set."
        )
    for label, want, got in (
        ("fudgeLJ", fudge.fudge_lj, written.get("fudgeLJ")),
        ("fudgeQQ", fudge_qq, written.get("fudgeQQ")),
    ):
        if want is None:
            continue
        # Compare against the decimal that was rendered, not against the exact
        # binary quotient: 1/1.2 has no finite decimal form, so demanding the
        # file reproduce it bit for bit would fail on a correct write.  The
        # rendering itself is checked for fidelity separately.
        rendered = float(_format_fudge(want))
        if abs(rendered - float(want)) > _FUDGE_RENDER_RTOL * max(1.0, abs(float(want))):
            raise SolventParameterError(
                f"Rendering {label}={want!r} as text lost too much precision "
                f"({rendered!r}); refusing to write {_where(top_path, context)}."
            )
        if got is None or abs(float(got) - rendered) > _FUDGE_VERIFY_ATOL:
            raise SolventParameterError(
                f"{_where(top_path, context)} still reads {label}="
                f"{'<missing>' if got is None else format(got, 'g')} after it was set to "
                f"{rendered:g}; the topology on disk does not match the AMBER parameters it "
                "came from."
            )
    if written.get("gen_pairs") == "yes" and parameterised_pairs:
        raise SolventParameterError(
            f"{_where(top_path, context)} sets gen-pairs = yes but {parameterised_pairs} "
            "[ pairs ] lines carry explicit parameters; GROMACS would use the explicit values for "
            "those and the generated, fudgeLJ-scaled values for the rest."
        )
    if written.get("gen_pairs") == "no" and bare_pairs:
        # grompp would fail on this rather than run it wrong, but failing here
        # names the cause instead of leaving a "No default Pair types" to decode.
        raise SolventParameterError(
            f"{_where(top_path, context)} sets gen-pairs = no but {bare_pairs} [ pairs ] lines "
            "carry no sigma/epsilon columns, so GROMACS has no 1-4 Lennard-Jones parameters for "
            "them at all."
        )

    if fudge.explicit_pairs:
        # The line count is not the pair count: GROMACS writes each
        # [ moleculetype ] once and multiplies it up through [ molecules ], so a
        # 183-lipid bilayer's pairs appear a few hundred times in the file and
        # tens of thousands of times in the simulation.  Comparing the *expanded*
        # count against the structure is what proves no molecule silently lost
        # its pairs to ParmEd's de-duplication.
        census = topology_pair_census(top_path)
        if census["expanded_total"] != fudge.explicit_pairs:
            raise SolventParameterError(
                f"{_where(top_path, context)} expands to {census['expanded_total']} 1-4 pairs but "
                f"{fudge.explicit_pairs} were derived from the AMBER parameters. The topology on "
                "disk does not describe the same set of 1-4 interactions as the prmtop it came "
                "from."
            )
        if census["expanded_bare"]:
            raise SolventParameterError(
                f"{_where(top_path, context)} expands to {census['expanded_bare']} 1-4 pairs with "
                "no explicit parameters."
            )
    return changed


def topology_pair_census(top_path: str) -> Dict[str, int]:
    """Count the ``[ pairs ]`` of *top_path*, both as written and as simulated.

    Keys: ``lines_bare`` / ``lines_parameterised`` are literal line counts;
    ``expanded_bare`` / ``expanded_parameterised`` / ``expanded_total`` multiply
    each ``[ moleculetype ]``'s rows by its ``[ molecules ]`` count, which is the
    number of 1-4 interactions GROMACS will actually evaluate.

    A ``[ moleculetype ]`` that ``[ molecules ]`` never mentions contributes to
    the line counts and not to the expanded ones -- which is correct, and is why
    the two are reported separately.
    """
    with open(top_path) as handle:
        lines = handle.readlines()

    per_molecule: Dict[str, List[int]] = {}
    counts: List[Tuple[str, int]] = []
    directive = ""
    expect_name = False
    current: Optional[str] = None

    for raw in lines:
        match = _DIRECTIVE_RE.match(raw)
        if match:
            directive = match.group(1).lower()
            expect_name = directive == "moleculetype"
            continue
        if raw.lstrip().startswith("#"):
            continue
        stripped = raw.split(";", 1)[0].strip()
        if not stripped:
            continue
        fields = stripped.split()
        if directive == "moleculetype" and expect_name:
            current = fields[0]
            per_molecule.setdefault(current, [0, 0])
            expect_name = False
        elif directive == "pairs" and current is not None and len(fields) >= 2:
            per_molecule[current][1 if len(fields) >= 5 else 0] += 1
        elif directive == "molecules" and len(fields) >= 2:
            try:
                counts.append((fields[0], int(fields[1])))
            except ValueError:
                continue

    census = {
        "lines_bare": sum(row[0] for row in per_molecule.values()),
        "lines_parameterised": sum(row[1] for row in per_molecule.values()),
        "expanded_bare": 0,
        "expanded_parameterised": 0,
    }
    for name, repeats in counts:
        row = per_molecule.get(name)
        if row is None:
            continue
        census["expanded_bare"] += row[0] * repeats
        census["expanded_parameterised"] += row[1] * repeats
    census["expanded_total"] = census["expanded_bare"] + census["expanded_parameterised"]
    return census


def topology_defaults(top_path: str) -> Optional[Dict[str, object]]:
    """Return the parsed ``[ defaults ]`` line of *top_path*, or ``None``.

    Keys: ``nbfunc``, ``comb_rule`` (ints), ``gen_pairs`` (str), ``fudgeLJ``,
    ``fudgeQQ`` (floats).  Missing trailing columns come back absent rather than
    defaulted, so a caller can tell "GROMACS will assume 1.0" from "the file
    says 1.0".
    """
    with open(top_path) as handle:
        lines = handle.readlines()
    index, fields = _locate_defaults(lines)
    if index is None:
        return None
    parsed: Dict[str, object] = {}
    for position, key in ((0, "nbfunc"), (1, "comb_rule")):
        if len(fields) > position:
            try:
                parsed[key] = int(fields[position])
            except ValueError:
                pass
    if len(fields) > 2:
        parsed["gen_pairs"] = fields[2].lower()
    for position, key in ((3, "fudgeLJ"), (4, "fudgeQQ")):
        if len(fields) > position:
            value = _as_float(fields[position])
            if value is not None:
                parsed[key] = value
    parsed["line"] = lines[index].rstrip("\n")
    return parsed


def _locate_defaults(lines: Sequence[str]) -> Tuple[Optional[int], List[str]]:
    """Return the index of the first data line under ``[ defaults ]``, and its fields."""
    in_defaults = False
    for i, raw in enumerate(lines):
        match = _DIRECTIVE_RE.match(raw)
        if match:
            in_defaults = match.group(1).lower() == "defaults"
            continue
        if not in_defaults or raw.lstrip().startswith("#"):
            continue
        stripped = raw.split(";", 1)[0].strip()
        if not stripped:
            continue
        return i, stripped.split()
    return None, []


def _count_pair_lines(lines: Sequence[str]) -> Tuple[int, int]:
    """Return (pairs without explicit parameters, pairs with them).

    ``fudgeLJ`` governs exactly the first group; the second carries its own V/W
    columns and ignores it.
    """
    bare = 0
    parameterised = 0
    directive = ""
    for raw in lines:
        match = _DIRECTIVE_RE.match(raw)
        if match:
            directive = match.group(1).lower()
            continue
        if directive != "pairs" or raw.lstrip().startswith("#"):
            continue
        fields = raw.split(";", 1)[0].split()
        if len(fields) < 2:
            continue
        if len(fields) >= 5:
            parameterised += 1
        else:
            bare += 1
    return bare, parameterised


def _format_fudge(value: Optional[float]) -> str:
    if value is None:
        return "n/a"
    # 1/1.2 needs more than the 2 significant figures '%g' would give it, and
    # GROMACS reads the field as a plain float.
    return f"{value:.8f}".rstrip("0").rstrip(".") or "0"


def _as_float(token: str) -> Optional[float]:
    try:
        return float(token)
    except (TypeError, ValueError):
        return None


def _where(path: str, context: str) -> str:
    name = os.path.basename(path)
    return f"{context} ({name})" if context else name


# --------------------------------------------------------------------------- #
# Probing
# --------------------------------------------------------------------------- #


def solvent_parameter_blocks(
    water_model: str,
    positive_ion: str,
    negative_ion: str,
    work_dir: str,
    amberhome: Optional[str] = None,
    timeout: int = 600,
) -> SolventBlocks:
    """Probe tleap for *water_model* + the two ions and return their blocks.

    Runs a one-shot tleap on a three-molecule unit (one water, one cation, one
    anion), converts the result with ParmEd and lifts out the atomtypes and
    moleculetypes.  Raises :class:`SolventParameterError` on any failure; the
    builder turns that into a fall back to the PACKMOL route.
    """
    water_key = str(water_model or "").strip().lower()
    leaprc = _WATER_LEAPRC.get(water_key)
    if leaprc is None:
        raise SolventParameterError(
            f"No tleap water model known for {water_model!r}; supported: "
            + ", ".join(sorted(_WATER_LEAPRC))
        )

    cation_name = str(positive_ion or "").strip()
    anion_name = str(negative_ion or "").strip()
    cation_unit = _ion_unit(cation_name, "positive")
    anion_unit = _ion_unit(anion_name, "negative")

    probe_dir = os.path.join(os.path.abspath(work_dir), "solvent_probe")
    shutil.rmtree(probe_dir, ignore_errors=True)
    os.makedirs(probe_dir, exist_ok=True)

    prmtop = os.path.join(probe_dir, "solvent_probe.prmtop")
    inpcrd = os.path.join(probe_dir, "solvent_probe.inpcrd")
    script = os.path.join(probe_dir, "solvent_probe.in")
    with open(script, "w") as handle:
        handle.write(
            "# Written by PRISM: one-shot probe for the water/ion parameters the\n"
            "# GROMACS membrane solvation path has to add to a dry topology.\n"
            f"source {leaprc}\n"
            # Order matters only in that it fixes the moleculetype order in the
            # ParmEd output; water last keeps the ion blocks adjacent.
            f"probe = combine {{ {cation_unit} {anion_unit} WAT }}\n"
            f"saveamberparm probe {os.path.basename(prmtop)} {os.path.basename(inpcrd)}\n"
            "quit\n"
        )

    _run_tleap(script, probe_dir, amberhome, timeout)
    if not os.path.isfile(prmtop):
        log = _read_tail(os.path.join(probe_dir, "leap.log"))
        raise SolventParameterError(
            f"tleap did not write {os.path.basename(prmtop)} for water model {water_model!r} "
            f"with ions {cation_unit}/{anion_unit}." + (f"\n{log}" if log else "")
        )

    probe_top = os.path.join(probe_dir, "solvent_probe.top")
    _parmed_convert(prmtop, inpcrd, probe_top)

    atomtypes, moleculetypes, charges, comb_rule = _parse_probe_topology(probe_top)

    # Rename to the names the rest of the fast path resolves against.
    rename = {
        "WAT": WATER_MOLECULETYPE,
        cation_unit: cation_name,
        anion_unit: anion_name,
    }
    moleculetypes, charges = _apply_renames(moleculetypes, charges, rename)

    water_block = dict(moleculetypes).get(WATER_MOLECULETYPE)
    if water_block is None:
        raise SolventParameterError(
            "The tleap probe produced no water moleculetype; cannot solvate."
        )
    water_sites = _count_atoms(water_block)
    _check_site_count(water_key, water_sites)

    return SolventBlocks(
        water_moleculetype=WATER_MOLECULETYPE,
        water_sites=water_sites,
        cation_moleculetype=cation_name,
        anion_moleculetype=anion_name,
        atomtypes=atomtypes,
        moleculetypes=moleculetypes,
        comb_rule=comb_rule,
        water_model=water_key,
        charges=charges,
        source=f"tleap {leaprc} ({cation_unit}/{anion_unit})",
    )


def _ion_unit(name: str, role: str) -> str:
    if not name:
        raise SolventParameterError(f"No {role} ion name was given.")
    unit = _ION_UNITS.get(name.strip().upper())
    if unit is None:
        raise SolventParameterError(
            f"Unknown {role} ion {name!r}; known names: " + ", ".join(sorted(_ION_UNITS))
        )
    return unit


def _run_tleap(script: str, work_dir: str, amberhome: Optional[str], timeout: int) -> None:
    env = dict(os.environ)
    executable = "tleap"
    if amberhome:
        env["AMBERHOME"] = amberhome
        candidate = os.path.join(amberhome, "bin", "tleap")
        if os.path.isfile(candidate):
            executable = candidate
        env["PATH"] = os.path.join(amberhome, "bin") + os.pathsep + env.get("PATH", "")
    if executable == "tleap" and shutil.which("tleap") is None:
        raise SolventParameterError(
            "tleap was not found on PATH, so the water/ion parameters cannot be probed."
        )
    try:
        proc = subprocess.run(
            [executable, "-f", os.path.basename(script)],
            cwd=work_dir,
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise SolventParameterError(f"tleap solvent probe timed out after {timeout} s.") from exc
    except OSError as exc:
        raise SolventParameterError(f"Could not launch tleap for the solvent probe: {exc}") from exc

    # tleap exits 0 even for a fatal parameter error, so the return code alone
    # is not evidence of success; the caller checks that the prmtop exists.
    if proc.returncode != 0:
        raise SolventParameterError(
            f"tleap solvent probe failed (exit {proc.returncode}):\n"
            + _tail(proc.stdout or "")
            + _tail(proc.stderr or "")
        )


def _parmed_convert(prmtop: str, inpcrd: str, out_top: str) -> None:
    try:
        import parmed as pmd
    except Exception as exc:  # pragma: no cover - environment dependent
        raise SolventParameterError(f"ParmEd is required to read the solvent probe: {exc}") from exc
    try:
        parm = pmd.load_file(prmtop, xyz=inpcrd)
    except Exception as exc:
        raise SolventParameterError(f"ParmEd could not read the solvent probe: {exc}") from exc
    try:
        # Not ``parm.save(..., format="gromacs")``: that leaves fudgeLJ at
        # ParmEd's 1.0 fallback whenever SCNB is not a single value, which for
        # the probe means "no dihedrals, so no SCNB at all".  The probe is inert
        # either way -- it has no [ pairs ] for a fudge factor to act on, and
        # ``merge_solvent_blocks`` copies only [ atomtypes ] and the
        # [ moleculetype ] blocks, never [ defaults ] -- but "inert" is checked
        # here instead of being taken on trust, and the same call is what the
        # solute topologies use, so the two conversions cannot drift apart.
        save_gromacs_topology(parm, out_top, overwrite=True, context="the solvent probe topology")
    except AmberFudgeError:
        # Wrong physics, not a tool fault: never softened into a fallback.
        raise
    except SolventParameterError:
        raise
    except Exception as exc:
        raise SolventParameterError(f"ParmEd could not convert the solvent probe: {exc}") from exc
    if not os.path.isfile(out_top):
        raise SolventParameterError("ParmEd returned without writing the solvent probe topology.")


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #


def _parse_probe_topology(
    path: str,
) -> Tuple[Tuple[Tuple[str, str], ...], Tuple[Tuple[str, str], ...], Dict[str, float], int]:
    """Split the probe topology into atomtypes, moleculetype blocks, charges."""
    with open(path) as handle:
        lines = handle.readlines()

    atomtypes: List[Tuple[str, str]] = []
    moleculetypes: List[Tuple[str, str]] = []
    charges: Dict[str, float] = {}
    comb_rule = 2

    directive = ""
    block_lines: List[str] = []
    block_name: Optional[str] = None
    block_charge = 0.0

    def flush() -> None:
        nonlocal block_lines, block_name, block_charge
        if block_name is not None:
            moleculetypes.append((block_name, "".join(block_lines).rstrip() + "\n"))
            charges[block_name] = block_charge
        block_lines = []
        block_name = None
        block_charge = 0.0

    for raw in lines:
        stripped = raw.split(";", 1)[0].strip()
        match = _DIRECTIVE_RE.match(raw)
        if match:
            new_directive = match.group(1).lower()
            if new_directive == "moleculetype":
                flush()
                block_lines = [raw]
                directive = new_directive
                continue
            if new_directive == "system":
                flush()
                directive = new_directive
                continue
            directive = new_directive
            if block_name is not None or block_lines:
                block_lines.append(raw)
            continue

        if block_lines:
            block_lines.append(raw)
            if not stripped:
                continue
            fields = stripped.split()
            if directive == "moleculetype" and block_name is None:
                block_name = fields[0]
            elif directive == "atoms" and len(fields) > 6:
                try:
                    block_charge += float(fields[6])
                except ValueError:
                    pass
            continue

        if not stripped:
            continue
        fields = stripped.split()
        if directive == "atomtypes":
            atomtypes.append((fields[0], raw.rstrip("\n") + "\n"))
        elif directive == "defaults" and len(fields) >= 2:
            try:
                comb_rule = int(fields[1])
            except ValueError:
                pass

    flush()
    if not atomtypes:
        raise SolventParameterError(f"{os.path.basename(path)} has no [ atomtypes ] section.")
    if not moleculetypes:
        raise SolventParameterError(f"{os.path.basename(path)} has no [ moleculetype ] blocks.")
    return tuple(atomtypes), tuple(moleculetypes), charges, comb_rule


def _apply_renames(
    moleculetypes: Sequence[Tuple[str, str]],
    charges: Dict[str, float],
    rename: Dict[str, str],
) -> Tuple[Tuple[Tuple[str, str], ...], Dict[str, float]]:
    """Rename moleculetype/residue/atom names, leaving atom *types* untouched.

    The type column is the link to ``[ atomtypes ]`` and must not move; the
    moleculetype name is what genion resolves, and the residue name is what the
    coordinate file and the purge key on.
    """
    out: List[Tuple[str, str]] = []
    new_charges: Dict[str, float] = {}
    for name, text in moleculetypes:
        target = rename.get(name, name)
        # Ion blocks are monatomic, so their atom name is renamed too; water
        # keeps O/H1/H2 (renaming those would break the element check in
        # solvate._solvent_rename_map).
        rename_atom = target != WATER_MOLECULETYPE
        out.append((target, _rewrite_block(text, name, target, rename_atom)))
        if name in charges:
            new_charges[target] = charges[name]
    return tuple(out), new_charges


def _rewrite_block(text: str, old: str, new: str, rename_atom: bool) -> str:
    """Rewrite one moleculetype block, replacing *old* names with *new*."""
    if old == new:
        return text
    directive = ""
    seen_name = False
    out: List[str] = []
    for raw in text.splitlines(keepends=True):
        match = _DIRECTIVE_RE.match(raw)
        if match:
            directive = match.group(1).lower()
            out.append(raw)
            continue
        stripped = raw.split(";", 1)[0].strip()
        if not stripped:
            out.append(raw)
            continue
        if directive == "moleculetype" and not seen_name:
            fields = stripped.split()
            nrexcl = fields[1] if len(fields) > 1 else "3"
            out.append(f"{new:<20s} {nrexcl}\n")
            seen_name = True
            continue
        if directive == "atoms":
            out.append(_rewrite_atom_line(raw, old, new, rename_atom))
            continue
        out.append(raw)
    return "".join(out)


def _rewrite_atom_line(raw: str, old: str, new: str, rename_atom: bool) -> str:
    """Replace the residue (and optionally atom) name columns of an atoms line."""
    body, sep, comment = raw.partition(";")
    fields = body.split()
    if len(fields) < 6:
        return raw
    if fields[3] == old:
        fields[3] = new
    if rename_atom and fields[4] == old:
        fields[4] = new
    rebuilt = (
        f"{fields[0]:>6s} {fields[1]:>10s} {fields[2]:>6s} {fields[3]:>6s} "
        f"{fields[4]:>6s} {fields[5]:>6s}"
    )
    if len(fields) > 6:
        rebuilt += " " + " ".join(fields[6:])
    return rebuilt + ("   " + sep + comment if sep else "\n")


def _count_atoms(block_text: str) -> int:
    directive = ""
    count = 0
    for raw in block_text.splitlines():
        match = _DIRECTIVE_RE.match(raw)
        if match:
            directive = match.group(1).lower()
            continue
        if raw.lstrip().startswith("#"):
            continue
        if directive == "atoms" and raw.split(";", 1)[0].strip():
            count += 1
    return count


def _check_site_count(water_key: str, sites: int) -> None:
    """Cross-check the probed water against the ``-cs`` box gmx solvate will use.

    Imported from :mod:`.solvate` rather than duplicated, so the two can never
    drift into disagreeing about what ``opc`` means.
    """
    try:
        from . import solvate as _solvate

        coord_file = getattr(_solvate, "_WATER_MODEL_COORD_FILE", {}).get(water_key)
        expected = getattr(_solvate, "_COORD_FILE_SITES", {}).get(coord_file)
    except Exception:  # pragma: no cover
        return
    if expected is None:
        return
    if sites != expected:
        raise SolventParameterError(
            f"tleap produced a {sites}-site water for {water_key!r}, but gmx solvate fills that "
            f"model from {coord_file} ({expected} sites). The topology and the coordinates would "
            "disagree about every water molecule. " + _PACKMOL_HINT
        )


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #


def assert_solvent_blocks(blocks: SolventBlocks) -> None:
    """Refuse blocks that would not survive ``gmx solvate``/``gmx genion``."""
    if not isinstance(blocks, SolventBlocks):
        raise SolventParameterError(
            f"Expected SolventBlocks, got {type(blocks).__name__}."
        )
    names = dict(blocks.moleculetypes)
    if blocks.water_moleculetype != WATER_MOLECULETYPE:
        raise SolventParameterError(
            f"Water moleculetype must be named {WATER_MOLECULETYPE!r} because that is the residue "
            f"name gmx solvate writes; got {blocks.water_moleculetype!r}."
        )
    for role, name in (
        ("water", blocks.water_moleculetype),
        ("positive ion", blocks.cation_moleculetype),
        ("negative ion", blocks.anion_moleculetype),
    ):
        if name not in names:
            raise SolventParameterError(
                f"No [ moleculetype ] block was probed for the {role} ({name!r}); "
                f"got: {', '.join(sorted(names)) or '<none>'}."
            )
    if blocks.water_sites not in (3, 4, 5):
        raise SolventParameterError(
            f"Implausible water site count {blocks.water_sites!r}."
        )
    if not blocks.atomtypes:
        raise SolventParameterError("The solvent probe returned no [ atomtypes ] entries.")
    if blocks.comb_rule != 2:
        raise SolventParameterError(
            f"The solvent probe used combination rule {blocks.comb_rule}, but AMBER-derived "
            "membrane topologies use rule 2 (sigma/epsilon)."
        )

    # Charge sanity: a mis-mapped ion unit (asking for Na+ and getting Ca2+)
    # would otherwise only show up as a system that will not neutralise.
    for role, name, expected in (
        ("water", blocks.water_moleculetype, 0.0),
        ("positive ion", blocks.cation_moleculetype, 1.0),
        ("negative ion", blocks.anion_moleculetype, -1.0),
    ):
        charge = blocks.charges.get(name)
        if charge is None:
            continue
        if role == "water":
            if abs(charge) > 1e-3:
                raise SolventParameterError(
                    f"Probed water carries a net charge of {charge:+.3f} e; it must be neutral."
                )
            continue
        # Only the sign is enforced: a divalent counter-ion is a legitimate
        # (if unusual) request, and genion handles it by charge.
        if charge * expected <= 0:
            raise SolventParameterError(
                f"Probed {role} {name!r} carries charge {charge:+.3f} e, which has the wrong sign."
            )


# --------------------------------------------------------------------------- #
# Merging
# --------------------------------------------------------------------------- #


def merge_solvent_blocks(top_path: str, blocks: SolventBlocks) -> None:
    """Splice *blocks* into the dry topology at *top_path*, in place.

    Adds the missing ``[ atomtypes ]`` entries and appends the water/ion
    ``[ moleculetype ]`` blocks immediately before ``[ system ]``.
    ``[ molecules ]`` is untouched -- ``gmx solvate`` and ``gmx genion`` own it.

    Idempotent: a topology that already carries all three moleculetypes is left
    alone, so a retry cannot duplicate them (which ``read_topology`` rejects).
    """
    if blocks is None:
        raise SolventParameterError("No solvent parameter blocks were provided.")
    top_path = os.path.abspath(top_path)
    if not os.path.isfile(top_path):
        raise SolventParameterError(f"Topology not found: {top_path}")

    with open(top_path) as handle:
        lines = handle.readlines()

    existing_types, existing_mols, system_at, atomtypes_end, comb_rule = _scan_topology(lines)

    wanted = {
        blocks.water_moleculetype,
        blocks.cation_moleculetype,
        blocks.anion_moleculetype,
    }
    if wanted <= existing_mols:
        return  # already merged

    clashes = wanted & existing_mols
    if clashes:
        raise SolventParameterError(
            f"{os.path.basename(top_path)} already defines "
            f"{', '.join(sorted(clashes))} but not the rest of the solvent blocks; "
            "refusing to produce a partially-merged topology. " + _PACKMOL_HINT
        )

    if comb_rule is not None and comb_rule != blocks.comb_rule:
        raise SolventParameterError(
            f"{os.path.basename(top_path)} uses combination rule {comb_rule} but the probed "
            f"solvent parameters use rule {blocks.comb_rule}; sigma/epsilon and C6/C12 are not "
            "interchangeable. " + _PACKMOL_HINT
        )
    if atomtypes_end is None:
        raise SolventParameterError(
            f"{os.path.basename(top_path)} has no [ atomtypes ] section to merge into; it is not "
            "a standalone topology. " + _PACKMOL_HINT
        )
    if system_at is None:
        raise SolventParameterError(
            f"{os.path.basename(top_path)} has no [ system ] section, so there is nowhere to "
            "insert the solvent moleculetypes."
        )

    new_types: List[str] = []
    for name, line in blocks.atomtypes:
        if name in existing_types:
            _assert_same_atomtype(name, existing_types[name], line, top_path)
            continue
        new_types.append(line)

    molecule_text = (
        "\n; --- solvent parameters added by PRISM (" + (blocks.source or "tleap probe") + ") ---\n"
        + "\n".join(text.rstrip() + "\n" for _name, text in blocks.moleculetypes)
        + "\n"
    )

    # Insert from the bottom up so the earlier index stays valid.
    lines[system_at:system_at] = [molecule_text]
    if new_types:
        lines[atomtypes_end:atomtypes_end] = [
            "; solvent atom types added by PRISM (" + (blocks.source or "tleap probe") + ")\n"
        ] + new_types

    with open(top_path, "w") as handle:
        handle.writelines(lines)


def _scan_topology(
    lines: Sequence[str],
) -> Tuple[Dict[str, str], set, Optional[int], Optional[int], Optional[int]]:
    """Return (atomtypes, moleculetype names, [system] index, atomtypes end, comb-rule)."""
    atomtypes: Dict[str, str] = {}
    moleculetypes: set = set()
    system_at: Optional[int] = None
    atomtypes_end: Optional[int] = None
    comb_rule: Optional[int] = None

    directive = ""
    expect_name = False
    for i, raw in enumerate(lines):
        match = _DIRECTIVE_RE.match(raw)
        if match:
            if directive == "atomtypes":
                atomtypes_end = i
            directive = match.group(1).lower()
            expect_name = directive == "moleculetype"
            if directive == "system" and system_at is None:
                system_at = i
            continue
        if raw.lstrip().startswith("#"):
            continue
        stripped = raw.split(";", 1)[0].strip()
        if not stripped:
            continue
        fields = stripped.split()
        if directive == "atomtypes":
            atomtypes[fields[0]] = stripped
        elif directive == "moleculetype" and expect_name:
            moleculetypes.add(fields[0])
            expect_name = False
        elif directive == "defaults" and comb_rule is None and len(fields) >= 2:
            try:
                comb_rule = int(fields[1])
            except ValueError:
                comb_rule = None
    if directive == "atomtypes":
        atomtypes_end = len(lines)
    return atomtypes, moleculetypes, system_at, atomtypes_end, comb_rule


def _assert_same_atomtype(name: str, existing: str, incoming: str, top_path: str) -> None:
    """Refuse to reuse a same-named atomtype that carries different parameters.

    Silently reusing the host's version would give the water or the ions the
    Lennard-Jones parameters of whatever lipid atom happened to share the name.
    """
    have = _atomtype_params(existing)
    want = _atomtype_params(incoming)
    if have is None or want is None:
        return
    for label, a, b in zip(("mass", "sigma", "epsilon"), have, want):
        if abs(a - b) > _PARAM_RTOL * max(1.0, abs(a), abs(b)):
            raise SolventParameterError(
                f"{os.path.basename(top_path)} already defines atom type {name!r} with a "
                f"different {label} ({a:g} vs the probed {b:g}). Merging would silently "
                "reparametrise one of them. " + _PACKMOL_HINT
            )


def _atomtype_params(line: str) -> Optional[Tuple[float, float, float]]:
    """Extract (mass, sigma, epsilon) from an ``[ atomtypes ]`` line."""
    fields = line.split(";", 1)[0].split()
    numbers: List[float] = []
    for token in fields[1:]:
        try:
            numbers.append(float(token))
        except ValueError:
            continue
    if len(numbers) < 3:
        return None
    # ParmEd writes: name at.num mass charge ptype sigma epsilon.  Mass is the
    # second number; sigma/epsilon are always the last two.
    mass = numbers[1] if len(numbers) >= 4 else numbers[0]
    return mass, numbers[-2], numbers[-1]


def _read_tail(path: str, lines: int = 20) -> str:
    try:
        with open(path) as handle:
            return _tail(handle.read(), lines)
    except OSError:
        return ""


def _tail(text: str, lines: int = 20) -> str:
    if not text:
        return ""
    return "\n".join(text.strip().splitlines()[-lines:])
