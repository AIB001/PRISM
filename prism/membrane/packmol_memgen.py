#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PACKMOL-Memgen backend: pack an all-atom lipid bilayer around a protein.

PACKMOL-Memgen ships with AmberTools and is fully CLI-driven. With
``--parametrize`` it emits an AMBER ``prmtop``/``crd`` pair (Lipid21 + the chosen
protein FF), which PRISM converts to GROMACS through ParmEd — reusing the same
``amb2gmx`` philosophy as the GAFF ligand path.

This module builds the command and (when AmberTools is installed) runs it. The
command construction is pure and unit-testable; the run is runtime-gated.

Two packing modes are supported:

* the default *wet* pack, in which PACKMOL places protein, lipids, water and
  ions and tleap parametrizes the finished system.  This is the validated path
  and nothing about it changes;
* an opt-in *dry* pack (``dry=True``), which packs only the protein and the
  lipids and leaves the water slab empty for ``gmx solvate``/``gmx genion`` to
  fill downstream.  Water is 99% of the molecules PACKMOL is asked to place on
  a realistic membrane system -- 131 635 of 133 015 on a class-C GPCR dimer --
  and PACKMOL places every one of them through GENCAN optimisation while
  ``gmx solvate`` stamps a pre-equilibrated box on a grid in seconds.  Removing
  them is what turns "killed before finishing" into a bounded build.

The dry pack is a *flag* change, not a hand-edited PACKMOL input: PRISM keeps
PACKMOL-Memgen's own tleap generation (leaprc selection, ``--ligand_param``,
``--leapline`` disulfides, box assignment) instead of reimplementing it.
"""

from __future__ import annotations

import importlib.util
import math
import os
import shutil
import struct
import subprocess
from dataclasses import dataclass
from numbers import Integral
from typing import Dict, List, Optional

# ORIENTATION_METHODS is imported rather than restated: a hand-copied set here
# silently rejected every new method config learned (``auto`` was the case that
# exposed it), and it rejected them deep inside command construction, after the
# solute had already been cleaned and oriented.
from .config import (
    MembraneConfig,
    ORIENTATION_METHODS,
    PACKMOL_PARAMETRIZED_LIPID_FFS,
    PDB_MIN_COORDINATE_WIDTH,
)


@dataclass
class PackmolMemgenResult:
    prmtop: Optional[str]
    inpcrd: Optional[str]
    packed_pdb: Optional[str]
    log: str
    success: bool
    #: True when this pack deliberately contains no water and no ions, so the
    #: caller must solvate it before the system is usable.  Defaulted so
    #: existing positional/keyword construction keeps working.
    dry: bool = False


#: Prefix for PACKMOL's own input/log files (``--packlog``).  PRISM sets it
#: explicitly rather than relying on the upstream default because it reads the
#: log back: packmol-memgen redirects the PACKMOL binary's stdout into this file
#: and re-emits only progress strings, so the packing verdict never reaches the
#: output PRISM captures.
PACKMOL_LOG_PREFIX = "packmol"

#: Where ``run()`` persists the captured tool output.  Distinct from
#: packmol-memgen's own ``packmol-memgen.log`` (hyphen), which it writes itself.
PRISM_LOG_NAME = "packmol_memgen.log"

#: Minimization engines in the order PRISM prefers them.  packmol-memgen
#: defaults to ``pmemd.cuda``, but its own help says "It's highly recommended
#: that the minimization steps are performed on cpus" because GPU minimization
#: of a freshly packed bilayer is a known cause of first-step crashes -- and a
#: present ``pmemd.cuda`` binary does not imply a usable GPU.
MINIMIZATION_ENGINES = ("pmemd", "sander", "pmemd.cuda")

_INSTALL_HINT = "Install: conda install -c conda-forge ambertools"

#: Flags that together make PACKMOL-Memgen pack a *dry* system.  All three are
#: required and none of them is optional:
#:
#: ``--solvent_parm``  overrides the density PACKMOL-Memgen uses to convert the
#:                     water-slab volume into a molecule count.  Its solvent
#:                     table is appended to (and therefore overrides) the
#:                     built-in one, and the count is ``int(volume * density *
#:                     N_A / MW)``, so a near-zero density makes it exactly 0.
#: ``--nocounter``     drops the ion blocks from the PACKMOL input *and* the
#:                     ``addionsrand`` lines from the generated tleap script.
#:                     Without it tleap aborts outright: ``addionsrand`` on a
#:                     solvent-free system is a fatal error, and no prmtop is
#:                     written at all.
#: ``--tight_box``     makes tleap ``set SYS box {x y z}`` from the packing
#:                     restraints instead of ``setBox SYS vdw``, which shrinks
#:                     the box onto the atoms and would delete the empty water
#:                     slab we still need.
DRY_PACK_OPTIONS = ("--nocounter", "--tight_box", "--solvent_parm")

#: Name of the one-line solvent table PRISM writes for a dry pack.
DRY_SOLVENT_PARM_NAME = "solvent_dry.parm"

#: Solvent PRISM leaves empty.  PACKMOL-Memgen's ``--solvents`` default (and
#: PRISM never overrides it), so this is the only entry that has to be shadowed.
DRY_SOLVENT_NAME = "WAT"

#: Density (g/cm^3) written for that entry.  It has to be small enough that
#: ``int(box_volume * density * N_A / MW)`` truncates to zero for any box PRISM
#: can build, and large enough that PACKMOL-Memgen's ``watnum / solvent_con``
#: division does not blow up.  At 1e-12 g/cm^3 the count stays 0 up to a
#: 3e13 A^3 box -- roughly a million times the volume of the largest membrane
#: system anyone builds -- while the divisor is still ~3e-14.
DRY_SOLVENT_DENSITY_G_CM3 = 1.0e-12

#: Molar mass PACKMOL-Memgen lists for ``WAT``; kept identical so only the
#: density differs from the shipped table.
DRY_SOLVENT_MOLAR_MASS = 18

#: Avogadro's number, and the water density PACKMOL-Memgen's own solvent table
#: uses.  Declared here because :func:`estimate_water_count` reproduces
#: PACKMOL-Memgen's molecule-count arithmetic and must not drift from it.
AVOGADRO = 6.02214076e23
WATER_DENSITY_G_CM3 = 0.997

#: PACKMOL-Memgen's ``--leaflet`` default (Angstrom, half the bilayer
#: thickness).  PRISM never passes ``--leaflet``, so this is what its box
#: sizing uses, and :func:`estimate_water_count` has to use the same value.
PACKMOL_LEAFLET_HALF_THICKNESS_A = 23.0

#: Residue labels that mean "bulk solvent" in an AMBER topology.  Matched
#: case-insensitively.
_BULK_WATER_RESIDUES = frozenset({"WAT", "HOH", "SOL", "TIP3", "T3P", "TP3", "T4P"})

#: Ion residue labels, matched CASE-SENSITIVELY and exactly.  A substring or
#: suffix test would be wrong: Lipid21 ships a residue literally named ``PH-``
#: (deprotonated phosphatidic acid head), so "ends with a minus sign" would
#: report every phosphatidic-acid lipid as a stray chloride.
_ION_RESIDUES = frozenset({
    "Na+", "Cl-", "K+", "Li+", "Rb+", "Cs+", "F-", "Br-", "IOD", "NH4", "H3O+",
})


def is_available() -> bool:
    """True if packmol-memgen is on PATH."""
    return shutil.which("packmol-memgen") is not None


def missing_requirements(config: Optional[MembraneConfig] = None) -> List[str]:
    """Return actionable messages for every prerequisite a real build needs.

    ``is_available()`` only answers whether the wrapper script exists.  The
    wrapper then resolves ``packmol``, ``tleap`` and the minimization engine
    from its own bundle or from ``$AMBERHOME/bin`` -- never from PATH -- and a
    missing one of those surfaces as an exit-1 line buried in the log *after*
    the solute has been cleaned, protonated and packed.  Reporting them
    together up front mirrors the GAFF generator's single "install AmberTools"
    message instead of discovering them one expensive run at a time.
    """
    if not is_available():
        return [f"packmol-memgen not found (provided by AmberTools). {_INSTALL_HINT}"]

    problems = []
    amberhome = _resolve_amberhome()
    if not amberhome:
        problems.append(
            "No AmberTools root found ($AMBERHOME is unset or does not contain dat/leap/cmd), "
            f"so tleap cannot parametrize the packed system. {_INSTALL_HINT}"
        )
    if _resolve_packmol(amberhome) is None and shutil.which("packmol") is None:
        problems.append(
            "packmol executable not found; packmol-memgen looks only in its own bundle and in "
            f"$AMBERHOME/bin. {_INSTALL_HINT}"
        )
    if amberhome and not os.path.exists(os.path.join(amberhome, "bin", "tleap")):
        problems.append(
            f"tleap not found in {os.path.join(amberhome, 'bin')}; the bilayer can be packed but "
            f"not parametrized. {_INSTALL_HINT}"
        )
    if getattr(config, "minimize", False) and _resolve_minimization_engine(amberhome) is None:
        problems.append(
            "No AMBER minimization engine (pmemd/sander) found in $AMBERHOME/bin, but membrane "
            f"'minimize' is enabled. {_INSTALL_HINT}, or set membrane 'minimize: false'."
        )
    # Only probe --help when a dry pack is actually requested: the probe costs
    # ~9 s because packmol-memgen imports its whole toolchain, and the default
    # (wet) path does not need the answer.  ``getattr`` because the 'solvate'
    # field is optional -- a directly constructed MembraneConfig from an older
    # caller simply does not have it.
    if str(getattr(config, "solvate", "packmol") or "").strip().lower() == "gromacs":
        try:
            unsupported = unsupported_dry_options()
        except RuntimeError as exc:  # --help unreadable; report, do not crash
            problems.append(str(exc))
        else:
            if unsupported:
                problems.append(
                    "Installed packmol-memgen does not advertise "
                    f"{', '.join(unsupported)}, so it cannot pack a lipids-only (dry) system. "
                    "Upgrade AmberTools or use membrane 'solvate: packmol'."
                )
    return problems


def _looks_like_amber_root(path) -> bool:
    """True when ``path`` is an AmberTools tree rather than an arbitrary prefix."""
    # dat/leap/cmd is what tleap needs and what distinguishes a real AmberTools
    # root from, say, /usr/local (which is what a standalone packmol install
    # would otherwise yield).
    return bool(path) and os.path.isdir(os.path.join(path, "dat", "leap", "cmd"))


def _resolve_amberhome() -> Optional[str]:
    """Return a validated AmberTools root, or None.

    packmol-memgen resolves tleap, reduce, memembed and the minimization engine
    from ``$AMBERHOME/bin`` and nothing else, so this value decides whether the
    build can be parametrized at all.  Derive it from packmol-memgen (or tleap)
    rather than from ``packmol``: packmol is commonly installed standalone in
    ``/usr/local/bin`` while AmberTools lives in a conda env, and that root
    would silently disable the legacy-lipid ``leaprc`` fix below and mis-resolve
    every other helper.
    """
    declared = os.environ.get("AMBERHOME") or None
    if _looks_like_amber_root(declared):
        return declared
    for executable in ("packmol-memgen", "tleap", "packmol"):
        path = shutil.which(executable)
        if not path:
            continue
        root = os.path.dirname(os.path.dirname(path))
        if _looks_like_amber_root(root):
            return root
    # Never erase an explicitly declared AMBERHOME just because it failed the
    # shape check; the user may know something this heuristic does not.
    return declared


def _resolve_packmol(amberhome) -> Optional[str]:
    """Return the ``packmol`` binary packmol-memgen itself would find."""
    try:
        spec = importlib.util.find_spec("packmol_memgen")
    except (ImportError, ValueError):
        spec = None
    for location in (spec.submodule_search_locations if spec else None) or ():
        bundled = os.path.join(location, "lib", "packmol", "packmol")
        if os.path.exists(bundled):
            return bundled
    if amberhome:
        candidate = os.path.join(amberhome, "bin", "packmol")
        if os.path.exists(candidate):
            return candidate
    return None


def _resolve_minimization_engine(amberhome) -> Optional[str]:
    """Return the AMBER engine to pass to ``--engine``, or None if none exists."""
    if not amberhome:
        return None
    for engine in MINIMIZATION_ENGINES:
        if os.path.exists(os.path.join(amberhome, "bin", engine)):
            return engine
    return None


#: PACKMOL-Memgen joins the frcmod and lib filenames of ``--ligand_param`` with
#: this character and splits the value back apart on it.
LIGAND_PARAM_SEPARATOR = ":"


def _packmol_ion_name(name: str) -> str:
    """Translate common GROMACS ion tokens to PACKMOL-Memgen names."""
    aliases = {
        "NA": "Na+",
        "NA+": "Na+",
        "CL": "Cl-",
        "CL-": "Cl-",
        "K": "K+",
        "K+": "K+",
    }
    token = str(name).strip()
    return aliases.get(token.upper(), token)


def _ligand_parameter_paths(config: MembraneConfig):
    """Return the absolute ``(frcmod, lib)`` paths to feed ``--ligand_param``.

    Absolutisation is not cosmetic here. ``run()`` launches packmol-memgen with
    ``cwd=work_dir``, and packmol-memgen passes each half of ``--ligand_param``
    straight into a tleap script as ``loadamberparams``/``loadoff`` arguments,
    which tleap resolves against that working directory. A caller-relative path
    would therefore silently point somewhere else -- the same trap the protein
    PDB is normalized for in ``run()``.

    ``expanduser`` runs first because membrane configs routinely carry
    ``~/params/lig.frcmod`` and ``abspath()`` would happily turn that into
    ``<cwd>/~/params/lig.frcmod``.

    This function only manipulates strings (no stat, no read), so
    ``build_command`` stays filesystem-free and unit-testable; existence is
    checked in ``run()``.
    """

    def resolve(value, field: str) -> str:
        path = os.path.abspath(os.path.expanduser(str(value).strip()))
        # PACKMOL-Memgen does `frcmod, lib = value.split(":")` and tleap
        # tokenizes its script lines on whitespace. Either character inside a
        # path does not fail loudly -- it silently resolves to a *different*,
        # wrong filename (or a truncated one), and the resulting topology would
        # be missing the ligand parameters. Refuse to emit a corrupt command.
        if LIGAND_PARAM_SEPARATOR in path:
            raise ValueError(
                f"Membrane {field} path may not contain {LIGAND_PARAM_SEPARATOR!r} because "
                f"PACKMOL-Memgen splits --ligand_param on it: {path!r}. "
                "Move the file to a path without that character."
            )
        if any(character.isspace() for character in path):
            raise ValueError(
                f"Membrane {field} path may not contain whitespace because tleap splits its "
                f"loadamberparams/loadoff lines on it: {path!r}. "
                "Move the file to a path without spaces."
            )
        return path

    return (
        resolve(config.ligand_frcmod, "ligand_frcmod"),
        resolve(config.ligand_lib, "ligand_lib"),
    )


def tleap_residue_numbers(pdb_path: str) -> Dict[tuple, Optional[tuple]]:
    """Map ``(chain, resSeq)`` to ``(tleap residue number, residue name)``.

    tleap does not answer to the residue numbers a PDB writes.  Measured against
    AmberTools' teLeap: it keeps the FIRST residue's sequence number and then
    increments by one per residue, ignoring the file's own gaps and the restart
    at a new chain.  A three-residue file numbered 41/50/60 becomes 41/42/43 --
    ``SYS.50`` does not exist, and a ``bond`` built from the file's own numbers
    either fails outright or lands on whatever residue holds that number.

    So the rule is ``first_resSeq + ordinal - 1``, and it is applied here rather
    than assumed to be the identity.  The patch route renumbers its solute 1..N
    before tleap sees it (:meth:`builder._write_tleap_ready_solute`), which makes
    the offset zero and the mapping positional; the PACKMOL route hands
    packmol-memgen the oriented structure with the depositor's numbering intact,
    where the offset is whatever the first residue happens to be -- 41 on 9OMO.
    One implementation covers both, because two would only agree until they did
    not.

    A key that appears more than once -- an insertion code, a chain whose
    numbering restarts, atoms of one residue split across the file -- maps to
    ``None`` so the caller refuses rather than guesses.
    """
    ordinals: Dict[tuple, Optional[tuple]] = {}
    ordinal = 0
    first_number: Optional[int] = None
    previous = None
    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            # A record too short to carry coordinates is not a residue tleap can
            # place, and it is not one prism.ptm.disulfides can detect an SG in
            # either (it parses the same columns).  Both sides therefore ignore
            # exactly the same records, which is what keeps the ordinals aligned
            # instead of shifting them by one per malformed line.  Indexing a
            # short line would otherwise raise a bare IndexError from the middle
            # of a build.
            if len(line) < PDB_MIN_COORDINATE_WIDTH:
                continue
            # chain + resSeq + insertion code: the full residue identity.
            key = (line[21], line[22:27])
            if key == previous:
                continue
            previous = key
            ordinal += 1
            try:
                short = (line[21].strip() or "A", int(line[22:26]))
            except ValueError:
                continue
            if first_number is None:
                first_number = short[1]
            number = first_number + ordinal - 1
            ordinals[short] = None if short in ordinals else (number, line[17:20].strip())
    return ordinals


def detect_disulfide_pairs(pdb_path: str):
    """Detect disulfides in *pdb_path* as ``(pairs, descriptions, unresolved)``.

    ``pairs`` are residue numbers in tleap's numbering of *this file*, ready for
    :func:`disulfide_leaplines`.

    This wraps :func:`prism.ptm.disulfides.detect_disulfides` -- PRISM's single
    disulfide detector, already used by the patch route -- rather than adding a
    second one.  The geometry is therefore identical on both routes by
    construction, which is the only way "both routes declare the same bonds"
    stays true after the next change to either.

    The residue sitting at each computed number must still be a cysteine.  That
    check is what makes an assumption violation loud: if anything shifted the
    ordinals, the bond would otherwise be declared against whatever residue
    happens to hold that number, which is silently wrong in a way no downstream
    stage can detect.
    """
    try:
        from ..ptm.disulfides import detect_disulfides
    except ImportError as exc:  # pragma: no cover - only on a broken install
        raise RuntimeError(
            f"prism.ptm.disulfides is not importable ({exc}), so the disulfides of the solute "
            "cannot be declared to tleap. Building without them would silently produce free "
            "cysteines with overlapping thiol hydrogens."
        ) from exc

    bonds = detect_disulfides(pdb_path)
    numbering = tleap_residue_numbers(pdb_path)
    pairs, described, unresolved = [], [], []
    for bond in bonds:
        numbers = []
        for sulfur in (bond.a, bond.b):
            entry = numbering.get((sulfur.chain, sulfur.resid))
            if entry is None or not str(entry[1]).upper().startswith("CY"):
                break
            numbers.append(entry[0])
        if len(numbers) != 2:
            unresolved.append(bond.describe())
            continue
        pairs.append((numbers[0], numbers[1]))
        described.append(
            f"{bond.a.chain}{bond.a.resid}-{bond.b.chain}{bond.b.resid} "
            f"(SG-SG {bond.distance:.2f} A; tleap {numbers[0]}-{numbers[1]})"
        )
    return pairs, described, unresolved


def disulfide_leaplines(pairs) -> List[str]:
    """Render ``--leapline`` bond commands for a solute's disulfide pairs.

    PACKMOL-Memgen cannot infer these: its tleap generator emits ``bond`` lines
    only for cardiolipin, while its protonation step renames an SS cysteine to
    CYX (dropping HG) without ever creating the S-S bond.  The result is a
    dangling sulfur carrying the CYX charge that assumes a bonded partner --
    silent, and present in nearly every GPCR/channel/transporter this module
    targets.

    ``pairs`` are 1-based residue numbers *as tleap numbers the packed solute*
    (protein first), which is what ``prism.ptm.disulfides.detect_disulfides``
    reports for the oriented solute PDB.  They are validated rather than
    formatted blindly: a wrong number bonds the wrong residues just as quietly
    as no bond at all.
    """
    if not pairs:
        return []
    lines = []
    seen = set()
    for pair in pairs:
        try:
            first, second = pair
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Membrane disulfide entry must be a pair of residue numbers: {pair!r}"
            ) from exc
        numbers = []
        for value in (first, second):
            if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
                raise ValueError(
                    f"Membrane disulfide residue numbers must be positive integers: {pair!r}"
                )
            numbers.append(int(value))
        if numbers[0] == numbers[1]:
            raise ValueError(f"Membrane disulfide pair bonds a residue to itself: {pair!r}")
        key = tuple(sorted(numbers))
        if key in seen:
            continue  # tleap rejects a duplicated bond outright
        seen.add(key)
        lines.append(f"bond SYS.{numbers[0]}.SG SYS.{numbers[1]}.SG")
    return lines


def write_dry_solvent_parm(work_dir: str) -> str:
    """Write the near-zero-density solvent table that empties the water slab.

    PACKMOL-Memgen reads its shipped ``solvent.parm`` and then *appends* the
    file given to ``--solvent_parm``, so a later line for the same solvent wins.
    It converts the slab volume into a molecule count with
    ``int(volume * density * N_A / MW)``, which truncates to zero at this
    density -- so the flag is a documented, single-purpose gate rather than a
    coincidence of the arithmetic.

    Only the density is changed.  The molar mass stays at PACKMOL-Memgen's own
    value because it appears in the divisor of a later expression as well;
    zeroing or omitting it would raise instead of packing nothing.

    Returns the absolute path, which is what must be handed to
    ``build_command(solvent_parm=...)``: the subprocess runs with
    ``cwd=work_dir`` and a relative path would resolve somewhere else.
    """
    work_dir = os.path.abspath(work_dir)
    os.makedirs(work_dir, exist_ok=True)
    path = os.path.join(work_dir, DRY_SOLVENT_PARM_NAME)
    with open(path, "w") as handle:
        handle.write(
            "# Written by PRISM for a lipids-only (dry) pack: PACKMOL places no water,\n"
            "# because 'gmx solvate' fills the slab afterwards in seconds.\n"
            "# Columns: solvent density[g/cm3] MW[g/mol] charge name source\n"
        )
        handle.write(
            f"{DRY_SOLVENT_NAME} {DRY_SOLVENT_DENSITY_G_CM3:g} "
            f"{DRY_SOLVENT_MOLAR_MASS} 0 Water prism-dry\n"
        )
    return path


def _dry_pack_options(solvent_parm: Optional[str]) -> List[str]:
    """Return the flags that turn a full pack into a lipids-only pack.

    ``solvent_parm`` is required rather than defaulted because the failure mode
    of omitting it is silent and severe: ``--nocounter --tight_box`` alone still
    packs every water, so the caller would get a fully solvated system that its
    own ``gmx solvate``/``gmx genion`` stage then solvates *again* -- water in
    the water, and a salt concentration nobody asked for.
    """
    if not solvent_parm or not str(solvent_parm).strip():
        raise ValueError(
            "A dry pack needs the near-zero-density solvent table; obtain it from "
            "write_dry_solvent_parm() rather than relying on --nocounter alone, which "
            "removes the ions but still packs every water molecule."
        )
    # cwd is the work directory during the run, so a caller-relative path would
    # resolve against the wrong root -- the same trap the protein PDB and the
    # ligand parameters are normalized for.
    path = os.path.abspath(os.path.expanduser(str(solvent_parm).strip()))
    return ["--nocounter", "--tight_box", "--solvent_parm", path]


def build_command(
    protein_pdb: str,
    config: MembraneConfig,
    parametrize: bool = True,
    protein_ff_option: Optional[str] = None,
    water_ff_option: Optional[str] = "--ffwat",
    engine: Optional[str] = None,
    packmol_path: Optional[str] = None,
    dry: bool = False,
    solvent_parm: Optional[str] = None,
    disulfide_pairs=None,
) -> List[str]:
    """Construct the packmol-memgen command line for the given config.

    The lipid argument uses PACKMOL-Memgen's ``A:B`` syntax for mixtures and the
    parallel ``--ratio`` flag for molar proportions.

    ``protein_ff_option``/``engine``/``packmol_path`` describe the *installed*
    release and are therefore inputs rather than guesses -- ``run()`` detects
    them and passes them in, keeping this function filesystem-free and
    unit-testable.  There is no safe default for ``protein_ff_option`` because
    AmberTools releases genuinely split between ``--fprot`` and ``--ffprot``;
    ``--ffwat`` has no such split.

    ``dry`` packs protein and lipids only and leaves the water slab empty for a
    downstream ``gmx solvate``/``gmx genion`` stage; it needs ``solvent_parm``
    from :func:`write_dry_solvent_parm`.  Everything else -- lipids, ratio,
    ``--dist``, ``--dist_wat``, orientation, disulfides, ligand parameters,
    force fields, iteration budgets, minimization -- is identical in both modes,
    so a dry pack is the same bilayer with the water left out.  In particular
    ``--dist_wat`` is still passed, because it is what sizes the empty slab.
    """
    config.raise_for_errors()
    lipids = ":".join(config.lipids)
    ratio = ":".join(str(r) for r in config.resolved_ratio())
    orient = (config.orient or "").lower()
    if orient not in ORIENTATION_METHODS:
        raise ValueError(f"Unsupported membrane orientation method: {config.orient!r}")

    cmd = [
        "packmol-memgen",
        "--pdb", protein_pdb,
        "--lipids", lipids,
        "--ratio", ratio,
        "--dist", f"{config.xy_distance:g}",
        "--dist_wat", f"{config.water_thickness:g}",
    ]
    if dry:
        # No salt flags at all on this path.  They are not merely redundant:
        # PACKMOL-Memgen logs an error line for "--saltcon without --salt", and
        # every ion it would place has to come from ``gmx genion`` instead --
        # counting ions against a box that is 99% empty would give a
        # concentration far from the one that was requested.
        cmd.extend(_dry_pack_options(solvent_parm))
    else:
        cmd.extend([
            "--salt",
            "--salt_c", _packmol_ion_name(config.positive_ion),
            "--salt_a", _packmol_ion_name(config.negative_ion),
            "--saltcon", f"{config.salt_concentration:g}",
        ])
    cmd.extend([
        # Pinned, not defaulted: PRISM parses this log back to recover PACKMOL's
        # packing verdict (see ``_packing_verdict_note``).
        "--packlog", PACKMOL_LOG_PREFIX,
    ])
    # packmol-memgen searches its own bundle and $AMBERHOME/bin for the packmol
    # binary, never PATH.  Hand it an explicit path when the caller found one
    # somewhere else, instead of letting it prepare the whole solute and then
    # exit 1 with "Packmol path not defined".
    if packmol_path:
        cmd.extend(["--packmol", packmol_path])
    # Orientation: if the protein is already membrane-aligned, tell packmol-memgen
    # not to re-orient; otherwise let its built-in memembed run.
    if orient != "memembed":
        cmd.append("--preoriented")
    # Outside the ``parametrize`` branch for the same reason as the ligand flags
    # below: tleap also runs behind --minimize, and a missing S-S bond there is
    # just as wrong.
    #
    # ``disulfide_pairs`` is a parameter rather than something read off the
    # filesystem here because this function must stay pure -- it is called with
    # paths that need not exist.  :func:`run` does the detection and passes the
    # result in; anything configured explicitly is merged with it, and
    # ``disulfide_leaplines`` de-duplicates.
    detected = list(disulfide_pairs or ())
    configured = list(getattr(config, "disulfide_pairs", None) or ())
    for leapline in disulfide_leaplines(detected + configured):
        cmd.extend(["--leapline", leapline])
    # Ligand flags are deliberately outside the ``parametrize`` branch:
    # --keepligs decides what gets packed, and --ligand_param also feeds the
    # tleap run behind --minimize ("in case of parametrizing OR minimizing the
    # system with non-canonical molecules"). Gating them on --parametrize would
    # drop the ligand from a pack-only or minimize-only build.
    if config.has_ligand():
        frcmod, lib = _ligand_parameter_paths(config)
        # MEMEMBED cleans the input PDB while orienting, which strips the
        # ligand out before packing. This is inert on the --preoriented path
        # (MEMEMBED never runs there), but it costs nothing and keeps
        # orient=memembed correct instead of silently ligand-free.
        cmd.append("--keepligs")
        cmd.extend(["--ligand_param", f"{frcmod}{LIGAND_PARAM_SEPARATOR}{lib}"])
        # --gaff2 is documented as meaningful only "if ligand parameters are
        # included" (packmol-memgen ignores it otherwise), so it stays nested
        # here rather than becoming a free-standing option. GAFF is the default.
        if config.gaff2:
            cmd.append("--gaff2")
    if parametrize:
        lipid_ff = (config.lipid_ff or "").lower()
        if lipid_ff not in PACKMOL_PARAMETRIZED_LIPID_FFS:
            raise ValueError(
                f"PACKMOL-Memgen AMBER parametrization does not support lipid_ff={config.lipid_ff!r}; "
                "use lipid17 or lipid21."
            )
        # PACKMOL-Memgen historically defaults to Lipid17.  Pass the requested
        # force field explicitly so PRISM's Lipid21 default is truthful.
        cmd.extend(["--fflip", lipid_ff])
        protein_ff, water_ff = config.packmol_forcefield_pair()
        # ff14SB/TIP3P are PACKMOL-Memgen's documented defaults.  Omitting
        # their flags avoids the historical --fprot/--ffprot spelling split;
        # non-default pairs are passed explicitly using flags detected from
        # the installed executable's help text.
        if (protein_ff, water_ff) != ("ff14SB", "tip3p"):
            if not protein_ff_option or not water_ff_option:
                raise ValueError(
                    f"Requesting {protein_ff}/{water_ff} needs the installed release's "
                    "force-field option spellings; obtain them from "
                    "_detect_forcefield_options() rather than guessing."
                )
            cmd.extend([protein_ff_option, protein_ff, water_ff_option, water_ff])
        cmd.append("--parametrize")
    # GENCAN iteration budgets.  The all-together stage dominates wall clock,
    # because each of its iterations evaluates overlaps across every molecule
    # in the box rather than within one species.
    if getattr(config, "packing_nloop", None) is not None:
        cmd.extend(["--nloop", str(int(config.packing_nloop))])
    if getattr(config, "packing_nloop_all", None) is not None:
        cmd.extend(["--nloop_all", str(int(config.packing_nloop_all))])
    if config.minimize:
        cmd.append("--minimize")
        # packmol-memgen's --engine defaults to pmemd.cuda.  Passing the engine
        # that actually exists avoids both a missing binary and GPU
        # minimization, which its own help warns can crash the first MD step.
        if engine:
            cmd.extend(["--engine", engine])
    return cmd


#: ``packmol-memgen --help`` costs ~9 s because it imports its whole toolchain,
#: and up to two independent checks consult it per build.  The installed
#: executable cannot change mid-process, so probe it at most once.
_HELP_TEXT_CACHE: Optional[str] = None


def _packmol_help_text(timeout: int = 60) -> str:
    """Return the installed executable's ``--help`` output (stdout + stderr)."""
    global _HELP_TEXT_CACHE
    if _HELP_TEXT_CACHE is not None:
        return _HELP_TEXT_CACHE
    try:
        proc = subprocess.run(
            ["packmol-memgen", "--help"],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"packmol-memgen --help did not answer within {timeout}s ({exc}); the executable "
            "is installed but unusably slow to start (a mounted Windows drive does this)."
        ) from exc
    except OSError as exc:
        raise RuntimeError(f"Could not inspect packmol-memgen --help: {exc}") from exc

    _HELP_TEXT_CACHE = (proc.stdout or "") + "\n" + (proc.stderr or "")
    return _HELP_TEXT_CACHE


def _detect_forcefield_options(timeout: int = 60):
    """Detect force-field option spellings used by the installed release.

    AmberTools releases in the wild use both ``--fprot`` and ``--ffprot``.
    Refuse to guess for a non-default force-field request because a silently
    ignored option would create a scientifically different system.
    """
    help_text = _packmol_help_text(timeout)
    protein_option = next((flag for flag in ("--fprot", "--ffprot") if flag in help_text), None)
    water_option = next((flag for flag in ("--ffwat", "--fwat") if flag in help_text), None)
    if not protein_option or not water_option:
        raise RuntimeError(
            "Installed packmol-memgen does not advertise compatible protein/water "
            "force-field options in --help; cannot safely request a non-default pair."
        )
    return protein_option, water_option


def _verify_ligand_options(cmd: List[str], timeout: int = 60) -> None:
    """Confirm the installed release advertises the ligand flags we emit.

    Unlike the ``--fprot``/``--ffprot`` split there is no spelling ambiguity to
    resolve here: every release inspected declares exactly ``--keepligs``,
    ``--ligand_param`` and ``--gaff2``, with no legacy aliases. So this is a
    capability check, not a spelling choice -- it turns argparse's
    "unrecognized arguments" (exit 2, after PRISM has already staged inputs)
    into an actionable message for an AmberTools too old to parametrize
    ligands. Same detect-don't-guess principle, applied where it pays.
    """
    ligand_flags = [flag for flag in ("--keepligs", "--ligand_param", "--gaff2") if flag in cmd]
    if not ligand_flags:
        return
    help_text = _packmol_help_text(timeout)
    missing = [flag for flag in ligand_flags if flag not in help_text]
    if missing:
        raise RuntimeError(
            "Installed packmol-memgen does not advertise "
            f"{', '.join(missing)} in --help; it is too old to build a membrane "
            "system with ligand parameters. Upgrade AmberTools."
        )


def unsupported_dry_options(timeout: int = 60) -> List[str]:
    """Return the dry-pack flags the installed release does not advertise.

    Same detect-don't-guess principle as ``_detect_forcefield_options``, applied
    where a *capability* rather than a spelling is in question.  Callers use it
    to refuse the fast path before hours of orientation and packing, instead of
    letting argparse exit 2 on "unrecognized arguments" afterwards.
    """
    help_text = _packmol_help_text(timeout)
    return [flag for flag in DRY_PACK_OPTIONS if flag not in help_text]


def supports_dry_pack(timeout: int = 60) -> bool:
    """True when the installed packmol-memgen can pack a lipids-only system."""
    if not is_available():
        return False
    try:
        return not unsupported_dry_options(timeout)
    except RuntimeError:
        # An unreadable --help means "unknown", and unknown must not enable a
        # non-default path: the wet pack is the validated one.
        return False


def _verify_dry_options(timeout: int = 60) -> None:
    """Raise unless every dry-pack flag exists in the installed release."""
    unsupported = unsupported_dry_options(timeout)
    if unsupported:
        raise RuntimeError(
            "Installed packmol-memgen does not advertise "
            f"{', '.join(unsupported)} in --help, so it cannot pack a lipids-only "
            "(dry) bilayer. Upgrade AmberTools, or build with "
            "--membrane-solvate packmol."
        )


def _require_readable_ligand_parameters(config: MembraneConfig) -> None:
    """Fail fast when the configured ligand frcmod/lib is missing or unreadable.

    tleap reports an absent parameter file as a generic parametrization failure
    buried hundreds of lines into packmol-memgen's log -- and only *after* the
    expensive packing stage has already run. Checking here (rather than in
    ``build_command``, which must stay pure and filesystem-free for unit tests)
    costs two stat calls and names the offending path directly.
    """
    for path, field in zip(_ligand_parameter_paths(config), ("ligand_frcmod", "ligand_lib")):
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"Membrane {field} not found: {path}. PACKMOL-Memgen needs both the AMBER "
                "frcmod and lib files to parametrize the ligand."
            )
        if not os.access(path, os.R_OK):
            raise PermissionError(f"Membrane {field} is not readable: {path}")


def _ensure_legacy_lipid_leaprc(config: MembraneConfig, work_dir: str, amberhome) -> Optional[str]:
    """Make a relocated legacy lipid ``leaprc`` available to tleap.

    AmberTools 24 archived older lipid force fields (Lipid17/14/11): their
    ``leaprc`` files moved from ``dat/leap/cmd/`` into ``dat/leap/cmd/oldff/``,
    which tleap does not search by default, so packmol-memgen's
    ``source leaprc.lipid17`` fails deep inside tleap even though the parameter
    (``.dat``) and library (``.lib``) files are still present.

    The fix is a temporary copy in the working directory: tleap searches the
    current directory first, and the archived ``leaprc`` refers to its ``.dat``
    and ``.lib`` files through tleap's search path rather than relatively, so a
    work-dir copy behaves exactly like the canonical one.  Restoring it into
    ``$AMBERHOME/dat/leap/cmd/`` instead would work too, but it permanently
    mutates a shared AmberTools install from inside a build -- nothing else in
    PRISM writes outside the output directory, and a truncated copy would leave
    a broken ``leaprc`` behind for every future tleap run on that machine.

    Returns the staged path, which the caller must remove after the run, or
    ``None`` when nothing was staged.
    """
    lipid_ff = (getattr(config, "lipid_ff", "") or "").strip().lower()
    if not lipid_ff or not amberhome:
        return None
    leaprc = f"leaprc.{lipid_ff}"
    cmd_dir = os.path.join(amberhome, "dat", "leap", "cmd")
    if os.path.isfile(os.path.join(cmd_dir, leaprc)):
        return None  # already on tleap's default search path
    legacy = os.path.join(cmd_dir, "oldff", leaprc)
    if not os.path.isfile(legacy):
        return None

    destination = os.path.join(work_dir, leaprc)
    if not os.path.exists(destination):
        try:
            shutil.copy2(legacy, destination)
            return destination
        except OSError:
            return None
    return None


def _persist_log(work_dir: str, log: str) -> None:
    """Write the captured tool output next to the build artifacts.

    This is the longest single step PRISM runs and its record used to be kept
    only when it failed -- so a successful-but-imperfect pack (see
    :func:`packing_quality_warning`), the tleap parametrization output and the
    resolved command line all vanished on the happy path.
    """
    try:
        with open(os.path.join(work_dir, PRISM_LOG_NAME), "w") as handle:
            handle.write(log)
    except OSError:
        pass  # bookkeeping must never mask the build's own outcome


def run(
    protein_pdb: str,
    config: MembraneConfig,
    work_dir: str,
    timeout: int = 7200,
    verbose: bool = True,
    dry: bool = False,
) -> PackmolMemgenResult:
    """Run packmol-memgen in ``work_dir``. Returns paths to the produced files.

    ``dry=True`` packs protein + lipids only; the returned prmtop/crd contain no
    water and no ions, and the caller MUST solvate and ionise them before the
    system is usable.  The result is only reported successful once the absence
    of solvent has been verified against PACKMOL's own input and the produced
    topology -- a dry flag that was silently ignored would otherwise hand the
    caller a fully solvated system to solvate a second time.

    Raises RuntimeError if packmol-memgen is not installed.
    """
    config.raise_for_errors()
    # The subprocess runs with ``cwd=work_dir``. Passing a caller-relative PDB
    # path unchanged would make it resolve relative to that new cwd (for
    # example ``out/GMX/.../out/GMX/.../protein.pdb``). Normalize both paths
    # before constructing the command so the standalone public API works with
    # relative output directories too.
    protein_pdb = os.path.abspath(protein_pdb)
    work_dir = os.path.abspath(work_dir)
    if not is_available():
        raise RuntimeError(
            "packmol-memgen not found. Install AmberTools (conda install -c conda-forge ambertools) "
            "or run PRISM from an environment that provides it (e.g. the AmberTools conda env)."
        )
    # Check the ligand parameters before anything expensive happens: a missing
    # frcmod/lib only surfaces as an opaque tleap error after packing.
    if config.has_ligand():
        _require_readable_ligand_parameters(config)
    os.makedirs(work_dir, exist_ok=True)

    # packmol-memgen resolves packmol, tleap, reduce and the minimization engine
    # from $AMBERHOME/bin (or its own bundle) and never from PATH, so the
    # environment handed to the subprocess decides what it can actually do.
    env = os.environ.copy()
    amberhome = _resolve_amberhome()
    if amberhome:
        env["AMBERHOME"] = amberhome
    # A packmol that is only on PATH is invisible to packmol-memgen; give it the
    # path explicitly rather than pointing AMBERHOME at packmol's prefix, which
    # would break every other helper it looks up there.
    packmol_override = None if _resolve_packmol(amberhome) else shutil.which("packmol")

    notes = []
    protein_ff, water_ff = config.packmol_forcefield_pair()
    forcefield_options = (None, None)
    if (protein_ff, water_ff) != ("ff14SB", "tip3p"):
        # The default pair emits no force-field flags at all, so there is no
        # spelling to resolve and no reason to pay for the --help probe.
        forcefield_options = _detect_forcefield_options()
    engine = _resolve_minimization_engine(amberhome) if config.minimize else None
    if config.minimize and engine is None:
        notes.append(
            "[PRISM] No AMBER minimization engine found in $AMBERHOME/bin; packmol-memgen will "
            "fall back to its own default (pmemd.cuda), which needs a working GPU."
        )
    solvent_parm = None
    if dry:
        # Refuse before orientation output is touched rather than after: an
        # AmberTools too old for these flags fails at argparse (exit 2) with a
        # message that says nothing about membranes.
        _verify_dry_options()
        solvent_parm = write_dry_solvent_parm(work_dir)
        notes.append(
            f"[PRISM] Lipids-only (dry) pack: PACKMOL places no water and no ions. "
            f"Solvent table {os.path.basename(solvent_parm)} sets the water density to "
            f"{DRY_SOLVENT_DENSITY_G_CM3:g} g/cm3 so the count truncates to zero. "
            "The caller must run gmx solvate/genion on the result."
        )
    # Disulfides.  packmol-memgen cannot infer these: its tleap generator emits
    # ``bond`` lines only for cardiolipin, and its protonation step renames an SS
    # cysteine to CYX (dropping HG) without ever creating the S-S bond, leaving a
    # dangling sulfur carrying a charge that assumes a bonded partner.  Until
    # this call existed the route asked for ``config.disulfide_pairs``, which was
    # never a field on MembraneConfig, so it declared none at all -- on a
    # structure like mGluR7 that is 19 missing bonds, free cysteines with thiol
    # hydrogens driven into occupied space, and a preflight energy of +1.8e16
    # kJ/mol.
    #
    # Detection runs on the very file packmol-memgen is about to read.  That is
    # what makes the ordinals valid: packmol-memgen copies the solute into the
    # packed system first and verbatim (measured -- the protein block of a packed
    # system is residue-for-residue identical to its input), so a residue's
    # position among the solute's residues is the same in both files.
    #
    # A protein PDB that is not there is not a build whose chemistry can come out
    # silently wrong -- packmol-memgen reports the missing input itself, and it
    # names it better than a disulfide detector would.  This is the one case
    # where skipping detection is not the silent omission this call exists to
    # fix, so it is spelled out rather than left to a bare except.
    disulfide_pairs, disulfides_described, unresolved = (
        detect_disulfide_pairs(protein_pdb) if os.path.isfile(protein_pdb) else ([], [], [])
    )
    if unresolved:
        raise RuntimeError(
            "Detected disulfide(s) in "
            f"{os.path.basename(protein_pdb)} whose residues could not be matched to tleap's "
            f"own residue numbering: {'; '.join(unresolved)}. Declaring them against a guessed "
            "number would bond the wrong residues silently, and not declaring them leaves free "
            "cysteines with overlapping thiol hydrogens. Renumber the structure so each chain "
            "numbers its residues uniquely (no insertion codes, no restarts within a chain)."
        )
    if disulfides_described:
        notes.append(
            f"[PRISM] Declaring {len(disulfides_described)} disulfide bond(s) to tleap: "
            + "; ".join(disulfides_described)
        )

    cmd = build_command(
        protein_pdb,
        config,
        parametrize=True,
        protein_ff_option=forcefield_options[0],
        water_ff_option=forcefield_options[1],
        engine=engine,
        packmol_path=packmol_override,
        dry=dry,
        solvent_parm=solvent_parm,
        disulfide_pairs=disulfide_pairs,
    )
    _verify_ligand_options(cmd)
    # Snapshot AFTER the solvent table is written so it is never mistaken for an
    # output of this run by the stale-artifact filters below.
    before_run = _file_signatures(work_dir)

    # AmberTools 24 relocated legacy lipid leaprc files (e.g. leaprc.lipid17)
    # into dat/leap/cmd/oldff/, which tleap does not search by default, so
    # packmol-memgen's `source leaprc.lipid17` fails deep inside tleap. Stage a
    # temporary copy in the working directory (tleap searches it first) that is
    # removed again afterwards.
    staged_leaprc = _ensure_legacy_lipid_leaprc(config, work_dir, amberhome)
    if staged_leaprc:
        notes.append(
            f"[PRISM] Staged {os.path.basename(staged_leaprc)} from AmberTools' oldff/ directory "
            "for the duration of this run; it is not on tleap's default search path."
        )

    if verbose:
        # Packing is the longest step of the whole build and used to run in
        # complete silence for up to packmol_timeout seconds.
        print(f"[membrane] packmol-memgen in {work_dir} (timeout {timeout:g}s; typically 10-60 min)")
        print("  " + " ".join(cmd))
        for note in notes:
            print(f"  {note}")

    # Guard the subprocess: a real pack+parametrize+minimize can exceed the
    # timeout, and the OS may raise (e.g. ENOMEM). Return a failed result rather
    # than propagating, so the builder degrades gracefully.
    try:
        proc = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, timeout=timeout, env=env)
    except subprocess.TimeoutExpired as exc:
        log = "\n".join(notes + [
            f"packmol-memgen timed out after {timeout}s: {exc}",
            "[PRISM] Packing cost grows with box volume and with how awkward the solute's "
            "surface is; a large complex needs longer. Raise the budget with the membrane "
            "YAML key 'packmol_timeout' (seconds). Building on a native Linux filesystem "
            "rather than a mounted Windows drive also helps substantially, because packmol "
            "rewrites the whole packed system on every iteration.",
        ])
        _persist_log(work_dir, log)
        return PackmolMemgenResult(
            prmtop=None, inpcrd=None, packed_pdb=None, log=log, success=False, dry=dry
        )
    except OSError as exc:
        log = "\n".join(notes + [f"packmol-memgen failed to launch: {exc}"])
        _persist_log(work_dir, log)
        return PackmolMemgenResult(
            prmtop=None, inpcrd=None, packed_pdb=None, log=log, success=False, dry=dry
        )
    finally:
        # The staged leaprc is only needed while tleap runs inside packmol-memgen;
        # remove it so the build directory keeps only real build artifacts.
        if staged_leaprc and os.path.isfile(staged_leaprc):
            try:
                os.remove(staged_leaprc)
            except OSError:
                pass
    log = "\n".join(notes + ["$ " + " ".join(cmd), (proc.stdout or "") + "\n" + (proc.stderr or "")])

    # Select a topology/coordinate pair with a common basename.  Independent
    # "first file" searches can combine leftovers from two previous builds (or
    # even PRISM's own topol.top) into an invalid pair.
    prmtop, inpcrd = _find_amber_pair(
        work_dir,
        prefer_minimized=config.minimize,
        changed_since=before_run,
    )
    packed = _find_first(
        work_dir,
        ("bilayer_",),
        prefix_match=True,
        suffix=".pdb",
        changed_since=before_run,
    )
    log += _packing_verdict_note(work_dir, packed, changed_since=before_run)

    success = proc.returncode == 0 and prmtop is not None and inpcrd is not None
    if proc.returncode != 0:
        # A failed process may leave a half-written pair. Do not expose those
        # files as usable AMBER inputs merely because their names match.
        prmtop = None
        inpcrd = None
    elif _packing_was_skipped(log):
        # PACKMOL-Memgen reuses an existing bilayer of the expected name instead
        # of packing ("Packed PDB ... found. Skipping PACKMOL") and still exits
        # 0, so it would happily parametrize the PREVIOUS run's system -- built
        # from a different protein, ligand or lipid composition -- and report
        # success.  The filename is derived from the input, so a re-run into a
        # dirty directory hits this every time.
        success = False
        log += (
            "\n[PRISM] PACKMOL-Memgen skipped packing and reused an existing bilayer in this "
            "directory, which would parametrize a previous build's system. Remove the output "
            "directory and re-run instead of reusing it.\n"
        )
    elif config.has_ligand() and packed:
        # Cheapest possible guard against the worst silent failure: a ligand
        # dropped during protonation/cleaning yields a perfectly valid bilayer
        # that simply has no ligand in it.
        missing = _ligand_absent_from(packed, config.ligand_resname)
        if missing:
            success = False
            log += (
                f"\n[PRISM] Residue {config.ligand_resname!r} does not appear in the packed system "
                f"{os.path.basename(packed)}; the ligand was dropped during packing.\n"
            )

    # Guard G1: the pack really is dry.  This runs last and only on the dry
    # path, because the whole fast path is built on the premise that the system
    # arrives solvent-free.  If a flag were silently ignored -- or a stale wet
    # topology were reused -- the downstream stage would solvate an already
    # solvated box and ionise it a second time, which no later check catches.
    if dry and success:
        violations = dry_pack_violations(work_dir, prmtop, changed_since=before_run)
        if violations:
            success = False
            prmtop = None
            inpcrd = None
            log += (
                "\n[PRISM] The lipids-only (dry) pack was requested but the result is not "
                "solvent-free:\n  - " + "\n  - ".join(violations) + "\n"
                "Refusing to hand a pre-solvated system to the GROMACS solvation stage. "
                "Build with --membrane-solvate packmol instead.\n"
            )

    _persist_log(work_dir, log)
    return PackmolMemgenResult(
        prmtop=prmtop, inpcrd=inpcrd, packed_pdb=packed, log=log, success=success, dry=dry
    )


def _packing_verdict_note(work_dir: str, packed_pdb, changed_since=None) -> str:
    """Return PACKMOL's own packing verdict as a log line, or ``""``.

    ``packing_quality_warning`` cannot see this in the captured output:
    packmol-memgen redirects the PACKMOL binary's stdout into ``packmol.log``
    and re-emits only a few progress strings, so "ENDED WITHOUT PERFECT
    PACKING" -- the module's single guard against a clash-ridden starting
    structure -- never reaches PRISM. Read it back from the file instead, and
    treat the ``_FORCED`` structure PACKMOL writes alongside it as the same
    signal.
    """
    verdict = "ENDED WITHOUT PERFECT PACKING"
    name = f"{PACKMOL_LOG_PREFIX}.log"
    signatures = _file_signatures(work_dir)
    # A log left by a previous build in the same directory says nothing about
    # this one.
    if name in signatures and (changed_since is None or signatures[name] != changed_since.get(name)):
        try:
            with open(os.path.join(work_dir, name), "r", errors="replace") as handle:
                if any(verdict in line.upper() for line in handle):
                    return f"\n[PRISM] PACKMOL reported: {verdict} (see {name}).\n"
        except OSError:
            pass
    if packed_pdb and os.path.exists(packed_pdb + "_FORCED"):
        return (
            f"\n[PRISM] PACKMOL wrote {os.path.basename(packed_pdb)}_FORCED, which means it "
            f"{verdict}.\n"
        )
    return ""


def packing_quality_warning(log: str) -> Optional[str]:
    """Describe an imperfect PACKMOL result, or None when packing converged.

    PACKMOL reports ``ENDED WITHOUT PERFECT PACKING`` and writes a ``_FORCED``
    structure when it exhausts its GENCAN budget with overlaps still present.
    It still exits 0, so nothing downstream notices. The residual clashes are
    normally absorbed by the minimization and the restrained equilibration, but
    the user has to be told -- a badly packed start is the usual cause of a
    first-step LINCS failure.
    """
    if "ENDED WITHOUT PERFECT PACKING" not in (log or "").upper():
        return None
    return (
        "PACKMOL exhausted its GENCAN budget with residual overlaps and wrote a forced "
        "structure. The minimization and restrained equilibration usually absorb this, but "
        "watch the first equilibration stage for LINCS warnings. To pack harder raise "
        "'packing_nloop_all'; to stop sooner (it often fails to converge regardless) lower it."
    )


#: PACKMOL-Memgen writes its generated PACKMOL script to ``<packlog>.inp``.
PACKMOL_INPUT_NAME = f"{PACKMOL_LOG_PREFIX}.inp"


def _packmol_input_structure_counts(path: str) -> Optional[Dict[str, int]]:
    """Return ``{structure filename: molecules requested}`` from a PACKMOL input.

    PACKMOL's input is a flat sequence of ``structure <file> ... end structure``
    blocks, each carrying one ``number N``.  Reading it back is the only way to
    see what PACKMOL was actually *asked* to place: packmol-memgen logs a
    summary, not the request, and a molecule count of zero is expressed by the
    block being absent rather than by ``number 0``.

    ``None`` means the file could not be read at all, which is different from
    "read it and found nothing" and must not be treated as a clean result.
    """
    counts: Dict[str, int] = {}
    current: Optional[str] = None
    try:
        with open(path, "r", errors="replace") as handle:
            for raw in handle:
                line = raw.strip()
                if line.startswith("structure "):
                    current = os.path.basename(line.split(None, 1)[1].strip())
                    counts.setdefault(current, 0)
                elif line.startswith("end structure"):
                    current = None
                elif current is not None and line.startswith("number "):
                    try:
                        counts[current] += int(float(line.split()[1]))
                    except (IndexError, ValueError):
                        # A malformed count is not evidence of absence; leave the
                        # block at its running total and let the topology check
                        # be the authority.
                        continue
    except OSError:
        return None
    return counts


def _is_bulk_solvent_structure(structure: str) -> bool:
    """True when a PACKMOL ``structure`` filename is water or an ion.

    PACKMOL-Memgen names these files after the residue: ``WAT.pdb`` for the
    solvent (PRISM never overrides ``--solvents``) and ``Na+.pdb`` / ``Cl-.pdb``
    / ``K+.pdb`` for the ions PRISM can request.  Lipid and protein blocks are
    named after their own residues and cannot collide with these.
    """
    name = os.path.splitext(str(structure).strip())[0]
    return name.upper() in _BULK_WATER_RESIDUES or name in _ION_RESIDUES


def _is_bulk_solvent_residue(label: str) -> bool:
    """True when an AMBER residue label is water or a free ion."""
    return str(label).upper() in _BULK_WATER_RESIDUES or str(label) in _ION_RESIDUES


def _amber_residue_label_counts(path: str) -> Optional[Dict[str, int]]:
    """Return a residue-label histogram from an AMBER topology, or ``None``.

    Parsed directly rather than through ParmEd: this module is deliberately
    dependency-free (it already reads POINTERS and the NetCDF header by hand),
    and loading a 350 000-atom prmtop through ParmEd to count residue names
    costs seconds and a large allocation for an answer that is one fixed-width
    section away.

    ``None`` (unreadable / no ``RESIDUE_LABEL`` section) is distinct from an
    empty histogram, so callers can fail closed instead of concluding "no
    water found".
    """
    counts: Dict[str, int] = {}
    try:
        with open(path, "r", encoding="ascii", errors="replace") as handle:
            for line in handle:
                if line.strip() != "%FLAG RESIDUE_LABEL":
                    continue
                for label_line in handle:
                    if label_line.startswith(("%FORMAT", "%COMMENT")):
                        continue
                    if label_line.startswith("%"):
                        return counts  # next section: the labels are complete
                    text = label_line.rstrip("\n")
                    # AMBER writes 20 labels per line in fixed 4-character
                    # fields; splitting on whitespace would merge an empty
                    # field into its neighbour and shift every later label.
                    for start in range(0, len(text), 4):
                        label = text[start:start + 4].strip()
                        if label:
                            counts[label] = counts.get(label, 0) + 1
                return counts
    except OSError:
        return None
    return None


def dry_pack_violations(work_dir: str, prmtop, changed_since=None) -> List[str]:
    """Return reasons a supposedly dry pack is not actually solvent-free.

    Two independent detectors, because they fail in different ways:

    * PACKMOL's own input says what was *requested*, and catches a
      ``--solvent_parm`` that never took effect;
    * the produced AMBER topology says what was *delivered*, and additionally
      catches solvent introduced by tleap (``addionsrand``) or a stale wet
      topology that packmol-memgen reused instead of regenerating.

    An empty list means verified dry.  If neither detector could be evaluated
    the result is a violation of its own: "unverifiable" must not read as
    "clean", because the whole downstream fast path assumes this guarantee.
    """
    problems: List[str] = []
    verified = False

    signatures = _file_signatures(work_dir)
    input_path = os.path.join(work_dir, PACKMOL_INPUT_NAME)
    # An input file left by a previous build says nothing about this one.
    if PACKMOL_INPUT_NAME in signatures and (
        changed_since is None
        or signatures[PACKMOL_INPUT_NAME] != changed_since.get(PACKMOL_INPUT_NAME)
    ):
        requested = _packmol_input_structure_counts(input_path)
        if requested is not None:
            verified = True
            for structure, number in sorted(requested.items()):
                if number > 0 and _is_bulk_solvent_structure(structure):
                    problems.append(
                        f"{PACKMOL_INPUT_NAME} asks PACKMOL to place {number} x {structure}"
                    )

    if prmtop:
        residues = _amber_residue_label_counts(prmtop)
        if residues is not None:
            verified = True
            for label, number in sorted(residues.items()):
                if _is_bulk_solvent_residue(label):
                    problems.append(
                        f"{os.path.basename(prmtop)} contains {number} {label!r} residues"
                    )

    if not verified:
        problems.append(
            f"neither {PACKMOL_INPUT_NAME} nor the AMBER topology could be read back, so the "
            "absence of water and ions could not be confirmed"
        )
    return problems


#: Atomic masses (Da) for the elements a membrane input can contain.  Only used
#: by :func:`estimate_water_count`, which needs the solute mass to reproduce
#: PACKMOL-Memgen's volume model; an unrecognised element contributes nothing,
#: exactly as in PACKMOL-Memgen's own accounting.
_ATOMIC_MASSES = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "NA": 22.990, "MG": 24.305, "P": 30.974, "S": 32.060, "CL": 35.450,
    "K": 39.098, "CA": 40.078, "MN": 54.938, "FE": 55.845, "CO": 58.933,
    "NI": 58.693, "CU": 63.546, "ZN": 65.380, "SE": 78.971, "BR": 79.904,
    "I": 126.904,
}


def _pdb_element(line: str) -> str:
    """Return the element symbol of a PDB coordinate record.

    Columns 77-78 are authoritative when present; otherwise fall back to the
    leading alphabetic characters of the atom name, which is what every PDB
    reader does for the pre-element-column files still in circulation.
    """
    element = line[76:78].strip()
    if element:
        return element.upper()
    name = line[12:16].strip()
    # A four-character atom name is left-justified with the element in column
    # 13; a shorter one is right-justified with the element in column 14.  In
    # both cases the element is the leading alphabetic run of the stripped name
    # once any leading digit (1HB, 2HB) is dropped.
    stripped = name.lstrip("0123456789")
    return "".join(character for character in stripped[:2] if character.isalpha()).upper()


def estimate_water_count(protein_pdb: str, config: MembraneConfig) -> int:
    """Estimate how many waters a *full* PACKMOL-Memgen pack would place.

    This exists so a caller can pick between the packing modes before paying
    for either: PACKMOL's cost is dominated by the water it is asked to place,
    and the number is not knowable from the config alone -- it is set by the
    solute's extramembrane extent, which for a class-C GPCR is a 18 nm
    Venus-flytrap span and for a small channel is nothing at all.

    The arithmetic deliberately mirrors PACKMOL-Memgen's own
    (``main.py``/``lib/utils.py``) rather than inventing a simpler model:

    * the box side is ``2 * max_xy_radius + 2 * --dist``;
    * z spans at least ``+/- (leaflet + --dist_wat)`` and otherwise the solute
      plus ``--dist_wat`` on each face;
    * the solute's own volume is removed from each water slab using the
      empirical protein density ``1.41 + 0.145 * exp(-MW/13000)`` g/cm^3
      (Protein Sci. 2004, 13:2825);
    * the count is the slab volume times water's number density.

    It is an estimate, not a prediction: PACKMOL-Memgen refines the solute
    volume with a Monte-Carlo profile, so expect agreement to ~10%, which is
    ample for a decision against a 30 000-molecule threshold.  Measured against
    the structure PACKMOL-Memgen itself packed (its ``PROT0.pdb``) this is
    within 1% on both reference systems: 10 013 vs 10 038 on the T4-lysozyme
    build and 130 333 vs 131 635 on the 9OMO GPCR dimer.

    One caveat, and it is a real one: with ``orient="memembed"``
    PACKMOL-Memgen re-orients the solute *itself* before measuring, so the file
    passed here is not the frame it will measure and the estimate becomes an
    upper bound (+37% on 9OMO).  That errs towards reporting more water than
    there will be, which is the harmless direction for a mode choice.

    Raises ``ValueError`` if the PDB has no readable coordinates -- callers
    selecting a mode automatically should treat that as "unknown" and keep the
    conservative default rather than guessing.
    """
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    masses: List[float] = []
    try:
        with open(protein_pdb, "r", errors="replace") as handle:
            for line in handle:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    continue  # a malformed record is not a reason to abort
                masses.append(_ATOMIC_MASSES.get(_pdb_element(line), 0.0))
    except OSError as exc:
        raise ValueError(f"Cannot estimate the water count from {protein_pdb}: {exc}") from exc
    if not xs:
        raise ValueError(f"No ATOM/HETATM records in {protein_pdb}; cannot estimate a water count.")

    # PACKMOL-Memgen centres the solute in xy before measuring, so the radius is
    # taken about the centroid rather than about the origin.
    center_x = sum(xs) / len(xs)
    center_y = sum(ys) / len(ys)
    max_radius = max(math.hypot(x - center_x, y - center_y) for x, y in zip(xs, ys))
    side = 2 * max_radius + 2 * float(config.xy_distance)

    leaflet = PACKMOL_LEAFLET_HALF_THICKNESS_A
    water_thickness = float(config.water_thickness)
    z_min = min(min(zs) - water_thickness, -(leaflet + water_thickness))
    z_max = max(max(zs) + water_thickness, leaflet + water_thickness)

    total_mass = sum(masses)
    if total_mass <= 0:
        raise ValueError(
            f"No recognisable elements in {protein_pdb}; cannot estimate a water count."
        )
    # g/cm^3 -> Da/A^3, which is the unit PACKMOL-Memgen divides masses by.
    density = (1.41 + 0.145 * math.exp(-total_mass / 13000.0)) * AVOGADRO / 1e24
    solute_upper = sum(m for m, z in zip(masses, zs) if z > leaflet) / density
    solute_lower = sum(m for m, z in zip(masses, zs) if z < -leaflet) / density

    upper_volume = max(side * side * (z_max - leaflet) - solute_upper, 0.0)
    lower_volume = max(side * side * (abs(z_min) - leaflet) - solute_lower, 0.0)
    number_density = WATER_DENSITY_G_CM3 * AVOGADRO / (DRY_SOLVENT_MOLAR_MASS * 1e24)
    # Truncated per slab, matching PACKMOL-Memgen's two independent int() casts.
    return int(upper_volume * number_density) + int(lower_volume * number_density)


def _packing_was_skipped(log: str) -> bool:
    """True when PACKMOL-Memgen reused an existing bilayer instead of packing."""
    lowered = (log or "").lower()
    return "skipping packmol" in lowered


def _ligand_absent_from(packed_pdb: str, resname: str) -> bool:
    """True when ``resname`` has no coordinate record in the packed system."""
    token = str(resname).strip().upper()
    if not token:
        return False
    try:
        with open(packed_pdb, "r", errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip().upper() == token:
                    return False
    except OSError:
        return False  # unreadable output is reported by the caller's own checks
    return True


def _find_first(
    directory: str,
    patterns,
    prefix_match: bool = False,
    suffix: str = "",
    changed_since=None,
) -> Optional[str]:
    """Find the first file in ``directory`` matching the given pattern(s)."""
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return None
    current_signatures = _file_signatures(directory) if changed_since is not None else None
    for name in names:
        if changed_since is not None and current_signatures.get(name) == changed_since.get(name):
            continue
        if prefix_match:
            if any(name.startswith(p) for p in patterns) and name.endswith(suffix):
                return os.path.join(directory, name)
        else:
            if any(name.endswith(p) for p in patterns):
                return os.path.join(directory, name)
    return None


def _file_signatures(directory: str):
    """Return lightweight signatures used to distinguish outputs from stale files."""
    try:
        names = os.listdir(directory)
    except OSError:
        return {}

    signatures = {}
    for name in names:
        path = os.path.join(directory, name)
        try:
            stat = os.stat(path)
        except OSError:
            continue
        if os.path.isfile(path):
            signatures[name] = (stat.st_mtime_ns, stat.st_size)
    return signatures


def _amber_topology_atom_count(path: str) -> Optional[int]:
    """Read NATOM from the POINTERS section of an AMBER topology."""
    try:
        with open(path, "r", encoding="ascii", errors="replace") as handle:
            for line in handle:
                if line.strip() != "%FLAG POINTERS":
                    continue
                for pointer_line in handle:
                    stripped = pointer_line.strip()
                    if not stripped or stripped.startswith("%FORMAT"):
                        continue
                    return int(stripped.split()[0])
    except (OSError, ValueError):
        pass
    return None


def _amber_coordinate_atom_count(path: str) -> Optional[int]:
    """Read NATOM from a text or classic-NetCDF AMBER coordinate file."""
    try:
        with open(path, "rb") as handle:
            magic = handle.read(4)
            if magic[:3] == b"CDF" and magic[3:4] in {b"\x01", b"\x02"}:
                # Classic NetCDF (CDF-1/CDF-2) begins with numrecs followed by
                # the dimension array. AMBER restart files expose NATOM as the
                # dimension named ``atom``. Only this small, stable header
                # fragment is needed; no optional NetCDF dependency is required.
                handle.read(4)
                tag, dimension_count = struct.unpack(">II", handle.read(8))
                if tag != 10:  # NC_DIMENSION
                    return None
                for _ in range(dimension_count):
                    name_length = struct.unpack(">I", handle.read(4))[0]
                    name = handle.read(name_length)
                    handle.read((-name_length) % 4)
                    dimension_length = struct.unpack(">I", handle.read(4))[0]
                    if name == b"atom":
                        return dimension_length
                return None

        with open(path, "r", encoding="ascii", errors="strict") as handle:
            next(handle)
            return int(next(handle).split()[0])
    except (OSError, StopIteration, UnicodeError, ValueError, IndexError, struct.error):
        return None


def _find_amber_pair(directory: str, prefer_minimized: bool = False, changed_since=None):
    """Return a matching AMBER topology/coordinate pair from ``directory``.

    When ``changed_since`` contains a pre-run signature snapshot, both files
    must have been created or updated by the current PACKMOL-Memgen invocation.
    This prevents a successful return code from reviving a stale pair left by
    an earlier build.
    """
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return None, None

    current_signatures = _file_signatures(directory)

    def changed(name):
        return name in current_signatures and (
            changed_since is None or current_signatures[name] != changed_since.get(name)
        )

    topology_names = [
        name
        for name in names
        if name != "topol.top" and name.endswith((".prmtop", ".top")) and changed(name)
    ]
    coordinate_suffixes = (".inpcrd", ".crd", ".rst7", ".restrt", ".ncrst")
    coordinate_names = [
        name for name in names if name.endswith(coordinate_suffixes) and changed(name)
    ]
    valid_topologies = []
    for topology_name in topology_names:
        atom_count = _amber_topology_atom_count(os.path.join(directory, topology_name))
        if atom_count is not None:
            valid_topologies.append((topology_name, atom_count))

    pairs = []
    for topology_name, topology_atoms in valid_topologies:
        topology_suffix = ".prmtop" if topology_name.endswith(".prmtop") else ".top"
        stem = topology_name[: -len(topology_suffix)]
        minimized_stem = stem[:-len("_lipid")] if stem.endswith("_lipid") else stem
        final_minimized_names = {
            f"{minimized_stem}_min.restrt",
            f"{minimized_stem}_min.ncrst",
            f"{minimized_stem}_min.rst7",
        }
        for coordinate_name in coordinate_names:
            direct_match = coordinate_name.startswith(stem + ".") or coordinate_name.startswith(stem + "_")
            final_minimized = coordinate_name in final_minimized_names
            direct_minimized = direct_match and "min" in coordinate_name[len(stem):].lower()
            if prefer_minimized:
                # Do not silently discard a requested minimization. In current
                # PACKMOL-Memgen, generic min.restrt is the restrained stage;
                # the final unrestrained result is <base>_min.restrt.
                if not (final_minimized or direct_minimized):
                    continue
            elif not direct_match or direct_minimized:
                continue

            coordinate_atoms = _amber_coordinate_atom_count(os.path.join(directory, coordinate_name))
            if coordinate_atoms != topology_atoms:
                continue

            if prefer_minimized:
                if final_minimized:
                    rank = 0  # final, unrestrained PACKMOL-Memgen minimization
                elif direct_minimized:
                    rank = 1
            else:
                rank = 0
            pairs.append((rank, topology_name, coordinate_name))

    if pairs:
        _, topology_name, coordinate_name = min(pairs)
        return os.path.join(directory, topology_name), os.path.join(directory, coordinate_name)

    return None, None
