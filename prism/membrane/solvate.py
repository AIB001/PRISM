#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GROMACS-side solvation and ionisation for the membrane fast path.

The default (slow) membrane path lets PACKMOL place every water and every ion
individually through GENCAN.  On a class-C GPCR that is ~130,000 molecules and
the pack never finishes.  ``gmx solvate`` does the same job in seconds by
stamping a pre-equilibrated water box on a grid, so the fast path packs only
the protein and the lipids and hands the aqueous phase to GROMACS::

    orient -> packmol-memgen (dry) -> ParmEd -> box normalise
           -> [ this module: solvate + purge + genion ]
           -> index -> restraints -> MDPs -> runscript

The reason this module exists rather than a two-line call to
:class:`~prism.utils.system.solvation.SolvationProcessorMixin` is that
``gmx solvate`` is *wrong* on a bilayer by default.  It inserts water wherever
vdW radii allow, and a freshly packed bilayer has metre-wide gaps between the
acyl tails; the result is water in the hydrophobic core, which is a
catastrophic and well known membrane-building failure.  Two things are needed:

1. a ``vdwradii.dat`` with inflated carbon radii (the canonical Lemkul recipe,
   C = 0.375 nm), which *reduces* core insertion but never eliminates it -- it
   is asymptotic, and even absurd radii leave tens of core waters while
   progressively under-filling the bulk; and
2. a deterministic **geometric purge** of everything that still lands between
   the two phosphate planes, which is the only thing that reaches zero.

Both are implemented here, together with the correctness guards that make the
result trustworthy: a fast wrong answer is worthless.  Every guard has a cheap
detector and a hard failure that names ``--membrane-solvate packmol`` as the
escape hatch.  There is deliberately **no silent fallback** to the slow path:
an automatic retry would burn hours without the user knowing the science
changed.

Ion placement does *not* reuse ``SolvationProcessorMixin._add_ions``.  That
method is correct for a soluble box and membrane-fatal here: ``genion -conc``
derives the ion count from the *total* box volume, but only ~40% of a membrane
box is aqueous, so a requested 0.15 M silently becomes ~0.27 M.  It also
rewrites the topology before placement and leaves it mutated after a fatal
error, and its post-hoc verification matches the substring ``'CA '`` against
alpha-carbon atom names, so it reports success unconditionally.  Rather than
"fix" the validated soluble path, this module computes ``-np``/``-nn``
explicitly from the aqueous water count and the topology's own net charge.
"""

from __future__ import annotations

import math
import os
import re
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

# Imported as a module, not by name: config owns the tunables this module reads
# back (the carbon radius and the bounds it is validated against), and going
# through the module keeps the single source of truth visible at each use site.
from . import config as _config

try:  # numpy drives the neighbour search; keep the import failure actionable
    import numpy as _np
except ImportError:  # pragma: no cover - numpy is a hard PRISM dependency
    _np = None

# The shared residue-name vocabulary.  ``config`` has no intra-package imports,
# so this cannot introduce a cycle no matter who imports this module first.
try:
    from . import config as membrane_config
except ImportError:  # pragma: no cover - only when imported outside the package
    from prism.membrane import config as membrane_config

try:
    from ..utils.system.base import SystemBuilderBase
    from ..utils.system.solvation import SolvationProcessorMixin
    from ..utils.system.topology import TopologyProcessorMixin
except ImportError:  # pragma: no cover - only when imported outside the package
    from prism.utils.system.base import SystemBuilderBase
    from prism.utils.system.solvation import SolvationProcessorMixin
    from prism.utils.system.topology import TopologyProcessorMixin

# Console output mirrors prism/membrane/builder.py: a missing colour module
# must never make this module unimportable.
try:
    from ..utils.colors import print_success, print_warning, print_info, number
except ImportError:  # pragma: no cover - exercised only on a broken install
    try:
        from prism.utils.colors import print_success, print_warning, print_info, number
    except ImportError:
        def print_success(x, **kwargs):
            print(f"✓ {x}")

        def print_warning(x, **kwargs):
            print(f"⚠ {x}")

        def print_info(x, **kwargs):
            print(f"ℹ {x}")

        def number(x):
            return x


__all__ = [
    "MembraneSolvationError",
    "LeafletPlanes",
    "SolvationOutcome",
    "MembraneSolvator",
    "leaflet_planes_from_parm",
    "leaflet_planes_from_gro",
    "read_topology",
    "MoleculeType",
    "TopologyView",
    "DEFAULT_CARBON_RADIUS_NM",
]


# --------------------------------------------------------------------------- #
# Constants.  Every magic number that a reviewer would otherwise have to guess
# at lives here with the measurement or convention that produced it.
# --------------------------------------------------------------------------- #

#: Residue/moleculetype name ``gmx solvate`` stamps down for water.  The name
#: is fixed by the ``-cs`` coordinate files GROMACS ships (spc216/tip4p/tip5p),
#: so the topology handed to :meth:`MembraneSolvator.run` *must* publish its
#: water under exactly this moleculetype name or nothing downstream lines up.
WATER_MOLECULETYPE = "SOL"

#: Bulk water number density at 300 K / 1 bar (molecules per nm^3).  Used only
#: to sanity-check that ``gmx solvate`` filled the aqueous slab (guard G4b).
WATER_NUMBER_DENSITY_NM3 = 33.4

#: Molarity of pure water (mol/L).  ``n_salt = conc * n_water / 55.56`` is the
#: same conversion ``gmx genion -conc`` uses internally, except genion derives
#: the water count from the *total* box volume, which is what makes it wrong
#: on a membrane.
WATER_MOLARITY = 55.56

#: How far a SUMMED system charge may sit from a whole number before PRISM
#: treats it as a broken topology rather than as rounding noise.
#:
#: Deliberately looser than the 1e-3 used for a single ion's charge below, and
#: the two must not be merged: one ion's charge is read from one line and is
#: exact to its printed precision, while this compares a sum over ~10^5 per-atom
#: charges.  Measured on a 6KQI membrane system, a genuinely +10 solute sums to
#: +10.002.  The failure this guards against -- a molecule counted the wrong
#: number of times, a truncated charge column -- moves the total by 0.5 e or
#: more, so 0.05 separates the two cleanly with two orders of magnitude of
#: headroom on each side.
CHARGE_ROUNDING_TOLERANCE = 0.05

#: Carbon van der Waals radius (nm) written into the membrane ``vdwradii.dat``.
#: 0.375 nm is the value the GROMACS membrane tutorials use; the stock value is
#: 0.17 nm.  Inflating it is necessary but *not* sufficient -- see the purge.
#:
#: Aliased from :mod:`prism.membrane.config` rather than restated: this is a
#: physical constant that config also validates against (its warning and fatal
#: bounds are expressed relative to it), and two independent 0.375 literals
#: would let one be tuned while the other silently kept the old value -- a
#: divergence that changes how much water ends up in the hydrophobic core and
#: that nothing downstream would report.
DEFAULT_CARBON_RADIUS_NM = _config.DEFAULT_VDWRADII_CARBON_NM

#: How far inside each phosphate plane the purge starts.  The headgroup region
#: legitimately holds water; only the acyl slab must be dry.  Trimming 0.3 nm
#: off each side keeps the interfacial hydration shell intact.
CORE_MARGIN_NM = 0.30

#: A water inside the acyl slab but within this distance of a solute heavy atom
#: is *kept*: it is first-shell hydration of a transmembrane pore or of a lipid
#: -facing cavity, not bulk water in the hydrophobic core.  0.45 nm is roughly
#: one water diameter plus a carbon radius.
PORE_RADIUS_NM = 0.45

#: The purge is a scalpel, not a sledgehammer.  Measured cost on real packs is
#: 8-9% of the inserted waters; anything approaching 20% means the leaflet
#: planes were detected wrongly and the run must stop rather than silently
#: dehydrate the box.
MAX_PURGE_FRACTION = 0.20

#: Minimum phosphorus atoms per leaflet for the plane fit to mean anything.
MIN_LEAFLET_P_ATOMS = 10

#: A real bilayer's phosphate planes are ~3.8 nm apart.  Below this the "two
#: leaflets" are an artefact of splitting one distribution in half.
MIN_LEAFLET_SEPARATION_NM = 2.0

#: Fractional tolerance on the bulk-water sanity check (guard G4b).  Wide,
#: because the estimate of the aqueous volume is coarse; it exists to catch a
#: vacuum box (~40% of the expected density), not to audit the density.
BULK_WATER_TOLERANCE = 0.25

#: Fractional tolerance on the delivered salt concentration (guard G6a).
SALT_CONCENTRATION_TOLERANCE = 0.10

#: Volume (nm^3) a single solute heavy atom removes from the aqueous phase.
#: Derived from partial specific volumes: lysozyme is 0.0173 nm^3 per heavy
#: atom, POPC ~0.024; the mean is close enough for a +/-25% band.
HEAVY_ATOM_VOLUME_NM3 = 0.020

#: Residual |net charge| (e) tolerated after ion placement.  Ion charges are
#: integral, so anything above this means the solute charge itself is
#: non-integral (a badly derived ligand charge set) and grompp will complain.
NET_CHARGE_TOLERANCE = 0.05

#: Names ``gmx solvate``/ParmEd/tleap use for water and monatomic ions.  Used
#: only by guard G1, to refuse to solvate a system that is already wet.
#:
#: This module's question is "is this a molecule of the aqueous phase, which I
#: may refuse or delete?", so it is water plus the *bath* ions and nothing else.
#: Structural metals are deliberately excluded by the shared vocabulary: a
#: buried catalytic iron is part of the solute, and a guard that counted it as
#: solvent would refuse to solvate a perfectly dry system.
_SOLVENT_MOLECULE_NAMES = membrane_config.SOLVENT_RESIDUE_NAMES

#: Water-model -> ``-cs`` coordinate file.  This duplicates the table inside
#: ``SolvationProcessorMixin._solvate``, deliberately: the site-count guard has
#: to know which pre-equilibrated box that method will stamp down, and the
#: table lives in a module the fast path is not allowed to modify.  Keep the
#: two in step -- a drift here shows up as a hard G3c failure, not as silence.
_WATER_MODEL_COORD_FILE = {
    "tip3p": "spc216.gro",
    "tip3p_original": "spc216.gro",
    "spc": "spc216.gro",
    "spce": "spc216.gro",
    "opc3": "spc216.gro",
    "tip4p": "tip4p.gro",
    "tip4pew": "tip4p.gro",
    "opc": "tip4p.gro",
    "tip5p": "tip5p.gro",
}
_COORD_FILE_SITES = {"spc216.gro": 3, "tip4p.gro": 4, "tip5p.gro": 5}

#: Sentinel distinguishing "caller did not pass a timeout" from "caller asked
#: for no timeout at all" (which ``SystemBuilderBase._run_command`` spells as
#: ``timeout=None``).
_TIMEOUT_DEFAULT = object()

_PACKMOL_HINT = (
    "Re-run with --membrane-solvate packmol to use the validated "
    "PACKMOL-Memgen solvation path instead."
)


class MembraneSolvationError(RuntimeError):
    """Raised when the fast solvation path cannot produce a trustworthy system.

    Always fatal by design: the alternative -- silently falling back to the
    slow path -- would spend hours re-doing the build with different science
    and no diagnostic.
    """


# --------------------------------------------------------------------------- #
# Leaflet geometry
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class LeafletPlanes:
    """Mean phosphate-plane z of each leaflet, **in nanometres**.

    Nanometres, not Angstrom, because everything downstream (the ``.gro``, the
    purge, the box) is in GROMACS units; :func:`leaflet_planes_from_parm` does
    the conversion once, at the only place an AMBER object is touched.
    """

    lower_z: float
    upper_z: float
    midplane_z: float
    n_lower: int
    n_upper: int

    @property
    def thickness(self) -> float:
        """Phosphate-to-phosphate distance (nm)."""
        return self.upper_z - self.lower_z

    @property
    def core_lower_z(self) -> float:
        """Lower bound of the slab the purge empties (nm)."""
        return self.lower_z + CORE_MARGIN_NM

    @property
    def core_upper_z(self) -> float:
        """Upper bound of the slab the purge empties (nm)."""
        return self.upper_z - CORE_MARGIN_NM

    def shifted(self, dz: float) -> "LeafletPlanes":
        """Return the same planes translated by *dz* nm along z.

        The box-normalisation step moves the whole system rigidly so the
        bilayer midplane sits at ``box_z / 2``; the planes have to follow, and
        a rigid translation is the only transformation applied, so this is
        exact rather than an approximation.
        """
        return LeafletPlanes(
            lower_z=self.lower_z + dz,
            upper_z=self.upper_z + dz,
            midplane_z=self.midplane_z + dz,
            n_lower=self.n_lower,
            n_upper=self.n_upper,
        )

    def describe(self) -> str:
        return (
            f"lower={self.lower_z:.3f} nm ({self.n_lower} P), "
            f"upper={self.upper_z:.3f} nm ({self.n_upper} P), "
            f"thickness={self.thickness:.3f} nm"
        )


def _planes_from_z(z_values: Sequence[float], source: str) -> LeafletPlanes:
    """Split a phosphate z distribution into two leaflets and fit both planes."""
    if _np is None:
        raise MembraneSolvationError("numpy is required for membrane leaflet detection.")
    z = _np.asarray(list(z_values), dtype=float)
    if z.size < 2 * MIN_LEAFLET_P_ATOMS:
        raise MembraneSolvationError(
            f"Only {z.size} lipid phosphorus atom(s) found in {source}; a bilayer needs at "
            f"least {2 * MIN_LEAFLET_P_ATOMS}. {_PACKMOL_HINT}"
        )
    # Median split, not a k-means/clustering step: a bilayer is by construction
    # two well separated populations, and the median is immune to a handful of
    # badly placed headgroups in a way a mean-based split is not.
    split = float(_np.median(z))
    lower = z[z < split]
    upper = z[z >= split]
    if lower.size < MIN_LEAFLET_P_ATOMS or upper.size < MIN_LEAFLET_P_ATOMS:
        raise MembraneSolvationError(
            f"Leaflet detection in {source} found {lower.size}/{upper.size} phosphorus atoms "
            f"per leaflet (minimum {MIN_LEAFLET_P_ATOMS} each). {_PACKMOL_HINT}"
        )
    lower_z = float(lower.mean())
    upper_z = float(upper.mean())
    if upper_z - lower_z < MIN_LEAFLET_SEPARATION_NM:
        raise MembraneSolvationError(
            f"Detected phosphate planes in {source} are only {upper_z - lower_z:.2f} nm apart "
            f"(expected >= {MIN_LEAFLET_SEPARATION_NM} nm for a bilayer); the leaflet split is "
            f"not meaningful. {_PACKMOL_HINT}"
        )
    return LeafletPlanes(
        lower_z=lower_z,
        upper_z=upper_z,
        midplane_z=0.5 * (lower_z + upper_z),
        n_lower=int(lower.size),
        n_upper=int(upper.size),
    )


def leaflet_planes_from_parm(parm, membrane_resnames: Set[str]) -> LeafletPlanes:
    """Locate both phosphate planes in a ParmEd structure (returns nm).

    Phosphorus is identified by ``atomic_number == 15``, never by atom name:
    AMBER Lipid21 splits one POPC into PA/PC/OL fragments whose phosphate atom
    is called ``P31``, CHARMM calls it ``P``, and a name-based rule would also
    sweep up ``PA``-prefixed carbons from other conventions.  The element is
    unambiguous and already in the prmtop.

    Only residues named in *membrane_resnames* contribute, so a phosphorylated
    serine (SEP/S1P) or a phosphate-bearing ligand cannot drag a plane away
    from the bilayer.
    """
    residues = getattr(parm, "residues", None)
    if not residues:
        raise MembraneSolvationError("ParmEd structure has no residues; cannot locate the bilayer.")
    wanted = {str(name).strip().upper() for name in (membrane_resnames or ())}
    if not wanted:
        raise MembraneSolvationError(
            "No membrane residue names supplied; cannot locate the bilayer leaflets."
        )
    z_values: List[float] = []
    for residue in residues:
        if str(residue.name).strip().upper() not in wanted:
            continue
        for atom in residue.atoms:
            if getattr(atom, "atomic_number", None) == 15:
                z_values.append(float(atom.xz) / 10.0)  # ParmEd is Angstrom, GROMACS is nm
    return _planes_from_z(z_values, "the packed AMBER system")


def leaflet_planes_from_gro(gro_path: str, membrane_resnames: Set[str]) -> Optional[LeafletPlanes]:
    """Best-effort plane fit straight from a ``.gro`` (returns nm, or None).

    A ``.gro`` carries no elements, so phosphorus has to be guessed from the
    atom name -- inside a lipid residue an atom name beginning with ``P`` is
    phosphorus in every convention PRISM supports, but the rule is a heuristic
    and so this returns ``None`` instead of raising when it finds nothing.  It
    is used only as a *cross-check* on the planes the caller passes in, never
    as the primary source: a stale (un-shifted) ``LeafletPlanes`` would make
    the purge delete the wrong slab, and that is worth one cheap re-derivation.
    """
    wanted = {str(name).strip().upper() for name in (membrane_resnames or ())}
    if not wanted:
        return None
    try:
        gro = _GroFile.read(gro_path)
    except (OSError, ValueError):
        return None
    z_values = [
        gro.coords[i][2]
        for i in range(gro.natoms)
        if gro.resnames[i].upper() in wanted and _leading_alpha(gro.atomnames[i]) == "P"
    ]
    try:
        return _planes_from_z(z_values, os.path.basename(gro_path))
    except MembraneSolvationError:
        return None


# --------------------------------------------------------------------------- #
# Minimal .gro reader/writer
#
# ParmEd could do this, but the purge has to preserve every byte of the atom
# lines it keeps (velocities, precision, the writer's own column conventions),
# and a round trip through a chemistry object silently reformats all of it.
# --------------------------------------------------------------------------- #


class _GroFile:
    """A ``.gro`` held as raw lines plus the fields the purge actually needs."""

    __slots__ = ("title", "atom_lines", "box_line", "path", "resids", "resnames", "atomnames", "coords")

    def __init__(self, title, atom_lines, box_line, path, resids, resnames, atomnames, coords):
        self.title = title
        self.atom_lines = atom_lines
        self.box_line = box_line
        self.path = path
        self.resids = resids
        self.resnames = resnames
        self.atomnames = atomnames
        self.coords = coords

    @property
    def natoms(self) -> int:
        return len(self.atom_lines)

    @classmethod
    def read(cls, path: str) -> "_GroFile":
        with open(path, "r") as fh:
            lines = fh.read().splitlines()
        if len(lines) < 3:
            raise ValueError(f"{path} is too short to be a .gro file.")
        try:
            declared = int(lines[1].strip())
        except ValueError as exc:
            raise ValueError(f"{path} line 2 is not an atom count.") from exc
        atom_lines = lines[2 : 2 + declared]
        if len(atom_lines) != declared:
            raise ValueError(
                f"{path} declares {declared} atoms but only {len(atom_lines)} atom lines follow."
            )
        if len(lines) < 2 + declared + 1:
            raise ValueError(f"{path} has no box line.")
        box_line = lines[2 + declared]

        resids: List[int] = []
        resnames: List[str] = []
        atomnames: List[str] = []
        coords: List[Tuple[float, float, float]] = []
        for line in atom_lines:
            resids.append(_int_or_zero(line[0:5]))
            resnames.append(line[5:10].strip())
            atomnames.append(line[10:15].strip())
            coords.append(_parse_gro_xyz(line))
        return cls(lines[0], atom_lines, box_line, path, resids, resnames, atomnames, coords)

    def box(self) -> Tuple[float, float, float]:
        parts = self.box_line.split()
        if len(parts) < 3:
            raise ValueError(f"{self.path} box line has fewer than three fields: {self.box_line!r}")
        return float(parts[0]), float(parts[1]), float(parts[2])

    def residue_slices(self) -> List[Tuple[int, int]]:
        """Half-open ``(start, stop)`` atom ranges, one per residue.

        A new residue starts wherever the residue number or the residue name
        changes.  ``.gro`` residue numbers wrap at 100000, but consecutive
        residues still differ, so the rule survives a 200k-atom membrane box.
        """
        slices: List[Tuple[int, int]] = []
        start = 0
        for i in range(1, self.natoms):
            if self.resids[i] != self.resids[i - 1] or self.resnames[i] != self.resnames[i - 1]:
                slices.append((start, i))
                start = i
        if self.natoms:
            slices.append((start, self.natoms))
        return slices

    def write(
        self,
        out_path: str,
        keep: Optional[Sequence[bool]] = None,
        rename: Optional[Dict[int, Tuple[str, str]]] = None,
    ) -> int:
        """Write the file back, optionally dropping atoms and renaming some.

        Everything from column 20 onward (coordinates, velocities, precision)
        is copied byte for byte; only the residue number, residue name, atom
        name and atom number fields are ever rewritten.  Both numbers are
        cosmetic to GROMACS -- it reads a ``.gro`` positionally -- but leaving
        gaps in them makes every downstream ``gmx`` message harder to read.
        """
        rename = rename or {}
        out_lines = [self.title]
        body: List[str] = []
        resnum = 0
        atomnum = 0
        for start, stop in self.residue_slices():
            if keep is not None and not keep[start]:
                continue
            resnum += 1
            for i in range(start, stop):
                atomnum += 1
                line = self.atom_lines[i]
                resname, atomname = rename.get(i, (self.resnames[i], self.atomnames[i]))
                body.append(
                    f"{resnum % 100000:5d}{resname[:5]:<5s}{atomname[:5]:>5s}"
                    f"{atomnum % 100000:5d}{line[20:]}"
                )
        out_lines.append(f"{len(body):5d}")
        out_lines.extend(body)
        out_lines.append(self.box_line)
        tmp_path = out_path + ".tmp"
        with open(tmp_path, "w") as fh:
            fh.write("\n".join(out_lines) + "\n")
        os.replace(tmp_path, out_path)  # atomic: a half-written .gro must never survive
        return len(body)


def _int_or_zero(text: str) -> int:
    try:
        return int(text.strip())
    except ValueError:
        return 0


def _parse_gro_xyz(line: str) -> Tuple[float, float, float]:
    """Extract x/y/z (nm) from a ``.gro`` atom line.

    ``.gro`` coordinate fields are fixed width (``%w.df`` with ``w = d + 5``),
    and the width is only discoverable from the line length because GROMACS
    supports ``-ndec``.  Whitespace splitting is the fallback, not the rule:
    adjacent negative velocities can run together and would silently shift the
    field indices.
    """
    tail = line[20:].rstrip()
    for width in (8, 9, 10, 11, 12, 13, 14, 15, 16):
        if len(tail) in (3 * width, 6 * width, 9 * width):
            try:
                return (
                    float(tail[0:width]),
                    float(tail[width : 2 * width]),
                    float(tail[2 * width : 3 * width]),
                )
            except ValueError:
                break
    parts = tail.split()
    if len(parts) < 3:
        raise ValueError(f"Cannot parse coordinates from .gro line: {line!r}")
    return float(parts[0]), float(parts[1]), float(parts[2])


def _leading_alpha(name: str) -> str:
    """First alphabetic character of an atom name, upper-cased ('' if none).

    PDB/gro atom names may be digit-prefixed (``1HB``), so the element guess
    has to skip leading digits rather than read column 0.
    """
    for ch in name:
        if ch.isalpha():
            return ch.upper()
    return ""


def _is_hydrogen_name(name: str) -> bool:
    return _leading_alpha(name) == "H"


# --------------------------------------------------------------------------- #
# Minimal GROMACS topology reader
#
# Only what the guards need: per-moleculetype atom count and total charge, and
# the [ molecules ] table.  ``#include`` is deliberately *not* followed: the
# membrane topology ParmEd emits is self-contained apart from the position
# restraint files, which carry no atoms and no charges.
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class MoleculeType:
    name: str
    natoms: int
    charge: float
    #: False when at least one ``[ atoms ]`` line omitted the charge column
    #: (legal in GROMACS -- it then falls back to the atomtype).  Net-charge
    #: arithmetic must refuse to guess in that case rather than under-count.
    charge_known: bool
    atom_names: Tuple[str, ...]
    residue_names: Tuple[str, ...]


@dataclass(frozen=True)
class TopologyView:
    moleculetypes: Dict[str, MoleculeType]
    molecules: Tuple[Tuple[str, int], ...]
    path: str

    def count_of(self, name: str) -> int:
        return sum(count for mol, count in self.molecules if mol == name)

    def natoms(self) -> int:
        total = 0
        for mol, count in self.molecules:
            mtype = self.moleculetypes.get(mol)
            if mtype is None:
                raise MembraneSolvationError(
                    f"[ molecules ] in {os.path.basename(self.path)} references moleculetype "
                    f"{mol!r}, which is not defined in the file."
                )
            total += mtype.natoms * count
        return total

    def net_charge(self) -> float:
        total = 0.0
        for mol, count in self.molecules:
            mtype = self.moleculetypes.get(mol)
            if mtype is None:
                raise MembraneSolvationError(
                    f"[ molecules ] references undefined moleculetype {mol!r}."
                )
            if not mtype.charge_known:
                raise MembraneSolvationError(
                    f"Moleculetype {mol!r} has [ atoms ] lines without an explicit charge column, "
                    "so the system net charge cannot be computed. " + _PACKMOL_HINT
                )
            total += mtype.charge * count
        return total


_DIRECTIVE_RE = re.compile(r"^\s*\[\s*([A-Za-z_0-9]+)\s*\]")


def read_topology(path: str) -> TopologyView:
    """Parse the parts of a GROMACS ``.top`` the fast-path guards depend on."""
    moleculetypes: Dict[str, MoleculeType] = {}
    molecules: List[Tuple[str, int]] = []

    directive = ""
    current_name: Optional[str] = None
    natoms = 0
    charge = 0.0
    charge_known = True
    atom_names: List[str] = []
    residue_names: List[str] = []

    def _flush():
        nonlocal current_name, natoms, charge, charge_known, atom_names, residue_names
        if current_name is not None:
            if current_name in moleculetypes:
                raise MembraneSolvationError(
                    f"Duplicate [ moleculetype ] {current_name!r} in {os.path.basename(path)}."
                )
            moleculetypes[current_name] = MoleculeType(
                name=current_name,
                natoms=natoms,
                charge=charge,
                charge_known=charge_known,
                atom_names=tuple(atom_names),
                residue_names=tuple(residue_names),
            )
        current_name = None
        natoms = 0
        charge = 0.0
        charge_known = True
        atom_names = []
        residue_names = []

    with open(path, "r") as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].rstrip()
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                # Preprocessor.  Restraint includes sit inside #ifdef POSRES and
                # contribute neither atoms nor charge, so skipping is exact.
                continue
            match = _DIRECTIVE_RE.match(line)
            if match:
                new_directive = match.group(1).lower()
                if new_directive == "moleculetype":
                    _flush()
                directive = new_directive
                continue
            fields = line.split()
            if directive == "moleculetype":
                if current_name is None:
                    current_name = fields[0]
            elif directive == "atoms" and current_name is not None:
                natoms += 1
                residue_names.append(fields[3] if len(fields) > 3 else "")
                atom_names.append(fields[4] if len(fields) > 4 else "")
                if len(fields) > 6:
                    try:
                        charge += float(fields[6])
                    except ValueError:
                        charge_known = False
                else:
                    charge_known = False
            elif directive == "molecules":
                if len(fields) >= 2:
                    try:
                        molecules.append((fields[0], int(fields[1])))
                    except ValueError:
                        raise MembraneSolvationError(
                            f"Malformed [ molecules ] entry in {os.path.basename(path)}: {line!r}"
                        )
    _flush()

    if not molecules:
        raise MembraneSolvationError(
            f"{os.path.basename(path)} has no [ molecules ] section; it is not a complete topology."
        )
    return TopologyView(moleculetypes=moleculetypes, molecules=tuple(molecules), path=path)


def _set_molecule_count(top_path: str, name: str, count: int) -> None:
    """Rewrite a single ``[ molecules ]`` entry in place.

    Refuses to touch a topology where the name appears more than once: the
    molecule order in that section is the coordinate order, so silently
    collapsing two blocks would scramble the system.
    """
    with open(top_path, "r") as fh:
        lines = fh.readlines()

    in_molecules = False
    hits: List[int] = []
    for i, raw in enumerate(lines):
        stripped = raw.split(";", 1)[0].strip()
        if not stripped:
            continue
        match = _DIRECTIVE_RE.match(raw)
        if match:
            in_molecules = match.group(1).lower() == "molecules"
            continue
        if in_molecules and stripped.split()[0] == name:
            hits.append(i)

    if len(hits) != 1:
        raise MembraneSolvationError(
            f"Expected exactly one [ molecules ] entry named {name!r} in "
            f"{os.path.basename(top_path)}, found {len(hits)}."
        )
    lines[hits[0]] = f"{name:<18} {count}\n"
    with open(top_path, "w") as fh:
        fh.writelines(lines)


# --------------------------------------------------------------------------- #
# Result
# --------------------------------------------------------------------------- #


@dataclass
class SolvationOutcome:
    """Everything the builder needs to report -- and to trust -- the fast path."""

    gro: str
    top: str
    n_water: int
    n_pos: int
    n_neg: int
    purged_core: int
    retained_pore: int
    effective_conc: float
    warnings: List[str] = field(default_factory=list)
    #: Water count ``gmx solvate`` inserted, before the core purge.
    inserted_water: int = 0
    #: The value guard G4b compared against, kept for the build log.
    expected_bulk_water: int = 0
    vdwradii_carbon: float = DEFAULT_CARBON_RADIUS_NM


# --------------------------------------------------------------------------- #
# The solvator
# --------------------------------------------------------------------------- #


class MembraneSolvator(SystemBuilderBase, SolvationProcessorMixin, TopologyProcessorMixin):
    """Solvate and ionise a dry, parametrised membrane system with GROMACS.

    Composition, not inheritance into ``MembraneBuilder``: the mixins need a
    ``config`` **dict** (``self.config["ions"]``) and an ``output_dir`` at
    construction time, whereas ``MembraneBuilder.config`` is a
    :class:`~prism.membrane.config.MembraneConfig` dataclass that is not
    subscriptable and whose output directory only arrives in ``build()``.  This
    class satisfies the mixin contract (``model_dir`` / ``overwrite`` /
    ``gmx_command`` / ``water_model_name`` / ``config['ions']``) and nothing
    else has to change.

    ``SolvationProcessorMixin._solvate`` is reused verbatim -- it already runs
    with ``cwd=self.model_dir``, which is exactly where ``gmx solvate`` looks
    for ``vdwradii.dat``, and it owns the water-model -> ``-cs`` mapping.
    ``_add_ions`` is *not* reused; see the module docstring.
    ``TopologyProcessorMixin`` is inherited for its grompp error-recovery
    helpers (``_fix_improper_from_grompp_errors``); its
    ``_update_topology_molecules`` is deliberately unused because its
    verification is vacuous on a membrane topology.
    """

    def __init__(
        self,
        model_dir: str,
        gmx_command: str,
        water_model: str,
        positive_ion: str,
        negative_ion: str,
        salt_concentration: float,
        overwrite: bool = False,
        timeout: int = 3600,
        verbose: bool = True,
    ) -> None:
        model_dir = os.path.abspath(model_dir)
        config = {
            "general": {"gmx_command": gmx_command},
            # The mixins read this shape; the fast path reads the attributes
            # below.  Both are populated so a future mixin reuse cannot find a
            # half-filled config.
            "ions": {
                "positive_ion": positive_ion,
                "negative_ion": negative_ion,
                "concentration": salt_concentration,
                "neutral": True,
            },
        }
        super().__init__(
            config=config,
            output_dir=os.path.dirname(model_dir) or ".",
            overwrite=overwrite,
            model_dirname=os.path.basename(model_dir),
        )
        self.water_model_name = str(water_model or "tip3p").strip().lower()
        self.positive_ion = str(positive_ion).strip()
        self.negative_ion = str(negative_ion).strip()
        self.salt_concentration = float(salt_concentration)
        self.timeout = int(timeout)
        self.verbose = bool(verbose)
        self._last_stdout = ""
        self._last_stderr = ""

    # -- plumbing ---------------------------------------------------------- #

    def _run_command(self, command, work_dir, input_str=None, env=None, timeout=_TIMEOUT_DEFAULT):
        """Capture the last tool output and apply this builder's timeout.

        Overridden rather than duplicated because ``_solvate`` (reused from the
        soluble path) discards stdout/stderr, and guard G4a has to read
        ``gmx solvate``'s warnings.  Delegating keeps the timeout/kill/report
        behaviour identical to every other PRISM subprocess.
        """
        effective = self.timeout if timeout is _TIMEOUT_DEFAULT else timeout
        stdout, stderr = super()._run_command(
            command, work_dir, input_str=input_str, env=env, timeout=effective
        )
        self._last_stdout, self._last_stderr = stdout or "", stderr or ""
        return stdout, stderr

    def _say(self, message: str) -> None:
        if self.verbose:
            print_info(message)

    # -- public entry point ------------------------------------------------ #

    def run(
        self,
        dry_gro: str,
        top: str,
        planes: LeafletPlanes,
        solute_resnames: Set[str],
        carbon_radius: float = DEFAULT_CARBON_RADIUS_NM,
        purge: bool = True,
        membrane_resnames: Optional[Set[str]] = None,
    ) -> SolvationOutcome:
        """Fill the aqueous slab, empty the hydrophobic core, add ions.

        Args:
            dry_gro: Coordinates of the water-free, box-normalised system.
            top: Topology for those coordinates.  Modified **in place** (that
                is what ``gmx solvate -p`` and ``gmx genion -p`` do), so the
                caller owns the backup of the pre-solvation state.
            planes: Leaflet planes for *these* coordinates, in nm.  If the
                caller normalised the box after fitting them, pass
                ``planes.shifted(dz)``, not the original object.
            solute_resnames: Protein (and embedded-ligand) residue names.  Only
                these protect a core water from the purge -- lipids must not,
                or the purge would keep exactly the waters it exists to remove.
            carbon_radius: Carbon vdW radius (nm) for the solvation-only
                ``vdwradii.dat``.
            purge: Set False only for research; the system will then contain
                water in the hydrophobic core.
            membrane_resnames: Optional; enables an independent re-derivation
                of the leaflet planes from *dry_gro* as a cross-check on the
                ``planes`` argument.  Extra keyword beyond the agreed design
                signature -- callers that omit it lose only the cross-check.

        Returns:
            A :class:`SolvationOutcome` whose ``gro`` is the final coordinate
            file (``solv_ions.gro``, or ``solv.gro`` when no ions were needed).
            The caller is responsible for publishing it under the deliverable
            name (``system.gro``).
        """
        if _np is None:
            raise MembraneSolvationError(
                "numpy is required for the GROMACS membrane solvation path. " + _PACKMOL_HINT
            )
        dry_gro = os.path.abspath(dry_gro)
        top = os.path.abspath(top)
        for required in (dry_gro, top):
            if not os.path.isfile(required):
                raise MembraneSolvationError(f"Missing input for membrane solvation: {required}")

        warnings: List[str] = []
        carbon_radius = float(carbon_radius)
        if not math.isfinite(carbon_radius) or carbon_radius <= 0:
            raise MembraneSolvationError(f"Carbon vdW radius must be positive, got {carbon_radius!r}.")

        # --- G3c / contract: the topology must already publish SOL and the ions
        view = read_topology(top)
        water_sites = self._check_solvent_moleculetypes(view)
        # --- G1: refuse to add water to a system that already has some
        self._require_dry(view)

        # --- plane sanity + independent cross-check against the coordinates
        self._validate_planes(planes)
        self._crosscheck_planes(dry_gro, planes, membrane_resnames, warnings)

        # --- solvate (G4a / G4b / G4c live inside)
        solv_gro, inserted, expected = self._solvate_membrane(
            dry_gro, top, planes, carbon_radius, warnings
        )

        # --- make the coordinate atom names agree with the topology's water.
        # GROMACS ignores .gro atom names, but grompp counts the mismatches and
        # reports "N non-matching atom names"; with ~200k waters that single
        # warning would be indistinguishable from a genuinely scrambled system
        # (guard G7).  Renaming positionally is safe precisely because the
        # names are ignored.
        rename = self._solvent_rename_map(solv_gro, view, water_sites, warnings)

        purged = 0
        retained = 0
        if purge:
            purged, retained = self._purge_core_water(
                solv_gro, top, planes, solute_resnames, water_sites, inserted, rename, warnings
            )
        else:
            warnings.append(
                "Core-water purge disabled: the system will contain water inside the hydrophobic "
                "core of the bilayer. This is only ever appropriate for methodological tests."
            )
            if rename:
                gro = _GroFile.read(solv_gro)
                gro.write(solv_gro, rename=rename)

        n_water = inserted - purged

        # --- ions
        n_pos, n_neg, effective_conc = self._ionise(solv_gro, top, n_water, warnings)
        final_gro = os.path.join(str(self.model_dir), "solv_ions.gro")
        if n_pos == 0 and n_neg == 0:
            # Nothing to place: publish the solvated coordinates under the name
            # the rest of the pipeline expects rather than branching downstream.
            shutil.copyfile(solv_gro, final_gro)

        # --- G6d: ions must not have landed in the hydrophobic core
        self._assert_no_core_ions(final_gro, top, planes, solute_resnames, n_pos, n_neg, warnings)

        outcome = SolvationOutcome(
            gro=final_gro,
            top=top,
            n_water=n_water - n_pos - n_neg,
            n_pos=n_pos,
            n_neg=n_neg,
            purged_core=purged,
            retained_pore=retained,
            effective_conc=effective_conc,
            warnings=warnings,
            inserted_water=inserted,
            expected_bulk_water=expected,
            vdwradii_carbon=carbon_radius,
        )
        if self.verbose:
            print_success(
                f"GROMACS solvation: {number(outcome.n_water)} water, "
                f"{number(n_pos)} {self.positive_ion}, {number(n_neg)} {self.negative_ion}, "
                f"{number(purged)} core waters purged ({number(retained)} pore waters kept), "
                f"effective salt {effective_conc:.3f} M"
            )
        return outcome

    # -- guards on the incoming topology ----------------------------------- #

    def _check_solvent_moleculetypes(self, view: TopologyView) -> int:
        """Guard G3c and the moleculetype-name contract.

        Returns the number of sites in the water moleculetype.

        ``gmx solvate`` appends ``SOL <n>`` to ``[ molecules ]`` using the
        residue name baked into the ``-cs`` file, and ``gmx genion -pname X``
        needs a moleculetype literally called ``X``.  Both are contracts the
        solvent-parameter blocks must already satisfy; checking them here turns
        two obscure downstream failures ("No such moleculetype NA", a grompp
        atom-count mismatch) into one actionable message before any work is done.
        """
        water = view.moleculetypes.get(WATER_MOLECULETYPE)
        if water is None:
            raise MembraneSolvationError(
                f"Topology {os.path.basename(view.path)} has no [ moleculetype ] named "
                f"{WATER_MOLECULETYPE!r}. gmx solvate writes its water under that name, so the "
                "solvent parameter blocks must be merged into the topology first. " + _PACKMOL_HINT
            )
        coord_file = _WATER_MODEL_COORD_FILE.get(self.water_model_name, "spc216.gro")
        expected_sites = _COORD_FILE_SITES[coord_file]
        if water.natoms != expected_sites:
            raise MembraneSolvationError(
                f"Water model mismatch: config asks for {self.water_model_name!r}, which "
                f"gmx solvate fills from {coord_file} ({expected_sites} sites), but the "
                f"{WATER_MOLECULETYPE!r} moleculetype in the topology has {water.natoms} atom(s). "
                "A 3-site block silently satisfying a 4-site request would give every water the "
                "wrong charge distribution. " + _PACKMOL_HINT
            )
        for role, name in (("positive", self.positive_ion), ("negative", self.negative_ion)):
            if name not in view.moleculetypes:
                available = ", ".join(sorted(view.moleculetypes)[:12]) or "<none>"
                raise MembraneSolvationError(
                    f"Topology has no [ moleculetype ] named {name!r} for the {role} ion; "
                    f"gmx genion -{'p' if role == 'positive' else 'n'}name would fail with "
                    f"'No such moleculetype'. Defined moleculetypes include: {available}. "
                    + _PACKMOL_HINT
                )
        return water.natoms

    @staticmethod
    def _require_dry(view: TopologyView) -> None:
        """Guard G1: the packed system must contain no water and no ions.

        If PACKMOL already placed the solvent, adding more would double-count
        it -- and the pack would have taken the hours the fast path exists to
        avoid, so a wet system here means the dry flags did not take effect.
        """
        wet = [
            (name, count)
            for name, count in view.molecules
            if count > 0 and name.strip().upper().rstrip("+-") in {
                n.rstrip("+-") for n in _SOLVENT_MOLECULE_NAMES
            }
        ]
        if wet:
            listed = ", ".join(f"{name} x{count}" for name, count in wet)
            raise MembraneSolvationError(
                "The packed system is not dry -- [ molecules ] already contains solvent/ions "
                f"({listed}). Solvating it again would double-count the aqueous phase. "
                "Check that packmol-memgen ran with --nocounter and a near-zero-density "
                "--solvent_parm. " + _PACKMOL_HINT
            )

    @staticmethod
    def _validate_planes(planes: LeafletPlanes) -> None:
        if not isinstance(planes, LeafletPlanes):
            raise MembraneSolvationError(
                f"Expected LeafletPlanes, got {type(planes).__name__}; the purge cannot run "
                "without the leaflet geometry."
            )
        if not (
            math.isfinite(planes.lower_z) and math.isfinite(planes.upper_z)
        ) or planes.thickness < MIN_LEAFLET_SEPARATION_NM:
            raise MembraneSolvationError(
                f"Implausible leaflet planes ({planes.describe()}); refusing to purge. "
                + _PACKMOL_HINT
            )

    def _crosscheck_planes(
        self,
        gro_path: str,
        planes: LeafletPlanes,
        membrane_resnames: Optional[Set[str]],
        warnings: List[str],
    ) -> None:
        """Re-derive the planes from the coordinates and compare.

        The single most damaging way this module can fail is to be handed
        planes fitted *before* the box normalisation shifted every atom: the
        purge would then empty a slab of bulk water and leave the core wet,
        and every other guard would pass.  One cheap re-fit removes that class
        of bug entirely.
        """
        if not membrane_resnames:
            return
        observed = leaflet_planes_from_gro(gro_path, membrane_resnames)
        if observed is None:
            warnings.append(
                "Could not re-derive the leaflet planes from the coordinates (no lipid "
                "phosphorus atom names recognised); the purge is running on the planes supplied "
                "by the caller without an independent cross-check."
            )
            return
        drift = max(abs(observed.lower_z - planes.lower_z), abs(observed.upper_z - planes.upper_z))
        if drift > 0.5:
            raise MembraneSolvationError(
                "Leaflet planes supplied to the solvator disagree with the coordinates by "
                f"{drift:.2f} nm (supplied {planes.describe()}; measured {observed.describe()}). "
                "This usually means the planes were not translated with the box normalisation. "
                + _PACKMOL_HINT
            )

    # -- solvation --------------------------------------------------------- #

    def _solvate_membrane(
        self,
        dry_gro: str,
        top: str,
        planes: LeafletPlanes,
        carbon_radius: float,
        warnings: List[str],
    ) -> Tuple[str, int, int]:
        """Run ``gmx solvate`` under inflated carbon radii.  Returns (gro, n, expected)."""
        model_dir = str(self.model_dir)
        solv_gro = os.path.join(model_dir, "solv.gro")
        # A leftover solv.gro from an interrupted run would make the reused
        # _solvate() return early without ever touching the topology, leaving a
        # coordinate file with waters and a [ molecules ] section without them.
        for stale in (solv_gro, os.path.join(model_dir, "solv_ions.gro")):
            if os.path.exists(stale):
                os.remove(stale)

        before = view_water_count(top)
        radii_path = self._write_vdwradii(carbon_radius)
        try:
            self._say(
                f"Solvating with inflated carbon radius C={carbon_radius:.3f} nm "
                f"(stock 0.17 nm) to keep water out of the acyl slab"
            )
            self._solvate(dry_gro, top)
            self._check_solvate_output(self._last_stdout, self._last_stderr, warnings)
        finally:
            # Guard G4c: vdwradii.dat is read from the working directory by
            # every later GROMACS tool too.  Leaving it behind would silently
            # change any subsequent solvation in this directory.
            if os.path.exists(radii_path):
                os.remove(radii_path)
        if os.path.exists(radii_path):  # pragma: no cover - unlink failed silently
            raise MembraneSolvationError(
                f"Could not remove {radii_path}; refusing to leave inflated vdW radii where a "
                "later GROMACS command would pick them up."
            )
        if not os.path.isfile(solv_gro):
            raise MembraneSolvationError("gmx solvate did not produce solv.gro.")

        inserted = view_water_count(top) - before
        if inserted <= 0:
            raise MembraneSolvationError(
                "gmx solvate added no water to the topology. The box is probably not larger than "
                "the packed system -- check the box-normalisation step. " + _PACKMOL_HINT
            )
        # Cross-check the topology against the coordinates: gmx solvate writes
        # both, and a disagreement means one of them is stale.
        gro = _GroFile.read(solv_gro)
        gro_waters = sum(1 for name in gro.resnames if name == WATER_MOLECULETYPE)
        sites = read_topology(top).moleculetypes[WATER_MOLECULETYPE].natoms
        if gro_waters != inserted * sites:
            raise MembraneSolvationError(
                f"gmx solvate wrote {gro_waters} water atoms into solv.gro but recorded "
                f"{inserted} x {sites} in the topology; the two files disagree."
            )

        expected = self._expected_bulk_waters(gro, planes)
        self._check_bulk_fill(inserted, expected, carbon_radius)
        self._say(
            f"gmx solvate inserted {number(inserted)} water molecules "
            f"(expected ~{number(expected)} for the aqueous volume)"
        )
        return solv_gro, inserted, expected

    def _write_vdwradii(self, carbon_radius: float) -> str:
        """Write a solvation-only ``vdwradii.dat`` with an inflated carbon radius.

        Only the ``??? C`` wildcard line is edited.  Nothing is *added*:
        GROMACS resolves radii by longest match, so a residue-scoped entry such
        as ``PA C 0.375`` silently overrides the wildcard for that residue, and
        the equivalent trick on hydrogen inflates SPC water's own hydrogens and
        under-fills the box by an order of magnitude with no error at all.
        """
        stock = self._stock_vdwradii()
        with open(stock, "r") as fh:
            lines = fh.readlines()

        edited = 0
        shadowed: List[str] = []
        out: List[str] = []
        for raw in lines:
            body = raw.split(";", 1)[0]
            fields = body.split()
            if len(fields) >= 3 and fields[1] == "C":
                if fields[0] in ("???", "*"):
                    out.append(f"{fields[0]:<4s} {fields[1]:<5s} {carbon_radius:.3f}\n")
                    edited += 1
                    continue
                # A residue-scoped carbon entry beats our wildcard for that
                # residue; report it rather than let it quietly defeat the
                # inflation on exactly the lipid we care about.
                shadowed.append(fields[0])
            out.append(raw)

        if edited != 1:
            raise MembraneSolvationError(
                f"Expected exactly one '??? C' wildcard line in {stock}, found {edited}. "
                "Refusing to solvate a bilayer with unknown carbon radii. " + _PACKMOL_HINT
            )
        if shadowed:
            raise MembraneSolvationError(
                f"{stock} contains residue-scoped carbon radii for {', '.join(sorted(set(shadowed)))}; "
                "GROMACS prefers the longest match, so those residues would be solvated with the "
                "stock radius and could take water into the hydrophobic core. " + _PACKMOL_HINT
            )

        radii_path = os.path.join(str(self.model_dir), "vdwradii.dat")
        with open(radii_path, "w") as fh:
            fh.writelines(out)
        return radii_path

    def _stock_vdwradii(self) -> str:
        """Locate the GROMACS-shipped ``vdwradii.dat``.

        Hand-writing a replacement is not an option: the stock file also
        carries the zero-radius entries for virtual sites (``SOL MW``, the
        vsite construction masses), and omitting them makes ``gmx solvate``
        treat every dummy as a real sphere.
        """
        candidates: List[str] = []
        for env_var, suffixes in (("GMXDATA", ("top/vdwradii.dat",)), ("GMXLIB", ("vdwradii.dat", "top/vdwradii.dat"))):
            root = os.environ.get(env_var)
            if root:
                candidates.extend(os.path.join(root, suffix) for suffix in suffixes)
        binary = shutil.which(self.gmx_command)
        if binary:
            prefix = os.path.dirname(os.path.dirname(os.path.realpath(binary)))
            candidates.append(os.path.join(prefix, "share", "gromacs", "top", "vdwradii.dat"))
            candidates.append(os.path.join(prefix, "share", "top", "vdwradii.dat"))
        for candidate in candidates:
            if os.path.isfile(candidate):
                return candidate
        raise MembraneSolvationError(
            "Could not find the GROMACS vdwradii.dat data file (looked in $GMXDATA, $GMXLIB and "
            f"next to '{self.gmx_command}'). Source GMXRC so $GMXDATA is set, or "
            + _PACKMOL_HINT
        )

    @staticmethod
    def _check_solvate_output(stdout: str, stderr: str, warnings: List[str]) -> None:
        """Guard G4a: the radii file must not have been silently mangled."""
        text = f"{stdout}\n{stderr}"
        if "warning double" in text.lower():
            raise MembraneSolvationError(
                "gmx solvate reported duplicate vdW radius entries ('Warning double ...'), which "
                "means the membrane vdwradii.dat is being shadowed by another entry and the "
                "inflated carbon radius may not have applied. " + _PACKMOL_HINT
            )
        for line in text.splitlines():
            if "Van der Waals radius" in line and "could not find" in line.lower():
                warnings.append(f"gmx solvate: {line.strip()}")

    def _expected_bulk_waters(self, gro: "_GroFile", planes: LeafletPlanes) -> int:
        """Coarse prediction of the bulk water count for guard G4b.

        Aqueous volume = box minus the acyl slab minus the volume the solute
        and headgroup heavy atoms outside that slab displace.  Deliberately
        crude: it exists to separate "filled" from "vacuum", and the tolerance
        is +/-25%.
        """
        box_x, box_y, box_z = gro.box()
        slab = max(0.0, min(planes.upper_z, box_z) - max(planes.lower_z, 0.0))
        aqueous = box_x * box_y * max(0.0, box_z - slab)
        heavy_outside = 0
        for i in range(gro.natoms):
            if gro.resnames[i] == WATER_MOLECULETYPE or _is_hydrogen_name(gro.atomnames[i]):
                continue
            z = gro.coords[i][2]
            if z < planes.lower_z or z > planes.upper_z:
                heavy_outside += 1
        aqueous -= heavy_outside * HEAVY_ATOM_VOLUME_NM3
        return int(max(0.0, aqueous) * WATER_NUMBER_DENSITY_NM3)

    @staticmethod
    def _check_bulk_fill(inserted: int, expected: int, carbon_radius: float) -> None:
        if expected <= 0:
            raise MembraneSolvationError(
                "Computed a non-positive aqueous volume for the solvated box; the box is not "
                "larger than the bilayer. " + _PACKMOL_HINT
            )
        ratio = inserted / float(expected)
        if abs(ratio - 1.0) > BULK_WATER_TOLERANCE:
            raise MembraneSolvationError(
                f"gmx solvate inserted {inserted} waters where ~{expected} were expected "
                f"({ratio:.2f}x). Under-filling usually means the carbon radius "
                f"({carbon_radius:.3f} nm) is too large for this system; over-filling means the "
                "leaflet planes or the box are wrong. Adjust membrane.vdwradii_carbon or "
                + _PACKMOL_HINT
            )

    # -- solvent atom-name alignment --------------------------------------- #

    def _solvent_rename_map(
        self, gro_path: str, view: TopologyView, water_sites: int, warnings: List[str]
    ) -> Dict[int, Tuple[str, str]]:
        """Map gro atom index -> (resname, atomname) taken from the topology.

        ``gmx solvate`` names its water after the ``-cs`` file (OW/HW1/HW2),
        while an AMBER-derived ``SOL`` block names it O/H1/H2.  GROMACS ignores
        ``.gro`` atom names ("atom names from topol.top will be used"), so
        renaming positionally changes nothing physical -- but it takes the
        grompp mismatch count from 3*N_water to zero, which is what makes the
        "N non-matching atom names" warning usable as a scrambled-system
        detector (guard G7).
        """
        water = view.moleculetypes[WATER_MOLECULETYPE]
        target_names = list(water.atom_names)
        target_resname = water.residue_names[0] if water.residue_names else WATER_MOLECULETYPE
        if len(target_names) != water_sites:
            return {}

        gro = _GroFile.read(gro_path)
        rename: Dict[int, Tuple[str, str]] = {}
        for start, stop in gro.residue_slices():
            if gro.resnames[start] != WATER_MOLECULETYPE:
                continue
            if stop - start != water_sites:
                raise MembraneSolvationError(
                    f"A {WATER_MOLECULETYPE} residue in {os.path.basename(gro_path)} has "
                    f"{stop - start} atoms, expected {water_sites}."
                )
            for offset, i in enumerate(range(start, stop)):
                current = gro.atomnames[i]
                target = target_names[offset]
                if current == target:
                    continue
                # Structural check before trusting positional correspondence:
                # oxygen must stay oxygen and hydrogen must stay hydrogen.  The
                # virtual site (MW vs EPW) is unconstrained -- it has no element.
                cur_el, tgt_el = _leading_alpha(current), _leading_alpha(target)
                if {cur_el, tgt_el} & {"O", "H"} and cur_el != tgt_el:
                    warnings.append(
                        f"Not renaming water atoms: site {offset + 1} is {current!r} in the "
                        f"coordinates but {target!r} in the topology, which are different "
                        "elements. grompp will report non-matching atom names."
                    )
                    return {}
                rename[i] = (target_resname, target)
            if gro.resnames[start] != target_resname:
                for i in range(start, stop):
                    rename.setdefault(i, (target_resname, gro.atomnames[i]))
        if rename:
            self._say(
                f"Aligning {number(len(rename))} water atom name(s) in the coordinates with the "
                "topology so grompp's mismatch count stays a meaningful signal"
            )
        return rename

    # -- the purge --------------------------------------------------------- #

    def _purge_core_water(
        self,
        gro_path: str,
        top_path: str,
        planes: LeafletPlanes,
        solute_resnames: Set[str],
        water_sites: int,
        inserted: int,
        rename: Dict[int, Tuple[str, str]],
        warnings: List[str],
    ) -> Tuple[int, int]:
        """Delete every water in the hydrophobic core.  Returns (purged, retained).

        This -- not ``vdwradii.dat`` -- is what reaches zero core waters.
        Inflating carbon is asymptotic: measured on real packs it takes the
        residue from ~2000 to ~600 and no radius takes it to 0, because a
        freshly packed bilayer has voids larger than any plausible radius.  The
        purge is geometric and therefore exact.
        """
        solute = {str(name).strip().upper() for name in (solute_resnames or ())}
        gro = _GroFile.read(gro_path)
        box = gro.box()
        coords = _np.asarray(gro.coords, dtype=float)

        solute_mask = _np.array(
            [
                gro.resnames[i].upper() in solute and not _is_hydrogen_name(gro.atomnames[i])
                for i in range(gro.natoms)
            ],
            dtype=bool,
        )
        # Only solute heavy atoms near the slab can shield a core water, and
        # restricting the reference set keeps the fallback neighbour search
        # cheap without changing a single answer.
        band = (coords[:, 2] > planes.core_lower_z - PORE_RADIUS_NM) & (
            coords[:, 2] < planes.core_upper_z + PORE_RADIUS_NM
        )
        reference = coords[solute_mask & band]

        # Candidate waters: oxygen inside the trimmed acyl slab.
        slices = gro.residue_slices()
        candidates: List[Tuple[int, int, int]] = []  # (start, stop, oxygen index)
        for start, stop in slices:
            if gro.resnames[start] != WATER_MOLECULETYPE:
                continue
            if stop - start != water_sites:
                raise MembraneSolvationError(
                    f"Water residue with {stop - start} atoms (expected {water_sites}); refusing "
                    "to purge a coordinate file whose water blocks are irregular."
                )
            oxygen = next(
                (i for i in range(start, stop) if _leading_alpha(gro.atomnames[i]) == "O"), start
            )
            z = coords[oxygen][2]
            if planes.core_lower_z < z < planes.core_upper_z:
                candidates.append((start, stop, oxygen))

        retained = 0
        purge_starts: List[int] = []
        if candidates:
            query = coords[[c[2] for c in candidates]]
            near = _any_within(query, reference, PORE_RADIUS_NM, box)
            for (start, _stop, _o), is_near in zip(candidates, near):
                if is_near:
                    retained += 1
                else:
                    purge_starts.append(start)

        purged = len(purge_starts)
        if purged == 0:
            if rename:
                gro.write(gro_path, rename=rename)
            self._say("No water landed in the hydrophobic core; nothing to purge")
            self._check_purge_outcome(gro_path, planes, solute, water_sites, box)
            return 0, retained

        # Guard G5b: a runaway purge means plane detection failed.
        fraction = purged / float(inserted)
        if fraction >= MAX_PURGE_FRACTION:
            raise MembraneSolvationError(
                f"The core-water purge would delete {purged} of {inserted} waters "
                f"({fraction:.1%}), far more than the ~8% a correctly detected bilayer needs. "
                f"The leaflet planes are almost certainly wrong ({planes.describe()}). "
                + _PACKMOL_HINT
            )

        keep = [True] * gro.natoms
        for start in purge_starts:
            for i in range(start, start + water_sites):
                keep[i] = False
        removed = gro.natoms - gro.write(gro_path, keep=keep, rename=rename)

        # Guard G5e: whole molecules only, and the topology must follow exactly.
        if removed != purged * water_sites:
            raise MembraneSolvationError(
                f"Purge removed {removed} atoms for {purged} waters of {water_sites} sites each; "
                "the deletion was not whole-molecule."
            )
        _set_molecule_count(top_path, WATER_MOLECULETYPE, inserted - purged)

        self._say(
            f"Purged {number(purged)} core water(s) ({fraction:.1%} of the inserted water); "
            f"kept {number(retained)} pore/interface water(s)"
        )
        if retained == 0 and reference.size:
            warnings.append(
                "No water was retained inside the transmembrane region. That is expected for a "
                "compact bundle, but a channel or transporter should keep a hydrated lumen -- "
                "inspect the pore before running production."
            )
        # Guard G5a: verify against the file that was actually written, so a
        # bug in the writer cannot pass on the strength of the in-memory state.
        self._check_purge_outcome(gro_path, planes, solute, water_sites, box)
        return purged, retained

    def _check_purge_outcome(
        self,
        gro_path: str,
        planes: LeafletPlanes,
        solute: Set[str],
        water_sites: int,
        box: Tuple[float, float, float],
    ) -> None:
        """Guard G5a: exactly zero bulk waters left in the hydrophobic core."""
        gro = _GroFile.read(gro_path)
        coords = _np.asarray(gro.coords, dtype=float)
        solute_mask = _np.array(
            [
                gro.resnames[i].upper() in solute and not _is_hydrogen_name(gro.atomnames[i])
                for i in range(gro.natoms)
            ],
            dtype=bool,
        )
        reference = coords[solute_mask]
        oxygens = [
            i
            for i in range(gro.natoms)
            if gro.resnames[i] == WATER_MOLECULETYPE
            and _leading_alpha(gro.atomnames[i]) == "O"
            and planes.core_lower_z < coords[i][2] < planes.core_upper_z
        ]
        if not oxygens:
            return
        near = _any_within(coords[oxygens], reference, PORE_RADIUS_NM, box)
        leftover = int((~near).sum())
        if leftover:
            raise MembraneSolvationError(
                f"{leftover} water molecule(s) remain in the hydrophobic core of the bilayer "
                f"(slab {planes.core_lower_z:.2f}-{planes.core_upper_z:.2f} nm) more than "
                f"{PORE_RADIUS_NM} nm from the solute. PACKMOL's reference value is zero; this "
                "system would be scientifically wrong. " + _PACKMOL_HINT
            )

    # -- ions -------------------------------------------------------------- #

    def _ionise(
        self, solv_gro: str, top_path: str, n_water: int, warnings: List[str]
    ) -> Tuple[int, int, float]:
        """Place counter-ions and salt.  Returns (n_pos, n_neg, effective_conc)."""
        view = read_topology(top_path)
        net_charge = view.net_charge()
        n_pos, n_neg = self._ion_counts(view, n_water, net_charge)
        if n_pos == 0 and n_neg == 0:
            # Only reachable for a neutral solute in a box too small to hold a
            # single formula unit at the requested molarity.  Say so rather
            # than let a silently salt-free system look deliberate.
            if self.salt_concentration > 0:
                warnings.append(
                    f"No ions placed: {n_water} water molecules cannot hold one formula unit at "
                    f"{self.salt_concentration:.3f} M. The system is salt-free."
                )
            self._say("System is neutral and no salt fits; skipping genion")
            return 0, 0, 0.0

        self._add_ions_explicit(solv_gro, top_path, n_pos, n_neg)

        final = read_topology(top_path)
        self._verify_ion_topology(final, n_water, n_pos, n_neg)

        residual = final.net_charge()
        if abs(residual) > NET_CHARGE_TOLERANCE:
            raise MembraneSolvationError(
                f"System net charge is {residual:+.3f} e after ion placement; expected ~0. "
                "grompp would warn and the Ewald neutralising background would absorb the "
                "difference. " + _PACKMOL_HINT
            )
        if abs(residual) > 0.01:
            warnings.append(
                f"Residual net charge {residual:+.4f} e comes from non-integral solute charges, "
                "not from the ion count; it cannot be removed by adding ions."
            )

        # Guard G6a: the delivered concentration, computed the way a chemist
        # would -- salt per water -- not from the total box volume the way
        # ``genion -conc`` does.
        n_water_final = final.count_of(WATER_MOLECULETYPE)
        n_salt = min(n_pos, n_neg)
        effective = WATER_MOLARITY * n_salt / n_water_final if n_water_final else 0.0
        if self.salt_concentration > 0:
            deviation = abs(effective - self.salt_concentration) / self.salt_concentration
            if deviation > SALT_CONCENTRATION_TOLERANCE:
                raise MembraneSolvationError(
                    f"Delivered salt concentration {effective:.3f} M differs from the requested "
                    f"{self.salt_concentration:.3f} M by {deviation:.0%}. "
                    f"({n_salt} salt pairs in {n_water_final} waters.)"
                )
        return n_pos, n_neg, effective

    @staticmethod
    def _round_to_integer_charge(charge: float, what: str) -> float:
        """Snap a summed charge to the whole number it physically is.

        Refuses rather than rounds when the value is not close to an integer:
        that means the topology itself is wrong (a truncated charge column, a
        molecule counted the wrong number of times), and quietly rounding it
        would hide the real fault behind a plausible ion count.
        """
        nearest = round(charge)
        if abs(charge - nearest) > CHARGE_ROUNDING_TOLERANCE:
            raise MembraneSolvationError(
                f"{what} is {charge:+.4f} e, which is not a whole number of elementary "
                f"charges (nearest is {nearest:+.0f}). A molecular system's charge must be "
                "integral, so the topology is inconsistent rather than merely imprecise. "
                + _PACKMOL_HINT
            )
        return float(nearest)

    def _ion_counts(self, view: TopologyView, n_water: int, net_charge: float) -> Tuple[int, int]:
        """Explicit ion counts from the aqueous water count and the net charge.

        ``genion -conc`` is not used: it converts a molarity into a count using
        the **total** box volume, and only ~40% of a membrane box is aqueous,
        so a requested 0.15 M is delivered as ~0.27 M.  Counting from the water
        molecules that actually exist is both correct and independent of the
        bilayer's share of the box.
        """
        q_pos = self._ion_charge(view, self.positive_ion, expect_positive=True)
        q_neg = self._ion_charge(view, self.negative_ion, expect_positive=False)

        # Stoichiometry of one neutral formula unit (NaCl -> 1:1, CaCl2 -> 1:2).
        gcd = math.gcd(int(round(abs(q_pos))), int(round(abs(q_neg)))) or 1
        stoich_pos = int(round(abs(q_neg))) // gcd
        stoich_neg = int(round(abs(q_pos))) // gcd

        units = int(round(self.salt_concentration * n_water / WATER_MOLARITY))
        n_pos = units * stoich_pos
        n_neg = units * stoich_neg

        # Neutralise on top of the salt, as genion -neutral would, but with the
        # solute charge read from the topology rather than from a tpr.
        #
        # Round the residual to the integer it physically is before deciding how
        # many ions to add.  Summing ~10^5 per-atom charges printed to 6 decimal
        # places accumulates a few thousandths, and ceil() promotes ANY positive
        # excess to a whole extra ion: a +10 solute read as +10.002 was given 11
        # chloride instead of 10, leaving the system at -1 e.  A 1e-9 epsilon
        # cannot absorb that -- the error is a million times larger -- and the
        # symptom is silent, because a system off by one elementary charge still
        # builds and still runs, with the Ewald background quietly soaking up
        # the difference.
        residual = net_charge + n_pos * q_pos + n_neg * q_neg
        residual = self._round_to_integer_charge(residual, "Residual system charge")
        if residual > 0:
            n_neg += int(math.ceil(residual / abs(q_neg) - CHARGE_ROUNDING_TOLERANCE))
        elif residual < 0:
            n_pos += int(math.ceil(-residual / abs(q_pos) - CHARGE_ROUNDING_TOLERANCE))

        if n_pos + n_neg > n_water:
            raise MembraneSolvationError(
                f"Ion placement would need {n_pos + n_neg} of only {n_water} water molecules."
            )
        self._say(
            f"Ion budget: solute net charge {net_charge:+.2f} e, {units} salt formula unit(s) "
            f"for {n_water} waters at {self.salt_concentration:.3f} M -> "
            f"{n_pos} {self.positive_ion} + {n_neg} {self.negative_ion}"
        )
        return n_pos, n_neg

    @staticmethod
    def _ion_charge(view: TopologyView, name: str, expect_positive: bool) -> float:
        mtype = view.moleculetypes.get(name)
        if mtype is None:  # already guarded, but this is arithmetic -- be explicit
            raise MembraneSolvationError(f"No [ moleculetype ] named {name!r} in the topology.")
        if not mtype.charge_known:
            raise MembraneSolvationError(f"Moleculetype {name!r} has no explicit charge column.")
        if mtype.natoms != 1:
            raise MembraneSolvationError(
                f"Ion moleculetype {name!r} has {mtype.natoms} atoms; only monatomic ions are "
                "supported by the fast path. " + _PACKMOL_HINT
            )
        charge = mtype.charge
        if expect_positive and charge <= 0:
            raise MembraneSolvationError(
                f"Positive ion {name!r} carries charge {charge:+.3f} in the topology. This is the "
                "signature of an atomtype collision (e.g. sodium written with the aromatic "
                "nitrogen type 'NA'), which grompp accepts without a single diagnostic. "
                + _PACKMOL_HINT
            )
        if not expect_positive and charge >= 0:
            raise MembraneSolvationError(
                f"Negative ion {name!r} carries charge {charge:+.3f} in the topology. "
                + _PACKMOL_HINT
            )
        if abs(charge - round(charge)) > 1e-3:
            raise MembraneSolvationError(
                f"Ion {name!r} has a non-integral charge {charge:+.4f}. " + _PACKMOL_HINT
            )
        return float(round(charge))

    def _add_ions_explicit(self, solv_gro: str, top_path: str, n_pos: int, n_neg: int) -> str:
        """grompp + genion with explicit ``-np``/``-nn``.

        ``genion -p`` rewrites the topology *before* it places anything and
        leaves it mutated if it then dies, with no backup of its own -- so a
        failed run would otherwise poison the topology for every retry.  The
        backup is taken first and restored on any failure.

        ``-rmin`` is deliberately never passed: the default 0.6 nm is enforced
        against every non-solvent atom, which is exactly what stops genion
        burying an ion in the acyl core, and overriding it removes that
        protection.
        """
        model_dir = str(self.model_dir)
        ions_gro = os.path.join(model_dir, "solv_ions.gro")
        ions_mdp = os.path.join(model_dir, "ions.mdp")
        ions_tpr = os.path.join(model_dir, "ions.tpr")
        backup = top_path + ".pre_genion"
        shutil.copyfile(top_path, backup)

        with open(ions_mdp, "w") as fh:
            fh.write("integrator = steep\nnsteps = 0\n")

        try:
            stdout, stderr = self._run_command(
                [
                    self.gmx_command, "grompp",
                    "-f", ions_mdp,
                    "-c", solv_gro,
                    "-p", top_path,
                    "-o", ions_tpr,
                    "-maxwarn", "5",
                ],
                model_dir,
            )
            # Guard G7 on the pre-genion system: there are no ions yet, so the
            # only legitimate mismatch count is zero.  A scrambled system
            # grompps to a valid tpr behind this single warning, which -maxwarn
            # then swallows -- so it has to be read, not merely allowed.
            self._check_atom_name_matching(stdout, stderr, allowed=0)

            self._run_command(
                [
                    self.gmx_command, "genion",
                    "-s", ions_tpr,
                    "-o", ions_gro,
                    "-p", top_path,
                    "-pname", self.positive_ion,
                    "-nname", self.negative_ion,
                    "-np", str(n_pos),
                    "-nn", str(n_neg),
                ],
                model_dir,
                input_str=f"{WATER_MOLECULETYPE}\n",
            )
        except Exception:
            shutil.copyfile(backup, top_path)
            raise
        finally:
            for temp in (ions_tpr, ions_mdp, os.path.join(model_dir, "mdout.mdp")):
                if os.path.exists(temp):
                    os.remove(temp)

        if not os.path.isfile(ions_gro):
            shutil.copyfile(backup, top_path)
            raise MembraneSolvationError("gmx genion did not produce solv_ions.gro.")
        self._align_ion_names(ions_gro, top_path, n_pos, n_neg)
        return ions_gro

    @staticmethod
    def _check_atom_name_matching(stdout: str, stderr: str, allowed: int) -> None:
        """Guard G7: topology and coordinates must describe the same system.

        A fully scrambled system produces a valid ``.tpr`` of the same size
        behind exactly one warning; only its count distinguishes it from a
        correct build.
        """
        pattern = re.compile(r"(\d+)\s+non-matching atom names")
        for match in pattern.finditer(f"{stdout}\n{stderr}"):
            count = int(match.group(1))
            if count > allowed:
                raise MembraneSolvationError(
                    f"grompp reported {count} non-matching atom names between the topology and "
                    f"the coordinates (at most {allowed} expected). The two files do not describe "
                    "the same system; grompp would still write a usable-looking .tpr. "
                    + _PACKMOL_HINT
                )

    def _align_ion_names(self, gro_path: str, top_path: str, n_pos: int, n_neg: int) -> None:
        """Give the placed ions the names the topology uses (cosmetic, best effort).

        ``genion`` names the ion residues after ``-pname``/``-nname``, which is
        not always byte-identical to the atom name inside the moleculetype.
        Aligning them keeps the *final* grompp's mismatch count at zero, so the
        builder's preflight can apply guard G7 strictly rather than having to
        tolerate ``N == n_ions``.  Never fatal: it changes nothing physical.
        """
        try:
            view = read_topology(top_path)
            gro = _GroFile.read(gro_path)
            targets = {}
            for name, expected in ((self.positive_ion, n_pos), (self.negative_ion, n_neg)):
                mtype = view.moleculetypes.get(name)
                if mtype is None or mtype.natoms != 1 or expected <= 0:
                    continue
                key = name.upper().rstrip("+-")
                targets[key] = (
                    mtype.residue_names[0] if mtype.residue_names else name,
                    mtype.atom_names[0] if mtype.atom_names else name,
                    expected,
                )
            rename: Dict[int, Tuple[str, str]] = {}
            seen: Dict[str, int] = {key: 0 for key in targets}
            for start, stop in gro.residue_slices():
                if stop - start != 1:
                    continue
                key = gro.resnames[start].upper().rstrip("+-")
                if key not in targets:
                    continue
                resname, atomname, _ = targets[key]
                seen[key] += 1
                if gro.resnames[start] != resname or gro.atomnames[start] != atomname:
                    rename[start] = (resname, atomname)
            # Only rewrite when every ion was accounted for; a count mismatch
            # means the identification was wrong and the file is left alone.
            if rename and all(seen[key] == targets[key][2] for key in targets):
                gro.write(gro_path, rename=rename)
        except Exception as exc:  # cosmetic only -- never fail the build for this
            self._say(f"Skipped ion atom-name alignment ({exc})")

    def _verify_ion_topology(
        self, view: TopologyView, n_water_before: int, n_pos: int, n_neg: int
    ) -> None:
        """Check genion's own topology edit instead of trusting it.

        The soluble path's verification falls back to searching the file for
        the substring ``'CA '``, which matches alpha-carbon atom names hundreds
        of times and therefore reports success unconditionally.  These are the
        counts genion must have written, and nothing else will do.
        """
        expected = {
            WATER_MOLECULETYPE: n_water_before - n_pos - n_neg,
            self.positive_ion: n_pos,
            self.negative_ion: n_neg,
        }
        for name, want in expected.items():
            got = view.count_of(name)
            if got != want:
                raise MembraneSolvationError(
                    f"After genion the topology lists {got} x {name}, expected {want}. "
                    f"The [ molecules ] section no longer matches the coordinates; the backup is "
                    f"at {os.path.basename(view.path)}.pre_genion. " + _PACKMOL_HINT
                )

    def _assert_no_core_ions(
        self,
        gro_path: str,
        top_path: str,
        planes: LeafletPlanes,
        solute_resnames: Set[str],
        n_pos: int,
        n_neg: int,
        warnings: List[str],
    ) -> None:
        """Guard G6d: no ion inside the hydrophobic core.

        Expected always to pass -- genion enforces ``-rmin 0.6`` against every
        non-solvent atom -- so a failure means that default was overridden
        somewhere, which is worth knowing loudly.  Also checks the coordinate
        file against the topology's atom total, which is the cheapest detector
        of a corrupted genion topology edit.
        """
        if not os.path.isfile(gro_path):
            return
        gro = _GroFile.read(gro_path)
        expected_atoms = read_topology(top_path).natoms()
        if gro.natoms != expected_atoms:
            raise MembraneSolvationError(
                f"{os.path.basename(gro_path)} holds {gro.natoms} atoms but the topology accounts "
                f"for {expected_atoms}. " + _PACKMOL_HINT
            )
        if n_pos == 0 and n_neg == 0:
            return
        coords = _np.asarray(gro.coords, dtype=float)
        keys = {self.positive_ion.upper().rstrip("+-"), self.negative_ion.upper().rstrip("+-")}
        buried = [
            i
            for i in range(gro.natoms)
            if gro.resnames[i].upper().rstrip("+-") in keys
            and planes.core_lower_z < coords[i][2] < planes.core_upper_z
        ]
        if buried:
            raise MembraneSolvationError(
                f"{len(buried)} ion(s) were placed inside the hydrophobic core of the bilayer. "
                "genion's default -rmin 0.6 should make this impossible. " + _PACKMOL_HINT
            )


def view_water_count(top_path: str) -> int:
    """Water molecules currently listed in a topology's ``[ molecules ]``."""
    return read_topology(top_path).count_of(WATER_MOLECULETYPE)


# --------------------------------------------------------------------------- #
# Neighbour search
# --------------------------------------------------------------------------- #


def _any_within(query, reference, cutoff: float, box: Tuple[float, float, float]):
    """Boolean mask: does each *query* point have a *reference* point within *cutoff*?

    Periodic in all three axes.  x/y periodicity is the one that matters -- a
    protein touching the box edge would otherwise look far from a water in its
    own periodic image and that water would be purged out of the pore -- and z
    periodicity is harmless because the two water slabs are contiguous across
    the boundary anyway.
    """
    n_query = len(query)
    if n_query == 0:
        return _np.zeros(0, dtype=bool)
    if len(reference) == 0:
        return _np.zeros(n_query, dtype=bool)

    box_arr = _np.asarray(box, dtype=float)
    # cKDTree's periodic mode requires every coordinate inside [0, L); wrapping
    # a copy keeps the guarantee without touching what gets written out.
    q = _np.mod(_np.asarray(query, dtype=float), box_arr)
    r = _np.mod(_np.asarray(reference, dtype=float), box_arr)
    try:
        from scipy.spatial import cKDTree

        tree = cKDTree(r, boxsize=box_arr)
        distances, _ = tree.query(q, k=1, distance_upper_bound=cutoff)
        return _np.isfinite(distances)
    except ImportError:
        pass

    # Fallback: chunked minimum-image brute force.  The reference set is
    # already restricted to the slab, so this is milliseconds, not minutes.
    result = _np.zeros(n_query, dtype=bool)
    chunk = max(1, int(4_000_000 // max(1, len(r))))
    for start in range(0, n_query, chunk):
        stop = min(start + chunk, n_query)
        delta = q[start:stop, None, :] - r[None, :, :]
        delta -= box_arr * _np.round(delta / box_arr)
        result[start:stop] = (_np.einsum("ijk,ijk->ij", delta, delta) < cutoff * cutoff).any(axis=1)
    return result
