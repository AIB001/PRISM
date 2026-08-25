"""Fast, deterministic protein-in-bilayer assembly from pre-equilibrated patches.

Motivation
----------
PRISM's original membrane route packs lipids from scratch with PACKMOL-Memgen.
That is unusable at scale: at the physically correct area-per-lipid of ~65 A^2 a
rigid POPC template's 95th-percentile radius (4.54 A) already equals the radius
of a 65 A^2 disk (4.55 A), so there is *zero* clearance for PACKMOL's 2.0 A
tolerance.  PACKMOL therefore never converges (it needs APL >= 90-95 A^2 to
satisfy its own tolerance), burns its whole iteration budget, is strictly
single-threaded, and emits rigid clones.  Measured: 17-30 min for a dry
protein + 650 lipid pack, and >2.9 h for a wet system.

Rigid-body geometric placement provably cannot reach lipid liquid density on its
own, which is exactly why CHARMM-GUI and VMD ship *pre-equilibrated* patches.
This module takes the same route: it tiles a bundled, genuinely equilibrated
bilayer patch, drops the oriented solute in without moving it, and deletes the
lipids that overlap it.  There is no optimiser and no randomness anywhere, so
the whole geometry step costs a few seconds and is bit-for-bit reproducible.

Bundled asset
-------------
Equilibrated AMBER Lipid21 bilayers from:

    Frankel, Patrick (2025). *Amber Lipid21 bilayer simulations.*
    Zenodo. https://doi.org/10.5281/zenodo.14776136
    Licensed CC-BY-4.0.

These are reproductions of the validation systems of the Lipid21 paper
(Dickson, Walker & Gould, *J. Chem. Theory Comput.* 2022, 18, 1726-1736,
doi:10.1021/acs.jctc.1c01217) -- 128 lipids + 5120 TIP3P waters (+128 K+ for
anionic lipids) each.  Residue names are already in native Lipid21 form
(PA / OL / PC / PE / PS / MY / SPM / SA / CHL) so the emitted coordinates load
directly under ``leaprc.lipid21``.  See ``prism/data/membrane_patches/NOTICE``
for the attribution the licence requires.

Ordering contract
-----------------
GROMACS requires every molecule's atoms to be contiguous, and ``gmx genconf
-nbox`` block-replicates, so naively multiplying ``[molecules]`` counts against
a tiled coordinate file is wrong.  Everything this module emits therefore obeys
one fixed contract, and :meth:`MembraneComposition.species` describes it
exactly:

1. **The solute comes first**, verbatim, in its input order, with its
   coordinates untouched.
2. **Then the lipids, grouped by species.**  Species appear in the order given
   by ``MembraneComposition.species``; that list is what a ``[molecules]``
   section must mirror.
3. **Within a species**, molecules are ordered lower leaflet first, then upper
   leaflet; within a leaflet by tile (row-major in x then y); within a tile by
   the molecule's index in the source patch.
4. **Within a molecule**, atoms keep the patch's original order, contiguous,
   with the Lipid21 residues in their canonical sequence.  A ``TER`` follows
   every molecule in PDB output.

No water and no ions are emitted -- this builds the *dry* system, and
:attr:`MembraneComposition.formal_charge` reports the charge the caller must
neutralise once it solvates.

Shrink-and-grow
---------------
Deleting the lipids that overlap the solute gets the *count* right but leaves
the survivors sitting on the lattice positions they occupied in a protein-free
bilayer, so the bilayer's inner boundary is ragged and stands ~10 A off the
protein: a vacuum annulus at exactly the place the first solvation shell of
lipid should be.  :func:`delete_clashing_lipids` therefore also implements the
algorithm of Wolf et al., *J. Comp. Chem.* **31**, 2169-2174 (2010) -- the one
``openmm.app.Modeller.addMembrane`` uses -- selected with ``shrink_factor``:

1. scale the solute down in x and y about the centre of its bilayer-spanning
   region (:func:`scale_solute`);
2. delete the lipids that overlap the *scaled* solute, so the survivors are
   left *inside* the real footprint rather than outside it;
3. remove the same number from each leaflet;
4. grow the solute back to full size in small increments, relaxing the bilayer
   at each one (:func:`grow_solute`, :func:`gromacs_relaxer`), which pushes the
   survivors outward into real van der Waals contact.

Steps 1-3 are pure geometry and cost milliseconds; step 4 needs a force field
and is driven through a caller-supplied ``relax`` callable so this module keeps
its "no optimiser, no randomness" character and stays importable without
GROMACS.

Public API
----------
``available_patches()``, ``patch_citation()``, ``load_patch()``,
``patch_geometry()``, ``tile_patch()``, ``required_tiling()``,
``read_solute()``, ``insert_solute()``, ``delete_clashing_lipids()``,
``build_dry_membrane_system()``, ``scale_solute()``,
``transmembrane_mask()``, ``transmembrane_cross_section()``,
``area_budget_shrink_factor()``, ``grow_solute()``, ``gromacs_relaxer()``,
and the dataclasses ``PatchSpec``, ``PatchGeometry``, ``LipidAssembly``,
``Solute``, ``DrySystem``, ``DeletionReport``, ``GrowthPlan``,
``GrowthReport``, ``MembraneComposition``.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "MembranePatchError",
    "PatchNotAvailableError",
    "PatchIntegrityError",
    "PatchSpec",
    "PatchGeometry",
    "LipidAssembly",
    "Solute",
    "DeletionReport",
    "GrowthPlan",
    "GrowthReport",
    "MembraneComposition",
    "DrySystem",
    "available_patches",
    "patch_citation",
    "patch_data_dir",
    "load_patch",
    "patch_geometry",
    "required_tiling",
    "tile_patch",
    "read_solute",
    "insert_solute",
    "delete_clashing_lipids",
    "build_dry_membrane_system",
    "transmembrane_mask",
    "transmembrane_cross_section",
    "scale_solute",
    "area_budget_shrink_factor",
    "grow_solute",
    "gromacs_relaxer",
    "LIPID21_RESIDUE_PATTERNS",
    "DEFAULT_CLASH_CUTOFF",
    "OPENMM_SHRINK_FACTOR",
    "SHRINK_FACTOR_BOUNDS",
    "DEFAULT_GROWTH_STEP",
]


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class MembranePatchError(RuntimeError):
    """Base class for every failure raised by this module."""


class PatchNotAvailableError(MembranePatchError):
    """Raised when no bundled patch exists for the requested lipid."""


class PatchIntegrityError(MembranePatchError):
    """Raised when a bundled patch fails its checksum or structural checks."""


# ---------------------------------------------------------------------------
# Lipid21 chemistry tables
#
# A Lipid21 lipid is split across several residues (head + tails).  A *molecule*
# is therefore a fixed sequence of residue names, and that sequence is what
# defines molecule boundaries in a coordinate-only file.
# ---------------------------------------------------------------------------

LIPID21_RESIDUE_PATTERNS: Dict[str, Tuple[str, ...]] = {
    "POPC": ("PA", "PC", "OL"),     # sn-1 palmitoyl, phosphatidylcholine, sn-2 oleoyl
    "POPE": ("PA", "PE", "OL"),
    "POPS": ("PA", "PS", "OL"),
    "DMPC": ("MY", "PC", "MY"),     # two myristoyl tails
    "PSM": ("PA", "SPM", "SA"),     # N-palmitoyl sphingomyelin
    "CHL": ("CHL",),                # cholesterol: a single residue, no phosphorus
}

#: Formal charge of each whole molecule, used to tell the caller how much
#: neutralising counter-charge the dry membrane still needs.
LIPID21_FORMAL_CHARGE: Dict[str, int] = {
    "POPC": 0,
    "POPE": 0,
    "POPS": -1,
    "DMPC": 0,
    "PSM": 0,
    "CHL": 0,
}

#: Atom name of the phosphorus used for leaflet identification.  Requirement:
#: leaflet assignment must come from lipid phosphorus, never from a guess.
PHOSPHORUS_ATOM_NAMES = frozenset({"P31", "P", "P1"})

#: Molecules with no phosphorus of their own.  They can ride along inside a
#: patch whose midplane is defined by a phosphorus-bearing partner, but a patch
#: made *only* of these cannot have its leaflets identified and is refused.
STEROL_SPECIES = frozenset({"CHL"})

#: Default clash criterion: delete a lipid if ANY of its heavy atoms is within
#: this distance of ANY solute heavy atom.  See ``docs`` in
#: :func:`delete_clashing_lipids` for the sweep that produced it.
DEFAULT_CLASH_CUTOFF = 2.6

#: The flat shrink factor used by ``openmm.app.Modeller.addMembrane`` and by
#: Wolf et al. (2010).  Exposed so a caller can reproduce them exactly, but it
#: is deliberately *not* the default here; see
#: :func:`area_budget_shrink_factor` for why a flat value misprices a protein
#: with a large transmembrane cross-section.
OPENMM_SHRINK_FACTOR = 0.5

#: Range the solved shrink factor is clamped to.  Below the lower bound the
#: survivors end up buried so deep in the real solute that the grow-back has to
#: drag them through it; above the upper bound the shrink does nothing that
#: plain deletion did not already do.
SHRINK_FACTOR_BOUNDS = (0.45, 0.95)

#: Largest lateral displacement, in Angstrom, that any bilayer-spanning solute
#: atom may make in a single grow-back increment.  See :class:`GrowthPlan`.
DEFAULT_GROWTH_STEP = 1.0

#: Effective heavy-atom radius used when rasterising the solute's cross-section.
#: Half of a C-C van der Waals contact; the absolute value only has to be
#: consistent with the ``target_apl`` the budget is solved against.
CROSS_SECTION_PROBE = 1.9

_MANIFEST_NAME = "manifest.json"
_PATCH_SUBDIR = os.path.join("data", "membrane_patches")

#: Environment variable pointing at a directory of extra/override patches.
PATCH_DIR_ENV = "PRISM_MEMBRANE_PATCH_DIR"


# ---------------------------------------------------------------------------
# Asset discovery
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PatchSpec:
    """Static description of one bundled patch."""

    key: str
    filename: str
    sha256: str
    n_atoms: int
    n_molecules: int
    composition: Tuple[Tuple[str, int], ...]
    box: Tuple[float, float, float]
    temperature_k: float
    source_file: str
    source_md5: str
    reference_apl: Optional[float] = None
    reference_note: str = ""
    aliases: Tuple[str, ...] = ()

    @property
    def lipid_species(self) -> Tuple[str, ...]:
        return tuple(name for name, _ in self.composition)


def patch_data_dir() -> Path:
    """Return the directory holding the bundled patches.

    Honours :data:`PATCH_DIR_ENV` so a user can point PRISM at extra patches
    (or newer ones) without reinstalling.  Falls back to the copy shipped
    inside the ``prism`` package.
    """
    override = os.environ.get(PATCH_DIR_ENV)
    if override:
        path = Path(override).expanduser()
        if not path.is_dir():
            raise PatchNotAvailableError(
                f"{PATCH_DIR_ENV} is set to {path!s} but that is not a directory."
            )
        return path
    return Path(__file__).resolve().parent.parent / _PATCH_SUBDIR


def _load_manifest() -> Tuple[Dict[str, PatchSpec], Dict[str, str], Dict[str, str]]:
    """Parse ``manifest.json``; return (specs, alias->key, provenance)."""
    directory = patch_data_dir()
    manifest_path = directory / _MANIFEST_NAME
    if not manifest_path.is_file():
        raise PatchNotAvailableError(
            "No membrane patch manifest found at {0}.\n"
            "The equilibrated bilayer patches ship inside the PRISM package; a "
            "missing manifest usually means PRISM was installed as a wheel built "
            "without package data.  Reinstall from source with 'pip install -e .', "
            "or set {1} to a directory containing manifest.json.".format(
                manifest_path, PATCH_DIR_ENV
            )
        )
    with open(manifest_path, "r", encoding="utf-8") as handle:
        raw = json.load(handle)

    specs: Dict[str, PatchSpec] = {}
    aliases: Dict[str, str] = {}
    for key, entry in sorted(raw.get("patches", {}).items()):
        spec = PatchSpec(
            key=key,
            filename=entry["filename"],
            sha256=entry["sha256"],
            n_atoms=int(entry["n_atoms"]),
            n_molecules=int(entry["n_molecules"]),
            composition=tuple((str(n), int(c)) for n, c in entry["composition"]),
            box=tuple(float(v) for v in entry["box"]),
            temperature_k=float(entry["temperature_k"]),
            source_file=entry["source_file"],
            source_md5=entry["source_md5"],
            reference_apl=(
                float(entry["reference_apl"])
                if entry.get("reference_apl") is not None
                else None
            ),
            reference_note=entry.get("reference_note", ""),
            aliases=tuple(entry.get("aliases", ())),
        )
        specs[key] = spec
        aliases[key.upper()] = key
        for alias in spec.aliases:
            aliases[alias.upper()] = key
    provenance = {
        k: v for k, v in raw.items() if k != "patches"
    }
    return specs, aliases, provenance


def available_patches() -> Dict[str, PatchSpec]:
    """Return every bundled patch, keyed by canonical name.

    Cheap: reads a small JSON manifest, never a coordinate file.
    """
    specs, _, _ = _load_manifest()
    return specs


def patch_citation() -> str:
    """Return the attribution text the CC-BY-4.0 licence requires."""
    _, _, provenance = _load_manifest()
    return provenance.get("attribution", "").strip()


def _resolve_key(lipid: str) -> Tuple[PatchSpec, Dict[str, PatchSpec]]:
    specs, aliases, _ = _load_manifest()
    key = aliases.get(str(lipid).strip().upper())
    if key is None:
        known = ", ".join(sorted(specs))
        alias_hint = ", ".join(
            sorted({a for s in specs.values() for a in s.aliases})
        )
        raise PatchNotAvailableError(
            "No pre-equilibrated membrane patch is bundled for lipid "
            f"'{lipid}'.\n"
            f"  Available patches: {known}\n"
            f"  Accepted aliases : {alias_hint}\n"
            "Options:\n"
            "  * pick one of the lipids above (POPC is the default and the best "
            "tested);\n"
            f"  * point {PATCH_DIR_ENV} at a directory with your own patch plus a "
            "manifest.json entry;\n"
            "  * or fall back to the PACKMOL-Memgen route, which can build "
            "arbitrary compositions but takes tens of minutes instead of seconds."
        )
    return specs[key], specs


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# Geometry containers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PatchGeometry:
    """Measured geometry of a patch or tiled patch.  All lengths in Angstrom."""

    box: Tuple[float, float, float]
    area: float
    n_upper: int
    n_lower: int
    area_per_lipid: float
    phosphate_thickness: float          # mean P-to-P separation of the leaflets
    core_half_width: float              # |z| below which no water was present
    species_counts: Dict[str, int]
    n_atoms: int
    n_molecules: int

    def describe(self) -> str:
        comp = ", ".join(f"{k}x{v}" for k, v in sorted(self.species_counts.items()))
        return (
            f"box {self.box[0]:.2f} x {self.box[1]:.2f} x {self.box[2]:.2f} A; "
            f"{self.n_molecules} molecules ({comp}); "
            f"{self.n_lower}/{self.n_upper} lipids lower/upper leaflet; "
            f"APL {self.area_per_lipid:.2f} A^2; "
            f"P-P {self.phosphate_thickness:.2f} A"
        )


@dataclass
class LipidAssembly:
    """A bilayer: flat per-atom arrays plus contiguous molecule offsets.

    Molecules are guaranteed contiguous: molecule ``m`` owns atoms
    ``mol_offsets[m]:mol_offsets[m + 1]``.
    """

    atom_names: np.ndarray      # (N,)   <U5
    res_names: np.ndarray       # (N,)   <U5
    res_local: np.ndarray       # (N,)   int32, residue index within the molecule
    elements: np.ndarray        # (N,)   <U2
    xyz: np.ndarray             # (N, 3) float64, Angstrom
    mol_offsets: np.ndarray     # (M+1,) int64
    mol_species: np.ndarray     # (M,)   <U8
    mol_leaflet: np.ndarray     # (M,)   int8, -1 lower / +1 upper
    box: np.ndarray             # (3,)   float64, Angstrom (orthorhombic)
    key: str = ""
    n_tiles: Tuple[int, int] = (1, 1)
    core_half_width: float = 0.0

    # -- basic accessors ---------------------------------------------------
    @property
    def n_atoms(self) -> int:
        return int(self.atom_names.shape[0])

    @property
    def n_molecules(self) -> int:
        return int(self.mol_species.shape[0])

    @property
    def heavy_mask(self) -> np.ndarray:
        return self.elements != "H"

    @property
    def atom_molecule(self) -> np.ndarray:
        """Per-atom molecule index."""
        counts = np.diff(self.mol_offsets)
        return np.repeat(np.arange(self.n_molecules, dtype=np.int64), counts)

    def species_counts(self) -> Dict[str, int]:
        names, counts = np.unique(self.mol_species, return_counts=True)
        return {str(n): int(c) for n, c in zip(names, counts)}

    def leaflet_counts(self) -> Dict[str, Tuple[int, int]]:
        """species -> (n_lower, n_upper)."""
        out: Dict[str, Tuple[int, int]] = {}
        for name in sorted(set(self.mol_species.tolist())):
            sel = self.mol_species == name
            out[str(name)] = (
                int(np.count_nonzero(sel & (self.mol_leaflet < 0))),
                int(np.count_nonzero(sel & (self.mol_leaflet > 0))),
            )
        return out

    def select_molecules(self, keep: np.ndarray) -> "LipidAssembly":
        """Return a new assembly containing only ``keep`` molecules.

        Molecule order is preserved, so the ordering contract survives.
        """
        keep = np.asarray(keep, dtype=bool)
        if keep.shape != (self.n_molecules,):
            raise ValueError("keep mask must have one entry per molecule")
        counts = np.diff(self.mol_offsets)
        atom_keep = np.repeat(keep, counts)
        new_counts = counts[keep]
        new_offsets = np.zeros(int(keep.sum()) + 1, dtype=np.int64)
        np.cumsum(new_counts, out=new_offsets[1:])
        return LipidAssembly(
            atom_names=self.atom_names[atom_keep],
            res_names=self.res_names[atom_keep],
            res_local=self.res_local[atom_keep],
            elements=self.elements[atom_keep],
            xyz=self.xyz[atom_keep],
            mol_offsets=new_offsets,
            mol_species=self.mol_species[keep],
            mol_leaflet=self.mol_leaflet[keep],
            box=self.box.copy(),
            key=self.key,
            n_tiles=self.n_tiles,
            core_half_width=self.core_half_width,
        )

    def reorder_molecules(self, order: np.ndarray) -> "LipidAssembly":
        """Return a new assembly with molecules permuted by ``order``."""
        order = np.asarray(order, dtype=np.int64)
        counts = np.diff(self.mol_offsets)
        starts = self.mol_offsets[:-1]
        atom_index = np.concatenate(
            [np.arange(starts[m], starts[m] + counts[m], dtype=np.int64) for m in order]
        ) if order.size else np.zeros(0, dtype=np.int64)
        new_counts = counts[order]
        new_offsets = np.zeros(order.size + 1, dtype=np.int64)
        np.cumsum(new_counts, out=new_offsets[1:])
        return LipidAssembly(
            atom_names=self.atom_names[atom_index],
            res_names=self.res_names[atom_index],
            res_local=self.res_local[atom_index],
            elements=self.elements[atom_index],
            xyz=self.xyz[atom_index],
            mol_offsets=new_offsets,
            mol_species=self.mol_species[order],
            mol_leaflet=self.mol_leaflet[order],
            box=self.box.copy(),
            key=self.key,
            n_tiles=self.n_tiles,
            core_half_width=self.core_half_width,
        )

    def molecule_centroids(self) -> np.ndarray:
        """(M, 3) centroid of every molecule (unweighted)."""
        counts = np.diff(self.mol_offsets)
        sums = np.add.reduceat(self.xyz, self.mol_offsets[:-1], axis=0)
        return sums / counts[:, None]


# ---------------------------------------------------------------------------
# Elements
# ---------------------------------------------------------------------------

def _elements_from_names(names: np.ndarray) -> np.ndarray:
    """Infer elements from Lipid21 atom names.

    Verified against every residue of every bundled patch: each hydrogen name
    begins with 'H' and every heavy-atom name begins with C, N, O, P, S or K,
    so the first-character rule is unambiguous here.  Anything else raises
    rather than being silently mis-typed as a heavy atom.
    """
    out = np.empty(names.shape, dtype="<U2")
    for i, raw in enumerate(names):
        name = str(raw).strip()
        if not name:
            raise PatchIntegrityError("empty atom name in patch")
        first = name[0]
        if first.isdigit():             # PDB-style '1HB'
            if len(name) < 2:
                raise PatchIntegrityError(f"uninterpretable atom name {name!r}")
            first = name[1]
        if first == "H":
            out[i] = "H"
        elif first in "CNOPSK":
            out[i] = first
        else:
            raise PatchIntegrityError(
                f"cannot infer element for atom name {name!r}; refusing to guess"
            )
    return out


def _elements_from_pdb(elem_field: Sequence[str], names: Sequence[str]) -> np.ndarray:
    """Prefer the PDB element column, fall back to the atom name."""
    out = np.empty(len(names), dtype="<U2")
    for i, (elem, name) in enumerate(zip(elem_field, names)):
        elem = (elem or "").strip()
        if elem:
            out[i] = elem[0].upper() if len(elem) == 1 else (
                elem[0].upper() + elem[1].lower()
            )
            continue
        text = str(name).strip()
        if not text:
            raise MembranePatchError("empty atom name in solute")
        first = text[0]
        if first.isdigit():
            first = text[1] if len(text) > 1 else "X"
        out[i] = first.upper()
    return out


# ---------------------------------------------------------------------------
# GRO reading / writing
# ---------------------------------------------------------------------------

def _read_gro(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read a (optionally gzipped) GRO file.

    Returns (res_ids, res_names, atom_names, xyz_angstrom, box_angstrom).
    Velocities, if present, are ignored.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    if len(lines) < 3:
        raise PatchIntegrityError(f"{path} is too short to be a GRO file")
    try:
        n_atoms = int(lines[1])
    except ValueError as exc:
        raise PatchIntegrityError(f"{path}: bad atom count line") from exc
    if len(lines) < n_atoms + 3:
        raise PatchIntegrityError(
            f"{path}: declares {n_atoms} atoms but has only {len(lines) - 3}"
        )

    res_ids = np.empty(n_atoms, dtype=np.int64)
    res_names = np.empty(n_atoms, dtype="<U5")
    atom_names = np.empty(n_atoms, dtype="<U5")
    xyz = np.empty((n_atoms, 3), dtype=np.float64)
    for i in range(n_atoms):
        line = lines[2 + i]
        res_ids[i] = int(line[0:5])
        res_names[i] = line[5:10].strip()
        atom_names[i] = line[10:15].strip()
        xyz[i, 0] = float(line[20:28])
        xyz[i, 1] = float(line[28:36])
        xyz[i, 2] = float(line[36:44])

    fields = lines[2 + n_atoms].split()
    if len(fields) < 3:
        raise PatchIntegrityError(f"{path}: bad box line")
    box = np.array([float(fields[0]), float(fields[1]), float(fields[2])])
    if len(fields) >= 9 and any(abs(float(f)) > 1e-6 for f in fields[3:9]):
        raise PatchIntegrityError(
            f"{path}: box is triclinic; only orthorhombic patches are supported"
        )
    return res_ids, res_names, atom_names, xyz * 10.0, box * 10.0


def _gro_lines(
    atom_names: np.ndarray,
    res_names: np.ndarray,
    res_seq: np.ndarray,
    xyz: np.ndarray,
    box: np.ndarray,
    title: str,
) -> List[str]:
    """Format GRO records.  Counters wrap modulo their field width, which is
    what GROMACS itself does; molecule identity is carried by *order*, not by
    these numbers."""
    out = [title, f"{xyz.shape[0]:5d}"]
    nm = xyz / 10.0
    for i in range(xyz.shape[0]):
        out.append(
            "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
            % (
                int(res_seq[i]) % 100000,
                str(res_names[i])[:5],
                str(atom_names[i])[:5],
                (i + 1) % 100000,
                nm[i, 0],
                nm[i, 1],
                nm[i, 2],
            )
        )
    out.append(
        "%10.5f%10.5f%10.5f" % (box[0] / 10.0, box[1] / 10.0, box[2] / 10.0)
    )
    return out


# ---------------------------------------------------------------------------
# Molecule segmentation
# ---------------------------------------------------------------------------

def _residue_boundaries(res_ids: np.ndarray, res_names: np.ndarray) -> np.ndarray:
    """Indices where a new residue starts, plus a terminating index."""
    n = res_ids.shape[0]
    if n == 0:
        return np.zeros(1, dtype=np.int64)
    changed = (res_ids[1:] != res_ids[:-1]) | (res_names[1:] != res_names[:-1])
    starts = np.concatenate(
        ([0], np.nonzero(changed)[0] + 1, [n])
    ).astype(np.int64)
    return starts


def _segment_molecules(
    res_ids: np.ndarray,
    res_names: np.ndarray,
    candidates: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split a residue stream into whole molecules by Lipid21 residue pattern.

    Longest pattern is tried first so segmentation is unambiguous, and the
    stream must be consumed exactly -- a leftover residue is an error rather
    than a silently dropped atom.

    Returns (mol_offsets, mol_species, res_local).
    """
    bounds = _residue_boundaries(res_ids, res_names)
    n_res = bounds.shape[0] - 1
    seq = [str(res_names[bounds[i]]) for i in range(n_res)]

    patterns = sorted(
        ((name, LIPID21_RESIDUE_PATTERNS[name]) for name in candidates),
        key=lambda item: (-len(item[1]), item[0]),
    )

    offsets: List[int] = [0]
    species: List[str] = []
    res_local = np.zeros(res_ids.shape[0], dtype=np.int32)

    i = 0
    while i < n_res:
        for name, pattern in patterns:
            k = len(pattern)
            if tuple(seq[i:i + k]) == pattern:
                for j in range(k):
                    res_local[bounds[i + j]:bounds[i + j + 1]] = j
                offsets.append(int(bounds[i + k]))
                species.append(name)
                i += k
                break
        else:
            raise PatchIntegrityError(
                "cannot segment patch into molecules: unexpected residue "
                f"sequence {seq[i:i + 4]!r} at residue {i}. "
                f"Known patterns: "
                + "; ".join(f"{n}={'-'.join(p)}" for n, p in patterns)
            )

    return (
        np.asarray(offsets, dtype=np.int64),
        np.asarray(species, dtype="<U8"),
        res_local,
    )


# ---------------------------------------------------------------------------
# Leaflet identification (phosphorus based)
# ---------------------------------------------------------------------------

def _wrap_symmetric(value: float, period: float) -> float:
    """Map ``value`` into ``(-period/2, period/2]``.

    Every z position in this module is only defined modulo the box, so it has to
    be reduced to one canonical branch before it can be compared against zero.
    Skipping this is not harmless: an midplane that is really at z = 0 comes back
    from :func:`_circular_mean` as z = Lz, and "correcting" by that much shifts
    the whole bilayer a full box off centre.
    """
    return value - period * round(value / period)


def _circular_mean(
    values: np.ndarray, period: float, *, min_resultant: float = 0.15
) -> float:
    """Mean of ``values`` on a circle of circumference ``period``.

    Needed because a deposited frame's bilayer can sit right on the periodic
    z boundary (two of the bundled patches do), where a plain arithmetic mean
    lands half a box away from the real midplane.  The result is reduced to
    ``(-period/2, period/2]``.

    Callers must pass a *single* slab (lipids only, water stripped).  A bimodal
    distribution with its modes half a box apart has a near-zero resultant and
    an essentially arbitrary mean direction, so a short resultant is treated as
    a failure rather than silently mis-centring the bilayer.
    """
    angle = values * (2.0 * math.pi / period)
    cos_mean = float(np.cos(angle).mean())
    sin_mean = float(np.sin(angle).mean())
    resultant = math.hypot(cos_mean, sin_mean)
    if resultant < min_resultant:
        raise MembranePatchError(
            f"cannot locate the bilayer midplane: the z distribution has a "
            f"circular resultant of only {resultant:.3f} (need "
            f"{min_resultant}), so its mean direction is ill-defined.  This "
            "happens when the slab fills the box or when solvent was not "
            "stripped before measuring."
        )
    mean = math.atan2(sin_mean, cos_mean)
    return _wrap_symmetric(mean * period / (2.0 * math.pi), period)


def _min_image(delta: np.ndarray, length: float) -> np.ndarray:
    return delta - length * np.round(delta / length)


def assign_leaflets(assembly: LipidAssembly) -> Tuple[np.ndarray, float, float]:
    """Assign every molecule to a leaflet using lipid phosphorus atoms.

    The midplane is derived *only* from phosphorus z-coordinates.  Molecules
    that have no phosphorus of their own (sterols) are then placed by their own
    centroid relative to that phosphorus-defined midplane; a patch with no
    phosphorus at all is refused rather than guessed at.

    Returns (leaflet, midplane_z, phosphate_thickness).
    """
    is_p = np.isin(assembly.atom_names, list(PHOSPHORUS_ATOM_NAMES))
    if not is_p.any():
        species = ", ".join(sorted(set(assembly.mol_species.tolist())))
        raise MembranePatchError(
            "Cannot identify leaflets: this patch contains no lipid phosphorus "
            f"atoms (species present: {species}).\n"
            "Leaflet assignment is defined here by phosphorus position, and "
            "refusing is safer than guessing from a sterol's ring orientation. "
            "Use a patch that contains at least one phosphocholine/-ethanolamine/"
            "-serine lipid, e.g. the POPC_CHOL patches for sterol-containing "
            "bilayers."
        )

    lz = float(assembly.box[2])
    atom_mol = assembly.atom_molecule

    # Coarse midplane from the whole lipid slab, robust to a bilayer that
    # straddles the periodic z boundary.
    coarse = _circular_mean(assembly.xyz[assembly.heavy_mask, 2], lz)

    # Refine using phosphorus only, unwrapped around the coarse estimate.
    pz = coarse + _min_image(assembly.xyz[is_p, 2] - coarse, lz)
    upper = pz > coarse
    if not upper.any() or upper.all():
        raise MembranePatchError(
            "Phosphorus atoms do not separate into two leaflets; this does not "
            "look like a bilayer patch."
        )
    z_up = float(pz[upper].mean())
    z_dn = float(pz[~upper].mean())
    midplane = _wrap_symmetric(0.5 * (z_up + z_dn), lz)
    thickness = abs(z_up - z_dn)

    # Molecules carrying phosphorus: leaflet straight from their own P.
    leaflet = np.zeros(assembly.n_molecules, dtype=np.int8)
    p_mol = atom_mol[is_p]
    p_z = midplane + _min_image(assembly.xyz[is_p, 2] - midplane, lz)
    for mol, z in zip(p_mol, p_z):
        leaflet[mol] = 1 if z > midplane else -1

    # Sterols and anything else without P: own centroid vs the P midplane.
    missing = np.nonzero(leaflet == 0)[0]
    if missing.size:
        centroids = assembly.molecule_centroids()
        cz = midplane + _min_image(centroids[missing, 2] - midplane, lz)
        leaflet[missing] = np.where(cz > midplane, 1, -1).astype(np.int8)
        unknown = {
            str(s)
            for s in assembly.mol_species[missing]
            if str(s) not in STEROL_SPECIES
        }
        if unknown:
            raise MembranePatchError(
                "Molecules without phosphorus and not registered as sterols "
                f"cannot be assigned to a leaflet: {sorted(unknown)}"
            )

    return leaflet, midplane, thickness


# ---------------------------------------------------------------------------
# Patch loading
# ---------------------------------------------------------------------------

def load_patch(lipid: str = "POPC", *, verify: bool = True) -> LipidAssembly:
    """Load the bundled pre-equilibrated patch for ``lipid``.

    The bundled asset is *dry* (lipids only) and canonical: molecules are whole,
    wrapped into the box in xy, and the bilayer midplane is already at z = 0.

    Parameters
    ----------
    lipid:
        Lipid name or alias, e.g. ``"POPC"``, ``"POPC:CHOL"``, ``"SM"``.
    verify:
        Check the file's SHA-256 against the manifest.  Leave this on; it is a
        few milliseconds and it catches a truncated or tampered asset before it
        turns into an inscrutable grompp failure.

    Raises
    ------
    PatchNotAvailableError
        If no patch is bundled for ``lipid`` (with a message listing what is).
    PatchIntegrityError
        If the asset fails its checksum or structural checks.
    """
    spec, _ = _resolve_key(lipid)
    path = patch_data_dir() / spec.filename
    if not path.is_file():
        raise PatchNotAvailableError(
            f"Patch '{spec.key}' is listed in the manifest but its file "
            f"{path} is missing.  Reinstall PRISM from source, or set "
            f"{PATCH_DIR_ENV}."
        )
    if verify:
        digest = _sha256(path)
        if digest != spec.sha256:
            raise PatchIntegrityError(
                f"Checksum mismatch for {path}.\n"
                f"  expected sha256 {spec.sha256}\n"
                f"  found    sha256 {digest}\n"
                "The bundled patch is corrupt; reinstall PRISM."
            )

    res_ids, res_names, atom_names, xyz, box = _read_gro(path)
    if xyz.shape[0] != spec.n_atoms:
        raise PatchIntegrityError(
            f"{path}: manifest says {spec.n_atoms} atoms, file has {xyz.shape[0]}"
        )

    candidates = [name for name, _ in spec.composition]
    unknown = [c for c in candidates if c not in LIPID21_RESIDUE_PATTERNS]
    if unknown:
        raise PatchIntegrityError(
            f"manifest lists species with no residue pattern: {unknown}"
        )
    offsets, species, res_local = _segment_molecules(res_ids, res_names, candidates)

    assembly = LipidAssembly(
        atom_names=atom_names,
        res_names=res_names,
        res_local=res_local,
        elements=_elements_from_names(atom_names),
        xyz=xyz,
        mol_offsets=offsets,
        mol_species=species,
        mol_leaflet=np.zeros(species.shape[0], dtype=np.int8),
        box=box,
        key=spec.key,
    )

    counts = assembly.species_counts()
    expected = {name: c for name, c in spec.composition}
    if counts != expected:
        raise PatchIntegrityError(
            f"{path}: composition {counts} does not match manifest {expected}"
        )

    leaflet, midplane, _ = assign_leaflets(assembly)
    assembly.mol_leaflet = leaflet
    if abs(midplane) > 1.0:
        raise PatchIntegrityError(
            f"{path}: bilayer midplane is at z = {midplane:.2f} A, expected 0; "
            "the bundled asset is not in canonical form."
        )
    assembly.core_half_width = _core_half_width(assembly)
    _check_molecules_whole(assembly)
    return assembly


def _check_molecules_whole(assembly: LipidAssembly) -> None:
    """Refuse a patch whose molecules are broken across the periodic boundary.

    Tiling translates whole molecules, so a molecule split across a boundary
    would be smeared across the tiled box.  Cheap guard, catastrophic failure
    if skipped.
    """
    offsets = assembly.mol_offsets
    lo = np.minimum.reduceat(assembly.xyz, offsets[:-1], axis=0)
    hi = np.maximum.reduceat(assembly.xyz, offsets[:-1], axis=0)
    extent = hi - lo
    limit = assembly.box * 0.5
    bad = np.nonzero((extent > limit).any(axis=1))[0]
    if bad.size:
        raise PatchIntegrityError(
            f"{bad.size} molecule(s) span more than half the box "
            f"(first offender: molecule {int(bad[0])}, extent "
            f"{extent[bad[0]].round(2).tolist()} A vs half-box "
            f"{limit.round(2).tolist()} A).  They are broken across the "
            "periodic boundary and cannot be tiled."
        )


def _core_half_width(assembly: LipidAssembly, midplane: float = 0.0) -> float:
    """Half-width of the phosphate-free slab around the midplane.

    Measured as the |z| of the innermost headgroup phosphorus, so the slab
    ``|z| < core_half_width`` provably contains no phosphate.  It is the
    working definition of "hydrophobic core" used by the clash and contact
    diagnostics.  (The *water*-free slab is wider still; water cannot reach
    past the phosphates, and the bundled patches are dry regardless.)
    """
    is_p = np.isin(assembly.atom_names, list(PHOSPHORUS_ATOM_NAMES))
    if not is_p.any():
        return 0.0
    dz = _min_image(assembly.xyz[is_p, 2] - midplane, float(assembly.box[2]))
    return float(np.abs(dz).min())


def patch_geometry(assembly: LipidAssembly) -> PatchGeometry:
    """Measure a patch or tiled patch.

    Area per lipid counts only phosphorus-bearing lipids per leaflet, which is
    the convention experimental APL values use; sterols are reported in
    ``species_counts`` but excluded from the APL denominator.
    """
    _, _, thickness = assign_leaflets(assembly)
    area = float(assembly.box[0] * assembly.box[1])
    lipid_mols = ~np.isin(assembly.mol_species, list(STEROL_SPECIES))
    n_up = int(np.count_nonzero(lipid_mols & (assembly.mol_leaflet > 0)))
    n_dn = int(np.count_nonzero(lipid_mols & (assembly.mol_leaflet < 0)))
    denom = max(1, (n_up + n_dn)) / 2.0
    return PatchGeometry(
        box=(float(assembly.box[0]), float(assembly.box[1]), float(assembly.box[2])),
        area=area,
        n_upper=n_up,
        n_lower=n_dn,
        area_per_lipid=area / denom,
        phosphate_thickness=thickness,
        core_half_width=assembly.core_half_width,
        species_counts=assembly.species_counts(),
        n_atoms=assembly.n_atoms,
        n_molecules=assembly.n_molecules,
    )


# ---------------------------------------------------------------------------
# Tiling
# ---------------------------------------------------------------------------

def required_tiling(
    assembly: LipidAssembly,
    target_x: float,
    target_y: float,
) -> Tuple[int, int]:
    """Smallest (nx, ny) whose tiled box covers ``target_x`` x ``target_y``."""
    if target_x <= 0 or target_y <= 0:
        raise ValueError("target box lengths must be positive")
    nx = max(1, int(math.ceil(target_x / float(assembly.box[0]) - 1e-9)))
    ny = max(1, int(math.ceil(target_y / float(assembly.box[1]) - 1e-9)))
    return nx, ny


def tile_patch(
    assembly: LipidAssembly,
    target_x: Optional[float] = None,
    target_y: Optional[float] = None,
    *,
    nx: Optional[int] = None,
    ny: Optional[int] = None,
    center_xy: Tuple[float, float] = (0.0, 0.0),
) -> LipidAssembly:
    """Replicate the patch in xy to cover a requested area.

    Only *integer* replication is performed.  Cutting a tiled patch back to an
    exact target box was considered and rejected: the cut plane is not a
    periodicity of the patch, so the two faces that meet across it are
    uncorrelated and butt into each other, which reintroduces exactly the hard
    overlaps and low-density seams this route exists to avoid.  Integer tiling
    keeps the equilibrated structure and its periodicity perfectly intact, at
    the cost of quantising the box to multiples of the patch edge.  Pass ``nx``
    / ``ny`` to take direct control.

    Both leaflets are replicated together by pure xy translation, so the
    leaflets stay consistent by construction and the per-leaflet counts simply
    scale by ``nx * ny`` (asserted below).

    Molecule ordering follows the module contract: species, then leaflet
    (lower, upper), then tile row-major in x then y, then source molecule index.

    Parameters
    ----------
    target_x, target_y:
        Box lengths to cover, in Angstrom.  Ignored if ``nx``/``ny`` are given.
    center_xy:
        The tiled patch is centred on this xy point.  Pass the solute's xy
        centre so the solute needs no translation.
    """
    if nx is None or ny is None:
        if target_x is None or target_y is None:
            raise ValueError("supply either target_x/target_y or nx/ny")
        auto_nx, auto_ny = required_tiling(assembly, target_x, target_y)
        nx = auto_nx if nx is None else nx
        ny = auto_ny if ny is None else ny
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be >= 1")

    lx, ly, lz = (float(v) for v in assembly.box)
    n_mol = assembly.n_molecules
    counts = np.diff(assembly.mol_offsets)
    starts = assembly.mol_offsets[:-1]

    # ---- wrap whole molecules into [0, L) in xy before replicating --------
    # The deposited frames are not wrapped at all (atoms sit several box
    # lengths away from the origin).  Replicating them unwrapped would leave
    # tiles overlapping each other and voids elsewhere.
    centroids = assembly.molecule_centroids()
    shift = np.zeros((n_mol, 3))
    shift[:, 0] = -np.floor(centroids[:, 0] / lx) * lx
    shift[:, 1] = -np.floor(centroids[:, 1] / ly) * ly
    per_atom_shift = np.repeat(shift, counts, axis=0)
    base_xyz = assembly.xyz + per_atom_shift

    # ---- deterministic molecule order ------------------------------------
    species_order = _canonical_species_order(assembly.mol_species)
    species_rank = np.array(
        [species_order.index(str(s)) for s in assembly.mol_species], dtype=np.int64
    )
    # leaflet: lower (-1) first -> rank 0, upper (+1) -> rank 1
    leaflet_rank = (assembly.mol_leaflet > 0).astype(np.int64)

    tiles = [(ix, iy) for ix in range(nx) for iy in range(ny)]

    out_species: List[np.ndarray] = []
    out_leaflet: List[np.ndarray] = []
    out_counts: List[np.ndarray] = []
    atom_index: List[np.ndarray] = []
    atom_shift: List[np.ndarray] = []

    for s_rank, name in enumerate(species_order):
        for l_rank in (0, 1):
            sel = np.nonzero((species_rank == s_rank) & (leaflet_rank == l_rank))[0]
            if sel.size == 0:
                continue
            sel = np.sort(sel)                      # source order within tile
            for ix, iy in tiles:
                delta = np.array([ix * lx, iy * ly, 0.0])
                for mol in sel:
                    atom_index.append(
                        np.arange(starts[mol], starts[mol] + counts[mol],
                                  dtype=np.int64)
                    )
                    atom_shift.append(np.broadcast_to(delta, (counts[mol], 3)))
                out_species.append(assembly.mol_species[sel])
                out_leaflet.append(assembly.mol_leaflet[sel])
                out_counts.append(counts[sel])

    idx = np.concatenate(atom_index)
    shifts = np.concatenate(atom_shift)
    new_counts = np.concatenate(out_counts)
    new_offsets = np.zeros(new_counts.size + 1, dtype=np.int64)
    np.cumsum(new_counts, out=new_offsets[1:])

    box = np.array([lx * nx, ly * ny, lz])
    # Centre the tiled slab on the requested xy point; z is already midplane-0.
    origin = np.array([
        center_xy[0] - 0.5 * box[0],
        center_xy[1] - 0.5 * box[1],
        0.0,
    ])

    tiled = LipidAssembly(
        atom_names=assembly.atom_names[idx],
        res_names=assembly.res_names[idx],
        res_local=assembly.res_local[idx],
        elements=assembly.elements[idx],
        xyz=base_xyz[idx] + shifts + origin,
        mol_offsets=new_offsets,
        mol_species=np.concatenate(out_species),
        mol_leaflet=np.concatenate(out_leaflet),
        box=box,
        key=assembly.key,
        n_tiles=(nx, ny),
        core_half_width=assembly.core_half_width,
    )

    # ---- leaflet consistency -------------------------------------------
    before = assembly.leaflet_counts()
    after = tiled.leaflet_counts()
    for name, (dn, up) in before.items():
        want = (dn * nx * ny, up * nx * ny)
        if after.get(name) != want:
            raise MembranePatchError(
                f"tiling broke leaflet consistency for {name}: expected {want}, "
                f"got {after.get(name)}"
            )
    if tiled.n_molecules != n_mol * nx * ny:
        raise MembranePatchError("tiling lost molecules")
    return tiled


def _canonical_species_order(species: np.ndarray) -> List[str]:
    """Stable species order: phosphorus-bearing lipids first (alphabetical),
    then sterols.  Deterministic and independent of file order."""
    names = sorted(set(str(s) for s in species))
    lipids = [n for n in names if n not in STEROL_SPECIES]
    sterols = [n for n in names if n in STEROL_SPECIES]
    return lipids + sterols


# ---------------------------------------------------------------------------
# Solute
# ---------------------------------------------------------------------------

@dataclass
class Solute:
    """An oriented solute, held verbatim so it can be emitted unmoved."""

    xyz: np.ndarray             # (N, 3) float64, Angstrom
    elements: np.ndarray        # (N,)   <U2
    atom_names: np.ndarray
    res_names: np.ndarray
    res_seq: np.ndarray
    chain_ids: np.ndarray
    records: List[str]          # verbatim ATOM/HETATM lines, no trailing newline
    ter_after: np.ndarray       # (N,) bool, True where a TER followed this atom
    path: str = ""

    @property
    def n_atoms(self) -> int:
        return int(self.xyz.shape[0])

    @property
    def heavy_mask(self) -> np.ndarray:
        return self.elements != "H"

    def extent(self) -> Tuple[np.ndarray, np.ndarray]:
        return self.xyz.min(axis=0), self.xyz.max(axis=0)

    def xy_center(self) -> Tuple[float, float]:
        lo, hi = self.extent()
        return (0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]))


def read_solute(path: os.PathLike | str) -> Solute:
    """Read an oriented solute from PDB.

    Coordinates are taken as-is: this module never rotates or translates the
    solute.  Original ATOM/HETATM lines are retained verbatim so PDB output
    reproduces chain ids, insertion codes, occupancies and B-factors exactly.
    """
    path = Path(path)
    if not path.is_file():
        raise MembranePatchError(f"solute file not found: {path}")

    records: List[str] = []
    names: List[str] = []
    res_names: List[str] = []
    res_seq: List[int] = []
    chains: List[str] = []
    elems: List[str] = []
    coords: List[Tuple[float, float, float]] = []
    ter_flags: List[bool] = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            tag = line[:6]
            if tag.startswith(("ATOM", "HETATM")):
                records.append(line)
                names.append(line[12:16].strip())
                res_names.append(line[17:20].strip())
                try:
                    res_seq.append(int(line[22:26]))
                except ValueError:
                    res_seq.append(len(res_seq) + 1)
                chains.append(line[21] if len(line) > 21 else " ")
                elems.append(line[76:78] if len(line) >= 78 else "")
                coords.append(
                    (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                )
                ter_flags.append(False)
            elif tag.startswith("TER") and ter_flags:
                ter_flags[-1] = True

    if not records:
        raise MembranePatchError(f"no ATOM/HETATM records in {path}")

    return Solute(
        xyz=np.asarray(coords, dtype=np.float64),
        elements=_elements_from_pdb(elems, names),
        atom_names=np.asarray(names, dtype="<U5"),
        res_names=np.asarray(res_names, dtype="<U5"),
        res_seq=np.asarray(res_seq, dtype=np.int64),
        chain_ids=np.asarray(chains, dtype="<U1"),
        records=records,
        ter_after=np.asarray(ter_flags, dtype=bool),
        path=str(path),
    )


# ---------------------------------------------------------------------------
# Assembly result containers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DeletionReport:
    """What :func:`delete_clashing_lipids` removed."""

    cutoff: float
    n_deleted: int
    deleted_by_species: Dict[str, int]
    deleted_by_leaflet: Tuple[int, int]      # (lower, upper)
    surviving_by_species: Dict[str, int]
    surviving_by_leaflet: Tuple[int, int]
    residual_contacts: int                   # solute-lipid heavy pairs < cutoff
    min_heavy_distance: float
    seconds: float

    #: 1.0 for plain deletion; < 1.0 when the lipids were tested against a
    #: solute shrunk in xy by that factor.  ``residual_contacts`` and
    #: ``min_heavy_distance`` are always measured against the reference the
    #: deletion actually used, i.e. the *scaled* solute when this is < 1.0.
    #: Distances to the full-size solute are in the two fields below and are
    #: expected to be short -- closing that gap is the point of the shrink.
    shrink_factor: float = 1.0

    #: Extra molecules removed purely to equalise the two leaflets, and the
    #: closest approach of a survivor to the *full-size* solute.
    n_leaflet_balanced: int = 0
    min_heavy_distance_full: float = float("nan")

    def describe(self) -> str:
        head = (
            f"cutoff {self.cutoff:.2f} A: deleted {self.n_deleted} lipids "
            f"({self.deleted_by_leaflet[0]} lower / {self.deleted_by_leaflet[1]} "
            f"upper); surviving "
            f"{self.surviving_by_leaflet[0]}/{self.surviving_by_leaflet[1]}; "
            f"residual contacts {self.residual_contacts}; "
            f"closest solute-lipid heavy pair {self.min_heavy_distance:.2f} A "
            f"[{self.seconds:.2f} s]"
        )
        if self.shrink_factor >= 1.0:
            return head
        return (
            f"shrink {self.shrink_factor:.3f} ({self.n_leaflet_balanced} extra "
            f"removed to balance leaflets); " + head
            + f"; closest approach to the full-size solute "
              f"{self.min_heavy_distance_full:.2f} A (grow-back pending)"
        )


@dataclass(frozen=True)
class GrowthPlan:
    """The schedule that turns a shrunken solute back into the real one.

    Holds both endpoints explicitly rather than a factor, so the interpolation
    is exact and independent of how the shrink was parameterised (a tapered
    shrink is not a single scalar).  Step ``i`` of ``n_steps`` places the solute
    at ``(1 - w) * scaled + w * target`` with ``w = i / (n_steps - 1)``, so step
    0 reproduces the geometry the deletion was performed against and the last
    step reproduces the solute exactly, bit for bit.
    """

    scaled_xyz: np.ndarray        # (N, 3) the solute the deletion was run against
    target_xyz: np.ndarray        # (N, 3) the real solute, verbatim
    n_steps: int
    shrink_factor: float
    center_xy: Tuple[float, float]
    band: float
    max_step_displacement: float
    travel: float                 # largest lateral displacement of any band atom

    def weight(self, step: int) -> float:
        if step < 0 or step >= self.n_steps:
            raise IndexError(f"step {step} outside 0..{self.n_steps - 1}")
        if self.n_steps == 1:
            return 1.0
        return step / float(self.n_steps - 1)

    def positions(self, step: int) -> np.ndarray:
        """Solute coordinates at ``step``.

        The last step returns ``target_xyz`` itself, not a rounded-off
        interpolation of it: the module's contract is that the solute is
        emitted exactly where it was read.
        """
        if step == self.n_steps - 1:
            return self.target_xyz.copy()
        w = self.weight(step)
        return (1.0 - w) * self.scaled_xyz + w * self.target_xyz

    def describe(self) -> str:
        return (
            f"grow shrink_factor {self.shrink_factor:.3f} -> 1.0 in "
            f"{self.n_steps} steps; largest lateral travel {self.travel:.1f} A "
            f"({self.max_step_displacement:.2f} A per step)"
        )


@dataclass(frozen=True)
class GrowthReport:
    """What :func:`grow_solute` did."""

    n_steps: int
    shrink_factor: float
    step_seconds: Tuple[float, ...]
    max_lipid_displacement: float        # over the whole growth, any heavy atom
    min_heavy_distance: float            # survivor-to-solute, at full size
    residual_hard_pairs: int             # heavy-atom pairs closer than 2.0 A
    seconds: float
    notes: Tuple[str, ...] = ()

    def describe(self) -> str:
        return (
            f"grew the solute back in {self.n_steps} steps "
            f"[{self.seconds:.1f} s]; lipids moved up to "
            f"{self.max_lipid_displacement:.1f} A; closest solute-lipid heavy "
            f"pair now {self.min_heavy_distance:.2f} A with "
            f"{self.residual_hard_pairs} pair(s) under 2.0 A"
        )


@dataclass(frozen=True)
class MembraneComposition:
    """Everything a caller needs to build a topology matching the coordinates.

    ``species`` is the authoritative ordering: a ``[molecules]`` section (or a
    tleap script) must list exactly these names with these counts, in this
    order, after the solute.
    """

    solute_atoms: int
    species: Tuple[Tuple[str, int], ...]
    residue_patterns: Dict[str, Tuple[str, ...]]
    per_leaflet: Dict[str, Tuple[int, int]]
    formal_charge: int
    box: Tuple[float, float, float]
    n_lipid_atoms: int
    n_atoms: int
    patch_key: str
    n_tiles: Tuple[int, int]
    citation: str

    def to_dict(self) -> dict:
        return {
            "solute_atoms": self.solute_atoms,
            "species": [list(item) for item in self.species],
            "residue_patterns": {
                k: list(v) for k, v in sorted(self.residue_patterns.items())
            },
            "per_leaflet": {k: list(v) for k, v in sorted(self.per_leaflet.items())},
            "formal_charge": self.formal_charge,
            "box": list(self.box),
            "n_lipid_atoms": self.n_lipid_atoms,
            "n_atoms": self.n_atoms,
            "patch_key": self.patch_key,
            "n_tiles": list(self.n_tiles),
            "ordering": (
                "solute first (verbatim), then lipids grouped by species in the "
                "order of 'species'; within a species lower leaflet then upper; "
                "molecules contiguous"
            ),
            "citation": self.citation,
        }

    def write_json(self, path: os.PathLike | str) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2, sort_keys=False)
            handle.write("\n")
        return path


@dataclass
class DrySystem:
    """Solute + bilayer, no water, no ions."""

    solute: Solute
    lipids: LipidAssembly
    box: np.ndarray
    deletion: Optional[DeletionReport] = None

    #: Set when the lipids were deleted against a shrunken solute.  While it is
    #: not ``None`` the bilayer is *not* finished: the survivors still sit
    #: inside the solute's real footprint and have to be pushed out by
    #: :func:`grow_solute` before the system is physical.
    growth: Optional[GrowthPlan] = None

    #: Set by :func:`grow_solute` once the grow-back has run.
    growth_report: Optional[GrowthReport] = None

    @property
    def needs_growth(self) -> bool:
        return self.growth is not None

    @property
    def n_atoms(self) -> int:
        return self.solute.n_atoms + self.lipids.n_atoms

    def composition(self) -> MembraneComposition:
        order = _canonical_species_order(self.lipids.mol_species)
        counts = self.lipids.species_counts()
        species = tuple((name, counts[name]) for name in order if counts.get(name))
        charge = sum(
            LIPID21_FORMAL_CHARGE.get(name, 0) * count for name, count in species
        )
        try:
            citation = patch_citation()
        except MembranePatchError:
            citation = ""
        return MembraneComposition(
            solute_atoms=self.solute.n_atoms,
            species=species,
            residue_patterns={
                name: LIPID21_RESIDUE_PATTERNS[name] for name, _ in species
            },
            per_leaflet=self.lipids.leaflet_counts(),
            formal_charge=charge,
            box=(float(self.box[0]), float(self.box[1]), float(self.box[2])),
            n_lipid_atoms=self.lipids.n_atoms,
            n_atoms=self.n_atoms,
            patch_key=self.lipids.key,
            n_tiles=self.lipids.n_tiles,
            citation=citation,
        )

    # -- writers ---------------------------------------------------------
    def _refuse_if_ungrown(self, what: str) -> None:
        """Refuse to emit a system whose grow-back has not been run.

        Between the shrink-deletion and :func:`grow_solute` the lipids are
        deliberately overlapping the solute -- on mGluR7 the closest heavy-atom
        pair is 0.27 A -- because they were placed against a shrunken copy of
        it.  Those coordinates are an intermediate, never an answer, and
        emitting them would hand the next stage a system no minimiser can
        rescue.  Cheap guard, silent catastrophe if skipped.
        """
        if self.growth is not None:
            raise MembranePatchError(
                f"Refusing to {what}: this system was built with "
                f"shrink_factor={self.growth.shrink_factor:.3f} and its lipids "
                "still overlap the solute (they were placed against a shrunken "
                "copy of it).  Run grow_solute() first, or rebuild with "
                "shrink_factor=None."
            )

    def write_pdb(self, path: os.PathLike | str) -> Path:
        """Write the merged dry system as PDB, one TER per molecule.

        Solute lines are echoed verbatim (only the atom serial is renumbered),
        guaranteeing its coordinates are untouched.  Atom serials wrap modulo
        100000 and residue numbers modulo 9999 because the PDB columns are that
        wide; molecule identity is carried by order and by the TER records, and
        consecutive residues always get different numbers.
        """
        self._refuse_if_ungrown("write a PDB")
        return self._write_pdb(path)

    def write_parametrization_pdb(self, path: os.PathLike | str) -> Path:
        """Write the system for a topology generator, NOT as a system.

        The grow-back is a chicken-and-egg: :func:`grow_solute` relaxes the
        bilayer with a force field, so a topology has to exist before the first
        increment -- but the only file the system exists as, at that moment, is
        the shrink-deletion intermediate whose lipids still overlap the solute.
        This writes exactly that file, once, so tleap can read it.

        It is sound because parametrization reads atom *identity and order*,
        never geometry: tleap assigns parameters from residue templates and
        infers no bond from a distance (PRISM declares every disulfide
        explicitly), so the topology it returns is the topology of the finished
        system, which is what makes it legitimate to keep using after the grow.
        The coordinates in this file are not, and must never reach a simulator
        -- which is why this is a separate name rather than a flag on
        :meth:`write_pdb`, and why :meth:`write_pdb` goes on refusing.  The
        REMARK says the same thing to whoever opens the file.

        Passing a finished (grown) system is allowed and simply writes it; the
        method's contract is "geometry not guaranteed", not "geometry wrong".
        """
        note = "PARAMETRIZATION INPUT ONLY -- do not simulate these coordinates."
        if self.growth is not None:
            note += (
                f"  Built with shrink_factor={self.growth.shrink_factor:.3f}; the "
                "lipids still overlap the solute and grow_solute() has not run."
            )
        return self._write_pdb(path, remarks=(note,))

    def _write_pdb(
        self, path: os.PathLike | str, remarks: Sequence[str] = ()
    ) -> Path:
        """The PDB writer itself, with no opinion on whether it should run."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        lines: List[str] = [
            "REMARK   1 PRISM fast membrane route: pre-equilibrated Lipid21 patch",
        ]
        for remark in remarks:
            for chunk in _wrap_remark(remark):
                lines.append(f"REMARK   1 {chunk}")
        for chunk in _wrap_remark(patch_citation()):
            lines.append(f"REMARK   1 {chunk}")
        lines.append(
            "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1"
            % (self.box[0], self.box[1], self.box[2], 90.0, 90.0, 90.0)
        )

        serial = 0
        for i, record in enumerate(self.solute.records):
            serial += 1
            padded = record.ljust(80)
            lines.append(padded[:6] + "%5d" % (serial % 100000) + padded[11:80].rstrip())
            if self.solute.ter_after[i]:
                serial += 1
                lines.append("TER   %5d" % (serial % 100000))
        if not self.solute.ter_after.any():
            serial += 1
            lines.append("TER   %5d" % (serial % 100000))

        lip = self.lipids
        res_seq = _global_residue_numbers(lip)
        offsets = lip.mol_offsets
        for m in range(lip.n_molecules):
            for a in range(offsets[m], offsets[m + 1]):
                serial += 1
                lines.append(
                    _pdb_atom_line(
                        serial=serial % 100000,
                        name=str(lip.atom_names[a]),
                        res_name=str(lip.res_names[a]),
                        chain="M",
                        res_seq=int(res_seq[a]),
                        xyz=lip.xyz[a],
                        element=str(lip.elements[a]),
                    )
                )
            serial += 1
            lines.append("TER   %5d" % (serial % 100000))
        lines.append("END")
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")
        return path

    def write_gro(self, path: os.PathLike | str, title: Optional[str] = None) -> Path:
        """Write the merged dry system as GRO."""
        self._refuse_if_ungrown("write a GRO")
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        lip = self.lipids
        sol = self.solute
        names = np.concatenate([sol.atom_names, lip.atom_names])
        res_names = np.concatenate([sol.res_names, lip.res_names])
        # Solute residue numbers are made contiguous from 1, then lipids follow.
        sol_res = _sequential_residue_numbers(sol.res_seq, sol.chain_ids)
        lip_res = _global_residue_numbers(lip) + int(sol_res.max())
        res_seq = np.concatenate([sol_res, lip_res])
        xyz = np.vstack([sol.xyz, lip.xyz])
        lines = _gro_lines(
            names, res_names, res_seq, xyz, self.box,
            title or f"PRISM dry membrane system ({lip.key}, "
                     f"{lip.n_tiles[0]}x{lip.n_tiles[1]} tiles)",
        )
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")
        return path


def _wrap_remark(text: str, width: int = 68) -> List[str]:
    out: List[str] = []
    for para in text.splitlines():
        para = para.strip()
        if not para:
            continue
        while len(para) > width:
            cut = para.rfind(" ", 0, width)
            cut = cut if cut > 0 else width
            out.append(para[:cut])
            para = para[cut:].lstrip()
        out.append(para)
    return out


def _sequential_residue_numbers(
    res_seq: np.ndarray, chain_ids: np.ndarray
) -> np.ndarray:
    """1-based residue counter that increments whenever (chain, resSeq) changes."""
    n = res_seq.shape[0]
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    changed = np.zeros(n, dtype=bool)
    changed[1:] = (res_seq[1:] != res_seq[:-1]) | (chain_ids[1:] != chain_ids[:-1])
    return np.cumsum(changed) + 1


def _global_residue_numbers(assembly: LipidAssembly) -> np.ndarray:
    """1-based residue counter over an assembly, restarting nowhere."""
    n = assembly.n_atoms
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    mol = assembly.atom_molecule
    changed = np.zeros(n, dtype=bool)
    changed[1:] = (assembly.res_local[1:] != assembly.res_local[:-1]) | (
        mol[1:] != mol[:-1]
    )
    return np.cumsum(changed) + 1


def _pdb_atom_line(
    *, serial: int, name: str, res_name: str, chain: str, res_seq: int,
    xyz: np.ndarray, element: str,
) -> str:
    # 4-character atom names start in column 13; shorter ones in column 14.
    field = name if len(name) >= 4 else " " + name
    return (
        "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
        % (
            serial,
            field[:4],
            res_name[:3],
            chain,
            ((res_seq - 1) % 9999) + 1,
            xyz[0],
            xyz[1],
            xyz[2],
            1.00,
            0.00,
            element[:2].rjust(2),
        )
    )


# ---------------------------------------------------------------------------
# Insertion
# ---------------------------------------------------------------------------

def insert_solute(
    lipids: LipidAssembly,
    solute: Solute,
    *,
    margin: float = 0.0,
    box_z: Optional[float] = None,
) -> DrySystem:
    """Place the oriented solute into the tiled bilayer without moving it.

    The solute's coordinates are used exactly as read: the bilayer patch is the
    thing that was positioned (midplane at z = 0, centred on the solute's xy
    centre by :func:`tile_patch`).  This function only validates the fit and
    merges.

    Parameters
    ----------
    margin:
        Minimum clearance required between the solute's xy bounding box and the
        tiled box edge.  A solute closer than this to the edge would interact
        with its own periodic image, and it also breaks the assumption that
        makes the clash search exact.
    box_z:
        Optional final box height.  Defaults to whatever spans both the solute
        and the bilayer with 10 A of headroom, so the caller has somewhere to
        put water.  The bilayer midplane stays at z = 0 either way.
    """
    lo, hi = solute.extent()
    box = lipids.box.astype(float).copy()
    # Clearance between the solute bbox and the periodic box faces, per side.
    #
    # Derived from the box and the solute extent alone, which makes it
    # independent of where the tiled slab's origin happens to sit.  Two earlier
    # formulations of this were wrong: pairing ``hi`` with the low face and
    # ``lo`` with the high face measures the box, not the gap (it returns
    # roughly a box length and so never trips), and locating the box centre from
    # ``lipids.xyz.min()`` is skewed because molecules are wrapped by centroid
    # and therefore overhang the face by up to half a molecule.
    #
    # This is exact whenever the slab is centred on the solute, which is what
    # :func:`tile_patch` does when passed ``center_xy=solute.xy_center()`` (the
    # module's own path), and it goes negative when the solute is wider than the
    # tiled box, which is the condition the guard below exists to catch.
    clearance = 0.5 * (box[:2] - (hi[:2] - lo[:2]))
    if (clearance < margin).any():
        raise MembranePatchError(
            "The solute does not fit inside the tiled bilayer with the "
            f"requested {margin:.1f} A margin: clearance is "
            f"{clearance.round(2).tolist()} A in x/y for a box of "
            f"{box[0]:.1f} x {box[1]:.1f} A.\n"
            "Tile one more patch along the short axis (increase nx/ny), or "
            "reduce the margin."
        )
    if (clearance < 0).any():
        raise MembranePatchError(
            "The solute sticks out of the tiled bilayer in xy; tile a larger "
            "patch."
        )

    if box_z is None:
        z_lo = min(float(lo[2]), float(lipids.xyz[:, 2].min()))
        z_hi = max(float(hi[2]), float(lipids.xyz[:, 2].max()))
        box[2] = (z_hi - z_lo) + 20.0
    else:
        box[2] = float(box_z)

    return DrySystem(solute=solute, lipids=lipids, box=box)


# ---------------------------------------------------------------------------
# Clash deletion
# ---------------------------------------------------------------------------

def _neighbour_pairs_within(
    query: np.ndarray,
    reference: np.ndarray,
    cutoff: float,
    box_xy: Tuple[float, float],
) -> np.ndarray:
    """Indices into ``query`` having a ``reference`` point within ``cutoff``.

    Periodic in x and y, unbounded in z (the bilayer is periodic laterally; z is
    not, and the solute may stick far out of the slab).  Uses ``scipy`` when it
    is importable and falls back to a plain-numpy cell list otherwise, so this
    module works on a base PRISM install where scipy is only an extra.
    """
    if query.size == 0 or reference.size == 0:
        return np.zeros(0, dtype=np.int64)
    lx, ly = float(box_xy[0]), float(box_xy[1])
    q = query.copy()
    r = reference.copy()
    q[:, 0] %= lx
    q[:, 1] %= ly
    r[:, 0] %= lx
    r[:, 1] %= ly

    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return _cell_list_within(q, r, cutoff, lx, ly)

    tree = cKDTree(r, boxsize=[lx, ly, 0])
    hits = tree.query_ball_point(q, cutoff, return_length=True)
    return np.nonzero(hits > 0)[0]


def _cell_list_within(
    query: np.ndarray,
    reference: np.ndarray,
    cutoff: float,
    lx: float,
    ly: float,
) -> np.ndarray:
    """Dependency-free periodic-in-xy neighbour search (scipy fallback)."""
    cell = max(cutoff, 1e-6)
    nx = max(1, int(math.floor(lx / cell)))
    ny = max(1, int(math.floor(ly / cell)))
    z0 = min(float(query[:, 2].min()), float(reference[:, 2].min()))
    nz_span = max(
        float(query[:, 2].max()), float(reference[:, 2].max())
    ) - z0
    nz = max(1, int(math.floor(nz_span / cell)) + 1)
    cx, cy = lx / nx, ly / ny

    def cell_of(points: np.ndarray) -> np.ndarray:
        i = np.clip((points[:, 0] / cx).astype(np.int64), 0, nx - 1)
        j = np.clip((points[:, 1] / cy).astype(np.int64), 0, ny - 1)
        k = np.clip(((points[:, 2] - z0) / cell).astype(np.int64), 0, nz - 1)
        return (i * ny + j) * nz + k

    ref_cell = cell_of(reference)
    order = np.argsort(ref_cell, kind="stable")
    ref_sorted = reference[order]
    cell_sorted = ref_cell[order]
    starts = np.searchsorted(cell_sorted, np.arange(nx * ny * nz), side="left")
    ends = np.searchsorted(cell_sorted, np.arange(nx * ny * nz), side="right")

    qi = np.clip((query[:, 0] / cx).astype(np.int64), 0, nx - 1)
    qj = np.clip((query[:, 1] / cy).astype(np.int64), 0, ny - 1)
    qk = np.clip(((query[:, 2] - z0) / cell).astype(np.int64), 0, nz - 1)

    cut2 = cutoff * cutoff
    hit = np.zeros(query.shape[0], dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                ii = (qi + di) % nx
                jj = (qj + dj) % ny
                kk = qk + dk
                valid = (kk >= 0) & (kk < nz)
                idx = np.nonzero(valid & ~hit)[0]
                if idx.size == 0:
                    continue
                cells = (ii[idx] * ny + jj[idx]) * nz + kk[idx]
                for q_i, c in zip(idx, cells):
                    s, e = starts[c], ends[c]
                    if s == e:
                        continue
                    d = ref_sorted[s:e] - query[q_i]
                    d[:, 0] -= lx * np.round(d[:, 0] / lx)
                    d[:, 1] -= ly * np.round(d[:, 1] / ly)
                    if (np.einsum("ij,ij->i", d, d) < cut2).any():
                        hit[q_i] = True
    return np.nonzero(hit)[0]


# ---------------------------------------------------------------------------
# Shrink-and-grow geometry
# ---------------------------------------------------------------------------

def transmembrane_mask(
    solute: Solute, band: float, *, heavy_only: bool = True
) -> np.ndarray:
    """Atoms of ``solute`` inside the bilayer slab ``|z| <= band``.

    ``band`` is meant to be a patch's ``core_half_width`` -- the |z| of the
    innermost headgroup phosphorus -- so the mask picks out exactly the part of
    the solute that a lipid tail can ever touch.  Everything the shrink does is
    keyed off this subset, because it is the only part of the solute the
    bilayer has an opinion about.
    """
    if band <= 0:
        raise ValueError("band must be positive")
    sel = np.abs(solute.xyz[:, 2]) <= float(band)
    if heavy_only:
        sel &= solute.heavy_mask
    return sel


def transmembrane_cross_section(
    solute: Solute,
    band: float,
    *,
    probe: float = CROSS_SECTION_PROBE,
    resolution: float = 0.5,
) -> float:
    """Lateral area, in A^2, that the solute takes away from the bilayer.

    Rasterises the xy projection of the bilayer-spanning heavy atoms, each
    inflated to ``probe``, onto a ``resolution`` grid and counts occupied
    cells.  A projection is the right measure and a convex hull is not: a
    V-shaped TM bundle, or the cleft between the two protomers of a dimer, is
    open to lipid and must not be charged to the protein.

    Deterministic, and cheap (one pass over a few thousand atoms).
    """
    sel = transmembrane_mask(solute, band)
    if not sel.any():
        return 0.0
    xy = solute.xyz[sel, :2]
    res = float(resolution)
    if res <= 0:
        raise ValueError("resolution must be positive")
    origin = xy.min(axis=0) - probe - res
    cells = np.floor((xy - origin) / res).astype(np.int64)

    r = int(math.ceil(probe / res))
    di, dj = np.mgrid[-r:r + 1, -r:r + 1]
    disc = (di ** 2 + dj ** 2) * res * res <= probe * probe
    off_i, off_j = di[disc], dj[disc]

    # Flatten (i, j) to one integer.  The stride has to span the *j* range, not
    # the i range, or two different cells collide and the area is undercounted.
    stride = int(cells[:, 1].max()) + 2 * r + 3
    flat = (cells[:, 0][:, None] + off_i[None, :]) * stride + (
        cells[:, 1][:, None] + off_j[None, :]
    )
    return float(np.unique(flat).size * res * res)


def scale_solute(
    solute: Solute,
    factor: float,
    *,
    center: Optional[Tuple[float, float]] = None,
    band: Optional[float] = None,
    taper: float = 1.0,
) -> Solute:
    """Return a copy of ``solute`` contracted in x and y by ``factor``.

    z is never touched, so the scaled copy still spans the bilayer in exactly
    the same slab and every leaflet sees the same contraction.

    ``center``
        Point the contraction is about.  Defaults to the xy centroid of the
        *bilayer-spanning* heavy atoms when ``band`` is given, and to the xy
        centre of the bounding box otherwise (which is what OpenMM uses).  The
        centroid is the better choice whenever the two differ: mGluR7's
        transmembrane centroid sits 19.7 A from its bounding-box centre because
        of the extracellular Venus-flytrap domain, and contracting about the
        box centre would translate the whole TM bundle sideways by that much on
        top of shrinking it, biasing which side of the protein keeps its
        lipids.

    ``band``, ``taper``
        With ``band`` set, the contraction is applied at full strength inside
        ``|z| <= band`` and eased smoothly (smoothstep) to no contraction at
        all by ``|z| = band * (1 + taper)``.  Outside the bilayer the shrink
        buys nothing -- there are no lipids there to make room for -- while it
        costs a great deal: a uniform 0.5x on this dimer drives each protomer's
        extracellular domain ~40 A laterally, straight through its partner,
        and that clash then dominates the very relaxation the grow-back relies
        on.  Pass ``band=None`` for the uniform, OpenMM-faithful behaviour.
    """
    factor = float(factor)
    if not 0.0 < factor <= 1.0:
        raise ValueError("shrink factor must lie in (0, 1]")
    if center is None:
        if band is not None:
            sel = transmembrane_mask(solute, band)
            if not sel.any():
                raise MembranePatchError(
                    f"No solute heavy atom lies inside the bilayer slab "
                    f"|z| <= {band:.2f} A, so there is no transmembrane region "
                    "to contract about.  Either the solute is not oriented "
                    "(its membrane-spanning axis must be z, with the midplane "
                    "at z = 0) or it is not a membrane protein at all."
                )
            center = tuple(solute.xyz[sel, :2].mean(axis=0))
        else:
            center = solute.xy_center()

    xyz = solute.xyz.copy()
    strength = np.ones(xyz.shape[0])
    if band is not None:
        if taper < 0:
            raise ValueError("taper must be >= 0")
        inner = float(band)
        outer = inner * (1.0 + float(taper))
        az = np.abs(xyz[:, 2])
        if outer > inner:
            t = np.clip((az - inner) / (outer - inner), 0.0, 1.0)
            strength = 1.0 - (t * t * (3.0 - 2.0 * t))     # smoothstep, C1
        else:
            strength = (az <= inner).astype(float)
    eff = 1.0 - strength * (1.0 - factor)                  # 1 outside, factor inside
    xyz[:, 0] = center[0] + eff * (xyz[:, 0] - center[0])
    xyz[:, 1] = center[1] + eff * (xyz[:, 1] - center[1])
    return replace(solute, xyz=xyz)


def _deleted_at(
    lipids: LipidAssembly,
    ref_xyz: np.ndarray,
    cutoff: float,
    *,
    heavy_only: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """(per-molecule nearest-reference distance, keep mask) at ``cutoff``."""
    box_xy = (float(lipids.box[0]), float(lipids.box[1]))
    hm = lipids.heavy_mask if heavy_only else np.ones(lipids.n_atoms, dtype=bool)
    per_mol = _per_molecule_min_distance(lipids.xyz[hm], lipids.atom_molecule[hm],
                                         lipids.n_molecules, ref_xyz, box_xy)
    return per_mol, per_mol >= cutoff


def _per_molecule_min_distance(
    xyz: np.ndarray,
    mol_of: np.ndarray,
    n_molecules: int,
    reference: np.ndarray,
    box_xy: Tuple[float, float],
) -> np.ndarray:
    """Closest approach of each molecule to ``reference``, periodic in xy."""
    out = np.full(n_molecules, np.inf)
    if xyz.size == 0 or reference.size == 0:
        return out
    lx, ly = float(box_xy[0]), float(box_xy[1])
    q = xyz.copy()
    q[:, 0] %= lx
    q[:, 1] %= ly
    r = reference.copy()
    r[:, 0] %= lx
    r[:, 1] %= ly
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        # scipy is only an extra, so this branch is the one a base install
        # takes.  The chunk is sized from the reference count to cap the
        # temporary at a few tens of MB: a fixed 2048-row chunk against a
        # 10,000-atom solute would allocate half a gigabyte per iteration.
        dist = np.empty(q.shape[0])
        chunk_rows = max(1, int(2_000_000 // max(1, r.shape[0])))
        for i in range(0, q.shape[0], chunk_rows):
            chunk = q[i:i + chunk_rows]
            d = chunk[:, None, :] - r[None, :, :]
            d[:, :, 0] -= lx * np.round(d[:, :, 0] / lx)
            d[:, :, 1] -= ly * np.round(d[:, :, 1] / ly)
            dist[i:i + chunk_rows] = np.sqrt(
                np.einsum("ijk,ijk->ij", d, d)
            ).min(axis=1)
    else:
        dist = cKDTree(r, boxsize=[lx, ly, 0]).query(q, k=1)[0]
    np.minimum.at(out, mol_of, dist)
    return out


def area_budget_shrink_factor(
    lipids: LipidAssembly,
    solute: Solute,
    cutoff: float,
    *,
    band: Optional[float] = None,
    target_apl: Optional[float] = None,
    taper: float = 1.0,
    heavy_only: bool = True,
    bounds: Tuple[float, float] = SHRINK_FACTOR_BOUNDS,
    tolerance: int = 1,
    max_iterations: int = 24,
) -> Tuple[float, dict]:
    """Solve for the shrink factor that deletes the *right number* of lipids.

    OpenMM shrinks by a flat 0.5 regardless of the protein.  That is a
    reasonable default only while the protein is small compared with the patch,
    because the lipids deleted scale roughly with the *scaled* footprint --
    about a quarter of the true one at 0.5 -- and the survivors then have to be
    squeezed into an area that was never meant to hold them.  Measured on the
    mGluR7 dimer (transmembrane cross-section 2405 A^2 in a 23205 A^2 box, a
    tenth of the whole bilayer), a flat 0.5 deletes 50 of 768 lipids where the
    area budget calls for 80, and the free area per lipid lands at 57.9 A^2
    against the patch's own 60.4 -- a 4% compression of a bilayer that was
    handed to us already equilibrated.

    So the quantity to hold fixed is not the factor, it is the area:

        (box area - solute cross-section) / (surviving lipids / 2) == target_apl

    with ``target_apl`` defaulting to the tiled patch's own area per lipid.
    That fixes the survivor count, and this function bisects the factor to hit
    it.  The result therefore scales with the protein automatically, and
    reduces to something close to 0.5 for a small one.

    Returns ``(factor, diagnostics)``.  Deterministic: bisection on a purely
    geometric count, with the bracket and every probe recorded in the
    diagnostics so the solve can be audited.

    The count is *usually* monotone in the factor but is not guaranteed to be:
    growing a solute that has a pore or an open cleft can release a lipid that
    was sitting inside it, and channels and transporters are exactly the
    clientele here.  Bisection still terminates and still returns a factor
    inside ``bounds``; ``diagnostics["achieved_deleted"]`` and
    ``diagnostics["clamped"]`` say how well it hit the budget, and callers
    should look at them rather than assume the target was met.
    """
    if band is None:
        band = lipids.core_half_width
    if band <= 0:
        raise MembranePatchError(
            "Cannot size the shrink: the patch reports no hydrophobic core "
            "half-width, so the bilayer-spanning slab is undefined."
        )
    area = float(lipids.box[0] * lipids.box[1])
    n_sterol = int(np.count_nonzero(np.isin(lipids.mol_species, list(STEROL_SPECIES))))
    n_lipid = lipids.n_molecules - n_sterol
    if target_apl is None:
        target_apl = area / max(1.0, n_lipid / 2.0)
    cross = transmembrane_cross_section(solute, band)
    wanted = 2.0 * (area - cross) / float(target_apl)
    # sterols ride along with their leaflet and are not counted in the APL
    wanted += n_sterol
    target_deleted = max(0, int(round(lipids.n_molecules - wanted)))

    lo, hi = float(bounds[0]), float(bounds[1])
    sel = solute.heavy_mask if heavy_only else np.ones(solute.n_atoms, dtype=bool)

    probes: Dict[float, int] = {}

    def deleted_at(factor: float) -> int:
        if factor in probes:
            return probes[factor]
        scaled = scale_solute(solute, factor, band=band, taper=taper)
        _, keep = _deleted_at(lipids, scaled.xyz[sel], cutoff,
                              heavy_only=heavy_only)
        probes[factor] = int(np.count_nonzero(~keep))
        return probes[factor]

    n_lo, n_hi = deleted_at(lo), deleted_at(hi)
    if target_deleted <= n_lo:
        best = lo
    elif target_deleted >= n_hi:
        best = hi
    else:
        for _ in range(max_iterations):
            mid = 0.5 * (lo + hi)
            n_mid = deleted_at(mid)
            if abs(n_mid - target_deleted) <= tolerance or hi - lo < 1e-3:
                lo = hi = mid
                break
            if n_mid < target_deleted:
                lo = mid
            else:
                hi = mid
        best = 0.5 * (lo + hi)
    best = min(max(best, float(bounds[0])), float(bounds[1]))
    achieved = deleted_at(best)
    return best, {
        "box_area": area,
        "cross_section": cross,
        "band": float(band),
        "target_apl": float(target_apl),
        "target_surviving": wanted,
        "target_deleted": target_deleted,
        "achieved_deleted": achieved,
        "clamped": bool(
            abs(best - bounds[0]) < 1e-9 or abs(best - bounds[1]) < 1e-9
        ) and abs(achieved - target_deleted) > tolerance,
        "probes": dict(sorted(probes.items())),
    }


def _growth_plan(
    solute: Solute,
    scaled: Solute,
    factor: float,
    band: float,
    max_step: float,
) -> GrowthPlan:
    sel = transmembrane_mask(solute, band)
    travel = float(
        np.linalg.norm(solute.xyz[sel, :2] - scaled.xyz[sel, :2], axis=1).max()
    ) if sel.any() else 0.0
    if max_step <= 0:
        raise ValueError("max_step_displacement must be positive")
    n_steps = max(2, int(math.ceil(travel / float(max_step))) + 1)
    centre = tuple(solute.xyz[sel, :2].mean(axis=0)) if sel.any() else solute.xy_center()
    return GrowthPlan(
        scaled_xyz=scaled.xyz.copy(),
        target_xyz=solute.xyz.copy(),
        n_steps=n_steps,
        shrink_factor=float(factor),
        center_xy=(float(centre[0]), float(centre[1])),
        band=float(band),
        max_step_displacement=float(max_step),
        travel=travel,
    )


def delete_clashing_lipids(
    system: DrySystem,
    cutoff: float = DEFAULT_CLASH_CUTOFF,
    *,
    heavy_only: bool = True,
    shrink_factor: Optional[float | str] = None,
    balance_leaflets: Optional[bool] = None,
    band: Optional[float] = None,
    taper: float = 1.0,
    max_step_displacement: float = DEFAULT_GROWTH_STEP,
) -> DrySystem:
    """Delete every lipid that overlaps the solute.

    Criterion: **delete a lipid if ANY of its heavy atoms lies within
    ``cutoff`` of ANY solute heavy atom.**  This was chosen by sweeping the
    alternatives; see the module report.  Any-atom/2.6 A leaves zero residual
    overlaps while deleting the fewest lipids.  Looser "N atoms within R" rules
    delete less but admit hard overlaps that no minimiser recovers from, and
    tighter cutoffs just delete more lipids for no gain.

    Hydrogens are ignored by default: the patch and the solute were equilibrated
    under different protocols, so heavy-atom overlap is the meaningful signal
    and H-H near-contacts relax away in the first few minimisation steps.

    Deterministic: a pure geometric predicate, evaluated once.  No iteration,
    no randomness.

    Shrink-and-grow
    ---------------
    ``shrink_factor`` switches on the Wolf et al. / OpenMM algorithm described
    in the module docstring.  It is ``None`` by default, which is exactly the
    behaviour above and the behaviour every existing caller and test relies on.

    * ``shrink_factor=None`` -- plain deletion.  The lipids that overlap the
      solute go; the rest stay exactly where the equilibrated patch put them,
      which leaves the bilayer's inner edge standing off the protein.
    * ``shrink_factor="auto"`` -- solve for the factor with
      :func:`area_budget_shrink_factor`, so the number of survivors matches the
      area actually left for them.  This is the recommended setting.
    * ``shrink_factor=<float>`` -- use that factor.  Pass
      :data:`OPENMM_SHRINK_FACTOR` to reproduce OpenMM exactly.

    With a shrink in force the returned system is **not yet physical**: its
    lipids overlap the real solute and ``.growth`` carries the plan that fixes
    that.  Feed it to :func:`grow_solute`.

    ``balance_leaflets`` defaults to ``True`` whenever a shrink is in force and
    ``False`` otherwise.  Equalising the leaflets costs molecules, so it is not
    imposed on the plain path where nobody asked for it; on the shrink path it
    is part of the published algorithm.  When the two leaflets differ, the
    surplus is made up from the survivors that are *closest* to the solute --
    those are the ones the grow-back would have had to shove hardest, so they
    are the cheapest to lose and dropping them relieves the tightest contacts.
    """
    start = time.perf_counter()
    lip = system.lipids
    sol = system.solute

    if band is None:
        band = lip.core_half_width
    if shrink_factor is not None and band <= 0:
        raise MembranePatchError(
            "Cannot shrink the solute: this patch reports no hydrophobic core "
            "half-width, so the bilayer-spanning slab that defines the shrink "
            "is undefined.  Pass an explicit 'band', or delete without a "
            "shrink factor."
        )
    budget: Optional[dict] = None
    scaled: Optional[Solute] = None
    factor = 1.0
    if shrink_factor is not None:
        if isinstance(shrink_factor, str):
            if shrink_factor.strip().lower() != "auto":
                raise ValueError(
                    f"shrink_factor must be a number, None, or 'auto', not "
                    f"{shrink_factor!r}"
                )
            factor, budget = area_budget_shrink_factor(
                lip, sol, cutoff, band=band, taper=taper, heavy_only=heavy_only
            )
        else:
            factor = float(shrink_factor)
            if not 0.0 < factor <= 1.0:
                raise ValueError("shrink_factor must lie in (0, 1]")
        # A factor of exactly 1 is not a shrink: scaling by it is the identity,
        # so building a scaled copy and a zero-travel growth plan would only
        # force the caller through a pointless pair of minimisations.
        if factor < 1.0:
            scaled = scale_solute(sol, factor, band=band, taper=taper)
    if balance_leaflets is None:
        balance_leaflets = scaled is not None

    reference = scaled if scaled is not None else sol

    sol_sel = reference.heavy_mask if heavy_only else np.ones(sol.n_atoms, dtype=bool)
    lip_sel = lip.heavy_mask if heavy_only else np.ones(lip.n_atoms, dtype=bool)
    sol_xyz = reference.xyz[sol_sel]
    lip_xyz = lip.xyz[lip_sel]
    lip_mol = lip.atom_molecule[lip_sel]

    box_xy = (float(lip.box[0]), float(lip.box[1]))
    hit_atoms = _neighbour_pairs_within(lip_xyz, sol_xyz, cutoff, box_xy)
    doomed = np.unique(lip_mol[hit_atoms])

    keep = np.ones(lip.n_molecules, dtype=bool)
    keep[doomed] = False

    n_balanced = 0
    if balance_leaflets:
        keep, n_balanced = _balance_leaflets(keep, lip, sol_xyz, box_xy)

    deleted_species: Dict[str, int] = {}
    for name in _canonical_species_order(lip.mol_species):
        n = int(np.count_nonzero((lip.mol_species == name) & ~keep))
        if n:
            deleted_species[name] = n
    del_lower = int(np.count_nonzero(~keep & (lip.mol_leaflet < 0)))
    del_upper = int(np.count_nonzero(~keep & (lip.mol_leaflet > 0)))

    survivors = lip.select_molecules(keep)

    # Verify the criterion actually held: no surviving heavy atom inside cutoff.
    surv_xyz = survivors.xyz[survivors.heavy_mask]
    residual = _neighbour_pairs_within(surv_xyz, sol_xyz, cutoff, box_xy).size
    min_dist = _min_distance(surv_xyz, sol_xyz, box_xy)

    report = DeletionReport(
        cutoff=float(cutoff),
        # counted from the final mask, not from ``doomed``, so it stays
        # consistent with deleted_by_species/-_by_leaflet once the leaflet
        # balancing has taken its extra molecules
        n_deleted=int(np.count_nonzero(~keep)),
        deleted_by_species=deleted_species,
        deleted_by_leaflet=(del_lower, del_upper),
        surviving_by_species=survivors.species_counts(),
        surviving_by_leaflet=(
            int(np.count_nonzero(survivors.mol_leaflet < 0)),
            int(np.count_nonzero(survivors.mol_leaflet > 0)),
        ),
        residual_contacts=int(residual),
        min_heavy_distance=float(min_dist),
        seconds=time.perf_counter() - start,
        shrink_factor=float(factor),
        n_leaflet_balanced=int(n_balanced),
        min_heavy_distance_full=(
            float(_min_distance(surv_xyz, sol.xyz[sol.heavy_mask], box_xy))
            if scaled is not None else float(min_dist)
        ),
    )
    plan = None
    if scaled is not None:
        plan = _growth_plan(sol, scaled, factor, float(band),
                            max_step_displacement)
    return DrySystem(
        solute=sol, lipids=survivors, box=system.box.copy(), deletion=report,
        growth=plan,
    )


def _balance_leaflets(
    keep: np.ndarray,
    lipids: LipidAssembly,
    reference_xyz: np.ndarray,
    box_xy: Tuple[float, float],
) -> Tuple[np.ndarray, int]:
    """Remove extra molecules so both leaflets lose the same number.

    An unequal deletion leaves the bilayer with a built-in area mismatch
    between its two halves, which relaxes by bending -- the leaflet holding
    more lipids buckles.  The published algorithm forbids it, and PRISM's plain
    deletion does not enforce it: on the mGluR7 dimer the shrink path deletes
    43 lower against 41 upper before balancing.

    Sterols are balanced together with their leaflet's lipids, because it is the
    leaflet's total area that has to match.
    """
    keep = keep.copy()
    lower = lipids.mol_leaflet < 0
    upper = lipids.mol_leaflet > 0
    n_del = {-1: int(np.count_nonzero(~keep & lower)),
             1: int(np.count_nonzero(~keep & upper))}
    surplus = max(n_del.values())
    if surplus == min(n_del.values()):
        return keep, 0

    hm = lipids.heavy_mask
    per_mol = _per_molecule_min_distance(
        lipids.xyz[hm], lipids.atom_molecule[hm], lipids.n_molecules,
        reference_xyz, box_xy,
    )
    added = 0
    for side, side_mask in ((-1, lower), (1, upper)):
        need = surplus - n_del[side]
        if need <= 0:
            continue
        candidates = np.nonzero(keep & side_mask)[0]
        if candidates.size < need:
            raise MembranePatchError(
                f"Cannot equalise the leaflets: the "
                f"{'lower' if side < 0 else 'upper'} leaflet needs {need} more "
                f"deletions but only has {candidates.size} molecules left.  The "
                "solute is not centred on the bilayer, or the clash cutoff is "
                "far too large."
            )
        # ties broken by molecule index so the choice is order-independent
        order = np.lexsort((candidates, per_mol[candidates]))[:need]
        keep[candidates[order]] = False
        added += int(need)
    return keep, added


def _min_distance(
    a: np.ndarray, b: np.ndarray, box_xy: Tuple[float, float]
) -> float:
    """Smallest periodic-in-xy distance between two point sets."""
    if a.size == 0 or b.size == 0:
        return float("inf")
    lx, ly = float(box_xy[0]), float(box_xy[1])
    aa = a.copy()
    bb = b.copy()
    aa[:, 0] %= lx
    aa[:, 1] %= ly
    bb[:, 0] %= lx
    bb[:, 1] %= ly
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        best = float("inf")
        for i in range(0, aa.shape[0], 4096):
            chunk = aa[i:i + 4096]
            d = chunk[:, None, :] - bb[None, :, :]
            d[:, :, 0] -= lx * np.round(d[:, :, 0] / lx)
            d[:, :, 1] -= ly * np.round(d[:, :, 1] / ly)
            best = min(best, float(np.sqrt((d * d).sum(axis=2)).min()))
        return best
    tree = cKDTree(bb, boxsize=[lx, ly, 0])
    return float(tree.query(aa, k=1)[0].min())


def contact_distances(
    system: DrySystem, *, band: Optional[float] = None
) -> np.ndarray:
    """Nearest-lipid distance for each solute surface atom facing the bilayer.

    This is the diagnostic for the defect that sank the earlier prototype: a
    vacuum annulus around the protein.  For every surviving lipid heavy atom in
    the hydrophobic band we find its nearest solute heavy atom; each solute atom
    that is nearest to at least one lipid atom is by definition on the
    lipid-facing surface, and we keep its true nearest-lipid distance.  The
    median of the result is the protein-lipid gap.  vdW contact is ~4 A.
    """
    lip = system.lipids
    sol = system.solute
    half = band if band is not None else max(lip.core_half_width, 10.0)

    lip_xyz = lip.xyz[lip.heavy_mask]
    lip_xyz = lip_xyz[np.abs(lip_xyz[:, 2]) <= half]
    sol_xyz = sol.xyz[sol.heavy_mask]
    sol_xyz = sol_xyz[np.abs(sol_xyz[:, 2]) <= half]
    if lip_xyz.size == 0 or sol_xyz.size == 0:
        return np.zeros(0)

    lx, ly = float(lip.box[0]), float(lip.box[1])
    lq = lip_xyz.copy()
    sq = sol_xyz.copy()
    lq[:, 0] %= lx
    lq[:, 1] %= ly
    sq[:, 0] %= lx
    sq[:, 1] %= ly
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:  # pragma: no cover - diagnostic only
        raise MembranePatchError(
            "contact_distances() needs scipy; install it or skip the diagnostic"
        ) from exc
    tree = cKDTree(sq, boxsize=[lx, ly, 0])
    dist, idx = tree.query(lq, k=1)
    best: Dict[int, float] = {}
    for d, i in zip(dist, idx):
        i = int(i)
        if d < best.get(i, float("inf")):
            best[i] = float(d)
    return np.sort(np.asarray(sorted(best.values())))


# ---------------------------------------------------------------------------
# Grow-back
# ---------------------------------------------------------------------------

#: A relaxer takes ``(step, n_steps, system)`` -- where ``system`` is a
#: :class:`DrySystem` holding the solute at that step's interpolated position
#: together with the current lipids -- and returns the relaxed lipid
#: coordinates, an ``(n_lipid_atoms, 3)`` array in Angstrom.  The solute must
#: come back exactly where it was handed over, which is why only the lipid
#: coordinates are returned at all.
Relaxer = Callable[[int, int, "DrySystem"], np.ndarray]


def grow_solute(
    system: DrySystem,
    relax: Relaxer,
    *,
    n_steps: Optional[int] = None,
    progress: Optional[Callable[[str], None]] = None,
) -> DrySystem:
    """Grow the shrunken solute back to full size, relaxing the bilayer.

    ``system`` must be a :class:`DrySystem` carrying a :class:`GrowthPlan`,
    i.e. one returned by :func:`delete_clashing_lipids` with a ``shrink_factor``.
    At each increment the solute is placed at the interpolated position and
    ``relax`` is called to let the lipids move out of its way; the last
    increment puts the solute back exactly where it was read, so the module's
    "the solute is emitted unmoved" contract survives the whole procedure.

    The solute is never handed to ``relax`` as something that may move.  It is
    the boundary condition, not a participant: a restrained solute would sag
    under the bilayer's lateral pressure and fight the interpolation, and the
    schedule would no longer mean what it says.

    Returns a new :class:`DrySystem` with the full-size solute, the relaxed
    lipids and ``growth=None`` -- the bilayer is finished.
    """
    plan = system.growth
    if plan is None:
        raise MembranePatchError(
            "This system was not built with a shrink, so there is nothing to "
            "grow.  Call delete_clashing_lipids(..., shrink_factor='auto') "
            "first."
        )
    if n_steps is not None:
        if n_steps < 2:
            raise ValueError("n_steps must be at least 2")
        plan = replace(plan, n_steps=int(n_steps))

    start = time.perf_counter()
    lipids = system.lipids
    initial = lipids.xyz.copy()
    step_seconds: List[float] = []
    notes: List[str] = []

    for step in range(plan.n_steps):
        t0 = time.perf_counter()
        stage = DrySystem(
            solute=replace(system.solute, xyz=plan.positions(step)),
            lipids=lipids,
            box=system.box.copy(),
        )
        new_xyz = np.asarray(relax(step, plan.n_steps, stage), dtype=np.float64)
        if new_xyz.shape != lipids.xyz.shape:
            raise MembranePatchError(
                f"relaxer returned {new_xyz.shape} coordinates for "
                f"{lipids.xyz.shape} lipid atoms at step {step}"
            )
        if not np.isfinite(new_xyz).all():
            raise MembranePatchError(
                f"relaxer returned non-finite lipid coordinates at grow step "
                f"{step} of {plan.n_steps}.  The bilayer blew up; reduce "
                "max_step_displacement or loosen the shrink factor."
            )
        lipids = replace(lipids, xyz=new_xyz)
        step_seconds.append(time.perf_counter() - t0)
        if progress is not None:
            progress(
                f"grow step {step + 1}/{plan.n_steps} "
                f"(weight {plan.weight(step):.3f}) "
                f"[{step_seconds[-1]:.1f} s]"
            )

    box_xy = (float(lipids.box[0]), float(lipids.box[1]))
    heavy = lipids.xyz[lipids.heavy_mask]
    sol_heavy = system.solute.xyz[system.solute.heavy_mask]
    min_dist = _min_distance(heavy, sol_heavy, box_xy)
    hard = _neighbour_pairs_within(heavy, sol_heavy, 2.0, box_xy).size
    if hard:
        notes.append(f"{hard} lipid heavy atom(s) still within 2.0 A of the solute")

    report = GrowthReport(
        n_steps=plan.n_steps,
        shrink_factor=plan.shrink_factor,
        step_seconds=tuple(step_seconds),
        max_lipid_displacement=float(
            np.linalg.norm(lipids.xyz - initial, axis=1).max()
        ) if lipids.n_atoms else 0.0,
        min_heavy_distance=float(min_dist),
        residual_hard_pairs=int(hard),
        seconds=time.perf_counter() - start,
        notes=tuple(notes),
    )
    return DrySystem(
        solute=system.solute, lipids=lipids, box=system.box.copy(),
        deletion=system.deletion, growth=None, growth_report=report,
    )


def _pbc_shift(system: DrySystem) -> np.ndarray:
    """Translation that centres the *solute* in the box.

    GROMACS folds every atom into the box and then rebuilds each molecule from
    its bond graph.  Lipids are repaired without trouble -- they are an order of
    magnitude smaller than half the box, and they are meant to overhang the
    faces, because :func:`tile_patch` wraps them by centroid.  The solute is the
    fragile one, so it is the only thing that gets centred here.

    Whether GROMACS can repair a broken solute depends on its longest *bonded*
    interaction, not on its overall size, so a protein wider than half the box
    is fine.  Not breaking it in the first place is still free, and it removes
    the failure mode entirely.

    z is handled differently from xy.  x and y are periodic for the bilayer, so
    a lipid that crosses a face there is physically unchanged.  z is not: the
    two leaflets are distinguishable, and a lipid folded across the z boundary
    comes back on the wrong side of the membrane.  The z shift is therefore
    taken from the solute *and* the lipids together, which always fits because
    :func:`insert_solute` sizes ``box_z`` from exactly that union plus 20 A.
    """
    lo, hi = system.solute.extent()
    box = np.asarray(system.box, dtype=float)
    span = hi - lo
    if (span > box + 1e-6).any():
        raise MembranePatchError(
            f"The solute spans {span.round(1).tolist()} A but its box is only "
            f"{box.round(1).tolist()} A, so no translation can keep it whole. "
            "Tile a larger patch (increase nx/ny), or raise box_z."
        )
    z_lo = min(float(lo[2]), float(system.lipids.xyz[:, 2].min()))
    z_hi = max(float(hi[2]), float(system.lipids.xyz[:, 2].max()))
    if z_hi - z_lo > box[2] + 1e-6:
        raise MembranePatchError(
            f"Solute and bilayer together span {z_hi - z_lo:.1f} A in z but "
            f"box_z is only {box[2]:.1f} A; raise box_z."
        )
    shift = 0.5 * (box - (hi + lo))
    shift[2] = 0.5 * (box[2] - (z_hi + z_lo))
    return shift


def gromacs_relaxer(
    *,
    topology: os.PathLike | str,
    mdp: os.PathLike | str,
    workdir: os.PathLike | str,
    index: Optional[os.PathLike | str] = None,
    freeze_group: str = "SOLU",
    gmx: str = "gmx",
    ntomp: Optional[int] = None,
    maxwarn: int = 20,
    final_mdp: Optional[os.PathLike | str] = None,
    keep_intermediates: bool = False,
    timeout: Optional[float] = None,
    env: Optional[Dict[str, str]] = None,
) -> Relaxer:
    """A :data:`Relaxer` that minimises the bilayer with ``gmx grompp``/``mdrun``.

    Why energy minimisation and not dynamics
    ----------------------------------------
    OpenMM grows the protein under Langevin dynamics at 10 K, taking 20 fs of
    dynamics per increment and therefore needing the increments to be tiny
    (0.25 A) -- of order a hundred of them.  Minimising *to convergence* at
    every increment instead leaves the bilayer fully relaxed at each stage
    rather than partly relaxed, so the increments can be several times larger
    for the same quality; on the mGluR7 dimer that is 15 steps instead of ~60.
    It also keeps the module deterministic, which a thermostat would not.

    Why freezing and not position restraints
    ----------------------------------------
    ``freezegrps`` makes the solute an immovable boundary, which is what the
    algorithm assumes.  Position restraints would let it sag under the
    bilayer's lateral pressure, and the interpolation schedule would stop
    describing where the solute actually is.  Freezing also needs no extra
    ``.itp`` and no force-constant choice.

    Parameters
    ----------
    topology:
        A ``.top`` whose ``[ molecules ]`` matches this system exactly -- the
        solute first, then the surviving lipids in
        :meth:`MembraneComposition.species` order.  Deleting lipids changes
        that count, so the topology has to be generated *after* the deletion.
    mdp:
        Minimisation parameters for the intermediate steps.  ``freezegrps`` and
        ``freezedim`` are appended (overriding anything already there), so the
        caller does not have to know about the freeze group.
    final_mdp:
        Optional, usually tighter, parameters for the last step, where the
        solute is at full size and the result is the one that gets kept.
    index:
        Index file defining ``freeze_group``.  One is written automatically
        from the solute/lipid split if this is ``None``.
    """
    # Resolved because ``_run`` launches gmx with cwd set to the stage
    # directory while passing it these paths; a relative workdir would then be
    # re-interpreted relative to the stage and every file would be "not found".
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    topology = Path(topology).resolve()
    mdp = Path(mdp).resolve()
    final_mdp = Path(final_mdp).resolve() if final_mdp is not None else None
    index = Path(index).resolve() if index is not None else None
    binary = shutil.which(gmx)
    if binary is None:
        raise MembranePatchError(
            f"'{gmx}' is not on PATH, so the bilayer cannot be relaxed while "
            "the solute grows back.  Either load a GROMACS module, pass "
            "gmx=<path>, or build without a shrink factor (which skips the "
            "grow-back and leaves the historical stand-off gap)."
        )
    if ntomp is None:
        ntomp = max(1, (os.cpu_count() or 1))

    def _write_mdp(source: Path, dest: Path) -> None:
        keep = [
            line for line in source.read_text().splitlines()
            if line.split("=")[0].strip().lower().replace("_", "-")
            not in ("freezegrps", "freezedim")
        ]
        keep += [
            "; appended by prism.membrane.patch.gromacs_relaxer",
            f"freezegrps = {freeze_group}",
            "freezedim  = Y Y Y",
        ]
        dest.write_text("\n".join(keep) + "\n")

    energies: List[Dict[str, float]] = []

    def relax(step: int, n_steps: int, system: DrySystem) -> np.ndarray:
        stage = workdir / f"grow{step:03d}"
        stage.mkdir(exist_ok=True)
        solute, lipids = system.solute, system.lipids
        shift = _pbc_shift(system)
        placed = DrySystem(
            solute=replace(solute, xyz=solute.xyz + shift),
            lipids=replace(lipids, xyz=lipids.xyz + shift),
            box=np.asarray(system.box, dtype=float),
        )
        conf = stage / "in.gro"
        placed.write_gro(conf, title=f"PRISM shrink-and-grow step {step}")

        ndx = index
        if ndx is None:
            ndx = stage / "freeze.ndx"
            _write_freeze_index(ndx, solute.n_atoms, lipids.n_atoms, freeze_group)

        run_mdp = stage / "em.mdp"
        _write_mdp(final_mdp if (final_mdp and step == n_steps - 1) else mdp, run_mdp)

        tpr = stage / "em.tpr"
        _run([binary, "grompp", "-f", str(run_mdp), "-c", str(conf),
              "-p", str(topology), "-n", str(ndx), "-o", str(tpr),
              "-po", str(stage / "mdout.mdp"), "-maxwarn", str(maxwarn)],
             stage, "grompp", env, timeout)
        _run([binary, "mdrun", "-s", str(tpr), "-deffnm", "em",
              "-ntmpi", "1", "-ntomp", str(ntomp), "-nb", "cpu"],
             stage, "mdrun", env, timeout)

        out = stage / "em.gro"
        if not out.is_file():
            raise MembranePatchError(
                f"gmx mdrun produced no {out.name} at grow step {step}; see "
                f"{stage / 'em.log'}"
            )
        _, _, _, xyz, _ = _read_gro(out)
        if xyz.shape[0] != solute.n_atoms + lipids.n_atoms:
            raise MembranePatchError(
                f"{out} has {xyz.shape[0]} atoms, expected "
                f"{solute.n_atoms + lipids.n_atoms} -- the topology's "
                "[ molecules ] section does not match this system."
            )
        energies.append(_parse_em_log(stage / "em.log"))
        relaxed = _rejoin_molecules(
            xyz[solute.n_atoms:] - shift, lipids.xyz, lipids.mol_offsets,
            placed.box,
        )
        if not keep_intermediates and step < n_steps - 1:
            shutil.rmtree(stage, ignore_errors=True)
        return relaxed

    relax.energies = energies               # type: ignore[attr-defined]
    return relax


def _write_freeze_index(
    path: Path, n_solute: int, n_lipid: int, freeze_group: str
) -> None:
    def block(name: str, lo: int, hi: int) -> str:
        out = [f"[ {name} ]"]
        row: List[str] = []
        for i in range(lo, hi + 1):
            row.append("%6d" % i)
            if len(row) == 15:
                out.append("".join(row))
                row = []
        if row:
            out.append("".join(row))
        return "\n".join(out)

    total = n_solute + n_lipid
    path.write_text(
        block(freeze_group, 1, n_solute) + "\n"
        + block("MEMB", n_solute + 1, total) + "\n"
        + block("System", 1, total) + "\n"
    )


def _rejoin_molecules(
    xyz: np.ndarray,
    previous: np.ndarray,
    offsets: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    """Undo GROMACS' wrapping: make molecules whole *and* put them back.

    ``mdrun`` folds every atom into the box, which does two separate things to
    a lipid that was straddling a face -- and :func:`tile_patch` wraps lipids by
    *centroid*, so molecules legitimately overhang the faces and straddling is
    the normal case, not an edge case:

    1. it splits the molecule across the boundary.  Re-joining is unambiguous
       by minimum image about the first atom, because a lipid is far smaller
       than half the box.
    2. it can move the whole re-joined molecule into the neighbouring periodic
       image, since the anchor atom itself may have been folded.  In xy that is
       physically a no-op, but it makes the emitted coordinates a scatter across
       images and it makes any displacement measured against the previous frame
       meaningless -- a no-op round trip reports moves of a whole box diagonal.
       In z it is not a no-op at all: an upper-leaflet lipid comes back below
       the lower leaflet.

    So each molecule is also re-imaged onto ``previous``: whichever box vector
    puts its centroid closest to where that molecule was before the relaxation.
    A lipid cannot travel half a box in one minimisation, so this is exact.
    """
    out = np.asarray(xyz, dtype=float).copy()
    previous = np.asarray(previous, dtype=float)
    box = np.asarray(box, dtype=float)
    starts = offsets[:-1]
    counts = np.diff(offsets)

    # (1) make every molecule whole about its own first atom
    anchors = np.repeat(out[starts], counts, axis=0)
    delta = out - anchors
    for dim in range(3):
        if box[dim] > 0:
            delta[:, dim] -= box[dim] * np.round(delta[:, dim] / box[dim])
    out = anchors + delta

    # (2) re-image each whole molecule onto where it was
    new_c = np.add.reduceat(out, starts, axis=0) / counts[:, None]
    old_c = np.add.reduceat(previous, starts, axis=0) / counts[:, None]
    offset = np.zeros_like(new_c)
    for dim in range(3):
        if box[dim] > 0:
            offset[:, dim] = -box[dim] * np.round(
                (new_c[:, dim] - old_c[:, dim]) / box[dim]
            )
    return out + np.repeat(offset, counts, axis=0)


def _parse_em_log(path: Path) -> Dict[str, float]:
    """Pull the converged energy and force out of an ``em.log``."""
    out: Dict[str, float] = {}
    if not path.is_file():
        return out
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("Potential Energy"):
            out["potential"] = float(line.split("=")[1])
        elif line.startswith("Maximum force"):
            out["fmax"] = float(line.split("=")[1].split()[0])
        elif "converged to Fmax" in line:
            out["converged"] = 1.0
        elif "did not converge" in line:
            out["converged"] = 0.0
    return out


def _run(cmd: List[str], cwd: Path, what: str,
         env: Optional[Dict[str, str]],
         timeout: Optional[float] = None) -> None:
    try:
        proc = subprocess.run(
            cmd, cwd=str(cwd), capture_output=True, text=True,
            env={**os.environ, **(env or {})}, timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise MembranePatchError(
            f"gmx {what} exceeded its {timeout:.0f} s budget in {cwd}.  Raise "
            "grow_timeout, or cut the minimisation's nsteps."
        ) from exc
    (cwd / f"{what}.log").write_text(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr or proc.stdout).splitlines()[-25:])
        raise MembranePatchError(
            f"gmx {what} failed (exit {proc.returncode}) in {cwd}:\n{tail}"
        )


# ---------------------------------------------------------------------------
# One-call convenience
# ---------------------------------------------------------------------------

def build_dry_membrane_system(
    solute: os.PathLike | str | Solute,
    *,
    lipid: str = "POPC",
    xy_margin: float = 12.0,
    nx: Optional[int] = None,
    ny: Optional[int] = None,
    clash_cutoff: float = DEFAULT_CLASH_CUTOFF,
    box_z: Optional[float] = None,
    shrink_factor: Optional[float | str] = None,
    balance_leaflets: Optional[bool] = None,
    max_step_displacement: float = DEFAULT_GROWTH_STEP,
) -> DrySystem:
    """Tile, insert and de-clash in one call.

    ``xy_margin`` is the clearance added on *each* side of the solute's xy
    bounding box when choosing the tiling, i.e. the minimum protein-to-image
    separation across the periodic boundary.

    ``shrink_factor`` selects the shrink-and-grow path; see
    :func:`delete_clashing_lipids`.  It defaults to ``None`` -- plain deletion,
    which needs no force field and no GROMACS.  With a shrink in force the
    returned system still has to be handed to :func:`grow_solute`.
    """
    sol = solute if isinstance(solute, Solute) else read_solute(solute)
    patch = load_patch(lipid)
    lo, hi = sol.extent()
    target_x = float(hi[0] - lo[0]) + 2.0 * xy_margin
    target_y = float(hi[1] - lo[1]) + 2.0 * xy_margin
    tiled = tile_patch(
        patch, target_x, target_y, nx=nx, ny=ny, center_xy=sol.xy_center()
    )
    system = insert_solute(tiled, sol, margin=0.0, box_z=box_z)
    return delete_clashing_lipids(
        system, clash_cutoff,
        shrink_factor=shrink_factor,
        balance_leaflets=balance_leaflets,
        max_step_displacement=max_step_displacement,
    )


# ---------------------------------------------------------------------------
# Asset derivation (provenance; run manually, not at import)
# ---------------------------------------------------------------------------

def derive_asset(
    source_gro: os.PathLike | str,
    species: Sequence[str],
    out_path: os.PathLike | str,
    *,
    drop_residues: Sequence[str] = ("WAT", "HOH", "SOL", "K+", "NA", "CL", "Na+", "Cl-"),
) -> dict:
    """Turn a deposited bilayer frame into a canonical dry PRISM patch.

    Deterministic pipeline, recorded here so the bundled asset can be
    regenerated and audited from the published source file:

    1. drop water and counter-ions (this route builds a *dry* system);
    2. segment whole molecules by Lipid21 residue pattern;
    3. wrap whole molecules into ``[0, L)`` in x and y by centroid;
    4. locate the bilayer midplane -- circular mean of lipid heavy-atom z
       (PBC-safe: two deposited frames sit on the z boundary), refined with
       phosphorus -- and shift it to z = 0, re-wrapping molecules in z;
    5. write a gzipped GRO with no velocities and a fixed gzip mtime, so the
       bytes are reproducible.

    Returns a manifest fragment describing the result.
    """
    source_gro = Path(source_gro)
    res_ids, res_names, atom_names, xyz, box = _read_gro(source_gro)

    drop = {d.upper() for d in drop_residues}
    keep_atoms = np.array([str(r).upper() not in drop for r in res_names])
    res_ids, res_names = res_ids[keep_atoms], res_names[keep_atoms]
    atom_names, xyz = atom_names[keep_atoms], xyz[keep_atoms]

    offsets, mol_species, res_local = _segment_molecules(res_ids, res_names, species)
    assembly = LipidAssembly(
        atom_names=atom_names,
        res_names=res_names,
        res_local=res_local,
        elements=_elements_from_names(atom_names),
        xyz=xyz,
        mol_offsets=offsets,
        mol_species=mol_species,
        mol_leaflet=np.zeros(mol_species.shape[0], dtype=np.int8),
        box=box,
    )
    _check_molecules_whole(assembly)

    counts = np.diff(assembly.mol_offsets)
    # (3) wrap in xy, (4) centre the midplane on z = 0 and wrap in z
    leaflet, midplane, _ = assign_leaflets(assembly)
    assembly.mol_leaflet = leaflet
    centroids = assembly.molecule_centroids()
    shift = np.zeros((assembly.n_molecules, 3))
    for dim in (0, 1):
        shift[:, dim] = -np.floor(centroids[:, dim] / box[dim]) * box[dim]
    cz = centroids[:, 2] - midplane
    shift[:, 2] = -midplane - box[2] * np.round(cz / box[2])
    assembly.xyz = assembly.xyz + np.repeat(shift, counts, axis=0)

    leaflet, midplane, thickness = assign_leaflets(assembly)
    assembly.mol_leaflet = leaflet
    if abs(midplane) > 1e-6:
        assembly.xyz[:, 2] -= midplane
        leaflet, midplane, thickness = assign_leaflets(assembly)
        assembly.mol_leaflet = leaflet
    assembly.core_half_width = _core_half_width(assembly)

    # Canonical molecule order == the module's ordering contract for one tile.
    ordered = tile_patch(assembly, nx=1, ny=1, center_xy=(box[0] / 2.0, box[1] / 2.0))

    res_seq = _global_residue_numbers(ordered)
    lines = _gro_lines(
        ordered.atom_names, ordered.res_names, res_seq, ordered.xyz, ordered.box,
        f"PRISM dry Lipid21 patch derived from {source_gro.name}",
    )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = ("\n".join(lines) + "\n").encode("utf-8")
    with open(out_path, "wb") as raw:
        # filename="" and mtime=0 keep the gzip header constant.  Without the
        # explicit empty filename, GzipFile copies fileobj.name into the header
        # and the "same" patch written to two paths differs in its bytes.
        with gzip.GzipFile(
            filename="", fileobj=raw, mode="wb", mtime=0, compresslevel=9
        ) as gz:
            gz.write(payload)

    geom = patch_geometry(ordered)
    with open(source_gro, "rb") as handle:
        src_md5 = hashlib.md5(handle.read()).hexdigest()
    return {
        "filename": out_path.name,
        "sha256": _sha256(out_path),
        "n_atoms": ordered.n_atoms,
        "n_molecules": ordered.n_molecules,
        "composition": [[k, v] for k, v in sorted(ordered.species_counts().items())],
        "box": [round(float(v), 5) for v in ordered.box],
        "source_file": source_gro.name,
        "source_md5": src_md5,
        "measured": {
            "area_per_lipid": round(geom.area_per_lipid, 3),
            "n_lower": geom.n_lower,
            "n_upper": geom.n_upper,
            "phosphate_thickness": round(geom.phosphate_thickness, 3),
            "core_half_width": round(geom.core_half_width, 3),
        },
    }
