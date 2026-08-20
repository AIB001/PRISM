#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Membrane system builder orchestration.

Pipeline (Tier-A, AMBER-consistent, reuses PRISM's amb2gmx philosophy):

    orient (OPM/PPM/preoriented)
      -> bilayer construction (PACKMOL-Memgen, or a bundled Lipid21 patch)
      -> ParmEd  AMBER prmtop/crd  ->  GROMACS .top/.gro
      -> water + ions
      -> index groups, position restraints, MDPs, run script
      -> membrane MDP protocol (semiisotropic, grouped thermostats)

The bilayer itself is built by one of two routes, selected by
``config.bilayer_method`` (see :data:`BILAYER_PACKMOL` / :data:`BILAYER_PATCH`):
PACKMOL-Memgen's rigid-body pack, or deterministic xy tiling of a bundled
pre-equilibrated patch followed by KD-tree deletion of the lipids the solute
displaces.  The patch route is dry by construction, so it runs with the GROMACS
solvation route below; it falls back to PACKMOL, with a message, whenever no
bundled patch covers the requested composition.  Index groups and position
restraints are generated after solvation on both routes, because inserting water
and ions shifts every atom number they refer to.

There are two routes through the middle of that pipeline, selected by
``config.solvate``:

``packmol`` (default)
    PACKMOL places protein, lipids, water and ions in one GENCAN optimisation.
    Validated, but water is ~99% of the molecules it has to place: a class-C
    GPCR homodimer needs ~131k waters and the pack does not finish.

``gromacs`` (opt-in fast path)
    PACKMOL packs only protein + lipids (``--nocounter`` plus a near-zero
    density solvent, see :mod:`.packmol_memgen`), tleap parametrizes that dry
    system, and water/ions are then placed by ``gmx solvate``/``gmx genion``
    (see :mod:`.solvate`).  This is ~6-8x faster end-to-end -- the remaining
    cost is lipid packing, not water -- and turns "never completes" into
    "under an hour" on the systems that motivated it.

    ``gmx solvate`` fills any void its vdW radii allow, including the gaps
    between acyl tails, so the fast path is only correct because it *also*
    inflates the carbon radius AND geometrically purges whatever water still
    lands in the hydrophobic core.  Radii alone are asymptotic (they leave
    hundreds of core waters on a freshly packed, not-yet-equilibrated bilayer),
    which is why the purge is mandatory rather than a refinement.

Everything that touches an external tool is runtime-gated with an actionable
capability report, so importing/constructing the builder never fails.  The fast
path's own modules are imported lazily for the same reason: an installation
without them still builds membranes the slow way.
"""

from __future__ import annotations

import glob
import hashlib
import inspect
import json
import os
import re
import shutil
import subprocess
import traceback
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

from . import packmol_memgen
from . import config as membrane_config
from .config import MembraneConfig, PACKMOL_PARAMETRIZED_LIPID_FFS
from .membrane_mdp import MEMBRANE_MDP_DIRNAME, write_membrane_mdps
from .orient import memembed_available, orient_protein
from .runscript import (
    INDEX_NAME,
    RUNSCRIPT_NAME,
    SYSTEM_GRO_NAME,
    TOPOL_NAME,
    write_membrane_runscript,
)

# The single AMBER->GROMACS topology writer.  ParmEd's own
# ``save(format="gromacs")`` does not carry AMBER's 1-4 scaling across faithfully
# (see prism.membrane.solvent_params), and the fix is deliberately kept in one
# module: a second copy here would be a second chance to get it wrong, and the
# two would drift.  Imported eagerly, unlike the rest of the fast path's modules,
# because every route -- PACKMOL included -- converts a topology, so there is no
# degraded mode in which builder should still convert without it.
from .solvent_params import SolventParameterError, save_gromacs_topology

# Console output goes through the shared helpers so a membrane run looks like
# every other PRISM run.  The fallback mirrors prism/utils/system/base.py: a
# missing colour module must never make the builder unimportable.
try:
    from ..utils.colors import print_success, print_warning, print_info, path, number
except ImportError:  # pragma: no cover - exercised only on a broken install
    try:
        from prism.utils.colors import print_success, print_warning, print_info, path, number
    except ImportError:
        def print_success(x, **kwargs):
            print(f"✓ {x}")

        def print_warning(x, **kwargs):
            print(f"⚠ {x}")

        def print_info(x, **kwargs):
            print(f"ℹ {x}")

        def path(x):
            return x

        def number(x):
            return x

#: Name of the membrane system directory.  This is the value fixed by the
#: canonical table in ``NAMING_CONVENTIONS.md`` section 1: ``GMX_PROLIG_`` plus
#: one mode token, exactly as for ``GMX_PROLIG_FEP`` / ``GMX_PROLIG_PMF`` /
#: ``GMX_PROLIG_MMPBSA``.  The directory holds the same shape of deliverable as
#: the soluble path (topology, coordinates, index, restraints, driver script).
#: The CLI help, the MCP tool docstring and ``prism.utils.system.base`` name it
#: in prose; keep them in step when changing it.  Every consumer that needs the
#: value imports it from here rather than repeating the literal.
MEMBRANE_SYSTEM_DIRNAME = "GMX_PROLIG_MEMB"

#: Basenames of the artifacts inside the system directory.  They are constants
#: rather than literals because three independent places have to agree on them
#: -- this module, the run-script generator and the MCP build validator -- and a
#: literal repeated per call site drifts the moment one of them is renamed.
#: The three the run script also consumes are aliased from :mod:`.runscript`
#: rather than restated, so the two modules cannot drift apart: that module
#: publishes them as the membrane file-name contract (``prism.sim`` otherwise
#: assumes the soluble ``solv_ions.gro``).
ORIENTED_PDB = "protein_oriented.pdb"
SYSTEM_GRO = SYSTEM_GRO_NAME
SYSTEM_TOP = TOPOL_NAME
SYSTEM_NDX = INDEX_NAME
#: Aliased from the module that actually writes the file, for the same reason
#: the three names above are: whoever renames it renames it once.
PACKMOL_LOG = packmol_memgen.PRISM_LOG_NAME
#: Log of a bundled-patch bilayer build.  Named for its own route so the file in
#: the system directory says which tool actually placed the lipids.
PATCH_LOG = "bilayer_patch.log"

#: Artifacts of the patch route.  The dry PDB is what tleap is handed, and it is
#: kept: it is the only record of the geometry the patch module produced before
#: parametrisation, and the difference between "the tiling was wrong" and "tleap
#: mis-parametrized a correct tiling" is not visible in anything else.
PATCH_DRY_PDB = "bilayer_patch_dry.pdb"
PATCH_PRMTOP = "bilayer_patch.prmtop"
PATCH_INPCRD = "bilayer_patch.inpcrd"
PATCH_LEAP_IN = "leap_patch.in"

#: Artifacts of the solute-preparation stage that precedes the patch route's
#: tiling (see :meth:`MembraneBuilder._prepare_patch_solute`).  The prepared
#: input and tleap's completed output are both kept: the difference between
#: them IS the set of atoms that used to appear only after the lipids had been
#: chosen, so it is the evidence that the ordering bug is gone.
PATCH_SOLUTE_INPUT_PDB = "solute_for_tleap.pdb"
PATCH_SOLUTE_COMPLETE_PDB = "solute_complete.pdb"
PATCH_SOLUTE_LEAP_IN = "leap_solute.in"
#: The completed solute as AMBER parameters, and the restrained minimisation
#: run over it.  Kept for the same reason as the PDBs: the difference between
#: ``solute_complete.pdb`` and ``solute_relaxed.pdb`` is the set of rebuilt
#: side chains that tleap placed on top of the rest of the protein.
PATCH_SOLUTE_PRMTOP = "solute_complete.prmtop"
PATCH_SOLUTE_INPCRD = "solute_complete.inpcrd"
#: The relaxed solute, in the order tleap's ``savepdb`` wrote the unit.
#:
#: That is NOT the order the deliverable is in.  A .gro carries no molecule
#: information, so ``system.gro`` has to list atoms in the order ``[ molecules ]``
#: expands -- whole molecules, i.e. connected components of the bond graph --
#: and on a disulfide-linked homodimer the two are wildly different: 9OMO's
#: first molecule is chain A's first segment bonded to chain B's first segment,
#: drawn from opposite ends of the file.  Measured, 11841 of this file's 24624
#: atoms sit at an index holding a different atom in ``system.gro``.
#:
#: Pairing this file with ``topol.top`` therefore yields a structure whose
#: "bonds" span the box -- 17.0 nm and 11.8 nm on 9OMO -- and mdrun aborts in
#: ``mshift``.  The name says which order it is in, and
#: :meth:`MembraneBuilder._solute_order_remarks` writes the same warning into
#: the file itself, because a name alone is a thing people skim.  It is kept as
#: evidence (the diff against ``solute_complete.pdb`` is the set of rebuilt side
#: chains) and as the structure the bilayer was built around, not as a
#: simulation input.
PATCH_SOLUTE_RELAXED_PDB = "solute_relaxed.leap_order.pdb"
PATCH_SOLUTE_MIN_IN = "min_solute.in"
PATCH_SOLUTE_MIN_OUT = "min_solute.out"
PATCH_SOLUTE_MIN_RST = "min_solute.rst"

#: Restrained-minimisation settings for the completed solute.
#:
#: 100 steps of steepest descent is not a relaxation of the structure; it is
#: exactly enough to undo the overlaps tleap's ideal-geometry rebuild creates.
#: Measured on 9OMO: the worst non-bonded contact goes from 0.097 A to 1.362 A
#: and every contact under 1 A disappears (121 of them), while the backbone
#: moves 0.042 A on average and 1.14 A at worst -- i.e. the fold, and the
#: bilayer frame the orientation step established, are preserved.
#: ``restraintmask`` holds the backbone rather than "the atoms that were
#: already there" because AMBER masks are expressed over names, not over
#: provenance, and holding the backbone is what protects the orientation.
#: The cutoff is 9 A because sander rejects anything smaller for a
#: non-periodic run ("unreasonably small cut for non-periodic run").
SOLUTE_MIN_CYCLES = 100
SOLUTE_MIN_CUTOFF_A = 9.0
SOLUTE_MIN_RESTRAINT_WT = 10.0

#: Separation (A) below which two non-bonded atoms of the same molecule are
#: overlapping rather than merely in contact.  A real structure has none: the
#: deposited 9OMO coordinates have zero, and tleap's completion of them has 114.
SEVERE_OVERLAP_A = 1.0

#: Backstop on how far the minimisation may move any atom.  Deliberately loose:
#: the atoms that move furthest are the ones tleap invented, and on 9OMO the
#: largest is a backbone carbonyl at a crystallographic gap that moves 7.0 A out
#: of its neighbour (mean over all atoms: 0.073 A).  The mis-ordering this used
#: to guard against is checked directly instead, by comparing the two
#: coordinate files tleap wrote from the same unit.  This only catches a
#: minimisation that has run away entirely.
SOLUTE_MIN_MAX_DISPLACEMENT_A = 25.0

#: Record of what a finished build was asked for, written into the system
#: directory.  Re-running into an existing directory compares the request
#: against this file rather than assuming that a directory holding a
#: ``system.gro`` holds *this* system.  See :meth:`MembraneBuilder._resume`.
BUILD_FINGERPRINT = "membrane_build.json"

#: Coordinates of the water-free system on the fast path, i.e. what ParmEd
#: writes before ``gmx solvate`` runs.  It is a separate name from
#: ``system.gro`` on purpose: the deliverable must never be a half-built system
#: that merely has the deliverable's name.
DRY_SYSTEM_GRO = "system_dry.gro"

#: PACKMOL-Memgen writes the protein it actually packed -- oriented by MEMEMBED,
#: translated into the box frame and renumbered -- as ``PROT<n>.pdb`` in its
#: working directory (upstream ``main.py``: ``outpdb=f"PROT{n}.pdb"``).  With
#: ``orient=memembed`` PRISM's own :data:`ORIENTED_PDB` is only the *input* it
#: handed over, still in the deposition frame, so this is where the oriented
#: coordinates really are.
PACKMOL_ORIENTED_PDB = "PROT0.pdb"

#: Where the pre-orientation input is parked when MEMEMBED's output cannot be
#: recovered.  A file named "oriented" holding un-oriented coordinates misleads
#: everyone who opens it, so it is renamed rather than left in place.
MEMEMBED_INPUT_PDB = "protein_memembed_input.pdb"

#: Solvation routes.  ``auto`` resolves to one of these before anything runs.
SOLVATE_PACKMOL = "packmol"
SOLVATE_GROMACS = "gromacs"

#: Bilayer-construction routes, selected by the configuration field named in
#: :data:`BILAYER_CONFIG_FIELD`.
#:
#: ``packmol`` (default)
#:     PACKMOL-Memgen packs every lipid as a rigid body in one GENCAN
#:     optimisation.  The validated route, and the ONLY route for a lipid with no
#:     bundled patch, so it must stay reachable.
#: ``patch``
#:     A bundled pre-equilibrated Lipid21 patch is tiled deterministically in xy,
#:     the oriented solute is inserted at z=0, and clashing lipids are deleted by
#:     KD-tree query (:mod:`prism.membrane.patch`).  Measured 54 s for a 315k-atom
#:     system against 17-30 min for PACKMOL's dry pack and 2.9 h for its wet pack.
#:     PACKMOL is not merely slower: at the physical area per lipid of 65 A^2 its
#:     rigid-body problem is geometrically unsolvable against its own 2.0 A
#:     tolerance, so on those systems it never converges at all.
#: ``auto``
#:     Use ``patch`` when a bundled patch covers the requested composition and
#:     the installation can run it; otherwise ``packmol``.
BILAYER_PACKMOL = "packmol"
BILAYER_PATCH = "patch"
BILAYER_AUTO = "auto"

#: Configuration field naming the bilayer route
#: (:attr:`~prism.membrane.config.MembraneConfig.bilayer_source`).  Read through
#: ``getattr`` so a configuration predating it resolves to the PACKMOL route,
#: which is what such a configuration can only have meant.
BILAYER_CONFIG_FIELD = "bilayer_source"

#: MDP stages without which the staged protocol cannot be run at all (``npt2``
#: is part of the protocol but a system is still usable if it is regenerated).
REQUIRED_MDP_STAGES = frozenset({"em", "nvt", "npt", "md"})

#: Shortest PDB coordinate record this module can read.  Defined in config.py so
#: this module and packmol_memgen's residue counter cannot disagree about which
#: records exist; re-exported here because that is where it has always been read
#: from.
PDB_MIN_COORDINATE_WIDTH = membrane_config.PDB_MIN_COORDINATE_WIDTH

#: GROMACS executables to look for, in the order :class:`GromacsEnvironment`
#: tries them.  Duplicated here rather than imported because constructing that
#: class prints and raises when GROMACS is absent, and an absent GROMACS must
#: only disable the fast path, not the builder.
GMX_COMMAND_CANDIDATES = ("gmx", "gmx_mpi", "gmx_d", "gmx_mpi_d")

#: Minimum gap (nm) left between the bilayer and its own periodic image in the
#: membrane plane.  Not a guess: the validated wet path delivers exactly this.
#: PACKMOL overshoots its own xy constraint box while packing (measured 85.84 A
#: of lipid inside a +-40.18 A constraint) and packmol-memgen compensates by
#: handing tleap a box 0.22 nm wider than the packed extent.  A smaller margin
#: would have head groups overlapping their image at step 0.
XY_IMAGE_MARGIN_NM = 0.22

#: G2a: how much of the requested water thickness must already be present in
#: the dry box before normalisation.  Below this the box was collapsed onto the
#: solute (tleap's ``setBox SYS vdw`` does exactly that), which is recoverable
#: -- :meth:`MembraneBuilder._normalize_box` rebuilds the box either way -- but
#: worth reporting because it means ``--tight_box`` did not take effect.
DRY_BOX_WATER_FRACTION = 1.8

#: G2a, patch route: how far (nm) the dry box may differ from the cell the
#: patch route explicitly handed tleap before it is worth reporting.  On this
#: route :data:`DRY_BOX_WATER_FRACTION` answers a question nobody asked --
#: ``insert_solute`` sizes z for the *contents* plus a fixed headroom and
#: ``_normalize_box`` sizes it for the water afterwards, so a dry box smaller
#: than the water slab is the design rather than a fault.  What can actually go
#: wrong is tleap not honouring ``set SYS box``, and that is what this measures.
#: One rounding step of .gro's 3-decimal-nm precision is 0.001 nm, so this is
#: two orders of magnitude above the noise floor and only fires on a real
#: disagreement.
TLEAP_BOX_DRIFT_TOLERANCE_NM = 0.05

#: How far (nm) a single atom may lie outside a box face in the membrane plane
#: on the patch route before the wrap is considered broken.
#:
#: ``patch.tile_patch`` wraps lipids into the tiled cell *by centroid*, so a
#: molecule straddling a face keeps its far end outside it.  That is correct --
#: it is what makes the tiling periodic -- and it means the coordinate extent is
#: legitimately wider than the period by up to one molecule.  Measured on
#: mGluR7/POPC the overhang is 0.82 nm per face in x and 0.81 nm in y, i.e. a
#: little under half an extended POPC.  This bound sits at roughly 2.5x that, so
#: it cannot fire on a correctly wrapped bilayer but still catches the failure
#: that matters: a box off by a whole tile (6 nm here), which would leave a
#: vacuum stripe exactly like the one this constant was introduced to police.
PATCH_MAX_XY_OVERHANG_NM = 2.0

#: G8: the energy gate applied to the finished system.  A correctly built
#: membrane system sits far from both limits (measured LJ(SR) per atom: 1.26
#: kJ/mol for the PACKMOL path, 3.63 for the fast path with the core purge),
#: while a scrambled topology/coordinate pairing reaches ~1e11.  The gate is
#: therefore wide enough never to fire on a merely imperfect pack.
PREFLIGHT_MAX_LJ_PER_ATOM_KJ = 50.0


class MembraneBuildError(RuntimeError):
    """Raised when a membrane build cannot produce a runnable system."""


@dataclass
class MembraneBuildResult:
    system_dir: str
    oriented_pdb: Optional[str] = None
    gro: Optional[str] = None
    top: Optional[str] = None
    ndx: Optional[str] = None  # index.ndx with SOLU/MEMB/SOLV groups (for grompp -n)
    posres: dict = field(default_factory=dict)  # stage macro -> position-restraint itp
    runscript: Optional[str] = None  # localrun.sh driving the staged protocol
    mdps: dict = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    success: bool = False
    note: str = ""
    #: Confirmations that are not actionable (orientation verified, TER records
    #: inserted).  Kept out of ``warnings`` so an empty ``warnings`` list stays
    #: usable as a health signal by CI and by the MCP layer.
    info: List[str] = field(default_factory=list)
    packmol_log: Optional[str] = None
    lipid_counts: dict = field(default_factory=dict)  # requested lipid -> residues packed
    #: Which solvation route actually ran (``packmol``/``gromacs``).  Recorded
    #: rather than inferred: ``solvate: auto`` and the automatic refusals both
    #: change it behind the user's back, and the two routes place water by
    #: completely different means, so a report that omits it is unreproducible.
    solvate_mode: Optional[str] = None
    #: Which bilayer-construction route actually ran (``packmol``/``patch``).
    #: Recorded for the same reason as ``solvate_mode``: ``auto`` and the
    #: automatic fallbacks change it behind the user's back, and the two routes
    #: place the lipids by completely different means.
    bilayer_method: Optional[str] = None
    #: Fast-path solvation tallies.  ``None`` on the PACKMOL route, where the
    #: numbers belong to PACKMOL and are not measured here.
    water_count: Optional[int] = None
    ion_counts: dict = field(default_factory=dict)  # ion moleculetype -> count
    purged_core_waters: Optional[int] = None
    retained_pore_waters: Optional[int] = None
    effective_salt_concentration: Optional[float] = None

    def raise_for_failure(self) -> "MembraneBuildResult":
        """Raise :class:`MembraneBuildError` unless the build is runnable.

        ``build()`` returns a dataclass so the MCP layer can report a partial
        result, which means a caller that only watches for exceptions would
        treat a failed build as a success.  Callers for which a partial
        scaffold is worthless -- the CLI and the Python API -- call this so the
        process exits non-zero like every other PRISM build mode.
        """
        if self.success:
            return self
        details = "".join(f"\n  - {message}" for message in self.warnings)
        raise MembraneBuildError((self.note or "Membrane build failed.") + details)


class MembraneBuilder:
    """Build an all-atom protein-in-bilayer system for GROMACS."""

    def __init__(self, config: MembraneConfig, verbose: bool = True, overwrite: bool = False):
        self.config = config
        self.verbose = verbose
        self.overwrite = overwrite
        #: The periodic cell the patch route handed tleap, in nm, or None on a
        #: route that never told tleap a box.  Recorded so ``_normalize_box``
        #: can check that the box it reads back is the one that was requested
        #: instead of applying the PACKMOL route's ``--tight_box`` heuristic to
        #: a build where PACKMOL never ran.
        self._tleap_box_nm: Optional[Tuple[float, float, float]] = None

    # ------------------------------------------------------------------ #
    def missing_capabilities(self) -> List[str]:
        """Return a list of missing-capability messages (empty if all present)."""
        missing = []
        if (self.config.lipid_ff or "").lower() not in PACKMOL_PARAMETRIZED_LIPID_FFS:
            missing.append(
                f"The automated AMBER backend cannot parametrize lipid_ff={self.config.lipid_ff!r}; "
                "use lipid17/lipid21 or provide an externally parametrized topology."
            )
        try:
            self.config.packmol_forcefield_pair()
        except ValueError as exc:
            missing.append(str(exc))
        if not packmol_memgen.is_available():
            missing.append(
                "packmol-memgen not found (provided by AmberTools). "
                "Install: conda install -c conda-forge ambertools"
            )
        if not self._parmed_available():
            missing.append(
                "ParmEd not importable (needed for AMBER->GROMACS conversion). "
                "Install: conda install -c conda-forge parmed  (ships with AmberTools)"
            )
        return missing

    def capability_report(self) -> List[str]:
        """Deprecated alias for :meth:`missing_capabilities`.

        The name reads as "a report" but the value is the list of things that
        are *absent*, so an installed environment reports ``[]``.  Kept because
        prism/builder/workflow.py and the module usage example call it.
        """
        return self.missing_capabilities()

    @staticmethod
    def _parmed_available() -> bool:
        try:
            import parmed  # noqa: F401

            return True
        except Exception:
            return False

    # ------------------------------------------------------------------ #
    #: Number of stages announced while building, for the progress lines.
    #: The fast path inserts box normalisation, solvent-parameter merging and
    #: GROMACS solvation between conversion and index generation.
    _BUILD_STAGES = 7
    _FAST_BUILD_STAGES = 10

    def build(self, protein_pdb: str, output_dir: str) -> MembraneBuildResult:
        """Run the full membrane build. Returns a :class:`MembraneBuildResult`.

        The pipeline is all-or-nothing: orientation, packing, conversion,
        restraints, MDPs and the run script are produced in a single pass, so it
        cannot be interleaved with the soluble step-by-step API.

        Raises on input that can never yield a system -- an invalid
        configuration, a structure outside the bilayer frame, or an output
        directory holding a previous build (see ``overwrite``).  External-tool
        faults are reported through the returned result instead; callers for
        which a partial scaffold is useless call
        :meth:`MembraneBuildResult.raise_for_failure`.
        """
        # Invalid scientific parameters must fail before creating an output
        # directory or launching PACKMOL-Memgen.  They are not missing-runtime
        # capabilities and must not be downgraded to a seemingly usable
        # "partial" scaffold.
        self.config.raise_for_errors()

        system_dir = os.path.join(output_dir, MEMBRANE_SYSTEM_DIRNAME)
        if self.overwrite:
            self._clear_previous_build(system_dir)
        else:
            resumed = self._resume(system_dir, output_dir, protein_pdb)
            if resumed is not None:
                self._report(resumed)
                return resumed
            self._refuse_partial_previous_build(system_dir)

        os.makedirs(system_dir, exist_ok=True)
        result = MembraneBuildResult(system_dir=system_dir)
        # NOT copied into result.warnings yet.  Some of these state which
        # bilayer route the build will take, and at this point that is a
        # composition-level guess: whether this INSTALLATION can run the patch
        # route is decided below, and copying the guess first made a build that
        # ran PACKMOL report that it had used a bundled patch.  See
        # _record_validation_warnings.
        validation_warnings = list(self.config.validation_warnings())

        # 0) Capability gate.
        missing = self.missing_capabilities()
        if missing:
            result.warnings.extend(validation_warnings)
            result.warnings.extend(missing)
            result.note = (
                "Membrane build prerequisites missing — emitted orientation + MDP scaffold only. "
                "Install the prerequisites listed in this report (or run from an AmberTools "
                "environment) to complete the bilayer build."
            )
            # Still produce what we can without external tools: orientation + MDPs.
            self._orient(protein_pdb, system_dir, result)
            self._emit_mdps(output_dir, result)
            self._report(result)
            return result

        # 0b) Choose the solvation route, and prepare everything the fast path
        #     needs, BEFORE packing.  Every reason to decline the fast path is
        #     therefore paid for in seconds instead of being discovered after a
        #     dry pack whose water can no longer be recovered.
        # ``resolved_solvate_mode`` already forces the GROMACS route when the
        # bilayer comes from a patch (that bilayer is dry, and there is no PACKMOL
        # run to place water in), so the two settings are not paired up again
        # here: a second implementation of that rule is one that can disagree.
        mode = self._resolve_solvate_mode(protein_pdb, result)
        solvent_blocks = None
        if mode == SOLVATE_GROMACS:
            solvent_blocks = self._probe_solvent_parameters(system_dir, result)
            if solvent_blocks is None:
                mode = SOLVATE_PACKMOL
        result.solvate_mode = mode
        fast = mode == SOLVATE_GROMACS
        # Chosen after the solvation route, because the patch route is only
        # available once GROMACS is placing the water.
        bilayer = self._resolve_bilayer_method(fast, result)
        result.bilayer_method = bilayer
        # Now that both routes are facts rather than intentions, the
        # configuration's own warnings can be reported without any of them
        # claiming a route this build did not take.
        self._record_validation_warnings(validation_warnings, bilayer, result)
        self._begin_stages(fast)

        # 1) Orient.
        self._progress(f"orienting the solute against the bilayer frame ({self.config.orient})")
        self._orient(protein_pdb, system_dir, result)

        # 2) Build the bilayer.  On the PACKMOL route the tool runs with its
        #    output captured, so the budget is announced: an unexplained silence
        #    of up to two hours is indistinguishable from a hang.
        try:
            try:
                pack = self._pack_bilayer(bilayer, result.oriented_pdb, system_dir, fast, result)
            except self._patch_errors() as exc:
                # The patch module's own failures (no such patch, a patch that
                # fails its checksum) are RuntimeErrors, so they used to escape
                # the handler below and abort the build -- including for
                # ``auto``, which promised only the best route available and is
                # entitled to fall back.  _patch_blockers now catches the
                # installation-level causes up front, so reaching here means the
                # patch failed once started; that is still recoverable, because
                # nothing downstream of packing has run yet.
                bilayer = self._fall_back_to_packmol(bilayer, exc, system_dir, result)
                result.bilayer_method = bilayer
                pack = self._pack_bilayer(bilayer, result.oriented_pdb, system_dir, fast, result)
        except (OSError, subprocess.SubprocessError) as exc:
            # packmol_memgen.run already turns its own launch/timeout faults
            # into a failed result; this catches what escapes around it (an
            # unwritable work dir, a killed helper).  Anything else is a bug in
            # PRISM and must not be relabelled as an external-tool failure.
            traceback.print_exc()
            result.warnings.append(f"Bilayer packing step failed: {exc}")
            self._adopt_memembed_orientation(system_dir, result)
            self._emit_mdps(output_dir, result)
            result.note = "Bilayer packing errored; emitted orientation + MDP scaffold only."
            self._report(result)
            return result

        # The log is the only record of how the bilayer was built, and the
        # packing-quality warning tells the user to read it -- so it is written
        # whether or not packing succeeded, not just on the failure path.  Each
        # route writes its own filename: a ``packmol_memgen.log`` that PACKMOL
        # never wrote would misdescribe the build to everyone who opens it.
        log_name = PATCH_LOG if bilayer == BILAYER_PATCH else PACKMOL_LOG
        result.packmol_log = self._write_packmol_log(
            system_dir, getattr(pack, "log", ""), name=log_name
        )
        # Likewise unconditional: MEMEMBED runs early enough that even a failed
        # pack usually leaves its oriented structure behind.
        self._adopt_memembed_orientation(system_dir, result)

        if bilayer != BILAYER_PATCH:
            # Parses PACKMOL's own convergence report, so it has nothing to say
            # about a bilayer that PACKMOL did not pack.
            quality = packmol_memgen.packing_quality_warning(getattr(pack, "log", ""))
            if quality:
                result.warnings.append(quality)

        if not getattr(pack, "success", False):
            builder_name = (
                "the bundled-patch bilayer route" if bilayer == BILAYER_PATCH else "packmol-memgen"
            )
            result.warnings.append(f"{builder_name} did not produce a parametrized system.")
            self._emit_mdps(output_dir, result)
            result.note = f"{builder_name} failed; see {result.packmol_log or log_name}."
            self._report(result)
            return result

        # 3) Convert AMBER -> GROMACS via ParmEd.
        self._progress("converting the AMBER topology to GROMACS (ParmEd)")
        gro = os.path.join(system_dir, SYSTEM_GRO)
        top = os.path.join(system_dir, SYSTEM_TOP)
        ndx = os.path.join(system_dir, SYSTEM_NDX)
        # On the fast path the converted coordinates are not yet the deliverable
        # -- they hold no water at all -- so they are written under their own
        # name and ``system.gro`` only appears once solvation has succeeded.
        converted_gro = os.path.join(system_dir, DRY_SYSTEM_GRO) if fast else gro
        converted = self._amber_to_gromacs(
            pack.prmtop, pack.inpcrd, top, converted_gro, result, dry=fast
        )
        if converted is not False:
            # 4-6) Water and ions.  Raises rather than degrading: a membrane
            #      system whose solvation failed a correctness guard must not be
            #      handed back looking like a build.
            if fast:
                gro = self._fast_solvate(
                    converted_gro,
                    top,
                    planes=converted if not isinstance(converted, bool) else None,
                    blocks=solvent_blocks,
                    system_dir=system_dir,
                    result=result,
                )
            result.gro, result.top = gro, top

            # 7) Thermostat-group index.  Generated from the FINAL coordinates
            #    on both routes: SOLU/MEMB/SOLV are absolute atom indices used
            #    as tc-grps by every MDP, so an index written before solvation
            #    would either be missing SOLV entirely (the dry system has no
            #    water) or silently describe a smaller system than the one
            #    grompp is handed.
            self._progress("writing the SOLU/MEMB/SOLV thermostat index")
            self._write_membrane_index_from_gro(gro, ndx, result)

            # 8) Position restraints.  The staged equilibration MDPs reference
            #    -DPOSRES/-DPOSRES_LIPID, so grompp fails outright unless the
            #    matching itp files exist and the topology includes them.  This
            #    runs on the final topology, after any solvent moleculetypes
            #    have been merged in, so the includes land in the blocks that
            #    survive to the end.
            self._progress("writing position restraints and including them from the topology")
            self._write_position_restraints(top, system_dir, result)

        # 9) Membrane MDP protocol.
        self._progress("writing the semiisotropic membrane MDP protocol")
        self._emit_mdps(output_dir, result)

        # The end-to-end gate: grompp the real deliverable and evaluate it at
        # step 0.  Only the fast path runs it, because only the fast path
        # assembles the topology and the coordinates from two different tools --
        # the failure mode it exists to catch (a scrambled correspondence that
        # grompps cleanly behind one swallowed warning) cannot arise when tleap
        # writes both.
        if fast and result.gro and result.top and result.ndx and result.mdps:
            self._preflight_energy(result.gro, result.top, result.ndx, output_dir, result)

        # 10) Run-script driving the staged protocol.  It must come last: it
        #     encodes the stage list from the MDPs and needs the index file.
        if result.top and result.gro and result.ndx and result.mdps:
            self._progress(f"writing {RUNSCRIPT_NAME}")
            try:
                result.runscript = write_membrane_runscript(
                    system_dir,
                    result.mdps,
                    gro=os.path.basename(result.gro),
                    top=os.path.basename(top),
                    ndx=os.path.basename(ndx),
                )
            except Exception as exc:  # a missing driver script is not fatal
                traceback.print_exc()
                result.warnings.append(f"Run-script generation failed: {exc}")

        result.success = self._is_complete(result)
        result.note = "Membrane system built." if result.success else "Conversion incomplete; see warnings."
        if result.success:
            # Written last, and only on success: the record says "this directory
            # holds a finished system built from that request", and a half-built
            # directory that claimed it would be reused by the next run.
            self._write_fingerprint(system_dir, protein_pdb, result)
        self._report(result)
        return result

    # ------------------------------------------------------------------ #
    @staticmethod
    def _is_complete(result: MembraneBuildResult) -> bool:
        """Whether *result* describes a system that can actually be run."""
        return bool(
            result.gro
            and result.top
            and result.ndx
            and result.runscript
            # Restraints are part of the deliverable, not a nicety: grompp
            # accepts a -DPOSRES define that nothing consumes, so a system
            # missing its restraint files runs the staged protocol completely
            # unrestrained without a single diagnostic.
            and result.posres
            and REQUIRED_MDP_STAGES.issubset(result.mdps)
            and all(os.path.isfile(mdp) for mdp in result.mdps.values())
        )

    # ------------------------------------------------------------------ #
    # Solvation route selection
    #
    # The fast path is opt-in and, once chosen, cannot be abandoned: a dry pack
    # has no water to fall back to.  So everything that can veto it is checked
    # here, before PACKMOL is launched, and a veto silently costs nothing.
    # After the pack, a failed guard is a hard error instead (see
    # ``_fast_solvate``) -- retrying on the slow path would burn hours without
    # the user knowing the science had changed underneath them.
    # ------------------------------------------------------------------ #

    def _resolve_solvate_mode(self, protein_pdb: str, result: MembraneBuildResult) -> str:
        """Return the solvation route this build will actually use.

        The composition-level half of the decision belongs to
        :class:`~prism.membrane.config.MembraneConfig`: it resolves ``auto``
        against the water estimate and downgrades a ``gromacs`` request that its
        own ``solvate_refusals()`` veto (a phosphorus-free bilayer, a lipid
        force field PACKMOL-Memgen cannot parametrize).  Those refusals are
        already in ``validation_warnings()``, which ``build`` collected.

        What is decided here is the half the configuration cannot see: whether
        this *installation* can run the fast path at all.
        """
        estimate = self._estimate_water_count(protein_pdb)
        requested = str(
            getattr(self.config, "solvate", SOLVATE_PACKMOL) or SOLVATE_PACKMOL
        ).strip().lower()
        resolver = getattr(self.config, "resolved_solvate_mode", None)
        if resolver is None:
            # A configuration without the field predates the fast path; it can
            # only have meant the route that existed then.
            mode = SOLVATE_GROMACS if requested == SOLVATE_GROMACS else SOLVATE_PACKMOL
        else:
            mode = str(resolver(estimate)).strip().lower()

        if mode != SOLVATE_GROMACS:
            self._advise_fast_path(requested, estimate, result)
            return SOLVATE_PACKMOL

        blockers = self._fast_path_blockers()
        if blockers:
            # Explicitly requested and explicitly declined: never silent, because
            # the consequence is hours of PACKMOL the user did not budget for.
            result.warnings.append(
                "membrane solvate=gromacs is not available for this build, falling back to the "
                f"{SOLVATE_PACKMOL} route (PACKMOL will place every water molecule, which is slow "
                "on large systems): " + "; ".join(blockers)
            )
            return SOLVATE_PACKMOL
        return SOLVATE_GROMACS

    # ------------------------------------------------------------------ #
    # Bilayer-construction route selection
    #
    # Decided BEFORE anything runs, exactly as the solvation route is, and for
    # the same reason: a half-built bilayer leaves a prmtop and a
    # ``bilayer_*.pdb`` in the system directory, and PACKMOL-Memgen reuses any
    # bilayer it finds there ("Packed PDB ... found. Skipping PACKMOL").  So a
    # fallback that happened *after* a failed patch build would hand PACKMOL the
    # patch route's leftovers to parametrize.  Every reason to decline the patch
    # route is therefore paid for in milliseconds up front, and a patch build
    # that fails once started is reported as a failure rather than retried.
    # ------------------------------------------------------------------ #

    #: Words that mark a configuration warning as a forward-looking claim about
    #: the bilayer route, e.g. "The bilayer will be tiled from a bundled
    #: pre-equilibrated patch ...".  Matched loosely and only used to decide
    #: whether a claim has been overtaken by events; if the wording in
    #: :mod:`~prism.membrane.config` changes so that nothing matches, the claim
    #: is reported as before and ``_resolve_bilayer_method``'s own message still
    #: states the route that actually ran alongside it.  The failure mode is
    #: therefore a redundant line, never a silently wrong one.
    _PATCH_CLAIM_MARKERS = ("patch", "bilayer")

    @classmethod
    def _claims_patch_route(cls, message: str) -> bool:
        lowered = str(message).lower()
        return (
            all(marker in lowered for marker in cls._PATCH_CLAIM_MARKERS)
            and "will be" in lowered
        )

    def _record_validation_warnings(
        self, warnings: Sequence[str], bilayer: str, result: MembraneBuildResult
    ) -> None:
        """Copy the configuration's warnings, minus any route claim now false.

        :meth:`MembraneConfig.validation_warnings` answers from the request
        alone, so when it says the bilayer will come from a bundled patch it is
        describing what the composition allows, not what this installation can
        do.  The installation-level half is decided here, and a build that fell
        back to PACKMOL must not carry a warning saying it used a patch -- that
        line is the provenance a reader trusts.  So the claim is only copied
        once the route has been resolved and still agrees with it.
        """
        for message in warnings:
            if bilayer != BILAYER_PATCH and self._claims_patch_route(message):
                # _resolve_bilayer_method has already recorded why, with the
                # installation-level reason this message could not know.
                continue
            result.warnings.append(message)

    def _requested_bilayer_method(self) -> str:
        """The bilayer route asked for, normalised. Defaults to PACKMOL."""
        value = getattr(self.config, BILAYER_CONFIG_FIELD, None)
        return str(value).strip().lower() if value else BILAYER_PACKMOL

    def _resolve_bilayer_method(self, fast: bool, result: MembraneBuildResult) -> str:
        """Return the bilayer route this build will actually use.

        The composition-level half of the decision belongs to
        :class:`~prism.membrane.config.MembraneConfig`: ``resolved_bilayer_source``
        resolves ``auto`` against ``bilayer_refusals()`` (no bundled patch for the
        composition, a phosphorus-free bilayer whose leaflets could not be
        identified), and ``validation_errors()`` has already refused an explicit
        ``patch`` request those rule out.  Those refusals are reported through
        ``validation_warnings()``, which ``build`` collected.

        What is decided here is the half the configuration cannot see: whether
        this *installation* can run the patch route at all.

        ``fast`` says whether GROMACS is placing the water.  The patch route
        hands back a dry system by construction, so the pairing with PACKMOL
        solvation would leave the system with no water at all; the configuration
        rejects that combination outright, so this is a backstop rather than the
        primary check.
        """
        requested = self._requested_bilayer_method()
        resolver = getattr(self.config, "resolved_bilayer_source", None)
        if resolver is None:
            # A configuration without the field predates the patch route; it can
            # only have meant the route that existed then.
            resolved = BILAYER_PATCH if requested == BILAYER_PATCH else BILAYER_PACKMOL
        else:
            resolved = str(resolver()).strip().lower()
        if resolved != BILAYER_PATCH:
            return BILAYER_PACKMOL

        blockers = self._patch_blockers(fast)
        if not blockers:
            return BILAYER_PATCH
        if requested == BILAYER_AUTO:
            # Nothing was promised, so this is information rather than a
            # correction: ``auto`` asked for the best available route and got it.
            result.info.append(
                f"Bilayer route: PACKMOL-Memgen (the {BILAYER_PATCH} route is not usable here: "
                + "; ".join(blockers)
                + ")."
            )
            return BILAYER_PACKMOL
        # Explicitly requested and explicitly declined: never silent, because the
        # consequence is minutes-to-hours of PACKMOL the user did not budget for,
        # and on a composition PACKMOL cannot converge on it is a build that never
        # finishes rather than a slow one.
        result.warnings.append(
            f"membrane {BILAYER_CONFIG_FIELD}={BILAYER_PATCH} is not available for this build, "
            f"falling back to the {BILAYER_PACKMOL} bilayer route (PACKMOL packs every lipid as a "
            "rigid body, which is far slower and may not converge at a physical area per lipid): "
            + "; ".join(blockers)
        )
        return BILAYER_PACKMOL

    def _patch_blockers(self, fast: bool) -> List[str]:
        """Reasons the bundled-patch bilayer route cannot be used here."""
        blockers: List[str] = []
        try:
            module = self._patch_module()
        except ImportError as exc:
            return [f"the patch-route module is unavailable ({exc})"]

        if self._first_attr(module, self._PATCH_BUILD_ENTRIES) is None:
            blockers.append(
                "prism.membrane.patch exposes none of "
                + "/".join(self._PATCH_BUILD_ENTRIES)
                + ", so there is no entry point to build a bilayer with"
            )
        # Appends its own reason when the answer is no.
        self._patch_assets_installed(module, blockers)
        if (self.config.orient or "").lower() == "memembed" and not memembed_available():
            # PRISM now runs MEMEMBED itself (orient.run_memembed) instead of
            # letting packmol-memgen do it out of sight, so this route is only
            # blocked when the binary is genuinely missing.  The old blanket
            # block was backwards: structures needing automatic orientation are
            # the ones absent from OPM -- new depositions -- so the fast route
            # was barred from precisely the cases that most need it.
            blockers.append(
                "orient=memembed needs the 'memembed' binary on PATH (it ships with AmberTools); "
                "install it, or orient the structure externally and use orient=preoriented"
            )
        if not fast:
            blockers.append(
                "the patch route produces a water-free system, so it needs the GROMACS solvation "
                f"route (set membrane solvate={SOLVATE_GROMACS}); this build is solvating with "
                f"{SOLVATE_PACKMOL}"
            )
        return blockers

    #: Entry point looked for on :mod:`prism.membrane.patch`, in order.  Probed
    #: rather than imported by name so a rename there degrades to PACKMOL with a
    #: message instead of raising AttributeError in the middle of a build.
    _PATCH_BUILD_ENTRIES = ("build_dry_membrane_system", "build_dry_system", "build")

    def _patch_assets_installed(self, module, blockers: List[str]) -> bool:
        """Whether this installation actually carries the bundled patch files.

        Only the installation-level question is asked here.  WHICH compositions
        have a patch is :meth:`MembraneConfig.bilayer_refusals`'s answer, and it
        is documented as the single source of truth for it -- a second
        implementation of the same rule is one that can disagree.  What the
        configuration cannot know is whether the manifest those names refer to was
        shipped, which is a real failure mode: the patches are package data, and a
        wheel built without it imports fine and has no patches.

        Still cheap, but no longer manifest-only.  The manifest is an index, and
        an index is not the thing it indexes: a wheel built without package data
        can ship ``manifest.json`` and none of the ``.gro.gz`` coordinate files
        it names, in which case this returned True, ``auto`` chose the patch
        route, and ``load_patch`` then raised ``PatchNotAvailableError`` /
        ``PatchIntegrityError`` in the middle of the build.  Those are
        RuntimeErrors, so the ``(OSError, SubprocessError)`` handler around
        packing did not catch them and ``auto`` -- which promised only "the best
        route available" -- aborted the build instead of using PACKMOL.

        So the coordinate file backing the requested patch is checked too, by
        existence, size and gzip header.  That is still milliseconds and never
        decompresses anything, and it is what the manifest cannot answer.
        """
        query = getattr(module, "available_patches", None)
        if query is None:
            return True  # nothing to interrogate; let the build report the truth
        try:
            specs = query()
        except Exception as exc:
            # A missing manifest lands here and its own message is the actionable
            # one (package data not installed, or the override directory is
            # wrong), so it is passed through rather than reworded.
            blockers.append(f"the bundled patch manifest could not be read ({exc})")
            return False
        if not specs:
            blockers.append("this installation carries no bundled bilayer patches")
            return False

        # WHICH composition has a patch is MembraneConfig.bilayer_refusals'
        # answer and is not re-derived here; this only asks whether the file
        # behind the name it chose was actually shipped.
        resolver = getattr(self.config, "bundled_patch", None)
        wanted = (resolver() if callable(resolver) else None) or ""
        if not wanted:
            return True
        spec = self._patch_spec(specs, str(wanted))
        if spec is None:
            blockers.append(
                f"the bundled patch manifest names no entry for {wanted!r} "
                f"(it has: {', '.join(sorted(specs))})"
            )
            return False
        try:
            directory = module.patch_data_dir()
        except Exception as exc:
            blockers.append(f"the bundled patch directory could not be located ({exc})")
            return False
        coordinates = os.path.join(str(directory), str(getattr(spec, "filename", "")))
        problem = self._unusable_patch_file(coordinates)
        if problem:
            blockers.append(
                f"the bundled {wanted} patch is indexed by the manifest but its coordinate file "
                f"{problem}"
            )
            return False
        return True

    @staticmethod
    def _patch_spec(specs, wanted: str):
        """The manifest entry for *wanted*, matched by key or by alias."""
        target = wanted.strip().upper()
        for key, spec in specs.items():
            if str(key).upper() == target:
                return spec
            if target in {str(a).upper() for a in getattr(spec, "aliases", ()) or ()}:
                return spec
        return None

    #: gzip's magic number.  Checked instead of decompressing: a truncated or
    #: git-lfs-pointer stand-in for a patch fails here in microseconds, and a
    #: file that passes is left for the patch module to verify by checksum.
    _GZIP_MAGIC = b"\x1f\x8b"

    @classmethod
    def _unusable_patch_file(cls, path: str) -> Optional[str]:
        """Why *path* cannot be a bundled patch, or None if it looks like one."""
        if not os.path.isfile(path):
            return f"is not installed ({path}); reinstall from source with 'pip install -e .'"
        try:
            if os.path.getsize(path) == 0:
                return f"is empty ({path})"
            if path.endswith(".gz"):
                with open(path, "rb") as handle:
                    if handle.read(2) != cls._GZIP_MAGIC:
                        return f"is not gzip data ({path}); it may be a git-lfs pointer"
        except OSError as exc:
            return f"could not be read ({path}: {exc})"
        return None

    @staticmethod
    def _first_attr(module, names):
        """First attribute of *module* named in *names* that is callable."""
        for name in names:
            candidate = getattr(module, name, None)
            if callable(candidate):
                return candidate
        return None

    @classmethod
    def _patch_errors(cls):
        """The patch module's own exception base, as a tuple for ``except``.

        Empty when the module or its error class cannot be imported, which
        makes the ``except`` clause match nothing -- the correct behaviour when
        there is no patch route to fail in the first place.
        """
        try:
            base = getattr(cls._patch_module(), "MembranePatchError", None)
        except ImportError:
            return ()
        return (base,) if isinstance(base, type) and issubclass(base, BaseException) else ()

    def _fall_back_to_packmol(self, bilayer: str, exc, system_dir: str, result) -> str:
        """Switch to PACKMOL after the patch route failed. Returns the new route.

        Only reachable before PACKMOL-Memgen has run, which is what makes the
        switch safe: packmol-memgen reuses any bilayer it finds in its working
        directory ("Packed PDB ... found. Skipping PACKMOL"), so the patch
        route's leftovers are cleared first rather than handed to it to
        parametrize.  An explicit ``bilayer_source=patch`` is still refused --
        falling back silently from a route the user named is how a build ends up
        being something other than what was asked for -- but it is refused with
        the patch module's own message rather than a traceback.
        """
        if self._requested_bilayer_method() != BILAYER_AUTO:
            raise MembraneBuildError(
                f"The {BILAYER_PATCH} bilayer route failed: {exc}\nIt was requested explicitly, so "
                f"PRISM is not silently substituting the {BILAYER_PACKMOL} route. Rebuild with "
                f"--membrane-{BILAYER_CONFIG_FIELD.replace('_', '-')} {BILAYER_PACKMOL} to use it."
            ) from exc
        for entry in self._matching_artifacts(system_dir, self._BLOCKING_ARTIFACT_GLOBS):
            try:
                os.remove(entry)
            except OSError as failure:
                raise MembraneBuildError(
                    f"The {BILAYER_PATCH} bilayer route failed ({exc}) and its leftovers could not "
                    f"be cleared before falling back ({entry}: {failure}); packmol-memgen would "
                    "reuse them instead of packing a new bilayer."
                ) from exc
        result.warnings.append(
            f"bilayer_source=auto selected the {BILAYER_PATCH} route, which then failed ({exc}); "
            f"rebuilt the bilayer with {BILAYER_PACKMOL} instead. PACKMOL packs every lipid as a "
            "rigid body, so this build is far slower and may not converge at a physical area per "
            "lipid."
        )
        return BILAYER_PACKMOL

    @staticmethod
    def _patch_module():
        """Import :mod:`prism.membrane.patch`, or raise ImportError.

        Imported lazily for the same reason as the fast path's own modules: an
        installation without it must still build membranes the PACKMOL way
        instead of making :mod:`prism.membrane.builder` unimportable.
        """
        from . import patch as patch_module

        return patch_module

    def _pack_bilayer(
        self,
        bilayer: str,
        oriented_pdb: str,
        system_dir: str,
        fast: bool,
        result: "MembraneBuildResult",
    ):
        """Build the bilayer by the selected route. Returns a pack result.

        Both routes return the same duck-typed shape -- ``prmtop``, ``inpcrd``,
        ``log``, ``success`` -- so everything downstream (ParmEd conversion,
        solvation, index, restraints, MDPs, run script) is identical either way
        and neither route is a special case of the other.
        """
        if bilayer != BILAYER_PATCH:
            self._progress(
                "packing the bilayer"
                + (" (protein + lipids only; water comes from GROMACS)" if fast else "")
                + " with PACKMOL-Memgen "
                f"(silent until it finishes; budget {self.config.packmol_timeout:.0f} s)",
            )
            return packmol_memgen.run(
                oriented_pdb,
                self.config,
                work_dir=system_dir,
                timeout=self.config.packmol_timeout,
                # Only passed when requested: the backend's ``dry`` support is
                # verified in _resolve_solvate_mode, and the slow path must keep
                # calling exactly the signature it always called.
                **({"dry": True} if fast else {}),
            )

        self._progress(
            "building the bilayer from a bundled pre-equilibrated Lipid21 patch "
            "(xy tiling, solute insertion, clash deletion; water comes from GROMACS)"
        )
        module = self._patch_module()
        build = self._first_attr(module, self._PATCH_BUILD_ENTRIES)
        if build is None:
            raise MembraneBuildError(
                "prism.membrane.patch exposes none of "
                + "/".join(self._PATCH_BUILD_ENTRIES)
                + f", so the {BILAYER_PATCH} bilayer route cannot be run. Rebuild with "
                f"--membrane-{BILAYER_CONFIG_FIELD.replace('_', '-')} {BILAYER_PACKMOL}."
            )

        # The patch identifier comes from the configuration, which owns the
        # composition -> patch mapping; it is resolved against the bundled
        # manifest's aliases by the patch module itself.  Falling back to the
        # lipid name keeps this working against a configuration that predates
        # ``bundled_patch()``.
        resolver = getattr(self.config, "bundled_patch", None)
        lipid = (resolver() if callable(resolver) else None) or sorted(self._requested_lipids())[0]

        # Optional arguments are passed only when the callable accepts them, the
        # same tolerance this module already applies to the solvator: the patch
        # module is a separate change and must be able to grow its signature
        # without this call site becoming wrong.  ``box_z`` is deliberately NOT
        # set: _normalize_box rebuilds the z box around the requested water
        # thickness after conversion, and two places sizing the same axis is how
        # they come to disagree.
        accepted = self._callable_parameters(build)
        extra = {}
        for name, value in (
            ("lipid", lipid),
            ("xy_margin", getattr(self.config, "patch_xy_margin", None)),
            ("clash_cutoff", getattr(self.config, "patch_clash_distance", None)),
            ("seam_clash_cutoff", getattr(self.config, "patch_seam_clash_distance", None)),
            # Shrink-and-grow (Wolf et al. / OpenMM).  Off unless asked for; see
            # _refuse_ungrown_system for why, and for what enabling it costs.
            ("shrink_factor", self._requested_shrink_factor()),
        ):
            if name in accepted and value is not None:
                extra[name] = value

        # The solute is completed by tleap BEFORE the lipids are chosen.  This
        # ordering is the whole point; see _prepare_patch_solute.
        solute_pdb, solute_log = self._prepare_patch_solute(oriented_pdb, system_dir, result)
        dry = build(solute_pdb, **extra)
        # A system that still carries a growth plan is an intermediate, not a
        # bilayer.  Checked here, before tleap spends four minutes on it.
        self._refuse_ungrown_system(dry)

        log = [
            f"PRISM fast membrane route: bundled pre-equilibrated {lipid} patch",
            f"solute:            {os.path.basename(oriented_pdb)}",
        ]
        log.extend(solute_log)
        composition = getattr(dry, "composition", None)
        if callable(composition):
            try:
                log.append(f"composition:       {composition()}")
            except Exception as exc:  # provenance is not worth failing a build
                log.append(f"composition:       unavailable ({exc})")
        deletion = getattr(dry, "deletion", None)
        if deletion is not None:
            log.append(f"clash deletion:    {deletion}")

        # The patch module returns geometry, not parameters, so this route runs
        # tleap itself.  PACKMOL-Memgen does the equivalent behind
        # ``--parametrize``; from here on the two routes are indistinguishable.
        dry_pdb = os.path.join(system_dir, PATCH_DRY_PDB)
        dry.write_pdb(dry_pdb)
        self._progress("parametrizing the dry patch system with tleap")
        prmtop, inpcrd, leap_log = self._parametrize_with_tleap(
            dry_pdb, self._patch_box(dry), system_dir, result
        )
        log.append(leap_log)
        return packmol_memgen.PackmolMemgenResult(
            prmtop=prmtop,
            inpcrd=inpcrd,
            packed_pdb=dry_pdb,
            log="\n".join(log),
            success=True,
            dry=True,
        )

    def _requested_shrink_factor(self):
        """``patch_shrink_factor`` from the configuration, or None.

        None means plain deletion, and that is the default -- deliberately.

        Shrink-and-grow (Wolf et al., the algorithm ``openmm.app.Modeller``
        uses) deletes lipids against a laterally shrunken copy of the solute and
        then grows the solute back, minimising the bilayer at each increment, so
        the survivors close in around it instead of being deleted out to the
        solute's full footprint.  Measured on mGluR7 it is a real improvement:
        the rim gap falls from 7.56 to 4.29 A, lipids in contact go from 13 to
        49, and the area per lipid *improves* from 64.02 to 60.62 A^2, toward
        the bundled patch's own 60.43.

        It is off by default on cost.  The grow-back is 14 GROMACS
        minimisations -- 13 short ones plus one full PME run -- and was measured
        at 821 s on 16 cores, against a ~616 s build and a stated ceiling of
        five minutes.  Turning it on by default would take a build that already
        sits at 2x the ceiling to roughly 4x it, for every user, including the
        ones who only want a system to look at.

        The case for on-by-default would be that the system is unusable without
        it, and after the box fix in :meth:`_xy_box` that is not so.  The rim
        gap is a local annulus that equilibration closes; the vacuum stripe the
        box bug left was a global fault that no equilibration recovers from,
        and it -- not the rim -- is what put the system at -5325 bar.  With the
        box correct, plain deletion delivers a bilayer that is physically sound
        at the start and merely takes longer to reach its equilibrium packing.
        That is an equilibration-time cost the user can choose to pay up front,
        which is exactly what an opt-in is.
        """
        value = getattr(self.config, "patch_shrink_factor", None)
        if value is None:
            return None
        if isinstance(value, str):
            text = value.strip().lower()
            if text in ("", "none", "off", "false", "no"):
                return None
            if text == "auto":
                return text  # patch.py solves for the factor itself
            try:
                # YAML quotes numbers often enough that passing "0.7" straight
                # through would fail in patch.py with a type complaint about a
                # value the user did spell correctly.
                return float(text)
            except ValueError:
                raise MembraneBuildError(
                    f"membrane patch_shrink_factor={value!r} is not a number, 'auto', or 'none'."
                ) from None
        return value

    def _refuse_ungrown_system(self, dry) -> None:
        """Refuse a system whose grow-back has not run. Returns None otherwise.

        ``patch.py`` builds the shrink-and-grow system correctly, but PRISM
        cannot yet finish it, and the honest thing is to say so before tleap
        spends four minutes producing a topology for coordinates that must not
        be shipped.  Two things have to be true before the grow-back can run
        here, and on this route only the first is:

        1. *Ordering.*  The relaxer minimises the bilayer against a GROMACS
           topology, so the topology has to exist before the grow-back rather
           than after it.  That much is fine -- ``grow_solute`` never moves the
           solute and tleap infers no bonds (PRISM declares every disulfide
           explicitly, from solute geometry that the grow-back leaves
           untouched), so a topology built before the grow is the same topology
           that describes the system after it.

        2. *Atom order.*  This is the blocker.  ``patch.py`` writes the relaxer's
           coordinates with the solute in the order it read it, while the
           topology PRISM hands the relaxer is in ``[ molecules ]`` order --
           whole molecules, i.e. connected components of the bond graph.  On a
           disulfide-linked homodimer those differ: measured on 9OMO, 11841 of
           the solute's 24624 atoms sit at an index holding a different atom.
           Pairing them gives the same multi-nanometre bonds Defect 3 is about,
           and mdrun aborts in ``mshift`` -- so the minimisation the grow-back
           depends on cannot run at all.

        Closing (2) means feeding ``patch.py`` a solute that is already in
        topology order, which changes the atom order of every patch build's
        deliverable, not just the ones that ask for a shrink.  That is a bigger
        change than this one, and it is not worth making blind for a feature
        that is off by default and costs 821 s when on.

        Raising is the right failure: the alternative is a system with 0.27 A
        lipid-solute overlaps, which is exactly what ``DrySystem.write_pdb``
        refuses to emit and what this check exists to explain rather than merely
        trip over.
        """
        if not getattr(dry, "growth", None):
            return
        factor = getattr(dry.growth, "shrink_factor", None)
        stated = f"{factor:.3f}" if isinstance(factor, (int, float)) else str(factor)
        raise MembraneBuildError(
            f"membrane patch_shrink_factor={stated} selected the shrink-and-grow lipid "
            "packing, but PRISM cannot complete it yet: the bilayer minimisation it relies on "
            "needs the solute in the same atom order as the topology, and on this route "
            "patch.py emits the solute in input order while the topology is in [ molecules ] "
            "order (on a disulfide-linked homodimer these differ in roughly half the solute's "
            "atoms, which makes the minimisation abort in mshift). The lipids are currently "
            f"placed against a shrunken solute and overlap the real one by "
            f"{self._closest_full_size_contact(dry)}, so this system must not be written. "
            "Rebuild without patch_shrink_factor (plain clash deletion): with the periodic box "
            "corrected that gives a physically sound bilayer, at a wider lipid-protein rim that "
            "equilibration closes."
        )

    @staticmethod
    def _closest_full_size_contact(dry) -> str:
        """How close the survivors come to the full-size solute, for the message."""
        distance = getattr(getattr(dry, "deletion", None), "min_heavy_distance_full", None)
        if not isinstance(distance, (int, float)):
            return "an unreported margin"
        return f"{distance:.2f} A at the closest heavy-atom pair"

    @staticmethod
    def _patch_box(dry) -> Tuple[float, float, float]:
        """The dry system's box in Angstrom, as tleap wants it."""
        box = getattr(dry, "box", None)
        if box is None or len(box) < 3:
            raise MembraneBuildError(
                "The patch route returned a system without a box, so tleap cannot be told what "
                "periodic cell to parametrize."
            )
        return float(box[0]), float(box[1]), float(box[2])

    # ------------------------------------------------------------------ #
    # Solute preparation for the patch route
    #
    # The patch route decides which lipids the solute displaces by measuring
    # the solute's coordinates.  tleap then completes that solute -- unresolved
    # side-chain tips, every hydrogen -- and those atoms appear AFTER the
    # lipids around them were chosen, at ideal template geometry, inside space
    # the deletion had just declared empty.  Measured on 9OMO before this fix:
    # 1130 heavy atoms and 2197 hydrogens added after deletion, minimum
    # protein-lipid separation 0.395 A, 6 pairs under 2.0 A, all of them
    # rebuilt aromatic side chains driven into acyl tails.
    #
    # Three fixes were available: re-run the deletion against the completed
    # structure and re-parametrize; complete the solute before the deletion; or
    # delete a second time and rebuild only the topology.  Completing first is
    # implemented here for three reasons.
    #
    #   * It is the only one where the deletion is evaluated against exactly
    #     the coordinates that end up in the topology, rather than against an
    #     approximation of them that is corrected afterwards.
    #   * It is the cheapest.  A re-deletion needs a SECOND tleap over the
    #     whole bilayer, because the topology has to match the new lipid set;
    #     completing first needs one tleap over the solute alone, which is a
    #     fifth of the atoms.  The topology-only variant is cheaper still and
    #     was rejected: it means hand-editing a prmtop and trusting the edit,
    #     with no tool left to disagree with.
    #   * It is checkable.  If the ordering worked, the bilayer's own tleap run
    #     adds zero atoms -- which _assert_solute_was_already_complete asserts
    #     rather than assumes, so the fix cannot rot back into the bug.
    #
    # Widening the clash cutoff instead was rejected on the measurements: the
    # atoms in question move several Angstrom, so a cutoff wide enough to cover
    # them only widens the vacuum annulus around the protein.
    # ------------------------------------------------------------------ #

    def _prepare_patch_solute(self, oriented_pdb: str, system_dir: str, result):
        """Complete the solute with tleap. Returns (pdb, log lines).

        The returned structure carries every atom the final topology will
        carry, in the frame the orientation step put it in, so the patch
        module's tiling and clash deletion see the finished solute.
        """
        self._progress(
            "completing the solute with tleap (missing heavy atoms, hydrogens, disulfides) "
            "before the lipids are chosen"
        )
        protein_ff, water_ff = self.config.packmol_forcefield_pair()
        prepared = os.path.join(system_dir, PATCH_SOLUTE_INPUT_PDB)
        complete = os.path.join(system_dir, PATCH_SOLUTE_COMPLETE_PDB)

        stripped, altlocs, protonation, labels = self._write_tleap_ready_solute(
            oriented_pdb, prepared
        )
        bonds, described = self._disulfide_bond_lines(prepared, rename=True, labels=labels)

        lines = [
            f"source leaprc.protein.{protein_ff}",
            f"source leaprc.water.{water_ff}",
        ]
        lines.extend(self._tleap_ligand_lines())
        prmtop = os.path.join(system_dir, PATCH_SOLUTE_PRMTOP)
        inpcrd = os.path.join(system_dir, PATCH_SOLUTE_INPCRD)
        lines.append(f"SYS = loadpdb {os.path.basename(prepared)}")
        lines.extend(bonds)
        # Parameters as well as coordinates: the completed solute has to be
        # relaxed before the lipids are chosen around it, and a minimiser needs
        # a topology.  It is the same unit in the same run, so the two files
        # cannot describe different structures.
        lines.append(
            f"saveAmberParm SYS {os.path.basename(prmtop)} {os.path.basename(inpcrd)}"
        )
        lines.append(f"savepdb SYS {os.path.basename(complete)}")
        lines.append("quit")

        output = self._run_tleap(
            lines,
            PATCH_SOLUTE_LEAP_IN,
            system_dir,
            produces=(complete, prmtop, inpcrd),
            doing="completing the solute before the bilayer is built around it",
        )
        heavy, hydrogens = self._leap_added_atoms(output)
        complete, relaxation = self._relax_completed_solute(
            complete, prmtop, inpcrd, system_dir, result
        )

        log = [
            f"solute prepared:   {os.path.basename(complete)} "
            f"(tleap completed it BEFORE clash deletion)",
            f"  hydrogens stripped from the input: {stripped} (tleap rebuilt them)",
            f"  heavy atoms tleap added:           {heavy}",
            f"  hydrogens tleap added:             {hydrogens}",
            f"  disulfide bonds declared:          {len(bonds)}",
        ]
        log.extend(f"    {entry}" for entry in described)
        log.extend(relaxation)

        result.info.append(
            f"Solute completed by tleap before the bilayer was built around it: "
            f"{heavy} heavy atom(s) and {hydrogens} hydrogen(s) added, {len(bonds)} disulfide(s) "
            f"declared. The clash deletion therefore measured the finished solute, not the input."
        )
        if stripped:
            result.info.append(
                f"Stripped {stripped} input hydrogen(s) before handing the solute to tleap, which "
                "rebuilds every hydrogen from its own templates; histidine protonation was carried "
                "over into the residue names first"
                + (f" ({protonation})." if protonation else ".")
            )
        if altlocs:
            result.warnings.append(
                f"The solute carries alternate conformations; kept the 'A' conformer and dropped "
                f"{altlocs} atom(s) from the others, because tleap would otherwise read the same "
                "atom twice."
            )
        return complete, log

    def _relax_completed_solute(self, complete, prmtop, inpcrd, system_dir, result):
        """Undo the overlaps tleap's rebuild leaves inside the solute.

        Completing a structure fixes where the lipids go, but it does not fix
        the protein: tleap places a missing side chain at ideal internal
        coordinates without looking at what is already there, so a side chain
        the crystallographer could not resolve is rebuilt straight through its
        neighbours.  Measured on 9OMO, whose deposited coordinates contain no
        non-bonded contact under 1 A at all, tleap's completion contains 114 --
        the worst being ``ILE449:CA -- ARG459:HH11`` at 0.097 A.  That is a
        protein-internal fault, so no amount of lipid deletion touches it, and
        it alone puts the finished system at +1.8e16 kJ/mol.

        A short restrained steepest descent removes them: the rebuilt tips have
        enormous forces on them and move first, while the backbone restraint
        holds the fold -- and, just as importantly, the membrane frame the
        orientation step established.  It runs on the solute ALONE, before the
        bilayer exists, which is both the cheapest place to do it (a fifth of
        the atoms, no solvent) and the only place where relaxing the protein
        cannot invalidate the lipid deletion, because the deletion has not
        happened yet.

        Skipped, not attempted, when the completed structure has no overlap to
        remove: that is the common case for a well-resolved input, and it is
        worth two seconds to find out instead of a minute and a half to fix
        nothing.
        """
        overlaps = self._severe_internal_overlaps(complete)
        if overlaps == 0:
            result.info.append(
                "Solute relaxation skipped: tleap's completed structure has no non-bonded contact "
                f"under {SEVERE_OVERLAP_A:.1f} A, so there is nothing to relieve."
            )
            return complete, ["  solute relaxation:                 not needed (no overlaps)"]
        # ``None`` means the check could not run (no SciPy).  Relaxing anyway is
        # the safe direction: it costs a minute, while skipping it costs the
        # build.
        found = "an unknown number of" if overlaps is None else str(overlaps)

        sander = self._sander_command()
        if sander is None:
            result.warnings.append(
                "sander was not found, so the solute's rebuilt side chains could not be relaxed. "
                f"tleap left {found} non-bonded contact(s) under {SEVERE_OVERLAP_A:.1f} A inside "
                "the protein, which will show up as an astronomical starting energy. Install "
                "AmberTools (conda install -c conda-forge ambertools)."
            )
            return complete, ["  solute relaxation:                 skipped (sander not found)"]

        control = os.path.join(system_dir, PATCH_SOLUTE_MIN_IN)
        with open(control, "w") as handle:
            handle.write(
                "PRISM: relieve the overlaps tleap's rebuilt side chains create\n"
                " &cntrl\n"
                f"  imin=1, maxcyc={SOLUTE_MIN_CYCLES:d}, ntmin=2,\n"
                "  ntb=0, igb=0, "
                f"cut={SOLUTE_MIN_CUTOFF_A:g},\n"
                f"  ntr=1, restraint_wt={SOLUTE_MIN_RESTRAINT_WT:g}, restraintmask='@CA,C,N,O',\n"
                "  ntpr=50, ntxo=1,\n"
                " /\n"
            )
        restart = os.path.join(system_dir, PATCH_SOLUTE_MIN_RST)
        for stale in (restart,):
            if os.path.isfile(stale):
                os.remove(stale)
        self._progress(
            f"relaxing {found} overlap(s) tleap's rebuilt side chains left in the solute "
            f"({SOLUTE_MIN_CYCLES} restrained steps)"
        )
        completed = self._run_gmx(
            [
                sander, "-O",
                "-i", os.path.basename(control),
                "-p", os.path.basename(prmtop),
                "-c", os.path.basename(inpcrd),
                "-ref", os.path.basename(inpcrd),
                "-o", PATCH_SOLUTE_MIN_OUT,
                "-r", os.path.basename(restart),
                "-x", os.devnull,
            ],
            system_dir,
            timeout=1800,
        )
        if completed is None or not os.path.isfile(restart):
            result.warnings.append(
                "The restrained minimisation of the solute did not run, so tleap's rebuilt side "
                f"chains still overlap the rest of the protein ({found} contact(s) under "
                f"{SEVERE_OVERLAP_A:.1f} A). See {PATCH_SOLUTE_MIN_OUT}."
            )
            return complete, ["  solute relaxation:                 failed (see " + PATCH_SOLUTE_MIN_OUT + ")"]

        relaxed = os.path.join(system_dir, PATCH_SOLUTE_RELAXED_PDB)
        moved, shifted = self._apply_restart_coordinates(
            complete, restart, inpcrd, relaxed,
            remarks=self._solute_order_remarks(self._count_pdb_atoms(complete)),
        )
        remaining = self._severe_internal_overlaps(relaxed)
        result.info.append(
            f"Solute relaxed before the bilayer was built around it: {found} non-bonded "
            f"contact(s) under {SEVERE_OVERLAP_A:.1f} A in tleap's rebuilt structure, "
            f"{remaining if remaining is not None else '?'} after {SOLUTE_MIN_CYCLES} restrained "
            f"steps ({shifted} atom(s) moved more than 2 A, largest {moved:.2f} A; these are the "
            "atoms tleap invented, and the backbone restraint held the fold and the bilayer frame)."
        )
        if remaining:
            result.warnings.append(
                f"{remaining} non-bonded contact(s) under {SEVERE_OVERLAP_A:.1f} A remain inside "
                f"the solute after relaxation; the starting energy will be high. Inspect "
                f"{PATCH_SOLUTE_RELAXED_PDB}."
            )
        return relaxed, [
            f"  solute relaxation:                 {found} -> "
            f"{remaining if remaining is not None else '?'} overlaps under "
            f"{SEVERE_OVERLAP_A:.1f} A, max displacement {moved:.2f} A ({shifted} atoms > 2 A)",
        ]

    @staticmethod
    def _sander_command() -> Optional[str]:
        """Path to sander, preferring the AmberTools tree tleap came from."""
        amberhome = MembraneBuilder._amberhome()
        if amberhome:
            candidate = os.path.join(amberhome, "bin", "sander")
            if os.path.isfile(candidate):
                return candidate
        return shutil.which("sander")

    @staticmethod
    def _severe_internal_overlaps(pdb_path: str) -> Optional[int]:
        """Non-bonded atom pairs closer than :data:`SEVERE_OVERLAP_A`, or None.

        Pairs in the same residue or in consecutive residues are excluded: those
        are the bonded and 1-3 contacts, whose separations are genuinely around
        1 A.  Everything else at that separation is an overlap -- two atoms two
        or more residues apart are never legitimately that close -- so this
        needs no connectivity and can run on a PDB.

        ``None`` when it could not be computed (no SciPy), which callers treat
        as "assume there is something to fix" rather than as zero.
        """
        try:
            import numpy as np
            from scipy.spatial import cKDTree
        except ImportError:
            return None
        coordinates, residues = [], []
        index = -1
        previous = None
        try:
            with open(pdb_path) as handle:
                for line in handle:
                    if not line.startswith(("ATOM", "HETATM")):
                        continue
                    key = (line[21], line[22:27])
                    if key != previous:
                        previous = key
                        index += 1
                    residues.append(index)
                    coordinates.append(
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    )
        except (OSError, ValueError):
            return None
        if len(coordinates) < 2:
            return 0
        points = np.asarray(coordinates, dtype=float)
        owner = np.asarray(residues, dtype=np.int64)
        pairs = cKDTree(points).query_pairs(r=SEVERE_OVERLAP_A, output_type="ndarray")
        if pairs.size == 0:
            return 0
        return int(np.count_nonzero(np.abs(owner[pairs[:, 0]] - owner[pairs[:, 1]]) > 1))

    @classmethod
    def _apply_restart_coordinates(
        cls, source_pdb: str, restart: str, reference: str, target_pdb: str,
        remarks: Sequence[str] = (),
    ):
        """Write *source_pdb* with *restart*'s coordinates. Returns (max, n>2A).

        *remarks* are written as REMARK records ahead of the structure. The
        reader in ``patch.py`` keeps only ATOM/HETATM lines, so they travel no
        further than this file -- which is the point: they are for the person
        who opens it.

        The restart and the PDB come from the same tleap unit, so their atom
        order is the same by construction, and substituting coordinates is
        exact in one pass where round-tripping through ParmEd would cost
        another parse of the topology.

        "By construction" is verified rather than assumed, but NOT by bounding
        how far an atom moved -- large motions here are legitimate and expected.
        The atoms that move furthest are precisely the ones tleap invented:
        measured on 9OMO, 14 atoms move more than 5 A and the largest is a
        backbone carbonyl at a crystallographic gap that tleap had placed by
        ideal geometry inside its neighbour, while the mean over all 24606
        atoms is 0.073 A.  A bound tight enough to catch a scrambled file would
        refuse the repair.

        What is checked instead is the correspondence itself: *reference* is the
        coordinate file tleap wrote alongside the PDB, so if those two agree
        atom for atom then the restart -- produced by sander from the topology
        tleap wrote in the same breath -- is in that order too.
        """
        coordinates = cls._read_amber_restart(restart)
        with open(source_pdb) as handle:
            lines = handle.readlines()
        atoms = [i for i, line in enumerate(lines) if line.startswith(("ATOM", "HETATM"))]
        if len(atoms) != len(coordinates):
            raise MembraneBuildError(
                f"{os.path.basename(restart)} holds {len(coordinates)} atoms but "
                f"{os.path.basename(source_pdb)} holds {len(atoms)}; the minimised coordinates "
                "cannot be matched to the structure."
            )

        original = [
            (float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54]))
            for i in atoms
        ]
        cls._assert_same_ordering(reference, original, source_pdb)

        largest, moved = 0.0, 0
        for (x, y, z), position, was in zip(coordinates, atoms, original):
            shift = ((x - was[0]) ** 2 + (y - was[1]) ** 2 + (z - was[2]) ** 2) ** 0.5
            largest = max(largest, shift)
            moved += shift > 2.0
            line = lines[position]
            lines[position] = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        if largest > SOLUTE_MIN_MAX_DISPLACEMENT_A:
            raise MembraneBuildError(
                f"The restrained minimisation moved an atom {largest:.1f} A, beyond the "
                f"{SOLUTE_MIN_MAX_DISPLACEMENT_A:.0f} A backstop. The solute would no longer be "
                "the structure the orientation step framed against the bilayer."
            )
        with open(target_pdb, "w") as handle:
            handle.writelines(f"REMARK   1 {text}\n" for text in remarks)
            handle.writelines(lines)
        return largest, moved

    @staticmethod
    def _count_pdb_atoms(pdb_path: str) -> Optional[int]:
        """ATOM/HETATM records in *pdb_path*, or None if it cannot be read.

        Only ever used to make a header message concrete, so an unreadable file
        drops the number rather than failing the build.
        """
        try:
            with open(pdb_path) as handle:
                return sum(
                    1 for line in handle if line.startswith(("ATOM", "HETATM"))
                )
        except OSError:
            return None

    @staticmethod
    def _solute_order_remarks(n_atoms: Optional[int] = None) -> List[str]:
        """The header that stops this file being used as a simulation input.

        Renaming the file says which order it is in; this says what that costs
        anyone who ignores the name.  Both are needed, because the failure it
        prevents is silent at every stage that could catch it -- grompp accepts
        the pair, and the first symptom is mdrun aborting in ``mshift`` over
        bonds several nanometres long.
        """
        counted = f"{n_atoms} atoms, " if n_atoms else ""
        return [
            "PRISM: relaxed solute, in tleap's savepdb order.",
            "",
            "THIS FILE DOES NOT MATCH topol.top. Do not pair them.",
            "",
            f"It holds the completed solute ({counted}rebuilt side chains relaxed by a",
            "short restrained minimisation) in the order tleap wrote the unit. The",
            "deliverable system.gro is in [ molecules ] order instead -- whole",
            "molecules, i.e. connected components of the bond graph -- because a .gro",
            "carries no molecule information and GROMACS pairs it with the topology",
            "by position alone. A disulfide-linked homodimer reorders heavily under",
            "that rule: the first molecule is one segment of each protomer, taken",
            "from opposite ends of this file.",
            "",
            "Pairing this file with topol.top produces bonds several nm long and",
            "mdrun aborts in mshift. For anything that needs the topology, take the",
            "coordinates from system.gro (or system_dry.gro), which is written in",
            "the order the topology expands to and is checked against it.",
            "",
            "Use this file to inspect what the bilayer was built around, and to diff",
            "against solute_complete.pdb to see which side chains tleap rebuilt.",
        ]

    #: How far the PDB and the coordinate file tleap wrote from the same unit
    #: may differ before they are not the same unit.  The PDB carries three
    #: decimals of an Angstrom and the restart seven, so a correct pair differs
    #: by half a thousandth (measured: 0.0008 A); a mis-ordered pair differs by
    #: Angstroms.
    _COORDINATE_AGREEMENT_A = 0.01

    @classmethod
    def _assert_same_ordering(cls, reference: str, coordinates, source_pdb: str) -> None:
        """Check that *reference* lists the same atoms in the same order."""
        try:
            expected = cls._read_amber_restart(reference)
        except (MembraneBuildError, OSError):
            return  # nothing to compare against; the count check above stands
        if len(expected) != len(coordinates):
            raise MembraneBuildError(
                f"{os.path.basename(reference)} holds {len(expected)} atoms but "
                f"{os.path.basename(source_pdb)} holds {len(coordinates)}."
            )
        for index, (want, got) in enumerate(zip(expected, coordinates), start=1):
            gap = max(abs(a - b) for a, b in zip(want, got))
            if gap > cls._COORDINATE_AGREEMENT_A:
                raise MembraneBuildError(
                    f"{os.path.basename(reference)} and {os.path.basename(source_pdb)} disagree at "
                    f"atom {index} by {gap:.3f} A, so they do not list the solute's atoms in the "
                    "same order and the minimised coordinates cannot be mapped back."
                )

    @staticmethod
    def _read_amber_restart(path: str) -> List[Tuple[float, float, float]]:
        """Coordinates from an AMBER ASCII restart (``ntxo=1``).

        Read as fixed 12-character fields rather than by splitting on
        whitespace: the format is ``6F12.7`` and adjacent values run together
        once a coordinate needs all twelve columns, at which point splitting
        silently returns fewer numbers than there are atoms.
        """
        with open(path) as handle:
            lines = handle.read().splitlines()
        if len(lines) < 2:
            raise MembraneBuildError(f"{os.path.basename(path)} is not an AMBER restart file.")
        try:
            count = int(lines[1].split()[0])
        except (IndexError, ValueError) as exc:
            raise MembraneBuildError(
                f"{os.path.basename(path)} has no atom count on line 2."
            ) from exc
        values: List[float] = []
        needed = 3 * count
        for line in lines[2:]:
            for start in range(0, len(line.rstrip()), 12):
                field = line[start:start + 12].strip()
                if not field:
                    continue
                try:
                    values.append(float(field))
                except ValueError:
                    continue
                if len(values) == needed:
                    break
            if len(values) == needed:
                break
        if len(values) != needed:
            raise MembraneBuildError(
                f"{os.path.basename(path)} holds {len(values) // 3} coordinates for {count} atoms."
            )
        return [tuple(values[i:i + 3]) for i in range(0, needed, 3)]

    #: Residue names whose protonation state is carried by WHICH hydrogens are
    #: present rather than by the residue name.  Everything else keeps its state
    #: in the name (ASH/GLH/LYN/CYM), which survives a hydrogen strip untouched.
    _HISTIDINE_NAMES = frozenset({"HIS", "HIE", "HID", "HIP", "HSD", "HSE", "HSP"})

    #: Which AMBER histidine a set of ring hydrogens means.  ND1-H only is the
    #: delta tautomer, NE2-H only the epsilon tautomer, both the cation.
    _HISTIDINE_STATES = {
        (True, True): "HIP",
        (True, False): "HID",
        (False, True): "HIE",
    }

    @classmethod
    def _is_hydrogen_record(cls, line: str) -> bool:
        """Whether a PDB ATOM/HETATM line describes a hydrogen.

        The element column decides whenever it is filled, which is what keeps
        an ``HG`` mercury ion from being read as a hydrogen -- the same
        misdetection this package already had to fix for CA and CD.  The
        name-based fallback is applied only to ATOM records, because the heavy
        atoms whose names begin with H are metals and metals are HETATM.
        """
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element:
            return element in ("H", "D")
        if not line.startswith("ATOM"):
            return False
        name = line[12:16].strip().upper()
        if not name:
            return False
        return name[0] == "H" or (name[0].isdigit() and name[1:2] == "H")

    def _write_tleap_ready_solute(self, source_pdb: str, target_pdb: str):
        """Write *source_pdb* without hydrogens. Returns (stripped, altlocs, note).

        tleap builds every hydrogen from its own templates, so an input's
        hydrogens are discarded in all but name -- and the names are the
        problem.  A deposited structure names them by its own convention, and
        tleap answers an unrecognised one with ``FATAL: Atom .R<HIE 58>.A<HD1
        18> does not have a type`` only after reading the entire system, so the
        build spends its whole budget to produce no topology at all.  Measured
        on 9OMO: 13 such atoms among 10120 hydrogens, and no prmtop.  Stripping
        them is what pdb4amber does for the same reason and costs nothing.

        What it would cost is the histidine protonation states.  A PDB that
        carries hydrogens states them ONLY through which hydrogens are present:
        strip blindly and every delta-protonated histidine silently becomes
        tleap's default HIE, moving a proton across the ring in every histidine
        in the system.  So the state is read off the hydrogens first and
        written into the residue name, which is where AMBER keeps it and which
        survives the strip.  Residues that already state their state in the
        name are left alone.
        """
        with open(source_pdb) as handle:
            lines = handle.readlines()

        # Pass 1: which ring hydrogens each histidine actually carries.
        ring_hydrogens: Dict[tuple, List[bool]] = {}
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[17:20].strip().upper() not in self._HISTIDINE_NAMES:
                continue
            name = line[12:16].strip().upper()
            if name not in ("HD1", "HE2"):
                continue
            state = ring_hydrogens.setdefault((line[21], line[22:27]), [False, False])
            state[0 if name == "HD1" else 1] = True

        renames = {
            key: self._HISTIDINE_STATES[tuple(state)]
            for key, state in ring_hydrogens.items()
            if tuple(state) in self._HISTIDINE_STATES
        }

        # Pass 2: rewrite.  Only coordinates and chain breaks matter to tleap;
        # CRYST1/ANISOU/CONECT and friends are dropped rather than passed
        # through, because a CONECT into an atom this pass just removed is a
        # record that describes a file that no longer exists.
        #
        # Residues are renumbered 1..N here.  That is not tidying: tleap keeps
        # the first residue's number and then increments, so on a structure
        # whose chains are each numbered 41..850 it answers to 41..1604 while
        # the file says 41..850 twice -- numbers that are neither unique nor
        # tleap's.  Renumbering makes tleap's rule the identity, which is the
        # only way the disulfide bond lines can be built from the file at all.
        out: List[str] = []
        stripped = altlocs = 0
        renamed: Dict[str, int] = {}
        labels: Dict[int, str] = {}
        number = 0
        previous = None
        for line in lines:
            if line.startswith("TER"):
                if out and not out[-1].startswith("TER"):
                    out.append("TER\n")
                continue
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < PDB_MIN_COORDINATE_WIDTH:
                # A truncated coordinate record cannot be renumbered (the chain
                # and residue columns are past its end) and cannot be placed.
                # Reporting it names the offending line; slicing it would raise
                # a bare IndexError from the middle of the build instead.
                raise MembraneBuildError(
                    f"Truncated {line[:6].strip()} record in the solute PDB: "
                    f"{len(line.rstrip())} characters, at least {PDB_MIN_COORDINATE_WIDTH} "
                    f"needed to reach the z column. Offending line: {line.rstrip()!r}. "
                    "Repair or remove it before building."
                )
            altloc = line[16:17].strip()
            if altloc and altloc.upper() != "A":
                altlocs += 1
                continue
            if self._is_hydrogen_record(line):
                stripped += 1
                continue
            key = (line[21], line[22:27])
            if key != previous:
                previous = key
                number += 1
                if number > 9999:
                    raise MembraneBuildError(
                        f"{os.path.basename(source_pdb)} holds more than 9999 residues, which "
                        "cannot be renumbered inside a PDB residue field. Build this solute with "
                        f"--membrane-{BILAYER_CONFIG_FIELD.replace('_', '-')} {BILAYER_PACKMOL}."
                    )
                # Kept so the build report can still name the disulfides the way
                # the depositor numbered them, which is how the user knows them.
                labels[number] = f"{line[21].strip() or 'A'}{line[22:27].strip()}"
            target = renames.get(key)
            if target and line[17:20].strip().upper() != target:
                renamed[target] = renamed.get(target, 0) + 1
                line = line[:17] + f"{target:>3s}" + line[20:]
            # A retained alternate conformer must lose its altLoc flag, or the
            # PDB still claims a partially occupied atom nothing else mentions.
            # The insertion code goes with the renumbering that replaces it.
            out.append(line[:16] + " " + line[17:22] + "%4d " % number + line[27:])
        if not out:
            raise MembraneBuildError(
                f"{os.path.basename(source_pdb)} holds no non-hydrogen ATOM/HETATM records, so "
                "there is no solute to build a bilayer around."
            )
        out.append("END\n")
        with open(target_pdb, "w") as handle:
            handle.writelines(out)

        note = ", ".join(f"{count} {name}" for name, count in sorted(renamed.items()))
        return stripped, altlocs, note, labels

    @staticmethod
    def _tleap_residue_ordinals(pdb_path: str) -> Dict[tuple, Optional[tuple]]:
        """Map ``(chain, resSeq)`` to ``(tleap residue number, residue name)``.

        :func:`packmol_memgen.disulfide_leaplines` emits ``bond SYS.<n>.SG``
        against the number tleap gives a residue, and tleap does not give it the
        number the PDB does.  Measured against AmberTools' teLeap: it keeps the
        FIRST residue's sequence number and then increments by one per residue,
        ignoring both the file's own gaps and the restart at a new chain.  On
        9OMO, whose two chains are each numbered 41..850, that makes chain A
        41..822 and chain B 823..1604 -- so the PDB's own numbers are neither
        unique nor what tleap answers to, and a ``bond`` built from them either
        fails outright or lands on whatever residue holds that number.

        :meth:`_write_tleap_ready_solute` therefore renumbers the solute 1..N
        before tleap ever sees it, which collapses tleap's rule to the identity
        and makes this mapping positional.  It is still computed from the file
        rather than assumed, because that is what makes the residue-name check
        in :meth:`_disulfide_bond_lines` able to catch a violated assumption.

        A key that appears more than once -- an insertion code, a chain whose
        numbering restarts, atoms of one residue split across the file -- maps
        to ``None`` so the caller refuses rather than guesses.

        The rule itself lives in :func:`packmol_memgen.tleap_residue_numbers`,
        which both bilayer routes call.  It is stated there in its general form
        (``first_resSeq + ordinal - 1``); this route's renumbering is the special
        case that collapses it to the identity, and keeping a second copy here
        would only mean the two agreed until one of them was edited.
        """
        return packmol_memgen.tleap_residue_numbers(pdb_path)

    def _disulfide_bond_lines(self, pdb_path: str, rename: bool = False, labels=None):
        """``bond`` lines for *pdb_path*'s disulfides, plus human descriptions.

        tleap infers no disulfide at all.  Left undeclared, the two cysteines
        are built as free cysteines carrying the HG thiol hydrogen that a
        disulfide does not have, and those hydrogens sit on top of each other:
        measured on 9OMO, ``CYS470:HG-CYS488:HG`` at 0.506 A, and 105 of the
        112 hard overlaps in the finished system were protein-internal pairs of
        exactly this kind.  A class-C GPCR's Venus-flytrap domain is full of
        them, so this is the dominant term rather than an edge case.

        Declaring the bond is only half of it: AMBER's bonded cystine is a
        different residue (CYX, with no HG), so the residues are renamed as
        well when ``rename`` is set.  A ``bond`` without the rename leaves a
        three-coordinate sulfur AND both thiol hydrogens in place, which is the
        bug rather than the fix.

        Pairs configured explicitly are honoured alongside the detected ones;
        :func:`packmol_memgen.disulfide_leaplines` de-duplicates and validates.
        """
        try:
            from ..ptm.disulfides import apply_disulfide_renaming, detect_disulfides
        except ImportError as exc:  # pragma: no cover - only on a broken install
            raise MembraneBuildError(
                f"prism.ptm.disulfides is not importable ({exc}), so the disulfides of the solute "
                "cannot be declared to tleap. Building without them would silently produce free "
                "cysteines with overlapping thiol hydrogens."
            ) from exc

        bonds = detect_disulfides(pdb_path)
        ordinals = self._tleap_residue_ordinals(pdb_path)
        pairs, described, unresolved = [], [], []
        for bond in bonds:
            numbers = []
            for sulfur in (bond.a, bond.b):
                entry = ordinals.get((sulfur.chain, sulfur.resid))
                # The residue found at that ordinal must still be a cysteine;
                # if it is not, the mapping is wrong and the bond would land on
                # whatever residue happens to hold that number.
                if entry is None or not str(entry[1]).upper().startswith("CY"):
                    break
                numbers.append(entry[0])
            if len(numbers) != 2:
                unresolved.append(bond.describe())
                continue
            pairs.append(tuple(numbers))
            # Named the way the depositor numbered them when that is still
            # known: "tleap residue 202" is not how anyone refers to a disulfide.
            where = (
                f"{labels.get(numbers[0], numbers[0])}-{labels.get(numbers[1], numbers[1])}"
                if labels
                else f"{numbers[0]}-{numbers[1]}"
            )
            described.append(f"{where} (SG-SG {bond.distance:.2f} A)")
        if unresolved:
            raise MembraneBuildError(
                "Detected disulfide(s) whose residues could not be matched to tleap's own residue "
                "numbering in "
                f"{os.path.basename(pdb_path)}: {'; '.join(unresolved)}. Declaring them against a "
                "guessed number would bond the wrong residues silently, and not declaring them "
                "leaves free cysteines with overlapping thiol hydrogens. Renumber the structure so "
                "each chain numbers its residues uniquely, or build with "
                f"--membrane-{BILAYER_CONFIG_FIELD.replace('_', '-')} {BILAYER_PACKMOL}."
            )

        if rename and bonds:
            # In place: the file was written by this class one step earlier and
            # is read fully before being rewritten.
            apply_disulfide_renaming(
                pdb_path, pdb_path, bonds, str(self.config.protein_forcefield)
            )

        configured = getattr(self.config, "disulfide_pairs", None)
        return packmol_memgen.disulfide_leaplines(pairs + list(configured or ())), described

    def _parametrize_with_tleap(self, dry_pdb: str, box, system_dir: str, result):
        """Run tleap over the dry patch system. Returns (prmtop, inpcrd, log).

        The script is the same shape PACKMOL-Memgen generates for its own
        ``--parametrize`` run -- lipid leaprc, protein leaprc, water leaprc, then
        ``loadpdb``/``set box``/``saveAmberParm`` -- so both routes hand ParmEd a
        topology built by the same tool from the same force fields.

        ``set SYS box`` is used rather than ``setBox SYS "vdw"`` because the patch
        route already knows the periodic cell it tiled, and letting tleap infer
        one from van der Waals radii is what collapses the box onto the contents.

        The solute reaching this point has already been completed by tleap (see
        :meth:`_prepare_patch_solute`), so this run must add nothing to it --
        which is asserted rather than assumed, because that invariant is the
        only thing making the clash deletion upstream describe the system that
        is actually parametrized here.
        """
        protein_ff, water_ff = self.config.packmol_forcefield_pair()
        prmtop = os.path.join(system_dir, PATCH_PRMTOP)
        inpcrd = os.path.join(system_dir, PATCH_INPCRD)

        lines = [
            f"source leaprc.{str(self.config.lipid_ff).strip().lower()}",
            f"source leaprc.protein.{protein_ff}",
            f"source leaprc.water.{water_ff}",
        ]
        lines.extend(self._tleap_ligand_lines())
        # The unit is named SYS because that is the name
        # ``packmol_memgen.disulfide_leaplines`` emits its bond commands against.
        lines.append(f"SYS = loadpdb {os.path.basename(dry_pdb)}")
        # Disulfides must be stated: tleap infers none, and a CYX whose partner
        # was never bonded is a dangling sulfur carrying a charge that assumes a
        # bond -- silent, and present in most GPCRs and channels.  They are
        # re-derived from THIS file rather than carried over from the solute
        # stage, so the residue numbers are the ones tleap is about to assign to
        # the very file it is about to read, with no assumption about how the
        # patch module ordered the system.
        bonds, described = self._disulfide_bond_lines(dry_pdb)
        lines.extend(bonds)
        lines.append("set SYS box { %.3f %.3f %.3f }" % box)
        lines.append(f"saveAmberParm SYS {os.path.basename(prmtop)} {os.path.basename(inpcrd)}")
        lines.append("quit")

        output = self._run_tleap(
            lines,
            PATCH_LEAP_IN,
            system_dir,
            produces=(prmtop, inpcrd),
            doing="parametrizing the patch-built bilayer",
        )
        self._assert_solute_was_already_complete(output)
        # Only now is the cell a fact about the deliverable rather than a
        # request, so this is where _normalize_box's check is armed.
        self._tleap_box_nm = tuple(value / 10.0 for value in box)
        if described:
            result.info.append(
                f"Disulfides declared to tleap for the bilayer system: {len(described)} "
                f"({'; '.join(described)})."
            )
        return prmtop, inpcrd, (
            f"tleap: wrote {os.path.basename(prmtop)}\n"
            f"disulfide bonds declared: {len(bonds)}\n" + self._tail(output, 15)
        )

    def _run_tleap(self, lines: Sequence[str], script_name: str, work_dir: str, produces, doing: str) -> str:
        """Write and run a tleap script; return its combined output.

        Shared by the two tleap runs this route makes -- completing the solute
        and parametrizing the bilayer -- so both get the same stale-output
        removal, the same "tleap exits 0 on a fatal error" handling and the same
        actionable message when AmberTools is absent.  A second copy of that
        handling is a second chance to accept a topology tleap never wrote.
        """
        for stale in produces:
            # tleap leaves the previous file in place when it fails, and a stale
            # file would be consumed as if this run had produced it.
            if os.path.isfile(stale):
                os.remove(stale)

        script = os.path.join(work_dir, script_name)
        with open(script, "w") as handle:
            handle.write("\n".join(lines) + "\n")

        amberhome = self._amberhome()
        tleap = os.path.join(amberhome, "bin", "tleap") if amberhome else "tleap"
        if not shutil.which(tleap) and not os.path.isfile(tleap):
            raise MembraneBuildError(
                f"tleap was not found ({tleap}), so the patch route cannot parametrize the bilayer "
                "it built. Install AmberTools (conda install -c conda-forge ambertools), or "
                f"rebuild with --membrane-{BILAYER_CONFIG_FIELD.replace('_', '-')} {BILAYER_PACKMOL}."
            )
        # _run_gmx is named for its first caller but is a plain "run this in that
        # directory, or tell me it would not start" helper; reused rather than
        # duplicated so both callers get the same timeout and failure handling.
        completed = self._run_gmx(
            [tleap, "-f", os.path.basename(script)], work_dir, timeout=3600
        )
        if completed is None:
            raise MembraneBuildError(f"tleap could not be started ({tleap}).")
        output = (completed.stdout or "") + (completed.stderr or "")
        missing = [os.path.basename(p) for p in produces if not os.path.isfile(p)]
        if missing:
            raise MembraneBuildError(
                f"tleap did not produce {', '.join(missing)} while {doing}. Its output ends:\n"
                + self._tail(output, 30)
            )
        # tleap exits 0 on a fatal error, so the files existing is the real test
        # above; this only surfaces the errors worth reading alongside them.
        fatals = [line for line in output.splitlines() if "FATAL" in line.upper()]
        if fatals:
            raise MembraneBuildError(
                f"tleap reported a fatal error while {doing}:\n" + "\n".join(fatals[:10])
            )
        return output

    #: tleap's own tally of what ``loadpdb`` had to invent, e.g.
    #: ``Leap added 3327 missing atoms according to residue templates:`` then
    #: ``1130 Heavy`` / ``2197 H / lone pairs``.  Parsed rather than inferred
    #: because it is the only place tleap states it, and the difference between
    #: "completed the solute" and "silently rebuilt it again" is exactly this
    #: number.
    _LEAP_ADDED_HEAVY_RE = re.compile(r"^\s*(\d+)\s+Heavy\s*$", re.MULTILINE)
    _LEAP_ADDED_H_RE = re.compile(r"^\s*(\d+)\s+H\s*/\s*lone pairs\s*$", re.MULTILINE)

    @classmethod
    def _leap_added_atoms(cls, output: str) -> Tuple[int, int]:
        """(heavy, hydrogen) atoms tleap reports adding, summed over its tallies."""
        heavy = sum(int(m) for m in cls._LEAP_ADDED_HEAVY_RE.findall(output or ""))
        hydrogens = sum(int(m) for m in cls._LEAP_ADDED_H_RE.findall(output or ""))
        return heavy, hydrogens

    @classmethod
    def _assert_solute_was_already_complete(cls, output: str) -> None:
        """The bilayer run must not add atoms the clash deletion never saw.

        This is the guard on the ordering fix.  Deletion decides which lipids
        the solute displaces from the solute's coordinates; if tleap then
        invents a side chain, that side chain was placed by ideal template
        geometry into a space the deletion had already declared empty, and the
        result is a hard overlap no minimiser recovers from (measured before the
        fix: 1130 heavy atoms added afterwards, 0.395 A minimum protein-lipid
        separation, 6 pairs under 2.0 A).  Since the solute is pre-completed by
        the same tleap with the same leaprc, the correct number here is zero,
        and anything else means the completion did not take.
        """
        heavy, _ = cls._leap_added_atoms(output)
        if heavy:
            raise MembraneBuildError(
                f"tleap added {heavy} missing heavy atom(s) to the assembled bilayer system, but "
                "the solute was completed before the lipids were chosen precisely so that it "
                "would add none. Those atoms were placed into space the clash deletion had "
                "already treated as empty, so the lipids around them are wrong. This is a PRISM "
                "bug rather than a bad input; rebuild with "
                f"--membrane-{BILAYER_CONFIG_FIELD.replace('_', '-')} {BILAYER_PACKMOL} meanwhile."
            )

    def _tleap_ligand_lines(self) -> List[str]:
        """``loadamberparams``/``loadoff`` lines for an embedded ligand, if any."""
        if not self.config.has_ligand():
            return []
        resolver = getattr(packmol_memgen, "_ligand_parameter_paths", None)
        if resolver is None:
            return []
        lines = []
        for entry in resolver(self.config) or ():
            # The backend hands back "frcmod:lib" pairs, the same strings it
            # passes to packmol-memgen's --ligand_param.
            for item in str(entry).split(packmol_memgen.LIGAND_PARAM_SEPARATOR):
                if not item:
                    continue
                verb = "loadoff" if os.path.splitext(item)[1].lower() in (".lib", ".off") else "loadamberparams"
                lines.append(f"{verb} {item}")
        return lines

    def _advise_fast_path(
        self, requested: str, estimate: Optional[int], result: MembraneBuildResult
    ) -> None:
        """Point at the fast path when PACKMOL is about to place a lot of water.

        Suppressed when the fast path was already asked for and declined: the
        reason has been reported, and repeating the suggestion would be advice
        the user has already followed.
        """
        if requested == SOLVATE_GROMACS:
            return
        refusals = getattr(self.config, "solvate_refusals", None)
        if refusals is not None and refusals():
            return
        threshold = int(getattr(self.config, "solvate_auto_threshold", 30000) or 30000)
        if estimate is None or estimate <= threshold:
            return
        result.info.append(
            f"~{estimate} water molecules projected; PACKMOL places each of them individually "
            "and will take many hours at this size. Consider --membrane-solvate gromacs "
            "(packs protein + lipids only, then gmx solvate/genion)."
        )

    def _estimate_water_count(self, protein_pdb: str) -> Optional[int]:
        """Projected water count for this configuration, or None if unknown."""
        estimator = getattr(packmol_memgen, "estimate_water_count", None)
        if estimator is None:
            return None
        try:
            return int(estimator(protein_pdb, self.config))
        except Exception:
            # A geometry estimate is advisory; a failure in it must not be able
            # to stop a build that would otherwise succeed.
            return None

    def _fast_path_blockers(self) -> List[str]:
        """Reasons the GROMACS solvation route cannot be used here."""
        blockers: List[str] = []
        try:
            self._fast_path_modules()
        except ImportError as exc:
            blockers.append(f"the fast-path solvation modules are unavailable ({exc})")
        if "dry" not in self._callable_parameters(packmol_memgen.run):
            blockers.append(
                "the installed packmol_memgen backend cannot pack without solvent "
                "(no 'dry' option on packmol_memgen.run)"
            )
        if self._gmx_command() is None:
            blockers.append(
                "no GROMACS executable found on PATH (tried "
                + ", ".join(GMX_COMMAND_CANDIDATES)
                + "), and gmx solvate/genion place the water on this route"
            )
        # Composition-level refusals (a phosphorus-free bilayer, an
        # unparametrizable lipid force field) are deliberately NOT repeated
        # here: MembraneConfig.resolved_solvate_mode has already applied them
        # and validation_warnings() has already reported them, and a second
        # implementation of the same rule is one that can disagree.
        return blockers

    @staticmethod
    def _callable_parameters(function) -> Set[str]:
        try:
            return set(inspect.signature(function).parameters)
        except (TypeError, ValueError):  # builtins / C callables
            return set()

    @staticmethod
    def _gmx_command() -> Optional[str]:
        for candidate in GMX_COMMAND_CANDIDATES:
            if shutil.which(candidate):
                return candidate
        return None

    @staticmethod
    def _fast_path_modules():
        """Import the fast path's own modules, or raise ImportError.

        Imported lazily and together so that an installation without them (or
        one where they fail to import) degrades to the PACKMOL route instead of
        making :mod:`prism.membrane.builder` unimportable -- the same rule the
        rest of this module applies to external tools.
        """
        from . import solvate as solvate_module
        from . import solvent_params as solvent_params_module

        return solvate_module, solvent_params_module

    def _probe_solvent_parameters(self, system_dir: str, result: MembraneBuildResult):
        """Harvest water/ion parameters for the requested models, or None.

        A dry prmtop contains no water and no ions, so the converted topology
        has neither the atom types nor the ``[ moleculetype ]`` blocks that
        ``gmx solvate``/``gmx genion`` are about to reference.  They are read
        out of a one-shot tleap probe rather than tabulated, because a table
        would answer an ``opc`` request with TIP3P and grompp would not object.
        """
        try:
            _, solvent_params = self._fast_path_modules()
            blocks = solvent_params.solvent_parameter_blocks(
                water_model=self.config.water_model,
                positive_ion=self.config.positive_ion,
                negative_ion=self.config.negative_ion,
                work_dir=system_dir,
                amberhome=self._amberhome(),
            )
            solvent_params.assert_solvent_blocks(blocks)
        except Exception as exc:
            traceback.print_exc()
            result.warnings.append(
                f"Could not obtain {self.config.water_model}/{self.config.positive_ion}/"
                f"{self.config.negative_ion} parameters from tleap ({exc}); falling back to the "
                f"{SOLVATE_PACKMOL} solvation route, which parametrizes its own solvent."
            )
            return None
        result.info.append(
            f"Solvent parameters probed from tleap: {getattr(blocks, 'water_moleculetype', 'SOL')} "
            f"({getattr(blocks, 'water_sites', '?')}-site), "
            f"{getattr(blocks, 'cation_moleculetype', '?')}/{getattr(blocks, 'anion_moleculetype', '?')}."
        )
        return blocks

    @staticmethod
    def _amberhome() -> Optional[str]:
        """AmberTools root as resolved by the packing backend, or None.

        Reuses that module's resolution rather than reading ``$AMBERHOME``: it
        validates the directory really is an AmberTools tree, which is the
        difference between finding tleap and finding an arbitrary prefix.
        """
        resolver = getattr(packmol_memgen, "_resolve_amberhome", None)
        if resolver is None:
            return None
        try:
            return resolver()
        except Exception:
            return None

    # ------------------------------------------------------------------ #
    # Re-running into an existing directory
    #
    # PACKMOL-Memgen reuses any bilayer it finds in its working directory
    # ("Packed PDB ... found. Skipping PACKMOL") and still exits 0, so a re-run
    # into a dirty directory either parametrizes the previous build's system or
    # fails after burning the whole packing budget.  Worse, the MDPs are
    # regenerated unconditionally from the *new* configuration while the stale
    # topology, coordinates and index survive -- localrun.sh's own preflight
    # then passes and runs the old system under the new protocol.  So a re-run
    # either resumes from a complete build or refuses, exactly as the soluble
    # path does with ``--overwrite``.
    # ------------------------------------------------------------------ #

    #: Artifacts whose presence makes a re-run wrong rather than merely noisy.
    #: ``PROT*.pdb``/``*_EMBED*.pdb`` are here for the same reason as the packed
    #: bilayer: packmol-memgen skips MEMEMBED when its output already exists, so
    #: a leftover would orient this build's protein from the previous one.
    _BLOCKING_ARTIFACT_GLOBS = (
        "bilayer_*.pdb", "*.prmtop", "*.crd", "*.rst7", "*.restrt",
        "PROT*.pdb", "*_EMBED*.pdb",
        SYSTEM_GRO, SYSTEM_TOP, SYSTEM_NDX, RUNSCRIPT_NAME, "posre_*.itp",
        # Fast-path intermediates.  A stale system_dry.gro or topol.top.pre_genion
        # describes a different system than the one this run packs, and the
        # vdwradii.dat left behind by an interrupted solvation would silently
        # re-inflate carbon radii for anything else run in this directory.
        DRY_SYSTEM_GRO, "solv.gro", "solv_ions.gro", "vdwradii.dat",
        # Patch-route intermediates.  A leftover tile or trimmed bilayer
        # describes a different composition than the one this run asks for, and
        # it is exactly as misleading as a stale packed bilayer.
        "patch_*.pdb", "patch_*.gro", "bilayer_patch*.pdb", "bilayer_patch*.gro",
        PATCH_DRY_PDB, PATCH_PRMTOP, PATCH_INPCRD,
        # The solute the patch route hands tleap.  A leftover describes the
        # PREVIOUS run's protein, and this run would tile a bilayer around it.
        PATCH_SOLUTE_INPUT_PDB, PATCH_SOLUTE_COMPLETE_PDB, PATCH_SOLUTE_RELAXED_PDB,
        PATCH_SOLUTE_PRMTOP, PATCH_SOLUTE_INPCRD, PATCH_SOLUTE_MIN_RST,
    )
    #: Everything ``overwrite`` clears: the above plus informational leftovers
    #: and the byproducts of the GROMACS solvation stage.
    _STALE_ARTIFACT_GLOBS = _BLOCKING_ARTIFACT_GLOBS + (
        PACKMOL_LOG, PATCH_LOG, MEMEMBED_INPUT_PDB, PATCH_LEAP_IN,
        PATCH_SOLUTE_LEAP_IN, PATCH_SOLUTE_MIN_IN, PATCH_SOLUTE_MIN_OUT, BUILD_FINGERPRINT,
        "leap.log", "leap_patch.log",
        "ions.tpr", "ions.mdp", "mdout.mdp", "solvent_params.itp",
        "topol.top.pre_genion", "preflight*", "#topol.top.*#",
    )

    @staticmethod
    def _matching_artifacts(system_dir: str, patterns) -> List[str]:
        found = set()
        for pattern in patterns:
            found.update(glob.glob(os.path.join(system_dir, pattern)))
        return sorted(entry for entry in found if os.path.isfile(entry))

    #: Configuration fields that determine WHICH system gets built, and so must
    #: match before a previous build may be handed back instead of a new one.
    #: An allowlist rather than the whole dataclass on purpose: a changed
    #: ``packmol_timeout`` or ``preflight_energy`` describes how a build is run,
    #: not what it produces, and forcing a rebuild for those would make
    #: ``--overwrite`` the only usable mode.  Everything here, by contrast,
    #: changes the deliverable -- including the ones that only reach the MDPs
    #: (temperature, production length, restraint force constants), because a
    #: resumed build returns the MDPs it found rather than regenerating them.
    _FINGERPRINT_FIELDS = (
        "lipids", "ratio", "lipid_ff", "orient", "pdb_id",
        "protein_forcefield", "water_model",
        "salt_concentration", "positive_ion", "negative_ion",
        "water_thickness", "xy_distance", "temperature", "production_ns",
        "minimize", "center_tolerance", "validate_orientation",
        "posres_fc_backbone", "posres_fc_sidechain", "posres_fc_lipid",
        "ligand_frcmod", "ligand_lib", "ligand_resname",
        "bilayer_source", "patch_clash_distance", "patch_seam_clash_distance",
        "patch_xy_margin", "solvate", "core_water_purge", "vdwradii_carbon",
    )

    def _fingerprint(self, protein_pdb: str) -> dict:
        """What this build was asked for, as a comparable JSON-safe dict.

        The solute is identified by content rather than by path: the failure
        this exists to stop is "same command, edited input", and a path compares
        equal across exactly that edit.
        """
        config = {}
        for name in self._FINGERPRINT_FIELDS:
            value = getattr(self.config, name, None)
            # json.dumps is the comparison, so the recorded value has to survive
            # a round trip -- tuples come back as lists and would never match.
            config[name] = json.loads(json.dumps(value, default=str))
        return {"version": 1, "config": config, "solute": self._file_identity(protein_pdb)}

    @staticmethod
    def _file_identity(path: str) -> dict:
        digest = hashlib.sha256()
        try:
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1 << 20), b""):
                    digest.update(chunk)
        except OSError as exc:
            return {"name": os.path.basename(path), "error": str(exc)}
        return {
            "name": os.path.basename(path),
            "bytes": os.path.getsize(path),
            "sha256": digest.hexdigest(),
        }

    def _write_fingerprint(self, system_dir: str, protein_pdb: str, result: MembraneBuildResult) -> None:
        """Record what this build was, so a re-run can tell whether it matches."""
        record = self._fingerprint(protein_pdb)
        record["routes"] = {
            "bilayer": result.bilayer_method,
            "solvate": result.solvate_mode,
        }
        try:
            with open(os.path.join(system_dir, BUILD_FINGERPRINT), "w") as handle:
                json.dump(record, handle, indent=2, sort_keys=True)
        except OSError as exc:
            # Losing the record costs the NEXT run its ability to resume; it does
            # not make this build wrong, so it is reported rather than fatal.
            result.warnings.append(
                f"Could not write {BUILD_FINGERPRINT} ({exc}); a re-run into this directory will "
                "refuse to reuse the system because it cannot verify what was built."
            )

    def _fingerprint_differences(self, system_dir: str, protein_pdb: str) -> Optional[List[str]]:
        """Fields where the previous build differs from this request.

        ``None`` means "cannot tell" -- no record, or an unreadable one -- which
        the caller treats exactly as a mismatch.  An empty list means the two
        requests agree.
        """
        path = os.path.join(system_dir, BUILD_FINGERPRINT)
        try:
            with open(path) as handle:
                previous = json.load(handle)
        except (OSError, ValueError):
            return None
        if not isinstance(previous, dict) or previous.get("version") != 1:
            return None
        current = self._fingerprint(protein_pdb)
        differences = []
        for name in self._FINGERPRINT_FIELDS:
            was = previous.get("config", {}).get(name)
            now = current["config"][name]
            if was != now:
                differences.append(f"{name}: built with {was!r}, requested {now!r}")
        was_solute = previous.get("solute") or {}
        now_solute = current["solute"]
        if was_solute.get("sha256") != now_solute.get("sha256"):
            differences.append(
                f"solute: built from {was_solute.get('name', '?')} "
                f"(sha256 {str(was_solute.get('sha256'))[:12]}...), requested "
                f"{now_solute.get('name')} (sha256 {str(now_solute.get('sha256'))[:12]}...)"
            )
        return differences

    def _resume(
        self, system_dir: str, output_dir: str, protein_pdb: str
    ) -> Optional[MembraneBuildResult]:
        """Return a populated result when *system_dir* already holds THIS build.

        Completeness alone is not enough to reuse a directory.  A membrane
        system is fully described by its configuration, and nothing in a
        finished system says which configuration produced it: re-running with
        ``lipids=POPE, bilayer_source=packmol, solvate=packmol, salt=0.5 M,
        T=350 K`` into a directory holding a POPC/patch/0.15 M/303 K build
        returned that build with ``success=True`` and no warnings at all.  So
        the request is compared against the record the build left behind, and a
        directory that cannot be shown to match is refused rather than reused.
        """
        gro = os.path.join(system_dir, SYSTEM_GRO)
        top = os.path.join(system_dir, SYSTEM_TOP)
        ndx = os.path.join(system_dir, SYSTEM_NDX)
        runscript = os.path.join(system_dir, RUNSCRIPT_NAME)
        if not all(os.path.isfile(entry) for entry in (gro, top, ndx, runscript)):
            return None
        posres = {os.path.basename(p): p for p in self._matching_artifacts(system_dir, ("posre_*.itp",))}
        if not posres:
            return None
        mdp_dir = os.path.join(output_dir, MEMBRANE_MDP_DIRNAME)
        mdps = {
            os.path.splitext(os.path.basename(f))[0]: f
            for f in sorted(glob.glob(os.path.join(mdp_dir, "*.mdp")))
        }
        if not REQUIRED_MDP_STAGES.issubset(mdps):
            return None

        differences = self._fingerprint_differences(system_dir, protein_pdb)
        if differences is None:
            raise MembraneBuildError(
                f"{system_dir} already holds a complete membrane system, but no readable "
                f"{BUILD_FINGERPRINT} saying what it was built from, so it cannot be shown to be "
                "the system now being asked for. Handing it back would answer this request with a "
                "different system's topology, coordinates and MDPs. Pass overwrite=True "
                "(CLI: --overwrite) to rebuild it, or build into a clean directory."
            )
        if differences:
            raise MembraneBuildError(
                f"{system_dir} already holds a membrane system, but it was built from a different "
                "request:\n"
                + "".join(f"  - {entry}\n" for entry in differences)
                + "Reusing it would return the earlier system under this request's name. Pass "
                "overwrite=True (CLI: --overwrite) to rebuild it, or build into a clean directory."
            )

        result = MembraneBuildResult(system_dir=system_dir, gro=gro, top=top, ndx=ndx)
        result.posres = posres
        result.runscript = runscript
        result.mdps = mdps
        oriented = os.path.join(system_dir, ORIENTED_PDB)
        if os.path.isfile(oriented):
            result.oriented_pdb = oriented
        # Either route's log, because either route may have written it.  Looking
        # only for packmol_memgen.log lost a resumed patch build the one
        # artifact that names the tool which actually placed its lipids.
        for name in (PACKMOL_LOG, PATCH_LOG):
            candidate = os.path.join(system_dir, name)
            if os.path.isfile(candidate):
                result.packmol_log = candidate
                break
        # Reported from the record rather than left as None: a resumed build
        # that cannot say how its bilayer was built is not reproducible, and
        # these two fields exist precisely because ``auto`` resolves them
        # behind the user's back.
        routes = self._previous_routes(system_dir)
        result.bilayer_method = routes.get("bilayer")
        result.solvate_mode = routes.get("solvate")
        result.success = True
        result.note = (
            f"Final model file {SYSTEM_GRO} already exists in {system_dir} and matches this "
            "request, skipping membrane build. Pass overwrite=True (CLI: --overwrite) to rebuild it."
        )
        return result

    @staticmethod
    def _previous_routes(system_dir: str) -> dict:
        try:
            with open(os.path.join(system_dir, BUILD_FINGERPRINT)) as handle:
                return json.load(handle).get("routes") or {}
        except (OSError, ValueError, AttributeError):
            return {}

    def _refuse_partial_previous_build(self, system_dir: str) -> None:
        """Refuse to build on top of a half-finished previous attempt."""
        leftovers = self._matching_artifacts(system_dir, self._BLOCKING_ARTIFACT_GLOBS)
        if not leftovers:
            return
        raise MembraneBuildError(
            f"{system_dir} already holds artifacts from an incomplete membrane build "
            f"({', '.join(os.path.basename(entry) for entry in leftovers)}). PACKMOL-Memgen would "
            "reuse the existing bilayer instead of packing a new one, and the MDPs would be "
            "regenerated around a topology this run did not build. Pass overwrite=True "
            "(CLI: --overwrite) to clear them, or build into a clean directory."
        )

    def _clear_previous_build(self, system_dir: str) -> None:
        """Delete a previous build's artifacts so this run really rebuilds."""
        for entry in self._matching_artifacts(system_dir, self._STALE_ARTIFACT_GLOBS):
            try:
                os.remove(entry)
            except OSError as exc:
                raise MembraneBuildError(
                    f"Cannot clear the previous membrane build: {entry} could not be removed ({exc})."
                ) from exc

    # ------------------------------------------------------------------ #
    def _orient(self, protein_pdb, system_dir, result: MembraneBuildResult):
        """Write the membrane-framed copy of the solute into *system_dir*.

        Failures here -- a structure outside the bilayer frame, an unreachable
        OPM entry, an orientation method PRISM does not automate -- can never be
        repaired by continuing, so they propagate instead of degrading into a
        one-line warning on a scaffold the caller may mistake for a build.
        """
        info = orient_protein(
            protein_pdb,
            os.path.join(system_dir, ORIENTED_PDB),
            self.config.orient,
            self.config.pdb_id,
            validate=self.config.validate_orientation,
            tolerance=self.config.center_tolerance,
        )
        result.oriented_pdb = info["oriented_pdb"]
        self._record_orientation_note(result, info.get("note"))

    #: Markers :mod:`prism.membrane.orient` uses when its note reports a
    #: problem.  Its note is never empty, so appending it unconditionally made
    #: every successful build report a WARNING and destroyed the signal value of
    #: ``warnings`` for CI and for the MCP payload.  (orient_protein returning a
    #: severity alongside the text would let this stop matching on strings.)
    _ORIENT_PROBLEM_MARKERS = ("NOT verified", "WARNING")

    def _adopt_memembed_orientation(self, system_dir: str, result: MembraneBuildResult) -> None:
        """Make ``protein_oriented.pdb`` honest after an ``orient=memembed`` run.

        With every other method PRISM orients the structure itself and writes
        the result to :data:`ORIENTED_PDB`.  ``memembed`` is different: the
        orientation happens *inside* packmol-memgen, so what PRISM writes under
        that name is the un-oriented input it handed over, while the oriented
        coordinates are left in ``PROT0.pdb``.  Measured on a class-C GPCR
        homodimer the two differ by ~200 A along z -- a file called "oriented"
        that is nothing of the kind, and the one file a user opens to check
        whether the orientation worked.

        So the oriented coordinates are copied back over it.  When they cannot
        be recovered (packing died before MEMEMBED ran) the input is renamed out
        of the way instead, because leaving a misnamed file is the failure this
        exists to prevent.
        """
        if (self.config.orient or "").lower() != "memembed":
            return
        target = os.path.join(system_dir, ORIENTED_PDB)
        source = os.path.join(system_dir, PACKMOL_ORIENTED_PDB)
        if os.path.isfile(source):
            try:
                shutil.copyfile(source, target)
            except OSError as exc:
                result.warnings.append(
                    f"Could not copy MEMEMBED's oriented structure over {ORIENTED_PDB} ({exc}); "
                    f"that file still holds the un-oriented input. The oriented coordinates are "
                    f"in {PACKMOL_ORIENTED_PDB}."
                )
                return
            result.oriented_pdb = target
            result.info.append(
                f"Orientation: {ORIENTED_PDB} now holds MEMEMBED's oriented coordinates "
                f"(copied from {PACKMOL_ORIENTED_PDB}); during packing it held the input frame, "
                "because MEMEMBED runs inside packmol-memgen rather than in PRISM."
            )
            return
        if not os.path.isfile(target):
            return
        fallback = os.path.join(system_dir, MEMEMBED_INPUT_PDB)
        try:
            os.replace(target, fallback)
        except OSError as exc:
            result.warnings.append(
                f"orient=memembed produced no {PACKMOL_ORIENTED_PDB}, and {ORIENTED_PDB} -- which "
                f"holds the UN-oriented input -- could not be renamed ({exc}). Do not treat it as "
                "an oriented structure."
            )
            return
        result.oriented_pdb = fallback
        result.warnings.append(
            f"orient=memembed produced no {PACKMOL_ORIENTED_PDB}, so no oriented structure was "
            f"recovered; the input was renamed to {MEMEMBED_INPUT_PDB} rather than left under a "
            "name that claims an orientation this build never performed."
        )

    @classmethod
    def _record_orientation_note(cls, result: MembraneBuildResult, note) -> None:
        if not note:
            return
        if any(marker in note for marker in cls._ORIENT_PROBLEM_MARKERS):
            result.warnings.append(note)
        else:
            result.info.append(f"Orientation: {note}")

    @staticmethod
    def _write_packmol_log(system_dir, log, name: str = PACKMOL_LOG) -> Optional[str]:
        log_path = os.path.join(system_dir, name)
        try:
            with open(log_path, "w") as fh:
                fh.write(log or "")
        except OSError:
            return None
        return log_path

    def _emit_mdps(self, output_dir, result: MembraneBuildResult):
        try:
            result.mdps = write_membrane_mdps(
                output_dir,
                temperature=self.config.temperature,
                ff_family=self.config.family(),
                production_ns=self.config.production_ns,
                posres_fc_backbone=self.config.posres_fc_backbone,
                posres_fc_sidechain=self.config.posres_fc_sidechain,
                posres_fc_lipid=self.config.posres_fc_lipid,
            )
        except Exception as exc:
            result.warnings.append(f"Membrane MDP generation failed: {exc}")

    @staticmethod
    def _append(result, field: str, message: str) -> None:
        """Append to ``result.<field>`` when the caller supplied that list.

        ``_amber_to_gromacs`` is also driven with lightweight stand-ins for the
        result object, so the conversion must not depend on the full dataclass
        being present just to record a note.
        """
        target = getattr(result, field, None)
        if isinstance(target, list):
            target.append(message)

    def _amber_to_gromacs(
        self, prmtop, inpcrd, out_top, out_gro, result: MembraneBuildResult, dry: bool = False
    ):
        """Convert an AMBER prmtop/inpcrd to GROMACS .top/.gro via ParmEd.

        Returns ``False`` when the conversion failed -- the pre-existing
        contract, which callers test directly.  On success it returns ``True``
        on the PACKMOL route and, on the fast route, the
        :class:`~prism.membrane.solvate.LeafletPlanes` measured from the same
        ParmEd structure: reading them here means the phosphorus positions are
        taken from the topology's own elements, before any coordinate file has
        had a chance to lose that information.

        The index file is deliberately NOT written here any more; it is written
        from the final coordinates once solvation is done (see ``build``).
        """
        try:
            import parmed as pmd

            before_top = self._file_signature(out_top)
            before_gro = self._file_signature(out_gro)
            parm = pmd.load_file(prmtop, xyz=inpcrd)
            if dry:
                # G1, before anything is written: a pack that was supposed to be
                # dry but is not would have its water counted twice -- once by
                # tleap and once by gmx solvate -- and nothing downstream can
                # detect that.
                self._assert_dry_pack(parm, os.path.dirname(os.path.abspath(out_top)))
            # NOT ``parm.save(out_top, format="gromacs")``: ParmEd does not carry
            # AMBER's per-dihedral 1-4 scaling across faithfully, and the failure
            # is silent -- grompp accepts the wrong factor without a murmur.  The
            # correction lives in solvent_params, which also reads the written
            # file back to prove the parameters landed.
            fudge = save_gromacs_topology(
                parm,
                out_top,
                overwrite=True,
                context="the dry membrane topology" if dry else "the membrane topology",
            )
            self._save_coordinates(parm, out_gro, out_top)
            self._append(result, "info", f"1-4 scaling: {fudge.describe()}")
            after_top = self._file_signature(out_top)
            after_gro = self._file_signature(out_gro)
            if after_top is None or after_gro is None:
                raise RuntimeError("ParmEd returned without writing both GROMACS output files.")
            if after_top == before_top or after_gro == before_gro:
                raise RuntimeError(
                    "ParmEd did not create or update both GROMACS outputs; refusing to reuse stale files."
                )
            if dry:
                return self._leaflet_planes(parm)
            return True
        except MembraneBuildError:
            # A correctness guard, not a tool fault: it must not be downgraded
            # into a warning on a scaffold the caller may mistake for a build.
            raise
        except SolventParameterError as exc:
            # Same reasoning, for the guards that now live in solvent_params.
            # These fire when the 1-4 scaling could not be established or could
            # not be verified in the bytes that were written -- i.e. exactly when
            # the topology on disk may be wrong while looking perfectly ordinary.
            # The broad ``except`` below would turn that into a warning and hand
            # back a directory the caller could mistake for a finished build, so
            # it is re-raised as a build error instead.
            raise MembraneBuildError(
                f"The GROMACS topology for {os.path.basename(prmtop)} could not be written with "
                f"the 1-4 scaling its AMBER parameters specify: {exc}"
            ) from exc
        except Exception as exc:
            # ParmEd raises a wide range of types; keep the catch broad but show
            # the traceback so a bug in PRISM is not mistaken for a tool fault.
            traceback.print_exc()
            result.warnings.append(
                f"ParmEd conversion failed ({exc}). As a fallback run: "
                f"parmed -i <(echo 'gromacs') or use acpype on {os.path.basename(prmtop)}."
            )
            return False

    # ------------------------------------------------------------------ #
    # Coordinate output
    #
    # A GROMACS .gro carries no molecule information: it is a flat list, and
    # GROMACS pairs it with the topology purely by position, expanding
    # ``[ molecules ]`` in order.  So the coordinates have to be written in the
    # order the topology's molecules imply, and that is NOT the order the atoms
    # sit in whenever a molecule's atoms are interleaved with another's.
    #
    # Which is exactly what a disulfide-linked homodimer is.  9OMO's protein
    # comes to tleap as eight segments (two chains, each broken by unresolved
    # loops); the intersubunit disulfide Cys136(A)-Cys136(B) bonds chain A's
    # first segment to chain B's first segment, so ParmEd's first molecule is
    # those two segments -- 2932 atoms drawn from opposite ends of the file --
    # while everything between them belongs to later molecules.
    #
    # ParmEd's own .gro writer knows this and tries to recover the order by
    # scanning, for every atom, the not-yet-written atoms for one with matching
    # type/name/residue name.  That is O(n^2) -- 250 s here -- and on this
    # system it does not even finish: it exhausts its candidates and raises
    # ``RuntimeError: Could not find <Atom N [1466]; In HIE 96>``, atom 1466
    # being precisely where the first molecule jumps to the other protomer.
    # Its ``combine='all'`` branch is linear but writes the structure's own
    # order, which here disagrees with the topology in 12239 of 109982 atoms.
    #
    # So the order is derived rather than guessed: molecules are the connected
    # components of the bond graph, enumerated by first appearance with their
    # atoms in ascending index order, which is the enumeration
    # ``Structure.split()`` performs and therefore the one the topology was
    # written from.  It is linear, and it is then CHECKED against the topology
    # that was actually written.
    # ------------------------------------------------------------------ #

    def _save_coordinates(self, parm, out_gro: str, out_top: Optional[str] = None) -> None:
        """Write *parm*'s coordinates in the order *out_top* expands to."""
        order = self._topology_atom_order(parm)
        if order is None:
            # No connectivity to group molecules by, so there is no better order
            # available than the one the object already has: let it write
            # itself, exactly as this did before the ordering was derived.
            parm.save(out_gro, overwrite=True)
            return
        self._write_gro_from_parm(parm, out_gro, order)
        if out_top:
            self._assert_gro_matches_topology(out_gro, out_top)

    @staticmethod
    def _topology_atom_order(parm) -> Optional[List[int]]:
        """Atom indices grouped into molecules, in ``[ molecules ]`` order.

        Molecules are the connected components of the bond graph.  They are
        emitted in order of the lowest atom index each contains, and within a
        molecule the atoms keep their relative order -- the enumeration
        ``Structure.split()`` produces, and so the one the topology lists.

        ``None`` when the object exposes no bond list at all -- not a structure
        whose molecules can be found -- so the caller can fall back.  An EMPTY
        bond list is a different answer (a box of monatomic ions) and is
        answered normally.
        """
        if getattr(parm, "bonds", None) is None:
            return None
        count = len(parm.atoms)
        parent = list(range(count))

        def find(node: int) -> int:
            while parent[node] != node:
                parent[node] = parent[parent[node]]  # path halving
                node = parent[node]
            return node

        for bond in parm.bonds:
            first, second = find(bond.atom1.idx), find(bond.atom2.idx)
            if first != second:
                # Union by index keeps the representative at the lowest atom,
                # so "order of first appearance" needs no extra bookkeeping.
                parent[max(first, second)] = min(first, second)

        groups: Dict[int, List[int]] = {}
        for index in range(count):
            groups.setdefault(find(index), []).append(index)
        order: List[int] = []
        for index in range(count):
            group = groups.pop(index, None)
            if group is not None:
                order.extend(group)
        return order

    @staticmethod
    def _write_gro_from_parm(parm, out_gro: str, order: Sequence[int]) -> None:
        """Write a .gro holding *parm*'s atoms in *order*.

        Velocities are deliberately not written: a system that has just been
        assembled has none, and the staged protocol generates them at the first
        NVT stage.
        """
        # ``parm.box`` is a numpy array, so it is compared to None rather than
        # tested for truth: an all-zero box is falsy but still an answer.
        raw = getattr(parm, "box", None)
        box = [] if raw is None else [float(value) for value in raw]
        if len(box) < 3:
            raise MembraneBuildError(
                "The converted structure has no periodic box, so its coordinates cannot be written "
                "for GROMACS."
            )
        if len(box) >= 6 and any(abs(float(angle) - 90.0) > 1e-3 for angle in box[3:6]):
            # Writing three lengths for a triclinic cell would silently turn it
            # into a rectangular one, and a membrane whose normal is no longer
            # z is not the system that was requested.
            raise MembraneBuildError(
                f"The converted structure has a triclinic box ({tuple(box[3:6])} degrees); this "
                "path builds rectangular membrane cells only."
            )
        coordinates = parm.coordinates
        atoms = parm.atoms
        out = ["Built by PRISM (membrane)\n", f"{len(order):5d}\n"]
        residue = 0
        previous = None
        for serial, index in enumerate(order, start=1):
            atom = atoms[index]
            if atom.residue.idx != previous:
                previous = atom.residue.idx
                residue += 1
            x, y, z = coordinates[index]
            out.append(
                f"{residue % 100000:5d}{atom.residue.name[:5]:<5}{atom.name[:5]:>5}"
                f"{serial % 100000:5d}{x / 10.0:8.3f}{y / 10.0:8.3f}{z / 10.0:8.3f}\n"
            )
        out.append("".join(f"{float(value) / 10.0:10.5f}" for value in box[:3]) + "\n")
        with open(out_gro, "w") as handle:
            handle.writelines(out)

    @classmethod
    def _assert_gro_matches_topology(cls, gro_path: str, top_path: str) -> None:
        """Refuse coordinates that do not line up with the topology just written.

        Cheap and exact: GROMACS pairs the two files by position, so expanding
        ``[ molecules ]`` gives the atom-and-residue-name sequence the topology
        expects, and it either equals the one in the .gro or the system is
        scrambled.  grompp would catch it later -- that is what caught it while
        this was being written -- but only after solvation has spent two minutes
        on a system that was already wrong, and it reports a count rather than
        the atom where the two diverge.
        """
        expected = cls._topology_atom_sequence(top_path)
        if expected is None:
            return  # unparsable topology layout; downstream grompp still checks
        try:
            with open(gro_path) as handle:
                lines = handle.read().splitlines()
            count = int(lines[1].strip())
            written = [(line[10:15].strip(), line[5:10].strip()) for line in lines[2:2 + count]]
        except (OSError, ValueError, IndexError):
            return
        if len(expected) != len(written):
            raise MembraneBuildError(
                f"{os.path.basename(top_path)} describes {len(expected)} atoms but "
                f"{os.path.basename(gro_path)} holds {len(written)}."
            )
        for position, (want, got) in enumerate(zip(expected, written), start=1):
            if want != got:
                raise MembraneBuildError(
                    "The GROMACS topology and coordinates describe different systems: at atom "
                    f"{position} the topology has {want[1]}:{want[0]} and the coordinates have "
                    f"{got[1]}:{got[0]}. Refusing to hand back a system whose atoms would be "
                    "assigned the wrong parameters."
                )

    @staticmethod
    def _topology_atom_sequence(top_path: str):
        """``[(atom name, residue name)]`` in ``[ molecules ]`` order, or None."""
        try:
            with open(top_path) as handle:
                text = handle.read()
        except OSError:
            return None
        if "#include" in text:
            # An assembled topology whose blocks live elsewhere cannot be
            # expanded from this file alone; leave it to grompp.
            return None
        types: Dict[str, List[Tuple[str, str]]] = {}
        molecules: List[Tuple[str, int]] = []
        section = current = None
        for raw in text.splitlines():
            line = raw.split(";")[0].strip()
            if not line:
                continue
            header = re.match(r"\[\s*(\w+)\s*\]", line)
            if header:
                section = header.group(1)
                if section == "moleculetype":
                    current = None
                continue
            if section == "moleculetype" and current is None:
                current = line.split()[0]
                types[current] = []
            elif section == "atoms" and current:
                fields = line.split()
                if len(fields) >= 5:
                    types[current].append((fields[4], fields[3]))
            elif section == "molecules":
                fields = line.split()
                if len(fields) >= 2 and fields[0] in types:
                    try:
                        molecules.append((fields[0], int(fields[1])))
                    except ValueError:
                        return None
                else:
                    return None
        if not molecules:
            return None
        sequence: List[Tuple[str, str]] = []
        for name, repeats in molecules:
            sequence.extend(types[name] * repeats)
        return sequence

    @staticmethod
    def _file_signature(file_path):
        """Lightweight signature used to distinguish this run from stale output."""
        try:
            stat = os.stat(file_path)
        except OSError:
            return None
        if not os.path.isfile(file_path):
            return None
        return stat.st_mtime_ns, stat.st_size

    # ------------------------------------------------------------------ #
    # Fast path: GROMACS solvation of a dry pack
    # ------------------------------------------------------------------ #

    def _assert_dry_pack(self, parm, work_dir: str) -> None:
        """G1: refuse to continue unless the pack really contains no solvent.

        The dry pack is produced by flags (``--nocounter`` plus a near-zero
        density solvent), and flags can be ignored by a release that spells them
        differently.  If they were, tleap would have parametrized a full water
        box that ``gmx solvate`` is about to fill *again*: the topology would
        claim twice the water it has room for, and no later stage looks.
        """
        counts: Dict[str, int] = {}
        for residue in getattr(parm, "residues", ()):
            name = (residue.name or "").strip().upper()
            if name in self._SOLV_RESNAMES:
                counts[name] = counts.get(name, 0) + 1
        if counts:
            raise MembraneBuildError(
                "The bilayer was packed for GROMACS solvation but the parametrized system "
                "already contains solvent ("
                + ", ".join(f"{name}={count}" for name, count in sorted(counts.items()))
                + "). Adding water on top of that would double-count it. Rebuild with "
                "--membrane-solvate packmol, which packs and parametrizes the solvent in one step."
            )
        # The PACKMOL input is the other half of the same statement, and it is
        # the half that says what was *asked* for rather than what came back.
        # packmol-memgen emits one "structure <species>.pdb" block per species
        # it wants placed, and omits the solvent block entirely when the count
        # is zero -- so a solvent block here means the dry lever did not take.
        packmol_inp = os.path.join(work_dir, f"{packmol_memgen.PACKMOL_LOG_PREFIX}.inp")
        try:
            with open(packmol_inp) as handle:
                requested = handle.read()
        except OSError:
            return  # absent input: nothing to contradict the residue histogram
        for match in re.findall(r"^\s*structure\s+(\S+)\.pdb\s*$", requested, re.MULTILINE):
            species = os.path.basename(match).strip().upper()
            if species in self._SOLV_RESNAMES:
                raise MembraneBuildError(
                    f"{os.path.basename(packmol_inp)} still asks PACKMOL to place {species}, so "
                    "this pack is not dry and gmx solvate would add a second solvent phase. "
                    "Rebuild with --membrane-solvate packmol."
                )

    def _leaflet_planes(self, parm):
        """Locate the two leaflets, or refuse the build.

        Unreachable in principle -- a composition without phosphorus is rejected
        before packing -- so reaching it means the packed bilayer is not the one
        that was requested.  There is no recovery: without the planes there is
        no core-water purge, and without the purge the fast path is unsafe.
        """
        solvate_module, _ = self._fast_path_modules()
        try:
            return solvate_module.leaflet_planes_from_parm(parm, self._membrane_resnames())
        except Exception as exc:
            raise MembraneBuildError(
                f"Could not locate the bilayer leaflets in the packed system ({exc}). The "
                "core-water purge that keeps gmx solvate out of the hydrophobic core depends on "
                "them, so this build cannot continue. Rebuild with --membrane-solvate packmol."
            ) from exc

    def _fast_solvate(
        self,
        dry_gro: str,
        top: str,
        planes,
        blocks,
        system_dir: str,
        result: MembraneBuildResult,
    ) -> str:
        """Turn the dry system into the solvated deliverable. Returns its .gro.

        Stages 4-6 of the fast path: normalise the box, merge the solvent's
        force-field blocks into the topology, then hand both to
        :class:`~prism.membrane.solvate.MembraneSolvator`, which runs
        ``gmx solvate`` with inflated carbon radii, purges whatever water still
        reached the hydrophobic core, and places ions by explicit count.
        """
        if planes is None:
            raise MembraneBuildError(
                "The bilayer leaflet planes were not measured during conversion, so the "
                "core-water purge cannot run. Rebuild with --membrane-solvate packmol."
            )
        solvate_module, solvent_params = self._fast_path_modules()

        # 4) Box.  Done here rather than with editconf because the bilayer must
        #    keep its position relative to the box faces: any recentring that is
        #    not a single rigid translation would shear the leaflets apart.
        self._progress("normalising the periodic box around the packed bilayer")
        normalized_gro, planes = self._normalize_box(dry_gro, planes, result)

        # 5) Solvent force field.
        self._progress("merging the water/ion parameters into the topology")
        try:
            solvent_params.merge_solvent_blocks(top, blocks)
        except Exception as exc:
            traceback.print_exc()
            raise MembraneBuildError(
                f"Could not add the water/ion parameters to {os.path.basename(top)} ({exc}). "
                "Rebuild with --membrane-solvate packmol."
            ) from exc

        # 6) Water and ions.
        self._progress("solvating with gmx solvate + genion (inflated C radii, core-water purge)")
        purge = bool(getattr(self.config, "core_water_purge", True))
        try:
            solvator = solvate_module.MembraneSolvator(
                model_dir=system_dir,
                gmx_command=self._gmx_command(),
                water_model=self.config.water_model,
                positive_ion=self.config.positive_ion,
                negative_ion=self.config.negative_ion,
                salt_concentration=self.config.salt_concentration,
                # Always overwrite: the builder has already refused or cleared
                # any previous build, so anything found here belongs to this run
                # and reusing it would mean reusing a half-finished stage.
                overwrite=True,
            )
            # Lets the solvator re-derive the leaflet planes from the
            # coordinates it was handed and cross-check them against ``planes``
            # -- the only independent test that the box normalisation above
            # moved the leaflets and the solute by the same amount.  Optional in
            # the agreed signature, so only passed when it is accepted.
            extra = (
                {"membrane_resnames": self._membrane_resnames()}
                if "membrane_resnames" in self._callable_parameters(solvator.run)
                else {}
            )
            outcome = solvator.run(
                normalized_gro,
                top,
                planes,
                self._solute_resnames(),
                carbon_radius=float(getattr(self.config, "vdwradii_carbon", 0.375)),
                purge=purge,
                **extra,
            )
        except MembraneBuildError:
            raise
        except Exception as exc:
            traceback.print_exc()
            raise MembraneBuildError(
                f"GROMACS solvation of the membrane system failed ({exc}). Rebuild with "
                "--membrane-solvate packmol to have PACKMOL place the water instead."
            ) from exc

        self._record_solvation(outcome, result)

        # The deliverable's name is part of the contract every consumer (the run
        # script, prism.sim, the MCP validator) relies on, and it is only earned
        # once solvation has succeeded.
        final_gro = os.path.join(system_dir, SYSTEM_GRO)
        self._adopt_artifact(getattr(outcome, "gro", None), final_gro, "solvated coordinates")
        self._adopt_artifact(getattr(outcome, "top", None), top, "solvated topology")
        return final_gro

    @staticmethod
    def _adopt_artifact(produced: Optional[str], expected: str, description: str) -> None:
        """Put *produced* at *expected* (a copy, so the intermediate survives)."""
        if not produced:
            raise MembraneBuildError(f"Solvation returned no {description}.")
        if not os.path.isfile(produced):
            raise MembraneBuildError(f"Solvation reported {description} at {produced}, which does not exist.")
        if os.path.abspath(produced) != os.path.abspath(expected):
            shutil.copyfile(produced, expected)

    def _record_solvation(self, outcome, result: MembraneBuildResult) -> None:
        """Copy the solvation tallies onto the result and report them."""
        result.water_count = getattr(outcome, "n_water", None)
        result.purged_core_waters = getattr(outcome, "purged_core", None)
        result.retained_pore_waters = getattr(outcome, "retained_pore", None)
        result.effective_salt_concentration = getattr(outcome, "effective_conc", None)
        positives = getattr(outcome, "n_pos", None)
        negatives = getattr(outcome, "n_neg", None)
        if positives is not None:
            result.ion_counts[str(self.config.positive_ion)] = positives
        if negatives is not None:
            result.ion_counts[str(self.config.negative_ion)] = negatives
        result.warnings.extend(getattr(outcome, "warnings", ()) or ())

        summary = (
            f"GROMACS solvation: {number(result.water_count)} waters, "
            + ", ".join(f"{name}={number(count)}" for name, count in sorted(result.ion_counts.items()))
            + (
                f", effective salt {result.effective_salt_concentration:.3f} M"
                if isinstance(result.effective_salt_concentration, float)
                else ""
            )
        )
        result.info.append(summary)
        if result.purged_core_waters is not None:
            result.info.append(
                f"Hydrophobic core: {number(result.purged_core_waters)} water molecule(s) deleted, "
                f"{number(result.retained_pore_waters or 0)} retained next to the solute "
                "(pore/interface water, which is real and must not be purged)."
            )
        # A channel or transporter that ends up with no pore water at all is
        # usually a sign the retention radius did not reach into the lumen; it
        # is not wrong enough to fail a build, but it is worth a look.
        if result.retained_pore_waters == 0:
            result.info.append(
                "No water was retained inside the membrane-spanning region; if this solute has a "
                "pore or lumen, inspect it before running."
            )

    # ------------------------------------------------------------------ #
    # Position restraints
    #
    # The staged equilibration protocol releases restraints over several
    # stages via the MDP ``define`` line, so a single topology has to carry
    # every restraint block, switched on by ``-DPOSRES`` / ``-DPOSRES_LIPID``
    # and scaled by the ``POSRES_FC_*`` macros.
    #
    # Restraints are per-``[ moleculetype ]`` and use atom numbers LOCAL to
    # that block, which is why this works off the emitted topology rather than
    # the ParmEd object: ParmEd names the blocks ``system1``, ``system2`` ...
    # in molecule order, and only their residue content says which is the
    # protein and which is a lipid.  Reading the file back keeps this
    # classification consistent with ``_write_membrane_index``.
    # ------------------------------------------------------------------ #

    #: Backbone atoms restrained at the stiffer POSRES_FC_BB constant.
    _BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
    #: Below this mass an atom is hydrogen; restraining hydrogens is pointless
    #: because they follow their heavy atom through the LINCS constraints.
    _HYDROGEN_MASS_CEILING = 2.0
    #: Phosphorus (30.97) identifies a lipid headgroup across every lipid type,
    #: without hardcoding per-lipid atom names.  Sulfur (32.06) sits safely
    #: outside this window.
    _PHOSPHORUS_MASS_RANGE = (30.5, 31.5)
    #: Sterols carry no phosphate.  Cholesterol's single hydroxyl oxygen sits at
    #: the same bilayer depth as a phospholipid's phosphate, and is what
    #: CHARMM-GUI restrains for these residues.  Without this fallback a
    #: POPC:CHL1 bilayer would have its cholesterol completely unrestrained.
    _OXYGEN_MASS_RANGE = (15.5, 16.5)

    def _write_position_restraints(self, top_path, system_dir, result: MembraneBuildResult):
        """Emit position-restraint itp files and include them from the topology."""
        try:
            with open(top_path) as fh:
                lines = fh.readlines()

            # Injecting a second time would ADD a [ position_restraints ] section
            # rather than replace it, multiplying every force constant by the
            # number of builds -- and grompp accepts that silently.
            already = [line for line in lines if "posre_solu.itp" in line or "posre_memb.itp" in line]
            if already:
                raise ValueError(
                    "topology already includes position-restraint files; refusing to inject "
                    "a second copy (that would multiply the restraint force constants)"
                )

            blocks = self._parse_moleculetype_blocks(lines)
            if not blocks:
                raise ValueError("no [ moleculetype ] sections found in the topology")

            membrane_resnames = self._membrane_resnames()
            solute_resnames = self._solute_resnames()
            protein_seen = 0
            lipid_seen = 0
            # Insertions are collected first and applied back-to-front so the
            # earlier blocks' line numbers stay valid.
            insertions = []
            written = {}

            for block in blocks:
                atoms = block["atoms"]
                if not atoms:
                    # Every [ moleculetype ] has an [ atoms ] section, so an
                    # empty one means the column layout defeated the parser.
                    # Continuing would quietly leave that molecule unrestrained.
                    raise ValueError(
                        "a [ moleculetype ] yielded no parsable [ atoms ] rows; "
                        "the topology layout was not understood"
                    )
                resnames = {atom["residue"] for atom in atoms}
                if resnames & solute_resnames:
                    protein_seen += 1
                    name = "posre_solu.itp" if protein_seen == 1 else f"posre_solu{protein_seen}.itp"
                    rows = self._protein_restraint_rows(atoms)
                    guard = "POSRES"
                elif resnames & membrane_resnames:
                    lipid_seen += 1
                    name = "posre_memb.itp" if lipid_seen == 1 else f"posre_memb{lipid_seen}.itp"
                    rows = self._lipid_restraint_rows(atoms)
                    guard = "POSRES_LIPID"
                else:
                    continue  # water and ions are never restrained
                if not rows:
                    # A solute or lipid block that yields no restrainable atom
                    # means the topology was parsed wrongly (e.g. a column
                    # layout we did not anticipate).  Failing here is what keeps
                    # an unrestrained run from being reported as restrained.
                    raise ValueError(
                        f"no restrainable atoms found for a {guard} block; "
                        "the topology's [ atoms ] layout was not understood"
                    )

                itp_path = os.path.join(system_dir, name)
                with open(itp_path, "w") as fh:
                    fh.write(self._restraint_itp_text(guard, rows))
                written[name] = itp_path
                insertions.append((block["end"], name, guard))

            if not protein_seen:
                # A membrane build always has a solute; finding none means the
                # classification failed rather than that none was wanted.
                raise ValueError(
                    "no protein/ligand [ moleculetype ] was recognised in the topology, "
                    "so the solute would be equilibrated without restraints"
                )

            for end_index, name, guard in sorted(insertions, reverse=True):
                lines[end_index:end_index] = [
                    "\n",
                    "; Position restraints (PRISM membrane protocol)\n",
                    f"#ifdef {guard}\n",
                    f'#include "{name}"\n',
                    "#endif\n",
                ]

            with open(top_path, "w") as fh:
                fh.writelines(lines)

            result.posres = written
            if self.verbose and written:
                print_success(f"Position restraints: {', '.join(sorted(written))}")
        except Exception as exc:
            # Without restraints the equilibration MDPs cannot be grompp'd, so
            # this is loud rather than silent -- but it must not destroy an
            # otherwise usable topology/coordinate pair.  The traceback stays
            # because a parser bug and a malformed topology look identical here.
            traceback.print_exc()
            result.warnings.append(
                f"Position-restraint generation failed ({exc}). grompp does NOT reject a "
                "-DPOSRES define that no topology consumes, so the staged equilibration would "
                "run completely unrestrained; this build is marked unsuccessful for that reason."
            )

    @staticmethod
    def _parse_moleculetype_blocks(lines):
        """Split a GROMACS topology into ``[ moleculetype ]`` blocks.

        Returns a list of ``{"end", "atoms"}`` where ``end`` is the line index
        at which an include may be appended and ``atoms`` holds the parsed
        ``[ atoms ]`` rows.
        """
        blocks = []
        current = None
        section = None
        for index, raw in enumerate(lines):
            line = raw.split(";", 1)[0].strip()
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip().lower()
                if section == "moleculetype":
                    if current is not None:
                        current["end"] = index
                        blocks.append(current)
                    current = {"end": len(lines), "atoms": []}
                elif section == "system" and current is not None:
                    current["end"] = index
                    blocks.append(current)
                    current = None
                continue
            if current is None or section != "atoms" or not line:
                continue
            fields = line.split()
            # nr type resnr residue atom cgnr charge mass
            if len(fields) < 8:
                continue
            try:
                number = int(fields[0])
                mass = float(fields[7])
            except ValueError:
                continue
            current["atoms"].append(
                {"nr": number, "residue": fields[3].upper(), "name": fields[4].upper(), "mass": mass}
            )
        if current is not None:
            blocks.append(current)
        return blocks

    def _protein_restraint_rows(self, atoms):
        """Backbone at POSRES_FC_BB, remaining heavy atoms at POSRES_FC_SC."""
        rows = []
        for atom in atoms:
            if atom["mass"] < self._HYDROGEN_MASS_CEILING:
                continue
            macro = "POSRES_FC_BB" if atom["name"] in self._BACKBONE_ATOMS else "POSRES_FC_SC"
            rows.append((atom["nr"], macro, macro, macro))
        return rows

    def _lipid_restraint_rows(self, atoms):
        """Restrain the headgroup phosphorus in z only.

        Restraining x/y as well would freeze lateral diffusion and prevent the
        bilayer from relaxing its area per lipid -- the very thing the
        semiisotropic barostat exists to allow.
        """
        low, high = self._PHOSPHORUS_MASS_RANGE
        anchors = [atom for atom in atoms if low < atom["mass"] < high]
        if not anchors:
            # A sterol (cholesterol) has no phosphate; anchor it by its hydroxyl
            # oxygen instead.  Falling through with an empty list would leave
            # that leaflet component silently free while everything else is
            # held -- the bilayer would deform around it during equilibration.
            oxygen_low, oxygen_high = self._OXYGEN_MASS_RANGE
            anchors = [atom for atom in atoms if oxygen_low < atom["mass"] < oxygen_high]
        return [(atom["nr"], "0", "0", "POSRES_FC_LIPID") for atom in anchors]

    @staticmethod
    def _restraint_itp_text(guard, rows):
        """Render a restraint itp whose force constants come from MDP macros."""
        defaults = {
            "POSRES": (("POSRES_FC_BB", "1000.0"), ("POSRES_FC_SC", "500.0")),
            "POSRES_LIPID": (("POSRES_FC_LIPID", "1000.0"),),
        }[guard]
        out = [
            "; Generated by PRISM (membrane protocol).\n",
            "; Force constants are supplied by the MDP 'define' line so a single\n",
            "; topology serves every stage of the staged restraint release.\n",
        ]
        for macro, fallback in defaults:
            out += [f"#ifndef {macro}\n", f"#define {macro} {fallback}\n", "#endif\n"]
        out += ["\n[ position_restraints ]\n", ";  atom  type            fx            fy            fz\n"]
        for number, fx, fy, fz in rows:
            out.append(f"{number:7d}     1 {fx:>13} {fy:>13} {fz:>13}\n")
        return "".join(out)

    def _solute_resnames(self):
        """Protein (and embedded-ligand) residue names.

        An embedded ligand is grouped with the protein everywhere -- thermostat
        group and position restraints alike -- because it is part of the solute
        being equilibrated, not part of the lipid or solvent bath.
        """
        names = set(self._PROTEIN_RESNAMES)
        if self.config.has_ligand():
            names.add(str(self.config.ligand_resname).strip().upper())
        return names

    def _membrane_resnames(self):
        """Residue names that belong to the bilayer, including Amber aliases."""
        names = {
            str(lipid).strip().upper()
            for lipid in self.config.lipids
            if isinstance(lipid, str) and lipid.strip()
        }
        for lipid in tuple(names):
            names.update(self._AMBER_LIPID_RESIDUES.get(lipid, ()))
        if "CHOL" in names:
            names.add("CHL1")
        if "CHL1" in names:
            names.add("CHOL")
        return names

    # Residue-name classification for the SOLU/MEMB/SOLV thermostat groups.
    #
    # This classifier asks "which thermostat group does this atom belong to?",
    # which is a stricter question than the one orient.py asks: it has to
    # partition *every* atom in the system, and it fails closed on anything it
    # cannot place.  Both sets are therefore composed from the shared vocabulary
    # in config.py rather than maintained here as a second copy -- the copies had
    # drifted, and the drift only ever surfaced as a hard ValueError at index
    # generation, after packing, parametrization and solvation had all run.
    #
    # SOLV is water + bath ions.  It now covers every ion spelling PRISM will
    # accept in ``positive_ion``/``negative_ion``, because that name propagates
    # verbatim into the residue name: ``positive_ion: "Li+"`` produced residues
    # literally named ``LI+``, which this classifier had no entry for and so
    # killed the build at the very last step.
    _SOLV_RESNAMES = membrane_config.SOLVENT_RESIDUE_NAMES
    # SOLU is the protein.  Deliberately narrower than orient's measurement set
    # by exactly ``UNK``: orient should still measure an unidentified polymer
    # residue, but assigning one to a thermostat group is a guess this classifier
    # must not make.  ``MSE`` (selenomethionine) is in here now -- its absence is
    # what made routine crystal structures fail at this step.
    _PROTEIN_RESNAMES = membrane_config.PROTEIN_RESIDUE_NAMES

    # Amber Lipid17/21 represents many phospholipids as several bonded
    # residues rather than retaining the input PACKMOL residue name.  For
    # example, one POPC molecule is written as PA (palmitoyl), PC (headgroup),
    # and OL (oleoyl).  These names come from AmberTools'
    # charmm-lipid-to-Amber mapping.  Restrict the aliases to lipids explicitly
    # requested by the user so an unrelated unknown residue is never silently
    # assigned to the membrane thermostat group.
    #
    # The unordered *membership* form is the only one this file needs; the patch
    # route's ordered/stoichiometric counterpart is a different table for a
    # different question, and config.py documents the exact relationship between
    # the two rather than pretending either one can be derived from the other by
    # accident.
    _AMBER_LIPID_RESIDUES = membrane_config.LIPID_FRAGMENT_RESIDUES

    def _lipid_residue_counts(self, residue_names: Sequence[str]) -> Dict[str, int]:
        """Count the residues attributable to each requested lipid.

        packmol-memgen turns each ratio into ``int(round(fraction * lipnum))``,
        so a minor component can round to zero and still produce a bilayer that
        passes every other check -- a POPC:CHL1 request silently delivered as
        pure POPC.  Fragment names are shared between lipids (POPC and POPE are
        both PA + head group + OL), so only fragments unique to one requested
        lipid can attribute a residue; a species with no unique fragment is left
        out rather than counted wrongly.

        These are residues, not molecules: Lipid21 splits one POPC into PA + PC
        + OL and one DPPC into PA + PA + PC, so a molecule tally would need a
        per-lipid head-group table.  Only the zero case is interpreted.

        Takes plain residue names rather than a ParmEd structure because the
        index is now built from the final coordinates, which are the only
        description of the system that includes the GROMACS-placed solvent.
        """
        if not residue_names:
            return {}
        names_upper = [str(name).strip().upper() for name in residue_names]
        aliases = {}
        for lipid in self._requested_lipids():
            names = {lipid} | set(self._AMBER_LIPID_RESIDUES.get(lipid, ()))
            if lipid in ("CHOL", "CHL1"):
                names |= {"CHOL", "CHL1"}
            aliases[lipid] = names
        counts = {}
        for lipid, names in aliases.items():
            shared = set().union(*(other for key, other in aliases.items() if key != lipid))
            exclusive = names - shared
            if not exclusive:
                continue
            counts[lipid] = sum(1 for name in names_upper if name in exclusive)
        return counts

    def _requested_lipids(self):
        return {
            str(lipid).strip().upper()
            for lipid in self.config.lipids
            if isinstance(lipid, str) and lipid.strip()
        }

    # ------------------------------------------------------------------ #
    # Coordinate-file handling
    #
    # Parsed here rather than through ParmEd/MDAnalysis because these are the
    # only two operations needed (a rigid translation and a residue-name scan),
    # both are exact in fixed-column .gro, and neither is worth a dependency
    # that the fast path would then have to be gated on.
    # ------------------------------------------------------------------ #

    @staticmethod
    def _gro_field_layout(record: str) -> Tuple[int, int]:
        """Return ``(field_width, decimals)`` of a .gro coordinate record.

        GROMACS itself derives the precision from the spacing of the first two
        decimal points rather than assuming ``%8.3f``: a high-precision file
        writes wider fields, and reading it at the wrong width silently yields
        nonsense coordinates instead of an error.
        """
        body = record[20:]
        first = body.find(".")
        second = body.find(".", first + 1)
        if first < 0 or second <= first:
            return 8, 3
        width = second - first
        decimals = width - 5
        if width < 6 or decimals < 1:
            return 8, 3
        return width, decimals

    @classmethod
    def _read_gro(cls, gro_path: str):
        """Parse a .gro into ``(title, prefixes, coords, tails, box, layout)``.

        ``prefixes`` and ``tails`` are kept verbatim so a rewrite touches only
        the coordinates: residue/atom naming and any velocities must survive
        byte-for-byte, or the topology stops matching the coordinates.
        """
        with open(gro_path) as handle:
            lines = handle.read().splitlines()
        if len(lines) < 3:
            raise ValueError(f"{os.path.basename(gro_path)} is too short to be a .gro file")
        try:
            natoms = int(lines[1].strip())
        except ValueError as exc:
            raise ValueError(
                f"{os.path.basename(gro_path)} has no atom count on line 2"
            ) from exc
        if len(lines) < natoms + 3:
            raise ValueError(
                f"{os.path.basename(gro_path)} declares {natoms} atoms but holds "
                f"{max(len(lines) - 3, 0)}"
            )
        records = lines[2:2 + natoms]
        box = [float(value) for value in lines[2 + natoms].split()]
        width, decimals = cls._gro_field_layout(records[0]) if records else (8, 3)
        prefixes, coords, tails = [], [], []
        for record in records:
            body = record[20:]
            try:
                xyz = tuple(float(body[i * width:(i + 1) * width]) for i in range(3))
            except ValueError as exc:
                raise ValueError(f"unparsable .gro coordinate record: {record!r}") from exc
            prefixes.append(record[:20])
            coords.append(xyz)
            tails.append(body[3 * width:])
        return lines[0], prefixes, coords, tails, box, (width, decimals)

    @staticmethod
    def _write_gro(gro_path, title, prefixes, coords, tails, box, layout) -> None:
        width, decimals = layout
        out = [f"{title}\n", f"{len(coords):5d}\n"]
        for prefix, (x, y, z), tail in zip(prefixes, coords, tails):
            out.append(
                f"{prefix}{x:{width}.{decimals}f}{y:{width}.{decimals}f}{z:{width}.{decimals}f}{tail}\n"
            )
        out.append("".join(f"{value:10.5f}" for value in box) + "\n")
        with open(gro_path, "w") as handle:
            handle.writelines(out)

    @staticmethod
    def _gro_residue_names(gro_path: str) -> Tuple[List[str], List[str]]:
        """Return ``(per-atom residue names, per-residue residue names)``."""
        with open(gro_path) as handle:
            lines = handle.read().splitlines()
        if len(lines) < 3:
            raise ValueError(f"{os.path.basename(gro_path)} is too short to be a .gro file")
        natoms = int(lines[1].strip())
        if len(lines) < natoms + 3:
            raise ValueError(
                f"{os.path.basename(gro_path)} declares {natoms} atoms but holds "
                f"{max(len(lines) - 3, 0)}"
            )
        atom_resnames: List[str] = []
        residue_names: List[str] = []
        previous_key = None
        for record in lines[2:2 + natoms]:
            # Columns 0-4 are the residue number and 5-9 the residue name; the
            # pair changes exactly at a residue boundary.  (.gro renumbers
            # residues modulo 100000, so the pair can repeat -- but never in two
            # adjacent residues, which is all this needs.)
            key = record[:10]
            name = record[5:10].strip().upper()
            atom_resnames.append(name)
            if key != previous_key:
                residue_names.append(name)
                previous_key = key
        return atom_resnames, residue_names

    def _normalize_box(self, gro_path: str, planes, result: MembraneBuildResult):
        """Rebuild the periodic box around the dry bilayer. Returns (gro, planes).

        tleap sizes the dry system's box for a system that has no water in it,
        so ``gmx solvate`` would either have nowhere to put the water slab or
        would fill a box of the wrong shape.  This sets the box PRISM's
        configuration asks for -- ``water_thickness`` on each face of the
        bilayer, an image gap of :data:`XY_IMAGE_MARGIN_NM` in the membrane
        plane -- and moves the contents into it.

        The move is a single rigid translation per axis and never a recentring
        (``editconf -c``/``-d``): the two leaflets and the solute must keep
        their relative geometry exactly, and the leaflet planes the core-water
        purge uses are absolute z values that are translated by the same amount.

        Contents are placed so that each face carries the requested water
        thickness, rather than by putting the bilayer midplane at ``box_z/2``:
        an asymmetric solute (a class-C GPCR's Venus-flytrap domain spans
        ~18 nm on one side only) would otherwise be pushed out of the box.
        """
        title, prefixes, coords, tails, box, layout = self._read_gro(gro_path)
        if len(box) < 3:
            raise MembraneBuildError(
                f"{os.path.basename(gro_path)} has no usable box vector; the dry system cannot be "
                "solvated. Rebuild with --membrane-solvate packmol."
            )
        if len(box) > 3 and any(abs(value) > 1e-6 for value in box[3:]):
            # A non-rectangular cell cannot be resized one axis at a time, and a
            # membrane never wants one: the bilayer normal must stay along z.
            raise MembraneBuildError(
                f"{os.path.basename(gro_path)} has a triclinic box, which this path cannot resize "
                "without shearing the bilayer. Rebuild with --membrane-solvate packmol."
            )
        if not coords:
            raise MembraneBuildError(f"{os.path.basename(gro_path)} contains no atoms.")

        lows = [min(c[axis] for c in coords) for axis in range(3)]
        highs = [max(c[axis] for c in coords) for axis in range(3)]
        extents = [high - low for low, high in zip(lows, highs)]
        water = float(self.config.water_thickness) / 10.0  # config is in Angstrom

        # G2a: report a dry box that is not the one the bilayer route asked for.
        # What "asked for" means differs by route, and applying one route's
        # question to the other is how this check came to fire on 100% of patch
        # builds while naming a tool that never ran.
        self._check_dry_box(box, extents, water, result)

        # xy comes from a different place on each route, because "the box the
        # route asked for" means a different thing on each.  See _xy_box.
        periodic_xy = (result.bilayer_method or "") == BILAYER_PATCH
        new_box, deltas_xy, xy_note = self._xy_box(
            box, lows, highs, extents, prefixes, coords, periodic_xy
        )
        new_box.append(extents[2] + 2.0 * water)
        deltas = (deltas_xy[0], deltas_xy[1], water - lows[2])
        moved = [
            (x + deltas[0], y + deltas[1], z + deltas[2]) for x, y, z in coords
        ]
        self._write_gro(gro_path, title, prefixes, moved, tails, new_box, layout)
        self._assert_box_normalized(
            gro_path, lows, highs, deltas, new_box, layout, periodic_xy
        )

        result.info.append(
            f"Box normalised to {new_box[0]:.2f} x {new_box[1]:.2f} x {new_box[2]:.2f} nm "
            f"({2 * water:.2f} nm of water along z, {xy_note})."
        )
        return gro_path, planes.shifted(deltas[2])

    def _xy_box(self, box, lows, highs, extents, prefixes, coords, periodic_xy):
        """The membrane-plane box and the xy translation into it.

        Returns ``([box_x, box_y], (dx, dy), note)``.

        The two routes disagree about what the coordinate extent means, and
        using one route's answer on the other is what put an 18.6 A vacuum
        stripe down the x face of every patch build.

        PACKMOL route -- the extent *is* the extent.  packmol-memgen packs
        molecules whole inside a constraint box and nothing is wrapped, so the
        widest atom marks the real edge of the system and the box has to clear
        it by :data:`XY_IMAGE_MARGIN_NM`.  ``max()`` against the box tleap
        reported is the right guard: it keeps a collapsed ``setBox SYS vdw``
        from cutting into the contents.

        Patch route -- the extent is meaningless and the period is exact.
        ``tile_patch`` lays down an integer number of copies of a
        pre-equilibrated patch and wraps each lipid into the cell *by centroid*,
        so a molecule that straddles a face keeps its far end outside it.  The
        cell is periodic *by construction*; the extent is simply that cell plus
        up to one lipid of legitimate overhang (measured: 20.255 nm of lipid
        across an 18.611 nm period in x, 14.092 across 12.469 in y).  Sizing the
        box to the extent therefore does not add clearance, it inserts vacuum:
        18.6 A in x and 18.4 A in y on mGluR7, a stripe with no lipid in it that
        the core-water purge then also empties of water.  Measured consequence:
        the outermost bins of the core-band density profile hold 38-118 atoms
        against a mean of 1124, and the system minimises at -5325 bar instead of
        +344 bar.  So xy is taken from the ``DrySystem`` box verbatim -- via the
        cell tleap was handed and echoed back -- and ``max()`` is not applied,
        because there is nothing here for it to protect against.

        The translation differs too.  Centring the *extent* is right when the
        extent is real; on the patch route it would drag the solute off the
        centre of the cell it was tiled around.  The raw frame is worse still --
        the mGluR7 solute sits at y = -92.1..4.0 nm, almost entirely outside its
        own box -- so xy cannot simply be left alone either.  What is centred
        instead is the solute, which is exactly ``patch._pbc_shift``'s rule and
        for its reason: the solute is the fragile molecule and must stay whole
        and inside, while the lipids are uniform, periodic and *meant* to hang
        over the faces.  ``gmx solvate`` was checked against this rather than
        assumed: handed atoms outside the box it writes them back unmoved and
        still excludes solvent from their in-box images under the minimum image
        convention, so the overhang costs nothing downstream.
        """
        if not periodic_xy:
            new_box = [
                max(box[0], extents[0] + XY_IMAGE_MARGIN_NM),
                max(box[1], extents[1] + XY_IMAGE_MARGIN_NM),
            ]
            deltas = tuple(
                (new_box[axis] - extents[axis]) / 2.0 - lows[axis] for axis in range(2)
            )
            return new_box, deltas, f"{XY_IMAGE_MARGIN_NM:.2f} nm image gap in xy"

        new_box = [float(box[0]), float(box[1])]
        anchor_lows, anchor_highs, anchor = self._solute_xy_span(
            prefixes, coords, lows, highs
        )
        deltas = tuple(
            new_box[axis] / 2.0 - (anchor_lows[axis] + anchor_highs[axis]) / 2.0
            for axis in range(2)
        )
        return (
            new_box,
            deltas,
            f"xy taken verbatim from the tiled patch period, which is exact by "
            f"construction, with the {anchor} centred in it",
        )

    def _solute_xy_span(self, prefixes, coords, lows, highs):
        """xy bounds of the non-lipid contents. Returns (lows, highs, label).

        Falls back to the whole system when the solute cannot be told apart,
        which keeps a composition this builder does not recognise centred on
        something rather than failing: an off-centre solute is a worse starting
        structure, not a wrong one, whereas the box size decided above is what
        actually has to be right.
        """
        membrane = {name.upper() for name in self._membrane_resnames()}
        solute = [
            xyz
            for prefix, xyz in zip(prefixes, coords)
            if prefix[5:10].strip().upper() not in membrane
        ]
        if not solute:
            return lows, highs, "contents"
        return (
            [min(c[axis] for c in solute) for axis in range(2)],
            [max(c[axis] for c in solute) for axis in range(2)],
            "solute",
        )

    def _check_dry_box(self, box, extents, water: float, result: MembraneBuildResult) -> None:
        """Report a dry box that is not the one this route asked tleap for.

        On the PACKMOL route the question is whether ``--tight_box`` took
        effect: packmol-memgen asks tleap for a cell sized around the water slab
        it packed, and ``setBox SYS vdw`` collapsing that onto the contents is a
        real, recoverable fault worth naming in the log.

        On the patch route that question is meaningless and the check was
        answering it anyway.  packmol-memgen never runs, so ``--tight_box`` is
        not a flag anybody passed; and the route deliberately does NOT size z
        for water -- ``insert_solute`` adds a fixed headroom to the contents and
        this method sizes the box for the water afterwards, on purpose, because
        two places sizing one axis is how they come to disagree.  So a dry z box
        smaller than the water slab was the design, and comparing a fixed 2.0 nm
        headroom against 1.8x a configurable 1.75 nm slab made the warning true
        by construction on every build.

        What can actually go wrong on the patch route is tleap not honouring the
        explicit ``set SYS box`` it was handed, so that is what is checked here
        -- against the cell that was requested, which the route knows exactly.
        """
        if (result.bilayer_method or "") != BILAYER_PATCH:
            if box[2] - extents[2] < DRY_BOX_WATER_FRACTION * water:
                result.warnings.append(
                    f"The dry system's box is only {box[2]:.2f} nm along z for a "
                    f"{extents[2]:.2f} nm solute+bilayer, i.e. tleap collapsed it onto the contents "
                    f"instead of leaving room for {2 * water:.2f} nm of water. PRISM is rebuilding "
                    "the box, so this is not fatal; it does mean packmol-memgen's --tight_box was "
                    "ignored."
                )
            return

        requested = self._tleap_box_nm
        if requested is None:
            return
        drift = [abs(box[axis] - requested[axis]) for axis in range(3)]
        if max(drift) > TLEAP_BOX_DRIFT_TOLERANCE_NM:
            result.warnings.append(
                "The dry system's box is "
                f"{box[0]:.2f} x {box[1]:.2f} x {box[2]:.2f} nm but the patch route handed tleap "
                f"{requested[0]:.2f} x {requested[1]:.2f} x {requested[2]:.2f} nm, so 'set SYS box' "
                "did not take effect and the cell the lipids were tiled into is not the cell that "
                "was parametrized. PRISM is rebuilding the box, so this is not fatal for z; in xy "
                "it means the bilayer's periodicity may not match its own image."
            )

    @classmethod
    def _assert_box_normalized(
        cls, gro_path, lows, highs, deltas, new_box, layout, periodic_xy=False
    ) -> None:
        """G2b/G2c: the written file really is the translated system, in the box.

        Re-read rather than trusted: the translation is arithmetic that cannot
        go wrong, but the fixed-column rewrite around it can (a coordinate that
        no longer fits its field runs into its neighbour and produces a
        plausible-looking file). Everything checked here is guaranteed by
        construction, so a failure is a PRISM bug, not a tool fault -- which is
        exactly why it must not be allowed to reach a user's trajectory.

        ``periodic_xy`` says the membrane plane is periodic by construction
        rather than by clearance -- the patch route, where lipids are wrapped by
        centroid and so overhang the faces on purpose.  The two xy checks below
        are both phrased in terms of the coordinate extent, which on that route
        does not measure the system's edge, and they would reject every correct
        patch build (``20.255 nm of contents`` against an exact ``18.611 nm``
        period reads as a box 1.6 nm too small).  They are replaced there by the
        check that is actually meaningful: that the overhang stays within one
        molecule of the face, so a wrap that has genuinely gone wrong -- a box
        short by a tile -- is still caught.  The rigid-translation check and
        every z check are route-independent and always run.
        """
        _, _, coords, _, box, _ = cls._read_gro(gro_path)
        tolerance = 1.5 * 10 ** -layout[1]  # one rounding step of the written precision
        for axis, name in enumerate("xyz"):
            observed_low = min(c[axis] for c in coords)
            observed_high = max(c[axis] for c in coords)
            if (
                abs(observed_low - (lows[axis] + deltas[axis])) > tolerance
                or abs(observed_high - (highs[axis] + deltas[axis])) > tolerance
            ):
                raise MembraneBuildError(
                    f"Box normalisation did not translate {name} rigidly "
                    f"(expected {lows[axis] + deltas[axis]:.3f}..{highs[axis] + deltas[axis]:.3f}, "
                    f"wrote {observed_low:.3f}..{observed_high:.3f}); refusing to solvate a "
                    "distorted system."
                )
            if periodic_xy and axis < 2:
                overhang = max(-observed_low, observed_high - box[axis], 0.0)
                if overhang > PATCH_MAX_XY_OVERHANG_NM:
                    raise MembraneBuildError(
                        f"After normalisation the contents hang {overhang:.3f} nm past the {name} "
                        f"face of a {box[axis]:.3f} nm box. The tiled patch wraps lipids by "
                        f"centroid, so up to about half a lipid of overhang is expected, but "
                        f"more than {PATCH_MAX_XY_OVERHANG_NM:.2f} nm means the box is not the "
                        "period the patch was tiled to -- solvating it would leave a vacuum "
                        "stripe along that face."
                    )
                continue
            if observed_low < -tolerance or observed_high >= box[axis]:
                raise MembraneBuildError(
                    f"After normalisation the system spans {observed_low:.3f}..{observed_high:.3f} nm "
                    f"in {name} but the box is {box[axis]:.3f} nm; atoms outside the box would be "
                    "wrapped into their own periodic image."
                )
            if axis < 2 and (observed_high - observed_low) + 0.20 > box[axis]:
                # In the membrane plane the bilayer is continuous across the
                # boundary, so too small a box means head groups overlapping
                # their own image at step 0 rather than merely a tight fit.
                raise MembraneBuildError(
                    f"The normalised box leaves under 0.2 nm between the bilayer and its periodic "
                    f"image in {name} ({box[axis]:.3f} nm box, "
                    f"{observed_high - observed_low:.3f} nm of contents)."
                )
        if any(abs(box[axis] - new_box[axis]) > 1e-4 for axis in range(3)):
            raise MembraneBuildError(
                f"The normalised box was written as {box[:3]} instead of {new_box}."
            )

    def _write_membrane_index_from_gro(
        self, gro_path: str, ndx_path: str, result: Optional[MembraneBuildResult] = None
    ) -> None:
        """Write the thermostat index from the system's final coordinates.

        Both routes end here.  The coordinates are the only artifact that
        describes the system grompp will actually be handed -- on the fast path
        the water and ions were placed by GROMACS and appear in no AMBER
        topology -- and SOLU/MEMB/SOLV are absolute atom indices, so an index
        derived from anything else would silently mis-couple the thermostats.
        """
        try:
            atom_resnames, residue_names = self._gro_residue_names(gro_path)
            self._write_membrane_index(atom_resnames, residue_names, ndx_path, result)
        except Exception as exc:
            # Non-fatal by the same reasoning as before -- the topology and
            # coordinates are still good -- but ``_is_complete`` refuses the
            # build, because every membrane MDP names these groups.
            traceback.print_exc()
            if result is not None:
                result.warnings.append(f"Index (SOLU/MEMB/SOLV) generation failed: {exc}")
            return
        if result is not None:
            result.ndx = ndx_path

    def _write_membrane_index(
        self,
        atom_resnames,
        residue_names=None,
        ndx_path: Optional[str] = None,
        result: Optional[MembraneBuildResult] = None,
    ):
        """Write a GROMACS index file with SOLU / MEMB / SOLV groups.

        SOLU = protein (incl. PTM residues); SOLV = water + ions; MEMB =
        explicitly configured lipids. Unknown residues fail closed instead of
        being silently coupled to the membrane thermostat.

        Takes a per-atom residue-name sequence and the per-residue sequence used
        for the lipid tally.  A ParmEd structure is accepted in the first
        position instead (``_write_membrane_index(parm, ndx_path)``) for callers
        that hold one: the classification is identical either way, and splitting
        it across two entry points is how the two would drift.
        """
        if hasattr(atom_resnames, "atoms"):
            parm = atom_resnames
            atom_resnames, residue_names, ndx_path, result = (
                [atom.residue.name for atom in parm.atoms],
                [residue.name for residue in getattr(parm, "residues", ())],
                residue_names,
                ndx_path if ndx_path is not None else result,
            )
        if not ndx_path:
            raise ValueError("_write_membrane_index needs an output path")

        solu, memb, solv = [], [], []
        membrane_resnames = self._membrane_resnames()
        solute_resnames = self._solute_resnames()
        unknown_resnames = set()
        for i, resname in enumerate(atom_resnames, start=1):
            rn = str(resname).strip().upper()
            if rn in solute_resnames:
                solu.append(i)
            elif rn in self._SOLV_RESNAMES:
                solv.append(i)
            elif rn in membrane_resnames:
                memb.append(i)
            else:
                unknown_resnames.add(rn or "<blank>")

        if unknown_resnames:
            raise ValueError(
                "Cannot classify residue name(s) for membrane thermostat groups: "
                + ", ".join(sorted(unknown_resnames))
                + ". Add supported protein/PTM names or list the lipid explicitly."
            )

        def _block(name, idxs):
            lines = [f"[ {name} ]"]
            for k in range(0, len(idxs), 15):
                lines.append(" ".join(str(x) for x in idxs[k:k + 15]))
            return "\n".join(lines) + "\n"

        groups = {"SOLU": solu, "MEMB": memb, "SOLV": solv}
        empty_groups = [name for name, idxs in groups.items() if not idxs]
        if empty_groups:
            raise ValueError(
                "Cannot generate membrane thermostat index: empty group(s) "
                f"{', '.join(empty_groups)}. Verify protein/lipid/solvent residue names."
            )
        # grompp matches tc-grps by name, so extra groups are free -- and
        # without a whole-system group the trjconv recipe localrun.sh prints
        # cannot be answered: with -n, GROMACS offers only this file's groups
        # and none of the three thermostat groups covers the box.
        groups["System"] = sorted(solu + memb + solv)
        # G9: the three thermostat groups must tile the system exactly.  They do
        # by construction above (every atom is classified or the build stops),
        # so a mismatch here means the atom list itself was truncated -- which is
        # what a half-read coordinate file looks like.
        if len(groups["System"]) != len(atom_resnames):
            raise ValueError(
                f"Thermostat groups cover {len(groups['System'])} atoms but the system has "
                f"{len(atom_resnames)}; every atom must belong to exactly one of SOLU/MEMB/SOLV."
            )

        with open(ndx_path, "w") as fh:
            for name, idxs in groups.items():
                fh.write(_block(name, idxs))

        lipid_counts = self._lipid_residue_counts(residue_names or ())
        absent = sorted(name for name, count in lipid_counts.items() if not count)
        if result is not None:
            result.lipid_counts = lipid_counts
            if absent:
                result.warnings.append(
                    "Requested lipid(s) absent from the packed bilayer: "
                    + ", ".join(absent)
                    + ". packmol-memgen rounds each species to a whole number of molecules, so a "
                    "small ratio component can round to zero; raise its ratio or the patch area."
                )

        if self.verbose:
            print_success(
                f"index.ndx groups: SOLU={number(len(solu))} MEMB={number(len(memb))} "
                f"SOLV={number(len(solv))} atoms"
            )
            if lipid_counts:
                print_info(
                    "Bilayer residues per requested lipid: "
                    + ", ".join(f"{name}={number(count)}" for name, count in sorted(lipid_counts.items()))
                )

    # ------------------------------------------------------------------ #
    # End-to-end gate
    # ------------------------------------------------------------------ #

    def _preflight_energy(self, gro, top, ndx, output_dir, result: MembraneBuildResult) -> None:
        """Grompp and evaluate the finished system at step 0.

        This is the only check that looks at the deliverable as GROMACS will:
        a topology and a coordinate file that disagree about which atom is
        which still grompp into a valid-looking .tpr behind a single warning,
        and the resulting run has a potential energy of order 1e11 kJ/mol
        instead of a negative one.  Nothing cheaper distinguishes "imperfectly
        packed" from "catastrophically wrong", and the two are indistinguishable
        by eye in the output files.

        Missing or unrunnable tools downgrade to a warning -- that is an
        environment problem.  A verdict the tools actually produced does not:
        the build is refused.
        """
        if not bool(getattr(self.config, "preflight_energy", True)):
            result.warnings.append(
                "Membrane preflight energy check disabled; the system was not verified to be "
                "grompp-able or to have a sane starting energy."
            )
            return
        gmx = self._gmx_command()
        # The em stage is the cheapest grompp that still exercises the real
        # thermostat groups and restraint defines.  ``output_dir`` is the
        # fallback source so the gate survives an MDP dict that came back short.
        mdp = result.mdps.get("em") or os.path.join(output_dir, MEMBRANE_MDP_DIRNAME, "em.mdp")
        if not gmx or not mdp or not os.path.isfile(mdp):
            result.warnings.append(
                "Skipped the membrane preflight energy check (no GROMACS executable or no em.mdp); "
                "run 'gmx grompp' on the system before trusting it."
            )
            return

        system_dir = os.path.dirname(os.path.abspath(gro))
        if self.verbose:
            print_info("Membrane preflight: grompp + zero-step mdrun on the finished system")
        grompp = [
            gmx, "grompp", "-f", mdp,
            "-c", os.path.basename(gro),
            # -r is required whenever the mdp carries a -DPOSRES define; passing
            # it unconditionally costs nothing and removes the coupling to which
            # stage the mdp came from.
            "-r", os.path.basename(gro),
            "-p", os.path.basename(top),
            "-n", os.path.basename(ndx),
            "-o", "preflight.tpr", "-po", "preflight_mdout.mdp",
            # High on purpose: the point is to READ the warnings, not to let
            # grompp's own threshold decide.  The atom-name mismatch this exists
            # to catch is one warning, which any -maxwarn >= 1 would swallow.
            "-maxwarn", "100",
        ]
        completed = self._run_gmx(grompp, system_dir, timeout=1800)
        if completed is None:
            result.warnings.append(
                "Could not run the membrane preflight check (GROMACS did not start); "
                "run 'gmx grompp' on the system before trusting it."
            )
            return
        output = (completed.stdout or "") + (completed.stderr or "")
        if completed.returncode != 0:
            raise MembraneBuildError(
                "The finished membrane system does not pass 'gmx grompp':\n"
                + self._tail(output)
            )
        self._assert_atom_names_match(output, result)
        self._assert_neutral(output)

        mdrun = [gmx, "mdrun", "-s", "preflight.tpr", "-deffnm", "preflight", "-nsteps", "0"]
        completed = self._run_gmx(mdrun, system_dir, timeout=3600)
        if completed is None or completed.returncode != 0:
            result.warnings.append(
                "The membrane system was grompp'd successfully but the zero-step energy "
                "evaluation did not run, so its starting energy is unverified."
            )
        else:
            self._assert_starting_energy(os.path.join(system_dir, "preflight.log"), gro, result)
        self._clean_preflight_artifacts(system_dir)

    #: Charge (in e) above which the system is genuinely unneutralised.  grompp
    #: emits its note for any residue at all, and summing tens of thousands of
    #: AMBER partial charges leaves float noise of ~1e-5 e on a system that is
    #: exactly neutral -- failing on that would reject correct builds.  Half an
    #: elementary charge cannot be noise: it is a missing ion.
    _NEUTRALITY_TOLERANCE_E = 0.5

    @classmethod
    def _assert_neutral(cls, grompp_output: str) -> None:
        """G6b: the system grompp saw must carry no net charge."""
        charges = [
            abs(float(match))
            for match in re.findall(
                r"non-zero total charge:\s*([-+0-9.eE]+)", grompp_output
            )
        ]
        if charges and max(charges) >= cls._NEUTRALITY_TOLERANCE_E:
            raise MembraneBuildError(
                f"The finished membrane system carries a net charge of {max(charges):.3f} e; "
                "neutralisation failed. Ewald would compensate it with a uniform background, "
                "which is not the system that was asked for."
            )

    @staticmethod
    def _clean_preflight_artifacts(system_dir: str) -> None:
        """Remove the zero-step run's byproducts, keeping only its log.

        ``preflight.gro`` next to ``system.gro`` is a coordinate file that looks
        like a deliverable and is not one; the log is the evidence worth keeping.
        """
        for name in (
            "preflight.tpr", "preflight.trr", "preflight.edr", "preflight.gro",
            "preflight.cpt", "preflight_prev.cpt", "preflight_mdout.mdp",
        ):
            try:
                os.remove(os.path.join(system_dir, name))
            except OSError:
                pass

    @staticmethod
    def _run_gmx(command: List[str], work_dir: str, timeout: int):
        """Run a GROMACS command, or return None if it could not be started."""
        try:
            return subprocess.run(
                command, cwd=work_dir, capture_output=True, text=True, timeout=timeout
            )
        except (OSError, subprocess.SubprocessError):
            traceback.print_exc()
            return None

    @staticmethod
    def _tail(text: str, lines: int = 20) -> str:
        return "\n".join((text or "").strip().splitlines()[-lines:])

    def _assert_atom_names_match(self, grompp_output: str, result: MembraneBuildResult) -> None:
        """G7: the topology and the coordinates must describe the same atoms.

        grompp reports a mismatch as a count in a single warning and then
        carries on, producing a .tpr of exactly the size a correct one would
        have.  The only benign source of mismatches on this path is genion
        renaming an ion (topology ``Na+`` vs coordinate ``Na``), which is
        bounded by the number of ions placed.
        """
        mismatched = sum(
            int(match) for match in re.findall(r"(\d+)\s+non-matching atom names", grompp_output)
        )
        if not mismatched:
            return
        allowed = sum(int(count) for count in result.ion_counts.values())
        if mismatched <= allowed:
            result.info.append(
                f"grompp reported {mismatched} non-matching atom name(s), within the "
                f"{allowed} ion(s) genion placed (it writes the element, the topology the "
                "moleculetype name)."
            )
            return
        raise MembraneBuildError(
            f"grompp reports {mismatched} non-matching atom names between {SYSTEM_TOP} and "
            f"{SYSTEM_GRO}, more than the {allowed} ion(s) that can legitimately differ. The "
            "topology and the coordinates describe different systems; refusing to hand back a "
            "system that would run at an astronomically wrong energy."
        )

    def _assert_starting_energy(self, log_path: str, gro: str, result: MembraneBuildResult) -> None:
        """G8: the step-0 energy must look like a solvated system's, not a clash."""
        try:
            with open(log_path) as handle:
                log_text = handle.read()
        except OSError:
            log_text = ""
        energies = self._parse_gmx_energies(log_text)
        potential = energies.get("Potential")
        if potential is None:
            # The energy-minimizer summary states the same number in prose when
            # the table did not survive (a converged-at-step-0 run prints only
            # the summary), so try it before giving up on the gate.
            summary = re.findall(r"Potential Energy\s*=\s*([-+0-9.eE]+)", log_text)
            potential = float(summary[-1]) if summary else None
        if potential is None:
            result.warnings.append(
                "Could not read the starting energy from the preflight run; the system is "
                "grompp-able but its energy is unverified."
            )
            return

        atom_count = 0
        try:
            with open(gro) as handle:
                handle.readline()
                atom_count = int(handle.readline().strip())
        except (OSError, ValueError):
            atom_count = 0
        lj = energies.get("LJ (SR)")
        lj_per_atom = lj / atom_count if (lj is not None and atom_count) else None

        if potential >= 0 or (lj_per_atom is not None and lj_per_atom > PREFLIGHT_MAX_LJ_PER_ATOM_KJ):
            raise MembraneBuildError(
                f"The finished membrane system has an unusable starting energy "
                f"(potential {potential:.3e} kJ/mol"
                + (f", LJ(SR) {lj_per_atom:.1f} kJ/mol per atom" if lj_per_atom is not None else "")
                + "). A correctly built system is strongly negative with a few kJ/mol of LJ per "
                "atom; this indicates overlapping atoms or a topology that does not match the "
                "coordinates."
            )
        result.info.append(
            f"Preflight energy: potential {potential:.3e} kJ/mol"
            + (f", LJ(SR) {lj_per_atom:.2f} kJ/mol per atom" if lj_per_atom is not None else "")
            + "."
        )

    @staticmethod
    def _parse_gmx_energies(log_text: str) -> Dict[str, float]:
        """Pull the last ``Energies (kJ/mol)`` table out of an mdrun log.

        The table is fixed 15-character columns of names over values, which is
        the only machine-readable form mdrun offers without a second ``gmx
        energy`` call that would need interactive input.
        """
        lines = log_text.splitlines()
        energies: Dict[str, float] = {}
        for index, line in enumerate(lines):
            if not line.strip().startswith("Energies (kJ/mol)"):
                continue
            block: Dict[str, float] = {}
            cursor = index + 1
            while cursor + 1 < len(lines):
                names = [
                    lines[cursor][i:i + 15].strip()
                    for i in range(0, len(lines[cursor].rstrip()), 15)
                ]
                raw_values = [
                    lines[cursor + 1][i:i + 15].strip()
                    for i in range(0, len(lines[cursor + 1].rstrip()), 15)
                ]
                if not names or not raw_values or len(raw_values) < len(names):
                    break
                try:
                    values = [float(value) for value in raw_values[:len(names)]]
                except ValueError:
                    break
                block.update(dict(zip(names, values)))
                cursor += 2
            if block:
                energies = block  # keep the LAST block: that is step 0's result
        return energies

    # ------------------------------------------------------------------ #
    def _begin_stages(self, fast: bool) -> None:
        """Reset the stage counter for a build on the given route."""
        self._stage_index = 0
        self._stage_total = self._FAST_BUILD_STAGES if fast else self._BUILD_STAGES

    def _progress(self, description: str) -> None:
        """Announce the next build stage.

        Deliberately not ``print_step``: the CLI workflow is already printing
        its own "Step n/N" line around this call, and nesting a second numbering
        under the same word reads as a contradiction.

        The counter advances rather than being passed in, because the two
        solvation routes run different stages and hand-numbering them is how the
        printed sequence ends up skipping or repeating a number.
        """
        self._stage_index = getattr(self, "_stage_index", 0) + 1
        if self.verbose:
            total = getattr(self, "_stage_total", self._BUILD_STAGES)
            print_info(f"Membrane stage {self._stage_index}/{total}: {description}")

    def _report(self, result: MembraneBuildResult) -> None:
        """Print the outcome of a build.

        prism/builder/workflow.py prints the artifact list again on the CLI
        path, but it never prints ``warnings``, so staying silent here would
        drop every diagnostic the build produced.  Warnings come first because
        the note asks the user to act on them.
        """
        if not self.verbose:
            return
        for message in result.warnings:
            print_warning(message)
        for message in result.info:
            print_info(message)
        if result.note:
            (print_success if result.success else print_warning)(result.note)
        if result.oriented_pdb:
            print(f"  oriented:   {path(result.oriented_pdb)}")
        if result.bilayer_method:
            # Which tool placed the lipids, for the same reason solvation is
            # reported: the automatic fallbacks change it silently, and it is not
            # visible anywhere in the output files.
            print(f"  bilayer:    {result.bilayer_method}")
        if result.solvate_mode:
            # Which tool placed the water is the single most consequential thing
            # about a membrane build that is not visible in the output files.
            detail = ""
            if result.water_count is not None:
                detail = f"  ({number(result.water_count)} waters"
                if result.purged_core_waters is not None:
                    detail += f", {number(result.purged_core_waters)} purged from the core"
                detail += ")"
            print(f"  solvation:  {result.solvate_mode}{detail}")
        if result.gro:
            print(f"  structure:  {path(result.gro)}")
        if result.top:
            print(f"  topology:   {path(result.top)}")
        if result.ndx:
            print(f"  index:      {path(result.ndx)}  (use: gmx grompp ... -n {SYSTEM_NDX})")
        if result.posres:
            print(f"  restraints: {', '.join(sorted(result.posres))}")
        if result.packmol_log:
            print(f"  pack log:   {path(result.packmol_log)}")
        if result.runscript:
            print(
                f"  run:        {path(result.runscript)}  "
                f"(cd {os.path.dirname(result.runscript)} && bash {RUNSCRIPT_NAME})"
            )
        if result.mdps:
            mdp_dir = os.path.dirname(next(iter(result.mdps.values())))
            print(f"  mdps:       {', '.join(sorted(result.mdps))} in {path(mdp_dir)}")
