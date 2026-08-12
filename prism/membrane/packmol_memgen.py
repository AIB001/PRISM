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
"""

from __future__ import annotations

import os
import shutil
import struct
import subprocess
from dataclasses import dataclass
from typing import List, Optional

from .config import MembraneConfig, PACKMOL_PARAMETRIZED_LIPID_FFS


@dataclass
class PackmolMemgenResult:
    prmtop: Optional[str]
    inpcrd: Optional[str]
    packed_pdb: Optional[str]
    log: str
    success: bool


def is_available() -> bool:
    """True if packmol-memgen is on PATH."""
    return shutil.which("packmol-memgen") is not None


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


def build_command(
    protein_pdb: str,
    config: MembraneConfig,
    parametrize: bool = True,
    protein_ff_option: str = "--fprot",
    water_ff_option: str = "--ffwat",
) -> List[str]:
    """Construct the packmol-memgen command line for the given config.

    The lipid argument uses PACKMOL-Memgen's ``A:B`` syntax for mixtures and the
    parallel ``--ratio`` flag for molar proportions.
    """
    config.raise_for_errors()
    lipids = ":".join(config.lipids)
    ratio = ":".join(str(r) for r in config.resolved_ratio())
    orient = (config.orient or "").lower()
    if orient not in {"preoriented", "opm", "ppm", "memembed"}:
        raise ValueError(f"Unsupported membrane orientation method: {config.orient!r}")

    cmd = [
        "packmol-memgen",
        "--pdb", protein_pdb,
        "--lipids", lipids,
        "--ratio", ratio,
        "--dist", f"{config.xy_distance:g}",
        "--dist_wat", f"{config.water_thickness:g}",
        "--salt",
        "--salt_c", _packmol_ion_name(config.positive_ion),
        "--salt_a", _packmol_ion_name(config.negative_ion),
        "--saltcon", f"{config.salt_concentration:g}",
    ]
    # Orientation: if the protein is already membrane-aligned, tell packmol-memgen
    # not to re-orient; otherwise let its built-in memembed run.
    if orient != "memembed":
        cmd.append("--preoriented")
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
            cmd.extend([protein_ff_option, protein_ff, water_ff_option, water_ff])
        cmd.append("--parametrize")
    if config.minimize:
        cmd.append("--minimize")
    return cmd


def _detect_forcefield_options(timeout: int = 15):
    """Detect force-field option spellings used by the installed release.

    AmberTools releases in the wild use both ``--fprot`` and ``--ffprot``.
    Refuse to guess for a non-default force-field request because a silently
    ignored option would create a scientifically different system.
    """
    try:
        proc = subprocess.run(
            ["packmol-memgen", "--help"],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"Could not inspect packmol-memgen --help: {exc}") from exc

    help_text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    protein_option = next((flag for flag in ("--fprot", "--ffprot") if flag in help_text), None)
    water_option = next((flag for flag in ("--ffwat", "--fwat") if flag in help_text), None)
    if not protein_option or not water_option:
        raise RuntimeError(
            "Installed packmol-memgen does not advertise compatible protein/water "
            "force-field options in --help; cannot safely request a non-default pair."
        )
    return protein_option, water_option


def _ensure_legacy_lipid_leaprc(config: MembraneConfig, work_dir: str, amberhome) -> Optional[str]:
    """Make a relocated legacy lipid ``leaprc`` available to tleap.

    AmberTools 24 archived older lipid force fields (Lipid17/14/11): their
    ``leaprc`` files moved from ``dat/leap/cmd/`` into ``dat/leap/cmd/oldff/``,
    which tleap does not search by default, so packmol-memgen's
    ``source leaprc.lipid17`` fails deep inside tleap even though the parameter
    (``.dat``) and library (``.lib``) files are still present.

    Preferred fix: restore the ``leaprc`` next to the current lipid leaprc files
    (``dat/leap/cmd/``) so tleap finds it on its default path. This is a one-time,
    idempotent copy of a file AmberTools already ships, and it leaves the build
    directory completely untouched. If that location is not writable (e.g. a
    read-only system Amber install), fall back to staging a temporary copy in the
    working directory (tleap searches the current directory first).

    Returns a path the caller must remove after the run (the temporary work-dir
    copy), or ``None`` when no cleanup is needed (already present, or restored to
    the canonical Amber directory).
    """
    lipid_ff = (getattr(config, "lipid_ff", "") or "").strip().lower()
    if not lipid_ff or not amberhome:
        return None
    leaprc = f"leaprc.{lipid_ff}"
    cmd_dir = os.path.join(amberhome, "dat", "leap", "cmd")
    canonical = os.path.join(cmd_dir, leaprc)
    if os.path.isfile(canonical):
        return None  # already on tleap's default search path
    legacy = os.path.join(cmd_dir, "oldff", leaprc)
    if not os.path.isfile(legacy):
        return None

    # Preferred: restore alongside the current lipid leaprc files (persistent).
    try:
        shutil.copy2(legacy, canonical)
        return None
    except OSError:
        pass

    # Fallback: temporary copy in the working directory, removed after the run.
    destination = os.path.join(work_dir, leaprc)
    if not os.path.exists(destination):
        try:
            shutil.copy2(legacy, destination)
            return destination
        except OSError:
            return None
    return None


def run(protein_pdb: str, config: MembraneConfig, work_dir: str, timeout: int = 7200) -> PackmolMemgenResult:
    """Run packmol-memgen in ``work_dir``. Returns paths to the produced files.

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
    os.makedirs(work_dir, exist_ok=True)
    protein_ff, water_ff = config.packmol_forcefield_pair()
    if (protein_ff, water_ff) == ("ff14SB", "tip3p"):
        forcefield_options = ("--fprot", "--ffwat")
    else:
        forcefield_options = _detect_forcefield_options()
    cmd = build_command(
        protein_pdb,
        config,
        parametrize=True,
        protein_ff_option=forcefield_options[0],
        water_ff_option=forcefield_options[1],
    )
    before_run = _file_signatures(work_dir)

    # packmol-memgen locates the `packmol` executable via $AMBERHOME/bin (it does
    # NOT search PATH). If AMBERHOME is unset but packmol is installed (e.g. in a
    # conda env), derive AMBERHOME from packmol's location so membrane builds work
    # out of the box rather than failing with "Packmol path defined but not found".
    env = os.environ.copy()
    if not env.get("AMBERHOME"):
        pk = shutil.which("packmol")
        if pk:
            env["AMBERHOME"] = os.path.dirname(os.path.dirname(pk))

    # AmberTools 24 relocated legacy lipid leaprc files (e.g. leaprc.lipid17)
    # into dat/leap/cmd/oldff/, which tleap does not search by default, so
    # packmol-memgen's `source leaprc.lipid17` fails deep inside tleap. Restore
    # it onto tleap's default search path (or, if that is read-only, stage a
    # temporary copy in the working directory that is removed afterwards).
    staged_leaprc = _ensure_legacy_lipid_leaprc(config, work_dir, env.get("AMBERHOME"))

    # Guard the subprocess: a real pack+parametrize+minimize can exceed the
    # timeout, and the OS may raise (e.g. ENOMEM). Return a failed result rather
    # than propagating, so the builder degrades gracefully.
    try:
        proc = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, timeout=timeout, env=env)
    except subprocess.TimeoutExpired as exc:
        return PackmolMemgenResult(prmtop=None, inpcrd=None, packed_pdb=None,
                                   log=f"packmol-memgen timed out after {timeout}s: {exc}", success=False)
    except OSError as exc:
        return PackmolMemgenResult(prmtop=None, inpcrd=None, packed_pdb=None,
                                   log=f"packmol-memgen failed to launch: {exc}", success=False)
    finally:
        # The staged leaprc is only needed while tleap runs inside packmol-memgen;
        # remove it so the build directory keeps only real build artifacts.
        if staged_leaprc and os.path.isfile(staged_leaprc):
            try:
                os.remove(staged_leaprc)
            except OSError:
                pass
    log = (proc.stdout or "") + "\n" + (proc.stderr or "")

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

    success = proc.returncode == 0 and prmtop is not None and inpcrd is not None
    if proc.returncode != 0:
        # A failed process may leave a half-written pair. Do not expose those
        # files as usable AMBER inputs merely because their names match.
        prmtop = None
        inpcrd = None
    return PackmolMemgenResult(prmtop=prmtop, inpcrd=inpcrd, packed_pdb=packed, log=log, success=success)


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
