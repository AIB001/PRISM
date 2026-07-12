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
import subprocess
from dataclasses import dataclass
from typing import List, Optional

from .config import MembraneConfig


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


def build_command(protein_pdb: str, config: MembraneConfig, parametrize: bool = True) -> List[str]:
    """Construct the packmol-memgen command line for the given config.

    The lipid argument uses PACKMOL-Memgen's ``A:B`` syntax for mixtures and the
    parallel ``--ratio`` flag for molar proportions.
    """
    lipids = ":".join(config.lipids)
    ratio = ":".join(str(r) for r in config.resolved_ratio())

    cmd = [
        "packmol-memgen",
        "--pdb", protein_pdb,
        "--lipids", lipids,
        "--ratio", ratio,
        "--dist", f"{config.xy_distance:g}",
        "--dist_wat", f"{config.water_thickness:g}",
        "--salt",
        "--salt_c", config.positive_ion,
        "--salt_a", config.negative_ion,
        "--saltcon", f"{config.salt_concentration:g}",
    ]
    # Orientation: if the protein is already membrane-aligned, tell packmol-memgen
    # not to re-orient; otherwise let its built-in memembed run.
    if config.orient in ("preoriented", "opm", "ppm", "memembed"):
        if config.orient != "memembed":
            cmd.append("--preoriented")
    if parametrize:
        cmd.append("--parametrize")
    if config.minimize:
        cmd.append("--minimize")
    return cmd


def run(protein_pdb: str, config: MembraneConfig, work_dir: str, timeout: int = 7200) -> PackmolMemgenResult:
    """Run packmol-memgen in ``work_dir``. Returns paths to the produced files.

    Raises RuntimeError if packmol-memgen is not installed.
    """
    if not is_available():
        raise RuntimeError(
            "packmol-memgen not found. Install AmberTools (conda install -c conda-forge ambertools) "
            "or run PRISM from an environment that provides it (e.g. the AmberTools conda env)."
        )
    os.makedirs(work_dir, exist_ok=True)
    cmd = build_command(protein_pdb, config, parametrize=True)

    # packmol-memgen locates the `packmol` executable via $AMBERHOME/bin (it does
    # NOT search PATH). If AMBERHOME is unset but packmol is installed (e.g. in a
    # conda env), derive AMBERHOME from packmol's location so membrane builds work
    # out of the box rather than failing with "Packmol path defined but not found".
    env = os.environ.copy()
    if not env.get("AMBERHOME"):
        pk = shutil.which("packmol")
        if pk:
            env["AMBERHOME"] = os.path.dirname(os.path.dirname(pk))

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
    log = (proc.stdout or "") + "\n" + (proc.stderr or "")

    prmtop = _find_first(work_dir, (".prmtop", ".top"))
    # When minimization is requested, packmol-memgen writes the minimized
    # coordinates to '<base>_min.restrt' (final) / 'min.restrt' (restrained stage).
    # Prefer those over the un-minimized tleap '.crd' so the stability-fixing
    # minimization is not silently discarded.
    inpcrd = None
    if config.minimize:
        inpcrd = (
            _find_first(work_dir, ("_min.restrt",), prefix_match=False)
            or _find_first(work_dir, ("min.restrt", ".restrt", ".ncrst"))
        )
    if inpcrd is None:
        inpcrd = _find_first(work_dir, (".inpcrd", ".crd", ".rst7"))
    packed = _find_first(work_dir, ("bilayer_",), prefix_match=True, suffix=".pdb")

    success = proc.returncode == 0 and prmtop is not None and inpcrd is not None
    return PackmolMemgenResult(prmtop=prmtop, inpcrd=inpcrd, packed_pdb=packed, log=log, success=success)


def _find_first(directory: str, patterns, prefix_match: bool = False, suffix: str = "") -> Optional[str]:
    """Find the first file in ``directory`` matching the given pattern(s)."""
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return None
    for name in names:
        if prefix_match:
            if any(name.startswith(p) for p in patterns) and name.endswith(suffix):
                return os.path.join(directory, name)
        else:
            if any(name.endswith(p) for p in patterns):
                return os.path.join(directory, name)
    return None
