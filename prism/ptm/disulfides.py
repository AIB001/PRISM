#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Disulfide-bond detection for GROMACS ``pdb2gmx``.

Disulfides are native to ``pdb2gmx``: when two cysteine ``SG`` atoms lie within
the ``specbond.dat`` reference distance (0.20 nm ± 10%) and the cysteines are
named so the force field recognises the bonded variant (``CYS2`` in CHARMM,
``CYX`` in AMBER), ``pdb2gmx`` forms the SS bond automatically. This module
detects candidate pairs geometrically (so the user gets an explicit, auditable
list) and works out the residue renaming / ``pdb2gmx`` flags required.

The PDB parser here is intentionally dependency-free (no mdtraj/Bio.PDB) so the
detection can run early in protein preparation.

This module deliberately prints nothing. It is a library layer, and only the
caller knows whether a finding is worth a line of console output and at what
weight; the staging layer renders the results (see :mod:`prism.ptm.stager`).
What the module does own is its *error* text, which follows the rest of PRISM in
naming the stage that needed a file rather than leaving the user with a bare
errno and a path.
"""

from __future__ import annotations

import math
import os
import shutil
from dataclasses import dataclass
from typing import List, Optional, Tuple

#: Lower bound on an SG-SG separation that can be a real disulfide. Anything
#: shorter is a clash or a duplicated altloc, not a bond.
_SG_SG_MIN = 1.5  # Angstrom
#: Upper bound. ``specbond.dat`` gives the reference SG-SG length as 0.20 nm and
#: ``pdb2gmx`` accepts it within 10% (1.8-2.2 A); this window is intentionally
#: wider so slightly distorted crystal geometries are still *reported* to the
#: user, while non-bonded pairs stay excluded. The consequence is that a pair
#: found between 2.2 and 2.5 A is a candidate PRISM lists but ``pdb2gmx`` may
#: still decline to bond on the force fields that rely on ``specbond.dat``
#: (CHARMM); on AMBER/OPLS the CYX rename forces the bond regardless.
_SG_SG_MAX = 2.5  # Angstrom (0.20 nm + ~25% — pragmatic upper bound)


@dataclass(frozen=True)
class CysSG:
    chain: str
    resid: int
    resname: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class Disulfide:
    a: CysSG
    b: CysSG
    distance: float

    def describe(self) -> str:
        return (
            f"{self.a.resname} {self.a.chain}{self.a.resid}"
            f" <-> {self.b.resname} {self.b.chain}{self.b.resid}"
            f" (SG-SG {self.distance:.2f} A)"
        )


def _require_input_pdb(pdb_path: str, purpose: str) -> None:
    """Fail with a message that names the stage, not just an errno and a path.

    ``open()`` on a missing file reports only the path, which inside a
    multi-stage build tells the user nothing about which step went wrong. The
    rest of PRISM always names the role of the file (``prism/core.py``:
    "Protein file not found: ..."), so disulfide handling does too.
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(
            f"PDB file for {purpose} not found: {pdb_path}. "
            "Disulfide handling reads the cleaned protein produced by the "
            "protein-preparation step; check that step completed, or pass an "
            "existing PDB."
        )


def _require_output_dir(out_path: str, purpose: str) -> None:
    """Check the destination directory before any work is done."""
    parent = os.path.dirname(os.path.abspath(out_path))
    if not os.path.isdir(parent):
        raise FileNotFoundError(
            f"Output directory for {purpose} does not exist: {parent}. "
            "Create it (or point the build at an existing output directory) "
            "before staging disulfides."
        )


def _parse_cys_sg(pdb_path: str) -> List[CysSG]:
    """Extract all cysteine SG atoms from a PDB file."""
    _require_input_pdb(pdb_path, "disulfide detection")
    out: List[CysSG] = []
    cys_names = {"CYS", "CYX", "CYM", "CYS2"}
    with open(pdb_path, "r") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "SG":
                continue
            resname = line[17:20].strip()
            if resname not in cys_names:
                continue
            chain = line[21].strip() or "A"
            try:
                resid = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            out.append(CysSG(chain, resid, resname, x, y, z))
    return out


def _dist(a: CysSG, b: CysSG) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def detect_disulfides(pdb_path: str) -> List[Disulfide]:
    """Detect disulfide bonds from cysteine SG-SG distances.

    Each SG is matched to at most one partner (its nearest in-range SG), which
    guards against the chained ``SG-SG-SG`` angle bug that crashes ``grompp``.

    A returned pair is a *candidate*: the accepted distance window is slightly
    wider than the one ``pdb2gmx`` applies to ``specbond.dat`` (see
    :data:`_SG_SG_MAX`). On force fields where the bond is formed from geometry
    rather than from the residue name (CHARMM), a pair near the upper bound can
    still be declined by ``pdb2gmx``, so callers should not report these as
    bonds that are guaranteed to exist in the final topology.
    """
    sgs = _parse_cys_sg(pdb_path)
    n = len(sgs)
    if n < 2:
        return []

    # candidate pairs within range, sorted by distance (greedy nearest-first)
    candidates: List[Tuple[float, int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            d = _dist(sgs[i], sgs[j])
            if _SG_SG_MIN <= d <= _SG_SG_MAX:
                candidates.append((d, i, j))
    candidates.sort()

    used = set()
    bonds: List[Disulfide] = []
    for d, i, j in candidates:
        if i in used or j in used:
            continue
        used.add(i)
        used.add(j)
        bonds.append(Disulfide(sgs[i], sgs[j], d))
    return bonds


def bonded_cys_resname(ff_family: str) -> Optional[str]:
    """Return the PDB residue name a *bonded* cysteine should carry, or None.

    * AMBER / OPLS use ``CYX`` (a valid 3-character residue name) for the
      disulfide-bonded cystine; renaming it in the input PDB is the established,
      non-interactive way to tell ``pdb2gmx`` the cysteine is in a disulfide.
    * CHARMM's bonded cysteine residue is ``CYS2`` *internally*, but that name is
      4 characters and is **not** a valid PDB residue name — writing it into the
      PDB corrupts the fixed-width columns. For CHARMM we therefore keep the
      residue as ``CYS`` and let ``pdb2gmx`` form the bond via ``specbond.dat``
      (distance-based), so this returns ``None`` (no rename).
    """
    fam = (ff_family or "").lower()
    if "charmm" in fam:
        return None
    # AMBER / OPLS use CYX (3 chars) for the disulfide-bonded cystine
    return "CYX"


def apply_disulfide_renaming(
    input_pdb: str,
    output_pdb: str,
    bonds: List[Disulfide],
    ff_family: str,
) -> int:
    """Rename the cysteines participating in disulfides to the bonded variant.

    Returns the number of residues renamed. For force fields whose bonded
    cysteine name is not a valid 3-character PDB residue name (e.g. CHARMM's
    ``CYS2``), no renaming is performed and ``pdb2gmx`` + ``specbond.dat`` forms
    the bond from the SG-SG geometry instead.
    """
    _require_input_pdb(input_pdb, "disulfide renaming")
    _require_output_dir(output_pdb, "disulfide renaming")

    target = bonded_cys_resname(ff_family)

    # No rename needed/possible (CHARMM, or an over-length name): copy through.
    if not bonds or target is None or len(target) > 3:
        if input_pdb != output_pdb:
            shutil.copy2(input_pdb, output_pdb)
        return 0

    to_rename = set()
    for b in bonds:
        to_rename.add((b.a.chain, b.a.resid))
        to_rename.add((b.b.chain, b.b.resid))

    renamed = 0
    with open(input_pdb, "r") as fh:
        lines = fh.readlines()

    new_lines = []
    seen = set()
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname in {"CYS", "CYX", "CYS2"}:
                chain = line[21].strip() or "A"
                try:
                    resid = int(line[22:26])
                except ValueError:
                    resid = None
                if resid is not None and (chain, resid) in to_rename:
                    line = line[:17] + f"{target:>3s}" + line[20:]
                    key = (chain, resid)
                    if key not in seen:
                        seen.add(key)
                        renamed += 1
        new_lines.append(line)

    with open(output_pdb, "w") as fh:
        fh.writelines(new_lines)
    return renamed
