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
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

# specbond.dat reference SG-SG length is 0.20 nm; allow a generous tolerance to
# catch slightly distorted crystal geometries while excluding non-bonded pairs.
_SG_SG_MIN = 1.5  # Angstrom
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


def _parse_cys_sg(pdb_path: str) -> List[CysSG]:
    """Extract all cysteine SG atoms from a PDB file."""
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
    target = bonded_cys_resname(ff_family)

    # No rename needed/possible (CHARMM, or an over-length name): copy through.
    if not bonds or target is None or len(target) > 3:
        if input_pdb != output_pdb:
            import shutil

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
