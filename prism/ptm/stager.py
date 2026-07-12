#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PTM staging: prepare a cleaned protein PDB for ``pdb2gmx`` with modifications.

The stager runs immediately after protein cleaning and before topology
generation. It:

1. Renames requested residues to their PTM code and normalises atom names so
   ``pdb2gmx`` matches the force-field ``.rtp`` entry.
2. Detects (or accepts) disulfide bonds and renames the bonded cysteines to the
   FF's cystine residue name (``CYS2`` / ``CYX``).
3. For the CHARMM36 family — whose bundled FF already contains the PTM residues
   — stages a local ``residuetypes.dat`` registering the PTM codes as
   ``Protein`` so the backbone is not split at the modified residue.
4. For the AMBER family, reports the route (``tleap`` + ``phosaa`` →
   ParmEd/ACPYPE) and gates on AmberTools availability.

The result object tells the system builder what changed and which ``pdb2gmx``
inputs/flags to use.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from .catalog import get_ptm
from .disulfides import Disulfide, apply_disulfide_renaming, bonded_cys_resname, detect_disulfides
from .spec import PTMConfig


def _ff_family(ff_name: Optional[str]) -> str:
    """Classify a protein force-field name into 'charmm' / 'amber' / 'opls'."""
    n = (ff_name or "").lower()
    if "charmm" in n:
        return "charmm"
    if "amber" in n or n.startswith("a99") or "ff19" in n or "ff14" in n:
        return "amber"
    if "opls" in n:
        return "opls"
    return "amber"


@dataclass
class PTMResult:
    """Outcome of the PTM staging step."""

    output_pdb: str
    family: str
    renamed_residues: List[str] = field(default_factory=list)
    disulfides: List[Disulfide] = field(default_factory=list)
    residuetypes_path: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    amber_route_required: bool = False

    def summary(self) -> str:
        lines = []
        if self.renamed_residues:
            lines.append(f"PTM residues: {', '.join(self.renamed_residues)}")
        if self.disulfides:
            lines.append(f"Disulfides ({len(self.disulfides)}): " + "; ".join(b.describe() for b in self.disulfides))
        if self.residuetypes_path:
            lines.append(f"Staged residuetypes.dat: {self.residuetypes_path}")
        for w in self.warnings:
            lines.append(f"WARNING: {w}")
        return "\n".join(lines) if lines else "No PTM changes applied."


class PTMStager:
    """Apply PTM and disulfide handling to a cleaned protein PDB."""

    def __init__(
        self,
        config: PTMConfig,
        forcefield_name: str,
        forcefield_dir: Optional[str] = None,
        gmx_command: str = "gmx",
        verbose: bool = True,
    ):
        self.config = config
        self.forcefield_name = forcefield_name
        self.family = _ff_family(forcefield_name)
        self.forcefield_dir = forcefield_dir
        self.gmx_command = gmx_command
        self.verbose = verbose

    # ------------------------------------------------------------------ #
    def apply(self, input_pdb: str, output_pdb: str, work_dir: Optional[str] = None) -> PTMResult:
        """Run the full PTM staging step, returning a :class:`PTMResult`."""
        work_dir = work_dir or os.path.dirname(os.path.abspath(output_pdb)) or "."
        result = PTMResult(output_pdb=output_pdb, family=self.family)

        # Surface configuration problems early (non-fatal).
        result.warnings.extend(self.config.validate())

        # 1) Apply requested PTM residue renames + atom-name normalisation.
        codes_used: List[str] = []
        self._rename_ptm_residues(input_pdb, output_pdb, result, codes_used)

        # 2) Disulfides (operates on the just-written file, in place).
        self._handle_disulfides(output_pdb, result)

        # 3) FF-specific staging.
        if self.family == "charmm" and codes_used:
            result.residuetypes_path = self._stage_residuetypes(work_dir, codes_used)
        elif self.family in ("amber", "opls") and codes_used:
            # phospho/methyl/acetyl on AMBER need the tleap->ACPYPE route, which
            # replaces the pdb2gmx protein-topology build for the modified chain.
            phospho = [c for c in codes_used if get_ptm(c) and get_ptm(c).category == "phosphorylation"]
            if phospho:
                result.amber_route_required = True
                # The residues were just renamed to AMBER names (e.g. SEP/TPO/PTR) that
                # amber*sb .rtp cannot build, and the tleap->ParmEd/ACPYPE route is not
                # wired into the default pdb2gmx path. Fail here with an actionable message
                # rather than letting pdb2gmx fail later with an "unknown residue" error.
                raise RuntimeError(
                    "AMBER phosphorylation PTM(s) "
                    f"{', '.join(sorted(set(phospho)))} require the tleap+"
                    f"{self.config.amber_phospho_ff} -> ParmEd/ACPYPE route, which is not "
                    "available in the default pdb2gmx build. Use --forcefield "
                    "charmm36-jul2022 for the fully automated phosphorylation path."
                )

        if self.verbose:
            print(result.summary())
        return result

    # ------------------------------------------------------------------ #
    def _rename_ptm_residues(self, input_pdb, output_pdb, result: PTMResult, codes_used: List[str]):
        """Rename requested residues to their PTM code; normalise atom names."""
        # Build a lookup of (chain, resid) -> (target_resname, atom_aliases)
        targets = {}
        for req in self.config.requests:
            d = get_ptm(req.code)
            if d is None:
                continue
            if self.family == "charmm":
                target = d.charmm_name_for_charge(req.charge)
            else:
                target = d.amber_name_for_charge(req.charge)
            if not target:
                result.warnings.append(
                    f"{req.code} not available in the {self.family} family at {req.chain}{req.resid}; skipped."
                )
                continue
            targets[(req.chain, req.resid)] = (target, d.atom_aliases)
            codes_used.append(target)
            result.renamed_residues.append(f"{req.chain}{req.resid}:{d.parent}->{target}")

        with open(input_pdb, "r") as fh:
            lines = fh.readlines()

        out = []
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21].strip() or "A"
                try:
                    resid = int(line[22:26])
                except ValueError:
                    out.append(line)
                    continue
                key = (chain, resid)
                if key in targets:
                    target, aliases = targets[key]
                    # PTM atoms are often HETATM; pdb2gmx wants ATOM for protein.
                    line = "ATOM  " + line[6:]
                    # normalise atom name if an alias is defined
                    atom_name = line[12:16].strip()
                    if atom_name in aliases:
                        new_name = aliases[atom_name]
                        line = line[:12] + f"{new_name:>4s}"[:4] + line[16:]
                    # rename the residue
                    line = line[:17] + f"{target:>3s}" + line[20:]
            out.append(line)

        with open(output_pdb, "w") as fh:
            fh.writelines(out)

    # ------------------------------------------------------------------ #
    def _handle_disulfides(self, pdb_path: str, result: PTMResult):
        mode = self.config.disulfides
        if mode in ("none", None, False):
            return
        bonds = detect_disulfides(pdb_path)
        if isinstance(mode, list):
            # restrict to user-specified residue positions
            wanted = {(str(c), int(i)) for c, i in mode}
            bonds = [
                b for b in bonds
                if (b.a.chain, b.a.resid) in wanted or (b.b.chain, b.b.resid) in wanted
            ]
        if not bonds:
            return
        n = apply_disulfide_renaming(pdb_path, pdb_path, bonds, self.forcefield_name)
        result.disulfides = bonds
        if self.verbose and bonds:
            target = bonded_cys_resname(self.forcefield_name)
            if n and target:
                print(f"  Detected {len(bonds)} disulfide bond(s); renamed {n} cysteines to '{target}'")
            else:
                print(f"  Detected {len(bonds)} disulfide bond(s); pdb2gmx + specbond.dat will form "
                      f"them from SG-SG geometry (cysteines kept as CYS for {self.family}).")

    # ------------------------------------------------------------------ #
    def _stage_residuetypes(self, work_dir: str, codes: List[str]) -> Optional[str]:
        """Write a local residuetypes.dat that registers PTM codes as Protein.

        ``pdb2gmx`` reads ``residuetypes.dat`` from the working directory first,
        so this must be a *complete* copy of the global file plus the PTM codes.
        """
        global_rt = self._find_global_residuetypes()
        dest = os.path.join(work_dir, "residuetypes.dat")
        existing = set()
        lines: List[str] = []

        if global_rt and os.path.exists(global_rt):
            with open(global_rt, "r") as fh:
                for line in fh:
                    lines.append(line.rstrip("\n"))
                    parts = line.split()
                    if parts:
                        existing.add(parts[0].upper())
        else:
            # Minimal fallback: standard residue classes pdb2gmx needs.
            lines.extend(self._fallback_residuetypes())
            for ln in lines:
                p = ln.split()
                if p:
                    existing.add(p[0].upper())

        added = []
        for code in sorted(set(codes)):
            if code.upper() not in existing:
                lines.append(f"{code}    Protein")
                added.append(code)

        with open(dest, "w") as fh:
            fh.write("\n".join(lines) + "\n")

        if self.verbose and added:
            print(f"  Registered PTM residues as Protein in {dest}: {', '.join(added)}")
        return dest

    def _find_global_residuetypes(self) -> Optional[str]:
        """Locate the GROMACS global residuetypes.dat."""
        # 1) GMXLIB / GMXDATA env hints
        for env in ("GMXLIB", "GMXDATA"):
            base = os.environ.get(env)
            if base:
                for cand in (
                    os.path.join(base, "residuetypes.dat"),
                    os.path.join(base, "top", "residuetypes.dat"),
                ):
                    if os.path.exists(cand):
                        return cand
        # 2) derive from the gmx binary location: <prefix>/share/gromacs/top
        gmx = shutil.which(self.gmx_command) or shutil.which("gmx")
        if gmx:
            prefix = Path(gmx).resolve().parent.parent
            cand = prefix / "share" / "gromacs" / "top" / "residuetypes.dat"
            if cand.exists():
                return str(cand)
        return None

    @staticmethod
    def _fallback_residuetypes() -> List[str]:
        std_aa = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "CYS2", "GLN", "GLU", "GLY", "HIS",
            "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "MET",
            "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "CYX", "ASH", "GLH", "LYN",
        ]
        lines = [f"{a}    Protein" for a in std_aa]
        lines += [
            "SOL    Water", "HOH    Water", "WAT    Water",
            "NA    Ion", "CL    Ion", "K    Ion", "MG    Ion", "CA    Ion", "ZN    Ion",
        ]
        return lines

    # ------------------------------------------------------------------ #
    def amber_tleap_available(self) -> bool:
        """Return True if an AmberTools tleap is reachable."""
        return shutil.which("tleap") is not None
