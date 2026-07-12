#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Membrane system builder orchestration.

Pipeline (Tier-A, AMBER-consistent, reuses PRISM's amb2gmx philosophy):

    orient (OPM/PPM/preoriented)
      -> PACKMOL-Memgen (Lipid21 bilayer + water + ions, --parametrize)
      -> ParmEd  AMBER prmtop/crd  ->  GROMACS .top/.gro
      -> membrane MDP protocol (semiisotropic, grouped thermostats)

Everything that touches an external tool is runtime-gated with an actionable
capability report, so importing/constructing the builder never fails.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from typing import List, Optional

from . import packmol_memgen
from .config import MembraneConfig
from .membrane_mdp import write_membrane_mdps
from .orient import orient_protein


@dataclass
class MembraneBuildResult:
    system_dir: str
    oriented_pdb: Optional[str] = None
    gro: Optional[str] = None
    top: Optional[str] = None
    ndx: Optional[str] = None  # index.ndx with SOLU/MEMB/SOLV groups (for grompp -n)
    mdps: dict = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    success: bool = False
    note: str = ""


class MembraneBuilder:
    """Build an all-atom protein-in-bilayer system for GROMACS."""

    def __init__(self, config: MembraneConfig, verbose: bool = True):
        self.config = config
        self.verbose = verbose

    # ------------------------------------------------------------------ #
    def capability_report(self) -> List[str]:
        """Return a list of missing-capability messages (empty if all present)."""
        missing = []
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

    @staticmethod
    def _parmed_available() -> bool:
        try:
            import parmed  # noqa: F401

            return True
        except Exception:
            return False

    # ------------------------------------------------------------------ #
    def build(self, protein_pdb: str, output_dir: str) -> MembraneBuildResult:
        """Run the full membrane build. Returns a :class:`MembraneBuildResult`."""
        system_dir = os.path.join(output_dir, "GMX_PROLIG_MEMB")
        os.makedirs(system_dir, exist_ok=True)
        result = MembraneBuildResult(system_dir=system_dir)
        result.warnings.extend(self.config.validate())

        # 0) Capability gate.
        missing = self.capability_report()
        if missing:
            result.warnings.extend(missing)
            result.note = (
                "Membrane build prerequisites missing — emitted orientation + MDP scaffold only. "
                "Install the tools above (or run from the AmberTools env) to complete the bilayer build."
            )
            # Still produce what we can without external tools: orientation + MDPs.
            self._orient_only(protein_pdb, system_dir, result)
            self._emit_mdps(output_dir, result)
            if self.verbose:
                self._print(result)
            return result

        # 1) Orient + 2) build bilayer. Guarded so an external-tool failure
        # (timeout, OSError, etc.) degrades to a scaffold rather than crashing.
        try:
            oriented = os.path.join(system_dir, "protein_oriented.pdb")
            info = orient_protein(protein_pdb, oriented, self.config.orient, self.config.pdb_id)
            result.oriented_pdb = info["oriented_pdb"]
            if info.get("note"):
                result.warnings.append(info["note"])

            pack = packmol_memgen.run(result.oriented_pdb, self.config, work_dir=system_dir)
        except Exception as exc:
            result.warnings.append(f"Membrane build step failed: {exc}")
            self._emit_mdps(output_dir, result)
            result.note = "Membrane build step errored; emitted orientation + MDP scaffold only."
            if self.verbose:
                self._print(result)
            return result

        if not pack.success:
            result.warnings.append("packmol-memgen did not produce a parametrized system (see log).")
            with open(os.path.join(system_dir, "packmol_memgen.log"), "w") as fh:
                fh.write(pack.log)
            self._emit_mdps(output_dir, result)
            result.note = "packmol-memgen failed; see packmol_memgen.log."
            if self.verbose:
                self._print(result)
            return result

        # 3) Convert AMBER -> GROMACS via ParmEd.
        gro = os.path.join(system_dir, "system.gro")
        top = os.path.join(system_dir, "topol.top")
        ndx = os.path.join(system_dir, "index.ndx")
        ok = self._amber_to_gromacs(pack.prmtop, pack.inpcrd, top, gro, result, ndx_path=ndx)
        if ok:
            result.gro, result.top = gro, top

        # 4) Membrane MDP protocol.
        self._emit_mdps(output_dir, result)

        result.success = bool(result.gro and result.top)
        result.note = "Membrane system built." if result.success else "Conversion incomplete; see warnings."
        if self.verbose:
            self._print(result)
        return result

    # ------------------------------------------------------------------ #
    def _orient_only(self, protein_pdb, system_dir, result: MembraneBuildResult):
        oriented = os.path.join(system_dir, "protein_oriented.pdb")
        try:
            info = orient_protein(protein_pdb, oriented, self.config.orient, self.config.pdb_id)
            result.oriented_pdb = info["oriented_pdb"]
            if info.get("note"):
                result.warnings.append(info["note"])
        except Exception as exc:
            result.warnings.append(f"Orientation step failed: {exc}")

    def _emit_mdps(self, output_dir, result: MembraneBuildResult):
        try:
            result.mdps = write_membrane_mdps(
                output_dir,
                temperature=self.config.temperature,
                ff_family=self.config.family(),
                production_ns=100.0,
            )
        except Exception as exc:
            result.warnings.append(f"Membrane MDP generation failed: {exc}")

    def _amber_to_gromacs(self, prmtop, inpcrd, out_top, out_gro, result: MembraneBuildResult, ndx_path=None) -> bool:
        """Convert an AMBER prmtop/inpcrd to GROMACS .top/.gro via ParmEd.

        Also writes an index.ndx with SOLU/MEMB/SOLV groups (referenced by the
        membrane MDP thermostat groups) classified from the ParmEd residues.
        """
        try:
            import parmed as pmd

            parm = pmd.load_file(prmtop, xyz=inpcrd)
            parm.save(out_top, format="gromacs", overwrite=True)
            parm.save(out_gro, overwrite=True)
            if ndx_path:
                try:
                    self._write_membrane_index(parm, ndx_path, result)
                    result.ndx = ndx_path
                except Exception as exc:  # index is non-fatal
                    result.warnings.append(f"Index (SOLU/MEMB/SOLV) generation failed: {exc}")
            return True
        except Exception as exc:
            result.warnings.append(
                f"ParmEd conversion failed ({exc}). As a fallback run: "
                f"parmed -i <(echo 'gromacs') or use acpype on {os.path.basename(prmtop)}."
            )
            return False

    # residue-name classification for the SOLU/MEMB/SOLV thermostat groups
    _SOLV_RESNAMES = {
        "WAT", "HOH", "SOL", "TIP3", "T3P", "TIP", "SPC",
        "NA", "NA+", "CL", "CL-", "K", "K+", "MG", "CA", "ZN", "SOD", "CLA", "POT",
    }
    _PROTEIN_RESNAMES = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
        "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM", "ASH", "GLH", "LYN",
        "SEP", "TPO", "PTR", "MLZ", "MLY", "M3L", "ALY", "2MR",
    }

    def _write_membrane_index(self, parm, ndx_path: str, result: "MembraneBuildResult"):
        """Write a GROMACS index file with SOLU / MEMB / SOLV groups.

        SOLU = protein (incl. PTM residues); SOLV = water + ions; MEMB =
        everything else (lipids). Atom indices are 1-based, in topology order.
        """
        solu, memb, solv = [], [], []
        for i, atom in enumerate(parm.atoms, start=1):
            rn = atom.residue.name.upper()
            if rn in self._PROTEIN_RESNAMES:
                solu.append(i)
            elif rn in self._SOLV_RESNAMES:
                solv.append(i)
            else:
                memb.append(i)

        def _block(name, idxs):
            lines = [f"[ {name} ]"]
            for k in range(0, len(idxs), 15):
                lines.append(" ".join(str(x) for x in idxs[k:k + 15]))
            return "\n".join(lines) + "\n"

        with open(ndx_path, "w") as fh:
            for name, idxs in (("SOLU", solu), ("MEMB", memb), ("SOLV", solv)):
                fh.write(_block(name, idxs))

        if not memb:
            result.warnings.append(
                "Index: no lipid (MEMB) atoms detected — verify lipid residue names; "
                "the membrane MDP tc-grps (SOLU MEMB SOLV) needs a non-empty MEMB group."
            )
        if self.verbose:
            print(f"  index.ndx: SOLU={len(solu)} MEMB={len(memb)} SOLV={len(solv)} atoms")

    def _print(self, result: MembraneBuildResult):
        print("\n[membrane] " + (result.note or ""))
        if result.oriented_pdb:
            print(f"  oriented: {result.oriented_pdb}")
        if result.gro:
            print(f"  structure: {result.gro}")
        if result.top:
            print(f"  topology:  {result.top}")
        if result.ndx:
            print(f"  index:     {result.ndx}  (use: gmx grompp ... -n index.ndx)")
        if result.mdps:
            print(f"  mdps: {', '.join(sorted(result.mdps))} in {os.path.dirname(next(iter(result.mdps.values())))}")
        for w in result.warnings:
            print(f"  WARNING: {w}")
