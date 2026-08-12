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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

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


def _pdb_residue_key(line: str) -> Optional[Tuple[str, int]]:
    """Return a normalized ``(chain, resid)`` key for a PDB atom record."""
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        resid = int(line[22:26])
    except (ValueError, IndexError):
        return None
    chain = line[21].strip() or "A"
    return chain, resid


def _known_source_names(definition) -> Set[str]:
    """Residue names that can legitimately represent one catalogued PTM."""
    names = {definition.code, definition.parent}
    for name in (
        definition.charmm_resname,
        definition.amber_resname,
        *definition.charmm_charge_variants.values(),
        *definition.amber_charge_variants.values(),
    ):
        if name:
            names.add(name)
    return {name.upper() for name in names}


def _normalized_target_atom_name(atom_name: str, aliases: Dict[str, str]) -> str:
    """Normalize one input atom name exactly as it will be written."""
    normalized = atom_name.strip().upper()
    return aliases.get(normalized, normalized)


def promote_requested_ptm_records(input_pdb: str, output_pdb: str, config: PTMConfig) -> int:
    """Promote requested PTM ``HETATM`` records before protein cleaning.

    ``ProteinCleaner`` groups ordinary HETATM records after the protein's TER
    record. Deposited modified amino acids normally use HETATM, so requested
    PTM records must be marked as protein before cleaning to keep their original
    place in the chain. The residue name and all coordinates remain unchanged.

    Returns the number of coordinate records promoted.
    """
    targets: Dict[Tuple[str, int], Set[str]] = {}
    for request in config.requests:
        definition = get_ptm(request.code)
        if definition is None:
            continue
        key = (str(request.chain).strip() or "A", int(request.resid))
        targets.setdefault(key, set()).update(_known_source_names(definition))

    with open(input_pdb, "r") as fh:
        lines = fh.readlines()

    promoted = 0
    output_lines = []
    for line in lines:
        key = _pdb_residue_key(line)
        if line.startswith("HETATM") and key in targets:
            resname = line[17:20].strip().upper()
            if resname in targets[key]:
                line = "ATOM  " + line[6:]
                promoted += 1
        output_lines.append(line)

    with open(output_pdb, "w") as fh:
        fh.writelines(output_lines)
    return promoted


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

        # Reject every invalid or unavailable request before an output file is
        # written.  A requested structural modification must never degrade into
        # a warning followed by an unchanged build.
        self._validate_requests()

        # 1) Apply requested PTM residue renames + atom-name normalisation.
        codes_used: List[str] = []
        self._rename_ptm_residues(input_pdb, output_pdb, result, codes_used)

        # 2) Disulfides (operates on the just-written file, in place).
        self._handle_disulfides(output_pdb, result)

        # 3) FF-specific staging.
        if self.family == "charmm" and codes_used:
            result.residuetypes_path = self._stage_residuetypes(work_dir, codes_used)

        if self.verbose:
            print(result.summary())
        return result

    def _validate_requests(self) -> None:
        """Fail before staging when a requested PTM cannot be honored."""
        problems = self.config.validate()
        if problems:
            details = "; ".join(problems)
            raise ValueError(f"Invalid PTM configuration: {details}")

        seen = set()
        resolved = []
        for request in self.config.requests:
            definition = get_ptm(request.code)
            # ``validate`` catches these; retain defensive checks for callers
            # supplying custom config-like objects.
            if definition is None:
                raise ValueError(
                    f"Unknown PTM code '{request.code}' at {request.chain}{request.resid}."
                )
            if not definition.validated:
                raise ValueError(
                    f"PTM '{request.code}' has no validated parameters: {definition.note}"
                )

            key = (str(request.chain).strip() or "A", int(request.resid))
            if key in seen:
                raise ValueError(
                    f"Duplicate PTM requests for {key[0]}{key[1]}; each residue can have only one PTM."
                )
            seen.add(key)

            if self.family == "charmm":
                target = definition.charmm_name_for_charge(request.charge)
            else:
                target = definition.amber_name_for_charge(request.charge)
            if not target:
                charge = (
                    f" charge {request.charge:+d}" if request.charge is not None else ""
                )
                raise ValueError(
                    f"PTM {request.code}{charge} is not available in the {self.family} "
                    f"force-field family at {request.chain}{request.resid}."
                )
            if len(target) > 4:
                raise ValueError(
                    f"PTM residue name '{target}' exceeds the 4-character residue-name "
                    "field supported by GROMACS/pdb2gmx."
                )
            resolved.append((request, definition, target))

        # The renamed AMBER variants (S1P/T1P/Y1P) are not catalog keys, so
        # route detection must use the original definitions rather than target
        # residue names.  Run this before writing a partially staged PDB.
        if self.family in ("amber", "opls"):
            phospho = sorted(
                {definition.code for _request, definition, _target in resolved
                 if definition.category == "phosphorylation"}
            )
            if phospho:
                family_label = self.family.upper()
                raise RuntimeError(
                    f"{family_label} phosphorylation PTM(s) {', '.join(phospho)} require the "
                    f"tleap+{self.config.amber_phospho_ff} -> ParmEd/ACPYPE route, which is not "
                    "available in the default pdb2gmx build. Use --forcefield "
                    "charmm36-jul2022 for the fully automated phosphorylation path."
                )

    # ------------------------------------------------------------------ #
    def _rename_ptm_residues(self, input_pdb, output_pdb, result: PTMResult, codes_used: List[str]):
        """Rename requested residues to their PTM code; normalise atom names."""
        # Build a lookup of (chain, resid) -> target metadata.
        targets = {}
        for req in self.config.requests:
            d = get_ptm(req.code)
            if d is None:
                raise ValueError(f"Unknown PTM code '{req.code}' at {req.chain}{req.resid}.")
            key = (str(req.chain).strip() or "A", int(req.resid))
            if key in targets:
                raise ValueError(
                    f"Duplicate PTM requests for {key[0]}{key[1]}; each residue can have only one PTM."
                )
            if self.family == "charmm":
                target = d.charmm_name_for_charge(req.charge)
            else:
                target = d.amber_name_for_charge(req.charge)
            if not target:
                charge = f" charge {req.charge:+d}" if req.charge is not None else ""
                raise ValueError(
                    f"PTM {req.code}{charge} is not available in the {self.family} "
                    f"force-field family at {req.chain}{req.resid}."
                )
            if len(target) > 4:
                raise ValueError(
                    f"PTM residue name '{target}' exceeds the 4-character residue-name "
                    "field supported by GROMACS/pdb2gmx."
                )
            targets[key] = (req, d, target, d.atom_aliases_for(self.family, target))

        with open(input_pdb, "r") as fh:
            lines = fh.readlines()

        # Validate every requested location before writing anything. Previously
        # a typo in chain/residue silently produced an unchanged PDB while the
        # result claimed success; a wrong parent residue was blindly relabelled.
        source_names: Dict[Tuple[str, int], Set[str]] = {}
        insertion_codes: Dict[Tuple[str, int], Set[str]] = {}
        normalized_atom_names: Dict[Tuple[str, int], Set[str]] = {}
        for line in lines:
            key = _pdb_residue_key(line)
            if key not in targets:
                continue
            source_names.setdefault(key, set()).add(line[17:20].strip().upper())
            insertion_codes.setdefault(key, set()).add(line[26].strip() if len(line) > 26 else "")
            aliases = targets[key][3]
            atom_name = _normalized_target_atom_name(line[12:16], aliases)
            normalized_atom_names.setdefault(key, set()).add(atom_name)

        for key, (req, definition, target, _aliases) in targets.items():
            location = f"{key[0]}{key[1]}"
            if key not in source_names:
                raise ValueError(f"PTM target {location} not found in input PDB.")
            if len(insertion_codes[key]) > 1:
                raise ValueError(
                    f"PTM target {location} is ambiguous because multiple insertion codes share that residue id."
                )
            unexpected = source_names[key] - _known_source_names(definition)
            if unexpected:
                names = ", ".join(sorted(unexpected))
                expected = ", ".join(sorted(_known_source_names(definition)))
                raise ValueError(
                    f"PTM target {location} has incompatible residue name(s) {names}; expected one of {expected}."
                )
            required_atoms = (
                set(definition.required_heavy_atoms_for(target))
                if self.family == "charmm"
                else set()
            )
            missing_atoms = required_atoms - normalized_atom_names.get(key, set())
            if missing_atoms:
                missing = ", ".join(sorted(missing_atoms))
                raise ValueError(
                    f"PTM target {location} is missing required {target} heavy atom(s): {missing}."
                )
            codes_used.append(target)
            result.renamed_residues.append(f"{req.chain}{req.resid}:{definition.parent}->{target}")

        out = []
        for line in lines:
            key = _pdb_residue_key(line)
            if key in targets:
                _req, _definition, target, aliases = targets[key]
                # PTM atoms are often HETATM; pdb2gmx wants ATOM for protein.
                line = "ATOM  " + line[6:]
                # Normalize case and aliases using the same rule as preflight.
                atom_name = line[12:16].strip()
                new_name = _normalized_target_atom_name(atom_name, aliases)
                if new_name != atom_name:
                    line = line[:12] + f"{new_name:>4s}"[:4] + line[16:]
                # rename the residue
                line = line[:17] + f"{target:<4s}" + line[21:]
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
