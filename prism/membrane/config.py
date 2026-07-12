#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Membrane-build configuration.

Holds the lipid composition, orientation source, and bilayer/solvent box
settings for an all-atom protein-in-membrane build. Parsed from CLI flags or a
YAML ``membrane:`` section.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional


# Common lipids understood by PACKMOL-Memgen (AMBER Lipid21 names). This is a
# convenience whitelist for validation/help; PACKMOL-Memgen supports many more.
COMMON_LIPIDS = {
    "POPC", "POPE", "POPS", "POPG", "POPA", "POPI",
    "DOPC", "DOPE", "DOPS", "DPPC", "DMPC", "DLPC",
    "CHL1", "CHOL", "PSM",  # cholesterol / sphingomyelin
}

# Lipid force fields and the protein-FF families they pair with.
LIPID_FF_FAMILY = {
    "lipid21": "amber",     # AMBER Lipid21 (packmol-memgen default)
    "slipids": "amber",     # Slipids-2020 (AMBER-compatible)
    "charmm36": "charmm",   # CHARMM36 lipids
}


@dataclass
class MembraneConfig:
    """Configuration for building a protein-membrane system."""

    enabled: bool = False
    #: Lipid composition (one or more lipid names).
    lipids: List[str] = field(default_factory=lambda: ["POPC"])
    #: Molar ratio per lipid (parallel to ``lipids``); defaults to equal parts.
    ratio: List[int] = field(default_factory=list)
    #: Orientation source: "opm" (fetch by PDB id), "ppm" (local/server),
    #: "preoriented" (input already membrane-aligned), or "memembed".
    orient: str = "preoriented"
    #: PDB id to fetch a pre-oriented structure from the OPM database.
    pdb_id: Optional[str] = None
    #: Lipid force field (selects the protein-FF family pairing).
    lipid_ff: str = "lipid21"
    #: Salt concentration (mol/L) for the aqueous phase.
    salt_concentration: float = 0.15
    positive_ion: str = "Na+"
    negative_ion: str = "Cl-"
    #: Water-layer thickness above/below the bilayer (Angstrom; --dist_wat).
    water_thickness: float = 17.5
    #: Minimum protein-to-box distance in the membrane plane (Angstrom; --dist).
    xy_distance: float = 17.0
    #: Run a short AMBER (sander/pmemd) minimization on the packed system
    #: before conversion (works around GROMACS-crash regressions).
    minimize: bool = True
    temperature: float = 303.15

    def family(self) -> str:
        """Protein force-field family implied by the lipid force field."""
        return LIPID_FF_FAMILY.get(self.lipid_ff.lower(), "amber")

    def resolved_ratio(self) -> List[int]:
        if self.ratio and len(self.ratio) == len(self.lipids):
            return self.ratio
        return [1] * len(self.lipids)

    def validate(self) -> List[str]:
        problems = []
        if not self.lipids:
            problems.append("No lipids specified for membrane build.")
        for lip in self.lipids:
            if lip.upper() not in COMMON_LIPIDS:
                problems.append(f"Lipid '{lip}' is not in the common-lipid whitelist (may still be valid for PACKMOL-Memgen).")
        if self.lipid_ff.lower() not in LIPID_FF_FAMILY:
            problems.append(f"Unknown lipid force field '{self.lipid_ff}'.")
        if self.orient == "opm" and not self.pdb_id:
            problems.append("orient=opm requires a --membrane-pdbid to fetch from the OPM database.")
        return problems


def parse_membrane_cli(
    membrane: bool,
    lipid: Optional[List[str]],
    lipid_ratio: Optional[List[int]],
    orient: Optional[str],
    pdb_id: Optional[str],
    lipid_ff: Optional[str],
    salt_concentration: float = 0.15,
    temperature: float = 303.15,
) -> MembraneConfig:
    """Build a :class:`MembraneConfig` from CLI arguments."""
    lipids = [l.upper() for l in (lipid or ["POPC"])]
    return MembraneConfig(
        enabled=bool(membrane),
        lipids=lipids,
        ratio=list(lipid_ratio or []),
        orient=(orient or "preoriented"),
        pdb_id=pdb_id,
        lipid_ff=(lipid_ff or "lipid21"),
        salt_concentration=salt_concentration,
        temperature=temperature,
    )


def parse_membrane_yaml(cfg: Optional[dict], temperature: float = 303.15) -> MembraneConfig:
    """Build a :class:`MembraneConfig` from a parsed YAML ``membrane:`` mapping."""
    if not cfg:
        return MembraneConfig(enabled=False)
    lipids = cfg.get("lipids") or cfg.get("lipid") or ["POPC"]
    if isinstance(lipids, str):
        lipids = [lipids]
    return MembraneConfig(
        enabled=bool(cfg.get("enabled", True)),
        lipids=[str(l).upper() for l in lipids],
        ratio=list(cfg.get("ratio", []) or []),
        orient=str(cfg.get("orient", "preoriented")),
        pdb_id=cfg.get("pdb_id"),
        lipid_ff=str(cfg.get("lipid_ff", "lipid21")),
        salt_concentration=float(cfg.get("salt_concentration", 0.15)),
        water_thickness=float(cfg.get("water_thickness", 17.5)),
        xy_distance=float(cfg.get("xy_distance", 17.0)),
        minimize=bool(cfg.get("minimize", True)),
        temperature=float(cfg.get("temperature", temperature)),
    )
