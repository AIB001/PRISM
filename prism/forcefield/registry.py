#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Canonical forcefield registry for PRISM.

This module centralizes backend metadata that was previously spread across
builders, simulation helpers, and FEP utilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


def _normalize_token(value: str) -> str:
    """Normalize forcefield names for alias matching."""
    return "".join(ch for ch in str(value).strip().lower() if ch.isalnum())


@dataclass(frozen=True)
class ForceFieldSpec:
    """Canonical metadata for a single ligand forcefield backend."""

    key: str
    display_name: str
    generator_class_name: str
    generator_output_subdir: str
    output_subdirs: tuple[str, ...]
    aliases: tuple[str, ...] = ()
    requires_external: bool = False
    supports_fep: bool = True
    generic_atom_names: bool = False
    dependencies: tuple[str, ...] = ()
    builder_output_subdir: Optional[str] = None

    def matches(self, value: str) -> bool:
        normalized = _normalize_token(value)
        if normalized == _normalize_token(self.key):
            return True
        return normalized in {_normalize_token(alias) for alias in self.aliases}

    @property
    def primary_output_subdir(self) -> str:
        """Preferred output location for PRISM builder workflows."""
        return self.builder_output_subdir or self.generator_output_subdir


_REGISTRY: tuple[ForceFieldSpec, ...] = (
    ForceFieldSpec(
        key="gaff",
        display_name="GAFF",
        generator_class_name="GAFFForceFieldGenerator",
        generator_output_subdir="LIG.amb2gmx",
        output_subdirs=("LIG.amb2gmx",),
        dependencies=("AmberTools (antechamber, parmchk2, tleap)",),
    ),
    ForceFieldSpec(
        key="gaff2",
        display_name="GAFF2",
        generator_class_name="GAFF2ForceFieldGenerator",
        generator_output_subdir="LIG.amb2gmx",
        output_subdirs=("LIG.amb2gmx",),
        dependencies=("AmberTools (antechamber, parmchk2, tleap)",),
    ),
    ForceFieldSpec(
        key="openff",
        display_name="OpenFF",
        generator_class_name="OpenFFForceFieldGenerator",
        generator_output_subdir="LIG.openff2gmx",
        output_subdirs=("LIG.openff2gmx",),
        generic_atom_names=True,
        dependencies=("openff-toolkit", "openff-interchange"),
    ),
    ForceFieldSpec(
        key="cgenff",
        display_name="CGenFF",
        generator_class_name="CGenFFForceFieldGenerator",
        generator_output_subdir="LIG.cgenff2gmx",
        output_subdirs=("LIG.cgenff2gmx",),
        requires_external=True,
        dependencies=("Web download from https://cgenff.com/",),
    ),
    ForceFieldSpec(
        key="charmm-gui",
        display_name="CHARMM-GUI",
        generator_class_name="CHARMMGUIForceFieldGenerator",
        generator_output_subdir="LIG.charmm2gmx",
        output_subdirs=("LIG.charmm2gmx",),
        aliases=("charmm_gui",),
        requires_external=True,
        dependencies=("CHARMM-GUI web server (http://charmm-gui.org)",),
    ),
    ForceFieldSpec(
        key="rtf",
        display_name="RTF",
        generator_class_name="RTFForceFieldGenerator",
        generator_output_subdir="LIG.rtf2gmx",
        output_subdirs=("LIG.rtf2gmx",),
        requires_external=True,
        dependencies=("RTF and PRM files from CGenFF",),
    ),
    ForceFieldSpec(
        key="opls",
        display_name="OPLS-AA",
        generator_class_name="OPLSAAForceFieldGenerator",
        generator_output_subdir="LIG.opls2gmx",
        output_subdirs=("LIG.opls2gmx",),
        aliases=("opls-aa",),
        requires_external=True,
        dependencies=("requests", "rdkit (optional, for alignment)"),
    ),
    ForceFieldSpec(
        key="mmff",
        display_name="MMFF",
        generator_class_name="MMFFForceFieldGenerator",
        generator_output_subdir="LIG.mmff2gmx",
        output_subdirs=("LIG.sp2gmx", "LIG.mmff2gmx"),
        requires_external=True,
        dependencies=("curl (command-line tool)",),
        builder_output_subdir="LIG.sp2gmx",
    ),
    ForceFieldSpec(
        key="match",
        display_name="MATCH",
        generator_class_name="MATCHForceFieldGenerator",
        generator_output_subdir="LIG.match2gmx",
        output_subdirs=("LIG.sp2gmx", "LIG.match2gmx"),
        requires_external=True,
        dependencies=("curl (command-line tool)",),
        builder_output_subdir="LIG.sp2gmx",
    ),
    ForceFieldSpec(
        key="hybrid",
        display_name="Both",
        generator_class_name="BothForceFieldGenerator",
        generator_output_subdir="LIG.both2gmx",
        output_subdirs=("LIG.sp2gmx", "LIG.both2gmx"),
        aliases=("both",),
        requires_external=True,
        dependencies=("curl (command-line tool)",),
        builder_output_subdir="LIG.sp2gmx",
    ),
)


def iter_forcefield_specs() -> tuple[ForceFieldSpec, ...]:
    """Return all registered forcefield specs."""
    return _REGISTRY


def get_forcefield_spec(name: str) -> ForceFieldSpec:
    """Return the canonical spec for a forcefield key or alias."""
    for spec in _REGISTRY:
        if spec.matches(name):
            return spec
    raise KeyError(f"Unknown forcefield: {name}")


def normalize_forcefield_name(name: str) -> str:
    """Resolve a forcefield alias to its canonical key."""
    return get_forcefield_spec(name).key


def get_forcefield_output_subdirs(name: str) -> tuple[str, ...]:
    """Return all recognized output subdirectories for a backend."""
    return get_forcefield_spec(name).output_subdirs


def get_primary_forcefield_output_subdir(name: str) -> str:
    """Return the preferred output subdirectory for a backend."""
    return get_forcefield_spec(name).primary_output_subdir


def iter_known_ligand_output_subdirs() -> tuple[str, ...]:
    """Return all known ligand output subdirs without duplicates."""
    ordered: list[str] = []
    for spec in _REGISTRY:
        for subdir in spec.output_subdirs:
            if subdir not in ordered:
                ordered.append(subdir)
    return tuple(ordered)


def iter_existing_ligand_output_dirs(base_dir: str | Path, forcefield: Optional[str] = None) -> tuple[Path, ...]:
    """Return existing ligand output directories under a base path."""
    base = Path(base_dir)
    candidate_subdirs = (
        get_forcefield_output_subdirs(forcefield) if forcefield is not None else iter_known_ligand_output_subdirs()
    )
    discovered: list[Path] = []
    seen: set[Path] = set()
    for subdir in candidate_subdirs:
        for candidate in sorted(base.glob(f"{subdir}*")):
            if candidate.is_dir() and candidate not in seen:
                discovered.append(candidate)
                seen.add(candidate)
    return tuple(discovered)


def resolve_ligand_artifact(base_dir: str | Path, filename: str, forcefield: Optional[str] = None) -> Optional[Path]:
    """
    Resolve a PRISM ligand artifact under a base directory.

    This checks both the direct path and the known backend output subdirectories.
    """
    base = Path(base_dir)
    direct = base / filename
    if direct.exists():
        return direct

    candidate_subdirs: Iterable[str]
    if forcefield is None:
        candidate_subdirs = iter_known_ligand_output_subdirs()
    else:
        candidate_subdirs = get_forcefield_output_subdirs(forcefield)

    for subdir in candidate_subdirs:
        candidate = base / subdir / filename
        if candidate.exists():
            return candidate
    return None


def resolve_ligand_artifact_dir(
    base_dir: str | Path, filename: str, forcefield: Optional[str] = None
) -> Optional[Path]:
    """Resolve the directory containing a ligand artifact."""
    resolved = resolve_ligand_artifact(base_dir, filename, forcefield=forcefield)
    if resolved is None:
        return None
    return resolved.parent


def discover_ligand_directory(base_dir: str | Path, forcefield: Optional[str] = None) -> Optional[Path]:
    """
    Discover the best ligand artifact directory below a base path.

    Prefers canonical PRISM layouts, then falls back to any immediate child
    containing a `.itp` file for older or non-standard fixtures.
    """
    resolved = resolve_ligand_artifact_dir(base_dir, "LIG.itp", forcefield=forcefield)
    if resolved is not None:
        return resolved

    base_path = Path(base_dir)
    itp_files = list(base_path.glob("*/*.itp"))
    if itp_files:
        return itp_files[0].parent
    return None
