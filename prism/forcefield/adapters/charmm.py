#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Shared CHARMM/CGenFF adapter utilities."""

from pathlib import Path
from typing import Optional


class CGenFFAdapter:
    """Adapter for CGenFF/CHARMM force field special handling.

    CGenFF has two main formats:
    1. CGenFF website format: Contains ``*_ffbonded.itp`` files with partial parameters
    2. CHARMM-GUI format: Contains complete ``charmm36.itp`` file with all parameters
    """

    CHARMM_FF_PATTERNS = ["charmm*.ff", "CHARMM*.ff"]

    @classmethod
    def detect(cls, lig_ff_path: Path) -> bool:
        """Detect whether a ligand force-field directory contains CHARMM assets."""
        for pattern in cls.CHARMM_FF_PATTERNS:
            if list(lig_ff_path.glob(pattern)):
                return True
        return False

    @classmethod
    def find_charmm_ff_dir(cls, lig_ff_path: Path) -> Optional[Path]:
        """Return the nested CHARMM force-field directory when present."""
        for pattern in cls.CHARMM_FF_PATTERNS:
            matches = list(lig_ff_path.glob(pattern))
            if matches:
                return matches[0]
        return None

    @classmethod
    def detect_format(cls, charmm_ff_dir: Path) -> str:
        """Detect whether the directory looks like CHARMM-GUI or CGenFF website output."""
        charmm36_itp = charmm_ff_dir / "charmm36.itp"
        if charmm36_itp.exists() and charmm36_itp.stat().st_size > 1000:
            return "charmm-gui"

        ffbonded_files = list(charmm_ff_dir.glob("*_ffbonded.itp"))
        if ffbonded_files and any(f.stat().st_size > 200 for f in ffbonded_files):
            return "cgenff-website"

        return "unknown"

    @classmethod
    def should_include_ligand_atomtypes(cls, ff_name: str, lig_ff_path: Path) -> bool:
        """Return False when CHARMM already provides the atomtypes in the main force field."""
        if ff_name.lower().startswith("charmm") and cls.detect(lig_ff_path):
            return False
        return True

    @classmethod
    def is_charmm_forcefield(cls, ff_name: str) -> bool:
        """Check whether the selected force-field family is CHARMM."""
        return ff_name.lower().startswith("charmm")
