#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CGenFF/CHARMM force field adapter.

This module centralizes all CGenFF/CHARMM-specific logic that was previously
scattered across multiple files (topology.py, workflow_fep.py, etc.).
"""

from pathlib import Path
from typing import Optional, List
import glob


class CGenFFAdapter:
    """Adapter for CGenFF/CHARMM force field special handling.

    CGenFF has two main formats:
    1. CGenFF website format: Contains *_ffbonded.itp files with partial parameters
    2. CHARMM-GUI format: Contains complete charmm36.itp file with all parameters

    This adapter provides unified detection and handling for both formats.
    """

    # Known CHARMM force field directory patterns
    CHARMM_FF_PATTERNS = ["charmm*.ff", "CHARMM*.ff"]

    @classmethod
    def detect(cls, lig_ff_path: Path) -> bool:
        """Detect if this is a CGenFF/CHARMM force field.

        Parameters:
        -----------
        lig_ff_path : Path
            Path to ligand force field directory

        Returns:
        --------
        bool : True if CGenFF/CHARMM format detected
        """
        # Check for charmm*.ff subdirectory
        for pattern in cls.CHARMM_FF_PATTERNS:
            if list(lig_ff_path.glob(pattern)):
                return True
        return False

    @classmethod
    def find_charmm_ff_dir(cls, lig_ff_path: Path) -> Optional[Path]:
        """Find the CHARMM force field subdirectory.

        Parameters:
        -----------
        lig_ff_path : Path
            Path to ligand force field directory

        Returns:
        --------
        Optional[Path] : Path to charmm*.ff directory, or None if not found
        """
        for pattern in cls.CHARMM_FF_PATTERNS:
            matches = list(lig_ff_path.glob(pattern))
            if matches:
                return matches[0]
        return None

    @classmethod
    def detect_format(cls, charmm_ff_dir: Path) -> str:
        """Detect CGenFF format type.

        Parameters:
        -----------
        charmm_ff_dir : Path
            Path to charmm*.ff directory

        Returns:
        --------
        str : Format type
            - "charmm-gui": Has complete charmm36.itp
            - "cgenff-website": Has *_ffbonded.itp files
            - "unknown": Neither format detected
        """
        # Check for CHARMM-GUI format (complete charmm36.itp)
        charmm36_itp = charmm_ff_dir / "charmm36.itp"
        if charmm36_itp.exists() and charmm36_itp.stat().st_size > 1000:
            return "charmm-gui"

        # Check for CGenFF website format (*_ffbonded.itp)
        ffbonded_files = list(charmm_ff_dir.glob("*_ffbonded.itp"))
        if ffbonded_files and any(f.stat().st_size > 200 for f in ffbonded_files):
            return "cgenff-website"

        return "unknown"

    @classmethod
    def should_include_ligand_atomtypes(cls, ff_name: str, lig_ff_path: Path) -> bool:
        """Determine if ligand atomtypes should be included in topology.

        For CGenFF/CHARMM, atomtypes are already in the main force field,
        so we skip including ligand-specific atomtypes.

        Parameters:
        -----------
        ff_name : str
            Main force field name (e.g., "charmm36-jul2022")
        lig_ff_path : Path
            Path to ligand force field directory

        Returns:
        --------
        bool : True if ligand atomtypes should be included
        """
        # If it's CHARMM and has charmm*.ff subdirectory, skip ligand atomtypes
        if ff_name.lower().startswith("charmm") and cls.detect(lig_ff_path):
            return False
        return True

    @classmethod
    def is_charmm_forcefield(cls, ff_name: str) -> bool:
        """Check if force field name indicates CHARMM.

        Parameters:
        -----------
        ff_name : str
            Force field name

        Returns:
        --------
        bool : True if CHARMM force field
        """
        return ff_name.lower().startswith("charmm")
