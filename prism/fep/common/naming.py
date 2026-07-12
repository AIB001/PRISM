#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM FEP System Naming Standardization

This module provides standardized naming conventions for FEP system output directories.
Ensures consistent directory naming across different force field combinations.
"""

import re
from pathlib import Path
from typing import Dict, Optional

from prism.forcefield.registry import normalize_forcefield_name


class FEPSystemNamer:
    """Standardization helper for FEP system directory names.

    Standard format: <protein_ff>-mut_<ligand_ff>
    Example: amber14sb_ol15-mut_gaff2
    """

    # Force-field name normalization map
    FORCEFIELD_NORMALIZATION = {
        # AMBER force fields
        "amber14sb_OL15": "amber14sb_ol15",
        "amber14sb": "amber14sb",
        "amber99sb-ildn": "amber99sb_ildn",
        "amber19sb": "amber19sb",
        "amber03": "amber03",
        # CHARMM force fields
        "charmm36_jul2022": "charmm36_jul2022",
        "charmm36m": "charmm36m",
        "charmm27": "charmm27",
        # OPLS force fields
        "oplsaa": "oplsaa",
        "oplsaam": "oplsaam",
        # Ligand force fields
        "gaff2": "gaff2",
        "gaff": "gaff",
        "openff": "openff",
        "cgenff": "cgenff",
        "opls": "opls",
    }

    @classmethod
    def generate_name(cls, protein_ff: str, ligand_ff: str) -> str:
        """Generate the standard FEP system directory name.

        Args:
            protein_ff: Protein force-field name (e.g. amber14sb_OL15)
            ligand_ff: Ligand force-field name (e.g. gaff2, opls)

        Returns:
            Normalized directory name, format: <protein_ff>-mut_<ligand_ff>
            e.g. amber14sb_ol15-mut_gaff2

        Examples:
            >>> FEPSystemNamer.generate_name("amber14sb_OL15", "gaff2")
            'amber14sb_ol15-mut_gaff2'
            >>> FEPSystemNamer.generate_name("charmm36_jul2022", "cgenff")
            'charmm36_jul2022-mut_cgenff'
        """
        # Normalize the force-field names
        protein = cls._normalize_forcefield_name(protein_ff)
        ligand = cls._normalize_forcefield_name(ligand_ff)

        # Build the standard directory name
        standard_name = f"{protein}-mut_{ligand}"

        return standard_name

    @classmethod
    def validate_name(cls, name: str) -> bool:
        """Validate that a directory name follows the FEP system naming convention.

        Convention: <protein_ff>-mut_<ligand_ff>
        - May contain only lowercase letters, digits, underscores and hyphens
        - Must contain the "-mut_" separator

        Args:
            name: Directory name

        Returns:
            Whether the name conforms to the convention

        Examples:
            >>> FEPSystemNamer.validate_name("amber14sb_ol15-mut_gaff2")
            True
            >>> FEPSystemNamer.validate_name("oplsaa")
            False
            >>> FEPSystemNamer.validate_name("amber14sb_ol15_gaff2")
            False
        """
        # Check the format: <protein_ff>-mut_<ligand_ff>
        pattern = r"^[a-z0-9_]+-mut_[a-z0-9_]+$"
        return bool(re.match(pattern, name))

    @classmethod
    def suggest_name(cls, protein_ff: str, ligand_ff: str) -> str:
        """Suggest a directory name (with hint information).

        Args:
            protein_ff: Protein force-field name
            ligand_ff: Ligand force-field name

        Returns:
            The suggested standard directory name
        """
        standard = cls.generate_name(protein_ff, ligand_ff)
        return standard

    @classmethod
    def parse_name(cls, name: str) -> Optional[Dict[str, str]]:
        """Parse a standard directory name and extract force-field information.

        Args:
            name: Standard directory name

        Returns:
            A dict with protein_ff and ligand_ff, or None if parsing fails

        Examples:
            >>> FEPSystemNamer.parse_name("amber14sb_ol15-mut_gaff2")
            {'protein_ff': 'amber14sb_ol15', 'ligand_ff': 'gaff2'}
        """
        if not cls.validate_name(name):
            return None

        parts = name.split("-mut_")
        if len(parts) == 2:
            return {"protein_ff": parts[0], "ligand_ff": parts[1]}
        return None

    @classmethod
    def _normalize_forcefield_name(cls, ff_name: str) -> str:
        """Normalize a force-field name.

        Conversion rules:
        1. Lowercase
        2. Remove redundant separators ('-' and '_')
        3. Apply the normalization map

        Args:
            ff_name: Original force-field name

        Returns:
            The normalized force-field name
        """
        # Check local mapping case-insensitively for protein force fields
        normalized_lower = ff_name.lower()
        for key, value in cls.FORCEFIELD_NORMALIZATION.items():
            if key.lower() == normalized_lower:
                return value

        # Then prefer the canonical ligand-forcefield registry when applicable.
        try:
            return normalize_forcefield_name(ff_name).replace("-", "_")
        except KeyError:
            pass

        # Normalize: lowercase and strip special characters
        normalized = ff_name.lower()
        # Remove all hyphens and underscores
        normalized = re.sub(r"[-_]+", "", normalized)

        return normalized

    @classmethod
    def get_default_output_dir(cls, protein_ff: str, ligand_ff: str, base_dir: str = ".") -> Path:
        """Get the default output directory path.

        Args:
            protein_ff: Protein force-field name
            ligand_ff: Ligand force-field name
            base_dir: Base directory (defaults to the current directory)

        Returns:
            The full output directory path
        """
        standard_name = cls.generate_name(protein_ff, ligand_ff)
        return Path(base_dir) / standard_name


def generate_fep_system_name(protein_ff: str, ligand_ff: str) -> str:
    """Convenience function to generate a FEP system directory name.

    Args:
        protein_ff: Protein force-field name
        ligand_ff: Ligand force-field name

    Returns:
        The normalized directory name

    Examples:
        >>> generate_fep_system_name("amber14sb_OL15", "gaff2")
        'amber14sb_ol15-mut_gaff2'
    """
    return FEPSystemNamer.generate_name(protein_ff, ligand_ff)


def validate_fep_system_name(name: str) -> bool:
    """Convenience function to validate a FEP system directory name.

    Args:
        name: Directory name

    Returns:
        Whether the name conforms to the convention
    """
    return FEPSystemNamer.validate_name(name)


def parse_fep_system_name(name: str) -> Optional[Dict[str, str]]:
    """Convenience function to parse a FEP system directory name.

    Args:
        name: Standard directory name

    Returns:
        A dict with force-field information, or None if parsing fails
    """
    return FEPSystemNamer.parse_name(name)
