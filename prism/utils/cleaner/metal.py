#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Metal detection utilities for protein cleaning.
"""


class MetalMixin:
    """Mixin for metal detection and processing."""

    def _is_metal_or_ion(self, residue_name: str, atom_name: str) -> bool:
        """
        Check if a residue/atom is a metal or ion.

        Parameters
        ----------
        residue_name : str
            Residue name from PDB
        atom_name : str
            Atom name from PDB

        Returns
        -------
        bool
            True if this is a metal or ion
        """
        # Check against known lists
        if residue_name in self.structural_metals or residue_name in self.non_structural_ions:
            return True

        # Also check atom name (sometimes metal is in atom name)
        if atom_name in self.structural_metals or atom_name in self.non_structural_ions:
            return True

        # Common metal ion patterns
        metal_patterns = ["ZN", "MG", "CA", "FE", "CU", "MN", "CO", "NI", "NA", "K", "CL"]
        for pattern in metal_patterns:
            if pattern in residue_name or pattern in atom_name:
                return True

        return False

    def _should_keep_metal_smart(self, residue_name: str, atom_name: str) -> bool:
        """
        Determine if a metal should be kept in 'smart' mode.

        Parameters
        ----------
        residue_name : str
            Residue name from PDB
        atom_name : str
            Atom name from PDB

        Returns
        -------
        bool
            True if this metal should be kept
        """
        # Check if it's a structural metal
        if residue_name in self.structural_metals or atom_name in self.structural_metals:
            return True

        # Check if it's explicitly a non-structural ion
        if residue_name in self.non_structural_ions or atom_name in self.non_structural_ions:
            return False

        # Unknown metal - be conservative and keep it
        if self.verbose:
            print(f"  Warning: Unknown metal/ion {residue_name}/{atom_name}, keeping by default")
        return True
