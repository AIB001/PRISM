#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MOL2 to XYZ Coordinate Converter for Gaussian Calculations

This module provides functionality to convert MOL2 files to XYZ format
for use with Gaussian quantum chemistry calculations.
"""

import os
from typing import List, Tuple, Optional

# Import color utilities
try:
    from ..utils.colors import print_success, print_info, print_warning, print_error
except ImportError:
    try:
        from prism.utils.colors import print_success, print_info, print_warning, print_error
    except ImportError:
        def print_success(x, **kwargs): print(f"[OK] {x}")
        def print_info(x, **kwargs): print(f"[INFO] {x}")
        def print_warning(x, **kwargs): print(f"[WARN] {x}")
        def print_error(x, **kwargs): print(f"[ERROR] {x}")


class CoordinateConverter:
    """
    Converter for molecular coordinate formats.

    Provides conversion between MOL2 and XYZ formats for Gaussian calculations.
    """

    # Element symbol mapping from atomic number
    ATOMIC_SYMBOLS = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
        9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
        23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
        30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
        37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
        44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
        51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
        58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
        65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
        72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
        79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
        86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U'
    }

    # Reverse mapping: symbol to atomic number
    SYMBOL_TO_NUMBER = {v: k for k, v in ATOMIC_SYMBOLS.items()}

    def __init__(self, verbose: bool = True):
        """
        Initialize the coordinate converter.

        Parameters
        ----------
        verbose : bool
            Whether to print progress information
        """
        self.verbose = verbose

    def parse_mol2(self, mol2_file: str) -> Tuple[List[Tuple[str, float, float, float]], List[float], str]:
        """
        Parse a MOL2 file to extract atomic coordinates and charges.

        Parameters
        ----------
        mol2_file : str
            Path to the MOL2 file

        Returns
        -------
        Tuple[List[Tuple], List[float], str]
            - List of (element, x, y, z) tuples
            - List of partial charges
            - Molecule name
        """
        atoms = []
        charges = []
        mol_name = os.path.splitext(os.path.basename(mol2_file))[0]
        in_atom_section = False

        with open(mol2_file, 'r') as f:
            for line in f:
                line_stripped = line.strip()

                # Check for MOLECULE section to get name
                if line_stripped.startswith("@<TRIPOS>MOLECULE"):
                    mol_name = next(f).strip()
                    continue

                # Check for ATOM section
                if line_stripped.startswith("@<TRIPOS>ATOM"):
                    in_atom_section = True
                    continue
                elif line_stripped.startswith("@<TRIPOS>"):
                    in_atom_section = False
                    continue

                # Parse atom lines
                if in_atom_section and line_stripped:
                    parts = line.split()
                    if len(parts) >= 6:
                        # MOL2 format: atom_id atom_name x y z atom_type [res_id res_name [charge]]
                        atom_name = parts[1]
                        x = float(parts[2])
                        y = float(parts[3])
                        z = float(parts[4])
                        atom_type = parts[5]

                        # Extract element from atom type (e.g., "C.ar" -> "C")
                        element = self._extract_element(atom_name, atom_type)

                        atoms.append((element, x, y, z))

                        # Extract charge if present
                        if len(parts) >= 9:
                            try:
                                charge = float(parts[8])
                                charges.append(charge)
                            except ValueError:
                                charges.append(0.0)
                        else:
                            charges.append(0.0)

        if not atoms:
            raise ValueError(f"No atoms found in MOL2 file: {mol2_file}")

        return atoms, charges, mol_name

    def _normalize_element_symbol(self, element: str) -> str:
        """
        Normalize element symbol to standard format (first letter uppercase, rest lowercase).

        Examples: "CL" -> "Cl", "BR" -> "Br", "c" -> "C", "na" -> "Na"

        Parameters
        ----------
        element : str
            Raw element symbol

        Returns
        -------
        str
            Normalized element symbol
        """
        element = element.strip()
        if len(element) == 0:
            return element
        elif len(element) == 1:
            return element.upper()
        else:
            return element[0].upper() + element[1:].lower()

    def _extract_element(self, atom_name: str, atom_type: str) -> str:
        """
        Extract element symbol from atom name or type.

        Parameters
        ----------
        atom_name : str
            Atom name from MOL2 (e.g., "C1", "N2", "CL01")
        atom_type : str
            Atom type from MOL2 (e.g., "C.ar", "N.am", "Cl", "Br")

        Returns
        -------
        str
            Element symbol (e.g., "C", "N", "Cl", "Br")
        """
        # First try to extract from atom type (before the dot)
        if '.' in atom_type:
            element = atom_type.split('.')[0]
        else:
            element = atom_type

        # Normalize element symbol (e.g., "CL" -> "Cl", "BR" -> "Br")
        element = self._normalize_element_symbol(element)

        # Check if valid element
        if element in self.SYMBOL_TO_NUMBER:
            return element

        # Try to extract from atom name (first 1-2 characters)
        # Try 2-character elements first (Cl, Br, Na, etc.)
        for length in [2, 1]:
            if len(atom_name) >= length:
                # Extract alphabetic characters only
                alpha_chars = ''.join(c for c in atom_name[:length+1] if c.isalpha())
                if len(alpha_chars) >= length:
                    candidate = self._normalize_element_symbol(alpha_chars[:length])
                    if candidate in self.SYMBOL_TO_NUMBER:
                        return candidate

        # Last resort: use first letter
        first_letter = ''.join(c for c in atom_name if c.isalpha())[:1].upper()
        if first_letter in self.SYMBOL_TO_NUMBER:
            return first_letter

        # Default to carbon if all else fails
        print_warning(f"Could not determine element for atom '{atom_name}' type '{atom_type}', defaulting to C")
        return 'C'

    def write_xyz(self, atoms: List[Tuple[str, float, float, float]],
                  xyz_file: str, comment: str = "") -> str:
        """
        Write atoms to XYZ format file.

        Parameters
        ----------
        atoms : List[Tuple]
            List of (element, x, y, z) tuples
        xyz_file : str
            Path to output XYZ file
        comment : str
            Comment line for XYZ file

        Returns
        -------
        str
            Path to the written XYZ file
        """
        with open(xyz_file, 'w') as f:
            # Number of atoms
            f.write(f"{len(atoms)}\n")

            # Comment line
            f.write(f"{comment}\n")

            # Atom coordinates
            for element, x, y, z in atoms:
                f.write(f"{element:2s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")

        return xyz_file

    def mol2_to_xyz(self, mol2_file: str, xyz_file: Optional[str] = None) -> str:
        """
        Convert MOL2 file to XYZ format.

        Parameters
        ----------
        mol2_file : str
            Path to input MOL2 file
        xyz_file : str, optional
            Path to output XYZ file. If not provided, uses same name with .xyz extension

        Returns
        -------
        str
            Path to the output XYZ file
        """
        if self.verbose:
            print_info(f"Converting MOL2 to XYZ: {mol2_file}")

        # Determine output file path
        if xyz_file is None:
            xyz_file = os.path.splitext(mol2_file)[0] + ".xyz"

        # Parse MOL2
        atoms, charges, mol_name = self.parse_mol2(mol2_file)

        if self.verbose:
            print_info(f"  Found {len(atoms)} atoms in molecule '{mol_name}'")

        # Write XYZ
        comment = f"Converted from {os.path.basename(mol2_file)} | {mol_name}"
        self.write_xyz(atoms, xyz_file, comment)

        if self.verbose:
            print_success(f"  XYZ file written: {xyz_file}")

        return xyz_file

    def get_atoms_from_mol2(self, mol2_file: str) -> List[Tuple[str, float, float, float]]:
        """
        Get atom list from MOL2 file.

        Parameters
        ----------
        mol2_file : str
            Path to MOL2 file

        Returns
        -------
        List[Tuple]
            List of (element, x, y, z) tuples
        """
        atoms, _, _ = self.parse_mol2(mol2_file)
        return atoms

    def get_charges_from_mol2(self, mol2_file: str) -> Tuple[List[float], int]:
        """
        Get partial charges and net charge from MOL2 file.

        Parameters
        ----------
        mol2_file : str
            Path to MOL2 file

        Returns
        -------
        Tuple[List[float], int]
            - List of partial charges
            - Net charge (rounded to integer)
        """
        _, charges, _ = self.parse_mol2(mol2_file)
        net_charge = int(round(sum(charges)))
        return charges, net_charge


def mol2_to_xyz(mol2_file: str, xyz_file: Optional[str] = None, verbose: bool = True) -> str:
    """
    Convert MOL2 file to XYZ format.

    Convenience function that wraps CoordinateConverter.mol2_to_xyz().

    Parameters
    ----------
    mol2_file : str
        Path to input MOL2 file
    xyz_file : str, optional
        Path to output XYZ file. If not provided, uses same name with .xyz extension
    verbose : bool
        Whether to print progress information

    Returns
    -------
    str
        Path to the output XYZ file

    Examples
    --------
    >>> from prism.gaussian import mol2_to_xyz
    >>> xyz_path = mol2_to_xyz("ligand.mol2")
    >>> print(f"XYZ file created: {xyz_path}")
    """
    converter = CoordinateConverter(verbose=verbose)
    return converter.mol2_to_xyz(mol2_file, xyz_file)
