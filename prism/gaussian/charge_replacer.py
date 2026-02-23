#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RESP Charge Replacer for GROMACS ITP Files

This module provides functionality to replace AM1-BCC charges in GROMACS ITP files
with RESP charges from Gaussian calculations.
"""

import shutil
from typing import Dict, List, Tuple, Optional

# Import color utilities
try:
    from ..utils.colors import print_success, print_info, print_warning, print_error
except ImportError:
    try:
        from prism.utils.colors import print_success, print_info, print_warning, print_error
    except ImportError:

        def print_success(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}[OK] {x}")

        def print_info(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}[INFO] {x}")

        def print_warning(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}[WARN] {x}")

        def print_error(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}[ERROR] {x}")


class RESPChargeReplacer:
    """
    Replace charges in GROMACS ITP files with RESP charges from MOL2 files.

    Supports two matching strategies:
    1. Name matching: Match atoms by atom name
    2. Sequential matching: Match atoms by order (fallback)
    """

    def __init__(self, resp_mol2_path: str, verbose: bool = True):
        """
        Initialize the RESP charge replacer.

        Parameters
        ----------
        resp_mol2_path : str
            Path to RESP MOL2 file with calculated charges
        verbose : bool
            Whether to print progress information
        """
        self.resp_mol2_path = resp_mol2_path
        self.verbose = verbose
        self.charges = {}  # atom_name -> charge
        self.charges_list = []  # [(atom_name, charge), ...] in order

        # Load charges on initialization
        self.load_resp_charges()

    def load_resp_charges(self) -> Dict[str, float]:
        """
        Load RESP charges from MOL2 file.

        Extracts atom names and charges from the @<TRIPOS>ATOM section.

        Returns
        -------
        Dict[str, float]
            Dictionary mapping atom names to charges
        """
        if self.verbose:
            print_info(f"Loading RESP charges from: {self.resp_mol2_path}")

        self.charges = {}
        self.charges_list = []
        in_atom_section = False

        with open(self.resp_mol2_path, "r") as f:
            for line in f:
                line_stripped = line.strip()

                if line_stripped.startswith("@<TRIPOS>ATOM"):
                    in_atom_section = True
                    continue
                elif line_stripped.startswith("@<TRIPOS>"):
                    in_atom_section = False
                    continue

                if in_atom_section and line_stripped:
                    parts = line.split()
                    if len(parts) >= 9:
                        # MOL2 ATOM format: id name x y z type res_id res_name charge
                        atom_name = parts[1]
                        try:
                            charge = float(parts[8])
                            self.charges[atom_name] = charge
                            self.charges_list.append((atom_name, charge))
                        except (ValueError, IndexError):
                            continue

        if self.verbose:
            print_info(f"  Loaded {len(self.charges)} atom charges")
            total_charge = sum(self.charges.values())
            print_info(f"  Total charge: {total_charge:.4f} (rounded: {round(total_charge)})")

        return self.charges

    def _parse_itp_atoms(self, itp_path: str) -> Tuple[List[dict], List[str], int, int]:
        """
        Parse atoms section from ITP file.

        Returns
        -------
        Tuple[List[dict], List[str], int, int]
            - List of atom dictionaries with parsed fields
            - List of original lines
            - Start line index of atoms section
            - End line index of atoms section
        """
        atoms = []
        lines = []
        in_atoms_section = False
        atoms_start = -1
        atoms_end = -1

        with open(itp_path, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            stripped = line.strip()

            # Detect section headers
            if stripped.startswith("[") and "]" in stripped:
                section_name = stripped.split("]")[0].split("[")[1].strip().lower()
                if section_name == "atoms":
                    in_atoms_section = True
                    atoms_start = i
                    continue
                elif in_atoms_section:
                    in_atoms_section = False
                    atoms_end = i
                    break

            # Parse atom lines in [ atoms ] section
            if in_atoms_section and stripped and not stripped.startswith(";"):
                # ITP atoms format: nr type resnr residue atom cgnr charge mass
                # Typical: 1  c3  1  LIG  C1  1  -0.108000  12.01
                parts = line.split()
                if len(parts) >= 7:
                    try:
                        atom = {
                            "line_idx": i,
                            "nr": int(parts[0]),
                            "type": parts[1],
                            "resnr": int(parts[2]),
                            "residue": parts[3],
                            "atom": parts[4],
                            "cgnr": int(parts[5]),
                            "charge": float(parts[6]),
                            "mass": float(parts[7]) if len(parts) > 7 else None,
                            "original_parts": parts,
                        }
                        atoms.append(atom)
                    except (ValueError, IndexError):
                        continue

        if atoms_end == -1:
            atoms_end = len(lines)

        return atoms, lines, atoms_start, atoms_end

    def _try_name_matching(self, itp_atoms: List[dict]) -> Optional[Dict[int, float]]:
        """
        Try to match charges by atom name.

        Returns
        -------
        Optional[Dict[int, float]]
            Dictionary mapping line index to new charge, or None if matching fails
        """
        matches = {}
        unmatched = []

        for atom in itp_atoms:
            atom_name = atom["atom"]

            if atom_name in self.charges:
                matches[atom["line_idx"]] = self.charges[atom_name]
            else:
                # Try variations
                found = False
                for variation in [atom_name.upper(), atom_name.lower(), atom_name.capitalize()]:
                    if variation in self.charges:
                        matches[atom["line_idx"]] = self.charges[variation]
                        found = True
                        break

                if not found:
                    unmatched.append(atom_name)

        # Check if name matching was successful
        if len(unmatched) == 0:
            return matches
        elif len(matches) > len(unmatched):
            # Partial match - warn but continue
            if self.verbose:
                print_warning(f"  Name matching: {len(matches)}/{len(itp_atoms)} atoms matched")
                print_warning(f"  Unmatched atoms: {unmatched[:5]}...")
            return None
        else:
            return None

    def _apply_sequential_matching(self, itp_atoms: List[dict]) -> Dict[int, float]:
        """
        Match charges sequentially by order.

        Returns
        -------
        Dict[int, float]
            Dictionary mapping line index to new charge
        """
        matches = {}

        if len(itp_atoms) != len(self.charges_list):
            print_warning(f"  Atom count mismatch: ITP has {len(itp_atoms)}, " f"RESP has {len(self.charges_list)}")

        # Match by order
        for i, atom in enumerate(itp_atoms):
            if i < len(self.charges_list):
                matches[atom["line_idx"]] = self.charges_list[i][1]

        return matches

    def replace_itp_charges(self, itp_path: str, output_path: Optional[str] = None, backup: bool = True) -> str:
        """
        Replace charges in ITP file with RESP charges.

        Parameters
        ----------
        itp_path : str
            Path to input ITP file
        output_path : str, optional
            Path for output file. If None, overwrites input file.
        backup : bool
            Whether to create backup of original file (default: True)

        Returns
        -------
        str
            Path to output file
        """
        if self.verbose:
            print_info(f"Replacing charges in: {itp_path}")

        # Set output path
        if output_path is None:
            output_path = itp_path

        # Create backup if requested and overwriting
        if backup and output_path == itp_path:
            backup_path = itp_path + ".bak"
            shutil.copy2(itp_path, backup_path)
            if self.verbose:
                print_info(f"  Backup created: {backup_path}")

        # Parse ITP
        itp_atoms, lines, atoms_start, atoms_end = self._parse_itp_atoms(itp_path)

        if not itp_atoms:
            raise ValueError(f"No atoms found in ITP file: {itp_path}")

        if self.verbose:
            print_info(f"  Found {len(itp_atoms)} atoms in ITP")

        # Try name matching first
        charge_map = self._try_name_matching(itp_atoms)

        if charge_map is None:
            # Fall back to sequential matching
            if self.verbose:
                print_info("  Using sequential matching (fallback)")
            charge_map = self._apply_sequential_matching(itp_atoms)

        # Apply charge replacements
        charges_replaced = 0
        for atom in itp_atoms:
            line_idx = atom["line_idx"]
            if line_idx in charge_map:
                new_charge = charge_map[line_idx]
                # Reconstruct the line with new charge
                parts = atom["original_parts"].copy()
                parts[6] = f"{new_charge:.6f}"

                # Preserve line format as much as possible
                new_line = self._format_atom_line(parts, lines[line_idx])
                lines[line_idx] = new_line

                charges_replaced += 1

        # Write output
        with open(output_path, "w") as f:
            f.writelines(lines)

        if self.verbose:
            print_success(f"  Replaced {charges_replaced} charges")
            print_success(f"  Output: {output_path}")

        return output_path

    def _format_atom_line(self, parts: List[str], original_line: str) -> str:
        """
        Format atom line preserving original spacing style.

        Parameters
        ----------
        parts : List[str]
            Parts of the atom line
        original_line : str
            Original line for reference

        Returns
        -------
        str
            Formatted atom line
        """
        # Try to preserve original column alignment
        line_ending = "\r\n" if original_line.endswith("\r\n") else "\n"
        # Standard ITP format with fixed column widths
        if len(parts) >= 8:
            return (
                f"{int(parts[0]):>6}  {parts[1]:<4}  {int(parts[2]):>4}  "
                f"{parts[3]:<4}  {parts[4]:<4}  {int(parts[5]):>4}  "
                f"{float(parts[6]):>10.6f}  {float(parts[7]):>8.4f}{line_ending}"
            )
        elif len(parts) >= 7:
            return (
                f"{int(parts[0]):>6}  {parts[1]:<4}  {int(parts[2]):>4}  "
                f"{parts[3]:<4}  {parts[4]:<4}  {int(parts[5]):>4}  "
                f"{float(parts[6]):>10.6f}{line_ending}"
            )
        else:
            # Fallback: simple space-separated
            return "  ".join(parts) + line_ending


def replace_itp_charges(
    itp_path: str, resp_mol2_path: str, output_path: Optional[str] = None, backup: bool = True, verbose: bool = True
) -> str:
    """
    Replace charges in ITP file with RESP charges from MOL2 file.

    Convenience function that wraps RESPChargeReplacer.

    Parameters
    ----------
    itp_path : str
        Path to input ITP file
    resp_mol2_path : str
        Path to RESP MOL2 file with calculated charges
    output_path : str, optional
        Path for output file. If None, overwrites input file.
    backup : bool
        Whether to create backup of original file (default: True)
    verbose : bool
        Whether to print progress information

    Returns
    -------
    str
        Path to output file

    Examples
    --------
    >>> from prism.gaussian import replace_itp_charges
    >>> replace_itp_charges("LIG.itp", "ligand_resp.mol2")
    """
    replacer = RESPChargeReplacer(resp_mol2_path, verbose=verbose)
    return replacer.replace_itp_charges(itp_path, output_path, backup)
