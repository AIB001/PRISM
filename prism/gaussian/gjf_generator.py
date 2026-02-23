#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gaussian Input File Generator for RESP Charge Calculations

This module generates Gaussian input files (.gjf) for:
1. Geometry optimization
2. Electrostatic potential (ESP) calculations for RESP fitting
"""

import os
from typing import List, Tuple

from .converter import CoordinateConverter
from .utils import get_charge_multiplicity, extract_optimized_coords_from_log

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


class GaussianInputGenerator:
    """
    Generator for Gaussian input files for RESP charge calculations.

    Supports both HF/6-31G* and B3LYP/6-31G* methods, which are commonly
    used for RESP charge fitting.
    """

    # Standard basis sets and methods for RESP
    METHODS = {"hf": "HF", "dft": "B3LYP"}

    BASIS_SET = "6-31G*"

    def __init__(self, method: str = "dft", nproc: int = 16, mem: str = "4GB"):
        """
        Initialize the Gaussian input file generator.

        Parameters
        ----------
        method : str
            Quantum chemistry method: 'hf' for HF/6-31G* or 'dft' for B3LYP/6-31G*
        nproc : int
            Number of processors to use (default: 16)
        mem : str
            Memory allocation (default: '4GB')
        """
        self.method = self.METHODS.get(method.lower(), "B3LYP")
        self.method_key = method.lower()
        self.basis_set = self.BASIS_SET
        self.nproc = nproc
        self.mem = mem

    def _write_link0(self, f, chk_file: str = None):
        """Write Link0 section (resource specifications)."""
        f.write(f"%nproc={self.nproc}\n")
        f.write(f"%mem={self.mem}\n")
        if chk_file:
            f.write(f"%chk={chk_file}\n")

    def _format_coords(self, atoms: List[Tuple[str, float, float, float]]) -> str:
        """Format atomic coordinates for Gaussian input."""
        lines = []
        for element, x, y, z in atoms:
            lines.append(f" {element:2s}    {x:14.8f}  {y:14.8f}  {z:14.8f}")
        return "\n".join(lines)

    def generate_opt_gjf(self, xyz_file: str, mol2_file: str, output_file: str) -> str:
        """
        Generate Gaussian input file for geometry optimization.

        Parameters
        ----------
        xyz_file : str
            Path to XYZ file with initial coordinates
        mol2_file : str
            Path to MOL2 file (for charge/multiplicity)
        output_file : str
            Path for output GJF file

        Returns
        -------
        str
            Path to generated GJF file
        """
        print_info(f"Generating optimization input: {output_file}")

        # Get charge and multiplicity
        charge, multiplicity = get_charge_multiplicity(mol2_file)
        print_info(f"  Charge: {charge}, Multiplicity: {multiplicity}")

        # Read coordinates from XYZ file
        atoms = self._read_xyz(xyz_file)

        # Generate checkpoint filename
        chk_file = os.path.splitext(os.path.basename(output_file))[0] + ".chk"

        with open(output_file, "w") as f:
            # Link0 section
            self._write_link0(f, chk_file)

            # Route section
            f.write(f"#p {self.method}/{self.basis_set} Opt\n")
            f.write("\n")

            # Title
            f.write(f"Geometry optimization for RESP - {self.method}/{self.basis_set}\n")
            f.write("\n")

            # Charge and multiplicity
            f.write(f"{charge} {multiplicity}\n")

            # Coordinates
            f.write(self._format_coords(atoms))
            f.write("\n\n")

        print_success(f"  Generated: {output_file}")
        return output_file

    def generate_esp_gjf(
        self, atoms: List[Tuple[str, float, float, float]], charge: int, multiplicity: int, output_file: str
    ) -> str:
        """
        Generate Gaussian input file for ESP calculation.

        Parameters
        ----------
        atoms : List[Tuple]
            List of (element, x, y, z) tuples
        charge : int
            Net molecular charge
        multiplicity : int
            Spin multiplicity
        output_file : str
            Path for output GJF file

        Returns
        -------
        str
            Path to generated GJF file
        """
        print_info(f"Generating ESP input: {output_file}")

        # Generate checkpoint filename
        chk_file = os.path.splitext(os.path.basename(output_file))[0] + ".chk"

        with open(output_file, "w") as f:
            # Link0 section
            self._write_link0(f, chk_file)

            # Route section - ESP calculation with Pop=MK for RESP fitting
            # IOp(6/33=2) requests ESP output for RESP
            # IOp(6/42=6) sets the density of points on layers
            f.write(f"#p {self.method}/{self.basis_set} Pop=MK IOp(6/33=2,6/42=6)\n")
            f.write("\n")

            # Title
            f.write(f"ESP calculation for RESP fitting - {self.method}/{self.basis_set}\n")
            f.write("\n")

            # Charge and multiplicity
            f.write(f"{charge} {multiplicity}\n")

            # Coordinates
            f.write(self._format_coords(atoms))
            f.write("\n\n")

        print_success(f"  Generated: {output_file}")
        return output_file

    def generate_esp_from_opt_log(self, log_file: str, output_file: str) -> str:
        """
        Generate ESP input file using optimized coordinates from a Gaussian log file.

        Parameters
        ----------
        log_file : str
            Path to optimization log file
        output_file : str
            Path for output GJF file

        Returns
        -------
        str
            Path to generated GJF file
        """
        print_info(f"Extracting optimized coordinates from: {log_file}")

        # Extract optimized coordinates
        atoms = extract_optimized_coords_from_log(log_file)
        if not atoms:
            raise ValueError(f"Could not extract coordinates from {log_file}")

        # Extract charge and multiplicity from log file
        charge, multiplicity = self._get_charge_mult_from_log(log_file)

        return self.generate_esp_gjf(atoms, charge, multiplicity, output_file)

    def generate_esp_direct(self, xyz_file: str, mol2_file: str, output_file: str) -> str:
        """
        Generate ESP input file directly from XYZ/MOL2 files (no optimization).

        Parameters
        ----------
        xyz_file : str
            Path to XYZ file with coordinates
        mol2_file : str
            Path to MOL2 file (for charge/multiplicity)
        output_file : str
            Path for output GJF file

        Returns
        -------
        str
            Path to generated GJF file
        """
        print_info(f"Generating direct ESP input: {output_file}")

        # Get charge and multiplicity
        charge, multiplicity = get_charge_multiplicity(mol2_file)

        # Read coordinates
        atoms = self._read_xyz(xyz_file)

        return self.generate_esp_gjf(atoms, charge, multiplicity, output_file)

    def _read_xyz(self, xyz_file: str) -> List[Tuple[str, float, float, float]]:
        """Read coordinates from XYZ file."""
        atoms = []
        with open(xyz_file, "r") as f:
            lines = f.readlines()

        # Skip first two lines (atom count and comment)
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0]
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                atoms.append((element, x, y, z))

        return atoms

    def _get_charge_mult_from_log(self, log_file: str) -> Tuple[int, int]:
        """Extract charge and multiplicity from Gaussian log file."""
        charge = 0
        multiplicity = 1

        with open(log_file, "r") as f:
            for line in f:
                if "Charge =" in line and "Multiplicity =" in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "Charge":
                            charge = int(parts[i + 2])
                        if part == "Multiplicity":
                            multiplicity = int(parts[i + 2])
                    break

        return charge, multiplicity

    def generate_all_inputs(self, mol2_file: str, output_dir: str, do_optimization: bool = False) -> dict:
        """
        Generate all necessary Gaussian input files.

        Parameters
        ----------
        mol2_file : str
            Path to input MOL2 file
        output_dir : str
            Output directory for generated files
        do_optimization : bool
            Whether to generate optimization input (default: False)

        Returns
        -------
        dict
            Dictionary with paths to generated files
        """
        os.makedirs(output_dir, exist_ok=True)

        files = {}

        # Convert MOL2 to XYZ
        converter = CoordinateConverter(verbose=False)
        xyz_file = os.path.join(output_dir, "ligand.xyz")
        converter.mol2_to_xyz(mol2_file, xyz_file)
        files["xyz"] = xyz_file

        if do_optimization:
            # Generate optimization input
            opt_gjf = os.path.join(output_dir, f"lig_opt_{self.method_key}.gjf")
            self.generate_opt_gjf(xyz_file, mol2_file, opt_gjf)
            files["opt_gjf"] = opt_gjf

            # ESP will be generated from optimized coordinates later
            esp_gjf = os.path.join(output_dir, f"lig_esp_{self.method_key}.gjf")
            files["esp_gjf"] = esp_gjf  # Placeholder path
        else:
            # Generate ESP input directly
            esp_gjf = os.path.join(output_dir, f"lig_esp_{self.method_key}.gjf")
            self.generate_esp_direct(xyz_file, mol2_file, esp_gjf)
            files["esp_gjf"] = esp_gjf

        return files
