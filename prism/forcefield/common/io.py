#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unified coordinate file I/O for PRISM forcefield module.

This module provides standardized functions for reading and writing
molecular coordinate files in various formats (GRO, MOL2, PDB).
It consolidates duplicate implementations from multiple forcefield
generators (GAFF, OpenFF, CGenFF, CHARMM-GUI, RTF, OPLS-AA).

All functions return/accept coordinates in Angstroms for internal processing,
with automatic conversion to/from nanometers for GRO files.
"""

from pathlib import Path
from typing import List, Tuple, Optional, Union
import numpy as np

from .units import angstrom_xyz_to_nm, nm_xyz_to_angstrom


def read_gro_coordinates(
    gro_file: Union[str, Path], atom_count: Optional[int] = None
) -> List[Tuple[int, str, np.ndarray]]:
    """Read coordinates from GROMACS GRO file.

    GRO files use nanometers as units. This function converts to Angstroms
    for consistency with other file formats.

    Parameters
    ----------
    gro_file : Union[str, Path]
        Path to GRO file
    atom_count : Optional[int]
        Expected number of atoms. If None, reads all atoms found.

    Returns
    -------
    List[Tuple[int, str, np.ndarray]]
        List of (atom_index, atom_name, coordinates) tuples.
        Coordinates are numpy arrays in Angstroms.

    Raises
    ------
    FileNotFoundError
        If GRO file doesn't exist
    ValueError
        If GRO file is malformed

    Examples
    --------
    >>> coords = read_gro_coordinates("ligand.gro")
    >>> for idx, name, xyz in coords:
    ...     print(f"Atom {idx} ({name}): {xyz}")
    """
    gro_path = Path(gro_file)
    if not gro_path.exists():
        raise FileNotFoundError(f"GRO file not found: {gro_path}")

    coords: List[Tuple[int, str, np.ndarray]] = []

    with gro_path.open("r", encoding="utf-8") as handle:
        # Skip title line
        handle.readline()

        # Read atom count
        count_line = handle.readline()
        if not count_line:
            raise ValueError("Malformed GRO file: missing atom count")

        total_atoms = int(count_line.strip())
        if atom_count is not None and atom_count > total_atoms:
            raise ValueError(f"Expected {atom_count} atoms but GRO only has {total_atoms}")

        # Read coordinates
        for i in range(total_atoms):
            line = handle.readline()
            if not line:
                raise ValueError(f"Unexpected end of GRO file at line {i+2}")

            # GRO format: fixed column positions
            # Columns: resid(5) resname(5) atomname(5) atomnr(5) x(8) y(8) z(8) vx(8) vy(8) vz(8)
            # We only need: atomnr, atomname, x, y, z
            try:
                atom_name = line[10:15].strip()
                atom_idx = int(line[15:20].strip())
                x_nm = float(line[20:28].strip())
                y_nm = float(line[28:36].strip())
                z_nm = float(line[36:44].strip())
            except (ValueError, IndexError) as exc:
                raise ValueError(f"Malformed GRO line {i+2}: {line.strip()}") from exc

            # Convert from nm to Angstroms
            xyz_nm = np.array([x_nm, y_nm, z_nm])
            xyz_angstrom = nm_xyz_to_angstrom(xyz_nm)

            coords.append((atom_idx, atom_name, xyz_angstrom))

            if atom_count is not None and len(coords) >= atom_count:
                break

    if atom_count is not None and len(coords) < atom_count:
        raise ValueError(f"Expected {atom_count} atoms but only found {len(coords)}")

    return coords


def read_mol2_coordinates(
    mol2_file: Union[str, Path], atom_count: Optional[int] = None
) -> List[Tuple[int, str, np.ndarray]]:
    """Read coordinates from MOL2 file.

    MOL2 files use Angstroms as units (no conversion needed).

    Parameters
    ----------
    mol2_file : Union[str, Path]
        Path to MOL2 file
    atom_count : Optional[int]
        Expected number of atoms. If None, reads all atoms found.

    Returns
    -------
    List[Tuple[int, str, np.ndarray]]
        List of (atom_index, atom_name, coordinates) tuples.
        Coordinates are numpy arrays in Angstroms.

    Raises
    ------
    FileNotFoundError
        If MOL2 file doesn't exist
    ValueError
        If MOL2 file has no @<TRIPOS>ATOM section

    Examples
    --------
    >>> coords = read_mol2_coordinates("ligand.mol2")
    >>> for idx, name, xyz in coords:
    ...     print(f"Atom {idx} ({name}): {xyz}")
    """
    mol2_path = Path(mol2_file)
    if not mol2_path.exists():
        raise FileNotFoundError(f"MOL2 file not found: {mol2_path}")

    coords: List[Tuple[int, str, np.ndarray]] = []

    with mol2_path.open("r", encoding="utf-8") as handle:
        in_atoms = False
        for line in handle:
            line = line.strip()
            if not line:
                continue

            # Look for ATOM section
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue

            # End of ATOM section
            if line.startswith("@<TRIPOS>") and in_atoms:
                break

            # Skip lines before ATOM section
            if not in_atoms:
                continue

            # Parse atom line
            parts = line.split()
            if len(parts) < 6:
                continue

            try:
                atom_idx = int(parts[0])
                atom_name = parts[1]
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
            except (ValueError, IndexError) as exc:
                raise ValueError(f"Malformed MOL2 atom line: {line}") from exc

            xyz = np.array([x, y, z])
            coords.append((atom_idx, atom_name, xyz))

            if atom_count is not None and len(coords) >= atom_count:
                break

    if not coords:
        raise ValueError("MOL2 file has no @<TRIPOS>ATOM section or no atoms found")

    if atom_count is not None and len(coords) < atom_count:
        raise ValueError(f"Expected {atom_count} atoms but only found {len(coords)}")

    return coords


def read_pdb_coordinates(
    pdb_file: Union[str, Path], atom_count: Optional[int] = None
) -> List[Tuple[int, str, np.ndarray]]:
    """Read coordinates from PDB file.

    PDB files use Angstroms as units (no conversion needed).

    Parameters
    ----------
    pdb_file : Union[str, Path]
        Path to PDB file
    atom_count : Optional[int]
        Expected number of atoms. If None, reads all atoms found.

    Returns
    -------
    List[Tuple[int, str, np.ndarray]]
        List of (atom_index, atom_name, coordinates) tuples.
        Coordinates are numpy arrays in Angstroms.

    Raises
    ------
    FileNotFoundError
        If PDB file doesn't exist
    ValueError
        If PDB file is malformed

    Examples
    --------
    >>> coords = read_pdb_coordinates("ligand.pdb")
    >>> for idx, name, xyz in coords:
    ...     print(f"Atom {idx} ({name}): {xyz}")
    """
    pdb_path = Path(pdb_file)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    coords: List[Tuple[int, str, np.ndarray]] = []

    with pdb_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            # Only read ATOM or HETATM records
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            # Parse PDB format (fixed columns)
            # Columns: 1-6=record, 7-11=serial, 13-16=atom, 17=altloc, 18-20=resname,
            #          22=chain, 23-26=resseq, 27=icode, 31-38=x, 39-46=y, 47-54=z
            try:
                serial = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except (ValueError, IndexError) as exc:
                raise ValueError(f"Malformed PDB atom line: {line}") from exc

            xyz = np.array([x, y, z])
            coords.append((serial, atom_name, xyz))

            if atom_count is not None and len(coords) >= atom_count:
                break

    if not coords:
        raise ValueError("PDB file has no ATOM/HETATM records")

    if atom_count is not None and len(coords) < atom_count:
        raise ValueError(f"Expected {atom_count} atoms but only found {len(coords)}")

    return coords


def write_gro_coordinates(
    output_file: Union[str, Path],
    coordinates: List[Tuple[int, str, np.ndarray]],
    residue_name: str = "LIG",
    title: str = "Generated by PRISM",
) -> None:
    """Write coordinates to GROMACS GRO file.

    GRO files use nanometers as units. This function converts from Angstroms.

    Parameters
    ----------
    output_file : Union[str, Path]
        Path to output GRO file
    coordinates : List[Tuple[int, str, np.ndarray]]
        List of (atom_index, atom_name, coordinates) tuples.
        Coordinates should be in Angstroms (will be converted to nm).
    residue_name : str
        Residue name to use in GRO file (default: "LIG")
    title : str
        Title line for GRO file

    Raises
    ------
    ValueError
        If coordinates list is empty

    Examples
    --------
    >>> coords = [(1, "CA", np.array([1.0, 2.0, 3.0]))]
    >>> write_gro_coordinates("output.gro", coords)
    """
    if not coordinates:
        raise ValueError("Cannot write GRO file with no coordinates")

    output_path = Path(output_file)

    with output_path.open("w", encoding="utf-8") as handle:
        # Write title
        handle.write(f"{title}\n")

        # Write atom count
        handle.write(f"{len(coordinates)}\n")

        # Write coordinates (convert Angstroms to nm)
        for atom_idx, atom_name, xyz_angstrom in coordinates:
            xyz_nm = angstrom_xyz_to_nm(xyz_angstrom)

            # GRO format: resid(5) resname(5) atomname(5) atomnr(5) x(8) y(8) z(8)
            # We use simplified format with default resid=1
            handle.write(
                f"{1:5d}{residue_name:5s}{atom_name:>5s}{atom_idx:5d}"
                f"{xyz_nm[0]:8.3f}{xyz_nm[1]:8.3f}{xyz_nm[2]:8.3f}\n"
            )

        # Write box vectors (use 1.0 nm cubic box as default)
        handle.write("   1.00000   1.00000   1.00000\n")


def detect_coordinate_file_format(file_path: Union[str, Path]) -> str:
    """Detect coordinate file format by extension and content.

    Parameters
    ----------
    file_path : Union[str, Path]
        Path to coordinate file

    Returns
    -------
    str
        File format: 'gro', 'mol2', 'pdb', or 'unknown'

    Examples
    --------
    >>> detect_coordinate_file_format("ligand.gro")
    'gro'
    >>> detect_coordinate_file_format("ligand.mol2")
    'mol2'
    """
    path = Path(file_path)

    # Check by extension first
    suffix = path.suffix.lower()
    if suffix == ".gro":
        return "gro"
    elif suffix == ".mol2":
        return "mol2"
    elif suffix == ".pdb":
        return "pdb"

    # Try to detect by content
    try:
        with path.open("r", encoding="utf-8") as handle:
            first_line = handle.readline().strip()

            # GRO files start with title (no special marker)
            # Check second line for atom count
            second_line = handle.readline().strip()
            if second_line.isdigit():
                return "gro"

            # MOL2 files have @<TRIPOS> markers
            if "@<TRIPOS>" in first_line or "@<TRIPOS>" in second_line:
                return "mol2"

            # PDB files start with ATOM or HETATM
            if first_line.startswith("ATOM") or first_line.startswith("HETATM"):
                return "pdb"

    except (IOError, UnicodeDecodeError):
        pass

    return "unknown"
