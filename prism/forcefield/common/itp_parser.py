#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unified ITP file parser for PRISM forcefield module.

This module provides a standardized parser for GROMACS ITP (Topology) files,
consolidating duplicate parsing logic from multiple forcefield generators.

The parser extracts standard sections: moleculetype, atoms, bonds, pairs,
angles, dihedrals, impropers, and atomtypes.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any, Optional, Sequence


@dataclass
class AtomRecord:
    """Atom record from ITP [atoms] section."""

    idx: int
    type_name: str
    resi: int
    resname: str
    atom_name: str
    charge: float
    mass: float


@dataclass
class BondRecord:
    """Bond record from ITP [bonds] section."""

    i: int
    j: int
    length_nm: float
    k_kj_per_nm2: float


@dataclass
class AngleRecord:
    """Angle record from ITP [angles] section."""

    i: int
    j: int
    k: int
    theta_deg: float
    k_kj_per_rad2: float


@dataclass
class DihedralRecord:
    """Dihedral record from ITP [dihedrals] section."""

    i: int
    j: int
    k: int
    l: int
    func: int
    phase_deg: float
    k_kj: float
    multiplicity: int


@dataclass
class AtomTypeRecord:
    """Atom type record from ITP [atomtypes] section."""

    name: str
    sigma_nm: float
    epsilon_kj: float


@dataclass
class MoleculeTypeRecord:
    """Moleculetype record from ITP [moleculetype] section."""

    name: str
    nrexcl: int


class ITPParser:
    """Unified GROMACS ITP file parser.

    This parser handles the subset of GROMACS ITP syntax used by PRISM
    forcefield generators. It extracts all standard sections and provides
    structured access to topology parameters.

    Attributes
    ----------
    path : Path
        Path to ITP file
    moleculetype : Optional[MoleculeTypeRecord]
        Molecule type information
    atoms : List[AtomRecord]
        Atom records
    bonds : List[BondRecord]
        Bond records
    pairs : List[BondRecord]
        Pair 1-4 interactions (stored as BondRecord for simplicity)
    angles : List[AngleRecord]
        Angle records
    dihedrals : List[DihedralRecord]
        Proper dihedral records
    impropers : List[DihedralRecord]
        Improper dihedral records
    atomtypes : List[AtomTypeRecord]
        Atom type records
    alias_map : Dict[str, int]
        Mapping from atom aliases to indices

    Examples
    --------
    >>> parser = ITPParser(Path("ligand.itp"))
    >>> parser.parse()
    >>> print(f"Found {len(parser.atoms)} atoms")
    >>> print(f"Found {len(parser.bonds)} bonds")
    """

    def __init__(self, path: Path):
        """Initialize parser.

        Parameters
        ----------
        path : Path
            Path to ITP file to parse
        """
        self.path = path
        self.moleculetype: Optional[MoleculeTypeRecord] = None
        self.atoms: List[AtomRecord] = []
        self.bonds: List[BondRecord] = []
        self.pairs: List[BondRecord] = []
        self.angles: List[AngleRecord] = []
        self.dihedrals: List[DihedralRecord] = []
        self.impropers: List[DihedralRecord] = []
        self.atomtypes: List[AtomTypeRecord] = []
        self.alias_map: Dict[str, int] = {}

    def parse(self) -> None:
        """Parse ITP file and extract all sections.

        Raises
        ------
        FileNotFoundError
            If ITP file doesn't exist
        ValueError
            If ITP file is malformed
        """
        if not self.path.exists():
            raise FileNotFoundError(f"ITP file not found: {self.path}")

        section: Optional[str] = None
        dihedral_target = "dihedrals"

        with self.path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()

                # Skip empty lines and pure comments
                if not line or line.startswith(";"):
                    # Check for improper/proper dihedral markers in comments
                    if section == "dihedrals" and line.strip():
                        lower = line.lower()
                        if "improper" in lower:
                            dihedral_target = "impropers"
                        elif "proper" in lower:
                            dihedral_target = "dihedrals"
                    continue

                # Remove inline comments
                line = line.split(";", 1)[0].strip()
                if not line:
                    continue

                # Skip standalone comment lines
                if line.startswith("#"):
                    continue

                # Section header
                if line.startswith("["):
                    section = line.strip("[]").split()[0].lower()
                    if section == "dihedrals":
                        dihedral_target = "dihedrals"
                    continue

                # Parse section content
                data = line.split()
                if section == "moleculetype":
                    self._parse_moleculetype(data)
                elif section == "atoms":
                    self._parse_atom(data)
                elif section == "bonds":
                    self._parse_bond(data)
                elif section == "pairs":
                    self._parse_pair(data)
                elif section == "angles":
                    self._parse_angle(data)
                elif section == "dihedrals":
                    record = self._parse_dihedral(data)
                    if dihedral_target == "impropers":
                        self.impropers.append(record)
                    else:
                        self.dihedrals.append(record)
                elif section == "atomtypes":
                    self._parse_atomtype(data)

    def _parse_moleculetype(self, data: Sequence[str]) -> None:
        """Parse [moleculetype] section line."""
        if len(data) < 2:
            return
        self.moleculetype = MoleculeTypeRecord(name=data[0], nrexcl=int(data[1]))

    def _parse_atom(self, data: Sequence[str]) -> None:
        """Parse [atoms] section line."""
        if not data or len(data) < 7:
            return

        # Handle atom aliases (e.g., "CA1" instead of index)
        alias_token: Optional[str] = None
        try:
            idx = int(data[0])
        except ValueError:
            alias_token = data[0]
            idx = len(self.atoms) + 1

        try:
            type_name = data[-7]
            resname = data[-5]
            resi = int(data[-6])
            atom_name = data[-4]
            charge = float(data[-2])
            mass = float(data[-1])
        except (ValueError, IndexError) as exc:
            raise ValueError(f"Invalid atom line: {' '.join(data)}") from exc

        record = AtomRecord(
            idx=idx,
            type_name=type_name,
            resi=resi,
            resname=resname,
            atom_name=atom_name,
            charge=charge,
            mass=mass,
        )
        self.atoms.append(record)

        if alias_token:
            self.alias_map[alias_token] = idx

    def _parse_bond(self, data: Sequence[str]) -> None:
        """Parse [bonds] section line."""
        if len(data) < 5:
            return

        i = self._resolve_index(data[0])
        j = self._resolve_index(data[1])

        try:
            length_nm = float(data[3])
            k_val = float(data[4])
        except ValueError as exc:
            raise ValueError(f"Invalid bond parameters: {' '.join(data)}") from exc

        self.bonds.append(BondRecord(i=i, j=j, length_nm=length_nm, k_kj_per_nm2=k_val))

    def _parse_pair(self, data: Sequence[str]) -> None:
        """Parse [pairs] section line."""
        if len(data) < 5:
            return

        i = self._resolve_index(data[0])
        j = self._resolve_index(data[1])

        try:
            length_nm = float(data[3])
            k_val = float(data[4])
        except ValueError as exc:
            raise ValueError(f"Invalid pair parameters: {' '.join(data)}") from exc

        self.pairs.append(BondRecord(i=i, j=j, length_nm=length_nm, k_kj_per_nm2=k_val))

    def _parse_angle(self, data: Sequence[str]) -> None:
        """Parse [angles] section line."""
        if len(data) < 6:
            return

        i = self._resolve_index(data[0])
        j = self._resolve_index(data[1])
        k = self._resolve_index(data[2])

        try:
            theta = float(data[4])
            k_val = float(data[5])
        except ValueError as exc:
            raise ValueError(f"Invalid angle parameters: {' '.join(data)}") from exc

        self.angles.append(AngleRecord(i=i, j=j, k=k, theta_deg=theta, k_kj_per_rad2=k_val))

    def _parse_dihedral(self, data: Sequence[str]) -> DihedralRecord:
        """Parse [dihedrals] section line."""
        if len(data) < 8:
            raise ValueError(f"Malformed dihedral line: {' '.join(data)}")

        i = self._resolve_index(data[0])
        j = self._resolve_index(data[1])
        k = self._resolve_index(data[2])
        l = self._resolve_index(data[3])

        try:
            func = int(data[4])
            phase = float(data[5])
            k_val = float(data[6])
            multiplicity = int(data[7])
        except ValueError as exc:
            raise ValueError(f"Invalid dihedral parameters: {' '.join(data)}") from exc

        return DihedralRecord(i=i, j=j, k=k, l=l, func=func, phase_deg=phase, k_kj=k_val, multiplicity=multiplicity)

    def _parse_atomtype(self, data: Sequence[str]) -> None:
        """Parse [atomtypes] section line."""
        if len(data) < 4:
            return

        name = data[0]

        try:
            # Handle different atomtype formats
            # Format 1: name mass charge ptype sigma epsilon
            # Format 2: name bond_type mass charge ptype sigma epsilon
            if len(data) >= 7:
                # Has bond_type column
                sigma_nm = float(data[-2])
                epsilon_kj = float(data[-1])
            else:
                # No bond_type column
                sigma_nm = float(data[-2])
                epsilon_kj = float(data[-1])

            self.atomtypes.append(AtomTypeRecord(name=name, sigma_nm=sigma_nm, epsilon_kj=epsilon_kj))
        except ValueError as exc:
            raise ValueError(f"Invalid atomtype line: {' '.join(data)}") from exc

    def _resolve_index(self, token: str) -> int:
        """Resolve atom index from token (integer or alias).

        Parameters
        ----------
        token : str
            Atom index as integer or alias string

        Returns
        -------
        int
            Resolved atom index

        Raises
        ------
        ValueError
            If token cannot be resolved to an index
        """
        # Try direct integer conversion
        try:
            return int(token)
        except ValueError:
            pass

        # Try alias lookup
        if token in self.alias_map:
            return self.alias_map[token]

        raise ValueError(f"Cannot resolve atom index: {token}")

    def get_section(self, section_name: str) -> List[Any]:
        """Get all records from a specific section.

        Parameters
        ----------
        section_name : str
            Name of section ('atoms', 'bonds', 'angles', 'dihedrals', etc.)

        Returns
        -------
        List[Any]
            List of records from the section

        Examples
        --------
        >>> parser = ITPParser(Path("ligand.itp"))
        >>> parser.parse()
        >>> atoms = parser.get_section('atoms')
        """
        section_map = {
            "atoms": self.atoms,
            "bonds": self.bonds,
            "pairs": self.pairs,
            "angles": self.angles,
            "dihedrals": self.dihedrals,
            "impropers": self.impropers,
            "atomtypes": self.atomtypes,
        }

        return section_map.get(section_name, [])

    def extract_section_text(self, section_name: str) -> str:
        """Extract raw text content from a specific section.

        Parameters
        ----------
        section_name : str
            Name of section to extract

        Returns
        -------
        str
            Raw text content of the section

        Raises
        ------
        ValueError
            If section not found

        Examples
        --------
        >>> parser = ITPParser(Path("ligand.itp"))
        >>> atoms_text = parser.extract_section_text('atoms')
        """
        if not self.path.exists():
            raise FileNotFoundError(f"ITP file not found: {self.path}")

        in_section = False
        section_lines = []

        with self.path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()

                # Start of section
                if stripped.startswith(f"[{section_name}]"):
                    in_section = True
                    continue

                # End of section (next section starts)
                if in_section and stripped.startswith("["):
                    break

                # Collect section content
                if in_section:
                    section_lines.append(line)

        if not section_lines:
            raise ValueError(f"Section [{section_name}] not found in {self.path}")

        return "".join(section_lines)

    @staticmethod
    def extract_section_from_content(content: str, section_name: str) -> str:
        """Extract a section from ITP content string.

        This is a static method that extracts a section from ITP file content
        without needing to instantiate the parser.

        Parameters
        ----------
        content : str
            ITP file content as string
        section_name : str
            Name of section to extract

        Returns
        -------
        str
            Raw text content of the section, or empty string if not found

        Examples
        --------
        >>> content = "[ atoms ]\\n...\\n[ bonds ]\\n..."
        >>> section = ITPParser.extract_section_from_content(content, "atoms")
        """
        import re

        pattern = rf"\[\s*{section_name}\s*\](.*?)(?=\[|$)"
        match = re.search(pattern, content, re.DOTALL)
        return match.group(1) if match else ""
