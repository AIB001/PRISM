#!/usr/bin/env python3
"""Utilities to convert Amber/GAFF GROMACS exports (ITP/TOP) into CHARMM RTF/PRM."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import math
import shutil

# Conversion constants
KJ_TO_KCAL = 0.239005736
NM_TO_ANG = 10.0
DEG_TO_RAD = math.pi / 180.0
SQ_NM_TO_SQ_ANG = NM_TO_ANG**2


@dataclass
class AtomRecord:
    idx: int
    type_name: str
    resi: int
    resname: str
    atom_name: str
    charge: float
    mass: float


@dataclass
class BondRecord:
    i: int
    j: int
    length_nm: float
    k_kj_per_nm2: float


@dataclass
class AngleRecord:
    i: int
    j: int
    k: int
    theta_deg: float
    k_kj_per_rad2: float


@dataclass
class DihedralRecord:
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
    name: str
    sigma_nm: float
    epsilon_kj: float


class ItpParser:
    """Lightweight parser for the subset of GROMACS ITP syntax produced by PRISM."""

    def __init__(self, path: Path):
        self.path = path
        self.atoms: List[AtomRecord] = []
        self.bonds: List[BondRecord] = []
        self.angles: List[AngleRecord] = []
        self.dihedrals: List[DihedralRecord] = []
        self.impropers: List[DihedralRecord] = []
        self.alias_map: Dict[str, int] = {}

    def parse(self) -> None:
        if not self.path.exists():
            raise FileNotFoundError(f"ITP file not found: {self.path}")

        section: Optional[str] = None
        dihedral_target = "dihedrals"

        with self.path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                stripped = raw_line.strip()
                if not stripped:
                    continue
                if stripped.startswith(";"):
                    lower = stripped.lower()
                    if section == "dihedrals" and "improper" in lower:
                        dihedral_target = "impropers"
                    elif section == "dihedrals" and "proper" in lower:
                        dihedral_target = "dihedrals"
                    continue

                line = raw_line.split(";", 1)[0].strip()
                if not line:
                    continue
                if line.startswith("#"):
                    continue
                if line.startswith("["):
                    section = line.strip("[]").split()[0].lower()
                    if section == "dihedrals":
                        dihedral_target = "dihedrals"
                    continue

                data = line.split()
                if section == "atoms":
                    self._parse_atom(data)
                elif section == "bonds":
                    self._parse_bond(data)
                elif section == "angles":
                    self._parse_angle(data)
                elif section == "dihedrals":
                    record = self._parse_dihedral(data)
                    if dihedral_target == "impropers":
                        self.impropers.append(record)
                    else:
                        self.dihedrals.append(record)

    def _parse_atom(self, data: Sequence[str]) -> None:
        if not data:
            return
        alias_token: Optional[str] = None
        try:
            idx = int(data[0])
        except ValueError:
            alias_token = data[0]
            idx = len(self.atoms) + 1

        if len(data) < 7:
            raise ValueError(f"Malformed atom line: {' '.join(data)}")

        try:
            mass = float(data[-1])
            charge = float(data[-2])
        except ValueError as exc:
            raise ValueError(f"Invalid charge/mass in atom line: {' '.join(data)}") from exc

        atom_name = data[-4]
        resname = data[-5]
        resi = int(data[-6])
        type_name = data[-7]
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
        if len(data) < 5:
            return
        i = self._resolve_index(data[0])
        j = self._resolve_index(data[1])
        length_nm = float(data[3])
        k_val = float(data[4])
        self.bonds.append(BondRecord(i=i, j=j, length_nm=length_nm, k_kj_per_nm2=k_val))

    def _parse_angle(self, data: Sequence[str]) -> None:
        if len(data) < 6:
            return
        i = self._resolve_index(data[0])
        j = self._resolve_index(data[1])
        k = self._resolve_index(data[2])
        theta = float(data[4])
        k_val = float(data[5])
        self.angles.append(AngleRecord(i=i, j=j, k=k, theta_deg=theta, k_kj_per_rad2=k_val))

    def _parse_dihedral(self, data: Sequence[str]) -> DihedralRecord:
        if len(data) < 8:
            raise ValueError(f"Malformed dihedral line: {' '.join(data)}")
        i = self._resolve_index(data[0])
        j = self._resolve_index(data[1])
        k = self._resolve_index(data[2])
        l = self._resolve_index(data[3])
        func = int(data[4])
        phase = float(data[5])
        k_val = float(data[6])
        multiplicity = int(data[7])
        return DihedralRecord(i=i, j=j, k=k, l=l, func=func, phase_deg=phase, k_kj=k_val, multiplicity=multiplicity)

    def _resolve_index(self, token: str) -> int:
        try:
            return int(token)
        except ValueError:
            if token in self.alias_map:
                return self.alias_map[token]
            raise


def parse_atomtypes(path: Path) -> Dict[str, AtomTypeRecord]:
    if not path.exists():
        raise FileNotFoundError(f"atomtypes_LIG.itp not found: {path}")

    records: Dict[str, AtomTypeRecord] = {}
    with path.open("r", encoding="utf-8") as handle:
        in_section = False
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith("["):
                name = line.strip("[]").split()[0].lower()
                in_section = name == "atomtypes"
                continue
            if not in_section:
                continue
            tokens = line.split()
            if len(tokens) < 7:
                continue
            name = tokens[0]
            sigma = float(tokens[5])
            epsilon = float(tokens[6])
            records[name] = AtomTypeRecord(name=name, sigma_nm=sigma, epsilon_kj=epsilon)
    return records


class AmberToCharmmConverter:
    """Convert Amber/GAFF ligand data from PRISM to CHARMM-compatible files."""

    def __init__(
        self,
        ligand_dir: Path,
        mol2_path: Path,
        output_prefix: Path,
        resname: str = "LIG",
        residue_name: Optional[str] = None,
        copy_mol2: bool = True,
    ) -> None:
        self.ligand_dir = Path(ligand_dir)
        self.mol2_path = Path(mol2_path)
        self.output_prefix = Path(output_prefix)
        self.resname = resname
        self.residue_name = residue_name or resname
        self.copy_mol2 = copy_mol2

        self.itp_path = self.ligand_dir / "LIG.itp"
        self.atomtypes_path = self.ligand_dir / "atomtypes_LIG.itp"
        self.gro_path = self.ligand_dir / "LIG.gro"

        if not self.mol2_path.exists():
            raise FileNotFoundError(f"Input MOL2 not found: {self.mol2_path}")

    def _make_unique_type(self, amber_type: str) -> str:
        """
        生成不超过4字符的GAFF原子类型名字（如 l_ca, l_sy）。
        """
        prefix = (self.residue_name[:1] or "L")[0].lower()
        core = amber_type.lower()
        if len(core) > 3:
            core = core[:3]
        return f"{prefix}{core}"

    def convert(self) -> Dict[str, Path]:
        parser = ItpParser(self.itp_path)
        parser.parse()
        atomtypes = parse_atomtypes(self.atomtypes_path)
        rtf_path = self.output_prefix.with_suffix(".rtf")
        prm_path = self.output_prefix.with_suffix(".prm")
        pdb_path = self.output_prefix.with_suffix(".pdb")
        mol2_out = self.output_prefix.parent / f"{self.output_prefix.name}_3D.mol2"

        self._write_rtf(parser.atoms, parser.bonds, rtf_path)
        self._write_prm(parser, atomtypes, prm_path)
        self._write_pdb(parser.atoms, pdb_path)

        outputs = {
            "rtf": rtf_path,
            "prm": prm_path,
            "pdb": pdb_path,
        }
        if self.copy_mol2:
            self._copy_mol2(mol2_out)
            outputs["mol2"] = mol2_out
        return outputs

    def _write_rtf(self, atoms: Sequence[AtomRecord], bonds: Sequence[BondRecord], path: Path) -> None:
        total_charge = sum(atom.charge for atom in atoms)
        lines: List[str] = []
        lines.append("* Generated by PRISM AmberToCharmmConverter")
        lines.append("*")
        lines.append(f"{len(atoms):5d} 1")
        lines.append("")
        lines.extend(self._format_mass_section(atoms))
        lines.append("")
        lines.append(f"RESI {self.residue_name:<8}{total_charge:12.3f}")
        lines.append("GROUP")
        for atom in atoms:
            unique_type = self._make_unique_type(atom.type_name)
            lines.append(f"ATOM {atom.atom_name:<4} {unique_type:<10} {atom.charge:8.4f} ! mass={atom.mass:6.3f}")
        lines.append("")
        lines.extend(self._format_bonds(atoms, bonds))
        lines.append("")
        lines.append("END\n")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines), encoding="utf-8")

    def _format_bonds(self, atoms: Sequence[AtomRecord], bonds: Sequence[BondRecord]) -> List[str]:
        name_lookup = {atom.idx: atom.atom_name for atom in atoms}
        bond_lines: List[str] = []
        for bond in bonds:
            a = name_lookup.get(bond.i)
            b = name_lookup.get(bond.j)
            if not a or not b:
                continue
            bond_lines.append(f"BOND {a:<4} {b:<4}")
        return bond_lines

    def _write_prm(self, parser: ItpParser, atomtypes: Dict[str, AtomTypeRecord], path: Path) -> None:
        lines: List[str] = []
        lines.append("* Parameters converted from GAFF/Amber by PRISM")
        lines.append("*")
        lines.append("")
        lines.extend(self._format_mass_section(parser.atoms))
        lines.append("")
        lines.append("BONDS")
        for bond in parser.bonds:
            name_i = self._make_unique_type(parser.atoms[bond.i - 1].type_name)
            name_j = self._make_unique_type(parser.atoms[bond.j - 1].type_name)
            k = self._convert_bond_k(bond.k_kj_per_nm2)
            req = bond.length_nm * NM_TO_ANG
            lines.append(f"{name_i:<10}{name_j:<10}{k:10.3f}{req:10.4f}")

        lines.append("")
        lines.append("ANGLES")
        for angle in parser.angles:
            ai = self._make_unique_type(parser.atoms[angle.i - 1].type_name)
            aj = self._make_unique_type(parser.atoms[angle.j - 1].type_name)
            ak = self._make_unique_type(parser.atoms[angle.k - 1].type_name)
            ktheta = self._convert_angle_k(angle.k_kj_per_rad2)
            th0 = angle.theta_deg
            lines.append(f"{ai:<10}{aj:<10}{ak:<10}{ktheta:10.3f}{th0:10.3f}")

        lines.append("")
        lines.append("DIHEDRALS")
        for dih in parser.dihedrals:
            ai = self._make_unique_type(parser.atoms[dih.i - 1].type_name)
            aj = self._make_unique_type(parser.atoms[dih.j - 1].type_name)
            ak = self._make_unique_type(parser.atoms[dih.k - 1].type_name)
            al = self._make_unique_type(parser.atoms[dih.l - 1].type_name)
            kchi = dih.k_kj * KJ_TO_KCAL
            multiplicity = dih.multiplicity if dih.multiplicity > 0 else 1
            lines.append(f"{ai:<10}{aj:<10}{ak:<10}{al:<10}{kchi:10.4f}{multiplicity:4d}{dih.phase_deg:10.3f}")

        lines.append("")
        lines.append("IMPROPERS")
        for imp in parser.impropers:
            ai = self._make_unique_type(parser.atoms[imp.i - 1].type_name)
            aj = self._make_unique_type(parser.atoms[imp.j - 1].type_name)
            ak = self._make_unique_type(parser.atoms[imp.k - 1].type_name)
            al = self._make_unique_type(parser.atoms[imp.l - 1].type_name)
            kpsi = imp.k_kj * KJ_TO_KCAL
            lines.append(f"{ai:<10}{aj:<10}{ak:<10}{al:<10}{kpsi:10.4f}{imp.phase_deg:10.3f}")

        lines.append("")
        lines.extend(self._format_nonbonded(atomtypes))
        lines.append("END\n")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines), encoding="utf-8")

    def _convert_bond_k(self, k_kj_per_nm2: float) -> float:
        return 0.5 * k_kj_per_nm2 * KJ_TO_KCAL / SQ_NM_TO_SQ_ANG

    def _convert_angle_k(self, k_kj_per_rad2: float) -> float:
        return 0.5 * k_kj_per_rad2 * KJ_TO_KCAL

    def _format_mass_section(self, atoms: Sequence[AtomRecord]) -> List[str]:
        masses: Dict[str, float] = {}
        elements: Dict[str, str] = {}
        for atom in atoms:
            unique_type = self._make_unique_type(atom.type_name)
            masses.setdefault(unique_type, atom.mass)
            elements.setdefault(unique_type, self._guess_element(atom.atom_name))
        entries: List[str] = []
        for idx, name in enumerate(sorted(masses.keys()), start=1):
            mass = masses[name]
            element = elements.get(name, name[0].upper())
            entries.append(f"MASS {idx:5d} {name:<10} {mass:10.4f} {element}")
        return entries

    def _format_nonbonded(self, atomtypes: Dict[str, AtomTypeRecord]) -> List[str]:
        lines = ["NONBONDED nbxmod 5 atom cdiel shift vatom vswitch vfswitch -1 0.0"]
        lines.append("!atom  ignored   epsilon   Rmin/2  ignored  eps,1-4 Rmin/2,1-4")
        for name, record in atomtypes.items():
            unique_type = self._make_unique_type(name)
            if record.sigma_nm == 0.0 and record.epsilon_kj == 0.0:
                rmin_over2 = 0.0
            else:
                rmin = (2 ** (1 / 6)) * record.sigma_nm * NM_TO_ANG
                rmin_over2 = rmin / 2.0
            eps = -abs(record.epsilon_kj * KJ_TO_KCAL)
            lines.append(
                f"{unique_type:<10}{0.0:10.2f}{eps:10.5f}{rmin_over2:10.4f}{0.0:10.2f}{eps:10.5f}{rmin_over2:10.4f}"
            )
        return lines

    def _write_pdb(self, atoms: Sequence[AtomRecord], path: Path) -> None:
        try:
            coordinates = self._read_gro_coordinates(len(atoms))
        except Exception:
            coordinates = self._read_mol2_coordinates(len(atoms))
        if len(coordinates) != len(atoms):
            raise ValueError("Mismatch between coordinate sources and atom records")

        lines: List[str] = []
        for atom, coords in zip(atoms, coordinates):
            x, y, z = coords
            element = self._guess_element(atom.atom_name)
            lines.append(
                f"HETATM{atom.idx:5d} {atom.atom_name:>4} {self.resname:<3} A{atom.resi:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"
            )
        lines.append("END")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines), encoding="utf-8")

    def _read_gro_coordinates(self, expected_atoms: int) -> List[Tuple[float, float, float]]:
        if not self.gro_path.exists():
            raise FileNotFoundError(f"GRO file not found: {self.gro_path}")
        coords: List[Tuple[float, float, float]] = []
        with self.gro_path.open("r", encoding="utf-8") as handle:
            header = handle.readline()
            count_line = handle.readline()
            if not count_line:
                raise ValueError("Malformed GRO file: missing atom count")
            total_atoms = int(count_line.strip())
            for _ in range(total_atoms):
                line = handle.readline()
                if not line:
                    break
                x = float(line[20:28]) * NM_TO_ANG
                y = float(line[28:36]) * NM_TO_ANG
                z = float(line[36:44]) * NM_TO_ANG
                coords.append((x, y, z))
        if len(coords) < expected_atoms:
            raise ValueError("Not enough coordinate entries in GRO file")
        return coords[:expected_atoms]

    def _read_mol2_coordinates(self, expected_atoms: int) -> List[Tuple[float, float, float]]:
        coords: List[Tuple[float, float, float]] = []
        with self.mol2_path.open("r", encoding="utf-8") as handle:
            in_atoms = False
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if line.startswith("@<TRIPOS>") and in_atoms:
                    break
                if not in_atoms:
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                coords.append((x, y, z))
        if len(coords) < expected_atoms:
            raise ValueError("Not enough coordinate entries in MOL2 file")
        return coords[:expected_atoms]

    def _guess_element(self, atom_name: str) -> str:
        stripped = atom_name.strip()
        if not stripped:
            return "C"
        first = stripped[0]
        second = stripped[1] if len(stripped) > 1 and stripped[1].islower() else ""
        if first.isdigit():
            first = stripped[1] if len(stripped) > 1 else "C"
            second = stripped[2] if len(stripped) > 2 and stripped[2].islower() else ""
        symbol = (first + second).upper()
        if len(symbol) == 2 and symbol[1] not in ("H", "B", "C", "N", "O", "F", "P", "S", "K", "V", "Y", "I", "W"):
            symbol = symbol[0]
        if symbol[0] not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            return "C"
        return symbol[:2].strip() or "C"

    def _copy_mol2(self, destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(self.mol2_path, destination)


def _build_cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert PRISM Amber/GAFF ligand outputs (LIG.amb2gmx) to CHARMM RTF/PRM files."
    )
    parser.add_argument("--ligand-dir", required=True, type=Path, help="Path to the LIG.amb2gmx directory")
    parser.add_argument(
        "--mol2",
        required=True,
        type=Path,
        help="Source MOL2 file (optionally copied as <prefix>_3D.mol2 with --copy-mol2)",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help="Output prefix (e.g., /path/to/38 -> writes 38.rtf/38.prm/38.pdb)",
    )
    parser.add_argument("--resname", default="LIG", help="Three/four-letter residue name for topology/PDB")
    parser.add_argument(
        "--residue-name",
        default=None,
        help="Residue label used inside RTF RESI line (defaults to --resname)",
    )
    parser.add_argument(
        "--no-copy-mol2",
        action="store_true",
        help="Skip copying the MOL2 (default is to write <prefix>_3D.mol2)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = _build_cli()
    args = parser.parse_args(argv)
    converter = AmberToCharmmConverter(
        ligand_dir=args.ligand_dir,
        mol2_path=args.mol2,
        output_prefix=args.output_prefix,
        resname=args.resname,
        residue_name=args.residue_name,
        copy_mol2=not args.no_copy_mol2,
    )
    outputs = converter.convert()
    for key, value in outputs.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()


__all__ = ["AmberToCharmmConverter"]
