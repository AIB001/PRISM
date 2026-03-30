#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Hybrid topology package builder for FEP calculations.

This module handles the assembly of hybrid ligand force field files (ITP, GRO, atomtypes)
from reference and mutant ligand components.
"""

import shutil
from pathlib import Path
from typing import Dict, Optional

from prism.fep.gromacs.itp_builder import ITPBuilder


class HybridPackageBuilder:
    """
    Builder for hybrid ligand force field packages.

    Assembles hybrid.itp, atomtypes_hybrid.itp, ff_hybrid.itp, hybrid.gro,
    and ligand_seed.pdb from reference and mutant ligand directories.
    """

    def __init__(
        self,
        hybrid_itp_filename: str = "hybrid.itp",
        hybrid_atomtypes_filename: str = "atomtypes_hybrid.itp",
        hybrid_forcefield_filename: str = "ff_hybrid.itp",
        hybrid_posre_filename: str = "posre_hybrid.itp",
        hybrid_gro_filename: str = "hybrid.gro",
        molecule_name: str = "HYB",
    ):
        """
        Initialize with customizable filenames.

        Parameters
        ----------
        hybrid_itp_filename : str
            Filename for hybrid ITP file (default: "hybrid.itp")
        hybrid_atomtypes_filename : str
            Filename for hybrid atomtypes file (default: "atomtypes_hybrid.itp")
        hybrid_forcefield_filename : str
            Filename for hybrid forcefield file (default: "ff_hybrid.itp")
        hybrid_posre_filename : str
            Filename for hybrid position restraint file (default: "posre_hybrid.itp")
        hybrid_gro_filename : str
            Filename for hybrid GRO file (default: "hybrid.gro")
        molecule_name : str
            Molecule name for hybrid topology (default: "HYB")
        """
        # Filenames for hybrid assets
        self.hybrid_itp_filename = hybrid_itp_filename
        self.hybrid_atomtypes_filename = hybrid_atomtypes_filename
        self.hybrid_forcefield_filename = hybrid_forcefield_filename
        self.hybrid_posre_filename = hybrid_posre_filename
        self.hybrid_gro_filename = hybrid_gro_filename
        self.molecule_name = molecule_name

    def prepare_hybrid_package_from_components(
        self,
        hybrid_dir: Path,
        hybrid_itp: Path,
        reference_ligand_dir: Path,
        mutant_ligand_dir: Optional[Path],
        seed_gro: Path,
    ) -> None:
        """
        Build complete hybrid ligand package from components.

        Parameters
        ----------
        hybrid_dir : Path
            Output directory for hybrid assets
        hybrid_itp : Path
            Path to hybrid ITP file (atom mapping)
        reference_ligand_dir : Path
            Directory containing reference ligand force field
        mutant_ligand_dir : Optional[Path]
            Directory containing mutant ligand force field (None for single-topology)
        seed_gro : Path
            Reference GRO file for coordinates
        """
        # Normalize and write hybrid ITP
        normalized_itp = self._normalize_hybrid_itp(hybrid_itp.read_text(), self.molecule_name)
        hybrid_itp_output = hybrid_dir / self.hybrid_itp_filename
        hybrid_itp_output.write_text(normalized_itp)

        # Add bonded sections if needed
        if self._itp_needs_bonded_sections(normalized_itp):
            # Auto-discover ligand ITP files (handle LIG.amb2gmx/ subdirectory)
            def _resolve_ligand_itp(ligand_dir: Path) -> str:
                """Resolve ligand ITP path, handling LIG.amb2gmx/ subdirectory."""
                # Try direct path
                direct_path = ligand_dir / "LIG.itp"
                if direct_path.exists():
                    return str(direct_path)

                # Try LIG.amb2gmx/LIG.itp (PRISM output structure)
                amb2gmx_path = ligand_dir / "LIG.amb2gmx" / "LIG.itp"
                if amb2gmx_path.exists():
                    return str(amb2gmx_path)

                raise FileNotFoundError(f"Cannot find LIG.itp in {ligand_dir}")

            ligand_a_itp = _resolve_ligand_itp(reference_ligand_dir)
            ligand_b_itp = _resolve_ligand_itp(mutant_ligand_dir) if mutant_ligand_dir else ligand_a_itp

            ITPBuilder.write_complete_hybrid_itp(
                output_path=str(hybrid_itp_output),
                hybrid_itp=str(hybrid_itp_output),
                ligand_a_itp=ligand_a_itp,
                ligand_b_itp=ligand_b_itp,
                molecule_name=self.molecule_name,
            )
            normalized_itp = hybrid_itp_output.read_text()

        # Build atomtypes file
        atomtypes = self._collect_atomtypes(reference_ligand_dir, mutant_ligand_dir)
        atomtypes_content = self._build_hybrid_atomtypes_itp(normalized_itp, atomtypes)
        (hybrid_dir / self.hybrid_atomtypes_filename).write_text(atomtypes_content)

        # Build force field file WITHOUT [ defaults ] block
        # Note: ff_hybrid.itp will be included in the main topology which
        # already has a [ defaults ] block from the protein force field.
        # Including [ defaults ] here would cause "Found a second defaults directive" error.
        ff_hybrid = f'#include "{self.hybrid_atomtypes_filename}"\n'
        (hybrid_dir / self.hybrid_forcefield_filename).write_text(ff_hybrid)

        # Build hybrid GRO file
        # Auto-discover mutant GRO file (handle LIG.amb2gmx/ subdirectory)
        mutant_gro = None
        if mutant_ligand_dir:
            # Try direct path
            direct_path = mutant_ligand_dir / "LIG.gro"
            if direct_path.exists():
                mutant_gro = direct_path
            else:
                # Try LIG.amb2gmx/LIG.gro (PRISM output structure)
                amb2gmx_path = mutant_ligand_dir / "LIG.amb2gmx" / "LIG.gro"
                if amb2gmx_path.exists():
                    mutant_gro = amb2gmx_path

        hybrid_gro_content = self._build_hybrid_gro(
            hybrid_itp_content=normalized_itp,
            reference_gro=seed_gro,
            mutant_gro=mutant_gro,
        )
        (hybrid_dir / self.hybrid_gro_filename).write_text(hybrid_gro_content)

        # Convert to PDB for visualization
        (hybrid_dir / "ligand_seed.pdb").write_text(self._gro_to_pdb(hybrid_dir / self.hybrid_gro_filename))

        # Copy position restraints if available
        # Auto-discover posre file (handle LIG.amb2gmx/ subdirectory)
        ref_posre = reference_ligand_dir / "posre_LIG.itp"
        if not ref_posre.exists():
            # Try LIG.amb2gmx/posre_LIG.itp (PRISM output structure)
            ref_posre = reference_ligand_dir / "LIG.amb2gmx" / "posre_LIG.itp"

        if ref_posre.exists():
            shutil.copy2(ref_posre, hybrid_dir / self.hybrid_posre_filename)

        # Write topology template
        self._write_hybrid_top_template(hybrid_dir / "hybrid.top")

    def _normalize_hybrid_itp(self, content: str, molecule_name: str) -> str:
        """
        Normalize hybrid ITP content.

        - Replace moleculetype name with molecule_name
        - Normalize atom lines to ensure consistent formatting
        """
        lines = []
        in_moleculetype = False
        in_atoms = False

        for line in content.splitlines():
            stripped = line.strip()

            if stripped.lower() == "[ moleculetype ]":
                in_moleculetype = True
                lines.append(line)
                continue

            if in_moleculetype and stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if len(parts) >= 2:
                    lines.append(f"{molecule_name}  {parts[1]}")
                    in_moleculetype = False
                    continue

            if stripped.lower() == "[ atoms ]":
                in_atoms = True
                lines.append(line)
                continue

            if in_atoms and stripped.startswith("["):
                in_atoms = False

            if in_atoms and stripped and not stripped.startswith(";"):
                lines.append(self._normalize_hybrid_atom_line(line))
                continue

            lines.append(line)

        return "\n".join(lines) + "\n"

    def _normalize_hybrid_atom_line(self, line: str) -> str:
        """Normalize hybrid atom line formatting."""
        return line.rstrip()

    def _itp_needs_bonded_sections(self, content: str) -> bool:
        """Check if ITP file is missing bonded sections."""
        lower = content.lower()
        return "[ bonds ]" not in lower and "[ angles ]" not in lower and "[ dihedrals ]" not in lower

    def _collect_atomtypes(self, reference_ligand_dir: Path, mutant_ligand_dir: Optional[Path]) -> Dict[str, str]:
        """
        Collect atomtype definitions from ligand directories.

        Returns
        -------
        Dict[str, str]
            Mapping of atomtype name to full definition line
        """
        atomtypes = {}
        for ligand_dir in [reference_ligand_dir, mutant_ligand_dir]:
            if ligand_dir is None:
                continue

            # Try to find atomtypes_LIG.itp (handle LIG.amb2gmx/ subdirectory)
            atomtypes_file = ligand_dir / "atomtypes_LIG.itp"
            if not atomtypes_file.exists():
                # Try LIG.amb2gmx/atomtypes_LIG.itp (PRISM output structure)
                atomtypes_file = ligand_dir / "LIG.amb2gmx" / "atomtypes_LIG.itp"

            if not atomtypes_file.exists():
                continue

            for atomtype_name, line in self._parse_atomtypes_file(atomtypes_file).items():
                atomtypes.setdefault(atomtype_name, line)
        return atomtypes

    def _parse_atomtypes_file(self, atomtypes_file: Path) -> Dict[str, str]:
        """Parse atomtypes ITP file."""
        atomtypes = {}
        in_section = False
        for line in atomtypes_file.read_text().splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ atomtypes ]":
                in_section = True
                continue
            if stripped.startswith("[") and in_section:
                break
            if not in_section or not stripped or stripped.startswith(";"):
                continue
            parts = stripped.split()
            atomtypes[parts[0]] = stripped
        return atomtypes

    def _build_hybrid_atomtypes_itp(self, hybrid_itp_content: str, source_atomtypes: Dict[str, str]) -> str:
        """
        Build atomtypes ITP file for hybrid ligand.

        Extracts used atomtypes from hybrid ITP and generates definitions,
        including dummy atomtypes (DUM_*) with zero LJ parameters.
        """
        used_types = self._parse_used_atomtypes(hybrid_itp_content)
        lines = [
            "[ atomtypes ]",
            "; merged atomtypes for hybrid ligand",
        ]

        written = set()
        dummy_requests = set()

        # Collect real and dummy atomtypes.
        # Real atomtypes are written here so unbound leg (ligand-only, no protein FF)
        # can find them. Bound leg topology inlines them separately from PRISM build,
        # so duplication is avoided by grompp's deduplication (same values = no error).
        real_requests = set()
        for atomtype in used_types:
            if atomtype.startswith("DUM_"):
                dummy_requests.add(atomtype)
            else:
                real_requests.add(atomtype)

        # Write real atomtypes first
        for atomtype in sorted(real_requests):
            if atomtype in written:
                continue
            base_line = source_atomtypes.get(atomtype)
            if base_line:
                lines.append(base_line)
                written.add(atomtype)

        # Write dummy atomtypes
        for dummy_atomtype in sorted(dummy_requests):
            if dummy_atomtype in written:
                continue
            base_type = dummy_atomtype[4:]
            base_line = source_atomtypes.get(base_type)
            bond_type = base_type
            if base_line:
                bond_type = base_line.split()[1]
            lines.append(
                f"{dummy_atomtype:<8s} {bond_type:<8s} 0.00000 0.00000 A 0.00000e+00 0.00000e+00 ; dummy from {base_type}"
            )
            written.add(dummy_atomtype)

        # Check for missing atomtypes
        missing_non_dummy = sorted(
            atomtype for atomtype in used_types if not atomtype.startswith("DUM_") and atomtype not in source_atomtypes
        )
        if missing_non_dummy:
            raise ValueError(f"Missing atomtype definitions for: {', '.join(missing_non_dummy)}")

        return "\n".join(lines) + "\n"

    def _extract_defaults_block(self, reference_ligand_dir: Path, mutant_ligand_dir: Optional[Path]) -> str:
        """
        Extract [ defaults ] block from ligand topologies.

        Falls back to inference if not found.
        """
        defaults_blocks = []
        for ligand_dir in [reference_ligand_dir, mutant_ligand_dir]:
            if ligand_dir is None:
                continue
            top_file = ligand_dir / "LIG.top"
            if not top_file.exists():
                continue
            block = self._parse_defaults_block(top_file.read_text())
            if block:
                defaults_blocks.append(block)

        if not defaults_blocks:
            return self._infer_defaults_block(reference_ligand_dir, mutant_ligand_dir)

        first = defaults_blocks[0]
        for block in defaults_blocks[1:]:
            if block != first:
                raise ValueError("Inconsistent [ defaults ] blocks between ligand topologies")
        return first

    def _infer_defaults_block(self, reference_ligand_dir: Path, mutant_ligand_dir: Optional[Path]) -> str:
        """Infer [ defaults ] block from atomtype naming conventions."""
        atomtypes = self._collect_atomtypes(reference_ligand_dir, mutant_ligand_dir)
        if atomtypes and all(name.startswith("opls_") or name.startswith("DUM_opls_") for name in atomtypes):
            return (
                "[ defaults ]\n"
                "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n"
                "1               3               yes             0.5     0.5\n"
            )
        return (
            "[ defaults ]\n"
            "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n"
            "1               2               yes             0.5     0.8333\n"
        )

    def _parse_defaults_block(self, top_content: str) -> str:
        """Parse [ defaults ] block from topology file."""
        lines = []
        in_section = False
        for line in top_content.splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ defaults ]":
                in_section = True
                lines.append(line)
                continue
            if in_section:
                if stripped.startswith("["):
                    break
                lines.append(line)
        return "\n".join(lines) + "\n" if lines else ""

    def _build_hybrid_gro(
        self,
        hybrid_itp_content: str,
        reference_gro: Path,
        mutant_gro: Optional[Path],
    ) -> str:
        """
        Build hybrid GRO file with coordinates from reference/mutant.

        Parameters
        ----------
        hybrid_itp_content : str
            Content of hybrid ITP file
        reference_gro : Path
            Reference ligand GRO file
        mutant_gro : Optional[Path]
            Mutant ligand GRO file (None for single-topology)

        Returns
        -------
        str
            Complete GRO file content
        """
        ref_coords, ref_box = self._parse_gro_atoms(reference_gro)
        mut_coords, _ = self._parse_gro_atoms(mutant_gro) if mutant_gro else ({}, None)
        hybrid_atoms, correspondence = self._parse_hybrid_atoms(hybrid_itp_content)

        lines = [
            "Hybrid ligand structure",
            f"{len(hybrid_atoms):5d}",
        ]

        for idx, atom in enumerate(hybrid_atoms, start=1):
            x, y, z = self._resolve_hybrid_atom_coordinates(atom, ref_coords, mut_coords, correspondence)
            lines.append(f"{1:5d}{self.molecule_name:<5s}{atom['name']:>5s}{idx:5d}{x:8.3f}{y:8.3f}{z:8.3f}")

        lines.append(self._normalize_box_line(ref_box))
        return "\n".join(lines) + "\n"

    def _parse_gro_atoms(self, gro_path: Optional[Path]) -> tuple[Dict[str, list[Dict[str, float]]], Optional[str]]:
        """
        Parse GRO file and extract atom coordinates.

        Returns
        -------
        tuple
            (atom_coords_by_name, box_line)
        """
        if gro_path is None or not gro_path.exists():
            return {}, None

        lines = gro_path.read_text().splitlines()
        if len(lines) < 3:
            return {}, None

        natoms = int(lines[1].strip())
        coords_by_name = {}

        for i in range(2, 2 + natoms):
            if i >= len(lines):
                break
            line = lines[i]
            if len(line) < 44:
                continue

            atom_name = line[10:15].strip()
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])

            coords_by_name.setdefault(atom_name, []).append({"x": x, "y": y, "z": z})

        box_line = lines[2 + natoms] if len(lines) > 2 + natoms else None
        return coords_by_name, box_line

    def _parse_hybrid_atoms(self, hybrid_itp_content: str) -> tuple[list[Dict[str, str]], Dict[str, str]]:
        """
        Parse [ atoms ] section from hybrid ITP.

        Also extracts FEbuilder correspondence comments.

        Returns
        -------
        tuple[list[Dict[str, str]], Dict[str, str]]
            (List of atom dictionaries, correspondence map)
            correspondence map: state A name -> state B name
        """
        atoms = []
        correspondence = {}
        in_atoms = False
        pending_correspondence = None

        for line in hybrid_itp_content.splitlines():
            stripped = line.strip()

            # Check for FEbuilder correspondence comment
            if stripped.startswith(";") and "name_b (state B):" in stripped:
                pending_correspondence = stripped.split("name_b (state B):", 1)[1].strip()
                continue

            if stripped.lower() == "[ atoms ]":
                in_atoms = True
                continue
            if in_atoms and stripped.startswith("["):
                break
            if not in_atoms or not stripped or stripped.startswith(";"):
                continue

            parts = stripped.split()
            if len(parts) < 5:
                continue

            atom_name = parts[4]
            atoms.append(
                {
                    "name": atom_name,
                    "type": parts[1],
                    "resname": parts[3],
                }
            )

            # Record correspondence if available
            if pending_correspondence:
                correspondence[atom_name] = pending_correspondence
                pending_correspondence = None

        return atoms, correspondence

    def _resolve_hybrid_atom_coordinates(
        self,
        atom: Dict[str, str],
        ref_coords: Dict[str, list[Dict[str, float]]],
        mut_coords: Dict[str, list[Dict[str, float]]],
        correspondence: Dict[str, str],
    ) -> tuple[float, float, float]:
        """
        Resolve coordinates for a hybrid atom.

        For dual-topology FEP:
        - Real atoms in state A get coordinates from reference
        - Real atoms in state B (DUM_ in state A) get coordinates from mutant
        - If FEbuilder correspondence exists, use corresponding atom coordinates from mutant

        Parameters
        ----------
        atom : Dict[str, str]
            Atom dictionary with keys: name, type, resname
        ref_coords : Dict[str, list[Dict[str, float]]]
            Reference ligand coordinates (by atom name)
        mut_coords : Dict[str, list[Dict[str, float]]]
            Mutant ligand coordinates (by atom name)
        correspondence : Dict[str, str]
            FEbuilder correspondence map (state A name -> state B name)

        Returns
        -------
        tuple[float, float, float]
            (x, y, z) coordinates
        """
        atom_name = atom["name"]
        atom_type = atom["type"]

        # Try reference first for state A atoms
        if atom_name in ref_coords and ref_coords[atom_name]:
            coord = ref_coords[atom_name][0]
            return coord["x"], coord["y"], coord["z"]

        # If dummy in reference, try mutant coordinates
        if atom_type.startswith("DUM_") and atom_name in mut_coords and mut_coords[atom_name]:
            coord = mut_coords[atom_name][0]
            return coord["x"], coord["y"], coord["z"]

        # If atom has FEbuilder correspondence, try using corresponding atom from mutant
        if atom_name in correspondence:
            corresponding_name_b = correspondence[atom_name]
            # Try with B suffix (e.g., H08 -> H08B)
            for suffix in ["", "B"]:
                candidate_name = f"{corresponding_name_b}{suffix}"
                if candidate_name in mut_coords and mut_coords[candidate_name]:
                    coord = mut_coords[candidate_name][0]
                    return coord["x"], coord["y"], coord["z"]

        # Default to origin (should not happen)
        return 0.0, 0.0, 0.0

    def _normalize_box_line(self, box_line: Optional[str]) -> str:
        """Normalize GRO box line."""
        if box_line is None:
            return "   5.00000   5.00000   5.00000"
        return box_line.strip()

    def _parse_used_atomtypes(self, hybrid_itp_content: str) -> set[str]:
        """
        Extract all atomtypes used in hybrid ITP [ atoms ] section.

        Returns
        -------
        set[str]
            Set of atomtype names
        """
        used_types = set()
        in_atoms = False

        for line in hybrid_itp_content.splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ atoms ]":
                in_atoms = True
                continue
            if in_atoms and stripped.startswith("["):
                break
            if not in_atoms or not stripped or stripped.startswith(";"):
                continue

            parts = stripped.split()
            if len(parts) >= 2:
                used_types.add(parts[1])  # typeA
            if len(parts) >= 9:
                used_types.add(parts[8])  # typeB (dual-topology)

        return used_types

    def _looks_like_float(self, value: str) -> bool:
        """Check if string looks like a float."""
        try:
            float(value)
            return True
        except ValueError:
            return False

    def _write_hybrid_top_template(self, output_path: Path) -> None:
        """Write hybrid topology template file."""
        content = f"""
; Hybrid ligand topology template
; This file is for reference only - actual topology is built by PRISM

#include "{self.hybrid_forcefield_filename}"
#include "{self.hybrid_itp_filename}"

[ system ]
Hybrid ligand

[ molecules ]
{self.molecule_name}  1
"""
        output_path.write_text(content.lstrip())

    def _gro_to_pdb(self, gro_file: Path) -> str:
        """
        Convert GRO file to PDB format.

        Simple conversion for visualization purposes.
        """
        lines = gro_file.read_text().splitlines()
        pdb_lines = ["REMARK   Generated from LIG.gro for PRISM-FEP scaffold"]
        atom_lines = lines[2:-1]

        for serial, line in enumerate(atom_lines, start=1):
            if len(line) < 44:
                continue
            resnum = int(line[:5].strip())
            resname = line[5:10].strip() or "LIG"
            atom_name = line[10:15].strip() or f"A{serial}"
            x = float(line[20:28].strip()) * 10.0
            y = float(line[28:36].strip()) * 10.0
            z = float(line[36:44].strip()) * 10.0
            element = "".join([c for c in atom_name if c.isalpha()])[:1].upper() or "C"
            pdb_lines.append(
                f"HETATM{serial:5d} {atom_name:>4s} {resname:>3s} A{resnum:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
            )

        pdb_lines.append("END")
        return "\n".join(pdb_lines) + "\n"
