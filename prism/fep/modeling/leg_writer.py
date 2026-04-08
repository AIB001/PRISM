#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP leg directory writer and topology manager.

This module handles writing bound/unbound leg directories, copying PRISM systems,
modifying topologies, and managing coordinate files.
"""

import os
import re
import shutil
import subprocess
from pathlib import Path


class LegWriter:
    """
    Writer for FEP leg directories (bound/unbound).

    Handles copying PRISM-built systems, replacing ligands with hybrid ligands,
    modifying topologies, and generating position restraints.
    """

    def __init__(self, output_dir: Path, molecule_name: str = "HYB"):
        """
        Initialize leg writer.

        Parameters
        ----------
        output_dir : Path
            Root output directory for FEP scaffold
        molecule_name : str
            Hybrid molecule name (default: "HYB")
        """
        self.output_dir = output_dir
        self.molecule_name = molecule_name

    def copy_prism_system_to_leg(self, prism_system_dir: Path, leg_dir: Path, leg_name: str) -> None:
        """
        Copy files from PRISM-built system to FEP leg directory.

        Parameters
        ----------
        prism_system_dir : Path
            Path to GMX_PROLIG_MD directory from PRISM
        leg_dir : Path
            Target FEP leg directory (bound or unbound)
        leg_name : str
            "bound" or "unbound"
        """
        print(f"\n{'='*70}")
        print(f"Copying {leg_name} leg from PRISM system")
        print(f"{'='*70}")

        # Check source files exist
        solv_ions_gro = prism_system_dir / "solv_ions.gro"
        topol_top = prism_system_dir / "topol.top"

        if not solv_ions_gro.exists():
            raise FileNotFoundError(f"solv_ions.gro not found in {prism_system_dir}")
        if not topol_top.exists():
            raise FileNotFoundError(f"topol.top not found in {prism_system_dir}")

        # Check if hybrid GRO exists
        hybrid_gro = self.output_dir / "common" / "hybrid" / "hybrid.gro"
        if hybrid_gro.exists():
            # Replace ligand atoms with hybrid ligand atoms
            self.replace_ligand_with_hybrid(solv_ions_gro, leg_dir / "input" / "conf.gro", hybrid_gro)
            print(f"  ✓ Replaced ligand with hybrid ligand → {leg_name}/input/conf.gro")
        else:
            # Copy coordinate file as-is
            shutil.copy2(solv_ions_gro, leg_dir / "input" / "conf.gro")
            print(f"  ✓ Copied solv_ions.gro → {leg_name}/input/conf.gro")

        # Copy all topology ITP files (chain topologies, posre files, etc.)
        self.copy_topology_itp_files(prism_system_dir, leg_dir)

        # Copy and modify topology to use hybrid ligand
        self.copy_and_modify_topology(topol_top, leg_dir / "topol.top")

        if hybrid_gro.exists():
            removed_sol = self.prune_overlapping_solvent(leg_dir / "input" / "conf.gro", leg_dir / "topol.top")
            if removed_sol:
                print(f"  ✓ Removed {removed_sol} overlapping SOL molecules after hybrid replacement")

        # Unbound leg topology from ligand-only builder may miss solvent/ion counts;
        # rebuild [ molecules ] from generated coordinate file to keep topology consistent.
        if leg_name == "unbound":
            self.sync_unbound_molecules_with_conf(leg_dir / "topol.top", leg_dir / "input" / "conf.gro")

        print(f"  ✓ Copied and modified topol.top → {leg_name}/topol.top")

        # Generate position restraint file for protein
        self.generate_posre_file(leg_dir, leg_name)

        print(f"  ✓ {leg_name.capitalize()} leg system ready")
        print(f"{'='*70}\n")

    def copy_topology_itp_files(self, source_dir: Path, target_dir: Path) -> None:
        """
        Copy all topology ITP files from source to target directory.

        This includes chain topology files (topol_Protein_*.itp) and
        position restraint files (posre_*.itp) that are referenced by topol.top.

        Force-field directories are intentionally not copied into each leg.
        Topologies keep the standard `<forcefield>.ff/...` include paths and rely
        on GROMACS search paths (working directory / GMXLIB / installation).
        That avoids vendoring redundant `.ff` trees into every repeat directory.

        Parameters
        ----------
        source_dir : Path
            Source directory (GMX_PROLIG_MD)
        target_dir : Path
            Target directory (bound/unbound leg)
        """
        # LIG* ITPs will be replaced with hybrid versions; skip them
        _lig_itps = {"LIG.itp", "posre_LIG.itp", "atomtypes_LIG.itp", "defaults_LIG.itp"}
        copied_count = 0
        for itp_file in source_dir.glob("*.itp"):
            if itp_file.name not in _lig_itps:
                shutil.copy2(itp_file, target_dir / itp_file.name)
                copied_count += 1
        if copied_count:
            print(f"  ✓ Copied {copied_count} topology ITP files")

    def generate_posre_file(self, leg_dir: Path, leg_name: str) -> None:
        """
        Generate position restraint file for protein atoms.

        For bound leg: extracts protein atoms from GRO, generates restraints with molecule-relative indices.
        For unbound leg: creates empty posre.itp to satisfy topology include.

        Parameters
        ----------
        leg_dir : Path
            FEP leg directory (bound or unbound)
        leg_name : str
            "bound" or "unbound"
        """
        posre_file = leg_dir / "posre.itp"
        conf_gro = leg_dir / "input" / "conf.gro"

        if not conf_gro.exists():
            print(f"  ⚠ Warning: {conf_gro} not found, skipping posre.itp generation")
            return

        # For unbound leg, create empty posre.itp
        if leg_name != "bound":
            with open(posre_file, "w") as f:
                f.write("; Position restraint file (empty for unbound leg)\n")
                f.write("[ position_restraints ]\n")
                f.write("; Empty - no position restraints for unbound leg\n")
            print(f"  ✓ Created empty posre.itp for {leg_name} leg")
            return

        # For bound leg: extract protein atoms manually from GRO, then generate posre
        try:
            protein_gro = leg_dir / "protein_only.gro"

            # Read conf.gro and extract protein atoms (skip hybrid ligand)
            lines = conf_gro.read_text().splitlines()
            if len(lines) < 3:
                raise ValueError("Invalid GRO file format")

            title = lines[0]
            n_atoms = int(lines[1].strip())
            atom_lines = lines[2 : 2 + n_atoms]
            box_line = lines[2 + n_atoms] if len(lines) > 2 + n_atoms else "   10.00000   10.00000   10.00000"

            # Filter protein atoms (residue name not HYB/LIG/SOL/ions)
            protein_atoms = []
            for line in atom_lines:
                if len(line) < 20:
                    continue
                resname = line[5:10].strip()
                if resname not in ["HYB", "LIG", "SOL", "NA", "CL", "K", "MG", "CA"]:
                    protein_atoms.append(line)

            if not protein_atoms:
                print(f"  ⚠ Warning: No protein atoms found in {conf_gro}")
                # Create empty posre
                with open(posre_file, "w") as f:
                    f.write("; Position restraint file (no protein found)\n")
                    f.write("[ position_restraints ]\n")
                return

            # Write protein-only GRO
            with open(protein_gro, "w") as f:
                f.write(f"{title} (protein only)\n")
                f.write(f"{len(protein_atoms)}\n")
                for line in protein_atoms:
                    f.write(line + "\n")
                f.write(box_line + "\n")

            # Generate posre from protein-only GRO (indices will be 1-based within protein)
            cmd_genrestr = [
                "gmx",
                "genrestr",
                "-f",
                str(protein_gro),
                "-o",
                str(posre_file),
                "-fc",
                "1000",
                "1000",
                "1000",
            ]
            subprocess.run(cmd_genrestr, input="0\n", text=True, capture_output=True, check=True, timeout=30)

            # Clean up temporary file
            if protein_gro.exists():
                protein_gro.unlink()

            print(f"  ✓ Generated posre.itp for {leg_name} leg (protein restraints)")

        except subprocess.CalledProcessError as e:
            print(f"  ⚠ Warning: Failed to generate posre.itp: {e}")
            # Create empty file to avoid topology errors
            with open(posre_file, "w") as f:
                f.write("; Position restraint file (fallback)\n")
                f.write("[ position_restraints ]\n")
                f.write("; Failed to generate - empty restraints\n")
        except Exception as e:
            print(f"  ⚠ Warning: Error generating posre.itp: {e}")

    def replace_ligand_with_hybrid(self, source_gro: Path, target_gro: Path, hybrid_gro: Path) -> None:
        """
        Replace ligand atoms in coordinate file with hybrid ligand atoms.

        Parameters
        ----------
        source_gro : Path
            Source GRO file (solv_ions.gro)
        target_gro : Path
            Target GRO file to write
        hybrid_gro : Path
            Hybrid ligand GRO file
        """
        # Read source GRO file
        with open(source_gro, "r") as f:
            source_lines = f.readlines()

        # Read hybrid ligand GRO file
        with open(hybrid_gro, "r") as f:
            hybrid_lines = f.readlines()

        # Parse hybrid ligand atoms
        hybrid_atoms = []
        for line in hybrid_lines[2:]:  # Skip title and atom count
            # Skip empty lines
            if len(line.strip()) == 0:
                continue
            # Check if this is a box line (contains only numbers, dots, and spaces)
            stripped = line.strip()
            if all(c in "0123456789. " for c in stripped):
                break  # Box line found
            # Parse GRO format
            hybrid_atoms.append(line)

        # Find ligand atoms in source file. Accept both the original LIG residue
        # name and the staged hybrid molecule name so an existing scaffold can be
        # refreshed in place after regenerating hybrid.gro.
        ligand_start = -1
        ligand_end = -1
        ligand_resnames = {"LIG", self.molecule_name}

        for i, line in enumerate(source_lines[2:], start=2):  # Skip title and atom count
            if len(line) > 10:
                line_head = line[:15]
                if any(resname in line_head for resname in ligand_resnames):
                    if ligand_start == -1:
                        ligand_start = i
                    ligand_end = i
                elif ligand_start != -1:
                    break

        # Replace ligand atoms with hybrid ligand atoms
        if ligand_start != -1:
            # Build new file content
            new_lines = source_lines[:ligand_start]

            # Get the starting atom index and residue number from the source file
            first_ligand_line = source_lines[ligand_start]
            try:
                res_num_str = first_ligand_line[0:5].strip()
                residue_number = int(res_num_str) if res_num_str else 164
                atom_num_str = first_ligand_line[15:20].strip()
                starting_atom_index = int(atom_num_str) if atom_num_str else 1
            except (ValueError, IndexError):
                residue_number = 164
                starting_atom_index = 1

            # editconf shifts the whole system box; hybrid.gro still has original ff coords
            dx, dy, dz = self._compute_ligand_shift(source_lines, ligand_start, ligand_end, hybrid_atoms)

            # Parse box size from source GRO (last line)
            box_line = source_lines[-1].strip().split()
            box_x = float(box_line[0]) if len(box_line) >= 1 else 10.0
            box_y = float(box_line[1]) if len(box_line) >= 2 else box_x
            box_z = float(box_line[2]) if len(box_line) >= 3 else box_x

            # GRO split: [resnum+resname, atomname, atomidx, x, y, z]
            for i, hybrid_line in enumerate(hybrid_atoms):
                parts = hybrid_line.split()
                if len(parts) >= 6:
                    try:
                        atom_idx = starting_atom_index + i
                        atom_name = parts[1]
                        x = float(parts[3]) + dx
                        y = float(parts[4]) + dy
                        z = float(parts[5]) + dz

                        # Apply PBC wrapping to keep coordinates inside the box [0, box)
                        while x < 0:
                            x += box_x
                        while x >= box_x:
                            x -= box_x
                        while y < 0:
                            y += box_y
                        while y >= box_y:
                            y -= box_y
                        while z < 0:
                            z += box_z
                        while z >= box_z:
                            z -= box_z

                        new_line = f"{residue_number:5d}{self.molecule_name:<5s}{atom_name:>5s}{atom_idx:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
                        new_lines.append(new_line)
                    except (ValueError, IndexError):
                        new_lines.append(hybrid_line)
                else:
                    new_lines.append(hybrid_line)

            new_lines.extend(source_lines[ligand_end + 1 :])

            # Update atom count in memory before writing
            num_atoms_old = int(source_lines[1].strip())
            num_atoms_new = num_atoms_old + len(hybrid_atoms) - (ligand_end - ligand_start + 1)
            new_lines[1] = f"{num_atoms_new}\n"

            with open(target_gro, "w") as f:
                f.writelines(new_lines)
        else:
            # No ligand found, copy as-is when writing to a different path.
            if source_gro.resolve() != target_gro.resolve():
                shutil.copy2(source_gro, target_gro)

    def _compute_ligand_shift(
        self,
        source_lines: list,
        ligand_start: int,
        ligand_end: int,
        hybrid_atoms: list,
    ) -> tuple:
        """
        Compute (dx, dy, dz) to align hybrid.gro coords to the editconf-centered system.

        solv_ions.gro has been shifted by gmx editconf; hybrid.gro still has the original
        forcefield GRO coords. We align via centroid of atoms present in both.
        Returns (0, 0, 0) if no common atoms found.
        """
        # Parse source ligand atoms using fixed-width GRO columns
        source_coords: dict = {}
        for line in source_lines[ligand_start : ligand_end + 1]:
            try:
                name = line[10:15].strip()
                source_coords[name] = (float(line[20:28]), float(line[28:36]), float(line[36:44]))
            except (ValueError, IndexError):
                pass

        # Parse hybrid atoms from split format [resnum+resname, atomname, atomidx, x, y, z]
        hybrid_coords: dict = {}
        for line in hybrid_atoms:
            parts = line.split()
            if len(parts) >= 6:
                try:
                    hybrid_coords[parts[1]] = (float(parts[3]), float(parts[4]), float(parts[5]))
                except (ValueError, IndexError):
                    pass

        common = [n for n in source_coords if n in hybrid_coords]
        if not common:
            print("  ⚠ Warning: no common atoms between source LIG and hybrid.gro; skipping coord shift")
            return (0.0, 0.0, 0.0)

        n = len(common)
        shift = tuple(sum(source_coords[a][i] - hybrid_coords[a][i] for a in common) / n for i in range(3))
        print(f"  ✓ Ligand coord shift: ({shift[0]:+.3f}, {shift[1]:+.3f}, {shift[2]:+.3f}) nm " f"[{n} common atoms]")
        return shift

    def prune_overlapping_solvent(self, conf_gro: Path, topol_top: Path, clash_cutoff_nm: float = 0.23) -> int:
        """Remove solvent molecules whose atoms overlap with hybrid heavy atoms."""
        lines = conf_gro.read_text().splitlines()
        if len(lines) < 4:
            return 0

        atom_lines = lines[2:-1]
        box_line = lines[-1]

        parsed = []
        for raw in atom_lines:
            if len(raw) < 44:
                parsed.append(None)
                continue
            try:
                parsed.append(
                    {
                        "line": raw,
                        "resnr": int(raw[0:5]),
                        "resname": raw[5:10].strip(),
                        "atomname": raw[10:15].strip(),
                        "atomnr": int(raw[15:20]),
                        "x": float(raw[20:28]),
                        "y": float(raw[28:36]),
                        "z": float(raw[36:44]),
                    }
                )
            except ValueError:
                parsed.append(None)

        hybrid_heavy = [
            a for a in parsed if a and a["resname"] == self.molecule_name and not a["atomname"].startswith("H")
        ]
        if not hybrid_heavy:
            return 0

        solvent_residues = set()
        cutoff2 = clash_cutoff_nm * clash_cutoff_nm
        for atom in parsed:
            if not atom or atom["resname"] != "SOL":
                continue
            for lig in hybrid_heavy:
                dx = atom["x"] - lig["x"]
                dy = atom["y"] - lig["y"]
                dz = atom["z"] - lig["z"]
                if dx * dx + dy * dy + dz * dz < cutoff2:
                    solvent_residues.add(atom["resnr"])
                    break

        if not solvent_residues:
            return 0

        kept = [a for a in parsed if not a or a["resnr"] not in solvent_residues]
        new_atom_lines = [a["line"] for a in kept if a]
        new_lines = [lines[0], str(len(new_atom_lines)), *new_atom_lines, box_line]
        conf_gro.write_text("\n".join(new_lines) + "\n")
        self._decrement_molecule_count(topol_top, "SOL", len(solvent_residues))
        return len(solvent_residues)

    def _decrement_molecule_count(self, topol_top: Path, molecule_name: str, removed_count: int) -> None:
        """Decrease a `[ molecules ]` entry after deleting solvent molecules."""
        if removed_count <= 0 or not topol_top.exists():
            return

        lines = topol_top.read_text().splitlines()
        in_molecules = False
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped.lower() == "[ molecules ]":
                in_molecules = True
                continue
            if in_molecules and stripped.startswith("["):
                break
            if not in_molecules or not stripped or stripped.startswith(";"):
                continue
            parts = stripped.split()
            if len(parts) >= 2 and parts[0] == molecule_name:
                new_count = int(parts[1]) - removed_count
                if new_count < 0:
                    raise ValueError(f"Cannot decrement {molecule_name} below zero in {topol_top}")
                lines[idx] = re.sub(rf"^{re.escape(parts[0])}\s+\d+", f"{parts[0]}        {new_count}", line)
                topol_top.write_text("\n".join(lines) + "\n")
                return

    def copy_and_modify_topology(self, source_top: Path, target_top: Path) -> None:
        """
        Copy topology file and modify to use hybrid ligand.

        Replaces ligand ITP include with hybrid ligand ITP and renames LIG to HYB in molecules section.
        Also adds hybrid atomtypes includes for unbound leg.
        """
        content = source_top.read_text()
        hybrid_prefix = self._hybrid_include_prefix(target_top.parent)
        hybrid_include = f'#include "{hybrid_prefix}/hybrid.itp"'
        hybrid_ff_include = f'#include "{hybrid_prefix}/ff_hybrid.itp"'

        # Check if we need to add hybrid atomtypes (for unbound leg)
        hybrid_ff_path = self.output_dir / "common" / "hybrid" / "ff_hybrid.itp"
        needs_hybrid_ff = hybrid_ff_path.exists()
        self._repair_hybrid_forcefield_include(hybrid_ff_path)

        # Replace ligand ITP include with hybrid ligand
        content = re.sub(
            r'#include\s+"[^"]*LIG\.itp"',
            hybrid_include,
            content,
        )

        # Deduplicate hybrid.itp includes (can happen if source topology was already modified)
        if content.count(hybrid_include) > 1:
            lines = content.splitlines()
            seen = False
            deduped = []
            for line in lines:
                if line.strip() == hybrid_include:
                    if not seen:
                        deduped.append(line)
                        seen = True
                else:
                    deduped.append(line)
            content = "\n".join(deduped) + "\n"

        # Add hybrid force field include if needed (after forcefield.itp)
        # Note: Must be after forcefield.itp to ensure [ defaults ] comes before [ atomtypes ]
        if needs_hybrid_ff and hybrid_ff_include not in content:
            # Find forcefield.itp include and insert after it
            lines = content.splitlines()
            for i, line in enumerate(lines):
                if "forcefield.itp" in line or "force_field.itp" in line:
                    # Insert after this line
                    lines.insert(i + 1, hybrid_ff_include)
                    break
            content = "\n".join(lines) + "\n"

        # Strip the ligand-specific [ atomtypes ] block that PRISM inserts for the
        # original ligand. The hybrid package already contributes a merged atomtype
        # include, so keeping both definitions triggers duplicate-atomtype warnings
        # in the bound leg grompp check.
        content = self._remove_top_level_atomtypes_block(content)

        # Ensure the hybrid topology include still exists after topology cleanup.
        if hybrid_include not in content:
            lines = content.splitlines()
            insert_idx = -1
            for i, line in enumerate(lines):
                if line.strip() == hybrid_ff_include:
                    insert_idx = i + 1
                    break
            if insert_idx == -1:
                for i, line in enumerate(lines):
                    if "forcefield.itp" in line or "force_field.itp" in line:
                        insert_idx = i + 1
                        break
            if insert_idx != -1:
                lines.insert(insert_idx, hybrid_include)
                content = "\n".join(lines) + "\n"

        # Some protein+ligand source topologies inline a ligand moleculetype name
        # immediately before the protein moleculetype entry. Once the hybrid ligand
        # is provided via `hybrid.itp`, keeping that stray HYB/LIG entry causes
        # `moleculetype ... is redefined` in grompp.
        content = self._remove_redundant_ligand_moleculetype_entries(content)

        content = self._replace_ligands_in_molecules_section(content)

        # Add posre.itp include if not present
        if "#ifdef POSRES" not in content:
            # Find [ moleculetype ] section for protein
            lines = content.splitlines()
            in_protein_moleculetype = False
            insert_idx = -1

            for i, line in enumerate(lines):
                if "[ moleculetype ]" in line.lower():
                    in_protein_moleculetype = True
                elif in_protein_moleculetype and "[ atoms ]" in line.lower():
                    # Insert before [ atoms ]
                    insert_idx = i
                    break

            if insert_idx > 0:
                posre_block = [
                    "",
                    "; Include Position restraint file",
                    "#ifdef POSRES",
                    '#include "posre.itp"',
                    "#endif",
                    "",
                ]
                for j, posre_line in enumerate(posre_block):
                    lines.insert(insert_idx + j, posre_line)
                content = "\n".join(lines) + "\n"

        target_top.write_text(content)

    def _repair_hybrid_forcefield_include(self, hybrid_ff_path: Path) -> None:
        """Backfill stale ff_hybrid.itp files that forgot the CGenFF bonded supplement include."""
        if not hybrid_ff_path.exists():
            return
        supplement = hybrid_ff_path.with_name("cgenff_bonded_hybrid.itp")
        if not supplement.exists():
            return
        content = hybrid_ff_path.read_text()
        include_line = '#include "cgenff_bonded_hybrid.itp"'
        if include_line in content:
            return
        lines = [line for line in content.splitlines() if line.strip()]
        if not lines:
            lines = ['#include "atomtypes_hybrid.itp"']
        lines.append(include_line)
        hybrid_ff_path.write_text("\n".join(lines) + "\n")

    def _remove_redundant_ligand_moleculetype_entries(self, content: str) -> str:
        """Drop stray HYB/LIG names embedded in a protein `[ moleculetype ]` section."""
        lines = content.splitlines()
        cleaned = []
        i = 0

        while i < len(lines):
            line = lines[i]
            cleaned.append(line)

            if line.strip().lower() != "[ moleculetype ]":
                i += 1
                continue

            i += 1
            section_lines = []
            while i < len(lines):
                current = lines[i]
                stripped = current.strip()
                if stripped.startswith("[") and stripped.lower() != "[ moleculetype ]":
                    break
                section_lines.append(current)
                i += 1

            entry_lines = []
            for offset, section_line in enumerate(section_lines):
                stripped = section_line.strip()
                if not stripped or stripped.startswith(";"):
                    continue
                parts = stripped.split()
                if len(parts) >= 2:
                    entry_lines.append((offset, parts[0]))

            remove_offsets = set()
            if len(entry_lines) > 1:
                non_ligand_entries = [name for _, name in entry_lines if name not in {self.molecule_name, "LIG"}]
                if non_ligand_entries:
                    remove_offsets = {offset for offset, name in entry_lines if name in {self.molecule_name, "LIG"}}

            for offset, section_line in enumerate(section_lines):
                if offset not in remove_offsets:
                    cleaned.append(section_line)

        return "\n".join(cleaned).rstrip() + "\n"

    def _remove_top_level_atomtypes_block(self, content: str) -> str:
        """Remove a standalone top-level [ atomtypes ] block from a copied topology."""
        lines = content.splitlines()
        cleaned = []
        in_atomtypes = False

        for line in lines:
            stripped = line.strip().lower()
            if stripped == "[ atomtypes ]":
                in_atomtypes = True
                continue
            if in_atomtypes and (
                stripped.startswith("#include") or (stripped.startswith("[") and stripped != "[ atomtypes ]")
            ):
                in_atomtypes = False
                cleaned.append(line)
                continue
            if in_atomtypes:
                continue
            cleaned.append(line)

        return "\n".join(cleaned).rstrip() + "\n"

    def _replace_ligands_in_molecules_section(self, content: str) -> str:
        """Replace LIG* entries with a single hybrid ligand entry inside [ molecules ]."""
        lines = content.splitlines()
        molecules_idx = -1
        section_end = len(lines)

        for i, line in enumerate(lines):
            if line.strip().lower() == "[ molecules ]":
                molecules_idx = i
                break

        if molecules_idx == -1:
            return content

        for i in range(molecules_idx + 1, len(lines)):
            stripped = lines[i].strip()
            if stripped.startswith("[") and stripped.lower() != "[ molecules ]":
                section_end = i
                break

        section_lines = lines[molecules_idx + 1 : section_end]
        filtered_section = []
        ligand_count = 0
        inserted_hybrid = False

        for line in section_lines:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                filtered_section.append(line)
                continue

            if re.match(r"^\s*LIG(?:_\d+)?\s+\d+", line):
                ligand_count += 1
                if not inserted_hybrid:
                    parts = stripped.split()
                    count = parts[1] if len(parts) > 1 else "1"
                    filtered_section.append(f"{self.molecule_name:<8} {count}")
                    inserted_hybrid = True
                continue

            filtered_section.append(line)

        if ligand_count > 0 and not inserted_hybrid:
            filtered_section.append(f"{self.molecule_name:<8} 1")

        rebuilt = lines[: molecules_idx + 1] + filtered_section + lines[section_end:]
        return "\n".join(rebuilt).rstrip() + "\n"

    def _hybrid_include_prefix(self, leg_dir: Path) -> str:
        """Return the relative include path from a leg directory to common/hybrid."""
        relative = os.path.relpath(self.output_dir / "common" / "hybrid", leg_dir)
        return relative.replace("\\", "/")

    def sync_unbound_molecules_with_conf(self, topol_path: Path, conf_gro_path: Path) -> None:
        """
        Synchronize [ molecules ] section with actual molecule counts in GRO file.

        For unbound leg, PRISM ligand-only builder may not add solvent/ions to topology.
        This method counts molecules in the GRO file and rebuilds the [ molecules ] section.

        Parameters
        ----------
        topol_path : Path
            Path to topology file
        conf_gro_path : Path
            Path to coordinate file (conf.gro)
        """
        if not conf_gro_path.exists():
            return

        # Count molecules in GRO file
        molecule_counts = {}
        lines = conf_gro_path.read_text().splitlines()
        if len(lines) < 3:
            return

        n_atoms = int(lines[1].strip())
        for i in range(2, 2 + n_atoms):
            if i >= len(lines):
                break
            line = lines[i]
            if len(line) < 10:
                continue
            resname = line[5:10].strip()
            molecule_counts[resname] = molecule_counts.get(resname, 0) + 1

        # Map residue names to molecule names
        resname_to_molname = {
            self.molecule_name: self.molecule_name,
            "LIG": self.molecule_name,
            "SOL": "SOL",
            "NA": "NA",
            "CL": "CL",
            "K": "K",
            "MG": "MG",
            "CA": "CA",
        }

        # Count molecules (not atoms)
        mol_counts = {}
        for resname, atom_count in molecule_counts.items():
            molname = resname_to_molname.get(resname, resname)
            if molname == "SOL":
                # Water has 3 atoms per molecule
                mol_counts[molname] = atom_count // 3
            elif molname in ["NA", "CL", "K", "MG", "CA"]:
                # Ions are single atoms
                mol_counts[molname] = atom_count
            elif molname == self.molecule_name:
                # Hybrid ligand - count as 1 molecule
                mol_counts[molname] = 1

        # Rebuild [ molecules ] section
        content = topol_path.read_text()
        lines = content.splitlines()

        # Find [ molecules ] section
        molecules_idx = -1
        for i, line in enumerate(lines):
            if line.strip().lower() == "[ molecules ]":
                molecules_idx = i
                break

        if molecules_idx == -1:
            return

        # Remove old [ molecules ] content
        new_lines = lines[: molecules_idx + 1]
        new_lines.append("; Compound        #mols")

        # Add molecules in standard order
        for molname in [self.molecule_name, "SOL", "NA", "CL", "K", "MG", "CA"]:
            if molname in mol_counts and mol_counts[molname] > 0:
                new_lines.append(f"{molname:<16s} {mol_counts[molname]}")

        topol_path.write_text("\n".join(new_lines) + "\n")

    def write_complete_topology(
        self,
        output_path: Path,
        forcefield: str,
        water_model: str,
        molecules: list,
        include_protein: bool = True,
    ) -> None:
        """
        Write complete topology file for FEP leg.

        Parameters
        ----------
        output_path : Path
            Output topology file path
        forcefield : str
            Force field name
        water_model : str
            Water model name
        molecules : list
            List of (molecule_name, count) tuples
        include_protein : bool
            Whether to include protein topology
        """
        hybrid_prefix = self._hybrid_include_prefix(output_path.parent)
        lines = [
            "; Topology for FEP calculation",
            "",
            f'#include "{forcefield}.ff/forcefield.itp"',
            f'#include "{hybrid_prefix}/ff_hybrid.itp"',
            "",
        ]

        if include_protein:
            lines.extend(
                [
                    "; Include protein topology",
                    '#include "protein.itp"',
                    "",
                    "; Include Position restraint file",
                    "#ifdef POSRES",
                    '#include "posre.itp"',
                    "#endif",
                    "",
                ]
            )

        lines.extend(
            [
                "; Include hybrid ligand topology",
                f'#include "{hybrid_prefix}/hybrid.itp"',
                "",
                f"; Include water topology",
                f'#include "{forcefield}.ff/{water_model}.itp"',
                "",
                "#ifdef POSRES_WATER",
                "; Position restraint for each water oxygen",
                "[ position_restraints ]",
                ";  i funct       fcx        fcy        fcz",
                "   1    1       1000       1000       1000",
                "#endif",
                "",
                "; Include topology for ions",
                f'#include "{forcefield}.ff/ions.itp"',
                "",
                "[ system ]",
                "FEP system",
                "",
                "[ molecules ]",
                "; Compound        #mols",
            ]
        )

        for mol_name, count in molecules:
            lines.append(f"{mol_name:<16s} {count}")

        output_path.write_text("\n".join(lines) + "\n")

    def extract_molecules_from_topology(self, topol_path: Path) -> list:
        """
        Extract [ molecules ] section from topology file.

        Returns
        -------
        list
            List of (molecule_name, count) tuples
        """
        if not topol_path.exists():
            return []

        content = topol_path.read_text()
        molecules = []
        in_molecules = False

        for line in content.splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ molecules ]":
                in_molecules = True
                continue
            if in_molecules:
                if stripped.startswith("["):
                    break
                if not stripped or stripped.startswith(";"):
                    continue
                parts = stripped.split()
                if len(parts) >= 2:
                    try:
                        molecules.append((parts[0], int(parts[1])))
                    except ValueError:
                        continue

        return molecules

    def write_complex_seed(self, receptor_pdb: Path, ligand_pdb: Path, output_pdb: Path) -> None:
        """
        Combine receptor and ligand PDB files into complex seed.

        Parameters
        ----------
        receptor_pdb : Path
            Receptor PDB file
        ligand_pdb : Path
            Ligand PDB file
        output_pdb : Path
            Output complex PDB file
        """
        receptor_records = self.extract_pdb_records(receptor_pdb)
        ligand_records = self.extract_pdb_records(ligand_pdb)

        merged = ["REMARK   PRISM-FEP bound-leg seed complex"]
        serial = 1
        for record in receptor_records + ligand_records:
            merged.append(record[:6] + f"{serial:5d}" + record[11:])
            serial += 1
        merged.append("END")
        output_pdb.write_text("\n".join(merged) + "\n")

    def extract_pdb_records(self, pdb_path: Path) -> list[str]:
        """Extract ATOM/HETATM records from PDB file."""
        records = []
        for line in pdb_path.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                records.append(line)
        return records

    def prepare_protein_for_system_builder(self, protein_path: Path, leg_dir: Path) -> Path:
        """
        Prepare protein PDB for PRISM system builder.

        For unbound leg, creates a dummy protein file.

        Parameters
        ----------
        protein_path : Path
            Original protein PDB path (None for unbound)
        leg_dir : Path
            Leg directory

        Returns
        -------
        Path
            Path to prepared protein file
        """
        if protein_path is None:
            return self.create_dummy_protein(leg_dir)
        return protein_path

    def create_dummy_protein(self, leg_dir: Path) -> Path:
        """Create dummy protein PDB for unbound leg."""
        dummy_pdb = leg_dir / "dummy_protein.pdb"
        with open(dummy_pdb, "w") as f:
            f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")
            f.write("END\n")
        return dummy_pdb

    def get_forcefield_indices(self, config: dict) -> tuple:
        """Get force field and water model indices from config."""
        ff_idx = config.get("forcefield_index", 14)  # amber14sb
        water_idx = config.get("water_model_index", 1)  # tip3p
        return ff_idx, water_idx

    def get_forcefield_info(self, _config: dict, ff_idx: int) -> dict:
        """Get force field information."""
        return {"name": "amber14sb", "index": ff_idx}

    def get_water_info(self, _conf, water_idx: int) -> dict:
        """Get water model information."""
        return {"name": "tip3p", "index": water_idx}
