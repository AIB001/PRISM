#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP leg directory writer and topology manager.

This module handles writing bound/unbound leg directories, copying PRISM systems,
modifying topologies, and managing coordinate files.
"""

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

        # Copy and modify topology to use hybrid ligand
        self.copy_and_modify_topology(topol_top, leg_dir / "topol.top")

        # Unbound leg topology from ligand-only builder may miss solvent/ion counts;
        # rebuild [ molecules ] from generated coordinate file to keep topology consistent.
        if leg_name == "unbound":
            self.sync_unbound_molecules_with_conf(leg_dir / "topol.top", leg_dir / "input" / "conf.gro")

        print(f"  ✓ Copied and modified topol.top → {leg_name}/topol.top")

        # Generate position restraint file for protein
        self.generate_posre_file(leg_dir, leg_name)

        print(f"  ✓ {leg_name.capitalize()} leg system ready")
        print(f"{'='*70}\n")

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

        # Find ligand atoms in source file (residue name LIG)
        ligand_start = -1
        ligand_end = -1

        for i, line in enumerate(source_lines[2:], start=2):  # Skip title and atom count
            # Check if this is a ligand line (residue name LIG)
            if len(line) > 10:
                # Try to find "LIG" in the line
                if "LIG" in line[:15]:  # Check first 15 characters for "LIG"
                    if ligand_start == -1:
                        ligand_start = i
                    ligand_end = i
                elif ligand_start != -1:
                    # We've moved past the ligand
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

            # Add hybrid ligand atoms with updated indices
            for i, hybrid_line in enumerate(hybrid_atoms):
                parts = hybrid_line.split()
                if len(parts) >= 7:
                    try:
                        atom_idx = starting_atom_index + i
                        atom_name = parts[2]
                        x = parts[4] if len(parts) > 4 else "0.000"
                        y = parts[5] if len(parts) > 5 else "0.000"
                        z = parts[6] if len(parts) > 6 else "0.000"

                        new_line = f"{residue_number:5d}{parts[1]:<5.5s}{atom_name:>5.5s}{atom_idx:5d}{x:>8.8s}{y:>8.8s}{z:>8.8s}\n"
                        new_lines.append(new_line)
                    except (ValueError, IndexError):
                        new_lines.append(hybrid_line)
                else:
                    new_lines.append(hybrid_line)

            new_lines.extend(source_lines[ligand_end + 1 :])

            # Write to target file
            with open(target_gro, "w") as f:
                f.writelines(new_lines)

            # Update atom count in line 2
            num_atoms_old = int(source_lines[1].strip())
            num_atoms_diff = len(hybrid_atoms) - (ligand_end - ligand_start + 1)
            num_atoms_new = num_atoms_old + num_atoms_diff

            # Read back the file to update the count
            with open(target_gro, "r") as f:
                target_content = f.readlines()

            target_content[1] = f"{num_atoms_new}\n"

            with open(target_gro, "w") as f:
                f.writelines(target_content)
        else:
            # No ligand found, copy as-is
            shutil.copy2(source_gro, target_gro)

    def copy_and_modify_topology(self, source_top: Path, target_top: Path) -> None:
        """
        Copy topology file and modify to use hybrid ligand.

        Replaces ligand ITP include with hybrid ligand ITP and renames LIG to HYB in molecules section.
        Also adds hybrid atomtypes includes for unbound leg.
        """
        content = source_top.read_text()

        # Check if we need to add hybrid atomtypes (for unbound leg)
        hybrid_ff_path = self.output_dir / "common" / "hybrid" / "ff_hybrid.itp"
        needs_hybrid_ff = hybrid_ff_path.exists()

        # Replace ligand ITP include with hybrid ligand
        content = re.sub(
            r'#include\s+"[^"]*LIG\.itp"',
            '#include "../common/hybrid/hybrid.itp"',
            content,
        )

        # Deduplicate hybrid.itp includes (can happen if source topology was already modified)
        hybrid_include = '#include "../common/hybrid/hybrid.itp"'
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
        if needs_hybrid_ff and '#include "../common/hybrid/ff_hybrid.itp"' not in content:
            # Find forcefield.itp include and insert after it
            lines = content.splitlines()
            for i, line in enumerate(lines):
                if "forcefield.itp" in line or "force_field.itp" in line:
                    # Insert after this line
                    lines.insert(i + 1, '#include "../common/hybrid/ff_hybrid.itp"')
                    break
            content = "\n".join(lines) + "\n"

        # Replace LIG with HYB in molecules section
        content = re.sub(
            r"^(\s*)LIG(\s+\d+.*?)$",
            rf"\1{self.molecule_name}\2",
            content,
            flags=re.MULTILINE,
        )

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
        lines = [
            "; Topology for FEP calculation",
            "",
            f'#include "{forcefield}.ff/forcefield.itp"',
            '#include "../common/hybrid/ff_hybrid.itp"',
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
                '#include "../common/hybrid/hybrid.itp"',
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
