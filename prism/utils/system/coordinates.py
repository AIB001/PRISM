#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coordinate file processing utilities for SystemBuilder.
"""

from pathlib import Path

try:
    from ..colors import print_success, number
except ImportError:
    from prism.utils.colors import print_success, number


class CoordinateProcessorMixin:
    """Mixin for coordinate-related operations."""

    def _combine_protein_multi_ligands(self, protein_gro: str, lig_ff_dirs: list) -> str:
        """Combines protein and multiple ligand GRO files into a complex."""
        print(f"\n  Step 4: Combining protein with {number(len(lig_ff_dirs))} ligand(s)...")
        complex_gro = self.model_dir / "pro_lig.gro"

        if complex_gro.exists() and not self.overwrite:
            print("Using existing complex.gro file.")
            return str(complex_gro)

        # Read protein coordinates
        with open(protein_gro, "r") as f_prot:
            prot_lines = f_prot.readlines()

        # Read all ligand coordinates
        all_lig_lines = []
        total_lig_atoms = 0

        for idx, lig_ff_dir in enumerate(lig_ff_dirs, 1):
            lig_gro_path = Path(lig_ff_dir) / "LIG.gro"
            if not lig_gro_path.exists():
                raise FileNotFoundError(f"Ligand {idx} GRO file not found: {lig_gro_path}")

            with open(lig_gro_path, "r") as f_lig:
                lig_lines = f_lig.readlines()

            # Extract coordinates (skip header, atom count, and box vectors)
            lig_coords = lig_lines[2:-1]

            # Renumber residue numbers and atom numbers for each ligand
            renumbered_coords = []
            for line in lig_coords:
                if len(line) >= 20:
                    # Update residue number (columns 0-4): each ligand gets its own residue number
                    new_res_num = idx
                    # Update atom number (columns 15-19): sequential numbering
                    current_atom_num = int(line[15:20].strip())
                    new_atom_num = total_lig_atoms + current_atom_num
                    # Reconstruct the line
                    new_line = f"{new_res_num:5d}" + line[5:15] + f"{new_atom_num:5d}" + line[20:]
                    renumbered_coords.append(new_line)
                else:
                    renumbered_coords.append(line)

            total_lig_atoms += len(lig_coords)
            all_lig_lines.extend(renumbered_coords)
            print(f"  Ligand {number(idx)}: {number(len(lig_coords))} atoms")

        # Calculate total atoms
        prot_atom_count = int(prot_lines[1].strip())
        total_atoms = prot_atom_count + total_lig_atoms

        # Write combined GRO file
        with open(complex_gro, "w") as f_out:
            f_out.write("Complex Protein-Ligands\n")
            f_out.write(f"{total_atoms:5d}\n")

            # Write all ligand coordinates first (to match topology molecule order)
            f_out.writelines(all_lig_lines)

            # Write protein coordinates, renumbering both residue and atom numbers
            prot_coords = prot_lines[2:-1]
            num_ligands = len(lig_ff_dirs)

            for line in prot_coords:
                if len(line) >= 20:
                    # Update residue number (columns 0-4): shift by number of ligands
                    current_res_num = int(line[:5].strip())
                    new_res_num = current_res_num + num_ligands
                    # Update atom number (columns 15-19): shift by total ligand atoms
                    current_atom_num = int(line[15:20].strip())
                    new_atom_num = current_atom_num + total_lig_atoms
                    # Reconstruct the line
                    new_line = f"{new_res_num:5d}" + line[5:15] + f"{new_atom_num:5d}" + line[20:]
                    f_out.write(new_line)
                else:
                    f_out.write(line)

            # Write box vectors from protein file
            f_out.write(prot_lines[-1])

        print_success(
            f"Combined system: {number(prot_atom_count)} protein atoms + {number(total_lig_atoms)} ligand atoms = {number(total_atoms)} total"
        )
        return str(complex_gro)

    def _create_box(self, complex_gro: str) -> str:
        """Creates the simulation box."""
        print("\n  Step 5: Creating simulation box...")
        boxed_gro = self.model_dir / "pro_lig_newbox.gro"

        if boxed_gro.exists() and not self.overwrite:
            print("Using existing boxed.gro file.")
            return str(boxed_gro)

        box_cfg = self.config["box"]

        if self.pmf_mode:
            # PMF mode: Create rectangular box, then extend ONLY Z direction
            print(f"  PMF Mode: Creating rectangular box with Z-axis extension for pulling")
            print(f"  Z extension: {self.box_extension[2]:.1f} nm")

            # First create rectangular (triclinic) box - NOT cubic!
            # This ensures X and Y dimensions fit the protein tightly
            temp_boxed = self.model_dir / "temp_boxed.gro"
            command = [
                self.gmx_command,
                "editconf",
                "-f",
                complex_gro,
                "-o",
                str(temp_boxed),
                "-bt",
                "triclinic",  # Use triclinic (rectangular) for PMF, not cubic
                "-d",
                str(box_cfg["distance"]),
            ]
            if box_cfg.get("center", True):
                command.append("-c")

            self._run_command(command, str(self.model_dir))

            # Read box vectors and extend ONLY Z direction (keep X, Y unchanged)
            self._extend_box_z(str(temp_boxed), str(boxed_gro))

            # Clean up temp file
            if temp_boxed.exists():
                temp_boxed.unlink()
        else:
            # Normal mode: Standard box creation
            command = [
                self.gmx_command,
                "editconf",
                "-f",
                complex_gro,
                "-o",
                str(boxed_gro),
                "-bt",
                box_cfg["shape"],
                "-d",
                str(box_cfg["distance"]),
            ]
            if box_cfg.get("center", True):
                command.append("-c")

            self._run_command(command, str(self.model_dir))

        return str(boxed_gro)

    def _extend_box_z(self, input_gro: str, output_gro: str):
        """
        Extend box in Z direction for PMF calculations.

        For SMD pulling:
        - After alignment, ligand is above protein in +Z direction
        - Ligand will be pulled further in +Z direction
        - Extra space is added at the +Z end of the box (no atom translation)
        - This ensures the pulling space is in the correct direction
        """
        with open(input_gro, "r") as f:
            lines = f.readlines()

        if len(lines) < 3:
            raise ValueError(f"Invalid GRO file: {input_gro}")

        # Parse box vectors from last line
        box_line = lines[-1].strip()
        box_parts = box_line.split()

        if len(box_parts) >= 3:
            # Box vectors: v1x v2y v3z [v1y v1z v2x v2z v3x v3y]
            v1x = float(box_parts[0])
            v2y = float(box_parts[1])
            v3z = float(box_parts[2])

            # For PMF mode: ONLY extend Z direction, keep X and Y unchanged
            new_v1x = v1x  # X unchanged
            new_v2y = v2y  # Y unchanged
            new_v3z = v3z + self.box_extension[2]  # Only Z extended

            print(f"    Original box: {v1x:.3f} x {v2y:.3f} x {v3z:.3f} nm")
            print(
                f"    Extended box: {new_v1x:.3f} x {new_v2y:.3f} x {new_v3z:.3f} nm (Z +{self.box_extension[2]:.1f} nm)"
            )
            print(f"    Extra space added at +Z end for pulling direction")

            # For SMD: Do NOT translate atoms - keep them at original positions
            # The extra space will be at the +Z end where the ligand will be pulled
            new_lines = []
            new_lines.append(lines[0])  # Title
            new_lines.append(lines[1])  # Atom count

            # Copy all atom lines without modification
            for line in lines[2:-1]:
                new_lines.append(line)

            # Write new box vectors (only Z changed)
            if len(box_parts) > 3:
                # Preserve off-diagonal elements if present
                new_box_line = f"   {new_v1x:.5f}   {new_v2y:.5f}   {new_v3z:.5f}"
                for i in range(3, len(box_parts)):
                    new_box_line += f"   {box_parts[i]}"
                new_box_line += "\n"
            else:
                new_box_line = f"   {new_v1x:.5f}   {new_v2y:.5f}   {new_v3z:.5f}\n"

            new_lines.append(new_box_line)

            # Write output file
            with open(output_gro, "w") as f:
                f.writelines(new_lines)

            print(f"    Box extended successfully")
        else:
            raise ValueError(f"Cannot parse box vectors from: {box_line}")
