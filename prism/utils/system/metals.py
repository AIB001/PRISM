#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Metal ion processing utilities for SystemBuilder.
"""

from prism.forcefield.common.units import angstrom_to_nm


class MetalProcessorMixin:
    """Mixin for metal ion processing operations."""

    def _extract_metals_from_pdb(self, pdb_file: str) -> list:
        """
        Extract metal ion information from PDB file.

        IMPORTANT: With the new cleaner.py implementation, metals are converted
        from HETATM to ATOM records, so we need to check BOTH record types.

        However, we must EXCLUDE protein residues like CYM, HID, HIE which are
        also ATOM records but are amino acids, not metals.

        Returns:
            List of dicts with keys: residue_name, atom_name, chain, resnum, coords
        """
        metals = []
        # True metal ions (must match cleaner.py definitions)
        metal_names = {
            "ZN",
            "MG",
            "CA",
            "FE",
            "CU",
            "MN",
            "CO",
            "NI",
            "CD",
            "HG",
            "ZN2",
            "MG2",
            "CA2",
            "FE2",
            "FE3",
            "CU2",
            "MN2",
        }

        # Amino acid residues that should NEVER be treated as metals
        # These include all standard + special protonation states
        protein_residues = {
            # Standard amino acids
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
            # Special protonation states
            "HID",
            "HIE",
            "HIP",
            "HSD",
            "HSE",
            "HSP",
            "HSPM",  # Histidine variants
            "CYM",
            "CYX",  # Cysteine variants (deprotonated, disulfide)
            "LYN",  # Lysine neutral
            "ASH",
            "GLH",  # Protonated acidic residues
            # Terminal variants
            "ACE",
            "NME",
            "NH2",
        }

        with open(pdb_file, "r") as f:
            for line in f:
                # Check both ATOM and HETATM (metals may be converted to ATOM by cleaner)
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    residue_name = line[17:20].strip().upper()
                    atom_name = line[12:16].strip().upper()

                    # CRITICAL: Skip if this is a protein residue
                    if residue_name in protein_residues:
                        continue

                    # Check if this is a metal (by residue name or atom name)
                    if residue_name in metal_names or atom_name in metal_names:
                        chain = line[21].strip()
                        resnum = line[22:26].strip()
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        metals.append(
                            {
                                "residue_name": residue_name,
                                "atom_name": atom_name,
                                "chain": chain,
                                "resnum": resnum,
                                "coords": (x, y, z),
                            }
                        )

        return metals

    def _add_metals_to_topology(self, pdb_file: str, gro_file: str, topol_file: str):
        """
        Add metal ions to topology and coordinates after pdb2gmx.

        pdb2gmx only processes ATOM records (protein), so we need to manually
        add HETATM records (metals) to both topology and coordinates.
        """
        metals = self._extract_metals_from_pdb(pdb_file)

        if not metals:
            return  # No metals to add

        print(f"\nAdding {len(metals)} metal ion(s) to topology and coordinates...")

        # Add metals to topology
        with open(topol_file, "r") as f:
            topol_lines = f.readlines()

        # Find [ molecules ] section
        molecules_idx = -1
        for i, line in enumerate(topol_lines):
            if line.strip() == "[ molecules ]":
                molecules_idx = i
                break

        if molecules_idx == -1:
            print("Warning: Could not find [ molecules ] section in topology")
            return

        # Count metals by residue name
        from collections import Counter

        metal_counts = Counter(m["residue_name"] for m in metals)

        # Insert metal molecule entries after [ molecules ] header
        insert_idx = molecules_idx + 1
        # Skip comment lines
        while insert_idx < len(topol_lines) and topol_lines[insert_idx].strip().startswith(";"):
            insert_idx += 1

        for metal_name, count in metal_counts.items():
            topol_lines.insert(insert_idx, f"{metal_name:<20} {count}\n")
            insert_idx += 1
            print(f"  Added {metal_name}: {count}")

        with open(topol_file, "w") as f:
            f.writelines(topol_lines)

        # Add metals to GRO file
        with open(gro_file, "r") as f:
            gro_lines = f.readlines()

        # Parse current atom count
        atom_count = int(gro_lines[1].strip())

        # Convert metal coordinates to GRO format and append before box line
        metal_gro_lines = []
        for i, metal in enumerate(metals, start=1):
            # GRO format: resnr resname atomname atomnr x y z
            resnum = metal["resnum"]
            resname = metal["residue_name"]
            atomname = metal["atom_name"]
            atom_num = atom_count + i
            x, y, z = metal["coords"]
            # Convert Angstroms to nm
            x_nm, y_nm, z_nm = angstrom_to_nm(x), angstrom_to_nm(y), angstrom_to_nm(z)

            # GRO format line
            line = f"{int(resnum):5d}{resname:5s}{atomname:>5s}{atom_num:5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n"
            metal_gro_lines.append(line)

        # Update atom count
        gro_lines[1] = f"{atom_count + len(metals):5d}\n"

        # Insert metal lines before box line (last line)
        gro_lines = gro_lines[:-1] + metal_gro_lines + [gro_lines[-1]]

        with open(gro_file, "w") as f:
            f.writelines(gro_lines)

        print(f"Successfully added {len(metals)} metal ion(s) to system")
