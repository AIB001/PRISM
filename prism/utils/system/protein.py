#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein processing utilities for SystemBuilder.
"""

import shutil
from pathlib import Path


class ProteinProcessorMixin:
    """Mixin for protein-related operations."""

    def _fix_protein(self, cleaned_protein: str) -> str:
        """Fix missing atoms in protein using pdbfixer."""
        print("\n  Step 1: Fixing protein with pdbfixer...")
        fixed_pdb = self.model_dir / "fixed_clean_protein.pdb"

        if fixed_pdb.exists() and not self.overwrite:
            print(f"Fixed protein PDB already exists at {fixed_pdb}, skipping pdbfixer.")
            return str(fixed_pdb)

        # Check if pdbfixer is available
        if not shutil.which("pdbfixer"):
            print("Warning: pdbfixer command not found. Skipping protein fixing step.")
            print("RECOMMENDED INSTALLATION:")
            print("  mamba install -c conda-forge pdbfixer")
            print("  # OR: conda install -c conda-forge pdbfixer")
            shutil.copy(cleaned_protein, fixed_pdb)
        else:
            # Run pdbfixer to fix missing atoms
            # Note: pdbfixer with --keep-heterogens=all will move metals to chain 'M'
            # This is OK - we extract them from chain 'M' after pdb2gmx
            command = [
                "pdbfixer",
                cleaned_protein,
                "--add-residues",
                "--output",
                str(fixed_pdb),
                "--add-atoms=heavy",
                "--keep-heterogens=all",
                f"--ph={self.config['simulation']['pH']}",
            ]
            self._run_command(command, str(self.model_dir.parent))

        # Convert AMBER histidine names to GROMACS-compatible names
        print("Converting AMBER histidine names to GROMACS format...")
        self._convert_amber_histidine(str(fixed_pdb))

        # Fix terminal atoms for GROMACS compatibility (OXT → O or OC1/OC2)
        print("Fixing terminal atoms for GROMACS compatibility...")
        from prism.utils.cleaner import fix_terminal_atoms

        # Get force field name if available
        ff_name = getattr(self, "ff_info", {}).get("name") if hasattr(self, "ff_info") and self.ff_info else None
        fix_terminal_atoms(str(fixed_pdb), str(fixed_pdb), force_field=ff_name, verbose=True)

        return str(fixed_pdb)

    def _convert_amber_histidine(self, pdb_file: str) -> None:
        """Convert AMBER histidine residue names to GROMACS-compatible names.

        AMBER force fields use: HSD, HSE, HSP
        Different GROMACS force fields use different naming:
        - AMBER force fields: HID, HIE, HIP
        - OPLS force fields: HISD, HISE, HISH
        - CHARMM force fields: HIS (generic name only - pdb2gmx handles protonation)

        This conversion is necessary for pdb2gmx to recognize the residues.

        Parameters
        ----------
        pdb_file : str
            Path to PDB file to modify
        """
        # Determine histidine mapping based on force field
        ff_info = getattr(self, "ff_info", None)
        ff_name = ff_info.get("name", "").lower() if ff_info else ""

        if "opls" in ff_name:
            # OPLS-AA force field uses HISD/HISE/HISH
            his_mapping = {
                "HSD": "HISD",  # δ-protonated histidine
                "HSE": "HISE",  # ε-protonated histidine
                "HSP": "HISH",  # doubly protonated histidine
            }
        elif "charmm" in ff_name:
            # CHARMM force fields use generic HIS only
            # pdb2gmx will handle protonation state automatically
            # Do NOT convert - keep HIS as HIS to avoid "Residue HIE not found" error
            print("  CHARMM force field detected - using generic HIS naming (pdb2gmx will handle protonation)")
            return
        else:
            # AMBER force fields use HID/HIE/HIP
            his_mapping = {
                "HSD": "HID",  # δ-protonated histidine
                "HSE": "HIE",  # ε-protonated histidine
                "HSP": "HIP",  # doubly protonated histidine
            }

        with open(pdb_file, "r") as f:
            lines = f.readlines()

        # Track conversion statistics
        conversion_stats = {name: 0 for name in his_mapping.keys()}
        new_lines = []

        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # PDB format: positions 17-20 (index 17:21) contain residue name
                # But we need to be careful with whitespace
                if len(line) >= 20:
                    resname = line[17:20].strip()  # Get residue name (columns 17-20)
                    if resname in his_mapping:
                        # Replace with GROMACS name
                        # Preserve formatting (right-justified in 3 columns)
                        new_resname = his_mapping[resname]
                        line = line[:17] + f"{new_resname:>3s}" + line[20:]
                        conversion_stats[resname] += 1
            new_lines.append(line)

        # Write back
        with open(pdb_file, "w") as f:
            f.writelines(new_lines)

        # Report conversions
        total_converted = sum(conversion_stats.values())
        if total_converted > 0:
            print(f"  Converted {total_converted} histidine residue(s):")
            for amber, gmx in his_mapping.items():
                if conversion_stats[amber] > 0:
                    print(f"    {amber} → {gmx}: {conversion_stats[amber]} residue(s)")
        else:
            print("  No AMBER histidine residues found (HSD/HSE/HSP)")

    def _apply_propka_renaming(self, pdb_file: str):
        """Apply PROPKA pKa-based histidine renaming to fixed PDB.

        This must run AFTER pdbfixer (which may revert HIS names back to HIS)
        and BEFORE _generate_topology / _rename_his_for_cmap.
        """
        from ..protonation import PropkaProtonator

        ph = self.config.get("protonation", {}).get("ph", self.config.get("simulation", {}).get("pH", 7.0))

        print(f"\n  Applying PROPKA pKa-based residue renaming (pH {ph})...")
        protonator = PropkaProtonator(ph=ph, verbose=True)
        stats = protonator.optimize_protein_protonation(pdb_file, pdb_file, ff_info=self.ff_info)

        renamed = stats.get("renamed", {})
        if renamed:
            print(f"  PROPKA: Renamed {len(renamed)} residue(s):")
            for (chain, resnum, resname), new_name in renamed.items():
                chain_label = chain or "-"
                print(f"    Chain {chain_label} Residue {resnum}: {resname} -> {new_name}")
        else:
            print("  PROPKA: No residues needed renaming")

    def _rename_his_for_cmap(self, pdb_file: str, ff_info: dict) -> bool:
        """Rename HIS → HISE in PDB for force fields with CMAP corrections.

        GROMACS pdb2gmx keeps the original PDB residue name in the topology
        even after selecting a specific protonation state (e.g. HIE).  For
        force fields that use CMAP corrections (like amber19sb), the CMAP
        type lookup is residue-name-dependent (e.g. XC-HIE, N-HIE), so
        leaving the name as HIS causes 'Unknown cmap torsion' errors.

        Renaming HIS → HISE in the PDB triggers the r2b mapping (HISE → HIE),
        so pdb2gmx writes the correct residue name and skips the interactive
        histidine protonation prompt.

        Returns True if any residues were renamed.
        """
        if not ff_info or "path" not in ff_info:
            return False

        ff_dir = Path(ff_info["path"])
        if not (ff_dir / "cmap.itp").exists():
            return False

        from ..protonation import default_histidine_name, get_ff_residue_names

        ff_residues = get_ff_residue_names(ff_info)
        if not ff_residues:
            return False

        target_name = default_histidine_name(ff_info)
        # Only auto-rename for Amber-style histidines (HIE exists in rtp).
        # For other force fields (e.g., CHARMM), keep interactive selection.
        if target_name != "HIE":
            return False

        # Read and rename HIS → force-field default (neutral histidine)
        # PDB format (0-indexed): positions 17-19 = residue name (3 chars).
        # Replace "HIS" with the target name so pdb2gmx uses the correct rtp block
        # and writes the residue name in the topology (avoids CMAP lookup errors).
        with open(pdb_file, "r") as f:
            lines = f.readlines()

        renamed_residues = set()
        new_lines = []
        for line in lines:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20] == "HIS" and line[20] == " ":
                resnum = line[22:26].strip()
                chain = line[21].strip()
                res_id = f"{chain}:{resnum}" if chain else resnum
                renamed_residues.add(res_id)
                line = line[:17] + f"{target_name:3s}" + line[20:]
            new_lines.append(line)

        if renamed_residues:
            with open(pdb_file, "w") as f:
                f.writelines(new_lines)
            print(f"  Renamed {len(renamed_residues)} HIS residue(s) → {target_name} for CMAP compatibility")
            return True

        return False
