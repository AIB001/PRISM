#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
I/O utilities for protein cleaning.
"""

import os
from typing import List, Tuple


class IOMixin:
    """Mixin for PDB file I/O operations."""

    def _read_pdb(self, pdb_file: str) -> Tuple[List[str], List[str]]:
        """
        Read PDB file and separate protein and heteroatom records.

        Returns
        -------
        tuple
            (protein_lines, hetatm_lines)
        """
        protein_lines = []
        hetatm_lines = []

        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    protein_lines.append(line)
                    self.stats["total_atoms"] += 1
                    self.stats["protein_atoms"] += 1
                elif line.startswith("HETATM"):
                    hetatm_lines.append(line)
                    self.stats["total_atoms"] += 1

        return protein_lines, hetatm_lines

    def _write_pdb(self, output_pdb: str, protein_lines: List[str], hetatm_lines: List[str], input_pdb: str):
        """
        Write cleaned PDB file with proper chain organization for pdb2gmx.

        Parameters
        ----------
        output_pdb : str
            Output file path
        protein_lines : list of str
            Protein ATOM records
        hetatm_lines : list of str
            Kept HETATM records
        input_pdb : str
            Original input file (for header information)
        """
        # Convert metal HETATM to ATOM for pdb2gmx compatibility
        metal_atom_lines, other_hetatm_lines = self._convert_metals_to_atom(hetatm_lines)

        # Group protein and metal atoms by chain for proper output
        chain_atoms = self._group_by_chain(protein_lines, metal_atom_lines)

        with open(output_pdb, "w") as out:
            # Write header
            out.write(f"REMARK   Cleaned by PRISM ProteinCleaner\n")
            out.write(f"REMARK   Original file: {os.path.basename(input_pdb)}\n")
            out.write(f"REMARK   Ion mode: {self.ion_mode}\n")
            out.write(f"REMARK   Distance cutoff: {self.distance_cutoff} A\n")
            out.write(f"REMARK   Keep crystal water: {self.keep_crystal_water}\n")
            out.write(f"REMARK   Remove crystallization artifacts: {self.remove_artifacts}\n")
            if self.stats.get("nterm_met_removed", 0) > 0:
                chains = ",".join(self.stats.get("nterm_met_chains", []))
                out.write(f"REMARK   Removed N-terminal MET on chain(s): {chains}\n")
            if metal_atom_lines:
                out.write(f"REMARK   Converted {len(metal_atom_lines)} metal HETATM to ATOM for pdb2gmx\n")
                out.write(f"REMARK   Metals integrated into their respective protein chains\n")

            # Write chains with protein atoms followed by their metals
            for chain_id in sorted(chain_atoms.keys()):
                protein_atoms, metal_atoms = chain_atoms[chain_id]

                # Write protein atoms
                for line in protein_atoms:
                    out.write(line)

                # Write metal atoms immediately after protein atoms in the same chain
                # This keeps them together for pdb2gmx
                for line in metal_atoms:
                    out.write(line)

                # TER record after each chain
                out.write("TER\n")

            # Write other HETATM records (cofactors, modified residues, etc.)
            for line in other_hetatm_lines:
                out.write(line)

            # Final END record
            out.write("END\n")
