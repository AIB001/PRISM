#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein processing utilities for protein cleaning.
"""

import numpy as np
from typing import List


class ProteinMixin:
    """Mixin for protein processing operations."""

    def _handle_nterm_met(self, protein_lines: List[str]) -> List[str]:
        """Optionally remove N-terminal MET residues per chain."""
        mode = (self.drop_nterm_met or "keep").lower()
        if mode not in {"keep", "drop", "auto"}:
            mode = "keep"

        # Identify first residue per chain
        first_res = {}
        for line in protein_lines:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]  # keep raw chain char (can be space)
            if chain in first_res:
                continue
            resname = line[17:20].strip()
            resnum = line[22:26]
            icode = line[26:27]
            first_res[chain] = (resname, resnum, icode)

        met_chains = {
            chain: (resnum, icode) for chain, (resname, resnum, icode) in first_res.items() if resname == "MET"
        }
        if not met_chains:
            return protein_lines

        ff_name = (self.forcefield_name or "").lower()
        should_drop = mode == "drop" or (mode == "auto" and "charmm36" in ff_name)

        if not should_drop:
            if self.verbose:
                chains = ",".join([c.strip() or "-" for c in met_chains.keys()])
                print(f"  Detected N-terminal MET on chain(s): {chains}")
                print("  Tip: set protein_preparation.nterm_met = 'drop' to remove them")
            return protein_lines

        # Drop N-terminal MET atoms
        new_lines = []
        removed_atoms = 0
        removed_chains = []
        for line in protein_lines:
            if not line.startswith("ATOM"):
                new_lines.append(line)
                continue
            chain = line[21]
            if chain not in met_chains:
                new_lines.append(line)
                continue
            resnum, icode = met_chains[chain]
            if line[22:26] == resnum and line[26:27] == icode and line[17:20].strip() == "MET":
                removed_atoms += 1
                if chain not in removed_chains:
                    removed_chains.append(chain)
                continue
            new_lines.append(line)

        if removed_atoms > 0:
            self.stats["nterm_met_removed"] = len(removed_chains)
            self.stats["nterm_met_atoms_removed"] = removed_atoms
            self.stats["nterm_met_chains"] = [c.strip() or "-" for c in removed_chains]
            self.stats["protein_atoms"] -= removed_atoms
            self.stats["total_atoms"] -= removed_atoms
            if self.verbose:
                chains = ",".join(self.stats["nterm_met_chains"])
                reason = "auto (CHARMM36)" if mode == "auto" else "user request"
                print(f"  Removed N-terminal MET on chain(s): {chains} ({reason})")

        return new_lines

    def _extract_protein_coords(self, protein_lines: List[str]) -> np.ndarray:
        """
        Extract heavy atom coordinates from protein ATOM records.

        Parameters
        ----------
        protein_lines : list of str
            ATOM records from PDB file

        Returns
        -------
        np.ndarray
            Coordinates array of shape (n_atoms, 3)
        """
        coords = []
        for line in protein_lines:
            # Skip hydrogen atoms
            atom_name = line[12:16].strip()
            if atom_name.startswith("H"):
                continue

            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue

        return np.array(coords) if coords else np.array([]).reshape(0, 3)

    def _group_by_chain(self, protein_lines: List[str], metal_lines: List[str]) -> dict:
        """
        Group protein and metal atoms by chain ID for proper PDB output.

        IMPORTANT: pdb2gmx requires each chain to have a consistent residue type.
        Since metals are type 'Ion' and proteins are type 'Protein', they cannot
        be in the same chain. We assign each metal to a NEW unique chain ID.

        Parameters
        ----------
        protein_lines : list of str
            Protein ATOM records
        metal_lines : list of str
            Metal ATOM records (converted from HETATM)

        Returns
        -------
        dict
            {chain_id: ([protein_atoms], [metal_atoms])}
        """
        import string

        chains = {}

        # Group protein atoms by chain
        for line in protein_lines:
            if len(line) > 21:
                chain_id = line[21]
                if chain_id not in chains:
                    chains[chain_id] = ([], [])
                chains[chain_id][0].append(line)

        # Assign metals to NEW unique chain IDs (not mixed with protein)
        # This is required because pdb2gmx doesn't allow mixed Protein/Ion types
        used_chains = set(chains.keys())
        available_chains = [c for c in string.ascii_uppercase + string.digits if c not in used_chains]

        metal_chain_idx = 0
        for line in metal_lines:
            if len(line) > 21:
                original_chain = line[21]

                # Assign a new unique chain ID for this metal
                if metal_chain_idx < len(available_chains):
                    new_chain_id = available_chains[metal_chain_idx]
                else:
                    # Fallback: use numbers or special characters
                    new_chain_id = str(metal_chain_idx % 10)

                metal_chain_idx += 1

                # Modify the line to use the new chain ID
                # PDB format: chain ID is at position 21
                modified_line = line[:21] + new_chain_id + line[22:]

                # Create new chain entry for this metal
                if new_chain_id not in chains:
                    chains[new_chain_id] = ([], [])
                chains[new_chain_id][1].append(modified_line)

                if self.verbose:
                    res_name = line[17:20].strip()
                    resnum = line[22:26].strip()
                    print(
                        f"  Reassigning {res_name} {resnum} from chain {original_chain} to chain {new_chain_id} (for pdb2gmx compatibility)"
                    )

        return chains
