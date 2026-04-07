#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HETATM processing utilities for protein cleaning.
"""

import numpy as np
from typing import List, Tuple


class HETATMMixin:
    """Mixin for HETATM record processing."""

    def _process_hetatm_records(self, hetatm_lines: List[str], protein_coords: np.ndarray) -> List[str]:
        """
        Process HETATM records according to ion handling mode.

        Parameters
        ----------
        hetatm_lines : list of str
            HETATM records from PDB file
        protein_coords : np.ndarray
            Protein heavy atom coordinates

        Returns
        -------
        list of str
            HETATM records to keep
        """
        kept_lines = []

        for line in hetatm_lines:
            residue_name = line[17:20].strip().upper()
            atom_name = line[12:16].strip().upper()

            # Handle water molecules
            if residue_name in self.WATER_NAMES:
                if self.keep_crystal_water:
                    kept_lines.append(line)
                    self.stats["water_kept"] += 1
                else:
                    self.stats["water_removed"] += 1
                continue

            # Handle crystallization artifacts
            if residue_name in self.CRYSTALLIZATION_ARTIFACTS:
                if self.remove_artifacts:
                    # Remove if flag is True (default)
                    self.stats["artifacts_removed"] += 1
                else:
                    # Keep if flag is False
                    kept_lines.append(line)
                continue

            # Check if this is a metal/ion
            is_metal = self._is_metal_or_ion(residue_name, atom_name)

            if not is_metal:
                # Not a metal/ion or artifact - keep it (might be a cofactor or modified residue)
                kept_lines.append(line)
                continue

            # Handle metal/ion based on mode
            should_keep = False

            if self.ion_mode == "keep_all":
                should_keep = True
            elif self.ion_mode == "smart":
                should_keep = self._should_keep_metal_smart(residue_name, atom_name)
            # elif self.ion_mode == 'remove_all': should_keep remains False

            # Distance filtering for kept metals (only if cutoff is set)
            if should_keep and self.distance_cutoff is not None and len(protein_coords) > 0:
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    metal_coord = np.array([x, y, z])

                    min_dist = self._calculate_min_distance(metal_coord, protein_coords)

                    if min_dist > self.distance_cutoff:
                        should_keep = False
                        self.stats["metals_removed_by_distance"] += 1
                        if self.verbose:
                            print(f"  Removed {residue_name} at distance {min_dist:.2f} Å from protein")
                except (ValueError, IndexError):
                    # If we can't parse coordinates, skip this atom
                    should_keep = False

            if should_keep:
                kept_lines.append(line)
                self.stats["metals_kept"] += 1
            else:
                self.stats["ions_removed"] += 1

        return kept_lines

    def _convert_metals_to_atom(self, hetatm_lines: List[str]) -> Tuple[List[str], List[str]]:
        """
        Convert metal ion HETATM records to ATOM records for pdb2gmx compatibility.

        This allows pdb2gmx to directly recognize and include metal ions in the topology,
        avoiding the need for manual addition later.

        Parameters
        ----------
        hetatm_lines : list of str
            HETATM records from cleaned PDB

        Returns
        -------
        tuple
            (metal_atom_lines, other_hetatm_lines)
            - metal_atom_lines: Metal ions converted to ATOM format
            - other_hetatm_lines: Other HETATM records (cofactors, etc.)
        """
        metal_atom_lines = []
        other_hetatm_lines = []

        for line in hetatm_lines:
            residue_name = line[17:20].strip().upper()
            atom_name = line[12:16].strip().upper()

            # Check if this is a metal ion that should be converted
            is_metal = (
                residue_name in self.structural_metals
                or atom_name in self.structural_metals
                or residue_name in self.non_structural_ions
                or atom_name in self.non_structural_ions
            )

            if is_metal:
                # Convert HETATM to ATOM
                # PDB format: columns 1-6 are record type
                atom_line = "ATOM  " + line[6:]
                metal_atom_lines.append(atom_line)

                if self.verbose:
                    chain = line[21] if len(line) > 21 else " "
                    resnum = line[22:26].strip() if len(line) > 26 else ""
                    print(f"  Converting {residue_name} (chain {chain}, res {resnum}) from HETATM to ATOM")
            else:
                # Keep as HETATM (cofactors, modified residues, etc.)
                other_hetatm_lines.append(line)

        return metal_atom_lines, other_hetatm_lines
