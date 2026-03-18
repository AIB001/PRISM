#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Core ProteinCleaner class for PDB file cleaning.
"""

import logging
from typing import List, Optional

from .geometry import GeometryMixin
from .metal import MetalMixin
from .io import IOMixin
from .protein import ProteinMixin
from .hetatm import HETATMMixin

logger = logging.getLogger(__name__)


class ProteinCleaner(
    IOMixin,
    ProteinMixin,
    HETATMMixin,
    MetalMixin,
    GeometryMixin,
):
    """
    Clean protein PDB files by removing artifacts and processing metals.

    This class removes common crystallization artifacts (glycerol, ethylene glycol, etc.),
    processes metal ions, and handles N-terminal methionine removal.
    """

    # Class constants for metal and artifact identification
    STRUCTURAL_METALS = {
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

    # Non-structural ions that can usually be safely removed
    NON_STRUCTURAL_IONS = {
        "NA",
        "CL",
        "K",
        "BR",
        "I",
        "F",
        "NA+",
        "CL-",
        "K+",
        "BR-",
        "I-",
        "SO4",
        "PO4",
        "NO3",
        "CO3",
        "SUL",
    }

    CRYSTALLIZATION_ARTIFACTS = {
        # Glycerol and polyols
        "GOL",
        "EDO",
        "MPD",
        "MRD",
        "PGO",
        "PG4",
        "BTB",
        # PEG oligomers
        "PEG",
        "PGE",
        "1PE",
        "P6G",
        "P33",
        "PE8",
        "2PE",
        "XPE",
        "12P",
        "33O",
        "P3E",
        # Sugars (unless covalently linked to protein)
        "NAG",
        "NDG",
        "BMA",
        "MAN",
        "GAL",
        "FUC",
        "GLC",
        "SIA",
        # Detergents
        "DMS",
        "BOG",
        "LMT",
        "C8E",
        "DAO",
        "SDS",
        "LDA",
        # Other common additives
        "ACT",
        "ACE",
        "FMT",
        "TRS",
        "EPE",
        "BCT",
        "BME",
    }

    # Water residue names
    WATER_NAMES = {"HOH", "WAT", "H2O", "TIP", "SPC", "SOL"}

    def __init__(
        self,
        ion_mode: str = "smart",
        distance_cutoff: float = 5.0,
        keep_crystal_water: bool = False,
        remove_artifacts: bool = True,
        keep_custom_metals: Optional[List[str]] = None,
        remove_custom_metals: Optional[List[str]] = None,
        verbose: bool = True,
        drop_nterm_met: Optional[str] = None,
        forcefield_name: Optional[str] = None,
    ):
        """
        Initialize the ProteinCleaner.

        Parameters
        ----------
        ion_mode : str
            Ion handling mode. Options:
            - 'keep_all': Keep all metal ions and heteroatoms (except water unless keep_crystal_water=True)
            - 'smart' (default): Keep structural metals, remove non-structural ions
            - 'remove_all': Remove all metal ions and heteroatoms
        distance_cutoff : float
            Maximum distance (in Angstroms) from protein for keeping metals.
            Metals farther than this from any protein heavy atom will be removed.
            Default: 5.0 Å
        keep_crystal_water : bool
            Whether to keep crystal water molecules. Default: False
            If True, water molecules are preserved in the output
        remove_artifacts : bool
            Whether to remove crystallization artifacts (GOL, EDO, PEG, NAG, etc.). Default: True
            If False, crystallization artifacts are kept in the output
        keep_custom_metals : list of str, optional
            Additional metal/ion names to always keep (e.g., ['MO', 'W'])
        remove_custom_metals : list of str, optional
            Additional metal/ion names to always remove
        verbose : bool
            Print detailed information during cleaning
        """
        if ion_mode not in ["keep_all", "smart", "remove_all"]:
            raise ValueError(f"Invalid ion_mode: {ion_mode}. Must be 'keep_all', 'smart', or 'remove_all'")

        self.ion_mode = ion_mode
        self.distance_cutoff = distance_cutoff
        self.keep_crystal_water = keep_crystal_water
        self.remove_artifacts = remove_artifacts
        self.verbose = verbose
        self.drop_nterm_met = drop_nterm_met
        self.forcefield_name = forcefield_name

        # Build custom metal sets
        self.structural_metals = self.STRUCTURAL_METALS.copy()
        if keep_custom_metals:
            self.structural_metals.update(m.upper() for m in keep_custom_metals)

        self.non_structural_ions = self.NON_STRUCTURAL_IONS.copy()
        if remove_custom_metals:
            self.non_structural_ions.update(m.upper() for m in remove_custom_metals)

        # Statistics
        self.stats = {
            "total_atoms": 0,
            "protein_atoms": 0,
            "water_removed": 0,
            "water_kept": 0,
            "ions_removed": 0,
            "metals_kept": 0,
            "metals_removed_by_distance": 0,
            "artifacts_removed": 0,
            "nterm_met_removed": 0,
            "nterm_met_atoms_removed": 0,
            "nterm_met_chains": [],
        }

    def clean_pdb(self, input_pdb: str, output_pdb: str) -> dict:
        """
        Clean a PDB file with intelligent ion handling.

        Parameters
        ----------
        input_pdb : str
            Path to input PDB file
        output_pdb : str
            Path to output cleaned PDB file

        Returns
        -------
        dict
            Statistics about the cleaning process
        """
        if not os.path.exists(input_pdb):
            raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")

        if self.verbose:
            print(f"\n=== Cleaning Protein Structure ===")
            print(f"Input: {input_pdb}")
            print(f"Output: {output_pdb}")
            print(f"Mode: {self.ion_mode}")
            print(f"Distance cutoff: {self.distance_cutoff} Å")
            print(f"Keep crystal water: {self.keep_crystal_water}")
            print(f"Remove crystallization artifacts: {self.remove_artifacts}")

        # Read and parse PDB file
        protein_lines, hetatm_lines = self._read_pdb(input_pdb)

        # Optionally drop N-terminal MET residues
        protein_lines = self._handle_nterm_met(protein_lines)

        # Extract protein coordinates for distance calculations
        protein_coords = self._extract_protein_coords(protein_lines)

        # Process HETATM records
        kept_hetatm_lines = []
        if self.ion_mode != "remove_all":
            kept_hetatm_lines = self._process_hetatm_records(hetatm_lines, protein_coords)

        # Write cleaned PDB
        self._write_pdb(output_pdb, protein_lines, kept_hetatm_lines, input_pdb)

        # Print statistics
        if self.verbose:
            self._print_statistics()

        return self.stats.copy()

    def _print_statistics(self):
        """Print cleaning statistics."""
        print("\n=== Cleaning Statistics ===")
        print(f"Total atoms read: {self.stats['total_atoms']}")
        print(f"Protein atoms: {self.stats['protein_atoms']}")
        if self.stats.get("nterm_met_removed", 0) > 0:
            chains = ",".join(self.stats.get("nterm_met_chains", []))
            print(f"N-terminal MET removed: {self.stats['nterm_met_removed']} chain(s) ({chains})")
        print(f"Water molecules removed: {self.stats['water_removed']}")
        if self.stats["water_kept"] > 0:
            print(f"Water molecules kept: {self.stats['water_kept']}")
        if self.stats["artifacts_removed"] > 0:
            print(f"Crystallization artifacts removed: {self.stats['artifacts_removed']} atoms")
        print(f"Ions/metals removed: {self.stats['ions_removed']}")
        print(f"Metals kept: {self.stats['metals_kept']}")
        if self.stats["metals_removed_by_distance"] > 0:
            print(f"Metals removed by distance filter: {self.stats['metals_removed_by_distance']}")

        total_kept = self.stats["protein_atoms"] + self.stats["metals_kept"] + self.stats["water_kept"]
        print(f"\nTotal atoms in output: {total_kept}")
        print("=" * 30)
