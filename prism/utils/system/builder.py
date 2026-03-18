#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Main SystemBuilder class that orchestrates the system building workflow.
"""

from .base import SystemBuilderBase
from .protein import ProteinProcessorMixin
from .topology import TopologyProcessorMixin
from .coordinates import CoordinateProcessorMixin
from .solvation import SolvationProcessorMixin
from .metals import MetalProcessorMixin

try:
    from ..colors import print_success, print_info
except ImportError:
    from prism.utils.colors import print_success, print_info


class SystemBuilder(
    SystemBuilderBase,
    ProteinProcessorMixin,
    TopologyProcessorMixin,
    CoordinateProcessorMixin,
    SolvationProcessorMixin,
    MetalProcessorMixin,
):
    """
    Builds a complete protein-ligand system for GROMACS simulations.

    This class encapsulates the GROMACS commands needed to process a protein,
    combine it with a ligand, create a simulation box, solvate it, and add ions.
    It is designed to be robust, especially in handling the output of GROMACS
    tools to ensure consistency between coordinate and topology files.
    """

    def build(
        self,
        cleaned_protein_path: str,
        lig_ff_dirs,
        ff_idx: int,
        water_idx: int,
        ff_info: dict = None,
        water_info: dict = None,
        nter: str = None,
        cter: str = None,
    ) -> str:
        """
        Runs the full system building workflow with support for multiple ligands.

        Args:
            cleaned_protein_path: Path to the cleaned protein PDB file.
            lig_ff_dirs: Path to directory containing ligand force field files (str),
                        or list of paths for multiple ligands, or None for protein-only.
            ff_idx: The index of the protein force field to use.
            water_idx: The index of the water model to use.
            ff_info: Dictionary containing force field info (name, dir, path)
            nter: N-terminal type for pdb2gmx (e.g., 'NH3+', 'NH2'). Default: auto-select
            cter: C-terminal type for pdb2gmx (e.g., 'COO-', 'COOH'). Default: auto-select
            water_info: Dictionary containing water model info (name, label, description)

        Returns:
            The path to the directory containing the final model files.
        """
        try:
            print_info("Building GROMACS Model")

            if (self.model_dir / "solv_ions.gro").exists() and not self.overwrite:
                print(f"Final model file solv_ions.gro already exists in {self.model_dir}, skipping model building.")
                return str(self.model_dir)

            # Store water model name for use in solvation
            self.water_model_name = water_info.get("name", "tip3p") if water_info else "tip3p"
            # Store force field info for terminal atom fixing
            self.ff_info = ff_info

            # Normalize lig_ff_dirs to always be a list (or None for protein-only)
            if lig_ff_dirs is None:
                lig_ff_list = None
            elif isinstance(lig_ff_dirs, str):
                # Legacy single ligand support
                lig_ff_list = [lig_ff_dirs]
            elif isinstance(lig_ff_dirs, list):
                lig_ff_list = lig_ff_dirs
            else:
                raise TypeError(f"lig_ff_dirs must be None, str, or list, got {type(lig_ff_dirs)}")

            fixed_pdb = self._fix_protein(cleaned_protein_path)

            # Apply PROPKA renaming AFTER pdbfixer (pdbfixer may revert HIS names)
            protonation_method = self.config.get("protonation", {}).get("method", "gromacs")
            if protonation_method == "propka":
                self._apply_propka_renaming(fixed_pdb)

            protein_gro, topol_top = self._generate_topology(
                fixed_pdb, ff_idx, water_idx, ff_info, water_info, nter, cter
            )

            # Process ligands if provided
            if lig_ff_list is not None and len(lig_ff_list) > 0:
                self._fix_topology_multi_ligands(topol_top, lig_ff_list)
                complex_gro = self._combine_protein_multi_ligands(protein_gro, lig_ff_list)
            else:
                # Protein-only system
                print("\nSkipping ligand processing (protein-only system)")
                complex_gro = protein_gro

            boxed_gro = self._create_box(complex_gro)
            solvated_gro = self._solvate(boxed_gro, topol_top)
            self._add_ions(solvated_gro, topol_top)

            print_success("System building completed")
            print(f"System files are in {self.model_dir}")
            return str(self.model_dir)

        except (RuntimeError, FileNotFoundError) as e:
            print(f"\nError during system build: {e}")
            import traceback

            traceback.print_exc()
            return ""
