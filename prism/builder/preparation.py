#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder - System Preparation Mixin

Handles protein cleaning/preparation, MDP file generation, model building,
local run script generation, and cleanup operations.
"""

import os
import glob
import shutil
from pathlib import Path

from ..utils.colors import (
    print_subheader,
)


class SystemPreparationMixin:
    """Mixin providing protein preparation, MDP generation, and script utilities."""

    def clean_protein(self, ion_mode=None, distance_cutoff=None, keep_crystal_water=None, remove_artifacts=None):
        """
        Clean the protein PDB file with intelligent metal ion handling.

        Parameters
        ----------
        ion_mode : str, optional
            Ion handling mode. Options:
            - 'keep_all': Keep all metal ions (except water unless keep_crystal_water=True)
            - 'smart' (default): Keep structural metals (Zn, Mg, Ca, Fe, etc.), remove non-structural ions (Na, Cl, etc.)
            - 'remove_all': Remove all metal ions
            If not specified, reads from config or defaults to 'smart'
        distance_cutoff : float, optional
            Maximum distance (Angstroms) from protein for keeping metals.
            Metals farther than this will be removed even if they are structural.
            If not specified, reads from config or defaults to 5.0 A
        keep_crystal_water : bool, optional
            Whether to keep crystal water molecules. Default: False
            If True, water molecules from the crystal structure are preserved.
            If not specified, reads from config or defaults to False
        remove_artifacts : bool, optional
            Whether to remove crystallization artifacts (GOL, EDO, PEG, NAG, etc.). Default: True
            If False, crystallization artifacts are kept in the output.
            If not specified, reads from config or defaults to True

        Returns
        -------
        str
            Path to cleaned (and optionally protonated) PDB file
        """
        from ..utils.cleaner import ProteinCleaner

        print_subheader("Cleaning Protein")

        cleaned_pdb = os.path.join(self.output_dir, f"{self.protein_name}_clean.pdb")

        # Check for existing cleaned file
        if os.path.exists(cleaned_pdb) and not self.overwrite:
            print(f"Using existing cleaned protein: {cleaned_pdb}")
            return cleaned_pdb

        # Get parameters from config if not explicitly provided
        if ion_mode is None:
            ion_mode = self.config.get("protein_preparation", {}).get("ion_handling_mode", "smart")
        if distance_cutoff is None:
            distance_cutoff = self.config.get("protein_preparation", {}).get("metal_distance_cutoff", 5.0)
        if keep_crystal_water is None:
            keep_crystal_water = self.config.get("protein_preparation", {}).get("keep_crystal_water", False)
        if remove_artifacts is None:
            remove_artifacts = self.config.get("protein_preparation", {}).get("remove_crystallization_artifacts", True)
        nterm_met = self.config.get("protein_preparation", {}).get("nterm_met", "keep")

        # Get custom metals from config
        keep_custom = self.config.get("protein_preparation", {}).get("keep_custom_metals", [])
        remove_custom = self.config.get("protein_preparation", {}).get("remove_custom_metals", [])

        # Initialize cleaner
        cleaner = ProteinCleaner(
            ion_mode=ion_mode,
            distance_cutoff=distance_cutoff,
            keep_crystal_water=keep_crystal_water,
            remove_artifacts=remove_artifacts,
            keep_custom_metals=keep_custom if keep_custom else None,
            remove_custom_metals=remove_custom if remove_custom else None,
            verbose=True,
            drop_nterm_met=nterm_met,
            forcefield_name=self.forcefield["name"] if self.forcefield else None,
        )

        # Clean the protein
        cleaner.clean_pdb(self.protein_path, cleaned_pdb)

        # Some prepared receptors keep non-canonical terminal cap atoms
        # (e.g. CAY/CY/OY or CAT/NT) on standard amino-acid residue names.
        # Strip them immediately so downstream pdbfixer/pdb2gmx sees a
        # canonical protein backbone.
        from prism.utils.cleaner import fix_terminal_atoms

        fix_terminal_atoms(
            cleaned_pdb,
            cleaned_pdb,
            force_field=self.forcefield["name"] if self.forcefield else None,
            verbose=True,
        )

        # Post-process: Fix terminal hydrogen names for AMBER compatibility
        self._fix_hydrogen_names(cleaned_pdb)

        print(f"Protein cleaned and saved to: {cleaned_pdb}")

        # Note: PROPKA HIS renaming (if --protonation propka) is applied later
        # in system.py's build() AFTER pdbfixer, to avoid pdbfixer reverting names.
        return cleaned_pdb

    def _fix_hydrogen_names(self, pdb_file: str) -> None:
        """
        Fix terminal hydrogen names for AMBER compatibility.

        AMBER force fields expect specific naming for N-terminal hydrogens:
        H1, H2, H3 instead of HN1, HN2, HN3
        """
        with open(pdb_file, "r") as f:
            lines = f.readlines()

        fixed_lines = []
        for line in lines:
            if line.startswith("ATOM"):
                # Fix terminal hydrogen names
                if "HN1" in line:
                    line = line.replace("HN1", "H1 ")
                elif "HN2" in line:
                    line = line.replace("HN2", "H2 ")
                elif "HN3" in line:
                    line = line.replace("HN3", "H3 ")
            fixed_lines.append(line)

        with open(pdb_file, "w") as f:
            f.writelines(fixed_lines)

    def build_model(self, cleaned_protein: str):
        """Build the GROMACS model with support for multiple ligands"""
        return self.system_builder.build(
            cleaned_protein,
            self.lig_ff_dirs,  # Now a list of force field directories
            self.forcefield_idx,
            self.water_model_idx,
            self.forcefield,  # Pass full force field info (name, dir, path)
            self.water_model,  # Pass full water model info
            nter=self.nter,  # N-terminal type
            cter=self.cter,  # C-terminal type
        )

    def generate_mdp_files(self):
        """Generate MDP files for MD simulations"""
        self.mdp_generator.generate_all()

    def save_config(self):
        """Save the current configuration to a file"""
        config_file = os.path.join(self.output_dir, "prism_config.yaml")
        self.config_manager.save_config(config_file)
        print(f"Configuration saved to: {config_file}")

    def cleanup(self):
        """Clean up temporary files"""
        print_subheader("Cleaning up temporary files")

        # Cleanup directories for GAFF, GAFF2, and OpenFF
        cleanup_dirs = ["forcefield", "temp_openff"]

        for dir_name in cleanup_dirs:
            cleanup_dir = os.path.join(self.output_dir, dir_name)
            if os.path.exists(cleanup_dir):
                if self.ligand_forcefield in ["gaff", "gaff2"]:
                    temp_patterns = [
                        "*.frcmod",
                        "*.prep",
                        "*.prmtop",
                        "*.rst7",
                        "*.log",
                        "*.in",
                        "ANTECHAMBER*",
                        "ATOMTYPE*",
                        "PREP*",
                        "NEWPDB*",
                        "sqm*",
                        "leap*",
                    ]

                    for pattern in temp_patterns:
                        for file_path in Path(cleanup_dir).glob(pattern):
                            try:
                                os.remove(file_path)
                            except Exception:
                                pass
                else:
                    # For OpenFF, remove the temporary directory entirely
                    try:
                        shutil.rmtree(cleanup_dir)
                    except Exception:
                        pass

        print("Cleanup completed")

    def generate_localrun_script(self):
        """Generate localrun.sh script for easy MD execution"""
        gmx_md_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
        if not os.path.exists(gmx_md_dir):
            print("Warning: GMX_PROLIG_MD directory not found, skipping localrun.sh generation")
            return None

        # Clean up Emacs backup files (#topol.top.1#, #topol.top.2#, etc.)
        backup_pattern = os.path.join(gmx_md_dir, "#*#")
        backup_files = glob.glob(backup_pattern)
        if backup_files:
            print(f"Cleaning up {len(backup_files)} Emacs backup file(s)...")
            for backup_file in backup_files:
                try:
                    os.remove(backup_file)
                    print(f"  Removed: {os.path.basename(backup_file)}")
                except Exception as e:
                    print(f"  Warning: Could not remove {backup_file}: {e}")

        script_path = os.path.join(gmx_md_dir, "localrun.sh")

        script_content = """#!/bin/bash

######################################################
# SIMULATION PART
######################################################

# Energy Minimization (EM)
mkdir -p em
if [ -f ./em/em.gro ]; then
    echo "EM already completed, skipping..."
elif [ -f ./em/em.tpr ]; then
    echo "EM tpr file found, continuing from checkpoint..."
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v -cpi ./em/em.cpt
else
    echo "Starting EM from scratch..."
    gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o ./em/em.tpr -maxwarn 999
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v
fi

# NVT Equilibration
mkdir -p nvt
if [ -f ./nvt/nvt.gro ]; then
    echo "NVT already completed, skipping..."
elif [ -f ./nvt/nvt.tpr ]; then
    echo "NVT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v -cpi ./nvt/nvt.cpt
else
    echo "Starting NVT from scratch..."
    gmx grompp -f ../mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro -p topol.top -o ./nvt/nvt.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v
fi

# NPT Equilibration
mkdir -p npt
if [ -f ./npt/npt.gro ]; then
    echo "NPT already completed, skipping..."
elif [ -f ./npt/npt.tpr ]; then
    echo "NPT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v -cpi ./npt/npt.cpt
else
    echo "Starting NPT from scratch..."
    gmx grompp -f ../mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v
fi

# Production MD
mkdir -p prod
if [ -f ./prod/md.gro ]; then
    echo "Production MD already completed, skipping..."
elif [ -f ./prod/md.tpr ]; then
    echo "Production MD tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v -cpi ./prod/md.cpt
else
    echo "Starting Production MD from scratch..."
    gmx grompp -f ../mdps/md.mdp -c ./npt/npt.gro -r ./npt/npt.gro -p topol.top -o ./prod/md.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v
fi
"""

        with open(script_path, "w") as f:
            f.write(script_content)

        # Make script executable
        os.chmod(script_path, 0o755)

        print(f"Local run script generated: {script_path}")
        return script_path
