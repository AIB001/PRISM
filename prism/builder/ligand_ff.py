#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder - Ligand Force Field Mixin

Handles ligand force field generation for all supported force field types
(GAFF, GAFF2, OpenFF, CGenFF, OPLS, MMFF, MATCH, Hybrid) including
RESP charge workflow (Gaussian or pre-computed).
"""

import os
import re
import shutil
import concurrent.futures

from ..utils.colors import (
    print_subheader,
    print_success,
    print_error,
    print_warning,
    print_info,
    success,
    path,
    number,
)
from ..forcefield.gaff import GAFFForceFieldGenerator
from ..forcefield.gaff2 import GAFF2ForceFieldGenerator
from ..forcefield.openff import OpenFFForceFieldGenerator
from ..forcefield.cgenff import CGenFFForceFieldGenerator
from ..forcefield.charmm_gui import CHARMMGUIForceFieldGenerator
from ..forcefield.opls_aa import OPLSAAForceFieldGenerator
from ..forcefield.swissparam import (
    BothForceFieldGenerator,
    MATCHForceFieldGenerator,
    HybridMMFFMATCHForceFieldGenerator,
)
from ..forcefield.rtf import RTFForceFieldGenerator


class LigandForceFieldMixin:
    """Mixin providing ligand force field generation capabilities."""

    def generate_ligand_forcefield(self):
        """Generate ligand force fields for all ligands using selected force field generator"""
        print_subheader(f"Generating Ligand Force Fields ({self.ligand_forcefield.upper()})")
        print(f"Processing {number(self.ligand_count)} ligand(s)...")

        # Determine force field suffix based on type
        if self.ligand_forcefield in ["gaff", "gaff2"]:
            ff_suffix = "amb2gmx"
        elif self.ligand_forcefield == "openff":
            ff_suffix = "openff2gmx"
        elif self.ligand_forcefield == "cgenff":
            ff_suffix = "cgenff2gmx"
        elif self.ligand_forcefield == "charmm-gui":
            ff_suffix = "charmm2gmx"
        elif self.ligand_forcefield == "opls":
            ff_suffix = "opls2gmx"
        elif self.ligand_forcefield in ["mmff", "match", "hybrid"]:
            ff_suffix = "sp2gmx"  # SwissParam
        elif self.ligand_forcefield == "rtf":
            ff_suffix = "rtf2gmx"  # CHARMM RTF
        else:
            ff_suffix = "lig2gmx"

        # Create base directory for ligand force fields
        # Single ligand: place directly in output_dir (backward compatibility)
        # Multi-ligand: place in Ligand_Forcefield/ subdirectory (new feature)
        if self.ligand_count == 1:
            ligand_ff_base_dir = self.output_dir
        else:
            ligand_ff_base_dir = os.path.join(self.output_dir, "Ligand_Forcefield")
            os.makedirs(ligand_ff_base_dir, exist_ok=True)

        def generate_single_ligand_ff(ligand_idx, ligand_path):
            """Generate force field for a single ligand"""
            ligand_num = ligand_idx + 1  # 1-based numbering
            print(
                f"\n{success('→')} Processing ligand {number(ligand_num)}/{number(self.ligand_count)}: {path(os.path.basename(ligand_path))}"
            )

            # Create temporary output directory for this ligand
            temp_output_dir = os.path.join(self.output_dir, f"_temp_lig_{ligand_num}")
            os.makedirs(temp_output_dir, exist_ok=True)

            # Check for existing or provided RESP file
            resp_mol2 = self._find_resp_file(ligand_path, ligand_idx)

            # Determine charge mode for GAFF/GAFF2
            # Use 'gas' (fast) if RESP charges will be applied, 'bcc' otherwise
            if resp_mol2 or self.gaussian_method:
                charge_mode = "gas"
                print(f"  Using 'gas' charge mode (RESP charges will be applied)")
            else:
                charge_mode = "bcc"

            # Generate force field based on type
            if self.ligand_forcefield == "gaff":
                generator = GAFFForceFieldGenerator(
                    ligand_path, temp_output_dir, overwrite=self.overwrite, charge_mode=charge_mode
                )
            elif self.ligand_forcefield == "gaff2":
                generator = GAFF2ForceFieldGenerator(
                    ligand_path, temp_output_dir, overwrite=self.overwrite, charge_mode=charge_mode
                )
            elif self.ligand_forcefield == "openff":
                generator = OpenFFForceFieldGenerator(
                    ligand_path,
                    temp_output_dir,
                    charge=self.config["ligand_forcefield"]["charge"],
                    overwrite=self.overwrite,
                )
            elif self.ligand_forcefield == "cgenff":
                # CGenFF website format (PDB+TOP files from cgenff.com)
                cgenff_dir = self.forcefield_paths[ligand_idx] if self.forcefield_paths else None
                generator = CGenFFForceFieldGenerator(
                    ligand_path, temp_output_dir, cgenff_dir=cgenff_dir, overwrite=self.overwrite
                )
            elif self.ligand_forcefield == "charmm-gui":
                # CHARMM-GUI format (ITP files from charmm-gui.org)
                charmm_gui_dir = self.forcefield_paths[ligand_idx] if self.forcefield_paths else None
                generator = CHARMMGUIForceFieldGenerator(
                    ligand_path, temp_output_dir, charmm_gui_dir=charmm_gui_dir, overwrite=self.overwrite
                )
            elif self.ligand_forcefield == "opls":
                generator = OPLSAAForceFieldGenerator(
                    ligand_path,
                    temp_output_dir,
                    charge=self.config["ligand_forcefield"]["charge"],
                    charge_model="cm1a",
                    align_to_input=True,
                    overwrite=self.overwrite,
                )
            elif self.ligand_forcefield == "mmff":
                generator = BothForceFieldGenerator(ligand_path, temp_output_dir, overwrite=self.overwrite)
            elif self.ligand_forcefield == "match":
                generator = MATCHForceFieldGenerator(ligand_path, temp_output_dir, overwrite=self.overwrite)
            elif self.ligand_forcefield == "hybrid":
                generator = HybridMMFFMATCHForceFieldGenerator(ligand_path, temp_output_dir, overwrite=self.overwrite)
            elif self.ligand_forcefield == "rtf":
                # RTF requires RTF+PRM+PDB files (from MATCH or CHARMM-GUI)
                # ligand_path can be:
                # 1. A directory containing .rtf, .prm, .pdb files
                # 2. A PDB file (will infer RTF directory from PDB file location)
                import glob

                # Determine RTF directory
                if os.path.isdir(ligand_path):
                    rtf_dir = ligand_path
                elif ligand_path.endswith(".pdb"):
                    # Extract directory from PDB file path
                    # Expected: /path/to/ligand.pdb
                    # Look for /path/to/ligand.rtf and /path/to/ligand.prm
                    pdb_base = ligand_path.replace(".pdb", "")
                    rtf_dir = os.path.dirname(ligand_path)
                elif self.forcefield_paths and ligand_idx < len(self.forcefield_paths):
                    rtf_dir = self.forcefield_paths[ligand_idx]
                else:
                    raise ValueError(
                        f"RTF force field requires a directory containing .rtf, .prm, .pdb files "
                        f"or a PDB file with corresponding RTF/PRM files in the same directory. "
                        f"Provided ligand_path: {ligand_path}"
                    )

                # Find RTF, PRM, PDB files
                rtf_files = glob.glob(os.path.join(rtf_dir, "*.rtf"))
                prm_files = glob.glob(os.path.join(rtf_dir, "*.prm"))
                pdb_files = glob.glob(os.path.join(rtf_dir, "*.pdb"))

                if not rtf_files or not prm_files or not pdb_files:
                    raise FileNotFoundError(
                        f"RTF force field requires .rtf, .prm, and .pdb files in {rtf_dir}. "
                        f"Found: {len(rtf_files)} RTF, {len(prm_files)} PRM, {len(pdb_files)} PDB files"
                    )

                # Use first matching files (usually one set per ligand)
                rtf_file = rtf_files[0]
                prm_file = prm_files[0]
                pdb_file = pdb_files[0]

                print_info(
                    f"  Using RTF files: {os.path.basename(rtf_file)}, {os.path.basename(prm_file)}, {os.path.basename(pdb_file)}"
                )

                generator = RTFForceFieldGenerator(
                    rtf_file=rtf_file,
                    prm_file=prm_file,
                    pdb_file=pdb_file,
                    output_dir=temp_output_dir,
                    overwrite=self.overwrite,
                )

            temp_ff_dir = generator.run()

            # Create final directory with proper naming
            # Single ligand: LIG.xxx2gmx (no suffix, backward compatibility)
            # Multi-ligand: LIG.xxx2gmx_1, LIG.xxx2gmx_2, etc.
            if self.ligand_count == 1:
                final_ff_dir = os.path.join(ligand_ff_base_dir, f"LIG.{ff_suffix}")
            else:
                final_ff_dir = os.path.join(ligand_ff_base_dir, f"LIG.{ff_suffix}_{ligand_num}")

            # Move files from temp directory to final location
            if os.path.exists(final_ff_dir) and self.overwrite:
                shutil.rmtree(final_ff_dir)
            shutil.move(temp_ff_dir, final_ff_dir)

            # Clean up temp directory
            if os.path.exists(temp_output_dir):
                shutil.rmtree(temp_output_dir)

            # Only rename if we have multiple ligands (ligand_count > 1)
            # Single ligand keeps the name "LIG" for backward compatibility
            if self.ligand_count > 1:
                self._rename_multi_ligand_files(final_ff_dir, ligand_num)

            # Verify required files exist
            required_files = ["LIG.itp", "LIG.gro"]
            for filename in required_files:
                filepath = os.path.join(final_ff_dir, filename)
                if not os.path.exists(filepath):
                    raise FileNotFoundError(f"Required file not found: {filepath}")

            # Handle RESP charges
            if resp_mol2:
                # Apply existing RESP charges
                self._apply_resp_charges(final_ff_dir, resp_mol2, ligand_num)
            elif self.gaussian_method and self.ligand_forcefield in ["gaff", "gaff2"]:
                # Generate Gaussian RESP workflow files
                resp_mol2 = self._handle_gaussian_workflow(ligand_path, final_ff_dir, ligand_idx)
                if resp_mol2:
                    self._apply_resp_charges(final_ff_dir, resp_mol2, ligand_num)

            print_success(f"  Ligand {number(ligand_num)} force field generated: {path(final_ff_dir)}")
            return final_ff_dir

        # Generate force fields in parallel using ThreadPoolExecutor
        print_info("Generating force fields in parallel...")
        self.lig_ff_dirs = []

        with concurrent.futures.ThreadPoolExecutor(max_workers=min(self.ligand_count, 4)) as executor:
            futures = {executor.submit(generate_single_ligand_ff, i, lp): i for i, lp in enumerate(self.ligand_paths)}

            for future in concurrent.futures.as_completed(futures):
                ligand_idx = futures[future]
                try:
                    ff_dir = future.result()
                    self.lig_ff_dirs.append((ligand_idx, ff_dir))
                except Exception as e:
                    print_error(f"Failed to generate force field for ligand {ligand_idx + 1}: {e}")
                    raise

        # Sort by ligand index to maintain order
        self.lig_ff_dirs.sort(key=lambda x: x[0])
        self.lig_ff_dirs = [ff_dir for _, ff_dir in self.lig_ff_dirs]

        print_success(f"All {number(self.ligand_count)} ligand force fields generated successfully")
        print(f"Force field directory: {path(ligand_ff_base_dir)}")

        return ligand_ff_base_dir

    def _rename_multi_ligand_files(self, final_ff_dir: str, ligand_num: int) -> None:
        """Rename LIG residue to LIG_N in ITP and GRO files for multi-ligand systems."""
        itp_file = os.path.join(final_ff_dir, "LIG.itp")
        gro_file = os.path.join(final_ff_dir, "LIG.gro")

        # Rename moleculetype in ITP file
        with open(itp_file, "r") as f:
            itp_lines = f.readlines()

        with open(itp_file, "w") as f:
            in_moleculetype = False
            moleculetype_replaced = False
            for line in itp_lines:
                # Detect [ moleculetype ] section
                if "[ moleculetype ]" in line:
                    in_moleculetype = True
                    f.write(line)
                    continue

                # Replace the moleculetype name in the next non-comment, non-empty line
                if in_moleculetype and not moleculetype_replaced:
                    if line.strip() and not line.strip().startswith(";"):
                        line = re.sub(r"\bLIG\b", f"LIG_{ligand_num}", line, count=1)
                        moleculetype_replaced = True
                        in_moleculetype = False

                # Replace residue names in [ atoms ] section
                if not line.strip().startswith(";") and " LIG " in line:
                    line = re.sub(r"(\s+)LIG(\s+)", rf"\1LIG_{ligand_num}\2", line)

                f.write(line)

        # Rename residue in GRO file - more robust replacement
        with open(gro_file, "r") as f:
            gro_lines = f.readlines()
        with open(gro_file, "w") as f:
            for i, line in enumerate(gro_lines):
                # Skip header (line 0) and atom count (line 1) and box vectors (last line)
                if i > 1 and i < len(gro_lines) - 1:
                    # GRO format: columns 0-4 (residue number), 5-9 (residue name)
                    if len(line) >= 20:  # Valid atom line
                        res_name = line[5:10].strip()
                        if res_name.startswith("LIG"):
                            new_res_name = f"LIG_{ligand_num}"
                            line = line[:5] + new_res_name.ljust(5)[:5] + line[10:]
                f.write(line)

    def _find_resp_file(self, ligand_path: str, ligand_idx: int):
        """
        Find RESP mol2 file for a ligand.

        Checks in order:
        1. Explicitly provided via --respfile
        2. Auto-detect ligand_resp.mol2 in same directory as ligand
        """
        # Check --respfile parameter
        if self.resp_files and ligand_idx < len(self.resp_files):
            resp_file = self.resp_files[ligand_idx]
            if os.path.exists(resp_file):
                print_info(f"  Using provided RESP file: {path(resp_file)}")
                return resp_file
            else:
                print_warning(f"  RESP file not found: {resp_file}")

        # Auto-detect ligand_resp.mol2 in same directory
        ligand_dir = os.path.dirname(ligand_path)
        auto_resp = os.path.join(ligand_dir, "ligand_resp.mol2")
        if os.path.exists(auto_resp):
            print_info(f"  Auto-detected RESP file: {path(auto_resp)}")
            return auto_resp

        return None

    def _handle_gaussian_workflow(self, ligand_path: str, ff_dir: str, ligand_idx: int):
        """
        Handle Gaussian RESP workflow for a ligand.

        If Gaussian (g16) is available, runs the calculation automatically.
        Otherwise, generates input files and a run script for manual execution.
        """
        try:
            from ..gaussian import GaussianRESPWorkflow, check_gaussian_available
        except ImportError:
            print_error("Gaussian module not available")
            return None

        print_info(f"  Setting up Gaussian RESP workflow...")
        if ligand_idx is not None:
            print_info(f"  Ligand index: {ligand_idx}")

        # Create Gaussian working directory
        gaussian_dir = os.path.join(ff_dir, "gaussian_resp")
        os.makedirs(gaussian_dir, exist_ok=True)

        # Initialize workflow
        workflow = GaussianRESPWorkflow(
            ligand_path=ligand_path,
            output_dir=gaussian_dir,
            method=self.gaussian_method,
            do_optimization=self.do_optimization,
            nproc=self.gaussian_nproc,
            mem=self.gaussian_mem,
            verbose=True,
        )

        # Check if Gaussian is available
        if check_gaussian_available():
            print_info(f"  Gaussian (g16) detected, running RESP calculation...")
            resp_mol2 = workflow.run()
            if resp_mol2:
                print_success(f"  RESP charges calculated: {path(resp_mol2)}")
                return resp_mol2
            else:
                print_warning(f"  Gaussian calculation failed")
                return None
        else:
            # Generate files for manual execution
            files = workflow.prepare_files()
            print_warning(f"  Gaussian (g16) not found in PATH")
            print_warning(f"  Generated Gaussian input files in: {path(gaussian_dir)}")
            print_warning(f"  To calculate RESP charges:")
            print_warning(f"    1. Run: bash {files['script']}")
            print_warning(f"    2. Re-run PRISM with: --respfile {os.path.join(gaussian_dir, 'ligand_resp.mol2')}")
            return None

    def _apply_resp_charges(self, ff_dir: str, resp_mol2: str, ligand_num: int) -> None:
        """Apply RESP charges to ITP file."""
        try:
            from ..gaussian import RESPChargeReplacer
        except ImportError:
            print_error("Gaussian module not available for charge replacement")
            return

        itp_path = os.path.join(ff_dir, "LIG.itp")
        if not os.path.exists(itp_path):
            print_error(f"  ITP file not found: {itp_path}")
            return

        print_info(f"  Applying RESP charges to ligand {number(ligand_num)}...")

        try:
            replacer = RESPChargeReplacer(resp_mol2, verbose=False)
            replacer.replace_itp_charges(itp_path, backup=True)
            print_success(f"  RESP charges applied from: {path(resp_mol2)}")
        except Exception as e:
            print_error(f"  Failed to apply RESP charges: {e}")
