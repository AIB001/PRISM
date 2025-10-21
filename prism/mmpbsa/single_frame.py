#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single-Frame MM/PBSA Calculator

This module performs MM/PBSA calculations on a single structure (e.g., from docking).
It converts a GRO file to XTC format and runs gmx_MMPBSA for binding free energy estimation.
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path


class SingleFrameMMPBSA:
    """
    Single-frame MM/PBSA calculator for protein-ligand complexes

    This class handles MM/PBSA calculations for a single structure, typically
    from molecular docking. It requires:
    1. A GRO file with the complex structure
    2. A PRISM-generated system directory with topology and force field parameters
    """

    def __init__(self, structure_file, system_dir, output_dir=None, overwrite=False):
        """
        Initialize Single-Frame MM/PBSA calculator

        Parameters
        ----------
        structure_file : str
            Path to GRO file containing the complex structure
        system_dir : str
            Path to PRISM-generated GMX_PROLIG_MD directory
        output_dir : str, optional
            Output directory for MM/PBSA results (default: same as structure_file directory)
        overwrite : bool
            Whether to overwrite existing output directory
        """
        self.structure_file = Path(structure_file).resolve()
        self.system_dir = Path(system_dir).resolve()

        # Default output directory to structure file's directory
        if output_dir is None:
            self.output_dir = self.structure_file.parent
        else:
            self.output_dir = Path(output_dir).resolve()

        self.overwrite = overwrite

        # Validate inputs
        self._validate_inputs()

        # Key files
        self.topol_file = self.system_dir / "topol.top"

        # Try to find ligand MOL2 file (required for GAFF parameters)
        ligand_mol2_locations = [
            self.system_dir / "forcefield" / "ligand.mol2",
            self.system_dir / "ligand.mol2",
            self.system_dir / ".." / "forcefield" / "ligand.mol2"
        ]
        self.ligand_mol2 = None
        for mol2_path in ligand_mol2_locations:
            if mol2_path.exists():
                self.ligand_mol2 = mol2_path.resolve()
                break

        # Try to find index file in multiple locations
        index_locations = [
            self.system_dir / "index.ndx",
            self.system_dir / "em" / "index.ndx",
            self.structure_file.parent / "index.ndx"
        ]
        self.index_file = None
        for index_path in index_locations:
            if index_path.exists():
                self.index_file = index_path
                break

        # Try to find TPR file in system directory (em/em.tpr)
        self.tpr_file = self.system_dir / "em" / "em.tpr"
        if not self.tpr_file.exists():
            # Try alternative location
            self.tpr_file = self.structure_file.parent / "em.tpr"
            if not self.tpr_file.exists():
                self.tpr_file = None

        # Find group numbers from index file
        self.protein_group = None
        self.ligand_group = None
        self.system_group = None
        self._parse_index_groups()

        print(f"\nInitialized Single-Frame MM/PBSA Calculator:")
        print(f"  Structure: {self.structure_file}")
        print(f"  System directory: {self.system_dir}")
        print(f"  Output directory: {self.output_dir}")
        print(f"  TPR file: {self.tpr_file if self.tpr_file else 'Not found, will use PDB'}")
        print(f"  Ligand MOL2: {self.ligand_mol2 if self.ligand_mol2 else 'Not found'}")
        print(f"  Protein group: {self.protein_group}")
        print(f"  Ligand group: {self.ligand_group}")
        print(f"  System group: {self.system_group}")

    def _validate_inputs(self):
        """Validate input files and directories"""
        if not self.structure_file.exists():
            raise FileNotFoundError(f"Structure file not found: {self.structure_file}")

        if not self.system_dir.exists():
            raise FileNotFoundError(f"System directory not found: {self.system_dir}")

        if not self.structure_file.suffix == '.gro':
            raise ValueError(f"Structure file must be a GRO file, got: {self.structure_file.suffix}")

        # Check for required files in system directory
        required_files = ['topol.top']
        for req_file in required_files:
            file_path = self.system_dir / req_file
            if not file_path.exists():
                raise FileNotFoundError(f"Required file not found in system directory: {req_file}")

    def _parse_index_groups(self):
        """Parse index file to find Protein, LIG, and System group numbers"""
        if self.index_file is None or not self.index_file.exists():
            if self.index_file is None:
                print(f"  Warning: Index file not found in any of the expected locations")
            else:
                print(f"  Warning: Index file not found at {self.index_file}")
            print(f"  Will generate index file during execution")
            return

        try:
            with open(self.index_file, 'r') as f:
                current_group_num = -1
                for line in f:
                    line = line.strip()

                    # Check for group header: [ group_name ]
                    if line.startswith('[') and line.endswith(']'):
                        current_group_num += 1
                        group_name = line[1:-1].strip()

                        # Match group names (case-insensitive)
                        group_name_lower = group_name.lower()
                        if group_name_lower == 'protein':
                            self.protein_group = current_group_num
                            print(f"  Found Protein group: {current_group_num}")
                        elif group_name_lower == 'lig' or group_name == 'LIG':
                            self.ligand_group = current_group_num
                            print(f"  Found LIG group: {current_group_num}")
                        elif group_name_lower == 'system':
                            self.system_group = current_group_num
                            print(f"  Found System group: {current_group_num}")

            # Validate that we found required groups
            if self.protein_group is None:
                raise ValueError("Could not find 'Protein' group in index file")
            if self.ligand_group is None:
                raise ValueError("Could not find 'LIG' group in index file")
            if self.system_group is None:
                print("  Warning: Could not find 'System' group, will use default (group 0)")
                self.system_group = 0

        except Exception as e:
            raise RuntimeError(f"Error parsing index file: {e}")

    def _run_command(self, cmd, cwd=None, input_text=None):
        """Run shell command and handle errors"""
        if isinstance(cmd, list):
            cmd_str = ' '.join(cmd)
        else:
            cmd_str = cmd

        print(f"\nExecuting: {cmd_str}")
        if cwd:
            print(f"Working directory: {cwd}")

        try:
            result = subprocess.run(
                cmd_str if isinstance(cmd, str) else cmd,
                shell=isinstance(cmd, str),
                cwd=cwd,
                input=input_text,
                capture_output=True,
                text=True,
                check=True
            )
            if result.stdout:
                print(result.stdout)
            return result
        except subprocess.CalledProcessError as e:
            print(f"\nError running command: {cmd_str}")
            print(f"Exit code: {e.returncode}")
            if e.stdout:
                print(f"STDOUT:\n{e.stdout}")
            if e.stderr:
                print(f"STDERR:\n{e.stderr}")
            raise RuntimeError(f"Command failed: {cmd_str}")

    def run(self):
        """Run the complete single-frame MM/PBSA workflow"""
        print(f"\n{'='*60}")
        print("Starting Single-Frame MM/PBSA Calculation")
        print(f"{'='*60}")

        try:
            # Step 1: Convert GRO to XTC
            print("\n[Step 1] Converting GRO to XTC format...")
            xtc_file = self._convert_gro_to_xtc()
            print(f"  ✓ Created: {xtc_file}")

            # Step 2: Prepare structure file (TPR or PDB)
            print("\n[Step 2] Preparing structure file...")
            if self.tpr_file and self.tpr_file.exists():
                structure_file = self.tpr_file
                print(f"  ✓ Using existing TPR: {structure_file}")
            else:
                structure_file = self._convert_gro_to_pdb()
                print(f"  ✓ Created PDB: {structure_file}")

            # Step 3: Generate/copy index file
            print("\n[Step 3] Preparing index file...")
            index_file = self._prepare_index_file()
            print(f"  ✓ Index file ready: {index_file}")

            # Step 4: Generate MM/PBSA input file
            print("\n[Step 4] Generating MM/PBSA input file...")
            mmpbsa_input = self._generate_mmpbsa_input()
            print(f"  ✓ Created: {mmpbsa_input}")

            # Step 5: Run gmx_MMPBSA
            print("\n[Step 5] Running gmx_MMPBSA...")
            self._run_gmx_mmpbsa(structure_file, xtc_file, index_file, mmpbsa_input)

            # Step 6: Parse results
            print("\n[Step 6] Parsing results...")
            results = self._parse_results()

            # Step 7: Clean up and organize results
            mmpbsa_dir = self._cleanup_and_organize()

            print(f"\n{'='*60}")
            print("Single-Frame MM/PBSA Calculation Completed Successfully!")
            print(f"{'='*60}")
            print(f"\nResults directory: {mmpbsa_dir}")
            self._print_results_summary(results, mmpbsa_dir)

            return results

        except Exception as e:
            print(f"\nError during MM/PBSA calculation: {e}")
            raise

    def _convert_gro_to_pdb(self):
        """Convert GRO file to PDB format (required by gmx_MMPBSA)"""
        pdb_file = self.output_dir / "complex.pdb"

        cmd = [
            "gmx", "editconf",
            "-f", str(self.structure_file),
            "-o", str(pdb_file)
        ]

        self._run_command(cmd)

        if not pdb_file.exists():
            raise FileNotFoundError(f"Failed to create PDB file: {pdb_file}")

        return pdb_file

    def _convert_gro_to_xtc(self):
        """Convert GRO file to single-frame XTC trajectory"""
        xtc_file = self.output_dir / "single_frame.xtc"

        cmd = [
            "gmx", "trjconv",
            "-f", str(self.structure_file),
            "-s", str(self.structure_file),
            "-o", str(xtc_file)
        ]

        # Select "System" (group 0)
        self._run_command(cmd, input_text="0\n")

        if not xtc_file.exists():
            raise FileNotFoundError(f"Failed to create XTC file: {xtc_file}")

        return xtc_file

    def _prepare_index_file(self):
        """Prepare index file for gmx_MMPBSA"""
        output_index = self.output_dir / "index.ndx"

        # Check if index file exists
        if self.index_file and self.index_file.exists():
            # Check if source and destination are the same
            if self.index_file.resolve() == output_index.resolve():
                print(f"  Using existing index file: {self.index_file}")
                return output_index
            else:
                print(f"  Copying existing index file from {self.index_file}")
                shutil.copy(self.index_file, output_index)
        else:
            # Generate index file from TPR (preferred) or structure
            print(f"  Generating index file from structure...")
            if self.tpr_file and self.tpr_file.exists():
                cmd = [
                    "gmx", "make_ndx",
                    "-f", str(self.tpr_file),
                    "-o", str(output_index)
                ]
            else:
                cmd = [
                    "gmx", "make_ndx",
                    "-f", str(self.structure_file),
                    "-o", str(output_index)
                ]
            # Quit immediately (just generate default groups)
            self._run_command(cmd, input_text="q\n")

            # Now parse the newly created index file to get group numbers
            self.index_file = output_index
            self._parse_index_groups()

        return output_index

    def _generate_mmpbsa_input(self):
        """Generate gmx_MMPBSA input file"""
        mmpbsa_input = self.output_dir / "mmpbsa.in"

        # Standard MM/PBSA parameters in Fortran namelist format
        # Based on gmx_MMPBSA best practices and official documentation
        content = """# gmx_MMPBSA input file for single-frame calculation
# Generated by PRISM
# Reference: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/

# General namelist variables
&general
  sys_name             = "Single-Frame-MMPBSA"
  startframe           = 1
  endframe             = 1
  interval             = 1
  forcefields          = "oldff/leaprc.ff99SB,leaprc.gaff"
  ions_parameters      = 1
  PBRadii              = 3          # mbondi3 radii (recommended for ff99SB)
  temperature          = 298.15
  qh_entropy           = 0
  interaction_entropy  = 0
  c2_entropy           = 0
  assign_chainID       = 0
  keep_files           = 2
  solvated_trajectory  = 1
  verbose              = 2
/

# Poisson-Boltzmann namelist variables
&pb
  ipb                  = 2          # Linear PB equation
  inp                  = 2          # Full nonpolar (SASA + SAV + WCA)
  indi                 = 1.0        # Solute dielectric constant
  exdi                 = 80.0       # Solvent dielectric constant (water)
  istrng               = 0.15       # Ionic strength 0.15 M (physiological)
  radiopt              = 1          # Use optimized radii
  prbrad               = 1.4        # Probe radius for SASA (Å)
  fillratio            = 4.0        # Grid fill ratio
  scale                = 2.0        # Grid spacing scale
  linit                = 1000       # Linear solver iterations
  eneopt               = 2          # Energy output option
/
"""

        with open(mmpbsa_input, 'w') as f:
            f.write(content)

        return mmpbsa_input

    def _run_gmx_mmpbsa(self, structure_file, xtc_file, index_file, mmpbsa_input):
        """Run gmx_MMPBSA calculation"""
        # gmx_MMPBSA will automatically generate AMBER topology from TPR
        # DO NOT provide -cp (topology) parameter - it causes atom count mismatch
        cmd = [
            "gmx_MMPBSA",
            "-O",  # Overwrite
            "-i", str(mmpbsa_input),
            "-cs", str(structure_file),  # TPR or PDB file for structure
            "-ci", str(index_file),
            "-cg", str(self.protein_group), str(self.ligand_group),  # Group numbers from index
            "-ct", str(xtc_file),
            "-o", str(self.output_dir / "FINAL_RESULTS_MMPBSA.dat"),
            "-eo", str(self.output_dir / "FINAL_RESULTS_MMPBSA.csv")
        ]

        # Add ligand MOL2 file if available (required for GAFF parameters)
        if self.ligand_mol2 and self.ligand_mol2.exists():
            cmd.extend(["-lm", str(self.ligand_mol2)])

        # Run gmx_MMPBSA - may exit with code 1 if gmx_MMPBSA_ana GUI fails
        # but the main calculation can still succeed, so check for output files
        try:
            self._run_command(cmd, cwd=str(self.output_dir))
        except RuntimeError as e:
            # Check if results files exist despite error
            results_dat = self.output_dir / "FINAL_RESULTS_MMPBSA.dat"
            results_csv = self.output_dir / "FINAL_RESULTS_MMPBSA.csv"

            if results_dat.exists() and results_csv.exists():
                print("\n  Warning: gmx_MMPBSA exited with error (likely GUI tool failure)")
                print("  However, results files were generated successfully - continuing...")
            else:
                # Real failure - re-raise
                raise

    def _parse_results(self):
        """Parse MM/PBSA results"""
        results_file = self.output_dir / "FINAL_RESULTS_MMPBSA.dat"

        if not results_file.exists():
            print(f"  Warning: Results file not found: {results_file}")
            return None

        results = {}

        try:
            with open(results_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    # Parse key-value pairs
                    if 'DELTA TOTAL' in line or 'ΔG' in line:
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                value = float(parts[-2])  # Value
                                error = float(parts[-1])  # Standard error
                                results['binding_energy'] = value
                                results['binding_energy_error'] = error
                            except (ValueError, IndexError):
                                pass

            return results

        except Exception as e:
            print(f"  Warning: Error parsing results: {e}")
            return None

    def _print_results_summary(self, results, mmpbsa_dir=None):
        """Print summary of MM/PBSA results"""
        if not results:
            print("\n  Warning: No results to display")
            return

        # Use mmpbsa_dir if provided, otherwise use output_dir
        result_dir = mmpbsa_dir if mmpbsa_dir else self.output_dir

        print("\n" + "="*60)
        print("MM/PBSA Results Summary")
        print("="*60)

        if 'binding_energy' in results:
            energy = results['binding_energy']
            error = results.get('binding_energy_error', 0.0)
            print(f"\nBinding Free Energy (ΔG):")
            print(f"  {energy:>10.2f} ± {error:.2f} kcal/mol")

        print("\nDetailed results available in:")
        print(f"  {result_dir / 'FINAL_RESULTS_MMPBSA.dat'}")
        print(f"  {result_dir / 'FINAL_RESULTS_MMPBSA.csv'}")
        print("="*60)

    def _cleanup_and_organize(self):
        """Clean up intermediate files and organize results into GMX_PROLIG_MMPBSA folder"""
        print("\n[Step 7] Organizing results and cleaning up intermediate files...")

        # Create GMX_PROLIG_MMPBSA directory
        mmpbsa_dir = self.output_dir / "GMX_PROLIG_MMPBSA"
        mmpbsa_dir.mkdir(exist_ok=True)

        # Files to keep (move to GMX_PROLIG_MMPBSA)
        important_files = [
            "FINAL_RESULTS_MMPBSA.dat",
            "FINAL_RESULTS_MMPBSA.csv",
            "gmx_MMPBSA.log",
            "mmpbsa.in",
            "RESULTS_gmx_MMPBSA.h5"
        ]

        # Move important files to GMX_PROLIG_MMPBSA directory
        for filename in important_files:
            src = self.output_dir / filename
            if src.exists():
                dst = mmpbsa_dir / filename
                shutil.move(str(src), str(dst))
                print(f"  ✓ Saved: {filename}")

        # Files to delete (intermediate/temporary files)
        patterns_to_delete = [
            "_GMXMMPBSA_*",
            "COM.prmtop",
            "REC.prmtop",
            "LIG.prmtop",
            "COM_traj_*.xtc",
            "leap.log",
            "single_frame.xtc",
            "index.ndx"
        ]

        import glob
        deleted_count = 0
        for pattern in patterns_to_delete:
            for filepath in glob.glob(str(self.output_dir / pattern)):
                try:
                    os.remove(filepath)
                    deleted_count += 1
                except Exception as e:
                    print(f"  Warning: Could not delete {filepath}: {e}")

        print(f"  ✓ Cleaned up {deleted_count} intermediate files")
        print(f"  ✓ Results saved to: {mmpbsa_dir}")

        return mmpbsa_dir
