#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SwissParam force field generator wrapper for PRISM
Supports MMFF-based, MATCH, and hybrid MMFF-based-MATCH force fields
Uses SwissParam command-line API (port 8443) for reliable force field generation
"""

import os
import tarfile
import shutil
import tempfile
import hashlib
import subprocess

# Import the base class
try:
    from .base import ForceFieldGeneratorBase
except ImportError:
    from base import ForceFieldGeneratorBase


class SwissParamForceFieldGenerator(ForceFieldGeneratorBase):
    """Base class for SwissParam-based force field generators using API"""

    def __init__(self, ligand_path, output_dir, approach, overwrite=False):
        """
        Initialize SwissParam force field generator

        Parameters:
        -----------
        ligand_path : str
            Path to the ligand file (MOL2/SDF/PDB)
        output_dir : str
            Directory where output files will be stored
        approach : str
            SwissParam approach: "MMFF-based", "MATCH", or "MMFF-based-MATCH"
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(ligand_path, output_dir, overwrite)

        # Validate approach
        valid_approaches = ["mmff-based", "match", "both"]
        if approach not in valid_approaches:
            raise ValueError(f"Invalid approach '{approach}'. Must be one of {valid_approaches}")

        self.approach = approach

        # Output directory for SwissParam files
        self.lig_ff_dir = os.path.join(self.output_dir, self.get_output_dir_name())

        # SwissParam API endpoint
        self.api_base = "https://www.swissparam.ch:8443"

        print(f"\nInitialized SwissParam Force Field Generator:")
        print(f"  Ligand: {self.ligand_path}")
        print(f"  Approach: {self.approach}")
        print(f"  Output directory: {self.output_dir}")

    def get_output_dir_name(self):
        """Get the output directory name for this force field"""
        if self.approach == "mmff-based":
            return "LIG.mmff2gmx"
        elif self.approach == "match":
            return "LIG.match2gmx"
        else:  # both
            return "LIG.both2gmx"

    def run(self):
        """Run the SwissParam force field generation workflow"""
        print(f"\n{'='*60}")
        print(f"Starting SwissParam Force Field Generation ({self.approach})")
        print(f"{'='*60}")

        try:
            # Check if output already exists
            if os.path.exists(self.lig_ff_dir) and not self.overwrite:
                if self.check_required_files(self.lig_ff_dir):
                    print(f"\nUsing cached SwissParam force field parameters from: {self.lig_ff_dir}")
                    print("All required files found:")
                    for f in sorted(os.listdir(self.lig_ff_dir)):
                        print(f"  - {f}")
                    return self.lig_ff_dir

            # Upload ligand to SwissParam and get results
            tar_data = self._submit_to_swissparam_api()

            # Extract and organize results
            sp_dir = self._extract_tarball(tar_data)

            # Find and copy relevant files
            files = self._find_swissparam_files(sp_dir)

            # Copy files to output directory
            self._copy_output_files(files, sp_dir)

            print(f"\n{'='*60}")
            print(f"SwissParam Force Field Generation Complete!")
            print(f"Output directory: {self.lig_ff_dir}")
            print(f"{'='*60}")

            return self.lig_ff_dir

        except Exception as e:
            print(f"\nERROR: SwissParam force field generation failed: {e}")
            raise

    def _submit_to_swissparam_api(self):
        """
        Submit ligand to SwissParam using command-line API and retrieve results

        Returns:
        --------
        bytes
            The tarball data containing force field files
        """
        print("\n[1] Submitting to SwissParam API and retrieving results...")

        # Map approach to API parameter
        approach_map = {"mmff-based": "mmff-based", "match": "match", "both": "both"}
        api_approach = approach_map[self.approach]

        # Submit job and retrieve tarball directly
        # API returns tarball directly after processing
        # Execute curl in the ligand file's directory using relative path
        ligand_dir = os.path.dirname(os.path.abspath(self.ligand_path))
        ligand_file = os.path.basename(self.ligand_path)

        submit_cmd = ["curl", "-F", f"myMol2=@{ligand_file}", f"{self.api_base}/startparam?approach={api_approach}"]

        try:
            print(f'    Running: curl -F "myMol2=@{ligand_file}" "{self.api_base}/startparam?approach={api_approach}"')
            print(f"    Working directory: {ligand_dir}")
            result = subprocess.run(submit_cmd, capture_output=True, timeout=120, cwd=ligand_dir)

            if result.returncode != 0:
                raise RuntimeError(f"curl failed: {result.stderr.decode()}")

            tar_data = result.stdout

            # Check if we got valid tarball data
            if len(tar_data) < 1000:
                error_msg = tar_data.decode("utf-8", errors="ignore").strip()
                if error_msg:
                    raise RuntimeError(f"SwissParam server error: {error_msg}")
                else:
                    raise RuntimeError(f"Incomplete data ({len(tar_data)} bytes)")

            print(f"    ✓ Retrieved tarball ({len(tar_data)} bytes)")
            return tar_data

        except subprocess.TimeoutExpired:
            raise RuntimeError("Timeout waiting for API (120s)")
        except Exception as e:
            error_str = str(e)
            if "could not be submitted to the queuing system" in error_str:
                raise RuntimeError(
                    f"SwissParam server error: {e}\n"
                    "    NOTE: This may be due to:\n"
                    "    - Too many requests in a short time (rate limiting)\n"
                    "    - Temporary server issues\n"
                    "    Please wait a few minutes before trying again."
                )
            raise RuntimeError(f"Failed: {e}")

    def _extract_tarball(self, tar_data):
        """Extract tarball to temporary directory"""
        print("\n[2] Extracting force field files...")

        temp_dir = tempfile.mkdtemp(prefix="swissparam_")

        # Save tarball
        tar_hash = hashlib.md5(tar_data).hexdigest()[:8]
        tar_file = os.path.join(temp_dir, f"swissparam_{tar_hash}.tar.gz")

        with open(tar_file, "wb") as f:
            f.write(tar_data)

        # Extract
        extract_dir = os.path.join(temp_dir, "extracted")
        os.makedirs(extract_dir, exist_ok=True)

        with tarfile.open(tar_file, "r:gz") as tar:
            tar.extractall(extract_dir)

        print(f"    ✓ Extracted to: {extract_dir}")

        return extract_dir

    def _find_swissparam_files(self, sp_dir):
        """
        Find relevant force field files in SwissParam output

        Returns dict with file paths for different file types
        """
        print("\n[3] Locating force field files...")

        files = {}

        # Walk through extracted directory
        for root, dirs, filenames in os.walk(sp_dir):
            for filename in filenames:
                filepath = os.path.join(root, filename)

                # GROMACS topology file (.itp)
                if filename.endswith(".itp"):
                    if "posre" not in filename.lower() and "itp" not in files:
                        files["itp"] = filepath
                        print(f"    ✓ Found .itp: {filename}")

                # CHARMM parameter files
                elif filename.endswith(".par"):
                    if "par" not in files:
                        files["par"] = filepath
                        print(f"    ✓ Found .par: {filename}")

                elif filename.endswith(".rtf"):
                    if "rtf" not in files:
                        files["rtf"] = filepath
                        print(f"    ✓ Found .rtf: {filename}")

                elif filename.endswith(".prm"):
                    if "prm" not in files:
                        files["prm"] = filepath
                        print(f"    ✓ Found .prm: {filename}")

                # Coordinate files
                elif filename.endswith(".pdb"):
                    if "pdb" not in files:
                        files["pdb"] = filepath
                        print(f"    ✓ Found .pdb: {filename}")

                elif filename.endswith(".crd"):
                    if "crd" not in files:
                        files["crd"] = filepath
                        print(f"    ✓ Found .crd: {filename}")

                elif filename.endswith(".gro"):
                    if "gro" not in files:
                        files["gro"] = filepath
                        print(f"    ✓ Found .gro: {filename}")

        # Check if we got any force field files
        if not any(files.get(k) for k in ["itp", "par", "rtf", "prm"]):
            raise RuntimeError(
                "No force field files found in SwissParam output. "
                "The server may have failed to generate parameters. "
                "Please check the ligand structure and try again."
            )

        return files

    def _copy_output_files(self, files, sp_dir):
        """Copy SwissParam files to output directory"""
        print("\n[4] Copying files to output directory...")
        if sp_dir:
            print(f"    Source directory: {sp_dir}")

        # Create output directory
        os.makedirs(self.lig_ff_dir, exist_ok=True)

        # Copy all found files
        for file_type, filepath in files.items():
            if filepath and os.path.exists(filepath):
                dest_name = f"ligand.{file_type}"
                dest_path = os.path.join(self.lig_ff_dir, dest_name)
                shutil.copy2(filepath, dest_path)
                print(f"    ✓ Copied {file_type}: {dest_name}")

        print(f"    ✓ All files copied to: {self.lig_ff_dir}")

    def check_required_files(self, output_dir):
        """
        Check if required force field files exist

        For SwissParam, we need at least topology files (.itp for GROMACS or .par/.rtf for CHARMM)
        """
        has_itp = any(f.endswith(".itp") for f in os.listdir(output_dir))
        has_charmm = any(f.endswith((".par", ".rtf")) for f in os.listdir(output_dir))

        return has_itp or has_charmm


# Create specific generator classes for each approach
class MMFFForceFieldGenerator(SwissParamForceFieldGenerator):
    """MMFF-based force field generator using SwissParam"""

    def __init__(self, ligand_path, output_dir, overwrite=False):
        super().__init__(ligand_path, output_dir, "mmff-based", overwrite)


class MATCHForceFieldGenerator(SwissParamForceFieldGenerator):
    """MATCH force field generator using SwissParam"""

    def __init__(self, ligand_path, output_dir, overwrite=False):
        super().__init__(ligand_path, output_dir, "match", overwrite)


class BothForceFieldGenerator(SwissParamForceFieldGenerator):
    """Both MMFF-based + MATCH force field generator using SwissParam"""

    def __init__(self, ligand_path, output_dir, overwrite=False):
        super().__init__(ligand_path, output_dir, "both", overwrite)


# Alias for backward compatibility
HybridForceFieldGenerator = BothForceFieldGenerator
HybridMMFFMATCHForceFieldGenerator = BothForceFieldGenerator


# Convenience function for external use
def generate_swissparam_ff(ligand_path, output_dir, approach="mmff-based", overwrite=False):
    """
    Generate SwissParam force field for a ligand

    Parameters:
    -----------
    ligand_path : str
        Path to ligand MOL2 file
    output_dir : str
        Output directory for force field files
    approach : str
        SwissParam approach: "mmff-based", "match", or "both"
    overwrite : bool
        Whether to overwrite existing files

    Returns:
    --------
    str
        Path to the directory containing force field files
    """
    generator = SwissParamForceFieldGenerator(ligand_path, output_dir, approach, overwrite)
    return generator.run()
