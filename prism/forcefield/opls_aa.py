#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OPLS-AA force field generator wrapper for PRISM
Uses LigParGen server to generate OPLS-AA parameters
"""

import os
import sys
import re
import zipfile
import shutil
import tempfile
from pathlib import Path
import importlib.util

# Import the base class
try:
    from .base import ForceFieldGeneratorBase
except ImportError:
    from base import ForceFieldGeneratorBase

# Import required dependencies
try:
    import requests
    import urllib3

    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
except ImportError as e:
    print(f"Warning: Missing requests dependency: {e}")
    print("OPLS-AA force field generator will not be available.")
    print("INSTALLATION:")
    print("  pip install requests")


class OPLSAAForceFieldGenerator(ForceFieldGeneratorBase):
    """OPLS-AA force field generator using LigParGen"""

    def __init__(self, ligand_path, output_dir, charge=0, charge_model="cm1a", align_to_input=True, overwrite=False):
        """
        Initialize OPLS-AA force field generator

        Parameters:
        -----------
        ligand_path : str
            Path to the ligand file (MOL2/SDF)
        output_dir : str
            Directory where output files will be stored
        charge : int
            Molecule charge (default: 0)
        charge_model : str
            Charge model: "cm1a" (default) or "cm5"
        align_to_input : bool
            Whether to align output coordinates to input structure (default: True)
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(ligand_path, output_dir, overwrite)

        self.charge = charge
        self.charge_model = charge_model
        self.align_to_input = align_to_input

        # Check RDKit availability for alignment
        self.rdkit_available = self._check_rdkit()
        if align_to_input and not self.rdkit_available:
            print("Warning: RDKit not available, coordinate alignment disabled")
            self.align_to_input = False

        # Detect file format
        self.file_format = self._detect_file_format()

        # Output directory for OPLS-AA files
        self.lig_ff_dir = os.path.join(self.output_dir, self.get_output_dir_name())

        print(f"\nInitialized OPLS-AA Force Field Generator:")
        print(f"  Ligand: {self.ligand_path}")
        print(f"  Charge: {self.charge}")
        print(f"  Charge model: {self.charge_model}")
        print(f"  Align to input: {self.align_to_input}")
        print(f"  Output directory: {self.output_dir}")

    def get_output_dir_name(self):
        """Get the output directory name for OPLS-AA"""
        return "LIG.opls2gmx"

    def run(self):
        """Run the OPLS-AA force field generation workflow"""
        print(f"\n{'='*60}")
        print("Starting OPLS-AA Force Field Generation")
        print(f"{'='*60}")

        try:
            # Check if output already exists
            if os.path.exists(self.lig_ff_dir) and not self.overwrite:
                if self.check_required_files(self.lig_ff_dir):
                    print(f"\nUsing cached OPLS-AA force field parameters from: {self.lig_ff_dir}")
                    print("All required files found:")
                    for f in sorted(os.listdir(self.lig_ff_dir)):
                        print(f"  - {f}")
                    print("\n(Use --overwrite to regenerate)")
                    return self.lig_ff_dir

            # Create output directory
            os.makedirs(self.lig_ff_dir, exist_ok=True)

            # Upload to LigParGen and download results
            lpg_dir = self._upload_to_ligpargen()

            # Process LigParGen output
            self._process_ligpargen_output(lpg_dir)

            # Cleanup temporary files
            self._cleanup_temp_files(lpg_dir)

            print(f"\n{'='*60}")
            print("Force field generation completed successfully!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {self.lig_ff_dir}")
            print("\nGenerated files:")
            for f in sorted(os.listdir(self.lig_ff_dir)):
                print(f"  - {f}")

            return self.lig_ff_dir

        except Exception as e:
            print(f"\nError during force field generation: {e}")
            import traceback

            traceback.print_exc()
            raise

    def _check_rdkit(self):
        """Check if RDKit is available"""
        return importlib.util.find_spec("rdkit.Chem") is not None

    def _detect_file_format(self):
        """Detect input file format"""
        file_ext = os.path.splitext(self.ligand_path)[1].lower()

        if file_ext == ".mol2":
            return "mol2"
        elif file_ext in [".sdf", ".sd", ".mol"]:
            return "mol"
        elif file_ext == ".pdb":
            return "pdb"
        else:
            print(f"Warning: Unknown file extension '{file_ext}', assuming MOL format")
            return "mol"

    def _prepare_upload_file(self, temp_dir):
        """
        Prepare file for upload to LigParGen (convert to PDB if needed)

        Parameters:
        -----------
        temp_dir : str
            Temporary directory for file conversion

        Returns:
        --------
        str : Path to file ready for upload
        """
        # If already PDB, use directly
        if self.ligand_path.lower().endswith(".pdb"):
            return self.ligand_path

        # Try to convert to PDB using RDKit
        if not self.rdkit_available:
            print("Warning: RDKit not available for format conversion, uploading original file")
            return self.ligand_path

        try:
            from rdkit import Chem

            print(f"Converting {self.file_format.upper()} to PDB format for LigParGen...")

            # Load molecule
            if self.file_format == "mol2":
                mol = Chem.MolFromMol2File(self.ligand_path, removeHs=False)
            else:  # sdf/mol
                mol = Chem.MolFromMolFile(self.ligand_path, removeHs=False)

            if mol is None:
                print("Warning: Could not load molecule for conversion, uploading original file")
                return self.ligand_path

            # Write PDB
            pdb_file = os.path.join(temp_dir, f"{Path(self.ligand_path).stem}.pdb")
            Chem.MolToPDBFile(mol, pdb_file)

            print(f"Converted to PDB: {os.path.basename(pdb_file)}")
            return pdb_file

        except Exception as e:
            print(f"Warning: Format conversion failed ({e}), uploading original file")
            return self.ligand_path

    def _upload_to_ligpargen(self, max_retries=3):
        """
        Upload ligand to LigParGen server and download results with retry mechanism

        Parameters:
        -----------
        max_retries : int
            Maximum number of retry attempts (default: 3)

        Returns:
        --------
        str : Path to extracted results directory
        """
        print("\n=== Uploading to LigParGen Server ===")

        # Create temporary working directory
        temp_dir = tempfile.mkdtemp(prefix="ligpargen_")

        last_error = None
        for attempt in range(max_retries):
            try:
                if attempt > 0:
                    print(f"\nRetry attempt {attempt + 1}/{max_retries}...")
                    import time

                    time.sleep(2)  # Wait 2 seconds before retry

                # Convert to PDB format if needed (LigParGen prefers PDB)
                upload_file = self._prepare_upload_file(temp_dir)

                # Prepare the upload
                url = "https://traken.chem.yale.edu/cgi-bin/results_lpg.py"

                with open(upload_file, "rb") as f:
                    files = {"molpdbfile": (os.path.basename(upload_file), f, "chemical/x-pdb")}
                    data = {
                        "chargetype": self.charge_model,
                        "dropcharge": str(self.charge),
                        "checkopt": "0",  # No optimization
                    }

                    print(
                        f"Uploading {os.path.basename(upload_file)} with charge={self.charge}, model={self.charge_model}"
                    )
                    response = requests.post(url, files=files, data=data, verify=False, timeout=60)

                if response.status_code != 200:
                    last_error = f"Error submitting form. Response status: {response.status_code}"
                    print(f"⚠ {last_error}")
                    continue

                # Check for errors
                if "Sorry, an error has been detected in your input data" in response.content.decode():
                    last_error = "LigParGen detected an error in the input data"
                    print(f"⚠ {last_error}")
                    if attempt < max_retries - 1:
                        continue
                    else:
                        raise ValueError(f"{last_error}. Please check your ligand file or try a different format.")

                print("✓ Form submitted successfully")

                # Parse response to find ZIP file
                zip_match = re.search(r'value="(/tmp/UNK_\w+\.zip)"', response.text)
                if not zip_match:
                    # Try alternative pattern
                    zip_match = re.search(r"(/tmp/[\w_]+\.zip)", response.text)
                    if not zip_match:
                        last_error = "No ZIP file found in LigParGen response"
                        print(f"⚠ {last_error}")
                        if attempt < max_retries - 1:
                            # Debug: save response for inspection
                            debug_file = os.path.join(temp_dir, f"ligpargen_response_attempt{attempt+1}.html")
                            with open(debug_file, "w") as f:
                                f.write(response.text)
                            print(f"  Debug: Response saved to {debug_file}")
                            continue
                        else:
                            # Check for common error messages
                            if "error" in response.text.lower():
                                error_match = re.search(r"<p[^>]*>(.*?error.*?)</p>", response.text, re.IGNORECASE)
                                if error_match:
                                    raise ValueError(f"LigParGen error: {error_match.group(1)}")
                            raise ValueError(
                                "No ZIP file found after all retry attempts. Check debug files for details."
                            )

                zip_path = zip_match.group(1)
                print(f"✓ Found result ZIP: {os.path.basename(zip_path)}")

                # Download the ZIP file
                zip_url = "https://traken.chem.yale.edu/cgi-bin/download_lpg.py"
                zip_data = {"fileout": zip_path}
                zip_response = requests.post(zip_url, data=zip_data, verify=False, timeout=60)

                if zip_response.status_code != 200:
                    last_error = f"Error downloading results. Status: {zip_response.status_code}"
                    print(f"⚠ {last_error}")
                    if attempt < max_retries - 1:
                        continue
                    else:
                        raise ValueError(last_error)

                # Save and extract ZIP file
                zip_file = os.path.join(temp_dir, os.path.basename(zip_path))
                with open(zip_file, "wb") as f:
                    f.write(zip_response.content)

                print(f"✓ Downloaded ZIP file: {os.path.basename(zip_file)}")

                # Extract the ZIP file
                extract_dir = os.path.join(temp_dir, os.path.basename(zip_path).replace(".zip", ""))
                with zipfile.ZipFile(zip_file, "r") as zip_ref:
                    zip_ref.extractall(extract_dir)

                print(f"✓ Extracted to: {extract_dir}")

                # Success - return the result
                return extract_dir

            except Exception as e:
                last_error = str(e)
                print(f"⚠ Error during upload: {e}")
                if attempt < max_retries - 1:
                    continue
                else:
                    # Cleanup on final failure
                    shutil.rmtree(temp_dir, ignore_errors=True)
                    raise

        # If we get here, all retries failed
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise ValueError(f"Failed after {max_retries} attempts. Last error: {last_error}")

    def _process_ligpargen_output(self, lpg_dir):
        """
        Process LigParGen output and convert to PRISM standard format

        Parameters:
        -----------
        lpg_dir : str
            Path to extracted LigParGen results
        """
        print("\n=== Processing LigParGen Output ===")

        # Find the GROMACS files in LigParGen output
        lpg_files = self._find_ligpargen_files(lpg_dir)

        # Optionally align coordinates to input structure
        if self.align_to_input and self.rdkit_available:
            self._align_coordinates(lpg_files["gro"])

        # Standardize to LIG naming
        self._standardize_files(lpg_files)

    def _find_ligpargen_files(self, lpg_dir):
        """Find relevant files in LigParGen output"""
        files = {}

        # LigParGen typically puts files in a subdirectory
        for root, dirs, filenames in os.walk(lpg_dir):
            for filename in filenames:
                filepath = os.path.join(root, filename)

                if filename.endswith(".gro"):
                    files["gro"] = filepath
                elif filename.endswith(".itp"):
                    files["itp"] = filepath
                elif filename.endswith(".top"):
                    files["top"] = filepath

        if not files.get("gro") or not files.get("itp"):
            raise ValueError(f"Could not find required GROMACS files in {lpg_dir}")

        print(f"Found LigParGen files:")
        for key, path in files.items():
            print(f"  {key}: {os.path.basename(path)}")

        return files

    def _align_coordinates(self, gro_file):
        """
        Align LigParGen output coordinates to input structure

        Parameters:
        -----------
        gro_file : str
            Path to LigParGen GRO file
        """
        print("Aligning coordinates to input structure...")

        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolAlign

            # Load input molecule
            if self.file_format == "mol2":
                mol_input = Chem.MolFromMol2File(self.ligand_path, removeHs=False)
            elif self.file_format == "pdb":
                mol_input = Chem.MolFromPDBFile(self.ligand_path, removeHs=False)
            else:
                mol_input = Chem.MolFromMolFile(self.ligand_path, removeHs=False)

            if mol_input is None:
                print("Warning: Could not load input molecule for alignment")
                return

            # Load LigParGen output (convert GRO to PDB first)
            temp_pdb = gro_file.replace(".gro", "_temp.pdb")
            self._convert_gro_to_pdb(gro_file, temp_pdb)

            mol_lpg = Chem.MolFromPDBFile(temp_pdb, removeHs=False)

            if mol_lpg is None:
                print("Warning: Could not load LigParGen output for alignment")
                os.remove(temp_pdb)
                return

            # Perform alignment
            rmsd = rdMolAlign.AlignMol(mol_lpg, mol_input)
            print(f"Alignment RMSD: {rmsd:.3f} Å")

            # Write aligned structure back to PDB, then convert to GRO
            Chem.MolToPDBFile(mol_lpg, temp_pdb)
            self._convert_pdb_to_gro(temp_pdb, gro_file)

            os.remove(temp_pdb)
            print("Coordinate alignment completed")

        except Exception as e:
            print(f"Warning: Coordinate alignment failed: {e}")

    def _convert_gro_to_pdb(self, gro_file, pdb_file):
        """Convert GRO to PDB format"""
        # Simple conversion - just reformatting
        with open(gro_file, "r") as f:
            lines = f.readlines()

        with open(pdb_file, "w") as f:
            atom_idx = 1
            for i, line in enumerate(lines):
                if i < 2 or i >= len(lines) - 1:  # Skip header and box line
                    continue

                # GRO format parsing
                if len(line) < 44:
                    continue

                resnum = int(line[0:5].strip())
                resname = line[5:10].strip()
                atomname = line[10:15].strip()
                atomnum = int(line[15:20].strip())
                x = float(line[20:28].strip()) * 10  # nm to Å
                y = float(line[28:36].strip()) * 10
                z = float(line[36:44].strip()) * 10

                # Write PDB format
                f.write(
                    f"ATOM  {atom_idx:5d}  {atomname:<4s}{resname:>3s}  {resnum:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atomname[0]:>2s}\n"
                )
                atom_idx += 1

            f.write("END\n")

    def _convert_pdb_to_gro(self, pdb_file, gro_file):
        """Convert PDB back to GRO format"""
        from rdkit import Chem

        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            return

        conf = mol.GetConformer()

        # Read original GRO to get formatting
        with open(gro_file, "r") as f:
            gro_lines = f.readlines()

        # Update coordinates
        with open(gro_file, "w") as f:
            f.write(gro_lines[0])  # Title
            f.write(gro_lines[1])  # Atom count

            for i, line in enumerate(gro_lines[2:-1]):  # Atom lines
                if i >= mol.GetNumAtoms():
                    break

                pos = conf.GetAtomPosition(i)
                x = pos.x / 10.0  # Å to nm
                y = pos.y / 10.0
                z = pos.z / 10.0

                # Preserve GRO formatting but update coordinates
                f.write(f"{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}\n")

            f.write(gro_lines[-1])  # Box vectors

    def _standardize_files(self, lpg_files):
        """
        Standardize LigParGen files to PRISM format (LIG naming)

        Parameters:
        -----------
        lpg_files : dict
            Dictionary of file paths
        """
        print("\n=== Standardizing to LIG naming ===")

        # Copy and standardize GRO file
        if "gro" in lpg_files:
            self._standardize_gro(lpg_files["gro"], os.path.join(self.lig_ff_dir, "LIG.gro"))

        # Process ITP file
        if "itp" in lpg_files:
            self._standardize_itp(lpg_files["itp"], os.path.join(self.lig_ff_dir, "LIG.itp"))

        # Process TOP file or create from ITP
        if "top" in lpg_files:
            self._standardize_top(lpg_files["top"], os.path.join(self.lig_ff_dir, "LIG.top"))
        else:
            self._create_top_from_itp(os.path.join(self.lig_ff_dir, "LIG.itp"))

        # Extract atomtypes to separate file
        self._extract_atomtypes(os.path.join(self.lig_ff_dir, "LIG.itp"))

        # Create position restraints file
        self._create_position_restraints(os.path.join(self.lig_ff_dir, "LIG.gro"))

        # Normalize total charge to nearest integer
        self._normalize_charges()

        print("Standardization complete")

    def _standardize_gro(self, source_gro, target_gro):
        """Standardize GRO file to use LIG residue name"""
        with open(source_gro, "r") as f:
            lines = f.readlines()

        with open(target_gro, "w") as f:
            f.write("LIG system\n")  # Title

            for i, line in enumerate(lines):
                if i == 0:
                    continue  # Skip original title
                elif i >= 2 and i < len(lines) - 1:  # Atom lines
                    if len(line) > 10:
                        line = line[:5] + "LIG".ljust(5) + line[10:]

                f.write(line)

    def _standardize_itp(self, source_itp, target_itp):
        """Standardize ITP file to use LIG molecule name"""
        with open(source_itp, "r") as f:
            content = f.read()

        # Replace molecule names with LIG (case-insensitive for all variants)
        # Common LigParGen molecule names: UNK, MOL, LIG, or ligand name
        content = re.sub(r"\bUNK\b", "LIG", content, flags=re.IGNORECASE)
        content = re.sub(r"\bunk\b", "LIG", content)
        content = re.sub(r"\bMOL\b", "LIG", content, flags=re.IGNORECASE)

        # Also replace in [ moleculetype ] section specifically
        content = re.sub(r"(\[\s*moleculetype\s*\]\s*;\s*\w+\s+\w+\s+)\S+", r"\1LIG", content)

        # Replace residue names in [ atoms ] section
        # Pattern: atom_nr  type  resnr  residue  atom  cgnr  charge  mass
        # We need to replace the residue name (4th column)
        # Use regex to replace while preserving formatting
        lines = content.split("\n")
        new_lines = []
        in_atoms = False

        for line in lines:
            if line.strip().startswith("[ atoms ]"):
                in_atoms = True
                new_lines.append(line)
            elif in_atoms and line.strip().startswith("["):
                in_atoms = False
                new_lines.append(line)
            elif in_atoms and not line.strip().startswith(";") and line.strip():
                # Use regex to replace residue name (4th field) while preserving spacing
                # Typical format: "     1  opls_800      1    UNK     C1      1    -0.120  12.011"
                # Replace any word in position 4 with LIG
                modified_line = re.sub(r"^(\s*\d+\s+\S+\s+\d+\s+)\S+(\s+)", r"\1LIG\2", line)
                new_lines.append(modified_line)
            else:
                new_lines.append(line)

        content = "\n".join(new_lines)

        # Add position restraints include if not present
        if "#ifdef POSRES" not in content:
            content += '\n#ifdef POSRES\n#include "posre_LIG.itp"\n#endif\n'

        with open(target_itp, "w") as f:
            f.write(content)

    def _standardize_top(self, source_top, target_top):
        """Standardize TOP file"""
        with open(source_top, "r") as f:
            content = f.read()

        # Replace molecule names
        content = re.sub(r"\bUNK\b", "LIG", content)
        content = re.sub(r"\bunk\b", "LIG", content)

        # Update includes
        content = re.sub(r'#include\s+"[^"]*\.itp"', '#include "LIG.itp"', content)

        with open(target_top, "w") as f:
            f.write(content)

    def _create_top_from_itp(self, itp_file):
        """Create TOP file from ITP"""
        if not os.path.exists(itp_file):
            raise FileNotFoundError(f"ITP file not found: {itp_file}")
        top_file = os.path.join(self.lig_ff_dir, "LIG.top")

        with open(top_file, "w") as f:
            f.write("; OPLS-AA topology for LIG\n\n")
            f.write("[ defaults ]\n")
            f.write("; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n")
            f.write("1         3          yes        0.5      0.5\n\n")
            f.write('#include "atomtypes_LIG.itp"\n')
            f.write('#include "LIG.itp"\n\n')
            f.write("[ system ]\n")
            f.write("LIG system\n\n")
            f.write("[ molecules ]\n")
            f.write("LIG              1\n")

    def _extract_atomtypes(self, itp_file):
        """Extract atomtypes to separate file"""
        atomtypes_file = os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")

        with open(itp_file, "r") as f:
            lines = f.readlines()

        # Find atomtypes section
        atomtypes_lines = []
        in_atomtypes = False
        new_itp_lines = []

        for line in lines:
            if line.strip().startswith("[ atomtypes ]"):
                in_atomtypes = True
                atomtypes_lines.append(line)
            elif in_atomtypes and line.strip().startswith("["):
                in_atomtypes = False
                new_itp_lines.append(line)
            elif in_atomtypes:
                atomtypes_lines.append(line)
            else:
                new_itp_lines.append(line)

        # Write atomtypes file
        if atomtypes_lines:
            with open(atomtypes_file, "w") as f:
                f.writelines(atomtypes_lines)

            # Update ITP file (remove atomtypes section)
            with open(itp_file, "w") as f:
                f.writelines(new_itp_lines)

    def _create_position_restraints(self, gro_file):
        """Create position restraints file"""
        posre_file = os.path.join(self.lig_ff_dir, "posre_LIG.itp")

        # Count atoms from GRO file
        with open(gro_file, "r") as f:
            lines = f.readlines()

        try:
            atom_count = int(lines[1].strip())
        except (IndexError, ValueError):
            atom_count = len(lines) - 3  # Estimate

        # Write position restraints
        with open(posre_file, "w") as f:
            f.write("[ position_restraints ]\n")
            f.write("; atom  type      fx      fy      fz\n")
            for i in range(1, atom_count + 1):
                f.write(f"{i:>6}     1  1000  1000  1000\n")

    def _cleanup_temp_files(self, temp_dir):
        """Clean up temporary files"""
        print("\n=== Cleaning up ===")
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)
            print("Temporary files removed")

    def _normalize_charges(self):
        """
        校正OPLS配体的总电荷到最接近的整数。

        LigParGen生成的电荷总和不一定是整数（例如-0.0001）。
        本方法根据误差大小采用两种策略：
        1. 小误差（< 0.01）：将误差分配给电荷绝对值最大的原子
        2. 大误差（≥ 0.01）：按比例缩放所有原子电荷

        这样可以保持总电荷为整数，同时最小化对电荷分布的影响。
        """
        itp_path = os.path.join(self.lig_ff_dir, "LIG.itp")
        if not os.path.exists(itp_path):
            return

        # 解析ITP文件中的原子电荷
        atoms = []
        with open(itp_path, "r") as f:
            in_atom_section = False
            for line in f:
                if line.startswith("[ atoms ]"):
                    in_atom_section = True
                    continue
                if in_atom_section:
                    if line.startswith("["):
                        # 进入新的section，退出
                        break
                    line = line.strip()
                    if not line or line.startswith(";"):
                        continue
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            nr = int(parts[0])
                            atom_type = parts[1]
                            resnr = int(parts[2])
                            residue = parts[3]
                            atom_name = parts[4]
                            cgnr = int(parts[5])
                            charge = float(parts[6])
                            mass = float(parts[7])
                            atoms.append(
                                {
                                    "nr": nr,
                                    "type": atom_type,
                                    "resnr": resnr,
                                    "residue": residue,
                                    "name": atom_name,
                                    "cgnr": cgnr,
                                    "charge": charge,
                                    "mass": mass,
                                    "line_start": line,
                                }
                            )
                        except (ValueError, IndexError) as e:
                            # 解析失败，跳过此行
                            continue

        if not atoms:
            print("  Warning: No atoms found in ITP file, skipping charge normalization")
            return

        # 计算当前总电荷
        current_total = sum(atom["charge"] for atom in atoms)

        # 确定期望电荷（四舍五入到最接近的整数）
        target_charge = round(current_total)

        # 如果总电荷已经是整数（误差小于1e-10），不需要校正
        if abs(current_total - target_charge) < 1e-10:
            print(f"  Total charge is already integer: {current_total:.6f}")
            return

        # 计算需要校正的量
        charge_error = target_charge - current_total
        print(f"  Normalizing charges: {current_total:.6f} → {target_charge:.6f} (Δ = {charge_error:+.6f})")

        # 根据误差大小选择策略
        if abs(charge_error) < 0.01:
            # 策略1：小误差，分配给电荷绝对值最大的原子
            # 这样可以最小化对整体电荷分布的影响
            max_charge_atom = max(atoms, key=lambda a: abs(a["charge"]))
            max_charge_atom["charge"] += charge_error
            print(
                f"  Applied small correction to atom {max_charge_atom['name']} "
                f"(charge: {max_charge_atom['charge'] - charge_error:.6f} → "
                f"{max_charge_atom['charge']:.6f})"
            )
        else:
            # 策略2：大误差，按比例缩放所有原子电荷
            # 这样保持电荷分布的相对关系
            scale_factor = target_charge / current_total
            for atom in atoms:
                atom["charge"] *= scale_factor
            print(f"  Applied proportional scaling (factor: {scale_factor:.6f}) to all atoms")

        # 重写ITP文件
        with open(itp_path, "r") as f:
            lines = f.readlines()

        # 替换原子行
        atom_dict = {atom["nr"]: atom for atom in atoms}
        new_lines = []
        in_atom_section = False

        for line in lines:
            if line.strip().startswith("[ atoms ]"):
                in_atom_section = True
                new_lines.append(line)
            elif in_atom_section:
                if line.strip().startswith("["):
                    # 进入新的section
                    in_atom_section = False
                    new_lines.append(line)
                elif line.strip() and not line.strip().startswith(";"):
                    # 这是原子行，尝试替换
                    parts = line.split()
                    if len(parts) >= 8:
                        try:
                            nr = int(parts[0])
                            if nr in atom_dict:
                                atom = atom_dict[nr]
                                # 重建行，保持原始格式
                                new_line = (
                                    f"{atom['nr']:>6}  "
                                    f"{atom['type']:<10}  "
                                    f"{atom['resnr']:>4}  "
                                    f"{atom['residue']:<5}  "
                                    f"{atom['name']:<5}  "
                                    f"{atom['cgnr']:>4}  "
                                    f"{atom['charge']:>8.6f}  "
                                    f"{atom['mass']:>8.3f}\n"
                                )
                                new_lines.append(new_line)
                                continue
                        except (ValueError, IndexError):
                            pass
                    # 保留原始行
                    new_lines.append(line)
            else:
                new_lines.append(line)

        # 写回文件
        with open(itp_path, "w") as f:
            f.writelines(new_lines)

        # 验证校正后的总电荷
        verified_total = sum(atom["charge"] for atom in atoms)
        print(f"  ✓ Verified total charge: {verified_total:.6f}")

        # 检查是否成功校正到整数
        if abs(verified_total - target_charge) > 1e-6:
            print(
                f"  ⚠ Warning: Charge normalization may have issues "
                f"(target: {target_charge:.6f}, actual: {verified_total:.6f})"
            )


if __name__ == "__main__":
    # Example usage
    import sys

    if len(sys.argv) < 2:
        print("Usage: python opls_aa.py <ligand_file.mol2>")
        sys.exit(1)

    ligand_file = sys.argv[1]
    output_dir = os.path.dirname(ligand_file) if os.path.dirname(ligand_file) else "."

    generator = OPLSAAForceFieldGenerator(
        ligand_path=ligand_file, output_dir=output_dir, charge=0, charge_model="cm1a", align_to_input=True
    )

    result_dir = generator.run()
    print(f"\nSuccess! Output in: {result_dir}")
