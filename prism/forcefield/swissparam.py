#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SwissParam force field generator wrapper for PRISM
Supports MMFF-based, MATCH, and hybrid MMFF-based-MATCH force fields
Uses SwissParam command-line API (port 8443) for reliable force field generation

Enhanced to parse RTF files and generate PRISM-compatible output formats
"""

import os
import tarfile
import shutil
import tempfile
import hashlib
import subprocess
from typing import Dict, Tuple

# Import the base class
try:
    from .base import ForceFieldGeneratorBase
    from .common.units import angstrom_to_nm
except ImportError:
    from base import ForceFieldGeneratorBase
    from common.units import angstrom_to_nm


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

        # Molecule type name for PRISM compatibility
        self.moleculetype_name = "LIG"

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

            # Generate PRISM format files if needed
            self._generate_prism_formats(files)

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

            # Small responses can still be valid gzip-compressed tarballs.
            # Only interpret small payloads as text when they are not gzip data.
            if len(tar_data) < 1000 and not tar_data.startswith(b"\x1f\x8b"):
                response_text = tar_data.decode("utf-8", errors="ignore").strip()

                # Check if this is a session number response
                if "Session number:" in response_text:
                    # Extract session number
                    import re

                    match = re.search(r"Session number:\s*(\d+)", response_text)
                    if match:
                        session_number = match.group(1)
                        print(f"    ⏳ Job queued, session number: {session_number}")

                        # Poll for completion
                        tar_data = self._poll_swissparam_session(session_number, ligand_dir=ligand_dir)
                    else:
                        raise RuntimeError(f"Could not parse session number from: {response_text}")
                else:
                    # Other error
                    if response_text:
                        raise RuntimeError(f"SwissParam server error: {response_text}")
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

    def _poll_swissparam_session(
        self, session_number: str, ligand_dir: str, timeout: int = 300, poll_interval: int = 10
    ) -> bytes:
        """
        Poll SwissParam session until completion and retrieve results.

        Parameters
        ----------
        session_number : str
            Session number returned by SwissParam
        ligand_dir : str
            Working directory for curl commands
        timeout : int
            Maximum time to wait in seconds (default: 300s = 5min)
        poll_interval : int
            Time between polls in seconds (default: 10s)

        Returns
        -------
        bytes
            Tarball data from completed session
        """
        import time

        print(f"\n[1b] Polling SwissParam session {session_number}...")
        print(f"    Timeout: {timeout}s, Poll interval: {poll_interval}s")

        start_time = time.time()
        attempts = 0

        while time.time() - start_time < timeout:
            attempts += 1

            # Check session status
            check_cmd = ["curl", "-s", f"https://www.swissparam.ch:8443/checksession?sessionNumber={session_number}"]

            try:
                result = subprocess.run(check_cmd, capture_output=True, timeout=30, cwd=ligand_dir)

                status_text = result.stdout.decode("utf-8", errors="ignore")

                # Check if calculation is finished
                if "Calculation is finished" in status_text:
                    print(f"    ✓ Session completed after {attempts} attempts ({time.time() - start_time:.1f}s)")

                    # Retrieve results
                    retrieve_cmd = [
                        "curl",
                        "-s",
                        f"https://www.swissparam.ch:8443/retrievesession?sessionNumber={session_number}",
                    ]

                    result = subprocess.run(retrieve_cmd, capture_output=True, timeout=120, cwd=ligand_dir)

                    tar_data = result.stdout

                    if len(tar_data) < 1000 and not tar_data.startswith(b"\x1f\x8b"):
                        raise RuntimeError(f"Retrieved incomplete data ({len(tar_data)} bytes)")

                    print(f"    ✓ Retrieved tarball ({len(tar_data)} bytes)")
                    return tar_data

                elif "still running" in status_text or "queued" in status_text:
                    print(f"    Attempt {attempts}: Job still running, waiting {poll_interval}s...")
                    time.sleep(poll_interval)

                else:
                    # Unknown status
                    print(f"    Status: {status_text[:200]}")
                    time.sleep(poll_interval)

            except subprocess.TimeoutExpired:
                print(f"    Attempt {attempts}: Timeout checking status, retrying...")
                time.sleep(poll_interval)

        raise RuntimeError(
            f"Timeout after {timeout}s waiting for session {session_number}\n"
            f"    You can manually check status with:\n"
            f"    curl 'https://www.swissparam.ch:8443/checksession?sessionNumber={session_number}'"
        )

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

        For PRISM compatibility, check for standard PRISM output files
        """
        required_files = ["LIG.gro", "LIG.itp", "LIG.top", "atomtypes_LIG.itp", "posre_LIG.itp"]

        # Check for PRISM standard files
        has_prism_files = all(os.path.exists(os.path.join(output_dir, f)) for f in required_files)

        if has_prism_files:
            return True

        # Fall back to checking for raw SwissParam files that can be converted
        has_itp = any(f.endswith(".itp") for f in os.listdir(output_dir))
        has_rtf = any(f.endswith(".rtf") for f in os.listdir(output_dir))
        has_pdb = any(f.endswith(".pdb") for f in os.listdir(output_dir))

        return has_itp or (has_rtf and has_pdb)

    def _generate_prism_formats(self, files):
        """
        Generate PRISM format files from downloaded SwissParam files

        Parameters:
        -----------
        files : dict
            Dictionary of found file paths
        """
        print("\n[5] Checking PRISM format files...")

        # Check if PRISM files already exist
        prism_files = ["LIG.gro", "LIG.itp", "LIG.top", "atomtypes_LIG.itp", "posre_LIG.itp"]
        missing_files = [f for f in prism_files if not os.path.exists(os.path.join(self.lig_ff_dir, f))]

        if not missing_files:
            print("    ✓ All PRISM format files already present")
            return

        print(f"    Generating missing PRISM files: {', '.join(missing_files)}")

        # Check if we have RTF files to parse
        rtf_file = files.get("rtf")
        pdb_file = files.get("pdb")

        if rtf_file and pdb_file:
            print("    Found RTF and PDB files - generating PRISM formats...")
            par_file = files.get("par")
            self._parse_rtf_and_generate_prism(rtf_file, pdb_file, par_file)
        elif files.get("itp") and files.get("gro"):
            print("    Using existing SwissParam ITP and GRO files")
            self._copy_swissparam_itp_to_prism(files)
        else:
            print("    ⚠ Warning: Cannot generate PRISM formats - missing required files")
            print(f"    Available files: {list(files.keys())}")

    def _parse_rtf_and_generate_prism(self, rtf_file, pdb_file, par_file=None):
        """
        Parse RTF file and generate PRISM format output

        Parameters:
        -----------
        rtf_file : str
            Path to RTF file
        pdb_file : str
            Path to PDB file
        par_file : str, optional
            Path to PAR file with parameters
        """
        print("\n[5a] Parsing RTF file...")

        try:
            # Parse RTF file
            rtf_data = self._parse_rtf_file(rtf_file)

            # Parse PAR file if available
            par_data = None
            if par_file:
                print("\n[5a'] Parsing PAR file for parameters...")
                par_data = self._parse_par_file(par_file)

            # Parse PDB coordinates
            coordinates = self._parse_pdb_coordinates(pdb_file)

            # Generate PRISM files
            print("\n[5b] Generating PRISM format files...")

            # Generate GRO file
            if not os.path.exists(os.path.join(self.lig_ff_dir, "LIG.gro")):
                self._generate_gro_file(coordinates, rtf_data)

            # Generate ITP file with parameters
            if not os.path.exists(os.path.join(self.lig_ff_dir, "LIG.itp")):
                self._generate_itp_file(rtf_data, par_data)

            # Generate atomtypes file
            if not os.path.exists(os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")):
                self._generate_atomtypes_file(rtf_data)

            # Generate position restraints
            if not os.path.exists(os.path.join(self.lig_ff_dir, "posre_LIG.itp")):
                self._generate_posre_file(rtf_data)

            # Generate topology file
            if not os.path.exists(os.path.join(self.lig_ff_dir, "LIG.top")):
                self._generate_top_file()

            print("    ✓ All PRISM format files generated successfully")

        except Exception as e:
            print(f"    ✗ Error generating PRISM formats: {e}")
            raise

    def _parse_par_file(self, par_file: str) -> Dict:
        """
        Parse PAR file for bond, angle, dihedral parameters

        Parameters:
        -----------
        par_file : str
            Path to PAR file

        Returns:
        --------
        dict : Parsed parameters with CHARMM units
        """
        print(f"      Reading PAR file: {os.path.basename(par_file)}")

        params = {"bonds": {}, "angles": {}, "dihedrals": {}}

        with open(par_file, "r") as f:
            lines = f.readlines()

        section = None
        for line in lines:
            line = line.strip()
            if not line or line.startswith("*"):
                continue

            # Section detection
            if line.startswith("BONDS"):
                section = "bonds"
                continue
            elif line.startswith("ANGLES"):
                section = "angles"
                continue
            elif line.startswith("DIHEDRALS"):
                section = "dihedrals"
                continue
            elif line.startswith("IMPROPER"):
                section = "impropers"
                continue

            # Parse parameters
            if section == "bonds":
                parts = line.split()
                if len(parts) >= 4:
                    type1, type2, k, b0 = parts[0], parts[1], float(parts[2]), float(parts[3])
                    key = tuple(sorted([type1, type2]))
                    params["bonds"][key] = {"k": k, "b0": b0}

            elif section == "angles":
                parts = line.split()
                if len(parts) >= 5:
                    type1, type2, type3, k, theta0 = parts[0], parts[1], parts[2], float(parts[3]), float(parts[4])
                    key = (type1, type2, type3)
                    params["angles"][key] = {"k": k, "theta0": theta0}

            elif section == "dihedrals":
                parts = line.split()
                if len(parts) >= 6:
                    type1, type2, type3, type4 = parts[0], parts[1], parts[2], parts[3]
                    k, delta = float(parts[4]), float(parts[5])
                    # Check for multiplicity (some formats have it)
                    mult = int(float(parts[6])) if len(parts) > 6 else 0
                    key = (type1, type2, type3, type4)
                    if key not in params["dihedrals"]:
                        params["dihedrals"][key] = []
                    params["dihedrals"][key].append({"k": k, "delta": delta, "mult": mult})

        print(
            f"        Found {len(params['bonds'])} bonds, {len(params['angles'])} angles, {len(params['dihedrals'])} dihedrals"
        )
        return params

    def _parse_rtf_file(self, rtf_file) -> Dict:
        """
        Parse RTF file for atoms, bonds, angles, dihedrals

        Parameters:
        -----------
        rtf_file : str
            Path to RTF file

        Returns:
        --------
        dict : Parsed RTF data
        """
        print(f"      Reading RTF file: {os.path.basename(rtf_file)}")

        rtf_data = {"atoms": {}, "bonds": [], "angles": [], "dihedrals": [], "atom_types": {}}

        with open(rtf_file, "r") as f:
            lines = f.readlines()

        in_atoms = False
        in_bonds = False
        in_angles = False
        in_dihedrals = False

        for line in lines:
            line = line.strip()
            if not line or line.startswith("*"):
                continue

            # Section detection
            if line.startswith("ATOM"):
                in_atoms = True
                in_bonds = in_angles = in_dihedrals = False
                parts = line.split()
                if len(parts) >= 4:
                    name = parts[1]
                    atom_type = parts[2]
                    charge = float(parts[3])
                    rtf_data["atoms"][name] = {"type": atom_type, "charge": charge}
                    rtf_data["atom_types"][name] = atom_type
                continue

            elif line.startswith("BOND"):
                in_bonds = True
                in_atoms = in_angles = in_dihedrals = False
                parts = line.split()[1:]
                for i in range(0, len(parts), 2):
                    if i + 1 < len(parts):
                        rtf_data["bonds"].append((parts[i], parts[i + 1]))
                continue

            elif line.startswith("ANGL") or line.startswith("ANGLE"):
                in_angles = True
                in_atoms = in_bonds = in_dihedrals = False
                parts = line.split()[1:]
                for i in range(0, len(parts), 3):
                    if i + 2 < len(parts):
                        rtf_data["angles"].append((parts[i], parts[i + 1], parts[i + 2]))
                continue

            elif line.startswith("IMPR") or line.startswith("IMPRO"):
                # Improper dihedrals - not parsed into bonds
                in_atoms = in_bonds = in_angles = in_dihedrals = False
                continue

            elif line.startswith("DIHE") or line.startswith("DIHEDRAL"):
                in_dihedrals = True
                in_atoms = in_bonds = in_angles = False
                parts = line.split()[1:]
                for i in range(0, len(parts), 4):
                    if i + 3 < len(parts):
                        rtf_data["dihedrals"].append((parts[i], parts[i + 1], parts[i + 2], parts[i + 3]))
                continue

            # Continue parsing multi-line sections
            if in_bonds and line:
                parts = line.split()
                for i in range(0, len(parts), 2):
                    if i + 1 < len(parts):
                        rtf_data["bonds"].append((parts[i], parts[i + 1]))

            elif in_angles and line:
                parts = line.split()
                for i in range(0, len(parts), 3):
                    if i + 2 < len(parts):
                        rtf_data["angles"].append((parts[i], parts[i + 1], parts[i + 2]))

            elif in_dihedrals and line:
                parts = line.split()
                for i in range(0, len(parts), 4):
                    if i + 3 < len(parts):
                        rtf_data["dihedrals"].append((parts[i], parts[i + 1], parts[i + 2], parts[i + 3]))

        print(f"      - Found {len(rtf_data['atoms'])} atoms")
        print(f"      - Found {len(rtf_data['bonds'])} bonds")
        print(f"      - Found {len(rtf_data['angles'])} angles")
        print(f"      - Found {len(rtf_data['dihedrals'])} dihedrals")

        return rtf_data

    def _parse_pdb_coordinates(self, pdb_file) -> Dict[str, Tuple[float, float, float]]:
        """
        Parse PDB file for coordinates

        Parameters:
        -----------
        pdb_file : str
            Path to PDB file

        Returns:
        --------
        dict : Atom name to coordinates mapping
        """
        print(f"      Reading PDB file: {os.path.basename(pdb_file)}")

        coordinates = {}
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        # Extract atom name (columns 13-16, 1-indexed)
                        atom_name = line[12:16].strip()

                        # Try standard PDB format first (columns 31-54, 1-indexed)
                        if len(line) >= 54:
                            x_str = line[30:38].strip()
                            y_str = line[38:46].strip()
                            z_str = line[46:54].strip()

                            # Check if we got valid numbers
                            if x_str and y_str and z_str:
                                try:
                                    x = float(x_str)
                                    y = float(y_str)
                                    z = float(z_str)
                                except ValueError:
                                    # Fall back to whitespace splitting
                                    raise ValueError("Invalid coordinate format")
                            else:
                                raise ValueError("Empty coordinate fields")
                        else:
                            raise ValueError("Line too short")

                        coordinates[atom_name] = (x, y, z)

                    except (ValueError, IndexError):
                        # Fallback: split by whitespace for non-standard formats
                        try:
                            parts = line.split()
                            if len(parts) >= 8:  # ATOM serial name x y z ...
                                atom_name = parts[2]  # Usually third column
                                x = float(parts[-6])  # 6th from end
                                y = float(parts[-5])  # 5th from end
                                z = float(parts[-4])  # 4th from end
                                coordinates[atom_name] = (x, y, z)
                        except (ValueError, IndexError):
                            # Skip malformed lines
                            continue

        print(f"      - Found coordinates for {len(coordinates)} atoms")
        return coordinates

    def _parse_prm_file(self, prm_file) -> Dict:
        """
        Parse PRM file for force field parameters

        Parameters:
        -----------
        prm_file : str
            Path to PRM file

        Returns:
        --------
        dict : Parsed PRM data
        """
        print(f"      Reading PRM file: {os.path.basename(prm_file)}")

        prm_data = {}
        # Basic PRM parsing - can be extended as needed
        with open(prm_file, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("*"):
                    # Parse mass, bond, angle, dihedral parameters as needed
                    pass

        return prm_data

    def _generate_gro_file(self, coordinates: Dict[str, Tuple[float, float, float]], rtf_data: Dict):
        """Generate GROMACS GRO file"""
        gro_file = os.path.join(self.lig_ff_dir, "LIG.gro")
        print(f"      Generating GRO file...")

        # Ensure output directory exists
        os.makedirs(self.lig_ff_dir, exist_ok=True)

        with open(gro_file, "w") as f:
            f.write("Generated by SwissParamForceFieldGenerator\n")
            f.write(f"{len(rtf_data['atoms']):5d}\n")

            atom_id = 0
            for atom_name in rtf_data["atoms"]:
                if atom_name in coordinates:
                    atom_id += 1
                    x, y, z = coordinates[atom_name]
                    # Convert Å to nm
                    x_nm = angstrom_to_nm(x)
                    y_nm = angstrom_to_nm(y)
                    z_nm = angstrom_to_nm(z)

                    f.write(f"{1:5d}{self.moleculetype_name:<5s}{atom_name:>5s}{atom_id:5d}")
                    f.write(f"{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n")

            # Box dimensions (dummy)
            f.write("   5.00000   5.00000   5.00000\n")

        print(f"      ✓ Generated: {os.path.basename(gro_file)}")

    def _parse_original_swissparam_itp(self) -> Dict[str, Dict]:
        """
        Parse the original SwissParam ligand.itp copied into the ligand FF directory.

        Returns
        -------
        Dict[str, Dict]
            Contains atomtypes and per-atom records keyed by atom name.
        """
        source_itp = os.path.join(self.lig_ff_dir, "ligand.itp")
        parsed = {"atomtypes": {}, "atoms_by_name": {}}

        if not os.path.exists(source_itp):
            return parsed

        section = None
        with open(source_itp, "r") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line or line.startswith(";"):
                    continue
                if line.startswith("[") and line.endswith("]"):
                    section = line.strip("[] ").lower()
                    continue

                parts = line.split()
                if section == "atomtypes" and len(parts) >= 7:
                    parsed["atomtypes"][parts[0]] = {
                        "mass": parts[2],
                        "charge": parts[3],
                        "ptype": parts[4],
                        "sigma": parts[5],
                        "epsilon": parts[6],
                    }
                elif section == "atoms" and len(parts) >= 8:
                    parsed["atoms_by_name"][parts[4]] = {
                        "type": parts[1],
                        "mass": parts[7],
                    }

        return parsed

    def _generate_itp_file(self, rtf_data: Dict, par_data: Dict = None):
        """Generate GROMACS ITP file with parameters"""
        itp_file = os.path.join(self.lig_ff_dir, "LIG.itp")
        print(f"      Generating ITP file...")

        # Ensure output directory exists
        os.makedirs(self.lig_ff_dir, exist_ok=True)
        original_itp = self._parse_original_swissparam_itp()

        with open(itp_file, "w") as f:
            f.write("; Generated by SwissParamForceFieldGenerator\n")
            f.write(f"; Approach: {self.approach}\n")
            f.write("\n")

            # Moleculetype
            f.write("[ moleculetype ]\n")
            f.write(f"{self.moleculetype_name}    3\n\n")

            # Atoms
            f.write("[ atoms ]\n")
            f.write(";   nr       type  resnr residue  atom   cgnr     charge       mass\n")

            atom_id = 0
            for atom_name, atom_data in rtf_data["atoms"].items():
                atom_id += 1
                atom_type = atom_data["type"]
                charge = atom_data["charge"]
                mass = original_itp["atoms_by_name"].get(atom_name, {}).get("mass", "12.01000")
                f.write(f"{atom_id:5d} {atom_type:>10s}      1    {self.moleculetype_name}   {atom_name:>4s}")
                f.write(f"{atom_id:5d}  {charge:10.6f}   {float(mass):10.5f}\n")

            name_to_id = {name: i + 1 for i, name in enumerate(rtf_data["atoms"].keys())}

            # Bonds with parameters
            if rtf_data["bonds"]:
                f.write("\n[ bonds ]\n")
                f.write(";  ai    aj funct   b0       kb\n")

                for atom1, atom2 in rtf_data["bonds"]:
                    if atom1 in name_to_id and atom2 in name_to_id:
                        id1 = name_to_id[atom1]
                        id2 = name_to_id[atom2]

                        # Look up parameters if available
                        if par_data and "atom_types" in rtf_data:
                            type1 = rtf_data["atom_types"].get(atom1)
                            type2 = rtf_data["atom_types"].get(atom2)
                            if type1 and type2:
                                key = tuple(sorted([type1, type2]))
                                if key in par_data.get("bonds", {}):
                                    # Convert CHARMM to GROMACS units
                                    # K: kcal/mol/Å² → kJ/mol/nm² = * 418.4
                                    # b0: Å → nm = * 0.1
                                    k_charmm = par_data["bonds"][key]["k"]
                                    b0_charmm = par_data["bonds"][key]["b0"]
                                    k_gromacs = k_charmm * 418.4
                                    b0_gromacs = b0_charmm * 0.1
                                    f.write(f"{id1:5d} {id2:5d}     1   {b0_gromacs:.6f}   {k_gromacs:.3f}\n")
                                else:
                                    # No parameters found, write without params
                                    f.write(f"{id1:5d} {id2:5d}     1\n")
                            else:
                                f.write(f"{id1:5d} {id2:5d}     1\n")
                        else:
                            f.write(f"{id1:5d} {id2:5d}     1\n")

            # Angles with parameters
            if rtf_data["angles"]:
                f.write("\n[ angles ]\n")
                f.write(";  ai    aj    ak funct   theta0   ktheta\n")

                for atom1, atom2, atom3 in rtf_data["angles"]:
                    if all(atom in name_to_id for atom in [atom1, atom2, atom3]):
                        id1 = name_to_id[atom1]
                        id2 = name_to_id[atom2]
                        id3 = name_to_id[atom3]

                        # Look up parameters if available
                        if par_data and "atom_types" in rtf_data:
                            type1 = rtf_data["atom_types"].get(atom1)
                            type2 = rtf_data["atom_types"].get(atom2)
                            type3 = rtf_data["atom_types"].get(atom3)
                            if type1 and type2 and type3:
                                key = (type1, type2, type3)
                                if key in par_data.get("angles", {}):
                                    # Convert CHARMM to GROMACS units
                                    # K: kcal/mol/rad² → kJ/mol/rad² = * 4.184
                                    # theta0: degrees (no conversion)
                                    k_charmm = par_data["angles"][key]["k"]
                                    theta0_charmm = par_data["angles"][key]["theta0"]
                                    k_gromacs = k_charmm * 4.184
                                    f.write(
                                        f"{id1:5d} {id2:5d} {id3:5d}     1   {theta0_charmm:.2f}   {k_gromacs:.3f}\n"
                                    )
                                else:
                                    f.write(f"{id1:5d} {id2:5d} {id3:5d}     1\n")
                            else:
                                f.write(f"{id1:5d} {id2:5d} {id3:5d}     1\n")
                        else:
                            f.write(f"{id1:5d} {id2:5d} {id3:5d}     1\n")

            # Dihedrals with parameters
            if rtf_data["dihedrals"]:
                f.write("\n[ dihedrals ]\n")
                f.write(";  ai    aj    ak    al funct   delta     kd\n")

                for atom1, atom2, atom3, atom4 in rtf_data["dihedrals"]:
                    if all(atom in name_to_id for atom in [atom1, atom2, atom3, atom4]):
                        id1 = name_to_id[atom1]
                        id2 = name_to_id[atom2]
                        id3 = name_to_id[atom3]
                        id4 = name_to_id[atom4]

                        # Look up parameters if available
                        if par_data and "atom_types" in rtf_data:
                            type1 = rtf_data["atom_types"].get(atom1)
                            type2 = rtf_data["atom_types"].get(atom2)
                            type3 = rtf_data["atom_types"].get(atom3)
                            type4 = rtf_data["atom_types"].get(atom4)
                            if type1 and type2 and type3 and type4:
                                key = (type1, type2, type3, type4)
                                if key in par_data.get("dihedrals", {}):
                                    # Get first dihedral parameter (may have multiple)
                                    dih_params = par_data["dihedrals"][key][0]
                                    # Convert CHARMM to GROMACS units
                                    # K: kcal/mol → kJ/mol = * 4.184
                                    # delta: degrees (no conversion)
                                    k_charmm = dih_params["k"]
                                    delta_charmm = dih_params["delta"]
                                    k_gromacs = k_charmm * 4.184
                                    # Use funct=9 for proper dihedral (Ryckaert-Bellemans)
                                    f.write(
                                        f"{id1:5d} {id2:5d} {id3:5d} {id4:5d}     9   {delta_charmm:.2f}   {k_gromacs:.3f}\n"
                                    )
                                else:
                                    f.write(f"{id1:5d} {id2:5d} {id3:5d} {id4:5d}     9\n")
                            else:
                                f.write(f"{id1:5d} {id2:5d} {id3:5d} {id4:5d}     9\n")
                        else:
                            f.write(f"{id1:5d} {id2:5d} {id3:5d} {id4:5d}     9\n")

        print(f"      ✓ Generated: {os.path.basename(itp_file)}")

    def _generate_atomtypes_file(self, rtf_data: Dict):
        """Generate atomtypes file"""
        atomtypes_file = os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")
        print(f"      Generating atomtypes file...")

        # Ensure output directory exists
        os.makedirs(self.lig_ff_dir, exist_ok=True)

        unique_types = set(atom_data["type"] for atom_data in rtf_data["atoms"].values())
        original_itp = self._parse_original_swissparam_itp()
        original_atomtypes = original_itp["atomtypes"]
        original_atoms_by_name = original_itp["atoms_by_name"]

        with open(atomtypes_file, "w") as f:
            f.write("; Generated by SwissParamForceFieldGenerator\n")
            f.write("[ atomtypes ]\n")
            f.write(";name   bond_type     mass     charge   ptype   sigma         epsilon\n")

            for atom_type in sorted(unique_types):
                mapped_params = None
                for atom_name, atom_data in rtf_data["atoms"].items():
                    if atom_data["type"] != atom_type:
                        continue
                    original_atom = original_atoms_by_name.get(atom_name)
                    if not original_atom:
                        continue
                    mapped_params = original_atomtypes.get(original_atom["type"])
                    if mapped_params:
                        break

                if mapped_params:
                    f.write(
                        f"{atom_type:>6s}  {atom_type:>6s}  "
                        f"{float(mapped_params['mass']):8.5f}  {float(mapped_params['charge']):8.5f}  "
                        f"{mapped_params['ptype']:>1s}   {float(mapped_params['sigma']):.6e}  {float(mapped_params['epsilon']):.6e}\n"
                    )
                else:
                    # Fallback for incomplete SwissParam outputs
                    f.write(f"{atom_type:>6s}  {atom_type:>6s}  12.01000  0.00000  A   3.399670e-01  4.577300e-01\n")

        print(f"      ✓ Generated: {os.path.basename(atomtypes_file)}")

    def _generate_posre_file(self, rtf_data: Dict):
        """Generate position restraints file"""
        posre_file = os.path.join(self.lig_ff_dir, "posre_LIG.itp")
        print(f"      Generating position restraints file...")

        # Ensure output directory exists
        os.makedirs(self.lig_ff_dir, exist_ok=True)

        with open(posre_file, "w") as f:
            f.write("; Generated by SwissParamForceFieldGenerator\n")
            f.write("[ position_restraints ]\n")
            f.write(";  i funct       fcx        fcy        fcz\n")

            atom_id = 0
            for atom_name in rtf_data["atoms"]:
                atom_id += 1
                f.write(f"{atom_id:4d}    1       1000       1000       1000\n")

        print(f"      ✓ Generated: {os.path.basename(posre_file)}")

    def _generate_top_file(self):
        """Generate topology file"""
        top_file = os.path.join(self.lig_ff_dir, "LIG.top")
        print(f"      Generating topology file...")

        # Ensure output directory exists
        os.makedirs(self.lig_ff_dir, exist_ok=True)

        with open(top_file, "w") as f:
            f.write("; Generated by SwissParamForceFieldGenerator\n")
            f.write(f"; Approach: {self.approach}\n")
            f.write("\n")

            f.write('#include "atomtypes_LIG.itp"\n')
            f.write('#include "LIG.itp"\n')
            f.write("\n")

            f.write("[ system ]\n")
            f.write("LIG\n")
            f.write("\n")

            f.write("[ molecules ]\n")
            f.write("LIG    1\n")

        print(f"      ✓ Generated: {os.path.basename(top_file)}")

    def _copy_swissparam_itp_to_prism(self, files):
        """
        Copy existing SwissParam ITP/GRO files to PRISM standard names

        Parameters:
        -----------
        files : dict
            Dictionary of found file paths
        """
        print("\n[5b] Copying SwissParam files to PRISM format...")

        # Copy ITP file
        if files.get("itp"):
            src = files["itp"]
            dst = os.path.join(self.lig_ff_dir, "LIG.itp")
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
                print(f"      ✓ Copied ITP file")

        # Copy GRO file
        if files.get("gro"):
            src = files["gro"]
            dst = os.path.join(self.lig_ff_dir, "LIG.gro")
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
                print(f"      ✓ Copied GRO file")

        # Generate other required files
        if not os.path.exists(os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")):
            self._generate_atomtypes_file({"atoms": {"DUMMY": {"type": "X", "charge": 0.0}}})

        if not os.path.exists(os.path.join(self.lig_ff_dir, "posre_LIG.itp")):
            # Generate minimal posre file
            posre_file = os.path.join(self.lig_ff_dir, "posre_LIG.itp")
            with open(posre_file, "w") as f:
                f.write("; Generated by SwissParamForceFieldGenerator\n")
                f.write("[ position_restraints ]\n")
                f.write(";  i funct       fcx        fcy        fcz\n")
                f.write("1    1       1000       1000       1000\n")

        if not os.path.exists(os.path.join(self.lig_ff_dir, "LIG.top")):
            self._generate_top_file()

        print("      ✓ All PRISM format files ready")


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
