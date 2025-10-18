#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SwissParam force field generator wrapper for PRISM
Supports MMFF-based, MATCH, and hybrid MMFF-based-MATCH force fields
"""

import os
import sys
import re
import time
import tarfile
import shutil
import tempfile
from pathlib import Path

# Import the base class
try:
    from .base import ForceFieldGeneratorBase
except ImportError:
    from base import ForceFieldGeneratorBase

# Import required dependencies
try:
    import mechanize
except ImportError as e:
    print(f"Warning: Missing mechanize dependency: {e}")
    print("SwissParam force field generators will not be available.")
    print("INSTALLATION:")
    print("  pip install mechanize")


class SwissParamForceFieldGenerator(ForceFieldGeneratorBase):
    """Base class for SwissParam-based force field generators"""

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
        valid_approaches = ["MMFF-based", "MATCH", "MMFF-based-MATCH"]
        if approach not in valid_approaches:
            raise ValueError(f"Invalid approach '{approach}'. Must be one of {valid_approaches}")

        self.approach = approach

        # Output directory for SwissParam files
        self.lig_ff_dir = os.path.join(self.output_dir, self.get_output_dir_name())

        print(f"\nInitialized SwissParam Force Field Generator:")
        print(f"  Ligand: {self.ligand_path}")
        print(f"  Approach: {self.approach}")
        print(f"  Output directory: {self.output_dir}")

    def get_output_dir_name(self):
        """Get the output directory name for this force field"""
        if self.approach == "MMFF-based":
            return "LIG.mmff2gmx"
        elif self.approach == "MATCH":
            return "LIG.match2gmx"
        else:  # MMFF-based-MATCH
            return "LIG.hybrid2gmx"

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
                    print("\n(Use --overwrite to regenerate)")
                    return self.lig_ff_dir

            # Create output directory
            os.makedirs(self.lig_ff_dir, exist_ok=True)

            # Upload to SwissParam and download results
            sp_dir = self._upload_to_swissparam()

            # Process SwissParam output
            self._process_swissparam_output(sp_dir)

            # Cleanup temporary files
            self._cleanup_temp_files(sp_dir)

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

    def _upload_to_swissparam(self):
        """
        Upload ligand to SwissParam server and download results

        Returns:
        --------
        str : Path to extracted results directory
        """
        print("\n=== Uploading to SwissParam Server ===")

        # Create temporary working directory
        temp_dir = tempfile.mkdtemp(prefix="swissparam_")

        try:
            # Initialize browser
            br = mechanize.Browser()
            br.set_handle_robots(False)
            br.set_handle_refresh(False)
            br.addheaders = [('User-agent', 'Mozilla/5.0')]

            # Open SwissParam
            url = "http://www.swissparam.ch/"
            print(f"Opening {url}")
            response = br.open(url)

            # Select the form
            br.form = list(br.forms())[0]

            # Upload file
            print(f"Uploading {os.path.basename(self.ligand_path)} with approach={self.approach}")
            br.form.add_file(open(self.ligand_path, 'rb'), 'text/plain',
                           os.path.basename(self.ligand_path), id="fileToUpload")

            # Select approach
            br.form.find_control("approach").value = [self.approach]

            # Submit form
            print("Submitting job to SwissParam...")
            response = br.submit()

            # Parse response to get job ID
            xml = response.read().strip()
            match = re.search(r'job=([^"&]+)', xml.decode())
            if not match:
                raise ValueError("Could not extract job ID from SwissParam response")

            jobid = match.group(1)
            print(f"Job submitted successfully. Job ID: {jobid}")

            # Wait for job to process
            print("Waiting for SwissParam to generate parameters...")
            time.sleep(20)  # Give SwissParam time to generate GROMACS files

            # Download results
            result_url = f"http://www.swissparam.ch/results/{jobid}/results.tar.gz"
            print(f"Downloading results from {result_url}")

            max_retries = 10
            tar_data = None
            for attempt in range(max_retries):
                try:
                    response = br.open(result_url)
                    tar_data = response.read()

                    if len(tar_data) > 100:  # Any reasonable tarball
                        print(f"Results downloaded (attempt {attempt + 1})")
                        break
                    else:
                        print(f"Result not ready yet (attempt {attempt + 1}/{max_retries}), waiting...")
                        time.sleep(3)
                except Exception as e:
                    if attempt < max_retries - 1:
                        print(f"Download failed (attempt {attempt + 1}/{max_retries}): {e}")
                        time.sleep(3)
                    else:
                        raise

            if tar_data is None or len(tar_data) == 0:
                raise ValueError("Failed to download results after multiple attempts")

            # Save tarball
            tar_file = os.path.join(temp_dir, f"swissparam_{jobid}.tar.gz")
            with open(tar_file, 'wb') as f:
                f.write(tar_data)

            print(f"Downloaded tarball: {os.path.basename(tar_file)}")

            # Extract tarball
            extract_dir = os.path.join(temp_dir, f"swissparam_{jobid}")
            os.makedirs(extract_dir, exist_ok=True)

            with tarfile.open(tar_file, 'r:gz') as tar:
                tar.extractall(path=extract_dir)

            print(f"Extracted to: {extract_dir}")

            return extract_dir

        except Exception as e:
            # Cleanup on error
            shutil.rmtree(temp_dir, ignore_errors=True)
            raise

    def _process_swissparam_output(self, sp_dir):
        """
        Process SwissParam output and convert to PRISM standard format

        Parameters:
        -----------
        sp_dir : str
            Path to extracted SwissParam results
        """
        print("\n=== Processing SwissParam Output ===")

        # Find SwissParam files
        sp_files = self._find_swissparam_files(sp_dir)

        # Convert and standardize files
        self._standardize_files(sp_files)

    def _find_swissparam_files(self, sp_dir):
        """Find relevant files in SwissParam output"""
        files = {}

        # DEBUG: Show all files in initial extraction
        print(f"\nDEBUG: Initial contents of {sp_dir}:")
        for root, dirs, filenames in os.walk(sp_dir):
            for filename in filenames:
                filepath = os.path.join(root, filename)
                rel_path = os.path.relpath(filepath, sp_dir)
                print(f"  {rel_path}")

        # Check for tautomer files first
        tautomer_files = []
        has_tautomers_file = False
        modified_mol2 = None

        for root, dirs, filenames in os.walk(sp_dir):
            for filename in filenames:
                filepath = os.path.join(root, filename)
                if filename.endswith('_t1.mol2') or filename.endswith('_t2.mol2'):
                    tautomer_files.append(filepath)
                    print(f"DEBUG: Found explicit tautomer: {filename}")
                elif filename == 'tautomers.smi' or filename == 'tautomersImage.txt':
                    has_tautomers_file = True
                    print(f"DEBUG: Found tautomers indicator file: {filename}")
                elif filename.endswith('.mol2') and not filename.endswith('.org') and not filename.endswith('_t1.mol2') and not filename.endswith('_t2.mol2'):
                    # This might be the modified molecule
                    modified_mol2 = filepath
                    print(f"DEBUG: Found modified MOL2: {filename}")

        if tautomer_files:
            print(f"SwissParam returned {len(tautomer_files)} explicit tautomer options")
            print(f"Auto-selecting first tautomer: {os.path.basename(tautomer_files[0])}")
            return self._resubmit_with_tautomer(tautomer_files[0])
        elif has_tautomers_file:
            # SwissParam detected tautomers - this is just informational
            # The GROMACS files should still be in the tarball, continue looking
            print("SwissParam detected tautomers (informational only), looking for GROMACS files...")

        # SwissParam output structure - look for GROMACS files
        for root, dirs, filenames in os.walk(sp_dir):
            for filename in filenames:
                filepath = os.path.join(root, filename)

                # SwissParam generates files with specific names
                if filename == 'swissparam.gro' or filename.endswith('.gro'):
                    if 'gro' not in files:  # Take first match
                        files['gro'] = filepath
                elif filename == 'swissparam.itp' or filename.endswith('.itp'):
                    if 'posre' not in filename.lower() and 'itp' not in files:
                        files['itp'] = filepath
                elif filename == 'swissparam.top' or (filename.endswith('.top') and 'swissparam' in filename):
                    if 'top' not in files:
                        files['top'] = filepath

        if not files.get('gro') or not files.get('itp'):
            raise ValueError(f"Could not find required GROMACS files in {sp_dir}")

        print(f"Found SwissParam files:")
        for key, path in files.items():
            print(f"  {key}: {os.path.basename(path)}")

        return files

    def _trigger_tautomer_continuation(self, jobid, sp_dir):
        """
        Trigger continuation when SwissParam detects tautomers

        Parameters:
        -----------
        jobid : str
            SwissParam job ID
        sp_dir : str
            Current extraction directory

        Returns:
        --------
        dict : Dictionary of GROMACS file paths
        """
        print("Triggering GROMACS file generation after tautomer detection...")

        try:
            # Initialize browser
            br = mechanize.Browser()
            br.set_handle_robots(False)
            br.set_handle_refresh(False)
            br.addheaders = [('User-agent', 'Mozilla/5.0')]

            # Access the results page (HTML)
            results_url = f"http://www.swissparam.ch/results/{jobid}"
            print(f"Accessing results page: {results_url}")
            response = br.open(results_url)

            # Try to find and submit the continuation form
            # The page might have a form to continue with the detected tautomer
            try:
                br.select_form(nr=0)  # Select first form
                response = br.submit()
                print("Submitted continuation request")
                time.sleep(5)  # Wait for processing
            except:
                print("No form found, proceeding with direct download...")

            # Now download the results again
            result_url = f"http://www.swissparam.ch/results/{jobid}/results.tar.gz"
            print(f"Downloading updated results...")

            max_retries = 20
            for attempt in range(max_retries):
                try:
                    response = br.open(result_url)
                    tar_data = response.read()

                    if len(tar_data) > 1000:
                        # Check if we have GROMACS files now
                        temp_tar = os.path.join(os.path.dirname(sp_dir), f"test_cont_{jobid}.tar.gz")
                        with open(temp_tar, 'wb') as f:
                            f.write(tar_data)

                        with tarfile.open(temp_tar, 'r:gz') as tar:
                            members = tar.getnames()
                            has_gro = any('gro' in m for m in members)
                            has_itp = any('itp' in m for m in members)

                            if has_gro and has_itp:
                                print(f"GROMACS files generated successfully!")
                                # Extract to new location
                                extract_dir = os.path.join(os.path.dirname(sp_dir), f"swissparam_cont_{jobid}")
                                os.makedirs(extract_dir, exist_ok=True)
                                tar.extractall(path=extract_dir)
                                os.remove(temp_tar)

                                # Find and return GROMACS files
                                files = {}
                                for root, dirs, filenames in os.walk(extract_dir):
                                    for filename in filenames:
                                        filepath = os.path.join(root, filename)
                                        if filename.endswith('.gro'):
                                            files['gro'] = filepath
                                        elif filename.endswith('.itp') and 'posre' not in filename.lower():
                                            files['itp'] = filepath
                                        elif filename.endswith('.top'):
                                            files['top'] = filepath

                                if files.get('gro') and files.get('itp'):
                                    return files
                            else:
                                print(f"Still no GROMACS files (attempt {attempt + 1}/{max_retries}), waiting...")
                                os.remove(temp_tar)
                                time.sleep(10)
                    else:
                        print(f"Waiting for results (attempt {attempt + 1}/{max_retries})...")
                        time.sleep(10)

                except Exception as e:
                    if attempt < max_retries - 1:
                        print(f"Download failed (attempt {attempt + 1}): {e}")
                        time.sleep(10)
                    else:
                        raise

            raise ValueError("Could not generate GROMACS files after tautomer continuation")

        except Exception as e:
            print(f"Error during tautomer continuation: {e}")
            raise

    def _resubmit_with_tautomer(self, tautomer_file):
        """
        Re-submit job with selected tautomer

        Parameters:
        -----------
        tautomer_file : str
            Path to tautomer MOL2 file

        Returns:
        --------
        dict : Dictionary of GROMACS file paths
        """
        print(f"Re-submitting with selected tautomer...")

        # Create new temp directory
        temp_dir = tempfile.mkdtemp(prefix="swissparam_tautomer_")

        try:
            # Initialize browser
            br = mechanize.Browser()
            br.set_handle_robots(False)
            br.set_handle_refresh(False)
            br.addheaders = [('User-agent', 'Mozilla/5.0')]

            # Open SwissParam
            url = "http://www.swissparam.ch/"
            br.open(url)
            br.form = list(br.forms())[0]

            # Upload tautomer file
            br.form.add_file(open(tautomer_file, 'rb'), 'text/plain',
                           os.path.basename(tautomer_file), id="fileToUpload")

            # Select approach
            br.form.find_control("approach").value = [self.approach]

            # Submit
            response = br.submit()

            # Parse job ID
            xml = response.read().strip()
            print(f"DEBUG: Tautomer submission response length: {len(xml)}")
            print(f"DEBUG: First 500 chars of response: {xml[:500].decode() if len(xml) > 0 else 'EMPTY'}")

            # Check if response contains an error message
            response_text = xml.decode()
            if 'ERROR:' in response_text or 'error' in response_text.lower():
                # Extract error message
                error_match = re.search(r'(MATCH ERROR:[^<]+)', response_text)
                if error_match:
                    error_msg = error_match.group(1).strip()
                    print(f"\n{'='*60}")
                    print(f"SwissParam Error: {error_msg}")
                    print(f"{'='*60}")
                    print("\nThis tautomer cannot be processed with MATCH force field.")
                    print("Possible solutions:")
                    print("  1. Try a different force field: --ligand-forcefield mmff")
                    print("  2. Try the hybrid approach: --ligand-forcefield hybrid")
                    print("  3. Use GAFF/GAFF2: --ligand-forcefield gaff2")
                    print("  4. Manually prepare the ligand structure to avoid tautomers")
                    raise ValueError(f"SwissParam MATCH force field failed: {error_msg}")
                else:
                    print(f"ERROR: SwissParam returned an error response")
                    print(f"Response preview:\n{response_text[:1000]}")
                    raise ValueError("SwissParam tautomer submission failed")

            # Try multiple patterns to extract job ID
            match = re.search(r'job=([^"&]+)', response_text)
            if not match:
                # Try alternative patterns
                match = re.search(r'/results/(\d+)', response_text)
            if not match:
                match = re.search(r'jobid["\s:=]+(\d+)', response_text, re.IGNORECASE)

            if not match:
                print(f"ERROR: Could not extract job ID from response")
                print(f"Full response:\n{response_text}")
                raise ValueError("Could not extract job ID from tautomer submission")

            jobid = match.group(1)
            print(f"Tautomer job submitted. Job ID: {jobid}")

            time.sleep(2)

            # Download results
            result_url = f"http://www.swissparam.ch/results/{jobid}/results.tar.gz"

            max_retries = 10
            for attempt in range(max_retries):
                try:
                    response = br.open(result_url)
                    tar_data = response.read()
                    if len(tar_data) > 0:
                        break
                    time.sleep(3)
                except Exception as e:
                    if attempt < max_retries - 1:
                        time.sleep(3)
                    else:
                        raise

            # Save and extract
            tar_file = os.path.join(temp_dir, f"swissparam_tautomer_{jobid}.tar.gz")
            with open(tar_file, 'wb') as f:
                f.write(tar_data)

            extract_dir = os.path.join(temp_dir, f"swissparam_tautomer_{jobid}")
            os.makedirs(extract_dir, exist_ok=True)

            with tarfile.open(tar_file, 'r:gz') as tar:
                tar.extractall(path=extract_dir)

            print(f"Extracted tautomer results")

            # DEBUG: Show all files in extracted directory
            print(f"\nDEBUG: Contents of {extract_dir}:")
            for root, dirs, filenames in os.walk(extract_dir):
                for filename in filenames:
                    filepath = os.path.join(root, filename)
                    rel_path = os.path.relpath(filepath, extract_dir)
                    print(f"  {rel_path}")

            # Find files in extracted directory
            files = {}
            for root, dirs, filenames in os.walk(extract_dir):
                for filename in filenames:
                    filepath = os.path.join(root, filename)

                    if filename == 'swissparam.gro' or filename.endswith('.gro'):
                        if 'gro' not in files:
                            files['gro'] = filepath
                            print(f"DEBUG: Found GRO file: {filename}")
                    elif filename == 'swissparam.itp' or filename.endswith('.itp'):
                        if 'posre' not in filename.lower() and 'itp' not in files:
                            files['itp'] = filepath
                            print(f"DEBUG: Found ITP file: {filename}")
                    elif filename == 'swissparam.top' or (filename.endswith('.top') and 'swissparam' in filename):
                        if 'top' not in files:
                            files['top'] = filepath
                            print(f"DEBUG: Found TOP file: {filename}")

            if not files.get('gro') or not files.get('itp'):
                print(f"\nDEBUG: Files found: {list(files.keys())}")
                raise ValueError(f"Could not find GROMACS files after tautomer submission")

            print(f"Found GROMACS files from tautomer:")
            for key, path in files.items():
                print(f"  {key}: {os.path.basename(path)}")

            return files

        except Exception as e:
            shutil.rmtree(temp_dir, ignore_errors=True)
            raise

    def _standardize_files(self, sp_files):
        """
        Standardize SwissParam files to PRISM format (LIG naming)

        Parameters:
        -----------
        sp_files : dict
            Dictionary of file paths
        """
        print("\n=== Standardizing to LIG naming ===")

        # Process GRO file
        if 'gro' in sp_files:
            self._standardize_gro(sp_files['gro'], os.path.join(self.lig_ff_dir, "LIG.gro"))

        # Process ITP file
        if 'itp' in sp_files:
            self._standardize_itp(sp_files['itp'], os.path.join(self.lig_ff_dir, "LIG.itp"))

        # Process TOP file or create from ITP
        if 'top' in sp_files:
            self._standardize_top(sp_files['top'], os.path.join(self.lig_ff_dir, "LIG.top"))
        else:
            self._create_top_from_itp(os.path.join(self.lig_ff_dir, "LIG.itp"))

        # Extract atomtypes to separate file
        self._extract_atomtypes(os.path.join(self.lig_ff_dir, "LIG.itp"))

        # Create position restraints file
        self._create_position_restraints(os.path.join(self.lig_ff_dir, "LIG.gro"))

        print("Standardization complete")

    def _standardize_gro(self, source_gro, target_gro):
        """Standardize GRO file to use LIG residue name"""
        with open(source_gro, 'r') as f:
            lines = f.readlines()

        with open(target_gro, 'w') as f:
            f.write("LIG system\n")  # Title

            for i, line in enumerate(lines):
                if i == 0:
                    continue  # Skip original title
                elif i >= 2 and i < len(lines) - 1:  # Atom lines
                    if len(line) > 10:
                        # Replace residue name (columns 5-10)
                        line = line[:5] + "LIG".ljust(5) + line[10:]

                f.write(line)

    def _standardize_itp(self, source_itp, target_itp):
        """Standardize ITP file to use LIG molecule name"""
        with open(source_itp, 'r') as f:
            content = f.read()

        # Replace common SwissParam molecule names with LIG (global)
        content = re.sub(r'\bUNL\b', 'LIG', content)  # SwissParam default
        content = re.sub(r'\bunl\b', 'LIG', content)
        content = re.sub(r'\bligand\b', 'LIG', content, flags=re.IGNORECASE)

        # Also replace the ligand name from file
        ligand_basename = Path(self.ligand_path).stem
        content = re.sub(rf'\b{re.escape(ligand_basename)}\b', 'LIG', content, flags=re.IGNORECASE)

        # Also replace in [ moleculetype ] section specifically
        content = re.sub(r'(\[\s*moleculetype\s*\]\s*;\s*\w+\s+\w+\s+)\S+', r'\1LIG', content)

        # Replace residue names in [ atoms ] section (4th column)
        lines = content.split('\n')
        new_lines = []
        in_atoms = False

        for line in lines:
            if line.strip().startswith('[ atoms ]'):
                in_atoms = True
                new_lines.append(line)
            elif in_atoms and line.strip().startswith('['):
                in_atoms = False
                new_lines.append(line)
            elif in_atoms and not line.strip().startswith(';') and line.strip():
                # Use regex to replace residue name (4th field) while preserving spacing
                modified_line = re.sub(
                    r'^(\s*\d+\s+\S+\s+\d+\s+)\S+(\s+)',
                    r'\1LIG\2',
                    line
                )
                new_lines.append(modified_line)
            else:
                new_lines.append(line)

        content = '\n'.join(new_lines)

        # Add position restraints include if not present
        if "#ifdef POSRES" not in content:
            content += "\n#ifdef POSRES\n#include \"posre_LIG.itp\"\n#endif\n"

        with open(target_itp, 'w') as f:
            f.write(content)

    def _standardize_top(self, source_top, target_top):
        """Standardize TOP file"""
        with open(source_top, 'r') as f:
            content = f.read()

        # Replace molecule names
        content = re.sub(r'\bUNL\b', 'LIG', content)
        content = re.sub(r'\bunl\b', 'LIG', content)
        content = re.sub(r'\bligand\b', 'LIG', content, flags=re.IGNORECASE)

        ligand_basename = Path(self.ligand_path).stem
        content = re.sub(rf'\b{re.escape(ligand_basename)}\b', 'LIG', content, flags=re.IGNORECASE)

        # Update includes
        content = re.sub(r'#include\s+"swissparam\.itp"', '#include "LIG.itp"', content)
        content = re.sub(r'#include\s+"[^"]*\.itp"', '#include "LIG.itp"', content)

        with open(target_top, 'w') as f:
            f.write(content)

    def _create_top_from_itp(self, itp_file):
        """Create TOP file from ITP"""
        top_file = os.path.join(self.lig_ff_dir, "LIG.top")

        with open(top_file, 'w') as f:
            f.write(f"; SwissParam ({self.approach}) topology for LIG\n\n")
            f.write("[ defaults ]\n")
            f.write("; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n")
            f.write("1         2          yes        1.0      1.0\n\n")
            f.write("#include \"atomtypes_LIG.itp\"\n")
            f.write("#include \"LIG.itp\"\n\n")
            f.write("[ system ]\n")
            f.write("LIG system\n\n")
            f.write("[ molecules ]\n")
            f.write("LIG              1\n")

    def _extract_atomtypes(self, itp_file):
        """Extract atomtypes to separate file"""
        atomtypes_file = os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")

        with open(itp_file, 'r') as f:
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

        # Write atomtypes file if atomtypes section exists
        if atomtypes_lines:
            with open(atomtypes_file, 'w') as f:
                f.writelines(atomtypes_lines)

            # Update ITP file (remove atomtypes section)
            with open(itp_file, 'w') as f:
                f.writelines(new_itp_lines)
        else:
            # If no atomtypes in ITP, create empty file to satisfy requirements
            with open(atomtypes_file, 'w') as f:
                f.write("; Atomtypes defined in force field\n")

    def _create_position_restraints(self, gro_file):
        """Create position restraints file"""
        posre_file = os.path.join(self.lig_ff_dir, "posre_LIG.itp")

        # Count atoms from GRO file
        with open(gro_file, 'r') as f:
            lines = f.readlines()

        try:
            atom_count = int(lines[1].strip())
        except (IndexError, ValueError):
            atom_count = len(lines) - 3  # Estimate

        # Write position restraints
        with open(posre_file, 'w') as f:
            f.write("[ position_restraints ]\n")
            f.write("; atom  type      fx      fy      fz\n")
            for i in range(1, atom_count + 1):
                f.write(f"{i:>6}     1  1000  1000  1000\n")

    def _cleanup_temp_files(self, temp_dir):
        """Clean up temporary files"""
        print("\n=== Cleaning up ===")
        if os.path.exists(temp_dir):
            parent_dir = os.path.dirname(temp_dir)
            shutil.rmtree(temp_dir, ignore_errors=True)
            print("Temporary files removed")


# Convenience classes for specific force field approaches
class MMFFForceFieldGenerator(SwissParamForceFieldGenerator):
    """MMFF-based force field generator using SwissParam"""

    def __init__(self, ligand_path, output_dir, overwrite=False):
        """
        Initialize MMFF-based force field generator

        Parameters:
        -----------
        ligand_path : str
            Path to the ligand file (MOL2/SDF/PDB)
        output_dir : str
            Directory where output files will be stored
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(ligand_path, output_dir, approach="MMFF-based", overwrite=overwrite)


class MATCHForceFieldGenerator(SwissParamForceFieldGenerator):
    """MATCH force field generator using SwissParam"""

    def __init__(self, ligand_path, output_dir, overwrite=False):
        """
        Initialize MATCH force field generator

        Parameters:
        -----------
        ligand_path : str
            Path to the ligand file (MOL2/SDF/PDB)
        output_dir : str
            Directory where output files will be stored
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(ligand_path, output_dir, approach="MATCH", overwrite=overwrite)


class HybridMMFFMATCHForceFieldGenerator(SwissParamForceFieldGenerator):
    """Hybrid MMFF-based-MATCH force field generator using SwissParam"""

    def __init__(self, ligand_path, output_dir, overwrite=False):
        """
        Initialize hybrid MMFF-based-MATCH force field generator

        Parameters:
        -----------
        ligand_path : str
            Path to the ligand file (MOL2/SDF/PDB)
        output_dir : str
            Directory where output files will be stored
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(ligand_path, output_dir, approach="MMFF-based-MATCH", overwrite=overwrite)


if __name__ == '__main__':
    # Example usage
    import sys

    if len(sys.argv) < 2:
        print("Usage: python swissparam.py <ligand_file.mol2> [approach]")
        print("  approach: MMFF-based (default), MATCH, or MMFF-based-MATCH")
        sys.exit(1)

    ligand_file = sys.argv[1]
    approach = sys.argv[2] if len(sys.argv) > 2 else "MMFF-based"
    output_dir = os.path.dirname(ligand_file) if os.path.dirname(ligand_file) else "."

    generator = SwissParamForceFieldGenerator(
        ligand_path=ligand_file,
        output_dir=output_dir,
        approach=approach
    )

    result_dir = generator.run()
    print(f"\nSuccess! Output in: {result_dir}")
