#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gaussian RESP Workflow Manager

This module manages the complete workflow for Gaussian RESP charge calculations,
including file preparation, script generation, and optional automatic execution.
"""

import os
from typing import Dict, Optional

from .converter import CoordinateConverter
from .gjf_generator import GaussianInputGenerator
from .utils import (
    check_gaussian_available,
    check_antechamber_available,
    get_charge_multiplicity,
)

# Import color utilities
try:
    from ..utils.colors import print_success, print_info, print_warning, print_error, print_header
except ImportError:
    try:
        from prism.utils.colors import print_success, print_info, print_warning, print_error, print_header
    except ImportError:
        def print_success(x, **kwargs): print(f"[OK] {x}")
        def print_info(x, **kwargs): print(f"[INFO] {x}")
        def print_warning(x, **kwargs): print(f"[WARN] {x}")
        def print_error(x, **kwargs): print(f"[ERROR] {x}")
        def print_header(x): print(f"\n{'='*60}\n{x}\n{'='*60}")


class GaussianRESPWorkflow:
    """
    Workflow manager for Gaussian RESP charge calculations.

    Handles the complete workflow from MOL2 input to RESP charges:
    1. Convert MOL2 to XYZ coordinates
    2. Generate Gaussian input files (opt and/or ESP)
    3. Generate run script for Gaussian calculations
    4. (Optional) Run Gaussian calculations
    5. Run RESP fitting via antechamber
    """

    def __init__(self, ligand_path: str, output_dir: str,
                 method: str = 'dft', do_optimization: bool = False,
                 nproc: int = 16, mem: str = '4GB', verbose: bool = True):
        """
        Initialize the Gaussian RESP workflow.

        Parameters
        ----------
        ligand_path : str
            Path to input ligand MOL2 file
        output_dir : str
            Output directory for Gaussian files
        method : str
            Quantum chemistry method: 'hf' or 'dft' (default: 'dft')
        do_optimization : bool
            Whether to perform geometry optimization (default: False)
        nproc : int
            Number of processors for Gaussian (default: 16)
        mem : str
            Memory allocation for Gaussian (default: '4GB')
        verbose : bool
            Whether to print progress information
        """
        self.ligand_path = os.path.abspath(ligand_path)
        self.output_dir = os.path.abspath(output_dir)
        self.method = method.lower()
        self.do_optimization = do_optimization
        self.nproc = nproc
        self.mem = mem
        self.verbose = verbose

        # Derived paths
        self.ligand_name = os.path.splitext(os.path.basename(ligand_path))[0]

        # Initialize generators
        self.converter = CoordinateConverter(verbose=verbose)
        self.gjf_generator = GaussianInputGenerator(method=method, nproc=nproc, mem=mem)

        # File paths (set during prepare_files)
        self.files = {}

    def prepare_files(self) -> Dict[str, str]:
        """
        Generate all necessary files for Gaussian RESP calculation.

        Returns
        -------
        Dict[str, str]
            Dictionary with paths to generated files:
            - xyz: XYZ coordinate file
            - opt_gjf: Optimization input file (if do_optimization=True)
            - esp_gjf: ESP input file
            - script: Run script path
        """
        if self.verbose:
            print_header("Preparing Gaussian RESP Files")
            print_info(f"  Ligand: {self.ligand_path}")
            print_info(f"  Method: {self.method.upper()}")
            print_info(f"  Optimization: {'Yes' if self.do_optimization else 'No'}")

        os.makedirs(self.output_dir, exist_ok=True)

        # Get charge and multiplicity (with verbose debug output)
        charge, multiplicity = get_charge_multiplicity(self.ligand_path, verbose=self.verbose)
        if self.verbose:
            print_info(f"  Charge: {charge}, Multiplicity: {multiplicity}")

        # 1. Convert MOL2 to XYZ
        xyz_file = os.path.join(self.output_dir, "ligand.xyz")
        self.converter.mol2_to_xyz(self.ligand_path, xyz_file)
        self.files['xyz'] = xyz_file

        # 2. Generate Gaussian input files
        if self.do_optimization:
            # Generate optimization input
            opt_gjf = os.path.join(self.output_dir, f"lig_opt_{self.method}.gjf")
            self.gjf_generator.generate_opt_gjf(xyz_file, self.ligand_path, opt_gjf)
            self.files['opt_gjf'] = opt_gjf

            # ESP file path (will be generated from optimized coords)
            esp_gjf = os.path.join(self.output_dir, f"lig_esp_{self.method}.gjf")
            self.files['esp_gjf'] = esp_gjf
        else:
            # Generate ESP input directly
            esp_gjf = os.path.join(self.output_dir, f"lig_esp_{self.method}.gjf")
            self.gjf_generator.generate_esp_direct(xyz_file, self.ligand_path, esp_gjf)
            self.files['esp_gjf'] = esp_gjf

        # 3. Generate run script
        script_path = self.generate_charge_script()
        self.files['script'] = script_path

        if self.verbose:
            print_success("Files generated successfully")
            for key, path in self.files.items():
                print_info(f"  {key}: {path}")

        return self.files

    def generate_charge_script(self) -> str:
        """
        Generate shell script for running Gaussian RESP workflow.

        Returns
        -------
        str
            Path to generated script
        """
        script_path = os.path.join(self.output_dir, "charge.sh")

        # Get method-specific filenames
        method = self.method
        isopt = "true" if self.do_optimization else "false"

        script_content = f'''#!/bin/bash
# PRISM Gaussian RESP Charge Calculation Script
# Ligand: {self.ligand_name}
# Generated by PRISM Gaussian module

set -e
cd "$(dirname "$0")"

METHOD="{method}"  # hf or dft
ISOPT="{isopt}"    # true or false

echo "========================================"
echo "PRISM Gaussian RESP Charge Calculation"
echo "Method: $METHOD | Optimization: $ISOPT"
echo "========================================"

# Check Gaussian
if ! command -v g16 &> /dev/null; then
    echo "ERROR: g16 not found in PATH"
    echo "Please load Gaussian environment first"
    exit 1
fi

# Check antechamber
if ! command -v antechamber &> /dev/null; then
    echo "ERROR: antechamber not found in PATH"
    echo "Please install AmberTools: conda install -c conda-forge ambertools"
    exit 1
fi

'''

        if self.do_optimization:
            script_content += f'''
# Step 1: Geometry optimization
echo ""
echo "Step 1: Geometry optimization..."
if [ -f "lig_opt_{method}.log" ] && grep -q "Normal termination" "lig_opt_{method}.log"; then
    echo "  Optimization already completed, skipping..."
else
    g16 lig_opt_{method}.gjf
    if ! grep -q "Normal termination" "lig_opt_{method}.log"; then
        echo "ERROR: Geometry optimization failed!"
        echo "Check lig_opt_{method}.log for details"
        exit 1
    fi
    echo "  Optimization completed successfully"
fi

# Step 2: Extract optimized coordinates and generate ESP input
echo ""
echo "Step 2: Extracting optimized coordinates..."
python3 - << 'PYTHON_SCRIPT'
import re
import sys

def extract_coords_from_log(log_file):
    """Extract optimized coordinates from Gaussian log file."""
    with open(log_file, 'r') as f:
        content = f.read()

    # Element mapping
    ELEMENTS = {{
        1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I'
    }}

    # Find the last Standard orientation
    pattern = r"Standard orientation:.*?-+\\n.*?-+\\n(.*?)-+\\n"
    matches = re.findall(pattern, content, re.DOTALL)

    if not matches:
        pattern = r"Input orientation:.*?-+\\n.*?-+\\n(.*?)-+\\n"
        matches = re.findall(pattern, content, re.DOTALL)

    if not matches:
        print("ERROR: Could not find coordinates in log file")
        sys.exit(1)

    coords = []
    for line in matches[-1].strip().split('\\n'):
        parts = line.split()
        if len(parts) >= 6:
            try:
                atomic_num = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                element = ELEMENTS.get(atomic_num, 'X')
                coords.append((element, x, y, z))
            except (ValueError, IndexError):
                continue

    return coords

def get_charge_mult_from_log(log_file):
    """Get charge and multiplicity from log file."""
    with open(log_file, 'r') as f:
        for line in f:
            if "Charge =" in line and "Multiplicity =" in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "Charge":
                        charge = int(parts[i + 2])
                    if part == "Multiplicity":
                        mult = int(parts[i + 2])
                return charge, mult
    return 0, 1

# Extract coordinates
coords = extract_coords_from_log("lig_opt_{method}.log")
charge, mult = get_charge_mult_from_log("lig_opt_{method}.log")

# Generate ESP input file
with open("lig_esp_{method}.gjf", 'w') as f:
    f.write("%nproc={self.nproc}\\n")
    f.write("%mem={self.mem}\\n")
    f.write("%chk=lig_esp_{method}.chk\\n")
    f.write("#p {self.gjf_generator.method}/6-31G* Pop=MK IOp(6/33=2,6/42=6)\\n")
    f.write("\\n")
    f.write("ESP calculation for RESP fitting\\n")
    f.write("\\n")
    f.write(f"{{charge}} {{mult}}\\n")
    for elem, x, y, z in coords:
        f.write(f" {{elem:2s}}    {{x:14.8f}}  {{y:14.8f}}  {{z:14.8f}}\\n")
    f.write("\\n\\n")

print(f"Generated ESP input with {{len(coords)}} atoms")
PYTHON_SCRIPT

echo "  ESP input file generated"

'''
        else:
            script_content += '''
echo ""
echo "Step 1: Skipping optimization (direct ESP calculation)"
'''

        # Common ESP and RESP steps
        script_content += f'''
# Step 3: ESP calculation
echo ""
echo "Step 3: ESP calculation..."
if [ -f "lig_esp_{method}.log" ] && grep -q "Normal termination" "lig_esp_{method}.log"; then
    echo "  ESP calculation already completed, skipping..."
else
    g16 lig_esp_{method}.gjf
    if ! grep -q "Normal termination" "lig_esp_{method}.log"; then
        echo "ERROR: ESP calculation failed!"
        echo "Check lig_esp_{method}.log for details"
        exit 1
    fi
    echo "  ESP calculation completed successfully"
fi

# Step 4: RESP fitting
echo ""
echo "Step 4: RESP fitting with antechamber..."
antechamber -i lig_esp_{method}.log -fi gout -o ligand_resp.mol2 -fo mol2 -c resp

if [ ! -f "ligand_resp.mol2" ]; then
    echo "ERROR: RESP fitting failed!"
    exit 1
fi

echo ""
echo "========================================"
echo "RESP charges calculated successfully!"
echo "========================================"
echo ""
echo "Output: $(pwd)/ligand_resp.mol2"
echo ""
echo "To use these charges in PRISM, re-run with:"
echo "  prism ... --respfile $(pwd)/ligand_resp.mol2"
echo ""
'''

        with open(script_path, 'w') as f:
            f.write(script_content)

        # Make executable
        os.chmod(script_path, 0o755)

        if self.verbose:
            print_success(f"Generated run script: {script_path}")

        return script_path

    def run_gaussian(self) -> Optional[str]:
        """
        Execute Gaussian calculations if g16 is available.

        Returns
        -------
        Optional[str]
            Path to ligand_resp.mol2 if successful, None otherwise
        """
        if not check_gaussian_available():
            if self.verbose:
                print_warning("Gaussian (g16) not available")
            return None

        if not check_antechamber_available():
            if self.verbose:
                print_warning("antechamber not available")
            return None

        if self.verbose:
            print_header("Running Gaussian RESP Calculation")

        import subprocess

        # Run the generated script
        script_path = self.files.get('script')
        if not script_path or not os.path.exists(script_path):
            print_error("Run script not found. Call prepare_files() first.")
            return None

        try:
            result = subprocess.run(
                ['bash', script_path],
                cwd=self.output_dir,
                capture_output=True,
                text=True
            )

            if result.returncode == 0:
                resp_mol2 = os.path.join(self.output_dir, "ligand_resp.mol2")
                if os.path.exists(resp_mol2):
                    if self.verbose:
                        print_success(f"RESP calculation completed: {resp_mol2}")
                    return resp_mol2

            if self.verbose:
                print_error("Gaussian calculation failed")
                if result.stdout:
                    print(result.stdout)
                if result.stderr:
                    print(result.stderr)

        except Exception as e:
            if self.verbose:
                print_error(f"Error running Gaussian: {e}")

        return None

    def run(self) -> Optional[str]:
        """
        Run the complete workflow.

        Returns
        -------
        Optional[str]
            Path to ligand_resp.mol2 if successful, None otherwise
        """
        # Prepare files
        self.prepare_files()

        # Try to run Gaussian
        if check_gaussian_available() and check_antechamber_available():
            return self.run_gaussian()
        else:
            if self.verbose:
                print_warning("Gaussian or antechamber not available")
                print_info(f"Files prepared in: {self.output_dir}")
                print_info(f"Run manually: bash {self.files.get('script', 'charge.sh')}")
            return None
