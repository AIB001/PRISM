#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GAFF force field generator wrapper for PRISM
"""

import os
import sys
import subprocess
import shutil
import re
from pathlib import Path

# Import the base class
try:
    from .base import ForceFieldGeneratorBase
except ImportError:
    from base import ForceFieldGeneratorBase

# Import color utilities
try:
    from ..utils.colors import print_success, print_info, print_warning, success, path, number
except ImportError:
    try:
        from prism.utils.colors import print_success, print_info, print_warning, success, path, number
    except ImportError:
        # Fallback if colors not available
        def print_success(x, **kwargs): print(f"✓ {x}")
        def print_info(x, **kwargs): print(f"ℹ {x}")
        def print_warning(x, **kwargs): print(f"⚠ {x}")
        def success(x): return f"✓ {x}"
        def path(x): return x
        def number(x): return x


class GAFFForceFieldGenerator(ForceFieldGeneratorBase):
    """GAFF force field generator wrapper"""
    
    def __init__(self, ligand_path, output_dir, overwrite=False):
        """Initialize GAFF force field generator"""
        super().__init__(ligand_path, output_dir, overwrite)
        
        # Check if RDKit is available
        self.rdkit_available = self._check_rdkit()
        
        # Detect file format
        self.file_format = self._detect_file_format()
        
        # Extract ligand information
        self.ligand_id = self._extract_ligand_id()
        
        # Auto-detect molecule net charge
        self.net_charge = self._auto_detect_net_charge()
        
        # Initialize list of possible molecule names
        self.possible_molecule_names = [
            self.ligand_id,
            self.ligand_name,
            self.ligand_id.upper(),
            self.ligand_name.upper(),
            "LIG"
        ]
        
        print(f"\n{'='*60}")
        print_info("Initialized GAFF Force Field Generator")
        print(f"{'='*60}")
        print(f"  Ligand: {path(self.ligand_path)}")
        print(f"  File format: {self.file_format.upper()}")
        print(f"  Ligand ID: {self.ligand_id}")
        print(f"  Net charge: {number(self.net_charge)}")
        print(f"  Output directory: {path(self.output_dir)}")
        print(f"  RDKit available: {success('Yes') if self.rdkit_available else '❌ No'}")
    
    def get_output_dir_name(self):
        """Get the output directory name for GAFF"""
        return "LIG.amb2gmx"
    
    def run(self):
        """Run the GAFF force field generation workflow"""
        print(f"\n{'='*60}")
        print_info("Starting GAFF Force Field Generation")
        print(f"{'='*60}")

        try:
            # Check if output already exists BEFORE trying to generate parameters
            # This allows using cached files even when antechamber/acpype are not available
            lig_dir = os.path.join(self.output_dir, "LIG.amb2gmx")
            if os.path.exists(lig_dir) and not self.overwrite:
                if self.check_required_files(lig_dir):
                    print_info(f"Using cached GAFF parameters from: {path(lig_dir)}")
                    print("All required files found:")
                    for f in sorted(os.listdir(lig_dir)):
                        print_success(f, prefix="  -")
                    print_warning("Use --overwrite to regenerate")
                    return lig_dir

            # Generate AMBER parameters
            print("  Generating AMBER parameters...")
            ff_dir = self.generate_amber_parameters()
            print_success("  AMBER parameters generated")

            # Find acpype output
            print("  Converting to GROMACS format...")
            ff_files = self.find_acpype_output(ff_dir)

            # Standardize to LIG naming
            lig_dir = self.standardize_to_LIG(ff_files)
            print_success("  Converted to GROMACS format")

            # Cleanup
            self.cleanup_temp_files()

            print(f"\n{'='*60}")
            print_success("GAFF force field generation completed!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {path(lig_dir)}")
            print("\nGenerated files:")
            for f in sorted(os.listdir(lig_dir)):
                print(f"  - {f}")

            return lig_dir

        except Exception as e:
            print(f"\nError during force field generation: {e}")
            raise
    
    # The following methods are copied from the original GAFF_ForceFieldGenerator.py
    # with minimal modifications to work within the wrapper structure
    
    def _check_rdkit(self):
        """Check if RDKit is available"""
        try:
            from rdkit import Chem
            return True
        except ImportError:
            return False
    
    def _detect_file_format(self):
        """Detect input file format (MOL2 or SDF)"""
        file_ext = os.path.splitext(self.ligand_path)[1].lower()
        
        if file_ext == '.mol2':
            return 'mol2'
        elif file_ext in ['.sdf', '.sd']:
            return 'sdf'
        else:
            print(f"Warning: Unknown file extension '{file_ext}', assuming MOL2 format")
            return 'mol2'
    
    def _extract_ligand_id(self):
        """Extract ligand ID from MOL2 or SDF file"""
        try:
            if self.file_format == 'mol2':
                with open(self.ligand_path, 'r') as f:
                    for line in f:
                        if line.startswith("@<TRIPOS>MOLECULE"):
                            ligand_id = next(f).strip()
                            return ligand_id
            elif self.file_format == 'sdf':
                if self.rdkit_available:
                    from rdkit import Chem
                    suppl = Chem.SDMolSupplier(self.ligand_path)
                    for mol in suppl:
                        if mol is not None:
                            if mol.HasProp('_Name'):
                                return mol.GetProp('_Name')
                            break
                else:
                    with open(self.ligand_path, 'r') as f:
                        first_line = f.readline().strip()
                        if first_line:
                            return first_line
        except Exception as e:
            print(f"Failed to extract ligand ID from {self.file_format.upper()} file: {e}")
        
        return self.ligand_name
    
    def _auto_detect_net_charge(self):
        """Automatically detect molecule net charge"""
        net_charge = 0
        
        try:
            if self.file_format == 'mol2':
                net_charge = self._detect_charge_from_mol2()
            elif self.file_format == 'sdf' and self.rdkit_available:
                net_charge = self._detect_charge_from_sdf()
            else:
                net_charge = self._detect_charge_with_antechamber()
        except Exception as e:
            print(f"Warning: Could not auto-detect charge ({e}), assuming neutral molecule")
            net_charge = 0
        
        if net_charge != 0:
            print(f"Auto-detected net charge: {net_charge}")
        
        return net_charge
    
    def _detect_charge_from_mol2(self):
        """Detect net charge from MOL2 file"""
        total_charge = 0.0
        in_atom_section = False
        
        with open(self.ligand_path, 'r') as f:
            for line in f:
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atom_section = True
                    continue
                elif line.startswith("@<TRIPOS>"):
                    in_atom_section = False
                elif in_atom_section and line.strip():
                    parts = line.split()
                    if len(parts) >= 9:
                        try:
                            charge = float(parts[8])
                            total_charge += charge
                        except ValueError:
                            pass
        
        net_charge = int(round(total_charge))
        return net_charge
    
    def _detect_charge_from_sdf(self):
        """Detect net charge from SDF file using RDKit"""
        from rdkit import Chem
        
        suppl = Chem.SDMolSupplier(self.ligand_path, removeHs=False)
        for mol in suppl:
            if mol is not None:
                net_charge = Chem.GetFormalCharge(mol)
                return net_charge
        
        return 0
    
    def _detect_charge_with_antechamber(self):
        """Use antechamber to detect charge"""
        temp_dir = os.path.join(self.output_dir, "temp_charge_detect")
        os.makedirs(temp_dir, exist_ok=True)
        
        try:
            temp_output = os.path.join(temp_dir, "temp.mol2")
            result = self.run_command([
                "antechamber",
                "-i", self.ligand_path,
                "-fi", self.file_format,
                "-o", temp_output,
                "-fo", "mol2",
                "-c", "gas",
                "-s", "0"
            ], check=False)
            
            if "Total number of electrons" in result:
                for line in result.split('\n'):
                    if "Net charge" in line:
                        match = re.search(r'Net charge:\s*([-\d]+)', line)
                        if match:
                            return int(match.group(1))
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
        
        return 0
    
    def run_command(self, cmd, cwd=None, check=True):
        """Run a shell command and return its output"""
        print(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")

        try:
            result = subprocess.run(
                cmd,
                cwd=cwd,
                check=check,
                text=True,
                capture_output=True
            )

            # Print full stdout if available
            if result.stdout.strip():
                print(f"STDOUT:\n{result.stdout}")

            # Print stderr even on success (antechamber often writes important info to stderr)
            if result.stderr.strip():
                print(f"STDERR:\n{result.stderr}")

            return result.stdout

        except subprocess.CalledProcessError as e:
            print(f"\n{'='*60}")
            print(f"ERROR executing command:")
            print(f"  Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
            print(f"  Return code: {e.returncode}")
            print(f"  Working directory: {cwd if cwd else 'current directory'}")
            print(f"{'='*60}")

            if e.stdout and e.stdout.strip():
                print(f"\nSTDOUT:\n{e.stdout}")

            if e.stderr and e.stderr.strip():
                print(f"\nSTDERR:\n{e.stderr}")
            else:
                print("\nNo error output captured from stderr")

            print(f"{'='*60}\n")

            if check:
                raise
            return e.stdout
    
    def check_command_exists(self, command):
        """Check if a command exists in the system"""
        try:
            subprocess.run(['which', command], check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError:
            return False
    
    def generate_amber_parameters(self):
        """Generate AMBER force field parameters using AmberTools"""
        print("\n=== Generating AMBER Parameters ===")

        # Check for required tools
        missing_tools = []
        for tool in ['antechamber', 'parmchk2', 'tleap', 'acpype']:
            if not self.check_command_exists(tool):
                missing_tools.append(tool)

        if missing_tools:
            raise RuntimeError(
                f"Missing required AmberTools components: {', '.join(missing_tools)}\n"
                f"Please install AmberTools via: conda install -c conda-forge ambertools\n"
                f"Or visit: http://ambermd.org/GetAmber.php"
            )

        # Check if sqm is available (required for BCC charges)
        sqm_available = self.check_command_exists('sqm')
        if not sqm_available:
            print("Warning: 'sqm' not found. AM1-BCC charge calculation may fail.")
            print("Will automatically fall back to gas phase charges if needed.")

        ff_dir = os.path.join(self.output_dir, "forcefield")
        os.makedirs(ff_dir, exist_ok=True)

        acpype_dirs = [d for d in os.listdir(ff_dir)
                      if (d.endswith('.acpype') or d.endswith('.amb2gmx'))
                      and os.path.isdir(os.path.join(ff_dir, d))]

        if acpype_dirs and not self.overwrite:
            print(f"Found existing acpype output, skipping generation")
            return ff_dir
        
        amber_mol2 = os.path.join(ff_dir, f"{self.ligand_name}.mol2")
        prep_file = os.path.join(ff_dir, f"{self.ligand_name}.prep")
        frcmod_file = os.path.join(ff_dir, f"{self.ligand_name}.frcmod")
        prmtop_file = os.path.join(ff_dir, f"{self.ligand_name}.prmtop")
        rst7_file = os.path.join(ff_dir, f"{self.ligand_name}.rst7")
        
        if self.file_format == 'mol2':
            self._process_mol2_format(amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file)
        else:
            # For SDF files, use direct antechamber conversion
            # This is the only method implemented in gaff2.py
            print("Processing SDF format...")
            self._process_sdf_format_direct(amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file)
        
        print("Converting to GROMACS format...")
        current_dir = os.getcwd()
        os.chdir(ff_dir)
        
        self.run_command([
            "acpype",
            "-p", os.path.basename(prmtop_file),
            "-x", os.path.basename(rst7_file),
            "-d"
        ])
        
        os.chdir(current_dir)
        
        return ff_dir
    
    # Include all the other necessary methods from the original file
    # (_process_mol2_format, _process_sdf_format_*, _create_topology_with_tleap, etc.)
    # These are omitted here for brevity but should be included in the actual implementation
    
    def _process_mol2_format(self, amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file):
        """Process MOL2 format files"""
        print("Processing MOL2 format...")

        # Get the ff_dir from amber_mol2 path
        ff_dir = os.path.dirname(amber_mol2)

        # Clean up any existing antechamber temp files before starting
        self._clean_antechamber_temp_files(ff_dir)

        # Copy input file to working directory to avoid path issues
        input_mol2 = os.path.join(ff_dir, f"input_{self.ligand_name}.mol2")
        shutil.copy2(self.ligand_path, input_mol2)

        print("Generating AM1-BCC charges...")
        cmd = [
            "antechamber",
            "-i", os.path.basename(input_mol2),
            "-fi", "mol2",
            "-o", os.path.basename(amber_mol2),
            "-fo", "mol2",
            "-c", "bcc",
            "-s", "2",
            "-at", "gaff"  # Explicitly specify GAFF atom types
        ]
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])

        try:
            self.run_command(cmd, cwd=ff_dir)
        except subprocess.CalledProcessError as e:
            print("\nAM1-BCC charge generation failed. Attempting fallback with gas phase charges...")
            # Try with simpler gas phase charges as fallback
            cmd_fallback = [
                "antechamber",
                "-i", os.path.basename(input_mol2),
                "-fi", "mol2",
                "-o", os.path.basename(amber_mol2),
                "-fo", "mol2",
                "-c", "gas",
                "-s", "2",
                "-at", "gaff"
            ]
            if self.net_charge != 0:
                cmd_fallback.extend(["-nc", str(self.net_charge)])
            self.run_command(cmd_fallback, cwd=ff_dir)

        print("Generating prep file...")
        # Generate prep file from the charged mol2 file
        # Important: We don't need to recalculate charges, just convert format
        cmd = [
            "antechamber",
            "-i", os.path.basename(amber_mol2),
            "-fi", "mol2",
            "-o", os.path.basename(prep_file),
            "-fo", "prepi",
            "-at", "gaff"
        ]
        # Add net charge if specified (required for proper file generation)
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])

        try:
            self.run_command(cmd, cwd=ff_dir)
        except subprocess.CalledProcessError as e:
            print("\nPrep file generation failed. Trying with explicit charge retention...")
            # Alternative: explicitly tell antechamber to read charges from mol2
            cmd_alt = [
                "antechamber",
                "-i", os.path.basename(amber_mol2),
                "-fi", "mol2",
                "-o", os.path.basename(prep_file),
                "-fo", "prepi",
                "-c", "rc",  # Read charges from input file
                "-at", "gaff"
            ]
            if self.net_charge != 0:
                cmd_alt.extend(["-nc", str(self.net_charge)])
            self.run_command(cmd_alt, cwd=ff_dir)

        print("Generating force field parameters...")
        self.run_command([
            "parmchk2",
            "-i", os.path.basename(prep_file),
            "-f", "prepi",
            "-o", os.path.basename(frcmod_file),
            "-a", "Y"  # Print all force field parameters
        ], cwd=ff_dir)

        print("Creating AMBER topology...")
        self._create_topology_with_tleap(amber_mol2, frcmod_file, prmtop_file, rst7_file)

        # Clean up intermediate input file
        try:
            os.remove(input_mol2)
        except:
            pass
    
    def _create_topology_with_tleap(self, amber_mol2, frcmod_file, prmtop_file, rst7_file):
        """Create AMBER topology using tleap"""
        ff_dir = os.path.dirname(amber_mol2)
        tleap_input = os.path.join(ff_dir, f"tleap_{self.ligand_name}.in")

        with open(tleap_input, 'w') as f:
            f.write(f"""source leaprc.protein.ff14SB
source leaprc.gaff

LIG = loadmol2 {os.path.basename(amber_mol2)}
loadamberparams {os.path.basename(frcmod_file)}

saveamberparm LIG {os.path.basename(prmtop_file)} {os.path.basename(rst7_file)}

quit
""")

        self.run_command(["tleap", "-f", os.path.basename(tleap_input)], cwd=ff_dir)
    
    def find_acpype_output(self, ff_dir):
        """Find the acpype output directory and files"""
        print("\n=== Locating acpype output ===")
        
        possible_dirs = []
        for item in os.listdir(ff_dir):
            item_path = os.path.join(ff_dir, item)
            if os.path.isdir(item_path) and (item.endswith('.acpype') or item.endswith('.amb2gmx')):
                possible_dirs.append(item_path)
                mol_name = item.split('.')[0]
                if mol_name not in self.possible_molecule_names:
                    self.possible_molecule_names.append(mol_name)
        
        if not possible_dirs:
            raise RuntimeError("No acpype output directory found")
        
        acpype_dir = possible_dirs[0]
        print(f"Found acpype directory: {acpype_dir}")
        
        files = {}
        for filename in os.listdir(acpype_dir):
            filepath = os.path.join(acpype_dir, filename)
            if filename.endswith('_GMX.itp'):
                files['itp'] = filepath
            elif filename.endswith('_GMX.gro'):
                files['gro'] = filepath
            elif filename.endswith('_GMX.top'):
                files['top'] = filepath
            elif filename.startswith('posre_') and filename.endswith('.itp'):
                files['posre'] = filepath
        
        print(f"Found files: {list(files.keys())}")
        return files
    
    def standardize_to_LIG(self, ff_files):
        """Standardize all files to use 'LIG' as molecule name"""
        print("\n=== Standardizing to LIG naming ===")
        
        lig_dir = os.path.join(self.output_dir, "LIG.amb2gmx")
        
        if os.path.exists(lig_dir) and not self.overwrite:
            print(f"LIG.amb2gmx already exists, skipping standardization")
            return lig_dir
        
        os.makedirs(lig_dir, exist_ok=True)
        
        std_files = {}
        
        if 'gro' in ff_files and ff_files['gro']:
            std_files['gro'] = os.path.join(lig_dir, "LIG.gro")
            self._standardize_file(ff_files['gro'], std_files['gro'])
        
        if 'top' in ff_files and ff_files['top']:
            std_files['top'] = os.path.join(lig_dir, "LIG.top")
            self._standardize_file(ff_files['top'], std_files['top'])
        
        std_files['itp'] = os.path.join(lig_dir, "LIG.itp")
        if 'itp' in ff_files and ff_files['itp']:
            self._standardize_file(ff_files['itp'], std_files['itp'], add_posre=True)
        elif 'top' in ff_files and ff_files['top']:
            self._extract_itp_from_top(ff_files['top'], std_files['itp'])
        
        if 'posre' in ff_files and ff_files['posre']:
            std_files['posre'] = os.path.join(lig_dir, "posre_LIG.itp")
            self._standardize_file(ff_files['posre'], std_files['posre'])
        
        if 'top' in ff_files and ff_files['top']:
            std_files['atomtypes'] = os.path.join(lig_dir, "atomtypes_LIG.itp")
            self._extract_atomtypes(ff_files['top'], std_files['atomtypes'])
        
        print(f"Standardization complete. Files in: {lig_dir}")
        return lig_dir
    
    def _standardize_file(self, source, target, add_posre=False):
        """Replace all molecule names with LIG in a file"""
        with open(source, 'r') as f:
            content = f.read()
        
        for mol_name in self.possible_molecule_names:
            if mol_name != "LIG":
                content = re.sub(r'\b{}\b'.format(re.escape(mol_name)), "LIG", content)
        
        if add_posre and "#ifdef POSRES" not in content:
            content += "\n#ifdef POSRES\n"
            content += '#include "posre_LIG.itp"\n'
            content += "#endif\n"
        
        with open(target, 'w') as f:
            f.write(content)
    
    def _extract_itp_from_top(self, top_file, itp_file):
        """Extract moleculetype section from TOP file to create ITP"""
        with open(top_file, 'r') as f:
            content = f.read()
        
        mol_match = re.search(r'\[ *moleculetype *\](.*?)(\[ *system|$)', content, re.DOTALL)
        
        if mol_match:
            mol_content = "[ moleculetype ]\n" + mol_match.group(1).strip()
            
            for mol_name in self.possible_molecule_names:
                if mol_name != "LIG":
                    mol_content = re.sub(r'\b{}\b'.format(re.escape(mol_name)), "LIG", mol_content)
            
            mol_content += "\n\n#ifdef POSRES\n"
            mol_content += '#include "posre_LIG.itp"\n'
            mol_content += "#endif\n"
            
            with open(itp_file, 'w') as f:
                f.write(mol_content)
        else:
            print("Warning: Could not extract moleculetype from TOP file")
    
    def _extract_atomtypes(self, top_file, atomtypes_file):
        """Extract atomtypes section from TOP file"""
        with open(top_file, 'r') as f:
            content = f.read()
        
        atomtypes_match = re.search(r'\[ *atomtypes *\](.*?)(?=\[ *[a-zA-Z]|\Z)', content, re.DOTALL)
        
        if atomtypes_match:
            atomtypes_content = "[ atomtypes ]\n" + atomtypes_match.group(1).strip() + "\n"
            with open(atomtypes_file, 'w') as f:
                f.write(atomtypes_content)
        else:
            print("Warning: Could not extract atomtypes from TOP file")
    
    def _clean_antechamber_temp_files(self, directory):
        """Clean up antechamber temporary files from a specific directory"""
        temp_patterns = ["ANTECHAMBER*", "ATOMTYPE.INF", "PREP.INF", "NEWPDB.PDB", "sqm.*"]
        for pattern in temp_patterns:
            for file_path in Path(directory).glob(pattern):
                try:
                    os.remove(file_path)
                except OSError:
                    pass

    def cleanup_temp_files(self):
        """Clean up temporary files"""
        print("\n=== Cleaning up ===")

        ff_dir = os.path.join(self.output_dir, "forcefield")
        if os.path.exists(ff_dir):
            patterns = ["*.frcmod", "*.prep", "*.prmtop", "*.rst7", "*.log", "*.in",
                       "ANTECHAMBER*", "ATOMTYPE*", "PREP*", "NEWPDB*", "sqm*", "input_*"]

            for pattern in patterns:
                for file_path in Path(ff_dir).glob(pattern):
                    try:
                        os.remove(file_path)
                    except OSError:
                        pass

        # Also clean up any antechamber temp files in current working directory
        # (These should not be created if cwd is set correctly, but clean them just in case)
        cwd = os.getcwd()
        self._clean_antechamber_temp_files(cwd)
    
    # Add remaining methods from original file as needed
    def _process_sdf_format_direct(self, amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file):
        """Process SDF format files with direct antechamber conversion"""
        print("Processing SDF format (direct method)...")
        
        cmd = [
            "antechamber",
            "-i", self.ligand_path,
            "-fi", "sdf",
            "-o", amber_mol2,
            "-fo", "mol2",
            "-c", "bcc",
            "-s", "2"
        ]
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])
        self.run_command(cmd)
        
        self._process_mol2_intermediate(amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file)
    
    def _process_mol2_intermediate(self, amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file, charge_method="wc"):
        """Process intermediate MOL2 file"""
        print("Generating prep file from charged MOL2...")
        cmd = [
            "antechamber",
            "-i", amber_mol2,
            "-fi", "mol2",
            "-o", prep_file,
            "-fo", "prepi",
            "-c", charge_method,
            "-s", "2"
        ]
        self.run_command(cmd)
        
        print("Generating force field parameters...")
        self.run_command([
            "parmchk2",
            "-i", prep_file,
            "-f", "prepi",
            "-o", frcmod_file
        ])
        
        print("Creating AMBER topology...")
        self._create_topology_with_tleap(amber_mol2, frcmod_file, prmtop_file, rst7_file)