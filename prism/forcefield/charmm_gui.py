#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CHARMM-GUI Force Field Generator for PRISM

This module processes CHARMM-GUI generated GROMACS files and converts them
to PRISM-standardized format (consistent with other force field generators).

CHARMM-GUI already generates GROMACS format files (ITP + structure files),
but they need to be reorganized to match PRISM's standard output structure.
"""

import os
import shutil
from pathlib import Path

try:
    from .base import ForceFieldGeneratorBase
except ImportError:
    from base import ForceFieldGeneratorBase


class CHARMMGUIForceFieldGenerator(ForceFieldGeneratorBase):
    """
    CHARMM-GUI force field generator - processes CHARMM-GUI output files

    CHARMM-GUI generates:
    - gromacs/LIG.itp (GROMACS topology)
    - gromacs/charmm36.itp (CHARMM36 force field)
    - gromacs/topol.top (topology file)
    - Various structure files (PDB, GRO, etc.)

    This converter reorganizes these files into PRISM-standard format:
    - LIG.charmm2gmx/LIG.itp
    - LIG.charmm2gmx/LIG.gro
    - LIG.charmm2gmx/atomtypes_LIG.itp
    - etc.
    """

    def __init__(self, ligand_path, output_dir, charmm_gui_dir=None, overwrite=False):
        """
        Initialize CHARMM-GUI Force Field Generator

        Parameters
        ----------
        ligand_path : str
            Path to ligand file (kept for API compatibility)
        output_dir : str
            Output directory for processed files
        charmm_gui_dir : str
            Path to directory containing CHARMM-GUI output (gromacs/ folder)
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(ligand_path, output_dir, overwrite)

        if charmm_gui_dir is None:
            raise ValueError(
                "charmm_gui_dir must be provided. " "Use --forcefield-path to specify the CHARMM-GUI output directory."
            )

        self.charmm_gui_dir = Path(charmm_gui_dir)

        if not self.charmm_gui_dir.exists():
            raise FileNotFoundError(f"CHARMM-GUI directory does not exist: {self.charmm_gui_dir}")

        # Look for gromacs/ subdirectory
        self.gromacs_dir = self.charmm_gui_dir / "gromacs"
        if not self.gromacs_dir.exists():
            # Try using the directory directly if it already contains gromacs files
            if (self.charmm_gui_dir / "LIG.itp").exists():
                self.gromacs_dir = self.charmm_gui_dir
            else:
                raise FileNotFoundError(
                    f"Cannot find gromacs/ directory in: {self.charmm_gui_dir}\n"
                    f"Expected structure: charmm_gui_dir/gromacs/LIG.itp"
                )

        self.moleculetype_name = "LIG"

        print(f"\nInitialized CHARMM-GUI Force Field Generator:")
        print(f"  CHARMM-GUI directory: {self.charmm_gui_dir}")
        print(f"  GROMACS directory: {self.gromacs_dir}")
        print(f"  Output directory: {self.output_dir}")

    def get_output_dir_name(self):
        """Get the output directory name for CHARMM-GUI"""
        return "LIG.charmm2gmx"

    def run(self):
        """Run the CHARMM-GUI conversion workflow"""
        print(f"\n{'='*60}")
        print("Starting CHARMM-GUI Force Field Processing")
        print(f"{'='*60}")

        try:
            # Check if output already exists
            lig_dir = os.path.join(self.output_dir, "LIG.charmm2gmx")
            if os.path.exists(lig_dir) and not self.overwrite:
                if self.check_required_files(lig_dir):
                    print(f"\nUsing cached CHARMM-GUI force field parameters from: {lig_dir}")
                    print("All required files found:")
                    for f in sorted(os.listdir(lig_dir)):
                        print(f"  - {f}")
                    print("\n(Use --overwrite to regenerate)")
                    return lig_dir

            # Step 1: Find input files
            print("\n[Step 1] Finding CHARMM-GUI input files...")
            self._find_input_files()

            # Step 2: Read CHARMM-GUI ITP file
            print("\n[Step 2] Reading CHARMM-GUI topology file...")
            itp_content = self._read_itp_file()
            num_atoms = self._count_atoms_in_itp(itp_content)
            print(f"  - Found {num_atoms} atoms")

            # Step 3: Generate or convert GRO file
            print("\n[Step 3] Generating GRO file...")
            gro_content = self._generate_or_convert_gro(num_atoms)
            print(f"  - Generated GRO file with {num_atoms} atoms")

            # Step 4: Extract atomtypes
            print("\n[Step 4] Extracting atom types...")
            atomtypes_content = self._extract_atomtypes(itp_content)
            num_atomtypes = len(
                [line for line in atomtypes_content.split("\n") if line.strip() and not line.startswith(";")]
            )
            print(f"  - Found {num_atomtypes} atom types")

            # Step 5: Create output directory
            print(f"\n[Step 5] Creating output directory: {lig_dir}")
            os.makedirs(lig_dir, exist_ok=True)

            # Step 6: Copy CHARMM36 force field if exists
            charmm36_src = self.gromacs_dir / "charmm36.itp"
            if charmm36_src.exists():
                charmm36_dst = Path(lig_dir) / "charmm36.ff"
                print(f"\n[Step 6.5] Setting up CHARMM36 force field...")
                if charmm36_dst.exists():
                    shutil.rmtree(charmm36_dst)
                charmm36_dst.mkdir(exist_ok=True)
                shutil.copy2(charmm36_src, charmm36_dst / "charmm36.itp")
                print(f"  - Copied CHARMM36 force field")

            # Step 7: Generate PRISM-format files
            print("\n[Step 7] Generating PRISM-formatted files...")

            # LIG.gro
            gro_file = os.path.join(lig_dir, "LIG.gro")
            with open(gro_file, "w") as f:
                f.write(gro_content)
            print(f"  - Generated: LIG.gro")

            # LIG.itp (reorganized)
            itp_file = os.path.join(lig_dir, "LIG.itp")
            with open(itp_file, "w") as f:
                f.write(itp_content)
            print(f"  - Generated: LIG.itp")

            # atomtypes_LIG.itp
            atomtypes_file = os.path.join(lig_dir, "atomtypes_LIG.itp")
            with open(atomtypes_file, "w") as f:
                f.write(atomtypes_content)
            print(f"  - Generated: atomtypes_LIG.itp")

            # posre_LIG.itp
            posre_content = self._generate_posre_itp(num_atoms)
            posre_file = os.path.join(lig_dir, "posre_LIG.itp")
            with open(posre_file, "w") as f:
                f.write(posre_content)
            print(f"  - Generated: posre_LIG.itp")

            # LIG.top (topology)
            top_content = self._generate_top_file()
            top_file = os.path.join(lig_dir, "LIG.top")
            with open(top_file, "w") as f:
                f.write(top_content)
            print(f"  - Generated: LIG.top")

            print(f"\n{'='*60}")
            print("CHARMM-GUI force field processing completed successfully!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {lig_dir}")
            print("\nGenerated files:")
            for f in sorted(os.listdir(lig_dir)):
                print(f"  - {f}")

            return lig_dir

        except Exception as e:
            print(f"\nError during CHARMM-GUI force field processing: {e}")
            raise

    def _find_input_files(self):
        """Find required CHARMM-GUI files"""
        # Check for LIG.itp
        itp_file = self.gromacs_dir / "LIG.itp"
        if not itp_file.exists():
            raise FileNotFoundError(
                f"LIG.itp not found in: {self.gromacs_dir}\n" f"Please ensure CHARMM-GUI output is complete."
            )

        self.itp_file = itp_file
        print(f"  Found LIG.itp: {itp_file}")

        # Check for structure files (PDB, GRO, or others)
        self.pdb_file = None
        self.gro_file = None

        # Priority: GRO > ligandrm.pdb > other PDB > other structure files
        # ligandrm.pdb is the CHARMM-GUI processed PDB with correct atom names
        for ext in [".gro", ".pdb", ".crd"]:
            structure_files = list(self.gromacs_dir.glob(f"*{ext}"))
            if structure_files:
                if ext == ".gro":
                    self.gro_file = structure_files[0]
                    print(f"  Found structure file: {structure_files[0].name}")
                    break
                elif ext == ".pdb":
                    # Prefer ligandrm.pdb over other PDB files
                    # Sort to ensure ligandrm.pdb comes first if it exists
                    structure_files.sort(key=lambda f: 1 if "ligandrm" in f.name.lower() else 0)
                    self.pdb_file = structure_files[-1]  # Get the last one (ligandrm.pdb if exists)
                    print(f"  Found PDB file: {self.pdb_file.name}")
                    break

        # Also check parent directories for structure files
        if not self.gro_file and not self.pdb_file:
            for ext in [".gro", ".pdb"]:
                structure_files = list(self.charmm_gui_dir.glob(f"*{ext}"))
                if structure_files:
                    if ext == ".gro":
                        self.gro_file = structure_files[0]
                        print(f"  Found structure file in parent dir: {structure_files[0].name}")
                        break
                    elif ext == ".pdb":
                        # Prefer ligandrm.pdb over other PDB files
                        # Sort to ensure ligandrm.pdb comes first if it exists
                        structure_files.sort(key=lambda f: 1 if "ligandrm" in f.name.lower() else 0)
                    self.pdb_file = structure_files[-1]  # Get the last one (ligandrm.pdb if exists)
                    print(f"  Found PDB file in parent dir: {self.pdb_file.name}")
                    break

    def _read_itp_file(self):
        """Read and process CHARMM-GUI ITP file"""
        with open(self.itp_file, "r") as f:
            content = f.read()
        return content

    def _count_atoms_in_itp(self, itp_content):
        """Count atoms in [atoms] section"""
        in_atoms = False
        count = 0
        for line in itp_content.split("\n"):
            if "[ atoms ]" in line:
                in_atoms = True
                continue
            if in_atoms and line.strip().startswith("["):
                break
            if in_atoms and line.strip() and not line.strip().startswith(";"):
                count += 1
        return count

    def _generate_or_convert_gro(self, num_atoms):
        """Generate or convert to GRO format"""
        # If GRO file exists, use it
        if self.gro_file and self.gro_file.exists():
            print(f"  Using existing GRO file: {self.gro_file.name}")
            with open(self.gro_file, "r") as f:
                return f.read()

        # If PDB file exists, convert to GRO
        if self.pdb_file and self.pdb_file.exists():
            print(f"  Converting PDB to GRO: {self.pdb_file.name}")
            return self._convert_pdb_to_gro(self.pdb_file)

        # If no structure file found, generate minimal GRO
        print(f"  Warning: No structure file found, generating minimal GRO")
        return self._generate_minimal_gro(num_atoms)

    def _convert_pdb_to_gro(self, pdb_file):
        """Convert PDB to GRO format"""
        gro_lines = []
        atom_count = 0
        x_coords, y_coords, z_coords = [], [], []

        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_count += 1
                    atom_name = line[12:16].strip()
                    x = float(line[30:38].strip()) / 10.0  # Å to nm
                    y = float(line[38:46].strip()) / 10.0
                    z = float(line[46:54].strip()) / 10.0

                    x_coords.append(x)
                    y_coords.append(y)
                    z_coords.append(z)

                    gro_line = f"{1:5d}{self.moleculetype_name:<5s}{atom_name:>5s}{atom_count:5d}"
                    gro_line += f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
                    gro_lines.append(gro_line)

        # Calculate box dimensions
        box_x = max(x_coords) - min(x_coords) + 1.0 if x_coords else 1.0
        box_y = max(y_coords) - min(y_coords) + 1.0 if y_coords else 1.0
        box_z = max(z_coords) - min(z_coords) + 1.0 if z_coords else 1.0

        # Write GRO content
        gro_content = f"{self.moleculetype_name} system\n"
        gro_content += f"{atom_count:5d}\n"
        gro_content += "".join(gro_lines)
        gro_content += f"{box_x:10.7f}{box_y:10.7f}{box_z:10.7f}\n"

        return gro_content

    def _generate_minimal_gro(self, num_atoms):
        """Generate minimal GRO file when no structure is available"""
        # Place all atoms at origin with box size 1nm
        gro_content = f"{self.moleculetype_name} system\n"
        gro_content += f"{num_atoms:5d}\n"
        for i in range(1, num_atoms + 1):
            gro_content += f"{1:5d}{self.moleculetype_name:<5s}ATOM{i:4d}{i:5d}"
            gro_content += f"0.0000.0000.000\n"
        gro_content += f"1.00000001.00000001.0000000\n"
        return gro_content

    def _extract_atomtypes(self, itp_content):
        """Extract atom types from ITP file"""
        content = "[ atomtypes ]\n"
        content += ";type, bondingtype, atomic_number, mass, charge, ptype, sigma, epsilon\n"

        # Parse atoms section to get unique atom types
        in_atoms = False
        atom_types = {}

        for line in itp_content.split("\n"):
            if "[ atoms ]" in line:
                in_atoms = True
                continue
            if in_atoms and line.strip().startswith("["):
                break
            if in_atoms and line.strip() and not line.strip().startswith(";"):
                parts = line.split()
                if len(parts) >= 8:
                    atom_type = parts[1]
                    mass = parts[7]
                    # Store mass for this atom type
                    if atom_type not in atom_types:
                        # Guess atomic number from mass
                        atomic_num = self._guess_atomic_number(float(mass))
                        atom_types[atom_type] = (atomic_num, mass)

        # Write atomtypes
        for atom_type, (atomic_num, mass) in sorted(atom_types.items()):
            content += f"{atom_type:10s}{atomic_num:5d}{mass:12s}  0.0000000000000000  A    0.0  0.0\n"

        return content

    def _guess_atomic_number(self, mass):
        """Guess atomic number from mass"""
        if mass < 1.5:
            return 1  # H
        elif mass < 13:
            return 6  # C
        elif mass < 15:
            return 7  # N
        elif mass < 17:
            return 8  # O
        elif mass < 19:
            return 9  # F
        elif mass < 24:
            return 11  # Na
        elif mass < 28:
            return 12  # Mg
        elif mass < 32:
            return 15  # P
        elif mass < 36:
            return 16  # S
        elif mass < 40:
            return 17  # Cl
        else:
            return 6  # Default to C

    def _generate_posre_itp(self, num_atoms):
        """Generate position restraint file"""
        content = "[ position_restraints ]\n"
        content += "; atom  type      fx      fy      fz\n"

        for i in range(1, num_atoms + 1):
            content += f"{i:6d}     1  1000  1000  1000\n"

        return content

    def _generate_top_file(self):
        """Generate LIG.top file"""
        content = "; Include CHARMM36 force field\n"

        # Check if charmm36.ff exists
        charmm36_dir = Path(self.output_dir) / "LIG.charmm2gmx" / "charmm36.ff"
        if charmm36_dir.exists():
            content += '#include "charmm36.ff/charmm36.itp"\n'
        else:
            # Fallback: include from original location
            charmm36_src = self.gromacs_dir / "charmm36.itp"
            if charmm36_src.exists():
                # Copy charmm36.itp to output
                import shutil as shutil_mod

                charmm36_dst = Path(self.output_dir) / "LIG.charmm2gmx" / "charmm36.itp"
                shutil_mod.copy2(charmm36_src, charmm36_dst)
                content += '#include "charmm36.itp"\n'

        content += "\n"
        content += '#include "atomtypes_LIG.itp"\n'
        content += '#include "LIG.itp"\n\n'

        content += "[ system ]\n"
        content += ";name\n"
        content += f"{self.moleculetype_name} system\n\n"

        content += "[ molecules ]\n"
        content += ";name\tnumber\n"
        content += f"{self.moleculetype_name}              1\n"

        return content
