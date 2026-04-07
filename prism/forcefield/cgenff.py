#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CGenFF force field generator wrapper for PRISM

CGenFF (CHARMM General Force Field) requires downloading ligand files
from the CGenFF web server: https://cgenff.com/

This module processes the downloaded CGenFF files and converts them
to PRISM-standardized GROMACS format.
"""

import os
import shutil
import re
from pathlib import Path
from typing import Dict, List, Tuple

# Import the base class
try:
    from .base import ForceFieldGeneratorBase
except ImportError:
    from base import ForceFieldGeneratorBase


class CGenFFForceFieldGenerator(ForceFieldGeneratorBase):
    """CGenFF force field generator - processes web-downloaded CGenFF files"""

    def __init__(self, ligand_path, output_dir, cgenff_dir=None, overwrite=False):
        """
        Initialize CGenFF force field generator

        Parameters
        ----------
        ligand_path : str
            Path to ligand file (not used for CGenFF, kept for API compatibility)
        output_dir : str
            Output directory for processed files
        cgenff_dir : str
            Path to directory containing downloaded CGenFF files (PDB and TOP)
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(ligand_path, output_dir, overwrite)

        if cgenff_dir is None:
            raise ValueError("cgenff_dir must be provided. Use --forcefield-path to specify the CGenFF directory.")

        self.cgenff_dir = Path(cgenff_dir)

        if not self.cgenff_dir.exists():
            raise FileNotFoundError(f"CGenFF directory does not exist: {self.cgenff_dir}")

        # Data storage
        self.atoms = []
        self.bonds = []
        self.pairs = []
        self.angles = []
        self.dihedrals_proper = []
        self.dihedrals_improper = []
        self.atomtypes = {}
        self.moleculetype_name = "LIG"

        # File paths
        self.pdb_file = None
        self.top_file = None

        print(f"\nInitialized CGenFF Force Field Generator:")
        print(f"  CGenFF directory: {self.cgenff_dir}")
        print(f"  Output directory: {self.output_dir}")

    def get_output_dir_name(self):
        """Get the output directory name for CGenFF"""
        return "LIG.cgenff2gmx"

    @property
    def capabilities(self):
        """CGenFF force field capabilities"""
        return {
            "supports_fep": True,
            "supports_pmf": True,
            "requires_external": True,  # Requires files from cgenff.com
            "has_protein_ff": False,
            "charge_method": "CGenFF",
        }

    def run(self):
        """Run the CGenFF force field conversion workflow"""
        print(f"\n{'='*60}")
        print("Starting CGenFF Force Field Processing")
        print(f"{'='*60}")

        try:
            # Check if output already exists
            lig_dir = os.path.join(self.output_dir, "LIG.cgenff2gmx")
            if os.path.exists(lig_dir) and not self.overwrite:
                if self.check_required_files(lig_dir):
                    print(f"\nUsing cached CGenFF force field parameters from: {lig_dir}")
                    print("All required files found:")
                    for f in sorted(os.listdir(lig_dir)):
                        print(f"  - {f}")
                    print("\n(Use --overwrite to regenerate)")
                    return lig_dir

            # Find input files
            print("\n[Step 1] Finding CGenFF input files...")
            self._find_input_files()

            # Parse PDB
            print("\n[Step 2] Converting PDB to GRO format...")
            gro_content, num_atoms = self._parse_pdb()
            print(f"  - Found {num_atoms} atoms")

            # Parse topology
            print("\n[Step 3] Parsing topology file...")
            self._parse_topology()
            print(f"  - Atoms: {len(self.atoms)}")
            print(f"  - Bonds: {len(self.bonds)}")
            print(f"  - Pairs: {len(self.pairs)}")
            print(f"  - Angles: {len(self.angles)}")
            print(f"  - Proper dihedrals: {len(self.dihedrals_proper)}")
            print(f"  - Improper dihedrals: {len(self.dihedrals_improper)}")

            # Clean up lone pairs (LP atoms) for halogens
            print("\n[Step 3.5] Cleaning lone pairs (LP atoms)...")
            self._remove_lone_pairs()

            # Create output directory
            print(f"\n[Step 4] Creating output directory: {lig_dir}")
            os.makedirs(lig_dir, exist_ok=True)

            # Copy charmm36.ff directory if it exists
            charmm_ff_src = self.cgenff_dir / "charmm36.ff"
            charmm_ff_dst = Path(lig_dir) / "charmm36.ff"
            if charmm_ff_src.exists():
                print(f"\n[Step 4.5] Copying CHARMM36 force field...")
                if charmm_ff_dst.exists():
                    shutil.rmtree(charmm_ff_dst)
                shutil.copytree(charmm_ff_src, charmm_ff_dst)
                print(f"  Copied CHARMM36 force field to output directory")

            # Generate files
            print("\n[Step 5] Generating PRISM-formatted files...")

            # LIG.gro
            gro_file = os.path.join(lig_dir, "LIG.gro")
            with open(gro_file, "w") as f:
                f.write(gro_content)
            print(f"  - Generated: LIG.gro")

            # atomtypes_LIG.itp
            atomtypes_content = self._generate_atomtypes_itp()
            atomtypes_file = os.path.join(lig_dir, "atomtypes_LIG.itp")
            with open(atomtypes_file, "w") as f:
                f.write(atomtypes_content)
            print(f"  - Generated: atomtypes_LIG.itp")

            # LIG.itp
            itp_content = self._generate_lig_itp()
            itp_file = os.path.join(lig_dir, "LIG.itp")
            with open(itp_file, "w") as f:
                f.write(itp_content)
            print(f"  - Generated: LIG.itp")

            # LIG.top
            top_content = self._generate_lig_top()
            top_file = os.path.join(lig_dir, "LIG.top")
            with open(top_file, "w") as f:
                f.write(top_content)
            print(f"  - Generated: LIG.top")

            # posre_LIG.itp
            posre_content = self._generate_posre_itp(num_atoms)
            posre_file = os.path.join(lig_dir, "posre_LIG.itp")
            with open(posre_file, "w") as f:
                f.write(posre_content)
            print(f"  - Generated: posre_LIG.itp")

            print(f"\n{'='*60}")
            print("CGenFF force field processing completed successfully!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {lig_dir}")
            print("\nGenerated files:")
            for f in sorted(os.listdir(lig_dir)):
                print(f"  - {f}")

            return lig_dir

        except Exception as e:
            print(f"\nError during CGenFF force field processing: {e}")
            raise

    def _find_input_files(self):
        """Find PDB and TOP files in CGenFF directory

        Prioritizes *_gmx.pdb files which are generated by CGenFF website
        and have atom names matching the topology files.
        """
        # Find topology file (TOP or RTF)
        top_files = list(self.cgenff_dir.glob("*_gmx.top"))
        if not top_files:
            top_files = list(self.cgenff_dir.glob("*.top"))
        if not top_files:
            # Try RTF files (CHARMM-GUI format)
            rtf_files = list(self.cgenff_dir.glob("*.rtf"))
            if not rtf_files:
                # Check lig/ subdirectory for RTF
                lig_dir = self.cgenff_dir / "lig"
                if lig_dir.exists():
                    rtf_files = list(lig_dir.glob("*.rtf"))
            if rtf_files:
                top_files = rtf_files
        if not top_files:
            # Last resort: try ITP files
            itp_files = list(self.cgenff_dir.glob("LIG.itp"))
            if itp_files:
                top_files = itp_files

        # Special handling for CHARMM-GUI LIG.top files that contain #include statements
        # These are full topology files, not molecule-only ITP files
        for top_file in top_files:
            with open(top_file, "r") as f:
                content = f.read()
                # Check if it contains #include "LIG.itp" and no [ atoms ] section
                if '#include "LIG.itp"' in content and "[ atoms ]" not in content:
                    # This is a full topology file, use LIG.itp instead
                    itp_file = self.cgenff_dir / "LIG.itp"
                    if itp_file.exists():
                        top_files = [itp_file]
                        print(f"  Detected CHARMM-GUI full topology file, using LIG.itp instead")
                        break

        if not top_files:
            raise FileNotFoundError(
                f"No topology file found in {self.cgenff_dir}\n"
                f"  Expected one of:\n"
                f"    - *_gmx.top from CGenFF website\n"
                f"    - *.rtf from CHARMM-GUI\n"
                f"    - lig/lig.rtf from CHARMM-GUI\n"
                f"    - LIG.itp from PRISM"
            )

        self.top_file = top_files[0]
        print(f"  Found topology file: {self.top_file.name}")

        # Find PDB file with matching atom names
        gmx_pdb_files = list(self.cgenff_dir.glob("*_gmx.pdb"))

        if gmx_pdb_files:
            # CGenFF website files - use first one
            self.pdb_file = gmx_pdb_files[0]
            print(f"  Found CGenFF-generated PDB file: {self.pdb_file.name}")
        else:
            # CHARMM-GUI or other: find PDB that matches topology
            pdb_files = list(self.cgenff_dir.glob("*.pdb"))
            if pdb_files:
                # Try to find PDB file with matching atom names
                matching_pdb = self._find_matching_pdb(pdb_files, self.top_file)
                if matching_pdb:
                    self.pdb_file = matching_pdb
                    print(f"  Found matching PDB file: {self.pdb_file.name}")
                else:
                    # No matching PDB found - use first but warn
                    self.pdb_file = pdb_files[0]
                    print(f"  Warning: No PDB file matches topology atom names")
                    print(f"  Using: {self.pdb_file.name}")
            else:
                # Last resort: try GRO files (CHARMM-GUI format)
                gro_files = list(self.cgenff_dir.glob("*.gro"))
                if gro_files:
                    self.pdb_file = gro_files[0]
                    print(f"  Found GRO file (will convert): {self.pdb_file.name}")
                else:
                    raise FileNotFoundError(
                        f"No coordinate file found in {self.cgenff_dir}\n"
                        f"  Expected: *_gmx.pdb, *.pdb, or *.gro (e.g., ligandrm.pdb)"
                    )
                print(f"  Note: This may cause atom mismatches")

        # Validate atom name consistency
        self._validate_pdb_top_consistency()

    def _find_matching_pdb(self, pdb_files, top_file):
        """Find PDB file whose atom names match TOP/ITP/RTF file

        Parameters
        ----------
        pdb_files : list
            List of PDB file paths to check
        top_file : Path
            TOP, ITP, or RTF file to match against

        Returns
        -------
        Path or None
            Matching PDB file, or None if no match found
        """
        # Read atom names from topology file (TOP, ITP, or RTF)
        top_atoms = set()

        # Check if it's an RTF file (has ATOM lines without [ atoms ] header)
        is_rtf = False
        with open(top_file, "r") as f:
            first_lines = [f.readline() for _ in range(20)]
            if any(line.startswith("ATOM") for line in first_lines):
                is_rtf = True

        if is_rtf:
            # Parse RTF format
            with open(top_file, "r") as f:
                for line in f:
                    if line.startswith("ATOM"):
                        parts = line.split()
                        if len(parts) >= 4:
                            atom_name = parts[1]
                            # Skip LP atoms
                            if "LP" not in atom_name.upper():
                                top_atoms.add(atom_name)
        else:
            # Parse TOP/ITP format ([ atoms ] section)
            with open(top_file, "r") as f:
                in_atoms = False
                for line in f:
                    if "[ atoms ]" in line.lower():
                        in_atoms = True
                        continue
                    if in_atoms:
                        if line.strip().startswith("["):
                            break
                        if line.strip().startswith(";"):
                            continue
                        parts = line.split()
                        if len(parts) >= 5:
                            try:
                                int(parts[0])  # Check if first column is atom index
                                atom_name = parts[4]
                                # Skip LP atoms
                                if "LP" not in atom_name.upper():
                                    top_atoms.add(atom_name)
                            except ValueError:
                                pass

        if not top_atoms:
            print(f"  Warning: Could not parse atom names from {top_file.name}")
            return None

        print(f"  Looking for PDB with {len(top_atoms)} matching atoms...")

        # Check each PDB file
        for pdb_file in pdb_files:
            pdb_atoms = set()
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atom_name = line[12:16].strip()
                        # Skip LP atoms
                        if "LP" not in atom_name.upper():
                            pdb_atoms.add(atom_name)

            # Check match
            if pdb_atoms == top_atoms:
                return pdb_file

            # Report mismatches for close matches
            if len(pdb_atoms & top_atoms) > len(top_atoms) * 0.8:  # 80% overlap
                missing = top_atoms - pdb_atoms
                extra = pdb_atoms - top_atoms
                print(f"    {pdb_file.name}: {len(pdb_atoms & top_atoms)}/{len(top_atoms)} match", end="")
                if missing:
                    miss_list = sorted(list(missing))[:3]
                    print(f" (missing: {miss_list}...)" if len(missing) > 3 else f" (missing: {miss_list})", end="")
                if extra:
                    extra_list = sorted(list(extra))[:3]
                    print(f" (extra: {extra_list}...)" if len(extra) > 3 else f" (extra: {extra_list})", end="")
                print()

        return None

    def _validate_pdb_top_consistency(self):
        """Validate that atom names in PDB match those in TOP file

        This is critical for CGenFF - the website generates *_gmx.pdb files
        with atom names that match the topology. Using different PDB files
        (e.g., original MOL2 converted to PDB) will cause parameter errors.
        """
        # Read atom names from PDB
        pdb_atom_names = []
        with open(self.pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()
                    # Skip LP atoms as they'll be removed
                    if "LP" not in atom_name.upper():
                        pdb_atom_names.append(atom_name)

        # Read atom names from TOP file [atoms] section
        top_atom_names = []
        with open(self.top_file, "r") as f:
            content = f.read()

        section = self._extract_section(content, "atoms")
        if section:
            for line in section.split("\n"):
                line = line.strip()
                if not line or line.startswith(";"):
                    continue

                parts = line.split()
                if len(parts) >= 5:
                    atom_name = parts[4]  # atom name is 5th column
                    # Skip LP atoms as they'll be removed
                    if "LP" not in atom_name.upper():
                        top_atom_names.append(atom_name)

        # Compare counts
        if len(pdb_atom_names) != len(top_atom_names):
            print(f"\n  WARNING: Atom count mismatch!")
            print(f"    PDB file: {len(pdb_atom_names)} atoms (excluding LP)")
            print(f"    TOP file: {len(top_atom_names)} atoms (excluding LP)")
            print(f"    This may indicate atom name inconsistency between files")

        # Sample first few atom names to check consistency
        sample_size = min(5, len(pdb_atom_names), len(top_atom_names))
        mismatch_found = False

        for i in range(sample_size):
            if pdb_atom_names[i] != top_atom_names[i]:
                if not mismatch_found:
                    print(f"\n  WARNING: Atom name mismatch detected!")
                    print(f"    Index | PDB      | TOP")
                    print(f"    ------|----------|----------")
                    mismatch_found = True
                print(f"    {i+1:5d} | {pdb_atom_names[i]:8s} | {top_atom_names[i]:8s}")

        if mismatch_found:
            print(f"\n  ERROR: Atom names in PDB and TOP files do not match!")
            print(f"  This will cause missing parameter errors in GROMACS.")
            print(f"\n  Solution:")
            print(f"    1. Ensure you're using the *_gmx.pdb file from CGenFF website")
            print(f"    2. Do NOT use the original MOL2 file converted to PDB")
            print(f"    3. The CGenFF website generates matching PDB and TOP files")
            raise ValueError("Atom name mismatch between PDB and TOP files")
        else:
            print(f"  Atom names validated: {len(pdb_atom_names)} atoms match between PDB and TOP")

    def _parse_pdb(self):
        """Parse coordinate file (PDB or GRO) and convert to GRO format (excluding LP atoms)"""
        # Detect file format
        is_gro = self.pdb_file.suffix.lower() == ".gro"

        if is_gro:
            return self._parse_gro_file()
        else:
            return self._parse_pdb_file()

    def _parse_pdb_file(self):
        """Parse PDB file and convert to GRO format (excluding LP atoms)"""
        gro_lines = []
        atom_count = 0

        # Box size calculation
        x_coords, y_coords, z_coords = [], [], []

        with open(self.pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Parse PDB line
                    atom_name = line[12:16].strip()

                    # Skip LP (lone pair) atoms - they will be removed from topology
                    if "LP" in atom_name.upper():
                        continue

                    atom_count += 1
                    x = float(line[30:38].strip()) / 10.0  # A to nm
                    y = float(line[38:46].strip()) / 10.0
                    z = float(line[46:54].strip()) / 10.0

                    x_coords.append(x)
                    y_coords.append(y)
                    z_coords.append(z)

                    # GRO format: resnum, resname, atomname, atomnum, x, y, z
                    gro_line = f"{1:5d}{self.moleculetype_name:<5s}{atom_name:>5s}{atom_count:5d}"
                    gro_line += f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
                    gro_lines.append(gro_line)

        # Calculate box dimensions (with 1.0 nm padding)
        box_x = max(x_coords) - min(x_coords) + 1.0
        box_y = max(y_coords) - min(y_coords) + 1.0
        box_z = max(z_coords) - min(z_coords) + 1.0

        # Write GRO file
        gro_content = f"{self.moleculetype_name} system\n"
        gro_content += f"{atom_count:5d}\n"
        gro_content += "".join(gro_lines)
        gro_content += f"{box_x:10.7f}{box_y:10.7f}{box_z:10.7f}\n"

        return gro_content, atom_count

    def _parse_gro_file(self):
        """Parse GRO file and convert to GRO format (excluding LP atoms)"""
        gro_lines = []
        atom_count = 0

        # Box size calculation
        x_coords, y_coords, z_coords = [], [], []

        with open(self.pdb_file, "r") as f:
            lines = f.readlines()

            # Skip header line (title)
            # Second line has atom count
            if len(lines) < 2:
                raise ValueError(f"GRO file too short: {self.pdb_file}")

            # Parse atoms (start from line 2, until box dimensions)
            for i, line in enumerate(lines[2:], start=2):
                # Check if we've reached the box dimensions line
                stripped = line.strip()
                if not stripped or len(stripped.split()) < 3:
                    # This is likely the box dimensions line
                    break

                parts = stripped.split()
                if len(parts) < 6:
                    continue

                # GRO format can have two forms:
                # 1. resnum and resname separate: "1 LIG N 1 ..."
                # 2. resnum and resname combined: "1LIG N 1 ..."
                if len(parts) >= 7:
                    # Try to parse as separate fields
                    try:
                        # First field might be resnum only
                        res_num = int(parts[0])
                        res_name = parts[1]
                        atom_name = parts[2]
                        atom_num = int(parts[3])
                        x = float(parts[4])
                        y = float(parts[5])
                        z = float(parts[6])
                    except ValueError:
                        # Combined format: extract resnum from start
                        first_field = parts[0]
                        # Find where numbers end and letters begin
                        import re

                        match = re.match(r"(\d+)([A-Za-z]+)", first_field)
                        if match:
                            res_num = int(match.group(1))
                            res_name = match.group(2)
                            atom_name = parts[1]
                            atom_num = int(parts[2])
                            x = float(parts[3])
                            y = float(parts[4])
                            z = float(parts[5])
                        else:
                            # Skip this line if we can't parse it
                            continue
                else:
                    # Assume combined format
                    first_field = parts[0]
                    import re

                    match = re.match(r"(\d+)([A-Za-z]+)", first_field)
                    if match:
                        res_num = int(match.group(1))
                        res_name = match.group(2)
                        atom_name = parts[1]
                        atom_num = int(parts[2])
                        x = float(parts[3])
                        y = float(parts[4])
                        z = float(parts[5]) if len(parts) > 5 else 0.0
                    else:
                        # Skip this line if we can't parse it
                        continue

                # Skip LP (lone pair) atoms
                if "LP" in atom_name.upper():
                    continue

                atom_count += 1
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)

                # Re-number atoms sequentially
                gro_line = f"{res_num:5d}{res_name:<5s}{atom_name:>5s}{atom_count:5d}"
                gro_line += f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
                gro_lines.append(gro_line)

        # Calculate box dimensions (with 1.0 nm padding)
        if x_coords:
            box_x = max(x_coords) - min(x_coords) + 1.0
            box_y = max(y_coords) - min(y_coords) + 1.0
            box_z = max(z_coords) - min(z_coords) + 1.0
        else:
            # Default box size if no atoms
            box_x = box_y = box_z = 1.0

        # Write GRO file
        gro_content = f"{self.moleculetype_name} system\n"
        gro_content += f"{atom_count:5d}\n"
        gro_content += "".join(gro_lines)
        gro_content += f"{box_x:10.7f}{box_y:10.7f}{box_z:10.7f}\n"

        return gro_content, atom_count

    def _parse_topology(self):
        """Parse topology file"""
        with open(self.top_file, "r") as f:
            content = f.read()

        # Parse sections
        self._parse_atoms_section(content)
        self._parse_bonds_section(content)
        self._parse_pairs_section(content)
        self._parse_angles_section(content)
        self._parse_dihedrals_section(content)

    def _parse_atoms_section(self, content):
        """Parse [atoms] section"""
        section = self._extract_section(content, "atoms")
        if not section:
            return

        for line in section.split("\n"):
            line = line.strip()
            if not line or line.startswith(";"):
                continue

            parts = line.split()
            if len(parts) >= 8:
                atom_id = parts[0]
                atom_type = parts[1]
                resnum = parts[2]
                resname = parts[3]
                atom_name = parts[4]
                cgnr = parts[5]
                charge = parts[6]
                mass = parts[7]

                # Store atom
                self.atoms.append(
                    {
                        "id": atom_id,
                        "type": atom_type,
                        "resnum": resnum,
                        "resname": self.moleculetype_name,  # Changed to LIG
                        "name": atom_name,
                        "cgnr": cgnr,
                        "charge": charge,
                        "mass": mass,
                    }
                )

                # Store atomtype if not exists
                if atom_type not in self.atomtypes:
                    self.atomtypes[atom_type] = {"name": atom_type, "mass": mass}

    def _read_bondtypes(self) -> Dict[Tuple[str, str, str], List[str]]:
        """Read bondtypes from bonded parameter files

        Looks for bondtypes in:
        1. cgenff_bonded_LIG.itp (PRISM-generated, ligand-specific)
        2. charmm36.ff/ffbonded.itp (CHARMM36 force field)

        Returns
        -------
        Dict[Tuple[str, str, str], List[str]]
            Mapping from (type1, type2, funct) -> [b0, kb]
            Includes both forward and reverse order for easy lookup
        """
        bondtypes = {}

        # Read from multiple sources in order of priority
        bonded_files = [
            self.cgenff_dir / "cgenff_bonded_LIG.itp",
            self.cgenff_dir / "charmm36.ff" / "ffbonded.itp",
        ]

        for bonded_file in bonded_files:
            if bonded_file.exists():
                print(f"  Reading bondtypes from: {bonded_file.name}")
                # Extract bondtypes section
                with open(bonded_file, "r") as f:
                    content = f.read()

                pattern = r"\[\s*bondtypes\s*\](.*?)(?=\[|$)"
                match = re.search(pattern, content, re.DOTALL)

                if match:
                    count = 0
                    for line in match.group(1).split("\n"):
                        line = line.strip()
                        if not line or line.startswith(";"):
                            continue
                        parts = line.split()
                        if len(parts) >= 5:
                            key = (parts[0], parts[1], parts[2])
                            # Only add if not already present (later sources don't override earlier ones)
                            if key not in bondtypes:
                                params = parts[3:5]
                                bondtypes[key] = params
                                # Add reverse order for matching (bonds are undirected)
                                bondtypes[(parts[1], parts[0], parts[2])] = params
                                count += 1

                    print(f"  Found {count} unique bondtypes in {bonded_file.name}")

        print(f"  Total bondtypes loaded: {len(bondtypes) // 2}")
        return bondtypes

    def _expand_bond_with_params(
        self, bond: List[str], atom_types: Dict[str, str], bondtypes: Dict[Tuple[str, str, str], List[str]]
    ) -> List[str]:
        """Expand bond entry with parameters from bondtypes

        Parameters
        ----------
        bond : List[str]
            Original bond entry [ai, aj, funct] or with params
        atom_types : Dict[str, str]
            Mapping from atom ID to atom type
        bondtypes : Dict[Tuple[str, str, str], List[str]]
            Bondtype parameters lookup

        Returns
        -------
        List[str]
            Bond entry with parameters [ai, aj, funct, b0, kb]
            Returns original if params not found or already present
        """
        # If already has parameters (5+ fields), return as-is
        if len(bond) >= 5:
            return bond

        ai, aj, funct = bond[0], bond[1], bond[2]

        # Get atom types for the two atoms
        type_a = atom_types.get(ai)
        type_b = atom_types.get(aj)

        if not type_a or not type_b:
            # Atom types not found, return original
            return bond

        # Look up bondtype parameters
        key = (type_a, type_b, funct)
        params = bondtypes.get(key)

        if params:
            return [ai, aj, funct] + params
        else:
            # Try reverse order (bonds are undirected)
            key_rev = (type_b, type_a, funct)
            params = bondtypes.get(key_rev)
            if params:
                return [ai, aj, funct] + params

        # No matching bondtype found, return original
        return bond

    def _parse_bonds_section(self, content):
        """Parse [bonds] section and expand with bondtype parameters"""
        section = self._extract_section(content, "bonds")
        if not section:
            return

        # Read bondtypes from bonded file
        bondtypes = self._read_bondtypes()

        # Build atom type mapping for parameter lookup
        atom_types = {atom["id"]: atom["type"] for atom in self.atoms}

        expanded_count = 0
        for line in section.split("\n"):
            line = line.strip()
            if not line or line.startswith(";"):
                continue

            parts = line.split()
            if len(parts) >= 3:
                # Expand with bondtype parameters if available
                original_len = len(parts)
                expanded = self._expand_bond_with_params(parts, atom_types, bondtypes)
                if len(expanded) > original_len:
                    expanded_count += 1
                self.bonds.append(expanded)

        if bondtypes and expanded_count > 0:
            print(f"  Expanded {expanded_count} bonds with parameters")

    def _parse_pairs_section(self, content):
        """Parse [pairs] section"""
        section = self._extract_section(content, "pairs")
        if not section:
            return

        for line in section.split("\n"):
            line = line.strip()
            if not line or line.startswith(";"):
                continue

            parts = line.split()
            if len(parts) >= 3:
                # Store all parameters
                self.pairs.append(parts)

    def _parse_angles_section(self, content):
        """Parse [angles] section"""
        section = self._extract_section(content, "angles")
        if not section:
            return

        for line in section.split("\n"):
            line = line.strip()
            if not line or line.startswith(";"):
                continue

            parts = line.split()
            if len(parts) >= 4:
                # Store all parameters
                self.angles.append(parts)

    def _parse_dihedrals_section(self, content):
        """Parse [dihedrals] section - both proper and improper"""
        # Find all dihedral sections
        pattern = r"\[\s*dihedrals\s*\](.*?)(?=\[|$)"
        matches = re.findall(pattern, content, re.DOTALL)

        for section in matches:
            for line in section.split("\n"):
                line = line.strip()
                if not line or line.startswith(";"):
                    continue

                parts = line.split()
                if len(parts) >= 5:
                    # Check function type (proper=9, improper=2)
                    funct = parts[4] if len(parts) > 4 else "9"
                    if funct == "2":
                        # Improper dihedral - store all parameters
                        self.dihedrals_improper.append(parts)
                    else:
                        # Proper dihedral - store all parameters
                        self.dihedrals_proper.append(parts)

    def _extract_section(self, content, section_name):
        """Extract a specific section from topology"""
        pattern = rf"\[\s*{section_name}\s*\](.*?)(?=\[|$)"
        match = re.search(pattern, content, re.DOTALL)
        return match.group(1) if match else None

    def _remove_lone_pairs(self):
        """
        Remove lone pair (LP) atoms from topology and transfer their charges to halogen atoms

        CGenFF generates LP atoms for halogens (F, Cl, Br, I) to represent lone pairs.
        These need to be removed for GROMACS compatibility:
        1. Identify LP atoms (names containing 'LP')
        2. Record LP charges and find bonded halogen atoms
        3. Remove LP atoms from all sections
        4. Transfer LP charges to halogen atoms
        """
        # Step 1: Identify LP atoms and build mapping
        lp_atoms = {}  # {lp_atom_id: {'charge': float, 'halogen_id': str}}
        halogen_elements = {"F", "Cl", "Br", "I", "CL"}

        for atom in self.atoms:
            atom_name = atom["name"].upper()
            # Check if atom name contains LP (e.g., LP1, LP2, LP)
            if "LP" in atom_name:
                lp_atoms[atom["id"]] = {"charge": float(atom["charge"]), "name": atom["name"], "halogen_id": None}

        if not lp_atoms:
            print("  No lone pair (LP) atoms found")
            return

        print(f"  Found {len(lp_atoms)} LP atoms to remove:")
        for lp_id, lp_data in lp_atoms.items():
            print(f"    - Atom {lp_id} ({lp_data['name']}): charge = {lp_data['charge']:.6f}")

        # Step 2: Find halogen atoms bonded to each LP
        # Look through bonds to find which halogen each LP is connected to
        for bond in self.bonds:
            atom1_id, atom2_id = bond[0], bond[1]

            # Check if one atom is LP and find its halogen partner
            if atom1_id in lp_atoms:
                # atom2 should be the halogen
                for atom in self.atoms:
                    if atom["id"] == atom2_id:
                        # Check if it's a halogen by name or type
                        atom_name_upper = atom["name"].upper()
                        is_halogen = any(hal in atom_name_upper for hal in halogen_elements)
                        if is_halogen or atom_name_upper[0] in halogen_elements:
                            lp_atoms[atom1_id]["halogen_id"] = atom2_id
                            print(f"    - LP {atom1_id} bonded to halogen {atom2_id} ({atom['name']})")
                        break
            elif atom2_id in lp_atoms:
                # atom1 should be the halogen
                for atom in self.atoms:
                    if atom["id"] == atom1_id:
                        atom_name_upper = atom["name"].upper()
                        is_halogen = any(hal in atom_name_upper for hal in halogen_elements)
                        if is_halogen or atom_name_upper[0] in halogen_elements:
                            lp_atoms[atom2_id]["halogen_id"] = atom1_id
                            print(f"    - LP {atom2_id} bonded to halogen {atom1_id} ({atom['name']})")
                        break

        # Step 3: Transfer charges from LP to halogen atoms
        charge_transfers = {}  # {halogen_id: total_lp_charge}
        for lp_id, lp_data in lp_atoms.items():
            if lp_data["halogen_id"]:
                halogen_id = lp_data["halogen_id"]
                if halogen_id not in charge_transfers:
                    charge_transfers[halogen_id] = 0.0
                charge_transfers[halogen_id] += lp_data["charge"]

        # Update halogen atom charges
        for atom in self.atoms:
            if atom["id"] in charge_transfers:
                old_charge = float(atom["charge"])
                new_charge = old_charge + charge_transfers[atom["id"]]
                atom["charge"] = f"{new_charge:.10f}"
                print(f"  Updated halogen {atom['id']} ({atom['name']}) charge: {old_charge:.6f} -> {new_charge:.6f}")

        # Step 4: Remove LP atoms from atoms list
        lp_atom_ids = set(lp_atoms.keys())
        self.atoms = [atom for atom in self.atoms if atom["id"] not in lp_atom_ids]
        print(f"  Removed {len(lp_atom_ids)} LP atoms from [atoms] section")

        # Step 5: Remove all interactions involving LP atoms
        # Helper function to check if any atom in a list is LP
        def contains_lp(atom_ids):
            return any(aid in lp_atom_ids for aid in atom_ids)

        # Remove bonds involving LP
        old_bonds = len(self.bonds)
        self.bonds = [bond for bond in self.bonds if not contains_lp(bond[:2])]
        print(f"  Removed {old_bonds - len(self.bonds)} bonds involving LP atoms")

        # Remove pairs involving LP
        old_pairs = len(self.pairs)
        self.pairs = [pair for pair in self.pairs if not contains_lp(pair[:2])]
        print(f"  Removed {old_pairs - len(self.pairs)} pairs involving LP atoms")

        # Remove angles involving LP
        old_angles = len(self.angles)
        self.angles = [angle for angle in self.angles if not contains_lp(angle[:3])]
        print(f"  Removed {old_angles - len(self.angles)} angles involving LP atoms")

        # Remove proper dihedrals involving LP
        old_proper = len(self.dihedrals_proper)
        self.dihedrals_proper = [dih for dih in self.dihedrals_proper if not contains_lp(dih[:4])]
        print(f"  Removed {old_proper - len(self.dihedrals_proper)} proper dihedrals involving LP atoms")

        # Remove improper dihedrals involving LP
        old_improper = len(self.dihedrals_improper)
        self.dihedrals_improper = [dih for dih in self.dihedrals_improper if not contains_lp(dih[:4])]
        print(f"  Removed {old_improper - len(self.dihedrals_improper)} improper dihedrals involving LP atoms")

        print(f"  Successfully removed all LP atoms and transferred charges to halogens")

    def _generate_atomtypes_itp(self):
        """Generate atomtypes_LIG.itp"""
        content = "[ atomtypes ]\n"
        content += ";type, bondingtype, atomic_number, mass, charge, ptype, sigma, epsilon\n"

        # Try to read from charmm36.ff/ffnonbonded.itp if available
        ffnonbonded = self.cgenff_dir / "charmm36.ff" / "ffnonbonded.itp"
        atomtype_params = {}

        if ffnonbonded.exists():
            with open(ffnonbonded, "r") as f:
                in_atomtypes = False
                for line in f:
                    if "[ atomtypes ]" in line:
                        in_atomtypes = True
                        continue
                    if in_atomtypes:
                        if line.strip().startswith("["):
                            break
                        if not line.strip() or line.strip().startswith(";"):
                            continue
                        parts = line.split()
                        if len(parts) >= 7:
                            atype = parts[0]
                            atomtype_params[atype] = line.strip()

        # Write atomtypes
        for atom_type in sorted(self.atomtypes.keys()):
            if atom_type in atomtype_params:
                content += atomtype_params[atom_type] + "\n"
            else:
                # Default values if not found
                mass = self.atomtypes[atom_type]["mass"]
                # Guess atomic number from mass
                atomic_num = self._guess_atomic_number(float(mass))
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
        elif mass < 20:
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

    def _generate_lig_itp(self):
        """Generate LIG.itp"""
        # NOTE: Do NOT include *_ffbonded.itp files here.
        # They must be included at the top-level topology (topol.top)
        # BEFORE any [ moleculetype ] sections, not inside LIG.itp.
        # This is because they contain [ bondtypes ], [ angletypes ],
        # [ dihedraltypes ] which must come before molecule definitions
        # in GROMACS topology files.
        # system.py handles this automatically when it detects charmm36.ff/.

        content = "[ moleculetype ]\n"
        content += f"{self.moleculetype_name}              3\n\n"

        # Atoms section
        content += "[ atoms ]\n"
        content += ";index, atom type, resnum, resname, name, cgnr, charge, mass\n"
        for atom in self.atoms:
            content += f"{atom['id']:6s} {atom['type']:14s} {atom['resnum']:4s} {atom['resname']:3s} "
            content += f"{atom['name']:5s} {atom['cgnr']:6s} {atom['charge']:>12s} {atom['mass']:>16s}\n"

        content += "\n"

        # Pairs section
        if self.pairs:
            content += "[ pairs ]\n"
            content += ";ai    aj   funct   parameters...\n"
            for pair in self.pairs:
                # Write all fields including any parameters
                content += "\t".join(f"{val:>6s}" if i < 3 else f"{val:>12s}" for i, val in enumerate(pair)) + "\n"
            content += "\n"

        # Bonds section
        if self.bonds:
            content += "[ bonds ]\n"
            content += ";ai    aj   funct   parameters...\n"
            for bond in self.bonds:
                # Write all fields including parameters
                content += "\t".join(f"{val:>6s}" if i < 3 else f"{val:>12s}" for i, val in enumerate(bond)) + "\n"
            content += "\n"

        # Angles section
        if self.angles:
            content += "[ angles ]\n"
            content += ";ai    aj    ak   funct   parameters...\n"
            for angle in self.angles:
                # Write all fields including parameters
                content += "\t".join(f"{val:>6s}" if i < 4 else f"{val:>12s}" for i, val in enumerate(angle)) + "\n"
            content += "\n"

        # Dihedrals section (proper)
        if self.dihedrals_proper:
            content += "[ dihedrals ]\n"
            content += ";ai    aj    ak    al   funct   parameters...\n"
            for dihedral in self.dihedrals_proper:
                # Write all fields including parameters
                content += "\t".join(f"{val:>6s}" if i < 5 else f"{val:>12s}" for i, val in enumerate(dihedral)) + "\n"
            content += "\n"

        # Dihedrals section (improper)
        if self.dihedrals_improper:
            content += "[ dihedrals ]\n"
            content += ";ai    aj    ak    al   funct   parameters...\n"
            for dihedral in self.dihedrals_improper:
                # Write all fields including parameters
                content += "\t".join(f"{val:>6s}" if i < 5 else f"{val:>12s}" for i, val in enumerate(dihedral)) + "\n"
            content += "\n"

        # Position restraints reference
        content += "#ifdef POSRES\n"
        content += '#include "posre_LIG.itp"\n'
        content += "#endif\n"

        return content

    def _generate_lig_top(self):
        """Generate LIG.top"""
        # Check if charmm36.ff directory exists
        charmm_ff_dir = self.cgenff_dir / "charmm36.ff"

        if charmm_ff_dir.exists():
            # Check for CHARMM-GUI format (complete charmm36.itp)
            charmm36_itp = charmm_ff_dir / "charmm36.itp"
            if charmm36_itp.exists():
                # CHARMM-GUI format: single file with all parameters
                content = "; Include CHARMM-GUI force field (complete)\n"
                content += f'#include "charmm36.ff/{charmm36_itp.name}"\n'
                content += "\n"
            else:
                # CGenFF website format: use forcefield.itp + ffbonded.itp
                content = "; Include CHARMM36 force field\n"
                content += '#include "charmm36.ff/forcefield.itp"\n'
                content += '#include "charmm36.ff/ffnonbonded.itp"\n'
                content += '#include "charmm36.ff/ffbonded.itp"\n'

                # Include ligand-specific bonded parameters if exists
                # Search for any *_ffbonded.itp files
                ffbonded_files = list(charmm_ff_dir.glob("*_ffbonded.itp"))
                for ffbonded_file in ffbonded_files:
                    # Check if file has content (not empty)
                    if ffbonded_file.stat().st_size > 200:  # More than just header
                        content += f'#include "charmm36.ff/{ffbonded_file.name}"\n'

                content += "\n"
        else:
            # Fallback: minimal topology
            content = "[ defaults ]\n"
            content += "; nbfunc\tcomb-rule\tgen-pairs\tfudgeLJ\tfudgeQQ\n"
            content += "     1\t     2\tyes   \t1.000000\t1.000000\n\n\n"
            content += '#include "atomtypes_LIG.itp"\n'

        content += '#include "LIG.itp"\n\n'

        content += "[ system ]\n"
        content += ";name\n"
        content += f"{self.moleculetype_name} system\n\n"

        content += "[ molecules ]\n"
        content += ";name\tnumber\n"
        content += f"{self.moleculetype_name}              1\n"

        return content

    def _generate_posre_itp(self, num_atoms):
        """Generate posre_LIG.itp"""
        content = "[ position_restraints ]\n"
        content += "; atom  type      fx      fy      fz\n"

        for i in range(1, num_atoms + 1):
            content += f"{i:6d}     1  1000  1000  1000\n"

        return content
