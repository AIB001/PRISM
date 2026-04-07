#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RTF (CHARMM Residue Topology File) force field generator for PRISM

This module processes CGenFF RTF+PRM files and converts them
to PRISM-standardized GROMACS format.
"""

import os
from pathlib import Path
from typing import Dict, Tuple

# Import the base class
try:
    from .base import ForceFieldGeneratorBase
except ImportError:
    from base import ForceFieldGeneratorBase


class RTFForceFieldGenerator(ForceFieldGeneratorBase):
    """RTF force field generator - processes CGenFF RTF+PRM files"""

    def __init__(self, rtf_file, prm_file, pdb_file, output_dir, overwrite=False):
        """
        Initialize RTF force field generator

        Parameters
        ----------
        rtf_file : str
            Path to RTF file (atom types, charges, bonds, angles)
        prm_file : str
            Path to PRM file (force field parameters)
        pdb_file : str
            Path to PDB file (coordinates)
        output_dir : str
            Output directory for processed files
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(pdb_file, output_dir, overwrite)  # Use pdb_file as ligand_path for API compatibility

        self.rtf_file = Path(rtf_file)
        self.prm_file = Path(prm_file)
        self.pdb_file = Path(pdb_file)

        if not self.rtf_file.exists():
            raise FileNotFoundError(f"RTF file does not exist: {self.rtf_file}")
        if not self.prm_file.exists():
            raise FileNotFoundError(f"PRM file does not exist: {self.prm_file}")
        if not self.pdb_file.exists():
            raise FileNotFoundError(f"PDB file does not exist: {self.pdb_file}")

        # Data storage
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.atomtypes = {}
        self.moleculetype_name = "LIG"

        # RTF data
        self.rtf_atoms = {}  # name -> {type, charge}
        self.rtf_masses = {}  # atom_type -> mass
        self.rtf_bonds = []  # [(atom1, atom2), ...]
        self.rtf_angles = []  # [(atom1, atom2, atom3), ...]
        self.rtf_dihedrals = []  # [(atom1, atom2, atom3, atom4), ...]
        self.auto_angles = False
        self.auto_dihedrals = False

        # PRM parameters (type-based)
        self.prm_bond_params = {}  # (type1, type2) -> (k_b, b0)
        self.prm_angle_params = {}  # (type1, type2, type3) -> (k_theta, theta0)
        self.prm_dihedral_params = {}  # (type1, type2, type3, type4) -> list of (k, phi0, multiplicity)

        print(f"\nInitialized RTF Force Field Generator:")
        print(f"  RTF file: {self.rtf_file}")
        print(f"  PRM file: {self.prm_file}")
        print(f"  PDB file: {self.pdb_file}")
        print(f"  Output directory: {self.output_dir}")

    def get_output_dir_name(self):
        """Get the output directory name for RTF"""
        return "LIG.rtf2gmx"

    def run(self):
        """Run the RTF force field conversion workflow"""
        print(f"\n{'='*60}")
        print("Starting RTF Force Field Processing")
        print(f"{'='*60}")

        try:
            # Check if output already exists
            lig_dir = os.path.join(self.output_dir, "LIG.rtf2gmx")
            if os.path.exists(lig_dir) and not self.overwrite:
                if self.check_required_files(lig_dir):
                    print(f"\nUsing cached RTF force field parameters from: {lig_dir}")
                    print("All required files found:")
                    for f in sorted(os.listdir(lig_dir)):
                        print(f"  - {f}")
                    print("\n(Use --overwrite to regenerate)")
                    return lig_dir

            # Parse RTF file
            print("\n[Step 1] Parsing RTF file...")
            self._parse_rtf()

            # Parse PRM file
            print("\n[Step 2] Parsing PRM file...")
            self._parse_prm()

            # Parse PDB file for coordinates
            print("\n[Step 3] Parsing PDB file...")
            coordinates = self._parse_pdb()

            # Generate GROMACS files
            print("\n[Step 4] Generating GROMACS files...")
            os.makedirs(lig_dir, exist_ok=True)

            # Generate GRO file
            self._generate_gro_file(lig_dir, coordinates)

            # Generate ITP file
            self._generate_itp_file(lig_dir)

            # Generate atomtypes file
            self._generate_atomtypes_file(lig_dir)

            # Generate position restraints
            self._generate_posre_file(lig_dir)

            # Generate topology file
            self._generate_top_file(lig_dir)

            print(f"\n{'='*60}")
            print("RTF force field processing completed successfully!")
            print(f"{'='*60}")

            print(f"Output files are in: {lig_dir}")
            print("\nGenerated files:")
            for f in sorted(os.listdir(lig_dir)):
                print(f"  - {f}")

            return lig_dir

        except Exception as e:
            print(f"\nError during RTF processing: {e}")
            raise

    def _parse_rtf(self):
        """Parse RTF file for atoms, bonds, angles, dihedrals"""
        print(f"  Reading RTF file: {self.rtf_file.name}")

        with open(self.rtf_file, "r") as f:
            lines = f.readlines()

        in_atoms = False
        in_bonds = False
        in_angles = False
        in_dihedrals = False
        in_impropers = False

        for line in lines:
            line = line.strip()
            if not line or line.startswith("*"):
                continue

            if line.startswith("MASS"):
                parts = line.split()
                if len(parts) >= 4:
                    atom_type = parts[2]
                    mass = float(parts[3])
                    self.rtf_masses[atom_type] = mass
                continue

            if line.upper().startswith("AUTO"):
                upper = line.upper()
                self.auto_angles = "ANGLE" in upper
                self.auto_dihedrals = "DIHE" in upper
                continue

            # Section detection
            if line.startswith("ATOM"):
                in_atoms = True
                in_bonds = in_angles = in_dihedrals = in_impropers = False
                # Parse atom line: ATOM name type charge
                parts = line.split()
                if len(parts) >= 4:
                    name = parts[1]
                    atom_type = parts[2]
                    charge = float(parts[3])
                    self.rtf_atoms[name] = {"type": atom_type, "charge": charge}
                continue

            elif line.startswith("BOND"):
                in_bonds = True
                in_atoms = in_angles = in_dihedrals = in_impropers = False
                # Parse bond line: BOND atom1 atom2 atom3 atom4 ...
                parts = line.split()[1:]  # Skip "BOND"
                for i in range(0, len(parts), 2):
                    if i + 1 < len(parts) and parts[i + 1] not in ["PATCH", "NONE", "IMPR", "FIRST", "LAST"]:
                        self.rtf_bonds.append((parts[i], parts[i + 1]))
                continue

            elif line.startswith("ANGL") or line.startswith("ANGLE"):
                in_angles = True
                in_atoms = in_bonds = in_dihedrals = in_impropers = False
                continue

            elif line.startswith("DIHE") or line.startswith("DIHEDRAL"):
                in_dihedrals = True
                in_atoms = in_bonds = in_angles = in_impropers = False
                continue

            elif line.startswith("IMPR") or line.startswith("IMPROPER"):
                in_impropers = True
                in_atoms = in_bonds = in_angles = in_dihedrals = False
                continue

            elif (
                line.startswith("PATCH")
                or line.startswith("NONE")
                or line.startswith("FIRST")
                or line.startswith("LAST")
            ):
                # Skip control statements
                in_atoms = in_bonds = in_angles = in_dihedrals = in_impropers = False
                continue

            # Continue parsing multi-line sections (but skip IMPROPERS)
            if in_bonds and line:
                parts = line.split()
                for i in range(0, len(parts), 2):
                    if i + 1 < len(parts) and parts[i + 1] not in ["PATCH", "NONE", "IMPR", "FIRST", "LAST"]:
                        self.rtf_bonds.append((parts[i], parts[i + 1]))

            elif in_angles and line:
                parts = line.split()
                for i in range(0, len(parts), 3):
                    if i + 2 < len(parts):
                        self.rtf_angles.append((parts[i], parts[i + 1], parts[i + 2]))

            elif in_dihedrals and line:
                parts = line.split()
                for i in range(0, len(parts), 4):
                    if i + 3 < len(parts):
                        self.rtf_dihedrals.append((parts[i], parts[i + 1], parts[i + 2], parts[i + 3]))

            # Skip IMPROPERS section for now (not used in basic ITP)

        print(f"  - Found {len(self.rtf_atoms)} atoms")
        print(f"  - Found {len(self.rtf_bonds)} bonds")
        if self.auto_angles and not self.rtf_angles:
            self.rtf_angles = self._autogenerate_angles()
        if self.auto_dihedrals and not self.rtf_dihedrals:
            self.rtf_dihedrals = self._autogenerate_dihedrals()
        print(f"  - Found {len(self.rtf_angles)} angles")
        print(f"  - Found {len(self.rtf_dihedrals)} dihedrals")

    def _build_bond_graph(self):
        """Return an undirected adjacency map from RTF bonds."""
        graph = {atom_name: set() for atom_name in self.rtf_atoms}
        for atom1, atom2 in self.rtf_bonds:
            if atom1 in graph and atom2 in graph:
                graph[atom1].add(atom2)
                graph[atom2].add(atom1)
        return graph

    def _autogenerate_angles(self):
        """Generate CHARMM AUTO ANGLES from the bond graph."""
        graph = self._build_bond_graph()
        angles = set()
        for center, neighbors in graph.items():
            ordered_neighbors = sorted(neighbors)
            for idx, atom1 in enumerate(ordered_neighbors):
                for atom3 in ordered_neighbors[idx + 1 :]:
                    angles.add((atom1, center, atom3))
        return sorted(angles)

    def _autogenerate_dihedrals(self):
        """Generate CHARMM AUTO DIHE torsions from the bond graph."""
        graph = self._build_bond_graph()
        dihedrals = set()
        for atom2, neighbors2 in graph.items():
            for atom3 in neighbors2:
                if atom2 >= atom3:
                    continue
                for atom1 in graph[atom2] - {atom3}:
                    for atom4 in graph[atom3] - {atom2}:
                        if atom1 == atom4:
                            continue
                        forward = (atom1, atom2, atom3, atom4)
                        reverse = tuple(reversed(forward))
                        dihedrals.add(min(forward, reverse))
        return sorted(dihedrals)

    def _parse_prm(self):
        """Parse PRM file for force field parameters"""
        print(f"  Reading PRM file: {self.prm_file.name}")

        with open(self.prm_file, "r") as f:
            lines = f.readlines()

        section = None
        for line in lines:
            line = line.strip()
            if not line or line.startswith("*"):
                continue

            # Detect sections
            if line.startswith("BONDS"):
                section = "bonds"
                continue
            elif line.startswith("ANGLES") or line.startswith("ANGL"):
                section = "angles"
                continue
            elif line.startswith("DIHEDRALS") or line.startswith("DIHE"):
                section = "dihedrals"
                continue
            elif line.startswith("NONBONDED") or line.startswith("IMPROPER"):
                section = None  # Skip these sections for now
                continue

            # Parse parameters based on section
            if section == "bonds":
                parts = line.split()
                if len(parts) >= 4:
                    type1, type2 = parts[0], parts[1]
                    k_b = float(parts[2])  # kcal/mol/Å^2
                    b0 = float(parts[3])  # Å
                    # Store both orderings for flexible lookup
                    self.prm_bond_params[(type1, type2)] = (k_b, b0)
                    self.prm_bond_params[(type2, type1)] = (k_b, b0)

            elif section == "angles":
                parts = line.split()
                if len(parts) >= 5:
                    type1, type2, type3 = parts[0], parts[1], parts[2]
                    k_theta = float(parts[3])  # kcal/mol/rad^2
                    theta0 = float(parts[4])  # degrees
                    self.prm_angle_params[(type1, type2, type3)] = (k_theta, theta0)
                    self.prm_angle_params[(type3, type2, type1)] = (k_theta, theta0)

            elif section == "dihedrals":
                parts = line.split()
                if len(parts) >= 7:
                    type1, type2, type3, type4 = parts[0], parts[1], parts[2], parts[3]
                    k = float(parts[4])  # kcal/mol
                    multiplicity = int(float(parts[5]))  # periodicity
                    phi0 = float(parts[6])  # degrees

                    key = (type1, type2, type3, type4)
                    if key not in self.prm_dihedral_params:
                        self.prm_dihedral_params[key] = []
                    self.prm_dihedral_params[key].append((k, phi0, multiplicity))

                    reverse_key = (type4, type3, type2, type1)
                    if reverse_key not in self.prm_dihedral_params:
                        self.prm_dihedral_params[reverse_key] = []
                    self.prm_dihedral_params[reverse_key].append((k, phi0, multiplicity))

        print(f"  - Found {len(self.prm_bond_params)//2} bond parameters")
        print(f"  - Found {len(self.prm_angle_params)} angle parameters")
        print(f"  - Found {len(self.prm_dihedral_params)} dihedral parameter sets")

    def _parse_pdb(self) -> Dict[str, Tuple[float, float, float]]:
        """Parse PDB file for coordinates"""
        print(f"  Reading PDB file: {self.pdb_file.name}")

        coordinates = {}
        with open(self.pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates[atom_name] = (x, y, z)

        print(f"  - Found coordinates for {len(coordinates)} atoms")
        return coordinates

    def _generate_gro_file(self, output_dir: str, coordinates: Dict[str, Tuple[float, float, float]]):
        """Generate GROMACS GRO file"""
        gro_file = os.path.join(output_dir, "LIG.gro")
        print(f"  Generating GRO file: {os.path.basename(gro_file)}")

        with open(gro_file, "w") as f:
            f.write("Generated by RTFForceFieldGenerator\n")
            f.write(f"{len(self.rtf_atoms):5d}\n")

            atom_id = 0
            for atom_name in self.rtf_atoms:
                if atom_name in coordinates:
                    atom_id += 1
                    x, y, z = coordinates[atom_name]
                    # Convert Å to nm
                    x_nm = x / 10.0
                    y_nm = y / 10.0
                    z_nm = z / 10.0

                    f.write(f"{1:5d}{self.moleculetype_name:<5s}{atom_name:>5s}{atom_id:5d}")
                    f.write(f"{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n")

            # Box dimensions (dummy)
            f.write("   5.00000   5.00000   5.00000\n")

    def _generate_itp_file(self, output_dir: str):
        """Generate GROMACS ITP file"""
        itp_file = os.path.join(output_dir, "LIG.itp")
        print(f"  Generating ITP file: {os.path.basename(itp_file)}")

        with open(itp_file, "w") as f:
            f.write("; Generated by RTFForceFieldGenerator\n")
            f.write(f"; RTF file: {self.rtf_file.name}\n")
            f.write(f"; PRM file: {self.prm_file.name}\n")
            f.write("\n")

            # Moleculetype
            f.write("[ moleculetype ]\n")
            f.write(f"{self.moleculetype_name}    3\n\n")

            # Atoms
            f.write("[ atoms ]\n")
            f.write(";   nr       type  resnr residue  atom   cgnr     charge       mass\n")

            atom_id = 0
            for atom_name, atom_data in self.rtf_atoms.items():
                atom_id += 1
                atom_type = atom_data["type"]
                charge = atom_data["charge"]
                mass = self.rtf_masses.get(atom_type, 12.01000)
                f.write(f"{atom_id:5d} {atom_type:>10s}      1    {self.moleculetype_name}   {atom_name:>4s}")
                f.write(f"{atom_id:5d}  {charge:10.6f}   {mass:8.5f}\n")

            # Bonds
            if self.rtf_bonds:
                f.write("\n[ bonds ]\n")
                f.write(";  ai    aj funct            c0            c1            c2            c3\n")

                # Create atom name to ID and type mappings
                name_to_id = {name: i + 1 for i, name in enumerate(self.rtf_atoms.keys())}
                name_to_type = {name: data["type"] for name, data in self.rtf_atoms.items()}

                for atom1, atom2 in self.rtf_bonds:
                    if atom1 in name_to_id and atom2 in name_to_id:
                        id1 = name_to_id[atom1]
                        id2 = name_to_id[atom2]
                        type1 = name_to_type[atom1]
                        type2 = name_to_type[atom2]

                        # Look up bond parameters
                        param_key = (type1, type2)
                        if param_key in self.prm_bond_params:
                            k_b, b0 = self.prm_bond_params[param_key]
                            # CHARMM harmonic bonds omit the 1/2 prefactor used by
                            # GROMACS funct=1, so the force constant doubles here.
                            k_gromacs = k_b * 836.8  # kcal/mol/A^2 -> kJ/mol/nm^2
                            b0_gromacs = b0 / 10.0  # nm
                            f.write(f"{id1:5d} {id2:5d}     1   {b0_gromacs:12.6f}  {k_gromacs:12.6f}\n")
                        else:
                            # Write placeholder if parameters not found
                            f.write(f"{id1:5d} {id2:5d}     1\n")

            # Angles
            if self.rtf_angles:
                f.write("\n[ angles ]\n")
                f.write(";  ai    aj    ak funct            c0            c1            c2            c3\n")

                name_to_id = {name: i + 1 for i, name in enumerate(self.rtf_atoms.keys())}
                name_to_type = {name: data["type"] for name, data in self.rtf_atoms.items()}

                for atom1, atom2, atom3 in self.rtf_angles:
                    if all(atom in name_to_id for atom in [atom1, atom2, atom3]):
                        id1 = name_to_id[atom1]
                        id2 = name_to_id[atom2]
                        id3 = name_to_id[atom3]
                        type1 = name_to_type[atom1]
                        type2 = name_to_type[atom2]
                        type3 = name_to_type[atom3]

                        # Look up angle parameters
                        param_key = (type1, type2, type3)
                        if param_key in self.prm_angle_params:
                            k_theta, theta0 = self.prm_angle_params[param_key]
                            # CHARMM harmonic angles omit the 1/2 prefactor used by
                            # GROMACS funct=1, so the force constant doubles here.
                            k_gromacs = k_theta * 8.368  # kcal/mol/rad^2 -> kJ/mol/rad^2
                            f.write(f"{id1:5d} {id2:5d} {id3:5d}     1   {theta0:12.6f}  {k_gromacs:12.6f}\n")
                        else:
                            # Write placeholder if parameters not found
                            f.write(f"{id1:5d} {id2:5d} {id3:5d}     1\n")

            # Dihedrals
            if self.rtf_dihedrals:
                f.write("\n[ dihedrals ]\n")
                f.write(
                    ";  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5\n"
                )

                name_to_id = {name: i + 1 for i, name in enumerate(self.rtf_atoms.keys())}
                name_to_type = {name: data["type"] for name, data in self.rtf_atoms.items()}

                for atom1, atom2, atom3, atom4 in self.rtf_dihedrals:
                    if all(atom in name_to_id for atom in [atom1, atom2, atom3, atom4]):
                        id1 = name_to_id[atom1]
                        id2 = name_to_id[atom2]
                        id3 = name_to_id[atom3]
                        id4 = name_to_id[atom4]
                        type1 = name_to_type[atom1]
                        type2 = name_to_type[atom2]
                        type3 = name_to_type[atom3]
                        type4 = name_to_type[atom4]

                        # Look up dihedral parameters
                        param_key = (type1, type2, type3, type4)
                        if param_key in self.prm_dihedral_params:
                            # Write each term in the Fourier series
                            for k, phi0, multiplicity in self.prm_dihedral_params[param_key]:
                                k_gromacs = k * 4.184  # kJ/mol
                                f.write(
                                    f"{id1:5d} {id2:5d} {id3:5d} {id4:5d}     9   {phi0:12.6f}  {k_gromacs:12.6f}  {multiplicity:5d}\n"
                                )
                        else:
                            # Write placeholder if parameters not found
                            f.write(f"{id1:5d} {id2:5d} {id3:5d} {id4:5d}     9\n")

    def _generate_atomtypes_file(self, output_dir: str):
        """Generate atomtypes file"""
        atomtypes_file = os.path.join(output_dir, "atomtypes_LIG.itp")
        print(f"  Generating atomtypes file: {os.path.basename(atomtypes_file)}")

        unique_types = set(atom_data["type"] for atom_data in self.rtf_atoms.values())

        with open(atomtypes_file, "w") as f:
            f.write("; Generated by RTFForceFieldGenerator\n")
            f.write("[ atomtypes ]\n")
            f.write(";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb\n")

            for atom_type in sorted(unique_types):
                # Default parameters (should be parsed from PRM file in full implementation)
                f.write(f"{atom_type:>6s}  {atom_type:>6s}  12.01000  0.00000  A   3.39967e-01  4.57730e-01\n")

    def _generate_posre_file(self, output_dir: str):
        """Generate position restraints file"""
        posre_file = os.path.join(output_dir, "posre_LIG.itp")
        print(f"  Generating position restraints: {os.path.basename(posre_file)}")

        with open(posre_file, "w") as f:
            f.write("; Generated by RTFForceFieldGenerator\n")
            f.write("[ position_restraints ]\n")
            f.write(";  i funct       fcx        fcy        fcz\n")

            atom_id = 0
            for atom_name in self.rtf_atoms:
                atom_id += 1
                f.write(f"{atom_id:4d}    1       1000       1000       1000\n")

    def _generate_top_file(self, output_dir: str):
        """Generate topology file"""
        top_file = os.path.join(output_dir, "LIG.top")
        print(f"  Generating topology file: {os.path.basename(top_file)}")

        with open(top_file, "w") as f:
            f.write("; Generated by RTFForceFieldGenerator\n")
            f.write(f"; RTF file: {self.rtf_file.name}\n")
            f.write(f"; PRM file: {self.prm_file.name}\n")
            f.write("\n")

            f.write('#include "atomtypes_LIG.itp"\n')
            f.write('#include "LIG.itp"\n')
            f.write("\n")

            f.write("[ system ]\n")
            f.write("LIG\n")
            f.write("\n")

            f.write("[ molecules ]\n")
            f.write("LIG    1\n")

    def check_required_files(self, lig_dir: str) -> bool:
        """Check if all required files exist"""
        required_files = ["LIG.itp", "LIG.gro", "atomtypes_LIG.itp", "posre_LIG.itp", "LIG.top"]
        return all(os.path.exists(os.path.join(lig_dir, f)) for f in required_files)


if __name__ == "__main__":
    # Example usage
    generator = RTFForceFieldGenerator(
        rtf_file="ligand.rtf", prm_file="ligand.prm", pdb_file="ligand.pdb", output_dir="output"
    )
    result = generator.run()
    print(f"Generated PRISM format files in: {result}")
