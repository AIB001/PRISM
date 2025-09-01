# # #!/usr/bin/env python3
# # # -*- coding: utf-8 -*-

# # """
# # OpenFF force field generator wrapper for PRISM
# # """

# # import os
# # import sys
# # import tempfile
# # import shutil
# # from pathlib import Path

# # # Import the base class
# # try:
# #     from .base import ForceFieldGeneratorBase
# # except ImportError:
# #     from base import ForceFieldGeneratorBase

# # # Import OpenFF specific dependencies
# # try:
# #     from openff.toolkit import Molecule, ForceField
# #     from openff.interchange import Interchange
# #     from openff.units import unit
# #     import numpy as np
# # except ImportError as e:
# #     print(f"Error: Missing OpenFF dependencies: {e}")
# #     print("Please install: pip install openff-toolkit openff-interchange")
# #     sys.exit(1)


# # class OpenFFForceFieldGenerator(ForceFieldGeneratorBase):
# #     """OpenFF force field generator wrapper"""
    
# #     def __init__(self, ligand_path, output_dir, charge=0, forcefield="openff-2.1.0", overwrite=False):
# #         """
# #         Initialize OpenFF force field generator
        
# #         Parameters:
# #         -----------
# #         ligand_path : str
# #             Path to the ligand file (SDF/MOL2)
# #         output_dir : str
# #             Directory where output files will be stored
# #         charge : int
# #             Molecule charge (default: 0)
# #         forcefield : str
# #             OpenFF force field version (default: "openff-2.1.0")
# #         overwrite : bool
# #             Whether to overwrite existing files
# #         """
# #         super().__init__(ligand_path, output_dir, overwrite)
        
# #         self.charge = charge
# #         self.forcefield_name = forcefield
# #         self.box_size = 1.0  # nm
        
# #         # Output directory for OpenFF files
# #         self.lig_ff_dir = os.path.join(self.output_dir, self.get_output_dir_name())
        
# #         print(f"\nInitialized OpenFF Force Field Generator:")
# #         print(f"  Ligand: {self.ligand_path}")
# #         print(f"  Charge: {self.charge}")
# #         print(f"  Force field: {self.forcefield_name}")
# #         print(f"  Output directory: {self.output_dir}")
    
# #     def get_output_dir_name(self):
# #         """Get the output directory name for OpenFF"""
# #         return "LIG.openff2gmx"
    
# #     def run(self):
# #         """Run the OpenFF force field generation workflow"""
# #         print(f"\n{'='*60}")
# #         print("Starting OpenFF Force Field Generation")
# #         print(f"{'='*60}")
        
# #         try:
# #             # Check if output already exists
# #             if os.path.exists(self.lig_ff_dir) and not self.overwrite:
# #                 if self.check_required_files(self.lig_ff_dir):
# #                     print(f"Output directory {self.lig_ff_dir} already exists with all required files.")
# #                     print("Skipping force field generation (use --overwrite to regenerate)")
# #                     return self.lig_ff_dir
            
# #             # Create output directory
# #             os.makedirs(self.lig_ff_dir, exist_ok=True)
            
# #             # Load molecule
# #             mol = self._load_molecule()
            
# #             # Generate force field parameters
# #             self._generate_parameters(mol)
            
# #             print(f"\n{'='*60}")
# #             print("Force field generation completed successfully!")
# #             print(f"{'='*60}")
# #             print(f"\nOutput files are in: {self.lig_ff_dir}")
# #             print("\nGenerated files:")
# #             for f in os.listdir(self.lig_ff_dir):
# #                 print(f"  - {f}")
            
# #             return self.lig_ff_dir
            
# #         except Exception as e:
# #             print(f"\nError during force field generation: {e}")
# #             import traceback
# #             traceback.print_exc()
# #             raise
    
# #     def _load_molecule(self):
# #         """Load molecule from file"""
# #         print("Loading molecule...")
        
# #         file_ext = Path(self.ligand_path).suffix.lower()
        
# #         if file_ext == '.sdf':
# #             return self._load_sdf()
# #         elif file_ext == '.mol2':
# #             return self._load_mol2()
# #         else:
# #             raise ValueError(f"Unsupported file format: {file_ext}")
    
# #     def _load_sdf(self):
# #         """Load molecule from SDF file"""
# #         try:
# #             mol = Molecule.from_file(self.ligand_path)
# #             print(f"Loaded molecule from SDF: {mol.n_atoms} atoms")
# #             return mol
# #         except Exception as e:
# #             print(f"Failed to load SDF directly: {e}")
# #             print("Trying alternative loading method...")
            
# #             try:
# #                 from rdkit import Chem
# #                 supplier = Chem.SDMolSupplier(self.ligand_path)
# #                 rdmol = supplier[0] if supplier else None
                
# #                 if not rdmol:
# #                     raise ValueError("No valid molecule found in SDF")
                
# #                 mol = Molecule.from_rdkit(rdmol)
# #                 print(f"Loaded molecule via RDKit: {mol.n_atoms} atoms")
# #                 return mol
# #             except ImportError:
# #                 raise RuntimeError("RDKit is required for this SDF file. Please install it.")
    
# #     def _load_mol2(self):
# #         """Load molecule from MOL2 file"""
# #         print("Converting MOL2 to SDF for OpenFF...")
        
# #         with tempfile.TemporaryDirectory() as temp_dir:
# #             temp_sdf = os.path.join(temp_dir, "temp.sdf")
            
# #             # Try to use antechamber for conversion
# #             try:
# #                 import subprocess
# #                 cmd = [
# #                     'antechamber',
# #                     '-i', self.ligand_path,
# #                     '-fi', 'mol2',
# #                     '-o', temp_sdf,
# #                     '-fo', 'sdf',
# #                     '-c', 'wc',  # Write charges
# #                     '-nc', str(self.charge)
# #                 ]
                
# #                 result = subprocess.run(cmd, capture_output=True, text=True)
                
# #                 if result.returncode == 0 and os.path.exists(temp_sdf):
# #                     mol = Molecule.from_file(temp_sdf)
# #                     print(f"Converted MOL2 to SDF: {mol.n_atoms} atoms")
# #                     return mol
# #                 else:
# #                     raise RuntimeError(f"Antechamber conversion failed: {result.stderr}")
                    
# #             except Exception as e:
# #                 print(f"Antechamber conversion failed: {e}")
                
# #                 # Try Open Babel as fallback
# #                 try:
# #                     import subprocess
# #                     cmd = [
# #                         'obabel',
# #                         '-i', 'mol2', self.ligand_path,
# #                         '-o', 'sdf', '-O', temp_sdf,
# #                         '-h'  # Add hydrogens
# #                     ]
                    
# #                     result = subprocess.run(cmd, capture_output=True, text=True)
                    
# #                     if result.returncode == 0 and os.path.exists(temp_sdf):
# #                         mol = Molecule.from_file(temp_sdf)
# #                         print(f"Converted MOL2 to SDF using Open Babel: {mol.n_atoms} atoms")
# #                         return mol
# #                     else:
# #                         raise RuntimeError("Open Babel conversion failed")
                        
# #                 except Exception:
# #                     raise RuntimeError("Could not convert MOL2 file. Please install antechamber or obabel.")
    
# #     def _generate_parameters(self, mol):
# #         """Generate OpenFF parameters and write GROMACS files"""
# #         print("Generating OpenFF parameters...")
        
# #         # Load force field
# #         ff = ForceField(f"{self.forcefield_name}.offxml")
        
# #         # Create interchange
# #         print("Creating OpenFF Interchange...")
# #         interchange = Interchange.from_smirnoff(ff, [mol])
        
# #         # Add box vectors
# #         interchange = self._add_box_vectors(interchange)
        
# #         # Write GROMACS files
# #         print("Writing GROMACS files...")
# #         self._write_gromacs_files(interchange)
    
# #     def _add_box_vectors(self, interchange):
# #         """Add box vectors to the system"""
# #         coords = interchange.positions.m_as(unit.angstrom)
# #         min_c, max_c = np.min(coords, axis=0), np.max(coords, axis=0)
# #         box_nm = (max_c - min_c + 2 * self.box_size * 10) * 0.1
# #         interchange.box = np.diag(box_nm) * unit.nanometer
# #         return interchange
    
# #     def _write_gromacs_files(self, interchange):
# #         """Write GROMACS format files"""
# #         output_prefix = os.path.join(self.lig_ff_dir, "LIG")
        
# #         # Write GRO file
# #         interchange.to_gro(f"{output_prefix}.gro")
        
# #         # Write TOP file
# #         interchange.to_top(f"{output_prefix}.top")
        
# #         # Process the topology files
# #         self._process_topology_files(output_prefix)
    
# #     def _process_topology_files(self, output_prefix):
# #         """Process and split topology files"""
# #         top_file = f"{output_prefix}.top"
        
# #         with open(top_file, 'r') as f:
# #             top_content = f.read()
        
# #         # Extract different sections
# #         atomtypes_content = self._extract_atomtypes_content(top_content)
# #         itp_content = self._extract_itp_content(top_content)
        
# #         # Modify ITP content to use LIG as resname
# #         modified_itp_content = self._modify_resname_to_lig(itp_content)
        
# #         # Count atoms for position restraints
# #         atom_count = self._count_atoms_in_itp(modified_itp_content)
        
# #         # Add position restraints include at the end of ITP
# #         while modified_itp_content.endswith("\n\n"):
# #             modified_itp_content = modified_itp_content[:-1]
# #         modified_itp_content += "\n\n#ifdef POSRES\n#include \"posre_LIG.itp\"\n#endif\n"
        
# #         # Write ITP file
# #         itp_file = f"{output_prefix}.itp"
# #         with open(itp_file, "w") as f:
# #             f.write(modified_itp_content)
        
# #         # Write atomtypes file
# #         atomtypes_file = os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")
# #         with open(atomtypes_file, "w") as f:
# #             f.write(atomtypes_content)
        
# #         # Write position restraints file
# #         posre_file = os.path.join(self.lig_ff_dir, "posre_LIG.itp")
# #         self._write_position_restraints(posre_file, atom_count)
        
# #         # Rebuild top file
# #         self._rebuild_top_file(top_content, itp_file, output_prefix)
        
# #         # Modify GRO file to use LIG residue name
# #         self._modify_gro_file(output_prefix)
    
# #     def _extract_atomtypes_content(self, top_content):
# #         """Extract atomtypes section from topology"""
# #         lines = top_content.split("\n")
# #         atomtypes_lines = []
# #         in_section = False
        
# #         for line in lines:
# #             stripped = line.strip()
# #             if stripped.startswith("[ atomtypes ]"):
# #                 in_section = True
# #                 atomtypes_lines.append(line)
# #             elif stripped.startswith("[") and in_section:
# #                 in_section = False
# #             elif in_section:
# #                 atomtypes_lines.append(line)
        
# #         return "\n".join(atomtypes_lines)
    
# #     def _extract_itp_content(self, top_content):
# #         """Extract moleculetype section from topology"""
# #         lines = top_content.split("\n")
# #         itp_lines = []
# #         in_section = False
        
# #         for line in lines:
# #             stripped = line.strip()
# #             if stripped.startswith("[ moleculetype ]"):
# #                 in_section = True
# #             elif stripped.startswith(("[ system ]", "[ molecules ]")):
# #                 in_section = False
            
# #             if in_section:
# #                 itp_lines.append(line)
        
# #         return "\n".join(itp_lines)
    
# #     def _modify_resname_to_lig(self, itp_content):
# #         """Modify residue names to LIG"""
# #         lines = itp_content.split("\n")
# #         modified_lines = []
# #         in_moleculetype = False
# #         in_atoms = False
        
# #         for line in lines:
# #             stripped = line.strip()
            
# #             # Handle moleculetype section
# #             if stripped.startswith("[ moleculetype ]"):
# #                 in_moleculetype = True
# #                 modified_lines.append(line)
# #             elif in_moleculetype and not stripped.startswith(";") and stripped:
# #                 # Replace molecule name with LIG
# #                 parts = line.split()
# #                 if len(parts) >= 2:
# #                     line = f"LIG              {parts[1]}"
# #                 in_moleculetype = False
# #                 modified_lines.append(line)
            
# #             # Handle atoms section
# #             elif stripped.startswith("[ atoms ]"):
# #                 in_atoms = True
# #                 modified_lines.append(line)
# #             elif in_atoms and stripped and not stripped.startswith(";"):
# #                 # Replace resname with LIG in atoms section
# #                 parts = line.split()
# #                 if len(parts) >= 8:
# #                     # Format: index, atom_type, resnum, resname, name, cgnr, charge, mass
# #                     parts[3] = "LIG"  # resname is the 4th field
# #                     # Reconstruct line with proper formatting
# #                     line = f"{parts[0]:>6} {parts[1]:<15} {parts[2]:>4} {parts[3]:<7} {parts[4]:<10} {parts[5]:>4} {parts[6]:>17} {parts[7]:>17}"
# #                 modified_lines.append(line)
# #             elif stripped.startswith("[") and in_atoms:
# #                 in_atoms = False
# #                 modified_lines.append(line)
# #             else:
# #                 modified_lines.append(line)
        
# #         return "\n".join(modified_lines)
    
# #     def _count_atoms_in_itp(self, itp_content):
# #         """Count atoms in ITP file"""
# #         lines = itp_content.split("\n")
# #         in_atoms = False
# #         atom_count = 0
        
# #         for line in lines:
# #             stripped = line.strip()
# #             if stripped.startswith("[ atoms ]"):
# #                 in_atoms = True
# #                 continue
# #             elif stripped.startswith("[") and in_atoms:
# #                 break
# #             elif in_atoms and stripped and not stripped.startswith(";"):
# #                 parts = stripped.split()
# #                 if len(parts) >= 8:
# #                     try:
# #                         atom_index = int(parts[0])
# #                         atom_count = max(atom_count, atom_index)
# #                     except (ValueError, IndexError):
# #                         pass
        
# #         return atom_count
    
# #     def _write_position_restraints(self, posre_file, atom_count):
# #         """Write position restraints file"""
# #         content = "[ position_restraints ]\n"
# #         content += "; atom  type      fx      fy      fz\n"
        
# #         for i in range(1, atom_count + 1):
# #             content += f"{i:>6}     1  1000  1000  1000\n"
        
# #         with open(posre_file, "w") as f:
# #             f.write(content)
    
# #     def _rebuild_top_file(self, top_content, itp_file, output_prefix):
# #         """Rebuild topology file with proper includes"""
# #         itp_name = os.path.basename(itp_file)
# #         lines = top_content.split("\n")
# #         new_lines = []
# #         skip_section = False
# #         skip_atomtypes = False
# #         added_includes = False
        
# #         for line in lines:
# #             stripped = line.strip()
            
# #             # Skip atomtypes section
# #             if stripped.startswith("[ atomtypes ]"):
# #                 skip_atomtypes = True
# #                 if not added_includes:
# #                     # Add includes after defaults section
# #                     new_lines.append("")
# #                     new_lines.append("#include \"atomtypes_LIG.itp\"")
# #                     new_lines.append(f"#include \"{itp_name}\"")
# #                     new_lines.append("")
# #                     added_includes = True
# #                 continue
# #             elif skip_atomtypes and stripped.startswith("[") and not stripped.startswith("[ atomtypes ]"):
# #                 skip_atomtypes = False
            
# #             # Skip moleculetype and related sections
# #             if stripped.startswith("[ moleculetype ]"):
# #                 skip_section = True
# #                 continue
# #             elif stripped.startswith(("[ system ]", "[ molecules ]")):
# #                 skip_section = False
# #                 new_lines.append(line)
# #             elif not skip_section and not skip_atomtypes:
# #                 new_lines.append(line)
        
# #         # Fix molecules section to use LIG
# #         final_lines = []
# #         in_molecules = False
# #         in_system = False
# #         for line in new_lines:
# #             if line.strip().startswith("[ system ]"):
# #                 in_system = True
# #                 final_lines.append(line)
# #             elif in_system and line.strip() and not line.strip().startswith(";"):
# #                 # Replace system name
# #                 final_lines.append("LIG system")
# #                 in_system = False
# #             elif line.strip().startswith("[ molecules ]"):
# #                 in_molecules = True
# #                 final_lines.append(line)
# #             elif in_molecules and line.strip() and not line.strip().startswith(";"):
# #                 # Replace molecule name with LIG
# #                 parts = line.split()
# #                 if len(parts) >= 2:
# #                     line = f"LIG              {parts[1]}"
# #                 final_lines.append(line)
# #                 in_molecules = False
# #             else:
# #                 final_lines.append(line)
        
# #         with open(f"{output_prefix}.top", "w") as f:
# #             f.write("\n".join(final_lines))
    
# #     def _modify_gro_file(self, output_prefix):
# #         """Modify GRO file to use LIG residue name"""
# #         gro_file = f"{output_prefix}.gro"
        
# #         with open(gro_file, 'r') as f:
# #             lines = f.readlines()
        
# #         # Modify residue names in GRO file
# #         modified_lines = []
# #         for i, line in enumerate(lines):
# #             if i == 0:
# #                 # First line is the title
# #                 modified_lines.append("LIG system\n")
# #             elif i >= 2 and i < len(lines) - 1:  # Skip header and box vectors
# #                 # GRO format: residue name is columns 6-10 (0-indexed: 5-10)
# #                 if len(line) > 10:
# #                     line = line[:5] + "LIG".ljust(5) + line[10:]
# #                 modified_lines.append(line)
# #             else:
# #                 modified_lines.append(line)
        
# #         with open(gro_file, 'w') as f:
# #             f.writelines(modified_lines)

# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# """
# OpenFF force field generator wrapper for PRISM
# """

# import os
# import sys
# import tempfile
# import shutil
# from pathlib import Path

# # Import the base class
# try:
#     from .base import ForceFieldGeneratorBase
# except ImportError:
#     from base import ForceFieldGeneratorBase

# # Import OpenFF specific dependencies
# try:
#     from openff.toolkit import Molecule, ForceField
#     from openff.interchange import Interchange
#     from openff.units import unit
#     import numpy as np
# except ImportError as e:
#     print(f"Error: Missing OpenFF dependencies: {e}")
#     print("Please install: pip install openff-toolkit openff-interchange")
#     sys.exit(1)


# class OpenFFForceFieldGenerator(ForceFieldGeneratorBase):
#     """OpenFF force field generator wrapper"""
    
#     def __init__(self, ligand_path, output_dir, charge=0, forcefield="openff-2.1.0", overwrite=False):
#         """
#         Initialize OpenFF force field generator
        
#         Parameters:
#         -----------
#         ligand_path : str
#             Path to the ligand file (SDF/MOL2)
#         output_dir : str
#             Directory where output files will be stored
#         charge : int
#             Molecule charge (default: 0)
#         forcefield : str
#             OpenFF force field version (default: "openff-2.1.0")
#         overwrite : bool
#             Whether to overwrite existing files
#         """
#         super().__init__(ligand_path, output_dir, overwrite)
        
#         self.charge = charge
#         self.forcefield_name = forcefield
#         self.box_size = 1.0  # nm
        
#         # Output directory for OpenFF files
#         self.lig_ff_dir = os.path.join(self.output_dir, self.get_output_dir_name())
        
#         print(f"\nInitialized OpenFF Force Field Generator:")
#         print(f"  Ligand: {self.ligand_path}")
#         print(f"  Charge: {self.charge}")
#         print(f"  Force field: {self.forcefield_name}")
#         print(f"  Output directory: {self.output_dir}")
    
#     def get_output_dir_name(self):
#         """Get the output directory name for OpenFF"""
#         return "LIG.openff2gmx"
    
#     def run(self):
#         """Run the OpenFF force field generation workflow"""
#         print(f"\n{'='*60}")
#         print("Starting OpenFF Force Field Generation")
#         print(f"{'='*60}")
        
#         try:
#             # Check if output already exists
#             if os.path.exists(self.lig_ff_dir) and not self.overwrite:
#                 if self.check_required_files(self.lig_ff_dir):
#                     print(f"Output directory {self.lig_ff_dir} already exists with all required files.")
#                     print("Skipping force field generation (use --overwrite to regenerate)")
#                     return self.lig_ff_dir
            
#             # Create output directory
#             os.makedirs(self.lig_ff_dir, exist_ok=True)
            
#             # Load molecule
#             mol = self._load_molecule()
            
#             # Generate force field parameters
#             self._generate_parameters(mol)
            
#             print(f"\n{'='*60}")
#             print("Force field generation completed successfully!")
#             print(f"{'='*60}")
#             print(f"\nOutput files are in: {self.lig_ff_dir}")
#             print("\nGenerated files:")
#             for f in os.listdir(self.lig_ff_dir):
#                 print(f"  - {f}")
            
#             return self.lig_ff_dir
            
#         except Exception as e:
#             print(f"\nError during force field generation: {e}")
#             import traceback
#             traceback.print_exc()
#             raise
    
#     def _load_molecule(self):
#         """Load molecule from file"""
#         print("Loading molecule...")
        
#         file_ext = Path(self.ligand_path).suffix.lower()
        
#         if file_ext == '.sdf':
#             return self._load_sdf()
#         elif file_ext == '.mol2':
#             return self._load_mol2()
#         else:
#             raise ValueError(f"Unsupported file format: {file_ext}")
    
#     def _load_sdf(self):
#         """Load molecule from SDF file"""
#         try:
#             mol = Molecule.from_file(self.ligand_path)
#             print(f"Loaded molecule from SDF: {mol.n_atoms} atoms")
#             return mol
#         except Exception as e:
#             print(f"Failed to load SDF directly: {e}")
#             print("Trying alternative loading method...")
            
#             try:
#                 from rdkit import Chem
#                 supplier = Chem.SDMolSupplier(self.ligand_path)
#                 rdmol = supplier[0] if supplier else None
                
#                 if not rdmol:
#                     raise ValueError("No valid molecule found in SDF")
                
#                 mol = Molecule.from_rdkit(rdmol)
#                 print(f"Loaded molecule via RDKit: {mol.n_atoms} atoms")
#                 return mol
#             except ImportError:
#                 raise RuntimeError("RDKit is required for this SDF file. Please install it.")
    
#     def _load_mol2(self):
#         """Load molecule from MOL2 file"""
#         print("Converting MOL2 to SDF for OpenFF...")
        
#         with tempfile.TemporaryDirectory() as temp_dir:
#             temp_sdf = os.path.join(temp_dir, "temp.sdf")
            
#             # Try to use antechamber for conversion
#             try:
#                 import subprocess
#                 cmd = [
#                     'antechamber',
#                     '-i', self.ligand_path,
#                     '-fi', 'mol2',
#                     '-o', temp_sdf,
#                     '-fo', 'sdf',
#                     '-c', 'wc',  # Write charges
#                     '-nc', str(self.charge)
#                 ]
                
#                 result = subprocess.run(cmd, capture_output=True, text=True)
                
#                 if result.returncode == 0 and os.path.exists(temp_sdf):
#                     mol = Molecule.from_file(temp_sdf)
#                     print(f"Converted MOL2 to SDF: {mol.n_atoms} atoms")
#                     return mol
#                 else:
#                     raise RuntimeError(f"Antechamber conversion failed: {result.stderr}")
                    
#             except Exception as e:
#                 print(f"Antechamber conversion failed: {e}")
                
#                 # Try Open Babel as fallback
#                 try:
#                     import subprocess
#                     cmd = [
#                         'obabel',
#                         '-i', 'mol2', self.ligand_path,
#                         '-o', 'sdf', '-O', temp_sdf,
#                         '-h'  # Add hydrogens
#                     ]
                    
#                     result = subprocess.run(cmd, capture_output=True, text=True)
                    
#                     if result.returncode == 0 and os.path.exists(temp_sdf):
#                         mol = Molecule.from_file(temp_sdf)
#                         print(f"Converted MOL2 to SDF using Open Babel: {mol.n_atoms} atoms")
#                         return mol
#                     else:
#                         raise RuntimeError("Open Babel conversion failed")
                        
#                 except Exception:
#                     raise RuntimeError("Could not convert MOL2 file. Please install antechamber or obabel.")
    
#     def _generate_parameters(self, mol):
#         """Generate OpenFF parameters and write GROMACS files"""
#         print("Generating OpenFF parameters...")
        
#         # Load force field
#         ff = ForceField(f"{self.forcefield_name}.offxml")
        
#         # Create interchange
#         print("Creating OpenFF Interchange...")
#         interchange = Interchange.from_smirnoff(ff, [mol])
        
#         # Add box vectors
#         interchange = self._add_box_vectors(interchange)
        
#         # Write GROMACS files
#         print("Writing GROMACS files...")
#         self._write_gromacs_files(interchange)
    
#     def _add_box_vectors(self, interchange):
#         """Add box vectors to the system"""
#         coords = interchange.positions.m_as(unit.angstrom)
#         min_c, max_c = np.min(coords, axis=0), np.max(coords, axis=0)
#         box_nm = (max_c - min_c + 2 * self.box_size * 10) * 0.1
#         interchange.box = np.diag(box_nm) * unit.nanometer
#         return interchange
    
#     def _write_gromacs_files(self, interchange):
#         """Write GROMACS format files"""
#         output_prefix = os.path.join(self.lig_ff_dir, "LIG")
        
#         # Write GRO file
#         interchange.to_gro(f"{output_prefix}.gro")
        
#         # Write TOP file
#         interchange.to_top(f"{output_prefix}.top")
        
#         # Process the topology files
#         self._process_topology_files(output_prefix)
    
#     def _process_topology_files(self, output_prefix):
#         """Process and split topology files"""
#         top_file = f"{output_prefix}.top"
        
#         with open(top_file, 'r') as f:
#             top_content = f.read()
        
#         # Extract different sections
#         atomtypes_content = self._extract_atomtypes_content(top_content)
#         itp_content = self._extract_itp_content(top_content)
        
#         # Modify ITP content to use LIG as resname
#         modified_itp_content = self._modify_resname_to_lig(itp_content)
        
#         # Count atoms for position restraints
#         atom_count = self._count_atoms_in_itp(modified_itp_content)
        
#         # Add position restraints include at the end of ITP
#         while modified_itp_content.endswith("\n\n"):
#             modified_itp_content = modified_itp_content[:-1]
#         modified_itp_content += "\n\n#ifdef POSRES\n#include \"posre_LIG.itp\"\n#endif\n"
        
#         # Write ITP file
#         itp_file = f"{output_prefix}.itp"
#         with open(itp_file, "w") as f:
#             f.write(modified_itp_content)
        
#         # Write atomtypes file
#         atomtypes_file = os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")
#         with open(atomtypes_file, "w") as f:
#             f.write(atomtypes_content)
        
#         # Write position restraints file
#         posre_file = os.path.join(self.lig_ff_dir, "posre_LIG.itp")
#         self._write_position_restraints(posre_file, atom_count)
        
#         # Rebuild top file
#         self._rebuild_top_file(top_content, itp_file, output_prefix)
        
#         # Modify GRO file to use LIG residue name
#         self._modify_gro_file(output_prefix)
    
#     def _extract_atomtypes_content(self, top_content):
#         """Extract atomtypes section from topology"""
#         lines = top_content.split("\n")
#         atomtypes_lines = []
#         in_section = False
        
#         for line in lines:
#             stripped = line.strip()
#             if stripped.startswith("[ atomtypes ]"):
#                 in_section = True
#                 atomtypes_lines.append(line)
#             elif stripped.startswith("[") and in_section:
#                 in_section = False
#             elif in_section:
#                 atomtypes_lines.append(line)
        
#         return "\n".join(atomtypes_lines)
    
#     def _extract_itp_content(self, top_content):
#         """Extract moleculetype section from topology"""
#         lines = top_content.split("\n")
#         itp_lines = []
#         in_section = False
        
#         for line in lines:
#             stripped = line.strip()
#             if stripped.startswith("[ moleculetype ]"):
#                 in_section = True
#             elif stripped.startswith(("[ system ]", "[ molecules ]")):
#                 in_section = False
            
#             if in_section:
#                 itp_lines.append(line)
        
#         return "\n".join(itp_lines)
    
#     def _modify_resname_to_lig(self, itp_content):
#         """Modify residue names to LIG - FIXED VERSION"""
#         lines = itp_content.split("\n")
#         modified_lines = []
#         current_section = None
        
#         for line in lines:
#             stripped = line.strip()
            
#             # Detect section changes
#             if stripped.startswith("[") and stripped.endswith("]"):
#                 current_section = stripped[1:-1].strip()
#                 modified_lines.append(line)
#                 continue
            
#             # Skip comments and empty lines
#             if not stripped or stripped.startswith(";"):
#                 modified_lines.append(line)
#                 continue
            
#             # Handle different sections
#             if current_section == "moleculetype":
#                 # Replace molecule name with LIG
#                 parts = line.split()
#                 if len(parts) >= 1:
#                     # Keep the original formatting but replace the first part with LIG
#                     remaining_parts = ' '.join(parts[1:]) if len(parts) > 1 else ""
#                     line = f"LIG              {remaining_parts}".rstrip()
#                 modified_lines.append(line)
            
#             elif current_section == "atoms":
#                 # Replace resname with LIG in atoms section
#                 parts = line.split()
#                 if len(parts) >= 8:
#                     # Format: index, atom_type, resnum, resname, name, cgnr, charge, mass
#                     parts[3] = "LIG"  # resname is the 4th field
#                     # Reconstruct line with proper formatting
#                     line = f"{parts[0]:>6} {parts[1]:<15} {parts[2]:>4} {parts[3]:<7} {parts[4]:<10} {parts[5]:>4} {parts[6]:>17} {parts[7]:>17}"
#                 modified_lines.append(line)
            
#             else:
#                 # For all other sections (bonds, angles, dihedrals, etc.), keep the line as is
#                 # This is crucial - no modifications should be made to other sections
#                 modified_lines.append(line)
        
#         return "\n".join(modified_lines)
    
#     def _count_atoms_in_itp(self, itp_content):
#         """Count atoms in ITP file"""
#         lines = itp_content.split("\n")
#         current_section = None
#         atom_count = 0
        
#         for line in lines:
#             stripped = line.strip()
            
#             # Detect section changes
#             if stripped.startswith("[") and stripped.endswith("]"):
#                 current_section = stripped[1:-1].strip()
#                 continue
            
#             # Only count in atoms section
#             if current_section == "atoms" and stripped and not stripped.startswith(";"):
#                 parts = stripped.split()
#                 if len(parts) >= 8:
#                     try:
#                         atom_index = int(parts[0])
#                         atom_count = max(atom_count, atom_index)
#                     except (ValueError, IndexError):
#                         pass
        
#         return atom_count
    
#     def _write_position_restraints(self, posre_file, atom_count):
#         """Write position restraints file"""
#         content = "[ position_restraints ]\n"
#         content += "; atom  type      fx      fy      fz\n"
        
#         for i in range(1, atom_count + 1):
#             content += f"{i:>6}     1  1000  1000  1000\n"
        
#         with open(posre_file, "w") as f:
#             f.write(content)
    
#     def _rebuild_top_file(self, top_content, itp_file, output_prefix):
#         """Rebuild topology file with proper includes"""
#         itp_name = os.path.basename(itp_file)
#         lines = top_content.split("\n")
#         new_lines = []
#         skip_section = False
#         skip_atomtypes = False
#         added_includes = False
        
#         for line in lines:
#             stripped = line.strip()
            
#             # Skip atomtypes section
#             if stripped.startswith("[ atomtypes ]"):
#                 skip_atomtypes = True
#                 if not added_includes:
#                     # Add includes after defaults section
#                     new_lines.append("")
#                     new_lines.append("#include \"atomtypes_LIG.itp\"")
#                     new_lines.append(f"#include \"{itp_name}\"")
#                     new_lines.append("")
#                     added_includes = True
#                 continue
#             elif skip_atomtypes and stripped.startswith("[") and not stripped.startswith("[ atomtypes ]"):
#                 skip_atomtypes = False
            
#             # Skip moleculetype and related sections
#             if stripped.startswith("[ moleculetype ]"):
#                 skip_section = True
#                 continue
#             elif stripped.startswith(("[ system ]", "[ molecules ]")):
#                 skip_section = False
#                 new_lines.append(line)
#             elif not skip_section and not skip_atomtypes:
#                 new_lines.append(line)
        
#         # Fix molecules section to use LIG
#         final_lines = []
#         in_molecules = False
#         in_system = False
#         for line in new_lines:
#             if line.strip().startswith("[ system ]"):
#                 in_system = True
#                 final_lines.append(line)
#             elif in_system and line.strip() and not line.strip().startswith(";"):
#                 # Replace system name
#                 final_lines.append("LIG system")
#                 in_system = False
#             elif line.strip().startswith("[ molecules ]"):
#                 in_molecules = True
#                 final_lines.append(line)
#             elif in_molecules and line.strip() and not line.strip().startswith(";"):
#                 # Replace molecule name with LIG
#                 parts = line.split()
#                 if len(parts) >= 2:
#                     line = f"LIG              {parts[1]}"
#                 final_lines.append(line)
#                 in_molecules = False
#             else:
#                 final_lines.append(line)
        
#         with open(f"{output_prefix}.top", "w") as f:
#             f.write("\n".join(final_lines))
    
#     def _modify_gro_file(self, output_prefix):
#         """Modify GRO file to use LIG residue name"""
#         gro_file = f"{output_prefix}.gro"
        
#         with open(gro_file, 'r') as f:
#             lines = f.readlines()
        
#         # Modify residue names in GRO file
#         modified_lines = []
#         for i, line in enumerate(lines):
#             if i == 0:
#                 # First line is the title
#                 modified_lines.append("LIG system\n")
#             elif i >= 2 and i < len(lines) - 1:  # Skip header and box vectors
#                 # GRO format: residue name is columns 6-10 (0-indexed: 5-10)
#                 if len(line) > 10:
#                     line = line[:5] + "LIG".ljust(5) + line[10:]
#                 modified_lines.append(line)
#             else:
#                 modified_lines.append(line)
        
#         with open(gro_file, 'w') as f:
#             f.writelines(modified_lines)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenFF force field generator wrapper for PRISM
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# Import the base class
try:
    from .base import ForceFieldGeneratorBase
except ImportError:
    from base import ForceFieldGeneratorBase

# Import OpenFF specific dependencies
try:
    from openff.toolkit import Molecule, ForceField
    from openff.interchange import Interchange
    from openff.units import unit
    import numpy as np
except ImportError as e:
    print(f"Error: Missing OpenFF dependencies: {e}")
    print("Please install: pip install openff-toolkit openff-interchange")
    sys.exit(1)


class OpenFFForceFieldGenerator(ForceFieldGeneratorBase):
    """OpenFF force field generator wrapper"""
    
    def __init__(self, ligand_path, output_dir, charge=0, forcefield="openff-2.1.0", overwrite=False):
        """
        Initialize OpenFF force field generator
        
        Parameters:
        -----------
        ligand_path : str
            Path to the ligand file (SDF/MOL2)
        output_dir : str
            Directory where output files will be stored
        charge : int
            Molecule charge (default: 0)
        forcefield : str
            OpenFF force field version (default: "openff-2.1.0")
        overwrite : bool
            Whether to overwrite existing files
        """
        super().__init__(ligand_path, output_dir, overwrite)
        
        self.charge = charge
        self.forcefield_name = forcefield
        self.box_size = 1.0  # nm
        
        # Output directory for OpenFF files
        self.lig_ff_dir = os.path.join(self.output_dir, self.get_output_dir_name())
        
        print(f"\nInitialized OpenFF Force Field Generator:")
        print(f"  Ligand: {self.ligand_path}")
        print(f"  Charge: {self.charge}")
        print(f"  Force field: {self.forcefield_name}")
        print(f"  Output directory: {self.output_dir}")
    
    def get_output_dir_name(self):
        """Get the output directory name for OpenFF"""
        return "LIG.openff2gmx"
    
    def run(self):
        """Run the OpenFF force field generation workflow"""
        print(f"\n{'='*60}")
        print("Starting OpenFF Force Field Generation")
        print(f"{'='*60}")
        
        try:
            # Check if output already exists
            if os.path.exists(self.lig_ff_dir) and not self.overwrite:
                if self.check_required_files(self.lig_ff_dir):
                    print(f"Output directory {self.lig_ff_dir} already exists with all required files.")
                    print("Skipping force field generation (use --overwrite to regenerate)")
                    return self.lig_ff_dir
            
            # Create output directory
            os.makedirs(self.lig_ff_dir, exist_ok=True)
            
            # Load molecule
            mol = self._load_molecule()
            
            # Generate force field parameters
            self._generate_parameters(mol)
            
            print(f"\n{'='*60}")
            print("Force field generation completed successfully!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {self.lig_ff_dir}")
            print("\nGenerated files:")
            for f in os.listdir(self.lig_ff_dir):
                print(f"  - {f}")
            
            return self.lig_ff_dir
            
        except Exception as e:
            print(f"\nError during force field generation: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    def _load_molecule(self):
        """Load molecule from file - IMPROVED VERSION"""
        print("Loading molecule...")
        
        file_ext = Path(self.ligand_path).suffix.lower()
        
        if file_ext == '.sdf':
            return self._load_sdf()
        elif file_ext == '.mol2':
            return self._load_mol2()
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
    
    def _load_sdf(self):
        """Load molecule from SDF file - IMPROVED VERSION"""
        try:
            mol = Molecule.from_file(self.ligand_path)
            print(f"Loaded molecule from SDF: {mol.n_atoms} atoms")
            return mol
        except Exception as e:
            print(f"Failed to load SDF directly: {e}")
            print("Trying RDKit loading method...")
            
            try:
                from rdkit import Chem
                supplier = Chem.SDMolSupplier(self.ligand_path)
                rdmol = supplier[0] if supplier else None
                
                if not rdmol:
                    raise ValueError("No valid molecule found in SDF")
                
                mol = Molecule.from_rdkit(rdmol)
                print(f"Loaded molecule via RDKit: {mol.n_atoms} atoms")
                return mol
            except ImportError:
                raise RuntimeError("RDKit is required for this SDF file. Please install it.")
            except Exception as rdkit_e:
                raise RuntimeError(f"Both OpenFF and RDKit failed to load SDF: {e}, {rdkit_e}")
    
    def _load_mol2(self):
        """Load molecule from MOL2 file - IMPROVED VERSION"""
        print("Converting MOL2 to SDF for OpenFF...")
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Copy input file to temp directory to avoid path issues
            temp_mol2 = os.path.join(temp_dir, "input.mol2")
            temp_sdf = os.path.join(temp_dir, "output.sdf")
            
            shutil.copy2(self.ligand_path, temp_mol2)
            
            # Try antechamber first with improved parameters
            try:
                import subprocess
                cmd = [
                    'antechamber',
                    '-i', 'input.mol2',   # Use relative path in temp dir
                    '-fi', 'mol2',
                    '-o', 'output.sdf',   # Use relative path in temp dir
                    '-fo', 'sdf',
                    '-c', 'bcc',          # Use BCC charge method instead of 'wc'
                    '-nc', str(self.charge)
                ]
                
                print(f"Running antechamber with command: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)
                
                if result.returncode == 0 and os.path.exists(temp_sdf):
                    mol = Molecule.from_file(temp_sdf)
                    print(f"Converted MOL2 to SDF using antechamber: {mol.n_atoms} atoms")
                    return mol
                else:
                    print(f"Antechamber failed with return code {result.returncode}")
                    print(f"STDERR: {result.stderr}")
                    print(f"STDOUT: {result.stdout}")
                    raise RuntimeError(f"Antechamber conversion failed")
                    
            except FileNotFoundError:
                print("Antechamber not found, trying Open Babel...")
            except Exception as e:
                print(f"Antechamber conversion failed: {e}")
            
            # Try Open Babel as fallback
            try:
                import subprocess
                cmd = [
                    'obabel',
                    '-i', 'mol2', 'input.mol2',
                    '-o', 'sdf', '-O', 'output.sdf',
                    '-h'  # Add hydrogens
                ]
                
                print(f"Running Open Babel with command: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)
                
                if result.returncode == 0 and os.path.exists(temp_sdf):
                    mol = Molecule.from_file(temp_sdf)
                    print(f"Converted MOL2 to SDF using Open Babel: {mol.n_atoms} atoms")
                    return mol
                else:
                    print(f"Open Babel failed with return code {result.returncode}")
                    print(f"STDERR: {result.stderr}")
                    raise RuntimeError("Open Babel conversion failed")
                    
            except FileNotFoundError:
                raise RuntimeError("Neither antechamber nor obabel found. Please install AmberTools or Open Babel.")
            except Exception as e:
                raise RuntimeError(f"Open Babel conversion failed: {e}")
            
            # If all methods fail
            raise RuntimeError("Could not convert MOL2 file. Please check the input file and ensure antechamber or obabel is installed.")
    
    def _generate_parameters(self, mol):
        """Generate OpenFF parameters and write GROMACS files"""
        print("Generating OpenFF parameters...")
        
        # Load force field
        ff = ForceField(f"{self.forcefield_name}.offxml")
        
        # Create interchange
        print("Creating OpenFF Interchange...")
        interchange = Interchange.from_smirnoff(ff, [mol])
        
        # Add box vectors
        interchange = self._add_box_vectors(interchange)
        
        # Write GROMACS files
        print("Writing GROMACS files...")
        self._write_gromacs_files(interchange)
    
    def _add_box_vectors(self, interchange):
        """Add box vectors to the system"""
        coords = interchange.positions.m_as(unit.angstrom)
        min_c, max_c = np.min(coords, axis=0), np.max(coords, axis=0)
        box_nm = (max_c - min_c + 2 * self.box_size * 10) * 0.1
        interchange.box = np.diag(box_nm) * unit.nanometer
        return interchange
    
    def _write_gromacs_files(self, interchange):
        """Write GROMACS format files"""
        output_prefix = os.path.join(self.lig_ff_dir, "LIG")
        
        # Write GRO file
        interchange.to_gro(f"{output_prefix}.gro")
        
        # Write TOP file
        interchange.to_top(f"{output_prefix}.top")
        
        # Process the topology files
        self._process_topology_files(output_prefix)
    
    def _process_topology_files(self, output_prefix):
        """Process and split topology files"""
        top_file = f"{output_prefix}.top"
        
        with open(top_file, 'r') as f:
            top_content = f.read()
        
        # Extract different sections
        atomtypes_content = self._extract_atomtypes_content(top_content)
        itp_content = self._extract_itp_content(top_content)
        
        # Modify ITP content to use LIG as resname
        modified_itp_content = self._modify_resname_to_lig(itp_content)
        
        # Count atoms for position restraints
        atom_count = self._count_atoms_in_itp(modified_itp_content)
        
        # Add position restraints include at the end of ITP
        while modified_itp_content.endswith("\n\n"):
            modified_itp_content = modified_itp_content[:-1]
        modified_itp_content += "\n\n#ifdef POSRES\n#include \"posre_LIG.itp\"\n#endif\n"
        
        # Write ITP file
        itp_file = f"{output_prefix}.itp"
        with open(itp_file, "w") as f:
            f.write(modified_itp_content)
        
        # Write atomtypes file
        atomtypes_file = os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")
        with open(atomtypes_file, "w") as f:
            f.write(atomtypes_content)
        
        # Write position restraints file
        posre_file = os.path.join(self.lig_ff_dir, "posre_LIG.itp")
        self._write_position_restraints(posre_file, atom_count)
        
        # Rebuild top file
        self._rebuild_top_file(top_content, itp_file, output_prefix)
        
        # Modify GRO file to use LIG residue name
        self._modify_gro_file(output_prefix)
    
    def _extract_atomtypes_content(self, top_content):
        """Extract atomtypes section from topology"""
        lines = top_content.split("\n")
        atomtypes_lines = []
        in_section = False
        
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("[ atomtypes ]"):
                in_section = True
                atomtypes_lines.append(line)
            elif stripped.startswith("[") and in_section:
                in_section = False
            elif in_section:
                atomtypes_lines.append(line)
        
        return "\n".join(atomtypes_lines)
    
    def _extract_itp_content(self, top_content):
        """Extract moleculetype section from topology"""
        lines = top_content.split("\n")
        itp_lines = []
        in_section = False
        
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("[ moleculetype ]"):
                in_section = True
            elif stripped.startswith(("[ system ]", "[ molecules ]")):
                in_section = False
            
            if in_section:
                itp_lines.append(line)
        
        return "\n".join(itp_lines)
    
    def _modify_resname_to_lig(self, itp_content):
        """Modify residue names to LIG - FIXED VERSION"""
        lines = itp_content.split("\n")
        modified_lines = []
        current_section = None
        
        for line in lines:
            stripped = line.strip()
            
            # Detect section changes
            if stripped.startswith("[") and stripped.endswith("]"):
                current_section = stripped[1:-1].strip()
                modified_lines.append(line)
                continue
            
            # Skip comments and empty lines
            if not stripped or stripped.startswith(";"):
                modified_lines.append(line)
                continue
            
            # Handle different sections
            if current_section == "moleculetype":
                # Replace molecule name with LIG
                parts = line.split()
                if len(parts) >= 1:
                    # Keep the original formatting but replace the first part with LIG
                    remaining_parts = ' '.join(parts[1:]) if len(parts) > 1 else ""
                    line = f"LIG              {remaining_parts}".rstrip()
                modified_lines.append(line)
            
            elif current_section == "atoms":
                # Replace resname with LIG in atoms section
                parts = line.split()
                if len(parts) >= 8:
                    # Format: index, atom_type, resnum, resname, name, cgnr, charge, mass
                    parts[3] = "LIG"  # resname is the 4th field
                    # Reconstruct line with proper formatting
                    line = f"{parts[0]:>6} {parts[1]:<15} {parts[2]:>4} {parts[3]:<7} {parts[4]:<10} {parts[5]:>4} {parts[6]:>17} {parts[7]:>17}"
                modified_lines.append(line)
            
            else:
                # For all other sections (bonds, angles, dihedrals, etc.), keep the line as is
                # This is crucial - no modifications should be made to other sections
                modified_lines.append(line)
        
        return "\n".join(modified_lines)
    
    def _count_atoms_in_itp(self, itp_content):
        """Count atoms in ITP file"""
        lines = itp_content.split("\n")
        current_section = None
        atom_count = 0
        
        for line in lines:
            stripped = line.strip()
            
            # Detect section changes
            if stripped.startswith("[") and stripped.endswith("]"):
                current_section = stripped[1:-1].strip()
                continue
            
            # Only count in atoms section
            if current_section == "atoms" and stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if len(parts) >= 8:
                    try:
                        atom_index = int(parts[0])
                        atom_count = max(atom_count, atom_index)
                    except (ValueError, IndexError):
                        pass
        
        return atom_count
    
    def _write_position_restraints(self, posre_file, atom_count):
        """Write position restraints file"""
        content = "[ position_restraints ]\n"
        content += "; atom  type      fx      fy      fz\n"
        
        for i in range(1, atom_count + 1):
            content += f"{i:>6}     1  1000  1000  1000\n"
        
        with open(posre_file, "w") as f:
            f.write(content)
    
    def _rebuild_top_file(self, top_content, itp_file, output_prefix):
        """Rebuild topology file with proper includes"""
        itp_name = os.path.basename(itp_file)
        lines = top_content.split("\n")
        new_lines = []
        skip_section = False
        skip_atomtypes = False
        added_includes = False
        
        for line in lines:
            stripped = line.strip()
            
            # Skip atomtypes section
            if stripped.startswith("[ atomtypes ]"):
                skip_atomtypes = True
                if not added_includes:
                    # Add includes after defaults section
                    new_lines.append("")
                    new_lines.append("#include \"atomtypes_LIG.itp\"")
                    new_lines.append(f"#include \"{itp_name}\"")
                    new_lines.append("")
                    added_includes = True
                continue
            elif skip_atomtypes and stripped.startswith("[") and not stripped.startswith("[ atomtypes ]"):
                skip_atomtypes = False
            
            # Skip moleculetype and related sections
            if stripped.startswith("[ moleculetype ]"):
                skip_section = True
                continue
            elif stripped.startswith(("[ system ]", "[ molecules ]")):
                skip_section = False
                new_lines.append(line)
            elif not skip_section and not skip_atomtypes:
                new_lines.append(line)
        
        # Fix molecules section to use LIG
        final_lines = []
        in_molecules = False
        in_system = False
        for line in new_lines:
            if line.strip().startswith("[ system ]"):
                in_system = True
                final_lines.append(line)
            elif in_system and line.strip() and not line.strip().startswith(";"):
                # Replace system name
                final_lines.append("LIG system")
                in_system = False
            elif line.strip().startswith("[ molecules ]"):
                in_molecules = True
                final_lines.append(line)
            elif in_molecules and line.strip() and not line.strip().startswith(";"):
                # Replace molecule name with LIG
                parts = line.split()
                if len(parts) >= 2:
                    line = f"LIG              {parts[1]}"
                final_lines.append(line)
                in_molecules = False
            else:
                final_lines.append(line)
        
        with open(f"{output_prefix}.top", "w") as f:
            f.write("\n".join(final_lines))
    
    def _modify_gro_file(self, output_prefix):
        """Modify GRO file to use LIG residue name"""
        gro_file = f"{output_prefix}.gro"
        
        with open(gro_file, 'r') as f:
            lines = f.readlines()
        
        # Modify residue names in GRO file
        modified_lines = []
        for i, line in enumerate(lines):
            if i == 0:
                # First line is the title
                modified_lines.append("LIG system\n")
            elif i >= 2 and i < len(lines) - 1:  # Skip header and box vectors
                # GRO format: residue name is columns 6-10 (0-indexed: 5-10)
                if len(line) > 10:
                    line = line[:5] + "LIG".ljust(5) + line[10:]
                modified_lines.append(line)
            else:
                modified_lines.append(line)
        
        with open(gro_file, 'w') as f:
            f.writelines(modified_lines)