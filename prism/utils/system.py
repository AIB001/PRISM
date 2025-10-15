# !filepath: prism/utils/system.py

# !/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM System Builder - Handles GROMACS system setup, solvation, and ion addition.
"""

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple


class SystemBuilder:
    """
    Builds a complete protein-ligand system for GROMACS simulations.

    This class encapsulates the GROMACS commands needed to process a protein,
    combine it with a ligand, create a simulation box, solvate it, and add ions.
    It is designed to be robust, especially in handling the output of GROMACS
    tools to ensure consistency between coordinate and topology files.
    """

    def __init__(self, config: Dict, output_dir: str, overwrite: bool = False):
        """
        Initializes the SystemBuilder.

        Args:
            config: A dictionary containing the simulation configuration.
            output_dir: The root directory for all output files.
            overwrite: If True, overwrite existing files.
        """
        self.config = config
        self.output_dir = Path(output_dir)
        self.overwrite = overwrite
        self.gmx_command = self.config.get("general", {}).get("gmx_command", "gmx")

        self.model_dir = self.output_dir / "GMX_PROLIG_MD"
        self.model_dir.mkdir(exist_ok=True)

    def _run_command(self, command: list, work_dir: str, input_str: Optional[str] = None) -> Tuple[str, str]:
        """
        Executes a shell command and handles errors.

        Args:
            command: The command to execute as a list of strings.
            work_dir: The directory in which to run the command.
            input_str: A string to be passed to the command's stdin.

        Returns:
            A tuple containing the stdout and stderr of the command.

        Raises:
            RuntimeError: If the command returns a non-zero exit code.
        """
        cmd_str = ' '.join(map(str, command))
        print(f"Executing in {work_dir}: {cmd_str}")

        process = subprocess.Popen(
            command,
            cwd=work_dir,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate(input=input_str)

        if process.returncode != 0:
            print("--- STDOUT ---")
            print(stdout)
            print("--- STDERR ---")
            print(stderr)
            raise RuntimeError(f"Command failed with exit code {process.returncode}: {cmd_str}")

        return stdout, stderr

    def build(self, cleaned_protein_path: str, lig_ff_dir: str, ff_idx: int, water_idx: int,
              ff_info: dict = None, water_info: dict = None) -> str:
        """
        Runs the full system building workflow.

        Args:
            cleaned_protein_path: Path to the cleaned protein PDB file.
            lig_ff_dir: Path to the directory containing ligand force field files.
            ff_idx: The index of the protein force field to use.
            water_idx: The index of the water model to use.
            ff_info: Dictionary containing force field info (name, dir, path)
            water_info: Dictionary containing water model info (name, label, description)

        Returns:
            The path to the directory containing the final model files.
        """
        try:
            print("\n=== Building GROMACS Model ===")

            if (self.model_dir / "solv_ions.gro").exists() and not self.overwrite:
                print(f"Final model file solv_ions.gro already exists in {self.model_dir}, skipping model building.")
                return str(self.model_dir)

            # Store water model name for use in solvation
            self.water_model_name = water_info.get('name', 'tip3p') if water_info else 'tip3p'

            fixed_pdb = self._fix_protein(cleaned_protein_path)
            protein_gro, topol_top = self._generate_topology(fixed_pdb, ff_idx, water_idx, ff_info, water_info)

            # Only process ligand if lig_ff_dir is provided (not None)
            if lig_ff_dir is not None:
                self._fix_topology(topol_top, lig_ff_dir)
                complex_gro = self._combine_protein_ligand(protein_gro, lig_ff_dir)
            else:
                # Protein-only system
                print("\nSkipping ligand processing (protein-only system)")
                complex_gro = protein_gro

            boxed_gro = self._create_box(complex_gro)
            solvated_gro = self._solvate(boxed_gro, topol_top)
            self._add_ions(solvated_gro, topol_top)

            print("\nSystem building completed successfully.")
            print(f"System files are in {self.model_dir}")
            return str(self.model_dir)

        except (RuntimeError, FileNotFoundError) as e:
            print(f"\nError during system build: {e}")
            import traceback
            traceback.print_exc()
            return ""

    def _fix_protein(self, cleaned_protein: str) -> str:
        """Fix missing atoms in protein using pdbfixer."""
        print("\nStep 1: Fixing protein with pdbfixer...")
        fixed_pdb = self.model_dir / "fixed_clean_protein.pdb"

        if fixed_pdb.exists() and not self.overwrite:
            print(f"Fixed protein PDB already exists at {fixed_pdb}, skipping pdbfixer.")
            return str(fixed_pdb)

        # Check if pdbfixer is available
        if not shutil.which("pdbfixer"):
            print("Warning: pdbfixer command not found. Skipping protein fixing step.")
            print("RECOMMENDED INSTALLATION:")
            print("  mamba install -c conda-forge pdbfixer")
            print("  # OR: conda install -c conda-forge pdbfixer")
            shutil.copy(cleaned_protein, fixed_pdb)
        else:
            # Run pdbfixer to fix missing atoms
            # Note: pdbfixer with --keep-heterogens=all will move metals to chain 'M'
            # This is OK - we extract them from chain 'M' after pdb2gmx
            command = [
                "pdbfixer",
                cleaned_protein,
                "--add-residues",
                "--output", str(fixed_pdb),
                "--add-atoms=heavy",
                "--keep-heterogens=all",
                f"--ph={self.config['simulation']['pH']}"
            ]
            self._run_command(command, str(self.model_dir.parent))

        # Fix terminal atoms for GROMACS compatibility (OXT → O)
        print("Fixing terminal atoms for GROMACS compatibility...")
        from prism.utils.cleaner import fix_terminal_atoms
        fix_terminal_atoms(str(fixed_pdb), str(fixed_pdb), verbose=True)

        return str(fixed_pdb)

    def _generate_topology(self, fixed_pdb: str, ff_idx: int, water_idx: int,
                          ff_info: dict = None, water_info: dict = None) -> Tuple[str, str]:
        """Generate topology using pdb2gmx, handling histidines."""
        print("\nStep 2: Generating topology with pdb2gmx...")

        protein_gro = self.model_dir / "pro.gro"
        topol_top = self.model_dir / "topol.top"

        if protein_gro.exists() and not self.overwrite:
            print(f"Protein topology already generated at {protein_gro}, skipping pdb2gmx.")
            return str(protein_gro), str(topol_top)

        # Build pdb2gmx command
        command = [
            self.gmx_command, "pdb2gmx",
            "-f", fixed_pdb,
            "-o", str(protein_gro),
            "-p", str(topol_top),
            "-ignh"
        ]

        # Use -ff and -water flags if force field info is provided
        # This is more reliable than interactive menu indexing
        if ff_info and 'dir' in ff_info:
            # Remove .ff extension from directory name for -ff flag
            ff_name = ff_info['dir'].replace('.ff', '')
            command.extend(["-ff", ff_name])
            print(f"Using force field: {ff_name} (from {ff_info.get('path', 'N/A')})")
        else:
            print(f"DEBUG: Force field index = {ff_idx} (interactive menu)")

        if water_info and 'name' in water_info:
            command.extend(["-water", water_info['name']])
            print(f"Using water model: {water_info['name']}")
        else:
            print(f"DEBUG: Water model index = {water_idx} (interactive menu)")

        # Handle histidine protonation states
        histidine_count = self._count_histidines(fixed_pdb)
        input_lines = []

        # If using flags, we don't need force field/water model selection
        if not (ff_info and water_info):
            input_lines = [str(ff_idx), str(water_idx)]

        if histidine_count > 0:
            print(f"Found {histidine_count} histidine residue(s), defaulting to HIE protonation state.")
            # HIE is typically the first and most common choice for neutral pH.
            input_lines.extend(["1"] * histidine_count)

        input_text = '\n'.join(input_lines)+'\n' if input_lines else None

        self._run_command(command, str(self.model_dir), input_str=input_text)

        # Check if pdb2gmx already handled metals (any ion chain topology exists)
        # PDBFixer may keep metals in original chains (D, F) or move to chain M
        # If pdb2gmx created ion chain topologies, metals are already included
        ion_chain_files = list(self.model_dir.glob("topol_Ion_chain_*.itp"))
        if ion_chain_files:
            print(f"\nMetals already processed by pdb2gmx ({len(ion_chain_files)} ion chain(s) found)")
        else:
            # pdb2gmx didn't process metals, add them manually
            self._add_metals_to_topology(fixed_pdb, str(protein_gro), str(topol_top))

        return str(protein_gro), str(topol_top)

    def _fix_topology(self, topol_path: str, lig_ff_dir: str):
        """Includes ligand parameters into the main topology file."""
        print("\nStep 3: Including ligand parameters in topology...")

        lig_ff_path = Path(lig_ff_dir)
        lig_itp_path = lig_ff_path / "LIG.itp"
        atomtypes_itp_path = lig_ff_path / "atomtypes_LIG.itp"

        if not lig_itp_path.exists():
            raise FileNotFoundError(f"Ligand ITP file not found: {lig_itp_path}")
        if not atomtypes_itp_path.exists():
            raise FileNotFoundError(f"Ligand atomtypes file not found: {atomtypes_itp_path}")

        with open(topol_path, 'r') as f:
            lines = f.readlines()

        # Check if already included to prevent duplication
        if any(f'#include "../{lig_ff_path.name}/LIG.itp"' in line for line in lines):
            print("Topology already includes ligand parameters.")
            return

        with open(atomtypes_itp_path, 'r') as f:
            atomtypes_content = f.read()

        new_lines = []
        molecules_added = False
        for line in lines:
            new_lines.append(line)
            # Insert atomtypes and ligand ITP after the main forcefield include
            if '#include' in line and '.ff/forcefield.itp' in line:
                new_lines.append("\n; Include ligand atomtypes\n")
                new_lines.append(atomtypes_content)
                new_lines.append("\n; Include ligand topology\n")
                new_lines.append(f'#include "../{lig_ff_path.name}/LIG.itp"\n')

            # Add ligand to the molecules section
            if line.strip() == "[ molecules ]":
                # This flag ensures we add the ligand after the [ molecules ] header
                molecules_added = True
            elif molecules_added and line.strip() and not line.strip().startswith(';'):
                # Insert LIG before the first molecule (e.g., Protein_chain_A)
                new_lines.insert(-1, "LIG                 1\n")
                molecules_added = False  # Prevent adding it again

        with open(topol_path, 'w') as f:
            f.writelines(new_lines)
        print("Topology file updated with ligand parameters.")

    def _combine_protein_ligand(self, protein_gro: str, lig_ff_dir: str) -> str:
        """Combines protein and ligand GRO files into a complex."""
        print("\nStep 4: Combining protein and ligand coordinates...")
        complex_gro = self.model_dir / "pro_lig.gro"
        lig_gro_path = Path(lig_ff_dir) / "LIG.gro"

        if complex_gro.exists() and not self.overwrite:
            print("Using existing complex.gro file.")
            return str(complex_gro)

        with open(protein_gro, "r") as f_prot, open(lig_gro_path, "r") as f_lig, open(complex_gro, "w") as f_out:
            prot_lines = f_prot.readlines()
            lig_lines = f_lig.readlines()

            total_atoms = int(prot_lines[1].strip())+int(lig_lines[1].strip())

            f_out.write("Complex Protein-Ligand\n")
            f_out.write(f"{total_atoms:5d}\n")
            
            # Write ligand coordinates first (to match topology molecule order)
            lig_coords = lig_lines[2:-1]
            f_out.writelines(lig_coords)
            
            # Write protein coordinates second, renumbering atoms sequentially
            prot_coords = prot_lines[2:-1]
            lig_atom_count = len(lig_coords)
            
            for line in prot_coords:
                # Update atom numbering to be sequential after ligand atoms
                parts = list(line)
                # Extract current atom number (positions 15-20 in GRO format)
                current_num = int(line[15:20].strip())
                new_num = current_num + lig_atom_count
                # Replace atom number in the line
                new_line = line[:15] + f"{new_num:5d}" + line[20:]
                f_out.write(new_line)
                
            f_out.write(prot_lines[-1])  # Box vectors

        return str(complex_gro)

    def _create_box(self, complex_gro: str) -> str:
        """Creates the simulation box."""
        print("\nStep 5: Creating simulation box...")
        boxed_gro = self.model_dir / "pro_lig_newbox.gro"

        if boxed_gro.exists() and not self.overwrite:
            print("Using existing boxed.gro file.")
            return str(boxed_gro)

        box_cfg = self.config['box']
        command = [
            self.gmx_command, "editconf",
            "-f", complex_gro,
            "-o", str(boxed_gro),
            "-bt", box_cfg['shape'],
            "-d", str(box_cfg['distance'])
        ]
        if box_cfg.get('center', True):
            command.append("-c")

        self._run_command(command, str(self.model_dir))
        return str(boxed_gro)

    def _solvate(self, boxed_gro: str, topol_top: str) -> str:
        """Solvates the system."""
        print("\nStep 6: Solvating the system...")
        solvated_gro = self.model_dir / "solv.gro"

        if solvated_gro.exists() and not self.overwrite:
            print("Using existing solv.gro file.")
            return str(solvated_gro)

        # Select appropriate water model coordinate file
        # Map water model names to their coordinate files
        water_coord_files = {
            'tip3p': 'spc216.gro',  # TIP3P uses SPC geometry
            'tip3p_original': 'spc216.gro',
            'spc': 'spc216.gro',
            'spce': 'spc216.gro',
            'tip4p': 'tip4p.gro',  # TIP4P has its own file with virtual sites
            'tip4pew': 'tip4p.gro',  # TIP4P/Ew uses same geometry
            'tip5p': 'tip5p.gro'  # TIP5P has its own file
        }

        water_model = getattr(self, 'water_model_name', 'tip3p').lower()
        water_file = water_coord_files.get(water_model, 'spc216.gro')

        print(f"Using water model: {water_model} (coordinate file: {water_file})")

        command = [
            self.gmx_command, "solvate",
            "-cp", boxed_gro,
            "-cs", water_file,
            "-o", str(solvated_gro),
            "-p", topol_top
        ]
        self._run_command(command, str(self.model_dir))
        return str(solvated_gro)

    def _add_ions(self, solvated_gro: str, topol_top: str):
        """
        Adds ions to neutralize the system and achieve a target salt concentration.
        This method uses a single, robust call to `gmx genion`.
        """
        print("\nStep 7: Adding ions...")
        ions_gro = self.model_dir / "solv_ions.gro"

        if ions_gro.exists() and not self.overwrite:
            print("Using existing solv_ions.gro file.")
            return

        ions_cfg = self.config['ions']
        temp_tpr = self.model_dir / "ions.tpr"

        # A minimal mdp is needed for grompp
        ions_mdp_path = self.model_dir / "ions.mdp"
        with open(ions_mdp_path, "w") as f:
            f.write("integrator=steep\nnsteps=0")

        grompp_cmd = [
            self.gmx_command, "grompp",
            "-f", str(ions_mdp_path),
            "-c", solvated_gro,
            "-p", topol_top,
            "-o", str(temp_tpr),
            "-maxwarn", "5"  # Allow some warnings
        ]
        self._run_command(grompp_cmd, str(self.model_dir))

        genion_cmd = [
            self.gmx_command, "genion",
            "-s", str(temp_tpr),
            "-o", str(ions_gro),
            "-p", topol_top,
            "-pname", ions_cfg['positive_ion'],
            "-nname", ions_cfg['negative_ion'],
        ]
        if ions_cfg.get('neutral', True):
            genion_cmd.append("-neutral")
        if ions_cfg.get('concentration', 0) > 0:
            genion_cmd.extend(["-conc", str(ions_cfg['concentration'])])

        # Provide "SOL" (the typical name for water) as input to stdin
        stdout, _ = self._run_command(genion_cmd, str(self.model_dir), input_str="SOL")

        # --- CRITICAL STEP: Parse genion output and fix topology ---
        self._update_topology_molecules(topol_top, stdout)

        # Clean up temporary files
        temp_tpr.unlink()
        ions_mdp_path.unlink()

    def _update_topology_molecules(self, topol_path: str, genion_stdout: str):
        """
        Parses the stdout of `gmx genion` and updates the `[ molecules ]`
        section of the topology file to ensure consistency.
        
        Modern GROMACS versions (2024+) automatically update the topology when using -p flag,
        so we verify the update was successful rather than parsing stdout.
        """
        print("Checking topology update from genion...")

        # First try to parse traditional output format for backward compatibility
        mol_counts = {}
        pattern = re.compile(r'^\s*([A-Z0-9_]+)\s+:\s+(\d+)')
        parsing = False
        
        for line in genion_stdout.splitlines():
            if 'Number of molecules:' in line:
                parsing = True
                continue
            if parsing:
                match = pattern.match(line)
                if match:
                    mol_name, count = match.groups()
                    mol_counts[mol_name] = int(count)
                elif not line.strip():
                    break

        # If parsing failed, check if GROMACS already updated the topology (modern versions)
        if not mol_counts:
            print("Modern GROMACS detected - topology should be auto-updated by genion -p flag.")
            # Verify the topology has ions by checking the [molecules] section
            with open(topol_path, "r") as f:
                content = f.read()
            
            if any(ion in content for ion in ['NA ', 'CL ', 'K ', 'BR ', 'CA ']):
                print("Topology successfully updated with ions.")
                return
            else:
                print("Warning: Could not verify ion addition in topology. Manual check recommended.")
                return

        with open(topol_path, "r") as f:
            top_lines = f.readlines()

        # Find the [ molecules ] section
        start_index, end_index = -1, -1
        for i, line in enumerate(top_lines):
            if line.strip() == "[ molecules ]":
                start_index = i
            elif start_index != -1 and line.strip().startswith('['):
                end_index = i
                break
        if start_index != -1 and end_index == -1:
            end_index = len(top_lines)

        if start_index == -1:
            raise RuntimeError("Could not find [ molecules ] section in topology file.")

        # Reconstruct the section
        new_molecules_section = [top_lines[start_index]]

        # Preserve non-ion/water molecules from the original file
        original_mols = {}
        for i in range(start_index+1, end_index):
            line = top_lines[i].strip()
            if not line or line.startswith(';'): continue
            parts = line.split()
            if len(parts) >= 2:
                mol_name = parts[0]
                if mol_name not in mol_counts:
                    original_mols[mol_name] = parts[1]

        # Write preserved molecules first
        for mol, count in original_mols.items():
            new_molecules_section.append(f"{mol:<18} {count}\n")

        # Then write the updated counts from genion
        for mol, count in mol_counts.items():
            if count > 0:
                new_molecules_section.append(f"{mol:<18} {count}\n")

        # Replace the old section with the new one
        final_lines = top_lines[:start_index]+new_molecules_section+top_lines[end_index:]

        with open(topol_path, "w") as f:
            f.writelines(final_lines)

        print("Topology `[ molecules ]` section has been successfully updated.")

    def _count_histidines(self, pdb_file: str) -> int:
        """Counts unique histidine residues in a PDB file."""
        his_residues = set()
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[17:20].strip() in ['HIS', 'HID', 'HIE', 'HIP']:
                    resnum = line[22:26].strip()
                    chain = line[21].strip()
                    res_id = f"{chain}:{resnum}" if chain else resnum
                    his_residues.add(res_id)
        return len(his_residues)

    def _extract_metals_from_pdb(self, pdb_file: str) -> list:
        """
        Extract metal ion information from PDB file.

        Returns:
            List of dicts with keys: residue_name, atom_name, chain, resnum, coords
        """
        metals = []
        metal_names = {'ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'CO', 'NI', 'CD', 'HG',
                      'ZN2', 'MG2', 'CA2', 'FE2', 'FE3', 'CU2', 'MN2'}

        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    residue_name = line[17:20].strip().upper()
                    atom_name = line[12:16].strip().upper()

                    # Check if this is a metal
                    if residue_name in metal_names or atom_name in metal_names:
                        chain = line[21].strip()
                        resnum = line[22:26].strip()
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        metals.append({
                            'residue_name': residue_name,
                            'atom_name': atom_name,
                            'chain': chain,
                            'resnum': resnum,
                            'coords': (x, y, z)
                        })

        return metals

    def _add_metals_to_topology(self, pdb_file: str, gro_file: str, topol_file: str):
        """
        Add metal ions to topology and coordinates after pdb2gmx.

        pdb2gmx only processes ATOM records (protein), so we need to manually
        add HETATM records (metals) to both topology and coordinates.
        """
        metals = self._extract_metals_from_pdb(pdb_file)

        if not metals:
            return  # No metals to add

        print(f"\nAdding {len(metals)} metal ion(s) to topology and coordinates...")

        # Add metals to topology
        with open(topol_file, 'r') as f:
            topol_lines = f.readlines()

        # Find [ molecules ] section
        molecules_idx = -1
        for i, line in enumerate(topol_lines):
            if line.strip() == '[ molecules ]':
                molecules_idx = i
                break

        if molecules_idx == -1:
            print("Warning: Could not find [ molecules ] section in topology")
            return

        # Count metals by residue name
        from collections import Counter
        metal_counts = Counter(m['residue_name'] for m in metals)

        # Insert metal molecule entries after [ molecules ] header
        insert_idx = molecules_idx + 1
        # Skip comment lines
        while insert_idx < len(topol_lines) and topol_lines[insert_idx].strip().startswith(';'):
            insert_idx += 1

        for metal_name, count in metal_counts.items():
            topol_lines.insert(insert_idx, f"{metal_name:<20} {count}\n")
            insert_idx += 1
            print(f"  Added {metal_name}: {count}")

        with open(topol_file, 'w') as f:
            f.writelines(topol_lines)

        # Add metals to GRO file
        with open(gro_file, 'r') as f:
            gro_lines = f.readlines()

        # Parse current atom count
        atom_count = int(gro_lines[1].strip())

        # Convert metal coordinates to GRO format and append before box line
        metal_gro_lines = []
        for i, metal in enumerate(metals, start=1):
            # GRO format: resnr resname atomname atomnr x y z
            resnum = metal['resnum']
            resname = metal['residue_name']
            atomname = metal['atom_name']
            atom_num = atom_count + i
            x, y, z = metal['coords']
            # Convert Angstroms to nm
            x_nm, y_nm, z_nm = x/10.0, y/10.0, z/10.0

            # GRO format line
            line = f"{int(resnum):5d}{resname:5s}{atomname:>5s}{atom_num:5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n"
            metal_gro_lines.append(line)

        # Update atom count
        gro_lines[1] = f"{atom_count + len(metals):5d}\n"

        # Insert metal lines before box line (last line)
        gro_lines = gro_lines[:-1] + metal_gro_lines + [gro_lines[-1]]

        with open(gro_file, 'w') as f:
            f.writelines(gro_lines)

        print(f"Successfully added {len(metals)} metal ion(s) to system")


if __name__ == '__main__':
    print("This script is intended to be used as a module within the PRISM package.")
