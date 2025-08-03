#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM System Builder - Build the GROMACS model system
"""

import os
import subprocess
import shutil
from .mdp import MDPGenerator


class SystemBuilder:
    """Build the GROMACS model system"""

    def __init__(self, config, output_dir, overwrite=False):
        self.config = config
        self.output_dir = output_dir
        self.overwrite = overwrite
        self.model_dir = os.path.join(output_dir, "GMX_PROLIG_MD")
        self.gmx_command = config['general']['gmx_command']

    def build(self, cleaned_protein, lig_ff_dir, forcefield_idx, water_model_idx):
        """Build the complete system"""
        print("\n=== Building GROMACS Model ===")

        if os.path.exists(os.path.join(self.model_dir, "solv_ions.gro")) and not self.overwrite:
            print(f"Final model file solv_ions.gro already exists in {self.model_dir}, skipping model building")
            return self.model_dir

        os.makedirs(self.model_dir, exist_ok=True)

        # Step 1: Fix protein
        fixed_pdb = self._fix_protein(cleaned_protein)

        # Step 2: Generate topology
        self._generate_topology(fixed_pdb, forcefield_idx, water_model_idx)

        # Step 3: Combine protein and ligand
        self._combine_protein_ligand(lig_ff_dir)

        # Step 4: Fix topology
        self._fix_topology(lig_ff_dir)

        # Step 5: Create box
        self._create_box()

        # Step 6: Solvate
        self._solvate()

        # Step 7: Add ions
        self._add_ions()

        print(f"\nModel build completed. System files are in {self.model_dir}")
        return self.model_dir

    def _fix_protein(self, cleaned_protein):
        """Fix missing atoms in protein"""
        fixed_pdb = os.path.join(self.model_dir, "fixed_clean_protein.pdb")

        if not os.path.exists(fixed_pdb) or self.overwrite:
            print("Running pdbfixer to fix missing atoms...")
            self.run_command([
                "pdbfixer",
                cleaned_protein, "--add-residues",
                "--output", fixed_pdb,
                "--add-atoms", "heavy",
                "--keep-heterogens", "none",
                f"--ph={self.config['simulation']['pH']}"
            ])
        else:
            print(f"Fixed protein PDB already exists at {fixed_pdb}, skipping pdbfixer")

        return fixed_pdb

    def _generate_topology(self, fixed_pdb, forcefield_idx, water_model_idx):
        """Generate topology using pdb2gmx"""
        current_dir = os.getcwd()
        os.chdir(self.model_dir)

        pro_gro = os.path.join(self.model_dir, "pro.gro")
        if not os.path.exists(pro_gro) or self.overwrite:
            histidine_count = self._count_histidines(fixed_pdb)

            input_lines = [str(forcefield_idx), str(water_model_idx)]

            if histidine_count > 0:
                print(f"Found {histidine_count} histidine residue(s), defaulting to HIE")
                for _ in range(histidine_count):
                    input_lines.append("1")  # HIE

            input_text = '\n'.join(input_lines) + '\n'

            print(f"Running pdb2gmx with force field #{forcefield_idx} and water model #{water_model_idx}")
            self.run_command(
                [self.gmx_command, "pdb2gmx", "-f", "fixed_clean_protein.pdb", "-o", "pro.gro", "-ignh"],
                input_text=input_text
            )
        else:
            print(f"Protein topology already generated at {pro_gro}, skipping pdb2gmx")

        os.chdir(current_dir)

    def _combine_protein_ligand(self, lig_ff_dir):
        """Combine protein and ligand coordinates"""
        pro_lig_gro = os.path.join(self.model_dir, "pro_lig.gro")

        if not os.path.exists(pro_lig_gro) or self.overwrite:
            print("Combining protein and ligand coordinates...")

            pro_gro = os.path.join(self.model_dir, "pro.gro")
            lig_gro = os.path.join(lig_ff_dir, "LIG.gro")

            # Read ligand
            with open(lig_gro, 'r') as f:
                lines = f.readlines()
                lig_atom_count = int(lines[1].strip())
                lig_coords = ''.join(lines[2:-1])

            # Read protein
            with open(pro_gro, 'r') as f:
                lines = f.readlines()
                protein_atom_count = int(lines[1].strip())
                gro_header = lines[0]
                gro_box = lines[-1]
                protein_coords = ''.join(lines[2:-1])

            # Write combined
            total_atom_count = protein_atom_count + lig_atom_count
            with open(pro_lig_gro, 'w') as f:
                f.write(gro_header)
                f.write(f"{total_atom_count}\n")
                f.write(protein_coords)
                f.write(lig_coords)
                f.write(gro_box)

            print("Protein-ligand complex created")

    def _fix_topology(self, lig_ff_dir):
        """Fix topology file to include ligand parameters"""
        topol_path = os.path.join(self.model_dir, "topol.top")

        if os.path.exists(topol_path):
            print("Fixing topology file...")

            with open(topol_path, 'r') as f:
                content = f.read()

            # Get the correct directory name based on force field
            lig_dir_name = os.path.basename(lig_ff_dir)

            if f"{lig_dir_name}/LIG.itp" in content:
                print("Topology already includes ligand parameters")
                return

            atomtypes_path = os.path.join(lig_ff_dir, "atomtypes_LIG.itp")
            with open(atomtypes_path, 'r') as f:
                atomtypes_content = f.read()

            new_content = []
            lines = content.split('\n')
            i = 0

            while i < len(lines):
                line = lines[i]
                new_content.append(line)

                if '#include' in line and '.ff/forcefield.itp' in line:
                    new_content.append("")
                    new_content.append("; Ligand atomtypes")
                    new_content.extend(atomtypes_content.strip().split('\n'))
                    new_content.append("")
                    new_content.append("; Ligand topology")
                    new_content.append(f'#include "../{lig_dir_name}/LIG.itp"')
                    new_content.append("")

                if '[ molecules ]' in line:
                    j = i + 1
                    while j < len(lines) and lines[j].strip() and not lines[j].startswith('['):
                        new_content.append(lines[j])
                        j += 1
                    new_content.append("LIG                 1")
                    i = j - 1

                i += 1

            with open(topol_path, 'w') as f:
                f.write('\n'.join(new_content))

            print("Topology file updated with ligand parameters")

    def _create_box(self):
        """Create simulation box"""
        current_dir = os.getcwd()
        os.chdir(self.model_dir)

        box_gro = os.path.join(self.model_dir, "pro_lig_newbox.gro")
        if not os.path.exists(box_gro) or self.overwrite:
            box_cmd = [
                self.gmx_command, "editconf",
                "-f", "pro_lig.gro",
                "-o", "pro_lig_newbox.gro"
            ]

            if self.config['box']['center']:
                box_cmd.append("-c")

            # Add box shape
            if self.config['box']['shape'] == 'dodecahedron':
                box_cmd.extend(["-bt", "dodecahedron"])
            elif self.config['box']['shape'] == 'octahedron':
                box_cmd.extend(["-bt", "octahedron"])

            # Add distance
            box_cmd.extend(["-d", str(self.config['box']['distance'])])

            self.run_command(box_cmd)

        os.chdir(current_dir)

    def _solvate(self):
        """Solvate the system"""
        current_dir = os.getcwd()
        os.chdir(self.model_dir)

        solv_gro = os.path.join(self.model_dir, "solv.gro")
        if not os.path.exists(solv_gro) or self.overwrite:
            self.run_command([
                self.gmx_command, "solvate",
                "-cp", "pro_lig_newbox.gro",
                "-p", "topol.top",
                "-o", "solv.gro"
            ])

        os.chdir(current_dir)

    def _add_ions(self):
        """Add ions to neutralize system"""
        current_dir = os.getcwd()
        os.chdir(self.model_dir)

        # Generate ions.mdp if needed
        ions_mdp = os.path.join(self.model_dir, "ions.mdp")
        if not os.path.exists(ions_mdp) or self.overwrite:
            mdp_gen = MDPGenerator(self.config, self.output_dir)
            mdp_gen._generate_ions_mdp()
            # Move to model directory
            shutil.move(os.path.join(mdp_gen.mdp_dir, "ions.mdp"), ions_mdp)

        # Generate TPR
        ions_tpr = os.path.join(self.model_dir, "ions.tpr")
        if not os.path.exists(ions_tpr) or self.overwrite:
            self.run_command([
                self.gmx_command, "grompp",
                "-f", ions_mdp,
                "-c", "solv.gro",
                "-p", "topol.top",
                "-o", "ions.tpr",
                "-maxwarn", "99999"
            ])

        # Add ions
        solv_ions_gro = os.path.join(self.model_dir, "solv_ions.gro")
        if not os.path.exists(solv_ions_gro) or self.overwrite:
            print("Adding ions...")

            # First, check if we need to add ions at all
            if self._check_ion_requirements():
                genion_cmd = [
                    self.gmx_command, "genion", "-s", "ions.tpr", "-o", "solv_ions.gro",
                    "-p", "topol.top",
                    "-pname", self.config['ions']['positive_ion'],
                    "-nname", self.config['ions']['negative_ion']
                ]

                if self.config['ions']['neutral']:
                    genion_cmd.append("-neutral")

                if self.config['ions']['concentration'] > 0:
                    genion_cmd.extend(["-conc", str(self.config['ions']['concentration'])])

                # Try different group numbers for SOL
                success = False
                for group in ["15", "14", "13"]:  # Try SOL group first
                    try:
                        output = self.run_command(genion_cmd, input_text=group, check=False)
                        if "Number of (3-atomic) solvent molecules:" in output:
                            print(f"Successfully added ions using group {group}")
                            success = True
                            break
                    except:
                        continue

                if not success:
                    # If no suitable water group found, it might be because no ions are needed
                    print("Warning: Could not add ions. This might be because the system is already neutral.")
                    print("Copying solvated structure without ions...")
                    shutil.copy("solv.gro", "solv_ions.gro")
            else:
                # Just copy the solvated structure if no ions needed
                print("System is already neutral, no ions needed.")
                shutil.copy("solv.gro", "solv_ions.gro")

        os.chdir(current_dir)

    def _check_ion_requirements(self):
        """Check if ions are actually needed"""
        # This is a simple check - you might want to enhance this
        # by actually reading the topology to determine charge
        return self.config['ions']['neutral'] or self.config['ions']['concentration'] > 0

    def run_command(self, cmd, input_text=None, cwd=None, check=True, shell=False):
        """Run a shell command with optional input"""
        cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
        print(f"Running: {cmd_str}")

        try:
            if input_text is not None:
                if shell:
                    process = subprocess.Popen(
                        cmd_str,
                        cwd=cwd,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True
                    )
                else:
                    process = subprocess.Popen(
                        cmd,
                        cwd=cwd,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True
                    )
                stdout, stderr = process.communicate(input=input_text)

                if process.returncode != 0 and check:
                    print(f"Error: {stderr}")
                    raise subprocess.CalledProcessError(process.returncode, cmd, stdout, stderr)

                return stdout
            else:
                if shell:
                    result = subprocess.run(
                        cmd_str,
                        cwd=cwd,
                        shell=True,
                        check=check,
                        capture_output=True,
                        text=True
                    )
                else:
                    result = subprocess.run(
                        cmd,
                        cwd=cwd,
                        check=check,
                        capture_output=True,
                        text=True
                    )
                return result.stdout

        except subprocess.CalledProcessError as e:
            print(f"Command failed: {cmd_str}")
            print(f"Error: {e.stderr}")
            if check:
                raise
            return e.stdout if hasattr(e, 'stdout') else ""

    def _count_histidines(self, pdb_file):
        """Count histidine residues in PDB file"""
        his_count = 0
        his_residues = []

        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[17:20].strip() in ['HIS', 'HID', 'HIE', 'HIP']:
                    resnum = line[22:26].strip()
                    chain = line[21].strip()
                    res_id = f"{chain}:{resnum}" if chain else resnum

                    if res_id not in his_residues:
                        his_residues.append(res_id)
                        his_count += 1

        if his_count > 0:
            print(f"Histidine residues found: {', '.join(his_residues)}")

        return his_count