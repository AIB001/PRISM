#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder PMF Mixin - PMF (Steered MD / Umbrella Sampling) methods

This mixin class provides all PMF-related methods for PRISMBuilder.
Extracted from the original monolithic builder.py for maintainability.
"""

import os

from ..utils.colors import (
    print_header,
    print_subheader,
    print_step,
    print_success,
    print_error,
    print_warning,
    print_info,
    path,
    number,
)


class PMFBuilderMixin:
    """Mixin class providing PMF workflow methods for PRISMBuilder."""

    def run_pmf(self):
        """Run the PMF workflow for steered MD system building"""
        print_header("PRISM PMF Builder Workflow")

        try:
            # Auto-compute box Z extension from pull distance if not user-specified
            pull_distance = self.config.get("pmf", {}).get("pull_distance", 3.6)
            if self.box_extension == (0.0, 0.0, 2.0):  # default, not user-specified
                z_ext = pull_distance + 1.0  # pull distance + 1 nm buffer
                self.box_extension = (0.0, 0.0, z_ext)
                self.system_builder.box_extension = self.box_extension
                print(
                    f"  Auto box Z extension: {z_ext:.1f} nm " f"(pull distance {pull_distance:.1f} nm + 1.0 nm buffer)"
                )

            # Save configuration for reference
            self.save_config()

            total_steps = 9

            # Step 1: Generate ligand force field
            print_step(1, total_steps, f"Generating ligand force field ({self.ligand_forcefield.upper()})")
            self.generate_ligand_forcefield()
            print_success(f"Ligand force field generated")

            # Step 2: Clean protein (minimal cleaning for PMF)
            print_step(2, total_steps, "Cleaning protein structure")
            cleaned_protein = self.clean_protein()
            print_success("Protein structure cleaned")

            # Step 3: Align structures for PMF (rotate to align pull vector with Z-axis)
            print_step(3, total_steps, "Aligning structures for PMF")
            aligned_protein, aligned_ligand = self.align_for_pmf(cleaned_protein)
            print_success("Structures aligned for PMF")

            # Step 3.5: Update LIG.gro coordinates with aligned ligand coordinates
            print_info("Updating ligand force field coordinates with aligned positions...")
            self.update_ligand_gro_coordinates(aligned_ligand)
            print_success("Ligand coordinates updated")

            # Step 4: Build model with extended Z-box
            print_step(4, total_steps, "Building GROMACS system (PMF mode)")
            model_dir = self.build_model(aligned_protein)
            if not model_dir:
                raise RuntimeError("Failed to build model")
            print_success("GROMACS PMF system built")

            # Step 5: Generate custom index file for SMD
            print_step(5, total_steps, "Generating SMD index file")
            index_path, pbcatom, n_near_atoms, n_near_res, n_ca = self.generate_smd_index()
            print_success(f"SMD index generated: {n_near_atoms} atoms in reference group, {n_ca} CA atoms frozen")

            # Step 6: Generate MDP files (including SMD parameters)
            print_step(6, total_steps, "Generating MD parameter files")
            self.generate_mdp_files()
            # Generate SMD-specific MDP file with custom groups
            pull_rate = self.config.get("pmf", {}).get("pull_rate", 0.01)
            pull_k = self.config.get("pmf", {}).get("pull_k", 1000.0)
            pull_distance = self.config.get("pmf", {}).get("pull_distance", 3.6)
            self.mdp_generator.generate_smd_mdp(
                pull_rate=pull_rate,
                pull_k=pull_k,
                pull_distance=pull_distance,
                ref_group="Protein_near_LIG",
                pull_group="LIG",
                pbcatom=pbcatom,
                freeze_group="CA_freeze",
            )
            print_success("MDP files generated (including smd.mdp)")

            # Step 7: Generate SMD run script
            print_step(7, total_steps, "Generating SMD run script")
            self.generate_smd_plot_script()
            script_path = self.generate_smd_script()
            print_success("SMD script generated")

            # Step 8: Generate umbrella sampling MDP
            print_step(8, total_steps, "Generating umbrella sampling MDP")
            umbrella_time_ns = self.config.get("pmf", {}).get("umbrella_time_ns", 10.0)
            self.mdp_generator.generate_umbrella_mdp(
                umbrella_time_ns=umbrella_time_ns,
                pull_k=pull_k,
                ref_group="Protein_near_LIG",
                pull_group="LIG",
                pbcatom=pbcatom,
                freeze_group="CA_freeze",
            )
            print_success("Umbrella sampling MDP generated")

            # Step 9: Generate umbrella sampling scripts
            print_step(9, total_steps, "Generating umbrella sampling scripts")
            self.generate_pmf_plot_script()
            self.generate_umbrella_script()
            print_success("Umbrella sampling scripts generated")

            print_header("PMF Workflow Complete!")
            print(f"\nOutput files are in: {path(self.output_dir)}")
            print(f"SMD system files are in: {path(os.path.join(self.output_dir, 'GMX_PROLIG_PMF'))}")
            print(f"MDP files are in: {path(self.mdp_generator.mdp_dir)}")
            print(f"Configuration saved in: {path(os.path.join(self.output_dir, 'prism_config.yaml'))}")
            print(f"\nProtein force field used: {number(self.forcefield['name'])}")
            print(f"Ligand force field used: {number(self.ligand_forcefield.upper())}")
            print(f"Water model used: {number(self.water_model['name'])}")
            print(f"Box extension (Z): {number(f'{self.box_extension[2]:.1f} nm')}")

            if script_path:
                gmx_smd_dir = os.path.join(self.output_dir, "GMX_PROLIG_PMF")
                print_header("Ready to Run PMF Simulations!")
                print(f"\nTo run the complete PMF workflow:")
                print(f"  1. Navigate to the SMD directory:")
                print(f"     {path(f'cd {gmx_smd_dir}')}")
                print(f"  2. Run SMD (EM -> NVT -> NPT -> pulling):")
                print(f"     {number('bash smd_run.sh')}")
                print(f"  3. Run umbrella sampling + WHAM:")
                print(f"     {number('bash umbrella_run.sh')}")

            return self.output_dir

        except Exception as e:
            print_error(f"PMF workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise

    def align_for_pmf(self, cleaned_protein):
        """
        Align protein and ligand structures for PMF calculations.

        Rotates the complex so the pull vector aligns with the Z-axis.

        Returns:
        --------
        tuple : (aligned_protein_path, aligned_ligand_path)
        """
        print_subheader("Aligning Structures for PMF")

        from ..pmf.alignment import PMFAligner

        # Create alignment work directory
        align_dir = os.path.join(self.output_dir, "PMF_ALIGNMENT")
        os.makedirs(align_dir, exist_ok=True)

        # Initialize aligner
        aligner = PMFAligner(pocket_cutoff=4.0, verbose=True)

        # Align protein and first ligand
        # For multiple ligands, we align based on the first one
        ligand_path = self.ligand_paths[0]

        alignment_results = aligner.align_for_pmf(
            protein_pdb=cleaned_protein, ligand_file=ligand_path, output_dir=align_dir, pullvec=self.pullvec
        )

        aligned_protein = alignment_results["aligned_protein"]
        aligned_ligand = alignment_results["aligned_ligand"]

        # Update ligand paths to use aligned ligand
        # The force field was already generated from original ligand, coordinates will be updated
        self.aligned_ligand_path = aligned_ligand

        print(f"  Pull vector aligned to Z-axis")
        print(f"  Aligned protein: {path(aligned_protein)}")
        print(f"  Aligned ligand: {path(aligned_ligand)}")

        return aligned_protein, aligned_ligand

    def update_ligand_gro_coordinates(self, aligned_ligand_file: str):
        """
        Update LIG.gro coordinates with aligned ligand coordinates.

        After alignment, the LIG.gro file still has the original coordinates.
        This method reads the aligned ligand MOL2/SDF and updates the GRO file
        with the new coordinates (converting from Angstroms to nm).

        Parameters:
        -----------
        aligned_ligand_file : str
            Path to the aligned ligand file (MOL2 or SDF)
        """
        from pathlib import Path as Pth

        # Parse aligned ligand to get new coordinates
        suffix = Pth(aligned_ligand_file).suffix.lower()

        aligned_coords = []  # List of (x, y, z) in Angstroms

        if suffix == ".mol2":
            with open(aligned_ligand_file, "r") as f:
                in_atom_section = False
                for line in f:
                    if line.startswith("@<TRIPOS>ATOM"):
                        in_atom_section = True
                        continue
                    elif line.startswith("@<TRIPOS>"):
                        in_atom_section = False
                        continue
                    if in_atom_section and line.strip():
                        parts = line.split()
                        if len(parts) >= 5:
                            try:
                                x = float(parts[2])
                                y = float(parts[3])
                                z = float(parts[4])
                                aligned_coords.append((x, y, z))
                            except ValueError:
                                continue
        elif suffix in [".sdf", ".sd"]:
            with open(aligned_ligand_file, "r") as f:
                lines = f.readlines()
            if len(lines) >= 4:
                try:
                    num_atoms = int(lines[3][:3].strip())
                    for i in range(4, 4 + num_atoms):
                        if i < len(lines):
                            line = lines[i]
                            x = float(line[0:10].strip())
                            y = float(line[10:20].strip())
                            z = float(line[20:30].strip())
                            aligned_coords.append((x, y, z))
                except (ValueError, IndexError):
                    pass

        if not aligned_coords:
            print_warning(f"Could not read coordinates from aligned ligand file: {aligned_ligand_file}")
            return

        print_info(f"  Read {len(aligned_coords)} atom coordinates from aligned ligand")

        # Update each LIG.gro file in the force field directories
        for ff_dir in self.lig_ff_dirs:
            gro_file = os.path.join(ff_dir, "LIG.gro")

            if not os.path.exists(gro_file):
                print_warning(f"  LIG.gro not found in {ff_dir}, skipping")
                continue

            # Read original GRO file
            with open(gro_file, "r") as f:
                gro_lines = f.readlines()

            if len(gro_lines) < 3:
                print_warning(f"  Invalid GRO file: {gro_file}")
                continue

            # GRO format:
            # Line 0: title
            # Line 1: number of atoms
            # Lines 2 to N+1: atom lines
            # Last line: box vectors

            title = gro_lines[0]
            num_atoms = int(gro_lines[1].strip())
            box_line = gro_lines[-1]

            if num_atoms != len(aligned_coords):
                print_warning(f"  Atom count mismatch: GRO has {num_atoms}, aligned has {len(aligned_coords)}")
                # Try to continue anyway with min atoms
                num_atoms = min(num_atoms, len(aligned_coords))

            # Update atom coordinates (convert Angstroms to nm)
            new_gro_lines = [title, gro_lines[1]]

            for i in range(num_atoms):
                old_line = gro_lines[2 + i]
                x_nm = aligned_coords[i][0] / 10.0  # Angstroms to nm
                y_nm = aligned_coords[i][1] / 10.0
                z_nm = aligned_coords[i][2] / 10.0

                # GRO format: residue number (5), residue name (5), atom name (5),
                # atom number (5), x (8.3f), y (8.3f), z (8.3f), [vx, vy, vz]
                # Total: first 20 chars are fixed, then coordinates

                # Keep the first 20 characters (residue info, atom name, atom number)
                prefix = old_line[:20]

                # Format new coordinates (8.3f format = 8 chars, 3 decimal places)
                new_line = f"{prefix}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n"
                new_gro_lines.append(new_line)

            # Add box line
            new_gro_lines.append(box_line)

            # Write updated GRO file
            with open(gro_file, "w") as f:
                f.writelines(new_gro_lines)

            print_success(f"  Updated coordinates in {path(gro_file)}")

    def generate_smd_index(self):
        """
        Generate custom index file for SMD simulations.

        Creates groups:
        - Protein_near_LIG: protein residues within 5A of ligand (reference/anchor group)
        - CA_freeze: all C-alpha atoms (freeze group during SMD)

        Returns:
        --------
        tuple: (index_path, pbcatom, n_near_atoms, n_near_residues, n_ca_atoms)
        """
        import numpy as np
        import subprocess

        gmx_smd_dir = os.path.join(self.output_dir, "GMX_PROLIG_PMF")
        gro_file = os.path.join(gmx_smd_dir, "solv_ions.gro")
        index_path = os.path.join(gmx_smd_dir, "index.ndx")

        # Step 1: Generate default index with gmx make_ndx
        print("  Generating default index groups...")
        make_ndx_cmd = [self.gromacs_env.gmx_command, "make_ndx", "-f", gro_file, "-o", index_path]
        process = subprocess.Popen(
            make_ndx_cmd,
            cwd=gmx_smd_dir,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        process.communicate(input="q\n")

        # Step 2: Parse GRO file
        print("  Parsing system coordinates...")
        with open(gro_file, "r") as f:
            lines = f.readlines()

        n_atoms = int(lines[1].strip())

        lig_indices = []
        lig_coords = []
        protein_residues = {}  # (resnum, resname) -> [(atom_idx, x, y, z), ...]
        ca_indices = []

        non_protein_resnames = {"SOL", "NA", "CL", "WAT", "HOH", "K", "MG", "Na+", "Cl-", "LIG"}

        for i in range(2, 2 + n_atoms):
            line = lines[i]
            resnum = int(line[0:5].strip())
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atom_idx = i - 1  # 1-based GROMACS numbering
            x = float(line[20:28]) * 10.0  # nm -> Angstrom
            y = float(line[28:36]) * 10.0
            z = float(line[36:44]) * 10.0

            if resname == "LIG":
                lig_indices.append(atom_idx)
                lig_coords.append([x, y, z])
            elif resname not in non_protein_resnames:
                reskey = (resnum, resname)
                if reskey not in protein_residues:
                    protein_residues[reskey] = []
                protein_residues[reskey].append((atom_idx, x, y, z))
                if atomname == "CA":
                    ca_indices.append(atom_idx)

        # Step 3: Find protein residues near ligand
        lig_coords_arr = np.array(lig_coords)

        # 5A cutoff for pull reference group
        print("  Finding protein residues near ligand (5 A cutoff for pull group)...")
        near_lig_reskeys_5A = set()
        for reskey, atoms in protein_residues.items():
            for atom_idx, x, y, z in atoms:
                diffs = lig_coords_arr - np.array([x, y, z])
                dists = np.linalg.norm(diffs, axis=1)
                if np.min(dists) <= 5.0:
                    near_lig_reskeys_5A.add(reskey)
                    break

        near_lig_indices = []
        for reskey in sorted(near_lig_reskeys_5A):
            for atom_idx, x, y, z in protein_residues[reskey]:
                near_lig_indices.append(atom_idx)
        near_lig_indices.sort()

        # 8A cutoff for freeze exclusion (pocket residues get unfrozen CA)
        print("  Finding pocket residues to unfreeze (8 A cutoff)...")
        pocket_reskeys_8A = set()
        for reskey, atoms in protein_residues.items():
            for atom_idx, x, y, z in atoms:
                diffs = lig_coords_arr - np.array([x, y, z])
                dists = np.linalg.norm(diffs, axis=1)
                if np.min(dists) <= 8.0:
                    pocket_reskeys_8A.add(reskey)
                    break

        # Build CA atom -> reskey mapping for exclusion
        ca_reskey_map = {}  # atom_idx -> reskey
        for reskey, atoms in protein_residues.items():
            for atom_idx, x, y, z in atoms:
                if atom_idx in ca_indices:
                    ca_reskey_map[atom_idx] = reskey

        # Exclude pocket CA atoms from freeze group
        ca_freeze_indices = [idx for idx in ca_indices if ca_reskey_map.get(idx) not in pocket_reskeys_8A]
        ca_unfrozen = len(ca_indices) - len(ca_freeze_indices)

        # Step 4: Find pbcatom (atom closest to center of Protein_near_LIG group)
        near_coords_list = []
        for reskey in near_lig_reskeys_5A:
            for atom_idx, x, y, z in protein_residues[reskey]:
                near_coords_list.append([atom_idx, x, y, z])
        near_coords_arr = np.array(near_coords_list)
        center = np.mean(near_coords_arr[:, 1:], axis=0)
        dists_to_center = np.linalg.norm(near_coords_arr[:, 1:] - center, axis=1)
        pbcatom = int(near_coords_arr[np.argmin(dists_to_center), 0])

        # Step 5: Append custom groups to index file
        print("  Writing custom groups to index file...")
        with open(index_path, "a") as f:

            def write_group(name, indices):
                f.write(f"\n[ {name} ]\n")
                for j, idx in enumerate(indices):
                    f.write(f"{idx:>6d}")
                    if (j + 1) % 15 == 0:
                        f.write("\n")
                if len(indices) % 15 != 0:
                    f.write("\n")

            write_group("Protein_near_LIG", near_lig_indices)
            write_group("CA_freeze", ca_freeze_indices)

        # Print summary
        near_res_info = sorted([f"{rn}{rnum}" for rnum, rn in near_lig_reskeys_5A])
        pocket_res_info = sorted([f"{rn}{rnum}" for rnum, rn in pocket_reskeys_8A])
        print(f"  Protein_near_LIG: {len(near_lig_indices)} atoms from {len(near_lig_reskeys_5A)} residues")
        print(f"    Residues: {', '.join(near_res_info)}")
        print(f"  Pocket residues (8 A, unfrozen): {len(pocket_reskeys_8A)} residues, {ca_unfrozen} CA atoms released")
        print(f"    Residues: {', '.join(pocket_res_info)}")
        print(f"  CA_freeze: {len(ca_freeze_indices)} C-alpha atoms (total {len(ca_indices)} - {ca_unfrozen} pocket)")
        print(f"  PBC reference atom: {pbcatom}")
        print(f"  Index file: {index_path}")

        return index_path, pbcatom, len(near_lig_indices), len(near_lig_reskeys_5A), len(ca_freeze_indices)

    def generate_smd_plot_script(self):
        """Generate Python script for plotting SMD force vs position"""
        gmx_smd_dir = os.path.join(self.output_dir, "GMX_PROLIG_PMF")
        smd_dir = os.path.join(gmx_smd_dir, "smd")
        os.makedirs(smd_dir, exist_ok=True)
        script_path = os.path.join(smd_dir, "plot_smd.py")

        script_content = '''#!/usr/bin/env python3
"""Plot SMD Force vs Position from GROMACS pull output files."""

import sys
import os
import numpy as np

def parse_xvg(filepath):
    """Parse a GROMACS .xvg file, return arrays of (time, values)."""
    times, values = [], []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                values.append(float(parts[1]))
    return np.array(times), np.array(values)

def main():
    if len(sys.argv) < 3:
        print("Usage: python plot_smd.py pullf.xvg pullx.xvg [output_dir]")
        sys.exit(1)

    pullf_file = sys.argv[1]
    pullx_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(pullf_file)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    time_f, force = parse_xvg(pullf_file)
    time_x, position = parse_xvg(pullx_file)

    # Align by time (should be identical, but be safe)
    if len(time_f) != len(time_x) or not np.allclose(time_f, time_x):
        common_times = np.intersect1d(time_f, time_x)
        mask_f = np.isin(time_f, common_times)
        mask_x = np.isin(time_x, common_times)
        force = force[mask_f]
        position = position[mask_x]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(position, force, color='#2166AC', linewidth=0.6, alpha=0.8)
    ax.set_xlabel('Position (nm)', fontsize=24)
    ax.set_ylabel('Force (kJ/mol/nm)', fontsize=24)
    ax.set_title('SMD Force vs Position', fontsize=28)
    ax.tick_params(axis='both', labelsize=22)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    output_path = os.path.join(output_dir, 'force_vs_position.png')
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_path}")

if __name__ == '__main__':
    main()
'''
        with open(script_path, "w") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        return script_path

    def generate_pmf_plot_script(self):
        """Generate Python script for plotting PMF profile from WHAM output"""
        gmx_smd_dir = os.path.join(self.output_dir, "GMX_PROLIG_PMF")
        umbrella_dir = os.path.join(gmx_smd_dir, "umbrella")
        os.makedirs(umbrella_dir, exist_ok=True)
        script_path = os.path.join(umbrella_dir, "plot_pmf.py")

        script_content = '''#!/usr/bin/env python3
"""Plot PMF profile from GROMACS WHAM output."""

import sys
import os
import numpy as np

def parse_xvg(filepath):
    """Parse a GROMACS .xvg file, return arrays of columns."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                data.append([float(x) for x in parts])
    return np.array(data)

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_pmf.py pmf.xvg [pmferror.xvg] [output_dir]")
        sys.exit(1)

    pmf_file = sys.argv[1]
    error_file = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].endswith('/') else None
    output_dir = sys.argv[-1] if sys.argv[-1].endswith('/') or os.path.isdir(sys.argv[-1]) else os.path.dirname(pmf_file)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Parse PMF data (columns: distance, PMF)
    pmf_data = parse_xvg(pmf_file)
    distance = pmf_data[:, 0]
    pmf = pmf_data[:, 1]

    # Parse error data if available
    has_error = False
    if error_file and os.path.isfile(error_file):
        try:
            err_data = parse_xvg(error_file)
            pmf_error = err_data[:, 1]
            has_error = True
        except Exception:
            pass

    # Find binding free energy (minimum PMF value)
    min_idx = np.argmin(pmf)
    min_pmf = pmf[min_idx]
    min_dist = distance[min_idx]

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(distance, pmf, color='#2166AC', linewidth=2.0, label='PMF')
    if has_error:
        ax.fill_between(distance, pmf - pmf_error, pmf + pmf_error,
                        color='#2166AC', alpha=0.2, label='Bootstrap error')

    # Mark minimum
    ax.axhline(y=min_pmf, color='#B2182B', linestyle='--', alpha=0.5, linewidth=1.0)
    ax.plot(min_dist, min_pmf, 'o', color='#B2182B', markersize=8)
    ax.annotate(f'  Min: {min_pmf:.2f} kcal/mol\\n  at {min_dist:.2f} nm',
                xy=(min_dist, min_pmf), fontsize=16,
                xytext=(min_dist + 0.3, min_pmf + 1.0))

    ax.set_xlabel('Distance (nm)', fontsize=24)
    ax.set_ylabel('PMF (kcal/mol)', fontsize=24)
    ax.set_title('Potential of Mean Force', fontsize=28)
    ax.tick_params(axis='both', labelsize=22)
    ax.grid(True, alpha=0.3)
    if has_error:
        ax.legend(fontsize=18, loc='best')
    fig.tight_layout()

    output_path = os.path.join(output_dir, 'pmf_profile.png')
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"PMF profile saved: {output_path}")
    print(f"Binding free energy: {min_pmf:.2f} kcal/mol at {min_dist:.2f} nm")

if __name__ == '__main__':
    main()
'''
        with open(script_path, "w") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        return script_path

    def generate_smd_script(self):
        """Generate smd_run.sh script for SMD simulation execution"""
        gmx_smd_dir = os.path.join(self.output_dir, "GMX_PROLIG_PMF")
        if not os.path.exists(gmx_smd_dir):
            print_warning("GMX_PROLIG_PMF directory not found, skipping smd_run.sh generation")
            return None

        script_path = os.path.join(gmx_smd_dir, "smd_run.sh")

        # Get pull parameters
        pull_rate = self.config.get("pmf", {}).get("pull_rate", 0.01)
        pull_k = self.config.get("pmf", {}).get("pull_k", 1000.0)
        pull_distance = self.config.get("pmf", {}).get("pull_distance", 3.6)

        script_content = f"""#!/bin/bash

######################################################
# PRISM SMD SIMULATION SCRIPT
# For Steered Molecular Dynamics (PMF calculations)
######################################################

set -e
export GMX_MAXBACKUP=-1

# Hardware settings - adjust as needed
NCPU=10
GPU_ID=0

echo "=============================================="
echo "     PRISM SMD Simulation Workflow"
echo "=============================================="
echo "Pull rate: {pull_rate} nm/ps"
echo "Pull constant: {pull_k} kJ/mol/nm^2"
echo "Pull distance: {pull_distance} nm"
echo ""

# Energy Minimization (EM)
echo ">>> Step 1: Energy Minimization"
mkdir -p em
if [ -f ./em/em.gro ]; then
    echo "EM already completed, skipping..."
elif [ -f ./em/em.tpr ]; then
    echo "EM tpr file found, continuing from checkpoint..."
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp $NCPU -gpu_id $GPU_ID -v -cpi ./em/em.cpt
else
    echo "Starting EM from scratch..."
    gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o ./em/em.tpr -maxwarn 999
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp $NCPU -gpu_id $GPU_ID -v
fi

# NVT Equilibration
echo ""
echo ">>> Step 2: NVT Equilibration"
mkdir -p nvt
if [ -f ./nvt/nvt.gro ]; then
    echo "NVT already completed, skipping..."
elif [ -f ./nvt/nvt.tpr ]; then
    echo "NVT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp $NCPU -nb gpu -bonded gpu -pme gpu -gpu_id $GPU_ID -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v -cpi ./nvt/nvt.cpt
else
    echo "Starting NVT from scratch..."
    gmx grompp -f ../mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro -p topol.top -o ./nvt/nvt.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp $NCPU -nb gpu -bonded gpu -pme gpu -gpu_id $GPU_ID -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v
fi

# NPT Equilibration
echo ""
echo ">>> Step 3: NPT Equilibration"
mkdir -p npt
if [ -f ./npt/npt.gro ]; then
    echo "NPT already completed, skipping..."
elif [ -f ./npt/npt.tpr ]; then
    echo "NPT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp $NCPU -nb gpu -bonded gpu -pme gpu -gpu_id $GPU_ID -s ./npt/npt.tpr -deffnm ./npt/npt -v -cpi ./npt/npt.cpt
else
    echo "Starting NPT from scratch..."
    gmx grompp -f ../mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp $NCPU -nb gpu -bonded gpu -pme gpu -gpu_id $GPU_ID -s ./npt/npt.tpr -deffnm ./npt/npt -v
fi

# Steered MD (SMD) - Pulling along Z-axis
echo ""
echo ">>> Step 4: Steered Molecular Dynamics (SMD)"
mkdir -p smd
if [ -f ./smd/smd.gro ]; then
    echo "SMD already completed, skipping..."
elif [ -f ./smd/smd.tpr ]; then
    echo "SMD tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp $NCPU -nb gpu -bonded gpu -pme gpu -gpu_id $GPU_ID -s ./smd/smd.tpr -deffnm ./smd/smd -v -cpi ./smd/smd.cpt -pf ./smd/pullf.xvg -px ./smd/pullx.xvg
else
    echo "Starting SMD from scratch..."
    # Index file with custom pull groups (Protein_near_LIG, LIG, CA_freeze)
    # was pre-generated during PRISM build
    echo "Using pre-built index file with custom pull groups..."

    gmx grompp -f ../mdps/smd.mdp -c ./npt/npt.gro -r ./npt/npt.gro -t ./npt/npt.cpt -p topol.top -n index.ndx -o ./smd/smd.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp $NCPU -nb gpu -bonded gpu -pme gpu -gpu_id $GPU_ID -s ./smd/smd.tpr -deffnm ./smd/smd -v -pf ./smd/pullf.xvg -px ./smd/pullx.xvg
fi

echo ""
echo "=============================================="
echo "     SMD Simulation Complete!"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - Trajectory: ./smd/smd.xtc"
echo "  - Pull force: ./smd/pullf.xvg"
echo "  - Pull position: ./smd/pullx.xvg"
echo ""
echo "Next steps for PMF calculation:"
echo "  1. Extract frames for umbrella sampling windows"
echo "  2. Run umbrella sampling at each window"
echo "  3. Analyze with WHAM: gmx wham -it tpr-files.dat -if pullf-files.dat"

# ============================================
# Auto-plot Force vs Position
# ============================================
if [ -f ./smd/pullf.xvg ] && [ -f ./smd/pullx.xvg ]; then
    echo ""
    echo ">>> Plotting Force vs Position..."
    python3 ./smd/plot_smd.py ./smd/pullf.xvg ./smd/pullx.xvg ./smd/
    if [ $? -eq 0 ]; then
        echo "Plot saved to ./smd/force_vs_position.png"
    else
        echo "Warning: Plotting failed. You can run manually:"
        echo "  python3 ./smd/plot_smd.py ./smd/pullf.xvg ./smd/pullx.xvg ./smd/"
    fi
fi
"""

        with open(script_path, "w") as f:
            f.write(script_content)

        # Make script executable
        os.chmod(script_path, 0o755)

        print(f"SMD run script generated: {script_path}")
        return script_path

    def generate_umbrella_script(self):
        """Generate umbrella_run.sh script for umbrella sampling + WHAM workflow"""
        gmx_smd_dir = os.path.join(self.output_dir, "GMX_PROLIG_PMF")
        if not os.path.exists(gmx_smd_dir):
            print_warning("GMX_PROLIG_PMF directory not found, skipping umbrella_run.sh generation")
            return None

        script_path = os.path.join(gmx_smd_dir, "umbrella_run.sh")

        # Get umbrella parameters
        umbrella_spacing = self.config.get("pmf", {}).get("umbrella_spacing", 0.12)
        wham_begin = self.config.get("pmf", {}).get("wham_begin", 1000)
        wham_bootstrap = self.config.get("pmf", {}).get("wham_bootstrap", 200)

        script_content = f"""#!/bin/bash

######################################################
# PRISM UMBRELLA SAMPLING + WHAM SCRIPT
# Run after smd_run.sh completes
######################################################

set -e
export GMX_MAXBACKUP=-1

# ============================================
# Configurable parameters - adjust as needed
# ============================================
WINDOW_SPACING={umbrella_spacing}    # nm between umbrella windows
NCPU=10                # CPU threads per window
GPU_ID=0               # GPU device ID
WHAM_BEGIN={wham_begin}          # ps to discard for WHAM equilibration
WHAM_BOOTSTRAP={wham_bootstrap}      # number of bootstrap iterations

echo "=============================================="
echo "     PRISM Umbrella Sampling Workflow"
echo "=============================================="
echo "Window spacing: ${{WINDOW_SPACING}} nm"
echo "CPUs per window: ${{NCPU}}"
echo ""

# ============================================
# Phase 1: Select umbrella windows from SMD
# ============================================
echo ">>> Phase 1: Selecting umbrella windows from SMD trajectory"

if [ ! -f ./smd/pullx.xvg ]; then
    echo "ERROR: ./smd/pullx.xvg not found. Run smd_run.sh first!"
    exit 1
fi

mkdir -p umbrella

python3 << 'PYEOF'
import numpy as np
import sys

spacing = {umbrella_spacing}
times, positions = [], []
with open('./smd/pullx.xvg', 'r') as f:
    for line in f:
        if line.startswith(('#', '@')):
            continue
        parts = line.split()
        if len(parts) >= 2:
            times.append(float(parts[0]))
            positions.append(float(parts[1]))

times = np.array(times)
positions = np.array(positions)

# Select windows at uniform spacing along the pull coordinate
pos_min, pos_max = positions.min(), positions.max()
target_positions = np.arange(pos_min, pos_max, spacing)

selected = []
for target in target_positions:
    idx = np.argmin(np.abs(positions - target))
    selected.append((len(selected), times[idx], positions[idx]))

with open('./umbrella/windows.txt', 'w') as f:
    f.write("# window_id  time(ps)  position(nm)\\n")
    for wid, t, p in selected:
        f.write(f"{{wid}}  {{t:.1f}}  {{p:.4f}}\\n")

print(f"Selected {{len(selected)}} umbrella windows from {{pos_min:.3f}} to {{pos_max:.3f}} nm")
PYEOF

N_WINDOWS=$(grep -v '^#' ./umbrella/windows.txt | wc -l)
echo "Selected $N_WINDOWS umbrella windows"
echo ""

# ============================================
# Phase 2: Extract frames from SMD trajectory
# ============================================
echo ">>> Phase 2: Extracting frames for each umbrella window"

N=0
while read -r WID TIME POS; do
    # Skip comment lines
    [[ "$WID" == \\#* ]] && continue

    CONF="./umbrella/conf${{N}}.gro"
    if [ -f "$CONF" ]; then
        echo "  Window $N: conf already exists, skipping..."
    else
        echo "  Window $N: extracting frame at t=${{TIME}} ps (position=${{POS}} nm)"
        echo 0 | gmx trjconv -s ./smd/smd.tpr -f ./smd/smd.xtc -dump "$TIME" -o "$CONF" > /dev/null 2>&1
    fi
    N=$((N + 1))
done < ./umbrella/windows.txt

echo "Extracted $N frames"
echo ""

# ============================================
# Phase 3: Run umbrella sampling windows
# ============================================
echo ">>> Phase 3: Running umbrella sampling windows"

N=0
while read -r WID TIME POS; do
    [[ "$WID" == \\#* ]] && continue

    WDIR="./umbrella/window_${{N}}"
    mkdir -p "$WDIR"

    if [ -f "${{WDIR}}/umbrella.gro" ]; then
        echo "  Window $N: already completed, skipping..."
    elif [ -f "${{WDIR}}/umbrella.tpr" ]; then
        echo "  Window $N: resuming from checkpoint..."
        gmx mdrun -ntmpi 1 -ntomp $NCPU -nb gpu -bonded gpu -pme gpu -gpu_id $GPU_ID \\
            -s "${{WDIR}}/umbrella.tpr" -deffnm "${{WDIR}}/umbrella" -v \\
            -cpi "${{WDIR}}/umbrella.cpt" \\
            -pf "${{WDIR}}/pullf.xvg" -px "${{WDIR}}/pullx.xvg"
    else
        echo "  Window $N: starting (position=${{POS}} nm)..."
        gmx grompp -f ../mdps/umbrella.mdp \\
            -c "./umbrella/conf${{N}}.gro" -r "./umbrella/conf${{N}}.gro" \\
            -p topol.top -n index.ndx \\
            -o "${{WDIR}}/umbrella.tpr" -maxwarn 999

        gmx mdrun -ntmpi 1 -ntomp $NCPU -nb gpu -bonded gpu -pme gpu -gpu_id $GPU_ID \\
            -s "${{WDIR}}/umbrella.tpr" -deffnm "${{WDIR}}/umbrella" -v \\
            -pf "${{WDIR}}/pullf.xvg" -px "${{WDIR}}/pullx.xvg"
    fi

    N=$((N + 1))
done < ./umbrella/windows.txt

echo ""
echo "All $N umbrella windows completed!"
echo ""

# ============================================
# Phase 4: WHAM analysis + PMF plot
# ============================================
echo ">>> Phase 4: WHAM analysis"

# Generate file lists for WHAM
echo "Generating tpr-files.dat and pullf-files.dat..."
> ./umbrella/tpr-files.dat
> ./umbrella/pullf-files.dat

for i in $(seq 0 $((N - 1))); do
    echo "./umbrella/window_${{i}}/umbrella.tpr" >> ./umbrella/tpr-files.dat
    echo "./umbrella/window_${{i}}/pullf.xvg" >> ./umbrella/pullf-files.dat
done

echo "Running WHAM..."
gmx wham -it ./umbrella/tpr-files.dat -if ./umbrella/pullf-files.dat \\
    -o ./umbrella/pmf.xvg -hist ./umbrella/histo.xvg \\
    -bsres ./umbrella/pmferror.xvg -nBootstrap $WHAM_BOOTSTRAP \\
    -bs-method b-hist -b $WHAM_BEGIN -unit kCal

echo ""
echo "=============================================="
echo "     WHAM Analysis Complete!"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - PMF profile: ./umbrella/pmf.xvg"
echo "  - Histograms:  ./umbrella/histo.xvg"
echo "  - PMF errors:  ./umbrella/pmferror.xvg"

# Plot PMF profile
if [ -f ./umbrella/pmf.xvg ]; then
    echo ""
    echo ">>> Plotting PMF profile..."
    python3 ./umbrella/plot_pmf.py ./umbrella/pmf.xvg ./umbrella/pmferror.xvg ./umbrella/
    if [ $? -eq 0 ]; then
        echo "Plot saved to ./umbrella/pmf_profile.png"
    else
        echo "Warning: Plotting failed. You can run manually:"
        echo "  python3 ./umbrella/plot_pmf.py ./umbrella/pmf.xvg ./umbrella/pmferror.xvg ./umbrella/"
    fi
fi
"""

        with open(script_path, "w") as f:
            f.write(script_content)

        os.chmod(script_path, 0o755)

        print(f"Umbrella sampling script generated: {script_path}")
        return script_path
