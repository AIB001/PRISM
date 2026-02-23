"""
REST2 Workflow - Orchestrates the full REST2 conversion pipeline.

Takes a GMX_PROLIG_MD directory (produced by PRISM's normal workflow)
and generates a GMX_PROLIG_REST2 directory ready for REST2 enhanced sampling.
"""

import os
import time

from .topology import (
    parse_gro,
    find_pocket_residues,
    merge_topology,
    parse_molecules_section,
    count_atoms_per_moltype,
    mark_solute_atoms,
)
from .builder import build_rest2


class REST2Workflow:
    """Orchestrate the full REST2 conversion pipeline.

    Converts a standard GROMACS MD directory (GMX_PROLIG_MD) into a REST2
    replica exchange directory (GMX_PROLIG_REST2).

    Args:
        md_dir: Path to GMX_PROLIG_MD directory
        output_dir: Path for output GMX_PROLIG_REST2 directory
        t_ref: Reference temperature in K (default: 310.0)
        t_max: Maximum effective temperature in K (default: 450.0)
        n_replicas: Number of REST2 replicas (default: 16)
        cutoff: Pocket detection cutoff in nm (default: 0.5)
        lig_name: Ligand residue name(s). Single string ('LIG') or list
                  (['LIG_1', 'LIG_2']) for multi-ligand systems.
        lig_moltype: Ligand moleculetype name(s). Defaults to same as lig_name.
        gmx: Path to gmx executable (default: 'gmx')
    """

    def __init__(
        self,
        md_dir,
        output_dir,
        t_ref=310.0,
        t_max=450.0,
        n_replicas=16,
        cutoff=0.5,
        lig_name="LIG",
        lig_moltype=None,
        gmx="gmx",
    ):
        self.md_dir = os.path.abspath(md_dir)
        self.output_dir = os.path.abspath(output_dir)
        self.t_ref = t_ref
        self.t_max = t_max
        self.n_replicas = n_replicas
        self.cutoff = cutoff
        self.gmx = gmx

        # Normalize lig_name / lig_moltype to lists for uniform handling
        if isinstance(lig_name, str):
            self.lig_names = [lig_name]
        else:
            self.lig_names = list(lig_name)

        if lig_moltype is None:
            self.lig_moltypes = list(self.lig_names)
        elif isinstance(lig_moltype, str):
            self.lig_moltypes = [lig_moltype]
        else:
            self.lig_moltypes = list(lig_moltype)

    def run(self):
        """Execute the full REST2 conversion pipeline.

        Returns:
            str: Path to the output GMX_PROLIG_REST2 directory
        """
        # Import color functions - use try/except for standalone usage
        try:
            from prism.utils.colors import print_header, print_step, print_success, print_error, print_info, path

            use_colors = True
        except ImportError:
            use_colors = False

        if use_colors:
            print_header("REST2 Enhanced Sampling Setup")
        else:
            print("=" * 60)
            print("  REST2 Enhanced Sampling Setup")
            print("=" * 60)

        print(f"  MD directory:   {self.md_dir}")
        print(f"  Output:         {self.output_dir}")
        print(f"  Replicas:       {self.n_replicas}")
        print(f"  T_ref:          {self.t_ref} K")
        print(f"  T_max:          {self.t_max} K")
        print(f"  Cutoff:         {self.cutoff} nm")
        if len(self.lig_names) == 1:
            print(f"  Ligand:         {self.lig_names[0]} (moltype: {self.lig_moltypes[0]})")
        else:
            print(f"  Ligands:        {', '.join(self.lig_names)} (moltypes: {', '.join(self.lig_moltypes)})")
        print()

        # Validate inputs
        topol_path = os.path.join(self.md_dir, "topol.top")
        gro_path = os.path.join(self.md_dir, "solv_ions.gro")

        if not os.path.isdir(self.md_dir):
            msg = f"MD directory not found: {self.md_dir}"
            if use_colors:
                print_error(msg)
            else:
                print(f"ERROR: {msg}")
            raise FileNotFoundError(msg)
        if not os.path.exists(topol_path):
            msg = f"topol.top not found in {self.md_dir}"
            if use_colors:
                print_error(msg)
            else:
                print(f"ERROR: {msg}")
            raise FileNotFoundError(msg)
        if not os.path.exists(gro_path):
            msg = f"solv_ions.gro not found in {self.md_dir}"
            if use_colors:
                print_error(msg)
            else:
                print(f"ERROR: {msg}")
            raise FileNotFoundError(msg)

        t_start = time.time()
        total_steps = 4

        # ---------------------------------------------------------------
        # Step 1: Parse GRO and topology
        # ---------------------------------------------------------------
        step_msg = "Parsing GRO and topology"
        if use_colors:
            print_step(1, total_steps, step_msg)
        else:
            print(f"Step 1/{total_steps}: {step_msg}...")

        atoms, box = parse_gro(gro_path)
        print(f"  Parsed {len(atoms)} atoms, box = {box}")

        molecules = parse_molecules_section(topol_path)
        print(f"  Molecules: {molecules}")

        # ---------------------------------------------------------------
        # Step 2: Merge topology with gmx grompp -pp
        # ---------------------------------------------------------------
        step_msg = "Merging topology (gmx grompp -pp)"
        if use_colors:
            print_step(2, total_steps, step_msg)
        else:
            print(f"\nStep 2/{total_steps}: {step_msg}...")

        os.makedirs(self.output_dir, exist_ok=True)
        processed_top = os.path.join(self.output_dir, "processed.top")
        merge_topology(topol_path, gro_path, processed_top, gmx=self.gmx)

        # Count atoms per moleculetype and build molecules_order
        atom_counts = count_atoms_per_moltype(processed_top)
        print(f"  Atoms per moleculetype: {atom_counts}")

        molecules_order = []
        for name, count in molecules:
            n_atoms = atom_counts.get(name, 0)
            molecules_order.append((name, count, n_atoms))

        # Find pocket residues (searches near all ligand resnames)
        pocket = find_pocket_residues(atoms, box, self.lig_names, molecules_order, self.cutoff)

        lig_label = ", ".join(self.lig_names)
        for chain, resids in sorted(pocket.items()):
            print(f"  {chain}: {len(resids)} residues within {self.cutoff} nm of {lig_label}")
            sorted_resids = sorted(resids)
            if len(sorted_resids) > 10:
                print(f"    Residues: {sorted_resids[:5]}...{sorted_resids[-5:]}")
            else:
                print(f"    Residues: {sorted_resids}")

        # ---------------------------------------------------------------
        # Step 3: Mark solute atoms
        # ---------------------------------------------------------------
        step_msg = "Marking solute atoms"
        if use_colors:
            print_step(3, total_steps, step_msg)
        else:
            print(f"\nStep 3/{total_steps}: {step_msg}...")

        marked_top = os.path.join(self.output_dir, "rest2_marked.top")
        mark_solute_atoms(processed_top, marked_top, self.lig_moltypes, pocket)

        # ---------------------------------------------------------------
        # Step 4: Build REST2 directory
        # ---------------------------------------------------------------
        step_msg = "Building REST2 directory (scaled topologies, MDPs, scripts)"
        if use_colors:
            print_step(4, total_steps, step_msg)
        else:
            print(f"\nStep 4/{total_steps}: {step_msg}...")

        build_rest2(
            md_dir=self.md_dir,
            output_dir=self.output_dir,
            marked_top_path=marked_top,
            original_top_path=topol_path,
            gro_path=gro_path,
            n_replicas=self.n_replicas,
            t_ref=self.t_ref,
            t_max=self.t_max,
            processed_top_path=processed_top,
        )

        elapsed = time.time() - t_start

        # ---------------------------------------------------------------
        # Completion
        # ---------------------------------------------------------------
        if use_colors:
            print_success(f"REST2 setup complete ({elapsed:.1f}s)")
            print_info(f"Output: {path(self.output_dir)}")
        else:
            print(f"\n{'=' * 60}")
            print(f"  REST2 setup complete in {elapsed:.1f} seconds")
            print(f"  Output: {self.output_dir}")
            print(f"{'=' * 60}")

        print(f"\n  Next steps:")
        print(f"  1. cd {self.output_dir}")
        print(f"  2. bash rest2_run.sh")

        return self.output_dir
