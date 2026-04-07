#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder - Workflow Mixin

Handles the main workflow orchestration: normal MD, REST2, and MM/PBSA modes.
"""

import os
import shutil

from ..utils.colors import (
    print_header,
    print_step,
    print_success,
    print_error,
    print_warning,
    path,
    number,
)


class WorkflowMixin:
    """Mixin providing workflow orchestration for normal MD, REST2, and MM/PBSA modes."""

    def run(self):
        """Run the complete workflow, dispatching to the appropriate mode."""
        if self.pmf_mode:
            return self.run_pmf()
        elif self.rest2_mode:
            return self.run_rest2()
        elif self.mmpbsa_mode:
            return self.run_mmpbsa()
        elif self.fep_mode:
            return self.run_fep()
        else:
            return self.run_normal()

    def run_normal(self):
        """Run the normal (non-PMF) workflow"""
        print_header("PRISM Builder Workflow")

        try:
            # Save configuration for reference
            self.save_config()

            # Step 1: Generate ligand force field
            print_step(1, 6, f"Generating ligand force field ({self.ligand_forcefield.upper()})")
            self.generate_ligand_forcefield()
            print_success(f"Ligand force field generated")

            # Step 2: Clean protein
            print_step(2, 6, "Cleaning protein structure")
            cleaned_protein = self.clean_protein()
            print_success("Protein structure cleaned")

            # Step 3: Build model
            print_step(3, 6, "Building GROMACS system")
            model_dir = self.build_model(cleaned_protein)
            if not model_dir:
                raise RuntimeError("Failed to build model")
            print_success("GROMACS system built")

            # Step 4: Generate MDP files
            print_step(4, 6, "Generating MD parameter files")
            self.generate_mdp_files()
            print_success("MDP files generated")

            # Step 5: Cleanup
            print_step(5, 6, "Cleaning up temporary files")
            self.cleanup()
            print_success("Cleanup completed")

            # Step 6: Generate local run script
            print_step(6, 6, "Generating simulation run script")
            script_path = self.generate_localrun_script()
            print_success("Run script generated")

            print_header("Workflow Complete!")
            print(f"\nOutput files are in: {path(self.output_dir)}")
            print(f"MD system files are in: {path(os.path.join(self.output_dir, 'GMX_PROLIG_MD'))}")
            print(f"MDP files are in: {path(self.mdp_generator.mdp_dir)}")
            print(f"Configuration saved in: {path(os.path.join(self.output_dir, 'prism_config.yaml'))}")
            print(f"\nProtein force field used: {number(self.forcefield['name'])}")
            print(f"Ligand force field used: {number(self.ligand_forcefield.upper())}")
            print(f"Water model used: {number(self.water_model['name'])}")

            if script_path:
                gmx_md_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
                print_header("Ready to Run MD Simulations!")
                print(f"\nTo run the MD workflow:")
                print(f"  1. Navigate to the MD directory:")
                print(f"     {path(f'cd {gmx_md_dir}')}")
                print(f"  2. Execute the simulation script:")
                print(f"     {number('bash localrun.sh')}")
                print(f"\nThe script will run:")
                print(f"  - Energy Minimization (EM)")
                print(f"  - NVT Equilibration")
                print(f"  - NPT Equilibration")
                print(f"  - Production MD")
                print(f"\nNote: Adjust GPU and thread settings in localrun.sh as needed")

            return self.output_dir

        except Exception as e:
            print_error(f"Workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise

    def run_rest2(self):
        """Run the REST2 workflow: build standard MD first, then convert to REST2"""
        print_header("PRISM REST2 Builder Workflow")

        try:
            # Phase 1: Build standard MD system via normal workflow
            print_step(1, 2, "Building standard MD system (normal workflow)")
            self.run_normal()

            # Phase 2: Convert GMX_PROLIG_MD -> GMX_PROLIG_REST2
            print_step(2, 2, "Converting to REST2 replica exchange setup")

            from ..rest2 import REST2Workflow

            md_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
            rest2_dir = os.path.join(self.output_dir, "GMX_PROLIG_REST2")

            # Ligand residue name: PRISM always standardizes to 'LIG' (single)
            # or 'LIG_1','LIG_2',... (multi-ligand)
            if self.ligand_count == 1:
                lig_name = "LIG"
            else:
                lig_name = [f"LIG_{i}" for i in range(1, self.ligand_count + 1)]

            workflow = REST2Workflow(
                md_dir=md_dir,
                output_dir=rest2_dir,
                t_ref=self.t_ref,
                t_max=self.t_max,
                n_replicas=self.n_replicas,
                cutoff=self.rest2_cutoff,
                lig_name=lig_name,
                gmx=self.gromacs_env.gmx_command,
            )

            workflow.run()

            print_header("REST2 Setup Complete!")
            print(f"\n  Standard MD files: {path(md_dir)}")
            print(f"  REST2 files:       {path(rest2_dir)}")
            print(f"\n  To run REST2 simulations:")
            print(f"  1. cd {rest2_dir}")
            print(f"  2. bash rest2_run.sh")

            return self.output_dir

        except Exception as e:
            print_error(f"REST2 workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise

    def run_mmpbsa(self):
        """Run MMPBSA workflow: build standard MD first, then set up for MMPBSA.

        Two sub-modes:
        - Single-frame (default): EM->NVT->NPT->gmx_MMPBSA on equilibrated structure.
        - Trajectory-based (--mmpbsa-traj N): EM->NVT->NPT->Production MD->gmx_MMPBSA.
        """
        try:
            from ..mmpbsa.generator import MMPBSAGenerator
        except ImportError:
            from prism.mmpbsa.generator import MMPBSAGenerator

        is_trajectory = self.mmpbsa_traj_ns is not None
        mode_label = "Trajectory" if is_trajectory else "Single-frame"
        print_header(f"PRISM MM/PBSA Builder Workflow ({mode_label})")

        try:
            # For trajectory mode, override production MD length
            if is_trajectory:
                self.config["simulation"]["production_time_ns"] = self.mmpbsa_traj_ns

            # Phase 1: Build standard MD system via normal workflow
            print_step(1, 2, "Building standard MD system (normal workflow)")
            self.run_normal()

            # Phase 2: Convert GMX_PROLIG_MD -> GMX_PROLIG_MMPBSA
            print_step(2, 2, "Setting up MM/PBSA calculation files")
            md_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
            mmpbsa_dir = os.path.join(self.output_dir, "GMX_PROLIG_MMPBSA")

            # Handle existing MMPBSA directory from previous runs
            if os.path.exists(mmpbsa_dir):
                print_warning(f"Removing existing {os.path.basename(mmpbsa_dir)}/ from previous run")
                shutil.rmtree(mmpbsa_dir)

            os.rename(md_dir, mmpbsa_dir)

            # Clean up files that are irrelevant for MMPBSA workflow
            localrun_path = os.path.join(mmpbsa_dir, "localrun.sh")
            if os.path.exists(localrun_path):
                os.remove(localrun_path)

            # In single-frame mode, remove the unnecessary md.mdp
            if not is_trajectory:
                md_mdp = os.path.join(self.output_dir, "mdps", "md.mdp")
                if os.path.exists(md_mdp):
                    os.remove(md_mdp)

            # Note: -lm (ligand mol2) is NOT needed for gmx_MMPBSA when
            # -cp topol.top is provided, because the GROMACS topology already
            # contains all ligand parameters. gmx_MMPBSA converts internally
            # via parmed. Omitting -lm avoids parmchk2 failures caused by
            # atom type mismatches in the original mol2 file.

            # Generate mmpbsa.in and mmpbsa_run.sh via MMPBSAGenerator
            generator = MMPBSAGenerator(
                mmpbsa_dir=mmpbsa_dir,
                protein_ff_name=self.forcefield["name"],
                ligand_forcefield=self.ligand_forcefield,
                temperature=self.config["simulation"]["temperature"],
                traj_ns=self.mmpbsa_traj_ns,
                trajectory_interval_ps=self.config.get("output", {}).get("trajectory_interval_ps", 500),
                gmx2amber=self.gmx2amber,
            )
            generator.generate_input()
            generator.generate_script()

            mode_desc = "trajectory" if is_trajectory else "single-frame"
            print_success(f"Generated mmpbsa.in ({mode_desc} mode)")
            print_success(f"Generated mmpbsa_run.sh ({mode_desc} mode)")

            # Print completion info
            print_header("MM/PBSA Setup Complete!")
            print(f"\n  MMPBSA files: {path(mmpbsa_dir)}")
            if is_trajectory:
                print(f"  Sub-mode: Trajectory-based ({self.mmpbsa_traj_ns} ns production MD)")
            else:
                print(f"  Sub-mode: Single-frame (no production MD)")
            print(f"\n  To run MM/PBSA workflow:")
            print(f"  1. cd {mmpbsa_dir}")
            print(f"  2. bash mmpbsa_run.sh")

            return self.output_dir

        except Exception as e:
            print_error(f"MMPBSA workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise
