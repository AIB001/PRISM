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

    #: Ligand force fields whose generator can emit the AMBER ligand library
    #: (``.frcmod`` + OFF ``.lib`` + ``savepdb``) that PACKMOL-Memgen consumes
    #: through ``--ligand_param``.  Every other backend produces GROMACS /
    #: CHARMM / OpenMM parameters only, so it cannot embed a ligand in a
    #: bilayer built by the AMBER toolchain.
    MEMBRANE_LIGAND_FORCEFIELDS = ("gaff", "gaff2")

    #: PDB records dropped from the protein when merging it with a ligand:
    #: file terminators (the ligand must come before them) and bookkeeping that
    #: becomes wrong once records are appended (MASTER is a record-count
    #: checksum, CONECT indices are not needed - tleap takes the ligand's
    #: connectivity from the ``.lib``).
    _MERGE_DROP_RECORDS = frozenset({"END", "ENDMDL", "CONECT", "MASTER"})

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
        elif getattr(self, "membrane_mode", False):
            return self.run_membrane()
        else:
            return self.run_normal()

    def run_membrane(self):
        """Build an all-atom protein-in-bilayer system, optionally embedding a ligand."""
        from ..membrane import MembraneBuilder
        from ..membrane.builder import MembraneBuildError

        print_header("PRISM Membrane Builder Workflow")
        try:
            embed_ligand = bool(self.ligand_paths)
            # Reject unsupported ligand inputs before anything is generated:
            # the parameterization below costs minutes of antechamber/sqm time,
            # so a late failure would only waste it and leave a half-built
            # output directory that looks like a real system.
            if embed_ligand:
                self._check_membrane_ligand_inputs()

            self.save_config()

            total_steps = 4 if embed_ligand else 3
            step = 0

            step += 1
            print_step(step, total_steps, "Cleaning protein structure")
            cleaned_protein = self.clean_protein()
            print_success("Protein structure cleaned")

            # Optional PTM / disulfide staging before membrane embedding.
            cleaned_protein = self.apply_ptm(cleaned_protein)

            build_structure = cleaned_protein
            ligand_resname = None
            if embed_ligand:
                step += 1
                print_step(
                    step,
                    total_steps,
                    f"Parameterizing ligand for bilayer embedding ({self.ligand_forcefield.upper()})",
                )
                # Same code path as run_normal(): keeps GAFF/GAFF2 selection,
                # RESP handling and caching identical between the soluble and
                # membrane workflows.
                self.generate_ligand_forcefield()

                ligand_resname = str(self.membrane_config.ligand_resname).strip().upper()
                frcmod, lib, ligand_pdb = self._membrane_ligand_amber_params(ligand_resname)
                if self.membrane_config.has_ligand():
                    # A pre-built pair from YAML would otherwise be replaced
                    # without a word, hiding which parameters were really used.
                    print_warning(
                        "Replacing the membrane YAML ligand_frcmod/ligand_lib pair with the "
                        "parameters PRISM just generated for the supplied ligand file "
                        f"({os.path.basename(self.ligand_paths[0])}); drop the ligand file to "
                        "keep the pre-built pair."
                    )
                self.membrane_config.ligand_frcmod = frcmod
                self.membrane_config.ligand_lib = lib
                # A half-supplied frcmod/lib pair (or any other broken field)
                # must fail here, not inside PACKMOL-Memgen.
                self.membrane_config.raise_for_errors()
                print_success(f"AMBER ligand parameters: {path(frcmod)}")
                print_success(f"AMBER ligand library:    {path(lib)}")

                # PACKMOL-Memgen packs the bilayer around whatever solute PDB it
                # is given, so the ligand has to travel inside that PDB.
                build_structure = self._merge_protein_ligand_pdb(
                    cleaned_protein,
                    ligand_pdb,
                    os.path.join(self.output_dir, "protein_ligand.pdb"),
                    expect_resname=ligand_resname,
                )
                print_success(
                    f"Merged protein + ligand ({ligand_resname}) structure: {path(build_structure)}"
                )

            step += 1
            print_step(step, total_steps, "Reporting membrane-build prerequisites")
            builder = MembraneBuilder(self.membrane_config, verbose=True)
            missing = builder.capability_report()
            if missing:
                print_warning("Some membrane-build tools are missing:")
                for m in missing:
                    print_warning(f"  - {m}")
                print_warning(
                    "PRISM will still emit the orientation + membrane MDP scaffold. "
                    "Install AmberTools (packmol-memgen + ParmEd) or run from that env to complete the bilayer."
                )

            step += 1
            print_step(step, total_steps, "Building membrane system (orient -> bilayer -> convert -> mdps)")
            result = builder.build(build_structure, self.output_dir)
            # Preserve the structured result for API/MCP callers.  The public
            # workflow still returns output_dir for compatibility with the
            # other builder modes.
            self.membrane_result = result

            print_header("Membrane Setup " + ("Complete!" if result.success else "(partial)"))
            print(f"\n  System dir: {path(result.system_dir)}")
            if result.gro:
                print(f"  Structure:  {path(result.gro)}")
                print(f"  Topology:   {path(result.top)}")
            if result.ndx:
                print(f"  Index:      {path(result.ndx)}")
            if result.posres:
                print(f"  Restraints: {', '.join(sorted(result.posres))}")
            if result.mdps:
                print(f"  Membrane MDPs: {', '.join(sorted(result.mdps))}")
            if embed_ligand:
                print(f"  Embedded ligand: residue {number(ligand_resname)} (packed with --keepligs/--ligand_param)")
                print(f"  Solute PDB: {path(build_structure)}")

            if result.success:
                print_header("Ready to Run Membrane MD Simulations!")
                print(f"\nThe membrane system in {path(result.system_dir)} contains:")
                print(f"  - Topology and coordinates ({os.path.basename(result.top)}, {os.path.basename(result.gro)})")
                print(f"  - Index file with SOLU/MEMB/SOLV thermostat groups ({os.path.basename(result.ndx)})")
                if result.posres:
                    print(
                        f"  - Position restraints for the staged equilibration "
                        f"({', '.join(sorted(result.posres))})"
                    )
                print(f"  - MDP protocol: {', '.join(sorted(result.mdps))}")
                print(f"\nTo run the staged membrane workflow:")
                print(f"  1. Navigate to the membrane system directory:")
                print(f"     {path(f'cd {result.system_dir}')}")
                print(f"  2. Execute the simulation script:")
                print(f"     {number('bash localrun.sh')}")
                print(f"\nNote: Adjust GPU and thread settings in localrun.sh as needed")
            else:
                print_warning(result.note)
                print(
                    "\n  A complete membrane build writes topology, coordinates, index.ndx, "
                    "position-restraint itp files and localrun.sh into "
                    f"{path(result.system_dir)}; the items listed above are what this run produced."
                )

            return self.output_dir

        except (ValueError, MembraneBuildError) as exc:
            # A rejected input is a usage error, not a crash: the message already
            # says what is wrong and how to fix it, and burying it under two
            # stack traces (one here, one from the re-raise) hides the only part
            # the user needs to read.
            print_error(f"Membrane build cannot proceed: {exc}")
            raise SystemExit(1) from None
        except Exception as e:
            print_error(f"Membrane workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise

    # ------------------------------------------------------------------ #
    # Ligand-in-bilayer helpers
    # ------------------------------------------------------------------ #

    def _check_membrane_ligand_inputs(self):
        """Reject ligand inputs a membrane build cannot honour.

        Every case here would otherwise end in a system that silently lacks the
        ligand (or an ``AttributeError`` deep inside the force-field backend),
        which is the worst possible outcome for a binding-site simulation.
        """
        if len(self.ligand_paths) > 1:
            names = ", ".join(os.path.basename(p) for p in self.ligand_paths)
            raise ValueError(
                f"Membrane mode can embed only one ligand, but {len(self.ligand_paths)} were "
                f"supplied ({names}). PACKMOL-Memgen's --ligand_param takes a single "
                "FRCMOD:LIB pair, so the remaining ligands could only be dropped silently. "
                "Build one ligand per membrane system, or pre-merge the ligands into a single "
                "AMBER library and pass it through the YAML membrane keys "
                "ligand_frcmod/ligand_lib (with ligand_resname)."
            )

        if self.ligand_forcefield not in self.MEMBRANE_LIGAND_FORCEFIELDS:
            supported = " or ".join(f"--ligand-forcefield {ff}" for ff in self.MEMBRANE_LIGAND_FORCEFIELDS)
            raise ValueError(
                f"Membrane mode cannot embed a ligand parameterized with "
                f"'{self.ligand_forcefield}'. PACKMOL-Memgen builds the bilayer with the AMBER "
                "toolchain and needs AMBER ligand parameters (a .frcmod plus an OFF .lib); only "
                f"the {'/'.join(ff.upper() for ff in self.MEMBRANE_LIGAND_FORCEFIELDS)} backends "
                f"produce those. Re-run with {supported}, or supply pre-built AMBER parameters "
                "through the YAML membrane keys ligand_frcmod/ligand_lib."
            )

        orient = str(getattr(self.membrane_config, "orient", "") or "").strip().lower()
        if orient == "opm":
            raise ValueError(
                "Membrane mode cannot embed a -lf/--ligand-file ligand with orient=opm: OPM "
                "orientation replaces the solute with the structure fetched for the given PDB "
                "id, which would discard the ligand PRISM merged into it. Orient the complex "
                "externally and use orient=preoriented (or orient=memembed, which orients the "
                "supplied structure in place)."
            )

    def _membrane_ligand_amber_params(self, resname):
        """Return ``(frcmod, lib, amber_pdb)`` for the single embedded ligand.

        The files are produced by :meth:`generate_ligand_forcefield` (inside the
        GAFF/GAFF2 generator's ``run()``) and published to the ligand force-field
        directory.  They are read back through the generator's own
        ``get_amber_ligand_params`` so the naming convention lives in exactly one
        place.  ``<RESNAME>_amber.pdb`` is ``savepdb`` output from the very tleap
        unit stored in the ``.lib``, so its atom names, atom order and residue
        name match the library by construction.
        """
        from ..forcefield.gaff import GAFFForceFieldGenerator
        from ..forcefield.gaff2 import GAFF2ForceFieldGenerator

        generator_class = {
            "gaff": GAFFForceFieldGenerator,
            "gaff2": GAFF2ForceFieldGenerator,
        }.get(self.ligand_forcefield)
        if generator_class is None:  # pragma: no cover - guarded earlier
            raise ValueError(
                f"No AMBER ligand-library backend for ligand force field '{self.ligand_forcefield}'."
            )

        # The generator resolves its products relative to
        # ``<output_dir>/<get_output_dir_name()>``; for a single ligand
        # generate_ligand_forcefield() publishes them to
        # ``<self.output_dir>/LIG.amb2gmx``.  Derive the root from the directory
        # actually produced so the two cannot drift apart.
        published_dir = self.lig_ff_dirs[0] if getattr(self, "lig_ff_dirs", None) else None
        generator_root = os.path.dirname(published_dir) if published_dir else self.output_dir
        generator = generator_class(self.ligand_paths[0], generator_root, overwrite=False)

        getter = getattr(generator, "get_amber_ligand_params", None)
        if getter is None:  # pragma: no cover - guarded earlier
            raise ValueError(
                f"The '{self.ligand_forcefield}' force-field backend does not produce AMBER "
                "ligand parameters, which PACKMOL-Memgen requires to embed a ligand."
            )

        params = getter(resname)
        if params is None:
            lig_dir = os.path.join(generator.output_dir, generator.get_output_dir_name())
            hint = ""
            default_resname = getattr(generator, "LIGAND_RESNAME_DEFAULT", "LIG")
            if resname != default_resname:
                fallback = getter(default_resname)
                if fallback:
                    hint = (
                        f" A library WAS produced under the default residue name "
                        f"'{default_resname}' ({os.path.basename(fallback[1])}): PRISM's "
                        f"{self.ligand_forcefield.upper()} backend always names the .lib unit "
                        f"'{default_resname}', and the unit name is what tleap matches against "
                        f"the PDB residue name. Set the membrane ligand residue name to "
                        f"'{default_resname}' (drop --membrane-ligand-resname / the YAML "
                        "membrane ligand_resname override)."
                    )
            raise RuntimeError(
                f"Ligand force field generation finished, but no AMBER ligand library for residue "
                f"'{resname}' was found in {lig_dir} (expected {resname}.frcmod, {resname}.lib and "
                f"{resname}_amber.pdb). PACKMOL-Memgen cannot embed the ligand without it, and "
                "PRISM refuses to fall through to a silently ligand-free bilayer. This usually "
                "means AmberTools (antechamber/parmchk2/tleap) is unavailable or the tleap step "
                "failed - see the '=== Generating AMBER ligand library ===' output above for the "
                "reason." + hint
            )
        return params

    @classmethod
    def _merge_protein_ligand_pdb(cls, protein_pdb, ligand_pdb, output_pdb, expect_resname=None):
        """Write a merged protein + ligand PDB for PACKMOL-Memgen.

        Layout::

            <protein records, END/ENDMDL/CONECT/MASTER removed>
            TER
            <ligand ATOM/HETATM records, verbatim>
            TER
            END

        The ``TER`` before **and** after the ligand block is load-bearing: tleap
        infers chain connectivity from record order, so without them it happily
        bonds the ligand to the neighbouring residue.  Atom serial numbers are
        deliberately left untouched - tleap takes the ligand's connectivity from
        the ``.lib`` unit, not from the PDB, so renumbering would buy nothing and
        risk breaking the correspondence with the library.

        Parameters
        ----------
        protein_pdb, ligand_pdb : str
            Input structures.  ``ligand_pdb`` must be the ``savepdb`` output that
            accompanies the AMBER ``.lib``.
        output_pdb : str
            Destination path.
        expect_resname : str, optional
            When given, every ligand residue name must equal it.  A mismatch is
            fatal: the membrane index groups and position restraints select the
            ligand by residue name, so the wrong name yields an unrestrained,
            wrongly-thermostatted ligand instead of an obvious failure.

        Returns
        -------
        str : ``output_pdb``.
        """

        def read_lines(pdb_path, label):
            try:
                with open(pdb_path, "r") as handle:
                    return handle.read().splitlines()
            except OSError as exc:
                raise ValueError(f"Cannot read {label} PDB {pdb_path!r}: {exc}") from exc

        merged = []
        protein_atoms = 0
        for line in read_lines(protein_pdb, "protein"):
            record = line[:6].strip().upper()
            if record in cls._MERGE_DROP_RECORDS:
                continue
            if record in ("ATOM", "HETATM"):
                protein_atoms += 1
            merged.append(line)
        if not protein_atoms:
            raise ValueError(f"No ATOM/HETATM records found in protein PDB: {protein_pdb}")

        ligand_atoms = [
            line
            for line in read_lines(ligand_pdb, "ligand")
            if line[:6].strip().upper() in ("ATOM", "HETATM")
        ]
        if not ligand_atoms:
            raise ValueError(f"No ATOM/HETATM records found in ligand PDB: {ligand_pdb}")

        if expect_resname:
            expected = str(expect_resname).strip().upper()
            found = sorted({line[17:20].strip().upper() for line in ligand_atoms})
            if found != [expected]:
                raise ValueError(
                    f"Ligand PDB {ligand_pdb} holds residue name(s) {', '.join(found) or '<blank>'} "
                    f"but the membrane build selects the ligand as '{expected}'. The index groups "
                    "and position restraints match on that name, so the mismatch would leave the "
                    "ligand unrestrained and outside the solute thermostat group."
                )

        # Trailing blank lines would push END away from the last record; drop
        # them before deciding whether a TER separator is still needed.
        while merged and not merged[-1].strip():
            merged.pop()
        if merged[-1][:6].strip().upper() != "TER":
            merged.append("TER")
        merged.extend(ligand_atoms)
        merged.append("TER")
        merged.append("END")

        with open(output_pdb, "w") as handle:
            handle.write("\n".join(merged) + "\n")
        return output_pdb

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

            # Step 2b: Apply PTM / disulfide staging (opt-in; no-op otherwise)
            cleaned_protein = self.apply_ptm(cleaned_protein)

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
