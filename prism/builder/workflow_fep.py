#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder - FEP Workflow Mixin

Handles the FEP (Free Energy Perturbation) workflow: building bound/unbound
MD systems, generating hybrid topology via atom mapping, and creating the
FEP scaffold.
"""

import os
import shutil
from pathlib import Path

from ..utils.colors import (
    print_header,
    print_step,
    print_success,
    print_error,
    path,
)
from ..forcefield.registry import (
    get_forcefield_output_subdirs,
    iter_existing_ligand_output_dirs,
)


class FEPWorkflowMixin:
    """Mixin providing FEP workflow for relative binding free energy calculations."""

    @staticmethod
    def _normalize_fep_output_dir(output_dir: str) -> str:
        """Collapse accidental ``.../GMX_PROLIG_FEP/GMX_PROLIG_FEP`` nesting."""
        normalized = os.path.abspath(output_dir)
        while (
            os.path.basename(normalized) == "GMX_PROLIG_FEP"
            and os.path.basename(os.path.dirname(normalized)) == "GMX_PROLIG_FEP"
        ):
            normalized = os.path.dirname(normalized)
        return normalized

    def run_fep(self):
        """Run the FEP workflow: build standard MD systems, generate hybrid topology, create FEP scaffold"""
        from ..fep.modeling import FEPScaffoldBuilder

        print_header("PRISM FEP Builder Workflow")

        if not self.mutant_ligand:
            raise ValueError("FEP mode requires --mutant ligand file")

        try:
            # Create FEP output directory first
            # Avoid double-nesting: if output_dir already ends with GMX_PROLIG_FEP, use it directly
            normalized_output_dir = self._normalize_fep_output_dir(self.output_dir)
            if os.path.basename(normalized_output_dir) == "GMX_PROLIG_FEP":
                fep_output = normalized_output_dir
            else:
                fep_output = os.path.join(normalized_output_dir, "GMX_PROLIG_FEP")

            os.makedirs(fep_output, exist_ok=True)

            root_md_dir = os.path.join(normalized_output_dir, "GMX_PROLIG_MD")
            if os.path.isdir(root_md_dir) and not os.listdir(root_md_dir):
                os.rmdir(root_md_dir)
                print(f"  ✓ Removed empty directory: {root_md_dir}")

            # Phase 1: Build bound system with reference ligand
            print_step(1, 5, "Building bound system with reference ligand")
            bound_output = os.path.join(fep_output, "_build/bound_md")

            self._build_standard_system(bound_output, use_protein=True)
            bound_system_dir = os.path.join(bound_output, "GMX_PROLIG_MD")

            # Save reference ligand FF directory from the freshly built bound system
            ref_ff_dir = self._resolve_generated_ligand_ff_dir(bound_output)
            print_success(f"Bound system built: {bound_system_dir}")
            print(f"  Reference ligand FF: {ref_ff_dir}")

            # Phase 2: Build unbound system (ligand only) with same box size as bound
            print_step(2, 5, "Building unbound system (ligand in water)")
            unbound_output = os.path.join(fep_output, "_build/unbound_md")

            # Read bound system box size
            bound_conf = os.path.join(bound_system_dir, "solv_ions.gro")
            box_size = self._read_box_size(bound_conf)
            print(f"  Using bound system box size: {box_size[0]:.3f} x {box_size[1]:.3f} x {box_size[2]:.3f} nm")

            self._build_standard_system(unbound_output, use_protein=False, box_size=box_size)
            unbound_system_dir = os.path.join(unbound_output, "GMX_PROLIG_MD")
            print_success(f"Unbound system built: {unbound_system_dir}")

            # Phase 3: Generate hybrid topology via service
            print_step(3, 5, "Generating hybrid topology via atom mapping")

            # ref_ff_dir was saved in Phase 1
            # Generate mutant ligand FF in _build directory
            mut_ff_output = os.path.join(fep_output, "_build/mutant_ligand_ff")
            self._generate_mutant_ligand_ff(self.mutant_ligand, mut_ff_output)
            mut_ff_dir = self._resolve_generated_ligand_ff_dir(mut_ff_output)

            # Build hybrid topology using HybridBuildService
            from ..fep.modeling.hybrid_service import HybridBuildService

            hybrid_output = os.path.join(fep_output, "common/hybrid")
            os.makedirs(hybrid_output, exist_ok=True)

            ref_visual_coord_source = str(
                HybridBuildService.resolve_coord_source(self.ligand_paths[0] if self.ligand_paths else "", ref_ff_dir)
            )
            mut_visual_coord_source = str(HybridBuildService.resolve_coord_source(self.mutant_ligand or "", mut_ff_dir))

            print(f"  Building hybrid topology:")
            print(f"    Reference FF: {ref_ff_dir}")
            print(f"    Mutant FF:    {mut_ff_dir}")
            print(f"    Output dir:   {hybrid_output}")

            hybrid_service = HybridBuildService(
                dist_cutoff=self.distance_cutoff,
                charge_cutoff=self.charge_cutoff,
                charge_common=self.charge_strategy,
                charge_reception=self.charge_reception,
            )

            hybrid_result = hybrid_service.build_from_forcefield_dirs(
                reference_ligand_dir=ref_ff_dir,
                mutant_ligand_dir=mut_ff_dir,
                hybrid_output_dir=hybrid_output,
                reference_coord_source=ref_visual_coord_source,
                mutant_coord_source=mut_visual_coord_source,
                molecule_name="HYB",
            )

            hybrid_itp = str(hybrid_result.hybrid_itp)
            hybrid_gro = str(hybrid_result.hybrid_gro)
            mapping = hybrid_result.mapping

            print(f"  Hybrid topology:")
            print(f"    ITP:  {hybrid_itp}")
            print(f"    GRO:  {hybrid_gro}")
            print(f"    Common: {len(mapping.common)}")
            print(f"    Transformed (ref): {len(mapping.transformed_a)}")
            print(f"    Transformed (mut): {len(mapping.transformed_b)}")

            # Export mapping HTML via MappingReportService
            from ..fep.visualize.reporting import MappingReportService

            report_service = MappingReportService(
                distance_cutoff=self.distance_cutoff,
                charge_cutoff=self.charge_cutoff,
                charge_strategy=self.charge_strategy,
                charge_reception=self.charge_reception,
            )

            report_service.generate_fep_mapping_html(
                output_dir=hybrid_output,
                ref_coord_source=ref_visual_coord_source,
                mut_coord_source=mut_visual_coord_source,
                mapping=mapping,
                ref_atoms=hybrid_result.ref_atoms,
                mut_atoms=hybrid_result.mut_atoms,
                fep_config=self.fep_config,
                config_file=self.config_file,
                ligand_forcefield=self.ligand_forcefield,
                forcefield_paths=self.forcefield_paths,
                distance_cutoff=self.distance_cutoff,
                charge_cutoff=self.charge_cutoff,
                charge_strategy=self.charge_strategy,
                charge_reception=self.charge_reception,
            )

            # Phase 4: Create FEP scaffold with complete systems
            print_step(4, 5, "Creating FEP scaffold with complete systems")

            fep_builder = FEPScaffoldBuilder(
                output_dir=fep_output,
                lambda_windows=self.lambda_windows,
                lambda_strategy=self.lambda_strategy,
                lambda_distribution=self.lambda_distribution,
                config=self.config,
                overwrite=False,
            )

            layout = fep_builder.build_from_components(
                receptor_pdb=self.protein_path,
                hybrid_itp=hybrid_itp,
                reference_ligand_dir=ref_ff_dir,
                mutant_ligand_dir=mut_ff_dir,
                bound_system_dir=bound_system_dir,
                unbound_system_dir=unbound_system_dir,
            )

            print_success(f"FEP scaffold created: {fep_output}")

            # Phase 5: Clean up intermediate build files
            print_step(5, 5, "Cleaning up intermediate build files")
            self._cleanup_fep_build_artifacts(fep_output)

            # Also clean up any files in output_dir root (outside GMX_PROLIG_FEP)
            cleanup_items = [
                os.path.join(normalized_output_dir, "mdps"),
            ]
            cleanup_items.extend(str(p) for p in iter_existing_ligand_output_dirs(normalized_output_dir))

            root_md_dir = os.path.join(normalized_output_dir, "GMX_PROLIG_MD")
            if os.path.isdir(root_md_dir) and not os.listdir(root_md_dir):
                cleanup_items.append(root_md_dir)

            for item in cleanup_items:
                if os.path.exists(item):
                    if os.path.isdir(item):
                        shutil.rmtree(item)
                    else:
                        os.remove(item)
                    print(f"  ✓ Removed {os.path.basename(item)}")

            print_success("Cleanup complete")

            print_header("FEP Workflow Complete!")
            print(f"\n  FEP scaffold:       {path(fep_output)}")
            print(f"\n  To run FEP simulations:")
            print(f"  1. cd {fep_output} && ./run_fep.sh bound")
            print(f"  2. cd {fep_output} && ./run_fep.sh unbound")
            print(f"  3. Analyze with gmx bar or alchemical-analysis")

            return self.output_dir

        except Exception as e:
            print_error(f"FEP workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise

    def _cleanup_fep_build_artifacts(self, fep_output: str) -> None:
        """Prune heavy intermediate system builds while keeping initial ligand FF outputs.

        Keep `_build/.../LIG.*` directories as the canonical backup of the initial
        reference/mutant ligand force fields. Remove nested `GMX_PROLIG_MD`
        directories generated only for scaffold assembly.
        """
        build_dir = Path(fep_output) / "_build"
        if not build_dir.exists():
            return

        for md_dir in sorted(build_dir.glob("**/GMX_PROLIG_MD")):
            if md_dir.is_dir():
                shutil.rmtree(md_dir)
                print(f"  ✓ Removed intermediate system directory: {md_dir}")

    def _read_box_size(self, gro_file: str) -> tuple:
        """Read box size from GRO file last line"""
        with open(gro_file, "r") as f:
            lines = f.readlines()
            box_line = lines[-1].strip().split()
            return tuple(float(x) for x in box_line[:3])

    def _build_standard_system(self, output_dir: str, use_protein: bool, box_size: tuple = None):
        """Build standard MD system (bound or unbound)"""
        os.makedirs(output_dir, exist_ok=True)

        original_fep = self.fep_mode
        original_output = self.output_dir
        original_protein = self.protein_path
        original_ligand_paths = self.ligand_paths
        original_lig_ff_dirs = list(self.lig_ff_dirs)

        # Convert ligand paths to absolute paths BEFORE changing output_dir
        # to avoid path resolution issues
        if isinstance(self.ligand_paths, str):
            self.ligand_paths = os.path.abspath(self.ligand_paths)
        elif isinstance(self.ligand_paths, list):
            self.ligand_paths = [os.path.abspath(p) for p in self.ligand_paths]

        # Also convert protein path to absolute
        if self.protein_path:
            self.protein_path = os.path.abspath(self.protein_path)

        self.fep_mode = False
        self.output_dir = output_dir

        try:
            # Update sub-component output directories
            self.system_builder.output_dir = Path(output_dir)
            self.system_builder.model_dir = Path(output_dir) / "GMX_PROLIG_MD"
            self.system_builder.model_dir.mkdir(exist_ok=True)
            self.mdp_generator.output_dir = output_dir

            if use_protein:
                # Normal protein+ligand system
                self.run_normal()
            else:
                # Ligand-only system: skip protein processing
                self._build_ligand_only_system(output_dir, box_size=box_size)
        finally:
            # Restore shared builder state even if the temporary build fails.
            self.fep_mode = original_fep
            self.output_dir = original_output
            self.protein_path = original_protein
            self.ligand_paths = original_ligand_paths
            self.lig_ff_dirs = original_lig_ff_dirs
            self.system_builder.output_dir = Path(original_output)
            self.system_builder.model_dir = Path(original_output) / "GMX_PROLIG_MD"
            self.mdp_generator.output_dir = original_output

    def _build_ligand_only_system(self, output_dir: str, box_size: tuple = None):
        """Build ligand-only system (no protein) for unbound FEP leg.

        Parameters
        ----------
        output_dir : str
            Output directory for the system.
        box_size : tuple, optional
            Box size (x, y, z) in nm. If provided, uses this exact box size.
            If None, uses box_distance parameter to create box around ligand.
        """
        print_step(1, 5, "Generating ligand force field for unbound system")

        original_output = self.output_dir
        self.output_dir = output_dir
        self.generate_ligand_forcefield()
        self.output_dir = original_output

        from ..fep.modeling import LigandOnlySystemBuilder

        ligand_builder = LigandOnlySystemBuilder(
            system_builder=self.system_builder,
            forcefield=self.forcefield,
            water_model=self.water_model,
            gromacs_env=self.gromacs_env,
            overwrite=self.overwrite,
            config=self.config,
            cgenff_supplement_builder=getattr(self.system_builder, "_build_cgenff_parameter_supplement", None),
        )
        ligand_builder.build(output_dir=output_dir, ligand_ff_dir=self.lig_ff_dirs[0], box_size=box_size)

    def _resolve_generated_ligand_ff_dir(self, output_dir: str) -> str:
        """Return the generated ligand force-field directory for the current ligand FF."""
        output_path = Path(output_dir)
        search_roots = [output_path]
        ligand_ff_dir = output_path / "Ligand_Forcefield"
        if ligand_ff_dir.exists():
            search_roots.append(ligand_ff_dir)

        ff_dirs = []
        seen = set()

        for root in search_roots:
            forcefield = getattr(self, "ligand_forcefield", None)
            preferred = iter_existing_ligand_output_dirs(root, forcefield=forcefield) if forcefield else ()
            for candidate in preferred:
                if candidate.is_dir() and candidate not in seen:
                    ff_dirs.append(candidate)
                    seen.add(candidate)
            for subdir in get_forcefield_output_subdirs(forcefield) if forcefield else ():
                for candidate in sorted(root.glob(f"{subdir}*")):
                    if candidate.is_dir() and candidate not in seen:
                        ff_dirs.append(candidate)
                        seen.add(candidate)
            for candidate in iter_existing_ligand_output_dirs(root):
                if candidate.is_dir() and candidate not in seen:
                    ff_dirs.append(candidate)
                    seen.add(candidate)
            for candidate in sorted(root.glob("LIG.*")):
                if candidate.is_dir() and candidate not in seen:
                    ff_dirs.append(candidate)
                    seen.add(candidate)

        if len(ff_dirs) == 0:
            raise FileNotFoundError(
                f"Expected at least one generated ligand FF dir in {output_dir}, found {len(ff_dirs)}"
            )
        # For multi-ligand case:
        # - FEP always builds reference ligand first (ligand index 0)
        # - Multi-ligand renaming names it LIG.sp2gmx_1 (1-based numbering)
        # - Look for _1 suffix first, then fall back to sorted first
        ff_dirs_names = [p.name for p in ff_dirs]
        # Check for directory ending with _1 (reference ligand)
        ref_candidates = [p for p in ff_dirs if p.name.endswith("_1")]
        if ref_candidates:
            # Found explicit _1, use it
            return str(ref_candidates[0])
        # No explicit _1, use sorted first
        return str(ff_dirs[0])

    def _generate_mutant_ligand_ff(self, mutant_ligand: str, output_dir: str):
        """Generate force field for mutant ligand"""
        # Reuse existing force field generation logic

        # Temporarily save current state
        original_ligands = self.ligand_paths
        original_output = self.output_dir
        original_lig_ff_dirs = self.lig_ff_dirs
        original_ff_paths = self.forcefield_paths

        # Set up for mutant ligand
        self.ligand_paths = [mutant_ligand]
        self.output_dir = output_dir
        # Use the second forcefield_path for mutant ligand (index 1)
        # In FEP mode: forcefield_paths = [ref_ff, mut_ff]
        if self.forcefield_paths and len(self.forcefield_paths) > 1:
            self.forcefield_paths = [self.forcefield_paths[1]]

        try:
            # Generate FF
            self.generate_ligand_forcefield()
        finally:
            # Restore shared builder state even if FF generation fails.
            self.ligand_paths = original_ligands
            self.output_dir = original_output
            self.lig_ff_dirs = original_lig_ff_dirs
            self.forcefield_paths = original_ff_paths
