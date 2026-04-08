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
from ..utils.system.topology import _parse_atomtypes, _should_vendor_forcefield
from ..forcefield.registry import (
    get_forcefield_output_subdirs,
    iter_existing_ligand_output_dirs,
    resolve_ligand_artifact,
)

try:
    from ..forcefield.adapters import CGenFFAdapter
except ImportError:
    # Fallback if adapters not available
    class CGenFFAdapter:
        @classmethod
        def should_include_ligand_atomtypes(cls, ff_name, lig_ff_path):
            return not (ff_name.startswith("charmm") and (lig_ff_path / "charmm36.ff").exists())

        @classmethod
        def is_charmm_forcefield(cls, ff_name):
            return ff_name.lower().startswith("charmm")

        @classmethod
        def detect(cls, lig_ff_path):
            return (lig_ff_path / "charmm36.ff").exists()


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

            # Resolve coordinate sources (MOL2 > PDB > GRO)
            def _resolve_coord_source(path: str, fallback: str) -> str:
                """If path is a directory (RTF), find PDB inside; else use as-is."""
                if path and os.path.isdir(path):
                    pdbs = list(Path(path).glob("*.pdb"))
                    return str(pdbs[0]) if pdbs else fallback
                return path if path and os.path.exists(path) else fallback

            ref_visual_coord_source = _resolve_coord_source(
                self.ligand_paths[0] if self.ligand_paths else "", ref_ff_dir
            )
            mut_visual_coord_source = _resolve_coord_source(self.mutant_ligand or "", mut_ff_dir)

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

        model_dir = Path(output_dir) / "GMX_PROLIG_MD"
        model_dir.mkdir(exist_ok=True, parents=True)

        ligand_ff_dir = Path(self.lig_ff_dirs[0])
        ligand_gro = Path(self._resolve_ligand_ff_artifact(str(ligand_ff_dir), "LIG.gro"))

        print_step(2, 5, "Creating ligand-only topology")

        # Copy force field to working directory (same as bound system)
        ff_name = self.forcefield["name"]
        water_model = self.water_model["name"]
        ff_idx = self.forcefield.get("index")
        ff_info = None
        use_ions_itp = "ions.itp"  # Default to old-style

        # If no index, try to find it by name
        if not ff_idx and self.gromacs_env:
            ff_name = self.forcefield.get("name")
            if ff_name:
                ff_name_lower = ff_name.lower()
                if ff_name_lower in self.gromacs_env.force_field_names:
                    ff_idx = self.gromacs_env.force_field_names[ff_name_lower]

        if ff_idx and self.gromacs_env:
            # Get force field info from environment
            if ff_idx in self.gromacs_env.force_fields:
                ff_info = self.gromacs_env.force_fields[ff_idx]

        # Determine which ion file to use (check before copying)
        if ff_info and "path" in ff_info:
            ff_path = Path(ff_info["path"])
            if ff_path.exists() and ff_path.is_dir():
                water_specific_ions = f"ions_{water_model}.itp"
                generic_ions = "ions.itp"

                # Check in force field directory
                if (ff_path / water_specific_ions).exists():
                    use_ions_itp = water_specific_ions
                    print(f"  Using water-specific ion file: {water_specific_ions}")
                elif (ff_path / generic_ions).exists():
                    use_ions_itp = generic_ions
                    print(f"  Using generic ion file: {generic_ions}")
                else:
                    # Neither exists - keep default and warn
                    print(f"  ⚠ Warning: Neither {water_specific_ions} nor {generic_ions} found in {ff_path}")
                    print(f"  Will attempt to use {use_ions_itp} (may fail at grompp)")

        if ff_info and "dir" in ff_info:
            ff_basename = ff_info["dir"]
            local_ff_dir = Path(model_dir) / ff_basename

            ff_path = Path(ff_info.get("path", ""))
            if not ff_path.exists() or not ff_path.is_dir():
                print(f"Warning: Force field path not found: {ff_path}")
                print(f"Using force field from system search paths: {ff_name}")
            else:
                # Validate force field completeness (support both old and new force field structures)
                required_files = ["forcefield.itp", "ffbonded.itp", "ffnonbonded.itp"]

                # Add the detected ion file to required files
                if use_ions_itp:
                    required_files.append(use_ions_itp)

                missing = [f for f in required_files if not (ff_path / f).exists()]
                if missing:
                    raise RuntimeError(
                        f"Force field is incomplete: {ff_info['name']}\n"
                        f"Missing files: {', '.join(missing)}\n"
                        f"Force field path: {ff_path}\n"
                        f"Source: {ff_info.get('source', 'unknown')}\n"
                        f"Please check your force field installation."
                    )

                if _should_vendor_forcefield(ff_info, model_dir):
                    # Copy if not exists, or validate and replace if incomplete
                    need_copy = not local_ff_dir.exists()
                    if not need_copy and self.overwrite:
                        missing_local = [f for f in required_files if not (local_ff_dir / f).exists()]
                        if missing_local:
                            print(f"Existing force field copy is incomplete, replacing...")
                            shutil.rmtree(local_ff_dir)
                            need_copy = True

                    if need_copy:
                        shutil.copytree(ff_path, local_ff_dir)
                        source_str = ff_info.get("source", ff_path)
                        print(f"Copied force field to working directory: {local_ff_dir}")
                        print(f"  Source: {source_str}")
                    else:
                        print(f"Using existing force field copy: {local_ff_dir}")
                else:
                    print(f"Using installed force field in place: {ff_path}")

        topol_path = model_dir / "topol.top"
        # Compute ligand include paths relative to the working model directory.
        # This must work for both legacy single-ligand layouts (output_dir/LIG.*)
        # and newer multi-ligand layouts (output_dir/Ligand_Forcefield/LIG.*).
        ligand_rel_dir = os.path.relpath(ligand_ff_dir, model_dir)
        cgenff_supplement = None
        main_atomtypes = set()
        if CGenFFAdapter.is_charmm_forcefield(ff_name) and CGenFFAdapter.detect(ligand_ff_dir):
            main_ff_dir = model_dir / f"{ff_name}.ff"
            if not main_ff_dir.exists() and ff_info and ff_info.get("path"):
                main_ff_dir = Path(ff_info["path"])
            ffnonbonded = main_ff_dir / "ffnonbonded.itp"
            if ffnonbonded.exists():
                main_atomtypes = _parse_atomtypes(ffnonbonded)
            charmm_ff_dir = CGenFFAdapter.find_charmm_ff_dir(ligand_ff_dir)
            cgenff_supplement = self.system_builder._build_cgenff_parameter_supplement(
                lig_ff_path=ligand_ff_dir,
                lig_itp_path=Path(self._resolve_ligand_ff_artifact(str(ligand_ff_dir), "LIG.itp")),
                charmm_ff_dir=charmm_ff_dir,
                main_bonded_files=[
                    p for p in [main_ff_dir / "ffbonded.itp", main_ff_dir / "ffmissingdihedrals.itp"] if p.exists()
                ],
                main_nonbonded_files=[p for p in [ffnonbonded] if p.exists()],
            )

        with open(topol_path, "w") as f:
            f.write("; Ligand-only topology for FEP unbound leg\n")
            f.write(f'#include "{ff_name}.ff/forcefield.itp"\n')
            atomtypes_lines = []
            atomtypes_itp_path = resolve_ligand_artifact(ligand_ff_dir, "atomtypes_LIG.itp")
            include_ligand_atomtypes = CGenFFAdapter.should_include_ligand_atomtypes(ff_name, ligand_ff_dir)
            if include_ligand_atomtypes and atomtypes_itp_path is not None and atomtypes_itp_path.exists():
                for line in atomtypes_itp_path.read_text().splitlines():
                    stripped = line.strip()
                    if not stripped or stripped.startswith(";") or stripped.startswith("["):
                        continue
                    atomtype = stripped.split()[0]
                    if atomtype not in main_atomtypes:
                        atomtypes_lines.append(line)
            if atomtypes_lines:
                f.write("\n; Include ligand-specific atomtypes\n")
                f.write("[ atomtypes ]\n")
                for line in atomtypes_lines:
                    f.write(f"{line}\n")
                f.write("\n")
            if cgenff_supplement is not None:
                f.write(f'#include "{ligand_rel_dir}/{cgenff_supplement.name}"\n')
            f.write(f'#include "{ligand_rel_dir}/LIG.itp"\n')
            f.write(f'#include "{ff_name}.ff/{water_model}.itp"\n')
            if use_ions_itp:
                f.write(f'#include "{ff_name}.ff/{use_ions_itp}"\n\n')
            else:
                # No ion file found - write comment and let grompp fail with clear error
                f.write(f"; WARNING: No ion file found for {ff_name} with water model {water_model}\n")
                f.write(f'#include "{ff_name}.ff/ions.itp"  ; This will likely fail\n\n')
            f.write("[ system ]\n")
            f.write("Ligand in water\n\n")
            f.write("[ molecules ]\n")
            # Read the actual moleculetype name from the ligand ITP file
            # This handles the multi-ligand renaming (LIG -> LIG_N) correctly
            ligand_itp = Path(self._resolve_ligand_ff_artifact(str(ligand_ff_dir), "LIG.itp"))
            ligand_mol_name = self._read_moleculetype_name(ligand_itp)
            f.write(f"{ligand_mol_name:<8} 1\n")

        shutil.copy(ligand_gro, model_dir / "lig.gro")

        print_step(3, 5, "Creating simulation box")
        boxed_gro = model_dir / "lig_newbox.gro"
        if box_size:
            print(f"  Using specified box size: {box_size[0]:.3f} x {box_size[1]:.3f} x {box_size[2]:.3f} nm")
            self.system_builder._run_command(
                [
                    self.system_builder.gmx_command,
                    "editconf",
                    "-f",
                    str(model_dir / "lig.gro"),
                    "-o",
                    str(boxed_gro),
                    "-bt",
                    "cubic",
                    "-box",
                    str(box_size[0]),
                    str(box_size[1]),
                    str(box_size[2]),
                    "-c",
                ],
                work_dir=str(model_dir),
            )
        else:
            box_distance = self.config.get("system", {}).get("box_distance", 1.5)
            self.system_builder._run_command(
                [
                    self.system_builder.gmx_command,
                    "editconf",
                    "-f",
                    str(model_dir / "lig.gro"),
                    "-o",
                    str(boxed_gro),
                    "-bt",
                    "cubic",
                    "-d",
                    str(box_distance),
                    "-c",
                ],
                work_dir=str(model_dir),
            )

        print_step(4, 5, "Solvating system")
        solvated_gro = self.system_builder._solvate(str(boxed_gro), str(topol_path))

        print_step(5, 5, "Adding ions")
        self.system_builder._add_ions(solvated_gro, str(topol_path))

        print_success(f"Ligand-only system built in {model_dir}")

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

    def _resolve_ligand_ff_artifact(self, ligand_ff_dir: str, filename: str) -> str:
        """Resolve an artifact within a generated ligand FF directory."""
        resolved = resolve_ligand_artifact(ligand_ff_dir, filename)
        if resolved is None:
            raise FileNotFoundError(f"Could not find {filename} under ligand FF directory: {ligand_ff_dir}")
        return str(resolved)

    def _read_moleculetype_name(self, itp_path: Path) -> str:
        """Read the moleculetype name from a ligand ITP file.

        Parameters
        ----------
        itp_path : Path
            Path to the LIG.itp file.

        Returns
        -------
        str
            The moleculetype name (e.g., LIG or LIG_1).

        Raises
        ------
        ValueError
            If the moleculetype section cannot be parsed.
        """
        in_moleculetype = False
        for raw_line in itp_path.read_text().splitlines():
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.lower() == "[ moleculetype ]":
                in_moleculetype = True
                continue
            if in_moleculetype:
                # First non-comment line after [ moleculetype ]
                # Format: "name  nrexcl"
                return line.split()[0]
        raise ValueError(f"Could not parse moleculetype from {itp_path}")

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
