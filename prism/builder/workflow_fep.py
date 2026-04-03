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
from typing import Optional, Tuple

from ..utils.colors import (
    print_header,
    print_step,
    print_success,
    print_error,
    print_warning,
    path,
)


class FEPWorkflowMixin:
    """Mixin providing FEP workflow for relative binding free energy calculations."""

    def run_fep(self):
        """Run the FEP workflow: build standard MD systems, generate hybrid topology, create FEP scaffold"""
        from ..fep.modeling import FEPScaffoldBuilder
        from ..fep.gromacs.itp_builder import ITPBuilder

        print_header("PRISM FEP Builder Workflow")

        if not self.mutant_ligand:
            raise ValueError("FEP mode requires --mutant ligand file")

        try:
            # Create FEP output directory first
            # Avoid double-nesting: if output_dir already ends with GMX_PROLIG_FEP, use it directly
            if os.path.basename(os.path.abspath(self.output_dir)) == "GMX_PROLIG_FEP":
                fep_output = os.path.abspath(self.output_dir)
            else:
                fep_output = os.path.join(self.output_dir, "GMX_PROLIG_FEP")
            os.makedirs(fep_output, exist_ok=True)

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

            # Phase 3: Generate hybrid topology
            print_step(3, 5, "Generating hybrid topology via atom mapping")

            # ref_ff_dir was saved in Phase 1
            # Generate mutant ligand FF in _build directory
            mut_ff_output = os.path.join(fep_output, "_build/mutant_ligand_ff")
            self._generate_mutant_ligand_ff(self.mutant_ligand, mut_ff_output)
            mut_ff_dir = self._resolve_generated_ligand_ff_dir(mut_ff_output)

            # Build hybrid topology using ITPBuilder
            hybrid_output = os.path.join(fep_output, "common/hybrid")
            os.makedirs(hybrid_output, exist_ok=True)

            # Read ligand atoms from ITP and GRO files
            from ..fep.io import read_ligand_from_prism
            from ..fep.core.hybrid_topology import HybridTopologyBuilder, LigandTopologyInput

            # Debug: Check if files exist
            ref_itp = os.path.join(ref_ff_dir, "LIG.itp")
            ref_gro = os.path.join(ref_ff_dir, "LIG.gro")
            mut_itp = os.path.join(mut_ff_dir, "LIG.itp")
            mut_gro = os.path.join(mut_ff_dir, "LIG.gro")

            print(f"\n  Checking force field files:")
            print(f"    ref_itp exists: {os.path.exists(ref_itp)} - {ref_itp}")
            print(f"    ref_gro exists: {os.path.exists(ref_gro)} - {ref_gro}")
            print(f"    mut_itp exists: {os.path.exists(mut_itp)} - {mut_itp}")
            print(f"    mut_gro exists: {os.path.exists(mut_gro)} - {mut_gro}")

            ref_mapping_coord_source = ref_gro
            mut_mapping_coord_source = mut_gro
            if not os.path.exists(ref_mapping_coord_source):
                raise FileNotFoundError(f"Generated reference ligand coordinates not found: {ref_mapping_coord_source}")
            if not os.path.exists(mut_mapping_coord_source):
                raise FileNotFoundError(f"Generated mutant ligand coordinates not found: {mut_mapping_coord_source}")

            ref_visual_coord_source = self.ligand_paths[0] if self.ligand_paths else ref_gro
            if not os.path.exists(ref_visual_coord_source):
                ref_visual_coord_source = ref_gro

            mut_visual_coord_source = self.mutant_ligand or mut_gro
            if not os.path.exists(mut_visual_coord_source):
                mut_visual_coord_source = mut_gro

            print("  Mapping coordinate sources:")
            print(f"    Reference coords: {ref_mapping_coord_source}")
            print(f"    Mutant coords:    {mut_mapping_coord_source}")

            ref_atoms = read_ligand_from_prism(itp_file=ref_itp, gro_file=ref_mapping_coord_source)

            mut_atoms = read_ligand_from_prism(itp_file=mut_itp, gro_file=mut_mapping_coord_source)

            # Perform atom mapping
            from ..fep.core.mapping import DistanceAtomMapper

            mapper = DistanceAtomMapper(
                dist_cutoff=self.distance_cutoff,
                charge_cutoff=self.charge_cutoff,
                charge_common=self.charge_strategy,
                charge_reception=self.charge_reception,
            )
            mapping = mapper.map(ref_atoms, mut_atoms)

            # Export an interactive mapping report into the final FEP output
            # so CLI users can inspect the hybridization result without running
            # a separate visualization helper.
            self._generate_fep_mapping_html(
                output_dir=hybrid_output,
                mapping=mapping,
                ref_atoms=ref_atoms,
                mut_atoms=mut_atoms,
                ref_coord_source=ref_visual_coord_source,
                mut_coord_source=mut_visual_coord_source,
            )

            # Build hybrid topology
            hybrid_builder = HybridTopologyBuilder(charge_strategy=self.charge_strategy)

            # Read ITP parameters for hybrid topology building
            ref_itp_data = ITPBuilder._parse_source_itp(Path(os.path.join(ref_ff_dir, "LIG.itp")).read_text())
            mut_itp_data = ITPBuilder._parse_source_itp(Path(os.path.join(mut_ff_dir, "LIG.itp")).read_text())

            params_a = LigandTopologyInput(
                masses={},
                bonds=ref_itp_data["sections"].get("bonds", []),
                pairs=ref_itp_data["sections"].get("pairs", []),
                angles=ref_itp_data["sections"].get("angles", []),
                dihedrals=ref_itp_data["sections"].get("dihedrals", []),
                impropers=ref_itp_data["sections"].get("impropers", []),
            )
            params_b = LigandTopologyInput(
                masses={},
                bonds=mut_itp_data["sections"].get("bonds", []),
                pairs=mut_itp_data["sections"].get("pairs", []),
                angles=mut_itp_data["sections"].get("angles", []),
                dihedrals=mut_itp_data["sections"].get("dihedrals", []),
                impropers=mut_itp_data["sections"].get("impropers", []),
            )

            hybrid_atoms = hybrid_builder.build(mapping, params_a, params_b, ref_atoms, mut_atoms)

            # Build hybrid bonded parameters from source ITPs
            # First write a temporary atoms-only ITP
            temp_itp = os.path.join(hybrid_output, "hybrid_atoms_temp.itp")
            ITPBuilder(hybrid_atoms, {}).write_itp(temp_itp, molecule_name="HYB")
            print(f"  Debug: temp_itp written: {os.path.exists(temp_itp)}")

            # Then build complete hybrid ITP with bonded terms
            hybrid_itp = os.path.join(hybrid_output, "hybrid.itp")
            hybrid_params = ITPBuilder.write_complete_hybrid_itp(
                output_path=hybrid_itp,
                hybrid_itp=temp_itp,
                ligand_a_itp=os.path.join(ref_ff_dir, "LIG.itp"),
                ligand_b_itp=os.path.join(mut_ff_dir, "LIG.itp"),
                molecule_name="HYB",
            )
            print(f"  Debug: hybrid.itp written: {os.path.exists(hybrid_itp)}")
            if os.path.exists(hybrid_itp):
                print(f"  Debug: hybrid.itp size: {os.path.getsize(hybrid_itp)} bytes")

            print_success(f"Hybrid topology: {hybrid_itp}")

            # Generate hybrid ligand structure file
            hybrid_gro = os.path.join(hybrid_output, "hybrid.gro")
            self._generate_hybrid_gro(hybrid_gro, hybrid_atoms, mapping, ref_ff_dir, mut_ff_dir)
            print_success(f"Hybrid structure: {hybrid_gro}")

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
                hybrid_gro=hybrid_gro,
                reference_ligand_dir=ref_ff_dir,
                mutant_ligand_dir=mut_ff_dir,
                bound_system_dir=bound_system_dir,
                unbound_system_dir=unbound_system_dir,
            )

            print_success(f"FEP scaffold created: {fep_output}")

            # Phase 5: Clean up intermediate build files
            print_step(5, 5, "Cleaning up intermediate build files")

            build_dir = os.path.join(fep_output, "_build")
            if os.path.exists(build_dir):
                shutil.rmtree(build_dir)
                print(f"  ✓ Removed build directory: {build_dir}")

            # Also clean up any files in output_dir root (outside GMX_PROLIG_FEP)
            cleanup_items = [
                os.path.join(self.output_dir, "mdps"),
            ]
            cleanup_items.extend(str(p) for p in Path(self.output_dir).glob("LIG.*") if p.is_dir())

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

    def _generate_fep_mapping_html(
        self,
        output_dir: str,
        mapping,
        ref_atoms: list,
        mut_atoms: list,
        ref_coord_source: str,
        mut_coord_source: str,
    ) -> None:
        """Write a default mapping HTML report for CLI-based FEP builds."""
        try:
            from ..fep.visualize.html import visualize_mapping_html
        except Exception as exc:
            print_warning(f"Skipping FEP mapping HTML generation: {exc}")
            return

        html_path = os.path.join(output_dir, "mapping.html")
        ref_pdb, ref_mol2 = self._prepare_mapping_visualization_inputs(ref_coord_source, output_dir, "ligand_a")
        mut_pdb, mut_mol2 = self._prepare_mapping_visualization_inputs(mut_coord_source, output_dir, "ligand_b")

        html_config = None
        if self.fep_config:
            try:
                from ..fep.config import FEPConfig

                html_config = FEPConfig(self.config_file, self.fep_config).get_html_config()
            except Exception as exc:
                print_warning(f"Could not load FEP HTML config: {exc}")

        ligand_a_name = Path(ref_coord_source).stem or "Ligand A"
        ligand_b_name = Path(mut_coord_source).stem or "Ligand B"

        try:
            visualize_mapping_html(
                mapping=mapping,
                pdb_a=ref_pdb,
                pdb_b=mut_pdb,
                mol2_a=ref_mol2,
                mol2_b=mut_mol2,
                atoms_a=ref_atoms,
                atoms_b=mut_atoms,
                output_path=html_path,
                title=f"FEP Mapping: {ligand_a_name} -> {ligand_b_name}",
                ligand_a_name=ligand_a_name,
                ligand_b_name=ligand_b_name,
                config=html_config,
            )
            print_success(f"Mapping HTML: {html_path}")
        except Exception as exc:
            print_warning(f"Failed to generate FEP mapping HTML: {exc}")

    def _prepare_mapping_visualization_inputs(
        self, coord_source: str, output_dir: str, stem: str
    ) -> Tuple[str, Optional[str]]:
        """Return a PDB path plus optional MOL2 template for mapping HTML generation."""
        source_path = Path(coord_source)
        suffix = source_path.suffix.lower()

        if suffix == ".pdb":
            mol2_candidates = [
                source_path.with_suffix(".mol2"),
                source_path.with_name(f"{source_path.stem}_3D.mol2"),
                source_path.with_name(f"{source_path.stem.lower()}_3D.mol2"),
            ]
            for mol2_candidate in mol2_candidates:
                if mol2_candidate.exists():
                    return str(source_path), str(mol2_candidate)
            return str(source_path), None

        if suffix == ".mol2":
            try:
                from rdkit import Chem

                mol = Chem.MolFromMol2File(str(source_path), sanitize=False, removeHs=False)
                if mol is None:
                    raise ValueError(f"RDKit could not parse {source_path}")

                pdb_path = Path(output_dir) / f"{stem}.pdb"
                Chem.MolToPDBFile(mol, str(pdb_path))
                return str(pdb_path), str(source_path)
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to prepare PDB for mapping visualization from {source_path}: {exc}"
                ) from exc

        raise RuntimeError(f"Unsupported coordinate source for mapping visualization: {coord_source}")

    def _generate_hybrid_gro(
        self, output_gro: str, hybrid_atoms: list, _mapping: "AtomMapping", ref_ff_dir: str, mut_ff_dir: str
    ) -> None:
        """
        Generate hybrid ligand GRO file with dummy atoms

        Parameters
        ----------
        output_gro : str
            Output path for hybrid GRO file
        hybrid_atoms : list
            List of HybridAtom objects from hybrid topology
        mapping : AtomMapping
            Atom mapping between reference and mutant ligands
        ref_ff_dir : str
            Reference ligand force field directory (contains LIG.gro)
        mut_ff_dir : str
            Mutant ligand force field directory (contains LIG.gro)
        """

        # Read reference and mutant ligand GRO files
        ref_gro_path = os.path.join(ref_ff_dir, "LIG.gro")
        mut_gro_path = os.path.join(mut_ff_dir, "LIG.gro")

        ref_coords = self._parse_gro(ref_gro_path)
        mut_coords = self._parse_gro(mut_gro_path)

        # Build coordinate lookup by atom name
        ref_coord_map = {atom["name"]: atom["coord"] for atom in ref_coords["atoms"]}
        mut_coord_map = {atom["name"]: atom["coord"] for atom in mut_coords["atoms"]}

        # Generate hybrid coordinates
        hybrid_atoms_data = []
        for hatom in hybrid_atoms:
            # Remove A/B suffixes from hybrid atom name to find original name
            base_name = hatom.name.rstrip("AB")

            if hatom.name.endswith("A"):
                # Transformed A atom (disappearing in state B)
                # Use coordinates from reference ligand
                coord = ref_coord_map.get(base_name, [0.0, 0.0, 0.0])
            elif hatom.name.endswith("B"):
                # Transformed B atom (appearing in state B)
                # Use coordinates from mutant ligand
                coord = mut_coord_map.get(base_name, [0.0, 0.0, 0.0])
            else:
                # Common atom
                # Use coordinates from reference ligand
                coord = ref_coord_map.get(hatom.name, [0.0, 0.0, 0.0])

            hybrid_atoms_data.append({"index": hatom.index, "name": hatom.name, "coord": coord})

        # Write hybrid GRO file
        with open(output_gro, "w") as f:
            f.write("Hybrid ligand structure generated by PRISM-FEP\n")
            f.write(f"{len(hybrid_atoms_data)}\n")

            # Get residue number from reference ligand
            if ref_coords["atoms"]:
                # Parse the first atom line to get residue number
                with open(ref_gro_path, "r") as ref_file:
                    ref_lines = ref_file.readlines()
                    first_atom_line = ref_lines[2]  # Skip title and atom count
                    res_num_str = first_atom_line[0:5].strip()
                    residue_number = int(res_num_str) if res_num_str else 1
            else:
                residue_number = 1

            for atom_data in hybrid_atoms_data:
                idx = atom_data["index"]
                name = atom_data["name"]
                x, y, z = atom_data["coord"]
                # GRO format: %5d%-5s%5s%5d%8.3f%8.3f%8.3f
                # residuenum (5 chars) + residuename (5 chars, left-justified)
                # + atomname (5 chars, right-justified) + atomnumber (5 chars) + x + y + z
                f.write(f"{residue_number:5d}LIG  {name:>5s}{idx:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")

            # Write box vectors (use default 1.0 nm, will be replaced by system box)
            f.write("   1.00000   1.00000   1.00000\n")

    def _parse_gro(self, gro_file: str) -> dict:
        """Parse GRO file and extract atoms and coordinates"""
        atoms = []
        with open(gro_file, "r") as f:
            lines = f.readlines()

        # Skip title line
        # Second line is atom count
        num_atoms = int(lines[1].strip())

        # Parse atom lines (format: %5d%-5s%5s%5d%8.3f%8.3f%8.3f)
        for i in range(num_atoms):
            line = lines[2 + i]
            # Skip box line (starts with space and has only numbers)
            if i == num_atoms - 1 and len(line.strip().split()) <= 3:
                break

            # More robust parsing using split
            parts = line.split()
            if len(parts) < 7:
                continue

            # GRO format: residue_number residue_name atom_name atom_number x y z
            # parts: [0]=residue_number, [1]=residue_name, [2]=atom_name, [3]=atom_number, [4]=x, [5]=y, [6]=z
            atom_name = parts[2].strip()
            try:
                x = float(parts[4])
                y = float(parts[5])
                z = float(parts[6])
                atoms.append({"name": atom_name, "coord": [x, y, z]})
            except (ValueError, IndexError):
                # Try fixed column format as fallback
                try:
                    atom_name = line[5:10].strip()
                    x = float(line[20:28].strip())
                    y = float(line[28:36].strip())
                    z = float(line[36:44].strip())
                    atoms.append({"name": atom_name, "coord": [x, y, z]})
                except (ValueError, IndexError):
                    continue

        return {"atoms": atoms}

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
        ligand_gro = ligand_ff_dir / "LIG.gro"
        if not ligand_gro.exists():
            raise FileNotFoundError(f"Ligand structure not found: {ligand_gro}")

        print_step(2, 5, "Creating ligand-only topology")

        # Copy force field to working directory (same as bound system)
        ff_name = self.forcefield["name"]
        water_model = self.water_model["name"]
        ff_idx = self.forcefield.get("index")
        ff_info = None
        use_ions_itp = "ions.itp"  # Default to old-style

        if ff_idx and self.gromacs_env:
            # Get force field info from environment
            if ff_idx in self.gromacs_env.force_fields:
                ff_info = self.gromacs_env.force_fields[ff_idx]

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

                # Check for ion files (support both old and new force field structures)
                water_specific_ions = f"ions_{water_model}.itp"
                if (ff_path / water_specific_ions).exists():
                    required_files.append(water_specific_ions)
                    use_ions_itp = water_specific_ions
                elif (ff_path / "ions.itp").exists():
                    # Old-style force field with generic ions.itp
                    required_files.append("ions.itp")
                    use_ions_itp = "ions.itp"
                else:
                    # No ion file found - will raise error below
                    use_ions_itp = None

                missing = [f for f in required_files if not (ff_path / f).exists()]
                if missing:
                    raise RuntimeError(
                        f"Force field is incomplete: {ff_info['name']}\n"
                        f"Missing files: {', '.join(missing)}\n"
                        f"Force field path: {ff_path}\n"
                        f"Source: {ff_info.get('source', 'unknown')}\n"
                        f"Please check your force field installation."
                    )

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

        topol_path = model_dir / "topol.top"
        ligand_rel_dir = f"../{ligand_ff_dir.name}"

        with open(topol_path, "w") as f:
            f.write("; Ligand-only topology for FEP unbound leg\n")
            f.write(f'#include "{ff_name}.ff/forcefield.itp"\n')
            f.write(f'#include "{ligand_rel_dir}/atomtypes_LIG.itp"\n')
            f.write(f'#include "{ligand_rel_dir}/LIG.itp"\n')
            f.write(f'#include "{ff_name}.ff/{water_model}.itp"\n')
            f.write(f'#include "{ff_name}.ff/{use_ions_itp}"\n\n')
            f.write("[ system ]\n")
            f.write("Ligand in water\n\n")
            f.write("[ molecules ]\n")
            f.write("LIG    1\n")

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
        ff_dirs = [p for p in Path(output_dir).glob("LIG.*") if p.is_dir()]
        if len(ff_dirs) != 1:
            raise FileNotFoundError(
                f"Expected exactly one generated ligand FF dir in {output_dir}, found {len(ff_dirs)}"
            )
        return str(ff_dirs[0])

    def _generate_mutant_ligand_ff(self, mutant_ligand: str, output_dir: str):
        """Generate force field for mutant ligand"""
        # Reuse existing force field generation logic

        # Temporarily save current state
        original_ligands = self.ligand_paths
        original_output = self.output_dir
        original_lig_ff_dirs = self.lig_ff_dirs

        # Set up for mutant ligand
        self.ligand_paths = [mutant_ligand]
        self.output_dir = output_dir

        try:
            # Generate FF
            self.generate_ligand_forcefield()
        finally:
            # Restore shared builder state even if FF generation fails.
            self.ligand_paths = original_ligands
            self.output_dir = original_output
            self.lig_ff_dirs = original_lig_ff_dirs
