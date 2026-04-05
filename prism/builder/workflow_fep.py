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
from ..utils.system.topology import _parse_atomtypes

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

            ref_visual_coord_source = self.ligand_paths[0] if self.ligand_paths else ref_gro
            if not os.path.exists(ref_visual_coord_source):
                ref_visual_coord_source = ref_gro

            mut_visual_coord_source = self.mutant_ligand or mut_gro
            if not os.path.exists(mut_visual_coord_source):
                mut_visual_coord_source = mut_gro

            # Use user-provided, receptor-aligned ligand coordinates for mapping whenever
            # available. Generated LIG.gro coordinates can differ from the aligned input
            # geometry and degrade the atom mapping and HTML report quality.
            ref_mapping_coord_source = ref_visual_coord_source if os.path.exists(ref_visual_coord_source) else ref_gro
            mut_mapping_coord_source = mut_visual_coord_source if os.path.exists(mut_visual_coord_source) else mut_gro
            if not os.path.exists(ref_mapping_coord_source):
                raise FileNotFoundError(f"Reference ligand coordinates not found: {ref_mapping_coord_source}")
            if not os.path.exists(mut_mapping_coord_source):
                raise FileNotFoundError(f"Mutant ligand coordinates not found: {mut_mapping_coord_source}")

            print("  Mapping coordinate sources:")
            print(f"    Reference coords: {ref_mapping_coord_source}")
            print(f"    Mutant coords:    {mut_mapping_coord_source}")

            ref_atoms = read_ligand_from_prism(itp_file=ref_itp, gro_file=ref_mapping_coord_source, state="a")

            mut_atoms = read_ligand_from_prism(itp_file=mut_itp, gro_file=mut_mapping_coord_source, state="b")

            # Perform atom mapping
            from ..fep.core.mapping import DistanceAtomMapper

            mapper = DistanceAtomMapper(
                dist_cutoff=self.distance_cutoff,
                charge_cutoff=self.charge_cutoff,
                charge_common=self.charge_strategy,
                charge_reception=self.charge_reception,
            )
            mapping = mapper.map(ref_atoms, mut_atoms)

            # Debug: Print mapping results
            print(f"  Debug: Mapping results:")
            print(f"    Common: {len(mapping.common)}")
            print(f"    Transformed A: {len(mapping.transformed_a)}")
            print(f"    Transformed B: {len(mapping.transformed_b)}")
            print(f"    Surrounding A: {len(mapping.surrounding_a)}")
            print(f"    Surrounding B: {len(mapping.surrounding_b)}")
            if mapping.transformed_a:
                print(f"    Transformed A atoms:")
                for atom in mapping.transformed_a:
                    print(f"      {atom.name}: type={atom.atom_type}, charge={atom.charge}")
            if mapping.transformed_b:
                print(f"    Transformed B atoms:")
                for atom in mapping.transformed_b:
                    print(f"      {atom.name}: type={atom.atom_type}, charge={atom.charge}")

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
            print(f"  Debug: hybrid_atoms count = {len(hybrid_atoms)}")
            transformed_atoms = [a for a in hybrid_atoms if a.name.endswith("A") or a.name.endswith("B")]
            print(f"  Debug: transformed atoms with A/B suffix: {len(transformed_atoms)}")
            for atom in transformed_atoms:
                print(f"    {atom.name}: A_type={atom.state_a_type}, B_type={atom.state_b_type}")

            # Build hybrid bonded parameters from source ITPs
            # First write a temporary atoms-only ITP
            temp_itp = os.path.join(hybrid_output, "hybrid_atoms_temp.itp")
            ITPBuilder(hybrid_atoms, {}).write_itp(temp_itp, molecule_name="HYB")
            print(f"  Debug: temp_itp written: {os.path.exists(temp_itp)}")

            # Then build complete hybrid ITP with bonded terms
            hybrid_itp = os.path.join(hybrid_output, "hybrid.itp")
            # Extract atom types for CHARMM force field dihedral handling
            from prism.fep.modeling.hybrid_package import HybridPackageBuilder

            hybrid_pkg_builder = HybridPackageBuilder()
            atom_types_a = hybrid_pkg_builder._parse_atom_types_from_itp(os.path.join(ref_ff_dir, "LIG.itp"))
            atom_types_b = hybrid_pkg_builder._parse_atom_types_from_itp(os.path.join(mut_ff_dir, "LIG.itp"))

            hybrid_params = ITPBuilder.write_complete_hybrid_itp(
                output_path=hybrid_itp,
                hybrid_itp=temp_itp,
                ligand_a_itp=os.path.join(ref_ff_dir, "LIG.itp"),
                ligand_b_itp=os.path.join(mut_ff_dir, "LIG.itp"),
                molecule_name="HYB",
                atom_types_a=atom_types_a,
                atom_types_b=atom_types_b,
            )
            print(f"  Debug: hybrid.itp written: {os.path.exists(hybrid_itp)}")
            if os.path.exists(hybrid_itp):
                print(f"  Debug: hybrid.itp size: {os.path.getsize(hybrid_itp)} bytes")

            print_success(f"Hybrid topology: {hybrid_itp}")

            # Generate hybrid ligand structure file
            hybrid_gro = os.path.join(hybrid_output, "hybrid.gro")
            self._generate_hybrid_gro(hybrid_gro, hybrid_atoms, mapping, ref_ff_dir, mut_ff_dir)
            print_success(f"Hybrid structure: {hybrid_gro}")

            # Export an interactive mapping report into the final FEP output
            # so CLI users can inspect the hybridization result without running
            # a separate visualization helper.
            # Note: Called after hybrid topology creation to ensure hybrid.itp and hybrid.gro exist
            self._generate_fep_mapping_html(
                output_dir=hybrid_output,
                ref_coord_source=ref_visual_coord_source,
                mut_coord_source=mut_visual_coord_source,
                mapping=mapping,
                ref_atoms=ref_atoms,
                mut_atoms=mut_atoms,
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
        ref_coord_source: str,
        mut_coord_source: str,
        mapping,
        ref_atoms,
        mut_atoms,
    ) -> None:
        """Write a default mapping HTML report for CLI-based FEP builds."""
        try:
            from ..fep.visualize.html import visualize_mapping_html
        except Exception as exc:
            print_warning(f"Skipping FEP mapping HTML generation: {exc}")
            return

        html_path = os.path.join(output_dir, "mapping.html")

        # Default: use the original mapped ligand atoms as the single source of
        # truth for classification. Do not recompute mapping from hybrid state
        # atoms here.
        from ..fep.io import read_mol2_atoms, write_ligand_to_pdb

        html_mapping = mapping
        html_atoms_a = ref_atoms
        html_atoms_b = mut_atoms

        ref_pdb = os.path.join(output_dir, "mapping_state_a.pdb")
        mut_pdb = os.path.join(output_dir, "mapping_state_b.pdb")
        write_ligand_to_pdb(html_atoms_a, ref_pdb, residue_name="LIG")
        write_ligand_to_pdb(html_atoms_b, mut_pdb, residue_name="LIG")

        print(f"  DEBUG MOL2 detection:")
        print(f"    ref_coord_source: {ref_coord_source}")
        print(f"    mut_coord_source: {mut_coord_source}")
        print(f"    ref_coord_source suffix: {Path(ref_coord_source).suffix}")
        print(f"    mut_coord_source suffix: {Path(mut_coord_source).suffix}")

        if Path(ref_coord_source).suffix == ".mol2":
            ref_mol2 = str(Path(ref_coord_source).absolute()) if Path(ref_coord_source).exists() else ref_coord_source
        elif Path(ref_coord_source).suffix == ".pdb":
            ref_mol2 = str(Path(ref_coord_source).with_suffix(".mol2"))
            if not Path(ref_mol2).exists():
                ref_mol2 = None
        else:
            ref_mol2 = None

        if Path(mut_coord_source).suffix == ".mol2":
            mut_mol2 = str(Path(mut_coord_source).absolute()) if Path(mut_coord_source).exists() else mut_coord_source
        elif Path(mut_coord_source).suffix == ".pdb":
            mut_mol2 = str(Path(mut_coord_source).with_suffix(".mol2"))
            if not Path(mut_mol2).exists():
                mut_mol2 = None
        else:
            mut_mol2 = None

        print(f"    ref_mol2: {ref_mol2}")
        print(f"    mut_mol2: {mut_mol2}")

        from ..fep.core.mapping import DistanceAtomMapper

        html_mapper = DistanceAtomMapper(
            dist_cutoff=self.distance_cutoff,
            charge_cutoff=self.charge_cutoff,
            charge_common=self.charge_strategy,
            charge_reception=self.charge_reception,
        )

        # For generic sequential force fields such as OPLS/OpenFF, the original
        # aligned MOL2 files are a better source for HTML mapping than the
        # parameterized topology output. If the MOL2-derived mapping is clearly
        # better, use it for visualization only.
        if ref_mol2 and mut_mol2:
            try:
                mol2_atoms_a = read_mol2_atoms(ref_mol2)
                mol2_atoms_b = read_mol2_atoms(mut_mol2)
                mol2_mapping = html_mapper.map(mol2_atoms_a, mol2_atoms_b)
                if len(mol2_mapping.common) > len(html_mapping.common):
                    html_mapping = mol2_mapping
                    html_atoms_a = mol2_atoms_a
                    html_atoms_b = mol2_atoms_b
                    write_ligand_to_pdb(html_atoms_a, ref_pdb, residue_name="LIG")
                    write_ligand_to_pdb(html_atoms_b, mut_pdb, residue_name="LIG")
                    print(
                        "  Using MOL2-derived mapping for HTML "
                        f"({len(mol2_mapping.common)} common vs {len(mapping.common)})"
                    )
            except Exception as exc:
                print_warning(f"Could not build MOL2-derived HTML mapping: {exc}")
        html_config = None
        if self.fep_config:
            try:
                from ..fep.config import FEPConfig

                html_config = FEPConfig(self.config_file, self.fep_config).get_html_config()
            except Exception as exc:
                print_warning(f"Could not load FEP HTML config: {exc}")

        if html_config is None:
            html_config = {}

        ff_cfg = dict(html_config.get("forcefield", {}))
        ff_type = ff_cfg.get("type", getattr(self, "ligand_forcefield", "unknown"))
        ff_cfg["type"] = self._describe_mapping_forcefield(ff_type)
        html_config["forcefield"] = ff_cfg

        ligand_a_name = Path(ref_coord_source).stem or "Ligand A"
        ligand_b_name = Path(mut_coord_source).stem or "Ligand B"

        try:
            visualize_mapping_html(
                mapping=html_mapping,
                pdb_a=ref_pdb,
                pdb_b=mut_pdb,
                mol2_a=ref_mol2,
                mol2_b=mut_mol2,
                atoms_a=html_atoms_a,
                atoms_b=html_atoms_b,
                output_path=html_path,
                title=f"FEP Mapping: {ligand_a_name} -> {ligand_b_name}",
                ligand_a_name=ligand_a_name,
                ligand_b_name=ligand_b_name,
                config=html_config,
            )
            print_success(f"Mapping HTML: {html_path}")
        except Exception as exc:
            print_warning(f"Failed to generate FEP mapping HTML: {exc}")

    def _describe_mapping_forcefield(self, ff_type: str) -> str:
        """Return a user-facing force-field label for mapping HTML."""
        base = str(ff_type).upper()
        lig_ff_lower = getattr(self, "ligand_forcefield", "").lower()

        # Handle both "cgenff" and "charmm-gui" force field types
        if lig_ff_lower not in ["cgenff", "charmm-gui"]:
            return base

        ff_paths = [Path(p) for p in (getattr(self, "forcefield_paths", None) or [])]
        if ff_paths and all((path / "charmm36.ff" / "charmm36.itp").exists() for path in ff_paths):
            return "CGENFF (CHARMM-GUI)"
        elif ff_paths and all((path / "gromacs" / "LIG.itp").exists() for path in ff_paths):
            return "CGENFF (CHARMM-GUI)"
        return "CGENFF (WEBSITE)"

    def _resolve_cgenff_mapping_topology(self, cgenff_path: str) -> str:
        """Return the source topology file that best represents the original CGenFF ligand."""
        path = Path(cgenff_path)
        if (path / "LIG.itp").exists():
            return str(path / "LIG.itp")
        gmx_tops = sorted(path.glob("*_gmx.top"))
        if gmx_tops:
            return str(gmx_tops[0])
        top_files = [candidate for candidate in sorted(path.glob("*.top")) if candidate.name != "LIG.top"]
        if top_files:
            return str(top_files[0])
        if (path / "LIG.top").exists() and (path / "LIG.itp").exists():
            return str(path / "LIG.itp")
        raise FileNotFoundError(f"No suitable CGenFF topology file found under {cgenff_path}")

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

        if suffix == ".gro":
            # Convert GRO to PDB using simple conversion (MDAnalysis writes incorrect elements)
            # Force custom conversion to ensure proper element symbols
            pdb_path = Path(output_dir) / f"{stem}.pdb"
            with open(source_path) as f:
                lines = f.readlines()
            with open(pdb_path, "w") as out:
                # Skip header (first 2 lines) and footer (last line)
                for i, line in enumerate(lines[2:-1]):
                    parts = line.split()
                    if len(parts) >= 5:
                        # GRO format: resid_resname atom_name atom_id x y z
                        # resid and resname are combined (e.g., "1HYB")
                        resid_resname = parts[0]
                        atom_name = parts[1]
                        atom_id = parts[2]
                        x = float(parts[3]) * 10  # nm to Å
                        y = float(parts[4]) * 10
                        z = float(parts[5]) * 10 if len(parts) > 5 else 0.0
                        # Split resid and resname
                        # resid is the numeric part, resname is the alphabetic part
                        resid = ""
                        resname = ""
                        for char in resid_resname:
                            if char.isdigit():
                                resid += char
                            else:
                                resname += char
                        # Extract element symbol BEFORE modifying resname
                        # Element is first character (uppercase), or first two if second is lowercase
                        atom_element = atom_name[0].upper()
                        if len(atom_name) > 1 and atom_name[1].islower():
                            atom_element = atom_name[:2].capitalize()
                        # Handle special cases
                        if atom_name.upper().startswith("CL"):
                            atom_element = "Cl"
                        elif atom_name.upper().startswith("BR"):
                            atom_element = "Br"
                        # Write PDB format
                        out.write(
                            f"ATOM  {atom_id:>5} {atom_name:<4} {resname:>3} {resid:>4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom_element:>2}\n"
                        )
                out.write("END\n")
            return str(pdb_path), None

        raise RuntimeError(f"Unsupported coordinate source for mapping visualization: {coord_source}")

    def _generate_state_specific_pdb(self, gro_file: str, output_dir: str, stem: str, state: str = "a") -> str:
        """
        Generate state-specific PDB file from hybrid GRO file.

        Filters atoms by state:
        - State "a": Includes common atoms + transformed atoms ending with digit+"A" (e.g., "00A", "01A")
        - State "b": Includes common atoms + transformed atoms ending with digit+"B" (e.g., "00B", "01B")

        Common atoms don't follow the digit+A/B pattern (e.g., C0A, C02, H11).

        Args:
            gro_file: Path to hybrid GRO file
            output_dir: Output directory
            stem: Output file stem
            state: "a" for state A, "b" for state B

        Returns:
            Path to generated PDB file
        """

        pdb_path = Path(output_dir) / f"{stem}.pdb"
        with open(gro_file) as f:
            lines = f.readlines()

        with open(pdb_path, "w") as out:
            # Skip header (first 2 lines) and footer (last line)
            for i, line in enumerate(lines[2:-1]):
                parts = line.split()
                if len(parts) >= 5:
                    # GRO format: resid_resname atom_name atom_id x y z
                    resid_resname = parts[0]
                    atom_name = parts[1]
                    atom_id = parts[2]
                    x = float(parts[3]) * 10  # nm to Å
                    y = float(parts[4]) * 10
                    z = float(parts[5]) * 10 if len(parts) > 5 else 0.0

                    # Filter by state
                    # Transformed atoms end with "A" or "B" (e.g., "00A", "01B", "C0EA", "O0MB")
                    # Common atoms don't end with "A" or "B", except for the special case "C0A"
                    # Note: "C0A" is a common atom even though it ends with "A"
                    is_transformed_a = atom_name.endswith("A") and atom_name != "C0A"
                    is_transformed_b = atom_name.endswith("B")
                    is_common = not (is_transformed_a or is_transformed_b)

                    if state == "a":
                        # State A: include common atoms + transformed A atoms
                        if not (is_common or is_transformed_a):
                            continue
                    else:  # state == "b"
                        # State B: include common atoms + transformed B atoms
                        if not (is_common or is_transformed_b):
                            continue

                    # Split resid and resname
                    resid = ""
                    resname = ""
                    for char in resid_resname:
                        if char.isdigit():
                            resid += char
                        else:
                            resname += char

                    # Extract element symbol
                    atom_element = atom_name[0].upper()
                    if len(atom_name) > 1 and atom_name[1].islower():
                        atom_element = atom_name[:2].capitalize()
                    # Handle special cases
                    if atom_name.upper().startswith("CL"):
                        atom_element = "Cl"
                    elif atom_name.upper().startswith("BR"):
                        atom_element = "Br"

                    # Write PDB format
                    out.write(
                        f"ATOM  {atom_id:>5} {atom_name:<4} {resname:>3} {resid:>4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom_element:>2}\n"
                    )
            out.write("END\n")

        return str(pdb_path)

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
            # Remove single A/B suffix from hybrid atom name to find original name
            # Note: Only remove ONE suffix character (A or B), not all trailing As/Bs
            if hatom.name.endswith("A"):
                base_name = hatom.name[:-1]  # Remove single 'A'
                # Transformed A atom (disappearing in state B)
                # Use coordinates from reference ligand
                coord = ref_coord_map.get(base_name, [0.0, 0.0, 0.0])
            elif hatom.name.endswith("B"):
                base_name = hatom.name[:-1]  # Remove single 'B'
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

        # Parse atom lines using GRO fixed-width columns.
        for i in range(num_atoms):
            line = lines[2 + i]
            try:
                atom_name = line[10:15].strip()
                x = float(line[20:28].strip())
                y = float(line[28:36].strip())
                z = float(line[36:44].strip())
                atoms.append({"name": atom_name, "coord": [x, y, z]})
            except (ValueError, IndexError):
                # Fall back to split parsing for non-standard whitespace-separated lines
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    atom_name = parts[1].strip()
                    x = float(parts[3])
                    y = float(parts[4])
                    z = float(parts[5])
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
                lig_itp_path=ligand_ff_dir / "LIG.itp",
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
            atomtypes_itp = ligand_ff_dir / "atomtypes_LIG.itp"
            include_ligand_atomtypes = CGenFFAdapter.should_include_ligand_atomtypes(ff_name, ligand_ff_dir)
            if include_ligand_atomtypes and atomtypes_itp.exists():
                for line in atomtypes_itp.read_text().splitlines():
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
            ligand_itp = ligand_ff_dir / "LIG.itp"
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
        # First search top-level (backward compatibility)
        ff_dirs = [p for p in Path(output_dir).glob("LIG.*") if p.is_dir()]
        if len(ff_dirs) != 1:
            # If not found, search Ligand_Forcefield subdirectory (new multi-ligand format)
            ligand_ff_dir = Path(output_dir) / "Ligand_Forcefield"
            if ligand_ff_dir.exists():
                ff_dirs = sorted([p for p in ligand_ff_dir.glob("LIG.*") if p.is_dir()], key=lambda p: p.name)
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
