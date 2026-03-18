#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Scaffold builders for PRISM-FEP modeling workflows.

This module does not run the full alchemical system build yet. It stages the
directory layout, hybrid ligand assets, seed coordinate files, topology
templates, and per-leg MDP files so the later full integration can plug into a
stable on-disk structure.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional

from prism.fep.gromacs.itp_builder import ITPBuilder
from prism.fep.gromacs.mdp_templates import write_fep_mdps


@dataclass
class FEPScaffoldLayout:
    """Filesystem layout for a staged FEP modeling workspace."""

    root: Path
    common_dir: Path
    hybrid_dir: Path
    protein_dir: Path
    bound_dir: Path
    unbound_dir: Path


class FEPScaffoldBuilder:
    """
    Stage a PRISM-style FEP workspace around an existing hybrid ligand topology.

    Inputs are intentionally minimal:
    - receptor PDB
    - hybrid ligand directory containing at least `LIG.itp` and `LIG.gro`
    - output directory

    The resulting layout separates `bound` and `unbound` legs and keeps shared
    assets in `common/`.
    """

    def __init__(
        self,
        output_dir: str,
        config: Optional[Dict] = None,
        lambda_windows: int = 32,
        lambda_strategy: str = "decoupled",
        lambda_distribution: str = "nonlinear",
        overwrite: bool = False,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.config = config or {}
        self.lambda_windows = lambda_windows
        self.lambda_strategy = lambda_strategy
        self.lambda_distribution = lambda_distribution
        self.overwrite = overwrite
        self.hybrid_itp_filename = "LIG.itp"
        self.hybrid_atomtypes_filename = "atomtypes_LIG.itp"
        self.hybrid_forcefield_filename = "ff_hybrid.itp"
        self.hybrid_posre_filename = "posre_LIG.itp"
        self.hybrid_gro_filename = "LIG.gro"
        self.molecule_name = "LIG"

    def build(self, receptor_pdb: str, hybrid_ligand_dir: str) -> FEPScaffoldLayout:
        """Create the FEP workspace scaffold."""
        receptor_path = Path(receptor_pdb)
        hybrid_dir = Path(hybrid_ligand_dir)
        self._validate_inputs(receptor_path, hybrid_dir)

        layout = self._prepare_layout()
        self._copy_common_assets(receptor_path, hybrid_dir, layout)

        ligand_seed_pdb = layout.common_dir / "hybrid" / "ligand_seed.pdb"
        self._write_bound_leg(layout, receptor_path, ligand_seed_pdb)
        self._write_unbound_leg(layout, ligand_seed_pdb)
        self._write_root_scripts(layout)
        self._write_manifest(layout, receptor_path, hybrid_dir)

        return layout

    def build_from_components(
        self,
        receptor_pdb: str,
        hybrid_itp: str,
        reference_ligand_dir: str,
        mutant_ligand_dir: Optional[str] = None,
        hybrid_gro: Optional[str] = None,
        bound_system_dir: Optional[str] = None,
        unbound_system_dir: Optional[str] = None,
    ) -> FEPScaffoldLayout:
        """
        Create the scaffold from an existing hybrid ITP plus ligand FF directories.

        This is the intended near-term bridge for cases where atom mapping and
        single-topology construction are already done and only the system/leg
        modeling scaffold is still missing.

        Parameters:
        -----------
        receptor_pdb : str
            Path to receptor PDB file
        hybrid_itp : str
            Path to hybrid ligand ITP file
        reference_ligand_dir : str
            Path to reference ligand force field directory
        mutant_ligand_dir : Optional[str]
            Path to mutant ligand force field directory
        hybrid_gro : Optional[str]
            Path to hybrid ligand GRO file (defaults to reference_ligand_dir/LIG.gro)
        bound_system_dir : Optional[str]
            Path to pre-built PRISM bound system (GMX_PROLIG_MD directory).
            If provided, copies complete system instead of creating placeholder.
        unbound_system_dir : Optional[str]
            Path to pre-built PRISM unbound system (GMX_PROLIG_MD directory).
            If provided, copies complete system instead of creating placeholder.
        """
        receptor_path = Path(receptor_pdb)
        hybrid_itp_path = Path(hybrid_itp)
        ref_dir = Path(reference_ligand_dir)
        mut_dir = Path(mutant_ligand_dir) if mutant_ligand_dir else None

        if not receptor_path.exists():
            raise FileNotFoundError(f"Receptor PDB not found: {receptor_path}")
        if not hybrid_itp_path.exists():
            raise FileNotFoundError(f"Hybrid ITP not found: {hybrid_itp_path}")
        if not ref_dir.exists():
            raise FileNotFoundError(f"Reference ligand FF dir not found: {ref_dir}")
        if mut_dir and not mut_dir.exists():
            raise FileNotFoundError(f"Mutant ligand FF dir not found: {mut_dir}")

        ref_gro = Path(hybrid_gro) if hybrid_gro else ref_dir / "LIG.gro"
        if not ref_gro.exists():
            raise FileNotFoundError(f"Hybrid/seed GRO not found: {ref_gro}")

        layout = self._prepare_layout()
        shutil.copy2(receptor_path, layout.protein_dir / "receptor.pdb")
        self._prepare_hybrid_package_from_components(layout, hybrid_itp_path, ref_dir, mut_dir, ref_gro)

        ligand_seed_pdb = layout.common_dir / "hybrid" / "ligand_seed.pdb"

        # Build bound leg
        if bound_system_dir:
            self._copy_prism_system_to_leg(Path(bound_system_dir), layout.bound_dir, "bound")
        else:
            self._write_bound_leg(layout, receptor_path, ligand_seed_pdb)

        # Build unbound leg
        if unbound_system_dir:
            self._copy_prism_system_to_leg(Path(unbound_system_dir), layout.unbound_dir, "unbound")
        else:
            self._write_unbound_leg(layout, ligand_seed_pdb)

        # Generate MDP files for both legs
        print(f"\n{'='*70}")
        print("Generating MDP files for FEP legs")
        print(f"{'='*70}")
        print(f"  Lambda strategy: {self.lambda_strategy}")
        print(f"  Lambda windows: {self.lambda_windows}")
        print(f"  Bound MDP dir: {layout.bound_dir / 'mdps'}")
        print(f"  Unbound MDP dir: {layout.unbound_dir / 'mdps'}")

        write_fep_mdps(
            output_dir=str(layout.bound_dir / "mdps"),
            lambda_strategy=self.lambda_strategy,
            lambda_distribution=self.lambda_distribution,
            lambda_windows=self.lambda_windows,
            config=self.config,
            leg_name="bound",
        )
        print(f"  ✓ Bound MDP files generated")

        write_fep_mdps(
            output_dir=str(layout.unbound_dir / "mdps"),
            lambda_strategy=self.lambda_strategy,
            lambda_distribution=self.lambda_distribution,
            lambda_windows=self.lambda_windows,
            config=self.config,
            leg_name="unbound",
        )
        print(f"  ✓ Unbound MDP files generated")
        print(f"{'='*70}\n")

        # Generate run scripts
        self._write_fep_run_script(layout.bound_dir, "bound")
        self._write_fep_run_script(layout.unbound_dir, "unbound")

        self._write_root_scripts(layout)
        self._write_manifest(
            layout,
            receptor_path,
            hybrid_itp_path.parent,
            extra_sources={
                "hybrid_itp": str(hybrid_itp_path.resolve()),
                "reference_ligand_dir": str(ref_dir.resolve()),
                "mutant_ligand_dir": str(mut_dir.resolve()) if mut_dir else None,
            },
        )

        return layout

    def _validate_inputs(self, receptor_path: Path, hybrid_dir: Path) -> None:
        if not receptor_path.exists():
            raise FileNotFoundError(f"Receptor PDB not found: {receptor_path}")
        if not hybrid_dir.exists():
            raise FileNotFoundError(f"Hybrid ligand directory not found: {hybrid_dir}")
        if not (hybrid_dir / "LIG.itp").exists():
            raise FileNotFoundError(f"Hybrid ligand ITP not found: {hybrid_dir / 'LIG.itp'}")
        if not (hybrid_dir / "LIG.gro").exists():
            raise FileNotFoundError(f"Hybrid ligand GRO not found: {hybrid_dir / 'LIG.gro'}")
        self.molecule_name = "LIG"

    def _prepare_layout(self) -> FEPScaffoldLayout:
        root = self.output_dir
        if root.exists() and self.overwrite:
            shutil.rmtree(root)

        layout = FEPScaffoldLayout(
            root=root,
            common_dir=root / "common",
            hybrid_dir=root / "common" / "hybrid",
            protein_dir=root / "common" / "protein",
            bound_dir=root / "bound",
            unbound_dir=root / "unbound",
        )

        for directory in [
            layout.common_dir,
            layout.hybrid_dir,
            layout.protein_dir,
            layout.bound_dir / "input",
            layout.bound_dir / "build",
            layout.bound_dir / "mdps",
            layout.unbound_dir / "input",
            layout.unbound_dir / "build",
            layout.unbound_dir / "mdps",
        ]:
            directory.mkdir(parents=True, exist_ok=True)

        return layout

    def _copy_common_assets(self, receptor_path: Path, hybrid_dir: Path, layout: FEPScaffoldLayout) -> None:
        shutil.copy2(receptor_path, layout.protein_dir / "receptor.pdb")

        for asset in self._iter_hybrid_assets(hybrid_dir):
            target = layout.hybrid_dir / asset.name
            shutil.copy2(asset, target)

        ligand_seed_pdb = layout.hybrid_dir / "ligand_seed.pdb"
        ligand_seed_pdb.write_text(self._gro_to_pdb(hybrid_dir / "LIG.gro"))

    def _iter_hybrid_assets(self, hybrid_dir: Path) -> Iterable[Path]:
        asset_names = [
            "LIG.itp",
            "LIG.gro",
            "LIG.top",
            "atomtypes_LIG.itp",
            "posre_LIG.itp",
        ]
        for name in asset_names:
            path = hybrid_dir / name
            if path.exists():
                yield path

    def _prepare_hybrid_package_from_components(
        self,
        layout: FEPScaffoldLayout,
        hybrid_itp: Path,
        reference_ligand_dir: Path,
        mutant_ligand_dir: Optional[Path],
        seed_gro: Path,
    ) -> None:
        self.hybrid_itp_filename = "hybrid.itp"
        self.hybrid_atomtypes_filename = "atomtypes_hybrid.itp"
        self.hybrid_forcefield_filename = "ff_hybrid.itp"
        self.hybrid_posre_filename = "posre_hybrid.itp"
        self.hybrid_gro_filename = "hybrid.gro"
        self.molecule_name = "HYB"

        normalized_itp = self._normalize_hybrid_itp(hybrid_itp.read_text(), self.molecule_name)
        hybrid_itp_output = layout.hybrid_dir / self.hybrid_itp_filename
        hybrid_itp_output.write_text(normalized_itp)

        if self._itp_needs_bonded_sections(normalized_itp):
            ITPBuilder.write_complete_hybrid_itp(
                output_path=str(hybrid_itp_output),
                hybrid_itp=str(hybrid_itp_output),
                ligand_a_itp=str(reference_ligand_dir / "LIG.itp"),
                ligand_b_itp=str(mutant_ligand_dir / "LIG.itp")
                if mutant_ligand_dir
                else str(reference_ligand_dir / "LIG.itp"),
                molecule_name=self.molecule_name,
            )
            normalized_itp = hybrid_itp_output.read_text()

        atomtypes = self._collect_atomtypes(reference_ligand_dir, mutant_ligand_dir)
        atomtypes_content = self._build_hybrid_atomtypes_itp(normalized_itp, atomtypes)
        (layout.hybrid_dir / self.hybrid_atomtypes_filename).write_text(atomtypes_content)

        defaults_block = self._extract_defaults_block(reference_ligand_dir, mutant_ligand_dir)
        if defaults_block:
            ff_hybrid = defaults_block + f'\n#include "{self.hybrid_atomtypes_filename}"\n'
            (layout.hybrid_dir / self.hybrid_forcefield_filename).write_text(ff_hybrid)

        hybrid_gro_content = self._build_hybrid_gro(
            hybrid_itp_content=normalized_itp,
            reference_gro=seed_gro,
            mutant_gro=(mutant_ligand_dir / "LIG.gro") if mutant_ligand_dir else None,
        )
        (layout.hybrid_dir / self.hybrid_gro_filename).write_text(hybrid_gro_content)
        (layout.hybrid_dir / "ligand_seed.pdb").write_text(
            self._gro_to_pdb(layout.hybrid_dir / self.hybrid_gro_filename)
        )

        ref_posre = reference_ligand_dir / "posre_LIG.itp"
        if ref_posre.exists():
            shutil.copy2(ref_posre, layout.hybrid_dir / self.hybrid_posre_filename)

        self._write_hybrid_top_template(layout.hybrid_dir / "hybrid.top")

    def _write_bound_leg(self, layout: FEPScaffoldLayout, receptor_path: Path, ligand_seed_pdb: Path) -> None:
        bound_input = layout.bound_dir / "input"
        shutil.copy2(layout.protein_dir / "receptor.pdb", bound_input / "receptor.pdb")
        shutil.copy2(ligand_seed_pdb, bound_input / "ligand_seed.pdb")
        self._write_complex_seed(
            bound_input / "receptor.pdb",
            bound_input / "ligand_seed.pdb",
            bound_input / "complex_seed.pdb",
        )

        # Note: This creates a placeholder. Use build_from_components() with PRISM-built systems for complete setup.
        self._create_placeholder_system("bound", layout, receptor_path)

        write_fep_mdps(
            output_dir=str(layout.bound_dir / "mdps"),
            lambda_strategy=self.lambda_strategy,
            lambda_distribution=self.lambda_distribution,
            lambda_windows=self.lambda_windows,
            config=self.config,
            leg_name="bound",
        )
        self._write_fep_run_script(layout.bound_dir, "bound")

    def _write_unbound_leg(self, layout: FEPScaffoldLayout, ligand_seed_pdb: Path) -> None:
        unbound_input = layout.unbound_dir / "input"
        shutil.copy2(ligand_seed_pdb, unbound_input / "ligand_seed.pdb")

        # Note: This creates a placeholder. Use build_from_components() with PRISM-built systems for complete setup.
        self._create_placeholder_system("unbound", layout, protein_path=None)

        write_fep_mdps(
            output_dir=str(layout.unbound_dir / "mdps"),
            lambda_strategy=self.lambda_strategy,
            lambda_distribution=self.lambda_distribution,
            lambda_windows=self.lambda_windows,
            config=self.config,
            leg_name="unbound",
        )
        self._write_fep_run_script(layout.unbound_dir, "unbound")

    def _write_manifest(
        self,
        layout: FEPScaffoldLayout,
        receptor_path: Path,
        hybrid_dir: Path,
        extra_sources: Optional[Dict[str, Optional[str]]] = None,
    ) -> None:
        manifest = {
            "scaffold_version": 1,
            "receptor_pdb": str(receptor_path.resolve()),
            "hybrid_ligand_dir": str(hybrid_dir.resolve()),
            "lambda_windows": self.lambda_windows,
            "lambda_strategy": self.lambda_strategy,
            "lambda_distribution": self.lambda_distribution,
            "hybrid_assets": {
                "molecule_name": self.molecule_name,
                "itp": self.hybrid_itp_filename,
                "atomtypes": self.hybrid_atomtypes_filename,
                "forcefield": self.hybrid_forcefield_filename,
                "gro": self.hybrid_gro_filename,
                "posre": self.hybrid_posre_filename,
            },
            "same_box_size": {
                "enabled": True,
                "policy": "build bound leg first and reuse the final box vectors for the unbound leg during full integration",
            },
            "reuse_plan": {
                "bound_leg": [
                    "Reuse PRISM protein preparation and topology generation",
                    "Reuse PRISM SystemBuilder for protein+hybrid-ligand box/solvation/ions",
                ],
                "unbound_leg": [
                    "Reuse the hybrid ligand assets staged in common/hybrid",
                    "Reuse PRISM boxing/solvation helpers with bound-leg box vectors",
                ],
            },
        }
        if extra_sources:
            manifest["sources"] = extra_sources
        (layout.root / "fep_scaffold.json").write_text(json.dumps(manifest, indent=2) + "\n")

    def _write_root_scripts(self, layout: FEPScaffoldLayout) -> None:
        script = """#!/usr/bin/env bash
set -euo pipefail

echo "PRISM-FEP scaffold created."
echo "Next intended order:"
echo "  1. Complete bound-leg system build"
echo "  2. Reuse bound-leg box vectors for unbound leg"
echo "  3. Run bound/localrun.sh"
echo "  4. Run unbound/localrun.sh"
"""
        path = layout.root / "run_all.sh"
        path.write_text(script)
        path.chmod(0o755)

    def _write_unbound_topology(self, topol_path: Path) -> None:
        content = """; PRISM-FEP unbound topology template
; This file is already close to a valid ligand-only topology if the force field
; include chain is complete inside common/hybrid/.
""".format(posre=self.hybrid_posre_filename)
        ff_path = topol_path.parent.parent / "common" / "hybrid" / self.hybrid_forcefield_filename
        atomtypes_path = topol_path.parent.parent / "common" / "hybrid" / self.hybrid_atomtypes_filename
        if ff_path.exists():
            content += f'#include "../common/hybrid/{self.hybrid_forcefield_filename}"\n'
        elif atomtypes_path.exists():
            content += f'#include "../common/hybrid/{self.hybrid_atomtypes_filename}"\n'
        content += f'#include "../common/hybrid/{self.hybrid_itp_filename}"\n\n'
        content += "[ system ]\nUnbound FEP leg\n\n"
        content += f"[ molecules ]\n{self.molecule_name} 1\n"
        topol_path.write_text(content)

    def _write_bound_topology_template(self, topol_path: Path) -> None:
        content = """; PRISM-FEP bound topology template
; Protein topology include will be supplied during full builder integration.
; The ligand include chain is already staged.

; TODO: include generated protein topology here, e.g.
; #include "topol_protein.itp"
"""
        atomtypes_path = topol_path.parent.parent / "common" / "hybrid" / self.hybrid_atomtypes_filename
        if atomtypes_path.exists():
            content += f'#include "../common/hybrid/{self.hybrid_atomtypes_filename}"\n'
        content += f'#include "../common/hybrid/{self.hybrid_itp_filename}"\n\n'
        content += "[ system ]\nBound FEP leg\n\n"
        content += f"[ molecules ]\n; Protein      1\n{self.molecule_name} 1\n"
        topol_path.write_text(content)

    def _write_hybrid_top_template(self, output_path: Path) -> None:
        content = "; PRISM-FEP hybrid ligand topology bundle\n"
        ff_path = output_path.parent / self.hybrid_forcefield_filename
        atomtypes_path = output_path.parent / self.hybrid_atomtypes_filename
        if ff_path.exists():
            content += f'#include "{self.hybrid_forcefield_filename}"\n'
        elif atomtypes_path.exists():
            content += f'#include "{self.hybrid_atomtypes_filename}"\n'
        content += f'#include "{self.hybrid_itp_filename}"\n\n'

        content += f"""
[ system ]
Hybrid ligand package

[ molecules ]
{self.molecule_name} 1
"""
        output_path.write_text(content)

    def _gro_to_pdb(self, gro_file: Path) -> str:
        lines = gro_file.read_text().splitlines()
        pdb_lines = ["REMARK   Generated from LIG.gro for PRISM-FEP scaffold"]
        atom_lines = lines[2:-1]
        for serial, line in enumerate(atom_lines, start=1):
            resnum = int(line[:5].strip())
            resname = line[5:10].strip() or "LIG"
            atom_name = line[10:15].strip() or f"A{serial}"
            atom_serial = int(line[15:20].strip())
            x = float(line[20:28].strip()) * 10.0
            y = float(line[28:36].strip()) * 10.0
            z = float(line[36:44].strip()) * 10.0
            element = "".join([c for c in atom_name if c.isalpha()])[:1].upper() or "C"
            pdb_lines.append(
                f"HETATM{serial:5d} {atom_name:>4s} {resname:>3s} A{resnum:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
            )
        pdb_lines.append("END")
        return "\n".join(pdb_lines) + "\n"

    def _write_complex_seed(self, receptor_pdb: Path, ligand_pdb: Path, output_pdb: Path) -> None:
        receptor_records = self._extract_pdb_records(receptor_pdb)
        ligand_records = self._extract_pdb_records(ligand_pdb)

        merged = ["REMARK   PRISM-FEP bound-leg seed complex"]
        serial = 1
        for record in receptor_records + ligand_records:
            merged.append(record[:6] + f"{serial:5d}" + record[11:])
            serial += 1
        merged.append("END")
        output_pdb.write_text("\n".join(merged) + "\n")

    def _extract_pdb_records(self, pdb_path: Path) -> list[str]:
        records = []
        for line in pdb_path.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                records.append(line)
        return records

    def _normalize_hybrid_itp(self, content: str, molecule_name: str) -> str:
        lines = content.splitlines()
        output = []
        in_moleculetype = False
        in_atoms = False
        replaced = False

        for line in lines:
            stripped = line.strip()
            if '#include "posre_LIG.itp"' in line:
                output.append(line.replace("posre_LIG.itp", self.hybrid_posre_filename))
                continue
            if stripped.lower() == "[ atoms ]":
                in_atoms = True
                output.append(line)
                continue
            if stripped.startswith("[") and in_atoms:
                in_atoms = False
            if in_atoms and stripped and not stripped.startswith(";"):
                output.append(self._normalize_hybrid_atom_line(stripped))
                continue
            output.append(line)
            if stripped.lower() == "[ moleculetype ]":
                in_moleculetype = True
                continue
            if in_moleculetype and stripped and not stripped.startswith(";") and not replaced:
                output[-1] = f"{molecule_name:16s}     3"
                replaced = True
                in_moleculetype = False

        return "\n".join(output) + "\n"

    def _normalize_hybrid_atom_line(self, line: str) -> str:
        parts = line.split()
        if len(parts) >= 10 and not self._looks_like_float(parts[8]) and len(parts) == 10:
            parts.append(parts[7])
        return " ".join(parts)

    def _itp_needs_bonded_sections(self, content: str) -> bool:
        lower = content.lower()
        return "[ bonds ]" not in lower and "[ angles ]" not in lower and "[ dihedrals ]" not in lower

    def _collect_atomtypes(self, reference_ligand_dir: Path, mutant_ligand_dir: Optional[Path]) -> Dict[str, str]:
        atomtypes = {}
        for ligand_dir in [reference_ligand_dir, mutant_ligand_dir]:
            if ligand_dir is None:
                continue
            atomtypes_file = ligand_dir / "atomtypes_LIG.itp"
            if not atomtypes_file.exists():
                continue
            for atomtype_name, line in self._parse_atomtypes_file(atomtypes_file).items():
                atomtypes.setdefault(atomtype_name, line)
        return atomtypes

    def _parse_atomtypes_file(self, atomtypes_file: Path) -> Dict[str, str]:
        atomtypes = {}
        in_section = False
        for line in atomtypes_file.read_text().splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ atomtypes ]":
                in_section = True
                continue
            if stripped.startswith("[") and in_section:
                break
            if not in_section or not stripped or stripped.startswith(";"):
                continue
            parts = stripped.split()
            atomtypes[parts[0]] = stripped
        return atomtypes

    def _build_hybrid_atomtypes_itp(self, hybrid_itp_content: str, source_atomtypes: Dict[str, str]) -> str:
        used_types = self._parse_used_atomtypes(hybrid_itp_content)
        lines = [
            "[ atomtypes ]",
            "; merged atomtypes for hybrid ligand",
        ]

        written = set()
        dummy_requests = set()
        for atomtype in used_types:
            if atomtype.startswith("DUM_"):
                dummy_requests.add(atomtype)
                continue
            source = source_atomtypes.get(atomtype)
            if source and atomtype not in written:
                lines.append(source)
                written.add(atomtype)

        for dummy_atomtype in sorted(dummy_requests):
            if dummy_atomtype in written:
                continue
            base_type = dummy_atomtype[4:]
            base_line = source_atomtypes.get(base_type)
            bond_type = base_type
            if base_line:
                bond_type = base_line.split()[1]
            lines.append(
                f"{dummy_atomtype:<8s} {bond_type:<8s} 0.00000 0.00000 A 0.00000e+00 0.00000e+00 ; dummy from {base_type}"
            )
            written.add(dummy_atomtype)

        missing_non_dummy = sorted(
            atomtype for atomtype in used_types if not atomtype.startswith("DUM_") and atomtype not in source_atomtypes
        )
        if missing_non_dummy:
            raise ValueError(f"Missing atomtype definitions for: {', '.join(missing_non_dummy)}")

        return "\n".join(lines) + "\n"

    def _extract_defaults_block(self, reference_ligand_dir: Path, mutant_ligand_dir: Optional[Path]) -> str:
        defaults_blocks = []
        for ligand_dir in [reference_ligand_dir, mutant_ligand_dir]:
            if ligand_dir is None:
                continue
            top_file = ligand_dir / "LIG.top"
            if not top_file.exists():
                continue
            block = self._parse_defaults_block(top_file.read_text())
            if block:
                defaults_blocks.append(block)

        if not defaults_blocks:
            return self._infer_defaults_block(reference_ligand_dir, mutant_ligand_dir)

        first = defaults_blocks[0]
        for block in defaults_blocks[1:]:
            if block != first:
                raise ValueError("Inconsistent [ defaults ] blocks between ligand topologies")
        return first

    def _infer_defaults_block(self, reference_ligand_dir: Path, mutant_ligand_dir: Optional[Path]) -> str:
        atomtypes = self._collect_atomtypes(reference_ligand_dir, mutant_ligand_dir)
        if atomtypes and all(name.startswith("opls_") or name.startswith("DUM_opls_") for name in atomtypes):
            return (
                "[ defaults ]\n"
                "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n"
                "1               3               yes             0.5     0.5\n"
            )
        return ""

    def _parse_defaults_block(self, top_content: str) -> str:
        lines = top_content.splitlines()
        block = []
        in_defaults = False
        for line in lines:
            stripped = line.strip()
            if stripped.lower() == "[ defaults ]":
                in_defaults = True
                block.append("[ defaults ]")
                continue
            if stripped.startswith("[") and in_defaults:
                break
            if in_defaults:
                block.append(line)
        if len(block) <= 1:
            return ""
        return "\n".join(block).rstrip() + "\n"

    def _build_hybrid_gro(
        self,
        hybrid_itp_content: str,
        reference_gro: Path,
        mutant_gro: Optional[Path],
    ) -> str:
        ref_atoms, ref_box = self._parse_gro_atoms(reference_gro)
        mut_atoms, _ = self._parse_gro_atoms(mutant_gro) if mutant_gro else ({}, None)
        hybrid_atoms = self._parse_hybrid_atoms(hybrid_itp_content)

        output = ["Hybrid ligand", f"{len(hybrid_atoms):5d}"]
        for index, atom in enumerate(hybrid_atoms, start=1):
            source_atom = self._resolve_hybrid_atom_coordinates(atom["atom_name"], ref_atoms, mut_atoms)
            output.append(
                f"{atom['resnr']:5d}{atom['residue']:<5s}{atom['atom_name']:>5s}{index:5d}"
                f"{source_atom['x']:8.3f}{source_atom['y']:8.3f}{source_atom['z']:8.3f}"
            )
        output.append(self._normalize_box_line(ref_box))
        return "\n".join(output) + "\n"

    def _parse_gro_atoms(self, gro_path: Optional[Path]) -> tuple[Dict[str, list[Dict[str, float]]], Optional[str]]:
        if gro_path is None or not gro_path.exists():
            return {}, None

        lines = gro_path.read_text().splitlines()
        atom_count = int(lines[1].strip())
        atoms = {}
        for line in lines[2 : 2 + atom_count]:
            atom_name = line[10:15].strip()
            atoms.setdefault(atom_name, []).append(
                {
                    "resnr": int(line[:5].strip()),
                    "resname": line[5:10].strip(),
                    "x": float(line[20:28].strip()),
                    "y": float(line[28:36].strip()),
                    "z": float(line[36:44].strip()),
                }
            )
        box = lines[2 + atom_count] if len(lines) > 2 + atom_count else None
        return atoms, box

    def _parse_hybrid_atoms(self, hybrid_itp_content: str) -> list[Dict[str, str]]:
        atoms = []
        in_atoms = False
        for line in hybrid_itp_content.splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ atoms ]":
                in_atoms = True
                continue
            if stripped.startswith("[") and in_atoms:
                break
            if not in_atoms or not stripped or stripped.startswith(";"):
                continue
            parts = stripped.split()
            if len(parts) < 5:
                continue
            atoms.append(
                {
                    "resnr": int(parts[2]),
                    "residue": parts[3],
                    "atom_name": parts[4],
                }
            )
        return atoms

    def _resolve_hybrid_atom_coordinates(
        self,
        atom_name: str,
        ref_atoms: Dict[str, list[Dict[str, float]]],
        mut_atoms: Dict[str, list[Dict[str, float]]],
    ) -> Dict[str, float]:
        candidates = []
        if atom_name.endswith("A"):
            base = atom_name[:-1]
            candidates = [(ref_atoms, atom_name), (ref_atoms, base), (mut_atoms, atom_name), (mut_atoms, base)]
        elif atom_name.endswith("B"):
            base = atom_name[:-1]
            candidates = [(mut_atoms, atom_name), (mut_atoms, base), (ref_atoms, atom_name), (ref_atoms, base)]
        else:
            candidates = [(ref_atoms, atom_name), (mut_atoms, atom_name)]

        for atom_pool, key in candidates:
            entries = atom_pool.get(key)
            if entries:
                return entries[0]
        raise ValueError(f"Could not resolve coordinates for hybrid atom: {atom_name}")

    def _normalize_box_line(self, box_line: Optional[str]) -> str:
        if box_line:
            try:
                values = [float(value) for value in box_line.split()[:3]]
                if len(values) == 3 and min(values) >= 2.2:
                    return box_line
            except ValueError:
                pass
        return "   3.00000   3.00000   3.00000"

    def _parse_used_atomtypes(self, hybrid_itp_content: str) -> set[str]:
        used_types = set()
        in_atoms = False
        for line in hybrid_itp_content.splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ atoms ]":
                in_atoms = True
                continue
            if stripped.startswith("[") and in_atoms:
                break
            if not in_atoms or not stripped or stripped.startswith(";"):
                continue
            parts = stripped.split()
            if len(parts) < 8:
                continue
            used_types.add(parts[1])
            if len(parts) >= 9 and not self._looks_like_float(parts[8]):
                used_types.add(parts[8])
        return used_types

    def _looks_like_float(self, value: str) -> bool:
        try:
            float(value)
        except ValueError:
            return False
        return True

    def build_from_prism_system(
        self,
        bound_system_dir: str,
        unbound_system_dir: Optional[str],
        hybrid_itp: str,
        reference_ligand_dir: str,
        mutant_ligand_dir: str,
    ) -> "FEPScaffoldLayout":
        """
        Build FEP scaffold from pre-built PRISM MD systems.

        This is the recommended workflow:
        1. Build bound system: prism protein.pdb ref_ligand.mol2 -o bound_md
        2. Build unbound system: prism dummy.pdb ref_ligand.mol2 -o unbound_md
        3. Generate hybrid topology: (mapping tool)
        4. Call this method to create FEP scaffold with complete systems

        Parameters:
        -----------
        bound_system_dir : str
            Path to PRISM-built bound system (GMX_PROLIG_MD directory)
        unbound_system_dir : Optional[str]
            Path to PRISM-built unbound system (or None to create placeholder)
        hybrid_itp : str
            Path to hybrid ligand ITP file
        reference_ligand_dir : str
            Path to reference ligand force field directory
        mutant_ligand_dir : str
            Path to mutant ligand force field directory

        Returns:
        --------
        FEPScaffoldLayout
            Layout object with paths to bound/unbound/common directories
        """
        layout = self._create_scaffold_structure()

        # Copy hybrid ligand files
        self._copy_hybrid_ligand_files(
            Path(hybrid_itp).parent, Path(reference_ligand_dir), Path(mutant_ligand_dir), layout
        )

        # Copy bound system
        self._copy_prism_system_to_leg(Path(bound_system_dir), layout.bound_dir, "bound")

        # Copy unbound system
        if unbound_system_dir:
            self._copy_prism_system_to_leg(Path(unbound_system_dir), layout.unbound_dir, "unbound")
        else:
            raise ValueError("unbound_system_dir is required - build it first with PRISMBuilder")

        # Generate MDP files for both legs
        print(f"\n{'='*70}")
        print("Generating MDP files for FEP legs")
        print(f"{'='*70}")
        print(f"  Lambda strategy: {self.lambda_strategy}")
        print(f"  Lambda windows: {self.lambda_windows}")
        print(f"  Bound MDP dir: {layout.bound_dir / 'mdps'}")
        print(f"  Unbound MDP dir: {layout.unbound_dir / 'mdps'}")

        write_fep_mdps(
            output_dir=str(layout.bound_dir / "mdps"),
            lambda_strategy=self.lambda_strategy,
            lambda_distribution=self.lambda_distribution,
            lambda_windows=self.lambda_windows,
            config=self.config,
            leg_name="bound",
        )
        print(f"  ✓ Bound MDP files generated")

        write_fep_mdps(
            output_dir=str(layout.unbound_dir / "mdps"),
            lambda_strategy=self.lambda_strategy,
            lambda_distribution=self.lambda_distribution,
            lambda_windows=self.lambda_windows,
            config=self.config,
            leg_name="unbound",
        )
        print(f"  ✓ Unbound MDP files generated")
        print(f"{'='*70}\n")

        # Generate run scripts
        self._write_fep_run_script(layout.bound_dir, "bound")
        self._write_fep_run_script(layout.unbound_dir, "unbound")

        return layout

    def _copy_prism_system_to_leg(self, prism_system_dir: Path, leg_dir: Path, leg_name: str) -> None:
        """
        Copy files from PRISM-built system to FEP leg directory.

        Parameters:
        -----------
        prism_system_dir : Path
            Path to GMX_PROLIG_MD directory from PRISM
        leg_dir : Path
            Target FEP leg directory (bound or unbound)
        leg_name : str
            "bound" or "unbound"
        """
        print(f"\n{'='*70}")
        print(f"Copying {leg_name} leg from PRISM system")
        print(f"{'='*70}")

        # Check source files exist
        solv_ions_gro = prism_system_dir / "solv_ions.gro"
        topol_top = prism_system_dir / "topol.top"

        if not solv_ions_gro.exists():
            raise FileNotFoundError(f"solv_ions.gro not found in {prism_system_dir}")
        if not topol_top.exists():
            raise FileNotFoundError(f"topol.top not found in {prism_system_dir}")

        # Copy coordinate file
        shutil.copy2(solv_ions_gro, leg_dir / "input" / "conf.gro")
        print(f"  ✓ Copied solv_ions.gro → {leg_name}/input/conf.gro")

        # Copy and modify topology to use hybrid ligand
        self._copy_and_modify_topology(topol_top, leg_dir / "topol.top")
        print(f"  ✓ Copied and modified topol.top → {leg_name}/topol.top")

        print(f"  ✓ {leg_name.capitalize()} leg system ready")
        print(f"{'='*70}\n")

    def _copy_and_modify_topology(self, source_top: Path, target_top: Path) -> None:
        """
        Copy topology file and modify to use hybrid ligand.

        Replaces ligand ITP include with hybrid ligand ITP.
        """
        content = source_top.read_text()

        # Replace ligand ITP include with hybrid ligand
        # Pattern: #include "path/to/LIG.itp" or similar
        import re

        content = re.sub(r'#include\s+"[^"]*LIG\.itp"', '#include "../common/hybrid/hybrid.itp"', content)

        target_top.write_text(content)

    def _prepare_protein_for_system_builder(self, protein_path: Path, leg_dir: Path) -> Path:
        """Prepare protein PDB for system builder"""
        cleaned_path = leg_dir / "build" / "cleaned_protein.pdb"
        cleaned_path.parent.mkdir(parents=True, exist_ok=True)

        # Use absolute path to avoid path resolution issues
        abs_protein_path = protein_path.resolve()
        shutil.copy2(abs_protein_path, cleaned_path)

        return cleaned_path.resolve()  # Return absolute path

    def _create_dummy_protein(self, leg_dir: Path) -> Path:
        """Create minimal dummy protein for ligand-only systems"""
        dummy_pdb = """REMARK   Dummy protein for ligand-only FEP system
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
END
"""
        dummy_path = leg_dir / "build" / "dummy_protein.pdb"
        dummy_path.parent.mkdir(parents=True, exist_ok=True)
        dummy_path.write_text(dummy_pdb)
        return dummy_path.resolve()  # Return absolute path

    def _get_forcefield_indices(self, config: dict) -> tuple:
        """Get force field and water model indices from config"""
        ff_name = config.get("forcefield", {}).get("protein", "amber14sb")
        water_name = config.get("forcefield", {}).get("water", "tip3p")
        ff_map = {"amber14sb": 0, "amber99sb-ildn": 1, "charmm36": 2}
        water_map = {"tip3p": 0, "tip4p": 1, "spc": 2}
        return ff_map.get(ff_name.lower(), 0), water_map.get(water_name.lower(), 0)

    def _get_forcefield_info(self, _config: dict, ff_idx: int) -> dict:
        """Get force field info dict"""
        ff_names = ["amber14sb", "amber99sb-ildn", "charmm36"]
        # SystemBuilder expects path to be set, use GROMACS built-in path
        return {"name": ff_names[ff_idx], "dir": f"{ff_names[ff_idx]}.ff", "path": f"{ff_names[ff_idx]}.ff"}

    def _get_water_info(self, _conf, water_idx: int) -> dict:
        """Get water model info dict"""
        water_names = ["tip3p", "tip4p", "spc"]
        return {"name": water_names[water_idx]}

    def _write_complete_topology(
        self, topol_path: Path, leg_name: str, layout: FEPScaffoldLayout, model_path: Path
    ) -> None:
        """Write complete topology including protein/water/ions"""
        content = f"; PRISM-FEP {leg_name} topology\n\n"

        # Include force field
        ff_path = layout.hybrid_dir / self.hybrid_forcefield_filename
        if ff_path.exists():
            content += f'#include "../common/hybrid/{self.hybrid_forcefield_filename}"\n'

        # Include protein topology (if bound leg)
        if leg_name == "bound":
            protein_itp = model_path / "topol_Protein_chain_A.itp"
            if protein_itp.exists():
                shutil.copy2(protein_itp, topol_path.parent / "topol_Protein_chain_A.itp")
                content += '#include "topol_Protein_chain_A.itp"\n'

        # Include hybrid ligand
        content += f'#include "../common/hybrid/{self.hybrid_itp_filename}"\n\n'

        # Include water and ions
        content += "; Include water topology\n"
        content += '#include "amber14sb.ff/tip3p.itp"\n\n'
        content += "#ifdef POSRES_WATER\n"
        content += "[ position_restraints ]\n"
        content += "   1    1       1000       1000       1000\n"
        content += "#endif\n\n"
        content += '#include "amber14sb.ff/ions.itp"\n\n'

        content += f"[ system ]\n{leg_name.capitalize()} FEP leg\n\n"
        content += "[ molecules ]\n"

        # Read molecule counts
        model_top = model_path / "topol.top"
        if model_top.exists():
            molecules = self._extract_molecules_from_topology(model_top)
            for mol_name, count in molecules:
                if mol_name == "LIG":
                    mol_name = self.molecule_name
                content += f"{mol_name:20s} {count}\n"

        topol_path.write_text(content)

    def _extract_molecules_from_topology(self, topol_path: Path) -> list:
        """Extract molecule list from topology"""
        molecules = []
        in_molecules = False
        for line in topol_path.read_text().splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ molecules ]":
                in_molecules = True
                continue
            if in_molecules and stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if len(parts) >= 2:
                    try:
                        molecules.append((parts[0], int(parts[1])))
                    except ValueError:
                        continue
        return molecules

    def _write_fep_run_script(self, leg_dir: Path, leg_name: str) -> None:
        """Generate real FEP run script"""
        script = f"""#!/usr/bin/env bash
set -euo pipefail

LEG_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
MDP_DIR="${{LEG_DIR}}/mdps"
INPUT_DIR="${{LEG_DIR}}/input"
BUILD_DIR="${{LEG_DIR}}/build"

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    PRISM-FEP {leg_name.upper()} LEG                                    ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

cd "${{LEG_DIR}}"

# Energy minimization
mkdir -p ${{BUILD_DIR}}

# Energy minimization
if [ -f ${{BUILD_DIR}}/em.gro ]; then
    echo "EM already completed, skipping..."
elif [ -f ${{BUILD_DIR}}/em.tpr ]; then
    echo "EM checkpoint found, resuming..."
    gmx mdrun -deffnm ${{BUILD_DIR}}/em -cpi ${{BUILD_DIR}}/em.cpt -v
else
    echo "Running energy minimization..."
    gmx grompp -f ${{MDP_DIR}}/em.mdp -c ${{INPUT_DIR}}/conf.gro -p topol.top -o ${{BUILD_DIR}}/em.tpr -maxwarn 2
    gmx mdrun -deffnm ${{BUILD_DIR}}/em -ntmpi 1 -ntomp 10 -v
fi

# NVT equilibration
if [ -f ${{BUILD_DIR}}/nvt.gro ]; then
    echo "NVT already completed, skipping..."
elif [ -f ${{BUILD_DIR}}/nvt.tpr ]; then
    echo "NVT checkpoint found, resuming..."
    gmx mdrun -deffnm ${{BUILD_DIR}}/nvt -cpi ${{BUILD_DIR}}/nvt.cpt -v
else
    echo "Running NVT equilibration..."
    gmx grompp -f ${{MDP_DIR}}/nvt.mdp -c ${{BUILD_DIR}}/em.gro -p topol.top -o ${{BUILD_DIR}}/nvt.tpr -maxwarn 2
    gmx mdrun -deffnm ${{BUILD_DIR}}/nvt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -v
fi

# NPT equilibration
if [ -f ${{BUILD_DIR}}/npt.gro ]; then
    echo "NPT already completed, skipping..."
elif [ -f ${{BUILD_DIR}}/npt.tpr ]; then
    echo "NPT checkpoint found, resuming..."
    gmx mdrun -deffnm ${{BUILD_DIR}}/npt -cpi ${{BUILD_DIR}}/npt.cpt -v
else
    echo "Running NPT equilibration..."
    gmx grompp -f ${{MDP_DIR}}/npt.mdp -c ${{BUILD_DIR}}/nvt.gro -t ${{BUILD_DIR}}/nvt.cpt -p topol.top -o ${{BUILD_DIR}}/npt.tpr -maxwarn 2
    gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -v
fi

# FEP production runs
echo "Running FEP production for all lambda windows..."
for lambda_mdp in ${{MDP_DIR}}/prod_*.mdp; do
    lambda_name=$(basename ${{lambda_mdp}} .mdp)
    if [ -f ${{BUILD_DIR}}/${{lambda_name}}.gro ]; then
        echo "  ${{lambda_name}}: already completed"
    elif [ -f ${{BUILD_DIR}}/${{lambda_name}}.tpr ]; then
        echo "  ${{lambda_name}}: resuming from checkpoint"
        gmx mdrun -deffnm ${{BUILD_DIR}}/${{lambda_name}} -cpi ${{BUILD_DIR}}/${{lambda_name}}.cpt -v
    else
        echo "  ${{lambda_name}}: starting"
        gmx grompp -f ${{lambda_mdp}} -c ${{BUILD_DIR}}/npt.gro -p topol.top -o ${{BUILD_DIR}}/${{lambda_name}}.tpr -maxwarn 2
        gmx mdrun -deffnm ${{BUILD_DIR}}/${{lambda_name}} -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -v
    fi
done

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    {leg_name.upper()} LEG COMPLETE                                   ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
"""
        path = leg_dir / "localrun.sh"
        path.write_text(script)
        path.chmod(0o755)
