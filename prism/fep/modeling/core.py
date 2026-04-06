#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Scaffold builders for PRISM-FEP modeling workflows.

This module orchestrates FEP scaffold building by delegating to specialized builders:
- HybridPackageBuilder: hybrid ligand force field assembly
- LegWriter: bound/unbound leg directory writing and topology management
- script_writer: shell script and SLURM submission script generation
"""

from __future__ import annotations

import json
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional

from prism.fep.gromacs.mdp_templates import write_fep_mdps
from prism.fep.modeling.hybrid_package import HybridPackageBuilder
from prism.fep.modeling.leg_writer import LegWriter
from prism.fep.modeling import script_writer


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
        self.output_dir = self._normalize_output_dir(Path(output_dir))
        self.config = config or {}
        self.lambda_windows = lambda_windows
        self.lambda_strategy = lambda_strategy
        self.lambda_distribution = lambda_distribution
        self.overwrite = overwrite
        self.replicas = max(
            1,
            int(
                self.config.get(
                    "replicas",
                    self.config.get("fep", {}).get("replicas", 1),
                )
            ),
        )

        # Initialize sub-builders
        self.hybrid_builder = HybridPackageBuilder()
        self.leg_writer = LegWriter(self.output_dir)
        self._sync_hybrid_assets()

    @staticmethod
    def _normalize_output_dir(output_dir: Path | str) -> Path:
        """Collapse accidental ``.../GMX_PROLIG_FEP/GMX_PROLIG_FEP`` nesting."""
        resolved = Path(output_dir).expanduser()
        while resolved.name == "GMX_PROLIG_FEP" and resolved.parent.name == "GMX_PROLIG_FEP":
            resolved = resolved.parent
        return resolved

    def _sync_hybrid_assets(self) -> None:
        """Mirror hybrid-package naming onto the scaffold builder and leg writer."""
        self.hybrid_itp_filename = self.hybrid_builder.hybrid_itp_filename
        self.hybrid_atomtypes_filename = self.hybrid_builder.hybrid_atomtypes_filename
        self.hybrid_forcefield_filename = self.hybrid_builder.hybrid_forcefield_filename
        self.hybrid_posre_filename = self.hybrid_builder.hybrid_posre_filename
        self.hybrid_gro_filename = self.hybrid_builder.hybrid_gro_filename
        self.molecule_name = self.hybrid_builder.molecule_name
        self.leg_writer.molecule_name = self.molecule_name

    def _discover_ligand_dir(self, base_dir: Path) -> Path:
        """
        Automatically discover ligand force field directory.

        Search patterns:
        1. base_dir/*/LIG.itp
        2. base_dir/*/*.itp
        3. base_dir/LIG.itp

        Parameters
        ----------
        base_dir : Path
            Base directory to search

        Returns
        -------
        Path
            Path to ligand directory

        Raises
        ------
        FileNotFoundError
            If no ligand ITP file found
        """
        # Pattern 1: */LIG.itp
        itp_files = list(base_dir.glob("*/LIG.itp"))
        if itp_files:
            return itp_files[0].parent

        # Pattern 2: */*.itp (any ITP file)
        itp_files = list(base_dir.glob("*/*.itp"))
        if itp_files:
            return itp_files[0].parent

        # Pattern 3: LIG.itp in base_dir
        if (base_dir / "LIG.itp").exists():
            return base_dir

        raise FileNotFoundError(
            f"Cannot find ligand ITP file in {base_dir}. " f"Expected patterns: */LIG.itp, */*.itp, or LIG.itp"
        )

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
        self._replicate_repeat_layout(layout)
        self._validate_scaffold_consistency(layout)
        script_writer.write_root_scripts(layout, self.config)
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
        auto_discover: bool = False,
    ) -> FEPScaffoldLayout:
        """
        Create the scaffold from an existing hybrid ITP plus ligand FF directories.

        This is the main entry point for building FEP scaffolds from PRISM-generated
        ligand force fields and hybrid topology.

        Parameters
        ----------
        receptor_pdb : str
            Path to receptor PDB file
        hybrid_itp : str
            Path to hybrid ITP file (atom mapping)
        reference_ligand_dir : str
            Directory containing reference ligand force field
        mutant_ligand_dir : Optional[str]
            Directory containing mutant ligand force field (None for single-topology)
        hybrid_gro : Optional[str]
            Path to hybrid GRO file (if None, uses reference_ligand_dir/LIG.gro)
        bound_system_dir : Optional[str]
            Path to PRISM-built bound system (GMX_PROLIG_MD directory)
        unbound_system_dir : Optional[str]
            Path to PRISM-built unbound system (GMX_PROLIG_MD directory)
        auto_discover : bool, optional
            Automatically discover ligand directories (default: False)

        Returns
        -------
        FEPScaffoldLayout
            The created scaffold layout
        """
        receptor_path = Path(receptor_pdb)
        hybrid_itp_path = Path(hybrid_itp)

        # Auto-discover ligand directories if requested
        if auto_discover:
            reference_ligand_path = self._discover_ligand_dir(Path(reference_ligand_dir))
            if mutant_ligand_dir:
                mutant_ligand_path = self._discover_ligand_dir(Path(mutant_ligand_dir))
            else:
                mutant_ligand_path = None
        else:
            reference_ligand_path = Path(reference_ligand_dir)
            mutant_ligand_path = Path(mutant_ligand_dir) if mutant_ligand_dir else None

        # Validate inputs
        if not receptor_path.exists():
            raise FileNotFoundError(f"Receptor PDB not found: {receptor_path}")
        if not hybrid_itp_path.exists():
            raise FileNotFoundError(f"Hybrid ITP not found: {hybrid_itp_path}")
        if not reference_ligand_path.exists():
            raise FileNotFoundError(f"Reference ligand directory not found: {reference_ligand_path}")
        if mutant_ligand_path and not mutant_ligand_path.exists():
            raise FileNotFoundError(f"Mutant ligand directory not found: {mutant_ligand_path}")

        # Prepare layout
        layout = self._prepare_layout()

        # Copy receptor
        shutil.copy2(receptor_path, layout.protein_dir / "receptor.pdb")

        # Determine seed GRO
        if hybrid_gro:
            seed_gro = Path(hybrid_gro)
        else:
            # Auto-discover LIG.gro (handle LIG.amb2gmx/ subdirectory)
            seed_gro = reference_ligand_path / "LIG.gro"
            if not seed_gro.exists():
                # Try LIG.amb2gmx/LIG.gro (PRISM output structure)
                seed_gro = reference_ligand_path / "LIG.amb2gmx" / "LIG.gro"
            if not seed_gro.exists():
                # Try any */LIG.gro pattern
                gro_files = list(reference_ligand_path.glob("*/LIG.gro"))
                if gro_files:
                    seed_gro = gro_files[0]

        if not seed_gro.exists():
            raise FileNotFoundError(f"Seed GRO file not found: {seed_gro}")

        # Build hybrid package
        self.hybrid_builder.prepare_hybrid_package_from_components(
            hybrid_dir=layout.hybrid_dir,
            hybrid_itp=hybrid_itp_path,
            reference_ligand_dir=reference_ligand_path,
            mutant_ligand_dir=mutant_ligand_path,
            seed_gro=seed_gro,
        )

        # Update naming to hybrid outputs for topology/script/template generation
        self._sync_hybrid_assets()

        ligand_seed_pdb = layout.hybrid_dir / "ligand_seed.pdb"

        # Write legs
        if bound_system_dir:
            self._write_bound_leg_from_prism(layout, Path(bound_system_dir))
        else:
            self._write_bound_leg(layout, receptor_path, ligand_seed_pdb)

        if unbound_system_dir:
            self._write_unbound_leg_from_prism(layout, Path(unbound_system_dir))
        else:
            self._write_unbound_leg(layout, ligand_seed_pdb)

        # Write scripts and manifest
        self._replicate_repeat_layout(layout)
        self._validate_scaffold_consistency(layout)
        script_writer.write_root_scripts(layout, self.config)
        self._write_manifest(
            layout,
            receptor_path,
            layout.hybrid_dir,
            extra_sources={
                "hybrid_itp": str(hybrid_itp_path.resolve()),
                "reference_ligand_dir": str(reference_ligand_path.resolve()),
                "mutant_ligand_dir": str(mutant_ligand_path.resolve()) if mutant_ligand_path else None,
            },
        )

        return layout

    def build_from_prism_system(
        self,
        receptor_pdb: str,
        hybrid_itp: str,
        reference_ligand_dir: str,
        mutant_ligand_dir: Optional[str],
        bound_prism_system: str,
        unbound_prism_system: str,
    ) -> FEPScaffoldLayout:
        """
        Build FEP scaffold from PRISM-built systems.

        This method integrates with PRISM-built MD systems, copying coordinate files
        and topologies from GMX_PROLIG_MD directories.

        Parameters
        ----------
        receptor_pdb : str
            Path to receptor PDB file
        hybrid_itp : str
            Path to hybrid ITP file
        reference_ligand_dir : str
            Directory containing reference ligand force field
        mutant_ligand_dir : Optional[str]
            Directory containing mutant ligand force field
        bound_prism_system : str
            Path to bound GMX_PROLIG_MD directory
        unbound_prism_system : str
            Path to unbound GMX_PROLIG_MD directory

        Returns
        -------
        FEPScaffoldLayout
            The created scaffold layout
        """
        return self.build_from_components(
            receptor_pdb=receptor_pdb,
            hybrid_itp=hybrid_itp,
            reference_ligand_dir=reference_ligand_dir,
            mutant_ligand_dir=mutant_ligand_dir,
            bound_system_dir=bound_prism_system,
            unbound_system_dir=unbound_prism_system,
        )

    def _validate_inputs(self, receptor_path: Path, hybrid_dir: Path) -> None:
        """Validate input files exist."""
        if not receptor_path.exists():
            raise FileNotFoundError(f"Receptor PDB not found: {receptor_path}")
        if not hybrid_dir.exists():
            raise FileNotFoundError(f"Hybrid ligand directory not found: {hybrid_dir}")

    def _prepare_layout(self) -> FEPScaffoldLayout:
        """Create directory structure for FEP scaffold."""
        root = self.output_dir
        # Note: Don't check if root.exists() here because previous steps
        # (bound/unbound system building, hybrid topology) may have already
        # created parts of the directory structure. Just ensure required
        # subdirectories exist.

        # Create directory structure
        common_dir = root / "common"
        hybrid_dir = common_dir / "hybrid"
        protein_dir = common_dir / "protein"
        bound_root = root / "bound"
        unbound_root = root / "unbound"
        bound_dir = bound_root / "repeat1"
        unbound_dir = unbound_root / "repeat1"

        for directory in [hybrid_dir, protein_dir, bound_dir / "input", unbound_dir / "input"]:
            directory.mkdir(parents=True, exist_ok=True)

        return FEPScaffoldLayout(
            root=root,
            common_dir=common_dir,
            hybrid_dir=hybrid_dir,
            protein_dir=protein_dir,
            bound_dir=bound_dir,
            unbound_dir=unbound_dir,
        )

    def _replicate_repeat_layout(self, layout: FEPScaffoldLayout) -> None:
        """Create repeat2..N using shared links for static assets."""
        if self.replicas <= 1:
            return

        for leg_name, seed_dir in (("bound", layout.bound_dir), ("unbound", layout.unbound_dir)):
            leg_root = layout.root / leg_name
            for replica_idx in range(2, self.replicas + 1):
                replica_dir = leg_root / f"repeat{replica_idx}"
                if replica_dir.exists():
                    continue
                replica_dir.mkdir(parents=True, exist_ok=True)
                for item in seed_dir.iterdir():
                    if item.name == "build" or item.name.startswith("window_"):
                        continue
                    target = replica_dir / item.name
                    if item.is_dir():
                        target.symlink_to(os.path.relpath(item, replica_dir), target_is_directory=True)
                    else:
                        target.symlink_to(os.path.relpath(item, replica_dir))

    def _copy_common_assets(self, receptor_path: Path, hybrid_dir: Path, layout: FEPScaffoldLayout) -> None:
        """Copy receptor and hybrid ligand assets to common directory."""
        shutil.copy2(receptor_path, layout.protein_dir / "receptor.pdb")

        for asset in self._iter_hybrid_assets(hybrid_dir):
            shutil.copy2(asset, layout.hybrid_dir / asset.name)

        # Generate ligand_seed.pdb from GRO file
        gro_file = hybrid_dir / self.hybrid_gro_filename
        if gro_file.exists():
            (layout.hybrid_dir / "ligand_seed.pdb").write_text(self.hybrid_builder._gro_to_pdb(gro_file))

    def _validate_scaffold_consistency(self, layout: FEPScaffoldLayout) -> None:
        """Fail fast when hybrid assets and staged leg coordinates disagree."""
        hybrid_itp_atoms = self._count_hybrid_itp_atoms(layout.hybrid_dir / self.hybrid_itp_filename)
        hybrid_gro_atoms = self._count_gro_atoms(layout.hybrid_dir / self.hybrid_gro_filename)

        if hybrid_itp_atoms != hybrid_gro_atoms:
            raise ValueError(
                "Hybrid asset mismatch: "
                f"`{layout.hybrid_dir / self.hybrid_itp_filename}` has {hybrid_itp_atoms} atoms, "
                f"but `{layout.hybrid_dir / self.hybrid_gro_filename}` has {hybrid_gro_atoms}."
            )

        for leg_dir in (layout.bound_dir, layout.unbound_dir):
            conf_gro = leg_dir / "input" / "conf.gro"
            ligand_atoms = self._count_leading_residue_atoms(conf_gro)
            if ligand_atoms != hybrid_itp_atoms:
                raise ValueError(
                    "Leg scaffold mismatch: "
                    f"`{conf_gro}` starts with {ligand_atoms} ligand atoms, "
                    f"but hybrid topology defines {hybrid_itp_atoms}."
                )

    @staticmethod
    def _count_hybrid_itp_atoms(itp_path: Path) -> int:
        if not itp_path.exists():
            raise FileNotFoundError(f"Hybrid ITP not found: {itp_path}")

        count = 0
        in_atoms = False
        for line in itp_path.read_text().splitlines():
            stripped = line.strip()
            if stripped.lower() == "[ atoms ]":
                in_atoms = True
                continue
            if in_atoms and stripped.startswith("["):
                break
            if in_atoms and stripped and not stripped.startswith(";"):
                count += 1
        return count

    @staticmethod
    def _count_gro_atoms(gro_path: Path) -> int:
        if not gro_path.exists():
            raise FileNotFoundError(f"GRO file not found: {gro_path}")
        lines = gro_path.read_text().splitlines()
        if len(lines) < 2:
            raise ValueError(f"Invalid GRO file: {gro_path}")
        return int(lines[1].strip())

    @staticmethod
    def _count_leading_residue_atoms(gro_path: Path) -> int:
        if not gro_path.exists():
            raise FileNotFoundError(f"Coordinate file not found: {gro_path}")

        lines = gro_path.read_text().splitlines()
        if len(lines) < 3:
            raise ValueError(f"Invalid GRO file: {gro_path}")

        natoms = int(lines[1].strip())
        atom_lines = lines[2 : 2 + natoms]
        if not atom_lines:
            return 0

        def residue_key(line: str) -> tuple[str, str]:
            resid = line[:5].strip()
            resname = line[5:10].strip()
            if resid or resname:
                return resid, resname
            parts = line.split()
            token = parts[0] if parts else ""
            match = re.match(r"^(\d+)([A-Za-z].*)$", token)
            if match:
                return match.group(1), match.group(2)
            return "", token

        first_key = residue_key(atom_lines[0])
        count = 0
        for line in atom_lines:
            if residue_key(line) != first_key:
                break
            count += 1
        return count

    def _iter_hybrid_assets(self, hybrid_dir: Path) -> Iterable[Path]:
        """Iterate over hybrid ligand asset files."""
        for name in [
            self.hybrid_itp_filename,
            self.hybrid_atomtypes_filename,
            self.hybrid_forcefield_filename,
            self.hybrid_gro_filename,
            self.hybrid_posre_filename,
        ]:
            path = hybrid_dir / name
            if path.exists():
                yield path

    def _write_bound_leg(self, layout: FEPScaffoldLayout, receptor_path: Path, ligand_seed_pdb: Path) -> None:
        """Write bound leg directory structure."""
        bound_input = layout.bound_dir / "input"
        shutil.copy2(layout.protein_dir / "receptor.pdb", bound_input / "receptor.pdb")
        shutil.copy2(ligand_seed_pdb, bound_input / "ligand_seed.pdb")
        self.leg_writer.write_complex_seed(
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

    def _write_unbound_leg(self, layout: FEPScaffoldLayout, ligand_seed_pdb: Path) -> None:
        """Write unbound leg directory structure."""
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

    def _write_bound_leg_from_prism(self, layout: FEPScaffoldLayout, prism_system_dir: Path) -> None:
        """Write bound leg from PRISM-built system."""
        self.leg_writer.copy_prism_system_to_leg(prism_system_dir, layout.bound_dir, "bound")

        write_fep_mdps(
            output_dir=str(layout.bound_dir / "mdps"),
            lambda_strategy=self.lambda_strategy,
            lambda_distribution=self.lambda_distribution,
            lambda_windows=self.lambda_windows,
            config=self.config,
            leg_name="bound",
        )

    def _write_unbound_leg_from_prism(self, layout: FEPScaffoldLayout, prism_system_dir: Path) -> None:
        """Write unbound leg from PRISM-built system."""
        self.leg_writer.copy_prism_system_to_leg(prism_system_dir, layout.unbound_dir, "unbound")

        write_fep_mdps(
            output_dir=str(layout.unbound_dir / "mdps"),
            lambda_strategy=self.lambda_strategy,
            lambda_distribution=self.lambda_distribution,
            lambda_windows=self.lambda_windows,
            config=self.config,
            leg_name="unbound",
        )

    def _create_placeholder_system(
        self, leg_name: str, layout: FEPScaffoldLayout, protein_path: Optional[Path] | None = None
    ) -> None:
        """Create placeholder topology and coordinate files.

        Parameters
        ----------
        leg_name : str
            Name of the leg ('bound' or 'unbound')
        layout : FEPScaffoldLayout
            Scaffold layout
        protein_path : Optional[Path]
            Not used in placeholder creation (kept for API compatibility)
        """
        _ = protein_path  # Mark as intentionally unused for compatibility
        """Create placeholder topology and coordinate files."""
        leg_dir = layout.bound_dir if leg_name == "bound" else layout.unbound_dir

        # Write placeholder topology
        if leg_name == "bound":
            self._write_bound_topology_template(leg_dir / "topol.top")
        else:
            self._write_unbound_topology(leg_dir / "topol.top")

        # Write placeholder coordinate file
        placeholder_gro = leg_dir / "input" / "conf.gro"
        if not placeholder_gro.exists():
            hybrid_gro = layout.hybrid_dir / self.hybrid_gro_filename
            if leg_name == "unbound" and hybrid_gro.exists():
                shutil.copy2(hybrid_gro, placeholder_gro)
            else:
                with open(placeholder_gro, "w") as f:
                    f.write("Placeholder system\n")
                    f.write("1\n")
                    f.write("    1HYB     C1    1   0.000   0.000   0.000\n")
                    f.write("   5.00000   5.00000   5.00000\n")

    def _write_manifest(
        self,
        layout: FEPScaffoldLayout,
        receptor_path: Path,
        hybrid_dir: Path,
        extra_sources: Optional[Dict[str, Optional[str]]] = None,
    ) -> None:
        """Write FEP scaffold manifest file."""
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

    def _write_unbound_topology(self, topol_path: Path) -> None:
        """Write unbound leg topology template."""
        hybrid_prefix = os.path.relpath(self.layout.common_hybrid_dir, topol_path.parent).replace("\\", "/")
        content = "; PRISM-FEP unbound topology template\n"
        defaults_path = topol_path.parent.parent / "common" / "hybrid" / "defaults_hybrid.itp"
        ff_path = topol_path.parent.parent / "common" / "hybrid" / self.hybrid_forcefield_filename
        atomtypes_path = topol_path.parent.parent / "common" / "hybrid" / self.hybrid_atomtypes_filename
        if defaults_path.exists():
            content += f'#include "{hybrid_prefix}/defaults_hybrid.itp"\n'
        if ff_path.exists():
            content += f'#include "{hybrid_prefix}/{self.hybrid_forcefield_filename}"\n'
        elif atomtypes_path.exists():
            content += f'#include "{hybrid_prefix}/{self.hybrid_atomtypes_filename}"\n'
        content += f'#include "{hybrid_prefix}/{self.hybrid_itp_filename}"\n\n'
        content += "[ system ]\nUnbound FEP leg\n\n"
        content += f"[ molecules ]\n{self.molecule_name} 1\n"
        topol_path.write_text(content)

    def _write_bound_topology_template(self, topol_path: Path) -> None:
        """Write bound leg topology template."""
        hybrid_prefix = os.path.relpath(self.layout.common_hybrid_dir, topol_path.parent).replace("\\", "/")
        content = "; PRISM-FEP bound topology template\n"
        defaults_path = topol_path.parent.parent / "common" / "hybrid" / "defaults_hybrid.itp"
        ff_path = topol_path.parent.parent / "common" / "hybrid" / self.hybrid_forcefield_filename
        atomtypes_path = topol_path.parent.parent / "common" / "hybrid" / self.hybrid_atomtypes_filename
        if defaults_path.exists():
            content += f'#include "{hybrid_prefix}/defaults_hybrid.itp"\n'
        if ff_path.exists():
            content += f'#include "{hybrid_prefix}/{self.hybrid_forcefield_filename}"\n'
        elif atomtypes_path.exists():
            content += f'#include "{hybrid_prefix}/{self.hybrid_atomtypes_filename}"\n'
        content += f'#include "{hybrid_prefix}/{self.hybrid_itp_filename}"\n\n'
        content += "[ system ]\nBound FEP leg\n\n"
        content += f"[ molecules ]\nProtein      1\n{self.molecule_name} 1\n"
        topol_path.write_text(content)
