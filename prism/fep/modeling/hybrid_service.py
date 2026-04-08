#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Hybrid topology orchestration services for FEP workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from prism.fep.common.io import read_ligand_from_prism
from prism.fep.core.hybrid_topology import HybridTopologyBuilder, LigandTopologyInput
from prism.fep.core.mapping import AtomMapping, DistanceAtomMapper
from prism.fep.gromacs.itp_builder import ITPBuilder
from prism.fep.modeling.hybrid_package import HybridPackageBuilder
from prism.forcefield.registry import discover_ligand_directory, resolve_ligand_artifact


@dataclass
class HybridBuildResult:
    """Artifacts and metadata produced by hybrid topology construction."""

    reference_itp: Path
    reference_gro: Path
    mutant_itp: Path
    mutant_gro: Path
    reference_coord_source: Path
    mutant_coord_source: Path
    mapping: AtomMapping
    ref_atoms: list
    mut_atoms: list
    hybrid_atoms: list
    hybrid_itp: Path
    hybrid_gro: Path
    hybrid_params: dict


class HybridBuildService:
    """Build hybrid topology artifacts from two PRISM ligand force-field directories."""

    def __init__(
        self,
        dist_cutoff: float = 0.6,
        charge_cutoff: float = 0.05,
        charge_common: str = "mean",
        charge_reception: str = "surround",
    ) -> None:
        self.dist_cutoff = dist_cutoff
        self.charge_cutoff = charge_cutoff
        self.charge_common = charge_common
        self.charge_reception = charge_reception

    @staticmethod
    def resolve_prism_ligand_dir(ligand_dir: str | Path) -> Path:
        """Resolve a PRISM ligand directory, including nested force-field layouts."""
        ligand_path = Path(ligand_dir)
        discovered_dir = discover_ligand_directory(ligand_path)
        if discovered_dir is not None:
            return discovered_dir
        return ligand_path

    @staticmethod
    def resolve_coord_source(path: str | Path | None, fallback: str | Path) -> Path:
        """Resolve an aligned coordinate source, preferring PDBs inside directory inputs."""
        fallback_path = Path(fallback)
        if path:
            source = Path(path)
            if source.is_dir():
                pdbs = list(source.glob("*.pdb"))
                if pdbs:
                    return pdbs[0]
            elif source.exists():
                return source
        return fallback_path

    @staticmethod
    def _resolve_required_artifact(ligand_dir: Path, filename: str) -> Path:
        resolved = resolve_ligand_artifact(ligand_dir, filename)
        if resolved is None:
            raise FileNotFoundError(f"Could not find {filename} under ligand FF directory: {ligand_dir}")
        return resolved

    @staticmethod
    def write_hybrid_gro(
        output_gro: str | Path,
        *,
        pkg_builder: HybridPackageBuilder,
        hybrid_itp: Path,
        reference_gro: Path,
        mutant_gro: Path,
        reference_coords_file: Path,
        mutant_coords_file: Path,
        reference_itp: Path,
        mutant_itp: Path,
    ) -> Path:
        """Write a hybrid GRO by reusing the canonical HybridPackageBuilder implementation."""
        output_path = Path(output_gro)
        gro_content = pkg_builder._build_hybrid_gro(
            hybrid_itp_content=hybrid_itp.read_text(),
            reference_gro=reference_gro,
            mutant_gro=mutant_gro,
            reference_coords_file=reference_coords_file,
            mutant_coords_file=mutant_coords_file,
            reference_itp=reference_itp,
            mutant_itp=mutant_itp,
        )
        output_path.write_text(gro_content)
        return output_path

    def build_from_forcefield_dirs(
        self,
        reference_ligand_dir: str | Path,
        mutant_ligand_dir: str | Path,
        hybrid_output_dir: str | Path,
        reference_coord_source: Optional[str | Path] = None,
        mutant_coord_source: Optional[str | Path] = None,
        molecule_name: str = "HYB",
    ) -> HybridBuildResult:
        """Build hybrid topology artifacts from two prepared ligand force-field directories."""
        ref_dir = self.resolve_prism_ligand_dir(reference_ligand_dir)
        mut_dir = self.resolve_prism_ligand_dir(mutant_ligand_dir)
        hybrid_dir = Path(hybrid_output_dir)
        hybrid_dir.mkdir(parents=True, exist_ok=True)

        ref_itp = self._resolve_required_artifact(ref_dir, "LIG.itp")
        ref_gro = self._resolve_required_artifact(ref_dir, "LIG.gro")
        mut_itp = self._resolve_required_artifact(mut_dir, "LIG.itp")
        mut_gro = self._resolve_required_artifact(mut_dir, "LIG.gro")

        ref_coord_file = self.resolve_coord_source(reference_coord_source, ref_gro)
        mut_coord_file = self.resolve_coord_source(mutant_coord_source, mut_gro)

        ref_atoms = read_ligand_from_prism(itp_file=str(ref_itp), gro_file=str(ref_coord_file), state="a")
        mut_atoms = read_ligand_from_prism(itp_file=str(mut_itp), gro_file=str(mut_coord_file), state="b")

        mapper = DistanceAtomMapper(
            dist_cutoff=self.dist_cutoff,
            charge_cutoff=self.charge_cutoff,
            charge_common=self.charge_common,
            charge_reception=self.charge_reception,
        )
        mapping = mapper.map(ref_atoms, mut_atoms)

        ref_itp_data = ITPBuilder._parse_source_itp(ref_itp.read_text())
        mut_itp_data = ITPBuilder._parse_source_itp(mut_itp.read_text())
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

        hybrid_builder = HybridTopologyBuilder(charge_strategy=self.charge_common)
        hybrid_atoms = hybrid_builder.build(mapping, params_a, params_b, ref_atoms, mut_atoms)

        temp_itp = hybrid_dir / "hybrid_atoms_temp.itp"
        ITPBuilder(hybrid_atoms, {}).write_itp(str(temp_itp), molecule_name=molecule_name)

        pkg_builder = HybridPackageBuilder()
        atom_types_a = pkg_builder._parse_atom_types_from_itp(str(ref_itp))
        atom_types_b = pkg_builder._parse_atom_types_from_itp(str(mut_itp))

        hybrid_itp = hybrid_dir / "hybrid.itp"
        hybrid_params = ITPBuilder.write_complete_hybrid_itp(
            output_path=str(hybrid_itp),
            hybrid_itp=str(temp_itp),
            ligand_a_itp=str(ref_itp),
            ligand_b_itp=str(mut_itp),
            molecule_name=molecule_name,
            ligand_a_structure=str(ref_coord_file),
            ligand_b_structure=str(mut_coord_file),
            atom_types_a=atom_types_a,
            atom_types_b=atom_types_b,
        )

        hybrid_gro = self.write_hybrid_gro(
            hybrid_dir / "hybrid.gro",
            pkg_builder=pkg_builder,
            hybrid_itp=hybrid_itp,
            reference_gro=ref_gro,
            mutant_gro=mut_gro,
            reference_coords_file=ref_coord_file,
            mutant_coords_file=mut_coord_file,
            reference_itp=ref_itp,
            mutant_itp=mut_itp,
        )

        return HybridBuildResult(
            reference_itp=ref_itp,
            reference_gro=ref_gro,
            mutant_itp=mut_itp,
            mutant_gro=mut_gro,
            reference_coord_source=ref_coord_file,
            mutant_coord_source=mut_coord_file,
            mapping=mapping,
            ref_atoms=ref_atoms,
            mut_atoms=mut_atoms,
            hybrid_atoms=hybrid_atoms,
            hybrid_itp=hybrid_itp,
            hybrid_gro=hybrid_gro,
            hybrid_params=hybrid_params,
        )
