#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Services for building ligand-only FEP preparation systems."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from prism.forcefield.adapters import CGenFFAdapter
from prism.forcefield.registry import resolve_ligand_artifact
from prism.utils.colors import print_step, print_success
from prism.utils.system.topology import _parse_atomtypes, _should_vendor_forcefield


@dataclass
class LigandOnlySystemResult:
    """Artifacts produced while building an unbound ligand-only system."""

    model_dir: Path
    topology: Path
    ligand_ff_dir: Path
    ligand_gro: Path
    boxed_gro: Path
    solvated_gro: str


class LigandOnlySystemBuilder:
    """Build the solvated ligand-only system used by the FEP unbound leg."""

    def __init__(
        self,
        *,
        system_builder,
        forcefield: dict,
        water_model: dict,
        gromacs_env=None,
        overwrite: bool = False,
        config: Optional[dict] = None,
        cgenff_supplement_builder: Optional[Callable[..., Optional[Path]]] = None,
    ) -> None:
        self.system_builder = system_builder
        self.forcefield = forcefield
        self.water_model = water_model
        self.gromacs_env = gromacs_env
        self.overwrite = overwrite
        self.config = config or {}
        self.cgenff_supplement_builder = cgenff_supplement_builder

    @staticmethod
    def _read_moleculetype_name(itp_path: Path) -> str:
        """Read the first moleculetype name defined in a ligand ITP."""
        in_moleculetype = False
        for raw_line in itp_path.read_text().splitlines():
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.lower() == "[ moleculetype ]":
                in_moleculetype = True
                continue
            if in_moleculetype:
                return line.split()[0]
        raise ValueError(f"Could not parse moleculetype from {itp_path}")

    @staticmethod
    def _resolve_required_artifact(ligand_ff_dir: Path, filename: str) -> Path:
        resolved = resolve_ligand_artifact(ligand_ff_dir, filename)
        if resolved is None:
            raise FileNotFoundError(f"Could not find {filename} under ligand FF directory: {ligand_ff_dir}")
        return resolved

    def _resolve_forcefield_info(self) -> tuple[str, str, Optional[int], Optional[dict], str]:
        """Resolve the active force-field metadata and ions include to use."""
        ff_name = self.forcefield["name"]
        water_model = self.water_model["name"]
        ff_idx = self.forcefield.get("index")
        ff_info = None
        use_ions_itp = "ions.itp"

        if not ff_idx and self.gromacs_env:
            ff_name_lower = ff_name.lower()
            if ff_name_lower in self.gromacs_env.force_field_names:
                ff_idx = self.gromacs_env.force_field_names[ff_name_lower]

        if ff_idx and self.gromacs_env and ff_idx in self.gromacs_env.force_fields:
            ff_info = self.gromacs_env.force_fields[ff_idx]

        if ff_info and "path" in ff_info:
            ff_path = Path(ff_info["path"])
            if ff_path.exists() and ff_path.is_dir():
                water_specific_ions = f"ions_{water_model}.itp"
                generic_ions = "ions.itp"
                if (ff_path / water_specific_ions).exists():
                    use_ions_itp = water_specific_ions
                    print(f"  Using water-specific ion file: {water_specific_ions}")
                elif (ff_path / generic_ions).exists():
                    use_ions_itp = generic_ions
                    print(f"  Using generic ion file: {generic_ions}")
                else:
                    print(f"  ⚠ Warning: Neither {water_specific_ions} nor {generic_ions} found in {ff_path}")
                    print(f"  Will attempt to use {use_ions_itp} (may fail at grompp)")

        return ff_name, water_model, ff_idx, ff_info, use_ions_itp

    def _vendor_forcefield_if_needed(
        self,
        *,
        model_dir: Path,
        ff_name: str,
        ff_info: Optional[dict],
        use_ions_itp: str,
    ) -> None:
        """Copy or validate the configured GROMACS force field for local use."""
        if not (ff_info and "dir" in ff_info):
            return

        ff_basename = ff_info["dir"]
        local_ff_dir = model_dir / ff_basename
        ff_path = Path(ff_info.get("path", ""))

        if not ff_path.exists() or not ff_path.is_dir():
            print(f"Warning: Force field path not found: {ff_path}")
            print(f"Using force field from system search paths: {ff_name}")
            return

        required_files = ["forcefield.itp", "ffbonded.itp", "ffnonbonded.itp", use_ions_itp]
        missing = [f for f in required_files if not (ff_path / f).exists()]
        if missing:
            raise RuntimeError(
                f"Force field is incomplete: {ff_info['name']}\n"
                f"Missing files: {', '.join(missing)}\n"
                f"Force field path: {ff_path}\n"
                f"Source: {ff_info.get('source', 'unknown')}\n"
                f"Please check your force field installation."
            )

        if not _should_vendor_forcefield(ff_info, model_dir):
            print(f"Using installed force field in place: {ff_path}")
            return

        need_copy = not local_ff_dir.exists()
        if not need_copy and self.overwrite:
            missing_local = [f for f in required_files if not (local_ff_dir / f).exists()]
            if missing_local:
                print("Existing force field copy is incomplete, replacing...")
                shutil.rmtree(local_ff_dir)
                need_copy = True

        if need_copy:
            shutil.copytree(ff_path, local_ff_dir)
            source_str = ff_info.get("source", ff_path)
            print(f"Copied force field to working directory: {local_ff_dir}")
            print(f"  Source: {source_str}")
        else:
            print(f"Using existing force field copy: {local_ff_dir}")

    def build(
        self, output_dir: str | Path, ligand_ff_dir: str | Path, box_size: tuple | None = None
    ) -> LigandOnlySystemResult:
        """Build a ligand-only solvated/ionized system for the unbound leg."""
        output_path = Path(output_dir)
        model_dir = output_path / "GMX_PROLIG_MD"
        model_dir.mkdir(exist_ok=True, parents=True)

        ligand_ff_path = Path(ligand_ff_dir)
        ligand_gro = self._resolve_required_artifact(ligand_ff_path, "LIG.gro")

        print_step(2, 5, "Creating ligand-only topology")

        ff_name, water_model, _ff_idx, ff_info, use_ions_itp = self._resolve_forcefield_info()
        self._vendor_forcefield_if_needed(
            model_dir=model_dir,
            ff_name=ff_name,
            ff_info=ff_info,
            use_ions_itp=use_ions_itp,
        )

        topol_path = model_dir / "topol.top"
        ligand_rel_dir = os.path.relpath(ligand_ff_path, model_dir)

        cgenff_supplement = None
        main_atomtypes = set()
        if CGenFFAdapter.is_charmm_forcefield(ff_name) and CGenFFAdapter.detect(ligand_ff_path):
            main_ff_dir = model_dir / f"{ff_name}.ff"
            if not main_ff_dir.exists() and ff_info and ff_info.get("path"):
                main_ff_dir = Path(ff_info["path"])
            ffnonbonded = main_ff_dir / "ffnonbonded.itp"
            if ffnonbonded.exists():
                main_atomtypes = _parse_atomtypes(ffnonbonded)
            if self.cgenff_supplement_builder is not None:
                charmm_ff_dir = CGenFFAdapter.find_charmm_ff_dir(ligand_ff_path)
                cgenff_supplement = self.cgenff_supplement_builder(
                    lig_ff_path=ligand_ff_path,
                    lig_itp_path=self._resolve_required_artifact(ligand_ff_path, "LIG.itp"),
                    charmm_ff_dir=charmm_ff_dir,
                    main_bonded_files=[
                        p for p in [main_ff_dir / "ffbonded.itp", main_ff_dir / "ffmissingdihedrals.itp"] if p.exists()
                    ],
                    main_nonbonded_files=[p for p in [ffnonbonded] if p.exists()],
                )

        atomtypes_lines = []
        atomtypes_itp_path = resolve_ligand_artifact(ligand_ff_path, "atomtypes_LIG.itp")
        include_ligand_atomtypes = CGenFFAdapter.should_include_ligand_atomtypes(ff_name, ligand_ff_path)
        if include_ligand_atomtypes and atomtypes_itp_path is not None and atomtypes_itp_path.exists():
            for line in atomtypes_itp_path.read_text().splitlines():
                stripped = line.strip()
                if not stripped or stripped.startswith(";") or stripped.startswith("["):
                    continue
                atomtype = stripped.split()[0]
                if atomtype not in main_atomtypes:
                    atomtypes_lines.append(line)

        ligand_itp = self._resolve_required_artifact(ligand_ff_path, "LIG.itp")
        ligand_mol_name = self._read_moleculetype_name(ligand_itp)

        with open(topol_path, "w") as handle:
            handle.write("; Ligand-only topology for FEP unbound leg\n")
            handle.write(f'#include "{ff_name}.ff/forcefield.itp"\n')
            if atomtypes_lines:
                handle.write("\n; Include ligand-specific atomtypes\n")
                handle.write("[ atomtypes ]\n")
                for line in atomtypes_lines:
                    handle.write(f"{line}\n")
                handle.write("\n")
            if cgenff_supplement is not None:
                handle.write(f'#include "{ligand_rel_dir}/{cgenff_supplement.name}"\n')
            handle.write(f'#include "{ligand_rel_dir}/LIG.itp"\n')
            handle.write(f'#include "{ff_name}.ff/{water_model}.itp"\n')
            if use_ions_itp:
                handle.write(f'#include "{ff_name}.ff/{use_ions_itp}"\n\n')
            else:
                handle.write(f"; WARNING: No ion file found for {ff_name} with water model {water_model}\n")
                handle.write(f'#include "{ff_name}.ff/ions.itp"  ; This will likely fail\n\n')
            handle.write("[ system ]\n")
            handle.write("Ligand in water\n\n")
            handle.write("[ molecules ]\n")
            handle.write(f"{ligand_mol_name:<8} 1\n")

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
        return LigandOnlySystemResult(
            model_dir=model_dir,
            topology=topol_path,
            ligand_ff_dir=ligand_ff_path,
            ligand_gro=ligand_gro,
            boxed_gro=boxed_gro,
            solvated_gro=solvated_gro,
        )
