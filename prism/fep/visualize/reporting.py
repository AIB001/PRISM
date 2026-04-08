#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mapping Report Service for FEP

This module provides a service layer for generating FEP mapping visualization reports.
Extracted from workflow_fep.py lines 377-503, 425-503, 505-583.
"""

import os
from pathlib import Path
from typing import Optional, Tuple

from prism.fep.common.coordinates import nm_xyz_to_angstrom


class MappingReportService:
    """Service for generating FEP mapping visualization reports."""

    def __init__(
        self,
        distance_cutoff: float = 0.6,
        charge_cutoff: float = 0.05,
        charge_strategy: str = "mean",
        charge_reception: str = "surround",
    ):
        """Initialize the mapping report service.

        Parameters
        ----------
        distance_cutoff : float, optional
            Distance cutoff for atom mapping (default: 0.6 nm)
        charge_cutoff : float, optional
            Charge cutoff for atom mapping (default: 0.05 e)
        charge_strategy : str, optional
            Strategy for handling common atom charges (default: "mean")
        charge_reception : str, optional
            Strategy for handling transformed atom charges (default: "surround")
        """
        self.distance_cutoff = distance_cutoff
        self.charge_cutoff = charge_cutoff
        self.charge_strategy = charge_strategy
        self.charge_reception = charge_reception

    def generate_fep_mapping_html(
        self,
        output_dir: str,
        ref_coord_source: str,
        mut_coord_source: str,
        mapping,
        ref_atoms,
        mut_atoms,
        fep_config=None,
        config_file=None,
        ligand_forcefield="unknown",
        forcefield_paths=None,
        distance_cutoff=None,
        charge_cutoff=None,
        charge_strategy=None,
        charge_reception=None,
    ) -> None:
        """Generate interactive mapping HTML report.

        Extracted from workflow_fep.py lines 377-503.

        Parameters
        ----------
        output_dir : str
            Output directory for mapping HTML
        ref_coord_source : str
            Path to reference coordinate source
        mut_coord_source : str
            Path to mutant coordinate source
        mapping : AtomMapping
            Atom mapping result
        ref_atoms : List
            Reference ligand atoms
        mut_atoms : List
            Mutant ligand atoms
        fep_config : Optional[str], optional
            Path to FEP config file
        config_file : Optional[str], optional
            Path to main config file
        ligand_forcefield : str, optional
            Force field type for display
        distance_cutoff : Optional[float], optional
            Override distance cutoff
        charge_cutoff : Optional[float], optional
            Override charge cutoff
        charge_strategy : Optional[str], optional
            Override charge strategy
        charge_reception : Optional[str], optional
            Override charge reception
        """
        try:
            from prism.fep.visualize.html import visualize_mapping_html
        except Exception as exc:
            print(f"Skipping FEP mapping HTML generation: {exc}")
            return

        html_path = os.path.join(output_dir, "mapping.html")

        from prism.fep.common.io import read_mol2_atoms, write_ligand_to_pdb

        html_mapping = mapping
        html_atoms_a = ref_atoms
        html_atoms_b = mut_atoms

        ref_pdb = os.path.join(output_dir, "mapping_state_a.pdb")
        mut_pdb = os.path.join(output_dir, "mapping_state_b.pdb")
        write_ligand_to_pdb(html_atoms_a, ref_pdb, residue_name="LIG")
        write_ligand_to_pdb(html_atoms_b, mut_pdb, residue_name="LIG")

        # Detect MOL2 files
        ref_mol2 = self._detect_mol2_file(ref_coord_source)
        mut_mol2 = self._detect_mol2_file(mut_coord_source)

        # For generic force fields (OPLS/OpenFF), try MOL2-derived mapping
        if ref_mol2 and mut_mol2:
            try:
                from prism.fep.core.mapping import DistanceAtomMapper

                html_mapper = DistanceAtomMapper(
                    dist_cutoff=distance_cutoff or self.distance_cutoff,
                    charge_cutoff=charge_cutoff or self.charge_cutoff,
                    charge_common=charge_strategy or self.charge_strategy,
                    charge_reception=charge_reception or self.charge_reception,
                )

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
                        f"  Using MOL2-derived mapping for HTML "
                        f"({len(mol2_mapping.common)} common vs {len(mapping.common)})"
                    )
            except Exception as exc:
                print(f"  Could not build MOL2-derived HTML mapping: {exc}")

        # Load HTML config
        html_config = {}
        if fep_config:
            try:
                from prism.fep.common.config import FEPConfig

                html_config = FEPConfig(config_file, fep_config).get_html_config()
            except Exception as exc:
                print(f"  Could not load FEP HTML config: {exc}")

        # Format force field type for display
        ff_cfg = dict(html_config.get("forcefield", {}))
        ff_type = ff_cfg.get("type", ligand_forcefield)
        ff_cfg["type"] = self._describe_mapping_forcefield(ff_type, forcefield_paths, ligand_forcefield)
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
            print(f"  ✓ Mapping HTML: {html_path}")
        except Exception as exc:
            print(f"  ⚠ Failed to generate FEP mapping HTML: {exc}")

    def prepare_visualization_inputs(self, coord_source: str, output_dir: str, stem: str) -> Tuple[str, Optional[str]]:
        """Return a PDB path plus optional MOL2 template for mapping HTML generation.

        Extracted from workflow_fep.py lines 425-503.

        Parameters
        ----------
        coord_source : str
            Path to coordinate source file
        output_dir : str
            Output directory for generated PDB
        stem : str
            Stem for output filename

        Returns
        -------
        Tuple[str, Optional[str]]
            (PDB path, MOL2 path or None)

        Raises
        ------
        RuntimeError
            If coordinate source format is not supported
        """
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
            # Convert GRO to PDB using simple conversion
            pdb_path = Path(output_dir) / f"{stem}.pdb"
            with open(source_path) as f:
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
                        z_raw = float(parts[5]) if len(parts) > 5 else 0.0
                        x, y, z = nm_xyz_to_angstrom(float(parts[3]), float(parts[4]), z_raw)

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
            return str(pdb_path), None

        raise RuntimeError(f"Unsupported coordinate source for mapping visualization: {coord_source}")

    def generate_state_specific_pdb(self, gro_file: str, output_dir: str, stem: str, state: str = "a") -> str:
        """Generate state-specific PDB file from hybrid GRO file.

        Extracted from workflow_fep.py lines 505-583.

        Filters atoms by state:
        - State "a": Includes common atoms + transformed atoms ending with "A"
        - State "b": Includes common atoms + transformed atoms ending with "B"

        Common atoms don't follow the digit+A/B pattern (e.g., C0A, C02, H11).

        Parameters
        ----------
        gro_file : str
            Path to hybrid GRO file
        output_dir : str
            Output directory
        stem : str
            Output file stem
        state : str, optional
            "a" for state A, "b" for state B (default: "a")

        Returns
        -------
        str
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
                    z_raw = float(parts[5]) if len(parts) > 5 else 0.0
                    x, y, z = nm_xyz_to_angstrom(float(parts[3]), float(parts[4]), z_raw)

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

    def _detect_mol2_file(self, coord_source: str) -> Optional[str]:
        """Detect MOL2 file for a given coordinate source.

        Parameters
        ----------
        coord_source : str
            Path to coordinate source file

        Returns
        -------
        Optional[str]
            Path to MOL2 file if found, None otherwise
        """
        path = Path(coord_source)

        if path.suffix == ".mol2":
            return str(path.absolute()) if path.exists() else coord_source
        elif path.suffix == ".pdb":
            mol2_path = path.with_suffix(".mol2")
            return str(mol2_path) if mol2_path.exists() else None

        return None

    def _describe_mapping_forcefield(self, ff_type: str, forcefield_paths=None, ligand_forcefield="") -> str:
        """Return a user-facing force-field label for mapping HTML.

        Restored from workflow_fep.py:395-409 to preserve CGenFF distinction.

        Parameters
        ----------
        ff_type : str
            Force field type
        forcefield_paths : Optional[List[str]], optional
            List of force field paths for CGenFF detection
        ligand_forcefield : str, optional
            Ligand force field name

        Returns
        -------
        str
            User-facing force field label with CGenFF source distinction
        """
        base = str(ff_type).upper()
        lig_ff_lower = (ligand_forcefield or "").lower()

        # Handle both "cgenff" and "charmm-gui" force field types
        if lig_ff_lower not in ["cgenff", "charmm-gui"]:
            return base

        # Distinguish CGenFF sources
        if forcefield_paths:
            ff_paths = [Path(p) for p in forcefield_paths]
            # Check for CHARMM-GUI style (has charmm36.ff/charmm36.itp)
            if ff_paths and all((path / "charmm36.ff" / "charmm36.itp").exists() for path in ff_paths):
                return "CGENFF (CHARMM-GUI)"
            # Check for CHARMM-GUI style (has gromacs/LIG.itp)
            elif ff_paths and all((path / "gromacs" / "LIG.itp").exists() for path in ff_paths):
                return "CGENFF (CHARMM-GUI)"

        # Default: website-generated CGenFF
        return "CGENFF (WEBSITE)"
