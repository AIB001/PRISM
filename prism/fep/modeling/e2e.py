#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
End-to-end FEP system building convenience functions.

This module provides high-level functions that automate the complete workflow
from ligand files to FEP simulation systems.
"""

from pathlib import Path
from typing import Optional

from prism.fep.modeling.core import FEPScaffoldBuilder
from prism.fep.modeling.hybrid_service import HybridBuildService


def build_fep_system_from_prism_ligands(
    receptor_pdb: str,
    reference_ligand_dir: str,
    mutant_ligand_dir: str,
    output_dir: str = "fep_output",
    hybrid_itp: Optional[str] = None,
    reference_coord_file: Optional[str] = None,
    mutant_coord_file: Optional[str] = None,
    lambda_windows: int = 32,
    lambda_strategy: str = "decoupled",
    dist_cutoff: float = 0.6,
    charge_cutoff: float = 0.05,
) -> Path:
    """
    Build complete FEP system from PRISM-formatted ligands.

    This is a convenience function that combines:
    1. Reading PRISM ligands
    2. Performing atom mapping
    3. Building FEP scaffold

    Parameters
    ----------
    receptor_pdb : str
        Path to receptor PDB file
    reference_ligand_dir : str
        Directory containing reference ligand PRISM files (LIG.itp, LIG.gro)
    mutant_ligand_dir : str
        Directory containing mutant ligand PRISM files (LIG.itp, LIG.gro)
    output_dir : str, optional
        Output directory for FEP system (default: "fep_output")
    hybrid_itp : Optional[str], optional
        Path to existing hybrid ITP file (if None, will perform mapping)
    lambda_windows : int, optional
        Number of lambda windows (default: 32)
    lambda_strategy : str, optional
        Lambda strategy (default: "decoupled")
    dist_cutoff : float, optional
        Distance cutoff for atom mapping (default: 0.6)
    charge_cutoff : float, optional
        Charge cutoff for atom mapping (default: 0.05)

    Returns
    -------
    Path
        Path to the created FEP output directory

    Examples
    --------
    >>> layout = build_fep_system_from_prism_ligands(
    ...     receptor_pdb="receptor.pdb",
    ...     reference_ligand_dir="ligand_a_ff",
    ...     mutant_ligand_dir="ligand_b_ff",
    ... )
    """
    hybrid_service = HybridBuildService(
        dist_cutoff=dist_cutoff,
        charge_cutoff=charge_cutoff,
        charge_common="mean",
        charge_reception="surround",
    )

    # Step 1: Build or validate hybrid topology
    if hybrid_itp is None:
        print(f"[1/3] Building hybrid topology from mapping...")

        temp_hybrid_dir = FEPScaffoldBuilder._normalize_output_dir(output_dir) / "temp_hybrid"
        temp_hybrid_dir.mkdir(parents=True, exist_ok=True)

        hybrid_result = hybrid_service.build_from_forcefield_dirs(
            reference_ligand_dir=reference_ligand_dir,
            mutant_ligand_dir=mutant_ligand_dir,
            hybrid_output_dir=str(temp_hybrid_dir),
            reference_coord_source=reference_coord_file,
            mutant_coord_source=mutant_coord_file,
            molecule_name="HYB",
        )

        # Use the generated hybrid ITP
        hybrid_itp = str(hybrid_result.hybrid_itp)
        print(f"  ✓ Hybrid ITP: {hybrid_itp} (generated from mapping)")
        print(f"  ✓ Common: {len(hybrid_result.mapping.common)}")
        print(f"  ✓ Transformed (ref): {len(hybrid_result.mapping.transformed_a)}")
        print(f"  ✓ Transformed (mut): {len(hybrid_result.mapping.transformed_b)}")
        print(f"[2/3] Hybrid topology ready")
    else:
        print(f"[1/3] Using existing hybrid ITP: {hybrid_itp}")
        print(f"[2/3] Hybrid topology ready")

    # Step 3: Build FEP scaffold
    print(f"[3/3] Building FEP scaffold...")
    builder = FEPScaffoldBuilder(
        output_dir=output_dir,
        lambda_windows=lambda_windows,
        lambda_strategy=lambda_strategy,
    )

    layout = builder.build_from_components(
        receptor_pdb=receptor_pdb,
        hybrid_itp=hybrid_itp,
        reference_ligand_dir=reference_ligand_dir,
        mutant_ligand_dir=mutant_ligand_dir,
    )

    print(f"  ✓ Output directory: {layout.root}")
    print(f"  ✓ Bound leg: {layout.bound_dir}")
    print(f"  ✓ Unbound leg: {layout.unbound_dir}")
    print(f"\n✅ FEP system built successfully!")

    return layout.root
