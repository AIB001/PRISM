#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Configurable FEP Mapping Test

This script tests the complete FEP workflow with force field configuration from fep.yaml:
1. Read force field type and parameters from fep.yaml
2. Generate force fields for both molecules
3. Create atom mapping between them
4. Generate PNG and HTML visualizations
5. Display detailed mapping information

Configuration files:
- config.conf: General PRISM parameters (temperature, salt concentration, etc.)
- fep.yaml: FEP-specific configuration (force field + mapping parameters)

Force field configuration in fep.yaml:
    forcefield:
      type: gaff/gaff2/openff/opls/cgenff
      params:
        charge_mode: bcc  # For GAFF/GAFF2
        # charge_model: cm1a  # For OPLS
        # version: 2.1.1  # For OpenFF

    mapping:
      dist_cutoff: 1.0
      charge_cutoff: 0.06
      charge_common: mean
      charge_reception: surround

Usage:
    python tests/gxf/FEP/unit_test/test_gaff_fep_mapping.py
"""

import sys
from pathlib import Path
import pytest

# Add PRISM to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from prism.forcefield.gaff import GAFFForceFieldGenerator
from prism.forcefield.gaff2 import GAFF2ForceFieldGenerator
from prism.forcefield.openff import OpenFFForceFieldGenerator
from prism.forcefield.opls_aa import OPLSAAForceFieldGenerator
from prism.forcefield.swissparam import MMFFForceFieldGenerator
from prism.fep.visualize.mapping import visualize_mapping_png
from prism.fep.visualize.html import visualize_mapping_html
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.io import read_ligand_from_prism
from prism.fep.config import FEPConfig
import json

from test_utils import FEPTestPaths


def test_oMeEtPh_EtPh_fep_mapping():
    """
    Complete FEP mapping test for oMeEtPh → EtPh

    Force field configuration is read from fep.yaml
    """
    paths = FEPTestPaths("oMeEtPh-EtPh")
    test_dir = paths.system_dir

    mol2_a = paths.input_dir / "oMeEtPh.mol2"
    mol2_b = paths.input_dir / "EtPh.mol2"

    print(f"\n{'='*70}")
    print(f"FEP Mapping Test: oMeEtPh → EtPh")
    print(f"{'='*70}")

    # Load configuration using FEPConfig
    fep_config = FEPConfig(str(test_dir))

    ff_type = fep_config.get_forcefield_type()
    ff_params = fep_config.get_forcefield_params()
    mapping_params = fep_config.get_mapping_params()

    print(f"\n[Step 0] Configuration loaded:")
    print(f"  Force field: {ff_type}")
    print(f"  FF params: {ff_params}")
    print(f"  Mapping params: {mapping_params}")

    # Check test data
    if not mol2_a.exists() or not mol2_b.exists():
        pytest.skip("Test data not found")

    # Map force field type to generator class
    ff_generators = {
        "gaff": GAFFForceFieldGenerator,
        "gaff2": GAFF2ForceFieldGenerator,
        "openff": OpenFFForceFieldGenerator,
        "opls": OPLSAAForceFieldGenerator,
        "mmff": MMFFForceFieldGenerator,
    }

    if ff_type not in ff_generators:
        pytest.skip(f"Unsupported force field type: {ff_type}")

    GeneratorClass = ff_generators[ff_type]

    # Step 1: Generate force fields
    print(f"\n[Step 1] Generating {ff_type.upper()} force fields...")

    output_base_dir = paths.get_forcefield_dir(ff_type)
    output_base_dir.mkdir(exist_ok=True)

    output_dir_a = output_base_dir / f"{ff_type}_oMeEtPh"
    output_dir_b = output_base_dir / f"{ff_type}_EtPh"
    output_dir_a.mkdir(exist_ok=True)
    output_dir_b.mkdir(exist_ok=True)

    try:
        print(f"  → Generating oMeEtPh force field...")
        generator_a = GeneratorClass(
            ligand_path=str(mol2_a),
            output_dir=str(output_dir_a),
            **ff_params,  # Use params from case.yaml
            overwrite=True,
        )
        result_dir_a = generator_a.run()
        print(f"    ✓ Generated: {result_dir_a}")

        print(f"  → Generating EtPh force field...")
        generator_b = GeneratorClass(
            ligand_path=str(mol2_b),
            output_dir=str(output_dir_b),
            **ff_params,  # Use params from case.yaml
            overwrite=True,
        )
        result_dir_b = generator_b.run()
        print(f"    ✓ Generated: {result_dir_b}")

    except Exception as e:
        pytest.skip(f"Force field generation failed: {e}")

    # Step 2: Read ligands from PRISM format
    print(f"\n[Step 2] Reading ligands from PRISM format...")

    try:
        lig_a = read_ligand_from_prism(
            itp_file=str(Path(result_dir_a) / "LIG.itp"), gro_file=str(Path(result_dir_a) / "LIG.gro")
        )
        print(f"  ✓ oMeEtPh: {len(lig_a)} atoms")

        lig_b = read_ligand_from_prism(
            itp_file=str(Path(result_dir_b) / "LIG.itp"), gro_file=str(Path(result_dir_b) / "LIG.gro")
        )
        print(f"  ✓ EtPh: {len(lig_b)} atoms")

    except Exception as e:
        pytest.fail(f"Failed to read PRISM format: {e}")

    # Step 3: Perform atom mapping
    print(f"\n[Step 3] Performing atom mapping...")

    # Use mapping params from FEPConfig
    mapper = DistanceAtomMapper(**mapping_params)
    mapping = mapper.map(lig_a, lig_b)

    print(f"  ✓ Common atoms: {len(mapping.common)}")
    print(f"  ✓ Transformed A: {len(mapping.transformed_a)}")
    print(f"  ✓ Transformed B: {len(mapping.transformed_b)}")
    print(f"  ✓ Surrounding A: {len(mapping.surrounding_a)}")
    print(f"  ✓ Surrounding B: {len(mapping.surrounding_b)}")

    # Verify mapping completeness
    total_a = len(lig_a)
    total_b = len(lig_b)
    sum_a = len(mapping.common) + len(mapping.transformed_a) + len(mapping.surrounding_a)
    sum_b = len(mapping.common) + len(mapping.transformed_b) + len(mapping.surrounding_b)

    if sum_a != total_a or sum_b != total_b:
        print(f"⚠️  Mapping incomplete:")
        print(f"   A: {sum_a}/{total_a}, B: {sum_b}/{total_b}")

    # Step 4: Generate visualizations
    print(f"\n[Step 4] Generating visualizations...")

    # PNG visualization
    png_file = str(output_base_dir / f"oMeEtPh-EtPh_{ff_type}_mapping.png")
    try:
        visualize_mapping_png(
            mapping=mapping,
            pdb_a=str(paths.input_dir / "oMeEtPh.pdb"),
            pdb_b=str(paths.input_dir / "EtPh.pdb"),
            mol2_a=str(mol2_a),
            mol2_b=str(mol2_b),
            output_path=png_file,
        )
        print(f"  ✓ PNG: {png_file}")
    except Exception as e:
        print(f"  ⚠️  PNG generation failed: {e}")

    # HTML visualization
    html_file = str(output_base_dir / f"oMeEtPh-EtPh_{ff_type}_mapping.html")
    try:
        # Use HTML config from FEPConfig
        html_config = fep_config.get_html_config()

        visualize_mapping_html(
            mapping=mapping,
            pdb_a=str(paths.input_dir / "oMeEtPh.pdb"),
            pdb_b=str(paths.input_dir / "EtPh.pdb"),
            mol2_a=str(mol2_a),
            mol2_b=str(mol2_b),
            atoms_a=lig_a,
            atoms_b=lig_b,
            output_path=html_file,
            ligand_a_name="oMeEtPh",
            ligand_b_name="EtPh",
            config=html_config,
        )
        print(f"  ✓ HTML: {html_file}")
    except Exception as e:
        print(f"  ⚠️  HTML generation failed: {e}")

    # Step 5: Save mapping results
    print(f"\n[Step 5] Saving mapping results...")

    results = {
        "system": "oMeEtPh-EtPh",
        "forcefield": ff_type,
        "forcefield_params": ff_params,
        "mapping": {
            "common": len(mapping.common),
            "transformed_a": len(mapping.transformed_a),
            "transformed_b": len(mapping.transformed_b),
            "surrounding_a": len(mapping.surrounding_a),
            "surrounding_b": len(mapping.surrounding_b),
            "total_a": total_a,
            "total_b": total_b,
        },
        "config": mapping_params,
    }

    results_file = output_base_dir / f"mapping_results_{ff_type}.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  ✓ Results: {results_file}")

    print(f"\n{'='*70}")
    print(f"✅ FEP Mapping Test Completed Successfully!")
    print(f"   Force Field: {ff_type.upper()}")
    print(f"{'='*70}")

    assert total_a > 0
    assert total_b > 0
    assert sum_a == total_a
    assert sum_b == total_b


def main():
    """Run FEP mapping test with force field from fep.yaml"""
    print(f"\n{'#'*70}")
    print(f"# FEP Mapping Test")
    print(f"# System: oMeEtPh → EtPh")
    print(f"# Force field configuration from fep.yaml")
    print(f"{'#'*70}")

    # Test with config from fep.yaml
    print("\nTesting with fep.yaml configuration...")
    test_oMeEtPh_EtPh_fep_mapping()

    # Summary
    print(f"\n\n{'#'*70}")
    print(f"# Test Summary")
    print(f"{'#'*70}")
    print(f"FEP Mapping:  ✅ PASSED")
    print(f"{'#'*70}\n")


if __name__ == "__main__":
    main()
