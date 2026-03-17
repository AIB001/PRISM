#!/usr/bin/env python3
"""
批量测试不同 charge_common 模式并生成 HTML
"""

from pathlib import Path
from prism.fep.config import FEPConfig
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.io import read_ligand_from_prism
from prism.fep.visualize.html import visualize_mapping_html
import shutil

test_dir = Path("/data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/oMeEtPh-EtPh")

# Test different charge_common modes
modes = ["ref", "mut", "mean", "none"]

for mode in modes:
    print(f"\n{'='*70}")
    print(f"Testing charge_common={mode}")
    print(f"{'='*70}")

    # Copy corresponding fep.yaml
    fep_file = test_dir / f"fep_{mode}.yaml"
    if not fep_file.exists():
        print(f"⚠️  {fep_file} not found, skipping")
        continue

    # Copy to fep.yaml (will be read by FEPConfig)
    shutil.copy(fep_file, test_dir / "fep.yaml")

    # Load config
    fep_config = FEPConfig(str(test_dir))

    # Read ligands
    lig_a = read_ligand_from_prism(
        str(test_dir / "gaff_test_output/gaff_oMeEtPh/LIG.amb2gmx/LIG.itp"),
        str(test_dir / "gaff_test_output/gaff_oMeEtPh/LIG.amb2gmx/LIG.gro"),
    )
    lig_b = read_ligand_from_prism(
        str(test_dir / "gaff_test_output/gaff_EtPh/LIG.amb2gmx/LIG.itp"),
        str(test_dir / "gaff_test_output/gaff_EtPh/LIG.amb2gmx/LIG.gro"),
    )

    # Mapping
    mapper = DistanceAtomMapper(**fep_config.get_mapping_params())
    mapping = mapper.map(lig_a, lig_b)

    # Note: mapping.map() now modifies lig_a and lig_b in-place
    # So lig_a and lig_b already have the correct charges

    # Generate HTML
    output_file = test_dir / f"gaff_test_output/oMeEtPh-EtPh_gaff_{mode}.html"
    visualize_mapping_html(
        mapping=mapping,
        pdb_a=str(test_dir / "oMeEtPh.pdb"),
        pdb_b=str(test_dir / "EtPh.pdb"),
        mol2_a=str(test_dir / "oMeEtPh.mol2"),
        mol2_b=str(test_dir / "EtPh.mol2"),
        atoms_a=lig_a,
        atoms_b=lig_b,
        output_path=str(output_file),
        ligand_a_name="oMeEtPh",
        ligand_b_name="EtPh",
        config=fep_config.get_html_config(),
    )
    print(f"✓ Generated: {output_file.name}")

print(f"\n{'='*70}")
print("All charge_common modes tested!")
print(f"{'='*70}")
