#!/usr/bin/env python3
"""
Test charge consistency after applying charge_common strategies
"""

from pathlib import Path
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.io import read_ligand_from_prism

test_dir = Path("/data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/oMeEtPh-EtPh")

print("Testing charge_common strategies...")
print("=" * 70)

# Test ref/mut/mean modes - charges should be identical
for mode in ["ref", "mut", "mean"]:
    print(f"\nTesting mode: {mode}")

    # Read fresh ligands for each test
    lig_a = read_ligand_from_prism(
        str(test_dir / "gaff_test_output/gaff_oMeEtPh/LIG.amb2gmx/LIG.itp"),
        str(test_dir / "gaff_test_output/gaff_oMeEtPh/LIG.amb2gmx/LIG.gro"),
    )
    lig_b = read_ligand_from_prism(
        str(test_dir / "gaff_test_output/gaff_EtPh/LIG.amb2gmx/LIG.itp"),
        str(test_dir / "gaff_test_output/gaff_EtPh/LIG.amb2gmx/LIG.gro"),
    )

    mapper = DistanceAtomMapper(dist_cutoff=1.0, charge_cutoff=0.06, charge_common=mode)
    mapping = mapper.map(lig_a, lig_b)

    max_diff = 0.0
    for atom_a, atom_b in mapping.common:
        diff = abs(atom_a.charge - atom_b.charge)
        if diff > max_diff:
            max_diff = diff

    print(f"  Common atoms: {len(mapping.common)}")
    print(f"  Surrounding atoms: {len(mapping.surrounding_a)}")
    print(f"  Max charge difference: {max_diff:.10f}")

    if max_diff < 1e-10:
        print(f"  ✓ PASS: All common atoms have identical charges")
    else:
        print(f"  ✗ FAIL: Common atoms have different charges!")

# Test 'none' mode - charges can be different
print(f"\nTesting mode: none")

# Read fresh ligands
lig_a = read_ligand_from_prism(
    str(test_dir / "gaff_test_output/gaff_oMeEtPh/LIG.amb2gmx/LIG.itp"),
    str(test_dir / "gaff_test_output/gaff_oMeEtPh/LIG.amb2gmx/LIG.gro"),
)
lig_b = read_ligand_from_prism(
    str(test_dir / "gaff_test_output/gaff_EtPh/LIG.amb2gmx/LIG.itp"),
    str(test_dir / "gaff_test_output/gaff_EtPh/LIG.amb2gmx/LIG.gro"),
)

mapper = DistanceAtomMapper(dist_cutoff=1.0, charge_cutoff=0.06, charge_common="none")
mapping = mapper.map(lig_a, lig_b)

different_count = 0
for atom_a, atom_b in mapping.common:
    if abs(atom_a.charge - atom_b.charge) > 1e-6:
        different_count += 1

print(f"  Common atoms: {len(mapping.common)}")
print(f"  Surrounding atoms: {len(mapping.surrounding_a)}")
print(f"  Atoms with different charges: {different_count}")

# In 'none' mode, we expect:
# 1. Fewer common atoms (only those with identical charges)
# 2. More surrounding atoms (those with any charge difference)
# 3. Common atoms should have identical charges (since only identical ones remain)
if len(mapping.common) < 10 and len(mapping.surrounding_a) > 5:
    print(f"  ✓ PASS: 'none' mode correctly classifies atoms with charge differences as surrounding")
else:
    print(f"  ✗ FAIL: Expected fewer common atoms and more surrounding atoms in 'none' mode!")

print("\n" + "=" * 70)
print("All tests completed!")
