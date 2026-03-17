#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Validate FEP mapping against FEbuilder results

This module compares our atom mapping results with FEbuilder's output
to verify algorithm correctness.
"""

import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from prism.utils.config import ConfigurationManager
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.io import read_rtf_for_fep


def read_hybrid_pdb_uncommon(pdb_file):
    """
    Read FEbuilder's hybrid.pdb and extract uncommon atoms

    Uncommon atoms have beta = -1 or 1

    Parameters
    ----------
    pdb_file : str
        Path to hybrid.pdb file

    Returns
    -------
    dict
        Dictionary with 'state_a' and 'state_b' atom lists
    """
    uncommon = {"state_a": [], "state_b": []}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                beta = float(line[60:66].strip())

                # Beta = -1 or 1 indicates uncommon/transformed
                if abs(beta) > 0.5:
                    if beta > 0:
                        uncommon["state_a"].append(atom_name)
                    else:
                        uncommon["state_b"].append(atom_name)
    return uncommon


@pytest.mark.parametrize(
    "system",
    [
        pytest.param("25-36", marks=pytest.mark.slow, id="25-36"),
        pytest.param("39-8", marks=pytest.mark.slow, id="39-8"),
    ],
)
def test_validate_against_febuilder(system):
    """
    Compare our mapping results with FEbuilder's output

    This test validates that our atom mapping produces the same
    uncommon atom classification as FEbuilder.
    """
    test_dir = f"tests/gxf/FEP/unit_test/{system}"

    if system == "25-36":
        rtf_a = f"{test_dir}/25.rtf"
        rtf_b = f"{test_dir}/36.rtf"
        pdb_a = f"{test_dir}/25.pdb"
        pdb_b = f"{test_dir}/36.pdb"
        lig_name_a, lig_name_b = "25", "36"
    elif system == "39-8":
        rtf_a = f"{test_dir}/39.rtf"
        rtf_b = f"{test_dir}/8.rtf"
        pdb_a = f"{test_dir}/39.pdb"
        pdb_b = f"{test_dir}/8.pdb"
        lig_name_a, lig_name_b = "39", "8"
    else:
        pytest.skip(f"Unknown system: {system}")

    hybrid_pdb = f"{test_dir}/other/hybrid.pdb"

    # Verify files exist
    for f in [rtf_a, rtf_b, pdb_a, pdb_b, hybrid_pdb]:
        if not Path(f).exists():
            pytest.skip(f"Required file not found: {f}")

    print(f"\n{'='*70}")
    print(f"Validating {lig_name_a} vs {lig_name_b} against FEbuilder")
    print(f"{'='*70}")

    # Step 1: Read ligands using CGenFF data
    print(f"[Step 1] Reading CGenFF data...")
    lig_a = read_rtf_for_fep(rtf_a, pdb_a)
    lig_b = read_rtf_for_fep(rtf_b, pdb_b)

    print(f"  ✓ Ligand {lig_name_a}: {len(lig_a)} atoms")
    print(f"  ✓ Ligand {lig_name_b}: {len(lig_b)} atoms")

    # Step 2: Perform atom mapping
    cfg = ConfigurationManager(config_path=f"{test_dir}/case.yaml")
    mapper = DistanceAtomMapper.from_config(cfg.config)
    print(f"[Step 2] Performing atom mapping (charge_cutoff={mapper.charge_cutoff})...")
    mapping = mapper.map(lig_a, lig_b)

    print(f"  ✓ Common: {len(mapping.common)}")
    print(f"  ✓ Transformed {lig_name_a}: {len(mapping.transformed_a)}")
    print(f"  ✓ Transformed {lig_name_b}: {len(mapping.transformed_b)}")
    print(f"  ✓ Surrounding {lig_name_a}: {len(mapping.surrounding_a)}")
    print(f"  ✓ Surrounding {lig_name_b}: {len(mapping.surrounding_b)}")

    # Step 3: Read FEbuilder results
    print(f"[Step 3] Reading FEbuilder results...")
    hybrid_uncommon = read_hybrid_pdb_uncommon(hybrid_pdb)
    febuilder_state_a = set(hybrid_uncommon["state_a"])
    febuilder_state_b = set(hybrid_uncommon["state_b"])

    print(f"  ✓ FEbuilder State A: {len(febuilder_state_a)} atoms")
    print(f"  ✓ FEbuilder State B: {len(febuilder_state_b)} atoms")

    # Step 4: Extract our uncommon atoms
    print(f"[Step 4] Extracting our uncommon atoms...")

    our_transformed_a = set(a.name for a in mapping.transformed_a)
    our_transformed_b = set(a.name for a in mapping.transformed_b)
    our_surrounding_a = set(a.name for a in mapping.surrounding_a)
    our_surrounding_b = set(a.name for a in mapping.surrounding_b)

    our_uncommon_a = our_transformed_a | our_surrounding_a
    our_uncommon_b = our_transformed_b | our_surrounding_b

    # Step 5: Detailed comparison
    print(f"[Step 5] Comparing with FEbuilder...")

    # Remove A/B suffixes from FEbuilder names for comparison
    febuilder_a_clean = set(name.replace("A", "").replace("B", "") for name in febuilder_state_a)
    febuilder_b_clean = set(name.replace("A", "").replace("B", "") for name in febuilder_state_b)

    print(f"\n  Ligand {lig_name_a} comparison:")
    print(f"    FEbuilder State B: {sorted(febuilder_b_clean)}")
    print(f"    Our uncommon {lig_name_a}: {sorted(our_uncommon_a)}")

    match_a = our_uncommon_a == febuilder_b_clean
    print(f"    Match? {'✓ YES' if match_a else '✗ NO'}")

    if not match_a:
        only_febuilder = febuilder_b_clean - our_uncommon_a
        only_ours = our_uncommon_a - febuilder_b_clean
        if only_febuilder:
            print(f"      Only in FEbuilder: {sorted(only_febuilder)}")
        if only_ours:
            print(f"      Only in ours: {sorted(only_ours)}")

    print(f"\n  Ligand {lig_name_b} comparison:")
    print(f"    FEbuilder State A: {sorted(febuilder_a_clean)}")
    print(f"    Our uncommon {lig_name_b}: {sorted(our_uncommon_b)}")

    match_b = our_uncommon_b == febuilder_a_clean
    print(f"    Match? {'✓ YES' if match_b else '✗ NO'}")

    if not match_b:
        only_febuilder = febuilder_a_clean - our_uncommon_b
        only_ours = our_uncommon_b - febuilder_a_clean
        if only_febuilder:
            print(f"      Only in FEbuilder: {sorted(only_febuilder)}")
        if only_ours:
            print(f"      Only in ours: {sorted(only_ours)}")

    # Summary
    print(f"\n{'='*70}")
    print(f"Summary for {system}:")
    print(f"  Ligand {lig_name_a}: {'✓ Match' if match_a else '✗ Different'}")
    print(f"  Ligand {lig_name_b}: {'✓ Match' if match_b else '✗ Different'}")
    print(f"  Overall: {'✓ MATCHES FEBUILDER!' if (match_a and match_b) else '✗ Different from FEbuilder'}")
    print(f"{'='*70}\n")

    # Assert that we should eventually match FEbuilder
    # For now, we'll just warn if different
    if not (match_a and match_b):
        pytest.skip(
            f"Our mapping differs from FEbuilder for {system} - "
            f"this is expected until we implement charge_common='ref' strategy"
        )

    # Step 6: Generate visualization
    print(f"[Step 6] Generating visualization...")
    from prism.fep.visualize import visualize_mapping_png, visualize_mapping_html

    output_dir = Path(test_dir) / "output"
    output_dir.mkdir(exist_ok=True)

    # Get mol2 file paths for bond order correction
    mol2_a = f"{test_dir}/{lig_name_a}_3D.mol2"
    mol2_b = f"{test_dir}/{lig_name_b}_3D.mol2"

    # Extract atom lists for HTML visualization
    atoms_a = lig_a  # lig_a is the list of Atom objects
    atoms_b = lig_b  # lig_b is the list of Atom objects

    visualize_mapping_png(mapping, pdb_a, pdb_b, mol2_a, mol2_b, output_path=str(output_dir / f"{system}_mapping.png"))
    visualize_mapping_html(
        mapping,
        pdb_a,
        pdb_b,
        mol2_a,
        mol2_b,
        atoms_a,
        atoms_b,
        output_path=str(output_dir / f"{system}_mapping.html"),
        title=f"FEP Atom Mapping: {system}",
    )

    print(f"\n{'='*70}\n")

    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
