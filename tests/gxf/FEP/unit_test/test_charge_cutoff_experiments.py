#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test different charge_cutoff values to match FEbuilder results

This experiment tests whether adjusting charge_cutoff can make
our mapping match FEbuilder without implementing charge_common strategy.
"""

import pytest
import sys
from pathlib import Path
from conftest import resolve_fep_case_dir

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.io import read_rtf_for_fep


def read_hybrid_pdb_uncommon(pdb_file):
    """Read FEbuilder's hybrid.pdb and extract uncommon atoms"""
    uncommon = {"state_a": [], "state_b": []}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                beta = float(line[60:66].strip())
                if abs(beta) > 0.5:
                    if beta > 0:
                        uncommon["state_a"].append(atom_name)
                    else:
                        uncommon["state_b"].append(atom_name)
    return uncommon


@pytest.mark.parametrize("charge_cutoff", [0.05, 0.1, 0.15, 0.2, 0.25])
def test_25_36_with_different_charge_cutoffs(charge_cutoff):
    """
    Test 25-36 system with different charge_cutoff values

    Goal: Find if increasing charge_cutoff can match FEbuilder without
    implementing charge_common='ref' strategy.
    """
    test_dir = resolve_fep_case_dir("25-36")

    # Read case configuration
    case_config_file = test_dir / "case.yaml"
    if not case_config_file.exists():
        pytest.skip(f"Case config not found: {case_config_file}")

    # File paths
    rtf_25 = test_dir / "25.rtf"
    rtf_36 = test_dir / "36.rtf"
    pdb_25 = test_dir / "25.pdb"
    pdb_36 = test_dir / "36.pdb"
    hybrid_pdb = test_dir / "other/hybrid.pdb"

    if not all(Path(f).exists() for f in [rtf_25, rtf_36, pdb_25, pdb_36, hybrid_pdb]):
        pytest.skip("Required files not found")

    # Read ligands
    lig25 = read_rtf_for_fep(str(rtf_25), str(pdb_25))
    lig36 = read_rtf_for_fep(str(rtf_36), str(pdb_36))

    # Perform mapping with specific charge_cutoff
    mapper = DistanceAtomMapper(dist_cutoff=0.6, charge_cutoff=charge_cutoff)
    mapping = mapper.map(lig25, lig36)

    # Read FEbuilder results
    hybrid_uncommon = read_hybrid_pdb_uncommon(hybrid_pdb)
    febuilder_state_b = set(hybrid_uncommon["state_b"])

    # Extract our uncommon atoms
    our_uncommon_25 = set(a.name for a in mapping.transformed_a) | set(a.name for a in mapping.surrounding_a)

    # Remove A/B suffixes from FEbuilder names
    febuilder_b_clean = set(name.replace("A", "").replace("B", "") for name in febuilder_state_b)

    # Check if matches
    matches = our_uncommon_25 == febuilder_b_clean

    print(f"\n{'='*70}")
    print(f"Testing 25-36 with charge_cutoff={charge_cutoff}")
    print(f"{'='*70}")
    print(f"  Common atoms: {len(mapping.common)}")
    print(f"  Our uncommon 25: {sorted(our_uncommon_25)}")
    print(f"  FEbuilder uncommon: {sorted(febuilder_b_clean)}")
    print(f"  Match? {'✓ YES' if matches else '✗ NO'}")

    if not matches:
        only_febuilder = febuilder_b_clean - our_uncommon_25
        only_ours = our_uncommon_25 - febuilder_b_clean
        if only_febuilder:
            print(f"    Only in FEbuilder: {sorted(only_febuilder)}")
        if only_ours:
            print(f"    Only in ours: {sorted(only_ours)}")
    print(f"{'='*70}")

    # For charge_cutoff=0.2, we expect to still have differences
    # because the issue is algorithm logic, not just the threshold
    if charge_cutoff >= 0.2:
        # Even with larger cutoff, C6/C7/C8 might still be classified wrong
        # because they need charge_common='ref' strategy
        print(f"\n  Note: Even with charge_cutoff={charge_cutoff}, " f"C6/C7/C8 may still be wrong")
        print(f"  This suggests we need charge_common='ref' strategy, " f"not just larger cutoff")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
