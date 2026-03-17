#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
End-to-end test: PRISM CGenFF → Atom Mapping

Uses CGenFF RTF files (atom types and charges) + PDB files (aligned coordinates).
Mapping parameters are read from case.yaml via ConfigurationManager.
"""

import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from prism.utils.config import ConfigurationManager
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.io import read_rtf_for_fep
from conftest import resolve_fep_case_dir


# Expected common atom counts per system (ground truth from FEbuilder)
EXPECTED_COMMON = {
    # 25-36 now resolves legacy `charge_reception=pert` through the current
    # `surround` behavior, giving 29 common atoms on the shipped hif2a data.
    "25-36": 29,
    "39-8": 31,
}


@pytest.mark.parametrize(
    "system",
    [
        pytest.param("25-36", marks=pytest.mark.slow, id="25-36"),
        pytest.param("39-8", marks=pytest.mark.slow, id="39-8"),
    ],
)
def test_cgenff_mapping_with_different_systems(system):
    """
    Test CGenFF mapping for different ligand pairs.

    Mapping parameters are read from case.yaml via ConfigurationManager.
    Run with: pytest -m slow -s
    """
    test_dir = resolve_fep_case_dir(system)

    if system == "25-36":
        rtf_a, pdb_a = test_dir / "25.rtf", test_dir / "25.pdb"
        rtf_b, pdb_b = test_dir / "36.rtf", test_dir / "36.pdb"
        lig_name_a, lig_name_b = "25", "36"
    elif system == "39-8":
        rtf_a, pdb_a = test_dir / "39.rtf", test_dir / "39.pdb"
        rtf_b, pdb_b = test_dir / "8.rtf", test_dir / "8.pdb"
        lig_name_a, lig_name_b = "39", "8"
    else:
        pytest.skip(f"Unknown system: {system}")

    if not all(Path(f).exists() for f in [rtf_a, rtf_b, pdb_a, pdb_b]):
        pytest.skip(f"Required files not found for {system}")

    cfg = ConfigurationManager(config_path=str(test_dir / "case.yaml"))
    mapper = DistanceAtomMapper.from_config(cfg.config)
    expected_common = EXPECTED_COMMON[system]

    print(f"\n{'='*70}")
    print(f"System: {lig_name_a} vs {lig_name_b}")
    print(f"Params: dist_cutoff={mapper.dist_cutoff}, charge_cutoff={mapper.charge_cutoff}")
    print(f"{'='*70}")

    lig_a = read_rtf_for_fep(str(rtf_a), str(pdb_a))
    lig_b = read_rtf_for_fep(str(rtf_b), str(pdb_b))
    print(f"  Ligand {lig_name_a}: {len(lig_a)} atoms")
    print(f"  Ligand {lig_name_b}: {len(lig_b)} atoms")

    mapping = mapper.map(lig_a, lig_b)

    print(f"\n  Common:           {len(mapping.common)}")
    print(f"  Transformed {lig_name_a}:  {[a.name for a in mapping.transformed_a]}")
    print(f"  Transformed {lig_name_b}:  {[a.name for a in mapping.transformed_b]}")
    print(f"  Surrounding {lig_name_a}: {[a.name for a in mapping.surrounding_a]}")
    print(f"  Surrounding {lig_name_b}: {[a.name for a in mapping.surrounding_b]}")

    assert len(mapping.common) == expected_common, f"Expected {expected_common} common atoms, got {len(mapping.common)}"
    print(f"\n  ✓ common atoms = {len(mapping.common)} (expected {expected_common})")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
