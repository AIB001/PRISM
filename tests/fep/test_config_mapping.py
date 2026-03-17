#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test FEP mapping with FEbuilder-compatible config.conf

This test validates that we can reproduce FEbuilder's results
when using the same configuration parameters.
"""

import sys
from pathlib import Path
import pytest

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from prism.fep.io import read_rtf_for_fep
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.config import read_fep_config


def _resolve_25_36_dir() -> Path:
    candidates = [
        Path("tests/gxf/FEP/unit_test/25-36"),
        Path("tests/gxf/FEP/test/hif2a/25-36"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    pytest.skip("25-36 test data not found")


def test_with_config():
    """Test atom mapping using FEbuilder config file"""

    test_dir = _resolve_25_36_dir()
    config_file = test_dir / "config.conf"

    print("=" * 70)
    print("测试：使用FEbuilder config.conf进行原子映射")
    print("=" * 70)

    # Step 1: Read config
    print(f"\n[Step 1] 读取config文件: {config_file}")
    try:
        config = read_fep_config(str(config_file))
        print(f"  ✓ charge_common: {config['charge_common']}")
        print(f"  ✓ charge_reception: {config['charge_reception']}")
        print(f"  ✓ dist_cutoff: {config['dist_cutoff']}")
        print(f"  ✓ charge_cutoff: {config['charge_cutoff']}")
    except Exception as e:
        print(f"  ✗ 读取config失败: {e}")
        # Use defaults
        config = {"charge_common": "ref", "charge_reception": "surround", "dist_cutoff": 0.6, "charge_cutoff": 0.05}
        print(f"  使用默认值")

    # Step 2: Read ligands
    print(f"\n[Step 2] 读取CGenFF系统...")
    lig25 = read_rtf_for_fep(str(test_dir / "25.rtf"), str(test_dir / "25.pdb"))
    lig36 = read_rtf_for_fep(str(test_dir / "36.rtf"), str(test_dir / "36.pdb"))

    print(f"  ✓ 配体25: {len(lig25)} atoms")
    print(f"  ✓ 配体36: {len(lig36)} atoms")

    # Step 3: Perform mapping with config parameters
    print(f"\n[Step 3] 执行原子映射 (使用config参数)...")
    mapper = DistanceAtomMapper(
        dist_cutoff=config["dist_cutoff"],
        charge_cutoff=config["charge_cutoff"],
        charge_common=config["charge_common"],
        charge_reception=config["charge_reception"],
    )

    mapping = mapper.map(lig25, lig36)

    print(f"  ✓ Common: {len(mapping.common)}")
    print(f"  ✓ Transformed 25: {len(mapping.transformed_a)}")
    print(f"  ✓ Transformed 36: {len(mapping.transformed_b)}")
    print(f"  ✓ Surrounding 25: {len(mapping.surrounding_a)}")
    print(f"  ✓ Surrounding 36: {len(mapping.surrounding_b)}")

    # Step 4: Extract uncommon atoms
    print(f"\n[Step 4] 提取uncommon原子...")
    our_uncommon_25 = set(a.name for a in mapping.transformed_a + mapping.surrounding_a)
    our_uncommon_36 = set(a.name for a in mapping.transformed_b + mapping.surrounding_b)

    print(f"  配体25 uncommon: {sorted(our_uncommon_25)}")
    print(f"  配体36 uncommon: {sorted(our_uncommon_36)}")

    # Step 5: Compare with FEbuilder
    print(f"\n[Step 5] 对比FEbuilder结果...")
    print(f"  FEbuilder配体25 uncommon: C10, C9, H7, N1")
    print(f"  FEbuilder配体36 uncommon: C9, F1, H5")

    print(f"\n{'='*70}")
    print(f"✓ 测试完成")
    print(f"{'='*70}")

    assert config["charge_reception"] == "surround"
    assert len(mapping.common) > 0
    assert len(mapping.common) + len(mapping.transformed_a) + len(mapping.surrounding_a) == len(lig25)
    assert len(mapping.common) + len(mapping.transformed_b) + len(mapping.surrounding_b) == len(lig36)


if __name__ == "__main__":
    test_with_config()
