#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
42-38 系统 GAFF2 FEP 测试
验证 GAFF2 力场与 PRISM FEP 流程
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from prism.builder.core import PRISMBuilder


def test_gaff2_fep():
    """测试 GAFF2 FEP 流程"""

    # Clean up previous builds to ensure fresh build
    import shutil
    import os

    output_dir = "gaff2_fixed"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        print(f"  清理旧构建: {output_dir}")

    print("\n" + "=" * 70)
    print("42-38 系统 GAFF2 FEP 测试")
    print("=" * 70)

    # GAFF2 配体 + AMBER 蛋白力场
    builder = PRISMBuilder(
        protein_path="input/receptor.pdb",
        ligand_paths=["input/42.pdb"],
        output_dir=output_dir,
        ligand_forcefield="gaff2",
        config_path="configs/case_gaff2.yaml",
        forcefield="amber99sb",
        water_model="tip3p",
        overwrite=True,
        fep_mode=True,
        mutant_ligand="input/38.pdb",
        fep_config="configs/fep_gaff2.yaml",
    )

    print(f"\n[Step 1] 运行 PRISMBuilder.run()...")
    result = builder.run()
    print(f"  ✓ Run completed: {result}")

    # 检查输出
    fep_dir = Path(output_dir) / "GMX_PROLIG_FEP"
    mapping_html = fep_dir / "common" / "hybrid" / "mapping.html"

    print(f"\n[Step 2] 检查 Mapping HTML...")
    if mapping_html.exists():
        print(f"  ✓ Mapping HTML 已生成: {mapping_html}")

        # 检查 unknown atoms
        html_content = mapping_html.read_text()
        unknown_count = html_content.count('"classification": "unknown"')
        print(f"  Unknown atoms: {unknown_count}")

        # 提取 mapping 统计
        import re

        match = re.search(r"Ligand A: common=(\d+), transformed=(\d+), surrounding=(\d+)", html_content)
        if match:
            common, transformed, surrounding = match.groups()
            print(f"  Mapping: common={common}, transformed={transformed}, surrounding={surrounding}")
    else:
        print(f"  ❌ Mapping HTML 不存在")
        return False

    print(f"\n[Step 3] 检查 FEP scaffold...")
    bound_dir = fep_dir / "bound" / "repeat1"
    unbound_dir = fep_dir / "unbound" / "repeat1"

    for leg_name, leg_dir in [("Bound", bound_dir), ("Unbound", unbound_dir)]:
        if not leg_dir.exists():
            print(f"  ❌ {leg_name} 目录不存在")
            return False
        print(f"  ✓ {leg_name} 目录存在")

    print(f"\n[Step 4] 检查 B-state 原子...")
    hybrid_itp = fep_dir / "common" / "hybrid" / "hybrid.itp"
    if hybrid_itp.exists():
        content = hybrid_itp.read_text()
        in_atoms = False
        total_atoms = 0
        b_state_atoms = 0

        for line in content.split("\n"):
            if "[ atoms ]" in line:
                in_atoms = True
                continue
            if in_atoms and line.strip().startswith("["):
                break
            if in_atoms and line.strip() and not line.strip().startswith(";"):
                parts = line.split()
                if len(parts) >= 8:
                    total_atoms += 1
                    if len(parts) >= 11:
                        b_state_atoms += 1

        print(f"  总原子数: {total_atoms}")
        print(f"  B-state 原子数: {b_state_atoms}")

        if b_state_atoms == 0:
            print(f"  ❌ 没有 B-state 原子!")
            return False
        else:
            print(f"  ✅ 有 {b_state_atoms} 个 B-state 原子")
    else:
        print(f"  ❌ hybrid.itp 不存在")
        return False

    print("\n" + "=" * 70)
    print("✅ GAFF2 FEP 测试完成")
    print("=" * 70)
    print(f"\n输出目录: {fep_dir}")
    print(f"Mapping HTML: {mapping_html}")
    print(f"\n下一步：运行 grompp 验证")

    return True


if __name__ == "__main__":
    success = test_gaff2_fep()
    sys.exit(0 if success else 1)
