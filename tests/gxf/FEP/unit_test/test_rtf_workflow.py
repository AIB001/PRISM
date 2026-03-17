#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
测试从 RTF+PRM 文件生成 PRISM 格式的工作流

这个测试展示了如何从 CGenFF 的 RTF+PRM 文件生成 PRISM 格式（ITP+GRO），
然后用于 FEP 原子映射和可视化。

注意：
- 输出目录使用 LIG.rtf2gmx/ 而不是 LIG.charmm2gmx/，避免覆盖现有文件
- 目前 PRISM 还没有 RTFForceFieldGenerator，这个测试展示了需求
"""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from prism.fep.io import read_rtf_for_fep, read_ligand_from_prism
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.visualize import visualize_mapping_png, visualize_mapping_html
from prism.utils.config import ConfigurationManager
from conftest import resolve_fep_case_dir


def test_rtf_to_prism_workflow_25_36():
    """
    测试 RTF+PRM → PRISM 格式的完整工作流（25-36 系统）

    工作流程：
    1. 方法 1（当前）：直接从 RTF+PDB 读取（快速测试）
    2. 方法 2（目标）：RTF+PRM → PRISM 格式 → FEP 映射

    注意：方法 2 需要实现 RTFForceFieldGenerator
    """

    test_dir = resolve_fep_case_dir("25-36")

    print("\n" + "=" * 70)
    print("RTF → PRISM 格式工作流测试：25-36 系统")
    print("=" * 70)

    # ============================================================
    # 方法 1: 直接从 RTF+PDB 读取（当前使用）
    # ============================================================
    print("\n[方法 1] 直接从 RTF+PDB 读取（快速测试）")

    rtf_25 = test_dir / "25.rtf"
    pdb_25 = test_dir / "25.pdb"
    rtf_36 = test_dir / "36.rtf"
    pdb_36 = test_dir / "36.pdb"

    if not all([rtf_25.exists(), pdb_25.exists(), rtf_36.exists(), pdb_36.exists()]):
        pytest.skip("RTF/PDB 文件不存在")

    lig_25 = read_rtf_for_fep(str(rtf_25), str(pdb_25))
    lig_36 = read_rtf_for_fep(str(rtf_36), str(pdb_36))

    print(f"  ✓ Ligand 25: {len(lig_25)} 个原子")
    print(f"  ✓ Ligand 36: {len(lig_36)} 个原子")

    # 执行映射
    config_file = test_dir / "case.yaml"
    cfg = ConfigurationManager(config_path=str(config_file))
    mapper = DistanceAtomMapper.from_config(cfg.config)
    mapping = mapper.map(lig_25, lig_36)

    print(f"  ✓ Common: {len(mapping.common)}")
    print(f"  ✓ Transformed 25: {len(mapping.transformed_a)}")
    print(f"  ✓ Transformed 36: {len(mapping.transformed_b)}")
    print(f"  ✓ Surrounding 25: {len(mapping.surrounding_a)}")
    print(f"  ✓ Surrounding 36: {len(mapping.surrounding_b)}")

    # 生成可视化
    output_dir = test_dir / "output"
    output_dir.mkdir(exist_ok=True)

    visualize_mapping_png(
        mapping=mapping,
        pdb_a=str(pdb_25),
        pdb_b=str(pdb_36),
        mol2_a=str(test_dir / "25_3D.mol2"),
        mol2_b=str(test_dir / "36_3D.mol2"),
        output_path=str(output_dir / "25-36_rtf_direct.png"),
    )

    visualize_mapping_html(
        mapping=mapping,
        pdb_a=str(pdb_25),
        pdb_b=str(pdb_36),
        mol2_a=str(test_dir / "25_3D.mol2"),
        mol2_b=str(test_dir / "36_3D.mol2"),
        atoms_a=lig_25,
        atoms_b=lig_36,
        output_path=str(output_dir / "25-36_rtf_direct.html"),
        title="RTF Direct Read: 25-36 System",
        ligand_a_name="Ligand 25",
        ligand_b_name="Ligand 36",
        config=cfg.config,
    )

    print("  ✓ PNG 可视化: 25-36_rtf_direct.png")
    print("  ✓ HTML 可视化: 25-36_rtf_direct.html")

    # ============================================================
    # 方法 2: RTF+PRM → PRISM 格式（目标，需要实现）
    # ============================================================
    print("\n[方法 2] RTF+PRM → PRISM 格式（需要实现 RTFForceFieldGenerator）")
    print("  目标流程：")
    print("    1. RTFForceFieldGenerator(rtf_file, prm_file, pdb_file)")
    print("       → 生成 LIG.rtf2gmx/LIG.itp + LIG.gro")
    print("    2. read_ligand_from_prism(itp_file, gro_file)")
    print("       → 读取 PRISM 格式")
    print("    3. DistanceAtomMapper.map(lig_a, lig_b)")
    print("       → 执行原子映射")
    print("    4. visualize_mapping_*(...)")
    print("       → 生成可视化")

    print("\n  优势：")
    print("    - 统一接口：与其他力场（GAFF/OpenFF/CHARMM-GUI）一致")
    print("    - 完整信息：包含键、角、二面角等拓扑信息")
    print("    - 可扩展：支持复杂的 CHARMM 力场参数")

    print("\n  实现建议：")
    print("    - 参考 CHARMMGUIForceFieldGenerator 的实现")
    print("    - 解析 RTF 文件获取原子类型、电荷、键、角等")
    print("    - 解析 PRM 文件获取力场参数")
    print("    - 生成标准的 GROMACS ITP 格式")
    print("    - 输出目录：LIG.rtf2gmx/ (避免与 charmm2gmx 冲突)")

    print("\n  文件对比：")
    print("    当前：25.rtf + 25.pdb → read_rtf_for_fep() → List[Atom]")
    print("    目标：25.rtf + 25.prm + 25.pdb → RTFForceFieldGenerator")
    print("          → LIG.rtf2gmx/LIG.itp + LIG.gro → read_ligand_from_prism()")

    print("\n" + "=" * 70)
    print("结论：")
    print("  - 方法 1 适合快速验证算法 ✓ 当前使用")
    print("  - 方法 2 是真正的 PRISM 集成 ⏳ 待实现")
    print("  - 不会覆盖现有的 LIG.charmm2gmx/ 目录")
    print("=" * 70 + "\n")


@pytest.mark.slow
def test_rtf_force_field_generator_25_36(tmp_path):
    """
    测试 RTFForceFieldGenerator

    这个测试展示了 RTFForceFieldGenerator 的实际使用
    """
    from prism.forcefield.rtf import RTFForceFieldGenerator

    test_dir = resolve_fep_case_dir("25-36")
    required = [
        test_dir / "25.rtf",
        test_dir / "25.prm",
        test_dir / "25.pdb",
        test_dir / "36.rtf",
        test_dir / "36.prm",
        test_dir / "36.pdb",
    ]
    if not all(path.exists() for path in required):
        pytest.skip("RTF/PRM/PDB files not found for 25-36")

    # Step 1: RTF+PRM → PRISM 格式
    generator_25 = RTFForceFieldGenerator(
        rtf_file=str(test_dir / "25.rtf"),
        prm_file=str(test_dir / "25.prm"),
        pdb_file=str(test_dir / "25.pdb"),
        output_dir=str(tmp_path / "lig25_output"),
    )
    lig_25_dir = generator_25.run()  # → LIG.rtf2gmx/

    generator_36 = RTFForceFieldGenerator(
        rtf_file=str(test_dir / "36.rtf"),
        prm_file=str(test_dir / "36.prm"),
        pdb_file=str(test_dir / "36.pdb"),
        output_dir=str(tmp_path / "lig36_output"),
    )
    lig_36_dir = generator_36.run()

    # Step 2: 读取 PRISM 格式
    lig_25 = read_ligand_from_prism(itp_file=f"{lig_25_dir}/LIG.itp", gro_file=f"{lig_25_dir}/LIG.gro")
    lig_36 = read_ligand_from_prism(itp_file=f"{lig_36_dir}/LIG.itp", gro_file=f"{lig_36_dir}/LIG.gro")

    # Step 3: 原子映射
    config_file = test_dir / "case.yaml"
    cfg = ConfigurationManager(config_path=str(config_file))
    mapper = DistanceAtomMapper.from_config(cfg.config)
    mapping = mapper.map(lig_25, lig_36)

    # Step 4: 可视化
    output_dir = tmp_path / "output"
    output_dir.mkdir(exist_ok=True)

    visualize_mapping_png(
        mapping=mapping,
        pdb_a=str(test_dir / "25.pdb"),
        pdb_b=str(test_dir / "36.pdb"),
        mol2_a=str(test_dir / "25_3D.mol2"),
        mol2_b=str(test_dir / "36_3D.mol2"),
        output_path=str(output_dir / "25-36_rtf_prism.png"),
    )

    print("✓ RTF → PRISM 格式工作流测试通过")


if __name__ == "__main__":
    test_rtf_to_prism_workflow_25_36()
