#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
测试 CHARMMGUIForceFieldGenerator

验证 CHARMM-GUI 生成的文件能否转换为 PRISM 标准格式
"""

import pytest
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from prism.forcefield.charmm_gui import CHARMMGUIForceFieldGenerator
from prism.fep.io import read_ligand_from_prism
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.visualize import visualize_mapping_png, visualize_mapping_html
from prism.utils.config import ConfigurationManager
from conftest import resolve_fep_case_dir


@pytest.mark.slow
def test_charmm_gui_to_prism_format_39_8():
    """
    测试 CHARMM-GUI → PRISM 格式转换（39-8 系统）

    流程：
    1. CHARMM-GUI 输出 (gromacs/LIG.itp)
    2. CHARMMGUIForceFieldGenerator 转换
    3. 生成 PRISM 格式 (LIG.charmm2gmx/LIG.itp + LIG.gro)
    4. FEP 模块读取 PRISM 格式
    5. 执行原子映射
    """

    print("\n" + "=" * 70)
    print("CHARMM-GUI → PRISM 格式转换测试：39-8 系统")
    print("=" * 70)

    test_dir = Path("tests/gxf/FEP/unit_test/39-8")

    # ============================================================
    # Step 1: 转换 Ligand 39
    # ============================================================
    print("\n[Step 1] 转换 Ligand 39: CHARMM-GUI → PRISM 格式")

    ligand_39_dir = test_dir / "39"
    output_39 = test_dir / "39"  # Will create 39/LIG.charmm2gmx/

    generator_39 = CHARMMGUIForceFieldGenerator(
        ligand_path=str(ligand_39_dir / "39.pdb"),
        output_dir=str(output_39),
        charmm_gui_dir=str(ligand_39_dir),
        overwrite=True,
    )

    prism_dir_39 = generator_39.run()
    itp_39 = Path(prism_dir_39) / "LIG.itp"
    gro_39 = Path(prism_dir_39) / "LIG.gro"

    assert itp_39.exists(), f"ITP 文件未生成: {itp_39}"
    assert gro_39.exists(), f"GRO 文件未生成: {gro_39}"
    print(f"  ✓ Ligand 39 PRISM 格式生成成功: {prism_dir_39}")

    # ============================================================
    # Step 2: 转换 Ligand 8
    # ============================================================
    print("\n[Step 2] 转换 Ligand 8: CHARMM-GUI → PRISM 格式")

    ligand_8_dir = test_dir / "8"
    output_8 = test_dir / "8"  # Will create 8/LIG.charmm2gmx/

    generator_8 = CHARMMGUIForceFieldGenerator(
        ligand_path=str(ligand_8_dir / "8.pdb"),
        output_dir=str(output_8),
        charmm_gui_dir=str(ligand_8_dir),
        overwrite=True,
    )

    prism_dir_8 = generator_8.run()
    itp_8 = Path(prism_dir_8) / "LIG.itp"
    gro_8 = Path(prism_dir_8) / "LIG.gro"

    assert itp_8.exists(), f"ITP 文件未生成: {itp_8}"
    assert gro_8.exists(), f"GRO 文件未生成: {gro_8}"
    print(f"  ✓ Ligand 8 PRISM 格式生成成功: {prism_dir_8}")

    # ============================================================
    # Step 3: FEP 模块读取 PRISM 格式文件
    # ============================================================
    print("\n[Step 3] FEP 模块读取 PRISM 格式文件")

    lig_39 = read_ligand_from_prism(itp_file=str(itp_39), gro_file=str(gro_39))
    print(f"  ✓ Ligand 39: {len(lig_39)} 个原子")

    lig_8 = read_ligand_from_prism(itp_file=str(itp_8), gro_file=str(gro_8))
    print(f"  ✓ Ligand 8: {len(lig_8)} 个原子")

    # ============================================================
    # Step 4: 执行原子映射
    # ============================================================
    print("\n[Step 4] 执行原子映射")

    config_file = test_dir / "case.yaml"
    cfg = ConfigurationManager(config_path=str(config_file))
    mapper = DistanceAtomMapper.from_config(cfg.config)
    mapping = mapper.map(lig_39, lig_8)

    print(f"  ✓ Common: {len(mapping.common)}")
    print(f"  ✓ Transformed 39: {len(mapping.transformed_a)}")
    print(f"  ✓ Transformed 8: {len(mapping.transformed_b)}")
    print(f"  ✓ Surrounding 39: {len(mapping.surrounding_a)}")
    print(f"  ✓ Surrounding 8: {len(mapping.surrounding_b)}")

    # ============================================================
    # Step 5: 生成可视化
    # ============================================================
    print("\n[Step 5] 生成可视化")

    output_dir = test_dir / "output"
    output_dir.mkdir(exist_ok=True)

    # Prefer drawing_3D.mol2 from ligand subdir (CHARMM-GUI output), fallback to *_3D.mol2
    mol2_39 = ligand_39_dir / "drawing_3D.mol2"
    if not mol2_39.exists():
        mol2_39 = test_dir / "39_3D.mol2"
    mol2_8 = ligand_8_dir / "drawing_3D.mol2"
    if not mol2_8.exists():
        mol2_8 = test_dir / "8_3D.mol2"

    pdb_39 = str(ligand_39_dir / "39.pdb")
    pdb_8 = str(ligand_8_dir / "8.pdb")

    png_output = output_dir / "39-8_charmm_gui_prism.png"
    visualize_mapping_png(
        mapping=mapping,
        pdb_a=pdb_39,
        pdb_b=pdb_8,
        mol2_a=str(mol2_39) if mol2_39.exists() else None,
        mol2_b=str(mol2_8) if mol2_8.exists() else None,
        output_path=str(png_output),
    )
    print(f"  ✓ PNG 可视化: {png_output}")

    html_output = output_dir / "39-8_charmm_gui_prism.html"
    visualize_mapping_html(
        mapping=mapping,
        pdb_a=pdb_39,
        pdb_b=pdb_8,
        mol2_a=str(mol2_39) if mol2_39.exists() else None,
        mol2_b=str(mol2_8) if mol2_8.exists() else None,
        atoms_a=lig_39,
        atoms_b=lig_8,
        output_path=str(html_output),
        title="CHARMM-GUI to PRISM to FEP: 39-8 System",
        ligand_a_name="Ligand 39",
        ligand_b_name="Ligand 8",
        config=cfg.config,  # Pass configuration for display
    )
    print(f"  ✓ HTML 可视化: {html_output}")

    # ============================================================
    # Step 6: 验证结果
    # ============================================================
    print("\n[Step 6] 验证结果")

    assert len(lig_39) > 0, "Ligand 39 应该有原子"
    assert len(lig_8) > 0, "Ligand 8 应该有原子"
    assert len(mapping.common) >= 25, f"Common 原子数应该 >= 25"

    print(f"\n{'='*70}")
    print("✓ CHARMM-GUI → PRISM 转换测试通过！")
    print(f"{'='*70}")
    print("\n文件流：")
    print(f"  CHARMM-GUI 输出 (gromacs/)")
    print(f"    ↓ CHARMMGUIForceFieldGenerator")
    print(f"  PRISM 格式 (LIG.charmm2gmx/)")
    print(f"    ↓ read_ligand_from_prism()")
    print(f"  Atom 列表")
    print(f"    ↓ DistanceAtomMapper")
    print(f"  AtomMapping")
    print(f"    ↓ visualize_mapping_*()")
    print(f"  可视化 (PNG+HTML)")
    print(f"{'='*70}\n")


def test_charmm_gui_file_structure():
    """
    测试 CHARMM-GUI 文件结构检测
    """
    print("\n" + "=" * 70)
    print("CHARMM-GUI 文件结构检测")
    print("=" * 70)

    test_cases = ["39-8", "25-36"]

    for case in test_cases:
        case_dir = resolve_fep_case_dir(case)
        print(f"\n测试案例: {case}")

        for ligand_name in ["39", "8", "25", "36"]:
            ligand_dir = case_dir / ligand_name
            if not ligand_dir.exists():
                continue

            print(f"\n  配体 {ligand_name}:")

            # 检查 gromacs 目录
            gromacs_dir = ligand_dir / "gromacs"
            if gromacs_dir.exists():
                print(f"    ✓ gromacs/ 目录存在")

                files = list(gromacs_dir.glob("*"))
                for f in sorted(files):
                    if f.is_file():
                        size = f.stat().st_size
                        print(f"      - {f.name} ({size} bytes)")

                # 检查必需文件
                lig_itp = gromacs_dir / "LIG.itp"
                if lig_itp.exists():
                    print(f"    ✓ LIG.itp 存在（CHARMM-GUI 拓扑）")
                else:
                    print(f"    ✗ LIG.itp 不存在")
            else:
                print(f"    ✗ gromacs/ 目录不存在")

            # 检查其他 CHARMM-GUI 文件
            rtf = ligand_dir / "lig.rtf"
            prm = ligand_dir / "lig.prm"

            if rtf.exists() or prm.exists():
                print(f"    ℹ 发现 CHARMM 文件:")
                if rtf.exists():
                    print(f"      - lig.rtf ({rtf.stat().st_size} bytes)")
                if prm.exists():
                    print(f"      - lig.prm ({prm.stat().st_size} bytes)")

    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    print("\n运行 CHARMM-GUI 转换测试...\n")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        test_charmm_gui_file_structure()

        try:
            test_charmm_gui_to_prism_format_39_8(tmp_path)
        except Exception as e:
            print(f"\n❌ 测试失败: {e}")
            import traceback

            traceback.print_exc()

    print("\n测试完成！")
