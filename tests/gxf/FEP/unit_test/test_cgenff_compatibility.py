#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
测试 CGenFFForceFieldGenerator 对 CHARMM-GUI 文件的兼容性

检查现有测试系统中的文件格式，验证是否需要扩展 CGenFFForceFieldGenerator
"""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from conftest import resolve_fep_case_dir


def test_charmm_gui_file_format_detection():
    """
    检测 CHARMM-GUI 生成的文件格式
    """
    print("\n" + "=" * 70)
    print("CHARMM-GUI 文件格式检测")
    print("=" * 70)

    test_cases = ["25-36", "39-8"]

    for case in test_cases:
        case_dir = resolve_fep_case_dir(case)
        print(f"\n测试案例: {case}")
        print(f"目录: {case_dir}")

        # 检查配体目录
        for ligand_name in ["25", "36", "39", "8"]:
            ligand_dir = case_dir / ligand_name
            if not ligand_dir.exists():
                continue

            print(f"\n  配体 {ligand_name}:")
            print(f"    目录: {ligand_dir}")

            # 列出所有文件
            files = list(ligand_dir.glob("*"))
            for f in sorted(files):
                if f.is_file():
                    print(f"      - {f.name} ({f.stat().st_size} bytes)")

            # 检查是否有 CGenFF 网站格式
            gmx_pdb = list(ligand_dir.glob("*_gmx.pdb"))
            gmx_top = list(ligand_dir.glob("*_gmx.top"))

            # 检查是否有 CHARMM-GUI 格式
            rtf = list(ligand_dir.glob("*.rtf"))
            prm = list(ligand_dir.glob("*.prm"))
            pdb = list(ligand_dir.glob("*.pdb"))

            if gmx_pdb and gmx_top:
                print(f"    ✓ CGenFF 网站格式: *_gmx.pdb + *_gmx.top")
            elif rtf and prm and pdb:
                print(f"    ✓ CHARMM-GUI 格式: .rtf + .prm + .pdb")
                print(f"      RTF: {rtf[0].name}")
                print(f"      PRM: {prm[0].name}")
                print(f"      PDB: {pdb[0].name}")
            else:
                print(f"    ⚠ 未知格式")

    print(f"\n{'='*70}\n")


@pytest.mark.slow
def test_cgenff_generator_with_charmm_gui_files():
    """
    测试 CGenFFForceFieldGenerator 能否处理 CHARMM-GUI 文件

    当前 CGenFFForceFieldGenerator 期望：
    - *_gmx.pdb (CGenFF 网站格式)
    - *_gmx.top (GROMACS topology 格式)

    CHARMM-GUI 提供：
    - *.pdb (CHARMM PDB 格式)
    - *.rtf (CHARMM topology 格式)
    - *.prm (CHARMM parameter 格式)

    这个测试检查兼容性。
    """

    print("\n" + "=" * 70)
    print("CGenFFForceFieldGenerator 兼容性测试")
    print("=" * 70)

    test_dir = resolve_fep_case_dir("25-36")
    ligand_25_dir = test_dir / "25"

    print(f"\n测试目录: {ligand_25_dir}")

    # 检查现有文件
    rtf_files = list(ligand_25_dir.glob("*.rtf"))
    prm_files = list(ligand_25_dir.glob("*.prm"))
    pdb_files = list(ligand_25_dir.glob("*.pdb"))

    print(f"\n现有文件:")
    print(f"  RTF 文件: {len(rtf_files)}")
    if rtf_files:
        print(f"    - {rtf_files[0].name}")
    print(f"  PRM 文件: {len(prm_files)}")
    if prm_files:
        print(f"    - {prm_files[0].name}")
    print(f"  PDB 文件: {len(pdb_files)}")
    if pdb_files:
        for f in pdb_files:
            print(f"    - {f.name}")

    # 检查是否有 CGenFF 网站格式
    gmx_pdb = list(ligand_25_dir.glob("*_gmx.pdb"))
    gmx_top = list(ligand_25_dir.glob("*_gmx.top"))

    print(f"\nCGenFF 网站格式:")
    print(f"  *_gmx.pdb: {len(gmx_pdb)}")
    print(f"  *_gmx.top: {len(gmx_top)}")

    if not (gmx_pdb and gmx_top):
        print(f"\n❌ 缺少 CGenFF 网站格式文件")
        print(f"\n当前状态:")
        print(f"  - 现有文件是 CHARMM-GUI 格式（RTF+PRM+PDB）")
        print(f"  - CGenFFForceFieldGenerator 期望 CGenFF 网站格式")
        print(f"    (*_gmx.pdb + *_gmx.top)")
        print(f"\n需要:")
        print(f"  1. 从 CGenFF 网站重新下载文件")
        print(f"     访问: https://cgenff.paramchem.org/")
        print(f"  2. 或者扩展 CGenFFForceFieldGenerator 支持 RTF+PRM")

        pytest.skip("需要 CGenFF 网站格式的 *_gmx.pdb 和 *_gmx.top 文件")

    print(f"\n{'='*70}\n")


def test_rtf_file_structure():
    """
    分析 RTF 文件的结构
    """
    print("\n" + "=" * 70)
    print("RTF 文件结构分析")
    print("=" * 70)

    test_dir = resolve_fep_case_dir("25-36")
    rtf_25 = test_dir / "25.rtf"

    if not rtf_25.exists():
        pytest.skip(f"RTF 文件不存在: {rtf_25}")

    print(f"\n文件: {rtf_25}")

    with open(rtf_25, "r") as f:
        lines = f.readlines()

    print(f"\n文件头:")
    for i, line in enumerate(lines[:10], 1):
        print(f"  {i}: {line.rstrip()}")

    print(f"\n文件格式特征:")

    # 检查关键字
    has_resi = any("RESI" in line for line in lines)
    has_atom = any("ATOM" in line for line in lines)
    has_mass = any("MASS" in line for line in lines)
    has_bond = any("BOND" in line for line in lines)
    has_angle = any("ANGLE" in line for line in lines)
    has_dihe = any("DIHEDRAL" in line for line in lines)

    print(f"  RESI (残基定义): {'✓' if has_resi else '✗'}")
    print(f"  ATOM (原子定义): {'✓' if has_atom else '✗'}")
    print(f"  MASS (质量): {'✓' if has_mass else '✗'}")
    print(f"  BOND (键): {'✓' if has_bond else '✗'}")
    print(f"  ANGLE (角): {'✓' if has_angle else '✗'}")
    print(f"  DIHEDRAL (二面角): {'✓' if has_dihe else '✗'}")

    print(f"\n这是标准的 CHARMM RTF 格式")
    print(f"与 GROMACS TOP 格式的区别:")
    print(f"  - CHARMM RTF: RESI, ATOM, MASS, BOND, ANGLE, DIHEDRAL")
    print(f"  - GROMACS TOP: [ moleculetype ], [ atoms ], [ bonds ], ...")

    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    print("\n运行 CHARMM-GUI 文件格式检测...\n")
    test_charmm_gui_file_format_detection()
    test_rtf_file_structure()
    print("\n检测完成！")
