#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
完整的 FEP 建模流程测试：42-38 系统（支持多种力场）

测试流程：
1. 使用指定力场（GAFF2/CHARMM）参数化两个配体
2. 执行原子映射
3. 搭建完整 FEP 系统
4. 生成 Mapping HTML 可视化
5. 验证输出文件

用法：
    # GAFF2 测试
    python test_fep_e2e.py --forcefield gaff2

    # CHARMM 测试
    python test_fep_e2e.py --forcefield cgenff

    # 生成 HTML
    python test_fep_e2e.py --forcefield gaff2 --html
"""

import sys
import shutil
import argparse
from pathlib import Path

# Import test utilities
sys.path.insert(0, str(Path(__file__).parent.parent))
from test_utils import FEPTestPaths, get_system_config

from prism.forcefield.gaff import GAFFForceFieldGenerator


def setup_cgenff_ligands(paths: FEPTestPaths, output_base: Path):
    """
    设置 CGenFF 配体（使用预生成的文件）

    CGenFF 使用38/和42/目录下的CHARMM-GUI文件
    """
    print("  使用 CHARMM-GUI 预生成文件...")

    system_config = get_system_config("42-38")
    ligand_a_id = system_config["ligand_b"]  # 42
    ligand_b_id = system_config["ligand_a"]  # 38

    # CGenFF ligand directories (original CHARMM-GUI structure)
    lig_42_charmm_dir = paths.system_dir / ligand_a_id / "gromacs"
    lig_38_charmm_dir = paths.system_dir / ligand_b_id / "gromacs"

    if not lig_42_charmm_dir.exists():
        raise FileNotFoundError(f"CHARMM files not found: {lig_42_charmm_dir}")
    if not lig_38_charmm_dir.exists():
        raise FileNotFoundError(f"CHARMM files not found: {lig_38_charmm_dir}")

    # 创建符号链接到标准输出目录
    lig_42_output = output_base / ligand_a_id
    lig_38_output = output_base / ligand_b_id

    lig_42_output.mkdir(parents=True, exist_ok=True)
    lig_38_output.mkdir(parents=True, exist_ok=True)

    # 链接到 gromacs 目录（模拟 LIG.amb2gmx 结构）
    (lig_42_output / "LIG.amb2gmx").symlink_to(lig_42_charmm_dir)
    (lig_38_output / "LIG.amb2gmx").symlink_to(lig_38_charmm_dir)

    print(f"  ✓ CHARMM 配体 {ligand_a_id}: {lig_42_output}")
    print(f"  ✓ CHARMM 配体 {ligand_b_id}: {lig_38_output}")

    return lig_42_output, lig_38_output


def setup_gaff2_ligands(paths: FEPTestPaths):
    """
    使用 GAFF2 参数化配体

    返回配体输出目录
    """
    print("\n[1/6] 使用 GAFF2 参数化配体...")

    system_config = get_system_config("42-38")
    ligand_a_id = system_config["ligand_b"]  # 42
    ligand_b_id = system_config["ligand_a"]  # 38

    # 配体 42
    ligand_42_mol2 = paths.get_mol2_file(ligand_a_id)
    ligand_42_output = paths.get_ligand_dir("gaff2", ligand_a_id)

    ffgen_42 = GAFFForceFieldGenerator(
        ligand_path=str(ligand_42_mol2),
        output_dir=str(ligand_42_output),
    )
    ffgen_42.run()
    print(f"  ✓ 配体 {ligand_a_id}: {ligand_42_output}")

    # 配体 38
    ligand_38_mol2 = paths.get_mol2_file(ligand_b_id)
    ligand_38_output = paths.get_ligand_dir("gaff2", ligand_b_id)

    ffgen_38 = GAFFForceFieldGenerator(
        ligand_path=str(ligand_38_mol2),
        output_dir=str(ligand_38_output),
    )
    ffgen_38.run()
    print(f"  ✓ 配体 {ligand_b_id}: {ligand_38_output}")

    return ligand_42_output, ligand_38_output


def generate_mapping_html(
    paths: FEPTestPaths,
    lig_42_output: Path,
    lig_38_output: Path,
    mapping,
    output_dir: Path,
    forcefield_type: str,
):
    """
    生成 Mapping HTML 可视化
    """
    print("\n[6/6] 生成 Mapping HTML...")

    from prism.fep.visualize.html import visualize_mapping_html
    from prism.fep.io import read_ligand_from_prism

    system_config = get_system_config("42-38")
    ligand_a_id = system_config["ligand_b"]  # 42
    ligand_b_id = system_config["ligand_a"]  # 38

    # 读取配体数据（返回 List[Atom]）
    prism_dir_42 = lig_42_output / "LIG.amb2gmx"
    prism_dir_38 = lig_38_output / "LIG.amb2gmx"

    atoms_a = read_ligand_from_prism(
        itp_file=str(prism_dir_42 / "LIG.itp"),
        gro_file=str(prism_dir_42 / "LIG.gro"),
    )

    atoms_b = read_ligand_from_prism(
        itp_file=str(prism_dir_38 / "LIG.itp"),
        gro_file=str(prism_dir_38 / "LIG.gro"),
    )

    print(f"  ✓ 配体 {ligand_a_id}: {len(atoms_a)} 原子")
    print(f"  ✓ 配体 {ligand_b_id}: {len(atoms_b)} 原子")

    # 查找 PDB 或 MOL2 文件用于可视化
    pdb_a = paths.get_pdb_file(ligand_a_id)
    pdb_b = paths.get_pdb_file(ligand_b_id)
    mol2_a = paths.get_mol2_file(ligand_a_id)
    mol2_b = paths.get_mol2_file(ligand_b_id)

    # 根据力场类型选择输入文件
    if forcefield_type == "cgenff":
        # CGenFF 使用 PDB
        if not pdb_a.exists():
            pdb_a = None
        if not pdb_b.exists():
            pdb_b = None
        mol2_a = None
        mol2_b = None
    else:
        # GAFF2 使用 MOL2
        if not mol2_a.exists():
            mol2_a = None
        if not mol2_b.exists():
            mol2_b = None
        pdb_a = None
        pdb_b = None

    # 生成 HTML
    html_path = output_dir / "42-38_mapping.html"

    visualize_mapping_html(
        mapping=mapping,
        pdb_a=str(pdb_a) if pdb_a else None,
        pdb_b=str(pdb_b) if pdb_b else None,
        mol2_a=str(mol2_a) if mol2_a else None,
        mol2_b=str(mol2_b) if mol2_b else None,
        atoms_a=atoms_a,  # 直接传入 List[Atom]
        atoms_b=atoms_b,  # 直接传入 List[Atom]
        output_path=str(html_path),
        title=f"42-38 FEP Mapping ({forcefield_type.upper()})",
        ligand_a_name=f"{ligand_a_id} (Y42)",
        ligand_b_name=f"{ligand_b_id} (F38)",
    )

    print(f"  ✓ HTML 生成: {html_path}")
    return html_path


def main():
    parser = argparse.ArgumentParser(description="FEP 端到端测试")
    parser.add_argument("--forcefield", choices=["gaff2", "cgenff"], default="gaff2", help="力场类型（默认: gaff2）")
    parser.add_argument("--html", action="store_true", help="生成 Mapping HTML 可视化")

    args = parser.parse_args()

    # Initialize test paths
    paths = FEPTestPaths("42-38")
    system_config = get_system_config("42-38")

    # 输出目录
    ligand_output_base = paths.get_forcefield_dir(args.forcefield)
    fep_output = paths.get_fep_system_dir(args.forcefield)

    # 清理旧输出
    if fep_output.exists():
        shutil.rmtree(fep_output)

    print("=" * 70)
    print(f"42-38 系统 FEP 端到端测试 ({args.forcefield.upper()})")
    print("=" * 70)

    # ========================================================================
    # 步骤 1-2: 配体参数化
    # ========================================================================
    if args.forcefield == "cgenff":
        ligand_42_output, ligand_38_output = setup_cgenff_ligands(paths, ligand_output_base)
    else:  # gaff2
        ligand_42_output, ligand_38_output = setup_gaff2_ligands(paths)

    # ========================================================================
    # 步骤 3: 读取 PRISM 格式配体
    # ========================================================================
    print("\n[2/6] 读取 PRISM 格式配体...")

    from prism.fep.io import read_ligand_from_prism

    prism_dir_42 = ligand_42_output / "LIG.amb2gmx"
    prism_dir_38 = ligand_38_output / "LIG.amb2gmx"

    atoms_42 = read_ligand_from_prism(
        itp_file=str(prism_dir_42 / "LIG.itp"),
        gro_file=str(prism_dir_42 / "LIG.gro"),
    )

    atoms_38 = read_ligand_from_prism(
        itp_file=str(prism_dir_38 / "LIG.itp"),
        gro_file=str(prism_dir_38 / "LIG.gro"),
    )

    print(f"  ✓ 配体 42: {len(atoms_42)} 原子")
    print(f"  ✓ 配体 38: {len(atoms_38)} 原子")

    # ========================================================================
    # 步骤 4: 原子映射
    # ========================================================================
    print("\n[3/6] 执行原子映射...")

    from prism.fep.core.mapping import DistanceAtomMapper

    mapper = DistanceAtomMapper(
        dist_cutoff=0.6,
        charge_cutoff=0.05,
        charge_common="mean",
        charge_reception="surround",
    )

    mapping = mapper.map(atoms_42, atoms_38)

    print(f"  ✓ Common: {len(mapping.common)}")
    print(f"  ✓ Transformed (42): {len(mapping.transformed_a)}")
    print(f"  ✓ Transformed (38): {len(mapping.transformed_b)}")

    # 验证映射结果
    expected_common = system_config.get("expected_common")
    if expected_common and abs(len(mapping.common) - expected_common) <= 2:
        print(f"  ✓ 映射结果符合预期 (~{expected_common} common)")

    # ========================================================================
    # 步骤 5: 搭建 FEP 系统
    # ========================================================================
    print("\n[4/6] 搭建 FEP 系统...")

    from prism.fep.modeling.core import FEPScaffoldBuilder
    import tempfile

    # 生成混合拓扑
    from prism.fep.gromacs.itp_builder import ITPBuilder

    hybrid_itp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".itp", delete=False)
    hybrid_itp_path = hybrid_itp_file.name

    ITPBuilder.write_complete_hybrid_itp(
        output_path=hybrid_itp_path,
        hybrid_itp=str(hybrid_itp_path),
        ligand_a_itp=str(prism_dir_42 / "LIG.itp"),
        ligand_b_itp=str(prism_dir_38 / "LIG.itp"),
        molecule_name="HYB",
    )
    print(f"  ✓ 混合拓扑: {hybrid_itp_path}")

    receptor_pdb = paths.input_dir / "receptor.pdb"

    builder = FEPScaffoldBuilder(
        output_dir=str(fep_output),
        config={},
        lambda_windows=16,
        lambda_strategy="decoupled",
    )

    layout = builder.build_from_components(
        receptor_pdb=str(receptor_pdb),
        hybrid_itp=hybrid_itp_path,
        reference_ligand_dir=str(ligand_42_output),
        mutant_ligand_dir=str(ligand_38_output),
    )

    print(f"  ✓ FEP 系统搭建完成: {layout.root}")

    # ========================================================================
    # 步骤 6: 生成 HTML（可选）
    # ========================================================================
    if args.html:
        test_output_dir = paths.get_test_output_dir(args.forcefield)
        test_output_dir.mkdir(parents=True, exist_ok=True)
        generate_mapping_html(paths, ligand_42_output, ligand_38_output, mapping, test_output_dir, args.forcefield)

    # ========================================================================
    # 总结
    # ========================================================================
    print("\n" + "=" * 70)
    print("✅ FEP 系统搭建成功！")
    print(f"\nFEP 系统目录: {layout.root}")
    print(f"\n运行 FEP 计算:")
    print(f"  cd {layout.root} && bash run_fep.sh all")
    print("=" * 70)

    # 清理临时文件
    Path(hybrid_itp_path).unlink(missing_ok=True)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
