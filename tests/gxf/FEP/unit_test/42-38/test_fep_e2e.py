#!/usr/bin/env python3
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
    python test_fep_e2e.py --forcefield charmm

    # 生成 HTML
    python test_fep_e2e.py --forcefield gaff2 --html
"""

import sys
import shutil
import argparse
from pathlib import Path

# 添加 PRISM 到路径
sys.path.insert(0, "/data2/gxf1212/work/PRISM")

from prism.forcefield.gaff import GAFFForceFieldGenerator


def setup_charmm_ligands(test_dir: Path, output_base: Path):
    """
    设置 CHARMM-GUI 配体（使用预生成的文件）

    CHARMM-GUI 文件结构：
    42/
    ├── 42.pdb, 42_modified.pdb, drawing_3D.mol2
    ├── gromacs/ (包含 LIG.itp, LIG.gro)
    └── lig/ (包含 stream files)
    """
    print("  使用 CHARMM-GUI 预生成文件...")

    # CHARMM 配体使用预生成的 PRISM 格式
    # 假设已有 gromacs/ 目录包含 LIG.itp 和 LIG.gro
    lig_42_charmm_dir = test_dir / "42" / "gromacs"
    lig_38_charmm_dir = test_dir / "38" / "gromacs"

    if not lig_42_charmm_dir.exists():
        raise FileNotFoundError(f"CHARMM files not found: {lig_42_charmm_dir}")
    if not lig_38_charmm_dir.exists():
        raise FileNotFoundError(f"CHARMM files not found: {lig_38_charmm_dir}")

    # 创建符号链接到标准输出目录
    lig_42_output = output_base / "42"
    lig_38_output = output_base / "38"

    lig_42_output.mkdir(parents=True, exist_ok=True)
    lig_38_output.mkdir(parents=True, exist_ok=True)

    # 链接到 gromacs 目录（模拟 LIG.amb2gmx 结构）
    (lig_42_output / "LIG.amb2gmx").symlink_to(lig_42_charmm_dir)
    (lig_38_output / "LIG.amb2gmx").symlink_to(lig_38_charmm_dir)

    print(f"  ✓ CHARMM 配体 42: {lig_42_output}")
    print(f"  ✓ CHARMM 配体 38: {lig_38_output}")

    return lig_42_output, lig_38_output


def setup_gaff2_ligands(test_dir: Path, output_base: Path):
    """
    使用 GAFF2 参数化配体

    返回配体输出目录
    """
    print("\n[1/6] 使用 GAFF2 参数化配体...")

    # 配体 42
    ligand_42_mol2 = test_dir / "42_3D.mol2"
    ligand_42_output = output_base / "42"

    ffgen_42 = GAFFForceFieldGenerator(
        ligand_path=str(ligand_42_mol2),
        output_dir=str(ligand_42_output),
    )
    ffgen_42.run()
    print(f"  ✓ 配体 42: {ligand_42_output}")

    # 配体 38
    ligand_38_mol2 = test_dir / "38_3D.mol2"
    ligand_38_output = output_base / "38"

    ffgen_38 = GAFFForceFieldGenerator(
        ligand_path=str(ligand_38_mol2),
        output_dir=str(ligand_38_output),
    )
    ffgen_38.run()
    print(f"  ✓ 配体 38: {ligand_38_output}")

    return ligand_42_output, ligand_38_output


def generate_mapping_html(
    test_dir: Path,
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

    print(f"  ✓ 配体 42: {len(atoms_a)} 原子")
    print(f"  ✓ 配体 38: {len(atoms_b)} 原子")

    # 查找 PDB 或 MOL2 文件用于可视化
    # 优先使用 PDB（CHARMM），其次 MOL2（GAFF2）
    pdb_a = test_dir / "42" / "42.pdb"
    pdb_b = test_dir / "38" / "38.pdb"
    mol2_a = test_dir / "42_3D.mol2"
    mol2_b = test_dir / "38_3D.mol2"

    # 根据力场类型选择输入文件
    if forcefield_type == "charmm":
        # CHARMM 使用 PDB
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
        ligand_a_name="42 (Y42)",
        ligand_b_name="38 (F38)",
    )

    print(f"  ✓ HTML 生成: {html_path}")
    return html_path


def main():
    parser = argparse.ArgumentParser(description="FEP 端到端测试")
    parser.add_argument("--forcefield", choices=["gaff2", "charmm"], default="gaff2", help="力场类型（默认: gaff2）")
    parser.add_argument("--html", action="store_true", help="生成 Mapping HTML 可视化")

    args = parser.parse_args()

    # 测试目录
    test_dir = Path("/data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/42-38")

    # 输出目录
    ligand_output_base = test_dir / f"{args.forcefield}_test"
    fep_output = test_dir / f"GMX_PROLIG_FEP_{args.forcefield}"

    # 清理旧输出
    if ligand_output_base.exists():
        shutil.rmtree(ligand_output_base)
    ligand_output_base.mkdir()

    if fep_output.exists():
        shutil.rmtree(fep_output)

    print("=" * 70)
    print(f"42-38 系统 FEP 端到端测试 ({args.forcefield.upper()})")
    print("=" * 70)

    # ========================================================================
    # 步骤 1-2: 配体参数化
    # ========================================================================
    if args.forcefield == "charmm":
        ligand_42_output, ligand_38_output = setup_charmm_ligands(test_dir, ligand_output_base)
    else:  # gaff2
        ligand_42_output, ligand_38_output = setup_gaff2_ligands(test_dir, ligand_output_base)

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
    )

    mapping = mapper.map(atoms_42, atoms_38)

    print(f"  ✓ Common atoms: {len(mapping.common)}")
    print(f"  ✓ Transformed (42): {len(mapping.transformed_a)}")
    print(f"  ✓ Transformed (38): {len(mapping.transformed_b)}")

    # ========================================================================
    # 步骤 5: 搭建 FEP 系统
    # ========================================================================
    print("\n[4/6] 搭建 FEP 系统...")

    # 生成混合拓扑
    from prism.fep.gromacs.itp_builder import ITPBuilder
    import tempfile

    hybrid_itp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".itp", delete=False)
    hybrid_itp_path = hybrid_itp_file.name

    ITPBuilder.write_complete_hybrid_itp(
        output_path=hybrid_itp_path,
        hybrid_itp=str(hybrid_itp_path),
        ligand_a_itp=str(prism_dir_42 / "LIG.itp"),
        ligand_b_itp=str(prism_dir_38 / "LIG.itp"),
        molecule_name="HYB",
    )

    print(f"  ✓ 混合拓扑生成: {hybrid_itp_path}")

    # 搭建 FEP 系统
    from prism.fep.modeling.core import FEPScaffoldBuilder

    receptor_pdb = test_dir / "receptor.pdb"

    # 读取力场对应的配置文件
    import yaml

    fep_config_path = test_dir / f"fep_{args.forcefield}.yaml"
    if not fep_config_path.exists():
        fep_config_path = test_dir / "fep.yaml"  # fallback
    with open(fep_config_path) as f:
        fep_config = yaml.safe_load(f)

    print(f"  ✓ 读取配置: {fep_config_path.name}")

    # 从 yaml 读取 lambda windows
    lambda_cfg = fep_config.get("lambda", {})
    lambda_windows = lambda_cfg.get("windows", 32)
    lambda_strategy = lambda_cfg.get("strategy", "decoupled")

    builder = FEPScaffoldBuilder(
        output_dir=str(fep_output),
        config=fep_config,
        lambda_windows=lambda_windows,
        lambda_strategy=lambda_strategy,
    )

    layout = builder.build_from_components(
        receptor_pdb=str(receptor_pdb),
        hybrid_itp=hybrid_itp_path,
        reference_ligand_dir=str(ligand_42_output),
        mutant_ligand_dir=str(ligand_38_output),
    )

    print(f"  ✓ FEP 系统搭建完成")
    print(f"    - 输出目录: {layout.root}")
    print(f"    - Bound leg: {layout.bound_dir}")
    print(f"    - Unbound leg: {layout.unbound_dir}")

    # ========================================================================
    # 步骤 6: 生成 HTML（可选）
    # ========================================================================
    html_path = None
    if args.html:
        html_path = generate_mapping_html(
            test_dir, ligand_42_output, ligand_38_output, mapping, fep_output, args.forcefield
        )

    # ========================================================================
    # 验证输出
    # ========================================================================
    print("\n" + "=" * 70)
    print("验证输出文件")
    print("=" * 70)

    import json

    checks = []

    # 混合拓扑
    checks.append(
        (
            "混合拓扑",
            [
                (layout.hybrid_dir / "hybrid.itp").exists(),
                (layout.hybrid_dir / "atomtypes_hybrid.itp").exists(),
                (layout.hybrid_dir / "ff_hybrid.itp").exists(),
            ],
        )
    )

    # Bound leg
    checks.append(
        (
            "Bound leg",
            [
                (layout.bound_dir / "topol.top").exists(),
                (layout.bound_dir / "input" / "conf.gro").exists(),
            ],
        )
    )

    # Unbound leg
    checks.append(
        (
            "Unbound leg",
            [
                (layout.unbound_dir / "topol.top").exists(),
                (layout.unbound_dir / "input" / "conf.gro").exists(),
            ],
        )
    )

    # MDP 文件
    bound_mdps = list((layout.bound_dir / "mdps").glob("prod_*.mdp"))
    checks.append((f"MDP 文件 ({len(bound_mdps)}/{lambda_windows})", [len(bound_mdps) == lambda_windows]))

    # Lambda schedule
    schedule_file = layout.bound_dir / "mdps" / "lambda_schedule.json"
    if schedule_file.exists():
        with open(schedule_file) as f:
            schedule = json.load(f)
        checks.append(
            (
                f"Lambda windows ({schedule.get('n_windows')}/{lambda_windows})",
                [schedule.get("n_windows") == lambda_windows],
            )
        )

    # 运行脚本
    checks.append(
        (
            "运行脚本",
            [
                (layout.root / "run_fep.sh").exists(),
            ],
        )
    )

    # HTML 文件
    if html_path and html_path.exists():
        checks.append(("Mapping HTML", [html_path.exists()]))

    # 打印结果
    all_pass = True
    for name, results in checks:
        if isinstance(results, list):
            if len(results) == 1 and isinstance(results[0], bool):
                status = "✅" if results[0] else "❌"
                print(f"{status} {name}")
                if not results[0]:
                    all_pass = False
            else:
                sub_pass = all(results)
                status = "✅" if sub_pass else "❌"
                print(f"{status} {name}")
                if not sub_pass:
                    all_pass = False

    print("\n" + "=" * 70)
    if all_pass:
        print("✅ 所有检查通过！FEP 系统搭建成功！")
        print(f"\n运行 FEP 计算:")
        print(f"  cd {layout.root} && bash run_fep.sh all")
        if html_path:
            print(f"\n查看 Mapping HTML:")
            print(f"  firefox {html_path}")
    else:
        print("❌ 有检查失败")
    print("=" * 70)

    # 清理临时文件
    Path(hybrid_itp_path).unlink(missing_ok=True)

    return 0 if all_pass else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
