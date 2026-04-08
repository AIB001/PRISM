#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
通用 FEP 系统完整流程测试
验证 PRISMBuilder.run_fep() 和 bound/unbound grompp

用法：
    # 从系统目录运行（自动检测配体文件）
    cd tests/gxf/FEP/unit_test/42-38 && python ../../test_run_fep.py
    cd tests/gxf/FEP/unit_test/25-36 && python ../../test_run_fep.py
    cd tests/gxf/FEP/unit_test/39-8 && python ../../test_run_fep.py

    # 指定力场
    python ../../test_run_fep.py --forcefield amber14sb_OL15
"""

import os
import sys
import argparse
import re
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from prism.builder.core import PRISMBuilder


def auto_detect_ligands(input_dir="input"):
    """自动检测配体文件

    优先级：
    1. 同名 *.mol2（通用；与 *.pdb 同 stem，优先用于 GAFF/OpenFF/OPLS）
    2. *_3d_pdb.mol2 (转换后的 MOL2，保留对齐坐标)
    3. *.pdb (对齐的 PDB)
    4. 其他 *.mol2
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        return None, None

    # 优先使用同名 mol2（例如 42.pdb -> 42.mol2）。
    # 这适用于 GAFF/OpenFF/OPLS 的常规流程，避免在目录同时存在
    # PDB 和 MOL2 时回退到 PDB，丢失已经准备好的 MOL2 键信息/坐标。
    plain_mol2 = sorted([f for f in input_path.glob("*.mol2") if "_3d" not in f.stem.lower() and "_3D" not in f.stem])
    if len(plain_mol2) >= 2:
        print("✓ 检测到同名 MOL2 文件")
        return _order_pair_by_case_name(input_path, plain_mol2)

    # 查找转换后的 MOL2 文件
    mol2_files = sorted(input_path.glob("*_3d_pdb.mol2"))
    if len(mol2_files) >= 2:
        print(f"✓ 检测到转换后的 MOL2 文件（保留对齐坐标）")
        return _order_pair_by_case_name(input_path, mol2_files)

    # 查找 PDB 文件
    pdb_files = sorted([f for f in input_path.glob("*.pdb") if "receptor" not in f.name.lower()])
    if len(pdb_files) >= 2:
        print(f"✓ 检测到 PDB 文件")
        return _order_pair_by_case_name(input_path, pdb_files)

    # 查找 MOL2 文件
    mol2_files = sorted(input_path.glob("*.mol2"))
    if len(mol2_files) >= 2:
        print(f"✓ 检测到 MOL2 文件")
        return _order_pair_by_case_name(input_path, mol2_files)

    return None, None


def _order_pair_by_case_name(input_path: Path, files):
    """Use directory name like ``42-38`` to preserve reference→mutant order."""
    case_name = input_path.resolve().parent.name
    if "-" in case_name:
        left, right = case_name.split("-", 1)
        stem_map = {f.stem.split("_")[0]: f.name for f in files}
        if left in stem_map and right in stem_map:
            return stem_map[left], stem_map[right]
    return files[0].name, files[1].name


def get_system_name():
    """从当前目录获取系统名称"""
    cwd = Path.cwd()
    return cwd.name


_LFF_CONFIG_MAP = {
    "gaff": ("case_gaff.yaml", "fep_gaff.yaml"),
    "gaff2": ("case_gaff2.yaml", "fep_gaff2.yaml"),
    "openff": ("case_openff.yaml", "fep_openff.yaml"),
    "opls": ("case_opls.yaml", "case_opls.yaml"),
    "oplsaa": ("case_opls.yaml", "case_opls.yaml"),
    "oplsaam": ("case_oplsaam.yaml", "case_oplsaam.yaml"),
    "cgenff": ("case_charmm.yaml", "fep_charmm.yaml"),
    "rtf": ("case_rtf.yaml", "fep_rtf.yaml"),
}

_LFF_FF_MAP = {
    "opls": "oplsaa",
    "oplsaa": "oplsaa",
    "oplsaam": "oplsaam",
}

_UNSUPPORTED_FEP_PROTEIN_FFS = {
    "oplsaam": (
        "OPLSAAM is currently unsupported for PRISM FEP system setup. "
        "The available oplsaam.ff dataset is missing required protein Ryckaert-Bellemans dihedral types, "
        "so bound-system grompp fails before FEP scaffold generation."
    ),
}


def _slug_case_token(value: str) -> str:
    token = re.sub(r"[^a-zA-Z0-9]+", "_", value.strip().lower())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "case"


def _default_case_dir(protein_ff: str, ligand_ff: str) -> str:
    return f"{_slug_case_token(protein_ff)}-mut_{_slug_case_token(ligand_ff)}"


def _resolve_configs(ligand_forcefield: str, forcefield: str):
    """Return (case_config, fep_config, protein_ff, water_model) for a ligand FF."""
    lff = ligand_forcefield.lower()
    protein_ff = forcefield or _LFF_FF_MAP.get(lff, ligand_forcefield)
    if protein_ff.lower() == "oplsaam":
        case_cfg, fep_cfg = ("case_oplsaam.yaml", "case_oplsaam.yaml")
    else:
        case_cfg, fep_cfg = _LFF_CONFIG_MAP.get(lff, ("case_gaff2.yaml", "fep_gaff2.yaml"))

    # Try to fall back to a gaff2 config pair when no dedicated config exists
    configs_dir = Path("configs")
    if not (configs_dir / case_cfg).exists():
        fallback = configs_dir / "case_gaff2.yaml"
        if fallback.exists():
            case_cfg = "case_gaff2.yaml"
    if not (configs_dir / fep_cfg).exists():
        fallback = configs_dir / "fep_gaff2.yaml"
        if fallback.exists():
            fep_cfg = "fep_gaff2.yaml"
        else:
            fep_cfg = case_cfg  # last resort

    water_model = "tip3p"
    return case_cfg, fep_cfg, protein_ff, water_model


def test_run_fep(
    forcefield="amber14sb_OL15", ligand_forcefield="gaff2", ref_ligand=None, mut_ligand=None, forcefield_paths=None
):
    """测试完整的 FEP 流程（run_fep）"""

    system_name = get_system_name()

    print("\n" + "=" * 70)
    print(f"{system_name} 系统 FEP 完整流程测试 ({forcefield} / {ligand_forcefield})")
    print("=" * 70)

    # RTF: ligand_path 传子目录（含 .rtf/.prm/.pdb），从 case_name 推断目录名
    if ligand_forcefield.lower() == "rtf":
        if ref_ligand is None or mut_ligand is None:
            # Find subdirectories in input/ that contain RTF files
            input_dirs = sorted([d.name for d in Path("input").iterdir() if d.is_dir() and any(d.glob("*.rtf"))])
            if len(input_dirs) >= 2:
                ref_ligand = f"input/{input_dirs[0]}"
                mut_ligand = f"input/{input_dirs[1]}"
            else:
                print(f"❌ RTF 模式需要手动指定 --ref-ligand 和 --mut-ligand（目录路径）")
                return False
        print(f"  Reference ligand dir: {ref_ligand}")
        print(f"  Mutant ligand dir:    {mut_ligand}")
    else:
        # 自动检测配体文件
        if ref_ligand is None or mut_ligand is None:
            ref_ligand, mut_ligand = auto_detect_ligands()

        if ref_ligand is None or mut_ligand is None:
            print(f"❌ 无法检测到配体文件")
            print(f"   请确保 input/ 目录包含两个配体文件：")
            print(f"   - *_3d_pdb.mol2 或 *.pdb 或 *.mol2")
            return False

    print(f"  Reference ligand: {ref_ligand}")
    print(f"  Mutant ligand:    {mut_ligand}")

    case_cfg, fep_cfg, protein_ff, water_model = _resolve_configs(ligand_forcefield, forcefield)
    unsupported_reason = _UNSUPPORTED_FEP_PROTEIN_FFS.get(protein_ff.lower())
    if unsupported_reason:
        print(f"  ✗ Unsupported protein force field for PRISM FEP: {protein_ff}")
        print(f"    {unsupported_reason}")
        return False

    print(f"  Case config:      configs/{case_cfg}")
    print(f"  FEP config:       configs/{fep_cfg}")
    print(f"  Protein FF:       {protein_ff}  Water: {water_model}")

    # 输出目录：默认使用 <protein_ff>-mut_<ligand_ff>
    output_dir = _default_case_dir(protein_ff, ligand_forcefield)

    # RTF already has "input/" prefix; others need it added
    if ligand_forcefield.lower() == "rtf":
        ref_ligand_path = ref_ligand
        mut_ligand_path = mut_ligand
    else:
        if os.path.isabs(ref_ligand) or str(ref_ligand).startswith("input/"):
            ref_ligand_path = ref_ligand
        else:
            ref_ligand_path = f"input/{ref_ligand}"
        if os.path.isabs(mut_ligand) or str(mut_ligand).startswith("input/"):
            mut_ligand_path = mut_ligand
        else:
            mut_ligand_path = f"input/{mut_ligand}"

    # Protein PDB: support both receptor.pdb and protein.pdb
    protein_path = "input/receptor.pdb"
    if not Path(protein_path).exists():
        protein_path = "input/protein.pdb"

    builder = PRISMBuilder(
        protein_path=protein_path,
        ligand_paths=[ref_ligand_path],
        output_dir=output_dir,
        ligand_forcefield=ligand_forcefield,
        config_path=f"configs/{case_cfg}",
        forcefield=protein_ff,
        water_model=water_model,
        overwrite=True,
        fep_mode=True,
        mutant_ligand=mut_ligand_path,
        fep_config=f"configs/{fep_cfg}",
        forcefield_path=forcefield_paths,
    )

    print(f"\n[Step 1] 运行 PRISMBuilder.run()...")
    result = builder.run()
    print(f"  ✓ Run completed: {result}")

    # 检查输出目录
    fep_dir = Path(output_dir) / "GMX_PROLIG_FEP"

    print(f"\n[Step 2] 检查 FEP 目录结构...")

    bound_dir = fep_dir / "bound" / "repeat1"
    unbound_dir = fep_dir / "unbound" / "repeat1"

    # 如果不存在，检查旧版本结构
    if not bound_dir.exists():
        bound_dir = fep_dir / "bound"
    if not unbound_dir.exists():
        unbound_dir = fep_dir / "unbound"

    for leg_name, leg_dir in [("Bound", bound_dir), ("Unbound", unbound_dir)]:
        if not leg_dir.exists():
            print(f"  ❌ {leg_name} 目录不存在: {leg_dir}")
            return False

        print(f"  ✓ {leg_name} 目录存在")

        # 检查关键文件
        conf_gro = leg_dir / "input" / "conf.gro"
        topol_top = leg_dir / "topol.top"

        if not conf_gro.exists():
            print(f"    ❌ {leg_name} conf.gro 不存在")
            return False
        else:
            # 检查不是占位符
            content = conf_gro.read_text()
            if "Placeholder" in content or len(content) < 100:
                print(f"    ⚠ {leg_name} conf.gro 是占位符或空文件")
            else:
                # 检查原子数量
                lines = content.split("\n")
                if len(lines) >= 2:
                    atom_count_line = lines[1].strip()
                    try:
                        atom_count = int(atom_count_line)
                        print(f"    ✓ {leg_name} conf.gro: {atom_count} 原子（真实系统）")
                    except ValueError:
                        print(f"    ⚠ {leg_name} conf.gro 格式异常")

        if not topol_top.exists():
            print(f"    ❌ {leg_name} topol.top 不存在")
            return False
        else:
            print(f"    ✓ {leg_name} topol.top 存在")

        print(f"    运行 {leg_name} grompp 验证...")
        tpr_name = "em_check.tpr"
        grompp_cmd = [
            "gmx",
            "grompp",
            "-f",
            "mdps/em.mdp",
            "-c",
            "input/conf.gro",
            "-r",
            "input/conf.gro",
            "-p",
            "topol.top",
            "-o",
            tpr_name,
            "-maxwarn",
            "2",
        ]
        result = subprocess.run(
            grompp_cmd,
            cwd=leg_dir,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"    ❌ {leg_name} grompp 失败")
            if result.stderr.strip():
                print(result.stderr.strip()[-1200:])
            if result.stdout.strip():
                print(result.stdout.strip()[-1200:])
            return False
        print(f"    ✓ {leg_name} grompp 通过: {leg_dir / tpr_name}")

    print(f"\n[Step 3] 检查 Mapping HTML...")
    mapping_html = fep_dir / "common" / "hybrid" / "mapping.html"
    if mapping_html.exists():
        print(f"    ✓ Mapping HTML 已生成: {mapping_html}")

        # 分析 mapping 质量
        content = mapping_html.read_text()

        atoms_a_match = re.search(r"const ATOMS_A = (\[.*?\]);", content, re.DOTALL)
        atoms_b_match = re.search(r"const ATOMS_B = (\[.*?\]);", content, re.DOTALL)

        if atoms_a_match and atoms_b_match:
            import json

            atoms_a = json.loads(atoms_a_match.group(1))
            atoms_b = json.loads(atoms_b_match.group(1))

            def count_classifications(atom_records):
                common = sum(1 for atom in atom_records if atom.get("classification") == "common")
                transformed = sum(1 for atom in atom_records if atom.get("classification") == "transformed")
                surrounding = sum(1 for atom in atom_records if atom.get("classification") == "surrounding")
                return common, transformed, surrounding

            common_a, trans_a, surr_a = count_classifications(atoms_a)
            common_b, trans_b, surr_b = count_classifications(atoms_b)

            total_a = common_a + trans_a + surr_a
            total_b = common_b + trans_b + surr_b

            print(f"\n    Mapping 统计:")
            print(f"      Ligand A: common={common_a}, transformed={trans_a}, surrounding={surr_a}, total={total_a}")
            print(f"      Ligand B: common={common_b}, transformed={trans_b}, surrounding={surr_b}, total={total_b}")

            # 评估质量
            common_count = max(common_a, common_b)
            if common_count >= 30:
                quality = "✅ 优秀"
            elif common_count >= 20:
                quality = "⚠️  一般"
            else:
                quality = "❌ 较差"

            print(f"      质量评估: {quality} ({common_count} common atoms)")
    else:
        print(f"    ❌ Mapping HTML 未生成: {mapping_html}")
        return False

    print("\n" + "=" * 70)
    print("✅ FEP 完整流程测试通过！")
    print("=" * 70)
    print(f"\n输出目录: {fep_dir.absolute()}")
    print(f"\nMapping HTML: {mapping_html.absolute()}")
    print(f"\n下一步：运行实际的 FEP 计算")
    print(f"  cd {fep_dir} && bash run_fep.sh all")

    return True


def main():
    parser = argparse.ArgumentParser(description="通用 FEP 系统完整流程测试")
    parser.add_argument("--forcefield", default="amber14sb_OL15", help="蛋白力场名称（默认: amber14sb_OL15）")
    parser.add_argument(
        "--ligand-forcefield",
        default="gaff2",
        dest="ligand_forcefield",
        help="配体力场名称（默认: gaff2）。支持: gaff, gaff2, openff, opls/oplsaa, cgenff",
    )
    parser.add_argument("--ref-ligand", default=None, help="Reference ligand 文件名（默认：自动检测）")
    parser.add_argument("--mut-ligand", default=None, help="Mutant ligand 文件名（默认：自动检测）")
    parser.add_argument(
        "--forcefield-path",
        "-ffp",
        action="append",
        dest="forcefield_paths",
        default=None,
        help="CGenFF 预生成目录路径（cgenff 力场时必须提供两个，ref 在前，mut 在后）",
    )

    args = parser.parse_args()
    success = test_run_fep(
        args.forcefield, args.ligand_forcefield, args.ref_ligand, args.mut_ligand, args.forcefield_paths
    )
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
