"""
完整的原子映射测试

测试所有系统的端到端工作流，验证 GRO 文件解析修复
"""

import pytest
from pathlib import Path
import numpy as np

from prism.fep.io import read_ligand_from_prism
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.visualize import visualize_mapping_png
from prism.fep.config import FEPConfig

from test_utils import FEPTestPaths


class TestCompleteMapping:
    """完整的原子映射测试"""

    @staticmethod
    def _require_files(*paths: Path):
        missing = [str(path) for path in paths if not path.exists()]
        if missing:
            pytest.skip(f"Required test fixtures not available: {', '.join(missing)}")

    def test_42_38_system(self):
        """测试 42-38 系统（手动生成的 GRO 文件）"""
        paths = FEPTestPaths("42-38")
        test_dir = paths.system_dir

        print("\n" + "=" * 80)
        print("测试 42-38 系统（手动生成 GRO）")
        print("=" * 80)

        self._require_files(
            test_dir / "42/gromacs/LIG.itp",
            test_dir / "42/gromacs/ligand.gro",
            test_dir / "38/gromacs/LIG.itp",
            test_dir / "38/gromacs/ligand.gro",
        )

        # Load configuration using FEPConfig
        fep_config = FEPConfig(str(test_dir))
        mapping_params = fep_config.get_mapping_params()
        print(f"从 fep.yaml 读取配置: {mapping_params}")

        # 读取配体
        lig_a = read_ligand_from_prism(
            itp_file=str(test_dir / "42/gromacs/LIG.itp"), gro_file=str(test_dir / "42/gromacs/ligand.gro")
        )

        lig_b = read_ligand_from_prism(
            itp_file=str(test_dir / "38/gromacs/LIG.itp"), gro_file=str(test_dir / "38/gromacs/ligand.gro")
        )

        print(f"配体 A: {len(lig_a)} 原子")
        print(f"配体 B: {len(lig_b)} 原子")

        # 验证坐标不是零
        assert not np.allclose(lig_a[0].coord, [0, 0, 0]), "配体 A 坐标不应该是零"
        assert not np.allclose(lig_b[0].coord, [0, 0, 0]), "配体 B 坐标不应该是零"

        # 执行映射
        mapper = DistanceAtomMapper(**mapping_params)
        mapping = mapper.map(lig_a, lig_b)

        print(
            f"映射结果: Common={len(mapping.common)}, "
            f"Transformed=A{len(mapping.transformed_a)}+B{len(mapping.transformed_b)}, "
            f"Surrounding=A{len(mapping.surrounding_a)}+B{len(mapping.surrounding_b)}"
        )

        # 验证映射结果
        total_a = len(lig_a)
        total_b = len(lig_b)
        sum_a = len(mapping.common) + len(mapping.transformed_a) + len(mapping.surrounding_a)
        sum_b = len(mapping.common) + len(mapping.transformed_b) + len(mapping.surrounding_b)

        assert sum_a == total_a, f"配体 A 原子数不匹配: {sum_a} != {total_a}"
        assert sum_b == total_b, f"配体 B 原子数不匹配: {sum_b} != {total_b}"

        # 氢原子统计
        h_common = sum(1 for a, b in mapping.common if a.element == "H")
        h_a = sum(1 for a in lig_a if a.element == "H")
        h_b = sum(1 for a in lig_b if a.element == "H")

        print(f"氢原子: A={h_a}, B={h_b}, Common={h_common}")
        print(f"氢原子映射率: {h_common/h_a*100:.1f}%")

        print("✅ 42-38 系统测试通过")

    def test_39_8_system(self):
        """测试 39-8 系统（PRISM 生成的 GRO 文件）"""
        test_dir = FEPTestPaths("39-8").system_dir

        print("\n" + "=" * 80)
        print("测试 39-8 系统（PRISM 生成 GRO）")
        print("=" * 80)

        self._require_files(
            test_dir / "39/LIG.charmm2gmx/LIG.itp",
            test_dir / "39/LIG.charmm2gmx/LIG.gro",
            test_dir / "8/LIG.charmm2gmx/LIG.itp",
            test_dir / "8/LIG.charmm2gmx/LIG.gro",
        )

        # Load configuration using FEPConfig
        fep_config = FEPConfig(str(test_dir))
        mapping_params = fep_config.get_mapping_params()
        print(f"从 fep.yaml 读取配置: {mapping_params}")

        # 读取配体
        lig_a = read_ligand_from_prism(
            itp_file=str(test_dir / "39/LIG.charmm2gmx/LIG.itp"), gro_file=str(test_dir / "39/LIG.charmm2gmx/LIG.gro")
        )

        lig_b = read_ligand_from_prism(
            itp_file=str(test_dir / "8/LIG.charmm2gmx/LIG.itp"), gro_file=str(test_dir / "8/LIG.charmm2gmx/LIG.gro")
        )

        print(f"配体 A: {len(lig_a)} 原子")
        print(f"配体 B: {len(lig_b)} 原子")

        # 验证坐标不是零
        assert not np.allclose(lig_a[0].coord, [0, 0, 0]), "配体 A 坐标不应该是零"
        assert not np.allclose(lig_b[0].coord, [0, 0, 0]), "配体 B 坐标不应该是零"

        print(f"配体 A 前3个原子坐标:")
        for atom in lig_a[:3]:
            print(f"  {atom.name}: ({atom.coord[0]:.3f}, {atom.coord[1]:.3f}, {atom.coord[2]:.3f}) Å")

        # 执行映射
        mapper = DistanceAtomMapper(**mapping_params)
        mapping = mapper.map(lig_a, lig_b)

        print(
            f"映射结果: Common={len(mapping.common)}, "
            f"Transformed=A{len(mapping.transformed_a)}+B{len(mapping.transformed_b)}, "
            f"Surrounding=A{len(mapping.surrounding_a)}+B{len(mapping.surrounding_b)}"
        )

        # 验证映射结果
        total_a = len(lig_a)
        total_b = len(lig_b)
        sum_a = len(mapping.common) + len(mapping.transformed_a) + len(mapping.surrounding_a)
        sum_b = len(mapping.common) + len(mapping.transformed_b) + len(mapping.surrounding_b)

        assert sum_a == total_a, f"配体 A 原子数不匹配: {sum_a} != {total_a}"
        assert sum_b == total_b, f"配体 B 原子数不匹配: {sum_b} != {total_b}"

        # 氢原子统计
        h_common = sum(1 for a, b in mapping.common if a.element == "H")
        h_a = sum(1 for a in lig_a if a.element == "H")
        h_b = sum(1 for a in lig_b if a.element == "H")

        print(f"氢原子: A={h_a}, B={h_b}, Common={h_common}")
        print(f"氢原子映射率: {h_common/h_a*100:.1f}%")

        # 分析 Surrounding 氢原子距离
        surrounding_h_a = [atom for atom in mapping.surrounding_a if atom.element == "H"]
        if surrounding_h_a:
            print(f"\nSurrounding 氢原子距离分析:")
            h_atoms_b = [a for a in lig_b if a.element == "H"]
            for atom in surrounding_h_a[:3]:  # 只显示前3个
                min_dist = min(np.linalg.norm(atom.coord - hb.coord) for hb in h_atoms_b)
                print(f"  {atom.name}: 最近距离 {min_dist:.3f} Å (阈值: 0.6 Å)")

        print("✅ 39-8 系统测试通过")

    def test_visualization_generation(self):
        """测试可视化生成（PNG + HTML）"""
        paths = FEPTestPaths("42-38")
        test_dir = paths.system_dir
        output_dir = paths.get_test_output_dir("cgenff")
        output_dir.mkdir(exist_ok=True)

        self._require_files(
            test_dir / "42/gromacs/LIG.itp",
            test_dir / "42/gromacs/ligand.gro",
            test_dir / "38/gromacs/LIG.itp",
            test_dir / "38/gromacs/ligand.gro",
        )

        print("\n" + "=" * 80)
        print("测试可视化生成（PNG + HTML）")
        print("=" * 80)

        # Load configuration using FEPConfig
        fep_config = FEPConfig(str(test_dir))
        mapping_params = fep_config.get_mapping_params()

        # 读取配体
        lig_a = read_ligand_from_prism(
            itp_file=str(test_dir / "42/gromacs/LIG.itp"), gro_file=str(test_dir / "42/gromacs/ligand.gro")
        )

        lig_b = read_ligand_from_prism(
            itp_file=str(test_dir / "38/gromacs/LIG.itp"), gro_file=str(test_dir / "38/gromacs/ligand.gro")
        )

        # 执行映射
        mapper = DistanceAtomMapper(**mapping_params)
        mapping = mapper.map(lig_a, lig_b)

        # 生成 PNG 可视化
        png_file = str(output_dir / "42-38_mapping.png")
        visualize_mapping_png(
            mapping=mapping,
            pdb_a=str(paths.get_pdb_file("42")),
            pdb_b=str(paths.get_pdb_file("38")),
            mol2_a=str(paths.get_mol2_file("42")),
            mol2_b=str(paths.get_mol2_file("38")),
            output_path=png_file,
        )

        assert Path(png_file).exists(), f"PNG 文件未生成: {png_file}"
        print(f"✅ PNG 可视化生成成功: {png_file}")

        # 生成 HTML 可视化
        from prism.fep.visualize import visualize_mapping_html

        html_file = str(output_dir / "42-38_mapping.html")

        # Use HTML config from FEPConfig
        html_config = fep_config.get_html_config()

        visualize_mapping_html(
            mapping=mapping,
            pdb_a=str(paths.get_pdb_file("42")),
            pdb_b=str(paths.get_pdb_file("38")),
            mol2_a=str(paths.get_mol2_file("42")),
            mol2_b=str(paths.get_mol2_file("38")),
            atoms_a=lig_a,
            atoms_b=lig_b,
            output_path=html_file,
            ligand_a_name="Ligand 42",
            ligand_b_name="Ligand 38",
            config=html_config,
        )

        assert Path(html_file).exists(), f"HTML 文件未生成: {html_file}"
        print(f"✅ HTML 可视化生成成功: {html_file}")

        print("✅ 可视化测试通过（PNG + HTML）")


if __name__ == "__main__":
    # 快速测试
    test = TestCompleteMapping()
    test.test_42_38_system()
    test.test_39_8_system()
    test.test_visualization_generation()

    print("\n" + "=" * 80)
    print("🎉 所有测试通过！GRO 文件解析修复成功！")
    print("=" * 80)
