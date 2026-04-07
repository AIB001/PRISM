---
name: fep-forcefield-integration
description: FEP module integration testing with various force fields
type: fep
---

# FEP 模块与各种力场的集成测试

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

## 核心职责

测试 PRISM 的所有力场生成器是否都能正确生成 PRISM 格式，并能够与 FEP 模块无缝集成。验证从"普通建模 → FEP 模式"的完整工作流。

## 命名规则 (CRITICAL)

- Case 目录统一使用 `<protein_ff>-mut_<ligand_ff>`（示例：`charmm36m-mut_mmff`）。
- 禁止 `.../GMX_PROLIG_FEP/GMX_PROLIG_FEP/` 重复嵌套目录。
- 禁止 `_pkgfix*`、`_final*`、`_new*` 这类临时后缀。

## 设计原则

**统一接口，解耦设计**：
```
原始配体文件 (MOL2/SDF/PDB/RTF+PRM)
    ↓ [力场生成器]
PRISM 格式 (ITP + GRO)
    ↓ [FEP 模块]
原子映射 + 可视化
```

FEP 模块**只关心 PRISM 格式**，不关心原始文件类型。

## 已支持的力场

| 力场 | 生成器类 | 输出目录 | 方法 | 状态 |
|------|----------|----------|------|------|
| GAFF | `GAFFForceFieldGenerator` | `LIG.amb2gmx` | 本地 (AmberTools) | ✅ |
| GAFF2 | `GAFF2ForceFieldGenerator` | `LIG.gaff2` | 本地 (AmberTools) | ✅ |
| OpenFF | `OpenFFForceFieldGenerator` | `LIG.openff2gmx` | 本地 | ✅ |
| OPLS-AA | `OPLSAAForceFieldGenerator` | `LIG.opls2gmx` | 服务器 (LigParGen) | ✅ |
| MMFF | `MMFFForceFieldGenerator` | `LIG.mmff2gmx` | 服务器 (SwissParam) | ✅ |
| MATCH | `MATCHForceFieldGenerator` | `LIG.match2gmx` | 服务器 (SwissParam) | ✅ |
| Both | `BothForceFieldGenerator` | `LIG.both2gmx` | 服务器 (SwissParam) | ✅ |
| CGenFF | `CGenFFForceFieldGenerator` | `LIG.cgenff2gmx` | 本地/网站 | ✅ |
| CHARMM-GUI | `CHARMMLegacyGenerator` | `LIG.charmm2gmx` | 本地 | ✅ |
| RTF | `RTFForceFieldGenerator` | `LIG.rtf2gmx` | 本地 (CHARMM RTF) | ✅ |

## PRISM 标准输出格式

所有力场生成器必须产生：

```
output_dir/LIG.<ff>2gmx/
├── LIG.gro              # 坐标文件（必需）
├── LIG.itp              # 分子拓扑（必需）
├── LIG.top              # 完整拓扑（必需）
├── atomtypes_LIG.itp    # 原子类型定义（必需）
└── posre_LIG.itp        # 位置限制（可选）
```

**FEP 模块只需要**：`LIG.itp` + `LIG.gro`

## 测试矩阵

### 测试案例：真实测试系统

| 力场 | 测试文件 | 测试案例 | 状态 | 测试命令 |
|------|----------|----------|------|----------|
| **CHARMM-GUI** | `test_charmm_gui_generator.py` | 39-8 + 42-38 | ✅ **已验证** (灰色原子已修复) | `python test_39_8_pdb_only.py`<br>`python test_42_38_gray_atom.py` |
| **RTF** | `test_run_fep.py` | p38-19-24 | ✅ **已验证** (2026-04-07) | `cd tests/gxf/FEP/unit_test/p38-19-24 && python ../test_run_fep.py --forcefield charmm36-jul2022 --ligand-forcefield rtf` |
| OpenFF | `test_openff_fep_mapping.py` | oMeEtPh-EtPh | ✅ **已验证** | `pytest test_openff_fep_mapping.py::TestOpenFFWithFEPMapping::test_end_to_end_fep_mapping` |
| GAFF | `test_gaff_fep_integration.py` | oMeEtPh-EtPh | ⚠️ 需要创建 | - |
| GAFF2 | `test_gaff2_fep_integration.py` | oMeEtPh-EtPh | ⚠️ 需要创建 | - |
| OPLS-AA | `test_opls_fep_integration.py` | oMeEtPh-EtPh | ⚠️ 需要创建 | - |
| CGenFF (网站) | `test_e2e_prism_fep.py` | 39-8 系统 | ⚠️ **被跳过** | 需要 CGenFF 网站文件 |

## 标准测试模板

每个力场的测试应该遵循以下模板：

```python
#!/usr/bin/env python3
"""
测试 <力场名> 与 FEP 模块的集成

完整流程：
1. 使用 <力场名> 生成器生成两个配体的 PRISM 格式
2. 使用 read_ligand_from_prism() 读取
3. 使用 DistanceAtomMapper 执行映射
4. 生成可视化（PNG + HTML）
"""

import pytest
from pathlib import Path
from prism.forcefield.<ff_module> import <FFGenerator>
from prism.fep.io import read_ligand_from_prism
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.visualize.html import visualize_mapping_html

def test_<ff>_fep_e2e(tmp_path):
    """端到端测试：<力场名> 生成 + FEP 映射"""

    # 测试数据
    mol2_a = "examples/ligands/oMeEtPh.mol2"
    mol2_b = "examples/ligands/EtPh.mol2"

    # Step 1: 生成配体 A 的力场
    output_a = tmp_path / "<ff>_oMeEtPh"
    generator_a = <FFGenerator>(
        ligand_path=mol2_a,
        output_dir=str(output_a),
        ...  # 力场特定参数
    )
    result_dir_a = generator_a.run()

    # Step 2: 生成配体 B 的力场
    output_b = tmp_path / "<ff>_EtPh"
    generator_b = <FFGenerator>(
        ligand_path=mol2_b,
        output_dir=str(output_b),
        ...  # 力场特定参数
    )
    result_dir_b = generator_b.run()

    # Step 3: 读取 PRISM 格式文件
    itp_a = Path(result_dir_a) / "LIG.itp"
    gro_a = Path(result_dir_a) / "LIG.gro"
    itp_b = Path(result_dir_b) / "LIG.itp"
    gro_b = Path(result_dir_b) / "LIG.gro"

    atoms_a = read_ligand_from_prism(str(itp_a), str(gro_a))
    atoms_b = read_ligand_from_prism(str(itp_b), str(gro_b))

    # Step 4: 执行原子映射
    mapper = DistanceAtomMapper(
        dist_cutoff=0.6,
        charge_cutoff=0.05,
        charge_common="mean"
    )
    mapping = mapper.map(atoms_a, atoms_b)

    # Step 5: 生成可视化
    output_html = tmp_path / "oMeEtPh-EtPh_<ff>_mapping.html"
    visualize_mapping_html(
        mapping=mapping,
        pdb_a="examples/ligands/oMeEtPh.pdb",
        pdb_b="examples/ligands/EtPh.pdb",
        mol2_a=mol2_a,
        mol2_b=mol2_b,
        atoms_a=atoms_a,
        atoms_b=atoms_b,
        output_path=str(output_html),
        title=f"<力场名> FEP Mapping: oMeEtPh → EtPh"
    )

    # 验证
    assert len(mapping.common) >= 15, "应该有足够的 common 原子"
    assert output_html.exists(), "HTML 文件应该生成"

    print(f"✓ {len(mapping.common)} common atoms")
    print(f"✓ Visualization: {output_html}")
```

## CLI 集成测试

测试 CLI 是否支持 `--fep` 模式：

```bash
# GAFF

## 关键验证点

### 1. PRISM 格式验证
- [ ] 所有力场都能生成 `LIG.itp` + `LIG.gro`
- [ ] ITP 文件包含 `[ atoms ]` 部分
- [ ] 原子类型和电荷正确

### 2. FEP 读取验证
- [ ] `read_ligand_from_prism()` 能读取所有力场的输出
- [ ] 坐标正确转换（nm → Å）
- [ ] 原子数量正确

### 3. 映射算法验证
- [ ] `DistanceAtomMapper` 对所有力场都能工作
- [ ] Common atoms 数量合理
- [ ] Transform/Surrounding 分类正确
- [ ] **所有原子都被分类（无灰色原子）**
- [ ] **电荷为0的原子数量在合理范围内**

### 4. 可视化验证
- [ ] HTML 报告生成成功
- [ ] PNG 可视化生成成功（或降级处理）
- [ ] 原子标签正确显示
- [ ] **未分类原子触发warning**
- [ ] **灰色原子显示warning提示**

### 5. CLI 集成验证
- [ ] `--fep --mutant` 参数对所有力场有效
- [ ] 输出目录结构正确
- [ ] 配置文件正确加载

### 6. 灰色原子检测（新增）
- [ ] 检查HTML中是否存在 `rgb(200, 200, 200)` 灰色原子
- [ ] 验证所有原子都有正确的classification（common/transformed/surrounding）
- [ ] 确认映射算法覆盖率：`total = common + transformed_a + transformed_b + surrounding_a + surrounding_b`
- [ ] 如果存在未分类原子，HTML应该显示warning框

**详细HTML检查方法**：
```bash
# 1. 检查灰色原子

### 7. 前端Warning机制（新增）
当出现以下情况时，HTML必须在**ligand可视化附近**显示warning框：

#### 触发条件
- 存在未分类原子（classification为unknown或缺失）
- 电荷为0的原子数量异常（>总原子数的10%）
- RDKit 2D对齐失败（PNG降级）
- 映射结果不合理（如common数量过少）

#### Warning样式
```html
<div class="warning-banner" style="background: #fff3cd; border: 1px solid #ffc107; padding: 15px; margin: 20px 0; border-radius: 8px;">
  <strong>⚠️ Warning:</strong> <span id="warning-message">...</span>
</div>
```

#### Warning位置
- 在分子可视化canvas附近（上方或下方）
- 使用醒目的颜色（橙色/红色背景）
- 不应该影响可视化，但必须显眼

## 力场特定注意事项

### GAFF/GAFF2
- 需要 AmberTools 安装
- 电荷模型：AM1-BCC（默认）或 RESP
- 原子类型：包含位置信息（ca, c3, ha, hc）

### OpenFF
- 需要 openff-toolkit 安装
- 支持多个 OpenFF 版本（2.0.0, 2.1.0, etc.）
- 原子类型通用（output_0, output_1）
- **特殊处理**：`ignore_atom_type=True`（因为类型不包含位置信息）

### OPLS-AA
- 需要 LigParGen 服务器访问
- 或使用本地 OPLS-AA 安装
- 电荷模型：CM1A 或 CM5

### MMFF/MATCH/Both (SwissParam)
- 需要 SwissParam 服务器访问
- 有请求频率限制
- 网络超时处理

### CGenFF
- 需要 CGenFF 网站文件（*_gmx.pdb + *_gmx.top）
- 或 CHARMM-GUI 文件（RTF + PRM）
- 参数化质量高，但流程复杂

### CHARMM-GUI
- 使用 CHARMM-GUI 生成的 PSF/PRM
- 需要转换到 GROMACS 格式
- 依赖 CHARMM 力场文件

## 测试优先级

### 优先级 1：核心力场（必须测试）
1. **GAFF** - 最常用，AmberTools 本地运行
2. **OpenFF** - 现代化力场，本地运行
3. **CGenFF** - CHARMM 生态，高质量参数

### 优先级 2：扩展力场（重要）
4. **OPLS-AA** - OPLS 生态
5. **GAFF2** - GAFF 改进版

### 优先级 3：服务器力场（可选）
6. **MMFF** - SwissParam 快速参数化
7. **MATCH** - SwissParam 高质量参数

## 相关代码路径

- 力场生成器基类：`prism/forcefield/base.py::ForceFieldGeneratorBase`
- FEP IO 模块：`prism/fep/io.py::read_ligand_from_prism()`
- 映射算法：`prism/fep/core/mapping.py::DistanceAtomMapper`
- CLI 参数：`prism/builder/cli.py` (line 343-372)
- 系统构建：`prism/utils/system/builder.py::SystemBuilder`

## 实施计划

### 阶段 1：创建标准测试（1-2 天）
- [ ] 创建 `test_gaff_fep_integration.py`
- [ ] 创建 `test_gaff2_fep_integration.py`
- [ ] 创建 `test_opls_fep_integration.py`

### 阶段 2：运行并修复（2-3 天）
- [ ] 运行所有测试
- [ ] 修复发现的 bug
- [ ] 完善错误处理

### 阶段 3：CLI 集成测试（1-2 天）
- [ ] 测试 `--fep --mutant` 参数
- [ ] 验证输出目录结构
- [ ] 测试配置文件加载

### 阶段 4：文档和示例（1 天）
- [ ] 更新 README
- [ ] 添加示例脚本
- [ ] 创建 troubleshooting 指南

## RTF 力场集成状态（2026-04-07）

### 成功标准
- ✅ **Build 阶段**：能够从 RTF 目录生成完整的 FEP 系统（bound/unbound legs）
- ✅ **Grompp 验证**：bound/unbound 的 `gmx grompp` 都能通过（maxwarn 10）
- ✅ **EM 收敛**：bound/unbound 的能量最小化都能收敛（Fmax < 200 kJ/mol/nm）
- ✅ **NVT 启动**：bound/unbound 的 NVT 都能稳定运行（2 ps smoke test 通过）

### 已修复的关键 Bug

**Bug 1: 混合拓扑键参数过度零化**
- **位置**: `prism/fep/gromacs/itp_builder.py:807-840`
- **问题**: 当 A/B 两态原子类型不同时，即使两态都有显式键参数，也会被错误地零化
- **影响**: 导致分子骨架在 lambda=0 时断裂，unbound leg 在 EM 中塌缩
- **修复**: 只在某一态缺少显式参数时才补零，保留两态都有的真实键参数
- **验证**: p38-19-24 系统 bound/unbound EM 都能正常收敛

**Bug 2: 坐标单位转换缺失**
- **位置**: `prism/fep/modeling/hybrid_package.py:746-791`
- **问题**: PDB/MOL2 文件使用 Å 单位，但代码未转换为 GROMACS 需要的 nm
- **影响**: 坐标放大 10 倍，导致"excluded atoms distance > 8 nm"错误
- **修复**: 在 `_parse_pdb_atoms()` 和 `_parse_mol2_atoms()` 中除以 10.0
- **验证**: hybrid.gro 中原子间距离恢复正常（< 2 nm）

### 测试系统
- **路径**: `tests/gxf/FEP/unit_test/p38-19-24/`
- **蛋白**: p38α MAPK (charmm36-jul2022)
- **配体**: 19 → 24 (RTF 格式)
- **输出**: `charmm36_jul2022-mut_rtf/GMX_PROLIG_FEP/`

### 验证命令
```bash
cd tests/gxf/FEP/unit_test/p38-19-24
python ../test_run_fep.py --forcefield charmm36-jul2022 --ligand-forcefield rtf

# 检查 EM 结果
grep "Maximum force" charmm36_jul2022-mut_rtf/GMX_PROLIG_FEP/bound/repeat1/build_em.log | tail -1
grep "Maximum force" charmm36_jul2022-mut_rtf/GMX_PROLIG_FEP/unbound/repeat1/build_em.log | tail -1

# 预期输出：
# bound: Fmax = 183.2 kJ/mol/nm (7778 steps)
# unbound: Fmax = 176.8 kJ/mol/nm (8690 steps)
```

### 当前限制
- RTF 输入必须是目录格式（包含 .rtf, .prm, .pdb 文件）
- 需要使用 `--ref-ligand input/19 --mut-ligand input/24` 指定目录路径
- 自动检测逻辑：扫描 `input/` 下包含 `*.rtf` 文件的子目录
