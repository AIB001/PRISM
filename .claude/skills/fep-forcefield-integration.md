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

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

prism protein.pdb oMeEtPh.mol2 -o output_gaff \
  -lff gaff \
  --fep \
  --mutant EtPh.mol2

# OpenFF

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

prism protein.pdb oMeEtPh.mol2 -o output_openff \
  -lff openff \
  --fep \
  --mutant EtPh.mol2

# OPLS-AA

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

prism protein.pdb oMeEtPh.mol2 -o output_opls \
  -lff opls \
  --fep \
  --mutant EtPh.mol2
```

**预期输出**：
```
output_<ff>/
├── reference/
│   └── LIG.<ff>2gmx/      # 参考配体的力场
├── mutant/
│   └── LIG.<ff>2gmx/      # 突变配体的力场
└── GMX_PROLIG_FEP/        # FEP 系统文件
    ├── bound/
    ├── unbound/
    └── common/
```

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

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep "rgb(200, 200, 200)" output/mapping.html  # 应该返回0

# 2. 检查unknown classification

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep '"classification": "unknown"' output/mapping.html  # 应该返回0

# 3. 检查原子数量

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep -o "const ATOMS_A = \[" output/mapping.html
# 然后数组长度应该等于 len(atoms_a)

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable


# 4. 检查原子属性完整性

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep -A3 "const ATOMS_A" output/mapping.html | head -10
# 每个原子应该有：

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# - name: 原子名称（如 "N", "C1", "H1"）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# - element: 元素符号（如 "N", "C", "H"）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# - charge: 非零电荷（如 -0.463, 0.07）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# - type: 原子类型（如 "NG1T1", "CG3C51"）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# - classification: "common" 或 "transformed" 或 "surrounding"

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# - fepColor: 颜色值（不应该有 rgb(200, 200, 200)）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable


# 5. 示例：正确格式的原子数据

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

{"id": "A1", "name": "N", "element": "N", "x": 217.73, "y": 63.84,
 "radius": 14, "charge": -0.463, "type": "NG1T1",
 "classification": "common", "fepColor": "rgb(204, 229, 77)",
 "elementColor": "#3050F8"}
```

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

## 现有测试文件状态（2026-03-29 更新）

| 测试文件 | 力场/系统 | 测试数据 | 状态 | 说明 |
|----------|----------|----------|------|------|
| `examples/test_charmm_gui_generator.py` | CHARMM-GUI | 39-8 系统 | ✅ **已验证** | **HTML 生成成功** (74KB) |
| `tests/test_openff_fep_mapping.py` | OpenFF | oMeEtPh-EtPh | ⚠️ 未验证 | 需要实际运行测试 |
| `examples/test_e2e_prism_fep.py` | CGenFF (网站) | 39-8 系统 | ⚠️ 被跳过 | 缺少 CGenFF 网站文件 (*_gmx.pdb/top) |
| `examples/test_cgenff_compatibility.py` | - | 39-8, 25-36 | ⚠️ 未验证 | 只检测文件格式，未测试 FEP 集成 |

**✅ 已验证成功**：
- CHARMM-GUI → PRISM → FEP 映射 → HTML 可视化
- 测试数据：`39-8/output/39-8_charmm_gui_prism.html` (74KB)
- 映射结果：25 common, 13+13 transformed

## 当前真实情况（2026-03-29）

### ❌ 所有测试都未通过！

**测试代码问题**：
1. `test_charmm_gui_generator.py` - 缺少 `resolve_fep_case_dir` 函数，无法运行
2. `test_e2e_prism_fep.py` - 被 pytest.skip 跳过（缺少文件）
3. `test_openff_fep_mapping.py` - 未实际运行验证

**测试数据状态**：
- ✅ CHARMM-GUI 数据存在（`39-8/39/gromacs/LIG.itp`）
- ✅ OpenFF 测试数据存在（`oMeEtPh-EtPh/*.mol2`）
- ❌ CGenFF 网站数据缺失（`*_gmx.pdb`, `*_gmx.top`）

### 🎯 下一步行动（优先级排序）

#### 立即可做（今天）

**1. 修复 CHARMM-GUI 测试**（推荐起点）
```bash
# 问题：缺少 resolve_fep_case_dir 函数

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 解决：直接使用硬编码路径

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 从项目根目录运行: cd /path/to/prism/
# 修复 test_charmm_gui_generator.py 中的导入错误

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 运行测试验证 HTML 生成

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

```

**2. 运行 OpenFF 测试**
```bash
cd /data2/gxf1212/work/PRISM
pytest tests/test_openff_fep_mapping.py::TestOpenFFWithFEPMapping::test_end_to_end_fep_mapping -v -s
# 验证是否能生成 HTML

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

```

#### 需要额外工作

**3. 创建 GAFF 集成测试**
- 使用 `oMeEtPh-EtPh` 测试数据
- 复用 CHARMM-GUI 测试模板
- 验证 GAFF → PRISM → FEP 流程

**4. 获取 CGenFF 网站文件**（可选）
- 访问 https://cgenff.paramchem.org/
- 上传 39 和 8 的配体结构
- 下载 `*_gmx.pdb` 和 `*_gmx.top`

## 快速验证脚本（实际可运行）

```bash
#!/bin/bash
# 快速验证 FEP + 力场集成的准备工作

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable


cd /data2/gxf1212/work/PRISM

echo "=== FEP + 力场集成现状检查 ==="

# 1. 检查力场生成器

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

echo "1. 检查力场生成器..."
python -c "
from prism.forcefield import list_available_generators
gens = list_available_generators()
for name in gens:
    print(f'  ✓ {name}')
" 2>&1 | grep -v "Warning\|PyMBAR\|JAX"

# 2. 检查测试数据

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

echo "2. 检查测试数据..."
echo "  CHARMM-GUI (39-8):"
ls 39-8/39/gromacs/LIG.itp 2>/dev/null && echo "    ✓ Ligand 39: LIG.itp 存在" || echo "    ✗ Ligand 39: 缺失"
ls 39-8/8/gromacs/LIG.itp 2>/dev/null && echo "    ✓ Ligand 8: LIG.itp 存在" || echo "    ✗ Ligand 8: 缺失"

echo "  OpenFF (oMeEtPh-EtPh):"
ls oMeEtPh-EtPh/oMeEtPh.mol2 2>/dev/null && echo "    ✓ oMeEtPh.mol2 存在" || echo "    ✗ oMeEtPh.mol2 缺失"
ls oMeEtPh-EtPh/EtPh.mol2 2>/dev/null && echo "    ✓ EtPh.mol2 存在" || echo "    ✗ EtPh.mol2 缺失"

# 3. 尝试运行测试（会失败，但能看到错误）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

echo "3. 尝试运行测试..."
echo "  OpenFF 测试:"
pytest tests/test_openff_fep_mapping.py::TestOpenFFWithFEPMapping::test_end_to_end_fep_mapping -v 2>&1 | head -5

echo "  CHARMM-GUI 测试:"
python examples/test_charmm_gui_generator.py 2>&1 | head -5

echo "=== 检查完成 ==="
```

## 常见问题

### Q1: 某个力场生成失败
**原因**：依赖未安装或网络问题
**解决**：
```bash
# GAFF/GAFF2

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

mamba install -c conda-forge ambertools

# OpenFF

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

mamba install -c conda-forge openff-toolkit openff-interchange

# 检查服务器连接

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

ping ligpargen.org
ping swissparam.ch
```

### Q2: HTML中出现灰色原子
**症状**：FEP classification模式下某些原子显示为灰色（未着色），classification为"unknown"

**根本原因**：
- **Bug位置**：`prism/fep/visualize/molecule.py`
- **问题**：RDKit的`Chem.AddHs()`添加隐式氢原子，创建了额外的原子（如"Atom32", "Atom33"）
- **影响**：这些额外原子的名字无法匹配到原始原子，导致匹配失败 → 灰色显示

**✅ 已修复** (2026-03-29):
1. 移除`pdb_to_mol()`中的`Chem.AddHs()`调用
2. 移除`prepare_mol_with_charges_and_labels()`fallback路径中的`Chem.AddHs()`调用
3. 对于CHARMM-GUI文件，跳过SMILES往返转换，直接使用PDB分子

**验证方法**：
```bash
# 运行测试确认无灰色原子

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 从项目根目录运行: cd /path/to/prism/
python test_39_8_pdb_only.py  # 39-8系统（复杂，2个突变点）
python test_42_38_gray_atom.py  # 42-38系统（简单，1个突变点）

# 检查HTML

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep "rgb(200, 200, 200)" output/39-8_pdb_only.html  # 应该返回0
grep '"classification": "unknown"' output/42-38_regression.html  # 应该返回0
```

**如果仍然出现灰色原子**（极少数情况）：
```python
# 检查映射覆盖率

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

total_a = len(mapping.common) + len(mapping.transformed_a) + len(mapping.surrounding_a)
total_b = len(mapping.common) + len(mapping.transformed_b) + len(mapping.surrounding_b)

if total_a != len(atoms_a):
    print(f"⚠️  配体A有 {len(atoms_a) - total_a} 个原子未分类")

if total_b != len(atoms_b):
    print(f"⚠️  配体B有 {len(atoms_b) - total_b} 个原子未分类")

# 检查电荷为0的原子

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

zero_charge = [a for a in atoms_a + atoms_b if abs(a.charge) < 0.001]
print(f"电荷为0的原子: {len(zero_charge)}")
```

**其他解决方案**：
- 检查输入文件（ITP/GRO）是否正确
- 确认原子类型和电荷参数完整
- 尝试调整 `dist_cutoff` 和 `charge_cutoff`

### Q3: PNG可视化失败
**症状**：`ValueError: Depict error: Substructure match with reference not found`

**根本原因**：
- RDKit无法找到两个分子的公共子结构
- 分子结构差异太大（如多突变点）
- MOL2/PDB文件的原子名称不匹配

**解决方案（已修复）**：
- 代码已添加降级处理（`prism/fep/visualize/mapping.py`）
- 当2D对齐失败时，使用默认坐标生成PNG
- 显示warning但不影响HTML生成

**预防措施**：
- 优先使用简单系统（如42-38单突变点）测试
- 检查MOL2和PDB文件的原子名称是否一致
- 如果PNG不是必需的，可以跳过，只生成HTML

### Q4: FEP mapping 结果不对
**原因**：原子类型不匹配或坐标问题
**调试**：
```python
# 检查原子类型

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

for atom in atoms_a:
    print(f"{atom.name}: type={atom.atom_type}, elem={atom.element}")

# 检查坐标

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

for atom in atoms_a[:3]:
    print(f"{atom.name}: coord={atom.coord}")

# 尝试不同的 cutoff

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

mapper = DistanceAtomMapper(dist_cutoff=1.0)  # 放宽阈值
```

### Q5: CLI `--fep` 模式不工作
**原因**：可能是 builder.py 未正确处理参数
**检查**：
```bash
# 查看帮助

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

prism --help | grep -A 5 fep

# 检查是否生成了正确的目录结构

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

ls -la output_*/reference/ output_*/mutant/
```

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
