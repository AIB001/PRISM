---
name: fep-mapping-viz
description: FEP hybrid topology atom mapping and visualization validation (joint inspection)
type: fep
---

# FEP Mapping 与可视化验证

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable


## 核心职责

检查 PRISM FEP 模块的 hybrid topology 构建和 atom mapping 算法，并通过 HTML 可视化相互验证结果一致性。确保 DistanceAtomMapper 算法与 FEbuilder 一致，且可视化正确反映分类结果。

## 关联性

**Mapping 算法 → 可视化**：
- DistanceAtomMapper 产生分类（common/transformed/surrounding）
- HTML 可视化展示这些分类
- **必须一致**：mapping 说 common，HTML 必须显示 common

**可视化 → Mapping 验证**：
- HTML 是算法结果的"眼睛"
- 通过可视化可快速发现算法错误
- common 数量左右不相等 → 算法有问题或显示有问题

## 关键检查点

### 1. Mapping 算法验证

**DistanceAtomMapper 7 步算法**（必须与 FEbuilder 一致）：
1. 距离筛选：`dist < dist_cutoff`
2. 元素匹配：`atom_a.element == atom_b.element`
3. 类型匹配：`atom_a.type == atom_b.type`（OpenFF/OPLS 除外）
4. 电荷匹配：`abs(atom_a.charge - atom_b.charge) < charge_cutoff`
5. 确定 common pairs
6. 确定 transformed（unique atoms）
7. 确定 surrounding（同位置但电荷/类型不同）

**配置来源优先级**（从高到低）：
1. **YAML 配置**：`fep.yaml` 的 `mapping` 部分
2. **FEPConfig 默认值**：`dist_cutoff=0.6, charge_cutoff=0.05, charge_common="mean"`
3. **DistanceAtomMapper 默认值**：与 FEPConfig 相同

**实际默认值**（代码中）：
- `dist_cutoff`: **0.6** Å（不是 1.0！）
- `charge_cutoff`: **0.05** e（不是 0.3！）
- `charge_common`: **"mean"**（不是 "ref"！）

**注意**：测试用例可能使用不同的值（如 `fep_ref.yaml` 中 `dist_cutoff=1.0`），但这是**测试配置**，不是代码默认值。

**代码路径**：
- 算法实现：`prism/fep/core/mapping.py::DistanceAtomMapper.map()` (line 179+)
- 配置读取：`prism/fep/config.py::FEPConfig.get_mapping_params()` (line 73-84)
- YAML 解析：`prism/fep/config.py::FEPConfig.__init__()` (line 37-56)
- FEbuilder 参考：`/home/gxf1212/data2/work/make_hybrid_top/FEbuilder/src/FEbuilder/setup/make_hybrid.py`

**检查命令**：
```bash
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh

# 1. 查看实际配置值（从 YAML）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

cat fep_ref.yaml
# 应该看到：dist_cutoff: 1.0, charge_cutoff: 0.06, charge_common: ref

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable


# 2. 验证配置读取

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

python -c "
from prism.fep.config import FEPConfig
cfg = FEPConfig('.', config_file='config_gaff.yaml', fep_file='fep_ref.yaml')
params = cfg.get_mapping_params()
print('Actual params:', params)
# 输出：dist_cutoff=1.0, charge_cutoff=0.06, charge_common='ref'

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

"

# 3. 运行 mapping 并检查分类数量

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

python tests/gxf/FEP/unit_test/test_gaff_fep_mapping.py

# 4. 检查日志输出

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 应该看到：

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# Common: 17, TransA: 4, TransB: 1, SurrA: 0, SurrB: 0  (REF 模式)

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# Common: 17, TransA: 4, TransB: 1, SurrA: 0, SurrB: 0  (MUT 模式)

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# Common: 17, TransA: 4, TransB: 1, SurrA: 0, SurrB: 0  (MEAN 模式)

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# Common: 4, TransA: 4, TransB: 1, SurrA: 13, SurrB: 13  (NONE 模式)

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

```

**常见问题**：
- ❌ "common 数量左右不一致" → 算法 bug 或坐标匹配问题
- ❌ "乙基末端甲基的 3 个氢对不上" → 检查 YAML 中的 `dist_cutoff` 值和坐标
- ❌ "C1/C2 本应 common 的侧链碳没被识别" → 检查原子类型匹配
- ❌ "为什么用了不同的 cutoff 值" → 检查是哪个 YAML 文件（fep_ref.yaml vs fep_mut.yaml）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 必须打开浏览器检查

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

mcp__plugin_playwright_playwright__browser_navigate http://localhost:8506/oMeEtPh-EtPh_gaff_ref.html
```

**HTML 关键检查项**：

#### a. 灰色原子检测（✅ 已修复 2026-03-29）
- ❌ 不应该出现灰色原子（`rgb(200, 200, 200)`）
- ❌ 不应该出现 `classification: "unknown"`
- ✅ 所有原子必须有正确的分类（common/transformed/surrounding）
- ✅ 所有原子必须有非零电荷（除非真的是中性原子）

**检查方法**：
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
# 数组长度应该等于 len(atoms_a)

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

```

**根本原因**（已修复）：
- **Bug位置**: `prism/fep/visualize/molecule.py`
- **问题**: RDKit的`Chem.AddHs()`添加隐式氢原子，创建额外原子（如"Atom32", "Atom33"）
- **修复**: 移除`Chem.AddHs()`调用，对于CHARMM-GUI文件直接使用PDB分子

**代码路径**：`prism/fep/visualize/molecule.py::pdb_to_mol()`, `prepare_mol_with_charges_and_labels()`

#### b. 原子标签（名字显示）
- ✅ 所有原子显示真实名字（如 "CD2", "H09A"）
- ❌ 不应该出现 "Atom15" 这种占位符
- ❌ 不应该只显示元素符号（如只显示 "C", "H"）

**代码路径**：`prism/fep/visualize/molecule.py::prepare_mol_with_charges_and_labels()`
- 坐标匹配：0.6 Å 阈值
- Fallback：base name 匹配（H07 → H07A）

#### b. 总电荷显示（**映射后**的电荷）
- **必须显示 mapping 修改后的电荷**，不是原始 GAFF 电荷
- REF 模式：State A 电荷，State B 被修改为 State A 电荷
- MUT 模式：State B 电荷，State A 被修改为 State B 电荷
- MEAN 模式：两边都是平均电荷
- 总电荷必须 ≈ 0（允许浮点误差）

**代码路径**：`prism/fep/visualize/html.py::_generate_canvas_html()`
```python
# 必须在 mapping 之后计算总电荷！

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

total_charge_a = sum(atom.charge for atom in atoms_a)  # ✅ 正确
total_charge_b = sum(atom.charge for atom in atoms_b)
```

**反面教材**（不要这样做）：
```python
# ❌ 错误：使用原始 GAFF 电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

total_charge_a = real_total_charge_a  # 这是未修改的电荷！
```

#### c. Common 数量一致性
- **左右两边 common 数量必须相等**
- HTML 底部 statistics 显示：`Common: 17`
- 两边的 `Atom Details & Statistics` 表格中 common 行数相等

**快速检查**：
```bash
# 从 HTML 提取 common 数量

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep -o "Common: [0-9]*" oMeEtPh-EtPh_gaff_ref.html | sort | uniq
# 应该只有一个数字，且两边相等

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

```

#### d. 按钮样式
- Coloring Mode 按钮必须**整个按钮区域高亮**，不能只有文字变蓝
- 左右不能有空白（margin: 0）

**代码路径**：`prism/fep/visualize/templates/styles.css`
```css
.toggle-buttons {
    margin: 0;  /* ✅ 正确：无左右空白 */
}
/* ❌ 错误：margin: 0 10px; 会有空白 */
```

#### e. 元素着色模式
- 碳：不能深黑（过重）
- 氢：不能浅灰（看不清）
- 白底白球：必须保留可见边框
- S（橙色）：字用黑色
- 后续 Br（暗红）、Cl（绿）、P（橙色）都要正确

**配色参考**：
- C: #909090（中灰）
- H: #E8E8E8（浅灰但可见）
- O: #FF0D0D（红）
- N: #3050F8（蓝）
- S: #FFFF30（橙），字黑色

### 3. 相互验证流程

**Step 1**: 生成 mapping 结果
```bash
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh
python tests/gxf/FEP/unit_test/regenerate_gaff_html.py --all-modes
```

**Step 2**: 检查算法输出
```bash
# 查看日志，确认分类数量

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 期望输出（REF/MEAN/MUT）：

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# Common: 17, TransA: 4, TransB: 1, SurrA: 0, SurrB: 0

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

```

**Step 3**: 用 playwright 检查 HTML
```bash
# 打开浏览器

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

mcp__plugin_playwright_playwright__browser_navigate file:///path/to/oMeEtPh-EtPh_gaff_ref.html

# 截图验证

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

mcp__plugin_playwright_playwright__browser_take_screenshot
```

**Step 4**: 交叉验证
- [ ] HTML 显示的 common 数量 = 算法输出的 common 数量
- [ ] HTML 左右 common 数量相等
- [ ] 总电荷 ≈ 0（mapping 修改后）
- [ ] 没有 "Atom15" 占位符
- [ ] 按钮样式正确（整个区域高亮）

### 4. 特殊力场处理

**OpenFF/OPLS-AA**：
- 原子类型通用（output_0, output_1），不包含位置信息
- **不能依赖 atom type 匹配**
- 优先级：距离 > 元素 > 电荷（忽略 type）

**GAFF/CGenFF**：
- 原子类型包含位置信息（ca, c3, ha, hc）
- 完整匹配：距离 > 元素 > **类型** > 电荷

**代码路径**：`prism/fep/core/mapping.py::DistanceAtomMapper._are_atoms_compatible()`
```python
# OpenFF/OPLS 特殊处理

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

if self.ignore_atom_type:  # 对 OpenFF/OPLS 设为 True
    # 跳过类型检查
else:
    if atom_a.type != atom_b.type:
        return False
```

## 典型问题案例

### 案例 1：灰色原子问题（✅ 已修复 2026-03-29）
**症状**：
- HTML中某些原子显示为灰色（`rgb(200, 200, 200)`）
- 这些原子的 `classification: "unknown"`, `charge: 0.0`
- 原子名称为 "Atom32", "Atom33" 等RDKit生成的名字

**根本原因**：
- **Bug位置**: `prism/fep/visualize/molecule.py`
- **问题1**: `pdb_to_mol()` 中调用 `Chem.AddHs()` 添加隐式氢原子
- **问题2**: `prepare_mol_with_charges_and_labels()` fallback路径中SMILES往返转换后再次调用 `Chem.AddHs()`
- **影响**: CHARMM-GUI的PDB文件已经包含所有显式氢原子，RDKit添加的隐式氢原子创建了额外原子，这些原子无法匹配到mapping数据

**修复**：
1. 移除 `pdb_to_mol()` 中的 `Chem.AddHs()` 调用
2. 移除 fallback 路径中的 `Chem.AddHs()` 调用
3. 对于CHARMM-GUI文件，跳过SMILES往返转换，直接使用PDB分子以保留原子名称

**验证**：
```bash
cd tests/gxf/FEP/unit_test
python test_42_38_gray_atom.py  # 简单系统（1个突变点）
python test_39_8_pdb_only.py    # 复杂系统（2个突变点）

# 检查HTML

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep "rgb(200, 200, 200)" output/42-38_regression.html  # 应该返回0
grep '"classification": "unknown"' output/39-8_pdb_only.html  # 应该返回0
```

**结果**：
- ✅ 42-38系统：0个灰色原子，38/38原子正确分类
- ✅ 39-8系统：0个灰色原子，38/38原子正确分类

### 案例 2：Common 数量不一致
**症状**：
- 算法输出：Common=17
- HTML 显示：左边 17 个，右边 15 个

**排查**：
1. 检查 `prepare_mol_with_charges_and_labels()` 是否所有原子都匹配成功
2. 检查是否有原子显示为 "Atom15"
3. 检查坐标匹配阈值（0.6 Å）

**代码位置**：`prism/fep/visualize/molecule.py:256-294`

### 案例 3：总电荷 ≠ 0（MUT 模式）
**症状**：
- MUT 模式总电荷显示 -0.000200

**原因**：
- 使用了原始 GAFF 电荷，而不是 mapping 修改后的电荷
- 代码：`total_charge_a = real_total_charge_a` ❌

**修复**：
```python
# 在 mapping 之后重新计算

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

total_charge_a = sum(atom.charge for atom in atoms_a)  # ✅
```

**代码位置**：`tests/gxf/FEP/unit_test/regenerate_gaff_html.py:98-107`

### 案例 4：原子显示为 "Atom15"
**症状**：
- HTML 中部分原子显示为 "Atom15", "Atom20" 而不是真实名字

**原因**：
- RDKit Mol 原子未能匹配到 hybrid topology atoms
- 坐标距离超过 0.6 Å 阈值
- base name fallback 也失败

**排查**：
1. 检查 mol2/SDF 坐标 vs gro 坐标
2. 检查元素是否匹配
3. 检查 base name 映射表

**代码位置**：`prism/fep/visualize/molecule.py:256-336`

### 案例 5：按钮样式问题
**症状**：
- Coloring Mode: element 按钮只有文字变蓝，左右有空白

**原因**：
- CSS 中 `margin: 0 10px` 导致左右 10px 空白

**修复**：
```css
.toggle-buttons {
    margin: 0;  /* 移除左右空白 */
}
```

**代码位置**：`prism/fep/visualize/templates/styles.css:447`

## 测试数据位置

- GAFF 测试：`tests/gxf/FEP/unit_test/oMeEtPh-EtPh/gaff_test_output/`
- 单配体 GAFF：`gaff_oMeEtPh/forcefield/BBE.amb2gmx/BBE_GMX.top`
- FEbuilder 处理：`common/hybrid/hybrid.itp`
- HTML 输出：`oMeEtPh-EtPh_gaff_{mode}.html`

## 快速验证脚本

```bash
cd /data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/oMeEtPh-EtPh

# 1. 重新生成所有模式的 HTML

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

python tests/gxf/FEP/unit_test/regenerate_gaff_html.py --all-modes

# 2. 检查总电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

for mode in ref mut mean none; do
    echo "=== $mode ==="
    grep "molecule-label molecule-label-a" gaff_test_output/oMeEtPh-EtPh_gaff_${mode}.html
    grep "molecule-label molecule-label-b" gaff_test_output/oMeEtPh-EtPh_gaff_${mode}.html
done

# 3. 检查 common 数量

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

for mode in ref mut mean none; do
    echo "=== $mode ==="
    grep -o "Common: [0-9]*" gaff_test_output/oMeEtPh-EtPh_gaff_${mode}.html | head -1
done

# 4. 用 playwright 检查（需要手动）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# mcp__plugin_playwright_playwright__browser_navigate file:///...

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

```

## 必须避免的错误

1. **不要只看本地文件**：必须用 playwright MCP 打开 HTML 检查真实渲染
2. **不要用原始电荷显示总电荷**：必须在 mapping 后重新计算
3. **不要忽略坐标匹配**：mol2 vs gro 坐标差异必须检查
4. **不要硬编码原子名字**：应该从 hybrid topology 读取
5. **不要混合不同力场的逻辑**：OpenFF 不能用类型匹配，GAFF 必须用

## 相关文档

- FEbuilder 算法：`/home/gxf1212/data2/work/make_hybrid_top/FEbuilder/src/FEbuilder/setup/make_hybrid.py`
- FEP 需求：`tests/gxf/FEP/plan/REQUIREMENTS.md`
- 可视化模板：`prism/fep/visualize/templates/`
- Mapping 配置：`tests/gxf/FEP/unit_test/oMeEtPh-EtPh/fep_*.yaml`
