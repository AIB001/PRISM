---
name: fep-visualization
description: FEP hybrid topology visualization validation
type: fep
---

# FEP 可视化验证

## 核心职责

通过 HTML 可视化验证 PRISM FEP 模块的 atom mapping 结果。确保可视化正确反映 mapping 算法的分类结果，并使用 playwright MCP 检查真实渲染。

## 快速验证

**使用自动化验证脚本**：
```bash
# 运行完整验证
python .claude/skills/fep-visualization-check.py <mapping.html>

# 示例
python .claude/skills/fep-visualization-check.py oplsaa/GMX_PROLIG_FEP/common/hybrid/mapping.html
```

**验证脚本检查项目**：
- ✅ 灰色原子检测（应为0）
- ✅ Unknown atoms检测（应为0）
- ✅ 原子标签检查（无AtomXX占位符）
- ✅ 总电荷显示（应≈0）
- ✅ Common数量一致性（左右相等）
- ✅ Common原子一一对应（完全匹配）
- ✅ 键级渲染（应有aromatic/double键）
- ✅ 原子属性完整性

## 手动验证（可选）

如果需要手动检查特定项目：

### 1. 灰色原子检测
```bash
grep "rgb(200, 200, 200)" output/mapping.html  # 应该返回0
grep '"classification": "unknown"' output/mapping.html  # 应该返回0
```

### 2. Common数量一致性
```bash
grep -o "Common: [0-9]*" output/mapping.html | head -1
# 左右两边应该显示相同的数字
```

### 3. 键级检查
```python
import json, re
from collections import Counter

html = Path("mapping.html").read_text()
bonds_a = json.loads(re.search(r'const BONDS_A = (\[.*?\]);', html, re.DOTALL).group(1))
bonds_b = json.loads(re.search(r'const BONDS_B = (\[.*?\]);', html, re.DOTALL).group(1))

bond_types_a = Counter([b.get('typeName', b.get('type', 1)) for b in bonds_a])
bond_types_b = Counter([b.get('typeName', b.get('type', 1)) for b in bonds_b])

print(f"配体A键级分布: {dict(bond_types_a)}")
print(f"配体B键级分布: {dict(bond_types_b)}")
```

### 4. 浏览器验证（必须）
```bash
# 使用 playwright MCP 检查真实渲染
mcp__plugin_playwright_playwright__browser_navigate file:///path/to/mapping.html
mcp__plugin_playwright_playwright__browser_take_screenshot
```

## 关键检查点

### a. 灰色原子检测（✅ 已修复 2026-03-29）
- ❌ 不应该出现灰色原子（`rgb(200, 200, 200)`）
- ❌ 不应该出现 `classification: "unknown"`
- ✅ 所有原子必须有正确的分类（common/transformed/surrounding）

**根本原因**（已修复）：
- **Bug位置**: `prism/fep/visualize/molecule.py`
- **问题**: RDKit的`Chem.AddHs()`添加隐式氢原子
- **修复**: 移除`Chem.AddHs()`调用

### b. 原子标签
- ✅ 所有原子显示真实名字（如 "CD2", "H09A"）
- ❌ 不应该出现 "Atom15" 这种占位符
- ❌ 不应该只显示元素符号

**代码路径**: `prism/fep/visualize/molecule.py::prepare_mol_with_charges_and_labels()`

### c. 总电荷显示（**映射后**的电荷）
- **必须显示 mapping 修改后的电荷**
- REF 模式：State A 电荷，State B 被修改为 State A 电荷
- MUT 模式：State B 电荷，State A 被修改为 State B 电荷
- MEAN 模式：两边都是平均电荷
- 总电荷必须 ≈ 0（允许浮点误差）

**代码路径**: `prism/fep/visualize/html.py::_generate_canvas_html()`

### d. Common 数量一致性
- **左右两边 common 数量必须相等**
- Common 原子必须一一对应（相同原子名）

**快速检查**：
```bash
python .claude/skills/fep-visualization-check.py mapping.html | grep "Common原子对应"
```

### e. 键级（Bond Order）渲染
- ✅ 必须有键级数据（BONDS_A, BONDS_B）
- ✅ 芳香环键应显示为 AROMATIC 类型
- ✅ 双键应显示为 DOUBLE 类型
- ✅ 单键应显示为 SINGLE 类型

**代码路径**: `prism/fep/visualize/mapping.py::_assign_bond_metadata()`

### f. 按钮样式
- Coloring Mode 按钮必须**整个按钮区域高亮**
- 左右不能有空白（margin: 0）

**代码路径**: `prism/fep/visualize/templates/styles.css`

### g. 元素着色模式
- 碳：不能深黑（过重）
- 氢：不能浅灰（看不清）
- 白底白球：必须保留可见边框

**配色参考**：
- C: #909090（中灰）
- H: #E8E8E8（浅灰但可见）
- O: #FF0D0D（红）
- N: #3050F8（蓝）
- S: #FFFF30（橙），字黑色

## 特殊力场处理

**OpenFF/OPLS-AA**：
- 原子类型通用（output_0, output_1）
- **不能依赖 atom type 匹配**
- 优先级：距离 > 元素 > 电荷（忽略 type）

**GAFF/CGenFF**：
- 原子类型包含位置信息（ca, c3, ha, hc）
- 完整匹配：距离 > 元素 > **类型** > 电荷

**代码路径**: `prism/fep/core/mapping.py::DistanceAtomMapper._are_atoms_compatible()`

## 典型问题案例

### 案例 1：灰色原子问题（✅ 已修复 2026-03-29）
**症状**：HTML中某些原子显示为灰色，classification为"unknown"

**修复**：移除 `pdb_to_mol()` 中的 `Chem.AddHs()` 调用

### 案例 2：Common 数量不一致
**症状**：算法输出 Common=17，HTML 显示左边 17 个，右边 15 个

**排查**：
1. 检查 `prepare_mol_with_charges_and_labels()` 是否所有原子都匹配成功
2. 检查是否有原子显示为 "Atom15"
3. 检查坐标匹配阈值（0.6 Å）

### 案例 3：总电荷 ≠ 0（MUT 模式）
**症状**：MUT 模式总电荷显示 -0.000200

**原因**：使用了原始 GAFF 电荷，而不是 mapping 修改后的电荷

**修复**：在 mapping 之后重新计算总电荷

### 案例 4：原子显示为 "Atom15"
**症状**：HTML 中部分原子显示为 "Atom15", "Atom20"

**原因**：RDKit Mol 原子未能匹配到 hybrid topology atoms

**排查**：
1. 检查 mol2/SDF 坐标 vs gro 坐标
2. 检查元素是否匹配
3. 检查 base name 映射表

### 案例 5：按钮样式问题
**症状**：Coloring Mode 按钮只有文字变蓝，左右有空白

**修复**：CSS 中设置 `margin: 0`

## 测试数据位置

- GAFF 测试：`tests/gxf/FEP/unit_test/oMeEtPh-EtPh/gaff_test_output/`
- HTML 输出：`oMeEtPh-EtPh_gaff_{mode}.html`

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
