---
name: gaff-debug
description: GAFF 力场参数与电荷调试
type: fep
---

# GAFF 力场参数与电荷调试

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable


## 核心职责

调试 GAFF 力场参数生成和电荷分配问题，追踪电荷从原始 mol2 → GAFF 处理 → hybrid topology → mapping 后的完整流程，验证单配体 vs hybrid topology 的一致性。

## 关键检查点

### 0. 配置验证（最重要！）

**第一步：总是先验证配置是否正确读取**

```bash
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh

# 1. 查看 YAML 配置文件

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

cat fep_ref.yaml
# 应该看到：

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# mapping:

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

#   dist_cutoff: 1.0

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

#   charge_cutoff: 0.06

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

#   charge_common: ref

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable


# 2. 验证 FEPConfig 是否正确读取

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
print('从 YAML 读取的参数:')
for k, v in params.items():
    print(f'  {k}: {v}')
"

# 3. 对比不同模式的配置

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
    cat fep_${mode}.yaml | grep -A 5 "mapping:"
done
```

**代码路径**：
- 配置读取：`prism/fep/config.py::FEPConfig.get_mapping_params()` (line 73-84)
- 默认值定义：`prism/fep/config.py::FEPConfig.get_mapping_params()` (line 79)

### 1. 电荷来源追踪

**完整流程**：
```
原始 mol2 电荷
    ↓
Antechamber (GAFF 参数化)
    ↓
单配体 GAFF topology (.top)
    ↓
FEbuilder 处理 (hybrid.itp)
    ↓
DistanceAtomMapper 修改 (charge_common 模式)
    ↓
HTML 显示的电荷
```

**关键文件**：
- 原始输入：`tests/gxf/FEP/unit_test/oMeEtPh-EtPh/oMeEtPh.mol2`
- State A GAFF：`gaff_test_output/gaff_oMeEtPh/forcefield/BBE.amb2gmx/BBE_GMX.top`
- State B GAFF：`gaff_test_output/gaff_EtPh/forcefield/BBE.amb2gmx/BBE_GMX.top`
- Hybrid topology：`gaff_test_output/GMX_PROLIG_FEP/common/hybrid/hybrid.itp`
- Mapping 后：查看 HTML 或 `atoms_a/atoms_b` 对象

**检查命令**：
```bash
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh

# 1. 原始 mol2 电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep -A 20 "@<TRIPOS>ATOM" oMeEtPh.mol2 | awk '{print $9, $10}'

# 2. GAFF 处理后电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep "^\s*\(.*\)\s*[0-9]" gaff_test_output/gaff_oMeEtPh/forcefield/BBE.amb2gmx/BBE_GMX.top | awk '{print $1, $NF}'

# 3. Hybrid topology 电荷（state A 和 B）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep "qtot" gaff_test_output/GMX_PROLIG_FEP/common/hybrid/hybrid.itp

# 4. HTML 显示的电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep -o '"charge": [0-9.-]*' gaff_test_output/oMeEtPh-EtPh_gaff_ref.html
```

### 2. 单配体 vs Hybrid 对比

**关键差异**：
- **单配体 GAFF**：只包含一个 ligand 的所有原子（21 个 for oMeEtPh）
- **Hybrid topology**：包含 state A + state B 的所有原子（25 个 total）
  - State A 真实原子：21 个
  - State B 真实原子：18 个
  - Dummy 原子：7 个（state A 中 state B 独有的）

**对比方法**（使用独立工具脚本）：
```bash
cd tests/gxf/FEP/unit_test
python regenerate_gaff_html.py --all-modes
```

**脚本内部逻辑**：
1. 读取 hybrid top
ology 原子（`read_ligand_from_prism(..., state="a"/"b")`）
2. 读取单配体 GAFF topology（`read_gaff_topology_atoms()`）
3. 用真实 GAFF 参数覆盖 hybrid 中非 dummy 原子的 type/charge
4. 在 mapping **之后**重新计算总电荷（line 108-110）

**代码位置**：`tests/gxf/FEP/unit_test/regenerate_gaff_html.py:26-84` (read_gaff_topology_atoms + merge_hybrid_with_real_params)

**注意**：这是独立工具脚本，不应该被 import！

### 3. Mapping 后电荷修改

**Charge common 模式**：
- **REF 模式**：State A 电荷保持，State B 改为 State A 电荷
- **MUT 模式**：State B 电荷保持，State A 改为 State B 电荷
- **MEAN 模式**：两边都改为平均电荷 `(charge_a + charge_b) / 2`
- **NONE 模式**：不修改电荷

**配置来源**：从 `fep.yaml` 的 `mapping.charge_common` 读取（默认 "mean"）

**代码位置**：
- Mapping 算法：`prism/fep/core/mapping.py::DistanceAtomMapper.map()` (line 179+)
- 电荷修改：`prism/fep/core/mapping.py::DistanceAtomMapper.map()` (line 285-301)
```python
# Step 3.5: Apply charge_common strategy

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

if self.charge_common == "ref":
    for atom_a, atom_b in common:
        atom_b.charge = atom_a.charge  # State B 改为 State A
elif self.charge_common == "mut":
    for atom_a, atom_b in common:
        atom_a.charge = atom_b.charge  # State A 改为 State B
elif self.charge_common == "mean":
    for atom_a, atom_b in common:
        ave = (atom_a.charge + atom_b.charge) / 2.0
        atom_a.charge = ave
        atom_b.charge = ave
elif self.charge_common == "none":
    pass  # 不修改
```

**验证方法**：
```bash
# 检查 MUT 模式：State A 的 aromatic carbon 电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 应该显示 State B 的电荷（-0.128），而不是 State A 的（-0.0733）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep -A 2 '"name": "CD2"' gaff_test_output/oMeEtPh-EtPh_gaff_mut.html | grep "charge"

# 期望输出（MUT 模式）：

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# State A CD2 charge: -0.128  (来自 EtPh)

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# State B CD2 charge: -0.128  (保持 EtPh)

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

```

### 4. 总电荷验证

**必须为 0**（映射后）：
```python
# ✅ 正确做法：在 mapping 后重新计算

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

atoms_a = copy.deepcopy(base_atoms_a)
atoms_b = copy.deepcopy(base_atoms_b)
mapping = mapper.map(atoms_a, atoms_b)  # 修改电荷

# 重新计算总电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

total_charge_a = sum(atom.charge for atom in atoms_a)
total_charge_b = sum(atom.charge for atom in atoms_b)

print(f"State A: {total_charge_a:.6f}")  # 应该 ≈ 0
print(f"State B: {total_charge_b:.6f}")  # 应该 ≈ 0
```

**代码位置**：`tests/gxf/FEP/unit_test/regenerate_gaff_html.py:98-107`

**❌ 错误做法**：
```python
# 不要使用原始 GAFF 电荷！

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

total_charge_a = real_total_charge_a  # 这是未修改的电荷
total_charge_b = real_total_charge_b
```

### 5. 原子类型映射

**GAFF 原子类型**（用于 mapping 匹配）：
- `ca`: 芳香碳（aromatic carbon）
- `c3`: sp3 碳（aliphatic carbon）
- `ha`: 芳香氢（aromatic hydrogen）
- `hc`: 脂肪氢（aliphatic hydrogen）
- `oh`, `os`: 醇/醚氧
- `n`, `na`, `nb`: 各种氮

**检查方法**：
```bash
# 查看单配体 GAFF 原子类型

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep -E "^\s+[0-9]+\s+[A-Z]" gaff_test_output/gaff_oMeEtPh/forcefield/BBE.amb2gmx/BBE_GMX.top | head -20

# 查看 hybrid topology 原子类型

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep "; type" gaff_test_output/GMX_PROLIG_FEP/common/hybrid/hybrid.itp
```

**代码位置**：
- GAFF 生成：`prism/forcefield/gaff.py::GAFFForceFieldGenerator._generate_parameters()`
- Hybrid 构建：`prism/fep/modeling/hybrid_package.py::_create_hybrid_topology()`

## 典型问题案例

### 案例 1：MUT 模式总电荷 ≠ 0

**症状**：
```
State A total charge: -0.000200  (应该是 0)
State B total charge: -0.000100  (应该是 0)
```

**原因**：
- HTML 显示的是原始 GAFF 电荷，不是 mapping 修改后的电荷
- 代码使用了 `real_total_charge_a` 而不是 `sum(atom.charge for atom in atoms_a)`

**修复**：
```python
# tests/gxf/FEP/unit_test/regenerate_gaff_html.py

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 在 mapping 之后重新计算

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

total_charge_a = sum(atom.charge for atom in atoms_a)
total_charge_b = sum(atom.charge for atom in atoms_b)
```

**验证**：
```bash
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh
python tests/gxf/FEP/unit_test/regenerate_gaff_html.py --all-modes

# 检查所有模式的总电荷

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
    grep "molecule-label" gaff_test_output/oMeEtPh-EtPh_gaff_${mode}.html | grep -o "Charge: [0-9.-]*"
done
```

### 案例 2：State A 显示 State B 电荷（MUT 模式）

**症状**：
```
MUT 模式 HTML 显示：
  State A CD2: -0.1280  (这是 EtPh 的电荷！)
  State B CD2: -0.1280  (正确，EtPh 的电荷)
```

**分析**：
- MUT 模式：State A 应该用 State B 电荷（-0.128）✅
- 但这些电荷来自哪里？
  - 如果来自 `merge_hybrid_with_real_params()` → 混合了 GAFF 电荷
  - 如果来自 `DistanceAtomMapper.map()` → 正确修改

**排查**：
```python
# 检查 merge_hybrid_with_real_params 是否正确 overlay

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# hybrid atoms 先读取，然后 real GAFF params overlay

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

# 应该只 overlay 电荷和类型，不改变其他属性

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable


# 检查 mapper 是否真的修改了 charges

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

for atom_a, atom_b in mapping.common:
    print(f"Before: A={atoms_a[atom_a].charge:.4f}, B={atoms_b[atom_b].charge:.4f}")
    # mapper.map() 修改
    print(f"After: A={atoms_a[atom_a].charge:.4f}, B={atoms_b[atom_b].charge:.4f}")
```

**代码位置**：
- Merge：`tests/gxf/FEP/unit_test/regenerate_gaff_html.py:58-84` (merge_hybrid_with_real_params)
- Mapping：`prism/fep/core/mapping.py:285-301`

### 案例 3：Hybrid topology 电荷不匹配

**症状**：
```
State A 单配体 GAFF: CD2 charge = -0.0733
Hybrid topology (state A): CD2 charge = -0.0734  (略有差异)
```

**原因**：
- FEbuilder 进行了电荷重新分配（charge redistribution）
- 确保总电荷 = 0

**检查**：
```bash
# 对比 hybrid.itp 和单配体 GAFF

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep "CD2" gaff_test_output/GMX_PROLIG_FEP/common/hybrid/hybrid.itp
grep "CD2" gaff_test_output/gaff_oMeEtPh/forcefield/BBE.amb2gmx/BBE_GMX.top
```

**这是正常的**：
- FEbuilder 会调整电荷以确保 neutrality
- 差异通常 < 0.001 e

### 案例 4：原子类型缺失

**症状**：
```
GROMACS error: Atom type X not found in force field
```

**原因**：
- GAFF 参数化失败
- ITP 文件缺少 `[ atomtypes ]` 声明

**检查**：
```bash
# 检查 ITP 文件是否包含 atomtypes

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

grep "\[ atomtypes \]" gaff_test_output/gaff_oMeEtPh/forcefield/BBE.amb2gmx/BBE_GMX.top

# 检查 GAFF 是否成功运行

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

ls gaff_test_output/gaff_oMeEtPr/forcefield/BBE.amb2gmx/
# 应该看到：BBE_GMX.top, BBE_GMX.itp, BBE_GMX.pdb

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

```

**代码位置**：
- GAFF 生成：`prism/forcefield/gaff.py::GAFFForceFieldGenerator.generate()`

## 快速验证脚本

```bash
cd /data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/oMeEtPh-EtPh

# 1. 重新生成 HTML（使用最新代码）

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

python tests/gxf/FEP/unit_test/regenerate_gaff_html.py --all-modes

# 2. 检查所有模式的总电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

echo "=== Total Charges ==="
for mode in ref mut mean none; do
    echo "Mode: $mode"
    grep "molecule-label" gaff_test_output/oMeEtPh-EtPh_gaff_${mode}.html | \
        grep -oP "Charge: \K[-0-9.]+" | \
        awk '{sum+=$1} END {print "  Sum:", sum}'
done

# 3. 检查 MUT 模式 aromatic carbon 电荷

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

echo "=== MUT Mode Aromatic Carbons ==="
grep -o '"name": "C[DE].*", "charge": [-0-9.]*' \
    gaff_test_output/oMeEtPh-EtPh_gaff_mut.html | \
    head -10

# 4. 对比 hybrid vs 单配体 GAFF

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

echo "=== Hybrid vs Single Ligand ==="
echo "State A CD2:"
grep "CD2" gaff_test_output/GMX_PROLIG_FEP/common/hybrid/hybrid.itp | head -1
grep "CD2" gaff_test_output/gaff_oMeEtPh/forcefield/BBE.amb2gmx/BBE_GMX.top | head -1
```

## 测试数据结构

```
tests/gxf/FEP/unit_test/oMeEtPh-EtPh/
├── oMeEtPh.mol2                    # 原始输入
├── EtPh.mol2
├── gaff_test_output/
│   ├── gaff_oMeEtPh/               # State A 单配体 GAFF
│   │   └── forcefield/BBE.amb2gmx/
│   │       └── BBE_GMX.top         # GAFF 电荷（未修改）
│   ├── gaff_EtPh/                  # State B 单配体 GAFF
│   │   └── forcefield/BBE.amb2gmx/
│   │       └── BBE_GMX.top
│   ├── GMX_PROLIG_FEP/
│   │   └── common/hybrid/
│   │       └── hybrid.itp          # FEbuilder 处理后的 hybrid
│   └── oMeEtPh-EtPh_gaff_*.html    # 最终可视化
└── fep_*.yaml                       # 配置文件
```

## 必须避免的错误

1. **不要混淆电荷来源**：
   - 原始 mol2 电荷 ≠ GAFF 电荷
   - 单配体 GAFF 电荷 ≠ Hybrid 电荷
   - Mapping 前电荷 ≠ Mapping 后电荷

2. **不要在 mapping 前计算总电荷**：
   ```python
   # ❌ 错误
   total_charge = sum(real_charges)  # 原始电荷
   mapper.map(atoms_a, atoms_b)
   visualize(total_charge=total_charge)  # 显示原始电荷！

   # ✅ 正确
   mapper.map(atoms_a, atoms_b)
   total_charge = sum(atom.charge for atom in atoms_a)  # mapping 后
   visualize(total_charge=total_charge)
   ```

3. **不要忽略浮点误差**：
   - 总电荷 ±0.001 是正常的
   - 但 -0.000200 说明有问题

4. **不要假设 hybrid 和单配体电荷相同**：
   - FEbuilder 会重新分配电荷
   - 差异 < 0.001 是正常的

## 相关代码

- GAFF 生成：`prism/forcefield/gaff.py::GAFFForceFieldGenerator`
- Hybrid 构建：`prism/fep/modeling/hybrid_package.py::_create_hybrid_topology()`
- Mapping 算法：`prism/fep/core/mapping.py::DistanceAtomMapper.map()`
- 可视化生成：`prism/fep/visualize/html.py::visualize_mapping_html()`
- 测试脚本：`tests/gxf/FEP/unit_test/regenerate_gaff_html.py`
