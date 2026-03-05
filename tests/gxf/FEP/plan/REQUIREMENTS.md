# PRISM-FEP 需求文档

## 1. 项目概述

### 1.1 目标
将 FEbuilder 的 FEP (Free Energy Perturbation) 双拓扑构建能力移植到 PRISM，支持 **GROMACS 的 A/B 状态拓扑**（`typeB/chargeB`）格式。

### 1.2 背景
- **FEbuilder**: 现有工具，支持 NAMD/CHARMM 格式的 FEP 系统构建
- **PRISM**: 支持 GROMACS 的蛋白-配体系统构建工具
- **需求**: 统一工作流，使用 GROMACS 进行 FEP 计算

### 1.3 核心价值
1. 复用 PRISM 现有的力场参数化流程 (GAFF, GAFF2, OpenFF, **CGenFF**, OPLS-AA)
2. 统一的分析工具链 (MD + FEP + PMF + REST2)
3. 更广泛的用户群体 (GROMACS 用户比 NAMD/CHARMM 更多)

### 1.4 主要功能
- **双拓扑构建**: 输入两个配体，输出 GROMACS A/B 双状态拓扑 (hybrid.itp)
- **原子映射**: 基于距离匹配，自动识别 common/transformed/surrounding 原子
- **力场支持**: 继承 PRISM: GAFF, GAFF2, OpenFF, CGenFF, OPLS-AA
- **系统构建**: 与蛋白质复合物组装、溶剂化、离子添加
- **FEP-MDP**: 自动生成 GROMACS FEP 计算所需的 MDP 文件
- **分析工具**: 解析 dhdl.xvg，计算 ΔΔG (支持 BAR/MBAR)

### 1.5 参考基线（已跑通体系）
- **位置**: `tests/gxf/FEP/ref`
- **拓扑来源**: pmx 生成的蛋白突变拓扑（示例 `system/newtop.top`）
- **力场**: `/home/gxf1212/data/local_programs/mutff/amber14sbmut.ff`
- **关键格式**: `[ atoms ]` 含 `typeB/chargeB/massB`，dummy 类型来自力场（`DUM_*`），不是固定 `dum`
- **运行模式**: 支持 λ-REMD (`gmx mdrun -multidir -replex`) 与非交换多窗口
- **分析基线**: `calculate.sh` 收集 `dhdl.xvg` 并 `gmx bar` 输出 `summary.dat`

---

## 2. 模块结构

```
prism/
├── fep/
│   ├── core/
│   │   ├── mapping.py             # 原子映射 (核心算法)
│   │   └── dual_topology.py       # 双拓扑构建 (GROMACS格式)
│   ├── gromacs/
│   │   ├── itp_builder.py         # ITP 文件生成
│   │   └── mdp_templates.py       # FEP-MDP 模板
│   ├── analysis/
│   │   ├── xvg_parser.py          # XVG 文件解析
│   │   └── estimators.py          # BAR/MBAR 分析
│   └── utils.py
├── builder/
│   ├── fep_builder.py             # FEP 构建器
│   └── cli.py                     # CLI 扩展
└── tests/
    └── fep/
        ├── test_mapping.py
        ├── test_dual_topology.py
        └── test_itp_builder.py
```

---

## 3. CLI 接口

```bash
# 基本用法
prism fep \
  --reference ligand.mol2 \
  --mutant mutant.mol2 \
  --receptor protein.pdb \
  --force-field gaff \
  --output fep_system

# 使用 CGenFF
prism fep \
  --reference ligand.mol2 \
  --mutant mutant.mol2 \
  --force-field cgenff \
  --forcefield-path toppar/ \
  --receptor protein.pdb \
  --protein-force-field charmm36 \
  --output fep_system

# 分析结果
prism fep-analyze --xvg dhdl.xvg --method bar
```

---

## 4. GROMACS vs NAMD 的关键差异

### 4.1 参数处理

**NAMD (FEbuilder)**:
- 需要参数合并 (merge_prm.py)
- 两个配体的键/角/二面角参数需要去重和合并
- 复杂的冲突处理逻辑

**GROMACS**:
- **不需要参数合并**
- 使用 typeB/chargeB 列，每个原子有两套状态
- bonds/angles/dihedrals 只需确保原子索引正确
- GROMACS 自动根据原子的 type/charge 选择状态

### 4.2 拓扑文件格式

**NAMD**: .rtf (topology) + .prm (parameters)，需要合并
**GROMACS**: .itp (include topology)，typeB/chargeB 格式，更简洁

### 4.3 实现影响

PRISM-FEP **不需要实现** FEbuilder 的 `merge_prm.py` 逻辑：
- ❌ 不需要复杂的参数值合并
- ✅ 只需要将原子索引映射到 hybrid_atoms
- ✅ GROMACS 会自动处理 A/B 状态的参数选择

---

## 5. 力场支持

### 5.1 支持的力场

| 力场 | CLI 参数 | 说明 | 示例文件 |
|------|----------|------|----------|
| GAFF | `--force-field gaff` | 通用小分子力场 | ligand.mol2 |
| GAFF2 | `--force-field gaff2` | GAFF 改进版 | ligand.mol2 |
| OpenFF | `--force-field openff` | 开放力场，基于 SMIRKS | ligand.sdf |
| **CGenFF** | `--force-field cgenff` | CHARMM 通用力场 | ligand.rtf + ligand.prm |
| OPLS-AA | `--force-field oplsaa` | OPLS 全原子力场 | ligand.mol2 |
| **amber14sbmut** | 外部力场 | pmx 蛋白突变拓扑用 | `amber14sbmut.ff` |

### 5.2 CGenFF 支持

**特点**:
- PRISM 已有 `prism/forcefield/cgenff.py` 模块
- 支持 CHARMM-GUI 输出格式 (`.rtf` + `.prm`)
- 适合含药物分子的系统

**测试数据**:
- 位置: `tests/gxf/FEP/test/hif2a/42-38/`
- 格式: CHARMM-GUI 输出
- 包含: `toppar/lig.rtf`, `toppar/lig.prm`

---

## 6. 详细设计文档

- **核心算法设计**: 见 `CORE_DESIGN.md`
  - 原子映射算法 (distance_based)
  - 双拓扑构建逻辑
  - HybridAtom 数据结构

- **GROMACS 输出格式**: 见 `OUTPUT_FORMAT.md`
  - ITP 文件格式 (typeB/chargeB)
  - MDP 模板文件
  - 虚拟原子处理

- **分析工具**: 见 `ANALYSIS.md`
  - XVG 文件解析
  - BAR/MBAR 方法
  - 自由能计算

---

## 7. 测试计划

### 7.1 测试数据
- **位置**: `tests/gxf/FEP/test/`
- **格式**: config.conf + 配体文件 + 受体
- **支持**: GAFF/OpenFF/CGenFF 格式
- **示例**: `hif2a/42-38/` (CHARMM-GUI 输出)
- **参考基线**: `tests/gxf/FEP/ref`（pmx 蛋白突变拓扑 + 可运行脚本）

### 7.2 测试脚本
```bash
# 批量测试
tests/gxf/FEP/scripts/run_fep_tests.sh

# 快速测试
tests/gxf/FEP/scripts/quick_test.sh <test_dir>
```

### 7.3 预期问题
| 问题 | 解决方案 |
|------|----------|
| CHARMM 格式 | 使用 PRISM 重新参数化或直接读取 .rtf/.prm |
| 虚拟原子 | 使用 dummy 类型（由力场提供，amber14sbmut 中为 `DUM_*`） |
| 力场兼容 | 统一使用一种力场，记录版本 |

---

## 9. 测试自动化详解

### 9.1 测试脚本设计

**批量测试脚本**: `tests/gxf/FEP/scripts/run_fep_tests.sh`

**设计思路**（伪代码）:
```bash
#!/bin/bash
# PRISM-FEP 批量测试脚本

# 1. 遍历所有测试案例
for test_dir in tests/gxf/FEP/test/*; do
    # 2. 解析 config.conf
    ref=$(grep "^ref=" "$test_dir/config.conf" | cut -d'=' -f2)
    mut=$(grep "^mut=" "$test_dir/config.conf" | cut -d'=' -f2)
    
    # 3. 查找 MOL2 文件
    ref_mol2=$(find "$test_dir" -name "*${ref}*.mol2")
    mut_mol2=$(find "$test_dir" -name "*${mut}*.mol2")
    
    # 4. 运行 PRISM-FEP
    prism fep \
        --reference "$ref_mol2" \
        --mutant "$mut_mol2" \
        --force-field gaff \
        --output "output/$(basename $test_dir)"
    
    # 5. 验证生成的拓扑
    cd "output/$(basename $test_dir)/complex/GMX_PROLIG_FEP"
    gmx grompp -f minim.mdp -c conf.gro -p topol.top -o test.tpr
    
    # 6. 记录结果
    if [ $? -eq 0 ]; then
        echo "✓ $test_dir"
    else
        echo "✗ $test_dir"
    fi
done

# 7. 生成汇总报告
generate_summary_report
```

**快速验证脚本**: `tests/gxf/FEP/scripts/quick_test.sh`
```bash
#!/bin/bash
# 快速验证单个测试案例

TEST_DIR="${1:?}"

# 解析配置并运行
prism fep \
    --reference "$(find $TEST_DIR -name "*.mol2" | head -1)" \
    --mutant "$(find $TEST_DIR -name "*.mol2" | tail -1)" \
    --force-field gaff \
    --output "test_output"
```

### 9.2 测试报告格式

**控制台输出**:
```
================================================================================
                    PRISM-FEP 批量测试报告
================================================================================
测试日期: 2025-01-15 14:30:25
测试目录: tests/gxf/FEP/test

[✓] RdRp/1p-cn-cch     remtp → remtp-1p-cn-cch
[✗] RdRp/2p-oh-och2ch3 remtp → remtp-2p-oh-och2ch3    (GROMACS验证失败)
[?] P38/42-38          lig42 → lig38                   (参数化超时)
[✓] hif2a/5-1         ref5 → mut1

测试完成: 4/4
成功: 2  失败: 1  超时: 1
详细报告: tests/gxf/FEP/output/test_report_20250115_143025.txt
================================================================================
```

**文本报告**: `tests/gxf/FEP/output/test_report.txt`
```text
================================================================================
测试案例: RdRp/1p-cn-cch
参考配体: remtp
突变配体: remtp-1p-cn-cch

结果: ✓ 通过
步骤:
  [1] 参数化: ✓ 成功 (耗时: 12.3s)
      - GAFF 力场参数生成完成
      - 生成的杂化拓扑包含 62 个原子

  [2] 系统构建: ✓ 成功 (耗时: 8.5s)
      - 蛋白-配体复合物构建完成
      - 加溶剂: TIP3P, 12.0 Å 盒子
      - 加离子: 0.15 M NaCl, 中和系统

  [3] 拓扑验证: ✓ 成功 (耗时: 2.1s)
      - GROMACS grompp 验证通过
      - 无警告或错误

输出目录: output/RdRp/1p-cn-cch/GMX_PROLIG_FEP
================================================================================
```

### 9.3 错误诊断

**常见错误及解决方案**:

| 错误信息 | 原因 | 解决方案 |
|---------|------|----------|
| Atom type not found | 力场参数缺失 | 检查 GAFF/OpenFF 参数化结果 |
| masses do not add up | dummy 类型或 massB 缺失 | 确保 dummy 类型来自 forcefield，且 `[ atoms ]` 含 `massB` |
| Inconsistent atom indices | 参数映射错误 | 检查 hybrid_atoms 索引映射 |
| grompp failed | 拓扑文件格式错误 | 验证 ITP 文件格式 |

---

## 10. 可视化报告

### 10.1 HTML 报告生成

**模块**: `prism/fep/analysis/report.py`

**功能**:
- 读取 FEP 计算结果（dhdl.xvg, 自由能值）
- 生成 HTML 报告，包含：
  - 测试案例信息和元数据
  - 原子映射可视化（2D 结构图）
  - 自由能变化曲线（ΔG vs λ）
  - λ 窗口能量分布直方图
  - 错误估计和置信区间
  - 与参考值对比（如果有）

**设计思路**:
1. 使用 Plotly.js 生成交互式图表
2. 使用 RDKit 生成 2D 分子结构图
3. 高亮显示 transformed atoms
4. 支持多个 λ 窗口的对比展示

**CLI 接口**:
```bash
prism fep-report --output-dir fep_system --report-type html
```

**输出文件**:
- `fep_report.html` - 主报告文件
- `plots/` - 图表数据目录
- `structures/` - 分子结构图

### 10.2 可视化元素

**1. 自由能变化曲线**:
- X 轴: λ (0.0 → 1.0)
- Y 轴: ΔG (kJ/mol)
- 显示各窗口的能量变化趋势
- 标注误差范围

**2. 原子映射图**:
- 使用 RDKit 生成 2D 分子结构
- 高亮显示 transformed atoms（红色）
- 显示 surrounding atoms（蓝色）
- 显示 common atoms（绿色）

**3. 能量分布**:
- 各 λ 窗口的 dH/dλ 分布直方图
- 显示采样充分性
- 标注平均值和标准差

**4. λ 窗口对比**:
- 多个子图展示不同 λ 窗口
- 显示能量收敛情况
- 标注异常窗口

### 10.3 报告模板

**HTML 结构**:
```html
<!DOCTYPE html>
<html>
<head>
    <title>PRISM-FEP Report: ligand_A → ligand_B</title>
    <script src="https://cdn.jsdelivr.net/npm/plotly.js@latest"></script>
    <style>
        /* 报告样式 */
        body { font-family: Arial, sans-serif; }
        .plot-container { width: 100%; height: 400px; }
        .info-box { background: #f0f0f0; padding: 15px; }
    </style>
</head>
<body>
    <h1>FEP Calculation Report</h1>
    
    <div class="info-box">
        <h2>System Information</h2>
        <p>Reference: {{ ref_name }}</p>
        <p>Mutant: {{ mut_name }}</p>
        <p>Force Field: {{ ff_name }}</p>
        <p>ΔG: {{ delta_g }} ± {{ error }} kJ/mol</p>
    </div>
    
    <h2>Free Energy Profile</h2>
    <div id="ddg-plot" class="plot-container"></div>
    
    <h2>Atom Mapping</h2>
    <div id="mapping-viz"></div>
    <img src="structures/mapping.png" alt="Atom Mapping">
    
    <h2>Energy Distribution</h2>
    <div id="lambda-plot" class="plot-container"></div>
    
    <script>
        // Plotly 图表生成代码
        // 使用 JSON 数据从 Python 传递
    </script>
</body>
</html>
```

### 10.4 数据收集

**收集的数据**:
1. 系统信息：配体名称、力场、参数
2. 计算结果：ΔG、误差、收敛性
3. 原子映射：common/transformed/surrounding 列表
4. 能量数据：各 λ 窗口的 dH/dλ 值
5. 验证结果：grompp 输出、警告信息

**数据格式**:
```python
{
    "system": {
        "reference": "ligand_A",
        "mutant": "ligand_B",
        "force_field": "gaff",
        "lambda_points": 11
    },
    "results": {
        "delta_g": -5.23,
        "error": 0.45,
        "converged": true
    },
    "mapping": {
        "common": 25,
        "transformed_a": 5,
        "transformed_b": 7,
        "surrounding": 3
    },
    "energy_data": {
        "lambdas": [0.0, 0.1, ..., 1.0],
        "dhdl": [10.2, 9.8, ..., -2.1],
        "errors": [0.5, 0.4, ..., 0.3]
    }
}
```

---

## 11. 参考资料

- FEbuilder: `/home/gxf1212/data2/work/make_hybrid_top/FEbuilder/`
- GROMACS FEP: https://manual.gromacs.org/current/fep.html
- pmx: https://github.com/deGrootLab/pmx
- PRISM CGenFF: `/data2/gxf1212/work/PRISM/prism/forcefield/cgenff.py`
- PRISM-Tutorial: `/home/gxf1212/data/work/PRISM-Tutorial`
- 参考体系与脚本: `tests/gxf/FEP/ref`
