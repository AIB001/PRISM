# PRISM-FEP 需求文档

## 1. 项目概述

### 1.1 目标
将 FEbuilder 的 FEP (Free Energy Perturbation) 双拓扑构建能力移植到 PRISM，支持 **GROMACS 的 A/B 状态拓扑**（`typeB/chargeB`）格式。

### 1.2 背景
| 角色 | 说明 |
|------|------|
| FEbuilder | 现有工具，支持 NAMD/CHARMM FEP 构建 |
| PRISM | GROMACS 的蛋白-配体建模工具 |
| 需求 | 统一工作流，使用 GROMACS 进行 FEP |

### 1.3 核心价值
| 价值点 | 说明 |
|--------|------|
| 复用建模能力 | 复用 PRISM 现有力场参数化流程 (GAFF/GAFF2/OpenFF/CGenFF/OPLS-AA) |
| 工具链统一 | MD + FEP + PMF + REST2 统一分析入口 |
| 用户覆盖 | 面向 GROMACS 用户群体 |

### 1.4 主要功能
| 功能 | 说明 |
|------|------|
| 双拓扑构建 | 输入两个配体，输出 GROMACS A/B 双状态拓扑 (hybrid.itp) |
| 原子映射 | 基于距离匹配，自动识别 common/transformed/surrounding 原子 |
| 力场支持 | GAFF/GAFF2/OpenFF/CGenFF/OPLS-AA |
| 系统构建 | 复合物组装、溶剂化、离子添加 |
| FEP-MDP | 自动生成 GROMACS FEP 计算所需的 MDP 文件 |
| 分析工具 | 解析 dhdl.xvg，计算 ΔΔG (支持 BAR/MBAR) |

### 1.5 参考基线与设计文档
| 项 | 内容 |
|----|------|
| **参考基线位置** | `tests/gxf/FEP/ref` |
| 拓扑格式 | GROMACS FEP 标准格式（示例 `system/newtop.top`） |
| 关键格式 | `[ atoms ]` 含 `typeB/chargeB/massB`，dummy 类型来自力场（`DUM_*`） |
| 参考力场 | `/home/gxf1212/data/local_programs/mutff/amber14sbmut.ff` |
| 运行模式 | λ-REMD (`gmx mdrun -multidir -replex`) 与非交换多窗口 |
| 分析基线 | `calculate.sh` 收集 `dhdl.xvg` + `gmx bar` 输出 `summary.dat` |

**详细设计文档**（本文档姊妹篇）：
| 文档 | 内容要点 |
|------|----------|
| `CORE_DESIGN.md` | 原子映射算法、双拓扑构建逻辑、HybridAtom 数据结构 |
| `OUTPUT_FORMAT.md` | ITP 格式、MDP 模板、dummy 原子处理 |
| `ANALYSIS.md` | XVG 解析、BAR/MBAR 方法、自由能计算 |

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
│   ├── core.py                    # PRISMBuilder 核心类（扩展支持 FEP）
│   ├── cli.py                     # 统一 CLI（添加 --fep 等参数）
│   └── fep_integration.py         # FEP 模式集成逻辑
└── tests/
    ├── fep/                       # 单元测试
    │   ├── test_mapping.py
    │   ├── test_dual_topology.py
    │   └── test_itp_builder.py
    └── gxf/FEP/
        ├── test/                  # 真实测试用例（case.yaml + 资源）
        │   ├── run_fep_tests.sh
        │   └── quick_test.sh
        └── ref/                   # 已跑通基线体系（pmx + amber14sbmut）
```

**架构说明**：
- **核心模块** (`fep/`)：FEP 算法实现，与 builder 解耦
- **集成模块** (`builder/fep_integration.py`)：连接 PRISMBuilder 和 FEP 模块
- **统一接口** (`builder/cli.py`)：在现有 CLI 基础上添加 `--fep` 等参数
- **向后兼容**：不影响现有 `prism` 命令的普通建模功能

---

## 3. 测试系统要求 (2024-03-12更新)

### 3.1 测试目录结构

```
tests/gxf/FEP/
├── unit_test/                      # 单元测试（快速验证）
│   ├── test_framework.py          # 框架测试（23个测试，全部通过）
│   ├── test_cgenff_mapping.py     # CGenFF端到端测试
│   ├── validate_vs_febuilder.py   # FEbuilder结果验证
│   ├── verify_cli_integration.py  # CLI集成验证
│   ├── test_charge_cutoff_experiments.py  # 参数实验
│   ├── verify_config_consistency.py       # 配置一致性验证
│   ├── 25-36/                     # 测试案例1
│   │   ├── case.yaml              # 测试配置（含test_requirements）
│   │   ├── config.conf            # FEbuilder兼容配置
│   │   ├── 25.rtf/25.pdb          # CGenFF文件
│   │   ├── 36.rtf/36.pdb
│   │   └── other/hybrid.pdb       # FEbuilder参考结果
│   └── 39-8/                      # 测试案例2
│       ├── case.yaml
│       └── ...
├── test/                          # 集成测试（真实案例）
│   ├── hif2a/                     # HIF2α体系
│   ├── RdRp/                      # RdRp体系
│   ├── P38/                       # P38体系
│   └── FEP+T4/                    # T4体系
└── ref/                           # 参考基线（pmx结果）
```

### 3.2 配置文件规范

#### case.yaml 格式
```yaml
case:
  group: hif2a
  name: 25-36
  path: hif2a/25-36

fundamental:
  ref: 25
  mut: 36
  rec: receptor.pdb
  ff_path: ../toppar

model:
  no_solvate: false
  charge_common: ref              # 与config.conf一致
  charge_reception: pert
  distance: 12

# 测试要求（unit_test专用）
test_requirements:
  mapping:
    dist_cutoff: 0.6
    charge_cutoff: 0.05           # 默认值
    charge_cutoff_alt: 0.2        # 备选值
    charge_common: ref

  expected:
    min_common_atoms: 25
    febuilder_match: true/false
    febuilder_uncommon_ref:        # FEbuilder的uncommon atoms
      ref_state: [...]
      mut_state: [...]

simulation:
  numofsteps: 500000
  repeats: 4

other:
  verbose: true
```

#### config.conf 格式（FEbuilder兼容）
```ini
[Model]
charge_common = ref
charge_reception = pert
distance = 12

[Other]
dist_cutoff = 0.6
charge_cutoff = 0.05
```

**关键要求**：
- ✅ `case.yaml` 和 `config.conf` 的基础参数必须一致
- ✅ `test_requirements` 部分仅存在于 `unit_test/` 目录
- ✅ `test/` 目录的配置保持与FEbuilder完全兼容

### 3.3 测试验证要求

#### 基础验证（所有案例）
- ✅ 能成功读取RTF+PDB文件
- ✅ 映射算法不报错
- ✅ 产生合理的分类结果

#### 高级验证（unit_test）
- ✅ 与FEbuilder结果逐原子对比
- ✅ 验证uncommon atoms分类
- ✅ 分析电荷差异原因
- ✅ 测试不同charge_cutoff值的效果

#### 配置一致性验证
- ✅ 使用 `verify_config_consistency.py` 检查所有测试案例
- ✅ 确保 `case.yaml` 和 `config.conf` 参数一致
- ✅ 验证基础参数（charge_common, charge_reception, distance）

### 3.4 测试运行控制

#### Pytest Markers
```python
@pytest.mark.slow  # 标记慢速测试（如25-36, 39-8）
@pytest.mark.integration  # 标记集成测试
```

#### 运行命令
```bash
# 快速测试（跳过慢速）
pytest tests/gxf/FEP/unit_test/ -m "not slow" -v

# 完整测试（包含慢速）
pytest tests/gxf/FEP/unit_test/ -m "slow" -v

# 特定系统测试
pytest tests/gxf/FEP/unit_test/test_cgenff_mapping.py::test_cgenff_mapping_with_different_systems[39-8] -v

# 参数实验
pytest tests/gxf/FEP/unit_test/test_charge_cutoff_experiments.py -v
```

### 3.5 当前测试状态 (2025-03-12)

#### ✅ 已验证案例 (完全匹配 FEbuilder)
| 案例 | Common | Surrounding (总和) | Transformed | charge_cutoff | 配置来源 |
|------|--------|-------------------|-------------|---------------|----------|
| 25-36 | 33 | 6 (lig25: C9,C10,H7 / lig36: C9,F1,H5) | 1 (N1) | 0.2 | case.yaml |
| 39-8 | 31 | 14 (lig39: 7个 / lig8: 7个) | 0 | 0.05 | case.yaml |

**验证方法**：
```bash
pytest tests/gxf/FEP/unit_test/validate_vs_febuilder.py -v -s -m slow
# 结果：两个系统完全匹配 FEbuilder 的 hybrid.pdb 标记
```

**配置架构**：
- `case.yaml` → `fep.mapping` → `ConfigurationManager` → `DistanceAtomMapper.from_config()`
- 测试文件无需任何默认值逻辑，全部由 PRISM 内部处理

#### 待添加案例
- [ ] FEP+T4 体系（4个案例）
- [ ] RdRp 体系
- [ ] P38 体系
- [ ] 更多HIF2α案例

### 3.6 关键发现

**charge_cutoff 参数建议**：
- 默认值：保持 **0.05**（不修改）
- 大多数情况：0.05 足够（如39-8）
- 特殊情况：可调整到 0.2（如25-36）
- 配置方式：通过 `case.yaml` 的 `test_requirements.mapping.charge_cutoff_alt` 指定

**重要结论**：
- ✅ 不需要实现复杂的 `charge_common='ref'` 策略
- ✅ 通过调整 `charge_cutoff` 即可匹配FEbuilder结果
- ✅ 保持默认参数不变，通过配置文件灵活调整

---

## 4. 三种原子类型分类标准

### 4.1 核心概念

FEbuilder将原子分为三种类型，每种在GROMACS FEP中有不同的处理方式：

| 类型 | 定义 | GROMACS ITP格式 | FEbuilder标记 | Dummy需求 |
|------|------|-----------------|--------------|-----------|
| **Common (公共)** | 类型相同，电荷相同 | typeA=typeB, chargeA=chargeB, massA=massB | beta=0 | 不需要 |
| **Surrounding (环境)** | 类型相同，电荷不同 | typeA=typeB, chargeA≠chargeB, massA=massB | beta≠0 (无A/B) | **不需要** |
| **Transformed (非公共)** | 类型不同，需要dummy | typeA≠typeB, 需dummy类型, chargeB=0 | beta≠0 (有A/B) | **需要** |

### 4.2 关键区别

**Surrounding原子**:
- ✅ 虽然在"uncommon"集合中，但**不需要dummy类型**
- ✅ 只需在ITP中使用typeB/chargeB列
- ✅ 原子类型保持不变，只是电荷值不同
- ✅ 例如：C9在25和36中都是CG2R66，但电荷不同

**Transformed原子**:
- ✅ 需要**完整的dummy类型替换**
- ✅ A状态使用真实类型，B状态使用dummy类型
- ✅ 或反之，取决于mutant/reference
- ✅ 例如：N1在25中存在，在36中不存在

**Common原子**:
- ✅ A/B状态完全相同，不需要任何特殊处理
- ✅ 在ITP中只需填写一次参数

### 4.3 FEbuilder hybrid.pdb标记规则

```
ATOM      1  N1   HYB L   1      18.205   3.518  -7.118  1.00  0.00      LIG N
                                                              ^^^^
                                                          beta=0: Common

ATOM     20 C15A  HYB L   1      26.554   5.241  -7.785  1.00 -1.00      LIG C
                                                              ^^^^^
                                                         beta=-1: State A (Transformed)

ATOM     39 C15B  HYB L   1      26.604   5.233  -7.747  1.00  1.00      LIG C
                                                              ^^^^^
                                                         beta=+1: State B (Transformed)
```

**标记说明**：
- **beta=0**: Common原子
- **beta=-1**: Transformed原子（State A存在，State B为dummy）
- **beta=+1**: Transformed原子（State B存在，State A为dummy）
- **无A/B后缀**: Surrounding原子（虽然beta≠0，但无A/B标记）

### 4.4 实现要求

**算法要求**：
- ✅ 我们的`AtomMapping`必须区分三种类型
- ✅ Surrounding原子不能归类为Transformed
- ✅ ITP生成时，Surrounding只需typeB/chargeB
- ✅ Transformed需要完整的dummy类型替换

**数据结构**：
```python
@dataclass
class AtomMapping:
    common: List[Tuple[Atom, Atom]]           # Common atoms
    transformed_a: List[Atom]                  # Transformed in A (need dummy in B)
    transformed_b: List[Atom]                  # Transformed in B (need dummy in A)
    surrounding_a: List[Atom]                  # Surrounding in A (need typeB/chargeB)
    surrounding_b: List[Atom]                  # Surrounding in B (need typeB/chargeB)
```

### 4.5 25-36系统实例分析

**Surrounding原子 (Ligand 36)**:
| 原子 | 36.rtf类型 | 电荷 | 说明 |
|------|-----------|------|------|
| C9 | CG2R66 | +0.286 | 类型相同，电荷不同，不需dummy |
| F1 | FGR1 | -0.196 | 类型不同但在距离内，surrounding处理 |
| H5 | HGR62 | +0.177 | 类型相同，电荷不同，不需dummy |

**Transformed原子 (Ligand 25)**:
| 原子 | 25.rtf类型 | 电荷 | 说明 |
|------|-----------|------|------|
| N1 | NG321 | -0.781 | 在36中不存在，需要dummy |

**验证结果**：
- ✅ **Surrounding**: 6个原子 (3+3)
- ✅ **Transformed**: 1个原子 (N1)
- ✅ **Common**: 33个原子
- ✅ **完全匹配FEbuilder分类**

---

## 5. CLI 接口

### 3.1 设计原则

PRISM-FEP 采用**统一接口设计**，在现有 `prism` 命令基础上扩展 FEP 功能，保持接口一致性和工作流集成性。

**核心优势**：
- ✅ 保持现有 `prism` 命令不变，用户无需学习新命令
- ✅ 通过 `--fep` 开关启用 FEP 模式，不影响普通建模用户
- ✅ 支持在同一工作流中集成 MD、FEP、PMF 等多种计算
- ✅ 配置文件可复用，MD 和 FEP 参数可统一管理

### 3.2 基本用法

```bash
# 普通 MD 建模（保持现有接口）
prism protein.pdb ligand.mol2 -o output_dir

# FEP 模式（添加 --fep 开关）
prism protein.pdb ligand.mol2 -o output_dir \
  --fep \
  --mutant mutant.mol2

# 使用 CGenFF 力场
prism protein.pdb ligand.mol2 -o output_dir \
  --fep \
  --mutant mutant.mol2 \
  --ligand-forcefield cgenff \
  --forcefield-path toppar/ \
  --forcefield charmm36

# 完整参数示例
prism protein.pdb ligand.mol2 -o fep_system \
  --fep \
  --mutant mutant.mol2 \
  --ligand-forcefield gaff2 \
  --forcefield amber14sb \
  --fep-config fep_setup.yaml \
  --config md_config.yaml
```

### 3.3 配置文件设计

**方案 1：分层 YAML（推荐）**
```yaml
# config.yaml
system:
  protein: protein.pdb
  ligand: ligand.mol2
  forcefield: amber14sb
  ligand_forcefield: gaff2

fep:
  enabled: true
  mutant: mutant.mol2
  distance_cutoff: 0.6
  charge_strategy: mean
  lambda_windows: 11
  soft_core:
    alpha: 0.5
    power: 1

md:
  temperature: 300
  pressure: 1.0
  time_step: 0.002
```

**方案 2：独立 FEP 配置**
```yaml
# fep_setup.yaml
mutant: mutant.mol2
distance_cutoff: 0.6
charge_strategy: mean  # ref | mut | mean
lambda_windows: 11
coulomb_lambdas: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
vdw_lambdas: [0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.0]
soft_core:
  alpha: 0.5
  power: 1
  sigma: 0.3
```

### 3.4 命令行参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--fep` | 启用 FEP 模式 | False |
| `--mutant` | 突变配体文件（FEP 必需） | - |
| `--fep-config` | FEP 配置文件路径 | fep_setup.yaml |
| `--distance-cutoff` | 原子映射距离阈值（Å） | 0.6 |
| `--charge-strategy` | 共用原子电荷策略 | mean |
| `--lambda-windows` | λ 窗口数量 | 11 |
| `--soft-core-alpha` | 软核心 α 参数 | 0.5 |

### 3.5 分析结果

```bash
# 分析 FEP 结果
prism analyze --xvg dhdl.xvg --method bar

# 生成 HTML 报告
prism report --output-dir fep_system --report-type html
```

---

## 4. GROMACS vs NAMD 的关键差异

| 维度 | NAMD (FEbuilder) | GROMACS | PRISM-FEP 结论 |
|------|------------------|---------|----------------|
| **参数处理** | 需要参数合并 (`merge_prm.py`)，两套参数去重合并，逻辑复杂 | **不需要参数合并**，依赖 typeB/chargeB 列，只需保持索引正确 | ✗ 不需要实现参数值合并<br>✓ 需要原子索引映射<br>✓ A/B 参数选择由 GROMACS 自动处理 |
| **拓扑文件格式** | `.rtf` + `.prm` (需要合并) | `.itp` (typeB/chargeB，pmx 还会含 massB) | 拓扑构建更简洁，只需原子索引映射 |
| **实现影响** | 需要复杂的冲突处理逻辑和参数值合并 | **更简洁**：GROMACS 根据原子的 type/charge 自动选择状态 | 无需 `merge_prm.py` 逻辑，简化开发工作量 |

---

## 5. 力场支持

### 5.1 PRISM 已有力场支持（FEP 直接复用）

**核心优势**：FEP **不需要新增力场支持**，直接复用 PRISM 现有力场生成器。所有力场均输出**统一的 GROMACS 格式**（.itp + .gro），FEP 只需读取并构建双拓扑。

| 力场 | CLI 参数 | PRISM 模块 | 输入格式 | 输出目录 | 依赖 |
|------|----------|-----------|----------|----------|------|
| **GAFF** | `--force-field gaff` | `prism.forcefield.gaff` | .mol2 | `LIG.amb2gmx/` | AmberTools |
| **GAFF2** | `--force-field gaff2` | `prism.forcefield.gaff2` | .mol2 | `LIG.amb2gmx/` | AmberTools |
| **OpenFF** | `--force-field openff` | `prism.forcefield.openff` | .sdf/.mol2 | `LIG.openff2gmx/` | openff-toolkit |
| **CGenFF** | `--force-field cgenff` | `prism.forcefield.cgenff` | .rtf + .prm | `LIG.cgenff2gmx/` | 网络下载 |
| **OPLS-AA** | `--force-field oplsaa` | `prism.forcefield.opls_aa` | .mol2 | `LIG.opls2gmx/` | LigParGen |
| **MMFF** | `--force-field mmff` | `prism.forcefield.swissparam` | .mol2 | `LIG.mmff2gmx/` | curl (SwissParam API) |
| **MATCH** | `--force-field match` | `prism.forcefield.swissparam` | .mol2 | `LIG.match2gmx/` | curl (SwissParam API) |
| **Both (MMFF+MATCH)** | `--force-field both` | `prism.forcefield.swissparam` | .mol2 | `LIG.both2gmx/` | curl (SwissParam API) |

### 5.2 CGenFF 特别说明

| 项 | 说明 |
|----|------|
| 模块 | `prism/forcefield/cgenff.py` |
| 格式 | CHARMM-GUI 输出 (`.rtf` + `.prm`) |
| 适用 | 含药物分子的系统 |
| 测试数据位置 | `tests/gxf/FEP/test/hif2a/42-38/` |
| 示例文件 | `toppar/lig.rtf`, `toppar/lig.prm` |

### 5.3 FEP 的力场处理流程

```
┌─────────────────────────────────────────────────────────┐
│  Step 1: 力场参数化（PRISM 现有模块，无需修改）            │
│  prism protein.pdb ligA.mol2 --fep --mutant ligB.mol2    │
│           --ligand-forcefield gaff                       │
└──────────────────────┬──────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────┐
│  Step 2: 读取 GROMACS 格式（FEP 新增）                    │
│  - 读取 ligA.amb2gmx/ligand.itp                          │
│  - 读取 ligB.amb2gmx/ligand.itp                          │
│  - 提取原子坐标、类型、电荷、参数                          │
└──────────────────────┬──────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────┐
│  Step 3: 原子映射 + 双拓扑构建（FEP 新增）                │
│  - mapping.py: 距离匹配，分类原子                         │
│  - dual_topology.py: 生成 hybrid.itp (含 typeB/chargeB)  │
└─────────────────────────────────────────────────────────┘
```

**关键点**：
- ✅ FEP **不需要处理力场参数化**，直接读取 PRISM 输出
- ✅ 所有力场输出**统一的 GROMACS 格式**，FEP 只需一套代码
- ✅ 新增力场支持（如未来新增）对 FEP **透明**，自动兼容
---

## 6. 测试计划与自动化

### 6.1 测试数据
| 项 | 说明 |
|----|------|
| 位置 | `tests/gxf/FEP/test/` |
| 格式 | 每个 case 必须包含 `case.yaml` + 配体文件 + 受体 |
| 说明 | `config.conf` 为历史输入，已生成 `case.yaml`，脚本以 `case.yaml` 为准 |
| 支持 | GAFF/OpenFF/CGenFF |
| 示例 | `hif2a/42-38/` (CHARMM-GUI 输出) |
| 参考基线 | `tests/gxf/FEP/ref`（pmx 蛋白突变拓扑 + 可运行脚本） |

**测试集合覆盖（现有样例）**：
| 组别 | Case | 关键输入/设置 | 备注 |
|------|------|---------------|------|
| FEP+T4 | Indole-FuPh 等 5 组 | `params=match`, `distance=8`, `temperature=302.15`, `concentration=0.05`, `disoff=1` | 配体含 `.rtf/.prm`，偏 CGenFF/CHARMM |
| RdRp | 1p-cn-cch | `recpsf` + `ff_path=toppar`, `no_solvate=True`, `charge_mode=mode2` | CHARMM/PSF 拓扑路径 |
| hif2a | 42-38, 39-8, 25-36 | `charge_mode=mode2`, `distance=12`, `no_solvate` 变化 | CGenFF/CHARMM，小分子配体 |
| P38 | 13-12 等 5 组 | `add_pdb=crw.pdb`, `distance=8`, `repeats=6` | 额外蛋白片段输入 |

**非自动化参考数据（仅供对照）**：
| 目录 | 说明 |
|------|------|
| `tests/gxf/FEP/test/CDI_FEP` | SwissParam/ParamChem 相关参考与旧流程脚本 |
| `tests/gxf/FEP/test/others` | 额外 CGenFF/CHARMM 参考拓扑与构建记录 |

### 6.2 测试脚本设计

**批量测试脚本**：`tests/gxf/FEP/test/run_fep_tests.sh`

**设计思路**（伪代码）：
```bash
#!/bin/bash
# PRISM-FEP 批量测试脚本

# 1. 遍历所有 case.yaml
for case_yaml in $(find tests/gxf/FEP/test -name case.yaml); do
    case_dir=$(dirname "$case_yaml")

    # 2. 使用 case.yaml 运行（主入口）
    prism --config "$case_yaml"

    # 3. 验证生成的拓扑（可选）
    # cd "$case_dir/output/GMX_PROLIG_FEP"
    # gmx grompp -f em.mdp -c conf.gro -p topol.top -o test.tpr
done

# 4. 生成汇总报告（后续实现）
generate_summary_report
```

**快速验证脚本**：`tests/gxf/FEP/test/quick_test.sh`
```bash
#!/bin/bash
# 快速验证单个测试案例
TEST_DIR="${1:?}"
prism --config "$TEST_DIR/case.yaml"
```

### 6.3 测试报告格式

**控制台输出**：
```
================================================================================
                    PRISM-FEP 批量测试报告
================================================================================
测试日期：2025-01-15 14：30：25
测试目录：tests/gxf/FEP/test

[✓] RdRp/1p-cn-cch     remtp → remtp-1p-cn-cch
[✗] RdRp/2p-oh-och2ch3 remtp → remtp-2p-oh-och2ch3    (GROMACS验证失败)
[?] P38/42-38          lig42 → lig38                   (参数化超时)
[✓] hif2a/5-1         ref5 → mut1

测试完成：4/4
成功：2  失败：1  超时：1
详细报告：tests/gxf/FEP/output/test_report_20250115_143025.txt
================================================================================
```

**文本报告**：`tests/gxf/FEP/output/test_report.txt`
```text
================================================================================
测试案例：RdRp/1p-cn-cch
参考配体：remtp
突变配体：remtp-1p-cn-cch

结果：✓ 通过
步骤：
  [1] 参数化：✓ 成功 (耗时：12.3s)
      - GAFF 力场参数生成完成
      - 生成的杂化拓扑包含 62 个原子

  [2] 系统构建：✓ 成功 (耗时：8.5s)
      - 蛋白-配体复合物构建完成
      - 加溶剂：TIP3P, 12.0 Å 盒子
      - 加离子：0.15 M NaCl, 中和系统

  [3] 拓扑验证：✓ 成功 (耗时：2.1s)
      - GROMACS grompp 验证通过
      - 无警告或错误

输出目录：output/RdRp/1p-cn-cch/GMX_PROLIG_FEP
================================================================================
```

### 6.4 预期问题与错误诊断

| 错误信息 | 原因 | 解决方案 |
|---------|------|----------|
| Atom type not found | 力场参数缺失 | 检查 GAFF/OpenFF 参数化结果 |
| masses do not add up | dummy 类型或 massB 缺失 | 确保 dummy 类型来自 forcefield，且 `[ atoms ]` 含 `massB` |
| Inconsistent atom indices | 参数映射错误 | 检查 hybrid_atoms 索引映射 |
| grompp failed | 拓扑文件格式错误 | 验证 ITP 文件格式 |
| CHARMM 格式不兼容 | 缺少 .rtf/.prm | 使用 PRISM 重新参数化或直接读取 |
| 虚拟原子问题 | dummy 类型未定义 | 使用 dummy 类型（由力场提供，amber14sbmut 中为 `DUM_*`） |
| 力场兼容性 | 不同力场混用 | 统一使用一种力场，记录版本 |

### 6.5 覆盖缺口与新增建议
| 覆盖点 | 建议 |
|--------|------|
| GAFF/GAFF2 | 新增 1-2 个仅 mol2/sdf 的小分子对（无 rtf/prm） |
| OpenFF | 新增 1 个 smirnoff/SMIRKS 的 ligand pair |
| OPLS-AA | 新增 1 个 OPLS 小分子对 |
| pmx 蛋白突变 | 新增 1 个直接读取 `newtop.top` 的测试入口 |
| 纯溶剂化 | 新增 1 例 `no_solvate=False` 的小分子对 |

---

## 7. 原子映射可视化 (已完成)

### 7.1 模块架构

**已完成的可视化模块**：
```
prism/fep/visualize/
├── __init__.py              # 导出 visualize_mapping_png, visualize_mapping_html
├── molecule.py              # 分子处理工具
│   ├── pdb_to_mol()                    # PDB → RDKit Mol 对象
│   ├── assign_bond_orders_from_mol2()  # 从 mol2 校正键级（RDKit）
│   └── prepare_mol_with_charges_and_labels()  # 添加电荷和标签
├── highlight.py             # 高亮颜色定义（FEbuilder 风格）
│   ├── COMMON_COLOR      = rgb(204, 229, 77)  # 绿色
│   ├── TRANSFORMED_COLOR = rgb(255, 77, 77)   # 红色
│   └── SURROUNDING_COLOR = rgb(77, 153, 255)  # 蓝色
├── mapping.py               # PNG 可视化（FEbuilder 风格）
│   └── visualize_mapping_png()  # MCS 对齐、键级校正、电荷标注
└── html.py                  # HTML 交互式可视化
    └── visualize_mapping_html()  # 独立 HTML 页面
```

### 7.2 PNG 可视化功能

**核心特性**：
| 特性 | 实现方式 | 状态 |
|------|----------|------|
| **MCS 对齐** | RDKit `rdFMCS` 最大公共子结构对齐 | ✅ |
| **键级校正** | 从 mol2 读取正确键级，使用 `AssignBondOrdersFromTemplate` | ✅ |
| **电荷标注** | 显示所有原子电荷（包括氢原子） | ✅ |
| **FEbuilder 风格** | 使用相同的 `MolDrawOptions` 和颜色方案 | ✅ |
| **高亮分类** | Common=绿, Transformed=红, Surrounding=蓝 | ✅ |
| **图例说明** | 显示各类型原子数量和处理方式 | ✅ |

**使用示例**：
```python
from prism.fep.visualize import visualize_mapping_png

visualize_mapping_png(
    mapping,                    # AtomMapping 对象
    pdb_a='25.pdb',             # 参考 PDB
    pdb_b='36.pdb',             # 突变 PDB
    mol2_a='25.mol2',           # 参考 mol2（键级校正）
    mol2_b='36.mol2',           # 突变 mol2
    output_path='25-36_mapping.png',
    edge=0.15,                  # 图像边距
    legends=('Ligand 25', 'Ligand 36')
)
```

### 7.3 HTML 交互式可视化（设计阶段）

**交互功能需求**：

| 功能 | 描述 | 优先级 | 实现建议 |
|------|------|--------|----------|
| **Hover 显示电荷** | 鼠标悬停显示原子电荷、类型、分类 | 高 | JavaScript + CSS tooltip |
| **对应原子高亮** | Hover 时高亮对面分子的对应原子 | 高 | 跨分子索引映射 |
| **电荷标签开关** | 切换是否在图像上直接显示电荷标签 | 中 | Checkbox 控制显示 |
| **导出图片** | 导出当前视图为 PNG（300 DPI） | 中 | Canvas API 转换 |
| **交互模式开关** | 配置项控制生成交互式 HTML 或静态 PNG | 高 | CLI 参数 `--interactive` |

**设计建议**：

1. **Hover Tooltip 内容**：
   ```
   原子: C1
   电荷: +0.286
   类型: CG2R66
   分类: Surrounding (lig25)
   对应: C9 (lig36, charge: +0.123)
   ```

2. **电荷标签显示策略**：
   - **默认关闭**：保持图像整洁
   - **开启时**：
     - 小分子（<30 原子）：显示所有电荷
     - 大分子（≥30 原子）：只显示差异原子（Surrounding、Transformed）

3. **HTML 独有功能**：
   | 功能 | 优势 |
   |------|------|
   动态筛选 | 点击图例显示/隐藏特定原子类型 |
   原子列表 | 侧边栏显示所有原子详细信息表 |
   搜索定位 | 搜索原子名称快速定位并高亮 |
   对比模式 | 并排显示参考和突变的详细信息 |
   数据导出 | 导出映射数据为 CSV/JSON 格式 |
   响应式设计 | 自适应屏幕尺寸，移动端友好 |

### 7.4 CLI 接口

```bash
# 生成 PNG 可视化（默认，FEbuilder 风格）
prism visualize-mapping --case 25-36 --output mapping.png

# 生成交互式 HTML
prism visualize-mapping --case 25-36 --output mapping.html --interactive

# 在 HTML 中显示电荷标签
prism visualize-mapping --case 25-36 --output mapping.html --show-charges

# 自定义图例
prism visualize-mapping --case 25-36 --legends "Reference" "Mutant"

# 调整图像边距
prism visualize-mapping --case 25-36 --edge 0.2
```

### 7.5 输出文件

| 输出类型 | 文件格式 | 用途 | 位置 |
|---------|---------|------|------|
| 静态图像 | PNG | 报告、论文、快速查看 | `output/mapping.png` |
| 交互式 | HTML | 深入分析、调试、演示 | `output/mapping.html` |
| 映射数据 | JSON | 数据交换、二次分析 | `output/mapping.json` |

### 7.6 验证状态

**已验证系统**：
| 系统 | Common | Surrounding | Transformed | 验证状态 |
|------|--------|-------------|-------------|----------|
| 25-36 | 33 | 6 (3+3) | 1 (N1) | ✅ 完全匹配 FEbuilder |
| 39-8 | 31 | 14 (7+7) | 0 | ✅ 完全匹配 FEbuilder |

**可视化验证**：
- ✅ PNG 输出与 FEbuilder 风格一致
- ✅ 键级校正正确（从 mol2 读取）
- ✅ 电荷显示完整（包括氢原子）
- ✅ 颜色高亮分类正确
- ⏳ HTML 交互功能（待实现）

---

## 8. FEP 结果分析报告（待实现）

### 8.1 功能需求

**模块**：`prism/fep/analysis/report.py`

生成完整的 FEP 计算报告，包含：
| 内容 | 说明 |
|------|------|
| 计算信息 | 案例与元数据 |
| 原子映射 | 嵌入 PNG/HTML 可视化 |
| ΔG 曲线 | ΔG vs λ |
| 能量分布 | λ 窗口直方图 |
| 误差评估 | 置信区间 |
| 收敛性分析 | HAMADA 统计量 |

### 8.2 设计思路

| 步骤 | 说明 |
|------|------|
| Plotly | 交互式图表（ΔG 曲线、能量分布） |
| RDKit | 2D 分子结构图（原子映射） |
| HTML | 统一报告模板 |
| 多窗口 | λ 窗口对比分析 |

### 8.3 CLI 接口

```bash
# 生成完整 HTML 报告
prism report --output-dir fep_system --report-type html

# 分析 FEP 结果
prism analyze --xvg dhdl.xvg --method bar

# 生成对比报告
prism report --compare system1 system2 --output comparison.html
```

### 8.4 可视化元素

| 可视化元素 | 说明 |
|------------|------|
| 自由能变化曲线 | λ (0→1) vs ΔG，标注误差 |
| 原子映射图 | RDKit 2D，高亮 transformed/surrounding/common |
| 能量分布 | 各 λ 窗口 dH/dλ 直方图 |
| λ 窗口对比 | 多子图展示收敛与异常 |

### 8.5 输出文件

| 输出 | 说明 |
|------|------|
| `fep_report.html` | 主报告文件 |
| `plots/` | 图表数据目录 |
| `structures/` | 分子结构图 |

### 8.6 数据收集

**收集的数据**：
1. 系统信息：配体名称、力场、参数
2. 计算结果：ΔG、误差、收敛性
3. 原子映射：common/transformed/surrounding 列表
4. 能量数据：各 λ 窗口的 dH/dλ 值
5. 验证结果：grompp 输出、警告信息

**数据格式**：
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

## 8. 参考资料

| 资源 | 位置 |
|------|------|
| FEbuilder | `/home/gxf1212/data2/work/make_hybrid_top/FEbuilder/` |
| GROMACS FEP | https://manual.gromacs.org/current/fep.html |
| pmx | https://github.com/deGrootLab/pmx |
| PRISM CGenFF | `/data2/gxf1212/work/PRISM/prism/forcefield/cgenff.py` |
| PRISM-Tutorial | `/home/gxf1212/data/work/PRISM-Tutorial` |
| 参考体系与脚本 | `tests/gxf/FEP/ref` |
| 测试用例集合 | `tests/gxf/FEP/test` |
