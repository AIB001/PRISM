# PRISM-FEP Module

PRISM的自由能微扰（Free Energy Perturbation, FEP）计算模块，用于相对结合自由能计算。

## 模块结构

```
prism/fep/
├── __init__.py      # 模块初始化和导出
├── atom.py             # 原子和原子映射类
├── hybrid_topology.py       # 混合拓扑生成
├── itp_builder.py          # ITP文件构建器
├── analysis.py           # FEP结果分析
├── fep_integration.py      # PRISM集成接口
└── config.py               # FEP配置管理
```

## CLI使用

### 基本用法

```bash
prism protein.pdb reference.mol2 -o output_dir \
  --fep \
  --mutant mutant.mol2
```

### 完整参数

```bash
prism protein.pdb reference.mol2 -o output_dir \
  --fep \
  --mutant mutant.mol2 \
  --distance-cutoff 0.6 \
  --charge-strategy mean \
  --lambda-windows 11 \
  --fep-config fep_config.yaml
```

## CLI参数说明

- `--fep`: 启用FEP模式
- `--mutant`: 突变体配体文件路径（必需）
- `--fep-config`: FEP配置文件路径（可选）
- `--distance-cutoff`: 原子映射距离阈值（默认: 0.6 Å）
- `--charge-strategy`: 公共原子电荷策略
  - `ref`: 使用参考配体电荷
  - `mut`: 使用突变体配体电荷
  - `mean`: 使用平均电荷（默认）
- `--lambda-windows`: Lambda窗口数量（默认: 11）

## Python API

### 基本使用

```python
from prism.builder.core import PRISMBuilder

builder = PRISMBuilder(
    protein_path="protein.pdb",
    ligand_paths=["reference.mol2"],
    output_dir="output",
    fep_mode=True,
    mutant_ligand="mutant.mol2",
    distance_cutoff=0.6,
    charge_strategy='mean',
    lambda_windows=11
)

builder.build()
```

### 使用配置文件

```python
builder = PRISMBuilder(
    protein_path="protein.pdb",
    ligand_paths=["reference.mol2"],
    output_dir="output",
    fep_mode=True,
    mutant_ligand="mutant.mol2",
    fep_config="fep_config.yaml"
)
```

## 配置文件格式

```yaml
fep:
  distance_cutoff: 0.6
  charge_strategy: mean
  lambda_windows: 11
  
  # 可选的高级设置
  soft_core:
    alpha: 0.5
    sigma: 0.3
  
  mdp_settings:
    nsteps: 5000000
    dt: 0.002
```

## 核心类

### Atom
表示单个原子及其属性。

### AtomMapping
管理参考配体和突变体配体之间的原子映射。

### HybridAtom
表示混合拓扑中的原子，包含A态和B态参数。

### ITPBuilder
构建GROMACS混合拓扑ITP文件。

### XVGParser
解析GROMACS XVG输出文件。

### FEstimator
使用MBAR/BAR方法估算自由能差。

## 工作流程

1. **原子映射**: 自动识别参考配体和突变体配体之间的公共原子
2. **混合拓扑生成**: 创建包含A态和B态参数的混合拓扑
3. **Lambda窗口设置**: 配置多个lambda值进行渐进式转换
4. **MD模拟**: 在每个lambda窗口运行平衡和生产模拟
5. **自由能分析**: 使用MBAR/BAR方法计算ΔΔG

## 集成状态

✓ CLI参数集成完成
✓ PRISMBuilder接口集成完成
✓ 配置系统集成完成
✓ 模块导入验证通过
✓ **原子映射算法已实现** ⭐
✓ **文件I/O模块已实现** ⭐
✓ **单元测试已完成** ⭐

## 已实现功能（2025-03-11）

### 1. 原子映射算法 ⭐
- `DistanceAtomMapper.map()` - 基于距离的原子匹配
- 支持可调距离阈值（默认0.6Å）
- 自动分类：common/transformed/surrounding
- 已通过真实配体测试（25-36）

### 2. 文件I/O模块 ⭐
- `read_ligand_from_prism()` - 从PRISM生成的ITP+GRO读取
- `read_mol2_atoms()` - 直接从MOL2读取（测试用）
- 自动坐标转换（nm → Å）
- 复用PRISM现有解析逻辑

### 3. 单元测试 ⭐
- 测试文件：`tests/fep/test_mapping_25_36.py`
- 覆盖：文件读取、映射算法、阈值对比
- 所有测试通过 ✅

### 使用示例

```python
from prism.fep import read_mol2_atoms, DistanceAtomMapper

# 读取配体
lig25 = read_mol2_atoms('tests/gxf/FEP/test/hif2a/25-36/25_3D.mol2')
lig36 = read_mol2_atoms('tests/gxf/FEP/test/hif2a/25-36/36_3D.mol2')

# 执行映射
mapper = DistanceAtomMapper(dist_cutoff=0.6)
mapping = mapper.map(lig25, lig36)

# 查看结果
print(f"Common: {len(mapping.common)}")
print(f"Transformed: {len(mapping.transformed_a) + len(mapping.transformed_b)}")
```

## 下一步开发

- [ ] 实现双拓扑构建（DualTopologyBuilder）
- [ ] 实现ITP文件生成（ITPBuilder）
- [ ] 实现MDP模板生成
- [ ] 集成PRISM力场系统（从ITP+GRO读取）
- [ ] 实现FEP分析工具（XVG解析、BAR/MBAR）
- [ ] 添加更多测试案例
- [ ] 性能优化和错误处理

## 参考

- GROMACS FEP教程: http://www.gromacs.org/Documentation/Tutorials
- Alchemical Free Energy Calculations: https://www.alchemistry.org/
