---
name: fep-analysis
description: FEP result analysis module debugging (multi-estimator TI/BAR/MBAR, automatic caching, convergence diagnostics)
type: fep
---

# FEP 结果分析模块调试

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

## 核心职责

调试 PRISM-FEP 分析模块，支持多估算器（TI/BAR/MBAR）对比分析、自动缓存、快速测试/完整分析模式。

## 快速检查清单

### ✅ 0. 数据准备

```bash
# 检查 XVG 文件是否存在
find <fep_dir> -name "*.xvg" | head -10

# 检查目录结构
# 应该看到: bound*/window_00/, unbound*/window_01/, ...
# 每个窗口下有: dhdl.xvg 或 prod.xvg
```

**参考数据位置**: `examples/fep_output/1st_RUN4/WT-M10-re/`（3次重复）

### ✅ 1. 快速测试（1-2分钟）

```bash
cd examples/fep_output/1st_RUN4/WT-M10-re
python test_multi_estimator.py
```

输出: `fep_multi_estimator_cache.json` + `fep_multi_estimator_report.html`

### ✅ 2. 完整分析（40-75分钟）

```bash
python test_multi_estimator.py --full
```

包含: time convergence + bootstrap (n=50) + 完整诊断

### ✅ 3. Python API

```python
from prism.fep.analysis import FEPMultiEstimatorAnalyzer

analyzer = FEPMultiEstimatorAnalyzer(
    bound_dirs=['bound1', 'bound2', 'bound3'],
    unbound_dirs=['unbound1', 'unbound2', 'unbound3'],
    estimators=['TI', 'BAR', 'MBAR'],
    skip_bootstrap=True,          # 快速模式
    skip_time_convergence=True,   # 快速模式
    cache_file='fep_multi_estimator_cache.json',
)

results = analyzer.analyze()

# 查看结果
for name, res in results.methods.items():
    print(f"{name}: ΔG = {res.delta_g:.3f} ± {res.delta_g_error:.3f}")

print(f"Agreement: {results.comparison['agreement']}")
print(f"ΔG range: {results.comparison['delta_g_range']:.3f}")
```

## 模块结构

```
prism/fep/analysis/
├── __init__.py
├── cli.py
├── plots.py
├── core/
│   ├── models.py          # FEResults, MultiEstimatorResults
│   ├── xvg_parser.py      # GROMACS XVG 解析
│   ├── estimators.py      # TI/BAR/MBAR 估算器
│   ├── profiles.py        # Lambda 剖面
│   └── convergence.py     # Time convergence + Bootstrap
├── analyzers/
│   ├── single.py          # FEPAnalyzer
│   └── multi.py           # FEPMultiEstimatorAnalyzer
└── reports/
    ├── html.py            # 单估算器报告
    ├── multi.py           # 多估算器报告（Tab切换）
    └── backend.py         # 多后端对比
```

## 核心功能

### 自动缓存

- **缓存文件**: JSON 格式，包含 metadata + methods + comparison
- **验证逻辑**: 自动检查 `estimators_used` 和 `temperature`
- **命中**: 配置匹配 → 秒级加载
- **未命中**: 重新运行并覆盖缓存

```bash
# 检查缓存内容
python3 -c "
import json
with open('fep_multi_estimator_cache.json', 'r') as f:
    data = json.load(f)
for method in ['TI', 'BAR', 'MBAR']:
    tc = data['methods'][method].get('time_convergence')
    bs = data['methods'][method].get('bootstrap')
    print(f'{method}: tc={bool(tc)}, bs={bool(bs)}')
"
```

### HTML 报告结构

1. **Summary Banner** - ΔΔG_bind ± error, 方法一致性警告
2. **Results Comparison** - TI/BAR/MBAR 对比表格, Box plots
3. **Estimator Tabs** (TI | BAR | MBAR):
   - Summary
   - Detailed Results (每次重复)
   - Lambda Profiles (ΔG(λ), ∂H/∂λ)
   - Convergence Diagnostics & Overlap:
     - Time-series convergence plot
     - Overlap matrix heatmap (仅 MBAR)
     - Bootstrap analysis (n=50, 80% sampling)

### 收敛性分析

**Time Convergence**: 将数据分段（10%, 20%, ..., 100%），计算累积 ΔG

**Bootstrap**: 重采样 n=50 次，每次抽取 80% frames，计算 ΔG 分布

**Overlap Matrix** (仅 MBAR):
- ✅ Good: min overlap > 0.3
- ⚠️ Acceptable: 0.1-0.3
- ❌ Poor: < 0.1（需要更多采样）

## 常见问题

### ❌ HTML 报告缺少收敛数据

**症状**: "Analysis was not performed"

**原因**: 快速模式（`skip_bootstrap=True, skip_time_convergence=True`）

**解决**:
```bash
rm fep_multi_estimator_cache.json
python test_multi_estimator.py --full
```

### ❌ 估算器结果不一致

**症状**: ΔG range > 1.0 kcal/mol, Agreement: poor

**检查**:
```python
mbar_result = results.methods['MBAR']
overlap_matrix = mbar_result.lambda_profiles['bound'][0]['overlap_matrix']
min_overlap = overlap_matrix.min()

if min_overlap < 0.1:
    print("Poor sampling - 增加采样时间或窗口数")
```

### ❌ 缓存不匹配

**症状**: "Cache mismatch: cache has ['TI', 'MBAR'], current has ['TI', 'BAR', 'MBAR']"

**自动处理**: 检测不匹配 → 自动重新运行

**手动强制**:
```bash
rm fep_multi_estimator_cache.json
python test_multi_estimator.py --no-cache
```

### ❌ ΔΔG 符号错误

**检查公式** (在 `analyzers/multi.py`):
```python
# 正确: ΔΔG = ΔG_bound - ΔG_unbound
# 如果 ΔΔG < 0: 结合更强（favorable）
# 如果 ΔΔG > 0: 结合更弱（unfavorable）
```

## XVG 文件解析

**实现**: `prism/fep/analysis/core/xvg_parser.py`

**支持格式**:
- `dhdl.xvg` (GROMACS 标准)
- `prod.xvg` (重命名版本)
- 自动处理 header (@, #, &)

**API**:
```python
from prism.fep.analysis.core.xvg_parser import parse_leg_data, find_xvg_file

# 查找文件
xvg_path = find_xvg_file(Path('bound1/window_00'))

# 解析 leg
datasets = parse_leg_data(
    leg_dir=Path('bound1'),
    leg_name='bound',
    estimator_name='MBAR',
    temperature=310.0,
    gmx_module=gmx,
    logger=logger,
)
```

## 绘图模块

**位置**: `prism/fep/analysis/plots.py`

**主要函数**:
- `build_lambda_profiles_html()` - ΔG(λ), ∂H/∂λ 曲线
- `build_time_convergence_html()` - Time convergence plot
- `build_bootstrap_html()` - Bootstrap distribution
- `build_overlap_matrix_html()` - Overlap heatmap

## 开发建议

### 已完成 ✅
- [x] 多估算器分析（TI, BAR, MBAR）
- [x] 自动缓存（JSON + 配置验证）
- [x] 快速/完整模式（skip_bootstrap, skip_time_convergence）
- [x] 进度条（tqdm）
- [x] HTML 报告（Tab 切换）
- [x] Time convergence + Bootstrap
- [x] Overlap matrix（MBAR only）
- [x] 多次重复支持

### 可能改进（按优先级）

1. **性能**: 并行运行估算器, 减少内存占用
2. **功能**: 更多估算器（FEP, US-tI）, 自定义 lambda 调度
3. **诊断**: 自动判断采样充分性, 推荐最优参数
4. **文档**: Jupyter notebook 教程, 视频演示

## 测试数据

**位置**: `examples/fep_output/1st_RUN4/WT-M10-re/`

**内容**: 3次重复（bound1/2/3, unbound1/2/3）, 32个lambda windows, prod.xvg格式

**运行**:
```bash
# 快速测试
python test_multi_estimator.py

# 完整分析
python test_multi_estimator.py --full

# 后台运行
python test_multi_estimator.py --full > run_full.log 2>&1 &
tail -f run_full.log
```

## 相关文档

- **模块**: `prism/fep/analysis/` (core, analyzers, reports)
- **测试**: `examples/fep_output/1st_RUN4/WT-M10-re/test_multi_estimator.py`
- **示例**: `examples/fep_output/1st_RUN4/WT-M10-re/fep_multi_estimator_report.html`

## 更新日志

**2026-03-29**:
- ✅ 更新模块结构（hierarchical: core/analyzers/reports）
- ✅ 添加自动缓存、快速/完整模式
- ✅ 更新 API 示例（FEPMultiEstimatorAnalyzer）
- ✅ 添加收敛诊断、overlap matrix
- ✅ 精简至 378 行（原 667 行）
