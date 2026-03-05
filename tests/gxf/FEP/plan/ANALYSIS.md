# PRISM-FEP 分析工具

本文档描述 FEP 结果分析工具和可视化报告生成。

## 1. XVG 文件解析模块 (`fep/analysis/xvg_parser.py`)

### 1.1 功能需求
- **输入**: GROMACS 输出的 `dhdl.xvg` 文件
- **输出**: (lambdas, dhdl, errors) 元组

### 1.2 设计思路

**XVG 文件格式**：
- 前几行是注释 (以 # 或 @ 开头)
- 数据部分: lambda dH/dl error
- 多列数据对应不同的 lambda 窗口

### 1.3 实现思路

```python
def parse_dhdl(xvg_file: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """解析 GROMACS dhdl.xvg 文件"""
    data = np.loadtxt(xvg_file, comments=['#', '@'])
    
    if data.shape[1] >= 3:
        return data[:, 0], data[:, 1], data[:, 2]  # 有误差列
    else:
        return data[:, 0], data[:, 1], np.zeros_like(data[:, 0])  # 无误差
```

---

## 2. 参考分析流程（当前可跑基线）

基线脚本见 `tests/gxf/FEP/ref/gxf_gmx_mdp_re/calculate.sh`，流程要点：
- 目录结构：`bound/` 与 `free/` 下各有 `repeat*/lambda_XX/`
- 每个 repeat 将 `lambda_XX/dhdl.xvg` 汇总到 `result/` 目录
- 使用 `gmx bar` 计算每个 repeat 的 ΔG：
  - `gmx bar -f lambda*.xvg -o bar.xvg -oi barint.xvg -oh histogram.xvg`
- 汇总输出：`summary.dat`（记录 bound/free 各 repeat 的 total 值）

该流程目前是“真实可跑”的默认基线，后续再扩展 alchemlyb/可视化即可。

---

## 3. 自由能估计器 (`fep/analysis/estimators.py`)

### 3.1 功能需求
- **输入**: dH/dl 数据列表
- **输出**: 自由能差 ΔG 和误差估计

### 3.2 设计思路

使用 alchemlyb 库进行自由能分析（可选增强）：
1. 解析 XVG 文件获取 dH/dl 数据
2. 使用 BAR 或 MBAR 方法计算自由能差
3. 提供误差估计

### 3.3 BAR 方法

**Bennett Acceptance Ratio (BAR)**：
- 适用于相邻 λ 窗口的自由能差计算
- 精度较高，计算量适中
- 最常用的方法

### 3.4 MBAR 方法

**Multistate Bennett Acceptance Ratio (MBAR)**：
- 同时使用所有 λ 窗口的数据
- 精度最高，计算量较大
- 适用于高质量数据

---

## 4. HTML 报告生成 (`fep/analysis/report.py`)

### 4.1 功能需求

生成 HTML 可视化报告（后续增强项），包含：
- 测试案例信息和元数据
- 原子映射可视化（2D 结构图）
- 自由能变化曲线（ΔG vs λ）
- λ 窗口能量分布直方图
- 错误估计和置信区间

### 4.2 设计思路

1. 使用 Plotly.js 生成交互式图表
2. 使用 RDKit 生成 2D 分子结构图
3. 高亮显示 transformed/surrounding/common atoms
4. 支持多个 λ 窗口的对比展示

### 4.3 CLI 接口

```bash
prism fep-report --output-dir fep_system --report-type html
```

### 4.4 可视化元素

**自由能变化曲线**：
- X 轴: λ (0.0 → 1.0)
- Y 轴: ΔG (kJ/mol)
- 标注误差范围

**原子映射图**：
- RDKit 2D 分子结构
- 高亮显示：
  - transformed atoms（红色）
  - surrounding atoms（蓝色）
  - common atoms（绿色）

**能量分布**：
- 各 λ 窗口的 dH/dλ 分布直方图
- 显示采样充分性

### 4.5 输出文件

- `fep_report.html` - 主报告
- `plots/` - 图表数据
- `structures/` - 分子结构图

---

## 5. 参考资料
- alchemlyb: https://alchemlyb.readthedocs.io/
- Plotly: https://plotly.com/python/
- RDKit: https://www.rdkit.org/
