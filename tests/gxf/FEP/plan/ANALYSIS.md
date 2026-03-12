# PRISM-FEP 分析工具

本文档描述 FEP 结果分析工具和可视化报告生成。

## 1. XVG 文件解析模块 (`fep/analysis/xvg_parser.py`)

### 1.1 功能需求
| 项 | 说明 |
|----|------|
| 输入 | GROMACS 输出的 `dhdl.xvg` 文件 |
| 输出 | `(lambdas, dhdl, errors)` 元组 |

### 1.2 设计思路

**XVG 文件格式**：
| 行类型 | 说明 |
|--------|------|
| 注释行 | 以 `#` 或 `@` 开头 |
| 数据行 | `lambda dH/dl error` |
| 多列 | 对应不同 λ 窗口 |

### 1.3 实现思路

```python
def parse_dhdl(xvg_file：str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]：
    """解析 GROMACS dhdl.xvg 文件"""
    data = np.loadtxt(xvg_file, comments=['#', '@'])
    
    if data.shape[1] >= 3：
        return data[：, 0], data[：, 1], data[：, 2]  # 有误差列
    else：
        return data[：, 0], data[：, 1], np.zeros_like(data[：, 0])  # 无误差
```

---

## 2. 参考分析流程（当前可跑基线）

基线脚本见 `tests/gxf/FEP/ref/gxf_gmx_mdp_re/calculate.sh`，流程要点：
| 步骤 | 说明 |
|------|------|
| 目录结构 | `bound/` 与 `free/` 下各有 `repeat*/lambda_XX/` |
| 数据汇总 | `lambda_XX/dhdl.xvg` 复制到 `result/` |
| BAR 计算 | `gmx bar -f lambda*.xvg -o bar.xvg -oi barint.xvg -oh histogram.xvg` |
| 汇总输出 | `summary.dat`（记录 bound/free 的 total 值） |

该流程目前是“真实可跑”的默认基线，后续再扩展 alchemlyb/可视化即可。

---

## 3. 自由能估计器 (`fep/analysis/estimators.py`)

### 3.1 功能需求
| 项 | 说明 |
|----|------|
| 输入 | dH/dl 数据列表 |
| 输出 | 自由能差 ΔG 与误差估计 |

### 3.2 设计思路

使用 alchemlyb 库进行自由能分析（可选增强）：
| 步骤 | 说明 |
|------|------|
| 解析 | 读取 XVG 获取 dH/dl |
| 计算 | BAR 或 MBAR |
| 误差 | 输出不确定度估计 |

### 3.3 BAR/MBAR 对比
| 方法 | 适用场景 | 特点 |
|------|----------|------|
| BAR | 相邻 λ 窗口 | 精度较高，计算量适中 |
| MBAR | 全部 λ 窗口 | 精度最高，计算量较大 |

---

## 4. 原子映射可视化 (`prism/fep/visualize/`)

### 4.1 模块结构

**已完成的可视化模块**：
```
prism/fep/visualize/
├── __init__.py              # 导出辅助函数
├── molecule.py              # 分子处理工具
│   ├── pdb_to_mol()                    # PDB → RDKit Mol
│   ├── assign_bond_orders_from_mol2()  # 从 mol2 校正键级
│   └── prepare_mol_with_charges_and_labels()  # 添加电荷和标签
├── highlight.py             # 高亮颜色定义
│   ├── COMMON_COLOR = rgb(204, 229, 77)      # 绿色
│   ├── TRANSFORMED_COLOR = rgb(255, 77, 77)   # 红色
│   └── SURROUNDING_COLOR = rgb(77, 153, 255)  # 蓝色
├── mapping.py               # PNG 可视化（FEbuilder 风格）
│   └── visualize_mapping_png()  # 生成 PNG 映射图
└── html.py                  # HTML 交互式可视化
    └── visualize_mapping_html()  # 生成 HTML 页面
```

### 4.2 PNG 可视化 (FEbuilder 风格)

**功能特性**：
| 特性 | 说明 |
|------|------|
| MCS 对齐 | 使用 RDKit `rdFMCS` 进行最大公共子结构对齐 |
| 键级校正 | 从 mol2 文件读取正确的键级信息 |
| 电荷标注 | 显示所有原子的电荷（包括氢原子） |
| 颜色高亮 | Common=绿色, Transformed=红色, Surrounding=蓝色 |
| 图例 | 显示各类型原子数量和说明 |
| 边距控制 | 可选调整图像边距 |

**使用示例**：
```python
from prism.fep.visualize import visualize_mapping_png

visualize_mapping_png(
    mapping,
    pdb_a='25.pdb',
    pdb_b='36.pdb',
    mol2_a='25.mol2',  # 用于键级校正
    mol2_b='36.mol2',
    output_path='25-36_mapping.png',
    edge=0.15,  # 边距
    legends=('Ligand 25', 'Ligand 36')
)
```

### 4.3 HTML 交互式可视化 (设计要求)

**交互功能需求**：

| 功能 | 描述 | 优先级 |
|------|------|--------|
| **Hover 显示电荷** | 鼠标悬停在原子上时，显示该原子的电荷信息 | 高 |
| **对应原子高亮** | Hover 时，如果对面分子有对应原子，同时高亮显示 | 高 |
| **电荷标签开关** | 提供开关控制是否在图像上直接显示电荷标签 | 中 |
| **导出图片** | 支持将当前视图导出为 PNG 图片 | 中 |
| **交互模式开关** | 配置项控制是否生成交互式 HTML 或静态 PNG | 高 |

**设计建议**：

1. **Hover Tooltip 内容**：
   - 当前原子名称（如 C1, H2）
   - 当前原子电荷（如 +0.123, -0.456）
   - 原子类型（如 CG2R66, NG321）
   - 分类信息（Common/Transformed/Surrounding）
   - 如果有对应原子：显示对应原子的名称和电荷

2. **电荷标签显示策略**：
   - **默认关闭**：避免图像过于拥挤
   - **开启时**：在小分子（<30 原子）上显示所有电荷
   - **大分子**：只显示差异电荷（Surrounding 和 Transformed）

3. **导出功能**：
   - 导出为高清 PNG（300 DPI）
   - 保留当前高亮和标签状态
   - 文件名自动添加时间戳

4. **静态 vs 交互式选择**：
   - 通过配置参数 `--interactive` 控制
   - 静态 PNG：用于报告、论文、快速查看
   - 交互式 HTML：用于深入分析、调试、演示

**HTML 特有功能（静态 PNG 不具备）**：

| 功能 | 优势 |
|------|------|
| **动态筛选** | 点击图例可显示/隐藏特定类型的原子 |
| **原子列表** | 侧边栏显示所有原子的详细信息表 |
| **搜索定位** | 搜索原子名称快速定位并高亮 |
| **对比模式** | 并排显示参考和突变分子的详细信息 |
| **数据导出** | 导出原子映射数据为 CSV/JSON 格式 |
| **响应式设计** | 自适应不同屏幕尺寸，移动端友好 |

### 4.4 CLI 接口

```bash
# 生成 PNG 可视化（默认）
prism visualize-mapping --case 25-36 --output mapping.png

# 生成交互式 HTML
prism visualize-mapping --case 25-36 --output mapping.html --interactive

# 在 HTML 中显示电荷标签
prism visualize-mapping --case 25-36 --output mapping.html --show-charges

# 自定义图例
prism visualize-mapping --case 25-36 --legends "Reference" "Mutant"
```

### 4.5 可视化输出

| 输出类型 | 文件格式 | 用途 | 位置 |
|---------|---------|------|------|
| 静态图像 | PNG | 报告、论文、快速查看 | `output/mapping.png` |
| 交互式 | HTML | 深入分析、调试、演示 | `output/mapping.html` |
| 映射数据 | JSON/CSV | 数据交换、二次分析 | `output/mapping.json` |

### 4.6 技术实现要点

**PNG 生成**：
- 使用 RDKit `MolDrawOptions` 设置绘制选项
- 通过 `rdFMCS` 计算最大公共子结构进行对齐
- 使用 `AssignBondOrdersFromTemplate` 从 mol2 校正键级
- 电荷信息通过 `atomNote` 属性添加

**HTML 生成**：
- 图像转为 base64 编码嵌入 HTML（单文件部署）
- 使用 JavaScript 实现 hover 交互和原子高亮
- CSS 实现响应式布局和动画效果
- 支持 Chart.js 或 D3.js 进行数据可视化

---

## 5. FEP 结果分析报告 (`fep/analysis/report.py`)

### 5.1 功能需求（待实现）

生成完整的 FEP 计算报告，包含：
| 内容 | 说明 |
|------|------|
| 计算信息 | 案例与元数据 |
| 原子映射 | 嵌入 PNG/HTML 可视化 |
| ΔG 曲线 | ΔG vs λ |
| 能量分布 | λ 窗口直方图 |
| 误差评估 | 置信区间 |
| 收敛性分析 | HAMADA 统计量 |

### 5.2 设计思路

| 步骤 | 说明 |
|------|------|
| Plotly | 交互式图表（ΔG 曲线、能量分布） |
| RDKit | 2D 分子结构图（原子映射） |
| HTML | 统一报告模板 |
| 多窗口 | λ 窗口对比分析 |

### 5.3 CLI 接口

```bash
# 生成完整 HTML 报告
prism report --output-dir fep_system --report-type html

# 分析 FEP 结果
prism analyze --xvg dhdl.xvg --method bar

# 生成对比报告
prism report --compare system1 system2 --output comparison.html
```

---

## 5. 参考资料
| 资源 | 地址 |
|------|------|
| alchemlyb | https://alchemlyb.readthedocs.io/ |
| Plotly | https://plotly.com/python/ |
| RDKit | https://www.rdkit.org/ |
