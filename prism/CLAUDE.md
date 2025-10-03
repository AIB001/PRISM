### rules
1.我们的对话采取中文对话，但是脚本的书写避免出现中文和emoji
2.update2_128.sh为更新到服务器的脚本，作文件整理时不可删除
3.对于文件的命名需要服务于项目，遵循一致性避免出现过于体现修改历史的命名
4.你不需要去测试脚本，我要放到服务器上来进行测试，我会把结果反馈给你你可以写测试脚本给我拿在服务器上测试，但是注意及时清除测试脚本，保证项目的干净

### 文件命名规范
1. **避免冗余后缀**: 文件名应简洁明确，避免如 `_only`, `_main`, `_final` 等体现修改历史的后缀
2. **功能导向命名**: 文件名应直接反映其主要功能，如 `pmf_remodel.py` 而非 `pmf_remodel_only.py`
3. **保持一致性**: 同类文件使用相同的命名模式和风格
4. **服务项目需求**: 文件名应便于用户理解和使用，避免过于技术化的历史标记

### 示例命名改进
- `pmf_remodel_only.py` → `pmf_remodel.py` (移除冗余的 `_only` 后缀)
- 保持核心功能文件的简洁命名，如 `test_pmf_remodeling.py`, `pmf_workflow_checker.py`
### importance
我们的脚本目前为本地环境，是通过本地修改完成后upload到服务器上进行测试的，因此只需要根据要求修改脚本，测试环节我会更新在服务器上自行测试

## PMF_tutorial
### tutorial格式
写成ipynb的格式，利用jupyter notebook向用户呈现我们的流程
### 要求
1.tutorial里避免出现中文和emoji
2.明显的I/O的显示
3.说明清楚我们的路径上的处理

## PMF Module API Design

### Two-API Architecture
PMF模块采用双API设计，只保留两个核心启动入口，避免API混乱:

#### 1. Remodeling API - 系统重建
**用途**: 将MD结果转换为PMF优化系统
**入口**: `PMFBuilder`

```python
from prism.pmf import PMFBuilder

# 初始化
builder = PMFBuilder(
    md_results_dir="./gromacssim",
    output_dir="./pmf_system",
    config=config_dict  # optional
)

# 执行重建
# equilibrate=True: 自动运行平衡 (本地测试)
# equilibrate=False: 仅生成脚本 (服务器部署)
results = builder.build(equilibrate=True)
```

#### 2. Runner API - PMF工作流
**用途**: 执行完整PMF计算流程
**入口**: `PMFRunner` 或 `run_pmf_workflow()`

```python
from prism.pmf import PMFRunner

# 方法1: 使用Runner类
runner = PMFRunner(config="pmf_config.yaml")
results = runner.run_complete_workflow(
    md_system_dir="./pmf_system",
    output_dir="./pmf_results"
)

# 方法2: 使用便捷函数
from prism.pmf import run_pmf_workflow
results = run_pmf_workflow(
    md_system_dir="./pmf_system",
    output_dir="./pmf_results",
    config="pmf_config.yaml"
)
```

### API使用规则

1. **禁止直接调用内部组件**
   - 不要直接使用 `PMFSystem`, `PMFWorkflow`, `SMDManager`, `UmbrellaManager` 等
   - 这些是内部实现，只通过两个主API调用

2. **只修改两个入口脚本**
   - 所有用户接口改动只在 `PMFBuilder` 和 `PMFRunner` 中进行
   - 保持入口简单清晰

3. **示例脚本规范**
   - examples/ 中的脚本只能使用这两个API
   - 不展示内部组件的直接使用

4. **配置文件优先**
   - 复杂参数通过YAML配置传递
   - API保持简洁

### 典型工作流

```python
# Step 1: Remodel MD system
from prism.pmf import PMFBuilder

builder = PMFBuilder("./md_results", "./pmf_system")
remodel_results = builder.build(equilibrate=True)

# Step 2: Run PMF calculations
from prism.pmf import run_pmf_workflow

pmf_results = run_pmf_workflow(
    md_system_dir="./pmf_system/GMX_PMF_SYSTEM",
    output_dir="./pmf_results"
)

print(f"Binding energy: {pmf_results['binding_energy']['value']:.2f} kcal/mol")
```

## PMF_code
### rebuild
主要建模部分已经实现高度的自动化
#### relaxation
这个模块的relaxation步骤用户可以选择local模式生成gromacs运行的脚本，也可以根据yaml的slurm文件配置来生成slurm脚本
