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

## PMF_code
### rebuild
主要建模部分已经实现高度的自动化
#### relaxation
这个模块的relaxation步骤用户可以选择local模式生成gromacs运行的脚本，也可以根据yaml的slurm文件配置来生成slurm脚本
