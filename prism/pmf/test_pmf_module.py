#!/usr/bin/env python3
  import sys
  import logging
  from pathlib import Path

  # 设置日志
  logging.basicConfig(level=logging.INFO,
                     format='%(asctime)s - %(levelname)s - %(message)s')

  def test_pmf_workflow():
      """测试PMF工作流"""

      try:
          from prism.pmf import PMFWorkflow, get_pmf_info

          # 显示模块信息
          info = get_pmf_info()
          print(f"PMF模块版本: {info['version']}")
          print(f"架构: {info['architecture']}")

          # 测试配置
          config = {
              'box': {
                  'box_distance': 1.2,
                  'pulling_distance': 2.0  # 测试用较小值
              },
              'smd': {
                  'pull_rate': 0.01,       # 测试用较快速率
                  'pull_k': 1000.0
              },
              'distance': {
                  'start': 0.3,
                  'end': 2.5
              },
              'equilibration': {
                  'em': {'nsteps': 10000},     # 测试用较少步数
                  'nvt': {'nsteps': 5000},     # 测试用较少步数  
                  'npt': {'nsteps': 10000}     # 测试用较少步数
              }
          }

          # 检查输入目录
          system_dir = "./gaff_model"  # 替换为你的实际路径
          if not Path(system_dir).exists():
              print(f"错误: 找不到系统目录 {system_dir}")
              print("请确保你有完整的PRISM MD结果")
              return False

          # 初始化工作流
          workflow = PMFWorkflow(system_dir, "./pmf_test", config)

          # 检查状态
          status = workflow.get_status()
          print(f"当前状态: {status['current_stage']}")
          print(f"下一步操作: {status['next_action']}")

          # 测试PMF系统构建
          print("\n=== 测试PMF系统构建 ===")
          build_results = workflow.build_pmf_system(equilibrate=False)  # 跳过平衡以加速测试

          print("PMF系统构建成功!")
          print(f"- 系统目录: {build_results['system_dir']}")
          print(f"- Z轴对齐: {build_results.get('alignment', {}).get('z_axis_aligned', 'N/A')}")

          # 测试SMD准备
          print("\n=== 测试SMD准备 ===")
          smd_results = workflow.prepare_smd(use_pmf_system=True)

          print("SMD准备成功!")
          print(f"- SMD目录: {smd_results['smd_dir']}")
          print(f"- 系统类型: {smd_results['system_type']}")
          print(f"- 运行脚本: {smd_results['run_script']}")

          # 检查生成的文件
          smd_dir = Path(smd_results['smd_dir'])
          required_files = ['smd.mdp', 'run_smd.sh', 'md.gro', 'topol.top']
          for file in required_files:
              if (smd_dir / file).exists():
                  print(f"✓ {file} 已生成")
              else:
                  print(f"✗ {file} 缺失")

          print("\n=== 测试完成 ===")
          print("如需运行SMD模拟，请执行:")
          print(f"cd {smd_results['smd_dir']} && bash run_smd.sh")

          return True

      except ImportError as e:
          print(f"导入错误: {e}")
          print("请确保PRISM已正确安装")
          return False

      except Exception as e:
          print(f"测试失败: {e}")
          logging.exception("详细错误信息:")
          return False

  def test_pmf_builder_direct():
      """直接测试PMF Builder"""

      try:
          from prism.pmf import pmf_builder

          print("=== 直接测试PMF Builder ===")

          builder = pmf_builder(
              md_results_dir="./gaff_model",  # 替换为实际路径
              output_dir="./pmf_builder_test",
              pulling_distance=2.0,
              box_distance=1.2
          )

          # 获取系统信息
          info = builder.get_system_info()
          print(f"构建器状态: {info.get('status', '未知')}")

          # 只测试结构提取（不做完整构建以加速测试）
          print("测试结构提取...")
          extraction_results = builder.extract_structures()

          print("结构提取成功!")
          print(f"- 蛋白质结构: {extraction_results['protein_structure']}")
          print(f"- 配体结构: {extraction_results['ligand_structure']}")
          print(f"- 复合物结构: {extraction_results['complex_structure']}")

          return True

      except Exception as e:
          print(f"PMF Builder测试失败: {e}")
          return False

  if __name__ == "__main__":
      print("开始PMF模块测试...")

      # 测试1: 完整工作流
      success1 = test_pmf_workflow()

      print("\n" + "="*50)

      # 测试2: 直接Builder
      success2 = test_pmf_builder_direct()

      if success1 and success2:
          print("\n所有测试通过! ✓")
      else:
          print("\n部分测试失败 ✗")
          sys.exit(1)