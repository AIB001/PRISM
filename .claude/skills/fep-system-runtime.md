---
name: fep-system-runtime
description: FEP 系统运行与执行调试
type: fep
---

# FEP 系统运行与执行调试

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

## 核心职责

检查 PRISM FEP 系统的运行时状态，验证脚本执行、GROMACS 命令、GPU 利用率、lambda 窗口执行策略、生产 MD 运行等运行时配置。

## ⚠️ 核心原则（CRITICAL）

### GPU 使用原则（最高优先级）

**为什么不使用 GPU？未经我允许的情况下要尽量使用 GPU！**

**必须使用 GPU**：
- ✅ **未经允许的情况下要尽量使用 GPU**
- ✅ **FEP 模板应该跟正常跑 MD 差不多**
- ✅ **本机有 GPU，必须充分利用**
- ✅ **默认启用 GPU，不要等待用户明确要求**
- ❌ **不要禁用 GPU 或使用 CPU-only 模式**
- ❌ **不要问"要不要用GPU"，直接用！**

**遇到性能问题时首先检查**：
1. GPU 是否启用？（`nvidia-smi` 检查利用率）
2. MDP 参数是否正确？（特别是 `calc-lambda-neighbors`）
3. 线程配置是否优化？（`-ntmpi 1 -ntomp 8`）
4. 不要怀疑 GROMACS 或常用参数值，优先认为是建模/topology/MDP参数的问题

**GPU 配置要求**：
```bash
# 必须包含以下参数
-nb gpu         # 非键相互作用使用 GPU
-bonded gpu     # 键相互作用使用 GPU
-pme gpu        # PME 使用 GPU
-update gpu     # 约束更新使用 GPU

# 线程配置（CRITICAL）
-ntmpi 1 -ntomp 8   # 单MPI进程，多OpenMP线程（GPU优化）
# ❌ 不要使用 -ntmpi 2 -ntomp 2（会降低GPU利用率）
# ❌ 不要使用 -ntmpi 4 -ntomp 2（这是CPU模式配置）
```

**性能基准**：
- 单 window GPU 利用率：**>80%**
- 单 window 性能：**>200 ns/day**
- 4 并行 windows GPU 利用率：**>70%**
- 4 并行 windows 性能：**>170 ns/day**

### 问题排查优先级

**优先认为是建模/topology/MDP参数的问题**：
- ✅ 优先检查：力场参数、拓扑文件、原子类型
- ✅ 优先检查：MDP 参数设置（特别是 `calc-lambda-neighbors`）
- ❌ 不要轻易认为是 GROMACS 软件问题
- ❌ 不要轻易认为是常用参数值的问题

**常见性能问题根源**（按优先级）：
1. `calc-lambda-neighbors = -1` → 计算所有 lambda 状态，性能暴跌（27-73倍差异）
2. GPU 未启用 → CPU 模式性能低下
3. 执行模式选择不当 → 未充分利用多 GPU
4. 线程配置错误 → 使用 `-ntmpi 2 -ntomp 2` 而非 `-ntmpi 1 -ntomp 8`

**实际性能数据**（calc-lambda-neighbors的影响）：
- ❌ `calc-lambda-neighbors = -1`：
  - GPU利用率：0-4%
  - 性能：3-8 ns/day
  - 原因：计算所有32个lambda状态之间的dH/dλ，计算量O(n²)
- ✅ `calc-lambda-neighbors = 1`：
  - 单window：GPU利用率89%，性能221 ns/day
  - 4并行windows：GPU利用率70-85%，性能170-191 ns/day
  - 原因：只计算相邻lambda状态，符合WHAM分析需求

## Lambda 窗口运行机制

### 动态窗口创建（重要）

**两阶段工作流**：

1. **Build 阶段**（搭建阶段）：
   - ✅ 生成 MDP 文件到 `mdps/` 目录
   - ✅ 生成 `run_fep.sh` 主脚本
   - ✅ 生成 `lambda_schedule.json`
   - ❌ **不创建** `window_XX/` 目录

2. **Run 阶段**（执行阶段）：
   - ✅ **动态创建** `window_XX/` 目录
   - ✅ 每个 lambda 窗口在执行时创建
   - ✅ 代码位置：`prism/fep/modeling/script_writer.py` 第214-219行

**验证命令**：
```bash
# Build 完成后（应该返回 0）
find bound -type d -name "window_*" | wc -l

# 运行后（应该返回 lambda_windows 数量）
bash bound/run_fep.sh bound
find bound -type d -name "window_*" | wc -l
```

## 执行模式选择（CRITICAL）

PRISM FEP 支持两种执行模式，通过配置简单切换：

### Standard 模式（默认）

**特点**：
- 每个 lambda 窗口独立运行
- 多个窗口并行在不同 GPU 上
- 适合：快速探索、资源有限

**资源语义**：
- 4 个 GPU = 4 个窗口并行（各自独立）
- 每个 GPU 运行 1-2 个窗口

**示例**（32 窗口，4 GPU）：
```bash
# GPU 0: windows 0-7
# GPU 1: windows 8-15
# GPU 2: windows 16-23
# GPU 3: windows 24-31

# 每个窗口 10 ns，总时间 ≈ 10 ns（并行）
```

**脚本**：`run_prod_standard.sh`

### Repex 模式（副本交换）

**特点**：
- Lambda 窗口间进行副本交换
- 提高采样效率
- 适合：精确计算、充分采样

**资源语义**：
- 4 个 GPU = 一个 multidir RE 作业共享 4 GPU
- 所有窗口一起做副本交换

**命令**：
```bash
gmx_mpi mdrun -multidir window_00,window_01,...,window_31 \
    -replex 1000 -nb gpu -bonded gpu -pme gpu
```

**脚本**：`run_prod_repex.sh`

### 模式切换

**配置文件**（`fep.yaml`）：
```yaml
execution:
  mode: standard   # 或 repex
```

**脚本生成**：
- 自动生成两种脚本：`run_prod_standard.sh`、`run_prod_repex.sh`
- 统一入口 `run_prod.sh` 根据 `execution.mode` 调用对应脚本

**关键区别**：
| 特性 | Standard | Repex |
|------|----------|-------|
| 窗口运行 | 独立并行 | 联合交换 |
| GPU 分配 | 每窗口独占 | 共享 multidir |
| 采样效率 | 标准 | 更高 |
| 实现复杂度 | 简单 | 中等 |

### 共享特性（两种模式相同）

**EM/NVT/NPT 平衡**：
- 两种模式都逐窗口独立跑
- 不使用副本交换

**分析层**：
- 最终每个 window 都产出 `dhdl.xvg`
- Analysis 不区分 standard/repex
- 完全共用分析逻辑

**建模层**：
- Mapping、hybrid topology、MDP 主体完全相同
- 只在 production 执行时分叉

## GPU 性能调优（CRITICAL）

### 性能瓶颈诊断

**症状：GPU 利用率极低（0-4%）**
```
nvidia-smi
# GPU-Util: 0-4%
# 性能：3-8 ns/day
```

**根本原因**：
- MDP 参数 `calc-lambda-neighbors = -1`
- 导致 GROMACS 计算所有 32 个 lambda 状态之间的 dH/dλ
- 计算量巨大：O(n²) 复杂度

**解决方案**：
```bash
# 修改 MDP 文件
calc-lambda-neighbors = 1  # 只计算相邻 lambda 状态

# 性能提升：
# 单 window: 221 ns/day, GPU 利用率 89%
# 4 并行: 170-191 ns/day, GPU 利用率 70-85%
# 提升倍数：27-73 倍
```

**为什么 calc-lambda-neighbors = 1 足够**：
- WHAM 分析只需要相邻 lambda 状态的 dH/dλ
- 这是 GROMACS 的最佳实践
- 不影响 FEP 计算精度

### 多 GPU 并行策略（Standard 模式）

**目标**：使用 4 个 GPU 并行运行多个 lambda 窗口

**配置原则**：
- 每个 GPU 运行 1-2 个窗口
- 每个窗口获得更多 CPU 资源
- 所有窗口并行运行，快速完成

**脚本示例**：
```bash
#!/bin/bash
# run_prod_standard.sh

NUM_GPUS=4
WINDOWS_PER_GPU=8

for gpu_id in $(seq 0 $((NUM_GPUS-1))); do
    start_window=$((gpu_id * WINDOWS_PER_GPU))
    end_window=$((start_window + WINDOWS_PER_GPU - 1))

    for window_idx in $(seq $start_window $end_window); do
        window_dir=$(printf "window_%02d" $window_idx)

        cd bound/$window_dir
        gmx mdrun -deffnm prod -nb gpu -bonded gpu -pme gpu \
            -ntmpi 1 -ntomp 8 -gpu_id $gpu_id &
        cd ../..
    done
done

wait
```

**性能预期**：
- 单窗口单 GPU：221 ns/day
- 4 并行窗口（4 GPU）：170-191 ns/day
- 总吞吐量：680-764 ns/day

### 副本交换配置（Repex 模式）

**目标**：所有窗口共享 4 GPU 进行副本交换

**配置原则**：
- 使用 `gmx_mpi mdrun -multidir`
- 所有窗口在一个 MPI 作业中
- 定期交换 lambda 状态

**脚本示例**：
```bash
#!/bin/bash
# run_prod_repex.sh

NUM_GPUS=4

# 收集所有窗口目录
window_dirs=""
for window_idx in $(seq 0 31); do
    window_dir=$(printf "window_%02d" $window_idx)
    window_dirs="${window_dirs}${window_dir},"
done
window_dirs=${window_dirs%,}  # 去掉最后的逗号

# 运行副本交换
cd bound
mpirun -np $NUM_GPUS gmx_mpi mdrun \
    -multidir ${window_dirs} \
    -deffnm prod \
    -replex 1000 \
    -nb gpu -bonded gpu -pme gpu \
    -ntomp 8
```

**性能预期**：
- GPU 利用率：70-85%
- 单窗口性能：170-191 ns/day
- 采样效率：比 standard 更高

## 运行前验证

### GPU 配置验证（必须！）

**在运行任何 FEP 计算前，必须验证 GPU 配置**：

```bash
cd gaff_test_output/GMX_PROLIG_FEP

# 1. 检查 MDP 文件中的关键参数
echo "=== 关键 MDP 参数 ==="
grep -E "calc-lambda-neighbors|free-energy" bound/mdps/prod_00.mdp

# 2. 检查执行模式
echo "=== 执行模式 ==="
grep "execution:" fep.yaml

# 3. 检查脚本中的 GPU 设置
echo "=== 脚本 GPU 设置 ==="
grep -E "gpu|GPU" bound/run_prod_*.sh | head -10

# 4. 验证 GPU 可用性
echo "=== GPU 可用性 ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
```

**常见问题**：
- ❌ `calc-lambda-neighbors = -1` → 改为 1
- ❌ 脚本中没有 `-nb gpu -bonded gpu -pme gpu` → 添加
- ❌ 使用 `-ntmpi 4 -ntomp 2`（CPU 模式）→ 改为 `-ntmpi 1 -ntomp 8`

### 快速验证脚本

```bash
cd gaff_test_output/GMX_PROLIG_FEP

# 1. 检查 GROMACS 可用
gmx --version

# 2. 检查 topology 语法
gmx check -f common/hybrid/hybrid.gro -p bound/topol.top

# 3. 尝试生成第一个窗口的 TPR
cd bound
gmx grompp -f mdps/em.mdp -c input/conf.gro -p topol.top -o em_test.tpr -maxwarn 2

# 4. 如果成功，检查 TPR
gmx check -s em_test.tpr

# 5. 清理测试文件
rm -f em_test.tpr *.cpt
```

**常见错误**：
- `Atom type XXX not found` → 力场参数缺失
- `No default bond type` → 键级参数缺失
- `number of coordinates in file` → 坐标数不匹配

## 脚本调用与问题解决（CRITICAL）

### 必须调用生成的脚本

**原则**：
- ✅ **必须调用 PRISM 生成的脚本**（`localrun.sh`, `run_fep.sh`, `run_prod_standard.sh` 等）
- ❌ **不要自己编写新的运行脚本**
- ✅ **从 MDP 生成角度解决问题**，而不是绕过生成的脚本

**为什么**：
- 生成的脚本已经包含了正确的 GPU 配置、线程配置、并行策略
- 生成的脚本处理了 resume 逻辑、错误处理、依赖关系
- 自己写脚本容易遗漏关键配置，导致性能问题

**正确做法**：
```bash
# ✅ 正确：调用生成的脚本
cd GMX_PROLIG_FEP/bound
bash localrun.sh                    # 单 leg 完整流程
bash run_prod_standard.sh           # 只运行 production 阶段

# ❌ 错误：自己写脚本
for window in window_*; do
    cd $window
    gmx mdrun -deffnm prod ...       # 容易遗漏关键参数
    cd ..
done
```

### 从 MDP 生成角度解决问题

**问题排查流程**：
1. **首先检查 MDP 文件**：
   ```bash
   grep -E "calc-lambda-neighbors|free-energy" bound/mdps/prod_00.mdp
   ```

2. **如果 MDP 参数错误**：
   - 修改 `prism/fep/gromacs/mdp_templates.py` 中的模板
   - 或修改 `prism/fep/modeling/mdp_generator.py` 中的生成逻辑
   - 重新 build FEP 系统

3. **不要在运行脚本中绕过问题**：
   - ❌ 不要用 `sed` 临时修改 MDP
   - ❌ 不要在运行脚本中硬编码参数
   - ✅ 从源头修复，重新生成正确的文件

### 多 GPU 并行策略详解

**用户明确要求**：
- "try using only 4 GPU so that each window get more CPUs"
- "每个window并行呢？这样的话，比如每个windows只跑10个纳秒，也能很快跑完"

**实现方式**：
1. **4 个 GPU 并行运行多个 lambda 窗口**
2. **每个 GPU 运行 1-2 个窗口**（不是每个窗口用4个GPU）
3. **每个窗口获得更多 CPU 资源**（通过 `-ntomp 8`）
4. **所有窗口并行执行**，快速完成

**GPU 分配策略**：
```bash
# 假设有 11 个 lambda windows，4 个 GPU
# GPU 0: window_00, window_01
# GPU 1: window_02, window_03
# GPU 2: window_04, window_05
# GPU 3: window_06, window_07
# (依此类推...)

# 每个 window 使用：
# - CUDA_VISIBLE_DEVICES=<gpu_id>  # 指定GPU
# -ntmpi 1 -ntomp 8                 # 单MPI进程，8个OpenMP线程
# -nb gpu -bonded gpu -pme gpu      # GPU加速
```

**为什么这样配置**：
- 每个 window 独占一个 GPU，避免 GPU 间通信开销
- `-ntmpi 1` 避免 MPI 进程间通信开销
- `-ntomp 8` 充分利用 CPU 核心，保持 GPU 忙碌
- 并行执行多个 windows，总吞吐量最大化

### 大致性能评估

**快速判断性能是否正常**：
```bash
# 检查 GPU 利用率
nvidia-smi
# 正常：70-89% GPU 利用率
# 异常：<10% GPU 利用率（检查 calc-lambda-neighbors）

# 检查单个 window 性能
grep "Performance" bound/window_00/prod.log
# 正常：>200 ns/day
# 异常：<50 ns/day（可能未使用 GPU）

# 检查 4 并行 windows 性能
# 总吞吐量应该是单 window 的 3-4 倍
# 正常：>600 ns/day（4个窗口并行）
# 异常：<200 ns/day（未真正并行）
```

**性能基准对照表**：
| 配置 | GPU利用率 | 性能 (ns/day) | 状态 |
|------|-----------|---------------|------|
| calc-lambda-neighbors = -1 | 0-4% | 3-8 | ❌ 性能灾难 |
| 单 window (正确配置) | 80-89% | 200-221 | ✅ 正常 |
| 4 并行 windows | 70-85% | 170-191 per window | ✅ 正常 |
| 4 并行总吞吐 | 70-85% | 680-764 total | ✅ 正常 |

### 脚本生成原则（核心开发原则）

**要尽量把脚本生成对，不是手动改！**

**开发优先级**：
1. ✅ **从源头修复**：修改 `prism/fep/modeling/script_writer.py` 或 `prism/fep/gromacs/mdp_templates.py`
2. ✅ **重新生成**：删除旧输出，重新 build，生成正确的脚本
3. ❌ **避免手动修改**：不要用 `sed` 或手动编辑生成的脚本
4. ❌ **避免临时补丁**：不要在运行脚本中硬编码修复

**为什么**：
- 手动改的脚本在下次 build 时会被覆盖
- 临时补丁无法解决根本问题，会反复出现
- 从源头修复可以确保所有用户都受益

**最终目标：用户能直接用**
```bash
# 用户应该只需要：
python tests/gxf/FEP/unit_test/42-38/test_run_fep.py --forcefield amber14sb
cd GMX_PROLIG_FEP/bound
bash localrun.sh  # 直接运行，无需修改

# 而不是：
# 1. 修改 MDP 文件
# 2. 修改运行脚本
# 3. 手动调整 GPU 参数
# 4. 然后才能运行
```

**问题修复流程**：
1. **发现问题**（如 GPU 利用率低）
2. **定位源头**（检查 MDP 模板、脚本生成代码）
3. **修复源头**（修改 Python 代码）
4. **验证修复**（重新 build，检查生成的脚本）
5. **提交代码**（确保所有用户都受益）

### GPU 回退策略（重要）

**原则：有 GPU 的系统不要回退到 CPU**

**智能回退逻辑**：
```bash
# 检测系统是否有 GPU
if nvidia-smi &> /dev/null; then
    # 有 GPU：不要回退到 CPU，GPU 失败就报错
    NUM_GPUS=$(nvidia-smi --list-gpus | wc -l)
    if [ $NUM_GPUS -ge 1 ]; then
        # GPU 可用，强制使用 GPU
        gmx mdrun -nb gpu -bonded gpu -pme gpu ... || {
            echo "Error: GPU run failed but GPU is available"
            echo "This indicates a configuration problem, not a hardware problem"
            exit 1
        }
    fi
else
    # 无 GPU：允许 CPU 回退
    gmx mdrun -nb gpu -bonded gpu -pme gpu ... || \
        gmx mdrun -ntomp 10 ...
fi
```

**系统特定配置**：
- **209 系统**（多 GPU 服务器）：
  - 有大量 GPU 可用
  - **不要回退到 CPU**，GPU 失败应该报错
  - 强制诊断问题（力场、topology、MDP 参数）

- **其他用户系统**：
  - 可能没有 GPU
  - **支持 CPU 回退**，确保兼容性
  - 但优先尝试 GPU

**脚本生成建议**：
```python
# prism/fep/modeling/script_writer.py
def detect_gpu_available():
    """检测系统是否有 GPU"""
    import subprocess
    try:
        subprocess.run(['nvidia-smi'], check=True,
                      capture_output=True, timeout=5)
        return True
    except:
        return False

def generate_mdrun_command(has_gpu, gpu_available_on_system):
    """根据系统配置生成 mdrun 命令"""
    if has_gpu and gpu_available_on_system:
        # 有 GPU 且系统支持：不要回退
        return """gmx mdrun -nb gpu -bonded gpu -pme gpu -ntmpi 1 -ntomp 8 || {
            echo "Error: GPU available but mdrun failed"
            echo "Check forcefield, topology, and MDP parameters"
            exit 1
        }"""
    else:
        # 无 GPU 或系统不支持：允许回退
        return """gmx mdrun -nb gpu -bonded gpu -pme gpu -ntmpi 1 -ntomp 8 || \
            gmx mdrun -ntmpi 1 -ntomp 10"""
```

## 运行测试

### 快速测试（可选）

```bash
cd gaff_test_output/GMX_PROLIG_FEP/bound

# 修改 MDP：只运行 100 步
sed -i 's/nsteps.*/nsteps = 100/' mdps/prod_00.mdp

# 生成第一个窗口的 TPR
window_dir="window_00"
mkdir -p $window_dir
gmx grompp -f mdps/prod_00.mdp -c input/conf.gro \
    -p topol.top -o $window_dir/prod.tpr -maxwarn 2

# 运行（GPU 模式）
cd $window_dir
gmx mdrun -deffnm prod -nb gpu -bonded gpu -pme gpu \
    -ntmpi 1 -ntomp 8 -gpu_id 0 -nsteps 100

# 检查输出
ls -lh prod.*  # 应该看到 .gro, .xtc, .edr, .dhdl.xvg

# 清理
cd ../..
rm -rf $window_dir
```

### 标准运行流程

**Standard 模式**：
```bash
cd gaff_test_output/GMX_PROLIG_FEP/bound
bash run_prod.sh standard  # 或直接 bash run_prod_standard.sh
```

**Repex 模式**：
```bash
cd gaff_test_output/GMX_PROLIG_FEP/bound
bash run_prod.sh repex    # 或直接 bash run_prod_repex.sh
```

**运行所有**：
```bash
cd gaff_test_output/GMX_PROLIG_FEP
bash run_fep.sh all
```

## 典型运行问题

### 案例 1：GPU 利用率极低（性能灾难）

**症状**：
```
nvidia-smi
# GPU-Util: 0-4%
# 性能：3-8 ns/day
```

**诊断**：
```bash
grep "calc-lambda-neighbors" bound/mdps/prod_00.mdp
# calc-lambda-neighbors = -1  ← 错误！
```

**解决**：
```bash
sed -i 's/calc-lambda-neighbors = -1/calc-lambda-neighbors = 1/g' \
    bound/mdps/prod_*.mdp unbound/mdps/prod_*.mdp
```

**性能提升**：27-73 倍

### 案例 2：脚本未使用 GPU

**症状**：脚本中缺少 GPU 参数

**诊断**：
```bash
grep -E "nb|bonded|pme" bound/run_prod_standard.sh
# 没有输出或只有注释
```

**解决**：脚本必须包含 `-nb gpu -bonded gpu -pme gpu`

### 案例 3：执行模式选择错误

**症状**：想要副本交换，但窗口独立运行

**诊断**：
```bash
grep "execution:" fep.yaml
# execution:
#   mode: standard  ← 错误！应该是 repex
```

**解决**：修改配置或直接调用 `bash run_prod_repex.sh`

### 案例 4：窗口目录未创建

**症状**：`find bound -type d -name "window_*"` 返回 0

**原因**：这是正常的！窗口目录在运行时才创建

**解决**：运行脚本会自动创建

## 运行时监控

### 监控 GPU 使用

```bash
# 在另一个终端监控 GPU
watch -n 1 nvidia-smi

# 检查 GROMACS 进程
ps aux | grep gmx

# Standard 模式：应该看到多个 gmx 进程
# Repex 模式：应该看到一个 gmx_mpi 进程
```

### 监控日志输出

```bash
# 实时查看日志
tail -f bound/window_00/prod.log

# 检查所有窗口进度
for dir in bound/window_*; do
    if [ -f "$dir/prod.log" ]; then
        progress=$(grep "Step" $dir/prod.log | tail -1)
        echo "$dir: $progress"
    fi
done

# Repex 模式：检查交换信息
grep "Replica exchange" bound/window_00/prod.log
```

### 检查输出文件

```bash
# 检查轨迹文件
ls -lh bound/window_*/*.xtc

# 检查能量文件
ls -lh bound/window_*/*.edr

# 检查 FEP 特有输出
ls -lh bound/window_*/dhdl.xvg

# Repex 模式：检查交换状态
ls -lh bound/window_*/replica_exchange.xvg
```

## 必须避免的错误

1. **不要禁用 GPU（最严重！）**：
   - ❌ 不使用 `-nb gpu -bonded gpu -pme gpu`
   - ✅ **未经允许的情况下要尽量使用 GPU**

2. **不要使用 calc-lambda-neighbors = -1**：
   - ❌ 计算所有状态，性能灾难
   - ✅ 只计算相邻状态（性能提升 27-73 倍）

3. **不要混淆执行模式**：
   - ❌ 想用 repex 但配置了 standard
   - ✅ 明确需求，选择正确模式

4. **不要硬编码资源参数**：
   - ❌ 脚本中写死 GPU 数量
   - ✅ 使用变量：`NUM_GPUS=${NUM_GPUS:-4}`

5. **不要跳过运行前验证**：
   - 必须先 `gmx check`
   - 必须先验证 GPU 配置
   - 必须检查执行模式

6. **不要忽略 bound/unbound 一致性**：
   - 两个 leg 的配置应该一致
   - 使用相同的执行模式

7. **不要轻易认为是软件问题**：
   - ❌ "GROMACS 参数有问题"
   - ✅ **优先认为是建模/topology/MDP 参数问题**

## 性能基准检查清单

**运行 FEP 前必须验证**：

### GPU 配置
- [ ] `calc-lambda-neighbors = 1`（不是 -1）
- [ ] 脚本包含 `-nb gpu -bonded gpu -pme gpu`
- [ ] 使用 `-ntmpi 1 -ntomp 8`（GPU 优化）

### 执行模式
- [ ] 配置文件明确 `execution.mode`
- [ ] 选择了正确的执行模式（standard/repex）
- [ ] 理解两种模式的资源语义差异

### 性能预期
- [ ] GPU 利用率预期：>80%
- [ ] 性能预期：>200 ns/day（单窗口）
- [ ] 多 GPU 并行策略合理

**如果性能不达标**：
1. 检查 `calc-lambda-neighbors`
2. 检查 GPU 参数是否启用
3. 检查执行模式是否正确
4. 检查 GPU 利用率（`nvidia-smi`）
5. **不要怀疑 GROMACS 或常用参数值**

## 相关代码

- 脚本生成：`prism/fep/modeling/script_writer.py::_write_fep_master_script()`
- Standard 脚本：`prism/fep/modeling/script_writer.py::_write_standard_script()`
- Repex 脚本：`prism/fep/modeling/script_writer.py::_write_repex_script()`
- GPU 检测：`prism/fep/modeling/script_writer.py::_detect_hardware()`
- MDP 生成：`prism/fep/modeling/mdp_generator.py::generate_production_mdps()`
- 窗口创建：`prism/fep/modeling/script_writer.py` 第214-219行
