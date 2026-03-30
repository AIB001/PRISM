---
name: fep-system-debug
description: FEP 体系搭建与运行调试
type: fep
---

# FEP 体系搭建与运行调试

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

## 核心职责

检查 PRISM FEP 模块搭建的完整体系是否能正常运行，验证 GROMACS 脚本、MDP 参数、bound/unbound leg 设置、lambda 调度等运行时配置。

## 标准化命名规范

**CRITICAL**: 所有代码必须遵循 PRISM 标准化命名规范（详见项目根目录 `NAMING_CONVENTIONS.md`）。

### 核心规则

1. **系统目录命名**：
   - ✅ `GMX_PROLIG_FEP/` （标准 FEP 系统目录）
   - ❌ `FEP_SYSTEM/`, `fep_output/`, `test_output/`

2. **力场输出目录**：
   - ✅ 使用 `ffgen.get_output_dir_name()` 获取
   - ✅ GAFF: `LIG.amb2gmx/`
   - ✅ OpenFF: `LIG.openff2gmx/`
   - ❌ 硬编码路径

3. **测试目录命名**：
   - ✅ `gaff_test/`, `openff_test/`, `test_42_38/`
   - ❌ `gaff2_e2e_test/`, `ligand_42/`, `temp_output/`

### 代码示例

```python
# ✅ 正确用法
from prism.forcefield.gaff import GAFFForceFieldGenerator

ffgen = GAFFForceFieldGenerator(
    ligand_path="ligand.mol2",
    output_dir="gaff_test",  # 清晰标识
)
ffgen.run()

# 获取标准化目录名
prism_dir = Path(ffgen.output_dir) / ffgen.get_output_dir_name()
# 结果: gaff_test/LIG.amb2gmx/

# ✅ FEP 系统使用标准命名
builder = FEPScaffoldBuilder(
    output_dir="GMX_PROLIG_FEP",  # 标准命名
)

# ❌ 错误用法
prism_dir = "gaff_test/LIG.amb2gmx"  # 硬编码
builder = FEPScaffoldBuilder(output_dir="FEP_SYSTEM")  # 不规范
```

### 测试文件检查清单

创建新测试时，确保：
- [ ] 使用 `{forcefield}_test/` 格式命名测试目录
- [ ] 使用 `ffgen.get_output_dir_name()` 获取配体目录
- [ ] FEP 系统使用 `GMX_PROLIG_FEP/` 标准命名
- [ ] 不使用描述性不足的名称（如 `temp/`, `output/`）

## 目录安排和复用机制

### Lambda 窗口生成机制（重要）

**两阶段工作流**：

1. **Build 阶段**（搭建阶段，`FEPScaffoldBuilder.build_from_components()`）：
   - ✅ 生成 MDP 文件到 `mdps/` 目录
   - ✅ 生成 `run_fep.sh` 主脚本
   - ✅ 生成 `lambda_schedule.json`（包含 lambda 值定义）
   - ❌ **不创建** `window_XX/` 目录
   - 示例：`lambda_windows=5` → 生成 `prod_00.mdp` 到 `prod_04.mdp`（5个文件）

2. **Run 阶段**（执行阶段，`bash run_fep.sh`）：
   - ✅ **动态创建** `window_XX/` 目录
   - ✅ 每个 lambda 窗口在执行时创建
   - ✅ 代码位置：`prism/fep/modeling/script_writer.py::write_fep_master_script()` 第214-219行
   - 示例：`for lambda_mdp in mdps/prod_*.mdp; do mkdir -p "window_${lambda_idx}"; done`

**常见误解**：
- ❌ "Lambda 窗口数量为 0 不是正常的" → **这是正常的！**
- ✅ Build 完成后，`window_XX/` 目录数量应该为 0
- ✅ 运行 `bash run_fep.sh` 后，才会创建 `window_XX/` 目录

**验证命令**：
```bash
# Build 完成后（应该返回 0）
find bound -type d -name "window_*" | wc -l

# 运行后（应该返回 lambda_windows 数量）
bash bound/run_fep.sh bound
find bound -type d -name "window_*" | wc -l
```

### Mapping 输出 → Building 输入的数据流

**数据流向**：
```
原始配体文件 (mol2/SDF)
    ↓
PRISM ForceFieldGenerator (参数化)
    ↓
PRISM 格式 (ITP+GRO)
    ↓
DistanceAtomMapper (原子映射)
    ↓
混合拓扑文件 (hybrid.itp + hybrid.gro)  ← Mapping 输出
    ↓
FEPScaffoldBuilder.build_from_components()  ← 自动复用
    ↓
完整 FEP 系统 (bound/unbound legs + MDPs + 脚本)
```

### Mapping 输出目录

**位置**: `common/hybrid/`

**关键文件**:
- `hybrid.itp` - 混合拓扑文件（核心，包含原子映射信息）
- `hybrid.gro` - 混合坐标文件
- `atomtypes_hybrid.itp` - 原子类型定义
- `ff_hybrid.itp` - 力场包含文件
- `posre_hybrid.itp` - 位置约束文件
- `ligand_seed.pdb` - 可视化种子文件

**生成代码**: `prism/fep/gromacs/itp_builder.py::ITPBuilder`

**示例输出结构**:
```
gaff_test_output/GMX_PROLIG_FEP/
└── common/
    └── hybrid/
        ├── hybrid.itp              ← 核心文件
        ├── hybrid.gro
        ├── atomtypes_hybrid.itp
        ├── ff_hybrid.itp
        ├── posre_hybrid.itp
        └── ligand_seed.pdb
```

### System Building 输入要求

**核心函数**: `prism/fep/modeling/core.py::FEPScaffoldBuilder.build_from_components()`

**输入参数**:
```python
def build_from_components(
    receptor_pdb: str,              # 受体 PDB 文件
    hybrid_itp: str,                # ← Mapping 输出的 hybrid.itp
    reference_ligand_dir: str,      # 参考配体力场目录 (包含 LIG.itp, LIG.gro)
    mutant_ligand_dir: Optional[str],  # 突变配体力场目录 (可选)
    hybrid_gro: Optional[str],      # 混合 GRO 文件 (可选，默认从 reference_ligand_dir 读取)
    bound_system_dir: Optional[str],    # PRISM bound 系统 (可选)
    unbound_system_dir: Optional[str],  # PRISM unbound 系统 (可选)
) -> FEPScaffoldLayout
```

**复用机制**:
- ✅ **自动复用**: `hybrid.itp` 直接作为 `hybrid_itp` 参数传入
- ✅ **自动复制**: `HybridPackageBuilder` 会将 mapping 输出复制到 `common/hybrid/`
- ✅ **无缝集成**: 无需手动移动或重命名文件

**自动化程度**: 85-90%

### 端到端工作流示例

**完整流程**（从配体到 FEP 系统）:
```python
# Step 1: 配体参数化
from prism.forcefield.gaff import GAFFForceFieldGenerator

generator_a = GAFFForceFieldGenerator(
    ligand_path="ligand_a.mol2",
    output_dir="ligand_a_ff",
)
prism_dir_a = generator_a.run()  # → ligand_a_ff/LIG.itp + LIG.gro

generator_b = GAFFForceFieldGenerator(
    ligand_path="ligand_b.mol2",
    output_dir="ligand_b_ff",
)
prism_dir_b = generator_b.run()  # → ligand_b_ff/LIG.itp + LIG.gro

# Step 2: 原子映射
from prism.fep.io import read_ligand_from_prism
from prism.fep.core.mapping import DistanceAtomMapper

lig_a = read_ligand_from_prism(
    itp_file="ligand_a_ff/LIG.itp",
    gro_file="ligand_a_ff/LIG.gro",
)

lig_b = read_ligand_from_prism(
    itp_file="ligand_b_ff/LIG.itp",
    gro_file="ligand_b_ff/LIG.gro",
)

mapper = DistanceAtomMapper(dist_cutoff=0.6, charge_cutoff=0.05)
mapping = mapper.map(lig_a, lig_b)  # → hybrid.itp

# Step 3: FEP 系统搭建（自动复用 mapping 输出）
from prism.fep.modeling.core import FEPScaffoldBuilder

builder = FEPScaffoldBuilder(
    output_dir="fep_output",
    lambda_windows=32,
)

layout = builder.build_from_components(
    receptor_pdb="receptor.pdb",
    hybrid_itp="hybrid.itp",              # ← Mapping 输出
    reference_ligand_dir="ligand_a_ff",   # ← PRISM 格式
    mutant_ligand_dir="ligand_b_ff",      # ← PRISM 格式
)

# 输出结构:
# fep_output/
# ├── common/hybrid/hybrid.itp  ← 从 mapping 输出复制
# ├── bound/
# └── unbound/
```

### 验证复用是否成功

**检查清单**:
```bash
# 1. 检查 mapping 输出是否存在
ls -lh common/hybrid/hybrid.itp

# 2. 检查 building 输入是否正确
python -c "
from pathlib import Path
hybrid_itp = Path('common/hybrid/hybrid.itp')
print(f'Hybrid ITP exists: {hybrid_itp.exists()}')
print(f'Size: {hybrid_itp.stat().st_size} bytes')
"

# 3. 检查 bound/unbound leg 是否正确生成
ls -d fep_output/bound fep_output/unbound

# 4. 检查关键文件
ls fep_output/common/hybrid/hybrid.itp
ls fep_output/bound/topol.top
ls fep_output/unbound/topol.top
```

## 关键检查点

### 0. 配置验证（第一步！）

**在调试任何问题前，先验证配置是否正确读取**

```bash
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh

# 1. 检查所有 YAML 配置文件
ls -la fep*.yaml config*.yaml

# 2. 验证 FEPConfig 读取
python -c "
from prism.fep.config import FEPConfig
cfg = FEPConfig('.', config_file='config_gaff.yaml', fep_file='fep_ref.yaml')
print('✓ FEPConfig loaded')
print('  Mapping params:', cfg.get_mapping_params())
print('  Simulation params:', cfg.get_simulation_params())
print('  Lambda params:', cfg.get_lambda_params())
"

# 3. 对比不同模式的配置差异
echo "=== Mapping 配置对比 ==="
for mode in ref mut mean none; do
    echo "--- $mode ---"
    cat fep_${mode}.yaml | grep -A 6 "mapping:" | head -7
done
```

**代码路径**：
- 配置管理：`prism/fep/config.py::FEPConfig` (整个类)
- Mapping 参数：`prism/fep/config.py::FEPConfig.get_mapping_params()` (line 73-84)
- Simulation 参数：`prism/fep/config.py::FEPConfig.get_simulation_params()` (line 93-113)
- Lambda 参数：`prism/fep/config.py::FEPConfig.get_lambda_params()` (line 126-150)

### 1. 参数传递优先级

**优先级顺序**（从高到低）：
1. CLI 参数（命令行直接传入）
2. `fep.yaml`（FEP 专用配置）
3. `config.yaml`（基础配置）
4. 代码默认值

**代码路径**：
- 配置读取：`prism/fep/io.py::FEPConfig.__init__()`
- 参数合并：`prism/fep/io.py::FEPConfig._merge_configs()`

**检查命令**：
```bash
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh

# 查看最终使用的配置
python -c "
from prism.fep.io import FEPConfig
cfg = FEPConfig('.', config_file='config_gaff.yaml', fep_file='fep_ref.yaml')
print('Mapping params:', cfg.get_mapping_params())
print('Build params:', cfg.get_build_params())
print('Run params:', cfg.get_run_params())
"
```

**常见问题**：
- ❌ "参数在多个层次重复出现" → 应该去重或模块化
- ❌ "CLI 参数未生效" → 检查优先级是否正确
- ❌ "fep.yaml 未被读取" → 检查文件路径和格式

### 2. Lambda 调度模式

**Decoupled vs Coupled**：
- **Decoupled**：bound 和 unbound leg 独立运行
  - Bound leg：蛋白-配体复合物
  - Unbound leg：配体在溶液中
  - 最终：ΔΔG = ΔG_bound - ΔG_unbound
- **Coupled**：同时计算（较少用）

**检查方法**：
```bash
# 查看 MDP 文件中的 lambda 设置
grep -h "free-energy" gaff_test_output/GMX_PROLIG_FEP/bound/*/mdp/*.mdp
grep -h "free-energy" gaff_test_output/GMX_PROLIG_FEP/unbound/*/mdp/*.mdp

# 查看 lambda 分布
grep "init-lambda-state" gaff_test_output/GMX_PROLIG_FEP/*/run.sh
```

**代码路径**：
- Lambda 生成：`prism/fep/modeling/hybrid_package.py::_generate_lambda_schedule()`
- MDP 生成：`prism/fep/modeling/hybrid_package.py::_write_mdp_files()`

### 3. 脚本生成检查

**应该生成的脚本**：
- `bound/run.sh` - bound leg 运行脚本
- `unbound/run.sh` - unbound leg 运行脚本
- 不应该有太多变体（local/parallel 应通过参数控制）

**检查清单**：
```bash
cd gaff_test_output/GMX_PROLIG_FEP

# 1. 检查脚本数量（不应该太多）
find . -name "run.sh" | wc -l
# 应该只有 2-4 个：bound+unbound，可能还有 replica 脚本

# 2. 检查脚本是否使用 GPU
grep -h "gpu" */run.sh | grep -v "^#"
# 应该看到：-nb gpu -bonded gpu -pme gpu

# 3. 检查并行设置
grep -h "ntomp\|ntmpi" */run.sh
# 应该合理利用 CPU 核心

# 4. 检查 bound/unbound 脚本一致性
diff bound/run.sh unbound/run.sh | head -20
# 主要差异应该是：工作目录、输入文件
```

**代码路径**：
- 脚本生成：`prism/fep/modeling/hybrid_package.py::_write_run_scripts()`

**常见问题**：
- ❌ "生成太多版本脚本" → 应该通过参数控制，不是生成不同脚本
- ❌ "bound/unbound 脚本各自一套" → 应该通用，只改必要部分
- ❌ "GPU 未启用" → 检查 GPU 检测逻辑

### 4. MDP 参数验证

**关键 MDP 参数**：
```bash
# 查看 free-energy 相关设置
grep -A 5 "free-energy" gaff_test_output/GMX_PROLIG_FEP/bound/lambda_0/mdp/md.mdp

# 应该看到：
; Free energy parameters
free-energy              = yes
init-lambda-state        = 0
delta-lambda             = 0
calc-lambda-neighbors     = 0  ; 或 1 (if using soft-core)
coupled-lambdas          = no  ; decoupled 模式
```

**检查清单**：
- [ ] `free-energy = yes`
- [ ] `init-lambda-state` 正确（0-1 分布）
- [ ] `calc-lambda-neighbors` 合理
- [ ] `coupled-lambdas = no` (decoupled 模式)
- [ ] 温度设置一致（通常 310 K）
- [ ] 时间步长合理（通常 0.002 ps = 2 fs）

**代码路径**：
- MDP 模板：`prism/fep/modeling/templates/fep.mdp.template`
- MDP 生成：`prism/fep/modeling/hybrid_package.py::_write_mdp_files()`

### 5. 目录结构验证

**标准目录结构**：
```
GMX_PROLIG_FEP/
├── bound/
│   ├── lambda_0/
│   │   ├── mdp/
│   │   │   └── md.mdp
│   │   ├── run.sh
│   │   └── # 其他文件
│   ├── lambda_1/
│   ├── ...
│   └── run.sh  (可选：批量运行脚本)
├── unbound/
│   ├── lambda_0/
│   ├── ...
│   └── run.sh
└── common/
    ├── hybrid/
    │   ├── hybrid.itp
    │   └── hybrid.gro
    └── topol.top
```

**检查命令**：
```bash
cd gaff_test_output/GMX_PROLIG_FEP

# 1. 检查目录存在
ls -la bound/ unbound/ common/

# 2. 检查 lambda 窗口数量
ls -d bound/lambda_* | wc -l
# 通常 11-21 个窗口

# 3. 检查关键文件
ls bound/lambda_0/mdp/md.mdp
ls common/hybrid/hybrid.itp
ls common/hybrid/hybrid.gro
ls common/topol.top

# 4. 检查文件完整性
for dir in bound/lambda_* unbound/lambda_*; do
    if [ ! -f "$dir/mdp/md.mdp" ]; then
        echo "Missing: $dir/mdp/md.mdp"
    fi
done
```

### 6. 运行前验证

**快速验证脚本**：
```bash
cd gaff_test_output/GMX_PROLIG_FEP

# 1. 检查 GROMACS 可用
gmx --version

# 2. 检查 topology 语法
gmx check -f common/hybrid/hybrid.gro -c common/hybrid/hybrid.gro -p common/topol.top

# 3. 尝试生成第一个 lambda 的 TPR
cd bound/lambda_0
gmx grompp -f mdp/md.mdp -c ../../common/hybrid/hybrid.gro -p ../../common/topol.top -n ../../common/index.ndx -o tpr.tpr

# 4. 如果成功，检查 TPR
gmx check -s tpr.tpr

# 5. 清理测试文件
rm -f tpr.tpr *.cpt
```

**常见错误**：
- `Atom type XXX not found` → 力场参数缺失
- `No default bond type` → 键级参数缺失
- `number of coordinates in file` → 坐标数不匹配

## 典型问题案例

### 案例 1：生成的脚本不反映配置

**症状**：
- 配置文件设置 `gpu_id = 1`
- 但脚本里使用 `gpu_id = 0`

**原因**：
- 参数传递断裂
- 脚本模板未正确使用配置

**排查**：
```python
# 检查配置读取
cfg = FEPConfig(...)
print("GPU ID:", cfg.get_run_params()['gpu_id'])

# 检查脚本生成
# prism/fep/modeling/hybrid_package.py::_write_run_scripts()
# 确保使用了 cfg.get_run_params()
```

**代码位置**：`prism/fep/modeling/hybrid_package.py::_write_run_scripts()`

### 案例 2：Bound/Unbound 脚本混乱

**症状**：
```
bound/
  ├── local_run.sh
  ├── parallel_run.sh
  ├── single_replica_run.sh
  └── ...
unbound/
  ├── local_run.sh
  ├── parallel_run.sh
  └── ...
```

**问题**：
- 太多脚本版本
- 应该通过参数控制，不是生成不同脚本

**修复**：
- 只生成一个通用 `run.sh`
- 通过 CLI 参数控制：`./run.sh --parallel --replica 1-3`

**代码位置**：`prism/fep/modeling/hybrid_package.py::_write_run_scripts()`

### 案例 3：MDP 参数未反映 lambda 调度

**症状**：
- 配置 11 个 lambda 窗口
- 但 MDP 中 `init-lambda-state` 都是 0

**原因**：
- MDP 模板未正确替换
- Lambda 循环未正确遍历

**检查**：
```bash
for i in bound/lambda_*; do
    lambda_id=$(basename $i | sed 's/lambda_//')
    echo "Lambda $lambda_id:"
    grep "init-lambda-state" $i/mdp/md.mdp
done
```

**代码位置**：`prism/fep/modeling/hybrid_package.py::_write_mdp_files()`

### 案例 4：GPU 未充分利用

**症状**：
- 脚本中只有 `-ntmpi 1 -ntomp 4`
- 机器有 10 核心 + GPU

**原因**：
- GPU 检测失败
- 默认值设置保守

**修复**：
```python
# 检测 GPU
def detect_gpu():
    try:
        result = subprocess.run(['nvidia-smi', '--query-gpu=name', '--format=csv,noheader'],
                              capture_output=True, text=True)
        if result.returncode == 0:
            return True
    except:
        pass
    return False

# 根据硬件设置并行参数
if has_gpu:
    ntmpi = 1  # GPU 模式通常用 1 MPI
    ntomp = 10 # 使用更多 OpenMP 线程
else:
    ntmpi = 4
    ntomp = 2
```

**代码位置**：`prism/fep/modeling/hybrid_package.py::_detect_hardware()`

### 案例 5: Mapping 输出无法被 Building 复用

**症状**：
- Mapping 生成了 `hybrid_mean_complete.itp`
- Building 报错找不到 `hybrid.itp`

**原因**：
- 文件命名约定不一致
- Mapping 输出包含模式名（mean/ref/mut/none）
- Building 期望固定文件名 `hybrid.itp`

**解决方案**：
```python
# 方案 1: 重命名文件（推荐）
from pathlib import Path
Path("hybrid_mean_complete.itp").rename("hybrid.itp")

# 方案 2: 在 build_from_components() 中指定正确路径
builder.build_from_components(
    ...
    hybrid_itp="hybrid_mean_complete.itp",  # 使用实际文件名
)

# 方案 3: 使用符号链接（快速测试）
ln -s hybrid_mean_complete.itp hybrid.itp
```

**检查方法**：
```bash
# 查看实际生成的文件名
ls -lh common/hybrid/*.itp

# 检查文件内容（第一行应该有分子名）
head -1 common/hybrid/hybrid*.itp

# 验证 mapping 是否正确
grep -c "common" common/hybrid/hybrid*.itp
```

**代码位置**：`prism/fep/modeling/hybrid_package.py::HybridPackageBuilder.__init__`

## 快速验证脚本

```bash
cd /data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/oMeEtPh-EtPh

echo "=== FEP System Debug ==="

# 1. 检查配置读取
echo "1. Configuration:"
python -c "
from prism.fep.io import FEPConfig
cfg = FEPConfig('.', config_file='config_gaff.yaml', fep_file='fep_ref.yaml')
print('  Lambda windows:', len(cfg.get_lambda_schedule()))
print('  GPU ID:', cfg.get_run_params().get('gpu_id', 'default'))
"

# 2. 检查目录结构
echo "2. Directory structure:"
if [ -d "gaff_test_output/GMX_PROLIG_FEP" ]; then
    cd gaff_test_output/GMX_PROLIG_FEP
    echo "  Bound lambdas: $(ls -d bound/lambda_* 2>/dev/null | wc -l)"
    echo "  Unbound lambdas: $(ls -d unbound/lambda_* 2>/dev/null | wc -l)"
    echo "  Run scripts: $(find . -name 'run.sh' | wc -l)"
else
    echo "  ERROR: GMX_PROLIG_FEP not found"
fi

# 3. 检查 MDP 参数
echo "3. MDP parameters:"
if [ -d "gaff_test_output/GMX_PROLIG_FEP/bound/lambda_0" ]; then
    grep "free-energy\|init-lambda" gaff_test_output/GMX_PROLIG_FEP/bound/lambda_0/mdp/md.mdp | head -5
fi

# 4. 检查 GPU 设置
echo "4. GPU settings:"
if [ -d "gaff_test_output/GMX_PROLIG_FEP" ]; then
    grep -h "gpu" gaff_test_output/GMX_PROLIG_FEP/*/run.sh 2>/dev/null | head -3
fi

echo "=== End ==="
```

## 运行测试（可选）

如果环境允许，可以运行短测试：
```bash
cd gaff_test_output/GMX_PROLIG_FEP/bound/lambda_0

# 修改 MDP：只运行 100 步
sed -i 's/nsteps.*/nsteps = 100/' mdp/md.mdp

# 生成 TPR
gmx grompp -f mdp/md.mdp -c ../../common/hybrid/hybrid.gro \
    -p ../../common/topol.top -n ../../common/index.ndx -o tpr.tpr

# 运行（如果 GPU 可用）
gmx mdrun -deffnm test -s tpr.tpr -ntmpi 1 -ntomp 4 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -nsteps 100

# 检查输出
ls -lh test.*  # 应该看到 .gro, .xtc, .edr

# 清理
rm -f test.* tpr.tpr #*.cpt
```

## 必须避免的错误

1. **不要生成太多脚本版本**：
   - ❌ `local_run.sh`, `parallel_run.sh`, `single_run.sh`...
   - ✅ 一个 `run.sh`，通过参数控制

2. **不要硬编码 GPU ID**：
   - ❌ `-gpu_id 0`
   - ✅ `-gpu_id ${GPU_ID:-0}` (可配置)

3. **不要忽略 bound/unbound 一致性**：
   - 两个 leg 的脚本应该几乎一样
   - 只改必要部分（输入文件、工作目录）

4. **不要跳过运行前验证**：
   - 必须先 `gmx check`
   - 必须先 `gmx grompp` 成功

## 相关代码

- 配置管理：`prism/fep/io.py::FEPConfig`
- 脚本生成：`prism/fep/modeling/hybrid_package.py::_write_run_scripts()`
- MDP 生成：`prism/fep/modeling/hybrid_package.py::_write_mdp_files()`
- Lambda 调度：`prism/fep/modeling/hybrid_package.py::_generate_lambda_schedule()`
- 目录构建：`prism/fep/modeling/hybrid_package.py::build_fep_system()`
