---
name: fep-system-setup
description: FEP system setup and scaffold generation debugging
type: fep
---

# FEP 体系搭建与脚手架生成调试

## 📝 Code Search and Navigation

**LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation:
- `mcp__cclsp__find_definition` - Find symbol definitions
- `mcp__cclsp__find_references` - Find symbol references
- `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
- Benefits: More accurate than grep, skips comments/strings, reduces token usage
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

## 核心职责

检查 PRISM FEP 模块搭建的完整体系，验证输入文件、力场参数、原子映射、混合拓扑生成、bound/unbound leg 脚手架结构。

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

## Bound 态坐标与复合物完整性检查（必须）

这部分是 **FEP 体系搭建检查的硬性要求**。只检查文件存在不够，必须检查配体是否真的被正确放回蛋白口袋中。

### 必查项

1. **`bound/repeat1/input/conf.gro` 中必须有 hybrid ligand**
   - 配体残基名应为 `HYB`，不是 `LIG`
   - `HYB` 原子数应与 `common/hybrid/hybrid.gro` 一致

2. **`common/hybrid/hybrid.gro` 不能是退化坐标**
   - 不能出现 “所有原子坐标都一样”
   - 不能全部是 `(0, 0, 0)`
   - 至少应有多个不同坐标点

3. **蛋白和配体必须在一起**
   - `bound/repeat1/input/conf.gro` 中 `HYB` 与蛋白原子之间应有合理空间关系
   - 不能出现配体整体漂到盒子另一边
   - 常用快速判断：
     - 配体-蛋白最近原子距离不是数纳米量级
     - 配体质心到蛋白 pocket 附近的距离应明显小于盒子边长

4. **EM 后仍要在一起**
   - 不能只看 `input/conf.gro`
   - 还要至少检查：
     - `bound/repeat1/build/em.gro`
   - 如果 `conf.gro` 正常但 `em.gro` 里 ligand 飞走，仍然算 setup 有问题

### 推荐检查命令

```bash
# 1. 检查 HYB 是否存在
python - <<'PY'
from pathlib import Path
p = Path("bound/repeat1/input/conf.gro")
lines = p.read_text().splitlines()
counts = {}
for line in lines[2:-1]:
    if len(line) >= 15:
        resname = line[5:10].strip()
        counts[resname] = counts.get(resname, 0) + 1
print("HYB atoms =", counts.get("HYB", 0))
PY

# 2. 检查 hybrid.gro 是否退化
python - <<'PY'
from pathlib import Path
p = Path("common/hybrid/hybrid.gro")
coords = []
for line in p.read_text().splitlines()[2:-1]:
    if len(line) >= 44:
        coords.append((float(line[20:28]), float(line[28:36]), float(line[36:44])))
print("atom count =", len(coords))
print("unique coords =", len(set(coords)))
PY
```

### 判定标准

- `HYB atoms > 0`
- `unique coords > 1`
- `bound/repeat1/input/conf.gro` 与 `bound/repeat1/build/em.gro` 中，配体都仍与蛋白接近
- 若发现：
  - `HYB` 不存在
  - `HYB` 全部坐标相同
  - `HYB` 全部是零坐标
  - 配体与蛋白明显错位
  
  则应直接判定 **体系搭建失败**，不能进入 runtime 测试

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
│   ├── input/
│   │   └── conf.gro          # 配置坐标文件
│   ├── mdps/                 # MDP 参数文件目录
│   │   ├── em.mdp
│   │   ├── nvt.mdp
│   │   ├── npt.mdp
│   │   ├── npt_short_00.mdp  # Lambda 专用 NPT (生产)
│   │   ├── npt_short_01.mdp
│   │   ├── ...
│   │   ├── prod_00.mdp       # Lambda 生产运行
│   │   ├── prod_01.mdp
│   │   └── ...
│   ├── topol.top             # 拓扑文件
│   └── run_fep.sh            # 运行脚本
├── unbound/                  # 同 bound 结构
│   └── ...
└── common/
    ├── hybrid/               # 混合拓扑
    │   ├── hybrid.itp
    │   └── hybrid.gro
    └── protein/              # 蛋白参考
        └── protein.pdb
```

**检查命令**：
```bash
cd gaff_test_output/GMX_PROLIG_FEP

# 1. 检查目录存在
ls -la bound/ unbound/ common/

# 2. 检查 MDP 文件数量（应该对应 lambda_windows 数量）
ls bound/mdps/prod_*.mdp | wc -l
# 例如 lambda_windows=32 → 应该有 32 个 prod_XX.mdp

# 3. 检查关键文件
ls bound/input/conf.gro
ls bound/topol.top
ls common/hybrid/hybrid.itp
ls common/hybrid/hybrid.gro

# 4. 检查文件内容
head -5 bound/mdps/prod_00.mdp  # 应该有 free-energy 参数
head -5 common/hybrid/hybrid.itp  # 应该有 [ moleculetype ]
```

### 6. 搭建后验证

**快速验证脚本**：
```bash
cd gaff_test_output/GMX_PROLIG_FEP

# 1. 检查 GROMACS 可用
gmx --version

# 2. 检查 topology 语法
gmx check -f common/hybrid/hybrid.gro -p bound/topol.top

# 3. 尝试生成第一个 lambda 的 TPR
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
# prism/fep/modeling/script_writer.py::_write_fep_master_script()
# 确保使用了 cfg.get_run_params()
```

**代码位置**：`prism/fep/modeling/script_writer.py::_write_fep_master_script()`

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
- 只生成一个通用 `run_fep.sh`
- 通过 CLI 参数控制：`bash run_fep.sh --mode standard --leg bound`

**代码位置**：`prism/fep/modeling/script_writer.py::_write_fep_master_script()`

### 案例 3：MDP 参数未反映 lambda 调度

**症状**：
- 配置 11 个 lambda 窗口
- 但 MDP 中 `init-lambda-state` 都是 0

**原因**：
- MDP 模板未正确替换
- Lambda 循环未正确遍历

**检查**：
```bash
for i in bound/mdps/prod_*.mdp; do
    lambda_id=$(basename $i | sed 's/prod_//' | sed 's/.mdp//')
    echo "Lambda $lambda_id:"
    grep "init-lambda-state" $i
done
```

**代码位置**：`prism/fep/modeling/mdp_generator.py::generate_production_mdps()`

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

**代码位置**：`prism/fep/modeling/script_writer.py::_detect_hardware()`

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

**代码位置**：`prism/fep/modeling/core.py::FEPScaffoldBuilder.build_from_components()`

## 快速验证脚本

```bash
cd /data2/gxf1212/work/PRISM/tests/gxf/FEP/unit_test/oMeEtPh-EtPh

echo "=== FEP System Setup Debug ==="

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
    echo "  Bound MDPs: $(ls bound/mdps/prod_*.mdp 2>/dev/null | wc -l)"
    echo "  Unbound MDPs: $(ls unbound/mdps/prod_*.mdp 2>/dev/null | wc -l)"
    echo "  Hybrid files: $(ls common/hybrid/*.* 2>/dev/null | wc -l)"
else
    echo "  ERROR: GMX_PROLIG_FEP not found"
fi

# 3. 检查 MDP 参数
echo "3. MDP parameters:"
if [ -f "gaff_test_output/GMX_PROLIG_FEP/bound/mdps/prod_00.mdp" ]; then
    grep "free-energy\|init-lambda" gaff_test_output/GMX_PROLIG_FEP/bound/mdps/prod_00.mdp | head -5
fi

echo "=== End ==="
```

## ITP 文件格式要求（重要）

### State B 参数写入规则

**CRITICAL**: GROMACS FEP 使用 typeB/chargeB/massB 列来编码 B 态参数。正确写入这些参数是 FEP 计算成功的关键。

#### 写入条件（`ITPBuilder._get_atoms_section()`）

**B 态参数写入逻辑**：
```python
# prism/fep/gromacs/itp_builder.py:132-136
if atom.state_b_type is not None and atom.state_b_charge is not None:
    if atom.state_b_type != atom.state_a_type or abs(atom.state_b_charge - atom.state_a_charge) > 1e-6:
        line += f" {atom.state_b_type:>10s} {atom.state_b_charge:10.6f}"
        mass_b = atom.mass_b if atom.mass_b is not None else atom.mass
        line += f" {mass_b:10.5f}"
```

**关键规则**：
- ✅ **写入 B 态**: `state_b_type != state_a_type` **或** `charge != 0.000001`
- ❌ **不写入 B 态**: `state_b_type == state_a_type` **且** `charge 相同`

#### 各类原子的 B 态参数

**1. Common atoms（共有原子）**：
```
; A态和B态参数相同
    1    CG2RC0      1     LIG     C1      1   0.345000  12.01100     CG2RC0   0.345000  12.01100
                              ↑ A态参数                  ↑ B态参数（相同）
```
- charge_common='mean' 时：A=B=平均电荷
- **有 B 态列**（但参数相同）

**2. Transformed A atoms（变换原子，仅A态存在）**：
```
; A态存在，B态为dummy
    38       FGA2      1     LIG    F3A      1  -0.230000   18.99800       DUM_FGA2   0.000000   18.99800
                              ↑ A态：FGA2              ↑ B态：DUM_（dummy）
```
- **必须有 B 态**：`typeB=DUM_FGA2`, `chargeB=0.0`
- 名称：`F3A`（加A后缀）

**3. Transformed B atoms（变换原子，仅B态存在）**：
```
; A态为dummy，B态存在
    39    DUM_HGA6      1     LIG    H1B      1   0.000000    1.00800       HGA6   0.109000   1.00800
                              ↑ A态：DUM_（dummy）        ↑ B态：HGA6
```
- **必须有 B 态**：`typeA=DUM_HGA6`, `chargeA=0.0`, `typeB=HGA6`
- 名称：`H1B`（加B后缀）

**4. Surrounding atoms（周围原子，电荷不同）**：
```
; A态和B态电荷不同
    20       HGA3      1     LIG    H10      1   0.090000    1.00800       HGA2   0.090000   1.00800
                              ↑ A态：HGA3, 0.09          ↑ B态：HGA2, 0.09（类型也可能不同）
```
- **必须有 B 态**：类型或电荷不同
- 用于电荷重分配

#### 验证 ITP 文件正确性

**检查清单**：

```bash
# 1. 检查原子总数（应该 = common + transformed + surrounding）
echo "=== 原子总数 ==="
grep -c "^\s*[0-9]" GMX_PROLIG_FEP/common/hybrid/hybrid.itp

# 2. 检查是否有 B 态参数
echo "=== B态参数统计 ==="
grep "^\s*[0-9]" GMX_PROLIG_FEP/common/hybrid/hybrid.itp | awk '{
    # 检查是否有第10列（typeB）
    has_typeb = (NF >= 10 && $10 !~ /^[0-9.]+$/)
    # 检查是否有第11列（chargeB）
    has_chargeb = (NF >= 11)
    # 检查是否有第12列（massB）
    has_massb = (NF >= 12)
    
    if (has_typeb || has_chargeb || has_massb) {
        print "Atom "$1": has B-state"
    }
}' | wc -l

# 3. 检查变换原子（应该有A/B后缀）
echo "=== 变换原子 ==="
grep "^\s*[0-9]" GMX_PROLIG_FEP/common/hybrid/hybrid.itp | awk '{print $5}' | grep -E "[AB]$" 

# 4. 检查dummy原子
echo "=== Dummy原子 ==="
grep "^\s*[0-9]" GMX_PROLIG_FEP/common/hybrid/hybrid.itp | awk '{
    if ($2 ~ /^DUM_/ || $10 ~ /^DUM_/) {
        print "Atom "$1": "$2" -> "$10
    }
}'

# 5. 检查F3→H1变换（42-38系统特例）
echo "=== F3/H1原子检查 ==="
grep "^\s*[0-9]" GMX_PROLIG_FEP/common/hybrid/hybrid.itp | grep -E "F3|H1" | awk '{
    printf "Atom %s (%s): A_type=%s, B_type=%s\n", $1, $5, $2, $10
}'
```

**预期输出**（42-38系统）：
```
=== 原子总数 ===
39  # 应该是39，不是38！

=== B态参数统计 ===
39  # 所有原子都应该检查

=== 变换原子 ===
F3A  # 变换A原子
H1B  # 变换B原子

=== Dummy原子 ===
Atom 38: FGA2 -> DUM_FGA2  # F3的B态是dummy
Atom 39: DUM_HGA6 -> HGA6  # H1的A态是dummy

=== F3/H1原子检查 ===
Atom 38 (F3A): A_type=FGA2, B_type=DUM_FGA2  # F→H变换
Atom 39 (H1B): A_type=DUM_HGA6, B_type=HGA6
```

#### 调试 State B 参数缺失

**问题症状**：
- 所有 dH/dλ = 0.0000000
- ΔΔG = 0.000 kcal/mol
- ITP 文件只有 38 个原子（应该是 39）

**诊断步骤**：

```python
# 1. 检查映射结果
from prism.fep.io import read_ligand_from_prism
from prism.fep.core.mapping import DistanceAtomMapper

ref_atoms = read_ligand_from_prism(
    itp_file="input/42/gromacs/LIG.itp",
    gro_file="input/42.pdb"
)

mut_atoms = read_ligand_from_prism(
    itp_file="input/38/gromacs/LIG.itp",
    gro_file="input/38.pdb"
)

mapper = DistanceAtomMapper(dist_cutoff=0.6, charge_cutoff=0.05)
mapping = mapper.map(ref_atoms, mut_atoms)

print(f"Transformed A: {len(mapping.transformed_a)}")  # 应该 >= 1
print(f"Transformed B: {len(mapping.transformed_b)}")  # 应该 >= 1

# 如果都是0，说明映射失败！

# 2. 检查混合拓扑构建
from prism.fep.core.hybrid_topology import HybridTopologyBuilder

builder = HybridTopologyBuilder(charge_strategy='mean')
hybrid_atoms = builder.build(mapping, params_a, params_b, ref_atoms, mut_atoms)

print(f"Hybrid atoms: {len(hybrid_atoms)}")  # 应该 = common + trans_a + trans_b

# 检查变换原子
transformed = [a for a in hybrid_atoms if a.name.endswith('A') or a.name.endswith('B')]
print(f"Transformed atoms: {len(transformed)}")
for atom in transformed:
    print(f"  {atom.name}: A_type={atom.state_a_type}, B_type={atom.state_b_type}")

# 3. 检查ITP写入
from prism.fep.gromacs.itp_builder import ITPBuilder

ITPBuilder(hybrid_atoms, {}).write_itp("test.itp", molecule_name="HYB")

# 检查写入的原子数
with open("test.itp") as f:
    content = f.read()
    atom_count = len([l for l in content.split('\n') if l.strip() and not l.startswith(';') and not l.startswith('[')])
    print(f"Atoms in ITP: {atom_count}")  # 应该 = len(hybrid_atoms)
```

**常见原因**：

1. **映射返回 0 个 transformed atoms**：
   - 原因：ITP 文件原子名称与 PDB 不一致
   - 解决：使用 CHARMM-GUI 格式（`input/XX/gromacs/LIG.itp`）

2. **HybridTopologyBuilder 没有创建变换原子**：
   - 原因：mapping.transformed_a 或 mapping.transformed_b 为空
   - 解决：检查映射参数（dist_cutoff, charge_cutoff）

3. **ITP 写入时丢失原子**：
   - 原因：_get_atoms_section() 过滤逻辑问题
   - 解决：确保变换原子的 state_b_type != state_a_type

#### CHARMM-GUI vs CGenFF 网站格式

**CRITICAL**: 原子名称一致性是映射成功的关键！

**CHARMM-GUI 格式**（推荐）：
```
PDB 文件:  ATOM     25  F3  LIG L   1  ...
ITP 文件:      25       FGA2      1      LIG     F3     25  ...
                           ↑ 原子名一致 ✅
```

**CGenFF 网站格式**（不推荐）：
```
PDB 文件:  ATOM     25  F3  LIG L   1  ...
ITP 文件:      25       FGA2      1    42_    F25     25  ...
                           ↑ 原子名不一致 ❌
```

**正确路径配置**：
```python
# ✅ 正确
cgenff_dir_42 = test_dir / "input" / "42" / "gromacs"
cgenff_dir_38 = test_dir / "input" / "38" / "gromacs"

# ❌ 错误
cgenff_dir_42 = test_dir / "input" / "42_cgenff"
cgenff_dir_38 = test_dir / "input" / "38_cgenff"
```

## 必须避免的错误

1. **不要生成太多脚本版本**：
   - ❌ `local_run.sh`, `parallel_run.sh`, `single_run.sh`...
   - ✅ 一个 `run_fep.sh`，通过参数控制

2. **不要硬编码 GPU ID**：
   - ❌ `-gpu_id 0`
   - ✅ `-gpu_id ${GPU_ID:-0}` (可配置)

3. **不要忽略 bound/unbound 一致性**：
   - 两个 leg 的脚本应该几乎一样
   - 只改必要部分（输入文件、工作目录）

4. **不要跳过搭建后验证**：
   - 必须先 `gmx check`
   - 必须先 `gmx grompp` 成功（至少测试 EM）

5. **不要混淆 build 和 run 阶段**：
   - Build 阶段不创建 `window_XX/` 目录
   - Run 阶段动态创建窗口目录

## 相关代码

- 配置管理：`prism/fep/config.py::FEPConfig`
- 脚本生成：`prism/fep/modeling/script_writer.py::_write_fep_master_script()`
- MDP 生成：`prism/fep/modeling/mdp_generator.py::generate_production_mdps()`
- Lambda 调度：`prism/fep/modeling/mdp_generator.py::_generate_lambda_schedule()`
- 目录构建：`prism/fep/modeling/core.py::FEPScaffoldBuilder.build_from_components()`
- 原子映射：`prism/fep/gromacs/itp_builder.py::ITPBuilder`
