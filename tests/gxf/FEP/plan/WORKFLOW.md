# PRISM-FEP Workflow Guide

Complete workflow for running Free Energy Perturbation calculations with PRISM.

## Table of Contents
1. [Quick Start](#quick-start)
2. [System Preparation](#system-preparation)
3. [FEP Building](#fep-building)
4. [Running Calculations](#running-calculations)
5. [Analysis](#analysis)

## Quick Start

### Minimal Example
```bash
# 1. Build hybrid topology and FEP system
prism protein.pdb ligand.mol2 -o output_dir -lff gaff2 --fep

# 2. Run bound leg
cd output_dir/GMX_PROLIG_FEP
./run_fep.sh bound

# 3. Run unbound leg
./run_fep.sh unbound

# 4. Analyze results
python -m prism.fep.analysis.cli --bound bound --unbound unbound
```

## System Preparation

### Input Requirements

#### Protein Structure
```bash
# Clean protein PDB
# - Remove waters, ions, ligands
# - Fix protonation states (HID/HIE/HIP for histidines)
# - Check for missing atoms/residues
```

#### Ligand Structures
```bash
# Reference ligand (state A)
ref_ligand.mol2  # or .pdb, .sdf

# Mutant ligand (state B)
mut_ligand.mol2   # or .pdb, .sdf
```

### Force Field Selection

```yaml
# Recommended: GAFF2 + AMBER14SB
forcefield: amber14sb
ligand_forcefield: gaff2

# Alternative options:
# - OpenFF + AMBER14SB
# - CHARMM36 + CGenFF
```

## FEP Building

### Step 1: Ligand Parameterization

PRISM automatically generates force field parameters for both ligands:

```bash
# GAFF2 (recommended)
prism -pf protein.pdb -lf ref.mol2 -lf mut.mol2 \
  -lff gaff2 -ff amber14sb -o output_dir
```

Output:
```
output_dir/
├── gaff2_ref/          # Reference ligand FF
│   └── LIG.amb2gmx/
│       ├── LIG.itp
│       ├── LIG.gro
│       └── atomtypes.itp
└── gaff2_mut/          # Mutant ligand FF
    └── LIG.amb2gmx/
        ├── LIG.itp
        ├── LIG.gro
        └── atomtypes.itp
```

### Step 2: Hybrid Topology Generation

```python
from prism.fep.modeling.core import FEPScaffoldBuilder

# Build FEP scaffold
builder = FEPScaffoldBuilder(
    output_dir="fep_output",
    lambda_windows=32,
    lambda_strategy="decoupled",
)

scaffold = builder.build_from_components(
    receptor_pdb="protein.pdb",
    hybrid_itp="hybrid.itp",
    reference_ligand_dir="gaff2_ref/LIG.amb2gmx",
    mutant_ligand_dir="gaff2_mut/LIG.amb2gmx",
)
```

Output structure:
```
fep_output/GMX_PROLIG_FEP/
├── common/
│   ├── protein/        # Receptor structure
│   └── hybrid/         # Hybrid topology files
│       ├── hybrid.itp
│       ├── atomtypes.itp
│       └── ligand_seed.pdb
├── bound/             # Bound leg (complex)
│   ├── input/
│   │   ├── conf.gro
│   │   ├── complex_seed.pdb
│   │   └── topol.top
│   ├── mdps/          # MDP parameter files
│   │   ├── em.mdp
│   │   ├── nvt.mdp
│   │   ├── npt.mdp
│   │   ├── npt_short_00.mdp
│   │   ├── npt_short_01.mdp
│   │   ├── ...
│   │   ├── prod_00.mdp
│   │   ├── prod_01.mdp
│   │   └── ...
│   └── build/         # Equilibration outputs
├── unbound/           # Unbound leg (ligand in water)
│   └── (same structure as bound/)
└── run_fep.sh         # Master execution script
```

### Step 3: Configuration

```yaml
# fep_config.yaml
fep:
  # Lambda schedule settings
  strategy: decoupled          # decoupled, coupled, or both
  lambda_windows: 32           # Total lambda windows
  distribution: nonlinear      # lambda point distribution
  
  # Replica runs for error estimation
  replicas: 3                  # Number of independent runs
  
  # Simulation settings
  production_time_ns: 10       # Production run length
  per_window_npt_time_ps: 100  # Per-window equilibration
  
  # Soft-core parameters
  soft_core_alpha: 0.5
  soft_core_sigma: 0.3

execution:
  mode: standard              # standard | repex
  # GPU configuration
  num_gpus: 4
  parallel_windows: 4          # Only used in standard mode
  omp_threads: 14              # OpenMP threads per GPU
  use_gpu_pme: true

simulation:
  # Standard PRISM settings
  temperature: 310
  pressure: 1.0
  dt: 0.002
```

## Running Calculations

### Master Script Usage

The master script (`run_fep.sh`) provides flexible execution control:

```bash
# Run specific legs
./run_fep.sh bound        # All bound replicas
./run_fep.sh unbound      # All unbound replicas
./run_fep.sh all          # Both legs

# Run specific replicas
./run_fep.sh bound1       # Single replica
./run_fep.sh bound1-3     # Range of replicas
./run_fep.sh bound1,unbound1  # Mixed targets
```

### Execution Stages

Each leg runs through these stages:

1. **Energy Minimization** (EM)
   - Steepest descent
   - Convergence: max force < 1000 kJ/mol/nm

2. **NVT Equilibration**
   - 500 ps (configurable)
   - Position restraints on heavy atoms
   - Temperature coupling (310 K)

3. **NPT Equilibration**
   - 500 ps (configurable)
   - Position restraints on heavy atoms
   - Pressure coupling (1 bar)

4. **Lambda Windows**
   - For each window (32 total):
     - Per-window NPT short (100 ps)
     - Production run (10 ns)

### Execution Modes

| Mode | 适用阶段 | 资源使用 | 说明 |
|------|----------|----------|------|
| `standard` | EM / NVT / NPT / Production | `parallel_windows` 个窗口并发；通常 1 个窗口占 1 张 GPU | 当前默认模式；适合普通多窗口并行 |
| `repex` | Production only | 一个 leg 的全部窗口组成 1 个 `gmx_mpi -multidir` 作业，共享 `num_gpus` 张 GPU | 仅生产阶段启用 λ 间副本交换；EM/NVT/NPT 仍逐窗口独立运行 |

**`standard` 示例**：`parallel_windows: 4`

```
Window 0 → GPU 0
Window 1 → GPU 1
Window 2 → GPU 2
Window 3 → GPU 3
Window 4 → GPU 0 (after window 0 completes)
...
```

**`repex` 示例**：`num_gpus: 4`

```bash
mpirun -oversubscribe -np 32 gmx_mpi mdrun \
  -deffnm prod -nb gpu -bonded gpu -pme gpu \
  -replex 1000 -multidir window_00 ... window_31
```

说明：
- `repex` 只改变运行脚本，不改变建模目录、MDP 模板或分析流程。
- `replex 1000`、`gmx_mpi`、`mpirun` 等高级细节默认写在脚本里，用户如需微调，直接修改生成脚本即可。

### Monitoring Progress

```bash
# Check running jobs
nvidia-smi                    # GPU utilization
tail -f bound/window_00/prod.log  # Log files

# Check completion
ls bound/window_*/prod.xvg   # Should show 32 files
```

对于 `repex` 模式，还应检查：

```bash
tail -f bound/prod.log
```

## Analysis

### Input Files

After completion, each lambda window produces:

```
bound/window_XX/
├── prod.xvg    # dH/dλ data (required for analysis)
├── prod.gro     # Final coordinates
├── prod.tpr     # Final checkpoint
└── prod.log     # Run log
```

### Running Analysis

```bash
# Basic analysis
python -m prism.fep.analysis.cli \
  --bound bound \
  --unbound unbound \
  --output results.html

# With replicas
python -m prism.fep.analysis.cli \
  --bound bound1,bound2,bound3 \
  --unbound unbound1,unbound2,unbound3 \
  --output results_replicas.html
```

### Output

Analysis produces:

- **ΔG**: Binding free energy with confidence intervals
- **Convergence**: Plot of ΔG vs simulation time
- **Overlap**: Pairwise histogram overlap matrix
- **Diagnostics**: Simulation quality metrics

## Best Practices

### 1. System Preparation
- Always check protein protonation states
- Verify ligand charges and protonation
- Remove crystallographic waters near binding site

### 2. Equilibration
- Never skip EM/NVT/NPT stages
- Check EM convergence (max force < 1000)
- Monitor system density during NPT

### 3. Lambda Windows
- 32 windows for decoupled strategy (recommended)
- 21 windows may suffice for simple transformations
- Check for "lambda gaps" in free energy profile

### 4. Replica Runs
- Use 3+ replicas for publication
- Start with 1 replica for testing
- Ensure replicas are independent (different velocities)

### 5. GPU Optimization
- Use 1 GPU per concurrent window
- Set `omp_threads` to match CPU cores per GPU
- Enable GPU PME for best performance

## Troubleshooting

### Common Issues

**EM not converging**
```
# Symptom: max force > 1000 kJ/mol/nm
# Solution: Increase em_nsteps or use steepest descent
```

**PME domain decomposition error**
```
# Symptom: "Domain decomposition was not successful"
# Solution: Increase box size or run longer equilibration
```

**Low GPU utilization**
```
# Symptom: GPU utilization < 50%
# Solution: Increase parallel_windows or check GPU settings
```

## Next Steps

- See [PARAMETERS.md](PARAMETERS.md) for detailed parameter reference
- See [GPU_OPTIMIZATION.md](GPU_OPTIMIZATION.md) for performance tuning
- See [ANALYSIS.md](ANALYSIS.md) for analysis details
