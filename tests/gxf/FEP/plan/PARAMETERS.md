# PRISM-FEP Parameter Reference

Complete reference for all FEP configuration parameters.

## Table of Contents
1. [Configuration File Format](#configuration-file-format)
2. [FEP Parameters](#fep-parameters)
3. [Execution Parameters](#execution-parameters)
4. [Simulation Parameters](#simulation-parameters)

## Configuration File Format

PRISM-FEP uses YAML configuration files for parameter management.

### File Locations

```
project_root/
├── fep_config.yaml           # Main FEP configuration
├── config_gaff.yaml          # Force field specific
└── prism/configs/default_config.yaml  # Default values
```

### Hierarchical Configuration

Parameters are merged in this priority order:
1. Command-line arguments (highest)
2. Project-specific config files
3. User config files
4. Default config files (lowest)

### Example Configuration

```yaml
# fep_config.yaml
fep:
  lambda_windows: 32
  replicas: 3

execution:
  num_gpus: 4
  parallel_windows: 4

simulation:
  temperature: 310
```

## FEP Parameters

### Lambda Strategy

```yaml
fep:
  strategy: decoupled    # Lambda transformation strategy
```

**Options**:
- `decoupled` (recommended): Separate electrostatics and VDW
  - Coulomb: 0 → 1 (windows 0-12)
  - VDW: 0 → 1 (windows 13-31)
  - Better convergence for large transformations
- `coupled`: Simultaneous electrostatics and VDW
  - Both: 0 → 1 (all windows)
  - Faster but may have hysteresis
- `both`: Run both strategies
  - Provides consistency check

**Default**: `decoupled`

### Lambda Windows

```yaml
fep:
  lambda_windows: 32          # Total number of lambda windows
  coul_windows: 12            # Coulomb windows (decoupled only)
  vdw_windows: 20             # VDW windows (decoupled only)
```

**Recommendations**:
- **Simple transformation** (small molecule change): 21 windows
- **Medium transformation** (functional group change): 32 windows (default)
- **Complex transformation** (ring change, large rearrangement): 43+ windows

**Distribution**:
```yaml
fep:
  distribution: nonlinear     # Lambda point distribution
```

**Options**:
- `linear`: Evenly spaced points
- `nonlinear`: Geometric spacing (recommended)
  - More points near λ=0 and λ=1
  - Better captures end-point singularities

### Replica Runs

```yaml
fep:
  replicas: 3                  # Number of independent runs
```

**Purpose**: Estimate statistical error and improve reliability

**Recommendations**:
- **Testing**: `replicas: 1`
- **Production**: `replicas: 3` (minimum for publication)
- **High precision**: `replicas: 5+`

**Directory Structure**:
```
GMX_PROLIG_FEP/
├── bound/                     # replicas: 1
├── bound1, bound2, bound3/    # replicas: 3
```

### Simulation Time

```yaml
fep:
  production_time_ns: 10       # Production run length per window (ns)
  per_window_npt_time_ps: 100  # Per-window NPT equilibration (ps)
```

**Guidelines**:
- **Quick testing**: 0.5-1 ns per window
- **Standard production**: 5-10 ns per window
- **High precision**: 20+ ns per window

**Per-window equilibration**:
- **Fast**: 50 ps
- **Standard**: 100 ps
- **Conservative**: 200 ps

### Soft-Core Parameters

```yaml
fep:
  soft_core_alpha: 0.5         # Soft-core alpha parameter
  soft_core_sigma: 0.3         # Soft-core sigma (nm)
```

**Purpose**: Prevent singularities when particles appear/disappear

**Default values** (recommended):
- `alpha: 0.5` (GROMACS default)
- `sigma: 0.3` nm (GROMACS default)

**Adjust if**:
- **Particle collapse**: Increase `alpha` or `sigma`
- **Poor convergence**: Decrease `alpha` or `sigma`

## Execution Parameters

### GPU Configuration

```yaml
execution:
  mode: standard              # standard | repex
  num_gpus: 4                  # Total number of available GPUs
  parallel_windows: 4          # Concurrent lambda windows in standard mode
  omp_threads: 14              # OpenMP threads per GPU
```

| 参数 | `standard` 模式 | `repex` 模式 |
|------|------------------|--------------|
| `mode` | `standard` | `repex` |
| `num_gpus` | 总 GPU 数；默认也作为 `parallel_windows` 的推断上限 | 总 GPU 数；整条 leg 的 multidir 作业共享这些 GPU |
| `parallel_windows` | 生效；控制同时跑多少个窗口 | 不生效；`repex` 下全部窗口进入一个 `gmx_mpi -multidir` 任务 |
| `omp_threads` | 每个窗口任务的 OpenMP 线程数 | 每个 MPI rank 的 OpenMP 线程数 |

**建议设置**:
- **1 GPU system**: `mode: standard`, `parallel_windows: 1`
- **4 GPU system, normal windows**: `mode: standard`, `parallel_windows: 4`
- **4 GPU system, lambda exchange**: `mode: repex`, `num_gpus: 4`

**Thread count**: 
```
omp_threads = CPU_cores_available / active_workers
```

说明：
- `repex` 需要的 MPI launcher、`replex 1000`、`gpu_id` 拼接等细节默认写入脚本，不作为 YAML 参数暴露。
- 如果用户需要修改这些高级执行细节，直接编辑生成的 `run_prod_repex.sh` 即可。

### GPU PME

```yaml
execution:
  use_gpu_pme: true            # Enable GPU PME calculation
```

**Recommendations**:
- **Large systems** (>50k atoms): `true` (faster)
- **Small systems** (<20k atoms): `false` (CPU may be faster)
- **Mixed precision**: Always use GPU PME

### GPU Bonded

```yaml
execution:
  use_gpu_bonded: true         # Enable GPU bonded interactions
```

**Recommendation**: Always `true` for modern GPUs

### Fallback to CPU

```yaml
execution:
  fallback_to_cpu: true         # Fall back to CPU if GPU fails
```

**Purpose**: Automatic error recovery

**Default**: `true`

## Simulation Parameters

### Basic Settings

```yaml
simulation:
  dt: 0.002                    # Time step (ps)
  temperature: 310              # Temperature (K)
  pressure: 1.0                 # Pressure (bar)
```

**Constraints**:
- `dt: 0.002` ps (2 fs) - maximum for hydrogen mass repartitioning
- `dt: 0.004` ps (4 fs) - if using LINCS with mass repartitioning

### Equilibration

```yaml
simulation:
  em_nsteps: 100000             # Energy minimization steps
  emtol: 1000.0                # EM convergence tolerance (kJ/mol/nm)
  equilibration_nvt_time_ps: 500   # NVT equilibration (ps)
  equilibration_npt_time_ps: 500   # NPT equilibration (ps)
```

**EM Settings**:
- `emtol < 1000`: Good convergence
- `emtol > 1000`: May need more minimization

**Equilibration time**:
- **Fast**: 200 ps each
- **Standard**: 500 ps each
- **Conservative**: 1000 ps each

### Thermostat and Barostat

```yaml
simulation:
  # Thermostat (V-rescale)
  tc-grps: "Protein Water_NA Water_ions Non-Protein"
  tau-t: 0.5
  ref-t: 310

  # Barostat (Parrinello-Rahman)
  tau-p: 2.0
  compressibility: 4.5e-5
  ref-p: 1.0
```

**Recommendations**:
- **NVT**: Use V-rescale thermostat
- **NPT**: Use Parrinello-Rahman barostat
- **Production**: Semi-isotropic pressure coupling for membrane

## Force Field Parameters

### Protein Force Field

```yaml
forcefield: amber14sb          # Protein force field
```

**Options**:
- `amber14sb`: AMBER ff14SB (recommended)
- `amber99sb-ildn`: AMBER ff99SB-ILDN
- `charmm36`: CHARMM36m

### Ligand Force Field

```yaml
ligand_forcefield: gaff2       # Ligand force field
```

**Options**:
- `gaff2`: GAFF2 (recommended)
- `gaff`: GAFF (legacy)
- `openff`: Open Force Field
- `oplsaa`: OPLS-AA
- `mmff94`: MMFF94

### Water Model

```yaml
water_model: tip3p              # Water model
```

**Options**:
- `tip3p`: TIP3P (default)
- `tip4p`: TIP4P-Ew
- `spce`: SPC/E
- `opc`: OPC (with amber19sb)

## Advanced Parameters

### Cutoffs

```yaml
simulation:
  rcoulomb: 1.0                # Coulomb cutoff (nm)
  rvdw: 1.0                    # VDW cutoff (nm)
  rlist: 1.0                   # Neighbor list cutoff (nm)
```

**Constraints**: `rvdw ≥ rcoulomb ≥ rlist`

### PME

```yaml
simulation:
  pme: yes                     # Enable PME
  pme-order: 4                 # PME interpolation order
  fourierspacing: 0.12         # PME grid spacing (nm)
```

**Recommendations**:
- `pme-order: 4` (default)
- `fourierspacing: 0.12` (default)

### Constraints

```yaml
simulation:
  constraints: h-bonds         # Bond constraints
  constraint_algorithm: lincs    # Constraint algorithm
```

**Options**:
- `h-bonds`: Bonds involving H only
- `all-bonds`: All bonds
- `none`: No constraints

## Parameter Validation

### Required Parameters

These parameters must be set:
```yaml
fep:
  lambda_windows: int > 0
  strategy: str

simulation:
  temperature: float > 0
  pressure: float > 0
```

### Conflicting Parameters

Avoid these conflicts:
```yaml
# BAD: Over-constrained system
simulation:
  dt: 0.005                    # Too large for H-bonds
  constraints: h-bonds

# GOOD: Consistent parameters
simulation:
  dt: 0.002                    # Appropriate
  constraints: h-bonds
```

## Parameter Tuning Guide

### Convergence Issues

**Problem**: Free energy not converged
```yaml
# Solution 1: Increase simulation time
fep:
  production_time_ns: 20        # Double from 10 ns

# Solution 2: More lambda windows
fep:
  lambda_windows: 43            # Increase from 32

# Solution 3: Better lambda distribution
fep:
  distribution: nonlinear       # Use geometric spacing
```

### Performance Issues

**Problem**: Simulation too slow
```yaml
# Solution 1: Optimize GPU usage
execution:
  parallel_windows: 4          # Increase concurrency
  use_gpu_pme: true

# Solution 2: Reduce per-window equilibration
fep:
  per_window_npt_time_ps: 50   # Reduce from 100 ps

# Solution 3: Fewer lambda windows
fep:
  lambda_windows: 21            # Reduce from 32
```

### Stability Issues

**Problem**: Simulation crashes or instability
```yaml
# Solution 1: Increase soft-core
fep:
  soft_core_alpha: 1.0         # Increase from 0.5
  soft_core_sigma: 0.3         # Keep at 0.3

# Solution 2: Longer equilibration
simulation:
  equilibration_npt_time_ps: 1000  # Increase from 500

# Solution 3: Smaller timestep
simulation:
  dt: 0.001                    # Reduce from 0.002
```

## Best Practices

### 1. Start Simple
```yaml
# Initial testing
fep:
  lambda_windows: 21
  replicas: 1
  production_time_ns: 1
```

### 2. Scale Up
```yaml
# Production
fep:
  lambda_windows: 32
  replicas: 3
  production_time_ns: 10
```

### 3. Validate
```yaml
# Check consistency
fep:
  strategy: both               # Run both decoupled and coupled
```

### 4. Document
```yaml
# Always keep a copy of your config
# Include in version control
# Document rationale for non-default choices
```

## Reference Tables

### Lambda Window Guidelines

| Transformation Type | Complexity | Recommended Windows |
|-------------------|------------|---------------------|
| Small molecule | Low | 21 |
| Functional group | Medium | 32 |
| Ring change | High | 43 |
| Large rearrangement | Very High | 50+ |

### Simulation Time Guidelines

| Goal | Time per Window | Total Time (32 windows) |
|------|-----------------|------------------------|
| Quick test | 0.5 ns | 16 ns |
| Standard | 5 ns | 160 ns |
| High quality | 10 ns | 320 ns |
| Publication | 20 ns | 640 ns |

### GPU Utilization Guidelines

| GPUs | parallel_windows | Estimated throughput |
|------|------------------|---------------------|
| 1 | 1 | 5-10 ns/day |
| 4 | 4 | 20-40 ns/day |
| 8 | 8 | 40-80 ns/day |
