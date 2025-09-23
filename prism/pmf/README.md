# PRISM PMF Module

The PMF (Potential of Mean Force) module provides comprehensive tools for calculating protein-ligand binding free energies through steered molecular dynamics (SMD) and umbrella sampling methodologies. This module integrates seamlessly with PRISM's architecture to deliver accurate binding energy calculations with automated workflows.

## Overview

The PMF module calculates binding free energies using a optimized workflow:
1. **File Preparation**: Copy essential files (topology, restraints) from MD results
2. **Structure Preprocessing**: Extract proper frame, remove periodic boundary conditions, center system
3. **SMD Simulation**: Pull the ligand away from the protein at constant velocity
4. **Umbrella Sampling**: Run multiple simulations at different protein-ligand distances with harmonic restraints
5. **WHAM Analysis**: Use the Weighted Histogram Analysis Method to calculate the potential of mean force
6. **Visualization**: Generate comprehensive plots and reports of binding energetics
7. **Optional System Rebuilding**: Advanced system modifications if needed (rarely required)

## Key Features

- **Automated Workflow**: Complete PMF calculations with minimal user intervention
- **Force Field Support**: Compatible with GAFF and OpenFF force fields used by PRISM
- **Adaptive Window Selection**: Intelligent selection of umbrella sampling windows based on SMD results
- **Error Analysis**: Bootstrap error estimation for statistical reliability
- **GPU Acceleration**: Full support for GPU-accelerated simulations
- **Comprehensive Visualization**: Automated generation of PMF curves, force profiles, and energy landscapes

## Quick Start

### Basic Usage

```python
import prism.pmf as pmf

# Recommended: PMF calculation with system rebuilding for optimal results
# Uses default settings with 2.0 nm pulling distance
system = pmf.pmf_system("./gaff_model", "./pmf_results")

# Customize pulling distance for your specific system
system = pmf.pmf_system(
    "./gaff_model", "./pmf_results",
    pulling_distance=3.5  # Custom 3.5 nm pulling distance
)

# Run complete automated workflow
results = system.run(mode='auto')
print(f"Binding energy: {results['binding_energy']['value']:.2f} kcal/mol")
```

### Step-by-Step Execution

```python
import prism.pmf as pmf

# Initialize PMF system with custom pulling distance
system = pmf.pmf_system(
    "./gaff_model", "./pmf_results",
    pulling_distance=2.5  # Customize pulling distance for your system
)

# Step 1: Prepare and rebuild system (automatic)
# - Copy essential files (topology, TPR, position restraints)
# - Preprocess structure (remove PBC, center system)
# - Create PMF-optimized box with GROMACS built-in algorithms
# - Rebuild system with extended box for pulling
smd_results = system.build(step='smd')
# Run manually: cd ./gaff_model/pmf_smd && bash run_smd.sh

# Step 2: Prepare umbrella sampling (after SMD completion)
umbrella_results = system.build_umbrella_step()
# Run manually: cd ./gaff_model/pmf_umbrella && bash run_all_umbrella.sh parallel

# Step 3: Perform WHAM analysis (after umbrella completion)
analysis_results = system.build_analysis_step()
```

### High-Level Runner API

```python
import prism.pmf as pmf

# Run complete workflow with runner
results = pmf.run_pmf_workflow(
    md_system_dir="./gaff_model",
    output_dir="./pmf_results",
    config="./my_pmf_config.yaml"
)

# Run individual steps  
smd_results = pmf.run_pmf_step("./gaff_model", "./pmf_results", "smd")
```

## Configuration

### Default Configuration

The module uses sensible defaults for production PMF calculations:

- **SMD**: 5 ns simulation with 0.005 nm/ps pull rate
- **Umbrella Sampling**: 22 ns per window (2 ns equilibration + 20 ns production)
- **Window Selection**: Adaptive spacing (0.1 nm close, 0.2 nm far)
- **Analysis**: 50 bootstrap iterations for error estimation

### Creating Custom Configurations

```python
# Create configuration templates
pmf.create_pmf_config("my_config.yaml", template="default")  # Standard
pmf.create_pmf_config("fast_config.yaml", template="fast")   # Quick testing
pmf.create_pmf_config("accurate_config.yaml", template="accurate")  # High accuracy
```

### Configuration Parameters

Key configuration sections:

```yaml
# System rebuilding settings
builder:
  rebuild_system: true
  pulling_distance: 2.0     # User-customizable pulling distance (nm)
  box_distance: 1.2         # Fixed dt distance for GROMACS standard box (nm)
  
# SMD simulation parameters  
smd:
  pull_rate: 0.005      # nm/ps
  pull_k: 1000.0        # kJ/mol/nm²
  nsteps: 2500000       # 5 ns simulation

# Umbrella sampling parameters
umbrella:
  sample_interval_near: 0.1     # nm, close distances
  sample_interval_far: 0.2      # nm, far distances
  production_time_ps: 20000     # 20 ns per window

# WHAM analysis parameters
analysis:
  begin_time_ps: 2000           # Skip first 2 ns
  bootstrap_iterations: 50      # Error analysis
  energy_unit: "kCal"          # Output unit
```

#### PMF Box Creation Algorithm

The PMF module uses GROMACS built-in algorithms for reliable box creation:

1. **Standard Box Creation**: Uses GROMACS `-d` option with fixed `box_distance` (default 1.2 nm)
2. **Box Extension**: Adds `2 × pulling_distance` to Z-axis dimension
3. **System Translation**: Moves system by `pulling_distance` in negative Z direction
4. **Result**: Optimal pulling space in +Z direction with safety margin in -Z direction

This ensures maximum available space for pulling simulations while maintaining system stability.

## System Requirements

### Input Requirements

1. **Completed MD System**: A PRISM-generated system directory containing:
   - `GMX_PROLIG_MD/solv_ions.gro` - System structure
   - `GMX_PROLIG_MD/topol.top` - System topology
   - MD trajectory files (optional, for system rebuilding)

2. **Software Dependencies**:
   - GROMACS (2018 or newer) with PMF tools
   - Python packages: numpy, matplotlib, scipy
   - Optional: GPU drivers for acceleration

### Output Structure

The PMF module generates results in a structured directory layout:

```
pmf_results/                     # Your specified output directory
├── analysis/                    # Final PMF analysis results
│   ├── pmf.xvg                 # PMF data file
│   ├── pmferror.xvg            # Error estimates
│   ├── pmf_curve.png           # PMF curve plot
│   ├── force_profile.png       # Force vs distance
│   ├── energy_landscape.png    # Energy surface visualization
│   └── pmf_analysis_report.txt # Summary report
└── pmf_config_used.yaml        # Configuration used for this calculation
```

**Working Files Location:**
- SMD and umbrella sampling simulations run alongside your MD directory
- Final analysis results are saved to your specified output directory
- This approach ensures efficient use of disk space and reliable access to all necessary files

## Advanced Usage

### Efficient File Management

PMF calculations efficiently manage files by working directly with your MD directory:

```python
import prism.pmf as pmf

# Default behavior - efficient file management
system = pmf.pmf_system("./gaff_model", "./pmf_results")
```

### Advanced Usage with Optional System Rebuilding

For special cases where you want to skip system rebuilding:

```python
# Advanced: Skip rebuilding (use only for special debugging cases)
system = pmf.pmf_system("./gaff_model", "./pmf_results", rebuild_system=False)
results = system.run(mode='auto')

# Standard recommended workflow (with rebuilding)
system = pmf.pmf_system("./gaff_model", "./pmf_results")  # rebuild_system=True by default

# Manual control over rebuilding process
smd_results = system.build(step='smd')          # Preprocesses files first
rebuild_results = system.rebuild_system()       # Rebuilds with PMF optimization
umbrella_results = system.build_umbrella_step() # Continue with umbrella sampling
analysis_results = system.build_analysis_step() # Final analysis
```

### Custom Analysis and Visualization

```python
# Analyze SMD results after completion
smd_analysis = system.analyze_smd()

# Detailed PMF analysis with custom plots
pmf_analysis = system.analyze_pmf()

# Export current configuration
system.export_config("./my_saved_config.yaml")
```

### Workflow Status and Monitoring

```python
# Check workflow status
status = system.get_status()
print(f"Current stage: {status['current_stage']}")
print(f"Next action: {status['next_action']}")

# Clean up specific components
system.clean(['smd'])  # Remove SMD results only
system.clean()         # Remove all PMF results
```

## Performance Optimization

### Hardware Recommendations

- **CPU**: Multi-core processor (8+ cores recommended)
- **GPU**: CUDA-compatible GPU for acceleration (optional but recommended)
- **Memory**: 16+ GB RAM for typical systems
- **Storage**: Fast SSD with sufficient space (50+ GB for complete workflow)

### Parallel Execution

```bash
# Umbrella sampling supports parallel execution
cd pmf_results/umbrella
bash run_all_umbrella.sh parallel  # Run windows in parallel
bash run_all_umbrella.sh sequential  # Run windows sequentially
```

### Time Estimates

Typical calculation times (protein-ligand complex with ~50,000 atoms):

- **SMD (5 ns)**: 1-6 hours (depending on hardware)
- **Umbrella Sampling**: Several hours to days (depends on number of windows and parallelization)
- **WHAM Analysis**: Minutes to hours (depends on bootstrap iterations)

## Troubleshooting

### Common Issues

1. **SMD Preparation Fails**
   ```
   Error: Could not locate GMX_PROLIG_MD directory
   ```
   - Ensure your input system was built with PRISM
   - Check that `solv_ions.gro` and `topol.top` exist

2. **Insufficient Pulling Distance**
   - Increase `z_extension` in configuration
   - Check that ligand fully unbinds during SMD

3. **Poor Window Selection**
   - Reduce `sample_interval_near` and `sample_interval_far`
   - Check SMD trajectory for proper unbinding pathway

4. **WHAM Analysis Fails**
   - Ensure at least 80% of umbrella windows completed successfully
   - Check umbrella sampling trajectories for proper convergence

### Performance Issues

- **Memory Errors**: Reduce system size or increase available memory
- **Slow Simulations**: Enable GPU acceleration or reduce simulation times
- **Disk Space**: Monitor storage usage, especially for trajectory files

### Getting Help

For technical issues:
1. Check log files in each step directory
2. Verify GROMACS installation and GPU drivers
3. Ensure all input files are properly formatted
4. Review configuration parameters for your system size

## API Reference

### Main Classes

#### `PMFSystem`
High-level interface for PMF calculations following PRISM patterns.

**Methods:**
- `build(step=None)` - Prepare PMF calculation components
- `run(mode='auto')` - Execute PMF workflow
- `rebuild_system(frame=None)` - Rebuild system with PMF optimization
- `get_status()` - Get workflow status
- `analyze_smd()` - Analyze SMD results
- `analyze_pmf()` - Analyze final PMF results
- `clean(components=None)` - Clean calculation results

#### `PMFBuilder`
System builder for PMF-optimized structures.

**Methods:**
- `build(equilibrate=True)` - Build PMF system with equilibration
- `extract_structures()` - Extract protein and ligand structures
- `rebuild_aligned_system()` - Reconstruct aligned system
- `run_equilibration()` - Run complete equilibration process

#### `PMFRunner`
High-level runner for complete PMF workflows.

**Methods:**
- `run_complete_workflow(md_system_dir, output_dir, steps=None)` - Run complete workflow
- `_run_builder_step()` - Execute system builder step
- `_run_smd_step()` - Execute SMD simulation step
- `_run_umbrella_step()` - Execute umbrella sampling step
- `_run_analysis_step()` - Execute WHAM analysis step

### Convenience Functions

- `pmf_system(system_dir, output_dir, config=None, rebuild_system=False, **kwargs)` - Create PMF system
- `run_pmf_workflow(md_system_dir, output_dir, config=None, steps=None)` - Run complete workflow
- `run_pmf_step(md_system_dir, output_dir, step, config=None)` - Run individual step
- `create_pmf_config(output_file, template="default")` - Create configuration file

## Integration with PRISM

The PMF module integrates seamlessly with PRISM workflows:

```python
import prism

# Build system with PRISM
system = prism.system("protein.pdb", "ligand.mol2")
build_results = system.build()

# Calculate PMF with integrated module
pmf_system = prism.pmf.pmf_system(
    build_results['output_dir'], 
    "./pmf_results"
)
pmf_results = pmf_system.run(mode='auto')

print(f"System built: {build_results['output_dir']}")
print(f"Binding energy: {pmf_results['binding_energy']['value']:.2f} kcal/mol")
```

## Citation

When using the PRISM PMF module in research, please cite:

```
PRISM: Protein Receptor Interaction Simulation Modeler
PMF Module for Binding Free Energy Calculations
```

## License

This module is part of the PRISM package and follows the same licensing terms.