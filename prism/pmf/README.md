# PRISM PMF Module

The PMF (Potential of Mean Force) module provides comprehensive tools for calculating protein-ligand binding free energies through steered molecular dynamics (SMD) and umbrella sampling methodologies. This module integrates seamlessly with PRISM's architecture to deliver accurate binding energy calculations with automated workflows.

## Overview

The PMF module calculates binding free energies by:
1. **System Optimization**: Rebuilding MD systems with PMF-optimized geometry and extended simulation boxes
2. **Steered MD (SMD)**: Pulling the ligand away from the protein at constant velocity to generate reaction coordinates
3. **Umbrella Sampling**: Running multiple simulations at different protein-ligand distances with harmonic restraints  
4. **WHAM Analysis**: Using the Weighted Histogram Analysis Method to calculate the potential of mean force
5. **Visualization**: Generating comprehensive plots and reports of binding energetics

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

# Using an existing MD system
system = pmf.pmf_system("./my_md_output", "./pmf_results")

# Run complete automated workflow
results = system.run(mode='auto')
print(f"Binding energy: {results['binding_energy']['value']:.2f} kcal/mol")
```

### Step-by-Step Execution

```python
import prism.pmf as pmf

# Initialize PMF system
system = pmf.pmf_system("./my_md_output", "./pmf_results")

# Step 1: Prepare SMD simulation
smd_results = system.build(step='smd')
# Run manually: cd ./pmf_results/smd && bash run_smd.sh

# Step 2: Prepare umbrella sampling (after SMD completion)
umbrella_results = system.build_umbrella_step()
# Run manually: cd ./pmf_results/umbrella && bash run_all_umbrella.sh parallel

# Step 3: Perform WHAM analysis (after umbrella completion)
analysis_results = system.build_analysis_step()
```

### High-Level Runner API

```python
import prism.pmf as pmf

# Run complete workflow with runner
results = pmf.run_pmf_workflow(
    md_system_dir="./my_md_output",
    output_dir="./pmf_results",
    config="./my_pmf_config.yaml"
)

# Run individual steps
smd_results = pmf.run_pmf_step("./system", "./output", "smd")
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
  z_extension: 2.5  # Extra box length for pulling (nm)
  
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

```
pmf_results/
├── smd/                          # SMD simulation
│   ├── smd.mdp                  # SMD parameters
│   ├── run_smd.sh               # Execution script
│   └── results/
│       ├── smd_pullf.xvg        # Force vs time
│       ├── smd_pullx.xvg        # Distance vs time
│       └── smd.gro              # Final structure
├── umbrella/                     # Umbrella sampling
│   ├── window_0/                # Individual windows
│   ├── window_1/
│   ├── ...
│   └── run_all_umbrella.sh      # Parallel execution script
├── analysis/                     # WHAM analysis
│   ├── pmf.xvg                  # PMF data
│   ├── pmferror.xvg             # Error estimates
│   ├── pmf_curve.png            # PMF plot
│   ├── force_profile.png        # Force vs distance
│   ├── energy_landscape.png     # Energy surface
│   └── pmf_analysis_report.txt  # Summary report
└── pmf_config_used.yaml         # Configuration used
```

## Advanced Usage

### System Rebuilding with PMF Optimization

For optimal PMF calculations, rebuild your system with PMF-specific settings:

```python
# Enable system rebuilding during PMF setup
system = pmf.pmf_system(
    "./my_md_output", 
    "./pmf_results",
    rebuild_system=True  # Rebuild with PMF optimization
)

# Or rebuild manually with specific settings
rebuild_results = system.rebuild_system(frame=-1)  # Use last MD frame
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
pmf_system = prism.pmf.pmf_system(build_results['output_dir'], "./pmf_results")
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