# PRISM PMF Module

## Overview

The PRISM PMF module provides comprehensive protein-ligand binding free energy calculations that are fully compatible with the PRISM project's architectural design patterns. This module refactors the functionality from the original pmf.py script into a modular, maintainable codebase that follows PRISM standards.

## Core Features

### 🏗️ PRISM Architecture Compatibility
- **Unified Interface**: Follows PRISM's `system.build()` design pattern
- **Modular Design**: Cleanly separated component managers
- **Configuration Management**: Fully compatible with PRISM configuration system
- **Error Handling**: Unified exception handling and logging

### 🔄 Dual Execution Modes
- **Automated Mode**: One-click execution of complete PMF calculation workflow
- **Step Control Mode**: Manual control over each calculation step
- **State Management**: Complete workflow state tracking and recovery

### 🏗️ PMF System Builder
- **Smart Reconstruction**: Extract and reconstruct systems from MD results
- **Geometry Optimization**: Automatic center-of-mass calculation and Z-axis alignment
- **Box Extension**: Automatic Z-axis extension to accommodate pulling process
- **Structure Extraction**: Automatic separation of Protein and LIG components
- **Complete Equilibration**: Automatic EM→NVT→NPT full equilibration process

### 📊 Complete Analysis Functionality
- **SMD Simulation**: Automated setup and execution
- **Umbrella Sampling**: Adaptive window selection and parallel execution
- **WHAM Analysis**: Comprehensive statistical analysis and error estimation
- **Visualization**: Professional-grade result plots and reports

## Quick Start

### Basic Usage (PRISM Style)

```python
import prism

# 1. Build system (using PRISM standard workflow)
system = prism.system("protein.pdb", "ligand.mol2")
build_results = system.build()

# 2. PMF calculation (integrated into PRISM workflow)
pmf = prism.pmf.pmf_system(build_results['output_dir'], "./pmf_results")

# Automated execution
results = pmf.run(mode='auto')
print(f"Binding energy: {results['binding_energy']['value']:.2f} kcal/mol")
```

### Using PMF System Builder

```python
import prism

# 1. Reconstruct PMF-optimized system from existing MD results
pmf = prism.pmf.pmf_system("./existing_md_results", "./pmf_optimized", 
                          rebuild_system=True)

# 2. Reconstruct system (auto extract, align, extend box, equilibrate)
rebuild_results = pmf.rebuild_system(frame=-1)  # Use last frame, auto equilibrate
print(f"Initial distance: {rebuild_results['alignment']['initial_distance']:.3f} nm")
print(f"System aligned to Z-axis, box extended, EM/NVT/NPT equilibration completed")

# 3. Run PMF calculation 
results = pmf.run(mode='auto')
print(f"Binding energy: {results['binding_energy']['value']:.2f} kcal/mol")
```

### Direct PMF Builder Usage

```python
import prism

# Independent PMF builder usage
builder = prism.pmf.pmf_builder("./md_results", "./pmf_system")

# Configure build parameters
config = {
    'reference_group': 'Protein',
    'moving_group': 'LIG', 
    'box': {'z_extension': 3.0}  # Larger Z-axis extension
}

# Build PMF-optimized system (with complete equilibration)
build_results = builder.build(equilibrate=True)
print(f"PMF system built and equilibrated: {build_results['system_dir']}")
print(f"Equilibration status: {build_results['status']}")  # pmf_ready_equilibrated

# Or step-by-step equilibration control
build_results = builder.build(equilibrate=False)  # Build only, no equilibration
print("=== Manual equilibration control ===")
em_results = builder.run_equilibration_step('em')    # Energy minimization
nvt_results = builder.run_equilibration_step('nvt')  # NVT equilibration
npt_results = builder.run_equilibration_step('npt')  # NPT equilibration

# Get final equilibrated system
final_system = builder.get_equilibrated_system()
print(f"Final equilibrated system: {final_system['structure']}")
```

### Equilibration Process Control

```python
import prism

# Detailed equilibration process control
builder = prism.pmf.pmf_builder("./md_results", "./pmf_system")

# 1. Build system geometry only
build_results = builder.build(equilibrate=False)
print("System geometry built, starting equilibration process...")

# 2. Check equilibration status
status = builder.get_equilibration_status()
print(f"Current status: {status['status']}")

# 3. Step-by-step equilibration
if status['em_ready']:
    em_results = builder.run_equilibration_step('em')
    print(f"✅ Energy minimization completed: {em_results['analysis']['final_energy']}")

if status['nvt_ready']:
    nvt_results = builder.run_equilibration_step('nvt') 
    print(f"✅ NVT equilibration completed: {nvt_results['analysis']['avg_temperature']} K")

if status['npt_ready']:
    npt_results = builder.run_equilibration_step('npt')
    print(f"✅ NPT equilibration completed: {npt_results['analysis']['avg_pressure']} bar")

# 4. Get final system
equilibrated_system = builder.get_equilibrated_system()
print(f"Equilibration complete! Final system: {equilibrated_system['structure']}")
```

### Step Control Mode

```python
import prism

# Build system using PRISM
system = prism.system("protein.pdb", "ligand.mol2")
system_output = system.build()

# Create PMF system
pmf = prism.pmf.PMFSystem(system_output['output_dir'], "./pmf_output")

# Step 1: Prepare SMD
print("=== Step 1: SMD Preparation ===")
smd_results = pmf.build(step='smd')
print(f"Execute SMD: {smd_results['instructions']}")

# Run SMD manually then continue
# bash {pmf_output}/smd/run_smd.sh

# Step 2: Prepare umbrella sampling  
print("=== Step 2: Umbrella Sampling Preparation ===")
umbrella_results = pmf.build(step='umbrella')
print(f"Execute umbrella sampling: {umbrella_results['instructions']}")

# Run umbrella sampling manually then continue
# bash {pmf_output}/umbrella/run_all_umbrella.sh parallel

# Step 3: Analysis
print("=== Step 3: WHAM Analysis ===")
analysis_results = pmf.build(step='analysis')
```

### Advanced Configuration

```python
import prism

# Custom configuration
config = {
    'reference_group': 'Protein',
    'moving_group': 'LIG',
    'smd': {
        'pull_rate': 0.01,      # Faster pulling velocity
        'pull_k': 2000.0        # Stronger force constant
    },
    'umbrella': {
        'sample_interval_near': 0.05,  # Dense sampling at close distances
        'production_time_ps': 30000    # Extended sampling time
    }
}

pmf = prism.pmf.pmf_system("./gaff_model", "./pmf_results", config=config)
results = pmf.run(mode='manual')  # Step control mode
```

## Architecture Design

### Module Structure

```
prism/pmf/
├── core.py          # High-level PMF interface (similar to prism.core)
├── workflow.py      # Workflow manager
├── pmf_builder.py   # PMF system builder (reconstruct from MD results)
├── equilibration.py # System equilibration manager (EM/NVT/NPT)
├── smd.py          # SMD simulation management
├── umbrella.py     # Umbrella sampling management
├── analyzer.py     # WHAM analysis and visualization
└── utils.py        # Utility functions
```

### Core Class Design

#### `PMFSystem` (High-level Interface)
```python
class PMFSystem:
    """High-level PMF interface, similar to PRISMSystem"""
    
    def __init__(self, system_dir, output_dir, config=None)
    def build(self, step=None)        # Prepare calculation components
    def build_step_by_step()          # Step-controlled building
    def run(self, mode='auto')        # Execute calculations
    def get_status()                  # Get status
    def analyze_smd()                 # SMD analysis
    def analyze_pmf()                 # PMF analysis
```

#### `PMFWorkflow` (Workflow Management)
```python
class PMFWorkflow:
    """PMF workflow manager, similar to PRISMBuilder"""
    
    def prepare_smd()                 # Prepare SMD
    def prepare_umbrella()            # Prepare umbrella sampling
    def prepare_analysis()            # Prepare analysis
    def run_analysis()                # Execute WHAM analysis
    def run_complete()                # Complete workflow
```

#### `PMFBuilder` (System Builder)
```python
class PMFBuilder:
    """PMF system builder, reconstruct optimized system from MD results"""
    
    def build(frame=None, equilibrate=True)  # Build PMF-optimized system
    def extract_structures()          # Extract protein and ligand structures
    def calculate_alignment_geometry() # Calculate center-of-mass and alignment vectors
    def rebuild_aligned_system()      # Reconstruct aligned and extended system
    def run_equilibration()           # Run complete equilibration process (EM→NVT→NPT)
    def run_equilibration_step(step)  # Run individual equilibration step
    def get_equilibration_status()    # Get equilibration status
    def get_equilibrated_system()     # Get final equilibrated system
    def get_system_info()             # Get system information
    def clean()                       # Clean temporary files
```

#### `PMFEquilibrationManager` (Equilibration Manager)
```python
class PMFEquilibrationManager:
    """PMF system equilibration manager, handles EM/NVT/NPT complete equilibration workflow"""
    
    def run_complete_equilibration()  # Run complete equilibration workflow
    def run_energy_minimization()     # Run energy minimization
    def run_nvt_equilibration()       # Run NVT equilibration
    def run_npt_equilibration()       # Run NPT equilibration
    def get_equilibration_status()    # Get equilibration status
    def get_final_equilibrated_system() # Get final equilibrated system
```

## Integration with PRISM

### Seamless Integration Example

```python
import prism

# Standard PRISM workflow
def complete_pmf_workflow(protein_file, ligand_file, output_base):
    """Complete PRISM + PMF workflow"""
    
    # 1. System building (PRISM standard)
    system = prism.system(protein_file, ligand_file, 
                         output_dir=f"{output_base}/system")
    build_results = system.build()
    print(f"✓ System building completed: {build_results['output_dir']}")
    
    # 2. PMF calculation (integrated)
    pmf = prism.pmf.pmf_system(
        build_results['output_dir'], 
        f"{output_base}/pmf"
    )
    pmf_results = pmf.run(mode='auto')
    print(f"✓ PMF calculation completed: {pmf_results['binding_energy']['value']:.2f} kcal/mol")
    
    return {
        'system': build_results,
        'pmf': pmf_results
    }

# Usage example
results = complete_pmf_workflow("protein.pdb", "ligand.mol2", "./analysis")
```

### Configuration System Compatibility

```python
# Using PRISM configuration file
pmf = prism.pmf.pmf_system("./system", "./pmf", config="pmf_config.yaml")

# Or using PRISM configuration dictionary
config = prism.config.load_config("prism_config.yaml")
pmf_config = config.get('pmf', {})
pmf = prism.pmf.pmf_system("./system", "./pmf", config=pmf_config)
```

## Output Structure

The PMF module follows PRISM's output organization pattern:

```
pmf_output/
├── smd/                    # SMD results
│   ├── run_smd.sh         # Execution script
│   ├── smd.mdp            # SMD parameters
│   ├── results/           # SMD results
│   └── analysis/          # SMD analysis
├── umbrella/              # Umbrella sampling
│   ├── run_all_umbrella.sh
│   ├── window_000/        # Sampling windows
│   ├── window_001/
│   └── ...
├── analysis/              # Final analysis
│   ├── pmf.xvg           # PMF data
│   ├── pmf_curve.png     # PMF plot
│   └── pmf_analysis_report.txt
└── workflow_state.yaml   # Workflow state
```

## State Management

### Status Monitoring

```python
pmf = prism.pmf.pmf_system("./system", "./pmf")

# Get detailed status
status = pmf.get_status()
print(f"Current stage: {status['current_stage']}")
print(f"Next action: {status['next_action']}")
print(f"File status: {status['files']}")
```

### Recovery and Cleanup

```python
# Clean specific components
pmf.clean(['smd'])              # Clean SMD results
pmf.clean(['umbrella', 'analysis'])  # Clean subsequent steps

# State recovery
pmf = prism.pmf.pmf_system("./system", "./pmf")  # Auto-recover state
```

## Error Handling

The PMF module adopts PRISM's error handling pattern:

```python
try:
    pmf = prism.pmf.pmf_system("./system", "./pmf")
    results = pmf.run(mode='auto')
except FileNotFoundError as e:
    print(f"System files missing: {e}")
except RuntimeError as e:
    print(f"Workflow error: {e}")
    # Get current status for debugging
    status = pmf.get_status()
    print(f"Stopped at: {status['current_stage']}")
```

## Advanced Features

### Custom Analysis

```python
# SMD results analysis
smd_analysis = pmf.analyze_smd()
print(f"SMD analysis completed: {smd_analysis['plots']}")

# Detailed PMF analysis
pmf_analysis = pmf.analyze_pmf()
print(f"PMF statistics: {pmf_analysis['pmf_statistics']}")
```

### Parameter Validation

```python
from prism.pmf import SMDBuilder

# Parameter validation
validation = SMDBuilder.validate_smd_parameters(config)
if not validation['valid']:
    print(f"Parameter issues: {validation['issues']}")
    print(f"Suggestions: {validation['warnings']}")
```

## Performance and Best Practices

### Recommended Workflow

1. **Development and Testing**: Use step control mode for easy debugging
2. **Production Calculations**: Use automated mode for efficiency
3. **Parameter Optimization**: Test with short times first, then long calculations
4. **Parallel Execution**: Use parallel mode for umbrella sampling

### System Requirements

- **GROMACS**: Version 2020+, with pull code support
- **Python Packages**: numpy, matplotlib, pyyaml
- **Computational Resources**: GPU recommended for umbrella sampling
- **Storage Space**: Sufficient space for trajectory storage

## Summary

The refactored PRISM PMF module provides the following advantages:

1. **Architectural Consistency**: Fully compliant with PRISM design patterns
2. **Code Quality**: Modular, maintainable, highly readable  
3. **Complete Functionality**: Complete workflow from SMD to WHAM
4. **Ease of Use**: High-level interface simplifies complex operations
5. **Flexibility**: Supports both automated and manual control modes
6. **Integration**: Seamless integration with PRISM ecosystem

This refactoring truly integrates PMF calculations into PRISM's core architecture, providing professional, reliable, and user-friendly PMF calculation capabilities!