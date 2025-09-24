# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PRISM (Protein Receptor Interaction Simulation Modeler) is a comprehensive Python tool for building protein-ligand systems for molecular dynamics simulations in GROMACS. It supports multiple force fields including GAFF (via AmberTools) and OpenFF (Open Force Field).

## Development Commands

### Installation and Setup

**Recommended Installation Method (2025)**: Use `mamba` for better dependency resolution

**Automatic Parallelization**: PRISM automatically enables OpenMP parallelization using all available CPU cores for trajectory analysis (MDTraj operations like DCD conversion). No additional configuration needed - just import PRISM and get accelerated performance.

```bash
# Install mamba (if not already installed)
conda install -c conda-forge mamba

# Install PRISM in development mode (basic analysis only)
pip install -e .

# Install with specific force field support (RECOMMENDED)
# Method 1: Using mamba (preferred for complex dependencies)
mamba install -c conda-forge openff-toolkit    # Includes openff-interchange
mamba install -c conda-forge ambertools        # For GAFF support

# Method 2: Using conda (alternative)
conda install -c conda-forge openff-toolkit ambertools

# Method 3: Using pip (works but not recommended for scientific dependencies)
pip install openff-toolkit openff-interchange

# Install PRISM with extras after dependencies
pip install -e .[openff]        # OpenFF support only
pip install -e .[gaff]          # GAFF support only
pip install -e .[all]           # All force fields

# Basic installation (analysis only, no force field generation)
pip install -r requirements.txt
```

**Installation Troubleshooting**: If you encounter conda dependency conflicts:

```bash
# Create isolated environment (HIGHLY RECOMMENDED)
mamba create -n prism-env -c conda-forge python=3.9 openff-toolkit ambertools
conda activate prism-env
pip install -e .

# Alternative: Use conda-forge channel priority
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install openff-toolkit ambertools

# For conda channel conflicts
conda install conda-forge::openff-toolkit conda-forge::ambertools

# Check installation
python -c "import prism as pm; print(pm.check_dependencies())"
```

### Running PRISM
```bash
# Basic usage with GAFF (default) - positional arguments
python prism/builder.py protein.pdb ligand.mol2 -o output_dir

# Using OpenFF force field
python prism/builder.py protein.pdb ligand.sdf -o output_dir --ligand-forcefield openff

# Alternative syntax using flags (more flexible argument order)
python prism/builder.py --protein-file protein.pdb --ligand-file ligand.mol2 -o output_dir

# Mixed order with flags (helpful for complex commands)
python prism/builder.py -o output_dir --ligand-forcefield openff --protein-file protein.pdb --ligand-file ligand.mol2

# With custom configuration
python prism/builder.py protein.pdb ligand.mol2 -o output_dir --config config.yaml

# List available force fields
python prism/builder.py --list-forcefields

# Export default configuration template
python prism/builder.py --export-config my_config.yaml
```

### Testing
```bash
# Test with example data (4xb4 complex)
cd test/4xb4
python ../../prism/builder.py 4xb4_protein.pdb 4xb4_ligand.mol2 -o test_output

# Run the generated simulation
cd test_output/GMX_PROLIG_MD
bash ../../../prism/sim/scripts/localrun.sh
```

### High-Level Python API
```python
import prism as pm

# Simple system building
system = pm.system("protein.pdb", "ligand.mol2")
output_dir = system.build()

# Step-by-step building with control
system = pm.system("protein.pdb", "ligand.sdf", ligand_forcefield="openff")
results = system.build_step_by_step()

# Running simulations
sim = pm.model("output/GMX_PROLIG_MD")
sim.run(engine="gmx")
```

## Architecture Overview

### Core Modules

1. **`prism.builder`** - Main entry point and workflow orchestration
   - `PRISMBuilder` class handles the complete workflow
   - Command-line interface with extensive options
   - Supports both GAFF and OpenFF force fields

2. **`prism.core`** - High-level Python API
   - `PRISMSystem` class provides simplified interface
   - Wraps `PRISMBuilder` functionality with user-friendly methods
   - Supports step-by-step building and configuration management

3. **`prism.forcefield/`** - Force field generation
   - `ForceFieldGeneratorBase` - Abstract base class
   - `GAFFForceFieldGenerator` - GAFF support via AmberTools/ACPYPE
   - `OpenFFForceFieldGenerator` - OpenFF support via openff-toolkit
   - Generates GROMACS-compatible .itp, .gro, and topology files

4. **`prism.utils/`** - Core utilities
   - `GromacsEnvironment` - GROMACS detection and management
   - `ConfigurationManager` - YAML configuration handling
   - `MDPGenerator` - MD parameter file generation
   - `SystemBuilder` - GROMACS system assembly

5. **`prism.sim/`** - Simulation execution
   - `SimulationModel` - High-level simulation interface
   - `GMXSimulator` - GROMACS-based simulations
   - `OpenMMSimulator` - OpenMM-based simulations
   - `localrun.sh` - Example GPU-accelerated simulation script

6. **`prism.pmf/`** - Advanced analysis
   - PMF calculation through steered MD and umbrella sampling
   - `SMDManager`, `UmbrellaManager` for specialized workflows
   - `PMFAnalyzer` for free energy calculations

### Configuration System

- **Default config**: `prism/configs/default_config.yaml`
- **Runtime config**: Generated in output directory as `prism_config.yaml`
- **Override hierarchy**: CLI args > config file > defaults
- **Key sections**: general, box, simulation, ions, constraints, energy_minimization, output, electrostatics, vdw, temperature_coupling, pressure_coupling

### Output Structure

Generated files follow this structure:
```
output_dir/
├── GMX_PROLIG_MD/           # Main GROMACS files
│   ├── solv_ions.gro        # Complete solvated system
│   ├── topol.top            # System topology
│   └── posre.itp            # Position restraints
├── LIG.amb2gmx/             # GAFF ligand files
│   ├── LIG.gro              # Ligand coordinates
│   ├── LIG.itp              # Ligand topology
│   ├── LIG.top              # Ligand topology file
│   ├── atomtypes_LIG.itp    # Atom type definitions
│   └── posre_LIG.itp        # Ligand position restraints
├── LIG.openff2gmx/          # OpenFF ligand files (alternative)
├── mdps/                    # MD parameter files
│   ├── em.mdp               # Energy minimization
│   ├── nvt.mdp              # NVT equilibration
│   ├── npt.mdp              # NPT equilibration
│   └── md.mdp               # Production MD
└── prism_config.yaml        # Used configuration
```

## Troubleshooting

### Common Issues

1. **Argument parsing errors**: Use correct order (protein ligand options) or alternative flags:
   ```bash
   # Correct: positional arguments first
   python prism/builder.py protein.pdb ligand.mol2 -o output_dir
   
   # Alternative: use flags for flexible order
   python prism/builder.py -o output_dir --protein-file protein.pdb --ligand-file ligand.mol2
   ```

2. **Ion addition failures**: The system now auto-detects water groups dynamically instead of using hardcoded values. If issues persist:
   - Check that the system contains water molecules
   - Verify GROMACS can find solvent groups in the topology
   - Check log output for specific error messages

3. **Force field errors**: Ensure dependencies are installed properly:

   **For GAFF support**:
   ```bash
   # Recommended installation
   mamba install -c conda-forge ambertools
   # Alternative
   conda install -c conda-forge ambertools
   # Last resort
   pip install ambertools  # May have dependency issues
   ```

   **For OpenFF support**:
   ```bash
   # Recommended installation (includes openff-interchange automatically)
   mamba install -c conda-forge openff-toolkit
   # Alternative
   conda install -c conda-forge openff-toolkit
   # Using pip (not recommended but works)
   pip install openff-toolkit openff-interchange
   ```

   **Common conda installation issues**:
   - Error: `ambertools does not exist`: Use `conda-forge::ambertools` or create new environment
   - Error: `PackagesNotFoundError`: Add conda-forge channel and set priority
   - Error: `Solving environment failed`: Use mamba instead of conda
   - Error: Channel conflicts: Create isolated environment with `mamba create -n newenv`

## Installation Troubleshooting Guide (2025)

### Common Installation Problems and Solutions

#### 1. **Conda/Mamba Channel Issues**

**Problem**: `PackagesNotFoundError: The following packages are not available from current channels`

**Solutions**:
```bash
# Solution A: Fix channel configuration
conda config --add channels conda-forge
conda config --set channel_priority strict

# Solution B: Use explicit channel syntax
conda install conda-forge::openff-toolkit conda-forge::ambertools

# Solution C: Create clean environment
mamba create -n prism-clean -c conda-forge python=3.9 openff-toolkit ambertools
conda activate prism-clean
pip install -e .
```

#### 2. **AmberTools "does not exist" Error**

**Problem**: `ambertools does not exist (perhaps a typo or a missing channel)`

**Solutions**:
```bash
# Recommended: Use mamba (better dependency solver)
mamba install -c conda-forge ambertools

# Alternative: Force conda-forge channel
conda install conda-forge::ambertools

# Check available versions
mamba search ambertools -c conda-forge
```

#### 3. **Environment Conflicts**

**Problem**: `Solving environment failed` or dependency conflicts

**Solutions**:
```bash
# Create isolated environment (BEST PRACTICE)
mamba create -n prism-env -c conda-forge python=3.9
conda activate prism-env
mamba install -c conda-forge openff-toolkit ambertools
pip install -e .

# Reset base environment if corrupted
conda clean --all
conda update conda
```

#### 4. **Platform-Specific Issues**

**macOS Apple Silicon**:
```bash
# Use native arm64 architecture
mamba create -n prism-env -c conda-forge python=3.9 openff-toolkit ambertools
# Avoid Rosetta/x86_64 unless absolutely necessary
```

**Windows Users**:
```bash
# Use Windows Subsystem for Linux (WSL)
# OpenFF does not officially support native Windows
```

**Linux**:
```bash
# Usually works without issues
mamba install -c conda-forge openff-toolkit ambertools
```

#### 5. **Mixed Installation Issues**

**Problem**: Mixing pip and conda installations causes conflicts

**Solutions**:
```bash
# Clean approach: Use conda/mamba for scientific packages
mamba install -c conda-forge openff-toolkit ambertools mdtraj
pip install -e .  # Only install PRISM via pip

# If you must use pip for everything:
pip install openff-toolkit openff-interchange mdtraj
pip install -e .[all]
```

#### 6. **Verification Commands**

After installation, verify everything works:
```bash
# Check PRISM installation
python -c "import prism as pm; print(pm.check_dependencies())"

# Test basic functionality
python -c "
import prism as pm
deps = pm.check_dependencies()
print('✓ GROMACS:', deps.get('gromacs', False))
print('✓ OpenFF:', deps.get('openff', False))
print('✓ AmberTools:', deps.get('antechamber', False))
print('✓ MDTraj:', deps.get('mdtraj', False))
"

# Test force field availability
python -c "
from prism.builder import PRISMBuilder
print('Available force fields:')
builder = PRISMBuilder()
print(builder.list_ligand_forcefields())
"
```

### Quick Fix Commands

For the most common issues, try these in order:

```bash
# 1. Quick fix for conda channel issues
mamba install -c conda-forge openff-toolkit ambertools

# 2. If mamba not available, install it first
conda install -c conda-forge mamba
mamba install -c conda-forge openff-toolkit ambertools

# 3. Nuclear option: fresh environment
mamba create -n prism-fresh -c conda-forge python=3.9 openff-toolkit ambertools mdtraj
conda activate prism-fresh
pip install -e .
```

4. **Memory errors**: Large systems may require more RAM during parameterization

### Debug Commands
```bash
# Check dependencies
python -c "import prism as pm; print(pm.check_dependencies())"

# List available force fields
python prism/builder.py --list-forcefields

# Test with minimal example
python prism/builder.py test/4xb4/4xb4_protein.pdb test/4xb4/4xb4_ligand.mol2 -o test_output
```

## Development Notes

- **Force field detection**: Uses GROMACS environment auto-detection
- **Dynamic ion addition**: Automatically detects water groups instead of hardcoded values
- **Error handling**: Comprehensive error messages with troubleshooting guidance
- **Flexibility**: Supports both file-based and programmatic configuration, multiple argument formats
- **Modularity**: Clean separation between force field generators, system building, and simulation
- **GPU support**: Optimized GROMACS commands for GPU acceleration in simulation scripts
- **Resume capability**: Simulation scripts check for existing files and can resume from checkpoints

## Entry Points

1. **Command line**: `python prism/builder.py` or `prism` (if installed)
2. **Python API**: `import prism as pm; pm.system()` or `pm.build_system()`
3. **Direct builder**: `from prism.builder import PRISMBuilder`
4. **Simulation**: `from prism.sim import model`

## Testing and Examples

- **Test data**: `test/4xb4/` contains complete protein-ligand example
- **Configuration examples**: `test/4xb4/prism_config.yaml` shows runtime configuration
- **MDP templates**: Generated dynamically based on configuration
- **Simulation scripts**: GPU-optimized examples in `prism/sim/scripts/`