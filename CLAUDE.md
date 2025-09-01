# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PRISM (Protein Receptor Interaction Simulation Modeler) is a comprehensive Python tool for building protein-ligand systems for molecular dynamics simulations in GROMACS. It supports multiple force fields including GAFF (via AmberTools) and OpenFF (Open Force Field).

## Development Commands

### Installation and Setup
```bash
# Install in development mode
pip install -e .

# Install with specific force field support
pip install -e .[gaff]          # GAFF support only
pip install -e .[openff]        # OpenFF support only  
pip install -e .[all]           # All force fields

# Basic installation
pip install -r requirements.txt
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

3. **Force field errors**: Ensure dependencies are installed:
   - **GAFF**: Requires AmberTools and ACPYPE
   - **OpenFF**: Requires openff-toolkit and openff-interchange

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