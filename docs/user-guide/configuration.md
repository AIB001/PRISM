# Configuration Guide

PRISM uses a flexible YAML-based configuration system that allows you to customize every aspect of your molecular dynamics simulations. This guide covers all configuration options and best practices.

## Configuration Overview

PRISM can be configured through:
1. **Command-line arguments** (highest priority)
2. **Configuration files** (YAML format)
3. **Default settings** (built-in defaults)

## Using Configuration Files

### Basic Usage

Create a custom configuration file and use it with PRISM:

```bash
# Use a custom configuration
prism protein.pdb ligand.mol2 -o output --config my_config.yaml

# Export the default configuration as a template
prism --export-config my_config.yaml
```

### Configuration File Structure

```yaml
# my_config.yaml
general:
  overwrite: false  # Whether to overwrite existing files
  
box:
  distance: 1.5     # Distance from protein to box edge (nm)
  shape: cubic      # Box shape: cubic, dodecahedron, octahedron
  center: true      # Center the protein in the box

simulation:
  temperature: 310  # Temperature in Kelvin
  pressure: 1.0     # Pressure in bar
  pH: 7.0          # pH for protonation states
  production_time_ns: 500  # Production run time
```

## Configuration Sections

### General Settings

Controls overall behavior of PRISM:

```yaml
general:
  overwrite: false  # Overwrite existing files
  # gmx_command is auto-detected but can be specified:
  # gmx_command: /usr/local/bin/gmx
```

### Box Settings

Define the simulation box parameters:

```yaml
box:
  distance: 1.5     # Minimum distance to box edge (nm)
  shape: cubic      # Options: cubic, dodecahedron, octahedron
  center: true      # Center protein in box
```

**Box Shape Recommendations:**
- `cubic`: Simple, good for most systems
- `dodecahedron`: ~29% less volume, more efficient
- `octahedron`: ~20% less volume than cubic

### Simulation Parameters

Core MD simulation settings:

```yaml
simulation:
  temperature: 310        # Body temperature (K)
  pressure: 1.0          # Atmospheric pressure (bar)
  pH: 7.0               # Physiological pH
  ligand_charge: 0      # Net charge of ligand
  production_time_ns: 500  # Production run duration
  dt: 0.002             # Time step (ps)
  equilibration_nvt_time_ps: 500  # NVT equilibration
  equilibration_npt_time_ps: 500  # NPT equilibration
```

**Temperature Guidelines:**
- 277 K: Near freezing (4°C)
- 298 K: Room temperature (25°C)
- 310 K: Body temperature (37°C)
- 323 K: Elevated temperature (50°C)

### Ion Settings

Control system neutralization and salt concentration:

```yaml
ions:
  neutral: true         # Neutralize system charge
  concentration: 0.15   # Salt concentration (M)
  positive_ion: NA      # Sodium ions
  negative_ion: CL      # Chloride ions
```

**Common Salt Concentrations:**
- 0.0 M: No added salt
- 0.15 M: Physiological conditions
- 0.5 M: High salt
- 1.0 M: Very high salt

### Energy Minimization

Parameters for energy minimization:

```yaml
energy_minimization:
  integrator: steep     # Steepest descent algorithm
  emtol: 200.0         # Convergence criterion (kJ/mol/nm)
  emstep: 0.01         # Initial step size
  nsteps: 10000        # Maximum steps
```

### Output Settings

Control output frequency and format:

```yaml
output:
  trajectory_interval_ps: 500   # Save coordinates every X ps
  energy_interval_ps: 10        # Save energies every X ps
  log_interval_ps: 10           # Log output frequency
  compressed_trajectory: true    # Use XTC compression
```

**Storage Considerations:**
- 500 ps interval: ~1 GB per μs for 50K atoms
- 100 ps interval: ~5 GB per μs for 50K atoms
- 10 ps interval: ~50 GB per μs for 50K atoms

### Advanced Settings

#### Electrostatics

```yaml
electrostatics:
  coulombtype: PME      # Particle Mesh Ewald
  rcoulomb: 1.0        # Coulomb cutoff (nm)
  pme_order: 4         # PME interpolation order
  fourierspacing: 0.16  # Grid spacing for FFT
```

#### Van der Waals

```yaml
vdw:
  rvdw: 1.0            # VdW cutoff (nm)
  dispcorr: EnerPres   # Long-range dispersion corrections
```

#### Temperature Coupling

```yaml
temperature_coupling:
  tcoupl: V-rescale    # Velocity rescaling thermostat
  tc_grps:             # Temperature coupling groups
    - Protein
    - Non-Protein
  tau_t:               # Coupling time constants (ps)
    - 0.1
    - 0.1
```

#### Pressure Coupling

```yaml
pressure_coupling:
  pcoupl: C-rescale    # C-rescale barostat (GROMACS 2021+)
  pcoupltype: isotropic  # Isotropic scaling
  tau_p: 1.0           # Coupling time constant (ps)
  compressibility: 4.5e-05  # Water compressibility
```

#### Constraints

```yaml
constraints:
  algorithm: lincs     # LINCS constraint algorithm
  type: h-bonds        # Constrain hydrogen bonds
  lincs_iter: 1        # LINCS iterations
  lincs_order: 4       # LINCS order
```

## Command-Line Overrides

Command-line arguments override configuration file settings:

```bash
# Override temperature
prism protein.pdb ligand.mol2 -o output --temperature 300

# Override multiple parameters
prism protein.pdb ligand.mol2 -o output \
  --forcefield amber14sb \
  --water tip4p \
  --temperature 298 \
  --production-ns 1000 \
  --box-distance 2.0
```

## Configuration Priority

Settings are applied in this order (highest to lowest priority):
1. Command-line arguments
2. Custom configuration file
3. Default configuration

Example:
```bash
# my_config.yaml has temperature: 310
# Command line overrides to 300
prism protein.pdb ligand.mol2 --config my_config.yaml --temperature 300
# Final temperature: 300 K
```

## Configuration Templates

### Membrane Protein Configuration

```yaml
# membrane_config.yaml
box:
  shape: cubic  # Better for membranes
  distance: 2.0  # Extra space for membrane

simulation:
  temperature: 310
  pressure: 1.0
  production_time_ns: 1000  # Longer for membrane systems

pressure_coupling:
  pcoupltype: semiisotropic  # For membranes
  tau_p: 5.0  # Slower coupling for membranes
```

### High-Throughput Screening

```yaml
# screening_config.yaml
general:
  overwrite: true  # Always overwrite

simulation:
  production_time_ns: 100  # Shorter runs
  equilibration_nvt_time_ps: 100
  equilibration_npt_time_ps: 100

output:
  trajectory_interval_ps: 1000  # Less frequent output
  compressed_trajectory: true
```

### Enhanced Sampling

```yaml
# enhanced_sampling.yaml
simulation:
  temperature: 323  # Higher temperature
  dt: 0.002
  production_time_ns: 2000  # Longer simulation

output:
  trajectory_interval_ps: 100  # More frequent sampling
```

## Python API Configuration

Using configuration in Python scripts:

```python
import prism

# Load from file
system = prism.PRISMSystem(
    "protein.pdb",
    "ligand.mol2",
    config="my_config.yaml"
)

# Or provide dictionary
config = {
    'simulation': {
        'temperature': 300,
        'production_time_ns': 1000
    },
    'box': {
        'distance': 2.0
    }
}

system = prism.PRISMSystem(
    "protein.pdb",
    "ligand.mol2",
    config=config
)

# Build system
system.build()
```

## Validation and Debugging

### Check Configuration

View the configuration used for a system:

```python
# The configuration is saved automatically
cat output/prism_config.yaml
```

### List Available Force Fields

```bash
prism --list-forcefields
```

### Validate Configuration

PRISM automatically validates configuration parameters:
- Temperature: 0-1000 K
- Pressure: 0.1-1000 bar
- Time step: 0.0001-0.004 ps
- Box distance: 0.5-10 nm

## Best Practices

1. **Start with defaults**: The default configuration works well for most systems
2. **Save your configuration**: Always save successful configurations for reproducibility
3. **Use appropriate time steps**: 
   - 2 fs (0.002 ps) with H-bond constraints
   - 1 fs (0.001 ps) without constraints
   - 4 fs (0.004 ps) with virtual sites
4. **Match temperature to experiment**: Use the same temperature as your experimental conditions
5. **Equilibrate properly**: Don't skip equilibration steps
6. **Monitor convergence**: Check that properties converge during equilibration

## Common Configuration Scenarios

### Drug-Target Binding

```yaml
simulation:
  temperature: 310  # Body temperature
  pH: 7.4          # Physiological pH
  production_time_ns: 500-1000
ions:
  concentration: 0.15  # Physiological
```

### Protein Folding

```yaml
simulation:
  temperature: 298  # or range for REMD
  production_time_ns: 1000-5000
box:
  distance: 2.0  # Extra space for unfolding
```

### Enzyme Catalysis

```yaml
simulation:
  temperature: 310
  pH: 7.0  # Or enzyme optimal pH
  dt: 0.001  # Smaller timestep for QM/MM
constraints:
  type: none  # For QM/MM regions
```

## Troubleshooting Configuration

### Common Issues

1. **"Invalid configuration key"**: Check for typos in YAML keys
2. **"Value out of range"**: Check parameter limits
3. **"File not found"**: Use absolute paths or check working directory
4. **"YAML parsing error"**: Check indentation and syntax

### Configuration Validation

```python
from prism.utils.config import ConfigurationManager

# Validate a configuration file
config_mgr = ConfigurationManager("my_config.yaml")
config_mgr.validate()  # Raises errors if invalid
```

## Next Steps

- Learn about [Input Files](input-files.md) preparation
- Understand [Force Fields](force-fields.md) selection
- Start [Building Systems](building-systems.md)