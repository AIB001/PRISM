# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PRISM (Protein Receptor Interaction Simulation Modeler) is a comprehensive Python tool for building protein-ligand systems for molecular dynamics simulations in GROMACS. It supports multiple force fields including GAFF (via AmberTools) and OpenFF (Open Force Field).

**Documentation Repository**: `/home/gxf1212/data/work/PRISM-Tutorial` - Official documentation and tutorials for PRISM

## Quick Start

### Installation
```bash
pip install -e .                    # Development mode
pip install -e .[all]               # All force fields
```

### Basic Usage
```bash
# CLI
python prism/builder.py protein.pdb ligand.mol2 -o output_dir

# Python API
import prism as pm
system = pm.system("protein.pdb", "ligand.mol2")
output_dir = system.build()
```

## Architecture Overview

**Core Modules**:
- `prism.builder` - Main entry point and workflow orchestration
- `prism.core` - High-level Python API with `PRISMSystem` class
- `prism.forcefield/` - Force field generation (GAFF, OpenFF)
- `prism.utils/` - GROMACS environment, configuration, system assembly
- `prism.sim/` - Simulation execution (GROMACS, OpenMM)
- `prism.pmf/` - PMF calculation via steered MD and umbrella sampling
- `prism.analysis/` - Trajectory analysis tools
  - `prism.analysis.calc/` - Calculation modules (RMSD, contacts, clustering, etc.)
  - `prism.analysis.plots/` - Hierarchical plotting modules organized by analysis type
  - `prism.analysis.contact/` - HTML visualization generation

**Configuration**: CLI args > config file > defaults (`prism/configs/default_config.yaml`)

**Output Structure**: `output_dir/GMX_PROLIG_MD/` (GROMACS files), `LIG.amb2gmx/` or `LIG.openff2gmx/` (ligand files), `mdps/` (MD parameters)

## Key Development Principles

- **Force field detection**: Auto-detect GROMACS environment
- **Dynamic ion addition**: Auto-detect water groups (not hardcoded)
- **Error handling**: Comprehensive error messages with troubleshooting guidance
- **Modularity**: Clean separation between force field generators, system building, and simulation
- **GPU support**: Optimized GROMACS commands for GPU acceleration
- **Resume capability**: Check for existing files and resume from checkpoints

## Common Issues

1. **Argument parsing**: Use positional args (protein ligand options) or flags (--protein-file --ligand-file)
2. **Ion addition failures**: Ensure system contains water molecules, check GROMACS can find solvent groups
3. **Force field errors**: Install required dependencies (GAFF needs AmberTools/ACPYPE, OpenFF needs openff-toolkit)
4. **Memory errors**: Large systems may require more RAM during parameterization

## Entry Points

- CLI: `python prism/builder.py` or `prism` (if installed)
- Python API: `import prism as pm; pm.system()` or `pm.build_system()`
- Direct builder: `from prism.builder import PRISMBuilder`
- Simulation: `from prism.sim import model`

## Testing

Test data in `test/4xb4/` contains complete protein-ligand example. MDP templates generated dynamically based on configuration.

Git commit message format: `[module] description;xxx`. Don't be too long!

## FEbuilder Ligand Conversion (AMBER → CHARMM)

### Quick Start
```bash
# Environment setup
conda activate /home/gxf1212/.conda/envs/prism
pip install -e .  # Enables prism package imports system-wide

# Batch convert all ligands to CHARMM/NAMD format
python test/modeling/FEbuilder_FF/AMBER/generate_charmm_ff.py
```

### Process & Implementation

**Conversion Pipeline**: `mol2` → `GAFF` (via antechamber/tleap) → `CHARMM RTF/PRM` (via `AmberToCharmmConverter`)

1. **GAFFForceFieldGenerator** populates `test/modeling/FEbuilder_FF/AMBER/_build/<name>/LIG.amb2gmx/` with GROMACS-format files
2. **AmberToCharmmConverter** converts GAFF parameters to CHARMM topology (RTF) and parameters (PRM) files
3. Output: `<name>.rtf`, `<name>.prm`, `<name>.pdb`, `<name>_3D.mol2` in the AMBER directory

### Important: Atom Type Naming for FEP

When converting multiple ligands for FEP (free energy perturbation) calculations, **each molecule gets a unique residue-name prefix** to avoid VMD psfgen type conflicts:
- Molecule `38.mol2` → types like `38_ca`, `38_c3`, `38_sy` (residue name "38")
- Molecule `42.mol2` → types like `42_ca`, `42_c3`, `42_sy` (residue name "42")

This prevents "duplicate type key" warnings when loading multiple force field topologies simultaneously in VMD.

**Customization**:
- Single molecule with custom residue name: `python generate_charmm_ff.py --resname XXX path/to/ligand.mol2`
- Force regeneration: add `--overwrite` flag (requires AmberTools installed)
