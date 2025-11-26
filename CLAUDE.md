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

## FEbuilder Ligand Conversion

- Environment: `conda activate /home/gxf1212/.conda/envs/prism` then `pip install -e .` (or `pip install .`) so the `prism` package and converters are importable system-wide.
- From repo root run `python test/modeling/FEbuilder_FF/AMBER/generate_charmm_ff.py` to batch convert every `*.mol2` in `test/modeling/FEbuilder_FF/AMBER/` (38/42 by default).
- The script calls `GAFFForceFieldGenerator` to populate `test/modeling/FEbuilder_FF/AMBER/_build/<name>/LIG.amb2gmx` and then feeds that into `AmberToCharmmConverter`, yielding `<name>.rtf/.prm/.pdb` plus `<name>_3D.mol2` inside the AMBER directory by default (override with `--charmm-dir` if needed).
- Regenerate individual molecules or switch residue name via `python .../generate_charmm_ff.py --resname XXX path/to/file.mol2`; add `--overwrite` if cached GAFF files should be replaced.
