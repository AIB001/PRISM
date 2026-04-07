# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PRISM (Protein Receptor Interaction Simulation Modeler) is a comprehensive Python tool for building protein-ligand systems for molecular dynamics simulations in GROMACS. It supports multiple force fields including GAFF (via AmberTools) and OpenFF (Open Force Field).

**Documentation Repository**: `/home/gxf1212/data/work/PRISM-Tutorial` - Official documentation and tutorials for PRISM

## Code Search and Navigation

**📝 LSP Tools Preferred**: For code navigation, use LSP tools over grep when available:
- **LSP tools** (`mcp__cclsp__` series): Semantic understanding, direct symbol location
  - `mcp__cclsp__find_definition` - Find symbol definitions accurately
  - `mcp__cclsp__find_references` - Find all symbol references
  - `mcp__cclsp__find_workspace_symbols` - Search workspace symbols
  - Benefits: More accurate, skips comments/strings, reduces token usage
- **Grep**: Text search requiring human interpretation
- **Availability**: LSP tools available in Claude Code tool context only
- **Fallback**: Use `Grep` tool when LSP tools are unavailable

## Quick Start

### Installation
```bash
pip install -e .                    # Development mode
pip install -e .[all]               # All force fields
```

### Basic Usage
```bash
# CLI - Single ligand
python prism/builder.py protein.pdb ligand.mol2 -o output_dir

# CLI - Multiple ligands
prism -pf protein.pdb -lf ligand1.mol2 -lf ligand2.mol2 -o output_dir -lff gaff -ff amber14sb

# CLI - CGenFF (requires one -ffp per ligand)
prism -pf protein.pdb -lf ligand1.mol2 -lf ligand2.mol2 -o output_dir \
  --ligand-forcefield cgenff -ffp /path/to/cgenff1 -ffp /path/to/cgenff2

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

**Output Structure**:
- Single ligand: `output_dir/LIG.amb2gmx/` or `LIG.openff2gmx/` (backward compatible)
- Multi-ligand: `output_dir/Ligand_Forcefield/LIG.amb2gmx_1/`, `LIG.amb2gmx_2/`, etc.
- GROMACS files: `output_dir/GMX_PROLIG_MD/`
- MD parameters: `output_dir/mdps/`

## Key Development Principles

- **Force field detection**: Auto-detect GROMACS environment
- **Dynamic ion addition**: Auto-detect water groups (not hardcoded)
- **Error handling**: Comprehensive error messages with troubleshooting guidance
- **Modularity**: Clean separation between force field generators, system building, and simulation
- **GPU support**: Optimized GROMACS commands for GPU acceleration
- **Resume capability**: Check for existing files and resume from checkpoints
- **Naming consistency**: Follow standardized naming conventions (see `NAMING_CONVENTIONS.md`)

## Standard Naming Conventions

**CRITICAL**: All new code MUST follow the standardized naming conventions defined in `NAMING_CONVENTIONS.md`.

### Key Rules

1. **System Directories**: Always use `GMX_PROLIG_*` prefix
   - ✅ `GMX_PROLIG_MD/`, `GMX_PROLIG_FEP/`, `GMX_PROLIG_PMF/`
   - ❌ `FEP_SYSTEM/`, `output/`, `my_system/`

2. **Force Field Directories**: Use `ffgen.get_output_dir_name()` method
   - ✅ `ff_dir = output_dir / ffgen.get_output_dir_name()`
   - ❌ `ff_dir = output_dir / "LIG.amb2gmx"`

3. **Test Directories**: Use descriptive, forcefield-specific names
   - ✅ `gaff_test/`, `openff_test/`, `test_42_38/`
   - ❌ `gaff2_e2e_test/`, `ligand_42/`, `test_output_final/`

4. **Never hardcode paths**:
   - ✅ `prism_dir = ligand_output / ffgen.get_output_dir_name()`
   - ❌ `prism_dir = "output/LIG.amb2gmx"`

5. **No nested duplicate system directories**:
   - ✅ `tests/gxf/FEP/unit_test/oMeEtPh-EtPh/charmm36m-mut_mmff/GMX_PROLIG_FEP/`
   - ❌ `tests/gxf/FEP/unit_test/oMeEtPh-EtPh/charmm36m-mut_mmff_pkgfix2/GMX_PROLIG_FEP/GMX_PROLIG_FEP/`

6. **Case directory naming must stay simple and readable**:
   - ✅ `<protein_ff>-mut_<ligand_ff>` (example: `charmm36m-mut_mmff`)
   - ❌ `*_pkgfix*`, `*_final*`, `*_new*`, or other ad-hoc suffixes

### Examples

```python
# ✅ Correct - Using standardized methods
ffgen = GAFFForceFieldGenerator(ligand_path="lig.mol2", output_dir="gaff_test")
ffgen.run()
prism_dir = Path(ffgen.output_dir) / ffgen.get_output_dir_name()

# ✅ Correct - Standard system naming
builder = FEPScaffoldBuilder(output_dir="GMX_PROLIG_FEP")

# ❌ Wrong - Hardcoded path
prism_dir = "gaff_test/LIG.amb2gmx"

# ❌ Wrong - Non-standard naming
builder = FEPScaffoldBuilder(output_dir="FEP_SYSTEM")
```

**See `NAMING_CONVENTIONS.md` for complete specification.**

## Common Issues

1. **Argument parsing**: Use positional args (protein ligand options) or flags (--protein-file --ligand-file)
2. **Multi-ligand CLI**: For multiple ligands, use `-pf` and `-lf` flags (not positional) for best compatibility
3. **CGenFF multi-ligand**: Requires one `--forcefield-path` (or `-ffp`) per ligand
4. **Ion addition failures**: Ensure system contains water molecules, check GROMACS can find solvent groups
5. **Force field errors**: Install required dependencies (GAFF needs AmberTools/ACPYPE, OpenFF needs openff-toolkit)
6. **Memory errors**: Large systems may require more RAM during parameterization

## Module-Specific Documentation

- **FEP Module**: See `prism/fep/CLAUDE.md` for FEP-specific guidelines and troubleshooting
- **PMF Module**: See `prism/pmf/CLAUDE.md` for PMF-specific guidelines

## Entry Points

- CLI: `python prism/builder.py` or `prism` (if installed)
- Python API: `import prism as pm; pm.system()` or `pm.build_system()`
- Direct builder: `from prism.builder import PRISMBuilder`
- Simulation: `from prism.sim import model`

## Testing

Test data in `test/4xb4/` contains complete protein-ligand example. MDP templates generated dynamically based on configuration.

### FEP System Testing

**Comprehensive FEP validation**: Use `tests/gxf/FEP/unit_test/test_run_fep.py` for complete FEP system testing:
- System building and directory structure validation
- Topology and coordinate file verification
- Grompp validation for both bound/unbound legs
- Mapping HTML quality checks
- Atom classification validation

**Usage** (works for all test systems):
```bash
# 42-38 system (HIF-2α)
cd tests/gxf/FEP/unit_test/42-38
python ../test_run_fep.py --forcefield amber14sb --ligand-forcefield gaff2

# 25-36 system (HIF-2α)
cd tests/gxf/FEP/unit_test/25-36
python ../test_run_fep.py --forcefield amber14sb --ligand-forcefield opls

# oMeEtPh-EtPh system (T4 lysozyme L99A)
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh
python ../test_run_fep.py --forcefield amber14sb --ligand-forcefield openff

# p38-19-24 system (p38α MAPK)
cd tests/gxf/FEP/unit_test/p38-19-24
python ../test_run_fep.py --forcefield charmm36-jul2022 --ligand-forcefield rtf
```

Git commit message format: `[module] description;xxx`. Don't be too long!

## Development Notes

- **LSP vs grep**: Prefer LSP tools (mcp__cclsp*) for symbol finding - more accurate, skips comments/strings
- **Commit messages**: Should reflect all recent changes in the commit
- **File creation**: Avoid creating files in root directory

### CPU/GPU资源管理

**CPU核心数要求**:
- **所有GROMACS任务必须使用14个CPU核心** (`-ntomp 14`)
- 56核心系统 ÷ 14核心/任务 = 最多4个并行任务
- **除非有SQM计算**，否则必须保证每个gmx任务使用14核心
- 更新脚本时使用: `find . -name "run_single_em_nvt.sh" -exec sed -i 's/-ntomp [0-9]*/-ntomp 14/g' {} \;`

**GPU分配**:
- 系统有8个GPU (GPU 0-7)，通常使用GPU 0-5
- 每个测试分配独立GPU，避免GPU过载
- 使用 `CUDA_VISIBLE_DEVICES=N` 指定GPU

**并行任务管理**:
- 最多同时运行4个FEP测试 (4 × 14 = 56核心)
- 如果超过4个并行任务，必须停止最不重要的一个
- 监控命令: `ps aux | grep "gmx mdrun" | grep -v grep | wc -l`

不用每次push；TODO完成了就清掉

尽量复用已有的测试脚本。
