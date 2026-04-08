# PRISM-FEP Module Instructions

This file provides guidance for Claude Code when working with FEP (Free Energy Perturbation) functionality in PRISM.

## Module Overview

PRISM-FEP implements free energy perturbation calculations for relative binding free energies between similar ligands.

**Key Features**:
- Automated atom mapping via DistanceAtomMapper (with configurable distance/charge cutoffs)
- Hybrid topology generation (A-state + B-state parameters)
- Lambda window setup for alchemical transformation (standard + repex modes)
- Multi-force field support: GAFF, CGenFF, OpenFF, OPLS-AA
- HTML visualization with bond order rendering and quality checks
- MOL2 coordinate file support with unified PDB/MOL2/GRO handling
- Repeat-exchange (lambda replica exchange) script generation

## Architecture

**Core Modules**:
- `prism/fep/core/` - Atom mapping and hybrid topology data structures
- `prism/fep/modeling/` - FEP system building and workflow orchestration
- `prism/fep/gromacs/` - GROMACS integration and MDP templates
- `prism/fep/visualize/` - HTML/PNG visualization with quality checks
- `prism/fep/io.py` - Unified PDB/MOL2/GRO coordinate handling
- `prism/fep/config.py` - YAML configuration management

**Entry Points**:
- CLI: `prism protein.pdb ref.mol2 --fep --mutant mut.mol2`
- Python API: `from prism.fep import FEPScaffoldBuilder`
- End-to-end: `python -m prism.fep.modeling.e2e`

## Standard Naming Conventions

**CRITICAL**: Follow PRISM-wide naming conventions (see `NAMING_CONVENTIONS.md`)

Key rules for FEP:
1. **System directories**: `GMX_PROLIG_FEP/`, `GMX_PROLIG_MD/`
2. **Force field directories**: Use `ffgen.get_output_dir_name()`
3. **Test directories**: `gaff_test/`, `test_42_38/`, NOT `test_output_final/`
4. **Never hardcode paths**: Always use dynamic methods
5. **Never create nested duplicate system directories**: prohibit `.../GMX_PROLIG_FEP/GMX_PROLIG_FEP/`
6. **Use simple case naming for force-field combinations**: `<protein_ff>-mut_<ligand_ff>` (example: `charmm36m-mut_mmff`), avoid ad-hoc suffixes
7. **Default FEP output naming**: when no explicit output directory is provided, default to `<protein_ff>-mut_<ligand_ff>` and place the scaffold at `<case_dir>/GMX_PROLIG_FEP/`
8. **42-38 cleanup convention**: move legacy/non-canonical outputs into `tests/gxf/FEP/unit_test/42-38/Archive/`

## Common Issues

## Critical Development Constraints

### Atom Name Mapping and Coordinate Handling (CRITICAL)

**Problem Background**:
- Force field conversion changes atom names (e.g., OPLS-AA: CZ→C, OpenFF: generic types)
- Mapping phase and visualization phase may use different coordinate sources
- This causes atom correspondence errors and incorrect visualizations

**Mandatory Requirements**:
- ✅ **Track atom name changes** during force field conversion
- ✅ **Coordinate source priority**: MOL2 > PDB > GRO (GRO as fallback only)
- ✅ **Sanity checks for common atoms**: All common atoms must have same element type
- ✅ **Chemical validation**: Carbon should not appear at chain ends, hydrogen should not connect to multiple bonds

**Forbidden Practices**:
- ❌ Using `mapping_state_a.pdb`/`mapping_state_b.pdb` as coordinate sources
- ❌ Using mapping files for coordinate lookup (always use original ligand GRO files)
- ❌ Relying on atom name-only matching without considering force field type

**Code Locations**:
- `prism/fep/io.py::_load_structure_coordinates()` - Unified coordinate handling
- `prism/fep/modeling/hybrid_package.py::_resolve_hybrid_atom_coordinates()` - Coordinate resolution
- `prism/builder/workflow_fep.py` (lines 254-272) - Mapping file generation with hybrid ITP atoms

**Verification**:
```python
# After hybrid GRO generation, verify:
unique_coords = len(set(all_coordinates))
total_atoms = len(all_coordinates)
zero_coords = sum(1 for c in all_coordinates if c == (0.0, 0.0, 0.0))
assert unique_coords == total_atoms, "All atoms must have unique coordinates"
assert zero_coords == 0, "No placeholder (0,0,0) coordinates allowed"
```

### Testing Parameter Passing Standards (CRITICAL)

**Forbidden Practices**:
- ❌ **Hardcoding parameters in test files** (e.g., `ligand_a_name=`, `dist_cutoff=`)
- ❌ **Bypassing YAML configuration** to pass parameters directly to DistanceAtomMapper
- ❌ **Manual parameter specification** that doesn't reflect real user workflow

**Correct Approach**:
- ✅ All test files must read configuration through `case.yaml`
- ✅ All FEP parameters must be passed via YAML files to simulate real user scenarios
- ✅ Configuration priority: CLI args > fep.yaml > config.yaml > code defaults

**Code Example**:
```python
# ❌ WRONG: Manual specification
mapper = DistanceAtomMapper(dist_cutoff=0.6, charge_cutoff=0.05)

# ✅ CORRECT: Read from YAML
config = FEPConfig(config_file="case.yaml")
mapper = DistanceAtomMapper(**config.mapping_params)
```

**Enforcement**:
- Test files should only contain `case.yaml` path and expected results
- All mapping parameters (dist_cutoff, charge_cutoff, etc.) must come from YAML
- This ensures tests reflect actual user workflow

### Script Generation Constraints (CRITICAL)

**Forbidden Practices**:
- ❌ **Generate excessive redundant scripts** (standard/repex/local/parallel, etc.)
- ❌ **Nest scripts inside leg directories** (e.g., `bound/scripts/run.sh`)
- ❌ **Hardcode temporary values** (e.g., "100 steps") into production code

**Correct Approach**:
- ✅ Generate at most 2-3 script variants (control via parameters)
- ✅ Place bound/unbound scripts at sibling level (not nested)
- ✅ Control parallel vs local via parameters, not separate scripts
- ✅ Implement short tests via configuration, not hardcoded values

**Script Structure**:
```
GMX_PROLIG_FEP/
├── bound/
│   ├── run_fep.sh          # Universal script (parameter-controlled)
│   └── window_00/ to window_N/
└── unbound/
    ├── run_fep.sh          # Universal script (parameter-controlled)
    └── window_00/ to window_N/
```

**Runtime Flexibility**:
- System building phase: Control total replica count
- Runtime phase: Support `--replica-range 1-2`, `--replica 1`, etc.
- Bound/unbound should use identical script structure
- Prefer using generated scripts over manual commands

### Visualization Quality Standards

**Mandatory Verification Metrics**:
- ✅ **Zero gray atoms**: `grep '"classification": "unknown"' mapping.html` returns 0
- ✅ **Total charge ≈ 0** for neutral molecules
- ✅ **Mapping coverage**: total = common + transformed + surrounding
- ✅ **Element consistency**: All common atoms have same element in both states
- ✅ **Chemical sanity**: No carbons at chain ends, no hydrogens with multiple bonds

**Forbidden Practices**:
- ❌ Purple backgrounds or flashy styling
- ❌ Deviating from mapping module style
- ❌ Debugging via images (must use log analysis instead)

**Frontend Display Requirements**:
- ✅ YAML parameters shown in collapsible section (non-default values highlighted)
- ✅ Important parameters have tooltips
- ✅ Clear title field sources (e.g., `Ligand A: 39`, not `Ligand A: Ligand 39`)
- ✅ FEP-specific parameters shown (force field type, cutoffs)
- ✅ General modeling parameters omitted (salt concentration, etc.)

### Force Field-Specific Constraints

**OpenFF/OPLS-AA Limitations**:
- ❌ **Cannot rely on atom type** for position distinction (use generic types)
- ✅ **Mapping must ignore atom type**, only use distance + element + charge
- ✅ Other force fields (GAFF/CGenFF) can use atom type normally

**RTF/CGenFF Requirements**:
- ✅ **Must integrate original CGenFF files** (toppar_c36_feb26)
- ✅ GROMACS official charmm36.ff lacks small molecule parameters
- ✅ RTF priority: MATCH > SwissParam > CGenFF website

**Protonation State Handling**:
- ✅ **Auto-map residue names based on RTP file** (not hardcoded Amber logic)
- ✅ Support alias mapping (HID→HSD, HIE→HSE, HIP→HSP)
- ✅ Warning on unmapped residues, preserve original name
- ✅ All protonation methods use same utils (PROPKA, etc.)

**Code Example**:
```python
# ❌ WRONG: Hardcoded Amber mapping
if residue_name == "HID":
    residue_name = "HSD"  # Amber-specific

# ✅ CORRECT: RTP-based auto-mapping
rtp_file = force_field_path / "aminoacids.rtp"
if residue_name in rtp_aliases:
    residue_name = rtp_aliases[residue_name]
elif residue_name not in available_rtp_names:
    logger.warning(f"Cannot map {residue_name}, keeping original")
```

### Engineering Standards (CRITICAL)

**Forbidden Practices**:
- ❌ **Delete existing HTML and test artifacts** (especially `*_test_output/*.html`)
- ❌ **Hardcode selection/resname** (e.g., `if 'HA' in u.select_atoms('all').resnames`)
- ❌ **Use Chinese punctuation** in code (must use English only)
- ❌ **Read mapping_state_a.pdb as coordinate source** (use original ligand GRO)

**Correct Practices**:
- ✅ All selections passed via config
- ✅ Coordinate source priority: original ligand GRO > mapping files
- ✅ Use hybrid.gro (not mapping files) as coordinate source
- ✅ Preserve all test artifacts as verification evidence

**Code Quality**:
- ✅ No hardcoded magic numbers (use configuration)
- ✅ Clear separation of concerns (no mixed responsibilities)
- ✅ Consistent naming conventions across modules
- ✅ Proper error handling and logging

### Hybrid Topology Coordinate Assignment

**Critical Issue**: Incorrect coordinate assignment causes EM crashes (6+ nm bond distances)

**Root Cause**:
- Hybrid topology atom names: HA, OA, C, H1, H2, C1, N, ... (specific)
- Mapping file atom names: H, O, C, H, H, C, N, ... (generic)
- Coordinate lookup fails → placeholder (0,0,0) → EM crash

**Solution Implemented**:
1. **workflow_fep.py (lines 254-272)**: Read atoms from hybrid ITP for mapping files
2. **hybrid_package.py (lines 202-208)**: Always use original ligand GRO as coordinate source
3. **hybrid_package.py (lines 975-982)**: Only support pure element names (HA, OA, CA, N, F)

**Verification After Build**:
```python
# Check hybrid.gro has all unique coordinates
coords = []
with open("hybrid.gro") as f:
    for line in f.readlines()[2:-1]:
        x, y, z = float(line[20:28]), float(line[28:36]), float(line[36:44])
        coords.append((x, y, z))

assert len(set(coords)) == len(coords), "All coordinates must be unique"
assert sum(1 for c in coords if c == (0.0, 0.0, 0.0)) == 0, "No placeholder coordinates"
```

### Documentation and Code Consistency

**Required Consistency**:
- ✅ Code behavior must match documentation
- ✅ Parameter names consistent across docs and code
- ✅ Examples reflect actual usage patterns
- ❌ No undocumented features or parameters

**Review Checklist**:
- [ ] All CLAUDE.md files reviewed together
- [ ] Code examples tested and verified
- [ ] Parameter documentation complete
- [ ] No contradictions between different docs

## Legacy Issues (Historical Reference)
**Symptom**: Some atoms appear gray with `classification: "unknown"`

**Solution**: Usually caused by RDKit `Chem.AddHs()` adding implicit hydrogens. For CHARMM-GUI systems, use PDB files directly (not MOL2) and skip SMILES round-trip.

**Code location**: `prism/fep/visualize/molecule.py::pdb_to_mol()`

### PNG Generation Fails with RDKit Error
**Symptom**: `ValueError: Depict error: Substructure match with reference not found`

**Cause**: Multi-point mutations (e.g., 39-8 system) have too different structures for RDKit alignment

**Solution**: HTML visualization works fine. PNG generates with warning using default coordinates.

### Mapping Results vs FEbuilder Discrepancies
**Issue**: DistanceAtomMapper shows more transformed atoms than FEbuilder

**Explanation**: DistanceAtomMapper uses systematic, conservative matching strategy:
- Multi-point mutations (39-8): May have additional transformed atoms due to real structural differences
- Single-point mutations (42-38): Should match FEbuilder results

**Verification**: Check atom type, charge, and position differences. Extra transformed atoms are usually legitimate.

### Multi-Force Field Support
**Supported Force Fields**:
- **GAFF/GAFF2**: Position-specific atom types (ca, c3, ha, hc) - full matching
- **CGenFF**: CHARMM General Force Field - full matching
- **RTF**: CHARMM Residue Topology File format - full matching
- **OpenFF**: Generic types (output_0, output_1) - distance + element + charge only
- **OPLS-AA**: Sequential types (opls_800, opls_801) - distance + element + charge only

**Implementation**: `DistanceAtomMapper._should_ignore_atom_type()` auto-detects force field type

### RTF (CHARMM) Force Field Support
**Status**: ✅ Fully supported (as of 2026-04-07)

**Key fixes applied**:
1. **Hybrid topology bond merging** (`prism/fep/gromacs/itp_builder.py`):
   - Fixed over-aggressive zeroing of bonded terms when atom types change between states
   - Now only zeros parameters when a term is missing in one state (not when both states have valid parameters)
   - Prevents artificial bond breaking that caused molecular collapse during EM

2. **Coordinate unit conversion** (`prism/fep/modeling/hybrid_package.py`):
   - Fixed Å→nm conversion for PDB/MOL2 coordinates (divide by 10.0)
   - GRO coordinates already in nm, no conversion needed

**Test system**: `tests/gxf/FEP/unit_test/p38-19-24/` (p38α MAPK, ligands 19→24)
- Protein FF: charmm36-jul2022
- Ligand FF: RTF
- Status: ✅ EM + NVT passing for both bound/unbound legs

### MOL2 Coordinate Support
**Feature**: Unified PDB/MOL2/GRO coordinate handling via `_load_structure_coordinates()`

**Benefits**:
- Use original ligand coordinates (MOL2/PDB) instead of GRO for better accuracy
- Automatic element symbol extraction from MOL2 atom names
- Fallback to GRO if MOL2 unavailable

**Code location**: `prism/fep/io.py::_load_structure_coordinates()`

## Testing

### Recommended Test Systems
- **42-38**: Single-point mutation ✅ Recommended
- **39-8**: Multi-point mutations (may trigger PNG alignment issues)
- **oMeEtPh-EtPh**: Terminal ethyl methyl transformation

### Test Checklist
- [ ] All atoms classified (no gray atoms in HTML)
- [ ] Total charges reasonable (≈0 for neutral molecules)
- [ ] Mapping coverage: `total = common + transformed + surrounding`
- [ ] PNG generates (or shows appropriate warning)
- [ ] HTML shows warning banner if problems detected
- [ ] Bond orders rendered correctly for aromatic rings
- [ ] Repeat-exchange scripts generated if repex mode enabled

### Repeat-Exchange Mode
**Feature**: Lambda replica exchange via `gmx_mpi -multidir`

**Configuration**:
```yaml
fep:
  mode: repex  # or 'standard'
```

**Generated scripts**:
- `run_prod_standard.sh` - Sequential lambda windows
- `run_prod_repex.sh` - Parallel replica exchange

**Code location**: `prism/fep/modeling/script_writer.py::write_production_scripts()`

### Bond Order Visualization
**Feature**: Professional bond order rendering for aromatic rings and conjugated systems

**Implementation**:
- Uses MOL2 file as bond order template
- RDKit MCS matching to assign bond orders
- Falls back to RDKit sanitization if MOL2 unavailable

**Code location**: `prism/fep/visualize/molecule.py::assign_bond_orders_from_mol2()`

### Test Scripts
```bash
# Unit tests
python tests/gxf/FEP/unit_test/test_42_38.py
python tests/gxf/FEP/unit_test/test_39_8.py

# End-to-end workflow
python -m prism.fep.modeling.e2e --help
```

## Visualization Quality Control

**Always verify HTML output after generation**:
1. Check for gray atoms (should be 0)
2. Verify mapping statistics make sense
3. Confirm total charges ≈ 0
4. Test button styling (FEP Classification, Element)
5. Verify atom labels show real names (not "Atom15")

**Verification commands**:
```bash
# Check for gray atoms
grep "rgb(200, 200, 200)" output/mapping.html  # Should return 0

# Check for unknown classification
grep '"classification": "unknown"' output/mapping.html  # Should return 0
```

## CHARMM-GUI Integration

**Special handling required**:
- PDB files already contain all explicit hydrogens
- Do NOT use RDKit `Chem.AddHs()` - creates extra atoms
- Skip SMILES round-trip to preserve atom names
- Atom types may differ from FEbuilder expectations

**Code path**: `prism/fep/visualize/molecule.py::pdb_to_mol()`

## Development Workflow

1. **Modify code** in `prism/fep/`
2. **Test** with 42-8 system first (simplest)
3. **Verify** with 39-8 system (complex, multi-point)
4. **Generate HTML** and check visualization quality
5. **Commit** with descriptive message: `[fep] brief description`

### FEP测试资源管理

**CPU/GPU资源要求**:
- **所有FEP测试必须使用14个CPU核心** (`-ntomp 14`)
- 56核心系统最多支持4个并行FEP任务 (4 × 14 = 56核心)
- **除非有SQM计算**，否则必须保证每个gmx任务使用14核心
- 使用 `find . -name "run_single_em_nvt.sh" -exec sed -i 's/-ntomp [0-9]*/-ntomp 14/g' {} \;` 更新脚本

**GPU分配策略**:
- 每个FEP测试分配独立GPU (GPU 0-5)
- 使用 `CUDA_VISIBLE_DEVICES=N` 指定GPU编号
- 避免GPU过载，确保每个测试有足够显存

**并行任务监控**:
- 最多同时运行4个FEP测试
- 监控命令: `ps aux | grep "gmx mdrun" | grep -v grep | wc -l`
- 如果超过4个并行任务，必须停止最不重要的一个
- 检查GPU状态: `nvidia-smi | grep "gmx"`

## Git Workflow

- Branch: `gxf` (FEP development branch)
- Main: `main` (integration target)
- Commit format: `[fep] description; details`
- Always verify HTML quality before committing

## Related Documentation

- `NAMING_CONVENTIONS.md` - Standardized naming conventions
- `README.md` - Module overview and API reference
- `CLAUDE.md` (root) - General PRISM instructions
- `.claude/skills/fep-mapping.md` - Atom mapping algorithm verification guide
- `.claude/skills/fep-visualization.md` - HTML visualization verification guide

## Module-Specific Documentation

- `prism/fep/core/CLAUDE.md` - Atom mapping and hybrid topology details
- `prism/fep/visualize/CLAUDE.md` - Visualization implementation and troubleshooting
- `prism/fep/modeling/CLAUDE.md` - FEP system building workflow
- `prism/fep/gromacs/CLAUDE.md` - GROMACS integration and MDP templates
