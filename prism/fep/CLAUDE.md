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

### Gray Atoms in HTML Visualization
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
- **GAFF/CGenFF**: Position-specific atom types (ca, c3, ha, hc) - full matching
- **OpenFF**: Generic types (output_0, output_1) - distance + element + charge only
- **OPLS-AA**: Sequential types (opls_800, opls_801) - distance + element + charge only

**Implementation**: `DistanceAtomMapper._should_ignore_atom_type()` auto-detects force field type

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
