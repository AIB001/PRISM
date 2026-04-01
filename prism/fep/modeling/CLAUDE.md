# PRISM-FEP Modeling Module

This file provides guidance for working with FEP system building and workflow orchestration.

## Module Overview

The modeling module orchestrates the complete FEP system building workflow:
- **Scaffold building**: Creates directory structure for bound/unbound legs
- **Hybrid topology generation**: Combines A-state and B-state parameters
- **Script generation**: Produces run scripts for standard and repex modes
- **End-to-end workflow**: Automates from ligands to simulation-ready systems

## Key Classes

### FEPScaffoldBuilder
**Purpose**: Create FEP directory structure and coordinate system building

**Directory structure created**:
```
GMX_PROLIG_FEP/
├── bound/              # Protein-ligand complex leg
│   ├── build/          # EM, NVT, NPT equilibration
│   ├── window_00/      # Lambda window directories
│   ├── window_01/
│   ├── ...
│   └── localrun.sh     # GPU-accelerated run script
├── unbound/            # Ligand in water leg
│   └── (same structure as bound)
└── common/
    └── hybrid/         # Hybrid topology files
        ├── hybrid.itp
        ├── atomtypes_hybrid.itp
        ├── defaults_hybrid.itp
        └── hybrid.gro
```

**Usage**:
```python
from prism.fep.modeling.core import FEPScaffoldBuilder

builder = FEPScaffoldBuilder(
    protein_path="receptor.pdb",
    ligand_ref_path="ref.mol2",
    ligand_mut_path="mut.mol2",
    output_dir="GMX_PROLIG_FEP"
)

scaffold = builder.build_from_ligands()
```

### HybridPackageBuilder
**Purpose**: Generate hybrid topology files (ITP, GRO)

**Key responsibilities**:
- Merge reference and mutant ligand parameters
- Create dummy atom types for transformed atoms
- Extract [defaults] block to separate file
- Handle [atomtypes] section for all force fields

**Output files**:
- `hybrid.itp`: Molecule definition with A/B state parameters
- `atomtypes_hybrid.itp`: Atom type definitions
- `defaults_hybrid.itp`: Nonbonded parameters (optional, extracted from hybrid.itp)
- `hybrid.gro`: Coordinates with hybrid topology

### ScriptWriter
**Purpose**: Generate GROMACS run scripts

**Modes**:
- **standard**: Sequential lambda windows
- **repex**: Lambda replica exchange via `gmx_mpi -multidir`

**Generated scripts**:
- `run_prod_standard.sh`: Standard mode production
- `run_prod_repex.sh`: Replica exchange production
- `localrun.sh`: GPU-accelerated launcher

## End-to-End Workflow

### build_fep_system_from_prism_ligands()
**Convenience function** for complete workflow automation:

```python
from prism.fep.modeling.e2e import build_fep_system_from_prism_ligands

output_dir = build_fep_system_from_prism_ligands(
    receptor_pdb="receptor.pdb",
    reference_ligand_dir="ligand_a_ff",
    mutant_ligand_dir="ligand_b_ff",
    output_dir="fep_output",
    lambda_windows=32,
    lambda_strategy="decoupled",
    dist_cutoff=0.6,
    charge_cutoff=0.05
)
```

**Steps**:
1. Read PRISM ligands (auto-discover LIG.amb2gmx/ subdirectory)
2. Perform atom mapping with DistanceAtomMapper
3. Build FEP scaffold with FEPScaffoldBuilder
4. Return path to output directory

### Command Line Interface
```bash
python -m prism.fep.modeling.e2e \
  --protein receptor.pdb \
  --reference ref.mol2 \
  --mutant mut.mol2 \
  --output fep_output \
  --lambda-windows 32 \
  --dist-cutoff 0.6
```

## Coordinate File Support

### MOL2 Support (New)
**Feature**: Use original ligand coordinates (MOL2/PDB) instead of GRO

**Benefits**:
- Better accuracy for ligand alignment
- Preserves bond order information
- Avoids GRO rounding issues

**Implementation**:
```python
# Auto-discover coordinate files
coord_files = ["ligand.mol2", "ligand.pdb", "processed.gro"]
for coord_file in coord_files:
    if Path(coord_file).exists():
        use this file
```

**Code location**: `prism/fep/io.py::_load_structure_coordinates()`

### Unified Coordinate Handling
**Function**: `_load_structure_coordinates()`

**Supports**:
- PDB files (standard)
- MOL2 files (new, with element extraction)
- GRO files (fallback)

**Element extraction from MOL2**:
```python
# Use atom name (not atom type) for element symbol
element = mol2_atom.name  # e.g., "C1" → "C"
```

## Lambda Strategies

### decoupled (Default)
**Separate A and B state decoupling**:
- λ = 0: Fully coupled state A
- λ = 1: Fully coupled state B
- Both states transform simultaneously

### separated
**Separate A→B and B→A transformations**:
- Forward leg: A state turns off, B state turns on
- Reverse leg: B state turns off, A state turns on

## Configuration

### YAML Configuration
```yaml
fep:
  mode: standard  # or 'repex'
  lambda_windows: 32
  lambda_strategy: decoupled
  mapping:
    dist_cutoff: 0.6
    charge_cutoff: 0.05
    charge_common: mean
    charge_reception: surround
```

### Programmatic Configuration
```python
from prism.fep.config import FEPConfig

config = FEPConfig(
    system_dir=".",
    config_file="config_gaff.yaml",
    fep_file="fep.yaml"
)

builder = FEPScaffoldBuilder(
    protein_path="receptor.pdb",
    ligand_ref_path="ref.mol2",
    ligand_mut_path="mut.mol2",
    config=config
)
```

## Common Issues

### Missing MOL2 Coordinates
**Problem**: GRO coordinates used instead of MOL2

**Solution**:
1. Check MOL2 file exists in ligand directory
2. Verify `builder.py` prioritizes MOL2 over GRO
3. Use `reference_coord_file` parameter in e2e workflow

**Code**: `builder/core.py::FEPScaffoldBuilder._get_ligand_coords()`

### [defaults] Block Error
**Problem**: "Found a second defaults directive" error

**Solution**: Extract [defaults] block to separate file

**Code**: `modeling/hybrid_package.py::_extract_defaults_block()`

### Wrong Hybrid GRO for Unbound Leg
**Problem**: Placeholder GRO used instead of real hybrid coordinates

**Solution**: Use `hybrid_package.gro` from HybridPackageBuilder

**Code**: `modeling/core.py::FEPScaffoldBuilder.build_from_ligands()`

## Testing

### Unit Tests
```bash
# Test scaffold building
python tests/gxf/FEP/unit_test/test_fep_scaffold.py

# Test e2e workflow
python tests/gxf/FEP/unit_test/test_fep_e2e.py
```

### Verification Checklist
- [ ] Directory structure created correctly
- [ ] Hybrid topology files present (hybrid.itp, atomtypes_hybrid.itp)
- [ ] Lambda window directories created (window_00 to window_N)
- [ ] Run scripts generated (standard and repex)
- [ ] MOL2 coordinates used (not GRO fallback)

## Code Locations

- **Scaffold builder**: `prism/fep/modeling/core.py`
- **Hybrid package**: `prism/fep/modeling/hybrid_package.py`
- **Script writer**: `prism/fep/modeling/script_writer.py`
- **E2E workflow**: `prism/fep/modeling/e2e.py`
- **Leg writer**: `prism/fep/modeling/leg_writer.py`

## Related Documentation

- `CLAUDE.md` (parent) - FEP module overview
- `../core/CLAUDE.md` - Atom mapping details
- `../gromacs/CLAUDE.md` - GROMACS integration
