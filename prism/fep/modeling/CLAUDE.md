# PRISM-FEP Modeling Module

This file provides guidance for working with FEP system building and workflow orchestration.

## Directory Naming Rule (CRITICAL)

- Use simple case names for force-field combinations: `<protein_ff>-mut_<ligand_ff>` (example: `charmm36m-mut_mmff`).
- Never create nested duplicate system directories like `.../GMX_PROLIG_FEP/GMX_PROLIG_FEP/`.
- Avoid ad-hoc suffixes (`_pkgfix*`, `_final*`, `_new*`) in directory names.
- Default output directory for FEP cases should follow `<protein_ff>-mut_<ligand_ff>` when users do not pass an explicit output path.
- For 42-38 maintenance, move legacy non-canonical outputs into `tests/gxf/FEP/unit_test/42-38/Archive/`.

## Module Overview

The modeling module orchestrates the complete FEP system building workflow:
- **Scaffold building**: Creates directory structure for bound/unbound legs
- **Hybrid topology generation**: Combines A-state and B-state parameters
- **Script generation**: Produces run scripts for standard and repex modes
- **End-to-end workflow**: Automates from ligands to simulation-ready systems

## Critical Development Constraints

### Hybrid Topology Coordinate Assignment (CRITICAL)

**Problem**: Incorrect coordinate assignment causes EM crashes with 6+ nm bond distances

**Root Cause**:
- Hybrid ITP atom names: HA, OA, C, H1, H2, C1, N, ... (specific)
- Mapping file atom names: H, O, C, H, H, C, N, ... (generic)
- Coordinate lookup fails → placeholder (0,0,0) → EM crash

**Mandatory Requirements**:
1. **workflow_fep.py (lines 254-272)**: Read atoms from hybrid ITP for mapping files
   ```python
   # ✅ CORRECT: Read from hybrid ITP
   hybrid_atoms_for_mapping_a = read_ligand_from_prism(
       itp_file=hybrid_itp,
       gro_file=hybrid_gro,
       state="a"
   )
   ```

2. **hybrid_package.py (lines 202-208)**: Always use original ligand GRO as coordinate source
   ```python
   # ✅ CORRECT: Always use original ligand GRO
   reference_coords_file = seed_gro    # NOT mapping_state_a.pdb
   mutant_coords_file = mutant_gro    # NOT mapping_state_b.pdb
   ```

3. **hybrid_package.py (lines 975-982)**: Only support pure element names
   ```python
   # ✅ CORRECT: Only pure element names (HA, OA, CA, N, F)
   if atom_name.isalpha() and len(atom_name) <= 3:
       element = self._extract_element_from_atom_name(atom_name)

   # ❌ WRONG: Don't support numbered names via element sequential matching
   # if (atom_name.isalpha() and len(atom_name) <= 3) or \
   #    (not atom_name.isalpha() and atom_name[0].isalpha()):
   ```

**Verification After Build**:
```python
# Check hybrid.gro coordinate uniqueness
coords = []
with open("hybrid.gro") as f:
    for line in f.readlines()[2:-1]:
        x, y, z = float(line[20:28]), float(line[28:36]), float(line[36:44])
        coords.append((x, y, z))

assert len(set(coords)) == len(coords), "All atoms must have unique coordinates"
assert sum(1 for c in coords if c == (0.0, 0.0, 0.0)) == 0, "No placeholder coordinates"
```

### Coordinate File Handling Standards

**Priority Order** (MUST follow):
1. **MOL2** - Original ligand coordinates with bond orders
2. **PDB** - Original ligand coordinates with full atom names
3. **GRO** - Fallback only (5-char atom name limit)

**Implementation** (hybrid_package.py):
```python
# ✅ CORRECT: Priority-based coordinate loading
def _load_structure_coordinates(coord_path):
    # Try MOL2 first (preserves full atom names)
    if coord_path.with_suffix(".mol2").exists():
        return self._parse_mol2_atoms(coord_path.with_suffix(".mol2"))

    # Try PDB second (preserves full atom names)
    if coord_path.with_suffix(".pdb").exists():
        return self._parse_pdb_atoms(coord_path.with_suffix(".pdb"))

    # Fall back to GRO (may have truncated names)
    return self._parse_gro_atoms(coord_path)
```

**Forbidden Practices**:
- ❌ Using `mapping_state_a.pdb`/`mapping_state_b.pdb` as coordinate sources
- ❌ Prioritizing mapping files over original ligand files
- ❌ Assuming GRO atom names match hybrid ITP atom names

### Script Generation Standards

**Forbidden Practices**:
- ❌ Generate excessive redundant scripts (standard/repex/local/parallel all separate)
- ❌ Nest scripts inside leg directories
- ❌ Hardcode temporary values (e.g., "100 steps")

**Correct Approach**:
- ✅ Generate universal scripts controlled by parameters
- ✅ Place bound/unbound scripts at sibling level
- ✅ Implement short tests via configuration

**Script Structure**:
```
GMX_PROLIG_FEP/
├── bound/
│   ├── run_fep.sh          # Universal script (handles all modes)
│   └── window_00/ to window_N/
└── unbound/
    ├── run_fep.sh          # Universal script (handles all modes)
    └── window_00/ to window_N/
```

**Runtime Flexibility**:
- Support `--replica-range 1-2` for partial runs
- Support `--replica 1` for single replica
- Support `--mode standard|repex` for mode switching
- Don't generate separate scripts for each mode

### Testing and Validation Requirements

**Test File Standards**:
- ✅ All parameters must come from `case.yaml`
- ❌ No hardcoded parameters in test files
- ✅ Must simulate real user workflow

**Validation After Build**:
```python
# Required validation checks
assert os.path.exists("hybrid.itp"), "hybrid.itp must exist"
assert os.path.exists("hybrid.gro"), "hybrid.gro must exist"
assert os.path.exists("topol.top"), "topol.top must exist"

# Check coordinate uniqueness
coords = parse_gro("hybrid.gro")
assert len(set(coords)) == len(coords), "All coordinates must be unique"

# Check for placeholder coordinates
zero_coords = sum(1 for c in coords if c == (0.0, 0.0, 0.0))
assert zero_coords == 0, "No placeholder (0,0,0) coordinates allowed"
```

### Force Field Integration Constraints

**OpenFF/OPLS-AA Limitations**:
- ❌ Cannot rely on atom type for position distinction
- ✅ Must ignore atom type in DistanceAtomMapper
- ✅ Only use distance + element + charge for matching

**RTF/CGenFF Requirements**:
- ✅ Must integrate original CGenFF files (toppar_c36_feb26)
- ✅ GROMACS official charmm36.ff lacks small molecule parameters
- ✅ RTF priority: MATCH > SwissParam > CGenFF website

**Protonation State Handling**:
- ✅ Auto-map residue names based on RTP file
- ✅ Support alias mapping (HID→HSD, HIE→HSE)
- ✅ Warning on unmapped residues
- ❌ Don't hardcode Amber-specific logic

### Configuration Management

**Parameter Priority** (MUST follow):
1. CLI arguments (highest)
2. fep.yaml (FEP-specific)
3. config.yaml (general)
4. Code defaults (lowest)

**Required Parameters**:
- dist_cutoff: Distance cutoff for mapping (default: 0.6 nm)
- charge_cutoff: Charge cutoff for mapping (default: 0.05)
- lambda_windows: Number of lambda windows (default: 11)
- lambda_strategy: decoupled or separated

**Lambda Strategy Clarification**:
- **decoupled**: Both A and B states transform simultaneously
- **separated**: Forward (A→B) and reverse (B→A) legs
- GROMACS only supports single topology (not dual topology)

### Engineering Standards

**Code Quality Requirements**:
- ✅ No hardcoded magic numbers
- ✅ Clear separation of concerns
- ✅ Consistent naming conventions
- ✅ Proper error handling and logging
- ❌ No Chinese punctuation in code
- ❌ No hardcoded selection/resname

**File Organization**:
- Keep modules focused (no God files >900 lines)
- Avoid unnecessary abstraction layers
- Delete unused code (don't keep "just in case")
- Don't create duplicate functionality

**Documentation Requirements**:
- ✅ All public functions have docstrings
- ✅ Complex logic has inline comments
- ✅ Examples reflect actual usage
- ❌ No undocumented features

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
