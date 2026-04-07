# PRISM-FEP Core Module

This file provides guidance for working with the core FEP data structures and atom mapping algorithms.

## Directory Naming Rule (CRITICAL)

- Use simple case names for force-field combinations: `<protein_ff>-mut_<ligand_ff>` (example: `charmm36m-mut_mmff`).
- Never create nested duplicate system directories like `.../GMX_PROLIG_FEP/GMX_PROLIG_FEP/`.
- Avoid ad-hoc suffixes (`_pkgfix*`, `_final*`, `_new*`) in directory names.
- Default output directory for FEP cases should follow `<protein_ff>-mut_<ligand_ff>` when users do not pass an explicit output path.
- For 42-38 maintenance, move legacy non-canonical outputs into `tests/gxf/FEP/unit_test/42-38/Archive/`.

## Module Overview

The core module implements the fundamental data structures and algorithms for FEP calculations:
- **Atom dataclass**: Represents atoms with position, charge, type, element
- **AtomMapping dataclass**: Result of atom mapping (common, transformed, surrounding)
- **DistanceAtomMapper**: Distance-based atom mapping algorithm

## Key Classes

### Atom
```python
@dataclass
class Atom:
    name: str          # Unique identifier (e.g., 'C1', 'H2')
    element: str       # Element symbol (C, H, N, O, etc.)
    coord: np.ndarray  # 3D coordinates [x, y, z] in Angstroms
    charge: float      # Partial charge from force field
    atom_type: str     # Force field atom type (e.g., 'ca', 'ha', 'na')
    index: int         # Atom index for referencing
```

**Usage**: Atom objects are the basic building blocks for FEP calculations.

### AtomMapping
```python
@dataclass
class AtomMapping:
    common: List[Tuple[Atom, Atom]]      # Shared atoms
    transformed_a: List[Atom]            # Unique to ligand A
    transformed_b: List[Atom]            # Unique to ligand B
    surrounding_a: List[Atom]            # Position-matched, divergent params in A
    surrounding_b: List[Atom]            # Position-matched, divergent params in B
```

**Usage**: This is the core result structure - all downstream processing depends on this classification.

### DistanceAtomMapper

**Parameters**:
- `dist_cutoff`: Distance threshold in Å (default: 0.6)
  - Smaller (0.3-0.5): Very similar molecules
  - Larger (0.8-1.0): Structurally divergent molecules
- `charge_cutoff`: Charge difference threshold (default: 0.05)
- `charge_common`: Strategy for common atom charges (ref/mut/mean/none)
- `charge_reception`: Strategy for surplus charges (unique/surround/none)
- `recharge_hydrogen`: Whether to perturb hydrogen charges (default: False)

**Algorithm** (7 steps):
1. Distance matching: `dist < dist_cutoff`
2. Element matching: `atom_a.element == atom_b.element`
3. Type matching: `atom_a.type == atom_b.type` (OpenFF/OPLS excepted)
4. Charge matching: `abs(atom_a.charge - atom_b.charge) < charge_cutoff`
5. Determine common pairs
6. Determine transformed (unique atoms)
7. Determine surrounding (position-matched but divergent parameters)

**Force Field Auto-Detection**:
```python
def _should_ignore_atom_type(self, atoms: List[Atom]) -> bool:
    """
    Generic/Sequential types (skip atom type check):
    - OpenFF: output_0, output_1, ...
    - OPLS-AA: opls_800, opls_801, ...

    Position-specific types (use atom type check):
    - GAFF: ca, c3, ha, hc, ...
    - CGenFF: CA, CB, CG, ...
    """
```

## Charge Handling Strategies

### charge_common
Controls how common atom charges are set in the hybrid topology:
- **"ref"**: Use reference ligand charges
- **"mut"**: Use mutant ligand charges
- **"mean"**: Use average of both charges (default)
- **"none"**: Keep original charges (no modification)

### charge_reception
Controls how surplus charges are distributed:
- **"unique"**: Only unique (non-hydrogen) atoms receive charge
- **"surround"**: Only surrounding atoms receive charge (default)
- **"surround_ext"**: Auto-extend to common atoms if charge/atom > 0.02
- **"none"**: No charge redistribution

**Legacy label normalization**: `"pert"` → `"surround"`

## Configuration

### From YAML
```yaml
fep:
  mapping:
    dist_cutoff: 0.6
    charge_cutoff: 0.05
    charge_common: mean
    charge_reception: surround
    recharge_hydrogen: false
```

### From Code
```python
from prism.fep.core.mapping import DistanceAtomMapper

mapper = DistanceAtomMapper(
    dist_cutoff=0.6,
    charge_cutoff=0.05,
    charge_common="mean",
    charge_reception="surround"
)

mapping = mapper.map(ligand_a_atoms, ligand_b_atoms)
```

### From Config Dict
```python
mapper = DistanceAtomMapper.from_config(config_dict)
```

## Common Issues

### Type Mismatch with OpenFF/OPLS
**Problem**: Atoms not matched despite correct positions

**Solution**: The mapper auto-detects OpenFF/OPLS and skips type checking. Check atom type names:
- OpenFF: `output_0`, `output_1`, ...
- OPLS: `opls_800`, `opls_801`, ...

### Too Many Surrounding Atoms
**Problem**: More surrounding atoms than expected

**Causes**:
- `charge_cutoff` too strict (try 0.1 instead of 0.05)
- Real parameter differences between force fields
- Protonation state differences

### Charge Imbalance
**Problem**: Total charge ≠ 0 after mapping

**Solution**:
1. Check `charge_common` strategy
2. Verify `charge_reception` is not "none"
3. Ensure `recharge_hydrogen=False` for neutral molecules

## Testing

### Unit Tests
```bash
python tests/gxf/FEP/unit_test/test_gaff_fep_mapping.py
```

### Expected Results
For oMeEtPh-EtPh system:
- **REF/MEAN/MUT**: Common=17, TransA=4, TransB=1, SurrA=0, SurrB=0
- **NONE**: Common=4, TransA=4, TransB=1, SurrA=13, SurrB=13

## Code Locations

- **Data structures**: `prism/fep/core/mapping.py:1-100`
- **DistanceAtomMapper**: `prism/fep/core/mapping.py:103-400`
- **Hybrid topology**: `prism/fep/core/hybrid_topology.py`

## Related Documentation

- `CLAUDE.md` (parent) - FEP module overview
- `../modeling/CLAUDE.md` - System building workflow
- `.claude/skills/fep-mapping.md` - Mapping algorithm verification guide
