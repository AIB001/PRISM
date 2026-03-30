# PRISM-FEP Module

Free Energy Perturbation (FEP) module for relative binding free energy calculations in GROMACS.

## Features

- **Automated atom mapping** via DistanceAtomMapper
- **Hybrid topology generation** (A-state + B-state parameters)
- **Lambda window setup** for alchemical transformation
- **HTML visualization** for mapping verification
- **CHARMM-GUI integration** for force field compatibility

## Quick Start

### CLI Usage

```bash
# Basic FEP calculation
prism protein.pdb reference.mol2 -o output --fep --mutant mutant.mol2

# With custom parameters
prism protein.pdb reference.mol2 -o output \
  --fep \
  --mutant mutant.mol2 \
  --distance-cutoff 0.6 \
  --charge-strategy mean \
  --lambda-windows 11
```

### Python API

```python
from prism.fep import FEPScaffoldBuilder

builder = FEPScaffoldBuilder(
    protein_path="protein.pdb",
    ligand_ref_path="reference.mol2",
    ligand_mut_path="mutant.mol2",
    output_dir="GMX_PROLIG_FEP"
)

scaffold = builder.build_from_ligands()
```

### End-to-End Workflow

```python
# Complete workflow from ligands to simulation-ready system
python -m prism.fep.modeling.e2e \
  --protein protein.pdb \
  --reference ref.mol2 \
  --mutant mut.mol2 \
  --output fep_output
```

## Architecture

```
prism/fep/
├── core/           # Core FEP workflow and mapping
├── modeling/       # Hybrid topology generation
├── gromacs/        # GROMACS integration and MDP templates
├── visualize/      # HTML/PNG visualization
├── io.py          # File I/O for hybrid topologies
└── config.py      # Configuration management
```

## Key Classes

- **DistanceAtomMapper** - Atom mapping based on distance, element, type, charge
- **FEPScaffoldBuilder** - Creates FEP directory structure (bound/unbound legs)
- **HybridPackageBuilder** - Generates hybrid topology files
- **FEPConfig** - Configuration management from YAML files

## Configuration

Create `fep.yaml` in your project directory:

```yaml
fep:
  distance_cutoff: 0.6      # Atom mapping distance threshold (Å)
  charge_cutoff: 0.05       # Charge difference threshold (e)
  charge_common: mean       # Common atom charge strategy (ref/mut/mean)
  charge_reception: surround # Transformed atom charge strategy
  lambda_windows: 11        # Number of lambda windows
```

## Output Structure

```
GMX_PROLIG_FEP/
├── bound/              # Protein-ligand complex leg
│   ├── build/          # EM, NVT, NPT equilibration
│   ├── window_XX/      # Lambda window directories
│   └── localrun.sh     # GPU-accelerated run script
├── unbound/            # Ligand in water leg
│   └── (same structure as bound)
└── common/
    └── hybrid/         # Hybrid topology files
        ├── hybrid.itp
        ├── atomtypes_hybrid.itp
        └── hybrid.gro
```

## Visualization

Generate HTML visualization to verify atom mapping:

```python
from prism.fep.visualize.html import visualize_mapping_html
from prism.fep.core.mapping import DistanceAtomMapper

# Create mapping
mapper = DistanceAtomMapper(dist_cutoff=0.6)
mapping = mapper.map(atoms_ref, atoms_mut)

# Generate HTML
visualize_mapping_html(
    mapping=mapping,
    pdb_a="ref.pdb",
    pdb_b="mut.pdb",
    atoms_a=atoms_ref,
    atoms_b=atoms_mut,
    output_path="mapping.html"
)
```

**Quality check**: No gray atoms, all atoms classified, total charge ≈ 0

## Testing

Recommended test systems:
- **42-38**: Single-point mutation ✅ (Recommended)
- **39-8**: Multi-point mutations
- **oMeEtPh-EtPh**: Terminal ethyl methyl transformation

```bash
# Run tests
python tests/gxf/FEP/unit_test/test_42_38.py
python tests/gxf/FEP/unit_test/test_39_8.py
```

## Documentation

- `CLAUDE.md` - Development guidelines and troubleshooting
- `NAMING_CONVENTIONS.md` - Standardized naming conventions
- `../../tests/gxf/FEP/unit_test/CLAUDE.md` - Test-specific guidelines

## Status

✅ Atom mapping implemented
✅ Hybrid topology generation working
✅ HTML visualization with quality checks
✅ CHARMM-GUI integration
✅ End-to-end workflow functional

## References

- GROMACS FEP Tutorial: http://www.gromacs.org/Documentation/Tutorials
- Alchemical Free Energy Calculations: https://www.alchemistry.org/
