# PRISM Force Field Generators

Unified interface for generating ligand force field parameters for GROMACS simulations.

## Available Force Fields

| Force Field | Generator Class | Output Directory | Method |
|------------|----------------|------------------|--------|
| **GAFF** | `GAFFForceFieldGenerator` | `LIG.amb2gmx` | Local (AmberTools) |
| **OpenFF** | `OpenFFForceFieldGenerator` | `LIG.openff2gmx` | Local |
| **OPLS-AA** | `OPLSAAForceFieldGenerator` | `LIG.opls2gmx` | Server (LigParGen) |
| **MMFF** | `MMFFForceFieldGenerator` | `LIG.mmff2gmx` | Server (SwissParam) |
| **MATCH** | `MATCHForceFieldGenerator` | `LIG.match2gmx` | Server (SwissParam) |
| **Hybrid** | `HybridMMFFMATCHForceFieldGenerator` | `LIG.hybrid2gmx` | Server (SwissParam) |

## Installation

```bash
# Base dependencies
pip install requests mechanize

# GAFF - AmberTools
mamba install -c conda-forge ambertools

# OpenFF
mamba install -c conda-forge openff-toolkit openff-interchange

# OPLS-AA - Optional RDKit for alignment
mamba install -c conda-forge rdkit
```

## Usage

All generators follow the same pattern:

```python
from prism.forcefield import <GeneratorClass>

generator = <GeneratorClass>(
    ligand_path="ligand.mol2",  # or .sdf
    output_dir="output",
    overwrite=False
)

result_dir = generator.run()
```

### Examples

**GAFF:**
```python
from prism.forcefield import GAFFForceFieldGenerator

generator = GAFFForceFieldGenerator("ligand.mol2", "output")
result_dir = generator.run()  # output/LIG.amb2gmx/
```

**OpenFF:**
```python
from prism.forcefield import OpenFFForceFieldGenerator

generator = OpenFFForceFieldGenerator(
    ligand_path="ligand.sdf",
    output_dir="output",
    charge=0,
    forcefield="openff-2.1.0"
)
result_dir = generator.run()  # output/LIG.openff2gmx/
```

**OPLS-AA:**
```python
from prism.forcefield import OPLSAAForceFieldGenerator

generator = OPLSAAForceFieldGenerator(
    ligand_path="ligand.mol2",
    output_dir="output",
    charge=0,
    charge_model="cm1a",  # or "cm5"
    align_to_input=True
)
result_dir = generator.run()  # output/LIG.opls2gmx/
```

**MMFF / MATCH / Hybrid:**
```python
from prism.forcefield import MMFFForceFieldGenerator, MATCHForceFieldGenerator, HybridMMFFMATCHForceFieldGenerator

# MMFF-based
generator = MMFFForceFieldGenerator("ligand.mol2", "output")
result_dir = generator.run()  # output/LIG.mmff2gmx/

# MATCH
generator = MATCHForceFieldGenerator("ligand.mol2", "output")
result_dir = generator.run()  # output/LIG.match2gmx/

# Hybrid
generator = HybridMMFFMATCHForceFieldGenerator("ligand.mol2", "output")
result_dir = generator.run()  # output/LIG.hybrid2gmx/
```

## Output Files

All generators produce standardized GROMACS files:

```
output/LIG.<ff>2gmx/
├── LIG.gro              # Coordinates
├── LIG.itp              # Molecule topology
├── LIG.top              # Complete topology
├── atomtypes_LIG.itp    # Atom types
└── posre_LIG.itp        # Position restraints
```

## Utility Functions

```python
from prism.forcefield import list_available_generators, get_generator_info

# List available generators
available = list_available_generators()

# Get detailed info
info = get_generator_info()
for name, details in info.items():
    print(f"{name}: {details['description']}")
```

## Input Formats

| Generator | MOL2 | SDF | PDB |
|-----------|------|-----|-----|
| GAFF, OpenFF, OPLS-AA | ✓ | ✓ | ✗ |
| MMFF, MATCH, Hybrid | ✓ | ✓ | ✓ |

## Troubleshooting

**Missing Dependencies:**
```bash
# Install missing packages as needed
pip install requests mechanize
mamba install -c conda-forge ambertools openff-toolkit rdkit
```

**Server Timeout:** Server-based generators (OPLS-AA, MMFF, MATCH, Hybrid) require internet connection. Large molecules may take longer to process.

**AmberTools:**
```bash
# Verify installation
which antechamber parmchk2 tleap acpype
```

## Command-Line Usage

```bash
python -m prism.forcefield.gaff ligand.mol2
python -m prism.forcefield.openff ligand.sdf
python -m prism.forcefield.opls_aa ligand.mol2
python -m prism.forcefield.swissparam ligand.mol2 MMFF-based
```

## References

- **GAFF**: Wang et al., J. Comput. Chem. 25: 1157-1174 (2004)
- **OpenFF**: [openforcefield.org](https://openforcefield.org)
- **OPLS-AA**: Jorgensen et al., JACS 118: 11225-11236 (1996) | [LigParGen](http://zarbi.chem.yale.edu/ligpargen/)
- **SwissParam**: Zoete et al., J. Comput. Chem. 32: 2359-2368 (2011) | [swissparam.ch](http://www.swissparam.ch)
