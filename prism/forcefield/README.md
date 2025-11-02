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
| **Both** | `BothForceFieldGenerator` | `LIG.both2gmx` | Server (SwissParam) |

## Installation

```bash
# Base dependencies
pip install requests

# GAFF - AmberTools
mamba install -c conda-forge ambertools

# OpenFF
mamba install -c conda-forge openff-toolkit openff-interchange

# OPLS-AA - Optional RDKit for alignment
mamba install -c conda-forge rdkit
```

## Usage

### High-Level API (Recommended)

```python
import prism as pm

# GAFF
system = pm.system("protein.pdb", "ligand.mol2", ligand_forcefield="gaff")
system.build()

# OpenFF
system = pm.system("protein.pdb", "ligand.sdf", ligand_forcefield="openff")
system.build()

# OPLS-AA
system = pm.system("protein.pdb", "ligand.mol2", ligand_forcefield="opls")
system.build()

# SwissParam - MMFF
system = pm.system("protein.pdb", "ligand.mol2", ligand_forcefield="mmff")
system.build()

# SwissParam - MATCH
system = pm.system("protein.pdb", "ligand.mol2", ligand_forcefield="match")
system.build()

# SwissParam - Both (MMFF + MATCH)
system = pm.system("protein.pdb", "ligand.mol2", ligand_forcefield="both")
system.build()
```

### Direct Generator Usage

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

**SwissParam - MMFF / MATCH / Both:**
```python
from prism.forcefield.swissparam import (
    MMFFForceFieldGenerator,
    MATCHForceFieldGenerator,
    BothForceFieldGenerator
)

# MMFF-based
generator = MMFFForceFieldGenerator("ligand.mol2", "output")
result_dir = generator.run()  # output/LIG.mmff2gmx/

# MATCH
generator = MATCHForceFieldGenerator("ligand.mol2", "output")
result_dir = generator.run()  # output/LIG.match2gmx/

# Both (MMFF + MATCH)
generator = BothForceFieldGenerator("ligand.mol2", "output")
result_dir = generator.run()  # output/LIG.both2gmx/
```

**SwissParam - Convenience Function:**
```python
from prism.forcefield.swissparam import generate_swissparam_ff

result_dir = generate_swissparam_ff(
    ligand_path="ligand.mol2",
    output_dir="output",
    approach="mmff-based",  # or "match", "both"
    overwrite=True
)
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

## SwissParam Special Notes

⚠️ **Request Limits**: SwissParam uses a web service with rate limiting
- Wait a few minutes between consecutive requests
- If you get "could not be submitted to the queuing system" error, wait and retry

⚠️ **Network Required**: Requires internet connection to access SwissParam server

**Supported Methods:**

| Method | Description | Use Case |
|--------|-------------|----------|
| MMFF | Based on MMFF94 force field | Fast parameterization for most drug-like molecules |
| MATCH | Based on MATCH parameters | More accurate charges and geometric parameters |
| Both | MMFF + MATCH combined | Most comprehensive parameter set |

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
| MMFF, MATCH, Both | ✓ | ✓ | ✓ |

## Troubleshooting

**Missing Dependencies:**
```bash
# Install missing packages as needed
pip install requests
mamba install -c conda-forge ambertools openff-toolkit rdkit
```

**Server Timeout:** Server-based generators (OPLS-AA, MMFF, MATCH, Both) require internet connection. Large molecules may take longer to process.

**SwissParam Rate Limiting:** If you encounter errors like "could not be submitted to the queuing system", you've likely hit the rate limit. Wait 5-10 minutes before retrying.

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
python -m prism.forcefield.swissparam ligand.mol2 mmff-based
```

## References

- **GAFF**: Wang et al., J. Comput. Chem. 25: 1157-1174 (2004)
- **OpenFF**: [openforcefield.org](https://openforcefield.org)
- **OPLS-AA**: Jorgensen et al., JACS 118: 11225-11236 (1996) | [LigParGen](http://zarbi.chem.yale.edu/ligpargen/)
- **SwissParam**: Zoete et al., J. Comput. Chem. 32: 2359-2368 (2011) | [swissparam.ch](http://www.swissparam.ch)
