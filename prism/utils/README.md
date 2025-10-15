# PRISM System Building Utilities

**For PRISM Team Members**

This document describes the key utilities for building GROMACS MD systems with PRISM, with focus on the new metal ion support and system building features.

## Quick Start

```python
from prism.builder import PRISMBuilder

# Build a protein-ligand system with metal ions
builder = PRISMBuilder(
    protein_path='protein.pdb',
    ligand_path='ligand.mol2',
    output_dir='output',
    forcefield='charmm36',
    water_model='tip4p'
)

# Clean protein (removes artifacts, preserves structural metals)
cleaned = builder.clean_protein()

# Generate ligand force field
lig_ff = builder.generate_ligand_ff()

# Build complete MD system
model_dir = builder.build_model(cleaned, lig_ff)
```

## Protein Cleaning (`cleaner.py`)

### Overview

The `ProteinCleaner` class intelligently processes PDB files to prepare them for MD simulations. It handles:

- **Metal ions** - Distinguishes structural metals (Zn, Mg, Ca) from buffer ions (Na, Cl)
- **Crystal water** - Optional retention of crystallographic water molecules (default: remove)
- **Crystallization artifacts** - Removes glycerol, PEG, detergents, sugars
- **Distance filtering** - Removes distant metals that may be crystallization artifacts
- **Flexible control** - Fine-tune what gets kept or removed for your specific research needs

### Three Ion Handling Modes

**1. `smart` (Default - Recommended)**
- ✅ Keeps: Structural metals (Zn, Mg, Ca, Fe, Cu, Mn, Co, Ni, Cd)
- ❌ Removes: Buffer ions (Na, Cl, K, SO4, PO4)
- ❌ Removes: Crystallization artifacts (GOL, EDO, PEG, NAG)

**2. `remove_all`**
- ❌ Removes: ALL metals and ions
- Use when: Metals are not important for your study

**3. `keep_all`**
- ✅ Keeps: ALL metals and ions
- Use when: You need complete control

### Usage Examples

#### Basic Usage (via PRISMBuilder)

```python
from prism.builder import PRISMBuilder

builder = PRISMBuilder(
    protein_path='1ABC.pdb',
    ligand_path='ligand.mol2',
    output_dir='output'
)

# Use default smart mode (removes water, keeps structural metals)
cleaned = builder.clean_protein()

# Specify all cleaning options
cleaned = builder.clean_protein(
    ion_mode='smart',              # keep structural metals
    distance_cutoff=5.0,           # remove metals >5Å from protein
    keep_crystal_water=False,      # remove water molecules
    remove_artifacts=True          # remove GOL, EDO, etc.
)
```

#### Crystal Water Processing Examples

```python
# Example 1: Protein with active site water (e.g., enzyme catalysis)
builder = PRISMBuilder(
    protein_path='enzyme_with_water.pdb',
    ligand_path='substrate.mol2',
    output_dir='output_with_water'
)

# Keep crystal water to preserve active site hydration
cleaned = builder.clean_protein(
    ion_mode='smart',
    keep_crystal_water=True,       # Keep all HOH molecules
    remove_artifacts=True
)

# Example 2: Clean structure without any heteroatoms except metals
builder = PRISMBuilder(
    protein_path='3HS4.pdb',
    ligand_path='AZM.sdf',
    output_dir='output_clean'
)

cleaned = builder.clean_protein(
    ion_mode='smart',              # Keep Zn (structural)
    keep_crystal_water=False,      # Remove all 404 HOH
    remove_artifacts=True          # Remove GOL
)

# Example 3: Study protein-water interactions, keep everything
cleaned = builder.clean_protein(
    ion_mode='keep_all',           # Keep all ions
    keep_crystal_water=True,       # Keep all water
    remove_artifacts=False         # Keep artifacts too
)
```

#### Advanced Usage (Direct ProteinCleaner)

```python
from prism.utils.cleaner import ProteinCleaner

cleaner = ProteinCleaner(
    ion_mode='smart',
    distance_cutoff=5.0,
    keep_crystal_water=False,
    remove_artifacts=True,
    keep_custom_metals=['MO', 'W'],     # Additional metals to keep
    remove_custom_metals=['CA'],        # Override: remove Ca
    verbose=True
)

cleaner.clean_pdb('input.pdb', 'output_clean.pdb')
```

### Crystal Water Processing

**Default Behavior:** Remove all crystal water molecules (`keep_crystal_water=False`)

Crystal water molecules from X-ray structures are typically removed because:
- They may not be relevant at physiological conditions
- MD simulations will add explicit solvent anyway
- They can cause numbering conflicts with simulation water

**When to Keep Crystal Water:** Set `keep_crystal_water=True` when:
- Water molecules are structurally important (e.g., buried in active site)
- You want to preserve specific protein-water interactions
- Studying hydration networks in binding sites

```python
# Keep all crystal water molecules
cleaned = builder.clean_protein(keep_crystal_water=True)

# Remove crystal water but keep structural metals
cleaned = builder.clean_protein(
    ion_mode='smart',
    keep_crystal_water=False  # Default
)
```

**Example:** 3HS4 (carbonic anhydrase) contains ~404 HOH molecules
- With `keep_crystal_water=False` (default): All removed, clean structure
- With `keep_crystal_water=True`: All 404 HOH retained for analysis

### What Gets Removed/Kept

**Structural Metals (kept in smart mode):**
- Zn²⁺ - Zinc fingers, enzyme active sites
- Mg²⁺ - ATP binding, nucleic acid stability
- Ca²⁺ - Calcium-binding proteins, structural stabilization
- Fe²⁺/Fe³⁺ - Heme groups, iron-sulfur clusters
- Cu²⁺ - Oxidoreductases
- Mn²⁺ - Superoxide dismutase

**Buffer Ions (removed in smart mode):**
- Na⁺, K⁺, Cl⁻ - Common crystallization buffer
- SO₄²⁻, PO₄³⁻ - Sulfate, phosphate
- Br⁻, I⁻, F⁻ - Halide ions

**Crystallization Artifacts (removed when `remove_artifacts=True`):**
- Polyols: GOL (glycerol), EDO (ethylene glycol), MPD
- PEG oligomers: PEG, PGE, 1PE, P6G, etc.
- Sugars: NAG, NDG, BMA, MAN, GAL, FUC
- Detergents: DMS, BOG, LMT, SDS
- Others: ACT (acetate), TRS (Tris), BME (β-mercaptoethanol)

**Water Molecules:**
- Crystal water (HOH) - Controlled by `keep_crystal_water` parameter
- Default: All removed for clean starting structure
- Optional: Keep all for studying hydration networks

## System Building (`system.py`)

### Overview

The `SystemBuilder` class handles the complete GROMACS workflow:
1. Fix protein structure (pdbfixer)
2. Generate topology (pdb2gmx)
3. Add ligand
4. Create simulation box
5. Solvate
6. Add neutralizing ions

### Key Features

**Metal Ion Support:**
- Automatically detects and includes metals processed by pdb2gmx
- Prevents duplication of metal atoms
- Handles metals in any chain (M, D2, F2, etc.)

**Force Field Selection:**
- Uses `-ff` and `-water` flags (more reliable than menu indices)
- Supports custom force field directories
- Handles CHARMM36, Amber, OPLS force fields

**Water Model Support:**
- Correct coordinate files for each water model
- TIP4P virtual site support (tip4p.gro)
- TIP3P, SPC, SPCE compatibility

**Ligand Integration:**
- Automatic topology merging
- Coordinate file combining
- Proper atom numbering

### Usage

```python
from prism.utils.system import SystemBuilder

config = {
    'simulation': {'pH': 7.0},
    'box': {'shape': 'cubic', 'distance': 1.5, 'center': True},
    'ions': {
        'neutral': True,
        'concentration': 0.15,
        'positive_ion': 'NA',
        'negative_ion': 'CL'
    }
}

builder = SystemBuilder(config, output_dir='GMX_PROLIG_MD')

# Build complete system
model_dir = builder.build(
    cleaned_protein_path='protein_clean.pdb',
    lig_ff_dir='LIG.openff2gmx',
    ff_idx=8,
    water_idx=6,
    ff_info={'name': 'charmm36-jul2022', 'dir': 'charmm36-jul2022.ff'},
    water_info={'name': 'tip4p'}
)
```

## Force Field Generation

### OpenFF (`openff.py`)

Modern quantum-mechanics-derived force field for small molecules.

```python
from prism.forcefield.openff import OpenFFGenerator

generator = OpenFFGenerator(
    ligand_path='ligand.mol2',
    output_dir='output',
    charge=0,
    forcefield='openff-2.1.0'
)

lig_ff_dir = generator.run()  # Returns path to LIG.openff2gmx/
```

**Features:**
- High-quality parameters from QM data
- Automatic charge calculation (AM1-BCC)
- MOL2/SDF input support
- Caching support (avoids regeneration)

### GAFF (`gaff.py`)

General Amber Force Field - classical approach for small molecules.

```python
from prism.forcefield.gaff import GAFFGenerator

generator = GAFFGenerator(
    ligand_path='ligand.mol2',
    output_dir='output',
    charge=0,
    multiplicity=1
)

lig_ff_dir = generator.run()  # Returns path to LIG.amb2gmx/
```

**Features:**
- Wide compatibility
- Fast parameter generation
- Requires AmberTools (antechamber, acpype)
- Caching support

## Current Limitations and Known Issues

### ⚠️ PDBFixer + Metal Chain Issue (In Progress)

**Problem:**
- PDBFixer may keep metals in the same chain as protein (e.g., chain A)
- GROMACS pdb2gmx requires ions in separate chains
- Error: "residues in chain do not have consistent type"

**Status:**
- Fix in development
- Will add post-pdbfixer chain separation step
- Metals will be moved to unique chain IDs before pdb2gmx

**Current Workaround:**
- Use structures where metals are already in separate chains
- Or manually separate metal chains in PDB before running PRISM

### ✅ Recently Fixed

1. **Force Field Selection** - Now uses `-ff` and `-water` flags instead of unreliable menu indices
2. **TIP4P Support** - Correct virtual site coordinate files
3. **Metal Duplication** - Checks for pdb2gmx-processed ion chains to prevent duplication
4. **Terminal Atoms** - Automatically fixes O/OXT positioning for GROMACS
5. **Ligand Caching** - Both OpenFF and GAFF respect cached files

---

**Last Updated:** 2025-10-15
**Status:** ✅ Production Ready (with noted limitations)
**Maintainer:** PRISM Development Team
