# Building Systems

This guide covers the complete process of building protein-ligand systems for molecular dynamics simulations using PRISM.

## Quick Start

### Basic System Building

The simplest way to build a system:

```bash
# Build with defaults (GAFF, amber99sb, tip3p)
prism protein.pdb ligand.mol2 -o my_system

# Build with OpenFF for ligand
prism protein.pdb ligand.sdf -o my_system --ligand-forcefield openff

# Build with custom force field
prism protein.pdb ligand.mol2 -o my_system --forcefield amber14sb --water tip4p
```

### Using Python API

```python
import prism

# Create system
system = prism.PRISMSystem(
    "protein.pdb",
    "ligand.mol2",
    output_dir="my_system"
)

# Build
system.build()

# Get output files
files = system.get_output_files()
print(f"System files: {files['gromacs_directory']}")
```

## The Building Process

PRISM performs these steps automatically:

1. **Force field generation** for ligand
2. **Protein cleaning** and preparation
3. **Complex assembly**
4. **Solvation** in water box
5. **Ion addition** for neutralization
6. **Topology generation**
7. **MDP file creation**

## Step-by-Step Building

### Step 1: Force Field Generation

PRISM generates ligand parameters automatically:

```python
from prism import PRISMBuilder

# Initialize builder
builder = PRISMBuilder(
    "protein.pdb",
    "ligand.mol2",
    "output",
    ligand_forcefield="gaff"  # or "openff"
)

# Generate ligand force field only
lig_ff_dir = builder.generate_ligand_forcefield()
print(f"Ligand parameters in: {lig_ff_dir}")
```

**What happens:**
- Charges calculated (AM1-BCC for GAFF)
- Atom types assigned
- Bonded parameters generated
- GROMACS-compatible files created

### Step 2: Protein Preparation

```python
# Clean protein
cleaned_protein = builder.clean_protein()
```

**Automatic cleaning:**
- Removes HETATM records
- Fixes terminal residue names
- Handles HIE/HID/HIP histidine naming

### Step 3: System Assembly

```python
# Build the complete model
model_dir = builder.build_model(cleaned_protein)
```

**Assembly process:**
1. Process protein with pdb2gmx
2. Combine protein and ligand
3. Define simulation box
4. Add solvent
5. Add ions

## Force Field Selection

### Protein Force Fields

Available options depend on your GROMACS installation:

```bash
# List available force fields
prism --list-forcefields
```

Common choices:

| Force Field | Use Case | Water Models |
|------------|----------|--------------|
| amber99sb | General purpose | tip3p, tip4p, spce |
| amber14sb | Improved backbone | tip3p, tip4p |
| amber99sb-ildn | Better side chains | tip3p, tip4p |
| charmm36 | Widely validated | tips3p |
| opls-aa | Good for small molecules | tip3p, tip4p, spce |

### Ligand Force Fields

PRISM supports two approaches:

#### GAFF (General AMBER Force Field)

```bash
prism protein.pdb ligand.mol2 -o output --ligand-forcefield gaff
```

**Advantages:**
- Well-tested
- Good for drug-like molecules
- Compatible with AMBER protein force fields

**Requirements:**
- AmberTools installed
- MOL2 format input

#### OpenFF (Open Force Field)

```bash
prism protein.pdb ligand.sdf -o output --ligand-forcefield openff
```

**Advantages:**
- Modern, systematically improved
- Better coverage of chemical space
- Actively developed

**Requirements:**
- openff-toolkit installed
- SDF or MOL2 format

### Force Field Compatibility

| Protein FF | Compatible Ligand FF | Recommended Water |
|-----------|---------------------|-------------------|
| amber99sb | GAFF, OpenFF | tip3p |
| amber14sb | GAFF, OpenFF | tip3p |
| charmm36 | CGenFF* | tips3p |
| opls-aa | OPLS-AA* | tip4p |

*Not directly supported by PRISM, requires manual setup

## Box Types and Solvation

### Box Shapes

Configure box shape in your config file or command line:

```bash
# Cubic box (default)
prism protein.pdb ligand.mol2 -o output --box-shape cubic

# Dodecahedron (more efficient)
prism protein.pdb ligand.mol2 -o output --box-shape dodecahedron

# Set distance to box edge
prism protein.pdb ligand.mol2 -o output --box-distance 1.5
```

### Box Size Considerations

```python
# Calculate box size needed
import mdtraj as md

# Load protein
traj = md.load("protein.pdb")

# Get dimensions
print(f"Protein size: {traj.xyz.max(axis=1) - traj.xyz.min(axis=1)} nm")

# Recommended box distance
min_distance = 1.0  # nm (minimum)
safe_distance = 1.5  # nm (recommended)
large_distance = 2.0  # nm (for unfolding studies)
```

### Solvation Options

```yaml
# In configuration file
box:
  distance: 1.5
  shape: dodecahedron
  center: true

ions:
  neutral: true
  concentration: 0.15  # Physiological salt
```

## Building Complex Systems

### Multiple Ligands

For systems with multiple ligands:

```python
# Build first ligand system
system1 = prism.PRISMSystem("protein.pdb", "ligand1.mol2", "temp1")
system1.build()

# Build second ligand system
system2 = prism.PRISMSystem("protein.pdb", "ligand2.mol2", "temp2")
system2.build()

# Combine manually (advanced)
# See Advanced Usage guide
```

### Membrane Proteins

Special considerations for membrane proteins:

```python
# 1. Orient protein in membrane first (use OPM, CHARMM-GUI)
# 2. Build membrane system (outside PRISM)
# 3. Add ligand to membrane system

# Example workflow:
# Use CHARMM-GUI for initial membrane setup
# Extract protein-ligand complex
# Use PRISM for ligand parameterization only
builder = PRISMBuilder("protein.pdb", "ligand.mol2", "output")
lig_params = builder.generate_ligand_forcefield()
# Then integrate with membrane system
```

### Protein-Protein Complexes

```python
# Prepare complex PDB with both proteins
# Each chain should be properly labeled

system = prism.PRISMSystem(
    "protein_complex.pdb",
    "ligand.mol2",
    "output",
    config={
        'box': {'distance': 2.0},  # Larger box for complexes
        'simulation': {'production_time_ns': 1000}
    }
)
system.build()
```

## System Validation

### Check Generated Files

After building, verify the output:

```python
# Check output structure
system = prism.PRISMSystem("protein.pdb", "ligand.mol2", "output")
system.build()

# List generated files
files = system.get_output_files()
for name, path in files.items():
    if os.path.exists(path):
        print(f"✓ {name}: {path}")
    else:
        print(f"✗ {name}: MISSING")
```

### Visualize System

```python
import mdtraj as md
import nglview as nv

# Load the solvated system
traj = md.load("output/GMX_PROLIG_MD/solv_ions.gro")

# Quick statistics
print(f"Total atoms: {traj.n_atoms}")
print(f"Box vectors: {traj.unitcell_vectors}")

# Count water molecules
water_atoms = traj.topology.select("water")
print(f"Water molecules: {len(water_atoms) // 3}")

# Count ions
na_ions = traj.topology.select("resname NA")
cl_ions = traj.topology.select("resname CL")
print(f"Na+ ions: {len(na_ions)}, Cl- ions: {len(cl_ions)}")

# Visualize
view = nv.show_mdtraj(traj)
view
```

### Validate Topology

```bash
# Check topology file
cd output/GMX_PROLIG_MD
grep -c "^SOL" topol.top  # Count solvent molecules
grep "qtot" topol.top | tail -1  # Check total charge
```

## Common Building Scenarios

### High-Throughput Screening

```python
import prism
from pathlib import Path

# List of ligands to process
ligands = Path("ligands").glob("*.mol2")

for ligand_file in ligands:
    ligand_name = ligand_file.stem
    
    system = prism.PRISMSystem(
        "protein.pdb",
        str(ligand_file),
        f"output/{ligand_name}",
        config="screening_config.yaml"
    )
    
    try:
        system.build()
        print(f"✓ Built system for {ligand_name}")
    except Exception as e:
        print(f"✗ Failed for {ligand_name}: {e}")
```

### Mutation Studies

```python
# For each mutant
mutants = ["WT", "A123V", "A123L", "A123F"]

for mutant in mutants:
    system = prism.PRISMSystem(
        f"protein_{mutant}.pdb",
        "ligand.mol2",
        f"output/{mutant}",
        config="mutation_study.yaml"
    )
    system.build()
```

### pH Studies

```python
# Build systems at different pH values
for pH in [5.0, 6.0, 7.0, 8.0]:
    # First prepare protein at specific pH
    prepare_at_ph("protein.pdb", f"protein_pH{pH}.pdb", pH)
    
    # Then build system
    system = prism.PRISMSystem(
        f"protein_pH{pH}.pdb",
        "ligand.mol2",
        f"output/pH_{pH}"
    )
    system.build()
```

## Build Options

### Command-Line Options

```bash
# Full command with all options
prism protein.pdb ligand.mol2 \
  --output my_system \
  --ligand-forcefield gaff \
  --forcefield amber14sb \
  --water tip3p \
  --box-shape dodecahedron \
  --box-distance 1.5 \
  --salt-concentration 0.15 \
  --temperature 310 \
  --pressure 1.0 \
  --production-ns 500 \
  --overwrite
```

### Python API Options

```python
system = prism.PRISMSystem(
    "protein.pdb",
    "ligand.mol2",
    output_dir="my_system",
    ligand_forcefield="gaff",
    forcefield="amber14sb",
    water_model="tip3p",
    overwrite=True
)

# Additional configuration
system.config['box']['distance'] = 2.0
system.config['simulation']['production_time_ns'] = 1000

# Build
system.build()
```

## Troubleshooting Build Issues

### Issue: "Force field not found"

```bash
# Check available force fields
prism --list-forcefields

# Use a known available force field
prism protein.pdb ligand.mol2 -o output --forcefield amber99sb
```

### Issue: "Cannot generate ligand parameters"

```python
# Check ligand structure
from rdkit import Chem

mol = Chem.MolFromMol2File("ligand.mol2")
if mol is None:
    print("Invalid MOL2 file")
else:
    # Check for common issues
    print(f"Atoms: {mol.GetNumAtoms()}")
    print(f"Bonds: {mol.GetNumBonds()}")
    
    # Visualize to check structure
    from rdkit.Chem import Draw
    Draw.MolToFile(mol, "ligand.png")
```

### Issue: "System has non-zero charge"

```bash
# Check topology
grep "qtot" output/GMX_PROLIG_MD/topol.top | tail -1

# If non-zero, check ion addition
# Increase salt concentration or check ligand charge
prism protein.pdb ligand.mol2 -o output --ligand-charge -1
```

### Issue: "Protein-ligand overlap"

This usually means the ligand wasn't properly positioned:

1. Check input files are in same coordinate frame
2. Ensure ligand is positioned where it should bind
3. Use docking to position ligand first if needed

## Performance Tips

### Speed Up Building

1. **Reuse force field parameters**:
```python
# Build once
system = prism.PRISMSystem("protein.pdb", "ligand.mol2", "output")
system.build()

# Reuse for similar systems
# Copy output/LIG.amb2gmx for other proteins
```

2. **Parallel processing**:
```python
from multiprocessing import Pool

def build_system(args):
    protein, ligand, output = args
    system = prism.PRISMSystem(protein, ligand, output)
    return system.build()

# Build multiple systems in parallel
with Pool(4) as pool:
    systems = [
        ("protein.pdb", "lig1.mol2", "out1"),
        ("protein.pdb", "lig2.mol2", "out2"),
        # ...
    ]
    pool.map(build_system, systems)
```

3. **Skip unnecessary steps**:
```bash
# If protein is already clean
prism protein_clean.pdb ligand.mol2 -o output

# If force field files exist
prism protein.pdb ligand.mol2 -o output  # Detects existing files
```

## Quality Assurance

### Pre-Build Checklist

- [ ] Protein structure is complete
- [ ] Ligand has correct protonation
- [ ] Ligand is positioned appropriately
- [ ] Force field choice is appropriate
- [ ] Box size is sufficient

### Post-Build Verification

- [ ] System has zero net charge
- [ ] No atom overlaps
- [ ] Correct number of waters/ions
- [ ] Topology file is complete
- [ ] MDP files are generated

### Visual Inspection

Always visualize the built system:

```bash
# Using VMD
vmd output/GMX_PROLIG_MD/solv_ions.gro

# Using PyMOL
pymol output/GMX_PROLIG_MD/solv_ions.gro
```

## Next Steps

- Choose appropriate [Force Fields](force-fields.md)
- Learn to [Run Simulations](running-simulations.md)
- Understand [Output Files](output-files.md)