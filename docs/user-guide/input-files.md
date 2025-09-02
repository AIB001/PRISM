# Input Files Guide

PRISM requires two primary input files: a protein structure and a ligand structure. This guide covers file formats, preparation methods, and best practices for preparing high-quality input files.

## Required Input Files

### Protein Structure (PDB)

PRISM expects protein structures in PDB format:

```bash
# Basic usage
prism protein.pdb ligand.mol2 -o output

# Alternative specification
prism --protein-file protein.pdb --ligand-file ligand.mol2 -o output
```

**PDB Requirements:**
- Standard PDB format (`.pdb` extension)
- Should contain only protein atoms (ATOM records)
- HETATM records are automatically removed
- Missing atoms/residues should be fixed beforehand

### Ligand Structure (MOL2/SDF)

PRISM supports two ligand formats:

- **MOL2** (`.mol2`): Tripos MOL2 format with charges
- **SDF** (`.sdf`, `.sd`): Structure Data File format

```bash
# Using MOL2 with GAFF
prism protein.pdb ligand.mol2 -o output

# Using SDF with OpenFF
prism protein.pdb ligand.sdf -o output --ligand-forcefield openff
```

## Protein Preparation

### 1. Download Structure

#### From PDB Database

```bash
# Download directly from PDB
wget https://files.rcsb.org/download/1ABC.pdb

# Or using PRISM utility (if available)
python -c "
from prism.utils import download_pdb
download_pdb('1ABC', 'protein.pdb')
"
```

#### From AlphaFold

```bash
# Download AlphaFold structure
wget https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.pdb
```

### 2. Clean the Structure

PRISM automatically cleans basic issues, but pre-cleaning is recommended:

#### Remove Non-Protein Atoms

```bash
# Remove water, ions, and ligands
grep "^ATOM" input.pdb > protein_clean.pdb

# Or using pymol
pymol -c -d "
load input.pdb;
remove resn HOH+WAT+NA+CL+SO4;
remove hetatm and not resn MSE;
save protein_clean.pdb;
quit
"
```

#### Using PDBFixer (Recommended)

```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Load and fix structure
fixer = PDBFixer(filename='input.pdb')

# Find and add missing residues
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Add missing hydrogens at pH 7.0
fixer.addMissingHydrogens(7.0)

# Save cleaned structure
with open('protein_clean.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
```

### 3. Handle Special Cases

#### Multiple Chains

```python
# Keep only chain A
from Bio import PDB

parser = PDB.PDBParser()
structure = parser.get_structure('protein', 'input.pdb')

# Extract chain A
chain_a = structure[0]['A']

# Save chain A
io = PDB.PDBIO()
io.set_structure(chain_a)
io.save('protein_chainA.pdb')
```

#### Mutations

```bash
# Using pymol for mutations
pymol -c -d "
load protein.pdb;
wizard mutagenesis;
mutation A/123/ALA;
save protein_mutant.pdb;
quit
"
```

#### Missing Loops

For missing loops, use modeling tools before PRISM:

```bash
# Using MODELLER (requires license)
python model_loops.py

# Or use online servers:
# - SWISS-MODEL: https://swissmodel.expasy.org
# - ModLoop: https://modbase.compbio.ucsf.edu/modloop/
```

## Ligand Preparation

### 1. MOL2 Format Preparation

#### From SMILES

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Create molecule from SMILES
smiles = "CC(=O)Nc1ccc(cc1)O"  # Paracetamol
mol = Chem.MolFromSmiles(smiles)

# Add hydrogens
mol = Chem.AddHs(mol)

# Generate 3D coordinates
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)

# Save as SDF (then convert to MOL2)
writer = Chem.SDWriter('ligand.sdf')
writer.write(mol)
writer.close()

# Convert to MOL2 using OpenBabel
import subprocess
subprocess.run(['obabel', '-isdf', 'ligand.sdf', 
                '-omol2', '-O', 'ligand.mol2'])
```

#### From PubChem

```bash
# Download SDF from PubChem (e.g., aspirin CID 2244)
wget https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF -O aspirin.sdf

# Convert to MOL2
obabel -isdf aspirin.sdf -omol2 -O aspirin.mol2
```

#### From PDB Ligand

```python
# Extract ligand from PDB complex
from pymol import cmd

cmd.load('complex.pdb')
cmd.select('ligand', 'resn LIG')  # Replace LIG with actual residue name
cmd.save('ligand.mol2', 'ligand')
```

### 2. SDF Format Preparation

#### Using RDKit

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Load molecule
mol = Chem.MolFromSmiles('CC(C)CC1=CC=CC=C1')

# Prepare 3D structure
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol)

# Save as SDF
writer = Chem.SDWriter('ligand.sdf')
writer.write(mol)
writer.close()
```

#### With Charges (for OpenFF)

```python
# OpenFF can calculate charges automatically
# But you can pre-calculate if needed

from openff.toolkit.topology import Molecule

# Load molecule
mol = Molecule.from_smiles('CC(=O)Nc1ccc(cc1)O')

# Generate conformer
mol.generate_conformers(n_conformers=1)

# Assign AM1-BCC charges
mol.assign_partial_charges('am1bcc')

# Save with charges
mol.to_file('ligand_charged.sdf', file_format='SDF')
```

### 3. Charge Assignment

#### For GAFF (MOL2)

PRISM automatically calculates AM1-BCC charges, but you can pre-calculate:

```bash
# Using antechamber
antechamber -i ligand.mol2 -fi mol2 \
            -o ligand_charged.mol2 -fo mol2 \
            -c bcc -nc 0  # nc = net charge
```

#### For OpenFF (SDF)

OpenFF handles charges automatically, but manual control is possible:

```python
from openff.toolkit.topology import Molecule

mol = Molecule.from_file('ligand.sdf')
mol.assign_partial_charges('am1bccelf10')  # Latest charge model
mol.to_file('ligand_charged.sdf')
```

## File Format Specifications

### PDB Format Requirements

```
ATOM      1  N   MET A   1      -8.901   4.127  -0.321  1.00 10.00           N
ATOM      2  CA  MET A   1      -8.608   4.135   1.107  1.00 10.00           C
```

**Column definitions:**
- 1-6: Record name (ATOM)
- 7-11: Atom serial number
- 13-16: Atom name
- 17: Alternate location
- 18-20: Residue name
- 22: Chain identifier
- 23-26: Residue sequence number
- 31-38: X coordinate
- 39-46: Y coordinate
- 47-54: Z coordinate
- 55-60: Occupancy
- 61-66: Temperature factor
- 77-78: Element symbol

### MOL2 Format Structure

```
@<TRIPOS>MOLECULE
ligand_name
   24    25     1     0     0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1         -0.7560    0.5450    0.0000 C.3     1 LIG        -0.0653
      2 C2          0.7560    0.5450    0.0000 C.ar    1 LIG        -0.0177

@<TRIPOS>BOND
     1     1     2 1
     2     2     3 ar
```

### SDF Format Structure

```
ligand_name
  RDKit          3D

 24 25  0  0  0  0  0  0  0  0999 V2000
   -0.7560    0.5450    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7560    0.5450    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
...
  1  2  1  0
  2  3  2  0
M  END
$$$$
```

## Quality Checks

### Protein Structure Validation

```python
# Check for missing atoms
from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser()
structure = parser.get_structure('protein', 'protein.pdb')

# Check for gaps in sequence
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    seq = pp.get_sequence()
    print(f"Chain sequence: {seq}")
    
# Check for clashes
from Bio.PDB import NeighborSearch

atoms = list(structure.get_atoms())
ns = NeighborSearch(atoms)
clashes = ns.search_all(1.5)  # Find atoms closer than 1.5 Å
print(f"Found {len(clashes)} clashes")
```

### Ligand Structure Validation

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

# Load ligand
mol = Chem.MolFromMol2File('ligand.mol2')

# Basic checks
print(f"Molecular weight: {Descriptors.MolWt(mol):.2f}")
print(f"LogP: {Descriptors.MolLogP(mol):.2f}")
print(f"Rotatable bonds: {Descriptors.NumRotatableBonds(mol)}")
print(f"H-bond donors: {Descriptors.NumHDonors(mol)}")
print(f"H-bond acceptors: {Descriptors.NumHAcceptors(mol)}")

# Check for problems
problems = Chem.DetectChemistryProblems(mol)
if problems:
    print("Chemistry problems detected:")
    for problem in problems:
        print(f"  - {problem.GetType()}: {problem.Message()}")
```

## Common Input File Issues

### Problem: Missing Residues

**Symptoms**: Gaps in protein structure
**Solution**: Use modeling tools or PDBFixer

```python
# Using PDBFixer
from pdbfixer import PDBFixer

fixer = PDBFixer(filename='protein.pdb')
fixer.findMissingResidues()
print(f"Missing residues: {fixer.missingResidues}")

# Add missing residues
fixer.findMissingAtoms()
fixer.addMissingAtoms()
```

### Problem: Incorrect Protonation

**Symptoms**: Wrong charges, simulation crashes
**Solution**: Set correct pH during preparation

```bash
# For protein (using pdb2pqr)
pdb2pqr30 --ph=7.0 --pdb-output protein_ph7.pdb protein.pdb

# For ligand (using antechamber)
antechamber -i ligand.mol2 -fi mol2 -o ligand_charged.mol2 -fo mol2 -c bcc
```

### Problem: Non-standard Residues

**Symptoms**: PRISM doesn't recognize residue
**Solution**: Convert or parameterize separately

```python
# Convert selenomethionine to methionine
sed 's/MSE/MET/g' protein.pdb > protein_fixed.pdb
```

### Problem: Alternative Conformations

**Symptoms**: Multiple positions for same atom
**Solution**: Keep only one conformation

```bash
# Remove alternative conformations (keep A)
grep -v '^ATOM.*[B-Z] [A-Z][A-Z][A-Z]' protein.pdb > protein_fixed.pdb
```

## Advanced Input Preparation

### Multi-Component Systems

For systems with multiple ligands or cofactors:

```python
# Prepare each component separately
prism protein.pdb ligand1.mol2 -o temp1 --ligand-forcefield gaff
prism protein.pdb ligand2.mol2 -o temp2 --ligand-forcefield gaff

# Then combine topologies manually
# (Advanced - see Advanced Usage guide)
```

### Covalent Ligands

For covalently bound ligands:

1. Keep ligand with protein in PDB
2. Define custom parameters
3. Use specialized tools (e.g., CGenFF for CHARMM)

### Metal Centers

For metalloproteins:

```python
# Keep metal ions in protein PDB
# PRISM will handle standard ions (Zn2+, Mg2+, etc.)
# For complex metal centers, use specialized force fields
```

## Best Practices

1. **Always visualize**: Check structures in PyMOL/ChimeraX before using
2. **Verify chemistry**: Ensure correct bond orders and charges
3. **Check stereochemistry**: Verify chirality is correct
4. **Document preparation**: Save all preparation steps for reproducibility
5. **Test with small systems**: Validate protocol with simple cases first

## Quality Control Checklist

Before running PRISM, verify:

- [ ] Protein has all residues (no gaps)
- [ ] All atoms present (no missing side chains)
- [ ] Correct protonation state
- [ ] No clashes (<1.5 Å between atoms)
- [ ] Ligand has correct charge
- [ ] Ligand has reasonable geometry
- [ ] File formats are correct
- [ ] Files are in the same directory or paths are correct

## Next Steps

- Learn about [Building Systems](building-systems.md)
- Understand [Force Fields](force-fields.md) selection
- Configure your [simulation parameters](configuration.md)