# PRISM Modeling Module

The PRISM modeling module provides advanced molecular modeling capabilities for protein-ligand system preparation and structure optimization.

## Features

- **Homology Modeling**: Build protein structures from sequences using template-based methods
- **Molecular Docking**: Dock ligands into protein binding sites with flexible protocols  
- **Conformer Generation**: Generate and optimize multiple conformers for flexible ligands
- **Structure Optimization**: Energy minimization and geometry optimization for proteins and ligands

## Usage Examples

### Homology Modeling

```python
import prism as pm

# Create homology model from sequence and template
modeler = pm.modeling.create_model(
    target_sequence="MKLLIV...",
    template_pdb="template.pdb"
)

# Build the model
model_pdb = modeler.build_homology_model()
```

### Molecular Docking

```python
import prism as pm

# Dock ligand into protein
docker = pm.modeling.dock_ligand(
    protein_pdb="protein.pdb",
    ligand_smiles="CC(=O)NC1=CC=C(C=C1)O",
    binding_site={'center': [10, 15, 20], 'radius': 15}
)

# Perform docking
results = docker.dock()
best_pose = docker.get_best_pose()
```

### Conformer Generation  

```python
import prism as pm

# Generate conformers for a ligand
conformer_gen = pm.modeling.generate_conformers(
    smiles="CC(=O)NC1=CC=C(C=C1)O",
    num_conformers=100
)

# Generate and optimize conformers
conformers_file = conformer_gen.generate_conformers()
best_conformer = conformer_gen.select_best_conformer()
```

### Structure Optimization

```python
import prism as pm

# Optimize protein or ligand structure
optimizer = pm.modeling.optimize_structure(
    input_file="protein.pdb",
    force_field="amber99sb-ildn"
)

# Perform optimization
optimized_structure = optimizer.optimize(method="energy_minimization")
validation = optimizer.validate_optimized_structure()
```

## Integration with PRISM Workflow

The modeling module integrates seamlessly with the main PRISM workflow:

```python
import prism as pm

# 1. Create homology model
modeler = pm.modeling.create_model("target_seq.fasta", "template.pdb")
protein_pdb = modeler.build_homology_model()

# 2. Generate ligand conformers  
conformer_gen = pm.modeling.generate_conformers("CC(=O)NC1=CC=C(C=C1)O")
ligand_sdf = conformer_gen.select_best_conformer()

# 3. Build MD system with modeled structures
system = pm.system(protein_pdb, ligand_sdf)
system.build()

# 4. Run simulation
sim = pm.model(system.get_output_files()['gromacs_directory'])
sim.run(engine="gmx")
```

## Dependencies

- **RDKit**: Required for ligand conformer generation and SMILES processing
- **BioPython**: Optional, for enhanced sequence analysis and PDB processing  
- **OpenMM**: Optional, for advanced structure optimization
- **GROMACS**: For structure optimization and energy minimization

## Module Structure

- `base.py`: Abstract base classes and common functionality
- `homology.py`: Homology modeling implementation
- `docking.py`: Molecular docking engine  
- `conformer.py`: Conformer generation and optimization
- `structure.py`: Structure optimization and refinement