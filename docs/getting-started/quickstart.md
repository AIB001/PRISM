# Quick Start Guide

This guide will help you get started with PRISM in under 5 minutes. We'll walk through the basic workflow of preparing a protein-ligand system for molecular dynamics simulation.

## Prerequisites Check

Before starting, ensure you have:
- ✅ PRISM installed and configured
- ✅ GROMACS properly set up
- ✅ Required force field dependencies (GAFF or OpenFF)
- ✅ Your protein structure file (PDB format)
- ✅ Your ligand structure file (MOL2 or SDF format)

## Basic Usage

### 1. Prepare Your Input Files

Ensure you have:
- **Protein file**: `protein.pdb` - Your protein structure
- **Ligand file**: `ligand.mol2` (for GAFF) or `ligand.sdf` (for OpenFF)

!!! tip "File Format Tips"
    - PDB files should have complete hydrogen atoms
    - MOL2 files work best with GAFF force field
    - SDF files are recommended for OpenFF force field
    - Ensure ligand has proper 3D coordinates and hydrogens

### 2. Run PRISM

#### Option A: Using GAFF Force Field (Default)

The simplest way to run PRISM with GAFF force field:

```bash
prism protein.pdb ligand.mol2 -o my_simulation
```

This command will:
- Process your protein with the AMBER force field
- Parameterize your ligand with GAFF
- Build a solvated system
- Generate all necessary simulation files

#### Option B: Using OpenFF Force Field

For more accurate ligand parameters using OpenFF:

```bash
prism protein.pdb ligand.sdf -o my_simulation --ligand-forcefield openff
```

#### Option C: Custom Configuration

For advanced users with specific requirements:

```bash
prism protein.pdb ligand.mol2 -o my_simulation --config my_config.yaml
```

### 3. Understanding the Output

After PRISM completes, you'll find the following structure in your output directory:

```
my_simulation/
├── LIG.amb2gmx/           # Ligand force field files (if using GAFF)
│   ├── LIG.gro           # Ligand coordinates
│   ├── LIG.itp           # Ligand topology
│   └── posre_LIG.itp     # Position restraints
├── LIG.openff2gmx/        # Ligand force field files (if using OpenFF)
├── GMX_PROLIG_MD/         # Main simulation directory
│   ├── solv_ions.gro     # Solvated system with ions
│   └── topol.top         # System topology
└── mdps/                  # Simulation protocols
    ├── em.mdp            # Energy minimization
    ├── nvt.mdp           # NVT equilibration
    ├── npt.mdp           # NPT equilibration
    └── md.mdp            # Production MD
```

## Running MD Simulations

### 1. Navigate to Simulation Directory

```bash
cd my_simulation/GMX_PROLIG_MD
```

### 2. Create Run Script

Create a file named `run_simulation.sh` with the following content:

```bash
#!/bin/bash

# PRISM Molecular Dynamics Simulation Script
# This script runs a complete MD simulation protocol

echo "========================================="
echo "Starting PRISM MD Simulation Protocol"
echo "========================================="

# Set simulation parameters
NTMPI=1        # Number of thread-MPI ranks
NTOMP=10       # Number of OpenMP threads per rank
GPU_ID=0       # GPU device ID

# Step 1: Energy Minimization
echo -e "\n[1/4] Running Energy Minimization..."
mkdir -p em
if [ -f ./em/em.gro ]; then
    echo "✓ EM already completed, skipping..."
else
    gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -r solv_ions.gro \
               -p topol.top -o ./em/em.tpr -maxwarn 999
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em \
              -ntmpi $NTMPI -ntomp $NTOMP -gpu_id $GPU_ID -v
    echo "✓ Energy minimization completed!"
fi

# Step 2: NVT Equilibration (Temperature)
echo -e "\n[2/4] Running NVT Equilibration..."
mkdir -p nvt
if [ -f ./nvt/nvt.gro ]; then
    echo "✓ NVT already completed, skipping..."
else
    gmx grompp -f ../mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro \
               -p topol.top -o ./nvt/nvt.tpr -maxwarn 999
    gmx mdrun -ntmpi $NTMPI -ntomp $NTOMP -nb gpu -bonded gpu \
              -pme gpu -gpu_id $GPU_ID -s ./nvt/nvt.tpr \
              -deffnm ./nvt/nvt -v
    echo "✓ NVT equilibration completed!"
fi

# Step 3: NPT Equilibration (Pressure)
echo -e "\n[3/4] Running NPT Equilibration..."
mkdir -p npt
if [ -f ./npt/npt.gro ]; then
    echo "✓ NPT already completed, skipping..."
else
    gmx grompp -f ../mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro \
               -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn 999
    gmx mdrun -ntmpi $NTMPI -ntomp $NTOMP -nb gpu -bonded gpu \
              -pme gpu -gpu_id $GPU_ID -s ./npt/npt.tpr \
              -deffnm ./npt/npt -v
    echo "✓ NPT equilibration completed!"
fi

# Step 4: Production MD
echo -e "\n[4/4] Running Production MD..."
mkdir -p prod
if [ -f ./prod/md.gro ]; then
    echo "✓ Production MD already completed!"
else
    gmx grompp -f ../mdps/md.mdp -c ./npt/npt.gro -r ./npt/npt.gro \
               -p topol.top -o ./prod/md.tpr -maxwarn 999
    gmx mdrun -ntmpi $NTMPI -ntomp $NTOMP -nb gpu -bonded gpu \
              -pme gpu -gpu_id $GPU_ID -s ./prod/md.tpr \
              -deffnm ./prod/md -v
    echo "✓ Production MD completed!"
fi

echo -e "\n========================================="
echo "Simulation Protocol Completed Successfully!"
echo "========================================="
echo "Results are in the 'prod' directory"
```

### 3. Run the Simulation

Make the script executable and run:

```bash
chmod +x run_simulation.sh
bash ./run_simulation.sh
```

!!! info "GPU Acceleration"
    The script is configured for GPU acceleration. If you don't have a GPU, remove the GPU-related flags:
    - Remove: `-nb gpu -bonded gpu -pme gpu -gpu_id 0`
    - The simulation will run on CPU only (slower but still functional)

### 4. Monitor Progress

The simulation will run through four stages:
1. **Energy Minimization** (~5 minutes)
2. **NVT Equilibration** (~10 minutes)
3. **NPT Equilibration** (~10 minutes)
4. **Production MD** (varies based on settings)

## Quick Examples

### Example 1: Basic Protein-Ligand System

```bash
# Prepare system with default settings
prism 1ABC.pdb drug.mol2 -o drug_simulation

# Run simulations
cd drug_simulation/GMX_PROLIG_MD
bash ../run_simulation.sh
```

### Example 2: Using Custom Water Model

```bash
# Create custom config file
cat > custom.yaml << EOF
water_model: tip4p
box_size: 1.0
EOF

# Run with custom config
prism protein.pdb ligand.mol2 -o output --config custom.yaml
```

### Example 3: High-Precision OpenFF Setup

```bash
# Use OpenFF for better ligand parameters
prism receptor.pdb compound.sdf -o precise_sim --ligand-forcefield openff

# Run with extended equilibration
cd precise_sim/GMX_PROLIG_MD
# Edit mdps/nvt.mdp and mdps/npt.mdp to increase nsteps if needed
bash ../run_simulation.sh
```

## Analyzing Results

After the simulation completes, you can analyze the trajectory:

```bash
# Check system stability (RMSD)
gmx rms -s ./em/em.tpr -f ./prod/md.xtc -o rmsd.xvg

# Analyze protein-ligand interactions
gmx distance -s ./prod/md.tpr -f ./prod/md.xtc -select 'com of group "Protein" plus com of group "LIG"' -o distance.xvg

# Extract snapshots
gmx trjconv -s ./prod/md.tpr -f ./prod/md.xtc -o trajectory.pdb -dt 100
```

## Common Options

| Option | Description | Example |
|--------|-------------|---------|
| `-o` | Output directory | `-o my_project` |
| `--ligand-forcefield` | Force field for ligand | `--ligand-forcefield openff` |
| `--config` | Custom configuration file | `--config custom.yaml` |
| `--water` | Water model | `--water tip4p` |
| `--box-size` | Box padding (nm) | `--box-size 1.2` |
| `--ions` | Ion concentration (M) | `--ions 0.15` |

## Troubleshooting Quick Tips

!!! warning "Common Issues"

    **Issue**: "Atom type not found" error
    - **Solution**: Ensure force field files are properly generated. Check the `LIG.amb2gmx/` or `LIG.openff2gmx/` directory.

    **Issue**: Simulation crashes during minimization
    - **Solution**: Check for clashes in input structures. Consider using a larger box size.

    **Issue**: "GPU not detected" message
    - **Solution**: Verify CUDA installation and GROMACS GPU support with `gmx mdrun -version`

    **Issue**: Simulation is very slow
    - **Solution**: Use GPU acceleration or increase the number of CPU cores with `-ntomp`

## Next Steps

Now that you've completed your first simulation:

1. **Explore Advanced Features**: Check the [User Guide](../user-guide/index.md) for detailed configuration options
2. **Optimize Performance**: Learn about [performance tuning](../user-guide/performance.md)
3. **Analyze Results**: See our [analysis tutorials](../tutorials/analysis.md)
4. **Customize Protocols**: Modify MDP files for your specific research needs

---

!!! success "Congratulations!"
    You've successfully set up and run your first protein-ligand MD simulation with PRISM! 🎉
    
    For more complex scenarios and advanced features, explore the rest of our documentation.