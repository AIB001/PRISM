# PRISM - Protein Receptor Interaction Simulation Modeler

PRISM is a comprehensive tool for building protein-ligand systems for molecular dynamics simulations in GROMACS. It supports multiple force fields for ligands including GAFF (General AMBER Force Field) and OpenFF (Open Force Field).

## Features

- Multiple Force Field Support
  - GAFF via AmberTools
  - OpenFF via openff-toolkit
- **Automatic System Building**: Complete workflow from PDB/MOL2/SDF to simulation-ready files
- **Flexible Configuration**: YAML-based configuration for easy customization
- **Smart File Processing**: Handles various input formats with automatic conversion
- **Position Restraints**: Automatic generation for equilibration protocols
- **Complete MDP Files**: Pre-configured protocols for minimization, equilibration, and production

## Installation

### Prerequisites

1. **GROMACS** (required)

   ```bash
   # Ubuntu/Debian, CUDA-toolkit is needed, Dependence installation
   nvcc --version # For check
   sudo apt-get install gcc
   sudo apt-get install g++
   sudo apt-get install cmake
   
   # Download GROMACS Package
   wget https://ftp.gromacs.org/gromacs/gromacs-2024.3.tar.gz
   tar xfz gromacs-2024.3.tar.gz
   
   # Prepare to build GROMACS with cmake
   cd gromacs-2024.3
   mkdir build
   cd build
   
   # Installation
   cmake .. -DGMX_MPI=ON \
   -DGMX_BUILD_OWN_FFTW=ON \
   -DGMX_GPU=CUDA \
   -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
   -DCUDA_INCLUDE_DIRS=/usr/local/cuda/include \
   -DCUDA_CUDART_LIBRARY=/usr/local/cuda/lib64 \
   -DCMAKE_INSTALL_PREFIX=/mnt/data/zf/gromacs-2024.3
   
   make -j${nproc}
   make check # Optional but recommended
   sudo make install
   source /mnt/data/zf/gromacs-2024.3/bin/GMXRC
   
   # install bioconda
   conda install -c bioconda
   ```

2. **Python 3.10** with required packages:

   ```bash
   # create environment 
   conda create -n prism python=3.10
   conda activate prism
   ```

3. **PDBFixer and pyyaml** (required)

   ```bash
   conda install -c conda-forge pdbfixer numpy scipy
   pip install pyyaml
   ```

### Force Field Specific Dependencies

#### For GAFF Support:

```bash
# AmberTools (required)
conda install conda-forge::ambertools
# ACPYPE (required)
pip install acpype
# Optional but recommended
conda install -c conda-forge rdkit
```

#### For OpenFF Support:

```bash
# OpenFF toolkit and dependencies
conda install -c conda-forge openff-toolkit openff-interchange

# RDKit (required for SDF handling)
conda install -c conda-forge rdkit

# OpenBabel (Optional bur Recommended)
conda install conda-forge::openbabel
```

### Installing PRISM

1. Clone or download the PRISM package

2. Install in development mode:

   ```bash
   cd PRISM
   pip install -e .
   ```

Or use directly without installation:

```bash
python /path/to/PRISM/prism/builder.py
```

## Quick Start

### Basic Usage

1. **Using GAFF (default)**:

   ```bash
   prism protein.pdb ligand.mol2 -o output_dir
   ```

2. **Using OpenFF**:

   ```bash
   prism protein.pdb ligand.sdf -o output_dir --ligand-forcefield openff
   ```

3. **With custom configuration**:

   ```bash
   prism protein.pdb ligand.mol2 -o output_dir --config my_config.yaml
   ```

### Running MD Simulations

After PRISM completes, you can run the simulations:

```bash
cd output_dir/GMX_PROLIG_MD
```

Create and save `localrun.sh`:

```bash
#!/bin/bash

######################################################
# SIMULATION PART
######################################################

# Energy Minimization (EM)
mkdir -p em
if [ -f ./em/em.gro ]; then
    echo "EM already completed, skipping..."
elif [ -f ./em/em.tpr ]; then
    echo "EM tpr file found, continuing from checkpoint..."
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v -cpi ./em/em.cpt
else
    echo "Starting EM from scratch..."
    gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o ./em/em.tpr -maxwarn 999
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v
fi

# NVT Equilibration
mkdir -p nvt
if [ -f ./nvt/nvt.gro ]; then
    echo "NVT already completed, skipping..."
elif [ -f ./nvt/nvt.tpr ]; then
    echo "NVT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v -cpi ./nvt/nvt.cpt
else
    echo "Starting NVT from scratch..."
    gmx grompp -f ../mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro -p topol.top -o ./nvt/nvt.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v
fi

# NPT Equilibration
mkdir -p npt
if [ -f ./npt/npt.gro ]; then
    echo "NPT already completed, skipping..."
elif [ -f ./npt/npt.tpr ]; then
    echo "NPT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v -cpi ./npt/npt.cpt
else
    echo "Starting NPT from scratch..."
    gmx grompp -f ../mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v
fi

# Production MD
mkdir -p prod
if [ -f ./prod/md.gro ]; then
    echo "Production MD already completed, skipping..."
elif [ -f ./prod/md.tpr ]; then
    echo "Production MD tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v -cpi ./prod/md.cpt
else
    echo "Starting Production MD from scratch..."
    gmx grompp -f ../mdps/md.mdp -c ./npt/npt.gro -r ./npt/npt.gro -p topol.top -o ./prod/md.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v
fi
```

run `bash localrun.sh`

## Configuration

PRISM uses YAML configuration files for customization. Key parameters include:

- **Force fields**: Choose from AMBER, CHARMM, or OPLS variants
- **Water models**: TIP3P, TIP4P, SPC, SPCE
- **Simulation parameters**: Temperature, pressure, time, etc.
- **Box settings**: Size, shape, solvation
- **Output controls**: Trajectory frequency, compression

See `configs/default.yaml` for a complete example.

## File Structure

```
PRISM/
├── prism/                    # Core modules
│   ├── __init__.py
│   ├── builder.py           # Main program
│   ├── forcefield/          # Force field generators
│   │   ├── __init__.py
│   │   ├── base.py         # Base class
│   │   ├── gaff.py         # GAFF wrapper
│   │   └── openff.py       # OpenFF wrapper
│   └── utils/              # Utilities
├── configs/                 # Example configurations
├── examples/               # Example input files
├── docs/                   # Documentation
└── README.md              # This file
```

## Output Files

PRISM generates a complete set of files ready for MD simulation:

- **Force field files** (`LIG.amb2gmx/` or `LIG.openff2gmx/`):
  - `LIG.gro`: Ligand coordinates
  - `LIG.itp`: Ligand topology
  - `atomtypes_LIG.itp`: Atom type definitions
  - `posre_LIG.itp`: Position restraints
- **System files** (`GMX_PROLIG_MD/`):
  - `solv_ions.gro`: Complete solvated system
  - `topol.top`: System topology
- **Protocol files** (`mdps/`):
  - `em.mdp`: Energy minimization
  - `nvt.mdp`: NVT equilibration
  - `npt.mdp`: NPT equilibration
  - `md.mdp`: Production run

## Troubleshooting

### Common Issues

1. **"Command not found" errors**: Ensure all dependencies are installed and in PATH
2. **Force field errors**: Check that AmberTools (for GAFF) or openff-toolkit (for OpenFF) is installed
3. **Memory errors**: Large systems may require more RAM, especially during parameterization

### Getting Help

- Check the log files in the output directory
- Ensure input files are properly formatted
- Verify all dependencies are correctly installed

## Citation

If you use PRISM in your research, please cite:

```LaTeX
@software{PRISM,
  author       = {Institute of Quantitative Biology, Zhejiang University; Theoretical Chemistry Institute, University of Wisconsin-Madison},
  title        = {PRISM: An Integrated Framework for High-Throughput Protein-Ligand Simulation Setup and Molecular Simulation-Based Drug Virtual Screening},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/aib001/PRISM}
}
```

## License

PRISM is released under the MIT License. Force field parameters are subject to their respective licenses.