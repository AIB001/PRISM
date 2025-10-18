# PRISM - Protein Receptor Interaction Simulation Modeler

PRISM is a comprehensive tool for building protein-ligand systems for molecular dynamics simulations in GROMACS. It supports multiple force fields for ligands including GAFF, GAFF2, OpenFF, CGenFF, OPLS-AA, and SwissParam (MMFF/MATCH/Hybrid).

## Features

- **Multiple Force Field Support**
  - **GAFF/GAFF2** via AmberTools
  - **OpenFF** via openff-toolkit
  - **CGenFF** (CHARMM General Force Field) via web download
  - **OPLS-AA** via LigParGen server
  - **SwissParam** (MMFF/MATCH/Hybrid) via SwissParam server
- **Automatic System Building**: Complete workflow from PDB/MOL2/SDF to simulation-ready files
- **Flexible Configuration**: YAML-based configuration for easy customization
- **Smart File Processing**: Handles various input formats with automatic conversion
- **Position Restraints**: Automatic generation for equilibration protocols
- **Complete MDP Files**: Pre-configured protocols for minimization, equilibration, and production
- **Metal Ion Support**: Automatic recognition and handling of metal ions (Zn, Ca, Mg, Fe, etc.)
- **Protonation States**: Support for special residue protonation states (CYM, HID, HIE, HIP, CYX, LYN, ASH, GLH)

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
   wget https://ftp.gromacs.org/gromacs/gromacs-2025.1.tar.gz
   tar xfz gromacs-2025.1.tar.gz
   
   # Prepare to build GROMACS with cmake
   cd gromacs-2025.1
   mkdir build
   cd build
   
   # Installation
   cmake .. -DGMX_MPI=ON \
   -DGMX_BUILD_OWN_FFTW=ON \
   -DGMX_GPU=CUDA \
   -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
   -DCUDA_INCLUDE_DIRS=/usr/local/cuda/include \
   -DCUDA_CUDART_LIBRARY=/usr/local/cuda/lib64 \
   -DCMAKE_INSTALL_PREFIX=~/gromacs-2025.1
   
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

#### For GAFF/GAFF2 Support:

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

# OpenBabel (Optional but Recommended)
conda install conda-forge::openbabel
```

#### For CGenFF Support:

CGenFF requires downloading force field files from the web server:

1. Visit https://cgenff.umaryland.edu/
2. Upload your ligand structure (MOL2/SDF)
3. Download the generated files (PDB and TOP files)
4. Place them in a directory
5. Use `--ligand-forcefield cgenff --forcefield-path /path/to/cgenff_files`

**Note**: For halogens (F, Cl, Br, I), PRISM automatically removes lone pair (LP) atoms and transfers their charges to halogen atoms.

#### For OPLS-AA Support:

OPLS-AA uses the LigParGen web server (requires internet connection):

```bash
# Required Python packages
pip install requests
```

Usage: `--ligand-forcefield opls`

**Note**: Internet connection required during force field generation.

#### For SwissParam Support (MMFF/MATCH/Hybrid):

SwissParam uses the SwissParam web server (requires internet connection):

```bash
# Required Python packages
pip install mechanize
```

Usage:
- `--ligand-forcefield mmff` (MMFF-based)
- `--ligand-forcefield match` (MATCH)
- `--ligand-forcefield hybrid` (Hybrid MMFF-based-MATCH)

**Note**: Internet connection required during force field generation.

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

2. **Using GAFF2 (improved GAFF)**:

   ```bash
   prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield gaff2
   ```

3. **Using OpenFF**:

   ```bash
   prism protein.pdb ligand.sdf -o output_dir --ligand-forcefield openff
   ```

4. **Using CGenFF**:

   ```bash
   # First download CGenFF files from https://cgenff.umaryland.edu/
   prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield cgenff --forcefield-path /path/to/cgenff_files
   ```

5. **Using OPLS-AA**:

   ```bash
   prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield opls
   ```

6. **Using SwissParam force fields**:

   ```bash
   # MMFF-based
   prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield mmff

   # MATCH
   prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield match

   # Hybrid MMFF-based-MATCH
   prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield hybrid
   ```

7. **With custom configuration**:

   ```bash
   prism protein.pdb ligand.mol2 -o output_dir --config my_config.yaml
   ```

8. **With specific protein force field**:

   ```bash
   prism protein.pdb ligand.mol2 -o output_dir --forcefield amber14sb --water tip4p
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

## Force Field Selection Guide

### When to Use Each Force Field

| Force Field | Best For | Pros | Cons | Internet Required |
|------------|----------|------|------|-------------------|
| **GAFF** | General small molecules | Widely used, well-tested | Older parameters | No |
| **GAFF2** | General small molecules | Improved over GAFF, better for pharmaceuticals | Newer, less tested | No |
| **OpenFF** | Drug-like molecules | Modern, data-driven parameters | Requires more dependencies | No |
| **CGenFF** | CHARMM compatibility | Consistent with CHARMM protein FF | Manual web download required | For download |
| **OPLS-AA** | All-atom simulations | Good for organic molecules | Internet required | Yes |
| **MMFF** | Quick parameterization | Fast, general purpose | Less accurate for MD | Yes |
| **MATCH** | CHARMM-style parameters | Good transferability | May fail for complex molecules | Yes |
| **Hybrid** | Balanced approach | Combines MMFF and MATCH | Internet required | Yes |

### Special Features

- **Metal Ions**: Automatically recognized and handled (Zn²⁺, Ca²⁺, Mg²⁺, Fe²⁺/³⁺, Cu²⁺, Mn²⁺, etc.)
- **Protonation States**: Support for CYM, HID, HIE, HIP, CYX, LYN, ASH, GLH
- **Halogen Handling**: CGenFF lone pairs (LP) automatically removed and charges transferred

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
│   │   ├── gaff.py         # GAFF force field
│   │   ├── gaff2.py        # GAFF2 force field
│   │   ├── openff.py       # OpenFF force field
│   │   ├── cgenff.py       # CGenFF force field
│   │   ├── opls_aa.py      # OPLS-AA force field
│   │   └── swissparam.py   # SwissParam force fields
│   └── utils/              # Utilities
│       ├── cleaner.py      # Protein cleaning & metal handling
│       ├── system.py       # System building
│       ├── config.py       # Configuration management
│       └── mdp.py          # MDP file generation
├── configs/                 # Example configurations
├── examples/               # Example input files
├── docs/                   # Documentation
└── README.md              # This file
```

## Output Files

PRISM generates a complete set of files ready for MD simulation:

- **Force field files** (directory name depends on force field):
  - `LIG.amb2gmx/` (GAFF/GAFF2)
  - `LIG.openff2gmx/` (OpenFF)
  - `LIG.cgenff2gmx/` (CGenFF)
  - `LIG.opls2gmx/` (OPLS-AA)
  - `LIG.mmff2gmx/`, `LIG.match2gmx/`, `LIG.hybrid2gmx/` (SwissParam)

  Each contains:
  - `LIG.gro`: Ligand coordinates
  - `LIG.itp`: Ligand topology
  - `atomtypes_LIG.itp`: Atom type definitions
  - `posre_LIG.itp`: Position restraints

- **System files** (`GMX_PROLIG_MD/`):
  - `solv_ions.gro`: Complete solvated system
  - `topol.top`: System topology
  - `protein_clean.pdb`: Cleaned protein structure
  - `topol_Protein.itp`: Protein topology
  - `topol_Ion*.itp`: Metal ion topologies (if present)

- **Protocol files** (`mdps/`):
  - `em.mdp`: Energy minimization
  - `nvt.mdp`: NVT equilibration
  - `npt.mdp`: NPT equilibration
  - `md.mdp`: Production run

- **Configuration backup**:
  - `prism_config.yaml`: Complete configuration used for this build

## Troubleshooting

### Common Issues

1. **"Command not found" errors**: Ensure all dependencies are installed and in PATH
   ```bash
   # Check GROMACS
   gmx --version

   # Check AmberTools (for GAFF/GAFF2)
   antechamber -h

   # Check Python packages
   python -c "import openff.toolkit"  # For OpenFF
   python -c "import requests"         # For OPLS-AA
   python -c "import mechanize"        # For SwissParam
   ```

2. **Force field errors**:
   - **GAFF/GAFF2**: Ensure AmberTools and ACPYPE are installed
   - **OpenFF**: Check openff-toolkit and openff-interchange
   - **CGenFF**: Verify you downloaded the correct files from the web server
   - **OPLS-AA**: Check internet connection and requests package
   - **SwissParam**: Check internet connection and mechanize package

3. **SwissParam errors**:
   - `MATCH ERROR: Could not Finish Charging`: Try `--ligand-forcefield mmff` or `--ligand-forcefield hybrid`
   - Tautomer issues: Use a different force field (GAFF2 recommended)

4. **CGenFF halogen issues**:
   - LP (lone pair) atoms are automatically removed
   - If you see LP-related errors, please report as a bug

5. **Metal ion issues**:
   - PRISM automatically detects common metal ions
   - Supported: Zn, Ca, Mg, Fe, Cu, Mn, Co, Ni, K, Na, Cl
   - If a metal is not recognized, it may be treated as a regular atom

6. **Protonation state issues**:
   - Ensure residues use standard names (CYM, HID, HIE, HIP, CYX, LYN, ASH, GLH)
   - GROMACS pdb2gmx should handle these automatically

7. **Memory errors**: Large systems may require more RAM, especially during parameterization

### Getting Help

- Check the log files in the output directory
- Ensure input files are properly formatted (PDB/MOL2/SDF)
- Verify all dependencies are correctly installed
- Check PRISM configuration in `prism_config.yaml`
- For internet-based force fields, verify network connectivity

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
