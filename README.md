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
   # Ubuntu/Debian
   sudo apt-get install gromacs
   
   # Using conda
   conda install -c bioconda gromacs
   ```

2. **Python 3.8+** with required packages:

   ```bash
   pip install pyyaml numpy
   ```

3. **PDBFixer** (required)

   ```bash
   conda install -c conda-forge pdbfixer
   ```

### Force Field Specific Dependencies

#### For GAFF Support:

```bash
# AmberTools (required)
conda install -c conda-forge ambertools

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

# Energy minimization
gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

# NVT equilibration
gmx grompp -f ../mdps/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT equilibration
gmx grompp -f ../mdps/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# Production MD
gmx grompp -f ../mdps/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

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
  author       = {Institute of Quantitative Biology, Zhejiang University},
  title        = {PRISM: An Integrated Framework for High-Throughput Protein-Ligand Simulation Setup and Molecular Simulation-Based Drug Virtual Screening},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/aib001/PRISM}
}
```

## License

PRISM is released under the MIT License. Force field parameters are subject to their respective licenses.

