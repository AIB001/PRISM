# Welcome to PRISM Documentation

<p align="center">
  <img src="assets/logo.png" alt="PRISM Logo" width="300">
</p>

<p align="center">
  <strong>Python-based Research Infrastructure for Structural Modeling</strong><br>
  <em>A comprehensive toolkit for molecular dynamics analysis and visualization</em>
</p>

<p align="center">
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/Python-3.8%2B-blue.svg" alt="Python Version"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License"></a>
  <a href="https://prism-docs.github.io"><img src="https://img.shields.io/badge/Docs-Latest-brightgreen.svg" alt="Documentation"></a>
  <a href="https://pypi.org/project/prism-md/"><img src="https://img.shields.io/pypi/v/prism-md.svg" alt="PyPI version"></a>
</p>

---

## Overview

**PRISM** is a powerful Python framework designed for molecular dynamics (MD) analysis, developed at the Institute of Quantitative Biology, Zhejiang University and Theoretical Chemistry Institute, University of Wisconsin-Madison. It provides researchers with an intuitive and efficient platform for analyzing MD trajectories, computing structural properties, and visualizing molecular systems.

## Quick Start

### Installation

Install PRISM using pip:

```bash
pip install prism-md
```

Or install from source:

```bash
git clone https://github.com/your-username/PRISM.git
cd PRISM
pip install -e .
```

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

## Key Features

<div class="feature-grid">
  <div class="feature-card">
    <h3>📈 Trajectory Analysis</h3>
    <p>Comprehensive tools for analyzing molecular dynamics trajectories including RMSD, RMSF, radius of gyration, and more.</p>
  </div>
  
  <div class="feature-card">
    <h3>🔬 Structural Analysis</h3>
    <p>Calculate radial distribution functions, hydrogen bonds, salt bridges, and other structural properties.</p>
  </div>
  
  <div class="feature-card">
    <h3>🎯 Atom Selection</h3>
    <p>Powerful and flexible atom selection language for precise analysis of molecular systems.</p>
  </div>
  
  <div class="feature-card">
    <h3>📊 Data Visualization</h3>
    <p>Built-in plotting functions for creating publication-quality figures and interactive visualizations.</p>
  </div>
  
  <div class="feature-card">
    <h3>⚙️ Format Support</h3>
    <p>Support for common MD formats including PDB, DCD, XTC, TRR, and more.</p>
  </div>
  
  <div class="feature-card">
    <h3>🔄 Workflow Integration</h3>
    <p>Seamlessly integrate with popular MD packages like GROMACS, AMBER, and NAMD.</p>
  </div>
</div>

## Core Modules

### Model Module 

### Analysis Module
The analysis module provides comprehensive tools for trajectory analysis:

- **Geometric Analysis**: RMSD, RMSF, radius of gyration, end-to-end distance
- **Structural Analysis**: RDF, hydrogen bonds, contacts, secondary structure
- **Dynamic Analysis**: MSD, diffusion coefficients, autocorrelation functions

### Visualization Module
Advanced visualization capabilities:

- **3D Rendering**: Interactive molecular visualization with customizable representations
- **2D Plots**: Publication-ready plots with matplotlib integration
- **Animation**: Create trajectory animations and movies

### Selection Module
Flexible atom selection system:

- **Basic Selections**: By name, type, residue, chain
- **Advanced Selections**: Distance-based, geometric selections
- **Custom Selections**: Define your own selection criteria

## System Requirements

- **Python**: 3.8 or higher
- **Memory**: 4 GB RAM minimum (8 GB recommended)
- **Operating System**: Linux, macOS, Windows 10+
- **Optional**: CUDA-capable GPU for acceleration

## Documentation Structure

- **[Getting Started](getting-started/index.md)**: Installation and basic usage
- **[User Guide](user-guide/index.md)**: Comprehensive feature documentation
- **[Tutorials](tutorials/index.md)**: Step-by-step guides and examples
- **[API Reference](api/index.md)**: Complete API documentation
- **[Examples](examples/index.md)**: Real-world use cases and scripts

## Contributing

We welcome contributions! Please see our [Contributing Guide](development/contributing.md) for details on how to get started.

### Development

```bash
# Clone repository
git clone https://github.com/your-username/PRISM.git
cd PRISM

# Install in development mode
pip install -e .[dev]

# Run tests
pytest tests/
```

## Citation

If you use PRISM in your research, please cite:

```bibtex
@software{prism2024,
  author = {Shi, Zhaoqi},
  title = {PRISM: Python-based Research Infrastructure for Structural Modeling},
  year = {2024},
  institution = {Theoretical Chemistry Institute, University of Wisconsin-Madison},
  url = {https://github.com/your-username/PRISM}
}
```

## Support

- **Documentation**: [https://prism-docs.github.io](https://prism-docs.github.io)
- **Issue Tracker**: [GitHub Issues](https://github.com/your-username/PRISM/issues)
- **Discussions**: [GitHub Discussions](https://github.com/your-username/PRISM/discussions)
- **Email**: zhaoqi.shi@wisc.edu

## License

PRISM is released under the [MIT License](about/license.md).

---

<p align="center">
  <strong>Developed at</strong><br>
  Theoretical Chemistry Institute<br>
  University of Wisconsin-Madison<br>
  <br>
  <strong>Author:</strong> Zhaoqi Shi
</p>