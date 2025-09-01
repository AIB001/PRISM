# Installation Guide

This guide provides detailed instructions for installing PRISM and its dependencies on various platforms.

## 📋 Prerequisites

PRISM requires several components to be installed before you can use it. Please ensure you have the following prerequisites installed in order.

### 1. GROMACS (Required)

GROMACS is essential for PRISM's molecular dynamics simulations. Follow these steps to install GROMACS with GPU support:

#### Check CUDA Toolkit
First, verify that you have CUDA installed:

```bash
nvcc --version
```

If CUDA is not installed, please install it from [NVIDIA's CUDA Toolkit page](https://developer.nvidia.com/cuda-downloads) before proceeding.

#### Install Build Dependencies

=== "Ubuntu/Debian"

    ```bash
    # Install essential build tools
    sudo apt-get update
    sudo apt-get install gcc g++ cmake
    sudo apt-get install build-essential
    ```

=== "CentOS/RHEL"

    ```bash
    # Install essential build tools
    sudo yum groupinstall "Development Tools"
    sudo yum install cmake gcc-c++
    ```

=== "macOS"

    ```bash
    # Install using Homebrew
    brew install cmake gcc
    ```

#### Download and Build GROMACS

```bash
# Download GROMACS 2024.3
wget https://ftp.gromacs.org/gromacs/gromacs-2024.3.tar.gz
tar xfz gromacs-2024.3.tar.gz

# Create build directory
cd gromacs-2024.3
mkdir build
cd build

# Configure with CMake (adjust paths as needed)
cmake .. -DGMX_MPI=ON \
         -DGMX_BUILD_OWN_FFTW=ON \
         -DGMX_GPU=CUDA \
         -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
         -DCUDA_INCLUDE_DIRS=/usr/local/cuda/include \
         -DCUDA_CUDART_LIBRARY=/usr/local/cuda/lib64 \
         -DCMAKE_INSTALL_PREFIX=/opt/gromacs-2024.3

# Build GROMACS (use all available cores)
make -j$(nproc)

# Run tests (optional but recommended)
make check

# Install GROMACS
sudo make install

# Source GROMACS environment
source /opt/gromacs-2024.3/bin/GMXRC
```

!!! tip "Add to your shell configuration"
    Add the following line to your `~/.bashrc` or `~/.zshrc` to automatically source GROMACS:
    ```bash
    source /opt/gromacs-2024.3/bin/GMXRC
    ```

#### Install Bioconda (Optional)

For additional bioinformatics tools:

```bash
# Install bioconda channel
conda config --add channels bioconda
conda config --add channels conda-forge
```

### 2. Python Environment (Required)

PRISM requires Python 3.10 or higher. We recommend using Conda for environment management:

```bash
# Create a new environment for PRISM
conda create -n prism python=3.10
conda activate prism
```

!!! info "Python Version"
    While Python 3.10 is recommended, PRISM supports Python 3.8 through 3.11. Python 3.10 offers the best balance of compatibility and performance.

### 3. Core Python Dependencies (Required)

Install essential Python packages:

```bash
# Activate PRISM environment
conda activate prism

# Install PDBFixer and core scientific packages
conda install -c conda-forge pdbfixer numpy scipy

# Install YAML support
pip install pyyaml
```

## 🧪 Force Field Specific Dependencies

Depending on which force fields you plan to use, install the appropriate packages:

### For GAFF Support

The General AMBER Force Field (GAFF) is commonly used for small molecules:

```bash
# AmberTools (required for GAFF)
conda install -c conda-forge ambertools

# ACPYPE for AMBER topology conversion (required)
pip install acpype

# RDKit for molecular structure handling (optional but recommended)
conda install -c conda-forge rdkit
```

### For OpenFF Support

The Open Force Field Initiative provides modern, reproducible force fields:

```bash
# OpenFF toolkit and dependencies
conda install -c conda-forge openff-toolkit openff-interchange

# RDKit for SDF file handling (required)
conda install -c conda-forge rdkit

# OpenBabel for file format conversion (optional but recommended)
conda install -c conda-forge openbabel
```

### For CHARMM Support

```bash
# CHARMM-GUI integration tools (optional)
pip install charmm-gui-api
```

## 📦 Installing PRISM

Once all prerequisites are installed, you can install PRISM itself:

### Method 1: Development Installation (Recommended)

This method allows you to modify PRISM and see changes immediately:

```bash
# Clone the PRISM repository
git clone https://github.com/your-username/PRISM.git
cd PRISM

# Install in development mode
pip install -e .
```

### Method 2: Direct Usage

You can use PRISM directly without installation:

```bash
# Run PRISM directly
python /path/to/PRISM/prism/builder.py
```

### Method 3: Standard Installation

For a standard installation:

```bash
# Clone the repository
git clone https://github.com/your-username/PRISM.git
cd PRISM

# Standard installation
pip install .
```

## ✅ Verification

After installation, verify that everything is working correctly:

### Test GROMACS Installation

```bash
# Check GROMACS version
gmx --version

# Should show version 2024.3 with CUDA support
```

### Test Python Environment

```python
# Activate PRISM environment
conda activate prism

# Test imports
python -c "import prism; print(f'PRISM version: {prism.__version__}')"
python -c "import pdbfixer; print('PDBFixer: OK')"
python -c "import yaml; print('YAML: OK')"
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"
```

### Test Force Field Support

```python
# Test GAFF support (if installed)
python -c "import acpype; print('ACPYPE: OK')"
python -c "import parmed; print('AmberTools: OK')"

# Test OpenFF support (if installed)
python -c "import openff.toolkit; print('OpenFF: OK')"
python -c "import openff.interchange; print('Interchange: OK')"
```

### Run PRISM Test Suite

```bash
# Navigate to PRISM directory
cd PRISM

# Run tests
pytest tests/

# Or run a quick test
python -m prism.tests.quick_test
```

## 🔧 Environment Variables

Set these environment variables for optimal performance:

```bash
# Add to ~/.bashrc or ~/.zshrc

# GROMACS path
export GMXRC=/opt/gromacs-2024.3/bin/GMXRC
source $GMXRC

# CUDA paths (adjust as needed)
export CUDA_HOME=/usr/local/cuda
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

# PRISM path (optional, for direct usage)
export PRISM_HOME=/path/to/PRISM
export PATH=$PRISM_HOME/bin:$PATH
```

## 🐛 Troubleshooting

### Common Installation Issues

??? failure "GROMACS build fails with CUDA errors"

    **Solution**: Ensure CUDA paths are correct:
    ```bash
    # Find CUDA installation
    which nvcc
    ls /usr/local/cuda*
    
    # Update CMake command with correct paths
    cmake .. -DCUDA_TOOLKIT_ROOT_DIR=/correct/cuda/path ...
    ```

??? failure "ImportError: No module named 'prism'"

    **Solution**: Ensure you're in the correct conda environment:
    ```bash
    conda activate prism
    which python  # Should point to conda environment
    pip list | grep prism
    ```

??? failure "PDBFixer installation fails"

    **Solution**: Install from conda-forge channel:
    ```bash
    conda install -c conda-forge pdbfixer
    # If still fails, try:
    conda install -c conda-forge openmm pdbfixer
    ```

??? failure "ACPYPE not found"

    **Solution**: Install with pip in the correct environment:
    ```bash
    conda activate prism
    pip install acpype --upgrade
    ```

### Platform-Specific Issues

#### GPU Support on Linux

If GPU is not detected:

```bash
# Check NVIDIA driver
nvidia-smi

# Check CUDA installation
nvcc --version
ldconfig -p | grep cuda

# Reinstall CUDA toolkit if necessary
```

#### macOS Specific

For macOS users without CUDA:

```bash
# Build GROMACS without GPU support
cmake .. -DGMX_MPI=ON \
         -DGMX_BUILD_OWN_FFTW=ON \
         -DGMX_GPU=OFF \
         -DCMAKE_INSTALL_PREFIX=/opt/gromacs-2024.3
```

## 📚 Additional Resources

- [GROMACS Installation Guide](http://manual.gromacs.org/documentation/current/install-guide/index.html)
- [AmberTools Documentation](https://ambermd.org/AmberTools.php)
- [OpenFF Documentation](https://docs.openforcefield.org/)
- [PRISM GitHub Repository](https://github.com/your-username/PRISM)

## 🤝 Getting Help

If you encounter issues during installation:

1. Check the [FAQ](../about/faq.md)
2. Search [GitHub Issues](https://github.com/AIB001/PRISM/issues)
3. Post on [GitHub Discussions](https://github.com/AIB001/PRISM/discussions)
4. Contact: [author email](mailto:zshi268@wisc.edu)

---

## ✨ Next Steps

After successful installation:

<div style="text-align: center; margin-top: 2rem;">
  <a href="quickstart/" class="md-button md-button--primary">Quick Start Tutorial</a>
  <a href="concepts/" class="md-button">Basic Concepts</a>
</div>