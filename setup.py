#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Setup script for PRISM
"""

from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Basic requirements
install_requires = [
    "numpy>=1.19.0",
    "pyyaml>=5.0",
    "MDAnalysis>=2.0.0",  # Required for analysis
    "matplotlib>=3.0.0",  # Required for plotting
    "seaborn>=0.11.0",    # Required for advanced plotting
    "pandas>=1.0.0",      # Required for data processing
    "scikit-learn>=1.0.0", # Required for clustering analysis
]

# Optional requirements for different force fields and advanced features
extras_require = {
    "gaff": [
        "acpype>=2021.0",
    ],
    "openff": [
        "openff-toolkit>=0.10.0",
        "openff-interchange>=0.2.0",
        "rdkit>=2021.0",
    ],
    "opls": [
        "requests>=2.25.0",
        "rdkit>=2021.0",  # Optional for coordinate alignment
    ],
    "swissparam": [
        # SwissParam now uses command-line API via curl (no Python dependencies)
    ],
    "protonation": [
        "meeko>=0.5.0",  # For hydrogen optimization and protonation
    ],
    "analysis": [
        "mdtraj>=1.9.0",      # Optional for enhanced trajectory analysis
        "scipy>=1.6.0",       # For statistical analysis
    ],
    "all": [
        "acpype>=2021.0",
        "openff-toolkit>=0.10.0",
        "openff-interchange>=0.2.0",
        "rdkit>=2021.0",
        "requests>=2.25.0",
        # SwissParam uses curl (no Python dependencies)
        "meeko>=0.5.0",
        "mdtraj>=1.9.0",
        "scipy>=1.6.0",
    ],
}

setup(
    name="prism-md",
    version="1.0.0",
    author="PRISM Development Team",
    author_email="",
    description="Protein Receptor Interaction Simulation Modeler - A tool for building protein-ligand MD systems",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/PRISM",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=install_requires,
    extras_require=extras_require,
    entry_points={
        "console_scripts": [
            "prism=prism.builder:main",
            "prism-builder=prism.builder:main",
        ],
    },
    include_package_data=True,
    package_data={
        "prism": ["configs/*.yaml"],
    },
    zip_safe=False,  # Important for proper package loading
)

print("Setup complete!")
print("\n" + "="*60)
print("PRISM INSTALLATION GUIDE")
print("="*60)
print("\nRECOMMENDED INSTALLATION (2025):")
print("  # 1. Install force field dependencies first (using mamba/conda)")
print("  mamba install -c conda-forge openff-toolkit ambertools")
print("  # OR using conda:")
print("  conda install -c conda-forge openff-toolkit ambertools")
print("\n  # 2. Then install PRISM:")
print("  pip install -e .                # Basic installation (includes core analysis)")
print("  pip install -e .[gaff]          # With GAFF support")
print("  pip install -e .[openff]        # With OpenFF support")
print("  pip install -e .[opls]          # With OPLS-AA support (LigParGen)")
print("  pip install -e .[swissparam]    # With SwissParam support (MMFF/MATCH/Hybrid)")
print("  pip install -e .[analysis]      # With enhanced analysis (MDTraj, SciPy)")
print("  pip install -e .[all]           # With all force fields and features")

print("\nALTERNATIVE (PIP ONLY - may have dependency issues):")
print("  pip install openff-toolkit openff-interchange")
print("  pip install -e .[all]")

print("\nFOR SCIDRAFT-STUDIO INTEGRATION (no force fields needed):")
print("  pip install -e .                # Basic installation sufficient")
print("  pip install -e .[analysis]      # With optional enhanced analysis")

print("\nTROUBLESHOOTING:")
print("  # If conda fails to find packages:")
print("  conda config --add channels conda-forge")
print("  conda config --set channel_priority strict")
print("  conda install conda-forge::openff-toolkit conda-forge::ambertools")
print("\n  # For complex dependency conflicts:")
print("  mamba create -n prism-env -c conda-forge python=3.9 openff-toolkit ambertools")
print("  conda activate prism-env")
print("  pip install -e .")

print("\n  # Check installation:")
print("  python -c \"import prism as pm; print(pm.check_dependencies())\"")

print("\nBASIC USAGE:")
print("  import prism as pm")
print("  system = pm.system('protein.pdb', 'ligand.mol2')")
print("  output_dir = system.build()")