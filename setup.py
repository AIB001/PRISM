# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

# """
# Setup script for PRISM
# """

# from setuptools import setup, find_packages
# import os

# # Read the README file
# with open("README.md", "r", encoding="utf-8") as fh:
#     long_description = fh.read()

# # Basic requirements
# install_requires = [
#     "numpy>=1.19.0",
#     "pyyaml>=5.0",
# ]

# # Optional requirements for different force fields
# extras_require = {
#     "gaff": [
#         "acpype>=2021.0",
#     ],
#     "openff": [
#         "openff-toolkit>=0.10.0",
#         "openff-interchange>=0.2.0",
#         "rdkit>=2021.0",
#     ],
#     "all": [
#         "acpype>=2021.0",
#         "openff-toolkit>=0.10.0",
#         "openff-interchange>=0.2.0",
#         "rdkit>=2021.0",
#     ],
# }

# setup(
#     name="prism-md",
#     version="1.0.0",
#     author="PRISM Development Team",
#     author_email="",
#     description="Protein Receptor Interaction Simulation Modeler - A tool for building protein-ligand MD systems",
#     long_description=long_description,
#     long_description_content_type="text/markdown",
#     url="https://github.com/yourusername/PRISM",
#     packages=find_packages(),
#     classifiers=[
#         "Development Status :: 4 - Beta",
#         "Intended Audience :: Science/Research",
#         "Topic :: Scientific/Engineering :: Chemistry",
#         "Topic :: Scientific/Engineering :: Bio-Informatics",
#         "License :: OSI Approved :: MIT License",
#         "Programming Language :: Python :: 3",
#         "Programming Language :: Python :: 3.8",
#         "Programming Language :: Python :: 3.9",
#         "Programming Language :: Python :: 3.10",
#         "Programming Language :: Python :: 3.11",
#     ],
#     python_requires=">=3.8",
#     install_requires=install_requires,
#     extras_require=extras_require,
#     entry_points={
#         "console_scripts": [
#             "prism=prism.builder:main",
#             "prism-builder=prism.builder:main",
#         ],
#     },
#     include_package_data=True,
#     package_data={
#         "prism": ["configs/*.yaml"],
#     },
# )

# # Create a requirements.txt file
# requirements_content = """# Core requirements
# numpy>=1.19.0
# pyyaml>=5.0

# # For GAFF support (optional)
# # Uncomment if using GAFF force field
# # acpype>=2021.0

# # For OpenFF support (optional)
# # Uncomment if using OpenFF force field
# # openff-toolkit>=0.10.0
# # openff-interchange>=0.2.0
# # rdkit>=2021.0

# # System dependencies (install separately)
# # - GROMACS
# # - AmberTools (for GAFF)
# # - PDBFixer
# """

# # Write requirements.txt
# with open("requirements.txt", "w") as f:
#     f.write(requirements_content)

# print("Setup complete!")
# print("\nTo install PRISM with all force fields:")
# print("  pip install -e .[all]")
# print("\nTo install with specific force field support:")
# print("  pip install -e .[gaff]    # For GAFF only")
# print("  pip install -e .[openff]  # For OpenFF only")

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Setup script for PRISM PMF System v2.0.0
Enhanced setup with complete CLI integration and optimization systems
"""

from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Enhanced requirements for v2.0.0 optimization systems
install_requires = [
    # Core dependencies
    "numpy>=1.19.0",
    "scipy>=1.7.0",
    "matplotlib>=3.3.0",
    "pandas>=1.3.0",
    "pyyaml>=5.0",
    
    # CLI and system management
    "click>=8.0.0",
    "rich>=10.0.0",
    "tabulate>=0.8.0",
    "psutil>=5.8.0",
    
    # Data management and analysis
    "h5py>=3.1.0",
    "joblib>=1.0.0",
    
    # Monitoring and logging
    "watchdog>=2.0.0",
    "coloredlogs>=15.0",
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
    "advanced": [
        # Advanced analysis and visualization
        "plotly>=5.0.0",
        "seaborn>=0.11.0",
        "scikit-learn>=1.0.0",
        
        # Performance optimization
        "numba>=0.55.0",
        "dask>=2021.0.0",
        
        # Web dashboard (optional)
        "flask>=2.0.0",
        "websockets>=10.0",
    ],
    "ml": [
        # Machine learning components
        "scikit-learn>=1.0.0",
        "tensorflow>=2.8.0",
        "torch>=1.11.0",
    ],
    "testing": [
        # Testing and quality assurance
        "pytest>=6.0.0",
        "pytest-cov>=3.0.0",
        "black>=22.0.0",
        "flake8>=4.0.0",
        "mypy>=0.950",
    ],
    "all": [
        # All optional dependencies
        "acpype>=2021.0",
        "openff-toolkit>=0.10.0",
        "openff-interchange>=0.2.0",
        "rdkit>=2021.0",
        "plotly>=5.0.0",
        "seaborn>=0.11.0",
        "scikit-learn>=1.0.0",
        "numba>=0.55.0",
        "dask>=2021.0.0",
        "flask>=2.0.0",
        "websockets>=10.0",
        "pytest>=6.0.0",
        "pytest-cov>=3.0.0",
    ],
}

setup(
    name="prism-md",
    version="2.0.0",
    author="PRISM Development Team",
    author_email="prism@research.org",
    description="Advanced Molecular Dynamics PMF Calculation Framework with 10 Optimization Systems",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/PRISM",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords="molecular dynamics, PMF, free energy, GROMACS, umbrella sampling, drug discovery",
    python_requires=">=3.8",
    install_requires=install_requires,
    extras_require=extras_require,
    
    # Complete CLI system integration
    entry_points={
        "console_scripts": [
            # Main CLI command with full optimization systems
            "prism=prism.cli.main:main",
            
            # Legacy compatibility commands
            "prism-builder=prism.builder:main",
            "prism-legacy=prism.builder:main",
            
            # Specialized commands for different workflows
            "prism-pmf=prism.pmf.workflow:cli_main",
            "prism-workflow=prism.cli.main:main",
            
            # Quick start scripts integration
            "quick-pmf=quick_start_pmf:main",
            "start-pmf=start_pmf_calculation:main",
            
            # System management utilities
            "prism-monitor=prism.cli.commands.monitor:main",
            "prism-config=prism.cli.commands.config:main",
            "prism-analysis=prism.cli.commands.analysis:main",
        ],
    },
    
    # Package data and resources
    include_package_data=True,
    package_data={
        "prism": [
            "configs/*.yaml",
            "configs/*.json", 
            "templates/*.mdp",
            "templates/*.sh",
            "utils/templates/*.yaml",
            "cli/templates/*",
        ],
        "prism.pmf": [
            "templates/*.mdp",
            "configs/*.yaml",
        ],
        "prism.testing": [
            "fixtures/*",
            "benchmarks/*",
        ],
    },
    
    # Additional metadata for v2.0.0
    project_urls={
        "Documentation": "https://prism-pmf.readthedocs.io/",
        "Source": "https://github.com/yourusername/PRISM",
        "Tracker": "https://github.com/yourusername/PRISM/issues",
        "Changelog": "https://github.com/yourusername/PRISM/blob/main/CHANGELOG.md",
    },
    
    zip_safe=False,  # Important for proper package loading
)

print("Setup complete!")
print("\nTo install PRISM:")
print("  pip install -e .                # Basic installation")
print("  pip install -e .[gaff]          # With GAFF support")
print("  pip install -e .[openff]        # With OpenFF support") 
print("  pip install -e .[all]           # With all force fields")
print("\nBasic usage:")
print("  import prism as pm")
print("  system = pm.system('protein.pdb', 'ligand.mol2')")
print("  output_dir = system.build()")