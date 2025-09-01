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
]

# Optional requirements for different force fields
extras_require = {
    "gaff": [
        "acpype>=2021.0",
    ],
    "openff": [
        "openff-toolkit>=0.10.0",
        "openff-interchange>=0.2.0",
        "rdkit>=2021.0",
    ],
    "all": [
        "acpype>=2021.0",
        "openff-toolkit>=0.10.0",
        "openff-interchange>=0.2.0",
        "rdkit>=2021.0",
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
print("\nTo install PRISM:")
print("  pip install -e .                # Basic installation")
print("  pip install -e .[gaff]          # With GAFF support")
print("  pip install -e .[openff]        # With OpenFF support") 
print("  pip install -e .[all]           # With all force fields")
print("\nBasic usage:")
print("  import prism as pm")
print("  system = pm.system('protein.pdb', 'ligand.mol2')")
print("  output_dir = system.build()")