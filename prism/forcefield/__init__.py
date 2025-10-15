#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Force field generators for PRISM

Available force field generators:
- GAFFForceFieldGenerator: GAFF force field using AmberTools
- OpenFFForceFieldGenerator: OpenFF force field
- OPLSAAForceFieldGenerator: OPLS-AA force field using LigParGen
- MMFFForceFieldGenerator: MMFF-based force field using SwissParam
- MATCHForceFieldGenerator: MATCH force field using SwissParam
- HybridMMFFMATCHForceFieldGenerator: Hybrid MMFF-based-MATCH force field
"""

# Try to import all available generators
__all__ = []

# Base class (always available)
try:
    from .base import ForceFieldGeneratorBase
    __all__.append("ForceFieldGeneratorBase")
except ImportError:
    pass

# GAFF force field
try:
    from .gaff import GAFFForceFieldGenerator
    __all__.append("GAFFForceFieldGenerator")
except ImportError:
    pass

# OpenFF force field
try:
    from .openff import OpenFFForceFieldGenerator
    __all__.append("OpenFFForceFieldGenerator")
except ImportError:
    pass

# OPLS-AA force field
try:
    from .opls_aa import OPLSAAForceFieldGenerator
    __all__.append("OPLSAAForceFieldGenerator")
except ImportError:
    pass

# SwissParam-based force fields
try:
    from .swissparam import (
        SwissParamForceFieldGenerator,
        MMFFForceFieldGenerator,
        MATCHForceFieldGenerator,
        HybridMMFFMATCHForceFieldGenerator
    )
    __all__.extend([
        "SwissParamForceFieldGenerator",
        "MMFFForceFieldGenerator",
        "MATCHForceFieldGenerator",
        "HybridMMFFMATCHForceFieldGenerator"
    ])
except ImportError:
    pass


def list_available_generators():
    """
    List all available force field generators

    Returns:
    --------
    list : Names of available generators
    """
    return __all__


def get_generator_info():
    """
    Get information about available force field generators

    Returns:
    --------
    dict : Dictionary with generator names and descriptions
    """
    info = {}

    if "GAFFForceFieldGenerator" in __all__:
        info["GAFF"] = {
            "class": "GAFFForceFieldGenerator",
            "description": "GAFF force field using AmberTools",
            "output_dir": "LIG.amb2gmx",
            "dependencies": ["AmberTools (antechamber, parmchk2, tleap)"]
        }

    if "OpenFFForceFieldGenerator" in __all__:
        info["OpenFF"] = {
            "class": "OpenFFForceFieldGenerator",
            "description": "OpenFF force field",
            "output_dir": "LIG.openff2gmx",
            "dependencies": ["openff-toolkit", "openff-interchange"]
        }

    if "OPLSAAForceFieldGenerator" in __all__:
        info["OPLS-AA"] = {
            "class": "OPLSAAForceFieldGenerator",
            "description": "OPLS-AA force field using LigParGen server",
            "output_dir": "LIG.opls2gmx",
            "dependencies": ["requests", "rdkit (optional, for alignment)"]
        }

    if "MMFFForceFieldGenerator" in __all__:
        info["MMFF"] = {
            "class": "MMFFForceFieldGenerator",
            "description": "MMFF-based force field using SwissParam server",
            "output_dir": "LIG.mmff2gmx",
            "dependencies": ["mechanize"]
        }

    if "MATCHForceFieldGenerator" in __all__:
        info["MATCH"] = {
            "class": "MATCHForceFieldGenerator",
            "description": "MATCH force field using SwissParam server",
            "output_dir": "LIG.match2gmx",
            "dependencies": ["mechanize"]
        }

    if "HybridMMFFMATCHForceFieldGenerator" in __all__:
        info["Hybrid MMFF-MATCH"] = {
            "class": "HybridMMFFMATCHForceFieldGenerator",
            "description": "Hybrid MMFF-based-MATCH force field using SwissParam",
            "output_dir": "LIG.hybrid2gmx",
            "dependencies": ["mechanize"]
        }

    return info
    