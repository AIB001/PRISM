#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Force field generators for PRISM

Available force field generators:
- GAFFForceFieldGenerator: GAFF force field using AmberTools
- GAFF2ForceFieldGenerator: GAFF2 force field using AmberTools (improved version)
- OpenFFForceFieldGenerator: OpenFF force field
- CGenFFForceFieldGenerator: CGenFF (CHARMM General Force Field) from web server
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

# GAFF2 force field
try:
    from .gaff2 import GAFF2ForceFieldGenerator
    __all__.append("GAFF2ForceFieldGenerator")
except ImportError:
    pass

# OpenFF force field
try:
    from .openff import OpenFFForceFieldGenerator
    __all__.append("OpenFFForceFieldGenerator")
except ImportError:
    pass

# CGenFF force field
try:
    from .cgenff import CGenFFForceFieldGenerator
    __all__.append("CGenFFForceFieldGenerator")
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
        BothForceFieldGenerator,
        HybridForceFieldGenerator  # Alias for backward compatibility
    )
    __all__.extend([
        "SwissParamForceFieldGenerator",
        "MMFFForceFieldGenerator",
        "MATCHForceFieldGenerator",
        "BothForceFieldGenerator",
        "HybridForceFieldGenerator"
    ])
except ImportError:
    pass

# Converters
try:
    from .converters import AmberToCharmmConverter
    __all__.append("AmberToCharmmConverter")
except ImportError:
    pass

# Gaussian RESP module (for charge replacement)
try:
    from ..gaussian import (
        GaussianRESPWorkflow,
        RESPChargeReplacer,
        replace_itp_charges
    )
    __all__.extend([
        "GaussianRESPWorkflow",
        "RESPChargeReplacer",
        "replace_itp_charges"
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

    if "GAFF2ForceFieldGenerator" in __all__:
        info["GAFF2"] = {
            "class": "GAFF2ForceFieldGenerator",
            "description": "GAFF2 force field using AmberTools (improved version with better torsion parameters)",
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

    if "CGenFFForceFieldGenerator" in __all__:
        info["CGenFF"] = {
            "class": "CGenFFForceFieldGenerator",
            "description": "CGenFF (CHARMM General Force Field) - requires web-downloaded files",
            "output_dir": "LIG.cgenff2gmx",
            "dependencies": ["Web download from https://cgenff.com/"]
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
            "description": "MMFF-based force field using SwissParam API",
            "output_dir": "LIG.mmff2gmx",
            "dependencies": ["curl (command-line tool)"]
        }

    if "MATCHForceFieldGenerator" in __all__:
        info["MATCH"] = {
            "class": "MATCHForceFieldGenerator",
            "description": "MATCH force field using SwissParam API",
            "output_dir": "LIG.match2gmx",
            "dependencies": ["curl (command-line tool)"]
        }

    if "BothForceFieldGenerator" in __all__:
        info["Both"] = {
            "class": "BothForceFieldGenerator",
            "description": "Both MMFF-based + MATCH force field using SwissParam API",
            "output_dir": "LIG.both2gmx",
            "dependencies": ["curl (command-line tool)"]
        }

    if "HybridForceFieldGenerator" in __all__:
        info["Hybrid (deprecated)"] = {
            "class": "HybridForceFieldGenerator",
            "description": "Alias for BothForceFieldGenerator (backward compatibility)",
            "output_dir": "LIG.both2gmx",
            "dependencies": ["curl (command-line tool)"]
        }

    return info
    
