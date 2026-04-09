#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Force field generators for PRISM.

Import failures are tracked in ``get_import_errors()`` so callers can diagnose
why a generator is unavailable instead of seeing it silently disappear from the
registry.
"""

from .registry import iter_forcefield_specs

# Try to import all available generators
__all__ = []
_IMPORT_ERRORS = {}


def _record_import_error(label, exc):
    """Track optional import failures for user-facing diagnostics."""
    _IMPORT_ERRORS[label] = f"{type(exc).__name__}: {exc}"


# Base class (always available)
try:
    from .base import ForceFieldGeneratorBase

    __all__.append("ForceFieldGeneratorBase")
except ImportError as exc:
    _record_import_error("ForceFieldGeneratorBase", exc)

# GAFF force field
try:
    from .gaff import GAFFForceFieldGenerator

    __all__.append("GAFFForceFieldGenerator")
except ImportError as exc:
    _record_import_error("GAFFForceFieldGenerator", exc)

# GAFF2 force field
try:
    from .gaff2 import GAFF2ForceFieldGenerator

    __all__.append("GAFF2ForceFieldGenerator")
except ImportError as exc:
    _record_import_error("GAFF2ForceFieldGenerator", exc)

# OpenFF force field
try:
    from .openff import OpenFFForceFieldGenerator

    __all__.append("OpenFFForceFieldGenerator")
except ImportError as exc:
    _record_import_error("OpenFFForceFieldGenerator", exc)

# CGenFF force field
try:
    from .cgenff import CGenFFForceFieldGenerator

    __all__.append("CGenFFForceFieldGenerator")
except ImportError as exc:
    _record_import_error("CGenFFForceFieldGenerator", exc)

# CHARMM-GUI force field
try:
    from .charmm_gui import CHARMMGUIForceFieldGenerator

    __all__.append("CHARMMGUIForceFieldGenerator")
except ImportError as exc:
    _record_import_error("CHARMMGUIForceFieldGenerator", exc)

# RTF force field
try:
    from .rtf import RTFForceFieldGenerator

    __all__.append("RTFForceFieldGenerator")
except ImportError as exc:
    _record_import_error("RTFForceFieldGenerator", exc)

# OPLS-AA force field
try:
    from .opls_aa import OPLSAAForceFieldGenerator

    __all__.append("OPLSAAForceFieldGenerator")
except ImportError as exc:
    _record_import_error("OPLSAAForceFieldGenerator", exc)

# SwissParam-based force fields
try:
    from .swissparam import (
        SwissParamForceFieldGenerator,
        MMFFForceFieldGenerator,
        MATCHForceFieldGenerator,
        BothForceFieldGenerator,
        HybridForceFieldGenerator,  # Alias for backward compatibility
    )

    __all__.extend(
        [
            "SwissParamForceFieldGenerator",
            "MMFFForceFieldGenerator",
            "MATCHForceFieldGenerator",
            "BothForceFieldGenerator",
            "HybridForceFieldGenerator",
        ]
    )
except ImportError as exc:
    _record_import_error("SwissParamForceFieldGenerator", exc)

# Converters
try:
    from .converters import AmberToCharmmConverter

    __all__.append("AmberToCharmmConverter")
except ImportError as exc:
    _record_import_error("AmberToCharmmConverter", exc)

# Gaussian RESP module (for charge replacement)
try:
    from ..gaussian import GaussianRESPWorkflow, RESPChargeReplacer, replace_itp_charges

    __all__.extend(["GaussianRESPWorkflow", "RESPChargeReplacer", "replace_itp_charges"])
except ImportError as exc:
    _record_import_error("GaussianRESPWorkflow", exc)


def list_available_generators():
    """
    List all available force field generators

    Returns:
    --------
    list : Names of available generators
    """
    return __all__


def get_import_errors():
    """Return optional generator import failures for diagnostics."""
    return dict(_IMPORT_ERRORS)


def get_generator_info():
    """
    Get information about available force field generators

    Returns:
    --------
    dict : Dictionary with generator names and descriptions
    """
    descriptions = {
        "gaff": "GAFF force field using AmberTools",
        "gaff2": "GAFF2 force field using AmberTools (improved version with better torsion parameters)",
        "openff": "OpenFF force field",
        "cgenff": "CGenFF (CHARMM General Force Field) - requires web-downloaded files",
        "charmm-gui": "CHARMM-GUI output converter - processes gromacs/ directory",
        "rtf": "RTF+PRM force field converter - processes CGenFF RTF/PRM files",
        "opls": "OPLS-AA force field using LigParGen server",
        "mmff": "MMFF-based force field using SwissParam API",
        "match": "MATCH force field using SwissParam API",
        "hybrid": "Both MMFF-based + MATCH force field using SwissParam API",
    }
    info = {}
    for spec in iter_forcefield_specs():
        if spec.generator_class_name not in __all__:
            continue
        info[spec.display_name] = {
            "class": spec.generator_class_name,
            "description": descriptions.get(spec.key, spec.display_name),
            "output_dir": spec.generator_output_subdir,
            "dependencies": list(spec.dependencies),
        }
    if "HybridForceFieldGenerator" in __all__:
        hybrid_spec = next(spec for spec in iter_forcefield_specs() if spec.key == "hybrid")
        info["Hybrid (deprecated)"] = {
            "class": "HybridForceFieldGenerator",
            "description": "Alias for BothForceFieldGenerator (backward compatibility)",
            "output_dir": hybrid_spec.generator_output_subdir,
            "dependencies": list(hybrid_spec.dependencies),
        }
    if _IMPORT_ERRORS:
        info["_unavailable"] = dict(sorted(_IMPORT_ERRORS.items()))
    return info
