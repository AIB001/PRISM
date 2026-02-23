#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRISM Analysis Resources

This module provides access to analysis resource files like shell scripts.
"""

from pathlib import Path


def get_namd_fep_script() -> Path:
    """
    Get path to NAMD FEP decomposition shell script.

    Returns
    -------
    Path
        Absolute path to mknamd_fep_decomp_convergence.sh
    """
    return Path(__file__).parent / "mknamd_fep_decomp_convergence.sh"


__all__ = ["get_namd_fep_script"]
