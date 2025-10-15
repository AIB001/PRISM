#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Resources module for PRISM analysis.

This module contains external scripts and resources used by analysis modules.
"""

import os
from pathlib import Path

# Get path to resources directory
RESOURCES_DIR = Path(__file__).parent

def get_resource_path(resource_name: str) -> str:
    """
    Get absolute path to a resource file.

    Parameters
    ----------
    resource_name : str
        Name of the resource file

    Returns
    -------
    str
        Absolute path to the resource file

    Raises
    ------
    FileNotFoundError
        If the resource file does not exist
    """
    resource_path = RESOURCES_DIR / resource_name
    if not resource_path.exists():
        raise FileNotFoundError(f"Resource not found: {resource_name}")
    return str(resource_path)

def get_namd_fep_script() -> str:
    """
    Get path to NAMD FEP decomposition bash script.

    Returns
    -------
    str
        Absolute path to mknamd_fep_decomp_convergence.sh
    """
    return get_resource_path('mknamd_fep_decomp_convergence.sh')
