#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Gaussian RESP Charge Calculation Module

This module provides tools for calculating high-precision RESP charges
using Gaussian quantum chemistry software, as an alternative to the
default AM1-BCC charges.

Features:
- MOL2 to XYZ coordinate conversion
- Gaussian input file generation (optimization and ESP)
- RESP charge extraction and ITP file updating
- Complete workflow management with script generation

Usage:
    # Direct API usage
    from prism.gaussian import GaussianRESPWorkflow

    workflow = GaussianRESPWorkflow(
        ligand_path="ligand.mol2",
        output_dir="gaussian_resp",
        method='dft',
        do_optimization=False
    )
    files = workflow.prepare_files()

    # Or use via PRISM CLI
    # prism protein.pdb ligand.mol2 -o output --gaussian dft
"""

from .converter import CoordinateConverter, mol2_to_xyz
from .gjf_generator import GaussianInputGenerator
from .charge_replacer import RESPChargeReplacer, replace_itp_charges
from .workflow import GaussianRESPWorkflow
from .utils import (
    check_gaussian_available,
    check_antechamber_available,
    calculate_charge_multiplicity,
    ATOMIC_SYMBOLS
)

__all__ = [
    # Converter
    "CoordinateConverter",
    "mol2_to_xyz",

    # GJF Generator
    "GaussianInputGenerator",

    # Charge Replacer
    "RESPChargeReplacer",
    "replace_itp_charges",

    # Workflow
    "GaussianRESPWorkflow",

    # Utilities
    "check_gaussian_available",
    "check_antechamber_available",
    "calculate_charge_multiplicity",
    "ATOMIC_SYMBOLS",
]
