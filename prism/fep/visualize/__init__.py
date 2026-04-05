#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP Visualization Module

Provides visualization tools for atom mapping results in FEP calculations.
Supports both PNG static images and interactive HTML visualizations.
"""

from .molecule import (
    pdb_to_mol,
    assign_bond_orders_from_mol2,
    assign_bond_orders_from_mol2_by_coords,
    prepare_mol_with_charges_and_labels,
)
from .highlight import create_highlight_info
from .mapping import visualize_mapping_png
from .html import visualize_mapping_html

__all__ = [
    "pdb_to_mol",
    "assign_bond_orders_from_mol2",
    "assign_bond_orders_from_mol2_by_coords",
    "prepare_mol_with_charges_and_labels",
    "create_highlight_info",
    "visualize_mapping_png",
    "visualize_mapping_html",
]
