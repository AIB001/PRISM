#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Coordinate conversion helpers for the FEP package.

Internal FEP coordinates are represented in nanometers (nm). File-format
boundaries use these helpers to convert to or from Angstrom-based formats such
as PDB and MOL2.
"""

from __future__ import annotations

from typing import Iterable, Tuple, Union

import numpy as np

Number = Union[int, float]

ANGSTROM_PER_NM = 10.0


def angstrom_to_nm(value: Union[Number, Iterable[Number], np.ndarray]):
    """Convert a scalar or vector from Angstrom to nanometer units."""
    if isinstance(value, np.ndarray):
        return value / ANGSTROM_PER_NM
    if isinstance(value, (int, float)):
        return float(value) / ANGSTROM_PER_NM
    return tuple(float(component) / ANGSTROM_PER_NM for component in value)


def nm_to_angstrom(value: Union[Number, Iterable[Number], np.ndarray]):
    """Convert a scalar or vector from nanometer to Angstrom units."""
    if isinstance(value, np.ndarray):
        return value * ANGSTROM_PER_NM
    if isinstance(value, (int, float)):
        return float(value) * ANGSTROM_PER_NM
    return tuple(float(component) * ANGSTROM_PER_NM for component in value)


def angstrom_xyz_to_nm(x: Number, y: Number, z: Number) -> Tuple[float, float, float]:
    """Convert three Angstrom coordinate components to nanometers."""
    coords = angstrom_to_nm((x, y, z))
    return coords[0], coords[1], coords[2]


def nm_xyz_to_angstrom(x: Number, y: Number, z: Number) -> Tuple[float, float, float]:
    """Convert three nanometer coordinate components to Angstrom."""
    coords = nm_to_angstrom((x, y, z))
    return coords[0], coords[1], coords[2]
