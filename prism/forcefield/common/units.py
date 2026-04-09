#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unit conversion utilities for PRISM forcefield module.

This module provides unified functions for converting between different
units commonly used in molecular simulations:
- Angstroms (Å) - used in PDB, MOL2 files
- Nanometers (nm) - used in GROMACS GRO files
"""

import numpy as np


def angstrom_to_nm(coord_angstrom: float) -> float:
    """Convert coordinate from Angstroms to nanometers.

    Parameters
    ----------
    coord_angstrom : float
        Coordinate value in Angstroms

    Returns
    -------
    float
        Coordinate value in nanometers

    Examples
    --------
    >>> angstrom_to_nm(10.0)
    1.0
    """
    return coord_angstrom / 10.0


def nm_to_angstrom(coord_nm: float) -> float:
    """Convert coordinate from nanometers to Angstroms.

    Parameters
    ----------
    coord_nm : float
        Coordinate value in nanometers

    Returns
    -------
    float
        Coordinate value in Angstroms

    Examples
    --------
    >>> nm_to_angstrom(1.0)
    10.0
    """
    return coord_nm * 10.0


def angstrom_xyz_to_nm(xyz_angstrom: np.ndarray) -> np.ndarray:
    """Convert XYZ coordinates from Angstroms to nanometers.

    Parameters
    ----------
    xyz_angstrom : np.ndarray
        XYZ coordinates in Angstroms, shape (3,) or (N, 3)

    Returns
    -------
    np.ndarray
        XYZ coordinates in nanometers, same shape as input

    Examples
    --------
    >>> angstrom_xyz_to_nm(np.array([10.0, 20.0, 30.0]))
    array([1., 2., 3.])
    """
    return xyz_angstrom / 10.0


def nm_xyz_to_angstrom(xyz_nm: np.ndarray) -> np.ndarray:
    """Convert XYZ coordinates from nanometers to Angstroms.

    Parameters
    ----------
    xyz_nm : np.ndarray
        XYZ coordinates in nanometers, shape (3,) or (N, 3)

    Returns
    -------
    np.ndarray
        XYZ coordinates in Angstroms, same shape as input

    Examples
    --------
    >>> nm_xyz_to_angstrom(np.array([1.0, 2.0, 3.0]))
    array([10., 20., 30.])
    """
    return xyz_nm * 10.0
