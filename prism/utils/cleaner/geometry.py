#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Geometry calculation utilities for protein cleaning.
"""

import numpy as np


class GeometryMixin:
    """Mixin for geometry calculations."""

    def _calculate_min_distance(self, point: np.ndarray, coords: np.ndarray) -> float:
        """
        Calculate minimum distance from a point to a set of coordinates.

        Parameters
        ----------
        point : np.ndarray
            Single coordinate (3,)
        coords : np.ndarray
            Array of coordinates (n, 3)

        Returns
        -------
        float
            Minimum distance in Angstroms
        """
        if len(coords) == 0:
            return float("inf")

        distances = np.sqrt(np.sum((coords - point) ** 2, axis=1))
        return float(np.min(distances))
