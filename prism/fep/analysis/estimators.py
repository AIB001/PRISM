#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Free Energy Estimators

BAR and MBAR estimators for calculating free energy differences from FEP data.
"""

from typing import Dict, Tuple


class FEstimator:
    """
    Free energy estimator for FEP calculations

    Supports both BAR (Bennett Acceptance Ratio) and MBAR
    (Multistate Bennett Acceptance Ratio) methods.

    Notes
    -----
    BAR is suitable for adjacent λ window pairs
    MBAR uses all λ windows simultaneously for better statistical efficiency
    """

    def __init__(self):
        self.delta_g = None
        self.error = None

    def bar(self, data: Dict) -> Tuple[float, float]:
        """
        Calculate free energy difference using BAR

        Parameters
        ----------
        data : Dict
            Dictionary containing:
            - 'forward_work': Work values from forward perturbations
            - 'reverse_work': Work values from reverse perturbations
            - 'temperatures': Temperature information

        Returns
        -------
        Tuple[float, float]
            (delta_g, error) - Free energy difference and error estimate

        Raises
        ------
        NotImplementedError
            To be implemented in next phase (or use alchemlyb/pymbar)
        """
        raise NotImplementedError("BAR calculation to be implemented (consider using alchemlyb)")

    def mbar(self, data: Dict) -> Tuple[float, float]:
        """
        Calculate free energy difference using MBAR

        Parameters
        ----------
        data : Dict
            Dictionary containing:
            - 'energy_matrix': Energy matrix for all λ windows
            - 'samples': Number of samples per window
            - 'temperatures': Temperature information

        Returns
        -------
        Tuple[float, float]
            (delta_g, error) - Free energy difference and error estimate

        Raises
        ------
        NotImplementedError
            To be implemented in next phase (or use pymbar/alchemlyb)
        """
        raise NotImplementedError("MBAR calculation to be implemented (consider using pymbar/alchemlyb)")

    def compute_convergence(self, data: Dict) -> Dict:
        """
        Compute convergence diagnostics for FEP calculations

        Parameters
        ----------
        data : Dict
            Time series data of free energy estimates

        Returns
        -------
        Dict
            Convergence metrics:
            - 'converged': Boolean indicating convergence
            - 'hysteresis': Cycle closure hysteresis
            - 'time_convergence': Time series convergence info
        """
        raise NotImplementedError("Convergence computation to be implemented")
