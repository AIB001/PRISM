#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP Analysis Data Models

This module contains data classes for FEP analysis results.
These models are shared across analyzers and reports.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class FEResults:
    """
    Container for FEP analysis results

    Attributes
    ----------
    delta_g : float
        Free energy difference (kcal/mol)
    delta_g_error : float
        Standard error of delta_g (kcal/mol)
    delta_g_bound : float
        Free energy difference for bound leg (kcal/mol)
    delta_g_unbound : float
        Free energy difference for unbound leg (kcal/mol)
    delta_g_components : Dict[str, float]
        Free energy decomposition by component (elec, vdw, etc.)
    convergence : Dict
        Convergence diagnostics
    metadata : Dict
        Analysis metadata (temperature, lambda windows, etc.)
    repeat_results : List[Dict[str, float]]
        Individual results for each repeat (if multiple repeats analyzed)
        Format: [{"repeat": 1, "bound": 54.07, "unbound": 49.84, "ddG": 4.23}, ...]
    n_repeats : int
        Number of repeats analyzed
    lambda_profiles : Optional[Dict]
        Lambda-dependent profiles for plotting (dg vs lambda, dhdl vs lambda)
    time_convergence : Optional[Dict]
        Time convergence data for report plotting
    bootstrap : Optional[Dict]
        Bootstrap resampling results for report plotting
    """

    delta_g: float = 0.0
    delta_g_error: float = 0.0
    delta_g_bound: float = 0.0
    delta_g_unbound: float = 0.0
    delta_g_components: Dict[str, float] = field(default_factory=dict)
    convergence: Dict = field(default_factory=dict)
    metadata: Dict = field(default_factory=dict)
    repeat_results: List[Dict[str, float]] = field(default_factory=list)
    n_repeats: int = 1
    lambda_profiles: Optional[Dict[str, Any]] = None
    time_convergence: Optional[Dict[str, Any]] = None
    bootstrap: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict:
        """Convert results to dictionary"""
        return {
            "delta_g": self.delta_g,
            "delta_g_error": self.delta_g_error,
            "delta_g_bound": self.delta_g_bound,
            "delta_g_unbound": self.delta_g_unbound,
            "delta_g_components": self.delta_g_components,
            "convergence": self.convergence,
            "metadata": self.metadata,
            "repeat_results": self.repeat_results,
            "n_repeats": self.n_repeats,
            "lambda_profiles": self.lambda_profiles,
            "time_convergence": self.time_convergence,
            "bootstrap": self.bootstrap,
        }


@dataclass
class MultiEstimatorResults:
    """
    Container for multi-estimator FEP analysis results

    Attributes
    ----------
    methods : Dict[str, FEResults]
        Dictionary of results from each estimator (e.g., {'TI': FEResults, 'BAR': FEResults, 'MBAR': FEResults})
    comparison : Dict[str, Any]
        Comparison metrics between estimators including:
        - delta_g_range: Range of ΔG values (max - min)
        - delta_g_std: Standard deviation of ΔG values
        - delta_g_mean: Mean ΔG across methods
        - agreement: 'good' (<0.5), 'moderate' (0.5-1.0), 'poor' (>1.0)
        - best_method: Recommended estimator (usually MBAR)
        - diverged: True if std > 1.0 kcal/mol
    metadata : Dict[str, Any]
        Analysis metadata including:
        - temperature: Simulation temperature
        - estimators_used: List of estimator names that succeeded
        - n_estimators: Number of estimators
        - backend: Analysis backend used
    """

    methods: Dict[str, "FEResults"] = field(default_factory=dict)
    comparison: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """
        Convert results to dictionary for JSON serialization

        Returns
        -------
        Dict
            Dictionary with 'methods', 'comparison', and 'metadata' keys
        """
        return {
            "methods": {name: results.to_dict() for name, results in self.methods.items()},
            "comparison": self.comparison,
            "metadata": self.metadata,
        }
