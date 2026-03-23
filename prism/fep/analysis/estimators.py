#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Free Energy Estimators

BAR and MBAR estimators for calculating free energy differences from FEP data.
"""

from pathlib import Path
from typing import Any, Dict, List, Tuple, Optional
import logging
import re
import subprocess

import pandas as pd


def compute_free_energy_alchemlyb(
    estimator_cls: Any, datasets: List[pd.DataFrame], logger: logging.Logger
) -> Tuple[float, float, Any, Optional[Any]]:
    """Fit an alchemlyb estimator and return ΔG/error in kcal/mol.

    Returns
    -------
    Tuple[float, float, Any, Optional[Any]]
        (delta_g_kcal, error_kcal, fitted_estimator, overlap_matrix)
    """
    combined = pd.concat(datasets)
    fitted = estimator_cls().fit(combined)

    delta_g_kj = float(fitted.delta_f_.iloc[0, -1])
    delta_g_kcal = delta_g_kj / 4.184

    if hasattr(fitted, "d_delta_f_"):
        error_kj = float(fitted.d_delta_f_.iloc[0, -1])
        error_kcal = error_kj / 4.184
    else:
        error_kcal = 0.0
        logger.warning("Error estimate not available for this estimator")

    # Extract overlap matrix for MBAR estimator
    overlap_matrix = None
    if hasattr(fitted, "overlap_matrix"):
        overlap_matrix = fitted.overlap_matrix_

    return delta_g_kcal, error_kcal, fitted, overlap_matrix


def compute_free_energy_gmx_bar(
    leg_dir: Path, leg_name: str, temperature: float, logger: logging.Logger
) -> Tuple[float, float]:
    """Run gmx bar and parse the reported ΔG/error in kcal/mol."""
    from prism.fep.analysis.xvg_parser import find_xvg_file

    xvg_files = []
    for window_dir in sorted(leg_dir.glob("window_*")):
        try:
            xvg_files.append(str(find_xvg_file(window_dir)))
        except FileNotFoundError:
            continue

    if not xvg_files:
        raise FileNotFoundError(f"No dhdl.xvg files found in {leg_dir}")

    logger.info(f"Found {len(xvg_files)} XVG files for {leg_name} leg")

    output_file = leg_dir / f"bar_{leg_name}.xvg"
    cmd = [
        "gmx",
        "bar",
        "-f",
        *xvg_files,
        "-o",
        str(output_file),
        "-temp",
        str(temperature),
        "-prec",
        "6",
    ]

    logger.info(f"Running: gmx bar -f ... -o {output_file}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout

        delta_g_kt = None
        error_kt = None
        for line in stdout.split("\n"):
            if line.strip().startswith("total"):
                match = re.search(r"DG\s+([+-]?\d+\.?\d*)\s*\+/-\s*([+-]?\d+\.?\d*)", line)
                if match:
                    delta_g_kt = float(match.group(1))
                    error_kt = float(match.group(2))
                    logger.info(f"Found total result: {delta_g_kt:.2f} ± {error_kt:.2f} kT")
                    break

        if delta_g_kt is None:
            for line in stdout.split("\n"):
                if "Delta G" in line or "estimate" in line:
                    match = re.search(r"([+-]?\d+\.?\d*)\s*\+/-\s*([+-]?\d+\.?\d*)", line)
                    if match:
                        delta_g_kt = float(match.group(1))
                        error_kt = float(match.group(2))
                        break

        if delta_g_kt is None and output_file.exists():
            with open(output_file) as handle:
                for line in handle:
                    if not line.startswith("#") and not line.startswith("@"):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                delta_g_kt = float(parts[0])
                                error_kt = float(parts[1])
                                break
                            except ValueError:
                                continue

        if delta_g_kt is None:
            raise ValueError("Could not parse gmx bar output")

        k_t_to_kj = (8.314 / 1000) * temperature
        delta_g_kcal = (delta_g_kt * k_t_to_kj) / 4.184
        error_kcal = ((error_kt or 0.0) * k_t_to_kj) / 4.184

        logger.info(f"gmx bar result: {delta_g_kcal:.2f} ± {error_kcal:.2f} kcal/mol")
        return delta_g_kcal, error_kcal
    except subprocess.CalledProcessError as exc:
        logger.error(f"gmx bar failed: {exc}")
        logger.error(f"stderr: {exc.stderr}")
        raise RuntimeError(f"gmx bar command failed: {exc}")


def summarize_repeat_results(
    bound_results: List[Tuple[float, float]],
    unbound_results: List[Tuple[float, float]],
) -> Dict[str, Any]:
    """Summarize per-repeat ΔG values and aggregate statistics."""
    repeat_results: List[Dict[str, float]] = []
    repeat_statistics: Dict[str, float] = {}
    n_repeats = 1

    if len(bound_results) != len(unbound_results):
        raise ValueError("Bound/unbound repeat counts do not match")

    delta_g_bound_values = [result[0] for result in bound_results]
    error_bound_values = [result[1] for result in bound_results]
    delta_g_unbound_values = [result[0] for result in unbound_results]
    error_unbound_values = [result[1] for result in unbound_results]

    delta_g_bound = float(pd.Series(delta_g_bound_values).mean())
    error_bound = float((pd.Series(error_bound_values).pow(2).sum() ** 0.5) / len(error_bound_values))
    delta_g_unbound = float(pd.Series(delta_g_unbound_values).mean())
    error_unbound = float((pd.Series(error_unbound_values).pow(2).sum() ** 0.5) / len(error_unbound_values))

    if len(bound_results) > 1:
        n_repeats = len(bound_results)
        for index, (bound_result, unbound_result) in enumerate(zip(bound_results, unbound_results), start=1):
            ddg = bound_result[0] - unbound_result[0]
            repeat_results.append(
                {
                    "repeat": index,
                    "bound": bound_result[0],
                    "unbound": unbound_result[0],
                    "ddG": ddg,
                }
            )

        ddg_values = pd.Series([result["ddG"] for result in repeat_results])
        bound_series = pd.Series(delta_g_bound_values)
        unbound_series = pd.Series(delta_g_unbound_values)
        repeat_statistics = {
            "bound_mean": float(bound_series.mean()),
            "bound_stderr": _compute_stderr(bound_series),
            "unbound_mean": float(unbound_series.mean()),
            "unbound_stderr": _compute_stderr(unbound_series),
            "ddG_mean": float(ddg_values.mean()),
            "ddG_stderr": _compute_stderr(ddg_values),
        }

    return {
        "delta_g_bound": delta_g_bound,
        "error_bound": error_bound,
        "delta_g_unbound": delta_g_unbound,
        "error_unbound": error_unbound,
        "repeat_results": repeat_results,
        "repeat_statistics": repeat_statistics,
        "n_repeats": n_repeats,
        "delta_g_bound_values": delta_g_bound_values,
        "delta_g_unbound_values": delta_g_unbound_values,
    }


def _compute_stderr(values: pd.Series) -> float:
    """Compute standard error for a series, returning 0.0 for singleton input."""
    if len(values) <= 1:
        return 0.0
    return float(values.std(ddof=1) / len(values) ** 0.5)


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
