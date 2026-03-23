#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Lambda-profile helpers for FEP analysis.
"""

from typing import Any, Dict, Optional


def extract_lambda_data(fitted_model: Any, estimator_name: str) -> Optional[Dict[str, Any]]:
    """Extract lambda-dependent profile data from a fitted estimator."""
    if fitted_model is None:
        return None

    lambda_data: Dict[str, Any] = {
        "lambdas": [],
        "dhdl": [],
        "cumulative_dg": [],
        "state_dg": [],
        "overlap_matrix": None,  # Only for MBAR
    }

    if estimator_name == "TI" and hasattr(fitted_model, "dhdl"):
        dhdl_df = fitted_model.dhdl
        n_states = len(dhdl_df)
        lambda_data["lambdas"] = list(range(n_states))
        lambda_data["dhdl"] = [float(v) for v in dhdl_df.sum(axis=1).values]
        if hasattr(fitted_model, "delta_f_"):
            lambda_data["cumulative_dg"] = [float(v) for v in fitted_model.delta_f_.iloc[0].values]
    elif hasattr(fitted_model, "delta_f_"):
        delta_f = fitted_model.delta_f_
        n_states = delta_f.shape[0]
        if n_states > 1:
            lambda_data["lambdas"] = [float(i) / (n_states - 1) for i in range(n_states)]
        else:
            lambda_data["lambdas"] = [0.0]

        for i in range(n_states):
            if i == 0:
                lambda_data["state_dg"].append(0.0)
            else:
                lambda_data["state_dg"].append(float(delta_f.iloc[0, i]))

        lambda_data["cumulative_dg"] = lambda_data["state_dg"].copy()
        lambda_data["dhdl"] = [0.0] * n_states

        # Extract overlap matrix for MBAR
        if hasattr(fitted_model, "overlap_matrix"):
            overlap_matrix = fitted_model.overlap_matrix
            if overlap_matrix is not None:
                # Convert to list for JSON serialization
                lambda_data["overlap_matrix"] = overlap_matrix.tolist()

    if not lambda_data["lambdas"]:
        return None

    return lambda_data


def build_lambda_profiles(fitted_bound: Any, fitted_unbound: Any, estimator_name: str) -> Dict[str, Any]:
    """Build bound/unbound lambda profiles for one or more repeats."""
    estimator = estimator_name.upper()

    if isinstance(fitted_bound, list):
        bound_profiles = [extract_lambda_data(fitted, estimator) for fitted in fitted_bound]
    else:
        bound_profiles = extract_lambda_data(fitted_bound, estimator)

    if isinstance(fitted_unbound, list):
        unbound_profiles = [extract_lambda_data(fitted, estimator) for fitted in fitted_unbound]
    else:
        unbound_profiles = extract_lambda_data(fitted_unbound, estimator)

    return {"bound": bound_profiles, "unbound": unbound_profiles}
