#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP Convergence Analysis

Time convergence and bootstrap error estimation for FEP calculations.
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Any

logger = logging.getLogger(__name__)


def compute_time_convergence(
    bound_data_list: List[List],
    unbound_data_list: List[List],
    estimator_class: Any,
    n_points: int = 10,
) -> Optional[Dict]:
    """
    Compute time convergence by fitting ΔG on increasing data fractions.

    Parameters
    ----------
    bound_data_list : list of lists of DataFrames (one per repeat)
    unbound_data_list : list of lists of DataFrames (one per repeat)
    estimator_class : alchemlyb estimator class (TI, BAR, or MBAR)
    n_points : number of convergence points to compute

    Returns
    -------
    dict with 'bound' and 'unbound' lists of {frac, dg, err} dicts,
    or None if analysis fails.
    """
    try:
        import pandas as pd
    except ImportError:
        return None

    if not bound_data_list:
        return None

    try:
        result: Dict[str, List] = {"bound": [], "unbound": []}

        for leg_name, data_list in [("bound", bound_data_list), ("unbound", unbound_data_list)]:
            all_dfs = []
            for repeat_data in data_list:
                all_dfs.extend(repeat_data)

            if not all_dfs:
                continue

            combined = pd.concat(all_dfs)
            n_states = len(all_dfs)
            rows_per_state = len(combined) // n_states if n_states > 0 else len(combined)

            fracs = [i / n_points for i in range(1, n_points + 1)]
            for frac in fracs:
                try:
                    n_take = max(1, int(rows_per_state * frac))
                    sub_dfs = []
                    idx_start = 0
                    for _ in range(n_states):
                        chunk = combined.iloc[idx_start : idx_start + n_take]
                        sub_dfs.append(chunk)
                        idx_start += rows_per_state
                    sub_combined = pd.concat(sub_dfs)
                    fitted = estimator_class().fit(sub_combined)
                    dg_kcal = float(fitted.delta_f_.iloc[0, -1]) / 4.184
                    err_kcal = float(fitted.d_delta_f_.iloc[0, -1]) / 4.184 if hasattr(fitted, "d_delta_f_") else 0.0
                    result[leg_name].append({"frac": frac, "dg": dg_kcal, "err": err_kcal})
                except Exception as e:
                    logger.debug(f"Time convergence point {frac:.1f} failed for {leg_name}: {e}")

        return result if (result["bound"] or result["unbound"]) else None
    except Exception as e:
        logger.warning(f"Time convergence analysis failed: {e}")
        return None


def compute_bootstrap(
    bound_data_list: List[List],
    unbound_data_list: List[List],
    estimator_class: Any,
    n_bootstrap: int = 50,
    fraction: float = 0.8,
) -> Optional[Dict]:
    """
    Bootstrap error estimation by resampling frames per lambda window.

    Parameters
    ----------
    bound_data_list : list of lists of DataFrames (one per repeat)
    unbound_data_list : list of lists of DataFrames (one per repeat)
    estimator_class : alchemlyb estimator class
    n_bootstrap : number of bootstrap iterations
    fraction : fraction of frames to sample each iteration

    Returns
    -------
    dict with ddG_values, mean, stderr, std, or None if analysis fails.
    """
    try:
        import pandas as pd
    except ImportError:
        return None

    if not bound_data_list:
        return None

    try:
        bound_dgs: List[float] = []
        unbound_dgs: List[float] = []

        for data_list, dgs_out in [(bound_data_list, bound_dgs), (unbound_data_list, unbound_dgs)]:
            all_dfs = []
            for repeat_data in data_list:
                all_dfs.extend(repeat_data)

            if not all_dfs:
                continue

            combined = pd.concat(all_dfs)
            n_states = len(all_dfs)
            rows_per_state = len(combined) // n_states if n_states > 0 else len(combined)

            for _ in range(n_bootstrap):
                try:
                    sub_dfs = []
                    idx_start = 0
                    for _ in range(n_states):
                        chunk = combined.iloc[idx_start : idx_start + rows_per_state]
                        n_take = max(1, int(len(chunk) * fraction))
                        sub_dfs.append(chunk.sample(n=n_take))
                        idx_start += rows_per_state
                    sub_combined = pd.concat(sub_dfs).sort_index()
                    fitted = estimator_class().fit(sub_combined)
                    dg_kcal = float(fitted.delta_f_.iloc[0, -1]) / 4.184
                    dgs_out.append(dg_kcal)
                except Exception as e:
                    logger.debug(f"Bootstrap iteration failed: {e}")

        if not bound_dgs or not unbound_dgs:
            return None

        n = min(len(bound_dgs), len(unbound_dgs))
        ddG_values = [bound_dgs[i] - unbound_dgs[i] for i in range(n)]

        return {
            "bound_dgs": bound_dgs,
            "unbound_dgs": unbound_dgs,
            "ddG_values": ddG_values,
            "ddG_mean": float(np.mean(ddG_values)),
            "ddG_std": float(np.std(ddG_values, ddof=1)),
            "ddG_stderr": float(np.std(ddG_values, ddof=1) / np.sqrt(n)),
            "n_bootstrap": n,
        }
    except Exception as e:
        logger.warning(f"Bootstrap analysis failed: {e}")
        return None
