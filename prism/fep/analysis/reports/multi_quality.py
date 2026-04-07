#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quality check functions for FEP multi-estimator reports

Provides advisory quality warnings based on FEP best practices.
"""

from typing import Optional, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from ..core.models import FEResults


def check_standard_error_quality(results: "FEResults") -> Optional[str]:
    """Check if standard error exceeds threshold (>1.0 kcal/mol)"""
    error = results.delta_g_error
    if error > 1.0:
        return (
            f"⚠️ <strong>High Standard Error:</strong> ΔΔG error ({error:.2f} kcal/mol) "
            f"exceeds recommended threshold (1.0 kcal/mol). Consider extending simulation "
            f"time or increasing number of repeats."
        )
    return None


def check_endpoint_catastrophe(results: "FEResults") -> Optional[str]:
    """Detect potential endpoint catastrophe in free energy profile"""
    lambda_profiles = getattr(results, "lambda_profiles", None)
    if not lambda_profiles:
        return None

    # Check bound leg
    bound_data = lambda_profiles.get("bound", {})
    if isinstance(bound_data, list):
        bound_data = bound_data[0] if bound_data else {}

    state_dg = bound_data.get("state_dg")
    if not state_dg or len(state_dg) < 3:
        return None

    # Calculate mean absolute dG (exclude endpoints)
    mean_abs_dg = np.mean(np.abs(state_dg[1:-1]))
    if mean_abs_dg < 0.1:  # Avoid division issues
        return None

    # Check if endpoints are abnormally large
    if abs(state_dg[0]) > 2 * mean_abs_dg or abs(state_dg[-1]) > 2 * mean_abs_dg:
        return (
            f"⚠️ <strong>Potential Endpoint Catastrophe:</strong> Large energy changes "
            f"detected at λ=0 or λ=1. Consider adding more lambda windows near endpoints."
        )

    return None


def check_time_convergence_quality(results: "FEResults") -> Optional[str]:
    """Check if simulation has converged over time"""
    time_conv = getattr(results, "time_convergence", None)
    if not time_conv:
        return None

    bound_conv = time_conv.get("bound", [])
    unbound_conv = time_conv.get("unbound", [])

    if not bound_conv or not unbound_conv:
        return None

    # Compare 50% vs 100% for both legs
    def check_leg(conv_data, leg_name):
        if len(conv_data) < 2:
            return None
        mid_idx = len(conv_data) // 2
        mid_point = conv_data[mid_idx]
        end_point = conv_data[-1]
        drift = abs(end_point["dg"] - mid_point["dg"])
        if drift > 0.5:
            return f"{leg_name}: {drift:.2f} kcal/mol drift"
        return None

    bound_issue = check_leg(bound_conv, "Bound")
    unbound_issue = check_leg(unbound_conv, "Unbound")

    issues = [x for x in [bound_issue, unbound_issue] if x]
    if issues:
        return (
            f"⚠️ <strong>Poor Time Convergence:</strong> {'; '.join(issues)}. "
            f"System may not be fully equilibrated. Consider extending simulation time."
        )

    return None


def check_overlap_matrix_quality(results: "FEResults") -> Optional[str]:
    """Check MBAR overlap matrix quality (MBAR only)"""
    lambda_profiles = getattr(results, "lambda_profiles", None)
    if not lambda_profiles:
        return None

    bound_data = lambda_profiles.get("bound", {})
    if isinstance(bound_data, list):
        bound_data = bound_data[0] if bound_data else {}

    overlap_matrix = bound_data.get("overlap_matrix")
    if overlap_matrix is None:
        return None

    overlap = np.array(overlap_matrix)
    n_states = overlap.shape[0]

    if n_states < 2:
        return None

    # Check minimum adjacent overlap
    min_overlap = 1.0
    for i in range(n_states - 1):
        min_overlap = min(min_overlap, overlap[i, i + 1], overlap[i + 1, i])

    if min_overlap < 0.03:
        return (
            f"❌ <strong>Poor Overlap Matrix:</strong> Minimum adjacent state overlap "
            f"({min_overlap:.3f}) is critically low (&lt;0.03). MBAR results may be unreliable. "
            f"Add more lambda windows or extend sampling time."
        )
    elif min_overlap < 0.1:
        return (
            f"⚠️ <strong>Marginal Overlap Matrix:</strong> Minimum adjacent state overlap "
            f"({min_overlap:.3f}) is below recommended threshold (0.1). Results may be inaccurate. "
            f"Consider adding lambda windows or extending sampling."
        )

    return None


def build_quality_warnings_section(estimator_name: str, results: "FEResults") -> str:
    """
    Build quality warnings section (conditional display)

    Returns empty string if no issues detected.
    """
    warnings = []

    # Run all quality checks
    checks = [
        check_standard_error_quality(results),
        check_endpoint_catastrophe(results),
        check_time_convergence_quality(results),
    ]

    # Add overlap check only for MBAR
    if estimator_name == "MBAR":
        checks.append(check_overlap_matrix_quality(results))

    # Collect non-None warnings
    warnings = [w for w in checks if w is not None]

    if not warnings:
        return ""

    # Build warning list
    warning_items = "\n".join([f"            <li>{w}</li>" for w in warnings])

    return f"""
        <div class="card" style="border-left: 4px solid #ff9800;">
            <h2>⚠️ Quality Assessment</h2>
            <div class="alert" style="background: #fff3e0; border-left: 4px solid #ff9800; padding: 15px; margin: 10px 0;">
                <p><strong>The following potential issues were detected:</strong></p>
                <ul style="margin: 10px 0; padding-left: 20px;">
{warning_items}
                </ul>
                <p style="margin-top: 15px; font-size: 0.9em; color: #666;">
                    <em><strong>Note:</strong> These are advisory warnings based on common best practices.
                    Results may still be valid depending on your specific system and requirements.</em>
                </p>
            </div>
        </div>
        """
