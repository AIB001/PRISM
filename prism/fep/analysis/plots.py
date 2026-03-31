#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Shared plotting functions for FEP analysis reports.

This module contains reusable plotting logic that can be used by both
single-estimator and multi-estimator report generators.
"""

from typing import Optional, Dict, Any, Tuple
import json as _json


def _axis_style(title: str) -> Dict[str, Any]:
    """Return a consistent Plotly axis style."""
    return {
        "title": {"text": f"<b>{title}</b>", "font": {"size": 16, "color": "#222"}},
        "showgrid": True,
        "gridcolor": "#e0e0e0",
        "zeroline": False,
        "tickfont": {"size": 14, "color": "#333"},
    }


def _base_layout(
    title: str,
    xaxis_title: str,
    yaxis_title: str,
    *,
    height: int,
    margin: Optional[Dict[str, int]] = None,
) -> Dict[str, Any]:
    """Return a shared Plotly layout with larger, bolder labels."""
    return {
        "title": {"text": f"<b>{title}</b>", "font": {"size": 16, "color": "#222"}},
        "xaxis": _axis_style(xaxis_title),
        "yaxis": _axis_style(yaxis_title),
        "paper_bgcolor": "rgba(0,0,0,0)",
        "plot_bgcolor": "rgba(255,255,255,0.8)",
        "margin": margin or {"l": 88, "r": 72, "t": 56, "b": 78},
        "hovermode": "x unified",
        "autosize": True,
        "height": height,
        "font": {
            "family": '-apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif',
            "size": 14,
            "color": "#333",
        },
    }


def build_lambda_plots_html(
    lambda_profiles: Optional[Dict[str, Any]],
    estimator_name: str,
    plot_suffix: str = "",
) -> Tuple[list, str]:
    """
    Build 4 separate lambda profile plot divs and Plotly JS.

    Returns 4 plot divs that can be arranged in a 2x2 grid:
    - dg-bound-plot{suffix}
    - dg-unbound-plot{suffix}
    - dhdl-bound-plot{suffix} (only for TI)
    - dhdl-unbound-plot{suffix} (only for TI)

    Parameters
    ----------
    lambda_profiles : dict or None
        As produced by FEPAnalyzer._build_lambda_profiles.
        Can contain single repeat or multiple repeats:
        - Single: {"bound": {...}, "unbound": {...}}
        - Multiple: {"bound": [repeat1, repeat2, ...], "unbound": [repeat1, repeat2, ...]}
    estimator_name : str
        Estimator name (TI, BAR, MBAR).
    plot_suffix : str, optional
        Suffix to add to plot IDs (for multi-estimator mode to avoid ID conflicts).

    Returns
    -------
    tuple[list, str]
        (list of 4 plot div HTML strings, combined script tag)
    """
    if not lambda_profiles:
        return (["<div></div>", "<div></div>", "<div></div>", "<div></div>"], "")

    bound_data = lambda_profiles.get("bound") or {}
    unbound_data = lambda_profiles.get("unbound") or {}

    # Support both single and multiple repeats
    bound_list = bound_data if isinstance(bound_data, list) else [bound_data]
    unbound_list = unbound_data if isinstance(unbound_data, list) else [unbound_data]

    is_ti = estimator_name.upper() == "TI"

    # Build unique plot IDs with suffix
    suffix = f"-{estimator_name}" if plot_suffix else ""
    bound_plot_id = f"dg-bound-plot{suffix}"
    unbound_plot_id = f"dg-unbound-plot{suffix}"
    dhdl_bound_plot_id = f"dhdl-bound-plot{suffix}"
    dhdl_unbound_plot_id = f"dhdl-unbound-plot{suffix}"

    # Build 4 separate plot divs
    plot_divs = [
        f'<div id="{bound_plot_id}" class="plot-container compact-plot" style="height:380px;min-height:380px;width:100%;overflow:hidden;"></div>',
        f'<div id="{unbound_plot_id}" class="plot-container compact-plot" style="height:380px;min-height:380px;width:100%;overflow:hidden;"></div>',
        f'<div id="{dhdl_bound_plot_id}" class="plot-container compact-plot" style="height:380px;min-height:380px;width:100%;overflow:hidden;"></div>'
        if is_ti
        else "<div></div>",
        f'<div id="{dhdl_unbound_plot_id}" class="plot-container compact-plot" style="height:380px;min-height:380px;width:100%;overflow:hidden;"></div>'
        if is_ti
        else "<div></div>",
    ]

    # Prepare data for Plotly
    bound_dg_traces = []
    unbound_dg_traces = []
    dhdl_bound_traces = []
    dhdl_unbound_traces = []

    for i, bound_profile in enumerate(bound_list):
        repeat_label = f"Repeat {i+1}" if len(bound_list) > 1 else "Bound"
        lambdas = bound_profile.get("lambdas", [])
        cumulative_dg = bound_profile.get("cumulative_dg", [])
        dhdl = bound_profile.get("dhdl", [])

        if lambdas and cumulative_dg:
            bound_dg_traces.append(
                {
                    "x": lambdas,
                    "y": cumulative_dg,
                    "mode": "lines+markers",
                    "name": repeat_label,
                    "line": {"width": 2},
                    "marker": {"size": 6},
                }
            )

        if is_ti and lambdas and dhdl:
            dhdl_bound_traces.append(
                {
                    "x": lambdas,
                    "y": dhdl,
                    "mode": "lines+markers",
                    "name": repeat_label,
                    "line": {"width": 2},
                    "marker": {"size": 6},
                }
            )

    for i, unbound_profile in enumerate(unbound_list):
        repeat_label = f"Repeat {i+1}" if len(unbound_list) > 1 else "Unbound"
        lambdas = unbound_profile.get("lambdas", [])
        cumulative_dg = unbound_profile.get("cumulative_dg", [])
        dhdl = unbound_profile.get("dhdl", [])

        if lambdas and cumulative_dg:
            unbound_dg_traces.append(
                {
                    "x": lambdas,
                    "y": cumulative_dg,
                    "mode": "lines+markers",
                    "name": repeat_label,
                    "line": {"width": 2},
                    "marker": {"size": 6},
                }
            )

        if is_ti and lambdas and dhdl:
            dhdl_unbound_traces.append(
                {
                    "x": lambdas,
                    "y": dhdl,
                    "mode": "lines+markers",
                    "name": repeat_label,
                    "line": {"width": 2},
                    "marker": {"size": 6},
                }
            )

    # Build Plotly layouts
    bound_layout = _base_layout(
        "Bound Leg ΔG vs λ",
        "λ",
        "ΔG (kcal/mol)",
        height=380,
        margin={"l": 90, "r": 70, "t": 56, "b": 86},
    )
    unbound_layout = _base_layout(
        "Unbound Leg ΔG vs λ",
        "λ",
        "ΔG (kcal/mol)",
        height=380,
        margin={"l": 90, "r": 70, "t": 56, "b": 86},
    )

    # Build script content
    script_content = f"""
        Plotly.newPlot('{bound_plot_id}', {_json.dumps(bound_dg_traces)}, {_json.dumps(bound_layout)}, {{responsive: true}});
        Plotly.newPlot('{unbound_plot_id}', {_json.dumps(unbound_dg_traces)}, {_json.dumps(unbound_layout)}, {{responsive: true}});
    """

    # Add dH/dλ plots for TI
    if is_ti:
        dhdl_bound_layout = {
            **_base_layout(
                "Bound Leg dH/dλ vs λ",
                "λ",
                "dH/dλ (kcal/mol)",
                height=380,
                margin={"l": 90, "r": 70, "t": 56, "b": 86},
            )
        }
        dhdl_unbound_layout = {
            **_base_layout(
                "Unbound Leg dH/dλ vs λ",
                "λ",
                "dH/dλ (kcal/mol)",
                height=380,
                margin={"l": 90, "r": 70, "t": 56, "b": 86},
            )
        }
        script_content += f"""
            Plotly.newPlot('{dhdl_bound_plot_id}', {_json.dumps(dhdl_bound_traces)}, {_json.dumps(dhdl_bound_layout)}, {{responsive: true}});
            Plotly.newPlot('{dhdl_unbound_plot_id}', {_json.dumps(dhdl_unbound_traces)}, {_json.dumps(dhdl_unbound_layout)}, {{responsive: true}});
        """

    return plot_divs, f"<script>{script_content}</script>"


def build_repeats_table_html(
    repeat_results: list,
    n_repeats: int,
    statistics: Optional[Dict] = None,
) -> str:
    """
    Build HTML table showing per-repeat ΔG values.

    Parameters
    ----------
    repeat_results : list
        List of dicts with keys: 'repeat', 'bound', 'unbound', 'ddG'.
    n_repeats : int
        Number of repeats.
    statistics : dict, optional
        Statistics dict with 'bound_mean', 'bound_stderr', etc.

    Returns
    -------
    str
        HTML table string.
    """
    rows = ""
    for r in repeat_results:
        rows += f"""
        <tr>
            <td>{r["repeat"]}</td>
            <td class="value">{r["bound"]:.2f}</td>
            <td class="value">{r["unbound"]:.2f}</td>
            <td class="value">{r["ddG"]:.2f}</td>
        </tr>
        """

    # Add statistics rows if available
    if statistics:
        rows += f"""
        <tr style="border-top: 2px solid #666; font-weight: bold; background: #e8f5e9;">
            <td>ave</td>
            <td class="value">{statistics.get("bound_mean", 0):.2f}</td>
            <td class="value">{statistics.get("unbound_mean", 0):.2f}</td>
            <td class="value">{statistics.get("ddG_mean", 0):.2f}</td>
        </tr>
        <tr style="font-weight: bold; background: #e8f5e9;">
            <td>stderr</td>
            <td class="value">{statistics.get("bound_stderr", 0):.2f}</td>
            <td class="value">{statistics.get("unbound_stderr", 0):.2f}</td>
            <td class="value">{statistics.get("ddG_stderr", 0):.2f}</td>
        </tr>
        """

    return f"""
    <table>
        <thead>
            <tr>
                <th>Repeats</th>
                <th>Bound ΔG (kcal/mol)</th>
                <th>Unbound ΔG (kcal/mol)</th>
                <th>ΔΔG (kcal/mol)</th>
            </tr>
        </thead>
        <tbody>
            {rows}
        </tbody>
    </table>
    <p style="margin-top: 15px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
        <strong>📐 Formula:</strong> ΔΔG = Bound - Unbound<br>
        <strong>📊 stderr:</strong> Standard error of the mean (σ/√n)<br>
        <strong>🔢 n_repeats:</strong> {n_repeats}
    </p>
    """


def build_repeats_outlier_plot_html(
    repeat_results: list,
    plot_suffix: str = "",
) -> str:
    """
    Build Plotly box plot showing distribution across repeats.

    Parameters
    ----------
    repeat_results : list
        List of dicts with keys: 'repeat', 'bound', 'unbound', 'ddG'.
    plot_suffix : str, optional
        Suffix for plot ID to avoid conflicts.

    Returns
    -------
    str
        HTML div and script for box plot.
    """
    if not repeat_results:
        return ""

    bound_values = [r["bound"] for r in repeat_results]
    unbound_values = [r["unbound"] for r in repeat_results]

    plot_id = f"repeats-outlier-{plot_suffix}" if plot_suffix else "repeats-outlier-plot"

    traces = []

    def create_box_scatter_group(values, name, color):
        return {
            "x": [name] * len(values),
            "y": values,
            "name": name,
            "type": "box",
            "marker": {"color": color},
            "boxpoints": "outliers",
            "jitter": 0.3,
            "pointpos": 0,
        }

    traces.append(create_box_scatter_group(bound_values, "Bound", "rgba(76,175,80,0.7)"))
    traces.append(create_box_scatter_group(unbound_values, "Unbound", "rgba(244,67,54,0.7)"))

    layout = {
        **_base_layout("Repeat Outlier Detection", "Measurement Type", "ΔG (kcal/mol)", height=340),
        "showlegend": False,
        "boxmode": "group",
    }

    return f"""
    <div id="{plot_id}" style="min-height:350px;"></div>
    <script>
    Plotly.newPlot('{plot_id}', {_json.dumps(traces)}, {_json.dumps(layout)}, {{responsive:true}});
    </script>
    <p style="margin-top: 10px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
        <strong>📊 Interpretation:</strong> Box plots show the distribution across repeats.
        Outliers (points outside 1.5×IQR) are displayed individually.
    </p>
    """


def build_overlap_matrix_html(
    lambda_profiles: Optional[Dict[str, Any]],
    estimator_name: str,
    plot_suffix: str = "",
) -> str:
    """
    Build overlap matrix visualization for bound and unbound legs (MBAR only).

    Parameters
    ----------
    lambda_profiles : dict or None
        Lambda profiles containing overlap_matrix for bound and unbound legs.
    estimator_name : str
        Estimator name (only MBAR has overlap matrix).
    plot_suffix : str, optional
        Suffix for plot ID.

    Returns
    -------
    str
        HTML div and script for overlap matrix heatmaps (bound + unbound).
    """
    if not lambda_profiles or estimator_name.upper() != "MBAR":
        return ""

    import numpy as np

    bound_data = lambda_profiles.get("bound", {})
    unbound_data = lambda_profiles.get("unbound", {})

    # Handle both single repeat and multiple repeats
    if isinstance(bound_data, list):
        bound_overlap = bound_data[0].get("overlap_matrix") if bound_data and bound_data[0] else None
    else:
        bound_overlap = bound_data.get("overlap_matrix")

    if isinstance(unbound_data, list):
        unbound_overlap = unbound_data[0].get("overlap_matrix") if unbound_data and unbound_data[0] else None
    else:
        unbound_overlap = unbound_data.get("overlap_matrix")

    if bound_overlap is None and unbound_overlap is None:
        return ""

    # Generate plots
    plots_html = []

    for leg_name, overlap_matrix in [("Bound", bound_overlap), ("Unbound", unbound_overlap)]:
        if overlap_matrix is None:
            continue

        overlap = np.array(overlap_matrix)
        n_states = overlap.shape[0]

        # Create Plotly heatmap
        heatmap_data = [
            {
                "z": overlap.tolist(),
                "x": list(range(n_states)),
                "y": list(range(n_states)),
                "colorscale": [
                    [0.0, "rgb(0,0,131)"],
                    [0.1, "rgb(0,0,255)"],
                    [0.3, "rgb(0,191,255)"],
                    [0.5, "rgb(0,255,255)"],
                    [0.7, "rgb(255,255,0)"],
                    [0.9, "rgb(255,127,0)"],
                    [1.0, "rgb(255,0,0)"],
                ],
                "colorbar": {
                    "title": {"text": "<b>Overlap</b>", "font": {"size": 15, "color": "#222"}},
                    "thickness": 15,
                },
                "type": "heatmap",
                "hovertemplate": "<b>%{x}</b> -> <b>%{y}</b><br>Overlap: %{z:.3f}<extra></extra>",
            }
        ]

        suffix_token = plot_suffix.lstrip("-")
        plot_id = (
            f"overlap-matrix-{leg_name.lower()}-{suffix_token}"
            if suffix_token
            else f"overlap-matrix-{leg_name.lower()}-plot"
        )

        layout = {
            **_base_layout(f"{leg_name} Leg - MBAR Overlap Matrix", "λ State", "λ State", height=520),
            "margin": {"l": 92, "r": 132, "t": 56, "b": 86},
            "hovermode": "closest",
            "yaxis": {
                **_axis_style("λ State"),
                "scaleanchor": "x",
                "constrain": "domain",
            },
            "xaxis": {
                **_axis_style("λ State"),
                "constrain": "domain",
            },
        }

        plots_html.append(f"""
        <div style="text-align:center;min-width:0;width:100%;">
            <div id="{plot_id}" class="plot-container overlap-plot" style="height:520px;min-height:520px;width:100%;"></div>
        </div>
        <script>
            Plotly.newPlot('{plot_id}', {_json.dumps(heatmap_data)}, {_json.dumps(layout)}, {{responsive:true}});
        </script>
        """)

    if not plots_html:
        return ""

    return f"""
    <div class="equal-grid two-col" style="margin-top:0;">
        {"".join(plots_html)}
    </div>
    <p style="margin-top: 10px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
        <strong>📊 Interpretation:</strong> Overlap matrix shows the phase space overlap between adjacent lambda windows.
        Values closer to 1.0 (red) indicate good overlap, while values near 0.0 (blue) suggest poor sampling.
        Ideally, all off-diagonal elements should be > 0.1 for reliable MBAR estimates.
    </p>
    """


def build_time_convergence_html(
    time_conv_data: Optional[Dict],
    plot_suffix: str = "",
) -> str:
    """
    Build time convergence plot HTML.

    Parameters
    ----------
    time_conv_data : dict or None
        Time convergence data with 'bound' and 'unbound' keys.
        Each value is a list of dicts with keys: 'frac', 'dg', 'err'.
    plot_suffix : str, optional
        Suffix for plot ID.

    Returns
    -------
    str
        HTML div and script for convergence plot.
    """
    if not time_conv_data:
        return """<div class="alert">
        <strong>ℹ️ Note:</strong> Time convergence analysis was not performed or data was insufficient.
    </div>"""

    bound_pts = time_conv_data.get("bound", [])
    unbound_pts = time_conv_data.get("unbound", [])

    if not bound_pts and not unbound_pts:
        return """<div class="alert">
        <strong>ℹ️ Note:</strong> Time convergence analysis yielded no data points.
    </div>"""

    suffix_token = plot_suffix.lstrip("-")
    plot_id = f"time-conv-{suffix_token}" if suffix_token else "time-conv-plot"

    traces = []
    if bound_pts:
        x = [p["frac"] for p in bound_pts]
        y = [p["dg"] for p in bound_pts]
        err = [p["err"] for p in bound_pts]
        traces.append(
            {
                "x": x,
                "y": y,
                "error_y": {"type": "data", "array": err, "visible": True},
                "name": "Bound",
                "mode": "lines+markers",
                "line": {"width": 3, "color": "rgb(76,175,80)"},
                "marker": {"size": 8, "symbol": "circle-open"},
            }
        )

    if unbound_pts:
        x = [p["frac"] for p in unbound_pts]
        y = [p["dg"] for p in unbound_pts]
        err = [p["err"] for p in unbound_pts]
        traces.append(
            {
                "x": x,
                "y": y,
                "error_y": {"type": "data", "array": err, "visible": True},
                "name": "Unbound",
                "mode": "lines+markers",
                "line": {"width": 3, "color": "rgb(244,67,54)"},
                "marker": {"size": 8, "symbol": "circle-open"},
            }
        )

    layout = _base_layout("Time Convergence of ΔG", "Fraction of simulation data", "ΔG (kcal/mol)", height=320)
    layout["margin"] = {"l": 88, "r": 118, "t": 56, "b": 78}
    layout["legend"] = {
        "x": 1.03,
        "y": 1.0,
        "xanchor": "left",
        "yanchor": "top",
        "bgcolor": "rgba(255,255,255,0.88)",
        "bordercolor": "rgba(0,0,0,0.08)",
        "borderwidth": 1,
    }

    return f"""
    <div id="{plot_id}" class="plot-container" style="min-height:320px;width:100%;"></div>
    <script>
    Plotly.newPlot('{plot_id}', {_json.dumps(traces)}, {_json.dumps(layout)}, {{responsive:true}});
    </script>
    <p style="margin-top: 10px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
        <strong>📈 Interpretation:</strong> A converged simulation shows a flat ΔG vs. simulation time.
        If ΔG is still changing at 100%, more sampling is needed.
    </p>
    """


def build_bootstrap_html(
    bootstrap_data: Optional[Dict],
    plot_suffix: str = "",
) -> str:
    """
    Build bootstrap analysis plot HTML with histogram and statistics table.

    Parameters
    ----------
    bootstrap_data : dict or None
        Bootstrap analysis data with keys: 'ddG_values', 'ddG_mean', 'ddG_std',
        'ddG_stderr', 'n_bootstrap'.
    plot_suffix : str, optional
        Suffix for plot ID.

    Returns
    -------
    str
        HTML div and script for bootstrap histogram and statistics table.
    """
    if not bootstrap_data:
        return """<div class="alert">
        <strong>ℹ️ Note:</strong> Bootstrap analysis was not performed or data was insufficient.
    </div>"""

    ddG_values = bootstrap_data.get("ddG_values", [])
    ddG_mean = bootstrap_data.get("ddG_mean", 0.0)
    ddG_std = bootstrap_data.get("ddG_std", 0.0)
    ddG_stderr = bootstrap_data.get("ddG_stderr", 0.0)
    n_bootstrap = bootstrap_data.get("n_bootstrap", 0)

    if not ddG_values:
        return """<div class="alert">
        <strong>ℹ️ Note:</strong> Bootstrap analysis yielded no data.
    </div>"""

    suffix_token = plot_suffix.lstrip("-")
    plot_id = f"bootstrap-hist-{suffix_token}" if suffix_token else "bootstrap-hist-plot"

    # Histogram trace
    hist_trace = {
        "x": ddG_values,
        "type": "histogram",
        "name": "ΔΔG bootstrap",
        "marker": {"color": "rgba(33,150,243,0.7)", "line": {"color": "rgb(33,150,243)", "width": 1}},
        "nbinsx": max(5, n_bootstrap // 5),
    }

    layout = {
        **_base_layout("Bootstrap ΔΔG Distribution", "ΔΔG (kcal/mol)", "Count", height=300),
        "shapes": [
            {
                "type": "line",
                "x0": ddG_mean,
                "x1": ddG_mean,
                "y0": 0,
                "y1": 1,
                "xref": "x",
                "yref": "paper",
                "line": {"color": "red", "width": 2, "dash": "dash"},
            }
        ],
    }

    return f"""
    <div id="{plot_id}" class="plot-container" style="min-height:300px;width:100%;"></div>
    <div style="margin-top:14px;min-width:0;width:100%;">
        <table>
            <thead><tr><th>Metric</th><th>Value</th></tr></thead>
            <tbody>
                <tr><td><strong>Mean ΔΔG</strong></td><td class="value">{ddG_mean:.3f} kcal/mol</td></tr>
                <tr><td><strong>Std Dev</strong></td><td class="value">{ddG_std:.3f} kcal/mol</td></tr>
                <tr><td><strong>Stderr</strong></td><td class="value">{ddG_stderr:.3f} kcal/mol</td></tr>
                <tr><td><strong>N Bootstrap</strong></td><td class="value">{n_bootstrap}</td></tr>
            </tbody>
        </table>
    </div>
    <script>
        Plotly.newPlot('{plot_id}', [{_json.dumps(hist_trace)}], {_json.dumps(layout)}, {{responsive:true}});
    </script>
    <p style="margin-top:10px;font-size:0.85em;color:#666;">
        Bootstrap resamples 80% of frames per window, {n_bootstrap} iterations.
        Red dashed line = mean.
    </p>
    """
