#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML Report Generator for FEP Analysis

Generates comprehensive, publication-ready HTML reports for FEP calculations.
Includes interactive plots, tables, and convergence diagnostics.

Modern UI design with static templates.
"""

from pathlib import Path
from typing import Optional, Dict, Any, Tuple, Union
from datetime import datetime


def _load_template(filename: str) -> str:
    """Load template file from templates directory."""
    template_dir = Path(__file__).parent / "templates"
    template_path = template_dir / filename

    if not template_path.exists():
        raise FileNotFoundError(f"Template file not found: {template_path}")

    with open(template_path, "r", encoding="utf-8") as f:
        return f.read()


class HTMLReportGenerator:
    """
    Generate HTML reports for FEP analysis results

    Creates comprehensive HTML reports with:
    - Summary of free energy results
    - Interactive plots (using Plotly)
    - Convergence diagnostics
    - Energy decomposition
    - Method parameters

    Parameters
    ----------
    results : FEResults
        Analysis results from FEPAnalyzer
    raw_data : dict, optional
        Raw data from analysis (contains fitted models for plotting)
    template_style : str, optional
        Style preference: 'modern' (default) or 'minimal'

    Examples
    --------
    >>> generator = HTMLReportGenerator(results)
    >>> html_path = generator.generate('report.html')
    """

    def __init__(self, results, raw_data=None, template_style: str = "modern", template_path: Optional[str] = None):
        """
        Initialize report generator

        Parameters
        ----------
        results : FEResults
            Analysis results
        raw_data : dict, optional
            Raw data for plotting
        template_style : str
            'modern' or 'minimal' style
        template_path : str, optional
            Custom template path (future feature)
        """
        self.results = results
        self.raw_data = raw_data or {}
        self.template_style = template_style
        self.template_path = template_path

    def generate(self, output_path: Union[str, Path]) -> Path:
        """
        Generate HTML report

        Parameters
        ----------
        output_path : str
            Path to output HTML file

        Returns
        -------
        Path
            Path to generated HTML file
        """
        html_content = self._generate_html()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            f.write(html_content)

        return output_path

    def _generate_html(self) -> str:
        """Generate complete HTML document"""

        # Load templates
        css_content = _load_template("styles.css")
        js_content = _load_template("script.js")

        # Extract result data
        results_dict = {
            "backend": getattr(self.results, "backend", "unknown"),
            "estimator": getattr(self.results, "estimator", "unknown"),
            "delta_g": getattr(self.results, "delta_g", 0.0),
            "delta_g_error": getattr(self.results, "delta_g_error", float("nan")),
            "delta_g_bound": getattr(self.results, "delta_g_bound", 0.0),
            "delta_g_unbound": getattr(self.results, "delta_g_unbound", 0.0),
        }

        results_dict["backend"] = self.results.metadata.get("backend", results_dict["backend"])
        results_dict["estimator"] = self.results.metadata.get("estimator", results_dict["estimator"])

        # Extract metadata
        temperature = self.results.metadata.get("temperature", 310.0)
        pressure = self.results.metadata.get("pressure", 1.0)
        n_lambda_windows = self.results.metadata.get("n_lambda_windows", 2)

        # Create lambda schedule for display (just show range)
        lambda_schedule = (
            [0.0, 1.0] if n_lambda_windows == 2 else [i / (n_lambda_windows - 1) for i in range(n_lambda_windows)]
        )

        # Build lambda profile plots from raw_data
        lambda_profiles = self.raw_data.get("lambda_profiles")
        estimator_name = results_dict.get("estimator", "MBAR")
        dhdl_divs, lambda_plot_script = self._build_lambda_plots_html(lambda_profiles, estimator_name)

        # Build convergence and bootstrap plots
        time_conv_html = self._build_time_convergence_html(self.raw_data.get("time_convergence"))
        bootstrap_html = self._build_bootstrap_html(self.raw_data.get("bootstrap"))
        repeats_html = self._build_repeats_table_html()
        repeats_outlier_html = self._build_repeats_outlier_plot_html()
        detailed_results_section = ""
        if repeats_html:
            detailed_results_section = f"""
        <div class="card">
            <h2>📋 Detailed Results</h2>
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;align-items:start;">
                <div>
                    <h3 style="color: #333; margin-bottom: 15px; font-size: 1.2em;">
                        📊 Multiple Repeats Statistics
                    </h3>
{repeats_html}
                </div>
                <div>
                    <h3 style="color: #333; margin-bottom: 15px; font-size: 1.2em;">
                        🎯 Outlier Detection
                    </h3>
{repeats_outlier_html}
                </div>
            </div>
        </div>
"""

        # Build HTML
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PRISM-FEP Analysis Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
{css_content}
    </style>
    <script>
{js_content}
    </script>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧪 FEP Analysis Report</h1>
            <div class="subtitle">
                Generated on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
            </div>
        </header>

        <div class="card">
            <h2>📊 Summary</h2>
            <div class="summary-grid">
                <div class="summary-item" style="grid-column: span 2; background: linear-gradient(135deg, #667eea15 0%, #764ba215 100%);">
                    <div class="label">ΔΔG Binding (kcal/mol)</div>
                    <div class="value">{results_dict["delta_g"]:.2f} ± {results_dict["delta_g_error"]:.2f}</div>
                    <div style="font-size: 0.8em; color: #666;">
                        ΔG_bound - ΔG_unbound = {results_dict["delta_g_bound"]:.2f} - {
            results_dict["delta_g_unbound"]:.2f}
                    </div>
                </div>
                <div class="summary-item">
                    <div class="label">Temperature</div>
                    <div class="value">{temperature:.1f} K</div>
                    <div style="font-size: 0.8em; color: #666;">{temperature - 273.15:.1f} °C</div>
                </div>
                <div class="summary-item">
                    <div class="label">Lambda Windows</div>
                    <div class="value">{n_lambda_windows}</div>
                    <div style="font-size: 0.8em; color: #666;">states</div>
                </div>
            </div>
            <div style="margin-top: 20px; padding: 15px; background: #f8f9fa; border-radius: 4px;">
                <strong>Method:</strong> {results_dict["backend"]} with {results_dict["estimator"]} estimator
            </div>
        </div>

        <div class="card">
            <h2>🔧 Simulation Parameters</h2>
            <table>
                <thead>
                    <tr>
                        <th>Parameter</th>
                        <th>Value</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>Temperature</strong></td>
                        <td class="value">{temperature:.1f} K ({temperature - 273.15:.1f} °C)</td>
                    </tr>
                    <tr>
                        <td><strong>Pressure</strong></td>
                        <td class="value">{pressure:.1f} bar</td>
                    </tr>
                    <tr>
                        <td><strong>Lambda Windows</strong></td>
                        <td class="value">{n_lambda_windows}</td>
                    </tr>
                    <tr>
                        <td><strong>Lambda Range</strong></td>
                        <td class="value">{lambda_schedule[0]:.3f} → {lambda_schedule[-1]:.3f}</td>
                    </tr>
                    <tr>
                        <td><strong>Estimator</strong></td>
                        <td class="value">{results_dict["estimator"]}</td>
                    </tr>
                </tbody>
            </table>
        </div>

{detailed_results_section}

        <div class="card">
            <h2>📈 Free Energy Profiles (λ-dependent)</h2>
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:16px;">
                <div class="plot-container" id="dg-bound-plot" style="min-height:300px;"></div>
                <div class="plot-container" id="dg-unbound-plot" style="min-height:300px;"></div>
            </div>
{dhdl_divs}
        </div>
{lambda_plot_script}

        <div class="card">
            <h2>📊 Convergence Diagnostics</h2>
{time_conv_html}
        </div>

        <div class="card">
            <h2>🔄 Bootstrap Analysis</h2>
{bootstrap_html}
        </div>

        <div class="card">
            <h2>📊 Overlap Matrix (MBAR)</h2>
            {self._build_overlap_matrix_html() or '<div class="alert"><strong>ℹ️ Note:</strong> Overlap matrix is only available for MBAR estimator.<br>Switch to MBAR backend to enable overlap matrix visualization.</div>'}
        </div>

        <div class="card">
            <h2>⚡ Energy Decomposition</h2>
            <div class="alert">
                <strong>📋 Note:</strong> Energy decomposition requires component-wise dH/dλ data.
                To enable energy decomposition, set <code>free-energy-components = yes</code> in MDP file.
            </div>
        </div>

        <div class="footer">
            <p>Generated by PRISM FEP Analyzer</p>
            <p>{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        </div>
    </div>
</body>
</html>
"""
        return html

    def _build_lambda_plots_html(
        self,
        lambda_profiles: Optional[Dict[str, Any]],
        estimator_name: str,
    ) -> Tuple[str, str]:
        """
        Build HTML div placeholders and Plotly JS for lambda-dependent plots.

        Parameters
        ----------
        lambda_profiles : dict or None
            As produced by FEPAnalyzer._build_lambda_profiles.
            Can contain single repeat or multiple repeats:
            - Single: {"bound": {...}, "unbound": {...}}
            - Multiple: {"bound": [repeat1, repeat2, ...], "unbound": [repeat1, repeat2, ...]}
        estimator_name : str
            Estimator name (TI, BAR, MBAR).

        Returns
        -------
        tuple[str, str]
            (dhdl_divs_html, script_tag_html)
        """
        import json as _json

        if not lambda_profiles:
            return "", ""

        bound_data = lambda_profiles.get("bound") or {}
        unbound_data = lambda_profiles.get("unbound") or {}

        # Support both single and multiple repeats
        # Check if data is a list (multiple repeats) or dict (single repeat)
        bound_list = bound_data if isinstance(bound_data, list) else [bound_data]
        unbound_list = unbound_data if isinstance(unbound_data, list) else [unbound_data]

        is_ti = estimator_name.upper() == "TI"

        # dH/dλ divs only for TI
        dhdl_divs = ""
        if is_ti:
            dhdl_divs = """
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-top:16px;">
                <div class="plot-container" id="dhdl-bound-plot" style="min-height:300px;"></div>
                <div class="plot-container" id="dhdl-unbound-plot" style="min-height:300px;"></div>
            </div>"""

        layout_common = {
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(255,255,255,0.8)",
            "margin": {"l": 60, "r": 30, "t": 45, "b": 50},
        }

        script_lines = ["<script>"]

        # Colors for multiple repeats
        repeat_colors = [
            "rgb(76,175,80)",  # Green
            "rgb(33,150,243)",  # Blue
            "rgb(255,152,0)",  # Orange
            "rgb(156,39,176)",  # Purple
            "rgb(0,150,136)",  # Teal
        ]

        def _dg_trace(data: dict, name: str, color: str) -> dict:
            return {
                "x": data.get("lambdas", []),
                "y": data.get("cumulative_dg", []),
                "name": name,
                "mode": "lines+markers",
                "line": {"width": 3, "color": color},
                "marker": {"size": 6},
            }

        def _dhdl_trace(data: dict, name: str, color: str) -> dict:
            return {
                "x": data.get("lambdas", []),
                "y": data.get("dhdl", []),
                "name": name,
                "mode": "lines+markers",
                "line": {"width": 2, "color": color},
                "marker": {"size": 5},
            }

        # Create traces for all repeats
        dg_bound_traces = []
        dg_unbound_traces = []
        dhdl_bound_traces = []
        dhdl_unbound_traces = []

        for i, bound in enumerate(bound_list):
            if not bound:
                continue
            repeat_name = f"Repeat {i + 1}" if len(bound_list) > 1 else "Bound Leg"
            color = repeat_colors[i % len(repeat_colors)]
            dg_bound_traces.append(_dg_trace(bound, repeat_name, color))
            if is_ti:
                dhdl_bound_traces.append(_dhdl_trace(bound, repeat_name, color))

        for i, unbound in enumerate(unbound_list):
            if not unbound:
                continue
            repeat_name = f"Repeat {i + 1}" if len(unbound_list) > 1 else "Unbound Leg"
            color = repeat_colors[i % len(repeat_colors)]
            dg_unbound_traces.append(_dg_trace(unbound, repeat_name, color))
            if is_ti:
                dhdl_unbound_traces.append(_dhdl_trace(unbound, repeat_name, color))

        bound_layout = {
            **layout_common,
            "title": {"text": "Bound Leg ΔG vs λ"},
            "xaxis": {"title": "λ"},
            "yaxis": {"title": "ΔG (kcal/mol)"},
        }
        unbound_layout = {
            **layout_common,
            "title": {"text": "Unbound Leg ΔG vs λ"},
            "xaxis": {"title": "λ"},
            "yaxis": {"title": "ΔG (kcal/mol)"},
        }
        bound_dhdl_layout = {
            **layout_common,
            "title": {"text": "Bound Leg dH/dλ vs λ"},
            "xaxis": {"title": "λ"},
            "yaxis": {"title": "dH/dλ (kcal/mol)"},
        }
        unbound_dhdl_layout = {
            **layout_common,
            "title": {"text": "Unbound Leg dH/dλ vs λ"},
            "xaxis": {"title": "λ"},
            "yaxis": {"title": "dH/dλ (kcal/mol)"},
        }

        # Plot ΔG vs λ
        if dg_bound_traces and dg_bound_traces[0]["x"]:
            bl = _json.dumps(bound_layout)
            script_lines.append(
                f"Plotly.newPlot('dg-bound-plot', {_json.dumps(dg_bound_traces)}, {bl}, {{responsive:true}});"
            )
        if dg_unbound_traces and dg_unbound_traces[0]["x"]:
            ul = _json.dumps(unbound_layout)
            script_lines.append(
                f"Plotly.newPlot('dg-unbound-plot', {_json.dumps(dg_unbound_traces)}, {ul}, {{responsive:true}});"
            )

        # Plot dH/dλ vs λ (TI only)
        if is_ti:
            if dhdl_bound_traces and dhdl_bound_traces[0]["x"]:
                bhl = _json.dumps(bound_dhdl_layout)
                script_lines.append(
                    f"Plotly.newPlot('dhdl-bound-plot', {_json.dumps(dhdl_bound_traces)}, {bhl}, {{responsive:true}});"
                )
            if dhdl_unbound_traces and dhdl_unbound_traces[0]["x"]:
                uhl = _json.dumps(unbound_dhdl_layout)
                script_lines.append(
                    f"Plotly.newPlot('dhdl-unbound-plot', {_json.dumps(dhdl_unbound_traces)}, {uhl}, {{responsive:true}});"
                )

        script_lines.append("</script>")
        return dhdl_divs, "\n".join(script_lines)

    def _build_repeats_table_html(self) -> str:
        """
        Build HTML table for multiple repeats statistics

        Returns
        -------
        str
            HTML table string or empty string if single repeat
        """
        # Defensive checks
        if not hasattr(self.results, "n_repeats") or self.results.n_repeats <= 1:
            return ""

        if not hasattr(self.results, "repeat_results") or not self.results.repeat_results:
            return ""

        # Build table rows for each repeat
        rows_html = ""
        for r in self.results.repeat_results:
            rows_html += f"""
                    <tr>
                        <td>{r["repeat"]}</td>
                        <td class="value">{r["bound"]:.2f}</td>
                        <td class="value">{r["unbound"]:.2f}</td>
                        <td class="value">{r["ddG"]:.2f}</td>
                    </tr>
        """

        # Add statistics rows (ave and stderr)
        stats = self.results.metadata.get("repeat_statistics", {})
        if stats:
            rows_html += f"""
                    <tr style="border-top: 2px solid #666; font-weight: bold; background: #e8f5e9;">
                        <td>ave</td>
                        <td class="value">{stats.get("bound_mean", 0):.2f}</td>
                        <td class="value">{stats.get("unbound_mean", 0):.2f}</td>
                        <td class="value">{stats.get("ddG_mean", 0):.2f}</td>
                    </tr>
                    <tr style="font-weight: bold; background: #e8f5e9;">
                        <td>stderr</td>
                        <td class="value">{stats.get("bound_stderr", 0):.2f}</td>
                        <td class="value">{stats.get("unbound_stderr", 0):.2f}</td>
                        <td class="value">{stats.get("ddG_stderr", 0):.2f}</td>
                    </tr>
        """

        # Build complete table HTML (without outer wrapper for side-by-side layout)
        table_html = f"""
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
{rows_html}
                    </tbody>
                </table>
                <p style="margin-top: 15px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
                    <strong>📐 Formula:</strong> ΔΔG = Bound - Unbound<br>
                    <strong>📊 stderr:</strong> Standard error of the mean (σ/√n)<br>
                    <strong>🔢 n_repeats:</strong> {self.results.n_repeats}
                </p>
        """

        return table_html

    def _build_repeats_outlier_plot_html(self) -> str:
        """
        Build HTML for outlier detection plot across multiple repeats.

        Creates a box plot with individual data points overlaid,
        similar to FEbuilder's dg-fepout.png visualization.

        Returns
        -------
        str
            HTML content with plotly scatter plot (empty string if single repeat).
        """
        # Defensive checks
        if not hasattr(self.results, "n_repeats") or self.results.n_repeats <= 1:
            return ""

        if not hasattr(self.results, "repeat_results") or not self.results.repeat_results:
            return ""

        import json as _json

        # Extract data from repeat_results (exclude ddG from outlier plot)
        bound_values = [r["bound"] for r in self.results.repeat_results]
        unbound_values = [r["unbound"] for r in self.results.repeat_results]

        # Build traces for box plots
        traces = []

        # Helper to create a box plot trace
        def create_box_scatter_group(values, name, color):
            box_trace = {
                "x": [name] * len(values),
                "y": values,
                "name": name,
                "type": "box",
                "marker": {"color": color},
                "boxpoints": "outliers",
                "jitter": 0.3,
                "pointpos": 0,
            }
            traces.append(box_trace)

        # Create box plots for bound and unbound only
        create_box_scatter_group(bound_values, "Bound", "rgba(76,175,80,0.7)")
        create_box_scatter_group(unbound_values, "Unbound", "rgba(244,67,54,0.7)")

        layout = {
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(255,255,255,0.8)",
            "margin": {"l": 70, "r": 30, "t": 45, "b": 60},
            "title": {"text": "Repeat Outlier Detection"},
            "xaxis": {"title": "Measurement Type"},
            "yaxis": {"title": "ΔG (kcal/mol)"},
            "showlegend": False,
            "boxmode": "group",
        }

        return f"""
            <div id="repeats-outlier-plot" style="min-height:350px;"></div>
            <script>
            Plotly.newPlot('repeats-outlier-plot', {_json.dumps(traces)}, {_json.dumps(layout)}, {{responsive:true}});
            </script>
            <p style="margin-top: 10px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
                <strong>📊 Interpretation:</strong> Box plots show the distribution across repeats.
                Outliers (points outside 1.5×IQR) are displayed individually.
            </p>
        """

    def _build_time_convergence_html(self, time_conv_data: Optional[Dict]) -> str:
        """
        Build HTML for time convergence analysis section.

        Parameters
        ----------
        time_conv_data : dict or None
            As returned by FEPAnalyzer._compute_time_convergence.
            Contains 'bound' and 'unbound' lists of {frac, dg, err} dicts.

        Returns
        -------
        str
            HTML content (empty string if no data).
        """
        if not time_conv_data:
            return """            <div class="alert">
                <strong>ℹ️ Note:</strong> Time convergence analysis was not performed or data was insufficient.
            </div>"""

        import json as _json

        bound_pts = time_conv_data.get("bound", [])
        unbound_pts = time_conv_data.get("unbound", [])

        if not bound_pts and not unbound_pts:
            return """            <div class="alert">
                <strong>ℹ️ Note:</strong> Time convergence analysis yielded no data points.
            </div>"""

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

        layout = {
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(255,255,255,0.8)",
            "margin": {"l": 70, "r": 30, "t": 45, "b": 60},
            "title": {"text": "Time Convergence of ΔG"},
            "xaxis": {"title": "Fraction of simulation data"},
            "yaxis": {"title": "ΔG (kcal/mol)"},
        }

        return f"""
            <div id="time-conv-plot" style="min-height:350px;"></div>
            <script>
            Plotly.newPlot('time-conv-plot', {_json.dumps(traces)}, {_json.dumps(layout)}, {{responsive:true}});
            </script>
            <p style="margin-top: 10px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
                <strong>📊 Interpretation:</strong> A converged simulation shows a flat ΔG vs. simulation time.
                If ΔG is still changing at 100%, more sampling is needed.
            </p>"""

    def _build_bootstrap_html(self, bootstrap_data: Optional[Dict]) -> str:
        """
        Build HTML for bootstrap error estimation section.

        Parameters
        ----------
        bootstrap_data : dict or None
            As returned by FEPAnalyzer._compute_bootstrap.

        Returns
        -------
        str
            HTML content.
        """
        if not bootstrap_data:
            return """            <div class="alert">
                <strong>ℹ️ Note:</strong> Bootstrap analysis was not performed or data was insufficient.
            </div>"""

        import json as _json

        ddG_values = bootstrap_data.get("ddG_values", [])
        ddG_mean = bootstrap_data.get("ddG_mean", 0.0)
        ddG_std = bootstrap_data.get("ddG_std", 0.0)
        ddG_stderr = bootstrap_data.get("ddG_stderr", 0.0)
        n_bootstrap = bootstrap_data.get("n_bootstrap", 0)

        if not ddG_values:
            return """            <div class="alert">
                <strong>ℹ️ Note:</strong> Bootstrap analysis yielded no data.
            </div>"""

        # Histogram trace
        hist_trace = {
            "x": ddG_values,
            "type": "histogram",
            "name": "ΔΔG bootstrap",
            "marker": {"color": "rgba(33,150,243,0.7)", "line": {"color": "rgb(33,150,243)", "width": 1}},
            "nbinsx": max(5, n_bootstrap // 5),
        }
        layout = {
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(255,255,255,0.8)",
            "margin": {"l": 70, "r": 30, "t": 45, "b": 60},
            "title": {"text": "Bootstrap ΔΔG Distribution"},
            "xaxis": {"title": "ΔΔG (kcal/mol)"},
            "yaxis": {"title": "Count"},
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
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;align-items:center;">
                <div id="bootstrap-hist-plot" style="min-height:300px;"></div>
                <div>
                    <table>
                        <thead><tr><th>Metric</th><th>Value</th></tr></thead>
                        <tbody>
                            <tr><td><strong>Mean ΔΔG</strong></td><td class="value">{ddG_mean:.3f} kcal/mol</td></tr>
                            <tr><td><strong>Std Dev</strong></td><td class="value">{ddG_std:.3f} kcal/mol</td></tr>
                            <tr><td><strong>Stderr</strong></td><td class="value">{ddG_stderr:.3f} kcal/mol</td></tr>
                            <tr><td><strong>N Bootstrap</strong></td><td class="value">{n_bootstrap}</td></tr>
                        </tbody>
                    </table>
                    <p style="margin-top:10px;font-size:0.85em;color:#666;">
                        Bootstrap resamples 80% of frames per window, {n_bootstrap} iterations.
                        Red dashed line = mean.
                    </p>
                </div>
            </div>
            <script>
            Plotly.newPlot('bootstrap-hist-plot', [{_json.dumps(hist_trace)}], {_json.dumps(layout)}, {{responsive:true}});
            </script>"""

    def _build_overlap_matrix_html(self) -> str:
        """
        Build HTML for MBAR overlap matrix visualization.

        Uses matplotlib to generate the overlap matrix plot and embeds it as base64 image.

        Returns
        -------
        str
            HTML content with overlap matrix plots (empty string if no overlap data).
        """
        bound_overlap = self.raw_data.get("overlap_matrix_bound")
        unbound_overlap = self.raw_data.get("overlap_matrix_unbound")

        if bound_overlap is None and unbound_overlap is None:
            return ""

        try:
            import io
            import base64
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from alchemlyb.visualisation import plot_mbar_overlap_matrix
        except ImportError:
            return """            <div class="alert">
                <strong>⚠️ Note:</strong> Overlap matrix visualization requires matplotlib.
                Install with: pip install matplotlib
            </div>"""

        # Generate plots and convert to base64
        images_html = []

        for leg_name, overlap_matrix in [("Bound", bound_overlap), ("Unbound", unbound_overlap)]:
            if overlap_matrix is None:
                continue

            try:
                fig, ax = plt.subplots(figsize=(6, 5))
                plot_mbar_overlap_matrix(overlap_matrix, ax=ax)
                ax.set_title(f"{leg_name} Leg Overlap Matrix")

                # Convert to base64
                buf = io.BytesIO()
                fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
                buf.seek(0)
                img_base64 = base64.b64encode(buf.read()).decode("utf-8")
                plt.close(fig)

                images_html.append(f"""
                    <div style="text-align:center;">
                        <img src="data:image/png;base64,{img_base64}"
                             alt="{leg_name} overlap matrix"
                             style="max-width:100%;height:auto;border-radius:4px;box-shadow:0 2px 4px rgba(0,0,0,0.1);">
                    </div>
                """)
            except Exception as exc:
                images_html.append(f"""
                    <div class="alert">
                        <strong>⚠️ Error:</strong> Could not generate {leg_name} overlap plot: {exc}
                    </div>
                """)

        if not images_html:
            return ""

        return f"""
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-top:20px;">
                {"".join(images_html)}
            </div>
            <p style="margin-top: 10px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
                <strong>📊 Interpretation:</strong> Overlap matrix shows the phase space overlap between adjacent lambda windows.
                Values closer to 1.0 (red) indicate good overlap, while values near 0.0 (blue) suggest poor sampling.
                Ideally, all off-diagonal elements should be > 0.1 for reliable MBAR estimates.
            </p>
        """


class MultiBackendReportGenerator:
    """
    Generate HTML reports comparing multiple backends

    Creates comprehensive HTML reports with:
    - Tab-based backend switching
    - Comparison plots
    - Per-backend detailed analysis

    Parameters
    ----------
    all_results : list of dict
        List of result dictionaries from multiple backends
    template_style : str, optional
        Style preference: 'modern' (default) or 'minimal'

    Examples
    --------
    >>> results = [
    ...     {'backend': 'alchemlyb', 'estimator': 'TI', 'delta_g': -10.5, ...},
    ...     {'backend': 'alchemlyb', 'estimator': 'BAR', 'delta_g': -11.2, ...},
    ... ]
    >>> generator = MultiBackendReportGenerator(results)
    >>> html_path = generator.generate('comparison_report.html')
    """

    def __init__(self, all_results: list, template_style: str = "modern"):
        """
        Initialize multi-backend report generator

        Parameters
        ----------
        all_results : list of dict
            Results from multiple backends
        template_style : str
            'modern' or 'minimal' style
        """
        self.all_results = all_results
        self.template_style = template_style

    def generate(self, output_path: Union[str, Path]) -> Path:
        """
        Generate HTML report

        Parameters
        ----------
        output_path : str
            Path to output HTML file

        Returns
        -------
        Path
            Path to generated HTML file
        """
        html_content = self._generate_html()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            f.write(html_content)

        return output_path

    def _generate_html(self) -> str:
        """Generate complete HTML document with backend tabs"""

        # Load templates
        css_content = _load_template("styles.css")
        js_content = _load_template("script.js")

        # Build comparison table
        comparison_rows = ""
        for result in self.all_results:
            error_display = (
                f"± {result['delta_g_error']:.2f}" if result["delta_g_error"] == result["delta_g_error"] else "N/A"
            )
            comparison_rows += f"""
                    <tr>
                        <td>{result["name"]}</td>
                        <td class="value">{result["delta_g"]:.2f}</td>
                        <td class="value">{error_display}</td>
                        <td class="value">{result["delta_g_bound"]:.2f}</td>
                        <td class="value">{result["delta_g_unbound"]:.2f}</td>
                    </tr>
            """

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>FEP Analysis - Multi-Backend Comparison</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
{css_content}
    </style>
    <script>
{js_content}
    </script>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧪 FEP Analysis - Multi-Backend Comparison</h1>
            <div class="subtitle">
                Generated on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
            </div>
        </header>

        <div class="card">
            <h2>📊 Results Comparison</h2>
            <table>
                <thead>
                    <tr>
                        <th>Backend / Estimator</th>
                        <th>ΔG (kcal/mol)</th>
                        <th>Error (kcal/mol)</th>
                        <th>Bound Leg</th>
                        <th>Unbound Leg</th>
                    </tr>
                </thead>
                <tbody>
{comparison_rows}
                </tbody>
            </table>
        </div>

        <div class="footer">
            <p>Generated by PRISM FEP Analyzer</p>
            <p>{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        </div>
    </div>
</body>
</html>
"""
        return html
