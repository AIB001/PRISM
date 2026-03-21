#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML Report Generator for FEP Analysis

Generates comprehensive, publication-ready HTML reports for FEP calculations.
Includes interactive plots, tables, and convergence diagnostics.

Modern UI design with static templates.
"""

from pathlib import Path
from typing import Optional, Dict, Any, Tuple
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

    def generate(self, output_path: str) -> Path:
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
        temperature = getattr(self.results, "temperature", 310.0)
        pressure = getattr(self.results, "pressure", 1.0)
        lambda_schedule = getattr(self.results, "lambda_schedule", [0.0, 1.0])

        # Build lambda profile plots from raw_data
        lambda_profiles = self.raw_data.get("lambda_profiles")
        estimator_name = results_dict.get("estimator", "MBAR")
        dhdl_divs, lambda_plot_script = self._build_lambda_plots_html(lambda_profiles, estimator_name)

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
                Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            </div>
        </header>

        <div class="card">
            <h2>📊 Summary</h2>
            <div class="summary-grid">
                <div class="summary-item">
                    <div class="label">Binding ΔG</div>
                    <div class="value">{results_dict['delta_g']:.2f} kcal/mol</div>
                    <div style="font-size: 0.8em; color: #666;">
                        {f"± {results_dict['delta_g_error']:.2f}" if results_dict['delta_g_error'] == results_dict['delta_g_error'] else "N/A"}
                    </div>
                </div>
                <div class="summary-item">
                    <div class="label">Bound Leg</div>
                    <div class="value">{results_dict['delta_g_bound']:.2f}</div>
                    <div style="font-size: 0.8em; color: #666;">kcal/mol</div>
                </div>
                <div class="summary-item">
                    <div class="label">Unbound Leg</div>
                    <div class="value">{results_dict['delta_g_unbound']:.2f}</div>
                    <div style="font-size: 0.8em; color: #666;">kcal/mol</div>
                </div>
                <div class="summary-item">
                    <div class="label">Temperature</div>
                    <div class="value">{temperature:.1f} K</div>
                    <div style="font-size: 0.8em; color: #666;">{temperature - 273.15:.1f} °C</div>
                </div>
            </div>
            <div style="margin-top: 20px; padding: 15px; background: #f8f9fa; border-radius: 4px;">
                <strong>Method:</strong> {results_dict['backend']} with {results_dict['estimator']} estimator
            </div>
        </div>

        <div class="card">
            <h2>📋 Detailed Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>Property</th>
                        <th>Value</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>Backend</strong></td>
                        <td class="value">{results_dict['backend']}</td>
                    </tr>
                    <tr>
                        <td><strong>Estimator</strong></td>
                        <td class="value">{results_dict['estimator']}</td>
                    </tr>
                    <tr>
                        <td><strong>Total Binding ΔG</strong></td>
                        <td class="value">{results_dict['delta_g']:.2f} kcal/mol</td>
                    </tr>
                    <tr>
                        <td><strong>Error Estimate</strong></td>
                        <td class="value">
                            {f"± {results_dict['delta_g_error']:.2f}" if results_dict['delta_g_error'] == results_dict['delta_g_error'] else "N/A"}
                        </td>
                    </tr>
                    <tr>
                        <td><strong>Bound Leg ΔG</strong></td>
                        <td class="value">{results_dict['delta_g_bound']:.2f} kcal/mol</td>
                    </tr>
                    <tr>
                        <td><strong>Unbound Leg ΔG</strong></td>
                        <td class="value">{results_dict['delta_g_unbound']:.2f} kcal/mol</td>
                    </tr>
                </tbody>
            </table>
        </div>

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
            <div class="alert">
                <strong>⚠️ Note:</strong> Convergence analysis requires trajectory data.
                For production analysis, ensure adequate sampling time and check:
                <ul style="margin: 10px 0 0 20px;">
                    <li>Time series of dH/dλ for stationarity</li>
                    <li>Block averaging for error convergence</li>
                    <li>Forward/reverse hysteresis for thermodynamic consistency</li>
                </ul>
            </div>
        </div>

        <div class="card">
            <h2>⚡ Energy Decomposition</h2>
            <div class="alert">
                <strong>📋 Note:</strong> Energy decomposition requires component-wise dH/dλ data.
                To enable energy decomposition, set <code>free-energy-components = yes</code> in MDP file.
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
                        <td class="value">{len(lambda_schedule)}</td>
                    </tr>
                    <tr>
                        <td><strong>Lambda Range</strong></td>
                        <td class="value">{lambda_schedule[0]:.3f} → {lambda_schedule[-1]:.3f}</td>
                    </tr>
                    <tr>
                        <td><strong>Estimator</strong></td>
                        <td class="value">{results_dict['estimator']}</td>
                    </tr>
                </tbody>
            </table>
        </div>

        <div class="footer">
            <p>Generated by PRISM FEP Analyzer</p>
            <p>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
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

        bound = lambda_profiles.get("bound") or {}
        unbound = lambda_profiles.get("unbound") or {}

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

        dg_bound_trace = _dg_trace(bound, "Bound Leg", "rgb(76,175,80)")
        dg_unbound_trace = _dg_trace(unbound, "Unbound Leg", "rgb(33,150,243)")

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

        if dg_bound_trace["x"]:
            bl = _json.dumps(bound_layout)
            script_lines.append(
                f"Plotly.newPlot('dg-bound-plot', {_json.dumps([dg_bound_trace])}, {bl}, {{responsive:true}});"
            )
        if dg_unbound_trace["x"]:
            ul = _json.dumps(unbound_layout)
            script_lines.append(
                f"Plotly.newPlot('dg-unbound-plot', {_json.dumps([dg_unbound_trace])}, {ul}, {{responsive:true}});"
            )

        if is_ti:
            dhdl_bound_trace = _dhdl_trace(bound, "Bound", "rgb(76,175,80)")
            dhdl_unbound_trace = _dhdl_trace(unbound, "Unbound", "rgb(33,150,243)")
            if dhdl_bound_trace["x"]:
                bhl = _json.dumps(bound_dhdl_layout)
                script_lines.append(
                    f"Plotly.newPlot('dhdl-bound-plot', {_json.dumps([dhdl_bound_trace])}, {bhl}, {{responsive:true}});"
                )
            if dhdl_unbound_trace["x"]:
                uhl = _json.dumps(unbound_dhdl_layout)
                script_lines.append(
                    f"Plotly.newPlot('dhdl-unbound-plot', {_json.dumps([dhdl_unbound_trace])}, {uhl}, {{responsive:true}});"
                )
        script_lines.append("</script>")
        return dhdl_divs, "\n".join(script_lines)


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

    def generate(self, output_path: str) -> Path:
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

        # Generate backend names
        backend_names = [f"{r['backend']}-{r['estimator']}" for r in self.all_results]

        # Build comparison table
        comparison_rows = ""
        for result in self.all_results:
            error_display = (
                f"± {result['delta_g_error']:.2f}" if result["delta_g_error"] == result["delta_g_error"] else "N/A"
            )
            comparison_rows += f"""
                    <tr>
                        <td>{result['name']}</td>
                        <td class="value">{result['delta_g']:.2f}</td>
                        <td class="value">{error_display}</td>
                        <td class="value">{result['delta_g_bound']:.2f}</td>
                        <td class="value">{result['delta_g_unbound']:.2f}</td>
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
                Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
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
            <p>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>
"""
        return html
