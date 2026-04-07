#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML Report Generator for FEP Analysis

Generates comprehensive, publication-ready HTML reports for FEP calculations.
Includes interactive plots, tables, and convergence diagnostics.

Modern UI design with static templates.
"""

from pathlib import Path
from typing import Optional, Dict, Any, Tuple, Union, TYPE_CHECKING
from datetime import datetime

from .. import plots

if TYPE_CHECKING:
    pass


def _load_template(filename: str) -> str:
    """Load template file from templates directory."""
    template_dir = Path(__file__).parent.parent / "templates"
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
            <div class="equal-grid two-col">
                <div style="min-width:0;width:100%;">
                    <h3 style="color: #333; margin-bottom: 15px; font-size: 1.2em;">
                        📊 Multiple Repeats Statistics
                    </h3>
{repeats_html}
                </div>
                <div style="min-width:0;width:100%;">
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
            <div class="equal-grid two-col">
                <div class="plot-container compact-plot" id="dg-bound-plot" style="height:380px;min-height:380px;width:100%;overflow:hidden;"></div>
                <div class="plot-container compact-plot" id="dg-unbound-plot" style="height:380px;min-height:380px;width:100%;overflow:hidden;"></div>
            </div>
{dhdl_divs}
        </div>
{lambda_plot_script}

        <div class="card">
            <h2>📊 Overlap Matrix (MBAR)</h2>
            {
                self._build_overlap_matrix_html()
                or '<div class="alert"><strong>ℹ️ Note:</strong> Overlap matrix is only available for MBAR estimator.<br>Switch to MBAR backend to enable overlap matrix visualization.</div>'
            }
        </div>

        <div class="card">
            <h2>📊 Diagnostics</h2>
            <div class="equal-grid two-col">
                <div style="min-width:0;width:100%;">
                    <h3 style="color: #333; margin-bottom: 15px; font-size: 1.2em;">
                        📈 Time Convergence
                    </h3>
{time_conv_html}
                </div>
                <div style="min-width:0;width:100%;">
                    <h3 style="color: #333; margin-bottom: 15px; font-size: 1.2em;">
                        🔄 Bootstrap Analysis
                    </h3>
{bootstrap_html}
                </div>
            </div>
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
        plot_suffix: str = "",
    ) -> Tuple[str, str]:
        """
        Build HTML div placeholders and Plotly JS for lambda-dependent plots.

        Delegates to shared plotting function in plots module.

        Parameters
        ----------
        lambda_profiles : dict or None
            As produced by FEPAnalyzer._build_lambda_profiles.
        estimator_name : str
            Estimator name (TI, BAR, MBAR).
        plot_suffix : str, optional
            Suffix to add to plot IDs.

        Returns
        -------
        tuple[str, str]
            (dhdl_divs_html, script_tag_html)
        """
        return plots.build_lambda_plots_html(lambda_profiles, estimator_name, plot_suffix)

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

        stats = self.results.metadata.get("repeat_statistics", {})
        return plots.build_repeats_table_html(
            self.results.repeat_results,
            self.results.n_repeats,
            stats,
        )

    def _build_repeats_outlier_plot_html(self) -> str:
        """
        Build HTML for outlier detection plot across multiple repeats.

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

        return plots.build_repeats_outlier_plot_html(self.results.repeat_results)

    def _build_time_convergence_html(self, time_conv_data: Optional[Dict]) -> str:
        """
        Build HTML for time convergence analysis section.

        Parameters
        ----------
        time_conv_data : dict or None
            As returned by FEPAnalyzer._compute_time_convergence.

        Returns
        -------
        str
            HTML content (empty string if no data).
        """
        return plots.build_time_convergence_html(time_conv_data)

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
        return plots.build_bootstrap_html(bootstrap_data)

    def _build_overlap_matrix_html(self) -> str:
        """
        Build HTML for MBAR overlap matrix visualization using Plotly.

        Returns
        -------
        str
            HTML content with overlap matrix plots (empty string if no overlap data).
        """
        lambda_profiles = getattr(self.results, "lambda_profiles", None)
        if not lambda_profiles:
            return ""

        estimator_name = getattr(self.results, "estimator", "MBAR")
        return plots.build_overlap_matrix_html(lambda_profiles, estimator_name)
