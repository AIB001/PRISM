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

from . import plots

if TYPE_CHECKING:
    from .analyzer import FEResults


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
            {
            self._build_overlap_matrix_html()
            or '<div class="alert"><strong>ℹ️ Note:</strong> Overlap matrix is only available for MBAR estimator.<br>Switch to MBAR backend to enable overlap matrix visualization.</div>'
        }
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


class MultiEstimatorReportGenerator:
    """
    Generate HTML reports comparing multiple estimators (TI/BAR/MBAR)

    Creates comprehensive HTML reports with:
    - Comparison table showing results from all estimators
    - Warning banner if results diverge significantly (>1 kcal/mol)
    - Tab-based switching between estimators
    - Complete analysis for each estimator (Summary, Lambda Profiles, etc.)
    - Overlap matrix only shown for MBAR (others show placeholder message)

    Parameters
    ----------
    multi_results : MultiEstimatorResults
        Results object containing data from multiple estimators

    Examples
    --------
    >>> from prism.fep.analysis import FEPMultiEstimatorAnalyzer, MultiEstimatorReportGenerator
    >>> analyzer = FEPMultiEstimatorAnalyzer(
    ...     bound_dir='fep/bound',
    ...     unbound_dir='fep/unbound',
    ...     estimators=['TI', 'BAR', 'MBAR']
    ... )
    >>> results = analyzer.analyze()
    >>> generator = MultiEstimatorReportGenerator(results)
    >>> html_path = generator.generate('comparison_report.html')
    """

    def __init__(self, multi_results):
        """
        Initialize multi-estimator report generator

        Parameters
        ----------
        multi_results : MultiEstimatorResults
            Results from multiple estimators
        """
        self.multi_results = multi_results
        self.template_style = "modern"

    def generate(self, output_path: Union[str, Path]) -> Path:
        """
        Generate HTML report with estimator comparison and tab switching

        If only one estimator is present, delegates to HTMLReportGenerator
        to avoid code duplication and ensure consistent formatting.

        Parameters
        ----------
        output_path : str or Path
            Path to output HTML file

        Returns
        -------
        Path
            Path to generated HTML file
        """
        # If only one estimator, use the standard HTMLReportGenerator
        # to avoid code duplication and ensure all plots are included
        if len(self.multi_results.methods) == 1:
            # Get the single estimator's results
            estimator_name = list(self.multi_results.methods.keys())[0]
            single_results = self.multi_results.methods[estimator_name]

            # Use HTMLReportGenerator for single estimator
            # This ensures all existing plots and features are included
            single_generator = HTMLReportGenerator(single_results)
            return single_generator.generate(output_path)

        # Multiple estimators: use multi-estimator format
        html_content = self._generate_html()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            f.write(html_content)

        return output_path

    def _generate_html(self) -> str:
        """Generate complete HTML document with comparison and tabs"""
        # Load templates
        css_content = _load_template("styles.css")
        js_content = _load_template("script.js")

        # Build warning banner (if diverged)
        warning_banner = self._build_warning_banner()

        # Build comparison table
        comparison_section = self._build_comparison_section()

        # Build estimator tabs
        estimator_tabs_content = self._build_estimator_tabs_content()

        # Build complete HTML
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>FEP Analysis - Multi-Estimator Comparison</title>
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
            <h1>🧪 FEP Analysis - Multi-Estimator Comparison</h1>
            <div class="subtitle">
                Generated on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
            </div>
        </header>

{warning_banner}

        {comparison_section}

        <div class="card">
            <h2>🔬 Detailed Analysis by Estimator</h2>

            <div class="tab-buttons">
                {self._build_tab_buttons()}
            </div>

{estimator_tabs_content}
        </div>

        <div class="footer">
            <p>Generated by PRISM FEP Analyzer (Multi-Estimator Mode)</p>
            <p>{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        </div>
    </div>
</body>
</html>"""
        return html

    def _build_tab_buttons(self) -> str:
        """Build tab button HTML for estimator switching"""
        buttons = []
        for estimator in self.multi_results.methods.keys():
            # Use format() to avoid f-string backslash issues
            buttons.append(
                '<button class="tab-btn" onclick="showEstimatorTab(\'{est}\')">{est}</button>'.format(est=estimator)
            )
        return "".join(buttons)

    def _build_warning_banner(self) -> str:
        """Build warning banner if results diverged"""
        if not self.multi_results.comparison.get("diverged", False):
            return ""

        return f"""
        <div class="card" style="background: #fff3cd; border-left: 4px solid #dc3545; margin-bottom: 20px;">
            <h3 style="color: #856404; margin-top: 0;">⚠️ Warning: Estimator Results Diverge</h3>
            <p style="margin-bottom: 10px;">
                Results from different estimators diverge significantly:
                <strong>range = {self.multi_results.comparison["delta_g_range"]:.2f} kcal/mol</strong>
            </p>
            <p style="margin-bottom: 0;">
                This may indicate:
            </p>
            <ul style="margin-top: 0; margin-bottom: 0;">
                <li>Insufficient sampling (check overlap matrix)</li>
                <li>Poor convergence (check time convergence)</li>
                <li>Inappropriate lambda spacing</li>
                <li>System far from equilibrium</li>
            </ul>
            <p style="margin-top: 10px; margin-bottom: 0;">
                <strong>💡 Recommendation:</strong> Use the estimator with the best sampling quality
                (usually MBAR with good overlap matrix) and consider extending simulation time.
            </p>
        </div>
        """

    def _build_comparison_section(self) -> str:
        """Build comparison section with table"""
        rows = ""
        for name, results in self.multi_results.methods.items():
            rows += f"""
                    <tr>
                        <td><strong>{name}</strong></td>
                        <td class="value">{results.delta_g:.2f}</td>
                        <td class="value">± {results.delta_g_error:.2f}</td>
                        <td class="value">{results.delta_g_bound:.2f}</td>
                        <td class="value">{results.delta_g_unbound:.2f}</td>
                    </tr>
            """

        # Add statistics row
        comp = self.multi_results.comparison
        stats_row = f"""
                    <tr style="border-top: 2px solid #666; font-weight: bold; background: #e8f5e9;">
                        <td>Statistics</td>
                        <td class="value">{comp["delta_g_mean"]:.2f}</td>
                        <td class="value">std = {comp["delta_g_std"]:.2f}</td>
                        <td class="value">range = {comp["delta_g_range"]:.2f}</td>
                        <td class="value">-</td>
                    </tr>
        """

        return f"""
        <div class="card">
            <h2>📊 Results Comparison</h2>
            <table>
                <thead>
                    <tr>
                        <th>Estimator</th>
                        <th>ΔΔG (kcal/mol)</th>
                        <th>Error (kcal/mol)</th>
                        <th>Bound ΔG (kcal/mol)</th>
                        <th>Unbound ΔG (kcal/mol)</th>
                    </tr>
                </thead>
                <tbody>
{rows}
{stats_row}
                </tbody>
            </table>
            <p style="margin-top: 15px; font-size: 0.9em; color: #666; background: #f8f9fa; padding: 10px; border-radius: 4px;">
                <strong>💡 Tip:</strong> Click estimator tabs below to see detailed analysis for each method.
                {
            "<br><strong>⚠️ Agreement:</strong> "
            + comp["agreement"]
            + " - Results "
            + (
                "agree well"
                if comp["agreement"] == "good"
                else "show moderate divergence"
                if comp["agreement"] == "moderate"
                else "diverge significantly"
            )
            + "."
        }
            </p>
        </div>
        """

    def _build_estimator_tabs_content(self) -> str:
        """Build tab content for each estimator"""
        tabs_content = ""

        # Get first estimator name as default
        estimator_names = list(self.multi_results.methods.keys())
        default_estimator = estimator_names[0] if estimator_names else "MBAR"

        for i, (estimator_name, results) in enumerate(self.multi_results.methods.items()):
            is_default = estimator_name == default_estimator
            tabs_content += self._build_single_estimator_tab(estimator_name, results, is_default)

        return tabs_content

    def _build_single_estimator_tab(self, estimator_name: str, results: "FEResults", is_default: bool) -> str:
        """Build content for one estimator tab (reusing HTMLReportGenerator methods)"""
        # Use HTMLReportGenerator's methods to build plots
        # This ensures all existing plots are included without code duplication
        temp_generator = HTMLReportGenerator(results)

        # Build lambda profile plots
        lambda_profiles = getattr(results, "lambda_profiles", None)
        lambda_divs, lambda_scripts = temp_generator._build_lambda_plots_html(
            lambda_profiles, estimator_name, plot_suffix=True
        )

        # Build overlap matrix section (real for MBAR, placeholder for others)
        if estimator_name == "MBAR":
            overlap_html = self._build_mbar_overlap_matrix_section(results)
        else:
            overlap_html = f"""
                <div class="card">
                    <h2>📊 Overlap Matrix ({estimator_name})</h2>
                    <div class="alert">
                        <strong>ℹ️ Info:</strong> Overlap matrix visualization is only available for MBAR estimator.
                        <br>Switch to the <strong>MBAR</strong> tab to view overlap matrix and assess sampling quality.
                        <br><br>
                        <em>Why MBAR only?</em><br>
                        MBAR uses a matrix of overlap weights between all lambda states, while TI and BAR
                        use local approximations. This makes MBAR more robust but also more computationally expensive.
                    </div>
                </div>
            """

        # Build repeat statistics table (if multiple repeats)
        repeats_html = self._build_repeats_table_for_estimator(results)

        # Build outlier plot (if multiple repeats)
        outlier_html = self._build_repeats_outlier_plot_for_estimator(results)

        return f"""
            <div id="estimator-{estimator_name}" class="tab-content" style="display: {"block" if is_default else "none"};">
                <div class="card">
                    <h2>📈 Lambda Profiles</h2>
                    <div style="display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-top:16px;">
                        <div class="plot-container" id="dg-bound-plot-{estimator_name}" style="min-height:300px;"></div>
                        <div class="plot-container" id="dg-unbound-plot-{estimator_name}" style="min-height:300px;"></div>
                    </div>
{lambda_divs}
                </div>
{lambda_scripts}
{outlier_html}
{overlap_html}
{repeats_html}
            </div>
        """

    def _build_mbar_overlap_matrix_section(self, results: "FEResults") -> str:
        """Build overlap matrix section for MBAR with actual Plotly heatmap"""
        import json as _json
        import numpy as np

        # Extract overlap matrix from lambda profiles
        lambda_profiles = getattr(results, "lambda_profiles", None)
        if not lambda_profiles:
            return self._build_overlap_matrix_placeholder("MBAR")

        # Get bound leg overlap matrix (use first repeat if multiple)
        bound_data = lambda_profiles.get("bound", {})
        if isinstance(bound_data, list):
            if not bound_data or not bound_data[0]:
                return self._build_overlap_matrix_placeholder("MBAR")
            overlap_matrix = bound_data[0].get("overlap_matrix")
        else:
            overlap_matrix = bound_data.get("overlap_matrix")

        if not overlap_matrix:
            return self._build_overlap_matrix_placeholder("MBAR")

        # Convert to numpy array
        overlap = np.array(overlap_matrix)
        n_states = overlap.shape[0]

        # Create Plotly heatmap
        heatmap_data = [
            {
                "z": overlap.tolist(),
                "x": [f"λ{i}" for i in range(n_states)],
                "y": [f"λ{i}" for i in range(n_states)],
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
                    "title": "Overlap",
                    "titleside": "right",
                    "thickness": 15,
                },
                "type": "heatmap",
                "reversescale": False,
            }
        ]

        layout = {
            "title": {
                "text": "MBAR Overlap Matrix - Phase Space Overlap Between Lambda States",
                "x": 0.5,
                "xanchor": "center",
                "font": {"size": 16},
            },
            "xaxis": {"title": "Lambda State", "side": "bottom"},
            "yaxis": {"title": "Lambda State", "side": "left"},
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(0,0,0,0)",
            "margin": {"l": 80, "r": 80, "t": 60, "b": 60},
            "width": 600,
            "height": 600,
        }

        plot_id = "mbar-overlap-matrix-plot"

        return f"""
                <div class="card">
                    <h2>📊 Overlap Matrix (MBAR)</h2>
                    <div id="{plot_id}" style="min-height:600px;"></div>
                    <p style="margin-top:15px;font-size:0.9em;color:#666;">
                        <strong>📊 Interpretation:</strong> Overlap matrix shows the phase space overlap between adjacent lambda windows.
                        <br>Values closer to 1.0 (red) indicate good overlap, while values near 0.0 (blue) suggest poor sampling.
                        <br>Ideally, all off-diagonal elements should be > 0.1 for reliable MBAR estimates.
                    </p>
                </div>
                <script>
                    Plotly.newPlot('{plot_id}', {_json.dumps(heatmap_data)}, {_json.dumps(layout)}, {{responsive:true}});
                </script>
        """

    def _build_overlap_matrix_placeholder(self, estimator_name: str) -> str:
        """Build placeholder message for non-MBAR estimators"""
        return f"""
                <div class="card">
                    <h2>📊 Overlap Matrix ({estimator_name})</h2>
                    <div class="alert">
                        <strong>ℹ️ Info:</strong> Overlap matrix visualization is only available for MBAR estimator.
                        <br>Switch to the <strong>MBAR</strong> tab to view overlap matrix and assess sampling quality.
                        <br><br>
                        <em>Why MBAR only?</em><br>
                        MBAR uses a matrix of overlap weights between all lambda states, while TI and BAR
                        use local approximations. This makes MBAR more robust but also more computationally expensive.
                    </div>
                </div>
        """

    def _build_repeats_outlier_plot_for_estimator(self, results: "FEResults") -> str:
        """Build outlier detection box plot for one estimator's results"""
        import json as _json

        if not hasattr(results, "n_repeats") or results.n_repeats <= 1:
            return ""

        if not hasattr(results, "repeat_results") or not results.repeat_results:
            return ""

        # Extract data for box plot
        bound_values = [r["bound"] for r in results.repeat_results]
        unbound_values = [r["unbound"] for r in results.repeat_results]

        # Build box plot traces
        box_traces = [
            {
                "x": ["Bound"] * len(bound_values),
                "y": bound_values,
                "name": "Bound",
                "type": "box",
                "marker": {"color": "rgba(76,175,80,0.7)"},
                "boxpoints": "outliers",
                "jitter": 0.3,
                "pointpos": 0,
            },
            {
                "x": ["Unbound"] * len(unbound_values),
                "y": unbound_values,
                "name": "Unbound",
                "type": "box",
                "marker": {"color": "rgba(244,67,54,0.7)"},
                "boxpoints": "outliers",
                "jitter": 0.3,
                "pointpos": 0,
            },
        ]

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

        plot_id = f"repeats-outlier-plot-{results.metadata.get('estimator', 'unknown')}"

        return f"""
                <div class="card">
                    <h3>📊 Multiple Repeats Statistics</h3>
                    <div id="{plot_id}" style="min-height:350px;"></div>
                    <p style="margin-top:15px;font-size:0.9em;color:#666;">
                        <strong>📊 Interpretation:</strong> Box plots show the distribution across repeats.
                        Points beyond whiskers are potential outliers.
                    </p>
                </div>
                <script>
                    Plotly.newPlot('{plot_id}', {_json.dumps(box_traces)}, {_json.dumps(layout)}, {{responsive:true}});
                </script>
        """

    def _build_repeats_table_for_estimator(self, results: "FEResults") -> str:
        """Build repeats table HTML for one estimator's results"""
        if not hasattr(results, "n_repeats") or results.n_repeats <= 1:
            return ""

        if not hasattr(results, "repeat_results") or not results.repeat_results:
            return ""

        # Build table rows for each repeat
        rows_html = ""
        for r in results.repeat_results:
            rows_html += f"""
                        <tr>
                            <td>{r["repeat"]}</td>
                            <td class="value">{r["bound"]:.2f}</td>
                            <td class="value">{r["unbound"]:.2f}</td>
                            <td class="value">{r["ddG"]:.2f}</td>
                        </tr>
            """

        # Add statistics rows (ave and stderr)
        stats = results.metadata.get("repeat_statistics", {})
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

        return f"""
                <div class="card">
                    <h3>📊 Multiple Repeats Statistics ({results.metadata.get("estimator", "")})</h3>
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
                        <strong>🔢 n_repeats:</strong> {results.n_repeats}
                    </p>
                </div>
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
