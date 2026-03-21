#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML Report Generator for FEP Analysis

Generates comprehensive, publication-ready HTML reports for FEP calculations.
Includes interactive plots, tables, and convergence diagnostics.
"""

from pathlib import Path
from typing import Optional
from datetime import datetime


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
    template_path : str, optional
        Path to custom HTML template (default: uses built-in template)

    Examples
    --------
    >>> generator = HTMLReportGenerator(results)
    >>> html_path = generator.generate('report.html')
    """

    def __init__(self, results, raw_data=None, template_path: Optional[str] = None):
        """
        Initialize report generator

        Parameters
        ----------
        results : FEResults
            Analysis results
        raw_data : dict, optional
            Raw data from analysis (contains fitted models for plotting)
        template_path : str, optional
            Custom template path
        """
        self.results = results
        self.raw_data = raw_data or {}
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

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PRISM-FEP Analysis Report</title>
    {self._get_css()}
    {self._get_javascript()}
</head>
<body>
    <div class="container">
        {self._generate_header()}
        {self._generate_summary()}
        {self._generate_results_table()}
        {self._generate_plots_section()}
        {self._generate_convergence_section()}
        {self._generate_decomposition_section()}
        {self._generate_methods_section()}
        {self._generate_footer()}
    </div>
</body>
</html>
"""
        return html

    def _get_css(self) -> str:
        """Get CSS styles"""
        return """
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
        }}

        .container {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}

        header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }}

        header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        header .subtitle {{
            font-size: 1.2em;
            opacity: 0.9;
        }}

        .card {{
            background: white;
            border-radius: 10px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}

        .card h2 {{
            color: #667eea;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #667eea;
        }}

        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}

        .summary-item {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}

        .summary-item .label {{
            font-size: 0.9em;
            color: #666;
            margin-bottom: 5px;
        }}

        .summary-item .value {{
            font-size: 1.8em;
            font-weight: bold;
            color: #333;
        }}

        .summary-item .error {{
            font-size: 0.9em;
            color: #999;
            margin-top: 5px;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }}

        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}

        th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #333;
        }}

        tr:hover {{
            background: #f8f9fa;
        }}

        .plot-container {{
            margin: 20px 0;
        }}

        footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 0.9em;
        }}

        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: 600;
        }}

        .badge-success {{
            background: #d4edda;
            color: #155724;
        }}

        .badge-warning {{
            background: #fff3cd;
            color: #856404;
        }}

        .badge-info {{
            background: #d1ecf1;
            color: #0c5460;
        }}
    </style>
"""

    def _get_javascript(self) -> str:
        """Get JavaScript for interactive plots"""
        return """
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script>
        // Plotly will be used for interactive plots
        function createPlot(divId, data, layout) {
            Plotly.newPlot(divId, data, layout, {responsive: true});
        }
    </script>
"""

    def _generate_header(self) -> str:
        """Generate report header"""
        return f"""
        <header>
            <h1>🧪 FEP Analysis Report</h1>
            <div class="subtitle">
                Free Energy Perturbation Calculation Results<br>
                Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            </div>
        </header>
"""

    def _generate_summary(self) -> str:
        """Generate summary section"""
        r = self.results
        return f"""
        <div class="card">
            <h2>Summary</h2>
            <div class="summary-grid">
                <div class="summary-item">
                    <div class="label">Binding Free Energy (ΔG)</div>
                    <div class="value">{r.delta_g:.2f}</div>
                    <div class="error">± {r.delta_g_error:.2f} kcal/mol</div>
                </div>
                <div class="summary-item">
                    <div class="label">Bound Leg (ΔG<sub>bound</sub>)</div>
                    <div class="value">{r.delta_g_bound:.2f}</div>
                    <div class="error">kcal/mol</div>
                </div>
                <div class="summary-item">
                    <div class="label">Unbound Leg (ΔG<sub>unbound</sub>)</div>
                    <div class="value">{r.delta_g_unbound:.2f}</div>
                    <div class="error">kcal/mol</div>
                </div>
                <div class="summary-item">
                    <div class="label">Estimator</div>
                    <div class="value"><span class="badge badge-info">{r.metadata.get('estimator', 'N/A')}</span></div>
                    <div class="error">T = {r.metadata.get('temperature', 310)} K</div>
                </div>
            </div>
        </div>
"""

    def _generate_results_table(self) -> str:
        """Generate detailed results table"""
        r = self.results
        return f"""
        <div class="card">
            <h2>Detailed Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>Quantity</th>
                        <th>Value</th>
                        <th>Unit</th>
                        <th>Description</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>ΔG<sub>bind</sub></strong></td>
                        <td>{r.delta_g:.2f} ± {r.delta_g_error:.2f}</td>
                        <td>kcal/mol</td>
                        <td>Binding free energy (unbound - bound)</td>
                    </tr>
                    <tr>
                        <td>ΔG<sub>bound</sub></td>
                        <td>{r.delta_g_bound:.2f}</td>
                        <td>kcal/mol</td>
                        <td>Free energy change for bound leg</td>
                    </tr>
                    <tr>
                        <td>ΔG<sub>unbound</sub></td>
                        <td>{r.delta_g_unbound:.2f}</td>
                        <td>kcal/mol</td>
                        <td>Free energy change for unbound leg</td>
                    </tr>
                    <tr>
                        <td>Temperature</td>
                        <td>{r.metadata.get('temperature', 310):.1f}</td>
                        <td>K</td>
                        <td>Simulation temperature</td>
                    </tr>
                    <tr>
                        <td>Lambda Windows</td>
                        <td>{r.metadata.get('n_lambda_windows', 'N/A')}</td>
                        <td>-</td>
                        <td>Number of λ windows</td>
                    </tr>
                </tbody>
            </table>
        </div>
"""

    def _generate_plots_section(self) -> str:
        """Generate plots section with Plotly visualizations"""
        # Check if we have fitted data for plotting
        if not self.raw_data:
            return self._generate_placeholder_plots()

        fitted_bound = self.raw_data.get("fitted_bound")
        fitted_unbound = self.raw_data.get("fitted_unbound")

        if fitted_bound is None and fitted_unbound is None:
            return self._generate_placeholder_plots()

        # Generate ΔG vs λ plot
        dg_plot_html = self._generate_dg_plot(fitted_bound, fitted_unbound)

        # Generate dH/dλ vs λ plot (for TI estimator)
        dhdl_plot_html = self._generate_dhdl_plot(fitted_bound, fitted_unbound)

        return f"""
        <div class="card">
            <h2>Free Energy Profile</h2>
            <div class="plot-container">
                <h3>ΔG vs λ (Cumulative Free Energy)</h3>
                {dg_plot_html}
            </div>
            <div class="plot-container">
                <h3>dH/dλ vs λ (Free Energy Gradient)</h3>
                {dhdl_plot_html}
            </div>
        </div>
"""

    def _generate_placeholder_plots(self) -> str:
        """Generate placeholder for plots section"""
        return """
        <div class="card">
            <h2>Free Energy Profile</h2>
            <div class="plot-container" id="dg-plot">
                <div style="text-align: center; padding: 40px; color: #999;">
                    <p>Interactive plots will be generated here</p>
                    <p><small>Requires plotting data from analysis</small></p>
                </div>
            </div>
        </div>
"""

    def _generate_dg_plot(self, fitted_bound, fitted_unbound) -> str:
        """Generate ΔG vs λ plot using Plotly"""
        import json

        plots_data = []

        # Plot bound leg
        if fitted_bound is not None:
            lambda_vals, dg_vals, dg_errors = self._extract_dg_data(fitted_bound)
            plots_data.append(
                {
                    "x": lambda_vals,
                    "y": [v / 4.184 for v in dg_vals],  # Convert kJ/mol to kcal/mol
                    "error_y": {"array": [e / 4.184 for e in dg_errors], "visible": True},
                    "mode": "lines+markers",
                    "name": "Bound Leg",
                    "line": {"color": "#667eea", "width": 2},
                    "marker": {"size": 6},
                }
            )

        # Plot unbound leg
        if fitted_unbound is not None:
            lambda_vals, dg_vals, dg_errors = self._extract_dg_data(fitted_unbound)
            plots_data.append(
                {
                    "x": lambda_vals,
                    "y": [v / 4.184 for v in dg_vals],  # Convert kJ/mol to kcal/mol
                    "error_y": {"array": [e / 4.184 for e in dg_errors], "visible": True},
                    "mode": "lines+markers",
                    "name": "Unbound Leg",
                    "line": {"color": "#764ba2", "width": 2},
                    "marker": {"size": 6},
                }
            )

        layout = {
            "xaxis": {"title": "λ (Lambda)", "tickformat": ".2f"},
            "yaxis": {"title": "ΔG (kcal/mol)", "zeroline": True},
            "hovermode": "closest",
            "showlegend": True,
            "height": 500,
        }

        return f"""
        <div id="dg-plot"></div>
        <script>
        var dgData = {json.dumps(plots_data)};
        var dgLayout = {json.dumps(layout)};
        Plotly.newPlot('dg-plot', dgData, dgLayout, {{responsive: true}});
        </script>
"""

    def _generate_dhdl_plot(self, fitted_bound, fitted_unbound) -> str:
        """Generate dH/dλ vs λ plot using Plotly (TI only)"""
        import json

        # Only for TI estimator
        if not hasattr(fitted_bound, "dhdl") if fitted_bound else True:
            return '<p style="color: #666; text-align: center;">dH/dλ plot only available for TI estimator</p>'

        plots_data = []

        # Plot bound leg dH/dλ
        if fitted_bound is not None and hasattr(fitted_bound, "dhdl"):
            lambda_vals, dhdl_vals = self._extract_dhdl_data(fitted_bound)
            # Plot vdw component (column index 2)
            if len(dhdl_vals) > 2:
                plots_data.append(
                    {
                        "x": lambda_vals,
                        "y": [v / 4.184 for v in dhdl_vals[2]],  # vdw component, kJ to kcal
                        "mode": "lines+markers",
                        "name": "Bound (vdW)",
                        "line": {"color": "#667eea", "width": 2},
                    }
                )

        # Plot unbound leg dH/dλ
        if fitted_unbound is not None and hasattr(fitted_unbound, "dhdl"):
            lambda_vals, dhdl_vals = self._extract_dhdl_data(fitted_unbound)
            if len(dhdl_vals) > 2:
                plots_data.append(
                    {
                        "x": lambda_vals,
                        "y": [v / 4.184 for v in dhdl_vals[2]],  # vdw component
                        "mode": "lines+markers",
                        "name": "Unbound (vdW)",
                        "line": {"color": "#764ba2", "width": 2},
                    }
                )

        if not plots_data:
            return '<p style="color: #666; text-align: center;">No dH/dλ data available</p>'

        layout = {
            "xaxis": {"title": "λ (Lambda)", "tickformat": ".2f"},
            "yaxis": {"title": "dH/dλ (kcal/mol)", "zeroline": True},
            "hovermode": "closest",
            "showlegend": True,
            "height": 500,
        }

        return f"""
        <div id="dhdl-plot"></div>
        <script>
        var dhdlData = {json.dumps(plots_data)};
        var dhdlLayout = {json.dumps(layout)};
        Plotly.newPlot('dhdl-plot', dhdlData, dhdlLayout, {{responsive: true}});
        </script>
"""

    def _extract_dg_data(self, fitted) -> tuple:
        """Extract lambda values, ΔG, and errors from fitted estimator"""
        # delta_f_ is a DataFrame with states as both index and columns
        # iloc[0] gives ΔG from state 0 to all states (in kT)

        # Get lambda values from index (MultiIndex with 4 components: mass, coul, vdw, bonded)
        # Use vdw-lambda as the primary lambda coordinate
        lambda_vals = []
        for idx in fitted.delta_f_.index:
            if isinstance(idx, tuple):
                # Use vdw-lambda (index 2) or coul-lambda (index 1)
                lambda_val = idx[2] if len(idx) > 2 else idx[0]
                lambda_vals.append(lambda_val)
            else:
                lambda_vals.append(idx)

        # Get ΔG values from first row (ΔG from state 0 to each state, in kT)
        dg_vals = fitted.delta_f_.iloc[0].values.tolist()

        # Get errors if available
        if hasattr(fitted, "d_delta_f_"):
            dg_errors = fitted.d_delta_f_.iloc[0].values.tolist()
        else:
            dg_errors = [0.0] * len(dg_vals)

        return lambda_vals, dg_vals, dg_errors

    def _extract_dhdl_data(self, fitted) -> tuple:
        """Extract lambda values and dH/dλ from fitted TI estimator"""
        # dhdl is a DataFrame with lambda states as index and components as columns
        # Columns: mass, coul, vdw, bonded

        lambda_vals = []
        for idx in fitted.dhdl.index:
            if isinstance(idx, tuple):
                lambda_val = idx[2] if len(idx) > 2 else idx[0]
                lambda_vals.append(lambda_val)
            else:
                lambda_vals.append(idx)

        # Get dH/dλ values for each component
        dhdl_vals = fitted.dhdl.values.T.tolist()  # Transpose to get [mass, coul, vdw, bonded]

        return lambda_vals, dhdl_vals

    def _generate_convergence_section(self) -> str:
        """Generate convergence analysis section"""
        conv = self.results.convergence
        status = conv.get("time_convergence", {}).get("status", "unknown")

        return f"""
        <div class="card">
            <h2>Convergence Analysis</h2>
            <table>
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>Total Samples</td>
                        <td>{conv.get('n_samples', 'N/A'):,}</td>
                        <td><span class="badge badge-success">OK</span></td>
                    </tr>
                    <tr>
                        <td>Time Convergence</td>
                        <td>{status}</td>
                        <td><span class="badge badge-info">Analysis</span></td>
                    </tr>
                    <tr>
                        <td>Hysteresis</td>
                        <td>{conv.get('hysteresis', 'N/A')}</td>
                        <td><span class="badge badge-info">Pending</span></td>
                    </tr>
                </tbody>
            </table>
        </div>
"""

    def _generate_decomposition_section(self) -> str:
        """Generate energy decomposition section"""
        components = self.results.delta_g_components

        if not components:
            return """
        <div class="card">
            <h2>Energy Decomposition</h2>
            <p style="color: #666; padding: 20px;">
                Energy decomposition not available. Enable component analysis
                during simulation to see contributions from electrostatic and
                van der Waals interactions.
            </p>
        </div>
"""

        rows = ""
        for comp_name, comp_value in components.items():
            rows += f"""
                    <tr>
                        <td><strong>{comp_name}</strong></td>
                        <td>{comp_value:.2f}</td>
                        <td>kcal/mol</td>
                    </tr>
"""

        return f"""
        <div class="card">
            <h2>Energy Decomposition</h2>
            <table>
                <thead>
                    <tr>
                        <th>Component</th>
                        <th>ΔG</th>
                        <th>Unit</th>
                    </tr>
                </thead>
                <tbody>
                    {rows}
                </tbody>
            </table>
        </div>
"""

    def _generate_methods_section(self) -> str:
        """Generate methods section"""
        estimator = self.results.metadata.get("estimator", "MBAR")
        method_descriptions = {
            "MBAR": "Multistate Bennett Acceptance Ratio - Uses all λ windows simultaneously for optimal statistical efficiency",
            "BAR": "Bennett Acceptance Ratio - Computes free energy differences between adjacent λ windows",
            "TI": "Thermodynamic Integration - Numerically integrates dH/dλ over λ",
        }

        return f"""
        <div class="card">
            <h2>Methods</h2>
            <h3>Free Energy Estimator</h3>
            <p><strong>{estimator}</strong></p>
            <p>{method_descriptions.get(estimator, 'Unknown estimator')}</p>

            <h3>Simulation Parameters</h3>
            <ul>
                <li>Temperature: {self.results.metadata.get('temperature', 310)} K</li>
                <li>Lambda windows: {self.results.metadata.get('n_lambda_windows', 'N/A')}</li>
                <li>Analysis performed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</li>
            </ul>
        </div>
"""

    def _generate_footer(self) -> str:
        """Generate report footer"""
        return f"""
        <footer>
            <p>Generated by <strong>PRISM-FEP</strong> (Protein Receptor Interaction Simulation Modeler)</p>
            <p>Version 0.1.0 | {datetime.now().year}</p>
        </footer>
"""
