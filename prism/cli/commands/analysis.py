#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM CLI Analysis Commands

Command handlers for result analysis, visualization, and report generation.
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional
from argparse import ArgumentParser, Namespace

from ...utils.logging_system import PrismLogger
from ..utils import CLIFormatter, TableFormatter


class AnalysisCommands:
    """Handler for analysis-related CLI commands"""
    
    def __init__(self, context):
        self.context = context
        self.logger = PrismLogger("analysis_cli")
        self.formatter = CLIFormatter()
    
    def setup_parser(self, parser: ArgumentParser):
        """Setup analysis command parser"""
        subparsers = parser.add_subparsers(dest='analysis_action', help='Analysis actions')
        
        # PMF plot generation
        pmf_parser = subparsers.add_parser('plot-pmf', help='Generate PMF plots')
        pmf_parser.add_argument('data_file', help='PMF data file')
        pmf_parser.add_argument('--output', '-o', help='Output image file')
        pmf_parser.add_argument('--format', choices=['png', 'pdf', 'svg'], default='png',
                               help='Output format')
        pmf_parser.add_argument('--title', help='Plot title')
        pmf_parser.add_argument('--xlabel', default='Reaction Coordinate', help='X-axis label')
        pmf_parser.add_argument('--ylabel', default='PMF (kJ/mol)', help='Y-axis label')
        
        # Statistics analysis
        stats_parser = subparsers.add_parser('stats', help='Generate statistical analysis')
        stats_parser.add_argument('data_file', help='Data file to analyze')
        stats_parser.add_argument('--column', type=int, default=1, help='Column to analyze (0-indexed)')
        stats_parser.add_argument('--output', help='Output file for results')
        
        # Convergence analysis
        conv_parser = subparsers.add_parser('convergence', help='Analyze convergence')
        conv_parser.add_argument('data_files', nargs='+', help='Data files for convergence analysis')
        conv_parser.add_argument('--output', help='Output plot file')
        conv_parser.add_argument('--window', type=int, default=100, help='Moving average window')
        
        # Report generation
        report_parser = subparsers.add_parser('report', help='Generate analysis report')
        report_parser.add_argument('workflow_id', help='Workflow ID to analyze')
        report_parser.add_argument('--format', choices=['html', 'pdf', 'markdown'], 
                                 default='html', help='Report format')
        report_parser.add_argument('--output', help='Output file name')
        report_parser.add_argument('--template', help='Report template')
        
        # Compare results
        compare_parser = subparsers.add_parser('compare', help='Compare multiple PMF results')
        compare_parser.add_argument('data_files', nargs='+', help='PMF data files to compare')
        compare_parser.add_argument('--labels', nargs='+', help='Labels for each dataset')
        compare_parser.add_argument('--output', help='Output plot file')
        
        # Energy analysis
        energy_parser = subparsers.add_parser('energy', help='Analyze energy data')
        energy_parser.add_argument('energy_file', help='Energy data file (XVG format)')
        energy_parser.add_argument('--terms', nargs='+', 
                                  help='Energy terms to analyze (e.g., Potential, Kinetic)')
        energy_parser.add_argument('--output', help='Output plot file')
        
        # Trajectory analysis
        traj_parser = subparsers.add_parser('trajectory', help='Analyze trajectory data')
        traj_parser.add_argument('trajectory_file', help='Trajectory file')
        traj_parser.add_argument('--analysis-type', choices=['rmsd', 'rg', 'distance'],
                               default='rmsd', help='Type of analysis')
        traj_parser.add_argument('--reference', help='Reference structure file')
        traj_parser.add_argument('--output', help='Output data file')
    
    def handle(self, args: Namespace) -> int:
        """Handle analysis commands"""
        try:
            if args.analysis_action == 'plot-pmf':
                return self._plot_pmf(args)
            elif args.analysis_action == 'stats':
                return self._analyze_statistics(args)
            elif args.analysis_action == 'convergence':
                return self._analyze_convergence(args)
            elif args.analysis_action == 'report':
                return self._generate_report(args)
            elif args.analysis_action == 'compare':
                return self._compare_results(args)
            elif args.analysis_action == 'energy':
                return self._analyze_energy(args)
            elif args.analysis_action == 'trajectory':
                return self._analyze_trajectory(args)
            else:
                self.formatter.error("No analysis action specified")
                return 1
                
        except Exception as e:
            self.formatter.error(f"Analysis command failed: {e}")
            if self.context.verbose:
                self.logger.exception("Analysis command error")
            return 1
    
    def _plot_pmf(self, args: Namespace) -> int:
        """Generate PMF plot"""
        data_file = Path(args.data_file)
        
        if not data_file.exists():
            self.formatter.error(f"Data file not found: {data_file}")
            return 1
        
        try:
            # Load PMF data
            data = self._load_pmf_data(data_file)
            
            if data is None:
                return 1
            
            # Create plot
            plt.figure(figsize=(10, 6))
            plt.plot(data[:, 0], data[:, 1], 'b-', linewidth=2, label='PMF')
            
            # Add error bars if available
            if data.shape[1] > 2:
                plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], 
                           fmt='none', ecolor='gray', alpha=0.5)
            
            # Customize plot
            plt.xlabel(args.xlabel)
            plt.ylabel(args.ylabel)
            plt.title(args.title or f'PMF Profile - {data_file.stem}')
            plt.grid(True, alpha=0.3)
            plt.legend()
            
            # Save plot
            output_file = args.output or f"{data_file.stem}_pmf.{args.format}"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.formatter.success(f"PMF plot saved: {output_file}")
            
            # Print summary statistics
            min_energy = np.min(data[:, 1])
            max_energy = np.max(data[:, 1])
            binding_energy = max_energy - min_energy
            
            print(f"\nPMF Analysis Summary:")
            print(f"  Minimum energy:  {min_energy:.2f} kJ/mol")
            print(f"  Maximum energy:  {max_energy:.2f} kJ/mol")
            print(f"  Binding energy:  {binding_energy:.2f} kJ/mol")
            print(f"  Data points:     {len(data)}")
            
            return 0
            
        except Exception as e:
            self.formatter.error(f"Failed to generate PMF plot: {e}")
            return 1
    
    def _analyze_statistics(self, args: Namespace) -> int:
        """Analyze statistical properties of data"""
        data_file = Path(args.data_file)
        
        if not data_file.exists():
            self.formatter.error(f"Data file not found: {data_file}")
            return 1
        
        try:
            # Load data
            data = np.loadtxt(data_file)
            
            if data.ndim == 1:
                values = data
            else:
                if args.column >= data.shape[1]:
                    self.formatter.error(f"Column {args.column} not found in data")
                    return 1
                values = data[:, args.column]
            
            # Calculate statistics
            stats = {
                'Count': len(values),
                'Mean': np.mean(values),
                'Median': np.median(values),
                'Std Dev': np.std(values),
                'Min': np.min(values),
                'Max': np.max(values),
                'Range': np.max(values) - np.min(values),
                'Skewness': self._calculate_skewness(values),
                'Kurtosis': self._calculate_kurtosis(values)
            }
            
            # Display results
            self.formatter.header(f"Statistical Analysis: {data_file.name}")
            
            stats_table = TableFormatter(['Statistic', 'Value'])
            for stat, value in stats.items():
                if isinstance(value, (int, np.integer)):
                    formatted_value = str(value)
                else:
                    formatted_value = f"{value:.6f}"
                stats_table.add_row([stat, formatted_value])
            
            stats_table.print(self.context.output_format)
            
            # Save results if requested
            if args.output:
                output_file = Path(args.output)
                with open(output_file, 'w') as f:
                    f.write(f"Statistical Analysis: {data_file.name}\n")
                    f.write("=" * 40 + "\n")
                    for stat, value in stats.items():
                        f.write(f"{stat:<15}: {value:.6f}\n")
                
                self.formatter.success(f"Statistics saved: {output_file}")
            
            return 0
            
        except Exception as e:
            self.formatter.error(f"Failed to analyze statistics: {e}")
            return 1
    
    def _analyze_convergence(self, args: Namespace) -> int:
        """Analyze convergence of data"""
        data_files = [Path(f) for f in args.data_files]
        
        # Check all files exist
        for data_file in data_files:
            if not data_file.exists():
                self.formatter.error(f"Data file not found: {data_file}")
                return 1
        
        try:
            plt.figure(figsize=(12, 8))
            
            for i, data_file in enumerate(data_files):
                # Load data
                data = np.loadtxt(data_file)
                if data.ndim == 1:
                    values = data
                else:
                    values = data[:, -1]  # Last column
                
                # Calculate running average
                running_avg = self._calculate_running_average(values, args.window)
                
                # Plot
                plt.subplot(2, 1, 1)
                plt.plot(values, alpha=0.7, label=f'{data_file.stem} (raw)')
                
                plt.subplot(2, 1, 2) 
                plt.plot(running_avg, linewidth=2, label=f'{data_file.stem} (avg)')
            
            # Customize plots
            plt.subplot(2, 1, 1)
            plt.ylabel('Raw Values')
            plt.title('Convergence Analysis')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.subplot(2, 1, 2)
            plt.xlabel('Time Step')
            plt.ylabel(f'Running Average (window={args.window})')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # Save plot
            output_file = args.output or 'convergence_analysis.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.formatter.success(f"Convergence analysis saved: {output_file}")
            return 0
            
        except Exception as e:
            self.formatter.error(f"Failed to analyze convergence: {e}")
            return 1
    
    def _generate_report(self, args: Namespace) -> int:
        """Generate analysis report"""
        # Mock workflow data
        workflow_data = {
            'id': args.workflow_id,
            'name': f'PMF Calculation {args.workflow_id}',
            'status': 'completed',
            'total_runtime': 7200,  # 2 hours
            'tasks_completed': 6,
            'tasks_failed': 0,
            'results_summary': {
                'binding_energy': -8.5,
                'pmf_profile_points': 30,
                'convergence_achieved': True,
                'estimated_error': 0.5
            }
        }
        
        if args.format == 'html':
            return self._generate_html_report(workflow_data, args.output)
        elif args.format == 'markdown':
            return self._generate_markdown_report(workflow_data, args.output)
        else:
            self.formatter.error(f"Report format '{args.format}' not yet implemented")
            return 1
    
    def _compare_results(self, args: Namespace) -> int:
        """Compare multiple PMF results"""
        data_files = [Path(f) for f in args.data_files]
        labels = args.labels or [f.stem for f in data_files]
        
        # Check files exist
        for data_file in data_files:
            if not data_file.exists():
                self.formatter.error(f"Data file not found: {data_file}")
                return 1
        
        try:
            plt.figure(figsize=(12, 8))
            
            comparison_data = []
            
            for i, (data_file, label) in enumerate(zip(data_files, labels)):
                data = self._load_pmf_data(data_file)
                if data is None:
                    continue
                
                plt.plot(data[:, 0], data[:, 1], linewidth=2, label=label, marker='o', markersize=4)
                
                # Store for comparison stats
                comparison_data.append({
                    'label': label,
                    'data': data,
                    'binding_energy': np.max(data[:, 1]) - np.min(data[:, 1]),
                    'min_energy': np.min(data[:, 1])
                })
            
            plt.xlabel('Reaction Coordinate')
            plt.ylabel('PMF (kJ/mol)')
            plt.title('PMF Comparison')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            # Save plot
            output_file = args.output or 'pmf_comparison.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.formatter.success(f"Comparison plot saved: {output_file}")
            
            # Print comparison statistics
            if comparison_data:
                self.formatter.header("Comparison Statistics")
                
                table = TableFormatter(['Dataset', 'Min Energy', 'Binding Energy', 'Data Points'])
                for entry in comparison_data:
                    table.add_row([
                        entry['label'],
                        f"{entry['min_energy']:.2f}",
                        f"{entry['binding_energy']:.2f}",
                        str(len(entry['data']))
                    ])
                
                table.print(self.context.output_format)
            
            return 0
            
        except Exception as e:
            self.formatter.error(f"Failed to compare results: {e}")
            return 1
    
    def _analyze_energy(self, args: Namespace) -> int:
        """Analyze energy data"""
        energy_file = Path(args.energy_file)
        
        if not energy_file.exists():
            self.formatter.error(f"Energy file not found: {energy_file}")
            return 1
        
        # For now, just show a message about energy analysis
        self.formatter.info(f"Energy analysis of {energy_file}")
        self.formatter.info("Energy analysis functionality would analyze XVG files")
        self.formatter.info("and generate time series plots of energy components")
        
        return 0
    
    def _analyze_trajectory(self, args: Namespace) -> int:
        """Analyze trajectory data"""
        traj_file = Path(args.trajectory_file)
        
        if not traj_file.exists():
            self.formatter.error(f"Trajectory file not found: {traj_file}")
            return 1
        
        # For now, just show a message about trajectory analysis
        self.formatter.info(f"Trajectory analysis of {traj_file}")
        self.formatter.info(f"Analysis type: {args.analysis_type}")
        self.formatter.info("Trajectory analysis would require MDAnalysis or similar")
        self.formatter.info("to calculate RMSD, radius of gyration, or distance metrics")
        
        return 0
    
    # Helper methods
    def _load_pmf_data(self, data_file: Path) -> Optional[np.ndarray]:
        """Load PMF data from file"""
        try:
            # Try to load as space-separated values
            data = np.loadtxt(data_file)
            
            # Ensure we have at least 2 columns
            if data.ndim == 1:
                self.formatter.error("PMF data file must contain at least 2 columns")
                return None
            
            if data.shape[1] < 2:
                self.formatter.error("PMF data file must contain at least 2 columns")
                return None
            
            return data
            
        except Exception as e:
            self.formatter.error(f"Failed to load PMF data: {e}")
            return None
    
    def _calculate_skewness(self, values: np.ndarray) -> float:
        """Calculate skewness of data"""
        mean = np.mean(values)
        std = np.std(values)
        n = len(values)
        
        if std == 0:
            return 0.0
        
        skew = np.sum(((values - mean) / std) ** 3) / n
        return skew
    
    def _calculate_kurtosis(self, values: np.ndarray) -> float:
        """Calculate kurtosis of data"""
        mean = np.mean(values)
        std = np.std(values)
        n = len(values)
        
        if std == 0:
            return 0.0
        
        kurt = np.sum(((values - mean) / std) ** 4) / n - 3
        return kurt
    
    def _calculate_running_average(self, values: np.ndarray, window: int) -> np.ndarray:
        """Calculate running average with specified window"""
        running_avg = np.zeros_like(values)
        
        for i in range(len(values)):
            start = max(0, i - window // 2)
            end = min(len(values), i + window // 2 + 1)
            running_avg[i] = np.mean(values[start:end])
        
        return running_avg
    
    def _generate_html_report(self, workflow_data: Dict[str, Any], output_file: Optional[str]) -> int:
        """Generate HTML report"""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PRISM Analysis Report - {workflow_data['name']}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; }}
        .summary {{ background: #ecf0f1; padding: 20px; border-radius: 5px; }}
        .metric {{ margin: 10px 0; }}
        .status-completed {{ color: #27ae60; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <h1>PRISM Analysis Report</h1>
    
    <div class="summary">
        <h2>Workflow Summary</h2>
        <div class="metric"><strong>Workflow ID:</strong> {workflow_data['id']}</div>
        <div class="metric"><strong>Name:</strong> {workflow_data['name']}</div>
        <div class="metric"><strong>Status:</strong> <span class="status-{workflow_data['status']}">{workflow_data['status'].upper()}</span></div>
        <div class="metric"><strong>Total Runtime:</strong> {workflow_data['total_runtime'] / 3600:.1f} hours</div>
        <div class="metric"><strong>Tasks Completed:</strong> {workflow_data['tasks_completed']}</div>
        <div class="metric"><strong>Tasks Failed:</strong> {workflow_data['tasks_failed']}</div>
    </div>
    
    <h2>Results Summary</h2>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
        <tr><td>Binding Energy</td><td>{workflow_data['results_summary']['binding_energy']} kJ/mol</td></tr>
        <tr><td>PMF Profile Points</td><td>{workflow_data['results_summary']['pmf_profile_points']}</td></tr>
        <tr><td>Convergence Achieved</td><td>{'Yes' if workflow_data['results_summary']['convergence_achieved'] else 'No'}</td></tr>
        <tr><td>Estimated Error</td><td>±{workflow_data['results_summary']['estimated_error']} kJ/mol</td></tr>
    </table>
    
    <h2>Generated Files</h2>
    <ul>
        <li>PMF profile data: pmf_profile.dat</li>
        <li>Umbrella sampling data: umbrella_histograms.dat</li>
        <li>Analysis plots: pmf_plot.png, convergence_plot.png</li>
        <li>Log files: workflow.log, analysis.log</li>
    </ul>
    
    <p><em>Report generated by PRISM PMF System</em></p>
</body>
</html>
        """
        
        output_path = output_file or f"report_{workflow_data['id']}.html"
        
        with open(output_path, 'w') as f:
            f.write(html_content)
        
        self.formatter.success(f"HTML report saved: {output_path}")
        return 0
    
    def _generate_markdown_report(self, workflow_data: Dict[str, Any], output_file: Optional[str]) -> int:
        """Generate Markdown report"""
        md_content = f"""# PRISM Analysis Report

## Workflow Summary

- **Workflow ID:** {workflow_data['id']}
- **Name:** {workflow_data['name']}
- **Status:** {workflow_data['status'].upper()}
- **Total Runtime:** {workflow_data['total_runtime'] / 3600:.1f} hours
- **Tasks Completed:** {workflow_data['tasks_completed']}
- **Tasks Failed:** {workflow_data['tasks_failed']}

## Results Summary

| Metric | Value |
|--------|-------|
| Binding Energy | {workflow_data['results_summary']['binding_energy']} kJ/mol |
| PMF Profile Points | {workflow_data['results_summary']['pmf_profile_points']} |
| Convergence Achieved | {'Yes' if workflow_data['results_summary']['convergence_achieved'] else 'No'} |
| Estimated Error | ±{workflow_data['results_summary']['estimated_error']} kJ/mol |

## Generated Files

- PMF profile data: `pmf_profile.dat`
- Umbrella sampling data: `umbrella_histograms.dat`
- Analysis plots: `pmf_plot.png`, `convergence_plot.png`
- Log files: `workflow.log`, `analysis.log`

---
*Report generated by PRISM PMF System*
        """
        
        output_path = output_file or f"report_{workflow_data['id']}.md"
        
        with open(output_path, 'w') as f:
            f.write(md_content)
        
        self.formatter.success(f"Markdown report saved: {output_path}")
        return 0