#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF Analysis module - WHAM analysis and visualization
"""

import os
import subprocess
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class PMFAnalyzer:
    """PMF analysis and WHAM runner"""
    
    def __init__(self, pmf_system):
        self.pmf_system = pmf_system
        self.analysis_dir = pmf_system.analysis_dir
        self.config = pmf_system.config['analysis']
    
    def run_wham(self, begin_time=None, bootstrap=None, energy_unit='kCal', **kwargs):
        """
        Run WHAM analysis on umbrella sampling results
        
        Parameters:
        -----------
        begin_time : int, optional
            Start time for analysis (ps)
        bootstrap : int, optional
            Number of bootstrap iterations
        energy_unit : str
            Energy unit for output
        
        Returns:
        --------
        Dict : WHAM analysis results
        """
        self.analysis_dir.mkdir(exist_ok=True)
        
        # Get parameters
        if begin_time is None:
            begin_time = self.config.get('begin_time_ps', 2000)
        if bootstrap is None:
            bootstrap = self.config.get('bootstrap', 50)
        
        # Prepare file lists
        wham_files = self._prepare_wham_files()
        
        # Run WHAM
        return self._run_wham_command(wham_files, begin_time, bootstrap, energy_unit)
    
    def visualize(self, pmf_file=None, error_file=None, **kwargs):
        """Generate PMF visualization plots"""
        if pmf_file is None:
            pmf_file = self.analysis_dir / "pmf.xvg"
        if error_file is None:
            error_file = self.analysis_dir / "pmferror.xvg"
        
        # Generate plots
        plots = {}
        
        # Main PMF plot
        plots['pmf_curve'] = self._plot_pmf_curve(pmf_file, error_file)
        
        # Force profile
        plots['force_profile'] = self._plot_force_profile(pmf_file)
        
        # Energy landscape
        plots['energy_landscape'] = self._plot_energy_landscape(pmf_file)
        
        return plots
    
    def _prepare_wham_files(self):
        """Prepare TPR and pullf file lists for WHAM"""
        umbrella_dir = self.pmf_system.umbrella_dir
        
        tpr_files = []
        pullf_files = []
        
        # Find all completed umbrella windows
        for window_dir in sorted(umbrella_dir.glob("window_*")):
            tpr_file = window_dir / "results" / "umbrella.tpr"
            pullf_file = window_dir / "results" / "umbrella_pullf.xvg"
            
            if tpr_file.exists() and pullf_file.exists():
                tpr_files.append(tpr_file)
                pullf_files.append(pullf_file)
        
        if not tpr_files:
            raise RuntimeError("No completed umbrella windows found")
        
        # Write file lists
        tpr_list = self.analysis_dir / "tpr-files.dat"
        pullf_list = self.analysis_dir / "pullf-files.dat"
        
        with open(tpr_list, 'w') as f:
            for tpr_file in tpr_files:
                f.write(f"{tpr_file.absolute()}\n")
        
        with open(pullf_list, 'w') as f:
            for pullf_file in pullf_files:
                f.write(f"{pullf_file.absolute()}\n")
        
        return {
            'tpr_list': tpr_list,
            'pullf_list': pullf_list,
            'n_windows': len(tpr_files)
        }
    
    def _run_wham_command(self, wham_files, begin_time, bootstrap, energy_unit):
        """Execute GROMACS WHAM command"""
        wham_cmd = [
            "gmx", "wham",
            "-it", str(wham_files['tpr_list']),
            "-if", str(wham_files['pullf_list']),
            "-b", str(begin_time),
            "-o", "pmf.xvg",
            "-unit", energy_unit,
            "-bsres", "pmferror.xvg",
            "-bsprof", "bsprofile.xvg",
            "-nBootstrap", str(bootstrap),
            "-bs-method", "b-hist",
            "-v"
        ]
        
        try:
            result = subprocess.run(
                wham_cmd,
                cwd=self.analysis_dir,
                capture_output=True,
                text=True,
                check=True,
                timeout=3600
            )
            
            logger.info("WHAM analysis completed successfully")
            
            # Calculate binding energy
            pmf_file = self.analysis_dir / "pmf.xvg"
            binding_energy = self._calculate_binding_energy(pmf_file)
            
            return {
                'pmf_file': str(pmf_file),
                'error_file': str(self.analysis_dir / "pmferror.xvg"),
                'profile_file': str(self.analysis_dir / "bsprofile.xvg"),
                'binding_energy': binding_energy,
                'n_windows': wham_files['n_windows'],
                'begin_time': begin_time,
                'bootstrap': bootstrap,
                'energy_unit': energy_unit
            }
            
        except subprocess.CalledProcessError as e:
            logger.error(f"WHAM analysis failed: {e.stderr}")
            raise RuntimeError(f"WHAM analysis failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            logger.error("WHAM analysis timed out")
            raise RuntimeError("WHAM analysis timed out")
    
    def _calculate_binding_energy(self, pmf_file):
        """Calculate binding energy from PMF data"""
        try:
            data = np.loadtxt(pmf_file, comments=['#', '@'])
            distances = data[:, 0]
            pmf_values = data[:, 1]
            
            min_pmf = np.min(pmf_values)
            max_pmf = np.max(pmf_values)
            binding_energy = max_pmf - min_pmf
            
            min_distance = distances[np.argmin(pmf_values)]
            max_distance = distances[np.argmax(pmf_values)]
            
            return {
                'value': binding_energy,
                'min_pmf': min_pmf,
                'max_pmf': max_pmf,
                'min_distance': min_distance,
                'max_distance': max_distance
            }
            
        except Exception as e:
            logger.warning(f"Could not calculate binding energy: {e}")
            return None
    
    def _plot_pmf_curve(self, pmf_file, error_file):
        """Create main PMF curve plot"""
        try:
            # Read PMF data
            pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
            distances = pmf_data[:, 0]
            pmf_values = pmf_data[:, 1]
            
            # Read error data if available
            errors = None
            if Path(error_file).exists():
                try:
                    error_data = np.loadtxt(error_file, comments=['#', '@'])
                    if error_data.shape[1] >= 2:
                        errors = error_data[:, 1]
                except:
                    pass
            
            # Create plot
            plt.figure(figsize=(12, 8))
            
            if errors is not None:
                plt.errorbar(distances, pmf_values, yerr=errors,
                           fmt='o-', capsize=3, linewidth=2, markersize=4,
                           color='blue', ecolor='lightblue', label='PMF ± Error')
            else:
                plt.plot(distances, pmf_values, 'o-', linewidth=2, markersize=4,
                        color='blue', label='PMF')
            
            plt.xlabel('Distance (nm)', fontsize=14)
            plt.ylabel('PMF (kcal/mol)', fontsize=14)
            plt.title('Potential of Mean Force', fontsize=16)
            plt.grid(True, alpha=0.3)
            plt.legend()
            
            # Add binding energy annotation
            min_idx = np.argmin(pmf_values)
            max_idx = np.argmax(pmf_values)
            binding_energy = pmf_values[max_idx] - pmf_values[min_idx]
            
            stats_text = f'Binding Energy: {binding_energy:.2f} kcal/mol\n'
            stats_text += f'Minimum: {distances[min_idx]:.3f} nm\n'
            stats_text += f'Range: {distances[0]:.2f} - {distances[-1]:.2f} nm'
            
            plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
                    fontsize=12, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
            
            plt.tight_layout()
            
            # Save plot
            plot_file = self.analysis_dir / "pmf_curve.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"PMF curve plot saved: {plot_file}")
            return str(plot_file)
            
        except Exception as e:
            logger.warning(f"Could not create PMF curve plot: {e}")
            return None
    
    def _plot_force_profile(self, pmf_file):
        """Create force profile plot from PMF gradient"""
        try:
            # Read PMF data
            pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
            distances = pmf_data[:, 0]
            pmf_values = pmf_data[:, 1]
            
            # Calculate force as negative gradient
            force_profile = -np.gradient(pmf_values, distances)
            
            # Create plot
            plt.figure(figsize=(12, 8))
            plt.plot(distances, force_profile, 'r-', linewidth=2, label='Force Profile')
            plt.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='Zero Force')
            
            plt.xlabel('Distance (nm)', fontsize=14)
            plt.ylabel('Force (kcal/mol/nm)', fontsize=14)
            plt.title('Force Profile (PMF Gradient)', fontsize=16)
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            
            # Save plot
            plot_file = self.analysis_dir / "force_profile.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Force profile plot saved: {plot_file}")
            return str(plot_file)
            
        except Exception as e:
            logger.warning(f"Could not create force profile plot: {e}")
            return None
    
    def _plot_energy_landscape(self, pmf_file):
        """Create energy landscape visualization"""
        try:
            # Read PMF data
            pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
            distances = pmf_data[:, 0]
            pmf_values = pmf_data[:, 1]
            
            # Create plot
            plt.figure(figsize=(12, 8))
            
            # Fill area under curve
            plt.fill_between(distances, pmf_values, alpha=0.3, color='lightblue', 
                           label='Energy Surface')
            plt.plot(distances, pmf_values, 'b-', linewidth=3, label='PMF')
            
            # Mark binding region (within 2 kcal/mol of minimum)
            min_idx = np.argmin(pmf_values)
            min_pmf = pmf_values[min_idx]
            threshold = min_pmf + 2.0
            
            bound_mask = pmf_values <= threshold
            if np.any(bound_mask):
                bound_distances = distances[bound_mask]
                plt.axvspan(bound_distances[0], bound_distances[-1], 
                          alpha=0.2, color='green', 
                          label='Favorable Binding (< 2 kcal/mol)')
            
            plt.xlabel('Distance (nm)', fontsize=14)
            plt.ylabel('PMF (kcal/mol)', fontsize=14)
            plt.title('Energy Landscape for Protein-Ligand Binding', fontsize=16)
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            
            # Save plot
            plot_file = self.analysis_dir / "energy_landscape.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Energy landscape plot saved: {plot_file}")
            return str(plot_file)
            
        except Exception as e:
            logger.warning(f"Could not create energy landscape plot: {e}")
            return None


class WhamRunner:
    """Standalone WHAM runner for advanced analysis"""
    
    def __init__(self, tpr_files, pullf_files, output_dir):
        self.tpr_files = tpr_files
        self.pullf_files = pullf_files
        self.output_dir = Path(output_dir)
    
    def run(self, **kwargs):
        """Run WHAM with custom parameters"""
        # Implementation for advanced WHAM analysis
        pass
    
    def bootstrap_analysis(self, n_iterations=100):
        """Run bootstrap error analysis"""
        # Implementation for bootstrap analysis
        pass