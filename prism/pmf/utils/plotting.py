#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF Plotting Utilities

Unified plotting functions for all PMF analysis stages:
- SMD: Force and distance profiles
- Umbrella: Sampling quality and convergence
- WHAM: PMF curves and energy landscapes
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional, Dict, Tuple, List
import logging

logger = logging.getLogger(__name__)


class PMFPlotter:
    """Unified plotting utility for PMF workflows"""

    def __init__(self, output_dir: Path, dpi: int = 300, style: str = 'default'):
        """
        Initialize plotter

        Args:
            output_dir: Directory to save plots
            dpi: Resolution for saved figures
            style: Matplotlib style to use
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.dpi = dpi
        self.style = style

        if style != 'default':
            plt.style.use(style)

    # ========================================================================
    # SMD Plotting Functions
    # ========================================================================

    def plot_smd_force_time(self, pullf_file: Path, output_name: str = 'force_vs_time.png') -> Optional[str]:
        """
        Plot SMD pulling force vs time

        Args:
            pullf_file: Path to smd_pullf.xvg
            output_name: Output filename

        Returns:
            Path to saved plot or None if failed
        """
        try:
            data = np.loadtxt(pullf_file, comments=['#', '@'])
            time = data[:, 0]  # ps
            force = data[:, 1]  # kJ/mol/nm

            plt.figure(figsize=(10, 6))
            plt.plot(time, force, linewidth=1.5, color='blue', alpha=0.8)
            plt.xlabel('Time (ps)', fontsize=12)
            plt.ylabel('Force (kJ/mol/nm)', fontsize=12)
            plt.title('SMD Pulling Force vs Time', fontsize=14)
            plt.grid(True, alpha=0.3)

            # Add statistics annotation
            stats_text = f'Max: {np.max(force):.2f}\nAvg: {np.mean(force):.2f}\nMin: {np.min(force):.2f}'
            plt.text(0.98, 0.98, stats_text, transform=plt.gca().transAxes,
                    fontsize=10, verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            plt.tight_layout()
            plot_path = self.output_dir / output_name
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logger.info(f"SMD force-time plot saved: {plot_path}")
            return str(plot_path)

        except Exception as e:
            logger.warning(f"Could not create SMD force-time plot: {e}")
            return None

    def plot_smd_distance_time(self, pullx_file: Path, output_name: str = 'distance_vs_time.png') -> Optional[str]:
        """
        Plot SMD pulling distance vs time

        Args:
            pullx_file: Path to smd_pullx.xvg
            output_name: Output filename

        Returns:
            Path to saved plot or None if failed
        """
        try:
            data = np.loadtxt(pullx_file, comments=['#', '@'])
            time = data[:, 0]  # ps
            distance = data[:, 1]  # nm

            plt.figure(figsize=(10, 6))
            plt.plot(time, distance, linewidth=1.5, color='green', alpha=0.8)
            plt.xlabel('Time (ps)', fontsize=12)
            plt.ylabel('Distance (nm)', fontsize=12)
            plt.title('SMD Pulling Distance vs Time', fontsize=14)
            plt.grid(True, alpha=0.3)

            # Add statistics
            total_dist = distance[-1] - distance[0]
            stats_text = f'Total: {total_dist:.3f} nm\nFinal: {distance[-1]:.3f} nm'
            plt.text(0.98, 0.02, stats_text, transform=plt.gca().transAxes,
                    fontsize=10, verticalalignment='bottom', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

            plt.tight_layout()
            plot_path = self.output_dir / output_name
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logger.info(f"SMD distance-time plot saved: {plot_path}")
            return str(plot_path)

        except Exception as e:
            logger.warning(f"Could not create SMD distance-time plot: {e}")
            return None

    def plot_smd_force_distance(self, pullf_file: Path, pullx_file: Path,
                               output_name: str = 'force_vs_distance.png') -> Optional[str]:
        """
        Plot SMD force vs distance

        Args:
            pullf_file: Path to smd_pullf.xvg
            pullx_file: Path to smd_pullx.xvg
            output_name: Output filename

        Returns:
            Path to saved plot or None if failed
        """
        try:
            pullf_data = np.loadtxt(pullf_file, comments=['#', '@'])
            pullx_data = np.loadtxt(pullx_file, comments=['#', '@'])

            time_f = pullf_data[:, 0]
            force = pullf_data[:, 1]

            time_x = pullx_data[:, 0]
            distance = pullx_data[:, 1]

            # Align data by interpolation
            from scipy.interpolate import interp1d
            f_interp = interp1d(time_f, force, kind='linear', fill_value='extrapolate')

            # Use distance time points as reference
            common_time = time_x[time_x <= time_f[-1]]
            force_aligned = f_interp(common_time)
            dist_aligned = distance[:len(common_time)]

            plt.figure(figsize=(10, 6))
            plt.plot(dist_aligned, force_aligned, linewidth=1.5, color='red', alpha=0.8)
            plt.xlabel('Distance (nm)', fontsize=12)
            plt.ylabel('Force (kJ/mol/nm)', fontsize=12)
            plt.title('SMD Force vs Distance', fontsize=14)
            plt.grid(True, alpha=0.3)

            plt.tight_layout()
            plot_path = self.output_dir / output_name
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logger.info(f"SMD force-distance plot saved: {plot_path}")
            return str(plot_path)

        except Exception as e:
            logger.warning(f"Could not create SMD force-distance plot: {e}")
            return None

    def plot_smd_work_profile(self, pullf_file: Path, pullx_file: Path,
                             output_name: str = 'work_profile.png') -> Optional[str]:
        """
        Plot SMD work profile (integral of force over distance)

        Args:
            pullf_file: Path to smd_pullf.xvg
            pullx_file: Path to smd_pullx.xvg
            output_name: Output filename

        Returns:
            Path to saved plot or None if failed
        """
        try:
            pullf_data = np.loadtxt(pullf_file, comments=['#', '@'])
            pullx_data = np.loadtxt(pullx_file, comments=['#', '@'])

            time_f = pullf_data[:, 0]
            force = pullf_data[:, 1]

            time_x = pullx_data[:, 0]
            distance = pullx_data[:, 1]

            # Align data
            from scipy.interpolate import interp1d
            from scipy.integrate import cumtrapz

            f_interp = interp1d(time_f, force, kind='linear', fill_value='extrapolate')
            common_time = time_x[time_x <= time_f[-1]]
            force_aligned = f_interp(common_time)
            dist_aligned = distance[:len(common_time)]

            # Calculate work as integral of force
            work = cumtrapz(force_aligned, dist_aligned, initial=0)

            plt.figure(figsize=(10, 6))
            plt.plot(dist_aligned, work, linewidth=2, color='purple', alpha=0.8)
            plt.xlabel('Distance (nm)', fontsize=12)
            plt.ylabel('Work (kJ/mol)', fontsize=12)
            plt.title('SMD Work Profile', fontsize=14)
            plt.grid(True, alpha=0.3)

            # Add final work value
            final_work = work[-1]
            plt.text(0.98, 0.98, f'Total Work: {final_work:.2f} kJ/mol',
                    transform=plt.gca().transAxes,
                    fontsize=12, verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.5))

            plt.tight_layout()
            plot_path = self.output_dir / output_name
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logger.info(f"SMD work profile plot saved: {plot_path}")
            return str(plot_path)

        except Exception as e:
            logger.warning(f"Could not create SMD work profile plot: {e}")
            return None

    def plot_all_smd(self, pullf_file: Path, pullx_file: Path) -> Dict[str, str]:
        """
        Generate all SMD plots

        Args:
            pullf_file: Path to smd_pullf.xvg
            pullx_file: Path to smd_pullx.xvg

        Returns:
            Dictionary of plot names and paths
        """
        plots = {}

        plots['force_time'] = self.plot_smd_force_time(pullf_file)
        plots['distance_time'] = self.plot_smd_distance_time(pullx_file)
        plots['force_distance'] = self.plot_smd_force_distance(pullf_file, pullx_file)
        plots['work_profile'] = self.plot_smd_work_profile(pullf_file, pullx_file)

        # Filter out None values
        plots = {k: v for k, v in plots.items() if v is not None}

        logger.info(f"Generated {len(plots)} SMD plots")
        return plots

    # ========================================================================
    # PMF/WHAM Plotting Functions
    # ========================================================================

    def plot_pmf_curve(self, pmf_file: Path, error_file: Optional[Path] = None,
                      output_name: str = 'pmf_curve.png') -> Optional[str]:
        """
        Plot PMF curve with optional error bars

        Args:
            pmf_file: Path to PMF profile file
            error_file: Path to PMF error file (optional)
            output_name: Output filename

        Returns:
            Path to saved plot or None if failed
        """
        try:
            pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
            distances = pmf_data[:, 0]
            pmf_values = pmf_data[:, 1]

            # Read error data if available
            errors = None
            if error_file and Path(error_file).exists():
                try:
                    error_data = np.loadtxt(error_file, comments=['#', '@'])
                    if error_data.shape[1] >= 2:
                        errors = error_data[:, 1]
                except Exception as e:
                    logger.warning(f"Could not read error file: {e}")

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
            plot_path = self.output_dir / output_name
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logger.info(f"PMF curve plot saved: {plot_path}")
            return str(plot_path)

        except Exception as e:
            logger.warning(f"Could not create PMF curve plot: {e}")
            return None

    def plot_force_profile(self, pmf_file: Path, output_name: str = 'force_profile.png') -> Optional[str]:
        """
        Plot force profile (negative gradient of PMF)

        Args:
            pmf_file: Path to PMF profile file
            output_name: Output filename

        Returns:
            Path to saved plot or None if failed
        """
        try:
            pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
            distances = pmf_data[:, 0]
            pmf_values = pmf_data[:, 1]

            # Calculate force as negative gradient
            force_profile = -np.gradient(pmf_values, distances)

            plt.figure(figsize=(12, 8))
            plt.plot(distances, force_profile, 'r-', linewidth=2, label='Force Profile')
            plt.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='Zero Force')

            plt.xlabel('Distance (nm)', fontsize=14)
            plt.ylabel('Force (kcal/mol/nm)', fontsize=14)
            plt.title('Force Profile (PMF Gradient)', fontsize=16)
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()

            plot_path = self.output_dir / output_name
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logger.info(f"Force profile plot saved: {plot_path}")
            return str(plot_path)

        except Exception as e:
            logger.warning(f"Could not create force profile plot: {e}")
            return None

    def plot_energy_landscape(self, pmf_file: Path, output_name: str = 'energy_landscape.png') -> Optional[str]:
        """
        Plot energy landscape visualization

        Args:
            pmf_file: Path to PMF profile file
            output_name: Output filename

        Returns:
            Path to saved plot or None if failed
        """
        try:
            pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
            distances = pmf_data[:, 0]
            pmf_values = pmf_data[:, 1]

            plt.figure(figsize=(12, 8))

            # Fill area under curve
            plt.fill_between(distances, pmf_values, alpha=0.3, color='lightblue',
                           label='Energy Surface')
            plt.plot(distances, pmf_values, 'b-', linewidth=3, label='PMF')

            # Mark binding region (within 2 kcal/mol of minimum)
            min_idx = np.argmin(pmf_values)
            min_value = pmf_values[min_idx]
            binding_threshold = min_value + 2.0

            binding_mask = pmf_values <= binding_threshold
            if np.any(binding_mask):
                plt.fill_between(distances, pmf_values, binding_threshold,
                               where=binding_mask, alpha=0.5, color='lightgreen',
                               label='Binding Region (≤2 kcal/mol)')

            # Mark minimum
            plt.plot(distances[min_idx], pmf_values[min_idx], 'r*',
                    markersize=20, label='Minimum')

            plt.xlabel('Distance (nm)', fontsize=14)
            plt.ylabel('PMF (kcal/mol)', fontsize=14)
            plt.title('Energy Landscape', fontsize=16)
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()

            plot_path = self.output_dir / output_name
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logger.info(f"Energy landscape plot saved: {plot_path}")
            return str(plot_path)

        except Exception as e:
            logger.warning(f"Could not create energy landscape plot: {e}")
            return None

    def plot_all_pmf(self, pmf_file: Path, error_file: Optional[Path] = None) -> Dict[str, str]:
        """
        Generate all PMF plots

        Args:
            pmf_file: Path to PMF profile file
            error_file: Path to PMF error file (optional)

        Returns:
            Dictionary of plot names and paths
        """
        plots = {}

        plots['pmf_curve'] = self.plot_pmf_curve(pmf_file, error_file)
        plots['force_profile'] = self.plot_force_profile(pmf_file)
        plots['energy_landscape'] = self.plot_energy_landscape(pmf_file)

        # Filter out None values
        plots = {k: v for k, v in plots.items() if v is not None}

        logger.info(f"Generated {len(plots)} PMF plots")
        return plots

    # ========================================================================
    # Umbrella Sampling Plotting Functions
    # ========================================================================

    def plot_umbrella_convergence(self, window_dirs: List[Path],
                                  output_name: str = 'umbrella_convergence.png') -> Optional[str]:
        """
        Plot umbrella sampling convergence

        Args:
            window_dirs: List of umbrella window directories
            output_name: Output filename

        Returns:
            Path to saved plot or None if failed
        """
        try:
            distances = []
            sampling_counts = []

            for window_dir in window_dirs:
                # Extract distance from directory name (e.g., "window_1.5nm")
                try:
                    dist_str = window_dir.name.split('_')[1].replace('nm', '')
                    distance = float(dist_str)

                    # Check for pullf file and count data points
                    pullf_file = window_dir / 'umbrella_pullf.xvg'
                    if pullf_file.exists():
                        data = np.loadtxt(pullf_file, comments=['#', '@'])
                        count = len(data)

                        distances.append(distance)
                        sampling_counts.append(count)
                except:
                    continue

            if not distances:
                logger.warning("No umbrella sampling data found for convergence plot")
                return None

            # Sort by distance
            sorted_indices = np.argsort(distances)
            distances = np.array(distances)[sorted_indices]
            sampling_counts = np.array(sampling_counts)[sorted_indices]

            plt.figure(figsize=(10, 6))
            plt.bar(distances, sampling_counts, width=0.05, color='skyblue', edgecolor='blue')
            plt.xlabel('Distance (nm)', fontsize=12)
            plt.ylabel('Sampling Count', fontsize=12)
            plt.title('Umbrella Sampling Coverage', fontsize=14)
            plt.grid(True, alpha=0.3, axis='y')

            plt.tight_layout()
            plot_path = self.output_dir / output_name
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logger.info(f"Umbrella convergence plot saved: {plot_path}")
            return str(plot_path)

        except Exception as e:
            logger.warning(f"Could not create umbrella convergence plot: {e}")
            return None


# Convenience functions for quick plotting

def plot_smd_results(pullf_file: Path, pullx_file: Path, output_dir: Path) -> Dict[str, str]:
    """
    Quick function to generate all SMD plots

    Args:
        pullf_file: Path to smd_pullf.xvg
        pullx_file: Path to smd_pullx.xvg
        output_dir: Directory to save plots

    Returns:
        Dictionary of plot names and paths
    """
    plotter = PMFPlotter(output_dir)
    return plotter.plot_all_smd(pullf_file, pullx_file)


def plot_pmf_results(pmf_file: Path, output_dir: Path, error_file: Optional[Path] = None) -> Dict[str, str]:
    """
    Quick function to generate all PMF plots

    Args:
        pmf_file: Path to PMF profile file
        output_dir: Directory to save plots
        error_file: Path to PMF error file (optional)

    Returns:
        Dictionary of plot names and paths
    """
    plotter = PMFPlotter(output_dir)
    return plotter.plot_all_pmf(pmf_file, error_file)
