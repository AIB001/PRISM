#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydrogen bond plotting utilities for PRISM analysis module.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple
from pathlib import Path

# Import publication style
from .publication_utils import apply_publication_style, get_standard_figsize


def plot_hbond_analysis(hbond_results: Dict,
                       output_path: str,
                       title: str = "") -> bool:
    """
    Plot hydrogen bond analysis results for multiple trajectories.

    Parameters
    ----------
    hbond_results : dict
        Dictionary with trajectory names as keys and hydrogen bond data as values
    output_path : str
        Path to save the plot
    title : str
        Plot title

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        print("🔗 Hydrogen bond 4-panel analysis: time series, distribution, averages, and statistics")
        apply_publication_style()

        fig, axes = plt.subplots(2, 2, figsize=get_standard_figsize('quad'))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

        # Plot 1: H-bond count time series
        ax1 = axes[0, 0]
        for i, (traj_name, hbond_data) in enumerate(hbond_results.items()):
            if isinstance(hbond_data, dict) and 'hbond_count' in hbond_data:
                hbond_count = hbond_data['hbond_count']
                if len(hbond_count) > 0:
                    time = np.arange(len(hbond_count)) * 0.5  # Assume 0.5ns interval
                    ax1.plot(time, hbond_count, color=colors[i % 3],
                            label=traj_name, linewidth=2, alpha=0.8)

        ax1.set_xlabel('Time (ns)', fontfamily='Times New Roman')
        ax1.set_ylabel('H-bonds', fontfamily='Times New Roman')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: H-bond count distribution
        ax2 = axes[0, 1]
        all_counts = []
        for hbond_data in hbond_results.values():
            if isinstance(hbond_data, dict) and 'hbond_count' in hbond_data:
                hbond_count = hbond_data['hbond_count']
                if len(hbond_count) > 0:
                    all_counts.extend(hbond_count)

        if all_counts:
            ax2.hist(all_counts, bins=15, alpha=0.7, color='lightgreen', edgecolor='black')
            ax2.set_xlabel('Hydrogen Bond Count', fontfamily='Times New Roman')
            ax2.set_ylabel('Frequency', fontfamily='Times New Roman')
            ax2.grid(True, alpha=0.3)

        # Plot 3: Average H-bond count per trajectory
        ax3 = axes[1, 0]
        traj_names = []
        avg_counts = []
        std_counts = []

        for traj_name, hbond_data in hbond_results.items():
            if isinstance(hbond_data, dict) and 'hbond_count' in hbond_data:
                hbond_count = hbond_data['hbond_count']
                if len(hbond_count) > 0:
                    traj_names.append(traj_name)
                    avg_counts.append(np.mean(hbond_count))
                    std_counts.append(np.std(hbond_count))

        if traj_names:
            bars = ax3.bar(traj_names, avg_counts, yerr=std_counts,
                          color=colors[:len(traj_names)], alpha=0.7, capsize=5)
            ax3.set_ylabel('Avg H-bonds', fontfamily='Times New Roman')
            ax3.grid(True, alpha=0.3, axis='y')

        # Plot 4: H-bond distance distribution (if available)
        ax4 = axes[1, 1]
        all_distances = []
        for hbond_data in hbond_results.values():
            if isinstance(hbond_data, dict) and 'distances' in hbond_data:
                distances = hbond_data['distances']
                if len(distances) > 0:
                    all_distances.extend(distances)

        if all_distances:
            ax4.hist(all_distances, bins=30, alpha=0.7, color='orange', edgecolor='black')
            ax4.set_xlabel('H-bond Distance (Å)', fontfamily='Times New Roman')
            ax4.set_ylabel('Frequency', fontfamily='Times New Roman')
            ax4.grid(True, alpha=0.3)
        else:
            ax4.text(0.5, 0.5, 'No Distance Data', ha='center', va='center',
                    transform=ax4.transAxes, fontsize=16, fontfamily='Times New Roman')
            ax4.set_xticks([])
            ax4.set_yticks([])

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in hydrogen bond plotting: {e}")
        return False


def plot_key_residue_hbonds(hbond_results: Dict,
                           key_residues: List[str],
                           output_path: str,
                           title: str = "Key Residue Hydrogen Bonds") -> bool:
    """
    Plot hydrogen bond analysis for specific key residues.

    Parameters
    ----------
    hbond_results : dict
        Dictionary with trajectory names as keys and hydrogen bond data as values
    key_residues : list
        List of key residue names to highlight
    output_path : str
        Path to save the plot
    title : str
        Plot title

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        print("🔗 Key residue hydrogen bond analysis: frequency and stability comparison")
        apply_publication_style()

        fig, axes = plt.subplots(1, 2, figsize=get_standard_figsize('horizontal'))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

        # Plot 1: Key residue H-bond frequencies
        ax1 = axes[0]
        residue_counts = {res: [] for res in key_residues}

        for traj_name, hbond_data in hbond_results.items():
            if isinstance(hbond_data, dict) and 'residue_hbonds' in hbond_data:
                for res in key_residues:
                    if res in hbond_data['residue_hbonds']:
                        residue_counts[res].append(hbond_data['residue_hbonds'][res])
                    else:
                        residue_counts[res].append(0)

        residue_names = list(residue_counts.keys())
        avg_counts = [np.mean(counts) if counts else 0 for counts in residue_counts.values()]
        std_counts = [np.std(counts) if counts and len(counts) > 1 else 0 for counts in residue_counts.values()]

        if residue_names:
            bars = ax1.bar(residue_names, avg_counts, yerr=std_counts,
                          color='skyblue', alpha=0.7, capsize=5)
            ax1.set_ylabel('Avg H-bonds', fontfamily='Times New Roman')
            ax1.set_xlabel('Key Residues', fontfamily='Times New Roman')
            ax1.grid(True, alpha=0.3, axis='y')
            plt.setp(ax1.get_xticklabels(), rotation=45)

        # Plot 2: H-bond stability over time
        ax2 = axes[1]
        for i, (traj_name, hbond_data) in enumerate(hbond_results.items()):
            if isinstance(hbond_data, dict) and 'hbond_stability' in hbond_data:
                stability = hbond_data['hbond_stability']
                if len(stability) > 0:
                    time = np.arange(len(stability)) * 0.5
                    ax2.plot(time, stability, color=colors[i % 3],
                            label=traj_name, linewidth=2, alpha=0.8)

        ax2.set_xlabel('Time (ns)', fontfamily='Times New Roman')
        ax2.set_ylabel('Stability', fontfamily='Times New Roman')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in key residue hydrogen bond plotting: {e}")
        return False


def plot_hbond_raincloud(hbond_data: Dict[str, Dict[str, float]],
                        output_path: str,
                        title: str = "Key Residue Hydrogen Bond Analysis") -> bool:
    """
    Create raincloud-style plots for key residue hydrogen bonds.

    Creates raincloud plots showing H-bond frequency distributions for key residues
    interacting with the LIG ligand, with violin, box, and scatter components.

    Parameters
    ----------
    hbond_data : dict
        Dictionary with residue names as keys and H-bond frequency data as nested dict
        Format: {'ASP623': {'Repeat1': 0.85, 'Repeat2': 0.92, 'Repeat3': 0.78}, ...}
    output_path : str
        Path to save the figure
    title : str
        Main plot title

    Returns
    -------
    bool
        True if successful
    """
    try:
        print("🎻 Key residue hydrogen bond violin plot showing distribution comparison across trajectories")
        apply_publication_style()

        # Define key residues for H-bond analysis
        key_residues = ['ASP618', 'ASP623', 'ASP760', 'ASN691', 'SER759', 'THR680', 'LYS551', 'ARG553', 'ARG555']

        # Filter to only include key residues that have data
        available_residues = [res for res in key_residues if res in hbond_data]

        if not available_residues:
            print("No key residue H-bond data available for raincloud plot")
            return False

        # Prepare data for raincloud plot
        all_data = []
        residue_labels = []

        for residue in available_residues:
            # Convert frequencies to percentages for consistency
            values = [hbond_data[residue].get(f'Repeat{i}', 0.0) * 100 for i in [1, 2, 3]]
            # Remove zero values to avoid empty distributions
            non_zero_values = [v for v in values if v > 0]
            if non_zero_values:
                all_data.append(non_zero_values)
                residue_labels.append(residue)
            else:
                # If all values are zero, add a small value for visualization
                all_data.append([0.1] * 3)
                residue_labels.append(residue)

        if not all_data:
            return False

        # Create raincloud plot
        fig, ax = plt.subplots(figsize=get_standard_figsize('single'))

        # Color scheme for the raincloud plot
        colors = [
            '#E8F4F8',  # Light blue
            '#F0E8F4',  # Light purple
            '#E8F8F0',  # Light green
            '#F8F0E8',  # Light orange
            '#F4E8E8',  # Light red
            '#F5F5DC',  # Beige
            '#E6E6FA',  # Lavender
            '#F0FFF0',  # Honeydew
            '#FFF5EE'   # Seashell
        ]

        edge_colors = [
            '#4A90A4',  # Dark blue
            '#8E44AD',  # Dark purple
            '#27AE60',  # Dark green
            '#E67E22',  # Dark orange
            '#E74C3C',  # Dark red
            '#D2B48C',  # Tan
            '#9370DB',  # Medium purple
            '#2E8B57',  # Sea green
            '#CD853F'   # Peru
        ]

        positions = np.arange(len(residue_labels))

        # 1. Create half violin plots (right side)
        for i, data in enumerate(all_data):
            # Create violin plot manually for half-violin effect
            if len(data) > 1:
                try:
                    from scipy import stats
                    density = stats.gaussian_kde(data)
                    y_data = np.linspace(min(data), max(data), 100)
                    x_data = density(y_data)

                    # Normalize and scale the violin width
                    x_data = x_data / np.max(x_data) * 0.3

                    # Plot right half of violin
                    ax.fill_betweenx(y_data, positions[i], positions[i] + x_data,
                                   color=colors[i % len(colors)], alpha=0.7, zorder=1)
                except:
                    # Fallback to simple violin if KDE fails
                    parts = ax.violinplot([data], positions=[positions[i]],
                                         showmeans=False, showmedians=False, showextrema=False)
                    for pc in parts['bodies']:
                        pc.set_facecolor(colors[i % len(colors)])
                        pc.set_alpha(0.7)
                        # Modify to make it half-violin
                        verts = pc.get_paths()[0].vertices
                        verts[:, 0] = np.where(verts[:, 0] < positions[i], positions[i], verts[:, 0])

        # 2. Create box plots
        for i, data in enumerate(all_data):
            pos = positions[i]
            if len(data) > 0:
                q1, median, q3 = np.percentile(data, [25, 50, 75])
                iqr = q3 - q1
                lower_whisker = max(min(data), q1 - 1.5 * iqr)
                upper_whisker = min(max(data), q3 + 1.5 * iqr)

                box_width = 0.1
                box = plt.Rectangle((pos - box_width / 2, q1), box_width, q3 - q1,
                                    facecolor='white', edgecolor=edge_colors[i % len(edge_colors)],
                                    linewidth=2, alpha=0.7, zorder=3)
                ax.add_patch(box)

                # Median line
                ax.plot([pos - box_width / 2, pos + box_width / 2], [median, median],
                        color=edge_colors[i % len(edge_colors)], linewidth=3, zorder=4)

                # Whiskers
                ax.plot([pos, pos], [q3, upper_whisker], color=edge_colors[i % len(edge_colors)], linewidth=2, zorder=3)
                ax.plot([pos, pos], [q1, lower_whisker], color=edge_colors[i % len(edge_colors)], linewidth=2, zorder=3)
                ax.plot([pos - 0.02, pos + 0.02], [upper_whisker, upper_whisker],
                        color=edge_colors[i % len(edge_colors)], linewidth=2, zorder=3)
                ax.plot([pos - 0.02, pos + 0.02], [lower_whisker, lower_whisker],
                        color=edge_colors[i % len(edge_colors)], linewidth=2, zorder=3)

        # 3. Create scatter plots (rain) - left side
        for i, data in enumerate(all_data):
            pos = positions[i]
            n_points = len(data)
            if n_points > 0:
                y_data = np.array(data)

                # Add jitter to x positions for scatter
                x_jitter = np.random.normal(pos - 0.15, 0.05, n_points)  # Left side with jitter
                x_jitter = np.clip(x_jitter, pos - 0.3, pos - 0.05)  # Keep within bounds

                ax.scatter(x_jitter, y_data, s=30, alpha=0.6,
                          color=edge_colors[i % len(edge_colors)], edgecolors='white',
                          linewidth=0.5, zorder=2)

        # 4. Formatting with proper tick alignment
        ax.set_xticks(positions)
        ax.set_xticklabels(residue_labels, rotation=45, ha='center', fontsize=18)
        ax.set_ylabel('Frequency (%)', fontsize=21, fontweight='bold')
        ax.set_xlabel('Key Residues', fontsize=21, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#FAFAFA')

        # Set y-axis limits
        all_values = [item for sublist in all_data for item in sublist]
        if all_values:
            ax.set_ylim(0, max(all_values) * 1.1)

        # Style axes
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_axisbelow(True)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in H-bond raincloud plot: {e}")
        return False


def plot_hbond_timeseries_multi_trajectory(timeseries_data: Dict[str, Dict[str, np.ndarray]],
                                           output_path: str,
                                           time_unit: str = 'ns',
                                           timestep_ps: float = 20.0,
                                           title: str = "") -> bool:
    """
    Plot hydrogen bond time series for multiple trajectories on the same plot.

    Parameters
    ----------
    timeseries_data : dict
        Nested dictionary: {trajectory_name: {residue_name: boolean_array}}
        Example: {'Repeat 1': {'ASP623': array([True, False, ...]), ...}, ...}
    output_path : str
        Path to save the plot
    time_unit : str
        Time unit for x-axis ('ns' or 'ps')
    timestep_ps : float
        Timestep in picoseconds (default: 20.0 ps)
    title : str
        Plot title

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        print("📈 H-bond time series: multi-trajectory comparison for key residues")
        apply_publication_style()

        # Extract all unique residues
        all_residues = set()
        for traj_data in timeseries_data.values():
            all_residues.update(traj_data.keys())

        if not all_residues:
            print("  ⚠ No H-bond time series data available")
            return False

        all_residues = sorted(all_residues)
        n_residues = len(all_residues)

        # Create subplots (rows for residues)
        n_cols = 2
        n_rows = (n_residues + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(12 * n_cols, 5 * n_rows))
        if n_rows == 1 and n_cols == 1:
            axes = np.array([[axes]])
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)

        # Colors for each trajectory
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
        traj_names = list(timeseries_data.keys())

        # Convert timestep to appropriate unit
        if time_unit == 'ns':
            time_conversion = timestep_ps / 1000.0  # ps to ns
            xlabel = 'Time (ns)'
        else:
            time_conversion = timestep_ps
            xlabel = 'Time (ps)'

        # Plot each residue
        for idx, residue in enumerate(all_residues):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row, col]

            # Plot each trajectory
            for traj_idx, traj_name in enumerate(traj_names):
                traj_data = timeseries_data[traj_name]

                if residue in traj_data:
                    hbond_array = traj_data[residue].astype(float)  # Convert bool to float
                    n_frames = len(hbond_array)
                    time_points = np.arange(n_frames) * time_conversion

                    # Plot as line (0/1 for H-bond absence/presence)
                    ax.plot(time_points, hbond_array, label=traj_name,
                           color=colors[traj_idx % len(colors)],
                           linewidth=1.5, alpha=0.8)

            ax.set_xlabel(xlabel, fontfamily='Times New Roman', fontweight='bold')
            ax.set_ylabel('H-bond Occupancy', fontfamily='Times New Roman', fontweight='bold')
            ax.text(0.02, 0.98, residue, transform=ax.transAxes,
                   fontfamily='Times New Roman', fontweight='bold', fontsize=14,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            ax.set_ylim(-0.1, 1.1)
            ax.set_yticks([0, 1])
            ax.set_yticklabels(['No', 'Yes'])
            ax.legend(loc='upper right', fontsize=10)
            ax.grid(True, alpha=0.3, axis='x')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        # Hide unused subplots
        for idx in range(n_residues, n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].set_visible(False)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        print(f"  ✓ Saved H-bond time series plot: {output_path}")
        return True

    except Exception as e:
        print(f"Error in H-bond time series plotting: {e}")
        import traceback
        traceback.print_exc()
        return False