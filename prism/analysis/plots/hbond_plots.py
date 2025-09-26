#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydrogen bond plotting utilities for PRISM analysis module.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
from pathlib import Path

# Import publication style
from .publication_utils import apply_publication_style


def plot_hbond_analysis(hbond_results: Dict,
                       output_path: str,
                       title: str = "Hydrogen Bond Analysis") -> bool:
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
        apply_publication_style()

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
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
        ax1.set_ylabel('Hydrogen Bond Count', fontfamily='Times New Roman')
        ax1.set_title('Hydrogen Bond Count vs Time', fontfamily='Times New Roman', fontweight='bold')
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
            ax2.set_title('H-bond Count Distribution', fontfamily='Times New Roman', fontweight='bold')
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
            ax3.set_ylabel('Average H-bond Count', fontfamily='Times New Roman')
            ax3.set_title('Average H-bond Count per Trajectory', fontfamily='Times New Roman', fontweight='bold')
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
            ax4.set_title('H-bond Distance Distribution', fontfamily='Times New Roman', fontweight='bold')
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
        apply_publication_style()

        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
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
            ax1.set_ylabel('Average H-bond Count', fontfamily='Times New Roman')
            ax1.set_xlabel('Key Residues', fontfamily='Times New Roman')
            ax1.set_title('Key Residue H-bond Frequencies', fontfamily='Times New Roman', fontweight='bold')
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
        ax2.set_ylabel('H-bond Stability', fontfamily='Times New Roman')
        ax2.set_title('H-bond Stability Over Time', fontfamily='Times New Roman', fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in key residue hydrogen bond plotting: {e}")
        return False