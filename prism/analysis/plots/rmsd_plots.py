#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD/RMSF plotting utilities for PRISM analysis module.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
from pathlib import Path

# Import publication style
from .publication_utils import apply_publication_style


def plot_rmsd_analysis(rmsd_results: Dict,
                      output_path: str,
                      title: str = "RMSD Analysis") -> bool:
    """
    Plot RMSD analysis results for multiple trajectories.

    Parameters
    ----------
    rmsd_results : dict
        Dictionary with trajectory names as keys and RMSD data as values
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

        # Plot 1: RMSD time series
        ax1 = axes[0, 0]
        for i, (traj_name, rmsd_data) in enumerate(rmsd_results.items()):
            if isinstance(rmsd_data, dict) and 'protein' in rmsd_data:
                protein_rmsd = rmsd_data['protein']
            else:
                protein_rmsd = rmsd_data

            if len(protein_rmsd) > 0:
                time = np.arange(len(protein_rmsd)) * 0.5  # Assume 0.5ns interval
                ax1.plot(time, protein_rmsd, color=colors[i % 3],
                        label=traj_name, linewidth=2, alpha=0.8)

        ax1.set_xlabel('Time (ns)', fontfamily='Times New Roman')
        ax1.set_ylabel('RMSD (Å)', fontfamily='Times New Roman')
        ax1.set_title('Protein RMSD vs Time', fontfamily='Times New Roman', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: RMSD distribution
        ax2 = axes[0, 1]
        all_rmsd = []
        for rmsd_data in rmsd_results.values():
            if isinstance(rmsd_data, dict) and 'protein' in rmsd_data:
                protein_rmsd = rmsd_data['protein']
            else:
                protein_rmsd = rmsd_data
            if len(protein_rmsd) > 0:
                all_rmsd.extend(protein_rmsd)

        if all_rmsd:
            ax2.hist(all_rmsd, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
            ax2.set_xlabel('RMSD (Å)', fontfamily='Times New Roman')
            ax2.set_ylabel('Frequency', fontfamily='Times New Roman')
            ax2.set_title('RMSD Distribution', fontfamily='Times New Roman', fontweight='bold')
            ax2.grid(True, alpha=0.3)

        # Plot 3: Average RMSD per trajectory
        ax3 = axes[1, 0]
        traj_names = []
        avg_rmsd = []
        std_rmsd = []

        for traj_name, rmsd_data in rmsd_results.items():
            if isinstance(rmsd_data, dict) and 'protein' in rmsd_data:
                protein_rmsd = rmsd_data['protein']
            else:
                protein_rmsd = rmsd_data

            if len(protein_rmsd) > 0:
                traj_names.append(traj_name)
                avg_rmsd.append(np.mean(protein_rmsd))
                std_rmsd.append(np.std(protein_rmsd))

        if traj_names:
            bars = ax3.bar(traj_names, avg_rmsd, yerr=std_rmsd,
                          color=colors[:len(traj_names)], alpha=0.7, capsize=5)
            ax3.set_ylabel('Average RMSD (Å)', fontfamily='Times New Roman')
            ax3.set_title('Average RMSD per Trajectory', fontfamily='Times New Roman', fontweight='bold')
            ax3.grid(True, alpha=0.3, axis='y')

        # Plot 4: Ligand RMSD (if available)
        ax4 = axes[1, 1]
        has_ligand = False

        for i, (traj_name, rmsd_data) in enumerate(rmsd_results.items()):
            if isinstance(rmsd_data, dict) and 'ligand' in rmsd_data:
                ligand_rmsd = rmsd_data['ligand']
                if len(ligand_rmsd) > 0:
                    time = np.arange(len(ligand_rmsd)) * 0.5
                    ax4.plot(time, ligand_rmsd, color=colors[i % 3],
                            label=traj_name, linewidth=2, alpha=0.8)
                    has_ligand = True

        if has_ligand:
            ax4.set_xlabel('Time (ns)', fontfamily='Times New Roman')
            ax4.set_ylabel('Ligand RMSD (Å)', fontfamily='Times New Roman')
            ax4.set_title('Ligand RMSD vs Time', fontfamily='Times New Roman', fontweight='bold')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        else:
            ax4.text(0.5, 0.5, 'No Ligand Data', ha='center', va='center',
                    transform=ax4.transAxes, fontsize=16, fontfamily='Times New Roman')
            ax4.set_xticks([])
            ax4.set_yticks([])

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in RMSD plotting: {e}")
        return False


def plot_rmsf_analysis(rmsf_results: Dict,
                      output_path: str,
                      title: str = "RMSF Analysis") -> bool:
    """
    Plot RMSF analysis results.

    Parameters
    ----------
    rmsf_results : dict
        Dictionary with trajectory names as keys and RMSF data as values
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

        # Plot 1: RMSF per residue
        ax1 = axes[0]
        for i, (traj_name, rmsf_data) in enumerate(rmsf_results.items()):
            if isinstance(rmsf_data, tuple) and len(rmsf_data) >= 2:
                residue_indices, rmsf_values = rmsf_data[:2]
            elif isinstance(rmsf_data, dict) and 'rmsf' in rmsf_data:
                rmsf_values = rmsf_data['rmsf']
                residue_indices = np.arange(len(rmsf_values))
            else:
                continue

            if len(rmsf_values) > 0:
                ax1.plot(residue_indices, rmsf_values, color=colors[i % 3],
                        label=traj_name, linewidth=2, alpha=0.8)

        ax1.set_xlabel('Residue Index', fontfamily='Times New Roman')
        ax1.set_ylabel('RMSF (Å)', fontfamily='Times New Roman')
        ax1.set_title('RMSF per Residue', fontfamily='Times New Roman', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: RMSF distribution
        ax2 = axes[1]
        all_rmsf = []
        for rmsf_data in rmsf_results.values():
            if isinstance(rmsf_data, tuple) and len(rmsf_data) >= 2:
                rmsf_values = rmsf_data[1]
            elif isinstance(rmsf_data, dict) and 'rmsf' in rmsf_data:
                rmsf_values = rmsf_data['rmsf']
            else:
                continue

            if len(rmsf_values) > 0:
                all_rmsf.extend(rmsf_values)

        if all_rmsf:
            ax2.hist(all_rmsf, bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
            ax2.set_xlabel('RMSF (Å)', fontfamily='Times New Roman')
            ax2.set_ylabel('Frequency', fontfamily='Times New Roman')
            ax2.set_title('RMSF Distribution', fontfamily='Times New Roman', fontweight='bold')
            ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in RMSF plotting: {e}")
        return False