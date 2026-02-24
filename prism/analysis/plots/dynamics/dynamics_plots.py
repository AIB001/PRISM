#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand dynamics plotting utilities for PRISM analysis module.

Provides publication-quality plotting functions for ligand dynamics analysis.
"""

import matplotlib.pyplot as plt
from typing import Dict, List, Optional

# Import publication style
from ..core.publication_utils import apply_publication_style, get_standard_figsize


def plot_ligand_dynamics_panel(
    dynamics_data_list: List[Dict], trajectory_names: List[str], output_path: str, title: str = ""
) -> bool:
    """
    Plot comprehensive ligand dynamics analysis in a multi-panel layout.

    Parameters
    ----------
    dynamics_data_list : list of dict
        List of dynamics data dictionaries from LigandDynamicsAnalyzer
    trajectory_names : list of str
        Names of trajectories
    output_path : str
        Path to save the plot
    title : str
        Overall plot title (optional)

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        print("📊 Ligand dynamics 4-panel analysis: RMSD, key distances, pocket volume, and distance distribution")
        apply_publication_style()

        fig, axes = plt.subplots(2, 2, figsize=get_standard_figsize("quad"))
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

        # Panel 1: Ligand RMSD time series
        ax1 = axes[0, 0]
        for i, (data, name) in enumerate(zip(dynamics_data_list, trajectory_names)):
            if data is not None and "ligand_rmsd" in data and "time" in data:
                time_ns = data["time"] / 1000.0  # Convert ps to ns
                ax1.plot(time_ns, data["ligand_rmsd"], color=colors[i % 3], label=name, linewidth=1.5, alpha=0.8)

        ax1.set_xlabel("Time (ns)", fontfamily="Times New Roman", fontsize=14)
        ax1.set_ylabel("Ligand RMSD (Å)", fontfamily="Times New Roman", fontsize=14)
        ax1.legend(frameon=True, fancybox=True, fontsize=12)
        ax1.grid(True, alpha=0.3)

        # Panel 2: Key distance time series (first key residue)
        ax2 = axes[0, 1]
        key_res = None
        for data in dynamics_data_list:
            if data and "distances" in data and data["distances"]:
                key_res = list(data["distances"].keys())[0]
                break

        if key_res:
            for i, (data, name) in enumerate(zip(dynamics_data_list, trajectory_names)):
                if data and "distances" in data and key_res in data["distances"]:
                    time_ns = data["time"] / 1000.0
                    distances_ang = data["distances"][key_res] * 10  # nm to Angstrom
                    ax2.plot(time_ns, distances_ang, color=colors[i % 3], label=name, linewidth=1.5, alpha=0.8)

            ax2.set_xlabel("Time (ns)", fontfamily="Times New Roman", fontsize=14)
            ax2.set_ylabel(f"Distance to {key_res} (Å)", fontfamily="Times New Roman", fontsize=14)
            ax2.legend(frameon=True, fancybox=True, fontsize=12)
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(
                0.5,
                0.5,
                "No Distance Data",
                ha="center",
                va="center",
                transform=ax2.transAxes,
                fontsize=16,
                fontfamily="Times New Roman",
            )

        # Panel 3: Binding pocket volume variation
        ax3 = axes[1, 0]
        has_pocket_data = False
        for i, (data, name) in enumerate(zip(dynamics_data_list, trajectory_names)):
            if data and "pocket_volumes" in data and len(data["pocket_volumes"]) > 0:
                time_ns = data["time"] / 1000.0
                # Normalize pocket volumes for better visualization
                volumes = data["pocket_volumes"] * 1000  # Scale for visibility
                ax3.plot(time_ns, volumes, color=colors[i % 3], label=name, linewidth=1.5, alpha=0.8)
                has_pocket_data = True

        if has_pocket_data:
            ax3.set_xlabel("Time (ns)", fontfamily="Times New Roman", fontsize=14)
            ax3.set_ylabel("Pocket Volume (nm³×10³)", fontfamily="Times New Roman", fontsize=14)
            ax3.legend(frameon=True, fancybox=True, fontsize=12)
            ax3.grid(True, alpha=0.3)
        else:
            ax3.text(
                0.5,
                0.5,
                "No Pocket Volume Data",
                ha="center",
                va="center",
                transform=ax3.transAxes,
                fontsize=16,
                fontfamily="Times New Roman",
            )

        # Panel 4: Distance distribution for all key residues
        ax4 = axes[1, 1]
        all_key_residues = set()
        for data in dynamics_data_list:
            if data and "distances" in data:
                all_key_residues.update(data["distances"].keys())

        if all_key_residues:
            # Plot top 3 key residues for clarity
            top_residues = sorted(list(all_key_residues))[:3]
            colors_res = ["#e74c3c", "#3498db", "#2ecc71"]

            for i, res in enumerate(top_residues):
                distances_all = []
                for data in dynamics_data_list:
                    if data and "distances" in data and res in data["distances"]:
                        distances_ang = data["distances"][res] * 10
                        distances_all.extend(distances_ang)

                if distances_all:
                    ax4.hist(
                        distances_all, bins=20, alpha=0.6, label=f"Lig-{res}", density=True, color=colors_res[i % 3]
                    )

            ax4.set_xlabel("Distance (Å)", fontfamily="Times New Roman", fontsize=14)
            ax4.set_ylabel("Probability Density", fontfamily="Times New Roman", fontsize=14)
            ax4.legend(frameon=True, fancybox=True, fontsize=12)
            ax4.grid(True, alpha=0.3)
        else:
            ax4.text(
                0.5,
                0.5,
                "No Distance Distribution",
                ha="center",
                va="center",
                transform=ax4.transAxes,
                fontsize=16,
                fontfamily="Times New Roman",
            )

        if title:
            plt.suptitle(title, fontsize=16, fontweight="bold")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        plt.close()

        return True

    except Exception as e:
        print(f"Error in ligand dynamics panel plotting: {e}")
        import traceback

        traceback.print_exc()
        return False


def plot_ligand_rmsd_comparison(
    dynamics_data_list: List[Dict], trajectory_names: List[str], output_path: str, title: str = ""
) -> bool:
    """
    Plot ligand RMSD comparison across multiple trajectories.

    Parameters
    ----------
    dynamics_data_list : list of dict
        List of dynamics data dictionaries
    trajectory_names : list of str
        Names of trajectories
    output_path : str
        Path to save the plot
    title : str
        Plot title

    Returns
    -------
    bool
        True if successful
    """
    try:
        print("📈 Ligand RMSD comparison plot showing stability across trajectories")
        apply_publication_style()

        fig, ax = plt.subplots(figsize=get_standard_figsize("single"))
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

        for i, (data, name) in enumerate(zip(dynamics_data_list, trajectory_names)):
            if data and "ligand_rmsd" in data and "time" in data:
                time_ns = data["time"] / 1000.0
                ax.plot(time_ns, data["ligand_rmsd"], color=colors[i % 3], label=name, linewidth=2, alpha=0.8)

        ax.set_xlabel("Time (ns)", fontfamily="Times New Roman", fontsize=16)
        ax.set_ylabel("Ligand RMSD (Å)", fontfamily="Times New Roman", fontsize=16)
        if title:
            ax.set_title(title, fontfamily="Times New Roman", fontsize=18, fontweight="bold")
        ax.legend(frameon=True, fancybox=True, fontsize=14)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        plt.close()

        return True

    except Exception as e:
        print(f"Error in ligand RMSD comparison plotting: {e}")
        return False


def plot_key_distances_timeseries(
    dynamics_data_list: List[Dict],
    trajectory_names: List[str],
    output_path: str,
    key_residues: Optional[List[str]] = None,
    title: str = "",
) -> bool:
    """
    Plot time series of distances to key residues.

    Parameters
    ----------
    dynamics_data_list : list of dict
        List of dynamics data dictionaries
    trajectory_names : list of str
        Names of trajectories
    output_path : str
        Path to save the plot
    key_residues : list of str, optional
        Specific residues to plot (if None, plot all available)
    title : str
        Plot title

    Returns
    -------
    bool
        True if successful
    """
    try:
        print("📏 Key residue distance time series showing ligand-protein interactions")
        apply_publication_style()

        # Collect all available residues
        all_residues = set()
        for data in dynamics_data_list:
            if data and "distances" in data:
                all_residues.update(data["distances"].keys())

        if not all_residues:
            print("No distance data available")
            return False

        # Select residues to plot
        residues_to_plot = key_residues if key_residues else sorted(list(all_residues))[:4]
        n_residues = len(residues_to_plot)

        # Create subplots
        fig, axes = plt.subplots(n_residues, 1, figsize=(12, 3 * n_residues))
        if n_residues == 1:
            axes = [axes]

        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

        for idx, res in enumerate(residues_to_plot):
            ax = axes[idx]

            for i, (data, name) in enumerate(zip(dynamics_data_list, trajectory_names)):
                if data and "distances" in data and res in data["distances"]:
                    time_ns = data["time"] / 1000.0
                    distances_ang = data["distances"][res] * 10
                    ax.plot(time_ns, distances_ang, color=colors[i % 3], label=name, linewidth=1.5, alpha=0.8)

            ax.set_ylabel(f"{res} Distance (Å)", fontfamily="Times New Roman", fontsize=14)
            ax.legend(frameon=True, fancybox=True, fontsize=12)
            ax.grid(True, alpha=0.3)

            if idx == n_residues - 1:
                ax.set_xlabel("Time (ns)", fontfamily="Times New Roman", fontsize=14)

        if title:
            plt.suptitle(title, fontsize=16, fontweight="bold")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        plt.close()

        return True

    except Exception as e:
        print(f"Error in key distances time series plotting: {e}")
        return False
