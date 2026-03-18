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
from ..core.publication_utils import apply_publication_style, PUBLICATION_FONTS, get_standard_figsize


def plot_rmsd_simple_timeseries(
    rmsd_results: Dict,
    output_path: str = None,
    title: str = "",
    ax: Optional[plt.Axes] = None,
    time_per_frame_ns: float = 0.5,
    ylabel: str = "Ligand RMSD (Å)",
    **kwargs,
):
    """
    Create simple RMSD time series plot matching the example format.

    Supports both standalone and panel modes for multi-panel figures.

    Parameters
    ----------
    rmsd_results : dict
        Dictionary with trajectory names as keys and RMSD data as values
    output_path : str, optional
        Path to save the plot. Required if ax is None (standalone mode).
    title : str
        Plot title (empty by default for publication style)
    ax : matplotlib.axes.Axes, optional
        If provided, plot on this axis (panel mode).
        If None, create new figure (standalone mode).
    time_per_frame_ns : float
        Time interval per frame in nanoseconds (default: 0.5)
    ylabel : str
        Y-axis label (default: "Ligand RMSD (Å)")
    **kwargs : dict
        Additional styling parameters (e.g., colors, linewidth)

    Returns
    -------
    bool or matplotlib.axes.Axes
        If standalone mode (ax=None): returns True/False for success
        If panel mode (ax provided): returns the axis object

    Examples
    --------
    # Standalone mode (backward compatible - original usage)
    >>> plot_rmsd_simple_timeseries(rmsd_data, "output.png")

    # Panel mode (new functionality for multi-panel figures)
    >>> fig, axes = plt.subplots(2, 2)
    >>> plot_rmsd_simple_timeseries(rmsd_data, ax=axes[0, 0])
    """
    try:
        print("📈 RMSD timeseries analysis showing structural stability over simulation time")

        # Determine mode: standalone (ax=None) or panel (ax provided)
        if ax is None:
            # Standalone mode - create new figure (original behavior)
            if output_path is None:
                raise ValueError("output_path is required when ax is None (standalone mode)")

            apply_publication_style()
            fig, ax = plt.subplots(figsize=get_standard_figsize("single"))
            own_figure = True
        else:
            # Panel mode - use provided axis
            own_figure = False

        # ========== Core plotting logic (shared by both modes) ==========

        # Use different colors for each trajectory
        colors = kwargs.get("colors", ["#4472C4", "#ED7D31", "#70AD47", "#FFC000", "#5B9BD5", "#C55A11"])
        linewidth = kwargs.get("linewidth", 1.5)
        alpha = kwargs.get("alpha", 0.85)

        # Plot all trajectories with different colors
        for i, (traj_name, rmsd_data) in enumerate(rmsd_results.items()):
            if isinstance(rmsd_data, dict) and "protein" in rmsd_data:
                protein_rmsd = rmsd_data["protein"]
            else:
                protein_rmsd = rmsd_data

            if len(protein_rmsd) > 0:
                time = np.arange(len(protein_rmsd)) * time_per_frame_ns
                color = colors[i % len(colors)]
                ax.plot(time, protein_rmsd, color=color, linewidth=linewidth, alpha=alpha, label=traj_name)

        # Format exactly like the example
        ax.set_xlabel("Time (ns)", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
        ax.tick_params(axis="both", labelsize=PUBLICATION_FONTS["tick_label"])

        # Remove grid and titles for clean look
        ax.grid(False)

        # Set clean axis limits
        all_rmsd = []
        for rmsd_data in rmsd_results.values():
            if isinstance(rmsd_data, dict) and "protein" in rmsd_data:
                all_rmsd.extend(rmsd_data["protein"])
            else:
                all_rmsd.extend(rmsd_data)

        if all_rmsd:
            ax.set_ylim(0, max(all_rmsd) * 1.05)

        # Add legend if multiple trajectories
        if len(rmsd_results) > 1:
            ax.legend(loc="best", fontsize=PUBLICATION_FONTS["legend"], framealpha=0.9)

        # Only add title if provided
        if title:
            ax.set_title(title, fontsize=PUBLICATION_FONTS["title"], fontweight="bold")

        # ========== Mode-specific post-processing ==========

        if own_figure:
            # Standalone mode: save and close figure
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
            plt.close()
            return True
        else:
            # Panel mode: return axis, caller handles saving
            return ax

    except Exception as e:
        print(f"Error in RMSD plotting: {e}")
        if own_figure:
            return False
        else:
            raise  # Panel mode: raise exception for caller to handle


def plot_rmsd_analysis(rmsd_results: Dict, output_path: str, title: str = "") -> bool:
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
        print("📈 RMSD 4-panel analysis: timeseries, distribution, averages, and ligand comparison")
        apply_publication_style()

        fig, axes = plt.subplots(2, 2, figsize=get_standard_figsize("quad"))
        if title:
            fig.suptitle(title, fontsize=16, fontweight="bold")
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

        # Plot 1: RMSD time series
        ax1 = axes[0, 0]
        for i, (traj_name, rmsd_data) in enumerate(rmsd_results.items()):
            if isinstance(rmsd_data, dict) and "protein" in rmsd_data:
                protein_rmsd = rmsd_data["protein"]
            else:
                protein_rmsd = rmsd_data

            if len(protein_rmsd) > 0:
                time = np.arange(len(protein_rmsd)) * 0.5  # Assume 0.5ns interval
                ax1.plot(time, protein_rmsd, color=colors[i % 3], label=traj_name, linewidth=2, alpha=0.8)

        ax1.set_xlabel("Time (ns)", fontfamily="Times New Roman")
        ax1.set_ylabel("RMSD (Å)", fontfamily="Times New Roman")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: RMSD distribution
        ax2 = axes[0, 1]
        all_rmsd = []
        for rmsd_data in rmsd_results.values():
            if isinstance(rmsd_data, dict) and "protein" in rmsd_data:
                protein_rmsd = rmsd_data["protein"]
            else:
                protein_rmsd = rmsd_data
            if len(protein_rmsd) > 0:
                all_rmsd.extend(protein_rmsd)

        if all_rmsd:
            ax2.hist(all_rmsd, bins=30, alpha=0.7, color="skyblue", edgecolor="black")
            ax2.set_xlabel("RMSD (Å)", fontfamily="Times New Roman")
            ax2.set_ylabel("Frequency", fontfamily="Times New Roman")
            ax2.grid(True, alpha=0.3)

        # Plot 3: Average RMSD per trajectory
        ax3 = axes[1, 0]
        traj_names = []
        avg_rmsd = []
        std_rmsd = []

        for traj_name, rmsd_data in rmsd_results.items():
            if isinstance(rmsd_data, dict) and "protein" in rmsd_data:
                protein_rmsd = rmsd_data["protein"]
            else:
                protein_rmsd = rmsd_data

            if len(protein_rmsd) > 0:
                traj_names.append(traj_name)
                avg_rmsd.append(np.mean(protein_rmsd))
                std_rmsd.append(np.std(protein_rmsd))

        if traj_names:
            bars = ax3.bar(traj_names, avg_rmsd, yerr=std_rmsd, color=colors[: len(traj_names)], alpha=0.7, capsize=5)
            ax3.set_ylabel("Average RMSD (Å)", fontfamily="Times New Roman")
            ax3.grid(True, alpha=0.3, axis="y")

        # Plot 4: Ligand RMSD (if available)
        ax4 = axes[1, 1]
        has_ligand = False

        for i, (traj_name, rmsd_data) in enumerate(rmsd_results.items()):
            if isinstance(rmsd_data, dict) and "ligand" in rmsd_data:
                ligand_rmsd = rmsd_data["ligand"]
                if len(ligand_rmsd) > 0:
                    time = np.arange(len(ligand_rmsd)) * 0.5
                    ax4.plot(time, ligand_rmsd, color=colors[i % 3], label=traj_name, linewidth=2, alpha=0.8)
                    has_ligand = True

        if has_ligand:
            ax4.set_xlabel("Time (ns)", fontfamily="Times New Roman")
            ax4.set_ylabel("Ligand RMSD (Å)", fontfamily="Times New Roman")
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        else:
            ax4.text(
                0.5,
                0.5,
                "No Ligand Data",
                ha="center",
                va="center",
                transform=ax4.transAxes,
                fontsize=16,
                fontfamily="Times New Roman",
            )
            ax4.set_xticks([])
            ax4.set_yticks([])

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        return True

    except Exception as e:
        print(f"Error in RMSD plotting: {e}")
        return False


def plot_rmsf_analysis(
    rmsf_results: Dict, output_path: str = None, title: str = "", ax: Optional[plt.Axes] = None, **kwargs
) -> bool:
    """
    Plot RMSF analysis results.

    Supports both standalone and panel modes. For panel mode with single axis,
    plots only the per-residue RMSF (not the distribution).

    Parameters
    ----------
    rmsf_results : dict
        Dictionary with trajectory names as keys and RMSF data as values
    output_path : str, optional
        Path to save the plot. Required if ax is None (standalone mode).
    title : str
        Plot title
    ax : matplotlib.axes.Axes, optional
        If provided, plot on this single axis (panel mode - per-residue only).
        If None, create 2-panel figure (standalone mode).
    **kwargs : dict
        Additional styling parameters

    Returns
    -------
    bool or matplotlib.axes.Axes
        If standalone mode (ax=None): returns True/False for success
        If panel mode (ax provided): returns the axis object

    Examples
    --------
    # Standalone mode (backward compatible - creates 2-panel figure)
    >>> plot_rmsf_analysis(rmsf_data, "output.png")

    # Panel mode (new - single panel with per-residue RMSF)
    >>> fig, axes = plt.subplots(2, 2)
    >>> plot_rmsf_analysis(rmsf_data, ax=axes[0, 0])
    """
    try:
        print("📈 RMSF analysis: per-residue fluctuations" + (" and distribution" if ax is None else ""))

        # Determine mode
        if ax is None:
            # Standalone mode - create 2-panel figure (original behavior)
            if output_path is None:
                raise ValueError("output_path is required when ax is None (standalone mode)")

            apply_publication_style()
            fig, axes = plt.subplots(1, 2, figsize=get_standard_figsize("horizontal"))
            ax1 = axes[0]
            ax2 = axes[1]
            own_figure = True
        else:
            # Panel mode - use single provided axis (per-residue RMSF only)
            ax1 = ax
            ax2 = None
            own_figure = False

        # ========== Core plotting logic ==========

        colors = kwargs.get("colors", ["#1f77b4", "#ff7f0e", "#2ca02c"])
        linewidth = kwargs.get("linewidth", 2)
        alpha = kwargs.get("alpha", 0.8)

        # Plot 1: RMSF per residue (always plotted)
        for i, (traj_name, rmsf_data) in enumerate(rmsf_results.items()):
            if isinstance(rmsf_data, tuple) and len(rmsf_data) >= 2:
                residue_indices, rmsf_values = rmsf_data[:2]
            elif isinstance(rmsf_data, dict) and "rmsf" in rmsf_data:
                rmsf_values = rmsf_data["rmsf"]
                residue_indices = np.arange(len(rmsf_values))
            else:
                continue

            if len(rmsf_values) > 0:
                ax1.plot(
                    residue_indices,
                    rmsf_values,
                    color=colors[i % len(colors)],
                    label=traj_name,
                    linewidth=linewidth,
                    alpha=alpha,
                )

        ax1.set_xlabel(
            "Residue Index", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold", fontfamily="Times New Roman"
        )
        ax1.set_ylabel(
            "RMSF (Å)", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold", fontfamily="Times New Roman"
        )
        ax1.legend(fontsize=PUBLICATION_FONTS["legend"])
        ax1.grid(True, alpha=0.3)
        ax1.tick_params(axis="both", labelsize=PUBLICATION_FONTS["tick_label"])

        # Plot 2: RMSF distribution (only in standalone mode)
        if ax2 is not None:
            all_rmsf = []
            for rmsf_data in rmsf_results.values():
                if isinstance(rmsf_data, tuple) and len(rmsf_data) >= 2:
                    rmsf_values = rmsf_data[1]
                elif isinstance(rmsf_data, dict) and "rmsf" in rmsf_data:
                    rmsf_values = rmsf_data["rmsf"]
                else:
                    continue

                if len(rmsf_values) > 0:
                    all_rmsf.extend(rmsf_values)

            if all_rmsf:
                ax2.hist(all_rmsf, bins=30, alpha=0.7, color="lightcoral", edgecolor="black")
                ax2.set_xlabel(
                    "RMSF (Å)",
                    fontsize=PUBLICATION_FONTS["axis_label"],
                    fontweight="bold",
                    fontfamily="Times New Roman",
                )
                ax2.set_ylabel(
                    "Frequency",
                    fontsize=PUBLICATION_FONTS["axis_label"],
                    fontweight="bold",
                    fontfamily="Times New Roman",
                )
                ax2.grid(True, alpha=0.3)
                ax2.tick_params(axis="both", labelsize=PUBLICATION_FONTS["tick_label"])

        # ========== Mode-specific post-processing ==========

        if own_figure:
            # Standalone mode: save and close
            if title:
                fig.suptitle(title, fontsize=PUBLICATION_FONTS["title"], fontweight="bold")
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            plt.close()
            return True
        else:
            # Panel mode: return axis
            if title:
                ax1.set_title(title, fontsize=PUBLICATION_FONTS["title"], fontweight="bold")
            return ax1

    except Exception as e:
        print(f"Error in RMSF plotting: {e}")
        if own_figure:
            return False
        else:
            raise


def plot_rmsf_with_auto_chains(
    rmsf_results: Dict[str, Tuple[np.ndarray, np.ndarray]],
    output_path: str,
    title: str = "",
    plot_type: str = "multi_panel",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    max_chains_per_row: int = 2,
) -> bool:
    """
    Create RMSF plots with automatic chain detection and flexible layout.

    Parameters
    ----------
    rmsf_results : dict
        Dictionary mapping chain names to (rmsf_values, residue_ids) tuples
        e.g., {"Chain PR1": (rmsf_array, residue_ids), ...}
    output_path : str
        Path to save the plot
    title : str
        Plot title (empty by default for publication style)
    plot_type : str
        Type of plot layout: "multi_panel", "single_panel", "separate_files"
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
    save_path : str, optional
        Deprecated: use output_path instead
    max_chains_per_row : int
        Maximum number of chains per row in multi-panel layout

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        apply_publication_style()

        if save_path is not None:
            output_path = save_path

        if not rmsf_results:
            print("No RMSF results provided for plotting")
            return False

        n_chains = len(rmsf_results)
        if n_chains == 0:
            print("No chain data available for RMSF plotting")
            return False

        # Sort chains by name for consistent ordering
        sorted_chains = sorted(rmsf_results.items(), key=lambda x: x[0])
        chain_names = [name for name, _ in sorted_chains]

        print(f"Creating RMSF plot for {n_chains} chains: {chain_names}")

        if plot_type == "multi_panel":
            return _plot_rmsf_multi_panel(sorted_chains, output_path, title, figsize, max_chains_per_row)
        elif plot_type == "single_panel":
            return _plot_rmsf_single_panel(sorted_chains, output_path, title, figsize)
        elif plot_type == "separate_files":
            return _plot_rmsf_separate_files(sorted_chains, output_path, title, figsize)
        else:
            print(f"Unknown plot_type: {plot_type}. Using multi_panel.")
            return _plot_rmsf_multi_panel(sorted_chains, output_path, title, figsize, max_chains_per_row)

    except Exception as e:
        print(f"Error in RMSF plotting with auto chains: {e}")
        return False


def _plot_rmsf_multi_panel(
    sorted_chains: List[Tuple[str, Tuple[np.ndarray, np.ndarray]]],
    output_path: str,
    title: str,
    figsize: Optional[Tuple[float, float]],
    max_chains_per_row: int,
) -> bool:
    """Create multi-panel RMSF plot with automatic chain type detection."""
    try:
        n_chains = len(sorted_chains)
        rows = (n_chains + max_chains_per_row - 1) // max_chains_per_row
        cols = min(n_chains, max_chains_per_row)

        # Auto-calculate figure size for publication quality
        if figsize is None:
            # Make panels similar in size to the original figure
            panel_width = 6  # Each panel width in inches
            panel_height = 4  # Each panel height in inches
            figsize = (panel_width * cols, panel_height * rows)

        fig, axes = plt.subplots(rows, cols, figsize=figsize)
        if rows == 1 and cols == 1:
            axes = [axes]
        elif rows == 1:
            axes = list(axes)
        else:
            axes = axes.flatten()

        # Enhanced chain type detection and colors
        protein_chains = [
            (name, data)
            for name, data in sorted_chains
            if "Protein" in name or ("Chain" in name and "Nucleic" not in name)
        ]
        nucleic_chains = [
            (name, data) for name, data in sorted_chains if "Nucleic" in name or "RNA" in name or "DNA" in name
        ]

        # Standard colors for publication
        protein_color = "#4472C4"  # Professional blue
        nucleic_color = "#E15759"  # Professional red/coral

        for i, (chain_name, (rmsf_values, residue_ids)) in enumerate(sorted_chains):
            ax = axes[i]

            # Determine chain type and appropriate labeling
            is_nucleic = any(keyword in chain_name.lower() for keyword in ["nucleic", "rna", "dna"])

            if is_nucleic:
                color = nucleic_color
                x_label = "Residue Number"
                # Clean up chain name for display
                display_name = chain_name.replace("Nucleic Chain", "Nucleic Chain").replace("Chain", "Chain")
                print(f"  📊 Plotting nucleic acid chain: {display_name} ({len(rmsf_values)} P atoms)")
            else:
                color = protein_color
                x_label = "Residue Number"
                # Clean up chain name for display
                display_name = chain_name.replace("Protein Chain", "Chain").replace("Chain", "Chain")
                print(f"  📊 Plotting protein chain: {display_name} ({len(rmsf_values)} CA atoms)")

            # Plot RMSF with filled area
            ax.plot(residue_ids, rmsf_values, color=color, linewidth=2, alpha=0.9)
            ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color=color)

            # Calculate and print statistics (removed from plot to avoid hiding content)
            mean_rmsf = np.mean(rmsf_values)
            std_rmsf = np.std(rmsf_values)
            print(f"  📊 {display_name} RMSF: Mean = {mean_rmsf:.1f} ± {std_rmsf:.2f} Å")

            # Add panel title (matching the original figure style)
            panel_title = display_name.replace("Protein ", "").replace("Nucleic ", "Nucleic ")
            ax.set_title(panel_title, fontsize=PUBLICATION_FONTS["title"], fontweight="bold", pad=10)

            # Formatting with appropriate x-axis label
            ax.set_xlabel(x_label, fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
            ax.set_ylabel("RMSF (Å)", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis="both", labelsize=PUBLICATION_FONTS["tick_label"])

            # Set y-axis to start from 0 for better comparison
            ax.set_ylim(bottom=0)

        # Hide unused subplots
        for i in range(n_chains, len(axes)):
            axes[i].set_visible(False)

        # Add main title if provided (but keep minimal per publication style)
        if title:
            fig.suptitle(title, fontsize=PUBLICATION_FONTS["title"] + 2, fontweight="bold", y=0.95)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        plt.close()

        print(f"✅ Multi-panel RMSF plot saved: {output_path}")
        return True

    except Exception as e:
        print(f"Error creating multi-panel RMSF plot: {e}")
        return False


def _plot_rmsf_single_panel(
    sorted_chains: List[Tuple[str, Tuple[np.ndarray, np.ndarray]]],
    output_path: str,
    title: str,
    figsize: Optional[Tuple[float, float]],
) -> bool:
    """Create single-panel RMSF plot with all chains."""
    try:
        # Auto-calculate figure size if not provided
        if figsize is None:
            figsize = (12, 8)

        fig, ax = plt.subplots(figsize=figsize)

        # Colors for different chain types
        n_chains = len(sorted_chains)
        colors = plt.cm.tab10(np.linspace(0, 1, n_chains))

        for i, (chain_name, (rmsf_values, residue_ids)) in enumerate(sorted_chains):
            color = colors[i]

            # Plot RMSF
            ax.plot(residue_ids, rmsf_values, color=color, linewidth=2, alpha=0.8, label=chain_name)

        # Formatting
        ax.set_xlabel("Residue Number", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
        ax.set_ylabel("RMSF (Å)", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
        if title:
            ax.set_title(title, fontsize=PUBLICATION_FONTS["title"], fontweight="bold")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=PUBLICATION_FONTS["legend"])
        ax.tick_params(axis="both", labelsize=PUBLICATION_FONTS["tick_label"])

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        plt.close()

        print(f"✅ Single-panel RMSF plot saved: {output_path}")
        return True

    except Exception as e:
        print(f"Error creating single-panel RMSF plot: {e}")
        return False


def _plot_rmsf_separate_files(
    sorted_chains: List[Tuple[str, Tuple[np.ndarray, np.ndarray]]],
    output_path: str,
    title: str,
    figsize: Optional[Tuple[float, float]],
) -> bool:
    """Create separate RMSF plots for each chain."""
    try:
        output_base = Path(output_path)
        success_count = 0

        # Auto-calculate figure size if not provided
        if figsize is None:
            figsize = (10, 6)

        for i, (chain_name, (rmsf_values, residue_ids)) in enumerate(sorted_chains):
            try:
                fig, ax = plt.subplots(figsize=figsize)

                # Determine color based on chain type
                if "Nucleic" in chain_name:
                    color = "red"
                else:
                    color = "blue"

                # Plot RMSF
                ax.plot(residue_ids, rmsf_values, color=color, linewidth=2, alpha=0.8)
                ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color=color)

                # Calculate statistics
                mean_rmsf = np.mean(rmsf_values)
                std_rmsf = np.std(rmsf_values)
                # Print statistics (removed from plot to avoid hiding content)
                print(f"  📊 {chain_name} RMSF: Mean = {mean_rmsf:.2f} ± {std_rmsf:.2f} Å")

                # Formatting
                ax.set_xlabel("Residue Number", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
                ax.set_ylabel("RMSF (Å)", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
                if title:
                    ax.set_title(f"{title} - {chain_name}", fontsize=PUBLICATION_FONTS["title"], fontweight="bold")
                ax.grid(True, alpha=0.3)
                ax.tick_params(axis="both", labelsize=PUBLICATION_FONTS["tick_label"])

                # Save individual file
                chain_output = (
                    output_base.parent
                    / f"{output_base.stem}_{chain_name.lower().replace(' ', '_')}{output_base.suffix}"
                )
                plt.tight_layout()
                plt.savefig(chain_output, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
                plt.close()

                print(f"✅ Individual RMSF plot saved: {chain_output}")
                success_count += 1

            except Exception as e:
                print(f"Error creating plot for {chain_name}: {e}")
                continue

        print(f"✅ Created {success_count}/{len(sorted_chains)} individual RMSF plots")
        return success_count > 0

    except Exception as e:
        print(f"Error creating separate RMSF files: {e}")
        return False


def plot_multi_trajectory_rmsf_comprehensive(
    multi_traj_data: Dict[str, Dict],
    output_path: str,
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    plot_individual_trajectories: bool = True,
    save_data: bool = True,
) -> bool:
    """
    Create comprehensive multi-trajectory RMSF plot showing all chains and trajectories.

    Parameters
    ----------
    multi_traj_data : Dict[str, Dict]
        Multi-trajectory RMSF data from calculate_multi_trajectory_rmsf()
    output_path : str
        Path to save the plot
    title : str, optional
        Plot title (empty by default for publication style)
    figsize : Tuple[float, float], optional
        Figure size (width, height)
    plot_individual_trajectories : bool
        Whether to show individual trajectory lines
    save_data : bool
        Whether to save data as CSV files

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        apply_publication_style()

        if not multi_traj_data:
            print("No multi-trajectory RMSF data provided for plotting")
            return False

        n_chains = len(multi_traj_data)
        if n_chains == 0:
            print("No chain data available for multi-trajectory RMSF plotting")
            return False

        # Auto-calculate figure size for publication quality
        if figsize is None:
            if n_chains == 1:
                figsize = (10, 6)
            elif n_chains == 2:
                figsize = (12, 5)
            else:
                figsize = (15, 8)

        # Create subplots
        if n_chains == 1:
            fig, axes = plt.subplots(1, 1, figsize=figsize)
            axes = [axes]
        elif n_chains == 2:
            fig, axes = plt.subplots(1, 2, figsize=figsize)
        else:
            # For more chains, use a grid layout
            rows = (n_chains + 1) // 2
            cols = 2
            fig, axes = plt.subplots(rows, cols, figsize=figsize)
            axes = axes.flatten()

        # Colors for different molecule types
        protein_color = "#4472C4"  # Professional blue
        nucleic_color = "#E15759"  # Professional red/coral
        individual_alpha = 0.3
        mean_alpha = 0.9

        print(f"📊 Creating comprehensive multi-trajectory RMSF plot for {n_chains} chains")

        chain_names = sorted(multi_traj_data.keys())

        for i, chain_name in enumerate(chain_names):
            ax = axes[i]
            chain_data = multi_traj_data[chain_name]

            # Determine chain type and colors
            is_nucleic = any(keyword in chain_name.lower() for keyword in ["nucleic", "rna", "dna"])

            if is_nucleic:
                main_color = nucleic_color
                x_label = "Residue Number"
                print(f"  📊 Plotting nucleic acid chain: {chain_name}")
            else:
                main_color = protein_color
                x_label = "Residue Number"
                print(f"  📊 Plotting protein chain: {chain_name}")

            # Get combined data
            combined_data = chain_data["combined"]
            mean_rmsf = combined_data["mean_rmsf"]
            std_rmsf = combined_data["std_rmsf"]
            residue_ids = combined_data["residue_ids"]

            # Plot individual trajectories if requested
            if plot_individual_trajectories and "trajectories" in chain_data:
                traj_names = sorted(chain_data["trajectories"].keys())
                for j, traj_name in enumerate(traj_names):
                    traj_rmsf, traj_residues = chain_data["trajectories"][traj_name]

                    # Ensure residue IDs match (they should if everything is working correctly)
                    if len(traj_rmsf) == len(residue_ids):
                        ax.plot(
                            residue_ids,
                            traj_rmsf,
                            color=main_color,
                            alpha=individual_alpha,
                            linewidth=1,
                            label=f"Trajectory {j+1}" if i == 0 else None,
                        )

            # Plot mean with error bars/band
            ax.plot(
                residue_ids,
                mean_rmsf,
                color=main_color,
                linewidth=3,
                alpha=mean_alpha,
                label="Mean" if i == 0 else None,
            )

            # Add standard deviation as filled area
            ax.fill_between(
                residue_ids,
                mean_rmsf - std_rmsf,
                mean_rmsf + std_rmsf,
                alpha=0.3,
                color=main_color,
                label="±1 SD" if i == 0 else None,
            )

            # Add statistics
            stats = chain_data["statistics"]
            overall_mean = stats["overall_mean"]
            overall_std = stats["overall_std"]
            n_trajectories = stats["n_trajectories"]

            # Print statistics (removed from plot to avoid hiding content)
            print(
                f"  📊 {chain_name} RMSF: Mean = {overall_mean:.1f} ± {overall_std:.2f} Å ({n_trajectories} trajectories)"
            )

            # Panel title
            clean_chain_name = chain_name.replace("Protein ", "").replace("Nucleic ", "Nucleic ")
            ax.set_title(clean_chain_name, fontsize=PUBLICATION_FONTS["title"], fontweight="bold", pad=10)

            # Formatting
            ax.set_xlabel(x_label, fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
            ax.set_ylabel("RMSF (Å)", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis="both", labelsize=PUBLICATION_FONTS["tick_label"])
            ax.set_ylim(bottom=0)

            # Add legend only to first subplot
            if i == 0 and plot_individual_trajectories:
                ax.legend(loc="upper right", fontsize=10, framealpha=0.9)

        # Hide unused subplots if any
        for i in range(n_chains, len(axes)):
            axes[i].set_visible(False)

        # Add main title if provided
        if title:
            fig.suptitle(title, fontsize=PUBLICATION_FONTS["title"] + 2, fontweight="bold", y=0.95)

        plt.tight_layout()

        # Save plot
        plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        plt.close()

        print(f"✅ Comprehensive multi-trajectory RMSF plot saved: {output_path}")

        # Save data to CSV if requested
        if save_data:
            output_path_obj = Path(output_path)
            csv_path = output_path_obj.with_suffix(".csv")

            # Create comprehensive data table
            all_data = []
            for chain_name, chain_data in multi_traj_data.items():
                combined_data = chain_data["combined"]
                stats = chain_data["statistics"]

                for i, (rmsf_val, res_id) in enumerate(zip(combined_data["mean_rmsf"], combined_data["residue_ids"])):
                    all_data.append(
                        {
                            "Chain": chain_name,
                            "Residue_ID": res_id,
                            "Mean_RMSF": rmsf_val,
                            "Std_RMSF": combined_data["std_rmsf"][i],
                            "Overall_Mean": stats["overall_mean"],
                            "Overall_Std": stats["overall_std"],
                            "N_Trajectories": stats["n_trajectories"],
                        }
                    )

            import pandas as pd

            df = pd.DataFrame(all_data)
            df.to_csv(csv_path, index=False)
            print(f"✅ RMSF data saved to CSV: {csv_path}")

        return True

    except Exception as e:
        print(f"Error creating comprehensive multi-trajectory RMSF plot: {e}")
        import traceback

        traceback.print_exc()
        return False


def plot_rmsd_time_series(
    rmsd_data: Dict[str, np.ndarray],
    times: Optional[np.ndarray] = None,
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    colors: Optional[List[str]] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot RMSD time series for multiple selections.

    Parameters
    ----------
    rmsd_data : dict
        Dictionary with selection names as keys and RMSD arrays as values
    times : np.ndarray, optional
        Time points. If None, use frame indices.
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    colors : list, optional
        Colors for each trace

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("single")
    fig, ax = plt.subplots(figsize=figsize)

    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(rmsd_data)))

    for i, (label, values) in enumerate(rmsd_data.items()):
        if times is None:
            plot_times = np.arange(len(values)) * 0.5  # Assume 0.5 ns per frame
        else:
            plot_times = times[: len(values)]

        ax.plot(plot_times, values, color=colors[i], label=label, linewidth=2, alpha=0.8)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (nm)")
    if title:
        ax.set_title(title)

    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"RMSD time series plot saved: {save_path}")

    return fig, ax


def plot_rmsf_per_residue(
    rmsf_values: np.ndarray,
    residue_ids: Optional[np.ndarray] = None,
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    highlight_threshold: Optional[float] = None,
    secondary_structure: Optional[Dict] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot RMSF values per residue.

    Parameters
    ----------
    rmsf_values : np.ndarray
        RMSF values for each residue
    residue_ids : np.ndarray, optional
        Residue IDs. If None, use sequential numbering.
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    highlight_threshold : float, optional
        Threshold for highlighting flexible regions
    secondary_structure : dict, optional
        Secondary structure annotations

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("single")
    fig, ax = plt.subplots(figsize=figsize)

    if residue_ids is None:
        residue_ids = np.arange(1, len(rmsf_values) + 1)

    # Main RMSF plot
    ax.plot(residue_ids, rmsf_values, "b-", linewidth=2, alpha=0.8)
    ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color="blue")

    # Highlight flexible regions if threshold provided
    if highlight_threshold is not None:
        flexible_mask = rmsf_values > highlight_threshold
        if np.any(flexible_mask):
            ax.fill_between(
                residue_ids,
                rmsf_values,
                highlight_threshold,
                where=flexible_mask,
                alpha=0.5,
                color="red",
                label=f"Flexible (RMSF > {highlight_threshold:.2f} nm)",
            )

    # Add secondary structure annotations if provided
    if secondary_structure:
        y_max = np.max(rmsf_values) * 1.1
        for ss_type, regions in secondary_structure.items():
            color = {"helix": "red", "sheet": "green", "loop": "gray"}.get(ss_type, "black")
            for start, end in regions:
                ax.axvspan(start, end, alpha=0.2, color=color, label=f"{ss_type.capitalize()}")

    ax.set_xlabel("Residue Number")
    ax.set_ylabel("RMSF (nm)")
    if title:
        ax.set_title(title)

    ax.grid(True, alpha=0.3)

    # Print statistics (removed from plot to avoid hiding content)
    mean_rmsf = np.mean(rmsf_values)
    std_rmsf = np.std(rmsf_values)
    print(f"  📊 RMSF statistics: Mean = {mean_rmsf:.3f} ± {std_rmsf:.3f} nm")

    if highlight_threshold is not None and np.any(flexible_mask):
        ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"RMSF plot saved: {save_path}")

    return fig, ax


def plot_rmsd_rmsf_combined(
    rmsd_data: Dict[str, np.ndarray],
    rmsf_data: Dict[str, np.ndarray],
    times: Optional[np.ndarray] = None,
    residue_ids: Optional[np.ndarray] = None,
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Create combined RMSD time series and RMSF plots.

    Parameters
    ----------
    rmsd_data : dict
        Dictionary with selection names as keys and RMSD arrays as values
    rmsf_data : dict
        Dictionary with selection names as keys and RMSF arrays as values
    times : np.ndarray, optional
        Time points for RMSD plot
    residue_ids : np.ndarray, optional
        Residue IDs for RMSF plot
    title : str
        Overall plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, (rmsd_ax, rmsf_ax)) objects
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    if title:
        fig.suptitle(title, fontsize=16, fontweight="bold")

    # RMSD subplot
    colors = plt.cm.tab10(np.linspace(0, 1, len(rmsd_data)))
    for i, (label, values) in enumerate(rmsd_data.items()):
        if times is None:
            plot_times = np.arange(len(values)) * 0.5
        else:
            plot_times = times[: len(values)]

        ax1.plot(plot_times, values, color=colors[i], label=label, linewidth=2, alpha=0.8)

    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel("RMSD (nm)")

    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # RMSF subplot
    for i, (label, values) in enumerate(rmsf_data.items()):
        if residue_ids is None:
            plot_residues = np.arange(1, len(values) + 1)
        else:
            plot_residues = residue_ids[: len(values)]

        ax2.plot(plot_residues, values, color=colors[i], label=label, linewidth=2, alpha=0.8)

    ax2.set_xlabel("Residue Number")
    ax2.set_ylabel("RMSF (nm)")

    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Combined RMSD/RMSF plot saved: {save_path}")


def plot_multi_chain_rmsf(
    chain_rmsf_data: Dict[str, Tuple[np.ndarray, np.ndarray]],
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    cols: int = 2,
) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Create multi-chain RMSF subplots with automatic layout.

    Parameters
    ----------
    chain_rmsf_data : dict
        Dictionary mapping chain names to (rmsf_values, residue_ids) tuples
    title : str
        Overall plot title
    figsize : tuple, optional
        Figure size. If None, calculated automatically.
    save_path : str, optional
        Path to save figure
    cols : int
        Number of columns in subplot grid

    Returns
    -------
    tuple
        (figure, axes_list) objects
    """
    n_chains = len(chain_rmsf_data)
    if n_chains == 0:
        raise ValueError("No chain data provided")

    # Calculate subplot layout
    rows = (n_chains + cols - 1) // cols

    # Auto-calculate figure size if not provided
    if figsize is None:
        figsize = (6 * cols, 4 * rows)

    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    if n_chains == 1:
        axes = [axes]
    elif rows == 1:
        axes = list(axes) if hasattr(axes, "__iter__") else [axes]
    else:
        axes = axes.flatten()

    # Colors for different chain types
    protein_colors = plt.cm.Blues(
        np.linspace(0.4, 0.9, sum(1 for name in chain_rmsf_data.keys() if "Protein" in name or "Chain" in name))
    )
    nucleic_colors = plt.cm.Reds(np.linspace(0.4, 0.9, sum(1 for name in chain_rmsf_data.keys() if "Nucleic" in name)))

    protein_idx, nucleic_idx = 0, 0

    for i, (chain_name, (rmsf_values, residue_ids)) in enumerate(chain_rmsf_data.items()):
        ax = axes[i]

        # Choose color based on chain type
        if "Nucleic" in chain_name:
            color = nucleic_colors[nucleic_idx] if nucleic_idx < len(nucleic_colors) else "red"
            nucleic_idx += 1
        else:
            color = protein_colors[protein_idx] if protein_idx < len(protein_colors) else "blue"
            protein_idx += 1

        # Plot RMSF
        ax.plot(residue_ids, rmsf_values, color=color, linewidth=2, alpha=0.8)
        ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color=color)

        # Formatting
        ax.set_xlabel("Residue Number")
        ax.set_ylabel("RMSF (Å)")

        ax.grid(True, alpha=0.3)

        # Print statistics (removed from plot to avoid hiding content)
        mean_rmsf = np.mean(rmsf_values)
        std_rmsf = np.std(rmsf_values)
        print(f"  📊 {chain_name} RMSF: Mean = {mean_rmsf:.2f} ± {std_rmsf:.2f} Å")

    # Hide unused subplots
    for i in range(n_chains, len(axes)):
        axes[i].set_visible(False)

    if title:
        fig.suptitle(title, fontsize=16, fontweight="bold")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-chain RMSF plot saved: {save_path}")

    return fig, axes


def plot_separate_rmsd(
    protein_rmsd: np.ndarray,
    ligand_rmsd: np.ndarray,
    times: Optional[np.ndarray] = None,
    protein_title: str = "",
    ligand_title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    protein_save_path: Optional[str] = None,
    ligand_save_path: Optional[str] = None,
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Create separate RMSD plots for protein and ligand.

    Parameters
    ----------
    protein_rmsd : np.ndarray
        Protein RMSD values
    ligand_rmsd : np.ndarray
        Ligand RMSD values
    times : np.ndarray, optional
        Time points. If None, use frame indices.
    protein_title : str
        Title for protein RMSD plot
    ligand_title : str
        Title for ligand RMSD plot
    figsize : tuple
        Figure size
    protein_save_path : str, optional
        Path to save protein RMSD plot
    ligand_save_path : str, optional
        Path to save ligand RMSD plot

    Returns
    -------
    tuple
        (figure, (protein_ax, ligand_ax)) objects
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)

    if times is None:
        times = np.arange(len(protein_rmsd))  # Frame numbers (no unit conversion needed)

    # Protein RMSD plot
    ax1.plot(times[: len(protein_rmsd)], protein_rmsd, color="blue", linewidth=2, alpha=0.8, label="Protein CA")
    ax1.fill_between(times[: len(protein_rmsd)], protein_rmsd, alpha=0.3, color="blue")
    ax1.set_xlabel("Frame")
    ax1.set_ylabel("RMSD (Å)")
    if protein_title:
        ax1.set_title(protein_title)

    ax1.grid(True, alpha=0.3)

    # Add protein statistics
    protein_mean = np.mean(protein_rmsd)
    protein_std = np.std(protein_rmsd)
    protein_stats = f"Mean: {protein_mean:.3f} ± {protein_std:.3f} Å"
    ax1.text(
        0.02,
        0.98,
        protein_stats,
        transform=ax1.transAxes,
        verticalalignment="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.8),
    )

    # Ligand RMSD plot
    ax2.plot(times[: len(ligand_rmsd)], ligand_rmsd, color="red", linewidth=2, alpha=0.8, label="Ligand")
    ax2.fill_between(times[: len(ligand_rmsd)], ligand_rmsd, alpha=0.3, color="red")
    ax2.set_xlabel("Frame")
    ax2.set_ylabel("RMSD (Å)")
    if ligand_title:
        ax2.set_title(ligand_title)

    ax2.grid(True, alpha=0.3)

    # Add ligand statistics
    ligand_mean = np.mean(ligand_rmsd)
    ligand_std = np.std(ligand_rmsd)
    ligand_stats = f"Mean: {ligand_mean:.2f} ± {ligand_std:.2f} Å"
    ax2.text(
        0.02,
        0.98,
        ligand_stats,
        transform=ax2.transAxes,
        verticalalignment="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="lightcoral", alpha=0.8),
    )

    plt.tight_layout()

    # Save individual plots if paths provided
    if protein_save_path:
        # Create individual protein plot
        fig_protein, ax_protein = plt.subplots(figsize=(10, 5))
        ax_protein.plot(times[: len(protein_rmsd)], protein_rmsd, color="blue", linewidth=2, alpha=0.8)
        ax_protein.fill_between(times[: len(protein_rmsd)], protein_rmsd, alpha=0.3, color="blue")
        ax_protein.set_xlabel("Frame")
        ax_protein.set_ylabel("RMSD (Å)")
        if protein_title:
            ax_protein.set_title(protein_title)

        ax_protein.grid(True, alpha=0.3)
        ax_protein.text(
            0.02,
            0.98,
            protein_stats,
            transform=ax_protein.transAxes,
            verticalalignment="top",
            fontsize=10,
            bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.8),
        )
        plt.tight_layout()
        fig_protein.savefig(protein_save_path, dpi=300, bbox_inches="tight")
        plt.close(fig_protein)
        logger.info(f"Protein RMSD plot saved: {protein_save_path}")

    if ligand_save_path:
        # Create individual ligand plot
        fig_ligand, ax_ligand = plt.subplots(figsize=(10, 5))
        ax_ligand.plot(times[: len(ligand_rmsd)], ligand_rmsd, color="red", linewidth=2, alpha=0.8)
        ax_ligand.fill_between(times[: len(ligand_rmsd)], ligand_rmsd, alpha=0.3, color="red")
        ax_ligand.set_xlabel("Frame")
        ax_ligand.set_ylabel("RMSD (Å)")
        if ligand_title:
            ax_ligand.set_title(ligand_title)

        ax_ligand.grid(True, alpha=0.3)
        ax_ligand.text(
            0.02,
            0.98,
            ligand_stats,
            transform=ax_ligand.transAxes,
            verticalalignment="top",
            fontsize=10,
            bbox=dict(boxstyle="round", facecolor="lightcoral", alpha=0.8),
        )
        plt.tight_layout()
        fig_ligand.savefig(ligand_save_path, dpi=300, bbox_inches="tight")
        plt.close(fig_ligand)
        logger.info(f"Ligand RMSD plot saved: {ligand_save_path}")

    return fig, (ax1, ax2)


def plot_multi_repeat_ligand_rmsd(
    rmsd_data: Dict[str, np.ndarray],
    times: Optional[np.ndarray] = None,
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    colors: Optional[List[str]] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot ligand RMSD for multiple repeat experiments.

    Parameters
    ----------
    rmsd_data : dict
        Dictionary mapping repeat names (e.g., 'Repeat 1', 'Repeat 2') to RMSD arrays
    times : np.ndarray, optional
        Time points. If None, use frame indices.
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    colors : list, optional
        Colors for each repeat

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("single")
    fig, ax = plt.subplots(figsize=figsize)

    if colors is None:
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    for i, (repeat_name, rmsd_values) in enumerate(rmsd_data.items()):
        if times is None:
            plot_times = np.arange(len(rmsd_values))  # Frame numbers (no unit conversion needed)
        else:
            plot_times = times[: len(rmsd_values)]

        color = colors[i % len(colors)]

        ax.plot(plot_times, rmsd_values, color=color, linestyle="-", linewidth=2.5, alpha=0.8, label=repeat_name)

    ax.set_xlabel("Frame")
    ax.set_ylabel("Ligand RMSD (Å)")
    if title:
        ax.set_title(title)

    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    # Add overall statistics
    all_values = np.concatenate([values for values in rmsd_data.values()])
    overall_mean = np.mean(all_values)
    overall_std = np.std(all_values)
    stats_text = f"Overall: {overall_mean:.2f} ± {overall_std:.2f} Å"

    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-repeat ligand RMSD plot saved: {save_path}")

    return fig, ax


def plot_multi_chain_rmsf_example_style(
    chain_data: Dict[str, Dict[str, np.ndarray]], output_path: str, title: str = "", separate_panels: bool = False
) -> bool:
    """
    Create 4-panel RMSF plot exactly matching the example multi-chain layout.

    Expected format: chain_data = {
        'PR1': {'rmsf': array, 'residue_ids': array, 'n_residues': int},
        'PR2': {'rmsf': array, 'residue_ids': array, 'n_residues': int},
        'PR3': {'rmsf': array, 'residue_ids': array, 'n_residues': int},
        'PR4': {'rmsf': array, 'residue_ids': array, 'n_residues': int}
    }

    Parameters
    ----------
    chain_data : dict
        Dictionary with chain names as keys and RMSF data as values
    output_path : str
        Path to save the figure
    title : str
        Main plot title (empty by default for publication style)
    separate_panels : bool
        If True, also save individual panel figures

    Returns
    -------
    bool
        True if successful
    """
    try:
        apply_publication_style()

        # Create 2x2 subplot layout matching example
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Define expected chains (PR1, PR2, PR3, PR4)
        chain_names = ["PR1", "PR2", "PR3", "PR4"]

        # Check if we have data
        if not chain_data:
            print("No chain data available for RMSF plotting")
            return False

        for i, chain_name in enumerate(chain_names):
            if i >= 4:  # Max 4 panels
                break

            ax = axes[i // 2, i % 2]

            if chain_name in chain_data:
                data = chain_data[chain_name]
                rmsf_values = data["rmsf"]
                residue_ids = data.get("residue_ids", np.arange(len(rmsf_values)))
                n_residues = len(rmsf_values)

                # Calculate and print statistics (removed from plot to avoid hiding content)
                mean_rmsf = np.mean(rmsf_values)
                std_rmsf = np.std(rmsf_values)
                print(f"  📊 {chain_name} RMSF: Mean = {mean_rmsf:.2f} ± {std_rmsf:.2f} Å")

                # Create filled area plot (exact match to example)
                ax.fill_between(residue_ids, 0, rmsf_values, alpha=0.6, color="lightblue")
                ax.plot(residue_ids, rmsf_values, color="steelblue", linewidth=1)

                # Panel title (shorter)

                # Set appropriate axis limits
                ax.set_ylim(0, max(rmsf_values) * 1.1)
            else:
                # Empty panel if no data
                ax.text(
                    0.5,
                    0.5,
                    f"No data for Chain {chain_name}",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=PUBLICATION_FONTS["axis_label"],
                    style="italic",
                )

            # Common formatting for all panels
            ax.set_xlabel("Residue Number", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
            ax.set_ylabel("RMSF (Å)", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis="both", labelsize=PUBLICATION_FONTS["tick_label"])

        # Save individual panels if requested
        if separate_panels:
            base_path = Path(output_path)
            for i, chain_name in enumerate(chain_names):
                if chain_name in chain_data:
                    separate_path = base_path.parent / f"{base_path.stem}_{chain_name.lower()}{base_path.suffix}"

                    # Create individual figure
                    fig_single, ax_single = plt.subplots(figsize=(8, 6))
                    data = chain_data[chain_name]
                    rmsf_values = data["rmsf"]
                    residue_ids = data.get("residue_ids", np.arange(len(rmsf_values)))
                    n_residues = len(rmsf_values)
                    mean_rmsf = np.mean(rmsf_values)
                    std_rmsf = np.std(rmsf_values)
                    print(f"  📊 {chain_name} (separate) RMSF: Mean = {mean_rmsf:.2f} ± {std_rmsf:.2f} Å")

                    ax_single.fill_between(residue_ids, 0, rmsf_values, alpha=0.6, color="lightblue")
                    ax_single.plot(residue_ids, rmsf_values, color="steelblue", linewidth=1)

                    ax_single.set_xlabel("Residue Number", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
                    ax_single.set_ylabel("RMSF (Å)", fontsize=PUBLICATION_FONTS["axis_label"], fontweight="bold")
                    ax_single.grid(True, alpha=0.3)
                    ax_single.set_ylim(0, max(rmsf_values) * 1.1)

                    plt.tight_layout()
                    fig_single.savefig(separate_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
                    plt.close(fig_single)

        # Overall title if provided (but empty by default)
        if title:
            fig.suptitle(title, fontsize=PUBLICATION_FONTS["title"], fontweight="bold", y=0.95)

        plt.tight_layout()
        if title:
            plt.subplots_adjust(top=0.90)  # Make room for title

        # Save main figure
        fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        plt.close(fig)

        return True

    except Exception as e:
        print(f"Error in multi-chain RMSF plotting: {e}")
        return False
