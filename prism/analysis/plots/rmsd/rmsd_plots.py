#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD plotting utilities for PRISM analysis module.
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Import publication style utilities
from ..core.publication_utils import (
    apply_publication_style,
    PUBLICATION_FONTS,
    get_standard_figsize,
    get_color_palette,
    setup_publication_figure,
    setup_publication_subplots,
    save_publication_figure,
    style_axes_for_publication,
    set_publication_labels,
    add_corner_annotation,
)

# Backward-compatible re-exports of RMSF functions (moved to rmsf_plots.py)
from . import rmsf_plots as _rmsf_plots

plot_rmsf_analysis = _rmsf_plots.plot_rmsf_analysis
plot_rmsf_with_auto_chains = _rmsf_plots.plot_rmsf_with_auto_chains
plot_multi_trajectory_rmsf_comprehensive = _rmsf_plots.plot_multi_trajectory_rmsf_comprehensive
plot_rmsf_per_residue = _rmsf_plots.plot_rmsf_per_residue
plot_rmsd_rmsf_combined = _rmsf_plots.plot_rmsd_rmsf_combined
plot_multi_chain_rmsf = _rmsf_plots.plot_multi_chain_rmsf
plot_multi_chain_rmsf_example_style = _rmsf_plots.plot_multi_chain_rmsf_example_style


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

            fig, ax = setup_publication_figure(panel_type="single")
            own_figure = True
        else:
            # Panel mode - use provided axis
            own_figure = False

        # ========== Core plotting logic (shared by both modes) ==========

        # Use color palette from publication_utils
        palette_name = kwargs.get("palette", "default")
        n_colors = len(rmsd_results)
        colors = kwargs.get("colors", get_color_palette(palette_name, n_colors))
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
            # Standalone mode: save and close figure using publication utils
            save_publication_figure(fig, output_path, dpi=300, format="png")
            plt.close(fig)
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

        fig, axes = setup_publication_subplots(2, 2, panel_type="quad")
        if title:
            fig.suptitle(title, fontsize=PUBLICATION_FONTS["title"], fontweight="bold")

        # Use color palette from publication_utils
        n_colors = len(rmsd_results)
        colors = get_color_palette("default", max(3, n_colors))

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

        set_publication_labels(ax1, "Time (ns)", "RMSD (Å)")
        ax1.legend()
        style_axes_for_publication(ax1, grid_alpha=0.3)

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
            set_publication_labels(ax2, "RMSD (Å)", "Frequency")
            style_axes_for_publication(ax2, grid_alpha=0.3)

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
            set_publication_labels(ax3, ylabel="Average RMSD (Å)")
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
            set_publication_labels(ax4, "Time (ns)", "Ligand RMSD (Å)")
            ax4.legend()
            style_axes_for_publication(ax4, grid_alpha=0.3)
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

        # Save using publication utils
        save_publication_figure(fig, output_path, dpi=300, format="png")
        plt.close(fig)

        return True

    except Exception as e:
        print(f"Error in RMSD plotting: {e}")
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
    fig, ax = setup_publication_subplots(figsize=figsize)

    if colors is None:
        # Use color palette from publication_utils
        colors = get_color_palette("default", len(rmsd_data))

    for i, (label, values) in enumerate(rmsd_data.items()):
        if times is None:
            plot_times = np.arange(len(values)) * 0.5  # Assume 0.5 ns per frame
        else:
            plot_times = times[: len(values)]

        ax.plot(plot_times, values, color=colors[i], label=label, linewidth=2, alpha=0.8)

    set_publication_labels(ax, "Time (ns)", "RMSD (nm)", title)
    style_axes_for_publication(ax, grid_alpha=0.3)
    ax.legend()

    if save_path:
        save_publication_figure(fig, save_path, dpi=300, format="png")
        logger.info(f"RMSD time series plot saved: {save_path}")

    return fig, ax


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
    if figsize is None:
        figsize = get_standard_figsize("vertical")
    fig, (ax1, ax2) = setup_publication_subplots(2, 1, figsize=figsize, panel_type="vertical")

    if times is None:
        times = np.arange(len(protein_rmsd))  # Frame numbers (no unit conversion needed)

    # Protein RMSD plot
    ax1.plot(times[: len(protein_rmsd)], protein_rmsd, color="blue", linewidth=2, alpha=0.8, label="Protein CA")
    ax1.fill_between(times[: len(protein_rmsd)], protein_rmsd, alpha=0.3, color="blue")
    set_publication_labels(ax1, "Frame", "RMSD (Å)", protein_title)
    style_axes_for_publication(ax1, grid_alpha=0.3)

    # Add protein statistics
    protein_mean = np.mean(protein_rmsd)
    protein_std = np.std(protein_rmsd)
    protein_stats = f"Mean: {protein_mean:.3f} ± {protein_std:.3f} Å"
    add_corner_annotation(ax1, protein_stats, facecolor="lightblue")

    # Ligand RMSD plot
    ax2.plot(times[: len(ligand_rmsd)], ligand_rmsd, color="red", linewidth=2, alpha=0.8, label="Ligand")
    ax2.fill_between(times[: len(ligand_rmsd)], ligand_rmsd, alpha=0.3, color="red")
    set_publication_labels(ax2, "Frame", "RMSD (Å)", ligand_title)
    style_axes_for_publication(ax2, grid_alpha=0.3)

    # Add ligand statistics
    ligand_mean = np.mean(ligand_rmsd)
    ligand_std = np.std(ligand_rmsd)
    ligand_stats = f"Mean: {ligand_mean:.2f} ± {ligand_std:.2f} Å"
    add_corner_annotation(ax2, ligand_stats, facecolor="lightcoral")

    plt.tight_layout()

    # Save individual plots if paths provided
    if protein_save_path:
        # Create individual protein plot
        fig_protein, ax_protein = setup_publication_subplots(figsize=(10, 5), panel_type="wide")
        ax_protein.plot(times[: len(protein_rmsd)], protein_rmsd, color="blue", linewidth=2, alpha=0.8)
        ax_protein.fill_between(times[: len(protein_rmsd)], protein_rmsd, alpha=0.3, color="blue")
        set_publication_labels(ax_protein, "Frame", "RMSD (Å)", protein_title)
        style_axes_for_publication(ax_protein, grid_alpha=0.3)
        add_corner_annotation(ax_protein, protein_stats, facecolor="lightblue")
        save_publication_figure(fig_protein, protein_save_path, dpi=300, format="png")
        plt.close(fig_protein)
        logger.info(f"Protein RMSD plot saved: {protein_save_path}")

    if ligand_save_path:
        # Create individual ligand plot
        fig_ligand, ax_ligand = setup_publication_subplots(figsize=(10, 5), panel_type="wide")
        ax_ligand.plot(times[: len(ligand_rmsd)], ligand_rmsd, color="red", linewidth=2, alpha=0.8)
        ax_ligand.fill_between(times[: len(ligand_rmsd)], ligand_rmsd, alpha=0.3, color="red")
        set_publication_labels(ax_ligand, "Frame", "RMSD (Å)", ligand_title)
        style_axes_for_publication(ax_ligand, grid_alpha=0.3)
        add_corner_annotation(ax_ligand, ligand_stats, facecolor="lightcoral")
        save_publication_figure(fig_ligand, ligand_save_path, dpi=300, format="png")
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
    fig, ax = setup_publication_subplots(figsize=figsize)

    if colors is None:
        # Use color palette from publication_utils
        n_colors = len(rmsd_data)
        colors = get_color_palette("default", max(5, n_colors))

    for i, (repeat_name, rmsd_values) in enumerate(rmsd_data.items()):
        if times is None:
            plot_times = np.arange(len(rmsd_values))  # Frame numbers (no unit conversion needed)
        else:
            plot_times = times[: len(rmsd_values)]

        color = colors[i % len(colors)]

        ax.plot(plot_times, rmsd_values, color=color, linestyle="-", linewidth=2.5, alpha=0.8, label=repeat_name)

    set_publication_labels(ax, "Frame", "Ligand RMSD (Å)", title)
    style_axes_for_publication(ax, grid_alpha=0.3)
    ax.legend(loc="best")

    # Add overall statistics
    all_values = np.concatenate([values for values in rmsd_data.values()])
    overall_mean = np.mean(all_values)
    overall_std = np.std(all_values)
    stats_text = f"Overall: {overall_mean:.2f} ± {overall_std:.2f} Å"

    add_corner_annotation(ax, stats_text, facecolor="wheat")

    if save_path:
        save_publication_figure(fig, save_path, dpi=300, format="png")
        logger.info(f"Multi-repeat ligand RMSD plot saved: {save_path}")

    return fig, ax
