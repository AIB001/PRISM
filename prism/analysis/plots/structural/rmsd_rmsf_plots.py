#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD and RMSF visualization functions.

Refactored from structural_plots.py for better modularity.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
import logging

# Import publication utilities
from ..publication_utils import PUBLICATION_FONTS, apply_publication_style, get_standard_figsize

# Import residue formatting utility
from pathlib import Path

# Apply global publication style on import
apply_publication_style()

logger = logging.getLogger(__name__)


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
