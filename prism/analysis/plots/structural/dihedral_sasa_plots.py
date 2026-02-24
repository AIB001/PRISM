#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dihedral angle and SASA visualization functions.

Refactored from structural_plots.py for better modularity.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
import logging

# Import publication utilities
from ..publication_utils import apply_publication_style, get_standard_figsize
# Import residue formatting utility

# Apply global publication style on import
apply_publication_style()

logger = logging.getLogger(__name__)


def plot_ramachandran(
    phi_angles: np.ndarray,
    psi_angles: np.ndarray,
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    density: bool = True,
    alpha: float = 0.6,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a Ramachandran plot (phi vs psi angles).

    Parameters
    ----------
    phi_angles : np.ndarray
        Phi dihedral angles in degrees
    psi_angles : np.ndarray
        Psi dihedral angles in degrees
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    density : bool
        Whether to show density contours
    alpha : float
        Point transparency

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("single")
    fig, ax = plt.subplots(figsize=figsize)

    # Flatten arrays if needed
    phi_flat = phi_angles.flatten() if phi_angles.ndim > 1 else phi_angles
    psi_flat = psi_angles.flatten() if psi_angles.ndim > 1 else psi_angles

    if density:
        # Create density plot
        ax.hist2d(phi_flat, psi_flat, bins=50, density=True, cmap="Blues", alpha=0.8)
        ax.scatter(phi_flat, psi_flat, alpha=alpha, s=1, color="red")
    else:
        # Simple scatter plot
        ax.scatter(phi_flat, psi_flat, alpha=alpha, s=10)

    # Add reference regions for common secondary structures
    # Alpha helix region
    ax.add_patch(
        plt.Rectangle((-120, -60), 60, 60, fill=False, edgecolor="green", linewidth=2, linestyle="--", label="α-helix")
    )

    # Beta sheet region
    ax.add_patch(
        plt.Rectangle((-180, 120), 120, 60, fill=False, edgecolor="blue", linewidth=2, linestyle="--", label="β-sheet")
    )

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("φ (degrees)")
    ax.set_ylabel("ψ (degrees)")
    if title:
        ax.set_title(title)

    ax.grid(True, alpha=0.3)
    ax.legend()

    # Add axis lines at 0
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Ramachandran plot saved: {save_path}")

    return fig, ax


def plot_dihedral_time_series(
    angles: np.ndarray,
    times: Optional[np.ndarray] = None,
    title: str = "",
    ylabel: str = "Angle (degrees)",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    residue_ids: Optional[List[int]] = None,
    max_traces: int = 10,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot dihedral angle time series.

    Parameters
    ----------
    angles : np.ndarray
        Dihedral angles array (n_frames, n_residues)
    times : np.ndarray, optional
        Time points. If None, use frame indices.
    title : str
        Plot title
    ylabel : str
        Y-axis label
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    residue_ids : list, optional
        Residue IDs for labeling
    max_traces : int
        Maximum number of traces to plot

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("single")
    fig, ax = plt.subplots(figsize=figsize)

    if times is None:
        times = np.arange(angles.shape[0])

    # Limit number of traces for readability
    n_residues = min(angles.shape[1], max_traces)
    colors = plt.cm.tab10(np.linspace(0, 1, n_residues))

    for i in range(n_residues):
        label = f"Res {residue_ids[i]}" if residue_ids else f"Residue {i+1}"
        ax.plot(times, angles[:, i], color=colors[i], alpha=0.8, linewidth=1.5, label=label)

    ax.set_xlabel("Time (ns)" if times is not None else "Frame")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    ax.grid(True, alpha=0.3)

    if n_residues <= 10:
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Dihedral time series plot saved: {save_path}")

    return fig, ax


def plot_sasa_comparison(
    sasa_data: Dict[str, np.ndarray],
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    plot_type: str = "box",
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot SASA comparison between different systems.

    Parameters
    ----------
    sasa_data : dict
        Dictionary with system names as keys and SASA arrays as values
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    plot_type : str
        Type of plot ('box', 'violin', 'time_series')

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("single")
    fig, ax = plt.subplots(figsize=figsize)

    if plot_type == "box":
        # Box plot
        data_list = [sasa_data[key] for key in sasa_data.keys()]
        labels = list(sasa_data.keys())

        bp = ax.boxplot(data_list, labels=labels, patch_artist=True)

        # Color boxes
        colors = plt.cm.Set3(np.linspace(0, 1, len(data_list)))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

    elif plot_type == "violin":
        # Use violin plot function
        data_dict = {key: list(values) for key, values in sasa_data.items()}
        return plot_violin_comparison(data_dict, title=title, ylabel="SASA (Ų)", figsize=figsize, save_path=save_path)

    elif plot_type == "time_series":
        # Time series plot
        colors = plt.cm.tab10(np.linspace(0, 1, len(sasa_data)))
        for i, (label, values) in enumerate(sasa_data.items()):
            times = np.arange(len(values)) * 0.5  # Assume 0.5 ns per frame
            ax.plot(times, values, color=colors[i], label=label, linewidth=2)

        ax.set_xlabel("Time (ns)")
        ax.legend()

    ax.set_ylabel("SASA (Ų)")
    if title:
        ax.set_title(title)

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"SASA comparison plot saved: {save_path}")

    return fig, ax


def plot_property_distribution(
    values: np.ndarray,
    title: str = "",
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    bins: int = 50,
    show_stats: bool = True,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot distribution of a property with statistics.

    Parameters
    ----------
    values : np.ndarray
        Property values
    title : str
        Plot title
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    bins : int
        Number of histogram bins
    show_stats : bool
        Whether to show statistics

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("single")
    fig, ax = plt.subplots(figsize=figsize)

    # Create histogram
    n, bins_edges, patches = ax.hist(values, bins=bins, alpha=0.7, color="skyblue", edgecolor="black")

    # Add statistical lines
    mean_val = np.mean(values)
    std_val = np.std(values)
    median_val = np.median(values)

    ax.axvline(mean_val, color="red", linestyle="--", linewidth=2, label=f"Mean: {mean_val:.2f}")
    ax.axvline(median_val, color="green", linestyle="--", linewidth=2, label=f"Median: {median_val:.2f}")

    # Add stats text
    if show_stats:
        # Print statistics (removed from plot to avoid hiding content)
        print(f"  📊 Distribution statistics:")
        print(f"    Mean: {mean_val:.2f}, Std: {std_val:.2f}")
        print(f"    Median: {median_val:.2f}, Min: {np.min(values):.2f}, Max: {np.max(values):.2f}")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Property distribution plot saved: {save_path}")

    return fig, ax
