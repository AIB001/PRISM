#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structural analysis plotting utilities.

Based on SciDraft-Studio implementation with PRISM integration.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def plot_violin_comparison(data_dict: Dict[str, List[float]],
                          title: str = "Distribution Comparison",
                          xlabel: str = "Category",
                          ylabel: str = "Value",
                          figsize: Tuple[float, float] = (10, 6),
                          save_path: Optional[str] = None,
                          colors: Optional[List[str]] = None,
                          show_points: bool = True,
                          show_box: bool = True,
                          **kwargs) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create violin plots for comparing distributions.

    Parameters
    ----------
    data_dict : dict
        Dictionary mapping labels to data lists
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
    colors : list, optional
        Colors for violins
    show_points : bool
        Whether to show individual points
    show_box : bool
        Whether to show box plot inside violin

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Set colors
    if colors is None:
        colors = plt.cm.Set3(np.linspace(0, 1, len(data_dict)))

    # Create violin plot
    parts = ax.violinplot([data_dict[key] for key in data_dict.keys()],
                         positions=range(len(data_dict)),
                         showmeans=True, showmedians=True)

    # Color the violins
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.7)

    # Add box plot if requested
    if show_box:
        bp = ax.boxplot([data_dict[key] for key in data_dict.keys()],
                       positions=range(len(data_dict)),
                       widths=0.1, patch_artist=True,
                       boxprops=dict(facecolor='white', alpha=0.5),
                       showfliers=False)

    # Add individual points if requested
    if show_points:
        for i, (label, values) in enumerate(data_dict.items()):
            y = values
            x = np.random.normal(i, 0.04, size=len(y))
            ax.scatter(x, y, alpha=0.3, s=20, color='black')

    # Style
    ax.set_xticks(range(len(data_dict)))
    ax.set_xticklabels(list(data_dict.keys()))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    # Add statistics text
    stats_text = "Mean ± SD:\n"
    for label, values in data_dict.items():
        stats_text += f"{label}: {np.mean(values):.2f} ± {np.std(values):.2f}\n"

    ax.text(0.02, 0.98, stats_text.strip(), transform=ax.transAxes,
            verticalalignment='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    # Save if requested
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Violin plot saved: {save_path}")

    return fig, ax


def plot_ramachandran(phi_angles: np.ndarray,
                     psi_angles: np.ndarray,
                     title: str = "Ramachandran Plot",
                     figsize: Tuple[float, float] = (8, 8),
                     save_path: Optional[str] = None,
                     density: bool = True,
                     alpha: float = 0.6) -> Tuple[plt.Figure, plt.Axes]:
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
    fig, ax = plt.subplots(figsize=figsize)

    # Flatten arrays if needed
    phi_flat = phi_angles.flatten() if phi_angles.ndim > 1 else phi_angles
    psi_flat = psi_angles.flatten() if psi_angles.ndim > 1 else psi_angles

    if density:
        # Create density plot
        ax.hist2d(phi_flat, psi_flat, bins=50, density=True,
                 cmap='Blues', alpha=0.8)
        ax.scatter(phi_flat, psi_flat, alpha=alpha, s=1, color='red')
    else:
        # Simple scatter plot
        ax.scatter(phi_flat, psi_flat, alpha=alpha, s=10)

    # Add reference regions for common secondary structures
    # Alpha helix region
    ax.add_patch(plt.Rectangle((-120, -60), 60, 60,
                              fill=False, edgecolor='green',
                              linewidth=2, linestyle='--',
                              label='α-helix'))

    # Beta sheet region
    ax.add_patch(plt.Rectangle((-180, 120), 120, 60,
                              fill=False, edgecolor='blue',
                              linewidth=2, linestyle='--',
                              label='β-sheet'))

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel('φ (degrees)')
    ax.set_ylabel('ψ (degrees)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Add axis lines at 0
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Ramachandran plot saved: {save_path}")

    return fig, ax


def plot_dihedral_time_series(angles: np.ndarray,
                             times: Optional[np.ndarray] = None,
                             title: str = "Dihedral Angle Time Series",
                             ylabel: str = "Angle (degrees)",
                             figsize: Tuple[float, float] = (12, 6),
                             save_path: Optional[str] = None,
                             residue_ids: Optional[List[int]] = None,
                             max_traces: int = 10) -> Tuple[plt.Figure, plt.Axes]:
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
    fig, ax = plt.subplots(figsize=figsize)

    if times is None:
        times = np.arange(angles.shape[0])

    # Limit number of traces for readability
    n_residues = min(angles.shape[1], max_traces)
    colors = plt.cm.tab10(np.linspace(0, 1, n_residues))

    for i in range(n_residues):
        label = f"Res {residue_ids[i]}" if residue_ids else f"Residue {i+1}"
        ax.plot(times, angles[:, i], color=colors[i], alpha=0.8,
               linewidth=1.5, label=label)

    ax.set_xlabel('Time (ns)' if times is not None else 'Frame')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    if n_residues <= 10:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Dihedral time series plot saved: {save_path}")

    return fig, ax


def plot_sasa_comparison(sasa_data: Dict[str, np.ndarray],
                        title: str = "SASA Comparison",
                        figsize: Tuple[float, float] = (10, 6),
                        save_path: Optional[str] = None,
                        plot_type: str = "box") -> Tuple[plt.Figure, plt.Axes]:
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
    fig, ax = plt.subplots(figsize=figsize)

    if plot_type == "box":
        # Box plot
        data_list = [sasa_data[key] for key in sasa_data.keys()]
        labels = list(sasa_data.keys())

        bp = ax.boxplot(data_list, labels=labels, patch_artist=True)

        # Color boxes
        colors = plt.cm.Set3(np.linspace(0, 1, len(data_list)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

    elif plot_type == "violin":
        # Use violin plot function
        data_dict = {key: list(values) for key, values in sasa_data.items()}
        return plot_violin_comparison(data_dict, title=title,
                                    ylabel="SASA (Ų)",
                                    figsize=figsize, save_path=save_path)

    elif plot_type == "time_series":
        # Time series plot
        colors = plt.cm.tab10(np.linspace(0, 1, len(sasa_data)))
        for i, (label, values) in enumerate(sasa_data.items()):
            times = np.arange(len(values)) * 0.5  # Assume 0.5 ns per frame
            ax.plot(times, values, color=colors[i], label=label, linewidth=2)

        ax.set_xlabel('Time (ns)')
        ax.legend()

    ax.set_ylabel('SASA (Ų)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"SASA comparison plot saved: {save_path}")

    return fig, ax


def plot_property_distribution(values: np.ndarray,
                              title: str = "Property Distribution",
                              xlabel: str = "Value",
                              ylabel: str = "Frequency",
                              figsize: Tuple[float, float] = (8, 6),
                              save_path: Optional[str] = None,
                              bins: int = 50,
                              show_stats: bool = True) -> Tuple[plt.Figure, plt.Axes]:
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
    fig, ax = plt.subplots(figsize=figsize)

    # Create histogram
    n, bins_edges, patches = ax.hist(values, bins=bins, alpha=0.7,
                                    color='skyblue', edgecolor='black')

    # Add statistical lines
    mean_val = np.mean(values)
    std_val = np.std(values)
    median_val = np.median(values)

    ax.axvline(mean_val, color='red', linestyle='--', linewidth=2,
              label=f'Mean: {mean_val:.2f}')
    ax.axvline(median_val, color='green', linestyle='--', linewidth=2,
              label=f'Median: {median_val:.2f}')

    # Add stats text
    if show_stats:
        stats_text = (f"Mean: {mean_val:.2f}\n"
                      f"Std: {std_val:.2f}\n"
                      f"Median: {median_val:.2f}\n"
                      f"Min: {np.min(values):.2f}\n"
                      f"Max: {np.max(values):.2f}")

        ax.text(0.75, 0.95, stats_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Property distribution plot saved: {save_path}")

    return fig, ax


def plot_rmsd_time_series(rmsd_data: Dict[str, np.ndarray],
                         times: Optional[np.ndarray] = None,
                         title: str = "RMSD Time Series",
                         figsize: Tuple[float, float] = (12, 6),
                         save_path: Optional[str] = None,
                         colors: Optional[List[str]] = None) -> Tuple[plt.Figure, plt.Axes]:
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
    fig, ax = plt.subplots(figsize=figsize)

    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(rmsd_data)))

    for i, (label, values) in enumerate(rmsd_data.items()):
        if times is None:
            plot_times = np.arange(len(values)) * 0.5  # Assume 0.5 ns per frame
        else:
            plot_times = times[:len(values)]

        ax.plot(plot_times, values, color=colors[i], label=label,
               linewidth=2, alpha=0.8)

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (nm)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"RMSD time series plot saved: {save_path}")

    return fig, ax


def plot_rmsf_per_residue(rmsf_values: np.ndarray,
                         residue_ids: Optional[np.ndarray] = None,
                         title: str = "RMSF per Residue",
                         figsize: Tuple[float, float] = (12, 6),
                         save_path: Optional[str] = None,
                         highlight_threshold: Optional[float] = None,
                         secondary_structure: Optional[Dict] = None) -> Tuple[plt.Figure, plt.Axes]:
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
    fig, ax = plt.subplots(figsize=figsize)

    if residue_ids is None:
        residue_ids = np.arange(1, len(rmsf_values) + 1)

    # Main RMSF plot
    ax.plot(residue_ids, rmsf_values, 'b-', linewidth=2, alpha=0.8)
    ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color='blue')

    # Highlight flexible regions if threshold provided
    if highlight_threshold is not None:
        flexible_mask = rmsf_values > highlight_threshold
        if np.any(flexible_mask):
            ax.fill_between(residue_ids, rmsf_values, highlight_threshold,
                          where=flexible_mask, alpha=0.5, color='red',
                          label=f'Flexible (RMSF > {highlight_threshold:.2f} nm)')

    # Add secondary structure annotations if provided
    if secondary_structure:
        y_max = np.max(rmsf_values) * 1.1
        for ss_type, regions in secondary_structure.items():
            color = {'helix': 'red', 'sheet': 'green', 'loop': 'gray'}.get(ss_type, 'black')
            for start, end in regions:
                ax.axvspan(start, end, alpha=0.2, color=color, label=f'{ss_type.capitalize()}')

    ax.set_xlabel('Residue Number')
    ax.set_ylabel('RMSF (nm)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    # Add statistics
    mean_rmsf = np.mean(rmsf_values)
    std_rmsf = np.std(rmsf_values)
    stats_text = f"Mean RMSF: {mean_rmsf:.3f} ± {std_rmsf:.3f} nm"

    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            verticalalignment='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    if highlight_threshold is not None and np.any(flexible_mask):
        ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"RMSF plot saved: {save_path}")

    return fig, ax


def plot_rmsd_rmsf_combined(rmsd_data: Dict[str, np.ndarray],
                           rmsf_data: Dict[str, np.ndarray],
                           times: Optional[np.ndarray] = None,
                           residue_ids: Optional[np.ndarray] = None,
                           title: str = "RMSD and RMSF Analysis",
                           figsize: Tuple[float, float] = (15, 8),
                           save_path: Optional[str] = None) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
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
    fig.suptitle(title, fontsize=16, fontweight='bold')

    # RMSD subplot
    colors = plt.cm.tab10(np.linspace(0, 1, len(rmsd_data)))
    for i, (label, values) in enumerate(rmsd_data.items()):
        if times is None:
            plot_times = np.arange(len(values)) * 0.5
        else:
            plot_times = times[:len(values)]

        ax1.plot(plot_times, values, color=colors[i], label=label,
                linewidth=2, alpha=0.8)

    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('RMSD (nm)')
    ax1.set_title('RMSD Time Series')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # RMSF subplot
    for i, (label, values) in enumerate(rmsf_data.items()):
        if residue_ids is None:
            plot_residues = np.arange(1, len(values) + 1)
        else:
            plot_residues = residue_ids[:len(values)]

        ax2.plot(plot_residues, values, color=colors[i], label=label,
                linewidth=2, alpha=0.8)

    ax2.set_xlabel('Residue Number')
    ax2.set_ylabel('RMSF (nm)')
    ax2.set_title('RMSF per Residue')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Combined RMSD/RMSF plot saved: {save_path}")

def plot_multi_chain_rmsf(chain_rmsf_data: Dict[str, Tuple[np.ndarray, np.ndarray]],
                           title: str = "Multi-Chain RMSF Analysis",
                           figsize: Optional[Tuple[float, float]] = None,
                           save_path: Optional[str] = None,
                           cols: int = 2) -> Tuple[plt.Figure, List[plt.Axes]]:
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
        axes = list(axes) if hasattr(axes, '__iter__') else [axes]
    else:
        axes = axes.flatten()

    # Colors for different chain types
    protein_colors = plt.cm.Blues(np.linspace(0.4, 0.9, sum(1 for name in chain_rmsf_data.keys() if 'Protein' in name or 'Chain' in name)))
    nucleic_colors = plt.cm.Reds(np.linspace(0.4, 0.9, sum(1 for name in chain_rmsf_data.keys() if 'Nucleic' in name)))

    protein_idx, nucleic_idx = 0, 0

    for i, (chain_name, (rmsf_values, residue_ids)) in enumerate(chain_rmsf_data.items()):
        ax = axes[i]

        # Choose color based on chain type
        if 'Nucleic' in chain_name:
            color = nucleic_colors[nucleic_idx] if nucleic_idx < len(nucleic_colors) else 'red'
            nucleic_idx += 1
        else:
            color = protein_colors[protein_idx] if protein_idx < len(protein_colors) else 'blue'
            protein_idx += 1

        # Plot RMSF
        ax.plot(residue_ids, rmsf_values, color=color, linewidth=2, alpha=0.8)
        ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color=color)

        # Formatting
        ax.set_xlabel('Residue/Nucleotide Number')
        ax.set_ylabel('RMSF (Å)')
        ax.set_title(f'{chain_name} (n={len(rmsf_values)})')
        ax.grid(True, alpha=0.3)

        # Add statistics
        mean_rmsf = np.mean(rmsf_values)
        std_rmsf = np.std(rmsf_values)
        stats_text = f"Mean: {mean_rmsf:.2f} ± {std_rmsf:.2f} Å"

        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Hide unused subplots
    for i in range(n_chains, len(axes)):
        axes[i].set_visible(False)

    plt.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Multi-chain RMSF plot saved: {save_path}")

    return fig, axes[:n_chains]


def plot_separate_rmsd(protein_rmsd: np.ndarray,
                       ligand_rmsd: np.ndarray,
                       times: Optional[np.ndarray] = None,
                       protein_title: str = "Protein RMSD (PR1 Chain)",
                       ligand_title: str = "Ligand RMSD (after protein alignment)",
                       figsize: Tuple[float, float] = (12, 8),
                       protein_save_path: Optional[str] = None,
                       ligand_save_path: Optional[str] = None) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
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
    ax1.plot(times[:len(protein_rmsd)], protein_rmsd, color='blue', linewidth=2, alpha=0.8, label='Protein CA')
    ax1.fill_between(times[:len(protein_rmsd)], protein_rmsd, alpha=0.3, color='blue')
    ax1.set_xlabel('Frame')
    ax1.set_ylabel('RMSD (Å)')
    ax1.set_title(protein_title)
    ax1.grid(True, alpha=0.3)

    # Add protein statistics
    protein_mean = np.mean(protein_rmsd)
    protein_std = np.std(protein_rmsd)
    protein_stats = f"Mean: {protein_mean:.3f} ± {protein_std:.3f} Å"
    ax1.text(0.02, 0.98, protein_stats, transform=ax1.transAxes,
             verticalalignment='top', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

    # Ligand RMSD plot
    ax2.plot(times[:len(ligand_rmsd)], ligand_rmsd, color='red', linewidth=2, alpha=0.8, label='Ligand')
    ax2.fill_between(times[:len(ligand_rmsd)], ligand_rmsd, alpha=0.3, color='red')
    ax2.set_xlabel('Frame')
    ax2.set_ylabel('RMSD (Å)')
    ax2.set_title(ligand_title)
    ax2.grid(True, alpha=0.3)

    # Add ligand statistics
    ligand_mean = np.mean(ligand_rmsd)
    ligand_std = np.std(ligand_rmsd)
    ligand_stats = f"Mean: {ligand_mean:.2f} ± {ligand_std:.2f} Å"
    ax2.text(0.02, 0.98, ligand_stats, transform=ax2.transAxes,
             verticalalignment='top', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))

    plt.tight_layout()

    # Save individual plots if paths provided
    if protein_save_path:
        # Create individual protein plot
        fig_protein, ax_protein = plt.subplots(figsize=(10, 5))
        ax_protein.plot(times[:len(protein_rmsd)], protein_rmsd, color='blue', linewidth=2, alpha=0.8)
        ax_protein.fill_between(times[:len(protein_rmsd)], protein_rmsd, alpha=0.3, color='blue')
        ax_protein.set_xlabel('Frame')
        ax_protein.set_ylabel('RMSD (Å)')
        ax_protein.set_title(protein_title)
        ax_protein.grid(True, alpha=0.3)
        ax_protein.text(0.02, 0.98, protein_stats, transform=ax_protein.transAxes,
                       verticalalignment='top', fontsize=10,
                       bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        plt.tight_layout()
        fig_protein.savefig(protein_save_path, dpi=300, bbox_inches='tight')
        plt.close(fig_protein)
        logger.info(f"Protein RMSD plot saved: {protein_save_path}")

    if ligand_save_path:
        # Create individual ligand plot
        fig_ligand, ax_ligand = plt.subplots(figsize=(10, 5))
        ax_ligand.plot(times[:len(ligand_rmsd)], ligand_rmsd, color='red', linewidth=2, alpha=0.8)
        ax_ligand.fill_between(times[:len(ligand_rmsd)], ligand_rmsd, alpha=0.3, color='red')
        ax_ligand.set_xlabel('Frame')
        ax_ligand.set_ylabel('RMSD (Å)')
        ax_ligand.set_title(ligand_title)
        ax_ligand.grid(True, alpha=0.3)
        ax_ligand.text(0.02, 0.98, ligand_stats, transform=ax_ligand.transAxes,
                      verticalalignment='top', fontsize=10,
                      bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
        plt.tight_layout()
        fig_ligand.savefig(ligand_save_path, dpi=300, bbox_inches='tight')
        plt.close(fig_ligand)
        logger.info(f"Ligand RMSD plot saved: {ligand_save_path}")

    return fig, (ax1, ax2)


def plot_multi_repeat_ligand_rmsd(rmsd_data: Dict[str, np.ndarray],
                                  times: Optional[np.ndarray] = None,
                                  title: str = "Ligand RMSD - Multiple Repeats",
                                  figsize: Tuple[float, float] = (12, 6),
                                  save_path: Optional[str] = None,
                                  colors: Optional[List[str]] = None) -> Tuple[plt.Figure, plt.Axes]:
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
    fig, ax = plt.subplots(figsize=figsize)

    if colors is None:
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    for i, (repeat_name, rmsd_values) in enumerate(rmsd_data.items()):
        if times is None:
            plot_times = np.arange(len(rmsd_values))  # Frame numbers (no unit conversion needed)
        else:
            plot_times = times[:len(rmsd_values)]

        color = colors[i % len(colors)]

        ax.plot(plot_times, rmsd_values, color=color, linestyle='-',
               linewidth=2.5, alpha=0.8, label=repeat_name)

    ax.set_xlabel('Frame')
    ax.set_ylabel('Ligand RMSD (Å)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')

    # Add overall statistics
    all_values = np.concatenate([values for values in rmsd_data.values()])
    overall_mean = np.mean(all_values)
    overall_std = np.std(all_values)
    stats_text = f"Overall: {overall_mean:.2f} ± {overall_std:.2f} Å"

    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            verticalalignment='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Multi-repeat ligand RMSD plot saved: {save_path}")

    return fig, ax


    return fig, (ax1, ax2)