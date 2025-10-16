#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Specialized visualization functions.

Refactored from structural_plots.py for better modularity.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Dict, List, Optional, Tuple
import logging

# Import publication utilities
from ..publication_utils import (
    get_publication_style, get_color_palette, fix_rotated_labels,
    add_statistical_annotations, style_axes_for_publication,
    save_publication_figure, validate_figure_size, PUBLICATION_FONTS, PUBLICATION_COLORS,
    apply_publication_style, setup_publication_figure, get_standard_figsize
)
# Import residue formatting utility
from ....utils.residue import format_residue_list
from pathlib import Path

# Apply global publication style on import
apply_publication_style()

logger = logging.getLogger(__name__)


def plot_violin_comparison(data_dict: Dict[str, List[float]],
                          title: str = "",
                          xlabel: str = "Category",
                          ylabel: str = "Value",
                          figsize: Optional[Tuple[float, float]] = None,
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
    if figsize is None:
        figsize = get_standard_figsize("single")
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

    ax.grid(True, alpha=0.3)

    # Print statistics (removed from plot to avoid hiding content)
    print("  📊 Property statistics (Mean ± SD):")
    for label, values in data_dict.items():
        print(f"    {label}: {np.mean(values):.2f} ± {np.std(values):.2f}")

    plt.tight_layout()

    # Save if requested
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Violin plot saved: {save_path}")

    return fig, ax




def plot_distance_time_series(
    distance_data: Dict[str, np.ndarray],
    times: Optional[np.ndarray] = None,
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Create 6-panel distance time series plot for key interactions.

    Parameters
    ----------
    distance_data : dict
        Dictionary mapping interaction names to distance time series
    times : np.ndarray, optional
        Time points. If None, use frame indices converted to ns
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes_list) objects
    """
    if not distance_data:
        raise ValueError("No distance data provided")

    # Apply publication style
    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        # Create 2x3 subplot grid
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

        interactions = list(distance_data.keys())[:6]  # Limit to 6
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

        for i, interaction in enumerate(interactions):
            ax = axes[i]
            distances = distance_data[interaction]

            if times is None:
                plot_times = np.arange(len(distances)) * 0.02  # 20 ps -> ns
            else:
                plot_times = times[:len(distances)]

            # Plot time series with error band
            ax.plot(plot_times, distances, color=colors[i],
                   linewidth=1.5, alpha=0.8)

            # Add running average
            window = min(50, len(distances)//10)
            if window > 1:
                running_avg = np.convolve(distances, np.ones(window)/window, mode='valid')
                avg_times = plot_times[window//2:len(running_avg)+window//2]
                ax.plot(avg_times, running_avg, 'r-', linewidth=2.5, alpha=0.9)

            # Add mean line
            mean_dist = np.mean(distances)
            ax.axhline(mean_dist, color='black', linestyle='--',
                      linewidth=1.5, alpha=0.7)

            # Print statistics (removed from plot to avoid hiding content)
            print(f"  📊 {interaction} distance: Mean = {mean_dist:.2f} ± {np.std(distances):.2f} Å")

            # Styling
            ax.set_xlabel('Time (ns)', fontsize=PUBLICATION_FONTS['axis_label'])
            ax.set_ylabel('Distance (Å)', fontsize=PUBLICATION_FONTS['axis_label'])

            ax.grid(True, alpha=0.3)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        # Hide unused subplots
        for i in range(len(interactions), len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"Distance time series plot saved: {save_path}")

    return fig, axes




def plot_magnesium_coordination(
    mg_data: Dict[str, Tuple[List[float], List[float]]],
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Create magnesium coordination analysis with dual scatter plots.

    Parameters
    ----------
    mg_data : dict
        Dictionary mapping Mg2+ sites to (distances_to_ligand, distances_to_protein) tuples
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes_list) objects
    """
    if not mg_data:
        raise ValueError("No magnesium coordination data provided")

    # Apply publication style
    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        fig, axes = plt.subplots(1, 2, figsize=figsize)

        colors = ['#FF6B6B', '#4ECDC4']  # Red and teal for Mg1 and Mg2

        mg_sites = list(mg_data.keys())

        for i, (site, (lig_dists, prot_dists)) in enumerate(mg_data.items()):
            color = colors[i % len(colors)]

            # Left panel: Mg-ligand distances vs time
            ax1 = axes[0]
            time_ns = np.arange(len(lig_dists)) * 0.02
            ax1.scatter(time_ns, lig_dists, color=color, alpha=0.6,
                       s=30, label=f'{site} (Ligand)', edgecolor='black', linewidth=0.5)

            # Right panel: Mg-protein distances vs time
            ax2 = axes[1]
            ax2.scatter(time_ns, prot_dists, color=color, alpha=0.6,
                       s=30, label=f'{site} (Protein)', edgecolor='black', linewidth=0.5)

        # Style left panel
        ax1.set_xlabel('Time (ns)', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax1.set_ylabel('Mg²⁺-Ligand Distance (Å)', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=PUBLICATION_FONTS['legend'])
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        # Add ideal coordination distance line
        ax1.axhline(2.1, color='green', linestyle='--', linewidth=2,
                   alpha=0.7, label='Ideal (2.1Å)')

        # Style right panel
        ax2.set_xlabel('Time (ns)', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax2.set_ylabel('Mg²⁺-Protein Distance (Å)', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=PUBLICATION_FONTS['legend'])
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

        # Add statistics
        for i, (site, (lig_dists, prot_dists)) in enumerate(mg_data.items()):
            lig_mean = np.mean(lig_dists)
            prot_mean = np.mean(prot_dists)

            ax1.text(0.02, 0.98 - i*0.1, f'{site}: {lig_mean:.2f}±{np.std(lig_dists):.2f}Å',
                    transform=ax1.transAxes, verticalalignment='top',
                    fontsize=PUBLICATION_FONTS['annotation'],
                    bbox=dict(boxstyle='round,pad=0.3', facecolor=colors[i], alpha=0.3))

            ax2.text(0.02, 0.98 - i*0.1, f'{site}: {prot_mean:.2f}±{np.std(prot_dists):.2f}Å',
                    transform=ax2.transAxes, verticalalignment='top',
                    fontsize=PUBLICATION_FONTS['annotation'],
                    bbox=dict(boxstyle='round,pad=0.3', facecolor=colors[i], alpha=0.3))

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"Magnesium coordination plot saved: {save_path}")

    return fig, axes




