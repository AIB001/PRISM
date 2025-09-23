#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparison plotting utilities for multi-system analysis.

Based on SciDraft-Studio implementation with PRISM integration.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple, Union
import logging

logger = logging.getLogger(__name__)


def plot_multi_system_comparison(data_dict: Dict[str, np.ndarray],
                                property_name: str = "Property",
                                figsize: Tuple[float, float] = (12, 8),
                                save_path: Optional[str] = None,
                                plot_types: List[str] = ["time_series", "distribution"],
                                colors: Optional[List[str]] = None) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Create comprehensive comparison plots for multiple systems.

    Parameters
    ----------
    data_dict : dict
        Dictionary with system names as keys and data arrays as values
    property_name : str
        Name of the property being plotted
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    plot_types : list
        Types of plots to include
    colors : list, optional
        Colors for each system

    Returns
    -------
    tuple
        (figure, list of axes) objects
    """
    n_plots = len(plot_types)
    fig, axes = plt.subplots(1, n_plots, figsize=figsize)

    if n_plots == 1:
        axes = [axes]

    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(data_dict)))

    for i, plot_type in enumerate(plot_types):
        ax = axes[i]

        if plot_type == "time_series":
            _plot_time_series_comparison(data_dict, ax, colors, property_name)
        elif plot_type == "distribution":
            _plot_distribution_comparison(data_dict, ax, colors, property_name)
        elif plot_type == "box":
            _plot_box_comparison(data_dict, ax, colors, property_name)
        elif plot_type == "violin":
            _plot_violin_comparison(data_dict, ax, colors, property_name)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Multi-system comparison plot saved: {save_path}")

    return fig, axes


def _plot_time_series_comparison(data_dict: Dict[str, np.ndarray],
                                ax: plt.Axes,
                                colors: List,
                                property_name: str):
    """Plot time series comparison."""
    for i, (label, values) in enumerate(data_dict.items()):
        times = np.arange(len(values)) * 0.5  # Assume 0.5 ns per frame
        ax.plot(times, values, color=colors[i], label=label, linewidth=2, alpha=0.8)

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(property_name)
    ax.set_title(f'{property_name} Time Series')
    ax.grid(True, alpha=0.3)
    ax.legend()


def _plot_distribution_comparison(data_dict: Dict[str, np.ndarray],
                                 ax: plt.Axes,
                                 colors: List,
                                 property_name: str):
    """Plot distribution comparison."""
    for i, (label, values) in enumerate(data_dict.items()):
        ax.hist(values, bins=30, alpha=0.6, color=colors[i], label=label,
               density=True, edgecolor='black', linewidth=0.5)

    ax.set_xlabel(property_name)
    ax.set_ylabel('Density')
    ax.set_title(f'{property_name} Distribution')
    ax.grid(True, alpha=0.3)
    ax.legend()


def _plot_box_comparison(data_dict: Dict[str, np.ndarray],
                        ax: plt.Axes,
                        colors: List,
                        property_name: str):
    """Plot box comparison."""
    data_list = [data_dict[key] for key in data_dict.keys()]
    labels = list(data_dict.keys())

    bp = ax.boxplot(data_list, labels=labels, patch_artist=True)

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel(property_name)
    ax.set_title(f'{property_name} Box Plot')
    ax.grid(True, alpha=0.3)


def _plot_violin_comparison(data_dict: Dict[str, np.ndarray],
                           ax: plt.Axes,
                           colors: List,
                           property_name: str):
    """Plot violin comparison."""
    data_list = [data_dict[key] for key in data_dict.keys()]
    labels = list(data_dict.keys())

    parts = ax.violinplot(data_list, positions=range(len(data_list)),
                         showmeans=True, showmedians=True)

    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.7)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylabel(property_name)
    ax.set_title(f'{property_name} Violin Plot')
    ax.grid(True, alpha=0.3)


def plot_correlation_matrix(data_dict: Dict[str, np.ndarray],
                           figsize: Tuple[float, float] = (10, 8),
                           save_path: Optional[str] = None,
                           method: str = "pearson") -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot correlation matrix between different properties or systems.

    Parameters
    ----------
    data_dict : dict
        Dictionary with property/system names as keys and data arrays as values
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    method : str
        Correlation method ('pearson', 'spearman', 'kendall')

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    import pandas as pd

    # Create DataFrame
    min_length = min(len(values) for values in data_dict.values())
    trimmed_data = {key: values[:min_length] for key, values in data_dict.items()}
    df = pd.DataFrame(trimmed_data)

    # Calculate correlation matrix
    corr_matrix = df.corr(method=method)

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Create heatmap
    im = ax.imshow(corr_matrix.values, cmap='coolwarm', aspect='auto',
                   vmin=-1, vmax=1)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Correlation Coefficient')

    # Set ticks and labels
    ax.set_xticks(range(len(corr_matrix.columns)))
    ax.set_yticks(range(len(corr_matrix.index)))
    ax.set_xticklabels(corr_matrix.columns, rotation=45, ha='right')
    ax.set_yticklabels(corr_matrix.index)

    # Add correlation values as text
    for i in range(len(corr_matrix.index)):
        for j in range(len(corr_matrix.columns)):
            text = ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}',
                          ha="center", va="center", color="black")

    ax.set_title(f'Correlation Matrix ({method.capitalize()})')
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Correlation matrix plot saved: {save_path}")

    return fig, ax


def plot_statistical_comparison(data_dict: Dict[str, np.ndarray],
                               figsize: Tuple[float, float] = (12, 8),
                               save_path: Optional[str] = None,
                               property_name: str = "Property") -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Create comprehensive statistical comparison plots.

    Parameters
    ----------
    data_dict : dict
        Dictionary with system names as keys and data arrays as values
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    property_name : str
        Name of the property

    Returns
    -------
    tuple
        (figure, list of axes) objects
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    colors = plt.cm.tab10(np.linspace(0, 1, len(data_dict)))

    # 1. Mean and standard deviation
    ax = axes[0]
    means = [np.mean(values) for values in data_dict.values()]
    stds = [np.std(values) for values in data_dict.values()]
    labels = list(data_dict.keys())

    x_pos = np.arange(len(labels))
    bars = ax.bar(x_pos, means, yerr=stds, capsize=5, color=colors, alpha=0.7)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel(f'{property_name} (Mean ± SD)')
    ax.set_title('Mean Values Comparison')
    ax.grid(True, alpha=0.3)

    # 2. Box plot
    ax = axes[1]
    data_list = [data_dict[key] for key in data_dict.keys()]
    bp = ax.boxplot(data_list, labels=labels, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel(property_name)
    ax.set_title('Distribution Comparison')
    ax.grid(True, alpha=0.3)

    # 3. Cumulative distribution
    ax = axes[2]
    for i, (label, values) in enumerate(data_dict.items()):
        sorted_values = np.sort(values)
        y = np.arange(1, len(sorted_values) + 1) / len(sorted_values)
        ax.plot(sorted_values, y, color=colors[i], label=label, linewidth=2)
    ax.set_xlabel(property_name)
    ax.set_ylabel('Cumulative Probability')
    ax.set_title('Cumulative Distribution')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # 4. Q-Q plot (comparing to normal distribution)
    ax = axes[3]
    from scipy import stats
    for i, (label, values) in enumerate(data_dict.items()):
        # Sample subset for readability
        sample_size = min(1000, len(values))
        sample = np.random.choice(values, sample_size, replace=False)
        stats.probplot(sample, dist="norm", plot=ax)
    ax.set_title('Q-Q Plot vs Normal Distribution')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Statistical comparison plot saved: {save_path}")

    return fig, axes


def plot_difference_analysis(data1: np.ndarray,
                            data2: np.ndarray,
                            label1: str = "System 1",
                            label2: str = "System 2",
                            property_name: str = "Property",
                            figsize: Tuple[float, float] = (12, 8),
                            save_path: Optional[str] = None) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Analyze and plot differences between two datasets.

    Parameters
    ----------
    data1, data2 : np.ndarray
        Data arrays to compare
    label1, label2 : str
        Labels for the datasets
    property_name : str
        Name of the property
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, list of axes) objects
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    # Ensure same length
    min_length = min(len(data1), len(data2))
    data1_trim = data1[:min_length]
    data2_trim = data2[:min_length]
    difference = data2_trim - data1_trim

    # 1. Time series comparison
    ax = axes[0]
    times = np.arange(min_length) * 0.5  # Assume 0.5 ns per frame
    ax.plot(times, data1_trim, label=label1, color='blue', alpha=0.7)
    ax.plot(times, data2_trim, label=label2, color='red', alpha=0.7)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(property_name)
    ax.set_title('Time Series Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. Difference time series
    ax = axes[1]
    ax.plot(times, difference, color='green', linewidth=1)
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax.fill_between(times, difference, 0, alpha=0.3, color='green')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(f'Δ{property_name} ({label2} - {label1})')
    ax.set_title('Difference Analysis')
    ax.grid(True, alpha=0.3)

    # 3. Distribution comparison
    ax = axes[2]
    ax.hist(data1_trim, bins=30, alpha=0.6, label=label1, color='blue', density=True)
    ax.hist(data2_trim, bins=30, alpha=0.6, label=label2, color='red', density=True)
    ax.set_xlabel(property_name)
    ax.set_ylabel('Density')
    ax.set_title('Distribution Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. Scatter plot
    ax = axes[3]
    ax.scatter(data1_trim, data2_trim, alpha=0.5, s=10)

    # Add diagonal line
    min_val = min(np.min(data1_trim), np.min(data2_trim))
    max_val = max(np.max(data1_trim), np.max(data2_trim))
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8)

    ax.set_xlabel(f'{property_name} ({label1})')
    ax.set_ylabel(f'{property_name} ({label2})')
    ax.set_title('Correlation Plot')
    ax.grid(True, alpha=0.3)

    # Add correlation coefficient
    corr_coef = np.corrcoef(data1_trim, data2_trim)[0, 1]
    ax.text(0.05, 0.95, f'r = {corr_coef:.3f}', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Difference analysis plot saved: {save_path}")

    return fig, axes