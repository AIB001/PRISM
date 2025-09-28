#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clustering visualization utilities for PRISM analysis module.

Provides publication-quality plots for trajectory clustering analysis including
PCA scatter plots, cluster timelines, population distributions, and RMSD matrices.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any
import logging

# Import publication utilities
from .publication_utils import (
    get_publication_style, get_color_palette, fix_rotated_labels,
    add_statistical_annotations, style_axes_for_publication,
    save_publication_figure, validate_figure_size, PUBLICATION_FONTS, PUBLICATION_COLORS,
    apply_publication_style, setup_publication_figure, get_standard_figsize
)

# Apply global publication style on import
apply_publication_style()

logger = logging.getLogger(__name__)


def plot_cluster_scatter(clustering_results: Dict[str, Any],
                        title: str = "",
                        figsize: Optional[Tuple[float, float]] = None,
                        save_path: Optional[str] = None,
                        show_centers: bool = True,
                        show_hulls: bool = False,
                        alpha: float = 0.7) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create 2D scatter plot of clusters in PCA space.

    Parameters
    ----------
    clustering_results : dict
        Results from ClusteringAnalyzer containing labels, coordinates_reduced, etc.
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    show_centers : bool
        Whether to show cluster centers
    show_hulls : bool
        Whether to show convex hulls around clusters
    alpha : float
        Point transparency

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if 'coordinates_reduced' not in clustering_results:
        raise ValueError("coordinates_reduced not found in clustering results")

    coordinates = clustering_results['coordinates_reduced']
    labels = clustering_results['labels']

    if figsize is None:
        figsize = get_standard_figsize("single")

    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        fig, ax = plt.subplots(figsize=figsize)

        # Get unique cluster labels (excluding noise if present)
        unique_labels = set(labels)
        if -1 in unique_labels:  # Remove noise label for DBSCAN
            unique_labels.remove(-1)

        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_labels)))

        # Plot noise points first (if any)
        if -1 in labels:
            noise_mask = labels == -1
            ax.scatter(coordinates[noise_mask, 0], coordinates[noise_mask, 1],
                      c='gray', alpha=alpha*0.5, s=30, label='Noise',
                      edgecolor='black', linewidth=0.5)

        # Plot clusters
        for i, (cluster_id, color) in enumerate(zip(sorted(unique_labels), colors)):
            cluster_mask = labels == cluster_id
            cluster_coords = coordinates[cluster_mask]

            ax.scatter(cluster_coords[:, 0], cluster_coords[:, 1],
                      c=[color], alpha=alpha, s=50,
                      label=f'Cluster {cluster_id}',
                      edgecolor='black', linewidth=0.5)

            # Show cluster centers
            if show_centers and 'cluster_centers' in clustering_results:
                centers = clustering_results['cluster_centers']
                if len(centers) > cluster_id:
                    ax.scatter(centers[cluster_id, 0], centers[cluster_id, 1],
                              c='red', marker='x', s=200, linewidths=3,
                              label=f'Center {cluster_id}' if i == 0 else "")

            # Show convex hulls
            if show_hulls and len(cluster_coords) > 2:
                try:
                    from scipy.spatial import ConvexHull
                    hull = ConvexHull(cluster_coords[:, :2])
                    for simplex in hull.simplices:
                        ax.plot(cluster_coords[simplex, 0], cluster_coords[simplex, 1],
                               'k--', alpha=0.3, linewidth=1)
                except ImportError:
                    logger.warning("scipy not available for convex hulls")

        # Styling
        ax.set_xlabel('Principal Component 1', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax.set_ylabel('Principal Component 2', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')


        # Add statistics
        n_clusters = len(unique_labels)
        if 'silhouette_score' in clustering_results:
            silhouette = clustering_results['silhouette_score']
            stats_text = f"Clusters: {n_clusters}\nSilhouette: {silhouette:.3f}"
        else:
            stats_text = f"Clusters: {n_clusters}"

        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=PUBLICATION_FONTS['annotation'],
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))

        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=PUBLICATION_FONTS['legend'])

        # Clean spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"Cluster scatter plot saved: {save_path}")

    return fig, ax


def plot_cluster_timeline(clustering_results: Dict[str, Any],
                         times: Optional[np.ndarray] = None,
                         title: str = "",
                         figsize: Optional[Tuple[float, float]] = None,
                         save_path: Optional[str] = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot cluster assignment over simulation time.

    Parameters
    ----------
    clustering_results : dict
        Results from ClusteringAnalyzer
    times : np.ndarray, optional
        Time points. If None, use frame indices
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    labels = clustering_results['labels']
    frame_indices = clustering_results.get('frame_indices', list(range(len(labels))))

    if times is None:
        times = np.array(frame_indices) * 0.02  # Assume 20 ps timesteps
    else:
        times = times[:len(labels)]

    if figsize is None:
        figsize = get_standard_figsize("single")

    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        fig, ax = plt.subplots(figsize=figsize)

        # Create color map for clusters
        unique_labels = sorted(set(labels))
        if -1 in unique_labels:
            unique_labels.remove(-1)
            unique_labels.append(-1)  # Put noise at end

        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_labels)))
        color_map = {label: colors[i] for i, label in enumerate(unique_labels)}
        if -1 in color_map:
            color_map[-1] = 'gray'  # Special color for noise

        # Plot timeline
        for i in range(len(times)):
            cluster = labels[i]
            color = color_map[cluster]
            ax.scatter(times[i], cluster, c=[color], s=20, alpha=0.8,
                      edgecolor='black', linewidth=0.3)

        # Add horizontal lines for each cluster
        for cluster in unique_labels:
            if cluster != -1:
                ax.axhline(y=cluster, color='gray', linestyle='--', alpha=0.3, linewidth=1)

        # Calculate cluster populations
        cluster_stats = {}
        for cluster in unique_labels:
            count = np.sum(labels == cluster)
            percentage = count / len(labels) * 100
            cluster_stats[cluster] = {'count': count, 'percentage': percentage}

        # Add population statistics
        stats_text = "Cluster populations:\n"
        for cluster in sorted(unique_labels):
            if cluster == -1:
                stats_text += f"Noise: {cluster_stats[cluster]['percentage']:.1f}%\n"
            else:
                stats_text += f"Cluster {cluster}: {cluster_stats[cluster]['percentage']:.1f}%\n"

        ax.text(0.02, 0.98, stats_text.strip(), transform=ax.transAxes,
                verticalalignment='top', fontsize=PUBLICATION_FONTS['annotation'],
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))

        # Styling
        ax.set_xlabel('Time (ns)', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax.set_ylabel('Cluster ID', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

        ax.grid(True, alpha=0.3, axis='x')

        # Set y-axis to show all clusters
        if unique_labels:
            y_min = min(unique_labels) - 0.5
            y_max = max(unique_labels) + 0.5
            ax.set_ylim(y_min, y_max)
            ax.set_yticks(unique_labels)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"Cluster timeline plot saved: {save_path}")

    return fig, ax


def plot_cluster_populations(clustering_results: Dict[str, Any],
                           title: str = "",
                           figsize: Optional[Tuple[float, float]] = None,
                           save_path: Optional[str] = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot cluster population distribution as bar chart.

    Parameters
    ----------
    clustering_results : dict
        Results from ClusteringAnalyzer
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    labels = clustering_results['labels']

    if figsize is None:
        figsize = get_standard_figsize("single")

    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        fig, ax = plt.subplots(figsize=figsize)

        # Calculate populations
        unique_labels, counts = np.unique(labels, return_counts=True)
        percentages = counts / len(labels) * 100

        # Prepare colors
        colors = []
        cluster_names = []
        for label in unique_labels:
            if label == -1:
                colors.append('gray')
                cluster_names.append('Noise')
            else:
                colors.append(plt.cm.Set3(label / max(1, max(unique_labels))))
                cluster_names.append(f'Cluster {label}')

        # Create bar chart
        bars = ax.bar(range(len(unique_labels)), percentages, color=colors,
                     alpha=0.8, edgecolor='black', linewidth=1.2)

        # Add percentage labels on bars
        for i, (bar, percentage) in enumerate(zip(bars, percentages)):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                   f'{percentage:.1f}%', ha='center', va='bottom',
                   fontsize=PUBLICATION_FONTS['value_text'], weight='bold')

        # Styling
        ax.set_xlabel('Clusters', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax.set_ylabel('Population (%)', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

        ax.set_xticks(range(len(unique_labels)))
        ax.set_xticklabels(cluster_names, fontsize=PUBLICATION_FONTS['tick_label'])
        ax.set_ylim(0, max(percentages) * 1.15)  # Add space for labels
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_axisbelow(True)

        # Clean spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"Cluster populations plot saved: {save_path}")

    return fig, ax


def plot_rmsd_matrix(rmsd_matrix: np.ndarray,
                    title: str = "",
                    figsize: Optional[Tuple[float, float]] = None,
                    save_path: Optional[str] = None,
                    cmap: str = 'viridis') -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot RMSD matrix as heatmap.

    Parameters
    ----------
    rmsd_matrix : np.ndarray
        Pairwise RMSD matrix
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    cmap : str
        Colormap name

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("single")

    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        fig, ax = plt.subplots(figsize=figsize)

        # Create heatmap
        im = ax.imshow(rmsd_matrix, cmap=cmap, aspect='auto')

        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label('RMSD (Å)', rotation=270, labelpad=25,
                      fontsize=PUBLICATION_FONTS['colorbar'], weight='bold')
        cbar.ax.tick_params(labelsize=PUBLICATION_FONTS['tick_label'])

        # Styling
        ax.set_xlabel('Frame Index', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax.set_ylabel('Frame Index', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')


        # Add statistics
        mean_rmsd = np.mean(rmsd_matrix)
        std_rmsd = np.std(rmsd_matrix)
        max_rmsd = np.max(rmsd_matrix)

        stats_text = f"Mean RMSD: {mean_rmsd:.2f} Å\nStd RMSD: {std_rmsd:.2f} Å\nMax RMSD: {max_rmsd:.2f} Å"
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=PUBLICATION_FONTS['annotation'],
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"RMSD matrix plot saved: {save_path}")

    return fig, ax


def plot_cluster_optimization(optimization_results: Dict[str, Any],
                            title: str = "",
                            figsize: Optional[Tuple[float, float]] = None,
                            save_path: Optional[str] = None) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Plot cluster optimization metrics (elbow method and silhouette scores).

    Parameters
    ----------
    optimization_results : dict
        Results from find_optimal_clusters
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, (inertia_ax, silhouette_ax)) objects
    """
    if figsize is None:
        figsize = get_standard_figsize("horizontal")

    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        cluster_range = optimization_results['cluster_range']
        silhouette_scores = optimization_results['silhouette_scores']
        inertias = optimization_results.get('inertias')
        optimal_clusters = optimization_results['optimal_clusters']

        # Plot 1: Elbow method (if inertias available)
        if inertias:
            ax1.plot(cluster_range, inertias, 'bo-', linewidth=2, markersize=8)
            ax1.axvline(x=optimal_clusters, color='red', linestyle='--', linewidth=2,
                       label=f'Optimal: {optimal_clusters}')
            ax1.set_xlabel('Number of Clusters', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
            ax1.set_ylabel('Inertia', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

            ax1.grid(True, alpha=0.3)
            ax1.legend(fontsize=PUBLICATION_FONTS['legend'])
        else:
            ax1.text(0.5, 0.5, 'Inertia not available\n(method-specific)',
                    transform=ax1.transAxes, ha='center', va='center',
                    fontsize=PUBLICATION_FONTS['annotation'])

        # Plot 2: Silhouette scores
        ax2.plot(cluster_range, silhouette_scores, 'ro-', linewidth=2, markersize=8)
        ax2.axvline(x=optimal_clusters, color='red', linestyle='--', linewidth=2,
                   label=f'Optimal: {optimal_clusters}')
        ax2.axhline(y=optimization_results['optimal_silhouette'], color='green',
                   linestyle=':', alpha=0.7, label=f'Max Score: {optimization_results["optimal_silhouette"]:.3f}')
        ax2.set_xlabel('Number of Clusters', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax2.set_ylabel('Silhouette Score', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=PUBLICATION_FONTS['legend'])

        # Clean spines for both plots
        for ax in [ax1, ax2]:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(1.2)
            ax.spines['bottom'].set_linewidth(1.2)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"Cluster optimization plot saved: {save_path}")

    return fig, (ax1, ax2)


def plot_multi_trajectory_clustering(clustering_results_dict: Dict[str, Dict[str, Any]],
                                   title: str = "",
                                   figsize: Optional[Tuple[float, float]] = None,
                                   save_path: Optional[str] = None) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Plot clustering results for multiple trajectories in a grid.

    Parameters
    ----------
    clustering_results_dict : dict
        Dictionary mapping trajectory names to clustering results
    title : str
        Overall plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes_list) objects
    """
    n_trajs = len(clustering_results_dict)
    cols = min(3, n_trajs)
    rows = (n_trajs + cols - 1) // cols

    if figsize is None:
        # 使用标准尺寸计算多panel尺寸
        base_width, base_height = get_standard_figsize('single')
        figsize = (base_width * cols * 0.8, base_height * rows * 0.8)

    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        fig, axes = plt.subplots(rows, cols, figsize=figsize)
        if n_trajs == 1:
            axes = [axes]
        elif rows == 1:
            axes = list(axes) if hasattr(axes, '__iter__') else [axes]
        else:
            axes = axes.flatten()

        for i, (traj_name, results) in enumerate(clustering_results_dict.items()):
            ax = axes[i]

            if 'coordinates_reduced' in results:
                coordinates = results['coordinates_reduced']
                labels = results['labels']

                # Create scatter plot
                unique_labels = set(labels)
                if -1 in unique_labels:
                    unique_labels.remove(-1)

                colors = plt.cm.Set3(np.linspace(0, 1, len(unique_labels)))

                # Plot noise points
                if -1 in labels:
                    noise_mask = labels == -1
                    ax.scatter(coordinates[noise_mask, 0], coordinates[noise_mask, 1],
                              c='gray', alpha=0.5, s=20, label='Noise')

                # Plot clusters
                for cluster_id, color in zip(sorted(unique_labels), colors):
                    cluster_mask = labels == cluster_id
                    cluster_coords = coordinates[cluster_mask]
                    ax.scatter(cluster_coords[:, 0], cluster_coords[:, 1],
                              c=[color], alpha=0.7, s=30, label=f'C{cluster_id}')

                # Show centers if available
                if 'cluster_centers' in results:
                    centers = results['cluster_centers']
                    ax.scatter(centers[:, 0], centers[:, 1], c='red',
                              marker='x', s=100, linewidths=2)


                ax.set_xlabel('PC1', fontsize=PUBLICATION_FONTS['axis_label'])
                ax.set_ylabel('PC2', fontsize=PUBLICATION_FONTS['axis_label'])
                ax.grid(True, alpha=0.3)

                # Add silhouette score
                if 'silhouette_score' in results:
                    silhouette = results['silhouette_score']
                    ax.text(0.02, 0.98, f'Silhouette: {silhouette:.3f}',
                           transform=ax.transAxes, verticalalignment='top',
                           fontsize=PUBLICATION_FONTS['annotation'],
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Hide unused subplots
        for i in range(n_trajs, len(axes)):
            axes[i].set_visible(False)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"Multi-trajectory clustering plot saved: {save_path}")

    return fig, axes