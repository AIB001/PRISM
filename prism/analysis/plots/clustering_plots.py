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


        # Print statistics (removed from plot to avoid hiding content)
        n_clusters = len(unique_labels)
        if 'silhouette_score' in clustering_results:
            silhouette = clustering_results['silhouette_score']
            print(f"  📊 Clustering statistics: {n_clusters} clusters, Silhouette score = {silhouette:.3f}")
        else:
            print(f"  📊 Clustering statistics: {n_clusters} clusters")

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

    # Debug logging
    logger.debug(f"Timeline plot: {len(labels)} labels, {len(frame_indices)} frame_indices")

    # Ensure frame_indices matches labels length
    # DO NOT regenerate - this would lose the actual frame numbers!
    # If mismatch occurs, it's a bug in clustering.py that should be fixed there
    if len(frame_indices) != len(labels):
        logger.error(f"Frame indices length mismatch: {len(frame_indices)} != {len(labels)}")
        logger.error(f"This is a bug in clustering analysis - frame_indices should match labels length")
        # Use what we have rather than creating wrong sequential indices
        if len(frame_indices) < len(labels):
            # Pad with extrapolated values if needed
            step = frame_indices[1] - frame_indices[0] if len(frame_indices) > 1 else 1
            while len(frame_indices) < len(labels):
                frame_indices.append(frame_indices[-1] + step)
        else:
            # Truncate if too long
            frame_indices = frame_indices[:len(labels)]

    if times is None:
        # Get timestep from results metadata, fallback to 0.5 ns for backward compatibility
        timestep_ns = clustering_results.get('timestep_ns', 0.5)
        times = np.array(frame_indices) * timestep_ns
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

        # Plot timeline as thick horizontal bar segments
        # Group consecutive frames with same cluster into continuous segments
        segments = []

        if len(labels) == 0:
            # No data to plot
            pass
        elif len(labels) == 1:
            # Single frame - create a single segment with some width
            time_width = 0.05 if len(times) == 1 else (times[-1] - times[0]) * 0.1
            segments.append({
                'cluster': labels[0],
                't_start': times[0] - time_width/2,
                't_end': times[0] + time_width/2,
                'start_idx': 0,
                'end_idx': 0
            })
        else:
            # Multiple frames - group into segments
            start_idx = 0

            for i in range(1, len(labels)):
                if labels[i] != labels[i-1]:
                    # Cluster changed, save previous segment
                    segments.append({
                        'cluster': labels[start_idx],
                        't_start': times[start_idx],
                        't_end': times[i-1],
                        'start_idx': start_idx,
                        'end_idx': i-1
                    })
                    start_idx = i

            # Add final segment
            segments.append({
                'cluster': labels[start_idx],
                't_start': times[start_idx],
                't_end': times[-1],
                'start_idx': start_idx,
                'end_idx': len(labels)-1
            })

        # Draw thick horizontal lines for each segment
        for seg in segments:
            cluster = seg['cluster']
            color = color_map[cluster]

            # Use thick horizontal lines to show cluster occupancy over time
            ax.hlines(y=cluster, xmin=seg['t_start'], xmax=seg['t_end'],
                     colors=color, linewidth=10, alpha=0.85, zorder=3)

            # Add small markers at segment boundaries for clarity
            if seg['start_idx'] == seg['end_idx']:
                # Single point segment - just add one marker
                ax.scatter([seg['t_start'] + (seg['t_end'] - seg['t_start'])/2], [cluster],
                          c=[color], s=60, alpha=0.9,
                          edgecolor='black', linewidth=0.5, zorder=4)
            else:
                # Multi-point segment - add markers at both ends
                ax.scatter([seg['t_start'], seg['t_end']], [cluster, cluster],
                          c=[color, color], s=40, alpha=0.9,
                          edgecolor='black', linewidth=0.5, zorder=4)

        # Add trajectory boundaries if combined trajectories
        if 'trajectory_metadata' in clustering_results:
            metadata = clustering_results['trajectory_metadata']
            frames_per_traj = metadata.get('frames_per_trajectory', [])

            if len(frames_per_traj) > 1:
                # Calculate boundary positions
                cumulative_frames = 0
                for i, n_frames in enumerate(frames_per_traj[:-1]):  # Don't add line after last trajectory
                    cumulative_frames += n_frames
                    boundary_time = times[cumulative_frames] if cumulative_frames < len(times) else times[-1]
                    ax.axvline(x=boundary_time, color='red', linestyle='--',
                              alpha=0.5, linewidth=2, zorder=1,
                              label='Trajectory boundary' if i == 0 else '')

        # Add horizontal lines for each cluster
        for cluster in unique_labels:
            if cluster != -1:
                ax.axhline(y=cluster, color='gray', linestyle='--', alpha=0.2, linewidth=1, zorder=1)

        # Calculate cluster populations
        cluster_stats = {}
        for cluster in unique_labels:
            count = np.sum(labels == cluster)
            percentage = count / len(labels) * 100
            cluster_stats[cluster] = {'count': count, 'percentage': percentage}

        # Print population statistics (removed from plot to avoid hiding content)
        print("  📊 Cluster populations:")
        for cluster in sorted(unique_labels):
            if cluster == -1:
                print(f"    Noise: {cluster_stats[cluster]['percentage']:.1f}% ({cluster_stats[cluster]['count']} frames)")
            else:
                print(f"    Cluster {cluster}: {cluster_stats[cluster]['percentage']:.1f}% ({cluster_stats[cluster]['count']} frames)")

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

        # Add legend if trajectory boundaries are shown
        if 'trajectory_metadata' in clustering_results:
            metadata = clustering_results['trajectory_metadata']
            if len(metadata.get('frames_per_trajectory', [])) > 1:
                ax.legend(loc='upper right', fontsize=PUBLICATION_FONTS['legend'])

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

        # Print RMSD statistics (removed from plot to avoid hiding content)
        print(f"  📊 RMSD matrix statistics: Mean = {mean_rmsd:.2f} Å, Std = {std_rmsd:.2f} Å, Max = {max_rmsd:.2f} Å")

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


def plot_cluster_transition_matrix(clustering_results: Dict[str, Any],
                                   title: str = "",
                                   figsize: Optional[Tuple[float, float]] = None,
                                   save_path: Optional[str] = None,
                                   normalize: bool = True) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot cluster transition probability matrix.

    Parameters
    ----------
    clustering_results : dict
        Results from ClusteringAnalyzer containing labels
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    normalize : bool
        If True, normalize rows to show transition probabilities

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    labels = clustering_results['labels']
    n_clusters = clustering_results.get('n_clusters', len(set(labels)))

    if figsize is None:
        figsize = get_standard_figsize("single")

    plt.style.use('default')
    with plt.rc_context(get_publication_style()):
        fig, ax = plt.subplots(figsize=figsize)

        # Calculate transition matrix
        transition_matrix = np.zeros((n_clusters, n_clusters))

        for t in range(len(labels) - 1):
            from_cluster = labels[t]
            to_cluster = labels[t + 1]
            if from_cluster >= 0 and to_cluster >= 0:  # Exclude noise (-1) for DBSCAN
                transition_matrix[from_cluster, to_cluster] += 1

        # Normalize if requested
        if normalize:
            row_sums = transition_matrix.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1  # Avoid division by zero
            transition_matrix = transition_matrix / row_sums

        # Create heatmap
        im = ax.imshow(transition_matrix, cmap='Blues', vmin=0, vmax=1 if normalize else None,
                      aspect='auto')

        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar_label = 'Transition Probability' if normalize else 'Transition Count'
        cbar.set_label(cbar_label, rotation=270, labelpad=25,
                      fontsize=PUBLICATION_FONTS['colorbar'], weight='bold')
        cbar.ax.tick_params(labelsize=PUBLICATION_FONTS['tick_label'])

        # Add text annotations
        for i in range(n_clusters):
            for j in range(n_clusters):
                value = transition_matrix[i, j]
                if normalize:
                    text = f'{value:.2f}'
                else:
                    text = f'{int(value)}'

                # Choose text color based on background
                text_color = 'white' if value > (0.5 if normalize else transition_matrix.max()/2) else 'black'
                ax.text(j, i, text, ha='center', va='center',
                       color=text_color, fontsize=PUBLICATION_FONTS['value_text'])

        # Styling
        ax.set_xlabel('To Cluster', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax.set_ylabel('From Cluster', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

        ax.set_xticks(range(n_clusters))
        ax.set_yticks(range(n_clusters))
        ax.set_xticklabels([f'C{i}' for i in range(n_clusters)],
                          fontsize=PUBLICATION_FONTS['tick_label'])
        ax.set_yticklabels([f'C{i}' for i in range(n_clusters)],
                          fontsize=PUBLICATION_FONTS['tick_label'])

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            logger.info(f"Cluster transition matrix saved: {save_path}")

    return fig, ax