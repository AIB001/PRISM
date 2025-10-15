#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extended clustering visualization utilities for PRISM analysis module.

Additional publication-quality plots for advanced clustering analysis including
RMSD distributions, lifetime analysis, inter-cluster distances, and contact fingerprints.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import mdtraj as md
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path
import logging

# Import publication utilities
from .publication_utils import (
    PUBLICATION_FONTS, PUBLICATION_COLORS,
    apply_publication_style, get_standard_figsize
)
from ...utils.residue import format_residue_list

logger = logging.getLogger(__name__)


def plot_cluster_rmsd_distributions(topology_file: str,
                                    trajectory_file: str,
                                    clustering_results: Dict[str, Any],
                                    align_selection: str = "protein and name CA",
                                    cluster_selection: str = "protein and name CA",
                                    save_path: Optional[str] = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot RMSD distribution within each cluster.

    Shows how tight/disperse each cluster is by calculating RMSD to cluster centroid.

    Parameters
    ----------
    topology_file : str
        Path to topology file
    trajectory_file : str
        Path to trajectory file
    clustering_results : dict
        Results from ClusteringAnalyzer
    align_selection : str
        Selection for alignment
    cluster_selection : str
        Selection for RMSD calculation
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    apply_publication_style()

    labels = clustering_results['labels']
    coordinates_reduced = clustering_results['coordinates_reduced']

    # Load trajectory
    traj = md.load(trajectory_file, top=topology_file)
    align_indices = traj.topology.select(align_selection)
    rmsd_indices = traj.topology.select(cluster_selection)

    # Align trajectory
    traj.superpose(traj, frame=0, atom_indices=align_indices)

    # Get unique clusters (excluding noise)
    unique_labels = sorted(set(labels))
    if -1 in unique_labels:
        unique_labels.remove(-1)

    # Calculate RMSD to centroid for each cluster
    cluster_rmsds = {}

    for cluster_id in unique_labels:
        cluster_mask = labels == cluster_id
        cluster_frames = np.where(cluster_mask)[0]

        if len(cluster_frames) > 0:
            # Find centroid frame (frame closest to geometric center in PCA space)
            cluster_coords = coordinates_reduced[cluster_mask]
            center = np.mean(cluster_coords, axis=0)
            distances = np.linalg.norm(cluster_coords - center, axis=1)
            centroid_idx = cluster_frames[np.argmin(distances)]

            # Calculate RMSD of all cluster frames to centroid
            rmsds = md.rmsd(traj[cluster_frames], traj[centroid_idx],
                           atom_indices=rmsd_indices)
            cluster_rmsds[cluster_id] = rmsds * 10  # Convert nm to Å

    # Create violin plot
    fig, ax = plt.subplots(figsize=get_standard_figsize('distribution'))

    # Prepare data for violin plot
    plot_data = []
    for cluster_id, rmsds in cluster_rmsds.items():
        for rmsd_val in rmsds:
            plot_data.append({
                'Cluster': f'C{cluster_id}',
                'RMSD': rmsd_val
            })

    df = pd.DataFrame(plot_data)

    # Create violin plot
    colors = PUBLICATION_COLORS['example'][:len(unique_labels)]
    sns.violinplot(data=df, x='Cluster', y='RMSD', hue='Cluster',
                   palette=colors, inner='box', legend=False, ax=ax)

    # Add statistics
    for i, cluster_id in enumerate(unique_labels):
        rmsds = cluster_rmsds[cluster_id]
        median = np.median(rmsds)
        mean = np.mean(rmsds)
        print(f"  Cluster {cluster_id}: median RMSD = {median:.2f} Å, mean = {mean:.2f} Å")

    ax.set_ylabel('RMSD to Centroid (Å)')
    ax.set_xlabel('Cluster')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_facecolor('#FAFAFA')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        logger.info(f"Cluster RMSD distribution saved: {save_path}")

    return fig, ax


def plot_cluster_lifetime_distribution(clustering_results: Dict[str, Any],
                                       save_path: Optional[str] = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot distribution of cluster residence times (lifetimes).

    Analyzes how long the system stays in each cluster before transitioning.

    Parameters
    ----------
    clustering_results : dict
        Results from ClusteringAnalyzer
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    apply_publication_style()

    labels = clustering_results['labels']
    timestep_ns = clustering_results.get('timestep_ns', 0.5)

    # Get unique clusters (excluding noise)
    unique_labels = sorted(set(labels))
    if -1 in unique_labels:
        unique_labels.remove(-1)

    # Calculate lifetimes for each cluster
    cluster_lifetimes = {cluster_id: [] for cluster_id in unique_labels}

    current_cluster = labels[0]
    current_lifetime = 1

    for i in range(1, len(labels)):
        if labels[i] == current_cluster:
            current_lifetime += 1
        else:
            # Transition occurred
            if current_cluster != -1:  # Exclude noise
                cluster_lifetimes[current_cluster].append(current_lifetime * timestep_ns)
            current_cluster = labels[i]
            current_lifetime = 1

    # Add final segment
    if current_cluster != -1:
        cluster_lifetimes[current_cluster].append(current_lifetime * timestep_ns)

    # Create histogram
    fig, ax = plt.subplots(figsize=get_standard_figsize('distribution'))

    colors = PUBLICATION_COLORS['example'][:len(unique_labels)]

    for i, cluster_id in enumerate(unique_labels):
        lifetimes = cluster_lifetimes[cluster_id]
        if len(lifetimes) > 0:
            ax.hist(lifetimes, bins=20, alpha=0.6, label=f'Cluster {cluster_id}',
                   color=colors[i], edgecolor='black', linewidth=0.5)

            median_lifetime = np.median(lifetimes)
            mean_lifetime = np.mean(lifetimes)
            print(f"  Cluster {cluster_id}: median lifetime = {median_lifetime:.1f} ns, "
                  f"mean = {mean_lifetime:.1f} ns, n_visits = {len(lifetimes)}")

    ax.set_xlabel('Residence Time (ns)')
    ax.set_ylabel('Frequency')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_facecolor('#FAFAFA')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        logger.info(f"Cluster lifetime distribution saved: {save_path}")

    return fig, ax


def plot_intercluster_distance_matrix(topology_file: str,
                                      trajectory_file: str,
                                      clustering_results: Dict[str, Any],
                                      align_selection: str = "protein and name CA",
                                      cluster_selection: str = "protein and name CA",
                                      save_path: Optional[str] = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot distance matrix between cluster centroids.

    Shows conformational similarity between different clusters.

    Parameters
    ----------
    topology_file : str
        Path to topology file
    trajectory_file : str
        Path to trajectory file
    clustering_results : dict
        Results from ClusteringAnalyzer
    align_selection : str
        Selection for alignment
    cluster_selection : str
        Selection for distance calculation
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    apply_publication_style()

    labels = clustering_results['labels']
    coordinates_reduced = clustering_results['coordinates_reduced']

    # Load trajectory
    traj = md.load(trajectory_file, top=topology_file)
    align_indices = traj.topology.select(align_selection)
    rmsd_indices = traj.topology.select(cluster_selection)

    # Align trajectory
    traj.superpose(traj, frame=0, atom_indices=align_indices)

    # Get unique clusters (excluding noise)
    unique_labels = sorted(set(labels))
    if -1 in unique_labels:
        unique_labels.remove(-1)

    n_clusters = len(unique_labels)

    # Find centroid frames
    centroid_frames = []
    for cluster_id in unique_labels:
        cluster_mask = labels == cluster_id
        cluster_frames = np.where(cluster_mask)[0]
        cluster_coords = coordinates_reduced[cluster_mask]
        center = np.mean(cluster_coords, axis=0)
        distances = np.linalg.norm(cluster_coords - center, axis=1)
        centroid_idx = cluster_frames[np.argmin(distances)]
        centroid_frames.append(centroid_idx)

    # Calculate pairwise RMSD between centroids
    distance_matrix = np.zeros((n_clusters, n_clusters))

    for i in range(n_clusters):
        for j in range(n_clusters):
            rmsd_val = md.rmsd(traj[centroid_frames[i]], traj[centroid_frames[j]],
                              atom_indices=rmsd_indices)[0]
            distance_matrix[i, j] = rmsd_val * 10  # Convert nm to Å

    # Create heatmap
    fig, ax = plt.subplots(figsize=get_standard_figsize('single'))

    im = ax.imshow(distance_matrix, cmap='YlOrRd', aspect='auto')

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('RMSD (Å)', rotation=270, labelpad=20,
                  fontsize=PUBLICATION_FONTS['colorbar'], fontweight='bold')
    cbar.ax.tick_params(labelsize=PUBLICATION_FONTS['tick_label'])

    # Add text annotations
    for i in range(n_clusters):
        for j in range(n_clusters):
            value = distance_matrix[i, j]
            text_color = 'white' if value > distance_matrix.max()/2 else 'black'
            ax.text(j, i, f'{value:.2f}', ha='center', va='center',
                   color=text_color, fontsize=PUBLICATION_FONTS['value_text'],
                   fontweight='bold')

    # Styling
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Cluster')
    ax.set_xticks(range(n_clusters))
    ax.set_yticks(range(n_clusters))
    ax.set_xticklabels([f'C{i}' for i in unique_labels])
    ax.set_yticklabels([f'C{i}' for i in unique_labels])

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        logger.info(f"Inter-cluster distance matrix saved: {save_path}")

    return fig, ax


def plot_cluster_contact_fingerprints(topology_file: str,
                                      trajectory_file: str,
                                      clustering_results: Dict[str, Any],
                                      key_residues: List[str],
                                      cutoff: float = 4.0,
                                      residue_format: str = "1letter",
                                      save_path: Optional[str] = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot contact probability fingerprints for each cluster.

    Shows how protein-ligand contacts differ between conformational clusters.

    Parameters
    ----------
    topology_file : str
        Path to topology file
    trajectory_file : str
        Path to trajectory file
    clustering_results : dict
        Results from ClusteringAnalyzer
    key_residues : list
        List of key residue names (3-letter code + number, e.g., ['ASP618', 'ASN691'])
    cutoff : float
        Distance cutoff in Angstroms
    residue_format : str
        Display format: "1letter" or "3letter"
    save_path : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    apply_publication_style()
    from ...analysis.calc.contacts import ContactAnalyzer
    from ...analysis.core.config import AnalysisConfig

    labels = clustering_results['labels']

    # Get unique clusters (excluding noise)
    unique_labels = sorted(set(labels))
    if -1 in unique_labels:
        unique_labels.remove(-1)

    # Calculate contacts for each cluster
    config = AnalysisConfig()
    analyzer = ContactAnalyzer(config)

    cluster_contacts = {}

    for cluster_id in unique_labels:
        cluster_mask = labels == cluster_id
        cluster_frames = np.where(cluster_mask)[0]

        # Load only cluster frames
        traj = md.load(trajectory_file, top=topology_file)
        cluster_traj = traj[cluster_frames]

        # Save temporary trajectory
        temp_traj = Path(f"/tmp/cluster_{cluster_id}_temp.xtc")
        cluster_traj.save_xtc(str(temp_traj))

        # Analyze contacts
        try:
            contacts = analyzer.analyze_key_residue_contacts(
                universe=topology_file,
                trajectory=str(temp_traj),
                key_residues=key_residues,
                cutoff=cutoff,
                step=1
            )
            cluster_contacts[cluster_id] = contacts

            # Clean up temp file
            temp_traj.unlink()

        except Exception as e:
            logger.warning(f"Contact analysis failed for cluster {cluster_id}: {e}")
            continue

    if not cluster_contacts:
        logger.error("No contact data available")
        return None, None

    # Create grouped bar plot
    fig, ax = plt.subplots(figsize=get_standard_figsize('wide'))

    residues = list(key_residues)
    n_residues = len(residues)
    n_clusters = len(cluster_contacts)

    x = np.arange(n_residues)
    width = 0.8 / n_clusters
    colors = PUBLICATION_COLORS['example'][:n_clusters]

    for i, cluster_id in enumerate(sorted(cluster_contacts.keys())):
        contacts = cluster_contacts[cluster_id]
        values = [contacts.get(res, 0.0) for res in residues]

        offset = (i - n_clusters/2 + 0.5) * width
        bars = ax.bar(x + offset, values, width, label=f'Cluster {cluster_id}',
                     color=colors[i], alpha=0.9, edgecolor='black', linewidth=0.5)

        # Print summary
        print(f"  Cluster {cluster_id} contact summary:")
        for res, val in zip(residues, values):
            if val > 50:  # Only print significant contacts
                print(f"    {res}: {val:.1f}%")

    # Format residue names
    residues_display = format_residue_list(residues, residue_format)

    # Styling
    ax.set_ylabel('Contact Probability (%)')
    ax.set_xlabel('Amino Acid Residues')
    ax.set_xticks(x)
    ax.set_xticklabels(residues_display, rotation=45, ha='center')
    ax.set_ylim(0, 110)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_facecolor('#FAFAFA')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        logger.info(f"Cluster contact fingerprints saved: {save_path}")

    return fig, ax


def save_cluster_representative_structures(topology_file: str,
                                          trajectory_file: str,
                                          clustering_results: Dict[str, Any],
                                          output_dir: str,
                                          align_selection: str = "protein and name CA"):
    """
    Save representative structures (centroids) for each cluster.

    Parameters
    ----------
    topology_file : str
        Path to topology file
    trajectory_file : str
        Path to trajectory file
    clustering_results : dict
        Results from ClusteringAnalyzer
    output_dir : str
        Directory to save centroid PDB files
    align_selection : str
        Selection for alignment

    Returns
    -------
    dict
        Mapping of cluster_id -> saved_file_path
    """
    labels = clustering_results['labels']
    coordinates_reduced = clustering_results['coordinates_reduced']

    # Load trajectory
    traj = md.load(trajectory_file, top=topology_file)
    align_indices = traj.topology.select(align_selection)
    traj.superpose(traj, frame=0, atom_indices=align_indices)

    # Get unique clusters (excluding noise)
    unique_labels = sorted(set(labels))
    if -1 in unique_labels:
        unique_labels.remove(-1)

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    for cluster_id in unique_labels:
        cluster_mask = labels == cluster_id
        cluster_frames = np.where(cluster_mask)[0]
        cluster_coords = coordinates_reduced[cluster_mask]

        # Find centroid frame
        center = np.mean(cluster_coords, axis=0)
        distances = np.linalg.norm(cluster_coords - center, axis=1)
        centroid_idx = cluster_frames[np.argmin(distances)]

        # Save centroid structure
        filename = f"centroid_cluster_{cluster_id}_frame_{centroid_idx}.pdb"
        filepath = output_path / filename
        traj[centroid_idx].save(str(filepath))

        saved_files[cluster_id] = str(filepath)
        logger.info(f"Saved centroid for cluster {cluster_id}: {filename}")

    return saved_files
