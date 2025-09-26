#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clustering analysis for MD trajectories using pure MDTraj implementation.

Completely rewritten to use MDTraj instead of MDAnalysis, following PRISM architecture requirements.
Implements GROMOS-style coordinate-based clustering and PCA-based methods.
"""

import numpy as np
import mdtraj as md
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import os
from pathlib import Path
from typing import List, Union, Optional, Tuple, Dict
import logging
import pickle

from ..config import AnalysisConfig

logger = logging.getLogger(__name__)


class ClusteringAnalyzer:
    """Pure MDTraj clustering analysis for protein-ligand trajectories"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path("./cache")
        self._cache_dir.mkdir(parents=True, exist_ok=True)

    def _load_and_align_trajectory(self,
                                  topology: str,
                                  trajectory: Union[str, List[str]],
                                  align_selection: str = "protein and name CA",
                                  start_frame: int = 0,
                                  end_frame: Optional[int] = None,
                                  step: int = 1) -> md.Trajectory:
        """
        Load and align trajectory using MDTraj.

        Parameters
        ----------
        topology : str
            Path to topology file
        trajectory : str or list of str
            Path(s) to trajectory file(s)
        align_selection : str
            MDTraj selection for alignment
        start_frame : int
            Starting frame
        end_frame : int, optional
            Ending frame
        step : int
            Step size

        Returns
        -------
        md.Trajectory
            Aligned trajectory
        """
        # Load trajectory - handle frame slicing separately
        # For large files, limit frames to avoid memory issues
        traj = md.load(trajectory, top=topology)

        # Apply frame slicing after loading
        if end_frame is None:
            end_frame = traj.n_frames

        # In debug/test mode, limit to max 200 frames for performance
        if step >= 30:  # Debug mode indicator
            max_frames = min(200, (end_frame - start_frame) // step)
            end_frame = min(end_frame, start_frame + max_frames * step)

        traj = traj[start_frame:end_frame:step]

        logger.info(f"Selected {traj.n_frames} frames from {start_frame}:{end_frame}:{step}")

        # Parse alignment selection and align trajectory
        if "segid" in align_selection:
            # Handle segid-based selection - convert to chainid for MDTraj
            # "segid PR1 and name CA" -> find main protein chain and use CA atoms
            align_indices = []

            # Find the largest protein chain (main protein chain)
            largest_chain = None
            max_protein_atoms = 0

            for chain in traj.topology.chains:
                try:
                    # Count protein atoms in this chain
                    chain_protein = traj.topology.select(f'chainid {chain.index} and protein')
                    if len(chain_protein) > max_protein_atoms:
                        max_protein_atoms = len(chain_protein)
                        largest_chain = chain
                except:
                    continue

            if largest_chain is not None and "name CA" in align_selection:
                # Select CA atoms from the main protein chain
                try:
                    align_indices = traj.topology.select(f'chainid {largest_chain.index} and name CA')
                except:
                    # Fallback: manual selection
                    ca_atoms = [atom.index for atom in largest_chain.atoms if atom.name == 'CA' and atom.residue.is_protein]
                    align_indices = np.array(ca_atoms)
            else:
                # Fallback to general protein selection
                align_indices = traj.topology.select("protein and name CA")

            align_indices = np.array(align_indices) if not isinstance(align_indices, np.ndarray) else align_indices
        else:
            # Standard MDTraj selection
            align_indices = traj.topology.select(align_selection)

        if len(align_indices) == 0:
            raise ValueError(f"No atoms found for alignment selection: {align_selection}")

        # Align trajectory to first frame
        traj.superpose(traj, frame=0, atom_indices=align_indices)

        logger.info(f"Loaded and aligned trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")
        logger.info(f"Alignment atoms: {len(align_indices)} atoms")

        return traj

    def _extract_coordinates(self,
                           traj: md.Trajectory,
                           cluster_selection: str = "protein") -> np.ndarray:
        """
        Extract coordinates for clustering using MDTraj.

        Parameters
        ----------
        traj : md.Trajectory
            Aligned trajectory
        cluster_selection : str
            Selection for clustering atoms

        Returns
        -------
        np.ndarray
            Coordinate matrix (n_frames, n_atoms * 3)
        """
        # Parse cluster selection
        if "segid" in cluster_selection:
            # Handle segid-based selection - find main protein chain
            cluster_indices = []

            # Find the largest protein chain (main protein chain)
            largest_chain = None
            max_protein_atoms = 0

            for chain in traj.topology.chains:
                try:
                    # Count protein atoms in this chain
                    chain_protein = traj.topology.select(f'chainid {chain.index} and protein')
                    if len(chain_protein) > max_protein_atoms:
                        max_protein_atoms = len(chain_protein)
                        largest_chain = chain
                except:
                    continue

            if largest_chain is not None:
                if "name CA" in cluster_selection:
                    # Only CA atoms from main protein chain
                    try:
                        cluster_indices = traj.topology.select(f'chainid {largest_chain.index} and name CA')
                    except:
                        ca_atoms = [atom.index for atom in largest_chain.atoms if atom.name == 'CA' and atom.residue.is_protein]
                        cluster_indices = np.array(ca_atoms)
                else:
                    # All atoms in the main protein chain
                    try:
                        cluster_indices = traj.topology.select(f'chainid {largest_chain.index} and protein')
                    except:
                        protein_atoms = [atom.index for atom in largest_chain.atoms if atom.residue.is_protein]
                        cluster_indices = np.array(protein_atoms)

            cluster_indices = np.array(cluster_indices) if not isinstance(cluster_indices, np.ndarray) else cluster_indices
        else:
            # Standard MDTraj selection
            if cluster_selection == "protein":
                cluster_indices = traj.topology.select("protein")
            else:
                cluster_indices = traj.topology.select(cluster_selection)

        if len(cluster_indices) == 0:
            raise ValueError(f"No atoms found for cluster selection: {cluster_selection}")

        # Extract coordinates
        cluster_coords = traj.xyz[:, cluster_indices, :]  # (n_frames, n_atoms, 3)
        coordinates = cluster_coords.reshape(traj.n_frames, -1)  # (n_frames, n_atoms * 3)

        logger.info(f"Extracted coordinates: {coordinates.shape[0]} frames, {len(cluster_indices)} atoms")
        return coordinates

    def perform_kmeans_clustering(self,
                                 universe: str,  # Now expects topology path
                                 trajectory: Union[str, List[str]],
                                 n_clusters: int = 5,
                                 align_selection: str = "protein and name CA",
                                 cluster_selection: str = "protein",
                                 start_frame: int = 0,
                                 end_frame: Optional[int] = None,
                                 step: int = 1,
                                 use_pca: bool = True,
                                 n_components: int = 10,
                                 cache_name: Optional[str] = None) -> Dict:
        """
        Perform K-means clustering on trajectory using pure MDTraj.

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str or list of str
            Path(s) to trajectory file(s)
        n_clusters : int
            Number of clusters
        align_selection : str
            Selection for alignment (MDTraj syntax)
        cluster_selection : str
            Selection for clustering (MDTraj syntax)
        start_frame : int
            Starting frame
        end_frame : int, optional
            Ending frame
        step : int
            Step size
        use_pca : bool
            Whether to use PCA for dimensionality reduction
        n_components : int
            Number of PCA components
        cache_name : str, optional
            Cache name

        Returns
        -------
        dict
            Clustering results with labels, centers, and metrics
        """
        try:
            # Create cache key
            if cache_name is None:
                cache_key = f"kmeans_{n_clusters}_{align_selection}_{cluster_selection}_{start_frame}_{end_frame}_{step}_{use_pca}_{n_components}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached clustering results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Load and align trajectory
            traj = self._load_and_align_trajectory(
                universe, trajectory, align_selection,
                start_frame, end_frame, step
            )

            # Extract coordinates for clustering
            coordinates = self._extract_coordinates(traj, cluster_selection)

            # Standardize coordinates
            scaler = StandardScaler()
            coordinates_scaled = scaler.fit_transform(coordinates)

            # Apply PCA if requested
            if use_pca and coordinates_scaled.shape[1] > n_components:
                pca = PCA(n_components=n_components)
                coordinates_reduced = pca.fit_transform(coordinates_scaled)
                logger.info(f"PCA explained variance ratio: {pca.explained_variance_ratio_.sum():.3f}")
            else:
                coordinates_reduced = coordinates_scaled
                pca = None

            # Perform K-means clustering
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            labels = kmeans.fit_predict(coordinates_reduced)

            # Calculate silhouette score
            silhouette = silhouette_score(coordinates_reduced, labels)

            # Prepare results
            results = {
                'labels': labels,
                'cluster_centers': kmeans.cluster_centers_,
                'silhouette_score': silhouette,
                'n_clusters': n_clusters,
                'coordinates_reduced': coordinates_reduced,
                'scaler': scaler,
                'pca': pca,
                'inertia': kmeans.inertia_,
                'frame_indices': list(range(start_frame, end_frame or traj.n_frames, step))
            }

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(results, f)

            logger.info(f"K-means clustering completed. Silhouette score: {silhouette:.3f}")
            return results

        except Exception as e:
            logger.error(f"Error in K-means clustering: {e}")
            raise

    def perform_dbscan_clustering(self,
                                 universe: str,
                                 trajectory: Union[str, List[str]],
                                 eps: float = 0.5,
                                 min_samples: int = 5,
                                 align_selection: str = "protein and name CA",
                                 cluster_selection: str = "protein",
                                 start_frame: int = 0,
                                 end_frame: Optional[int] = None,
                                 step: int = 1,
                                 use_pca: bool = True,
                                 n_components: int = 10,
                                 cache_name: Optional[str] = None) -> Dict:
        """
        Perform DBSCAN clustering on trajectory using pure MDTraj.

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str or list of str
            Path(s) to trajectory file(s)
        eps : float
            DBSCAN eps parameter
        min_samples : int
            DBSCAN min_samples parameter
        align_selection : str
            Selection for alignment (MDTraj syntax)
        cluster_selection : str
            Selection for clustering (MDTraj syntax)
        start_frame : int
            Starting frame
        end_frame : int, optional
            Ending frame
        step : int
            Step size
        use_pca : bool
            Whether to use PCA
        n_components : int
            Number of PCA components
        cache_name : str, optional
            Cache name

        Returns
        -------
        dict
            Clustering results
        """
        try:
            # Load and align trajectory
            traj = self._load_and_align_trajectory(
                universe, trajectory, align_selection,
                start_frame, end_frame, step
            )

            # Extract coordinates for clustering
            coordinates = self._extract_coordinates(traj, cluster_selection)

            # Standardize coordinates
            scaler = StandardScaler()
            coordinates_scaled = scaler.fit_transform(coordinates)

            # Apply PCA if requested
            if use_pca and coordinates_scaled.shape[1] > n_components:
                pca = PCA(n_components=n_components)
                coordinates_reduced = pca.fit_transform(coordinates_scaled)
                logger.info(f"PCA explained variance ratio: {pca.explained_variance_ratio_.sum():.3f}")
            else:
                coordinates_reduced = coordinates_scaled
                pca = None

            # Perform DBSCAN clustering
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            labels = dbscan.fit_predict(coordinates_reduced)

            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            n_noise = list(labels).count(-1)

            # Calculate silhouette score (if there are clusters)
            if n_clusters > 1:
                # Exclude noise points for silhouette score
                mask = labels != -1
                if mask.sum() > 1:
                    silhouette = silhouette_score(coordinates_reduced[mask], labels[mask])
                else:
                    silhouette = -1
            else:
                silhouette = -1

            results = {
                'labels': labels,
                'n_clusters': n_clusters,
                'n_noise': n_noise,
                'silhouette_score': silhouette,
                'eps': eps,
                'min_samples': min_samples,
                'coordinates_reduced': coordinates_reduced,
                'scaler': scaler,
                'pca': pca,
                'frame_indices': list(range(start_frame, end_frame or traj.n_frames, step))
            }

            logger.info(f"DBSCAN clustering completed. Clusters: {n_clusters}, Noise: {n_noise}")
            return results

        except Exception as e:
            logger.error(f"Error in DBSCAN clustering: {e}")
            raise

    def find_optimal_clusters(self,
                             universe: str,
                             trajectory: Union[str, List[str]],
                             max_clusters: int = 10,
                             method: str = "kmeans",
                             align_selection: str = "protein and name CA",
                             cluster_selection: str = "protein",
                             start_frame: int = 0,
                             end_frame: Optional[int] = None,
                             step: int = 1,
                             use_pca: bool = True,
                             n_components: int = 10) -> Dict:
        """
        Find optimal number of clusters using various metrics with pure MDTraj.

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str or list of str
            Path(s) to trajectory file(s)
        max_clusters : int
            Maximum number of clusters to test
        method : str
            Clustering method ('kmeans' or 'agglomerative')
        align_selection : str
            Selection for alignment (MDTraj syntax)
        cluster_selection : str
            Selection for clustering (MDTraj syntax)
        start_frame : int
            Starting frame
        end_frame : int, optional
            Ending frame
        step : int
            Step size
        use_pca : bool
            Whether to use PCA
        n_components : int
            Number of PCA components

        Returns
        -------
        dict
            Optimization results with metrics for different cluster numbers
        """
        try:
            # Load and align trajectory
            traj = self._load_and_align_trajectory(
                universe, trajectory, align_selection,
                start_frame, end_frame, step
            )

            # Extract coordinates for clustering
            coordinates = self._extract_coordinates(traj, cluster_selection)

            # Standardize coordinates
            scaler = StandardScaler()
            coordinates_scaled = scaler.fit_transform(coordinates)

            # Apply PCA if requested
            if use_pca and coordinates_scaled.shape[1] > n_components:
                pca = PCA(n_components=n_components)
                coordinates_reduced = pca.fit_transform(coordinates_scaled)
            else:
                coordinates_reduced = coordinates_scaled

            cluster_range = range(2, max_clusters + 1)
            inertias = []
            silhouette_scores = []

            for n_clusters in cluster_range:
                if method == "kmeans":
                    clusterer = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
                elif method == "agglomerative":
                    clusterer = AgglomerativeClustering(n_clusters=n_clusters)
                else:
                    raise ValueError(f"Unknown clustering method: {method}")

                labels = clusterer.fit_predict(coordinates_reduced)
                silhouette = silhouette_score(coordinates_reduced, labels)
                silhouette_scores.append(silhouette)

                if hasattr(clusterer, 'inertia_'):
                    inertias.append(clusterer.inertia_)

                logger.info(f"Clusters: {n_clusters}, Silhouette: {silhouette:.3f}")

            # Find optimal number of clusters (highest silhouette score)
            optimal_idx = np.argmax(silhouette_scores)
            optimal_clusters = cluster_range[optimal_idx]

            results = {
                'cluster_range': list(cluster_range),
                'silhouette_scores': silhouette_scores,
                'inertias': inertias if inertias else None,
                'optimal_clusters': optimal_clusters,
                'optimal_silhouette': silhouette_scores[optimal_idx],
                'method': method
            }

            logger.info(f"Optimal clusters: {optimal_clusters} (silhouette: {silhouette_scores[optimal_idx]:.3f})")
            return results

        except Exception as e:
            logger.error(f"Error in cluster optimization: {e}")
            raise

    def calculate_rmsd_matrix(self,
                            universe: str,
                            trajectory: Union[str, List[str]],
                            align_selection: str = "protein and name CA",
                            cluster_selection: str = "protein and name CA",
                            start_frame: int = 0,
                            end_frame: Optional[int] = None,
                            step: int = 1) -> np.ndarray:
        """
        Calculate pairwise RMSD matrix using pure MDTraj (GROMOS-style).

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str or list of str
            Path(s) to trajectory file(s)
        align_selection : str
            Selection for alignment
        cluster_selection : str
            Selection for RMSD calculation
        start_frame : int
            Starting frame
        end_frame : int, optional
            Ending frame
        step : int
            Step size

        Returns
        -------
        np.ndarray
            Pairwise RMSD matrix
        """
        # Load and align trajectory
        traj = self._load_and_align_trajectory(
            universe, trajectory, align_selection,
            start_frame, end_frame, step
        )

        # Parse cluster selection for RMSD atoms
        if "segid" in cluster_selection:
            parts = cluster_selection.split("and")
            # Assume CA atoms from main protein chain
            rmsd_indices = []
            for chain in traj.topology.chains:
                if chain.index == 0:  # Main protein chain
                    ca_atoms = [atom.index for atom in chain.atoms if atom.name == 'CA']
                    rmsd_indices.extend(ca_atoms)
                    break
            rmsd_indices = np.array(rmsd_indices)
        else:
            rmsd_indices = traj.topology.select(cluster_selection)

        # Calculate pairwise RMSD matrix
        rmsd_matrix = np.zeros((traj.n_frames, traj.n_frames))

        for i in range(traj.n_frames):
            for j in range(i, traj.n_frames):
                # Calculate RMSD between frames i and j
                rmsd_val = md.rmsd(traj[i:i+1], traj[j:j+1], atom_indices=rmsd_indices)[0]
                rmsd_matrix[i, j] = rmsd_val
                rmsd_matrix[j, i] = rmsd_val  # Symmetric matrix

        logger.info(f"Calculated RMSD matrix: {rmsd_matrix.shape}")
        return rmsd_matrix