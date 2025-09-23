#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clustering analysis for MD trajectories using various algorithms.

Based on SciDraft-Studio implementation with PRISM integration.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
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
    """Clustering analysis for protein-ligand trajectories"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path("./analysis_results/clustering")
        self._cache_dir.mkdir(parents=True, exist_ok=True)

    def _prepare_alignment(self,
                          universe: mda.Universe,
                          align_selection: str = "protein and name CA",
                          cluster_selection: str = "protein",
                          start_frame: int = 0,
                          end_frame: Optional[int] = None,
                          step: int = 1) -> np.ndarray:
        """
        Prepare trajectory for clustering by aligning and extracting coordinates.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe object
        align_selection : str
            Selection for alignment
        cluster_selection : str
            Selection for clustering
        start_frame : int
            Starting frame
        end_frame : int, optional
            Ending frame
        step : int
            Step size

        Returns
        -------
        np.ndarray
            Coordinate matrix (n_frames, n_atoms * 3)
        """
        if end_frame is None:
            end_frame = len(universe.trajectory)

        # Align trajectory
        align.AlignTraj(universe, universe, select=align_selection).run(
            start=start_frame, stop=end_frame, step=step
        )

        # Extract coordinates
        cluster_atoms = universe.select_atoms(cluster_selection)
        n_frames = len(range(start_frame, end_frame, step))
        coordinates = np.zeros((n_frames, len(cluster_atoms) * 3))

        for i, ts in enumerate(universe.trajectory[start_frame:end_frame:step]):
            coordinates[i] = cluster_atoms.positions.flatten()

        logger.info(f"Prepared {n_frames} frames with {len(cluster_atoms)} atoms")
        return coordinates

    def perform_kmeans_clustering(self,
                                 universe: Union[mda.Universe, str],
                                 trajectory: Optional[Union[str, List[str]]] = None,
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
        Perform K-means clustering on trajectory.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or topology path
        trajectory : str, list of str, optional
            Trajectory file(s)
        n_clusters : int
            Number of clusters
        align_selection : str
            Selection for alignment
        cluster_selection : str
            Selection for clustering
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
            # Handle universe input
            if isinstance(universe, str):
                if trajectory is None:
                    raise ValueError("Trajectory path required when universe is a string")
                universe = mda.Universe(universe, trajectory)

            # Create cache key
            if cache_name is None:
                cache_key = f"kmeans_{n_clusters}_{align_selection}_{cluster_selection}_{start_frame}_{end_frame}_{step}_{use_pca}_{n_components}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached clustering results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Prepare coordinates
            coordinates = self._prepare_alignment(
                universe, align_selection, cluster_selection,
                start_frame, end_frame, step
            )

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
                'frame_indices': list(range(start_frame, end_frame, step))
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
                                 universe: Union[mda.Universe, str],
                                 trajectory: Optional[Union[str, List[str]]] = None,
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
        Perform DBSCAN clustering on trajectory.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or topology path
        trajectory : str, list of str, optional
            Trajectory file(s)
        eps : float
            DBSCAN eps parameter
        min_samples : int
            DBSCAN min_samples parameter
        align_selection : str
            Selection for alignment
        cluster_selection : str
            Selection for clustering
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
            # Handle universe input
            if isinstance(universe, str):
                if trajectory is None:
                    raise ValueError("Trajectory path required when universe is a string")
                universe = mda.Universe(universe, trajectory)

            # Prepare coordinates
            coordinates = self._prepare_alignment(
                universe, align_selection, cluster_selection,
                start_frame, end_frame, step
            )

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
                'frame_indices': list(range(start_frame, end_frame, step))
            }

            logger.info(f"DBSCAN clustering completed. Clusters: {n_clusters}, Noise: {n_noise}")
            return results

        except Exception as e:
            logger.error(f"Error in DBSCAN clustering: {e}")
            raise

    def find_optimal_clusters(self,
                             universe: Union[mda.Universe, str],
                             trajectory: Optional[Union[str, List[str]]] = None,
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
        Find optimal number of clusters using various metrics.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or topology path
        trajectory : str, list of str, optional
            Trajectory file(s)
        max_clusters : int
            Maximum number of clusters to test
        method : str
            Clustering method ('kmeans' or 'agglomerative')
        align_selection : str
            Selection for alignment
        cluster_selection : str
            Selection for clustering
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
            # Handle universe input
            if isinstance(universe, str):
                if trajectory is None:
                    raise ValueError("Trajectory path required when universe is a string")
                universe = mda.Universe(universe, trajectory)

            # Prepare coordinates
            coordinates = self._prepare_alignment(
                universe, align_selection, cluster_selection,
                start_frame, end_frame, step
            )

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