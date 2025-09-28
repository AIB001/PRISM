#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structural analysis for MD trajectories.

Based on SciDraft-Studio implementation with PRISM integration.
"""

import numpy as np
import MDAnalysis as mda
from typing import List, Union, Optional, Tuple, Dict, Callable
import logging
from pathlib import Path
import pickle

from ..core.config import AnalysisConfig

logger = logging.getLogger(__name__)


class StructuralAnalyzer:
    """Structural properties analysis for protein-ligand systems"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path("./cache")
        self._cache_dir.mkdir(parents=True, exist_ok=True)

    def calculate_radius_of_gyration(self,
                                   universe: Union[mda.Universe, str],
                                   trajectory: Optional[Union[str, List[str]]] = None,
                                   selection: str = "protein",
                                   start_frame: int = 0,
                                   end_frame: Optional[int] = None,
                                   step: int = 1,
                                   cache_name: Optional[str] = None) -> np.ndarray:
        """
        Calculate radius of gyration for a selection over trajectory.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for analysis. Default is "protein".
        start_frame : int, optional
            Starting frame for analysis. Default is 0.
        end_frame : int, optional
            Ending frame for analysis. If None, use all frames.
        step : int, optional
            Step size for frame analysis. Default is 1.
        cache_name : str, optional
            Custom cache name for storing results.

        Returns
        -------
        np.ndarray
            Radius of gyration values for each frame in Angstroms.
        """
        try:
            # Handle universe input
            if isinstance(universe, str):
                if trajectory is None:
                    raise ValueError("Trajectory path required when universe is a string")
                universe = mda.Universe(universe, trajectory)

            # Set up frame range
            if end_frame is None:
                end_frame = len(universe.trajectory)

            # Create cache key
            if cache_name is None:
                cache_key = f"rgyr_{selection}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached radius of gyration results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Get atoms selection
            atoms = universe.select_atoms(selection)
            if len(atoms) == 0:
                raise ValueError(f"No atoms found for selection: {selection}")

            # Calculate radius of gyration for each frame
            rgyr_values = []
            for ts in universe.trajectory[start_frame:end_frame:step]:
                rgyr = atoms.radius_of_gyration()
                rgyr_values.append(rgyr)

            rgyr_values = np.array(rgyr_values)

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(rgyr_values, f)

            logger.info(f"Radius of gyration calculation completed. Mean: {np.mean(rgyr_values):.2f} Å")
            return rgyr_values

        except Exception as e:
            logger.error(f"Error calculating radius of gyration: {e}")
            raise

    def calculate_end_to_end_distance(self,
                                    universe: Union[mda.Universe, str],
                                    trajectory: Optional[Union[str, List[str]]] = None,
                                    selection: str = "protein",
                                    start_frame: int = 0,
                                    end_frame: Optional[int] = None,
                                    step: int = 1,
                                    cache_name: Optional[str] = None) -> np.ndarray:
        """
        Calculate end-to-end distance (N-terminus to C-terminus).

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for protein. Default is "protein".
        start_frame : int, optional
            Starting frame for analysis
        end_frame : int, optional
            Ending frame for analysis
        step : int, optional
            Step size for frame analysis
        cache_name : str, optional
            Custom cache name

        Returns
        -------
        np.ndarray
            End-to-end distances for each frame in Angstroms.
        """
        try:
            # Handle universe input
            if isinstance(universe, str):
                if trajectory is None:
                    raise ValueError("Trajectory path required when universe is a string")
                universe = mda.Universe(universe, trajectory)

            # Set up frame range
            if end_frame is None:
                end_frame = len(universe.trajectory)

            # Create cache key
            if cache_name is None:
                cache_key = f"end2end_{selection}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached end-to-end distance results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Get first and last CA atoms
            try:
                first_ca = universe.select_atoms(f"({selection}) and name CA")[0]
                last_ca = universe.select_atoms(f"({selection}) and name CA")[-1]
            except IndexError:
                raise ValueError("Could not find CA atoms for end-to-end distance calculation")

            # Calculate end-to-end distance for each frame
            distances = []
            for ts in universe.trajectory[start_frame:end_frame:step]:
                distance = np.linalg.norm(first_ca.position - last_ca.position)
                distances.append(distance)

            distances = np.array(distances)

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(distances, f)

            logger.info(f"End-to-end distance calculation completed. Mean: {np.mean(distances):.2f} Å")
            return distances

        except Exception as e:
            logger.error(f"Error calculating end-to-end distance: {e}")
            raise

    def calculate_custom_property(self,
                                universe: Union[mda.Universe, str],
                                property_func: Callable,
                                trajectory: Optional[Union[str, List[str]]] = None,
                                selection: str = "protein",
                                start_frame: int = 0,
                                end_frame: Optional[int] = None,
                                step: int = 1,
                                cache_name: Optional[str] = None) -> np.ndarray:
        """
        Calculate a custom property over trajectory frames.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        property_func : callable
            Function that takes atoms selection and returns a scalar value
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for analysis
        start_frame : int, optional
            Starting frame for analysis
        end_frame : int, optional
            Ending frame for analysis
        step : int, optional
            Step size for frame analysis
        cache_name : str, optional
            Custom cache name

        Returns
        -------
        np.ndarray
            Property values for each frame

        Examples
        --------
        >>> def custom_distance(atoms):
        ...     return np.linalg.norm(atoms[0].position - atoms[-1].position)
        >>> analyzer = StructuralAnalyzer(config)
        >>> distances = analyzer.calculate_custom_property(
        ...     universe, custom_distance, selection="protein and name CA"
        ... )
        """
        try:
            # Handle universe input
            if isinstance(universe, str):
                if trajectory is None:
                    raise ValueError("Trajectory path required when universe is a string")
                universe = mda.Universe(universe, trajectory)

            # Set up frame range
            if end_frame is None:
                end_frame = len(universe.trajectory)

            # Create cache key
            if cache_name is None:
                func_name = getattr(property_func, '__name__', 'custom')
                cache_key = f"custom_{func_name}_{selection}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached custom property results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Get atoms selection
            atoms = universe.select_atoms(selection)
            if len(atoms) == 0:
                raise ValueError(f"No atoms found for selection: {selection}")

            # Calculate property for each frame
            values = []
            for ts in universe.trajectory[start_frame:end_frame:step]:
                value = property_func(atoms)
                values.append(value)

            values = np.array(values)

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(values, f)

            logger.info(f"Custom property calculation completed. Mean: {np.mean(values):.2f}")
            return values

        except Exception as e:
            logger.error(f"Error calculating custom property: {e}")
            raise

    def compare_property_distributions(self,
                                     universes: List[Union[mda.Universe, str]],
                                     property_func: Optional[Callable] = None,
                                     trajectories: Optional[List[str]] = None,
                                     labels: Optional[List[str]] = None,
                                     selection: str = "protein",
                                     start_frame: int = 0,
                                     end_frame: Optional[int] = None,
                                     step: int = 1) -> Dict[str, List[float]]:
        """
        Compare property distributions across multiple systems.

        Parameters
        ----------
        universes : list of MDAnalysis.Universe or str
            List of universe objects or topology paths
        property_func : callable, optional
            Function to calculate property. If None, uses radius of gyration.
        trajectories : list of str, optional
            List of trajectory paths
        labels : list of str, optional
            Labels for each system
        selection : str, optional
            Selection string for analysis
        start_frame : int, optional
            Starting frame for analysis
        end_frame : int, optional
            Ending frame for analysis
        step : int, optional
            Step size for frame analysis

        Returns
        -------
        dict
            Dictionary mapping labels to property value lists

        Examples
        --------
        >>> def gyration_radius(atoms):
        ...     return atoms.radius_of_gyration()
        >>> results = analyzer.compare_property_distributions(
        ...     universes=[u1, u2, u3],
        ...     property_func=gyration_radius,
        ...     labels=['WT', 'Mutant1', 'Mutant2']
        ... )
        """
        if labels is None:
            labels = [f"System {i + 1}" for i in range(len(universes))]

        if len(labels) != len(universes):
            raise ValueError("Number of labels must match number of universes")

        results = {}

        for i, (univ, label) in enumerate(zip(universes, labels)):
            # Create universe if needed
            if isinstance(univ, str):
                if trajectories and i < len(trajectories) and trajectories[i]:
                    univ = mda.Universe(univ, trajectories[i])
                else:
                    raise ValueError(f"Trajectory needed for universe {i}")

            # Get selection
            atoms = univ.select_atoms(selection)

            # Calculate property for each frame
            values = []
            frame_range = range(start_frame,
                              end_frame if end_frame else len(univ.trajectory),
                              step)

            for frame_idx in frame_range:
                univ.trajectory[frame_idx]
                if property_func:
                    value = property_func(atoms)
                else:
                    # Default to radius of gyration
                    value = atoms.radius_of_gyration()
                values.append(value)

            results[label] = values

        logger.info(f"Property distribution comparison completed for {len(universes)} systems")
        return results