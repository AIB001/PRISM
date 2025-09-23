#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD calculation for MD trajectories.

Based on SciDraft-Studio implementation with PRISM integration.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from typing import List, Union, Optional, Tuple, Dict
import os
import logging
from pathlib import Path
import pickle
from scipy.signal import savgol_filter

from ..config import AnalysisConfig

logger = logging.getLogger(__name__)


class RMSDAnalyzer:
    """RMSD analysis for protein-ligand systems"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path("./analysis_results/rmsd")
        self._cache_dir.mkdir(parents=True, exist_ok=True)

    def calculate_rmsd(self,
                      universe: Union[mda.Universe, str],
                      trajectory: Optional[Union[str, List[str]]] = None,
                      selection: str = "protein and name CA",
                      reference: Optional[mda.Universe] = None,
                      ref_frame: int = 0,
                      start_frame: int = 0,
                      end_frame: Optional[int] = None,
                      step: int = 1,
                      cache_name: Optional[str] = None) -> np.ndarray:
        """
        Calculate RMSD for a trajectory.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for RMSD calculation. Default is "protein and name CA".
        reference : MDAnalysis.Universe, optional
            Reference structure. If None, use first frame.
        ref_frame : int, optional
            Reference frame index. Default is 0.
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
            RMSD values for each frame in Angstroms.
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
                cache_key = f"rmsd_{selection}_{ref_frame}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached RMSD results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Set reference
            if reference is None:
                reference = universe
                universe.trajectory[ref_frame]
                ref_coords = universe.select_atoms(selection).positions.copy()
            else:
                reference.trajectory[ref_frame]
                ref_coords = reference.select_atoms(selection).positions.copy()

            # Align trajectory if needed
            align.AlignTraj(universe, reference,
                          select=selection,
                          filename=None).run(start=start_frame, stop=end_frame, step=step)

            # Calculate RMSD
            rmsd_analysis = rms.RMSD(universe.select_atoms(selection),
                                   reference.select_atoms(selection),
                                   ref_frame=ref_frame)
            rmsd_analysis.run(start=start_frame, stop=end_frame, step=step)

            rmsd_values = rmsd_analysis.results.rmsd[:, 2]  # Extract RMSD column

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(rmsd_values, f)

            logger.info(f"RMSD calculation completed. Mean: {np.mean(rmsd_values):.2f} Å")
            return rmsd_values

        except Exception as e:
            logger.error(f"Error calculating RMSD: {e}")
            raise

    def calculate_rmsf(self,
                      universe: Union[mda.Universe, str],
                      trajectory: Optional[Union[str, List[str]]] = None,
                      selection: str = "protein and name CA",
                      start_frame: int = 0,
                      end_frame: Optional[int] = None,
                      step: int = 1,
                      cache_name: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate Root Mean Square Fluctuation (RMSF).

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for RMSF calculation
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
        tuple of np.ndarray
            (rmsf_values, residue_ids) - RMSF per residue and corresponding residue IDs
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
                cache_key = f"rmsf_{selection}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached RMSF results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Align trajectory first
            average = align.AverageStructure(universe, select=selection).run(
                start=start_frame, stop=end_frame, step=step)
            ref = average.results.universe

            aligner = align.AlignTraj(universe, ref, select=selection, in_memory=True)
            aligner.run(start=start_frame, stop=end_frame, step=step)

            # Calculate RMSF
            atoms = universe.select_atoms(selection)
            rmsf_analysis = rms.RMSF(atoms)
            rmsf_analysis.run(start=start_frame, stop=end_frame, step=step)

            rmsf_values = rmsf_analysis.results.rmsf
            residue_ids = np.array([atom.resid for atom in atoms])

            results = (rmsf_values, residue_ids)

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(results, f)

            logger.info(f"RMSF calculation completed. Mean: {np.mean(rmsf_values):.2f} Å")
            return results

        except Exception as e:
            logger.error(f"Error calculating RMSF: {e}")
            raise

    def smooth_rmsd(self, rmsd_values: np.ndarray,
                   window_length: Optional[int] = None,
                   polyorder: int = 3) -> np.ndarray:
        """
        Apply Savitzky-Golay smoothing to RMSD data.

        Parameters
        ----------
        rmsd_values : np.ndarray
            Raw RMSD values
        window_length : int, optional
            Window length for smoothing. If None, use config value.
        polyorder : int, optional
            Polynomial order for smoothing

        Returns
        -------
        np.ndarray
            Smoothed RMSD values
        """
        if window_length is None:
            window_length = self.config.smooth_window

        if len(rmsd_values) < window_length:
            logger.warning("Data too short for smoothing, returning original")
            return rmsd_values

        # Ensure window_length is odd
        if window_length % 2 == 0:
            window_length += 1

        return savgol_filter(rmsd_values, window_length, polyorder)