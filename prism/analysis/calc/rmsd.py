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
from ...utils.topology import detect_all_chains

logger = logging.getLogger(__name__)


class RMSDAnalyzer:
    """RMSD analysis for protein-ligand systems"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path("./cache")
        self._cache_dir.mkdir(parents=True, exist_ok=True)

    def calculate_rmsd(self,
                      universe: Union[mda.Universe, str],
                      trajectory: Optional[Union[str, List[str]]] = None,
                      selection: str = "protein and name CA",
                      align_selection: Optional[str] = None,
                      calculate_selection: Optional[str] = None,
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
            Deprecated: use align_selection and calculate_selection instead.
        align_selection : str, optional
            Selection string for trajectory alignment. If None, uses selection.
            Default behavior aligns on "protein and name CA".
        calculate_selection : str, optional
            Selection string for RMSD calculation. If None, uses align_selection.
            Allows calculating RMSD on different atoms than alignment (e.g., ligand after protein alignment).
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

            # Handle selection parameters
            if align_selection is None:
                align_selection = selection
            if calculate_selection is None:
                calculate_selection = align_selection

            # Create cache key
            if cache_name is None:
                cache_key = f"rmsd_align_{align_selection}_calc_{calculate_selection}_{ref_frame}_{start_frame}_{end_frame}_{step}"
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
                ref_coords = universe.select_atoms(calculate_selection).positions.copy()
            else:
                reference.trajectory[ref_frame]
                ref_coords = reference.select_atoms(calculate_selection).positions.copy()

            # Align trajectory if needed with parallel processing
            logger.info(f"Aligning trajectory on '{align_selection}' using {os.environ.get('OMP_NUM_THREADS', 'auto')} CPU threads")
            align.AlignTraj(universe, reference,
                          select=align_selection,
                          filename=None).run(start=start_frame, stop=end_frame, step=step)

            # Calculate RMSD on specified selection
            logger.info(f"Calculating RMSD for '{calculate_selection}'")
            rmsd_analysis = rms.RMSD(universe.select_atoms(calculate_selection),
                                   reference.select_atoms(calculate_selection),
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
                      align_selection: Optional[str] = None,
                      calculate_selection: Optional[str] = None,
                      start_frame: int = 0,
                      end_frame: Optional[int] = None,
                      step: int = 1,
                      cache_name: Optional[str] = None) -> Union[Tuple[np.ndarray, np.ndarray], Dict[str, Tuple[np.ndarray, np.ndarray]]]:
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
            Deprecated: use align_selection and calculate_selection instead.
        align_selection : str, optional
            Selection string for trajectory alignment. If None, uses selection.
            Default behavior aligns on "protein and name CA".
        calculate_selection : str, optional
            Selection string for RMSF calculation. If None, auto-detects protein chains.
            When None, returns dict with chain-separated RMSF results.
            When specified, returns single tuple for that selection.
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
        Union[Tuple[np.ndarray, np.ndarray], Dict[str, Tuple[np.ndarray, np.ndarray]]]
            If calculate_selection specified: (rmsf_values, residue_ids)
            If calculate_selection is None: {"Chain X": (rmsf_values, residue_ids), ...}
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

            # Handle selection parameters
            if align_selection is None:
                align_selection = selection

            # Auto-detect chains if calculate_selection is None
            if calculate_selection is None:
                all_chains = detect_all_chains(universe)
                # Combine protein and nucleic chains
                combined_chains = {}
                combined_chains.update(all_chains["protein"])
                combined_chains.update(all_chains["nucleic"])

                if len(combined_chains) == 0:
                    # Fallback to original selection
                    calculate_selection = align_selection
                    single_selection = True
                else:
                    single_selection = False
                    protein_chains = combined_chains
            else:
                single_selection = True

            # Create cache key
            if cache_name is None:
                if single_selection:
                    calc_key = calculate_selection if calculate_selection else align_selection
                    cache_key = f"rmsf_align_{align_selection}_calc_{calc_key}_{start_frame}_{end_frame}_{step}"
                else:
                    cache_key = f"rmsf_align_{align_selection}_calc_auto_chains_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached RMSF results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Align trajectory first with parallel processing
            logger.info(f"Computing average structure using {os.environ.get('OMP_NUM_THREADS', 'auto')} CPU threads")
            average = align.AverageStructure(universe, select=align_selection).run(
                start=start_frame, stop=end_frame, step=step)
            ref = average.results.universe

            logger.info(f"Aligning trajectory on '{align_selection}' to average structure")
            aligner = align.AlignTraj(universe, ref, select=align_selection, in_memory=True)
            aligner.run(start=start_frame, stop=end_frame, step=step)

            # Calculate RMSF
            if single_selection:
                # Single selection RMSF
                calc_sel = calculate_selection if calculate_selection else align_selection
                logger.info(f"Calculating RMSF for '{calc_sel}' with parallel processing")
                atoms = universe.select_atoms(calc_sel)
                rmsf_analysis = rms.RMSF(atoms)
                rmsf_analysis.run(start=start_frame, stop=end_frame, step=step)

                rmsf_values = rmsf_analysis.results.rmsf
                residue_ids = np.array([atom.resid for atom in atoms])
                results = (rmsf_values, residue_ids)
            else:
                # Multi-chain RMSF
                logger.info("Calculating RMSF for multiple chains (protein + nucleic) with parallel processing")
                results = {}
                for chain_name, chain_selection in protein_chains.items():
                    logger.info(f"  Processing {chain_name}: {chain_selection}")
                    atoms = universe.select_atoms(chain_selection)
                    if len(atoms) == 0:
                        logger.warning(f"No atoms found for {chain_name} with selection '{chain_selection}'")
                        continue

                    rmsf_analysis = rms.RMSF(atoms)
                    rmsf_analysis.run(start=start_frame, stop=end_frame, step=step)

                    rmsf_values = rmsf_analysis.results.rmsf
                    residue_ids = np.array([atom.resid for atom in atoms])
                    results[chain_name] = (rmsf_values, residue_ids)

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(results, f)

            if single_selection:
                rmsf_mean = np.mean(results[0])
                logger.info(f"RMSF calculation completed. Mean: {rmsf_mean:.2f} Å")
            else:
                logger.info(f"RMSF calculation completed for {len(results)} chains")
                for chain_name, (rmsf_vals, _) in results.items():
                    logger.info(f"  {chain_name}: Mean RMSF = {np.mean(rmsf_vals):.2f} Å")

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