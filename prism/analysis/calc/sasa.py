#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solvent Accessible Surface Area (SASA) analysis for MD trajectories.

Uses MDAnalysis capabilities for SASA calculation.
"""

import numpy as np
import MDAnalysis as mda
from typing import List, Union, Optional, Tuple, Dict
import logging
from pathlib import Path
import pickle

try:
    from MDAnalysis.analysis.sasa import SASAAnalysis

    SASA_AVAILABLE = True
except ImportError:
    SASA_AVAILABLE = False

from ..core.config import AnalysisConfig

logger = logging.getLogger(__name__)


class SASAAnalyzer:
    """Solvent Accessible Surface Area analysis for protein-ligand systems"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path("./cache")
        self._cache_dir.mkdir(parents=True, exist_ok=True)

        if not SASA_AVAILABLE:
            logger.warning("MDAnalysis SASA module not available. " "Some functionality may be limited.")

    def calculate_sasa(
        self,
        universe: Union[mda.Universe, str],
        trajectory: Optional[Union[str, List[str]]] = None,
        selection: str = "protein",
        probe_radius: float = 1.4,
        start_frame: int = 0,
        end_frame: Optional[int] = None,
        step: int = 1,
        cache_name: Optional[str] = None,
    ) -> np.ndarray:
        """
        Calculate total SASA for a selection over trajectory.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for analysis. Default is "protein".
        probe_radius : float, optional
            Probe radius for SASA calculation in Angstroms. Default is 1.4.
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
            Total SASA values for each frame in Angstroms^2.
        """
        if not SASA_AVAILABLE:
            raise ImportError(
                "MDAnalysis SASA module not available. " "Please update MDAnalysis or use alternative method."
            )

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
                cache_key = f"sasa_{selection}_{probe_radius}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached SASA results from {cache_file}")
                with open(cache_file, "rb") as f:
                    return pickle.load(f)

            # Get atoms selection
            atoms = universe.select_atoms(selection)
            if len(atoms) == 0:
                raise ValueError(f"No atoms found for selection: {selection}")

            # Run SASA analysis
            sasa_analysis = SASAAnalysis(atoms, probe_radius=probe_radius)
            sasa_analysis.run(start=start_frame, stop=end_frame, step=step)

            # Extract total SASA values
            sasa_values = sasa_analysis.results.sasa

            # Cache results
            with open(cache_file, "wb") as f:
                pickle.dump(sasa_values, f)

            logger.info(f"SASA calculation completed. Mean: {np.mean(sasa_values):.2f} Å²")
            return sasa_values

        except Exception as e:
            logger.error(f"Error calculating SASA: {e}")
            raise

    def calculate_sasa_per_residue(
        self,
        universe: Union[mda.Universe, str],
        trajectory: Optional[Union[str, List[str]]] = None,
        selection: str = "protein",
        probe_radius: float = 1.4,
        start_frame: int = 0,
        end_frame: Optional[int] = None,
        step: int = 1,
        cache_name: Optional[str] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate SASA per residue over trajectory.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for analysis
        probe_radius : float, optional
            Probe radius for SASA calculation
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
            (sasa_per_residue, residue_ids) - SASA per residue per frame and residue IDs
        """
        if not SASA_AVAILABLE:
            raise ImportError("MDAnalysis SASA module not available.")

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
                cache_key = f"sasa_residue_{selection}_{probe_radius}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached per-residue SASA results from {cache_file}")
                with open(cache_file, "rb") as f:
                    return pickle.load(f)

            # Get atoms selection
            atoms = universe.select_atoms(selection)
            if len(atoms) == 0:
                raise ValueError(f"No atoms found for selection: {selection}")

            # Group by residue
            residues = atoms.residues
            residue_ids = [res.resid for res in residues]

            # Calculate SASA for each residue
            n_frames = len(range(start_frame, end_frame, step))
            sasa_per_residue = np.zeros((n_frames, len(residues)))

            for i, residue in enumerate(residues):
                res_atoms = residue.atoms
                sasa_analysis = SASAAnalysis(res_atoms, probe_radius=probe_radius)
                sasa_analysis.run(start=start_frame, stop=end_frame, step=step)
                sasa_per_residue[:, i] = sasa_analysis.results.sasa

            results = (sasa_per_residue, np.array(residue_ids))

            # Cache results
            with open(cache_file, "wb") as f:
                pickle.dump(results, f)

            logger.info(f"Per-residue SASA calculation completed for {len(residues)} residues")
            return results

        except Exception as e:
            logger.error(f"Error calculating per-residue SASA: {e}")
            raise

    def calculate_sasa_difference(
        self,
        universe1: Union[mda.Universe, str],
        universe2: Union[mda.Universe, str],
        trajectory1: Optional[Union[str, List[str]]] = None,
        trajectory2: Optional[Union[str, List[str]]] = None,
        selection: str = "protein",
        probe_radius: float = 1.4,
        start_frame: int = 0,
        end_frame: Optional[int] = None,
        step: int = 1,
    ) -> Dict[str, np.ndarray]:
        """
        Calculate SASA difference between two systems (e.g., apo vs holo).

        Parameters
        ----------
        universe1 : MDAnalysis.Universe or str
            First universe (e.g., apo state)
        universe2 : MDAnalysis.Universe or str
            Second universe (e.g., holo state)
        trajectory1 : str, list of str, optional
            Trajectory for first system
        trajectory2 : str, list of str, optional
            Trajectory for second system
        selection : str, optional
            Selection string for analysis
        probe_radius : float, optional
            Probe radius for SASA calculation
        start_frame : int, optional
            Starting frame for analysis
        end_frame : int, optional
            Ending frame for analysis
        step : int, optional
            Step size for frame analysis

        Returns
        -------
        dict
            Dictionary with 'system1', 'system2', and 'difference' SASA values
        """
        try:
            # Calculate SASA for both systems
            sasa1 = self.calculate_sasa(
                universe1,
                trajectory1,
                selection,
                probe_radius,
                start_frame,
                end_frame,
                step,
                cache_name=f"sasa_sys1_{hash(str(universe1))}",
            )

            sasa2 = self.calculate_sasa(
                universe2,
                trajectory2,
                selection,
                probe_radius,
                start_frame,
                end_frame,
                step,
                cache_name=f"sasa_sys2_{hash(str(universe2))}",
            )

            # Calculate difference
            min_length = min(len(sasa1), len(sasa2))
            sasa1_trimmed = sasa1[:min_length]
            sasa2_trimmed = sasa2[:min_length]
            difference = sasa2_trimmed - sasa1_trimmed

            results = {
                "system1": sasa1_trimmed,
                "system2": sasa2_trimmed,
                "difference": difference,
                "mean_difference": np.mean(difference),
                "std_difference": np.std(difference),
            }

            logger.info(
                f"SASA difference calculation completed. "
                f"Mean difference: {results['mean_difference']:.2f} ± "
                f"{results['std_difference']:.2f} Å²"
            )

            return results

        except Exception as e:
            logger.error(f"Error calculating SASA difference: {e}")
            raise

    def calculate_buried_surface_area(
        self,
        universe: Union[mda.Universe, str],
        trajectory: Optional[Union[str, List[str]]] = None,
        selection1: str = "protein",
        selection2: str = "resname LIG",
        probe_radius: float = 1.4,
        start_frame: int = 0,
        end_frame: Optional[int] = None,
        step: int = 1,
        cache_name: Optional[str] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Calculate buried surface area between two selections (e.g., protein-ligand).

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection1 : str, optional
            First selection (e.g., protein)
        selection2 : str, optional
            Second selection (e.g., ligand)
        probe_radius : float, optional
            Probe radius for SASA calculation
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
        dict
            Dictionary with SASA values and buried surface area
        """
        try:
            # Calculate SASA for individual components
            sasa1 = self.calculate_sasa(
                universe,
                trajectory,
                selection1,
                probe_radius,
                start_frame,
                end_frame,
                step,
                cache_name=f"bsa_comp1_{cache_name}" if cache_name else None,
            )

            sasa2 = self.calculate_sasa(
                universe,
                trajectory,
                selection2,
                probe_radius,
                start_frame,
                end_frame,
                step,
                cache_name=f"bsa_comp2_{cache_name}" if cache_name else None,
            )

            # Calculate SASA for the complex
            complex_selection = f"({selection1}) or ({selection2})"
            sasa_complex = self.calculate_sasa(
                universe,
                trajectory,
                complex_selection,
                probe_radius,
                start_frame,
                end_frame,
                step,
                cache_name=f"bsa_complex_{cache_name}" if cache_name else None,
            )

            # Calculate buried surface area
            # BSA = (SASA1 + SASA2 - SASA_complex) / 2
            min_length = min(len(sasa1), len(sasa2), len(sasa_complex))
            bsa = (sasa1[:min_length] + sasa2[:min_length] - sasa_complex[:min_length]) / 2

            results = {
                "sasa_component1": sasa1[:min_length],
                "sasa_component2": sasa2[:min_length],
                "sasa_complex": sasa_complex[:min_length],
                "buried_surface_area": bsa,
                "mean_bsa": np.mean(bsa),
                "std_bsa": np.std(bsa),
            }

            logger.info(
                f"Buried surface area calculation completed. "
                f"Mean BSA: {results['mean_bsa']:.2f} ± "
                f"{results['std_bsa']:.2f} Å²"
            )

            return results

        except Exception as e:
            logger.error(f"Error calculating buried surface area: {e}")
            raise
