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
import mdtraj as md

from ..core.config import AnalysisConfig
from ...utils.topology import detect_all_chains
from ...utils.mdtraj_topology import (
    detect_all_chains_mdtraj, get_primary_protein_chain_mdtraj,
    create_chain_selections_mdtraj, validate_chain_selections_mdtraj
)

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
                      universe: Union[mda.Universe, str, md.Trajectory],
                      trajectory: Optional[Union[str, List[str]]] = None,
                      selection: str = "protein and name CA",
                      align_selection: Optional[str] = None,
                      calculate_selection: Optional[str] = None,
                      start_frame: int = 0,
                      end_frame: Optional[int] = None,
                      step: int = 1,
                      cache_name: Optional[str] = None,
                      use_mdtraj: bool = True) -> Union[Tuple[np.ndarray, np.ndarray], Dict[str, Tuple[np.ndarray, np.ndarray]]]:
        """
        Calculate Root Mean Square Fluctuation (RMSF) with automatic chain detection.

        This method calculates RMSF for molecular dynamics trajectories with automatic
        detection of protein and nucleic acid chains. Uses proper atom selections:
        CA atoms for protein chains and P atoms for nucleic acid chains.

        Parameters
        ----------
        universe : MDAnalysis.Universe, mdtraj.Trajectory, or str
            Universe/trajectory object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for RMSF calculation. Default is "protein and name CA".
            Deprecated: use align_selection and calculate_selection instead.
        align_selection : str, optional
            Selection string for trajectory alignment. If None, uses selection.
            Default behavior aligns on primary protein chain CA atoms.
        calculate_selection : str, optional
            Selection string for RMSF calculation. If None, auto-detects all chains.
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
        use_mdtraj : bool, optional
            Whether to use MDTraj for chain detection (recommended). Default True.

        Returns
        -------
        Union[Tuple[np.ndarray, np.ndarray], Dict[str, Tuple[np.ndarray, np.ndarray]]]
            If calculate_selection specified: (rmsf_values, residue_ids)
            If calculate_selection is None: {"Chain Name": (rmsf_values, residue_ids), ...}
            where Chain Name includes "Protein Chain X" or "Nucleic Chain Y"
        """
        try:
            # Handle input - support both MDAnalysis and MDTraj
            mdtraj_trajectory = None

            if use_mdtraj:
                # Try to use MDTraj for chain detection and analysis
                if isinstance(universe, md.Trajectory):
                    mdtraj_trajectory = universe
                elif isinstance(universe, str):
                    if trajectory is None:
                        raise ValueError("Trajectory path required when universe is a string")
                    # Load with MDTraj
                    logger.info("Loading trajectory with MDTraj for enhanced chain detection")
                    mdtraj_trajectory = md.load(trajectory, top=universe)
                else:
                    # Convert MDAnalysis to MDTraj if needed
                    logger.info("Converting MDAnalysis universe for MDTraj analysis")
                    # For now, fallback to MDAnalysis workflow
                    use_mdtraj = False

            if not use_mdtraj or mdtraj_trajectory is None:
                # Fallback to original MDAnalysis workflow
                if isinstance(universe, str):
                    if trajectory is None:
                        raise ValueError("Trajectory path required when universe is a string")
                    universe = mda.Universe(universe, trajectory)
                elif isinstance(universe, md.Trajectory):
                    raise ValueError("MDTraj trajectory provided but use_mdtraj=False")

            # Set up frame range
            if use_mdtraj and mdtraj_trajectory is not None:
                if end_frame is None:
                    end_frame = len(mdtraj_trajectory)
                # Apply frame slicing
                if start_frame > 0 or end_frame < len(mdtraj_trajectory) or step > 1:
                    mdtraj_trajectory = mdtraj_trajectory[start_frame:end_frame:step]
                    logger.info(f"Using frames {start_frame}:{end_frame}:{step} ({len(mdtraj_trajectory)} frames)")
            else:
                if end_frame is None:
                    end_frame = len(universe.trajectory)

            # Handle selection parameters and chain detection
            if align_selection is None:
                align_selection = selection

            # Auto-detect chains if calculate_selection is None
            if calculate_selection is None:
                if use_mdtraj and mdtraj_trajectory is not None:
                    # Use MDTraj chain detection
                    logger.info("Auto-detecting chains using MDTraj")
                    all_chains = detect_all_chains_mdtraj(mdtraj_trajectory)

                    # Get primary protein chain for alignment
                    primary_protein = get_primary_protein_chain_mdtraj(mdtraj_trajectory)
                    if primary_protein:
                        align_selection = primary_protein
                        logger.info(f"Using primary protein chain for alignment: {align_selection}")

                    # Combine protein and nucleic chains
                    combined_chains = {}
                    combined_chains.update(all_chains["protein"])
                    combined_chains.update(all_chains["nucleic"])
                else:
                    # Use MDAnalysis chain detection (fallback)
                    all_chains = detect_all_chains(universe)
                    combined_chains = {}
                    combined_chains.update(all_chains["protein"])
                    combined_chains.update(all_chains["nucleic"])

                if len(combined_chains) == 0:
                    # Fallback to original selection
                    calculate_selection = align_selection
                    single_selection = True
                    logger.warning("No chains auto-detected, using fallback selection")
                else:
                    single_selection = False
                    protein_chains = combined_chains
                    logger.info(f"Auto-detected {len(combined_chains)} chains for RMSF analysis")
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

            # Calculate RMSF using appropriate method
            if use_mdtraj and mdtraj_trajectory is not None:
                # Use MDTraj for RMSF calculation with OpenMP parallelization
                logger.info("Using MDTraj for RMSF calculation with automatic alignment")

                # Configure OpenMP threads for optimal performance
                original_omp_threads = os.environ.get('OMP_NUM_THREADS')
                cpu_count = os.cpu_count()
                recommended_threads = max(1, cpu_count - 1)  # Leave one core for system tasks
                os.environ['OMP_NUM_THREADS'] = str(recommended_threads)
                logger.info(f"Configured OpenMP with {recommended_threads} threads (available CPUs: {cpu_count})")

                try:
                    if single_selection:
                        # Single selection RMSF with MDTraj
                        calc_sel = calculate_selection if calculate_selection else align_selection
                        results = self._calculate_rmsf_mdtraj_single(mdtraj_trajectory, calc_sel, align_selection)
                    else:
                        # Multi-chain RMSF with MDTraj
                        results = self._calculate_rmsf_mdtraj_multi(mdtraj_trajectory, protein_chains, align_selection)
                finally:
                    # Restore original OMP_NUM_THREADS setting
                    if original_omp_threads is not None:
                        os.environ['OMP_NUM_THREADS'] = original_omp_threads
                    else:
                        os.environ.pop('OMP_NUM_THREADS', None)
            else:
                # Use MDAnalysis for RMSF calculation (fallback)
                logger.info("Using MDAnalysis for RMSF calculation")

                # Align trajectory first with parallel processing
                logger.info(f"Computing average structure using {os.environ.get('OMP_NUM_THREADS', 'auto')} CPU threads")
                average = align.AverageStructure(universe, select=align_selection).run(
                    start=start_frame, stop=end_frame, step=step)
                ref = average.results.universe

                logger.info(f"Aligning trajectory on '{align_selection}' to average structure")
                aligner = align.AlignTraj(universe, ref, select=align_selection, in_memory=True)
                aligner.run(start=start_frame, stop=end_frame, step=step)

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

    def calculate_rmsd_mdtraj(self,
                              topology_file: str,
                              trajectory_files: Union[str, List[str]],
                              align_selection: Optional[str] = None,
                              calculate_selection: Optional[str] = None,
                              auto_detect_chains: bool = True,
                              start_frame: int = 0,
                              end_frame: Optional[int] = None,
                              step: int = 1,
                              ref_frame: int = 0,
                              cache_name: Optional[str] = None) -> Union[np.ndarray, Dict[str, np.ndarray]]:
        """
        Calculate RMSD using MDTraj with automatic chain detection.

        Parameters
        ----------
        topology_file : str
            Path to topology file
        trajectory_files : str or list of str
            Path(s) to trajectory file(s)
        align_selection : str, optional
            Selection string for trajectory alignment. If None and auto_detect_chains is True,
            uses primary protein chain CA atoms (preferentially chain P).
        calculate_selection : str, optional
            Selection string for RMSD calculation. If None and auto_detect_chains is True,
            calculates both protein and ligand RMSD separately.
        auto_detect_chains : bool
            Whether to automatically detect chains. Default is True.
        start_frame : int, optional
            Starting frame for analysis
        end_frame : int, optional
            Ending frame for analysis
        step : int, optional
            Step size for frame analysis
        ref_frame : int, optional
            Reference frame index. Default is 0.
        cache_name : str, optional
            Custom cache name

        Returns
        -------
        Union[np.ndarray, Dict[str, np.ndarray]]
            If calculate_selection specified: RMSD time series (frames,)
            If auto_detect_chains is True: {"protein": rmsd_array, "ligand": rmsd_array}
            Time units are in nanoseconds, RMSD values in Angstroms.
        """
        try:
            # Load trajectory using MDTraj
            logger.info("Loading trajectory with MDTraj")
            traj = md.load(trajectory_files, top=topology_file)

            # Apply frame slicing
            traj = traj[start_frame:end_frame:step]
            logger.info(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

            # Handle selection parameters with chain P preference
            if align_selection is None and auto_detect_chains:
                # Try to find chain P first, then fallback to primary protein chain
                from ...utils.ligand import identify_ligand_residue

                # Check for chain P specifically
                try:
                    chain_p_atoms = traj.topology.select("protein and chainid P and name CA")
                    if len(chain_p_atoms) > 0:
                        align_selection = "protein and chainid P and name CA"
                        logger.info(f"Using chain P for alignment: {len(chain_p_atoms)} CA atoms")
                    else:
                        # Try chain 0 (often the main protein chain)
                        chain_0_atoms = traj.topology.select("protein and chainid 0 and name CA")
                        if len(chain_0_atoms) > 0:
                            align_selection = "protein and chainid 0 and name CA"
                            logger.info(f"Using chain 0 for alignment: {len(chain_0_atoms)} CA atoms")
                        else:
                            # Fallback to all protein CA atoms
                            align_selection = "protein and name CA"
                            logger.info("Using all protein CA atoms for alignment")
                except:
                    align_selection = "protein and name CA"
                    logger.info("Fallback to all protein CA atoms for alignment")

            # Handle calculate_selection for auto-detection
            single_selection = False
            if calculate_selection is None and auto_detect_chains:
                # Auto-detect protein and ligand for separate RMSD calculations
                selections = {}

                # Protein RMSD (same as alignment selection, but can be different)
                if "chainid P" in align_selection:
                    selections["protein"] = "protein and chainid P and name CA"
                elif "chainid 0" in align_selection:
                    selections["protein"] = "protein and chainid 0 and name CA"
                else:
                    selections["protein"] = "protein and name CA"

                # Ligand RMSD - try to identify ligand
                try:
                    ligand = identify_ligand_residue(traj)
                    if ligand:
                        ligand_selection = f"resname {ligand.name}"
                        ligand_atoms = traj.topology.select(ligand_selection)
                        if len(ligand_atoms) > 0:
                            selections["ligand"] = ligand_selection
                            logger.info(f"Found ligand {ligand.name}: {len(ligand_atoms)} atoms")
                        else:
                            logger.warning(f"Ligand {ligand.name} found but no atoms selected")
                    else:
                        logger.warning("No ligand found for RMSD calculation")
                except Exception as e:
                    logger.warning(f"Ligand detection failed: {e}")

                if len(selections) == 0:
                    # Fallback to single protein selection
                    calculate_selection = align_selection
                    single_selection = True
                    logger.warning("No selections found, falling back to alignment selection")
                else:
                    logger.info(f"Using auto-detected selections: {list(selections.keys())}")
                    protein_ligand_selections = selections
            else:
                single_selection = True
                if calculate_selection is None:
                    calculate_selection = align_selection if align_selection else "protein and name CA"

            # Create cache key
            if cache_name is None:
                if single_selection:
                    calc_key = str(calculate_selection if calculate_selection else "protein")
                    align_key = str(align_selection if align_selection else "protein")
                    cache_key = f"rmsd_mdtraj_align_{align_key}_calc_{calc_key}_{ref_frame}_{start_frame}_{end_frame}_{step}"
                else:
                    sel_names = "_".join(sorted(protein_ligand_selections.keys()))
                    align_key = str(align_selection if align_selection else "auto")
                    cache_key = f"rmsd_mdtraj_align_{align_key}_calc_{sel_names}_{ref_frame}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "").replace("{", "").replace("}", "")
            else:
                cache_key = str(cache_name)

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached RMSD results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Align trajectory using MDTraj
            if align_selection:
                logger.info(f"Aligning trajectory on '{align_selection}' using MDTraj")
                align_indices = traj.topology.select(align_selection)

                if len(align_indices) > 0:
                    # Superpose to reference frame
                    traj.superpose(traj, frame=ref_frame, atom_indices=align_indices)
                    logger.info(f"Alignment completed on {len(align_indices)} atoms to frame {ref_frame}")
                else:
                    logger.warning(f"No atoms found for alignment selection: {align_selection}")

            # Calculate RMSD
            if single_selection:
                # Single selection RMSD
                calc_sel = calculate_selection if calculate_selection else "protein and name CA"
                logger.info(f"Calculating RMSD for '{calc_sel}'")

                calc_indices = traj.topology.select(calc_sel)
                if len(calc_indices) == 0:
                    raise ValueError(f"No atoms found for calculation selection: {calc_sel}")

                # Calculate RMSD using MDTraj relative to reference frame
                rmsd_values = md.rmsd(traj, traj, frame=ref_frame, atom_indices=calc_indices)

                results = rmsd_values * 10.0  # Convert nm to Angstroms
                logger.info(f"Single selection RMSD completed. Mean: {np.mean(results):.2f} Å")
            else:
                # Multi-selection RMSD (protein and ligand)
                logger.info("Calculating RMSD for protein and ligand using MDTraj")
                results = {}

                for selection_name, selection_string in protein_ligand_selections.items():
                    logger.info(f"  Processing {selection_name}: {selection_string}")

                    calc_indices = traj.topology.select(selection_string)
                    if len(calc_indices) == 0:
                        logger.warning(f"No atoms found for {selection_name} with selection '{selection_string}'")
                        continue

                    # Calculate RMSD for this selection
                    rmsd_values = md.rmsd(traj, traj, frame=ref_frame, atom_indices=calc_indices)
                    rmsd_values = rmsd_values * 10.0  # Convert nm to Angstroms

                    results[selection_name] = rmsd_values
                    logger.info(f"  {selection_name}: Mean RMSD = {np.mean(rmsd_values):.2f} Å")

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(results, f)

            logger.info("RMSD calculation completed successfully")
            return results

        except Exception as e:
            logger.error(f"Error calculating RMSD with MDTraj: {e}")
            raise

    def calculate_rmsf_mdtraj(self,
                              topology_file: str,
                              trajectory_files: Union[str, List[str]],
                              align_selection: Optional[str] = None,
                              calculate_selection: Optional[str] = None,
                              auto_detect_chains: bool = True,
                              start_frame: int = 0,
                              end_frame: Optional[int] = None,
                              step: int = 1,
                              cache_name: Optional[str] = None) -> Union[Tuple[np.ndarray, np.ndarray], Dict[str, Tuple[np.ndarray, np.ndarray]]]:
        """
        Calculate RMSF using MDTraj with automatic chain detection.

        Parameters
        ----------
        topology_file : str
            Path to topology file
        trajectory_files : str or list of str
            Path(s) to trajectory file(s)
        align_selection : str, optional
            Selection string for trajectory alignment. If None and auto_detect_chains is True,
            uses primary protein chain CA atoms.
        calculate_selection : str, optional
            Selection string for RMSF calculation. If None and auto_detect_chains is True,
            auto-detects all chains.
        auto_detect_chains : bool
            Whether to automatically detect chains. Default is True.
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
            If auto_detect_chains is True: {"Chain X": (rmsf_values, residue_ids), ...}
        """
        try:
            # Load trajectory using MDTraj
            logger.info("Loading trajectory with MDTraj")
            traj = md.load(trajectory_files, top=topology_file)

            # Apply frame slicing
            traj = traj[start_frame:end_frame:step]

            logger.info(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

            # Handle selection parameters
            if align_selection is None and auto_detect_chains:
                # Auto-detect primary protein chain for alignment
                align_selection = get_primary_protein_chain_mdtraj(traj)
                if align_selection is None:
                    logger.warning("No protein chain found for alignment, using all CA atoms")
                    align_selection = "name CA"

            # Handle calculate_selection for auto-detection
            single_selection = False
            if calculate_selection is None and auto_detect_chains:
                # Auto-detect all chains
                chain_selections = create_chain_selections_mdtraj(traj, include_protein=True, include_nucleic=True)

                if len(chain_selections) == 0:
                    # Fallback to CA atoms
                    logger.warning("No chains detected, falling back to CA atoms")
                    calculate_selection = "name CA"
                    single_selection = True
                else:
                    # Validate chain selections
                    validation_results = validate_chain_selections_mdtraj(traj, chain_selections)
                    valid_chains = {name: sel for name, sel in chain_selections.items() if validation_results[name]}

                    if len(valid_chains) == 0:
                        logger.warning("No valid chains found, falling back to CA atoms")
                        calculate_selection = "name CA"
                        single_selection = True
                    else:
                        logger.info(f"Using {len(valid_chains)} auto-detected chains for RMSF calculation")
                        protein_chains = valid_chains
            else:
                single_selection = True
                if calculate_selection is None:
                    calculate_selection = align_selection if align_selection else "name CA"

            # Create cache key
            if cache_name is None:
                if single_selection:
                    calc_key = calculate_selection if calculate_selection else "name CA"
                    align_key = align_selection if align_selection else "name CA"
                    cache_key = f"rmsf_mdtraj_align_{align_key}_calc_{calc_key}_{start_frame}_{end_frame}_{step}"
                else:
                    chain_names = "_".join(sorted(protein_chains.keys()))
                    align_key = align_selection if align_selection else "auto"
                    cache_key = f"rmsf_mdtraj_align_{align_key}_calc_{chain_names}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached RMSF results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Align trajectory using MDTraj
            if align_selection:
                logger.info(f"Aligning trajectory on '{align_selection}' using MDTraj")
                align_indices = traj.topology.select(align_selection)

                if len(align_indices) > 0:
                    # Superpose to first frame
                    traj.superpose(traj, frame=0, atom_indices=align_indices)
                    logger.info(f"Alignment completed on {len(align_indices)} atoms")
                else:
                    logger.warning(f"No atoms found for alignment selection: {align_selection}")

            # Calculate RMSF
            if single_selection:
                # Single selection RMSF
                calc_sel = calculate_selection if calculate_selection else "name CA"
                logger.info(f"Calculating RMSF for '{calc_sel}'")

                calc_indices = traj.topology.select(calc_sel)
                if len(calc_indices) == 0:
                    raise ValueError(f"No atoms found for calculation selection: {calc_sel}")

                # Calculate RMSF using MDTraj
                rmsf_values = md.rmsf(traj, traj, frame=0, atom_indices=calc_indices)

                # Get residue IDs
                residue_ids = []
                for idx in calc_indices:
                    atom = traj.topology.atom(idx)
                    residue_ids.append(atom.residue.index)
                residue_ids = np.array(residue_ids)

                results = (rmsf_values, residue_ids)
                logger.info(f"Single selection RMSF completed. Mean: {np.mean(rmsf_values):.2f} Å")
            else:
                # Multi-chain RMSF
                logger.info("Calculating RMSF for multiple chains using MDTraj")
                results = {}

                for chain_name, chain_selection in protein_chains.items():
                    logger.info(f"  Processing {chain_name}: {chain_selection}")

                    calc_indices = traj.topology.select(chain_selection)
                    if len(calc_indices) == 0:
                        logger.warning(f"No atoms found for {chain_name} with selection '{chain_selection}'")
                        continue

                    # Calculate RMSF for this chain
                    rmsf_values = md.rmsf(traj, traj, frame=0, atom_indices=calc_indices)

                    # Get residue IDs for this chain
                    residue_ids = []
                    for idx in calc_indices:
                        atom = traj.topology.atom(idx)
                        residue_ids.append(atom.residue.index)
                    residue_ids = np.array(residue_ids)

                    results[chain_name] = (rmsf_values, residue_ids)
                    logger.info(f"  {chain_name}: {len(rmsf_values)} atoms, Mean RMSF = {np.mean(rmsf_values):.2f} Å")

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(results, f)

            logger.info("RMSF calculation completed successfully")
            return results

        except Exception as e:
            logger.error(f"Error calculating RMSF with MDTraj: {e}")
            raise

    def _calculate_rmsf_mdtraj_single(self, trajectory: md.Trajectory,
                                     calculate_selection: str,
                                     align_selection: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate RMSF for a single selection using MDTraj with OpenMP parallelization.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            MDTraj trajectory object
        calculate_selection : str
            Selection string for RMSF calculation
        align_selection : str
            Selection string for alignment

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            RMSF values and residue/nucleotide IDs
        """
        try:
            # Get atom indices for alignment and calculation
            align_indices = trajectory.topology.select(align_selection)
            calc_indices = trajectory.topology.select(calculate_selection)

            if len(align_indices) == 0:
                raise ValueError(f"No atoms found for alignment selection: {align_selection}")
            if len(calc_indices) == 0:
                raise ValueError(f"No atoms found for calculation selection: {calculate_selection}")

            logger.info(f"Aligning trajectory using {len(align_indices)} atoms")
            logger.info(f"Calculating RMSF for {len(calc_indices)} atoms with {os.environ.get('OMP_NUM_THREADS', 'auto')} CPU threads")

            # Align trajectory
            aligned_traj = trajectory.superpose(trajectory, atom_indices=align_indices)

            # Pre-center trajectory for better performance
            logger.info("Pre-centering trajectory for optimized RMSF calculation")
            centered_traj = aligned_traj.center_coordinates()

            # Calculate RMSF for selected atoms using MDTraj's built-in function with parallelization
            logger.info("Calculating RMSF using MDTraj's optimized parallel implementation")
            rmsf_values = md.rmsf(centered_traj, centered_traj, frame=0, atom_indices=calc_indices, parallel=True)

            # Convert from nm to Angstrom (MDTraj uses nm)
            rmsf_values *= 10.0

            # Calculate per-atom RMSF (average over x,y,z coordinates)
            rmsf_per_atom = np.sqrt(np.mean(rmsf_values.reshape(-1, 3)**2, axis=1))

            # Get residue IDs for the calculation atoms
            residue_ids = np.array([trajectory.topology.atom(i).residue.resSeq for i in calc_indices])

            logger.info(f"RMSF calculation completed. Mean RMSF: {np.mean(rmsf_per_atom):.2f} Å")

            return rmsf_per_atom, residue_ids

        except Exception as e:
            logger.error(f"Error in single RMSF calculation with MDTraj: {e}")
            raise

    def _calculate_rmsf_mdtraj_multi(self, trajectory: md.Trajectory,
                                   chain_selections: Dict[str, str],
                                   align_selection: str) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """
        Calculate RMSF for multiple chains using MDTraj.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            MDTraj trajectory object
        chain_selections : Dict[str, str]
            Dictionary mapping chain names to selection strings
        align_selection : str
            Selection string for alignment

        Returns
        -------
        Dict[str, Tuple[np.ndarray, np.ndarray]]
            Dictionary mapping chain names to (rmsf_values, residue_ids) tuples
        """
        try:
            # Get alignment indices
            align_indices = trajectory.topology.select(align_selection)
            if len(align_indices) == 0:
                raise ValueError(f"No atoms found for alignment selection: {align_selection}")

            logger.info(f"Aligning trajectory using {len(align_indices)} atoms")

            # Align trajectory once
            aligned_traj = trajectory.superpose(trajectory, atom_indices=align_indices)

            results = {}

            for chain_name, chain_selection in chain_selections.items():
                try:
                    logger.info(f"Processing {chain_name}: {chain_selection}")

                    # Get atom indices for this chain
                    calc_indices = trajectory.topology.select(chain_selection)

                    if len(calc_indices) == 0:
                        logger.warning(f"No atoms found for {chain_name} with selection '{chain_selection}'")
                        continue

                    # Calculate RMSF for this chain
                    mean_positions = np.mean(aligned_traj.xyz[:, calc_indices], axis=0)
                    rmsf_values = np.sqrt(np.mean((aligned_traj.xyz[:, calc_indices] - mean_positions)**2, axis=0))
                    # Convert from nm to Angstrom
                    rmsf_values *= 10.0

                    # Calculate per-atom RMSF (average over x,y,z coordinates)
                    rmsf_per_atom = np.sqrt(np.mean(rmsf_values**2, axis=1))

                    # Get residue IDs
                    residue_ids = np.array([trajectory.topology.atom(i).residue.resSeq for i in calc_indices])

                    results[chain_name] = (rmsf_per_atom, residue_ids)
                    logger.info(f"  {chain_name}: {len(rmsf_per_atom)} atoms, Mean RMSF = {np.mean(rmsf_per_atom):.2f} Å")

                except Exception as e:
                    logger.error(f"Error processing {chain_name}: {e}")
                    continue

            return results

        except Exception as e:
            logger.error(f"Error in multi-chain RMSF calculation with MDTraj: {e}")
            raise

    def calculate_multi_trajectory_rmsf(self,
                                      topology_file: str,
                                      trajectory_files: List[str],
                                      auto_detect_chains: bool = True,
                                      align_selection: Optional[str] = None,
                                      start_frame: int = 0,
                                      end_frame: Optional[int] = None,
                                      step: int = 1,
                                      cache_name: Optional[str] = None) -> Dict[str, Dict]:
        """
        Calculate RMSF across multiple trajectories with automatic chain detection.

        This method analyzes multiple trajectory replicates and provides combined
        statistics including mean, standard deviation, and individual trajectory data.

        Parameters
        ----------
        topology_file : str
            Path to topology file
        trajectory_files : List[str]
            List of trajectory file paths
        auto_detect_chains : bool
            Whether to automatically detect chains. Default is True.
        align_selection : str, optional
            Selection string for trajectory alignment. If None, uses primary protein chain.
        start_frame : int, optional
            Starting frame for analysis
        end_frame : int, optional
            Ending frame for analysis
        step : int, optional
            Step size for frame analysis
        cache_name : str, optional
            Custom cache name for combined results

        Returns
        -------
        Dict[str, Dict]
            Multi-trajectory RMSF data structure:
            {
                "Chain X": {
                    "trajectories": {
                        "traj1": (rmsf_values, residue_ids),
                        "traj2": (rmsf_values, residue_ids),
                        ...
                    },
                    "combined": {
                        "mean_rmsf": np.array,
                        "std_rmsf": np.array,
                        "residue_ids": np.array
                    },
                    "statistics": {
                        "overall_mean": float,
                        "overall_std": float,
                        "trajectory_means": List[float]
                    }
                }
            }
        """
        try:
            logger.info(f"Starting multi-trajectory RMSF analysis for {len(trajectory_files)} trajectories")

            # Cache key for combined results
            if cache_name is None:
                traj_names = [Path(f).stem for f in trajectory_files]
                cache_key = f"multi_traj_rmsf_{len(trajectory_files)}trajs_{start_frame}_{end_frame}_{step}"
                cache_key += f"_{'_'.join(traj_names[:3])}"  # Include first 3 traj names
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached multi-trajectory RMSF results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Calculate RMSF for each trajectory
            trajectory_results = {}
            all_chain_names = set()

            for i, traj_file in enumerate(trajectory_files):
                if not Path(traj_file).exists():
                    logger.warning(f"Trajectory file not found: {traj_file}")
                    continue

                traj_name = f"trajectory_{i+1}"
                logger.info(f"Processing {traj_name}: {Path(traj_file).name}")

                try:
                    # Calculate RMSF for this trajectory using our enhanced method
                    rmsf_data = self.calculate_rmsf(
                        universe=topology_file,
                        trajectory=traj_file,
                        calculate_selection=None if auto_detect_chains else align_selection,
                        align_selection=align_selection,
                        use_mdtraj=True,
                        start_frame=start_frame,
                        end_frame=end_frame,
                        step=step,
                        cache_name=f"{cache_key}_{traj_name}"
                    )

                    if isinstance(rmsf_data, dict):
                        trajectory_results[traj_name] = rmsf_data
                        all_chain_names.update(rmsf_data.keys())
                        logger.info(f"  ✓ {traj_name}: {len(rmsf_data)} chains detected")
                    else:
                        # Single selection result - wrap in dict
                        trajectory_results[traj_name] = {"default_chain": rmsf_data}
                        all_chain_names.add("default_chain")
                        logger.info(f"  ✓ {traj_name}: single chain analysis")

                except Exception as e:
                    logger.error(f"Error processing {traj_name}: {e}")
                    continue

            if not trajectory_results:
                logger.error("No trajectories processed successfully")
                return {}

            logger.info(f"Processed {len(trajectory_results)} trajectories, found chains: {sorted(all_chain_names)}")

            # Combine results across trajectories
            combined_results = {}

            for chain_name in sorted(all_chain_names):
                logger.info(f"Combining data for {chain_name}")

                # Collect data from all trajectories for this chain
                chain_data = {
                    "trajectories": {},
                    "combined": {},
                    "statistics": {}
                }

                traj_rmsf_values = []
                traj_means = []
                reference_residue_ids = None

                for traj_name, traj_results in trajectory_results.items():
                    if chain_name in traj_results:
                        rmsf_values, residue_ids = traj_results[chain_name]

                        # Store individual trajectory data
                        chain_data["trajectories"][traj_name] = (rmsf_values, residue_ids)

                        # Collect for averaging
                        traj_rmsf_values.append(rmsf_values)
                        traj_means.append(np.mean(rmsf_values))

                        # Use first trajectory's residue IDs as reference
                        if reference_residue_ids is None:
                            reference_residue_ids = residue_ids
                        elif not np.array_equal(residue_ids, reference_residue_ids):
                            logger.warning(f"Residue ID mismatch for {chain_name} in {traj_name}")

                if traj_rmsf_values:
                    # Calculate combined statistics
                    traj_rmsf_array = np.array(traj_rmsf_values)
                    mean_rmsf = np.mean(traj_rmsf_array, axis=0)
                    std_rmsf = np.std(traj_rmsf_array, axis=0)

                    chain_data["combined"] = {
                        "mean_rmsf": mean_rmsf,
                        "std_rmsf": std_rmsf,
                        "residue_ids": reference_residue_ids
                    }

                    chain_data["statistics"] = {
                        "overall_mean": float(np.mean(mean_rmsf)),
                        "overall_std": float(np.mean(std_rmsf)),
                        "trajectory_means": traj_means,
                        "n_trajectories": len(traj_rmsf_values),
                        "n_residues": len(mean_rmsf)
                    }

                    combined_results[chain_name] = chain_data

                    logger.info(f"  ✓ {chain_name}: {len(traj_rmsf_values)} trajectories, "
                              f"mean RMSF = {chain_data['statistics']['overall_mean']:.2f} ± "
                              f"{chain_data['statistics']['overall_std']:.2f} Å")

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(combined_results, f)

            logger.info(f"Multi-trajectory RMSF analysis completed for {len(combined_results)} chains")
            return combined_results

        except Exception as e:
            logger.error(f"Error in multi-trajectory RMSF calculation: {e}")
            raise