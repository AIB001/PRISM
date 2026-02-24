#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dihedral angle analysis for MD trajectories.

Includes backbone dihedrals (phi, psi, chi), custom dihedrals, and Ramachandran analysis.
"""

import numpy as np
import MDAnalysis as mda
from typing import List, Union, Optional, Dict
import logging
from pathlib import Path
import pickle
import importlib.util

DIHEDRAL_AVAILABLE = importlib.util.find_spec("MDAnalysis.analysis.dihedrals") is not None

from ..core.config import AnalysisConfig

logger = logging.getLogger(__name__)


class DihedralAnalyzer:
    """Dihedral angle analysis for protein-ligand systems"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path("./cache")
        self._cache_dir.mkdir(parents=True, exist_ok=True)

        if not DIHEDRAL_AVAILABLE:
            logger.warning("MDAnalysis dihedral module not available. " "Some functionality may be limited.")

    def calculate_backbone_dihedrals(
        self,
        universe: Union[mda.Universe, str],
        trajectory: Optional[Union[str, List[str]]] = None,
        selection: str = "protein",
        start_frame: int = 0,
        end_frame: Optional[int] = None,
        step: int = 1,
        cache_name: Optional[str] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Calculate backbone dihedral angles (phi, psi) for protein residues.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for protein. Default is "protein".
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
        dict
            Dictionary with 'phi', 'psi' arrays and 'residue_ids'
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
                cache_key = f"backbone_{selection}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached backbone dihedral results from {cache_file}")
                with open(cache_file, "rb") as f:
                    return pickle.load(f)

            # Get protein atoms
            protein = universe.select_atoms(selection)
            residues = protein.residues

            # Initialize arrays
            n_frames = len(range(start_frame, end_frame, step))
            phi_angles = []
            psi_angles = []
            residue_ids = []

            # Calculate phi and psi for each residue
            for residue in residues:
                try:
                    # Skip first residue (no phi) and last residue (no psi)
                    if residue.resindex == 0:
                        continue
                    if residue.resindex == len(residues) - 1:
                        continue

                    # Get atoms for phi angle (C(i-1) - N(i) - CA(i) - C(i))
                    prev_residue = residues[residue.resindex - 1]
                    try:
                        c_prev = prev_residue.atoms.select_atoms("name C")[0]
                        n_curr = residue.atoms.select_atoms("name N")[0]
                        ca_curr = residue.atoms.select_atoms("name CA")[0]
                        c_curr = residue.atoms.select_atoms("name C")[0]
                    except IndexError:
                        continue  # Skip if atoms are missing

                    # Get atoms for psi angle (N(i) - CA(i) - C(i) - N(i+1))
                    if residue.resindex < len(residues) - 1:
                        next_residue = residues[residue.resindex + 1]
                        try:
                            n_next = next_residue.atoms.select_atoms("name N")[0]
                        except IndexError:
                            continue  # Skip if next N is missing
                    else:
                        continue

                    # Calculate dihedral angles for each frame
                    phi_values = []
                    psi_values = []

                    for ts in universe.trajectory[start_frame:end_frame:step]:
                        # Calculate phi dihedral
                        phi = self._calculate_dihedral(
                            c_prev.position, n_curr.position, ca_curr.position, c_curr.position
                        )
                        phi_values.append(phi)

                        # Calculate psi dihedral
                        psi = self._calculate_dihedral(
                            n_curr.position, ca_curr.position, c_curr.position, n_next.position
                        )
                        psi_values.append(psi)

                    phi_angles.append(phi_values)
                    psi_angles.append(psi_values)
                    residue_ids.append(residue.resid)

                except Exception as e:
                    logger.debug(f"Skipping residue {residue.resid}: {e}")
                    continue

            # Convert to numpy arrays
            phi_angles = np.array(phi_angles).T  # Shape: (n_frames, n_residues)
            psi_angles = np.array(psi_angles).T
            residue_ids = np.array(residue_ids)

            results = {"phi": phi_angles, "psi": psi_angles, "residue_ids": residue_ids}

            # Cache results
            with open(cache_file, "wb") as f:
                pickle.dump(results, f)

            logger.info(f"Backbone dihedral calculation completed for {len(residue_ids)} residues")
            return results

        except Exception as e:
            logger.error(f"Error calculating backbone dihedrals: {e}")
            raise

    def calculate_custom_dihedral(
        self,
        universe: Union[mda.Universe, str],
        atom_indices: List[int],
        trajectory: Optional[Union[str, List[str]]] = None,
        start_frame: int = 0,
        end_frame: Optional[int] = None,
        step: int = 1,
        cache_name: Optional[str] = None,
    ) -> np.ndarray:
        """
        Calculate a custom dihedral angle defined by four atom indices.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        atom_indices : list of int
            Four atom indices defining the dihedral angle
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
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
            Dihedral angles for each frame in degrees
        """
        if len(atom_indices) != 4:
            raise ValueError("Exactly four atom indices required for dihedral calculation")

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
                cache_key = f"custom_dihedral_{'-'.join(map(str, atom_indices))}_{start_frame}_{end_frame}_{step}"
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached custom dihedral results from {cache_file}")
                with open(cache_file, "rb") as f:
                    return pickle.load(f)

            # Get atoms
            atoms = [universe.atoms[idx] for idx in atom_indices]

            # Calculate dihedral for each frame
            dihedral_values = []
            for ts in universe.trajectory[start_frame:end_frame:step]:
                positions = [atom.position for atom in atoms]
                dihedral = self._calculate_dihedral(*positions)
                dihedral_values.append(dihedral)

            dihedral_values = np.array(dihedral_values)

            # Cache results
            with open(cache_file, "wb") as f:
                pickle.dump(dihedral_values, f)

            logger.info(f"Custom dihedral calculation completed. " f"Mean: {np.mean(dihedral_values):.2f}°")
            return dihedral_values

        except Exception as e:
            logger.error(f"Error calculating custom dihedral: {e}")
            raise

    def calculate_chi_angles(
        self,
        universe: Union[mda.Universe, str],
        trajectory: Optional[Union[str, List[str]]] = None,
        selection: str = "protein",
        start_frame: int = 0,
        end_frame: Optional[int] = None,
        step: int = 1,
        cache_name: Optional[str] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Calculate chi1 dihedral angles for protein side chains.

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for protein
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
            Dictionary with 'chi1' angles and 'residue_ids'
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
                cache_key = f"chi1_{selection}_{start_frame}_{end_frame}_{step}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached chi1 results from {cache_file}")
                with open(cache_file, "rb") as f:
                    return pickle.load(f)

            # Get protein atoms
            protein = universe.select_atoms(selection)
            residues = protein.residues

            chi1_angles = []
            residue_ids = []

            # Calculate chi1 for each residue with side chains
            for residue in residues:
                try:
                    # Chi1 is N-CA-CB-CG (for most residues)
                    n_atom = residue.atoms.select_atoms("name N")[0]
                    ca_atom = residue.atoms.select_atoms("name CA")[0]
                    cb_atoms = residue.atoms.select_atoms("name CB")

                    if len(cb_atoms) == 0:
                        continue  # Skip glycine (no CB)

                    cb_atom = cb_atoms[0]

                    # Get CG atom (varies by residue type)
                    cg_atoms = residue.atoms.select_atoms("name CG or name CG1 or name SG or name OG or name OG1")
                    if len(cg_atoms) == 0:
                        continue  # Skip if no suitable CG atom

                    cg_atom = cg_atoms[0]

                    # Calculate chi1 for each frame
                    chi1_values = []
                    for ts in universe.trajectory[start_frame:end_frame:step]:
                        chi1 = self._calculate_dihedral(
                            n_atom.position, ca_atom.position, cb_atom.position, cg_atom.position
                        )
                        chi1_values.append(chi1)

                    chi1_angles.append(chi1_values)
                    residue_ids.append(residue.resid)

                except Exception as e:
                    logger.debug(f"Skipping residue {residue.resid}: {e}")
                    continue

            # Convert to numpy arrays
            chi1_angles = np.array(chi1_angles).T  # Shape: (n_frames, n_residues)
            residue_ids = np.array(residue_ids)

            results = {"chi1": chi1_angles, "residue_ids": residue_ids}

            # Cache results
            with open(cache_file, "wb") as f:
                pickle.dump(results, f)

            logger.info(f"Chi1 calculation completed for {len(residue_ids)} residues")
            return results

        except Exception as e:
            logger.error(f"Error calculating chi1 angles: {e}")
            raise

    def _calculate_dihedral(self, pos1: np.ndarray, pos2: np.ndarray, pos3: np.ndarray, pos4: np.ndarray) -> float:
        """
        Calculate dihedral angle between four positions.

        Parameters
        ----------
        pos1, pos2, pos3, pos4 : np.ndarray
            Positions of the four atoms

        Returns
        -------
        float
            Dihedral angle in degrees
        """
        # Calculate vectors
        b1 = pos2 - pos1
        b2 = pos3 - pos2
        b3 = pos4 - pos3

        # Calculate cross products
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        # Normalize
        n1_norm = n1 / np.linalg.norm(n1)
        n2_norm = n2 / np.linalg.norm(n2)

        # Calculate dihedral angle
        cos_angle = np.dot(n1_norm, n2_norm)
        # Clamp to avoid numerical errors
        cos_angle = np.clip(cos_angle, -1.0, 1.0)

        # Calculate angle
        angle = np.arccos(cos_angle)

        # Determine sign
        sign = np.dot(np.cross(n1_norm, n2_norm), b2 / np.linalg.norm(b2))
        if sign < 0:
            angle = -angle

        # Convert to degrees
        return np.degrees(angle)

    def ramachandran_analysis(
        self,
        universe: Union[mda.Universe, str],
        trajectory: Optional[Union[str, List[str]]] = None,
        selection: str = "protein",
        start_frame: int = 0,
        end_frame: Optional[int] = None,
        step: int = 1,
        cache_name: Optional[str] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Perform Ramachandran plot analysis (phi vs psi angles).

        Parameters
        ----------
        universe : MDAnalysis.Universe or str
            Universe object or path to topology file
        trajectory : str, list of str, optional
            Path(s) to trajectory file(s)
        selection : str, optional
            Selection string for protein
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
            Dictionary with 'phi', 'psi' arrays, 'residue_ids', and statistics
        """
        try:
            # Calculate backbone dihedrals
            backbone_results = self.calculate_backbone_dihedrals(
                universe, trajectory, selection, start_frame, end_frame, step, cache_name
            )

            phi_angles = backbone_results["phi"]
            psi_angles = backbone_results["psi"]
            residue_ids = backbone_results["residue_ids"]

            # Flatten arrays for Ramachandran plot
            phi_flat = phi_angles.flatten()
            psi_flat = psi_angles.flatten()

            # Calculate statistics
            phi_mean = np.mean(phi_flat)
            phi_std = np.std(phi_flat)
            psi_mean = np.mean(psi_flat)
            psi_std = np.std(psi_flat)

            results = {
                "phi": phi_angles,
                "psi": psi_angles,
                "phi_flat": phi_flat,
                "psi_flat": psi_flat,
                "residue_ids": residue_ids,
                "phi_mean": phi_mean,
                "phi_std": phi_std,
                "psi_mean": psi_mean,
                "psi_std": psi_std,
                "n_frames": phi_angles.shape[0],
                "n_residues": phi_angles.shape[1],
            }

            logger.info(
                f"Ramachandran analysis completed. "
                f"Phi: {phi_mean:.1f}° ± {phi_std:.1f}°, "
                f"Psi: {psi_mean:.1f}° ± {psi_std:.1f}°"
            )

            return results

        except Exception as e:
            logger.error(f"Error in Ramachandran analysis: {e}")
            raise
