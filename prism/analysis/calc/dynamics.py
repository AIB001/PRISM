#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand dynamics analysis for MD trajectories.

Implements ligand-specific dynamics analysis including RMSD, distance tracking,
and binding pocket volume analysis using MDTraj.
"""

import numpy as np
import mdtraj as md
from typing import List, Union, Optional, Tuple, Dict
import logging
from pathlib import Path

from ..core.config import AnalysisConfig
from ..core.parallel import default_processor

logger = logging.getLogger(__name__)


class LigandDynamicsAnalyzer:
    """Ligand dynamics analysis for protein-ligand systems"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path("./cache")
        self._cache_dir.mkdir(parents=True, exist_ok=True)

    def analyze_ligand_dynamics(self,
                              topology: str,
                              trajectory: Union[str, List[str]],
                              ligand_name: str = "UNK",
                              alignment_selection: str = "protein and name CA",
                              key_residues: Optional[List[str]] = None,
                              step: int = 1) -> Dict:
        """
        Analyze ligand-specific dynamics including RMSD, distances, and pocket volume.

        Parameters
        ----------
        topology : str
            Path to topology file
        trajectory : str or list of str
            Path(s) to trajectory file(s)
        ligand_name : str
            Ligand residue name
        alignment_selection : str
            Selection for trajectory alignment
        key_residues : list of str
            Key residues for distance analysis (e.g., ['D623', 'N691', 'S759'])
        step : int
            Frame step for analysis

        Returns
        -------
        dict
            Dictionary containing dynamics analysis results
        """
        try:
            # Load trajectory with parallel processing configuration
            with default_processor.configure_omp_threads():
                traj = md.load(str(trajectory), top=str(topology))
            if step > 1:
                traj = traj[::step]
            logger.info(f"Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")

            # Align trajectory - convert MDAnalysis syntax to MDTraj syntax
            mdtraj_alignment_sel = alignment_selection.replace("segid PR1", "chainid 0").replace(" and name CA", " and name CA")
            alignment_atoms = traj.topology.select(mdtraj_alignment_sel)
            if len(alignment_atoms) > 0:
                with default_processor.configure_omp_threads():
                    traj.superpose(traj, frame=0, atom_indices=alignment_atoms)
                logger.info(f"Aligned trajectory using {len(alignment_atoms)} atoms")

            # Find ligand atoms
            ligand_atoms = traj.topology.select(f'resname {ligand_name}')
            if len(ligand_atoms) == 0:
                logger.warning(f"No ligand atoms found for {ligand_name}")
                return {}

            logger.info(f"Found ligand: {ligand_name} with {len(ligand_atoms)} atoms")

            # Calculate ligand RMSD
            ligand_rmsd = self._calculate_ligand_rmsd(traj, ligand_atoms)

            # Calculate key distances if residues provided
            distance_data = {}
            if key_residues:
                distance_data = self._calculate_key_distances(traj, ligand_atoms, key_residues)

            # Calculate binding pocket volume
            pocket_volumes = self._calculate_pocket_volume(traj, key_residues or [])

            return {
                'time': traj.time,
                'ligand_rmsd': ligand_rmsd,
                'distances': distance_data,
                'pocket_volumes': pocket_volumes,
                'ligand_name': ligand_name,
                'n_frames': traj.n_frames
            }

        except Exception as e:
            logger.error(f"Error in ligand dynamics analysis: {e}")
            return {}

    def _calculate_ligand_rmsd(self, traj: md.Trajectory, ligand_atoms: np.ndarray) -> np.ndarray:
        """Calculate ligand RMSD relative to first frame"""
        try:
            ligand_coords = traj.xyz[:, ligand_atoms, :]
            ref_coords = ligand_coords[0]

            rmsd_values = []
            for frame_coords in ligand_coords:
                rmsd = np.sqrt(np.mean((frame_coords - ref_coords)**2)) * 10  # Convert to Angstrom
                rmsd_values.append(rmsd)

            return np.array(rmsd_values)
        except Exception as e:
            logger.error(f"Error calculating ligand RMSD: {e}")
            return np.array([])

    def _calculate_key_distances(self, traj: md.Trajectory,
                               ligand_atoms: np.ndarray,
                               key_residues: List[str]) -> Dict[str, np.ndarray]:
        """Calculate distances from ligand center of mass to key residues"""
        distance_data = {}

        try:
            # Calculate ligand center of mass with parallel processing
            with default_processor.configure_omp_threads():
                ligand_com = md.compute_center_of_mass(traj.atom_slice(ligand_atoms))

            for res in key_residues:
                try:
                    # Try different selection formats
                    res_atoms = []

                    # Try by residue name
                    test_atoms = traj.topology.select(f'resname {res}')
                    if len(test_atoms) > 0:
                        res_atoms = test_atoms
                    else:
                        # Try by residue number (extract number part)
                        res_num = ''.join(filter(str.isdigit, res))
                        if res_num:
                            test_atoms = traj.topology.select(f'resid {res_num}')
                            if len(test_atoms) > 0:
                                res_atoms = test_atoms

                    if len(res_atoms) > 0:
                        # Calculate residue center of mass with parallel processing
                        with default_processor.configure_omp_threads():
                            res_com = md.compute_center_of_mass(traj.atom_slice(res_atoms))
                        distances = np.linalg.norm(ligand_com - res_com, axis=1)
                        distance_data[res] = distances
                        logger.debug(f"Distance to {res}: {np.mean(distances):.2f} ± {np.std(distances):.2f} nm")

                except Exception as e:
                    logger.warning(f"Could not calculate distance to {res}: {e}")

        except Exception as e:
            logger.error(f"Error calculating key distances: {e}")

        return distance_data

    def _calculate_pocket_volume(self, traj: md.Trajectory,
                               pocket_residues: List[str]) -> np.ndarray:
        """Calculate approximate binding pocket volume variation"""
        try:
            pocket_atoms = []

            for res in pocket_residues:
                try:
                    # Try different selection methods
                    res_atoms = traj.topology.select(f'resname {res}')
                    if len(res_atoms) == 0:
                        res_num = ''.join(filter(str.isdigit, res))
                        if res_num:
                            res_atoms = traj.topology.select(f'resid {res_num}')

                    if len(res_atoms) > 0:
                        pocket_atoms.extend(res_atoms)
                except:
                    continue

            if len(pocket_atoms) > 3:
                pocket_atoms = np.array(pocket_atoms)
                pocket_coords = traj.xyz[:, pocket_atoms, :]

                # Calculate approximate pocket volume using standard deviation
                pocket_volumes = []
                for frame in range(traj.n_frames):
                    coords = pocket_coords[frame]
                    center = np.mean(coords, axis=0)
                    distances = np.linalg.norm(coords - center, axis=1)
                    volume_approx = np.std(distances) ** 3  # Rough volume estimate
                    pocket_volumes.append(volume_approx)

                return np.array(pocket_volumes)
            else:
                logger.warning("Insufficient pocket atoms for volume calculation")
                return np.zeros(traj.n_frames)

        except Exception as e:
            logger.error(f"Error calculating pocket volume: {e}")
            return np.zeros(getattr(traj, 'n_frames', 0))