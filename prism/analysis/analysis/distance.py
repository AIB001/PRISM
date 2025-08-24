import logging
import numpy as np
from typing import Dict, List, Tuple, Optional

logger = logging.getLogger(__name__)

class DistanceAnalyzer:
    """Analyzer for distance calculations"""
    
    def __init__(self, trajectory_manager, config):
        self.traj = trajectory_manager.traj
        self.ligand_atoms = trajectory_manager.ligand_atoms
        self.protein_residues = trajectory_manager.protein_residues
        self.config = config
    
    def calculate_distance_analysis(self) -> Dict[str, Dict[str, float]]:
        """Calculate distance statistics for all residues"""
        if self.traj is None:
            raise ValueError("No trajectory provided")
        
        distance_stats = {}
        
        for residue in self.protein_residues:
            residue_atoms = [
                a.index for a in residue.atoms 
                if a.element.symbol != 'H'
            ]
            
            if not residue_atoms:
                continue
            
            try:
                distances = self._calculate_min_distances_mdtraj(
                    self.traj, self.ligand_atoms, residue_atoms
                )
                distances_angstrom = distances * 10.0  # nm to angstrom
                
                key = f"{residue.name} {residue.resSeq}"
                distance_stats[key] = {
                    'mean': float(np.mean(distances_angstrom)),
                    'min': float(np.min(distances_angstrom)),
                    'max': float(np.max(distances_angstrom)),
                    'std': float(np.std(distances_angstrom)),
                    'contact_proportion': float(np.mean(distances_angstrom < (self.config.distance_cutoff_nm * 10.0)))
                }
            except Exception as e:
                logger.warning(f"Failed to analyze distance for residue {residue.name} {residue.resSeq}: {e}")
                continue
        
        return distance_stats
    
    @staticmethod
    def _calculate_min_distances_mdtraj(traj, atoms1, atoms2):
        """Calculate minimum distances between two atom groups"""
        try:
            distances = []
            for frame_idx in range(traj.n_frames):
                coords1 = traj.xyz[frame_idx, atoms1]
                coords2 = traj.xyz[frame_idx, atoms2]
                
                diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
                dist_matrix = np.sqrt(np.sum(diff**2, axis=2))
                min_dist = np.min(dist_matrix)
                distances.append(min_dist)
            
            return np.array(distances)
        except Exception as e:
            logger.error(f"Error calculating distances: {e}")
            raise
    
    @staticmethod
    def select_top_contacts(contact_proportions: Dict[str, float], top_n: int = 10) -> List[Tuple[str, float]]:
        """Select top N contacting residues"""
        sorted_contacts = sorted(
            contact_proportions.items(),
            key=lambda x: x[1],
            reverse=True
        )
        return sorted_contacts[:top_n]