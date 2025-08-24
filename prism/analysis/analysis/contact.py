import logging
import numpy as np
from typing import Dict, List, Tuple, Optional
import mdtraj as md
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)

class ContactAnalyzer:
    """Analyzer for protein-ligand contacts"""
    
    def __init__(self, trajectory_manager, config):
        self.traj = trajectory_manager.traj
        self.ligand_atoms = trajectory_manager.ligand_atoms
        self.protein_residues = trajectory_manager.protein_residues
        self.protein_atoms = trajectory_manager.protein_atoms
        self.config = config
    
    def calculate_contact_proportions(self) -> Dict[str, float]:
        """Calculate contact proportions for all residues"""
        if self.traj is None:
            raise ValueError("No trajectory provided")
        
        contact_proportions = {}
        
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
                contact_states = self._detect_contacts(distances, self.config)
                contact_proportion = np.mean(contact_states)
                
                key = f"{residue.name} {residue.resSeq}"
                contact_proportions[key] = contact_proportion
            except Exception as e:
                logger.warning(f"Failed to analyze residue {residue.name} {residue.resSeq}: {e}")
                continue
        
        return contact_proportions
    
    def analyze_residue_contact(self, residue_id: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Analyze contact for a specific residue"""
        if self.traj is None:
            raise ValueError("No trajectory provided")
        
        try:
            res_name, res_seq = residue_id.split()
            res_seq = int(res_seq)
        except ValueError:
            logger.error(f"Invalid residue identifier: {residue_id}")
            return None, None
        
        target_res = next(
            (r for r in self.protein_residues 
             if r.name == res_name and r.resSeq == res_seq),
            None
        )
        
        if not target_res:
            logger.error(f"Residue not found: {residue_id}")
            return None, None
        
        residue_atoms = [
            a.index for a in target_res.atoms 
            if a.element.symbol != 'H'
        ]
        
        try:
            distances = self._calculate_min_distances_mdtraj(
                self.traj, self.ligand_atoms, residue_atoms
            )
            contact_states = self._detect_contacts(distances, self.config)
            return distances, contact_states
        except Exception as e:
            logger.error(f"Failed to analyze residue contact: {e}")
            return None, None
    
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
    def _detect_contacts(distances: np.ndarray, config) -> np.ndarray:
        """Detect contact states using hysteresis"""
        distances = np.asarray(distances).flatten()
        contact_states = np.zeros(len(distances), dtype=int)
        in_contact = False
        
        enter_threshold = config.contact_enter_threshold_nm
        exit_threshold = config.contact_exit_threshold_nm
        
        for i, dist in enumerate(distances):
            if not in_contact and dist < enter_threshold:
                in_contact = True
            elif in_contact and dist > exit_threshold:
                in_contact = False
            contact_states[i] = 1 if in_contact else 0
        
        return contact_states
    
    def analyze_atomic_contacts(self, distance_cutoff: float) -> Tuple[np.ndarray, Dict, Dict]:
        """Analyze atomic-level contacts for visualization"""
        contact_counts = np.zeros(len(self.ligand_atoms))
        residue_contacts = {}  # {(ligand_idx, residue_info): count}
        residue_atom_details = {}  # {residue_info: {ligand_idx: count}}
        
        for frame_idx in range(self.traj.n_frames):
            lig_coords = self.traj.xyz[frame_idx, self.ligand_atoms]
            prot_coords = self.traj.xyz[frame_idx, self.protein_atoms]
            distances = cdist(lig_coords, prot_coords)
            
            contact_mask = distances < distance_cutoff
            lig_indices, prot_indices = np.where(contact_mask)
            
            for lig_idx, prot_idx in zip(lig_indices, prot_indices):
                prot_atom = self.protein_atoms[prot_idx]
                atom_obj = self.traj.topology.atom(prot_atom)
                residue = atom_obj.residue
                chain_id = residue.chain.index
                residue_info = f"{residue.name}{residue.resSeq}({chr(65+chain_id)})"
                
                key = (lig_idx, residue_info)
                residue_contacts[key] = residue_contacts.get(key, 0) + 1
                contact_counts[lig_idx] += 1
                
                if residue_info not in residue_atom_details:
                    residue_atom_details[residue_info] = {}
                if lig_idx not in residue_atom_details[residue_info]:
                    residue_atom_details[residue_info][lig_idx] = 0
                residue_atom_details[residue_info][lig_idx] += 1
        
        contact_frequencies = contact_counts / self.traj.n_frames
        
        return contact_frequencies, residue_contacts, residue_atom_details