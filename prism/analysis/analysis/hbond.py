import logging
import numpy as np
from typing import Dict, List, Tuple, Optional
import mdtraj as md

logger = logging.getLogger(__name__)

class HydrogenBondAnalyzer:
    """Analyzer for hydrogen bonds"""
    
    def __init__(self, trajectory_manager, config):
        self.traj = trajectory_manager.traj
        self.protein_donors = trajectory_manager.protein_donors
        self.protein_acceptors = trajectory_manager.protein_acceptors
        self.ligand_donors = trajectory_manager.ligand_donors
        self.ligand_acceptors = trajectory_manager.ligand_acceptors
        self.config = config
    
    def analyze_hydrogen_bonds(self, parallel: bool = True) -> Tuple[List[Tuple[str, float, float, float]], Dict]:
        """Analyze hydrogen bonds between protein and ligand"""
        if self.traj is None:
            raise ValueError("No trajectory provided")
        
        hbond_pairs, hbond_triplets = self._prepare_hbond_pairs()
        
        if not hbond_pairs:
            logger.warning("No hydrogen bond pairs found")
            return [], {}
        
        total_frames = self.traj.n_frames
        hbond_stats = {}
        
        hbond_stats = self._analyze_hbonds_serial(hbond_pairs, hbond_triplets, hbond_stats)
        
        return self._calculate_hbond_frequencies(hbond_stats, total_frames)
    
    def _prepare_hbond_pairs(self) -> Tuple[List, List]:
        """Prepare hydrogen bond donor-acceptor pairs"""
        hbond_pairs = []
        hbond_triplets = []
        
        self._add_hbond_pairs(
            self.protein_donors, self.ligand_acceptors, hbond_pairs, hbond_triplets, 
            donor_is_protein=True
        )
        
        self._add_hbond_pairs(
            self.ligand_donors, self.protein_acceptors, hbond_pairs, hbond_triplets,
            donor_is_protein=False
        )
        
        return hbond_pairs, hbond_triplets
    
    def _add_hbond_pairs(self, donors: List[int], acceptors: List[int], 
                        hbond_pairs: List, hbond_triplets: List, donor_is_protein: bool):
        """Add hydrogen bond pairs for analysis"""
        frame0 = self.traj.xyz[0]
        
        for donor_idx in donors:
            donor_atom = self.traj.topology.atom(donor_idx)
            residue_atoms = list(donor_atom.residue.atoms)
            
            hydrogens = []
            for atom in residue_atoms:
                if atom.element.symbol == 'H':
                    dist = np.linalg.norm(frame0[donor_idx] - frame0[atom.index])
                    if dist < self.config.bond_length_threshold_nm:
                        hydrogens.append(atom.index)
            
            if hydrogens:
                hydrogen_idx = hydrogens[0]
                for acceptor_idx in acceptors:
                    hbond_pairs.append([donor_idx, acceptor_idx])
                    hbond_triplets.append([hydrogen_idx, donor_idx, acceptor_idx])
    
    def _analyze_hbonds_serial(self, hbond_pairs: List, hbond_triplets: List, hbond_stats: Dict) -> Dict:
        """Analyze hydrogen bonds serially"""
        for frame_idx in range(self.traj.n_frames):
            frame_stats = self._compute_hbonds_for_frame(frame_idx, hbond_pairs, hbond_triplets)
            for key, data in frame_stats.items():
                if key not in hbond_stats:
                    hbond_stats[key] = {'count': 0, 'distance': [], 'angle': []}
                hbond_stats[key]['count'] += data['count']
                hbond_stats[key]['distance'].extend(data['distance'])
                hbond_stats[key]['angle'].extend(data['angle'])
        return hbond_stats
    
    def _compute_hbonds_for_frame(self, frame_idx: int, hbond_pairs: List, hbond_triplets: List) -> Dict:
        """Compute hydrogen bonds for a single frame"""
        frame_stats = {}
        try:
            dists = md.compute_distances(self.traj[frame_idx], hbond_pairs, periodic=False)
            angles = md.compute_angles(self.traj[frame_idx], hbond_triplets, periodic=False)
            
            for idx, pair in enumerate(hbond_pairs):
                dist = dists[0][idx]
                angle = np.degrees(angles[0][idx])
                
                if (dist < self.config.hbond_distance_cutoff_nm and 
                    angle > self.config.hbond_angle_cutoff_deg):
                    
                    key = self._create_hbond_key(pair)
                    if key not in frame_stats:
                        frame_stats[key] = {'count': 0, 'distance': [], 'angle': []}
                    
                    frame_stats[key]['count'] += 1
                    frame_stats[key]['distance'].append(dist * 10.0)  # nm to angstrom
                    frame_stats[key]['angle'].append(angle)
        except Exception as e:
            logger.warning(f"Error processing frame {frame_idx}: {e}")
        return frame_stats
    
    def _create_hbond_key(self, pair: List[int]) -> str:
        """Create hydrogen bond identifier string"""
        donor_idx, acceptor_idx = pair
        donor_atom = self.traj.topology.atom(donor_idx)
        acceptor_atom = self.traj.topology.atom(acceptor_idx)
        
        if donor_atom.residue.is_protein:
            return (f"{donor_atom.residue.name} {donor_atom.residue.resSeq} "
                   f"({donor_atom.name}) → Ligand")
        else:
            return (f"Ligand → {acceptor_atom.residue.name} "
                   f"{acceptor_atom.residue.resSeq} ({acceptor_atom.name})")
    
    def _calculate_hbond_frequencies(self, hbond_stats: Dict, total_frames: int) -> Tuple[List[Tuple[str, float, float, float]], Dict]:
        """Calculate hydrogen bond frequencies"""
        hbond_frequencies = []
        for key, data in hbond_stats.items():
            freq = data['count'] / total_frames
            avg_dist = np.mean(data['distance']) if data['distance'] else 0
            avg_angle = np.mean(data['angle']) if data['angle'] else 0
            hbond_frequencies.append((key, freq, avg_dist, avg_angle))
        
        hbond_frequencies.sort(key=lambda x: x[1], reverse=True)
        return hbond_frequencies, hbond_stats