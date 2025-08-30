import logging
import numpy as np
import mdtraj as md
from typing import List, Tuple, Dict
from .config import AnalysisConfig

logger = logging.getLogger(__name__)

class HBondAnalyzer:
    
    def __init__(self, config: AnalysisConfig):
        self.config = config
        self.protein_donors = []
        self.protein_acceptors = []
        self.ligand_donors = []
        self.protein_residues = []
        self.ligand_acceptors = []
        self.hbond_stats = {}
          
    def set_trajectory_data(self, traj, ligand_residue, protein_residues=None):  
        self.traj = traj
        self.ligand_residue = ligand_residue

        if protein_residues is not None:
            self.protein_residues = protein_residues
        else:
            self.protein_residues = []
            standard_aa = {
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
            }
            for residue in traj.topology.residues:
                if residue.name in standard_aa:
                    self.protein_residues.append(residue)
    
        try:
            self._prepare_hbond_analysis()
            logger.warning(f"Hydrogen bond analysis prepared with {len(self.protein_residues)} protein residues")
        except Exception as e:
            logger.warning(f"Failed to prepare hydrogen bond analysis: {e}")
            logger.warning(f"Ligand residue: {self.ligand_residue}")
            logger.warning(f"Protein residues count: {len(self.protein_residues)}")
            logger.warning(f"Trajectory frames: {traj.n_frames if hasattr(traj, 'n_frames') else 'unknown'}")
        
    def _prepare_hbond_analysis(self):
        try:
            if not hasattr(self, 'traj') or self.traj is None:
                raise ValueError("No trajectory data available")
            
            frame0 = self.traj.xyz[0]

            if not self.protein_residues:
                logger.warning("No protein residues found for hydrogen bond analysis")
                return

            self._find_donors_acceptors(frame0, self.protein_residues, is_protein=True)

            if self.ligand_residue:
                self._find_donors_acceptors(frame0, [self.ligand_residue], is_protein=False)
            else:
                logger.warning("No ligand residue found for hydrogen bond analysis")
        
            logger.warning(f"Found - Protein donors: {len(self.protein_donors)}, acceptors: {len(self.protein_acceptors)}")
            logger.warning(f"Found - Ligand donors: {len(self.ligand_donors)}, acceptors: {len(self.ligand_acceptors)}")

            total_donors = len(self.protein_donors) + len(self.ligand_donors)
            total_acceptors = len(self.protein_acceptors) + len(self.ligand_acceptors)
        
            if total_donors == 0 or total_acceptors == 0:
                logger.warning(f"Insufficient donors ({total_donors}) or acceptors ({total_acceptors}) for hydrogen bond analysis")
            
        except Exception as e:
            logger.error(f"Failed to prepare hydrogen bond analysis: {e}")
            import traceback
            logger.error(traceback.format_exc())
            logger.warning("Hydrogen bond analysis will be skipped")
    
    def _find_donors_acceptors(self, frame0: np.ndarray, residues: List, is_protein: bool):
        donor_list = self.protein_donors if is_protein else self.ligand_donors
        acceptor_list = self.protein_acceptors if is_protein else self.ligand_acceptors
        
        for residue in residues:
            residue_atoms = list(residue.atoms)
            
            for atom in residue_atoms:
                try:
                    if atom.element.symbol in ['N', 'O']:
                        is_donor = False
                        for h_atom in residue_atoms:
                            if h_atom.element.symbol == 'H':
                                dist = np.linalg.norm(
                                    frame0[atom.index] - frame0[h_atom.index]
                                )
                                if dist < self.config.bond_length_threshold_nm:
                                    is_donor = True
                                    break
                        
                        if is_donor:
                            donor_list.append(atom.index)
                        acceptor_list.append(atom.index)
                except AttributeError:
                    continue
    
    def analyze_hydrogen_bonds(self, universe) -> List[Tuple[str, float, float, float]]:
        if not hasattr(self, 'traj') or self.traj is None:
            logger.warning("No MDTraj trajectory available for hydrogen bond analysis")
            return []
        
        hbond_pairs, hbond_triplets = self._prepare_hbond_pairs()
        
        if not hbond_pairs:
            logger.warning("No hydrogen bond pairs found")
            return []
        
        total_frames = self.traj.n_frames
        hbond_stats = {}
        
        hbond_stats = self._analyze_hbonds_serial(hbond_pairs, hbond_triplets, hbond_stats)
        
        return self._calculate_hbond_frequencies(hbond_stats, total_frames)
    
    def _prepare_hbond_pairs(self) -> Tuple[List, List]:
        hbond_pairs = []
        hbond_triplets = []
        
        self._add_hbond_pairs(self.protein_donors, self.ligand_acceptors, 
                             hbond_pairs, hbond_triplets, donor_is_protein=True)
        
        self._add_hbond_pairs(self.ligand_donors, self.protein_acceptors,
                             hbond_pairs, hbond_triplets, donor_is_protein=False)
        
        return hbond_pairs, hbond_triplets
    
    def _add_hbond_pairs(self, donors: List[int], acceptors: List[int], 
                        hbond_pairs: List, hbond_triplets: List, donor_is_protein: bool):
            
        frame0 = self.traj.xyz[0]
        
        for donor_idx in donors:
            try:
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
            except (AttributeError, IndexError):
                continue
    
    def _analyze_hbonds_serial(self, hbond_pairs: List, hbond_triplets: List,
                             hbond_stats: Dict) -> Dict:
        logger.warning(f"Analyzing hydrogen bonds for {self.traj.n_frames} frames...")
        
        frame_step = max(1, self.traj.n_frames // 500)
        
        for frame_idx in range(0, self.traj.n_frames, frame_step):
            try:
                frame_stats = self._compute_hbonds_for_frame(frame_idx, hbond_pairs, hbond_triplets)
                for key, data in frame_stats.items():
                    if key not in hbond_stats:
                        hbond_stats[key] = {'count': 0, 'distance': [], 'angle': []}
                    hbond_stats[key]['count'] += data['count']
                    hbond_stats[key]['distance'].extend(data['distance'])
                    hbond_stats[key]['angle'].extend(data['angle'])
            except Exception as e:
                logger.warning(f"Error processing frame {frame_idx}: {e}")
                continue
                
        return hbond_stats
    
    def _compute_hbonds_for_frame(self, frame_idx: int, hbond_pairs: List, 
                                hbond_triplets: List) -> Dict:
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
                    frame_stats[key]['distance'].append(dist * 10.0)
                    frame_stats[key]['angle'].append(angle)
        except Exception as e:
            logger.warning(f"Error processing frame {frame_idx}: {e}")
        return frame_stats
    
    def _create_hbond_key(self, pair: List[int]) -> str:
        try:
            donor_idx, acceptor_idx = pair
            donor_atom = self.traj.topology.atom(donor_idx)
            acceptor_atom = self.traj.topology.atom(acceptor_idx)
            
            if donor_idx in self.protein_donors:
                return (f"{donor_atom.residue.name} {donor_atom.residue.resSeq} "
                       f"({donor_atom.name}) -> Ligand")
            else:
                return (f"Ligand -> {acceptor_atom.residue.name} "
                       f"{acceptor_atom.residue.resSeq} ({acceptor_atom.name})")
        except:
            return f"Unknown_HBond_{donor_idx}_{acceptor_idx}"
    
    def _calculate_hbond_frequencies(self, hbond_stats: Dict, total_frames: int) -> List[Tuple[str, float, float, float]]:
        hbond_frequencies = []
        for key, data in hbond_stats.items():
            freq = data['count'] / total_frames
            avg_dist = np.mean(data['distance']) if data['distance'] else 0
            avg_angle = np.mean(data['angle']) if data['angle'] else 0
            hbond_frequencies.append((key, freq, avg_dist, avg_angle))
        
        hbond_frequencies.sort(key=lambda x: x[1], reverse=True)
        self.hbond_stats = hbond_stats
        return hbond_frequencies
