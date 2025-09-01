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
        self.traj = None
        self.ligand_residue = None
          
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
        
    def _prepare_hbond_analysis(self):
        try:
            if not hasattr(self, 'traj') or self.traj is None:
                raise ValueError("No trajectory data available")
            
            frame0 = self.traj.xyz[0]

            if not self.protein_residues:
                logger.warning("No protein residues found for hydrogen bond analysis")
                return

            # OPTIMIZATION: Only find donors/acceptors near the ligand
            self._find_nearby_donors_acceptors(frame0)
        
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
    
    def _find_nearby_donors_acceptors(self, frame0: np.ndarray):
        """OPTIMIZED: Only find donors/acceptors within reasonable distance of ligand"""
        # Get ligand center
        ligand_atoms = list(self.ligand_residue.atoms) if self.ligand_residue else []
        if not ligand_atoms:
            # Fallback to original method if no ligand
            self._find_donors_acceptors(frame0, self.protein_residues, is_protein=True)
            if self.ligand_residue:
                self._find_donors_acceptors(frame0, [self.ligand_residue], is_protein=False)
            return
            
        ligand_positions = [frame0[atom.index] for atom in ligand_atoms]
        ligand_center = np.mean(ligand_positions, axis=0)
        
        # OPTIMIZATION: Only consider residues within 1.0 nm of ligand
        DISTANCE_CUTOFF = 1.0  # nm
        
        nearby_residues = []
        for residue in self.protein_residues:
            residue_atoms = list(residue.atoms)
            if residue_atoms:
                residue_center = np.mean([frame0[atom.index] for atom in residue_atoms], axis=0)
                dist = np.linalg.norm(residue_center - ligand_center)
                if dist < DISTANCE_CUTOFF:
                    nearby_residues.append(residue)
        
        logger.warning(f"Filtering: {len(nearby_residues)} nearby residues out of {len(self.protein_residues)} total")
        
        # Find donors/acceptors only in nearby residues
        self._find_donors_acceptors(frame0, nearby_residues, is_protein=True)
        
        # Find ligand donors/acceptors
        if self.ligand_residue:
            self._find_donors_acceptors(frame0, [self.ligand_residue], is_protein=False)
    
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
        """OPTIMIZED hydrogen bond analysis with frame skipping"""
        if not hasattr(self, 'traj') or self.traj is None:
            logger.warning("No MDTraj trajectory available for hydrogen bond analysis")
            return []
        
        hbond_pairs, hbond_triplets = self._prepare_hbond_pairs()
        
        if not hbond_pairs:
            logger.warning("No hydrogen bond pairs found")
            return []
        
        # OPTIMIZATION: Limit number of pairs analyzed
        MAX_PAIRS = 500  # Limit to 500 most likely pairs
        if len(hbond_pairs) > MAX_PAIRS:
            logger.warning(f"Limiting analysis to {MAX_PAIRS} pairs (from {len(hbond_pairs)})")
            # Prioritize pairs by initial distance
            frame0 = self.traj.xyz[0]
            pair_distances = []
            for i, pair in enumerate(hbond_pairs):
                dist = np.linalg.norm(frame0[pair[0]] - frame0[pair[1]])
                pair_distances.append((i, dist))
            
            # Sort by distance and keep closest pairs
            pair_distances.sort(key=lambda x: x[1])
            selected_indices = [idx for idx, _ in pair_distances[:MAX_PAIRS]]
            
            hbond_pairs = [hbond_pairs[i] for i in selected_indices]
            hbond_triplets = [hbond_triplets[i] for i in selected_indices]
        
        total_frames = self.traj.n_frames
        hbond_stats = {}
        
        # Use optimized analysis with frame skipping
        hbond_stats = self._analyze_hbonds_fast(hbond_pairs, hbond_triplets, hbond_stats)
        
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
        
        # OPTIMIZATION: Pre-filter by distance in first frame
        MAX_INITIAL_DISTANCE = 0.5  # nm
        
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
                        # OPTIMIZATION: Check initial distance
                        initial_dist = np.linalg.norm(frame0[donor_idx] - frame0[acceptor_idx])
                        if initial_dist < MAX_INITIAL_DISTANCE:
                            hbond_pairs.append([donor_idx, acceptor_idx])
                            hbond_triplets.append([hydrogen_idx, donor_idx, acceptor_idx])
            except (AttributeError, IndexError):
                continue
    
    def _analyze_hbonds_fast(self, hbond_pairs: List, hbond_triplets: List,
                             hbond_stats: Dict) -> Dict:
        """FAST analysis using frame skipping and vectorization"""
        n_frames = self.traj.n_frames
        n_pairs = len(hbond_pairs)
        
        if n_pairs == 0:
            return hbond_stats
        
        logger.warning(f"Analyzing {n_pairs} potential H-bond pairs across {n_frames} frames...")
        
        # OPTIMIZATION: Skip frames for faster analysis
        FRAME_SKIP = max(1, n_frames // 100)  # Analyze ~100 frames max
        frames_to_analyze = list(range(0, n_frames, FRAME_SKIP))
        n_analyzed_frames = len(frames_to_analyze)
        
        logger.warning(f"  Analyzing every {FRAME_SKIP} frames ({n_analyzed_frames} total frames)")
        
        distance_cutoff = self.config.hbond_distance_cutoff_nm
        angle_cutoff = self.config.hbond_angle_cutoff_deg
        
        # Initialize statistics
        for idx, pair in enumerate(hbond_pairs):
            key = self._create_hbond_key(pair)
            hbond_stats[key] = {'count': 0, 'distance': [], 'angle': []}
        
        # OPTIMIZATION: Process all pairs at once for each frame
        for frame_idx in frames_to_analyze:
            # Get single frame
            frame = self.traj[frame_idx]
            
            # Compute all distances and angles at once
            try:
                distances = md.compute_distances(frame, hbond_pairs, periodic=False)[0]
                angles = np.degrees(md.compute_angles(frame, hbond_triplets, periodic=False)[0])
                
                # Find H-bonds
                is_hbond = (distances < distance_cutoff) & (angles > angle_cutoff)
                
                # Update statistics only for detected H-bonds
                for idx in np.where(is_hbond)[0]:
                    key = self._create_hbond_key(hbond_pairs[idx])
                    hbond_stats[key]['count'] += FRAME_SKIP  # Account for skipped frames
                    hbond_stats[key]['distance'].append(distances[idx] * 10.0)
                    hbond_stats[key]['angle'].append(angles[idx])
                    
            except Exception as e:
                logger.warning(f"Error processing frame {frame_idx}: {e}")
                continue
        
        # Filter out H-bonds with very low occurrence
        min_occurrence = max(1, int(0.01 * n_frames))  # At least 1% of frames
        filtered_stats = {}
        for key, data in hbond_stats.items():
            if data['count'] >= min_occurrence:
                filtered_stats[key] = data
        
        logger.warning(f"Found {len(filtered_stats)} significant H-bonds")
        
        return filtered_stats
    
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