import os
import logging
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
import mdtraj as md
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")

from .config import AnalysisConfig

logger = logging.getLogger(__name__)

class TrajectoryManager:
    
    def __init__(self, topology_file: str, trajectory_file: str, config: AnalysisConfig):
        self.config = config
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.traj = None
        
        self.ligand_atoms = None
        self.protein_atoms = None
        self.protein_residues = None
        self.ligand_residue = None
        
        self.protein_donors = []
        self.protein_acceptors = []
        self.ligand_donors = []
        self.ligand_acceptors = []
        
        self.load_trajectory()
        self.prepare_system()
    
    def load_trajectory(self):
        try:
            if not os.path.exists(self.topology_file):
                raise FileNotFoundError(f"Topology file not found: {self.topology_file}")
            if not os.path.exists(self.trajectory_file):
                raise FileNotFoundError(f"Trajectory file not found: {self.trajectory_file}")
            
            self.traj = md.load(self.trajectory_file, top=self.topology_file)
            if self.traj.n_frames == 0:
                raise ValueError("Trajectory has no frames")
            
            return self.traj
        except Exception as e:
            logger.error(f"Failed to load trajectory: {e}")
            raise
    
    def prepare_system(self):
        try:
            self._select_atoms()
            self._prepare_hbond_groups()
            
            logger.info(f"Found {len(self.ligand_atoms)} ligand atoms, {len(self.protein_atoms)} protein atoms")
        except Exception as e:
            logger.error(f"Failed to prepare system: {e}")
            raise
    
    def _select_atoms(self):
        self.ligand_atoms = self.traj.topology.select(
            f"resname {self.config.ligand_resname} and element != 'H'"
        )
        self.protein_atoms = self.traj.topology.select("protein and element != 'H'")
        
        self.protein_residues = [
            res for res in self.traj.topology.residues 
            if res.is_protein and res.name != 'HOH'
        ]
        
        self.ligand_residue = next(
            (res for res in self.traj.topology.residues 
             if res.name == self.config.ligand_resname), None
        )
    
    def _prepare_hbond_groups(self):
        try:
            self.protein_donors = []
            self.protein_acceptors = []
            self.ligand_donors = []
            self.ligand_acceptors = []
            
            frame0 = self.traj.xyz[0]
            self._find_donors_acceptors(frame0, self.protein_residues, is_protein=True)
            if self.ligand_residue:
                self._find_donors_acceptors(frame0, [self.ligand_residue], is_protein=False)
            
            logger.info(f"Protein donors: {len(self.protein_donors)}, acceptors: {len(self.protein_acceptors)}")
            logger.info(f"Ligand donors: {len(self.ligand_donors)}, acceptors: {len(self.ligand_acceptors)}")
        except Exception as e:
            logger.error(f"Failed to prepare hydrogen bond analysis: {e}")
            raise
    
    def _find_donors_acceptors(self, frame0: np.ndarray, residues: List, is_protein: bool):
        if is_protein:
            donor_list = self.protein_donors
            acceptor_list = self.protein_acceptors
        else:
            donor_list = self.ligand_donors
            acceptor_list = self.ligand_acceptors
        
        for residue in residues:
            residue_atoms = list(residue.atoms)
            
            for atom in residue_atoms:
                if atom.element.symbol in ['N', 'O']:
                    is_donor = False
                    for h_atom in residue_atoms:
                        if h_atom.element.symbol == 'H':
                            dist = np.linalg.norm(
                                frame0[atom.index] - frame0[h_atom.index]
                            )
                            if dist < self.config.bond_length_threshold_nm:
                                donor_list.append(atom.index)
                                is_donor = True
                                break
                    
                    acceptor_list.append(atom.index)
    
    def get_system_info(self) -> Dict[str, Any]:
        if self.traj is None:
            return {}
        
        return {
            'trajectory_info': {
                'n_frames': int(self.traj.n_frames),
                'n_atoms': int(self.traj.n_atoms),
                'time_step_ps': float(self.traj.timestep),
                'total_time_ns': float(self.traj.time[-1] / 1000),
                'box_dimensions': self.traj.unitcell_lengths[0].tolist() if self.traj.unitcell_lengths is not None else None
            },
            'ligand_info': {
                'residue_name': self.config.ligand_resname,
                'n_atoms': len(self.ligand_atoms),
                'atom_indices': self.ligand_atoms.tolist()
            },
            'protein_info': {
                'n_atoms': len(self.protein_atoms),
                'n_residues': len(self.protein_residues),
                'residue_names': [res.name for res in self.protein_residues],
                'residue_numbers': [res.resSeq for res in self.protein_residues]
            },
            'hbond_info': {
                'protein_donors': len(self.protein_donors),
                'protein_acceptors': len(self.protein_acceptors),
                'ligand_donors': len(self.ligand_donors),
                'ligand_acceptors': len(self.ligand_acceptors)
            }
        }
    
    def validate_system(self) -> bool:
        if self.traj is None:
            logger.error("No trajectory loaded")
            return False
        
        if len(self.ligand_atoms) == 0:
            logger.error(f"No ligand atoms found for resname '{self.config.ligand_resname}'")
            return False
        
        if len(self.protein_atoms) == 0:
            logger.error("No protein atoms found")
            return False
        
        return True