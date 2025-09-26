import os
import logging
import numpy as np
import MDAnalysis as mda
from typing import Tuple
from .config import AnalysisConfig
# DistanceCalculator removed - using MDTraj directly

logger = logging.getLogger(__name__)

class MultiSystemAnalyzer:
    """Multi-system analysis module with fixed implementation"""
    
    def __init__(self, config: AnalysisConfig):
        self.config = config
  
    def calculate_unique_contacts(self, universe, sel1: str, sel2: str, 
                                 cutoff_nm: float) -> np.ndarray:
        """Calculate cumulative unique contacts over trajectory"""
        try:
            group1 = universe.select_atoms(sel1)
            group2 = universe.select_atoms(sel2)
            
            if len(group1) == 0 or len(group2) == 0:
                logger.warning(f"Empty selection: {sel1}={len(group1)}, {sel2}={len(group2)}")
                return np.array([])
            
            seen_contacts = set()
            cumulative_contacts = []
            cutoff_angstrom = cutoff_nm * 10.0

            for ts in universe.trajectory:
                try:
                    # TODO: Replace with MDTraj implementation
                    # dist_matrix = DistanceCalculator.calculate_distance_array_mda(
                        group1, group2, universe.dimensions
                    )
                    contact_indices = np.where(dist_matrix < cutoff_angstrom)

                    # Use atom-residue contacts instead of just atom-atom
                    frame_contacts = {
                        (group1[i].index, group2[j].resid)
                        for i, j in zip(*contact_indices)
                    }

                    seen_contacts.update(frame_contacts)
                    cumulative_contacts.append(len(seen_contacts))
                except Exception as e:
                    logger.warning(f"Error in frame {ts.frame}: {e}")
                    cumulative_contacts.append(len(seen_contacts))

            return np.array(cumulative_contacts)
        except Exception as e:
            logger.error(f"Error calculating unique contacts: {e}")
            return np.array([])
    
    def track_initial_contacts(self, universe, sel1: str, sel2: str, cutoff_nm: float) -> np.ndarray:
        """Track how many initial contacts are retained over time"""
        try:
            group1 = universe.select_atoms(sel1)
            group2 = universe.select_atoms(sel2)
            
            if len(group1) == 0 or len(group2) == 0:
                logger.warning(f"Empty selection: {sel1}={len(group1)}, {sel2}={len(group2)}")
                return np.array([])
            
            # Get initial contacts from first frame
            universe.trajectory[0]
            cutoff_angstrom = cutoff_nm *10.0
            
            # TODO: Replace with MDTraj implementation
            # dist_matrix = DistanceCalculator.calculate_distance_array_mda(
                group1, group2, universe.dimensions
            )
            contact_indices = np.where(dist_matrix < cutoff_angstrom)
            initial_contacts = {
                (group1[i].index, group2[j].resid)
                for i, j in zip(*contact_indices)
            }
            
            if len(initial_contacts) == 0:
                logger.warning("No initial contacts found")
                return np.array([0] * len(universe.trajectory))
            
            remaining_contacts = []
            for ts in universe.trajectory:
                try:
                    # TODO: Replace with MDTraj implementation
                    # dist_matrix = DistanceCalculator.calculate_distance_array_mda(
                        group1, group2, universe.dimensions
                    )
                    contact_indices = np.where(dist_matrix < cutoff_angstrom)
                    current_contacts = {
                        (group1[i].index, group2[j].resid)
                        for i, j in zip(*contact_indices)
                    }
                    remaining = sum(1 for contact in initial_contacts if contact in current_contacts)
                    remaining_contacts.append(remaining)
                except Exception as e:
                    logger.warning(f"Error in frame {ts.frame}: {e}")
                    remaining_contacts.append(0)
            
            return np.array(remaining_contacts)
        except Exception as e:
            logger.error(f"Error tracking initial contacts: {e}")
            return np.array([])
    
    def extract_molecular_weight(self, prod_dir: str, name: str, ligand_resname: str) -> float:
        """Extract molecular weight from topology file"""
        tpr_path = os.path.join(prod_dir, 'md.tpr')
        if not os.path.exists(tpr_path):
            return 0.0
        
        try:
            u = mda.Universe(tpr_path)
            ligand = u.select_atoms(f"resname {ligand_resname}")
            if len(ligand) == 0:
                # Try to find non-standard residues
                ligand = u.select_atoms("not (protein or resname SOL HOH WAT NA CL K MG CA ZN)")
                if len(ligand) == 0:
                    logger.warning(f"No ligand found for {name}")
                    return 0.0
            
            return float(ligand.masses.sum()) if len(ligand) > 0 else 0.0
        except Exception as e:
            logger.warning(f"Failed to extract MW for {name}: {e}")
            return 0.0
    
    def process_single_system(self, name: str, traj_path: str, prod_dir: str, top_path: str, 
                              analysis_type: str, selection1: str, selection2: str, 
                              cutoff_nm: float, timestep: float, ligand_resname: str) -> Tuple[str, np.ndarray, np.ndarray, float]:
        """Process a single system for multi-system analysis"""
        try:
            logger.warning(f"Processing system: {name}")
            u = mda.Universe(top_path, traj_path)
            
            if analysis_type == "overview":
                contacts = self.calculate_unique_contacts(u, selection1, selection2, cutoff_nm)
            else:  # decay
                contacts = self.track_initial_contacts(u, selection1, selection2, cutoff_nm)
            
            times = np.arange(len(contacts)) * timestep
            weight = self.extract_molecular_weight(prod_dir, name, ligand_resname)

            return name, contacts, times, weight
        except Exception as e:
            logger.error(f"Analysis failed for {name}: {e}")
            return name, np.array([]), np.array([]), 0.0
