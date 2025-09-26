import logging
import numpy as np
from ..config import AnalysisConfig, convert_numpy_types
import mdtraj as md

logger = logging.getLogger(__name__)



class ContactAnalyzer:
    """Contact analysis module"""
    
    def __init__(self, config: AnalysisConfig):
        self.config = config
    
    def identify_ligand_residue(self, traj):
        """Automatically identify ligand residue using PRISM utilities"""
        from ...utils.ligand import identify_ligand_residue
        return identify_ligand_residue(traj)
    
    def get_heavy_atoms(self, residue):
        """Get heavy atom indices from residue"""
        heavy_atoms = []
        for atom in residue.atoms:
            if atom.element.symbol != 'H':
                heavy_atoms.append(atom.index)
        return heavy_atoms
    
    def get_protein_heavy_atoms(self, traj):
        """Get protein heavy atom indices"""
        protein_atoms = []
        standard_aa = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        for atom in traj.topology.atoms:
            if (atom.residue.name in standard_aa and 
                atom.element.symbol in ['C', 'N', 'O', 'S']):
                protein_atoms.append(atom.index)
        
        return protein_atoms
    
    def calculate_contact_proportions(self, traj):
        """Calculate contact proportions using new method"""
        ligand_residue = self.identify_ligand_residue(traj)
        if not ligand_residue:
            raise ValueError("Could not identify ligand residue")
        
        logger.warning(f"Analyzing ligand: {ligand_residue.name}{ligand_residue.resSeq}")
        
        ligand_atoms = self.get_heavy_atoms(ligand_residue)
        protein_atoms = self.get_protein_heavy_atoms(traj)
        
        atom_pairs = []
        pair_info = []
        
        for lig_atom in ligand_atoms:
            for prot_atom in protein_atoms:
                atom_pairs.append([lig_atom, prot_atom])
                residue = traj.topology.atom(prot_atom).residue
                residue_id = f"{residue.name}{residue.resSeq}"
                pair_info.append((lig_atom, prot_atom, residue_id))
        
        n_frames = traj.n_frames
        traj_sample = traj
        frame_indices = np.arange(n_frames)
        
        distances = md.compute_distances(traj_sample, atom_pairs, opt=True)        
        contact_threshold = self.config.contact_enter_threshold_nm
        contacts = distances < contact_threshold
        
        contact_counts = {}
        residue_contacts = {}
        ligand_atom_contacts = {}
        residue_ligand_atoms = {}
        contact_distances = {}
        
        for i, (lig_atom, prot_atom, residue_id) in enumerate(pair_info):
            contact_frames = int(np.sum(contacts[:, i]))
            
            if contact_frames > 0:
                contact_mask = contacts[:, i]
                avg_distance = float(np.mean(distances[contact_mask, i]))
                
                key = (lig_atom, residue_id)
                if key not in contact_counts:
                    contact_counts[key] = 0
                    contact_distances[key] = []
                contact_counts[key] += contact_frames
                contact_distances[key].append(avg_distance)
                
                if residue_id not in residue_contacts:
                    residue_contacts[residue_id] = 0
                    residue_ligand_atoms[residue_id] = {}
                residue_contacts[residue_id] += contact_frames
                
                if lig_atom not in residue_ligand_atoms[residue_id]:
                    residue_ligand_atoms[residue_id][lig_atom] = 0
                residue_ligand_atoms[residue_id][lig_atom] += contact_frames
                
                if lig_atom not in ligand_atom_contacts:
                    ligand_atom_contacts[lig_atom] = 0
                ligand_atom_contacts[lig_atom] += contact_frames
        
        n_frames_analyzed = len(frame_indices)
        contact_frequencies = {}
        avg_contact_distances = {}
        
        for key, count in contact_counts.items():
            freq = float(count / n_frames_analyzed)
            contact_frequencies[key] = freq
            avg_contact_distances[key] = float(np.mean(contact_distances[key]))
        
        residue_proportions = {}
        residue_avg_distances = {}
        residue_best_ligand_atoms = {}
        
        for residue_id, count in residue_contacts.items():
            if residue_id in residue_ligand_atoms:
                best_lig_atom = max(residue_ligand_atoms[residue_id].items(), key=lambda x: x[1])[0]
                residue_best_ligand_atoms[residue_id] = best_lig_atom
            
            ligand_atoms_contacting = set()
            residue_distances = []
            for (lig_atom, res_id), freq in contact_frequencies.items():
                if res_id == residue_id:
                    ligand_atoms_contacting.add(lig_atom)
                    if (lig_atom, res_id) in avg_contact_distances:
                        residue_distances.append(avg_contact_distances[(lig_atom, res_id)])
            
            if ligand_atoms_contacting:
                max_possible = n_frames_analyzed * len(ligand_atoms_contacting)
                proportion = float(count / max_possible)
                residue_proportions[residue_id] = proportion
                if residue_distances:
                    residue_avg_distances[residue_id] = float(np.mean(residue_distances))
        
        logger.warning(f"Found {len(contact_frequencies)} significant atom-residue contacts")
        logger.warning(f"Found {len(residue_proportions)} contacting residues")
        
        result = {
            'contact_frequencies': convert_numpy_types(contact_frequencies),
            'residue_proportions': convert_numpy_types(residue_proportions),
            'residue_avg_distances': convert_numpy_types(residue_avg_distances),
            'residue_best_ligand_atoms': convert_numpy_types(residue_best_ligand_atoms),
            'ligand_atom_contacts': convert_numpy_types(ligand_atom_contacts),
            'ligand_atoms': convert_numpy_types(ligand_atoms),
            'ligand_residue': ligand_residue,
            'total_frames': int(n_frames_analyzed)
        }
        
        return result

    def detect_contacts(self, distances):
        """Detect contact states from distance array"""
        contact_states = np.zeros(len(distances), dtype=int)
        contact_states[distances < self.config.contact_enter_threshold_nm] = 1
        return contact_states

    def calculate_contacts(self, universe, trajectory=None, selection1="", selection2="", cutoff=4.0, step=1, cache_name=None):
        """
        Calculate contacts between two selections using pure MDTraj.

        Parameters
        ----------
        universe : str
            Path to topology file (PDB/GRO)
        trajectory : str
            Path to trajectory file (XTC/DCD/TRR)
        selection1 : str
            First selection (MDTraj selection syntax)
        selection2 : str
            Second selection (MDTraj selection syntax)
        cutoff : float
            Contact cutoff distance in Angstrom
        step : int
            Step size for trajectory frames
        cache_name : str
            Cache name (for compatibility, not used)

        Returns
        -------
        np.ndarray
            Boolean array indicating contacts for each frame
        """
        try:
            # Load trajectory with MDTraj
            if trajectory:
                traj = md.load(trajectory, top=universe)
            else:
                traj = md.load(universe)  # For PDB-only cases

            # Apply frame slicing
            if step > 1:
                traj = traj[::step]

            # Convert cutoff from Angstrom to nm for MDTraj
            cutoff_nm = cutoff / 10.0

            # Select atoms using MDTraj selection syntax
            try:
                atoms1 = traj.topology.select(selection1)
                atoms2 = traj.topology.select(selection2)
            except Exception as e:
                logger.error(f"Selection error: {e}")
                logger.error(f"selection1: '{selection1}', selection2: '{selection2}'")
                return np.array([])

            if len(atoms1) == 0 or len(atoms2) == 0:
                logger.warning(f"Empty selection: group1={len(atoms1)}, group2={len(atoms2)}")
                return np.array([])

            logger.info(f"Contact analysis: {len(atoms1)} atoms vs {len(atoms2)} atoms")
            logger.info(f"Analyzing {traj.n_frames} frames with cutoff {cutoff} Å")

            # Create atom pairs for distance calculation
            atom_pairs = []
            for atom1 in atoms1:
                for atom2 in atoms2:
                    atom_pairs.append([atom1, atom2])

            if len(atom_pairs) == 0:
                logger.warning("No atom pairs generated")
                return np.array([])

            # Calculate distances for all pairs across all frames
            distances = md.compute_distances(traj, atom_pairs)

            # Check for contacts: any pair within cutoff counts as contact for that frame
            contacts_per_frame = np.any(distances < cutoff_nm, axis=1)

            logger.info(f"Contact analysis complete: {np.sum(contacts_per_frame)} contact frames out of {traj.n_frames}")

            return contacts_per_frame

        except Exception as e:
            logger.error(f"Contact calculation failed: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return np.array([])
