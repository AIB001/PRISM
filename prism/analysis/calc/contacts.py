import logging
import numpy as np
import pickle
from pathlib import Path
from ..core.config import AnalysisConfig, convert_numpy_types
import mdtraj as md
from ..core.parallel import default_processor
from ...utils.residue import normalize_residue_to_3letter

logger = logging.getLogger(__name__)



class ContactAnalyzer:
    """Contact analysis module"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._cache_dir = Path(config.cache_dir)
        self._cache_dir.mkdir(parents=True, exist_ok=True)
    
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
        
        # Configure OpenMP threads for parallel distance computation
        with default_processor.configure_omp_threads():
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

            # Calculate distances for all pairs across all frames with parallel processing
            with default_processor.configure_omp_threads():
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

    def analyze_key_residue_contacts(self, universe, trajectory, key_residues=None, cutoff=4.0, step=1, top_n=15, cache_name=None):
        """
        Analyze contact probabilities for all protein residues with ligand, then return top N.

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str
            Path to trajectory file
        key_residues : list, optional
            If provided, analyze only these specific residues. If None, analyze all protein residues
        cutoff : float
            Contact cutoff distance in Angstroms
        step : int
            Frame step size
        top_n : int
            Number of top residues to return (by contact probability)
        cache_name : str, optional
            Custom cache name for storing results

        Returns
        -------
        dict
            Dictionary with residue names as keys and contact probabilities as values
            Format: {'ASP618': 85.2, 'ASN691': 67.8, ...} (percentages 0-100)
            Sorted by contact probability (highest first), limited to top_n residues
        """
        try:
            # Create cache key
            if cache_name is None:
                traj_name = Path(trajectory).stem if trajectory else "topology_only"
                key_res_str = "all" if key_residues is None else "_".join(key_residues[:5])  # First 5 residues
                cache_key = f"key_residue_contacts_{traj_name}_{cutoff}A_{step}step_{top_n}top_{key_res_str}"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached key residue contact results from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Load trajectory
            if trajectory:
                traj = md.load(trajectory, top=universe)
            else:
                traj = md.load(universe)

            if step > 1:
                traj = traj[::step]

            cutoff_nm = cutoff / 10.0

            # Identify ligand
            ligand_selection = "resname LIG"
            ligand_atoms = traj.topology.select(ligand_selection)

            if len(ligand_atoms) == 0:
                logger.warning("No ligand atoms found with selection 'resname LIG'")
                return {}

            # If specific residues provided, use them; otherwise analyze all protein residues
            if key_residues is not None:
                # Normalize residue names to 3-letter format for MDTraj compatibility
                residues_to_analyze = [normalize_residue_to_3letter(res) for res in key_residues]
                logger.info(f"Analyzing {len(residues_to_analyze)} specified residues vs {len(ligand_atoms)} ligand atoms")
            else:
                # Get all protein residues
                logger.info("Analyzing all protein residues vs ligand")
                protein_residues = []
                for residue in traj.topology.residues:
                    if residue.is_protein:
                        # Format: 'ASP618' (residue name + residue id)
                        res_name = f"{residue.name}{residue.resSeq}"
                        protein_residues.append(res_name)

                residues_to_analyze = protein_residues
                logger.info(f"Found {len(residues_to_analyze)} protein residues to analyze")

            residue_contact_probs = {}

            for residue in residues_to_analyze:
                try:
                    # Try different residue selection formats
                    residue_selections = [
                        f"resname {residue[:3]} and resid {residue[3:]}",  # ASP618 -> resname ASP and resid 618
                        f"residue {residue[3:]}",  # ASP618 -> residue 618
                        f"resSeq {residue[3:]} and resname {residue[:3]}"  # ASP618 -> resSeq 618 and resname ASP
                    ]

                    residue_atoms = []
                    for sel in residue_selections:
                        try:
                            atoms = traj.topology.select(sel)
                            if len(atoms) > 0:
                                residue_atoms = atoms
                                break
                        except:
                            continue

                    if len(residue_atoms) == 0:
                        logger.warning(f"Could not find residue {residue}")
                        continue

                    # Create atom pairs for distance calculation
                    atom_pairs = []
                    for lig_atom in ligand_atoms:
                        for res_atom in residue_atoms:
                            atom_pairs.append([lig_atom, res_atom])

                    if len(atom_pairs) == 0:
                        continue

                    # Calculate distances with parallel processing
                    with default_processor.configure_omp_threads():
                        distances = md.compute_distances(traj, atom_pairs)

                    # Check for contacts: any pair within cutoff counts as contact for that frame
                    contacts_per_frame = np.any(distances < cutoff_nm, axis=1)

                    # Calculate contact probability as percentage
                    contact_prob = float(np.mean(contacts_per_frame) * 100)
                    residue_contact_probs[residue] = min(100.0, max(0.0, contact_prob))

                    logger.info(f"{residue}: {contact_prob:.1f}% contact probability")

                except Exception as e:
                    logger.warning(f"Error analyzing residue {residue}: {e}")
                    continue

            logger.info(f"Contact analysis complete: {len(residue_contact_probs)} residues analyzed")

            # If no specific residues were provided (analyzing all), select top N by contact probability
            if key_residues is None and len(residue_contact_probs) > top_n:
                # Sort by contact probability (descending) and take top N
                sorted_residues = sorted(residue_contact_probs.items(), key=lambda x: x[1], reverse=True)
                top_residues = dict(sorted_residues[:top_n])

                logger.info(f"Selected top {top_n} residues with highest contact probabilities:")
                for i, (residue, prob) in enumerate(sorted_residues[:top_n], 1):
                    logger.info(f"  {i:2d}. {residue}: {prob:.1f}%")

                # Cache results
                with open(cache_file, 'wb') as f:
                    pickle.dump(top_residues, f)
                logger.info(f"Cached key residue contact results to {cache_file}")

                return top_residues
            else:
                # Cache results
                with open(cache_file, 'wb') as f:
                    pickle.dump(residue_contact_probs, f)
                logger.info(f"Cached key residue contact results to {cache_file}")

                return residue_contact_probs

        except Exception as e:
            logger.error(f"Key residue contact analysis failed: {e}")
            return {}

    def analyze_contact_numbers_timeseries(self, universe, trajectory, cutoff=4.0, step=1, cache_name=None):
        """
        Analyze contact numbers (atom pairs within cutoff) over time for LIG-protein.

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str
            Path to trajectory file
        cutoff : float
            Contact cutoff distance in Angstroms
        step : int
            Frame step size
        cache_name : str, optional
            Custom cache name for storing results

        Returns
        -------
        dict
            Dictionary containing:
            - 'contact_numbers': np.array of contact numbers per frame
            - 'times': np.array of time points (ns)
            - 'total_pairs': total number of atom pairs analyzed
        """
        try:
            # Create cache key
            if cache_name is None:
                traj_name = Path(trajectory).stem if trajectory else "topology_only"
                cache_key = f"contact_numbers_timeseries_{traj_name}_{cutoff}A_{step}step"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"Loading cached contact numbers timeseries from {cache_file}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)

            # Load trajectory
            if trajectory:
                traj = md.load(trajectory, top=universe)
            else:
                traj = md.load(universe)

            if step > 1:
                traj = traj[::step]

            cutoff_nm = cutoff / 10.0

            # Get ligand and protein atoms
            ligand_atoms = traj.topology.select("resname LIG")
            protein_atoms = traj.topology.select("protein")

            if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
                logger.warning(f"No ligand ({len(ligand_atoms)}) or protein ({len(protein_atoms)}) atoms found")
                return {}

            # Create all atom pairs
            atom_pairs = []
            for lig_atom in ligand_atoms:
                for prot_atom in protein_atoms:
                    atom_pairs.append([lig_atom, prot_atom])

            logger.info(f"Analyzing contact numbers: {len(ligand_atoms)} ligand × {len(protein_atoms)} protein = {len(atom_pairs)} pairs")

            # Calculate distances for all pairs across all frames with parallel processing
            with default_processor.configure_omp_threads():
                distances = md.compute_distances(traj, atom_pairs)

            # Count contacts per frame
            contacts_per_frame = np.sum(distances < cutoff_nm, axis=1)

            # Calculate time points (assuming 0.5 ns intervals)
            times = np.arange(len(contacts_per_frame)) * 0.5 * step

            logger.info(f"Contact numbers range: {np.min(contacts_per_frame)}-{np.max(contacts_per_frame)} contacts")

            results = {
                'contact_numbers': contacts_per_frame,
                'times': times,
                'total_pairs': len(atom_pairs),
                'n_frames': len(contacts_per_frame)
            }

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(results, f)
            logger.info(f"Cached contact numbers timeseries to {cache_file}")

            return results

        except Exception as e:
            logger.error(f"Contact numbers timeseries analysis failed: {e}")
            return {}

    def analyze_residue_contact_numbers(self, universe, trajectory, key_residues=None, cutoff=4.0, step=1):
        """
        Analyze contact numbers between ligand and specific residues with normalization.

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str
            Path to trajectory file
        key_residues : list, optional
            List of residue names to analyze
        cutoff : float
            Contact cutoff distance in Angstroms
        step : int
            Frame step size

        Returns
        -------
        dict
            Dictionary with residue names as keys and contact data as values
            Format: {'ASP618': {'contact_numbers': array, 'normalized': array, 'n_atoms': int}, ...}
        """
        try:
            # Load trajectory
            if trajectory:
                traj = md.load(trajectory, top=universe)
            else:
                traj = md.load(universe)

            if step > 1:
                traj = traj[::step]

            cutoff_nm = cutoff / 10.0

            # Default key residues
            if key_residues is None:
                key_residues = ['ASP618', 'ASP623', 'ASP760', 'ASN691', 'SER759',
                               'THR680', 'LYS551', 'ARG553', 'ARG555']
            else:
                # Normalize residue names to 3-letter format for MDTraj compatibility
                key_residues = [normalize_residue_to_3letter(res) for res in key_residues]

            # Get ligand atoms
            ligand_atoms = traj.topology.select("resname LIG")
            if len(ligand_atoms) == 0:
                logger.warning("No ligand atoms found")
                return {}

            logger.info(f"Analyzing residue contact numbers: {len(ligand_atoms)} ligand atoms")

            residue_contact_data = {}

            for residue in key_residues:
                try:
                    # Try different residue selection formats
                    residue_selections = [
                        f"resname {residue[:3]} and resid {residue[3:]}",
                        f"residue {residue[3:]}",
                        f"resSeq {residue[3:]} and resname {residue[:3]}"
                    ]

                    residue_atoms = []
                    for sel in residue_selections:
                        try:
                            atoms = traj.topology.select(sel)
                            if len(atoms) > 0:
                                residue_atoms = atoms
                                break
                        except:
                            continue

                    if len(residue_atoms) == 0:
                        logger.warning(f"Could not find residue {residue}")
                        continue

                    # Create atom pairs
                    atom_pairs = []
                    for lig_atom in ligand_atoms:
                        for res_atom in residue_atoms:
                            atom_pairs.append([lig_atom, res_atom])

                    if len(atom_pairs) == 0:
                        continue

                    # Calculate distances and contact numbers with parallel processing
                    with default_processor.configure_omp_threads():
                        distances = md.compute_distances(traj, atom_pairs)
                    contact_numbers = np.sum(distances < cutoff_nm, axis=1)

                    # Normalize by residue atom count
                    normalized_contacts = contact_numbers / len(residue_atoms)

                    residue_contact_data[residue] = {
                        'contact_numbers': contact_numbers,
                        'normalized': normalized_contacts,
                        'n_atoms': len(residue_atoms),
                        'n_pairs': len(atom_pairs)
                    }

                    logger.info(f"{residue}: {len(residue_atoms)} atoms, avg contacts={np.mean(contact_numbers):.1f}, normalized={np.mean(normalized_contacts):.2f}")

                except Exception as e:
                    logger.warning(f"Error analyzing residue {residue}: {e}")
                    continue

            logger.info(f"Residue contact numbers analysis complete: {len(residue_contact_data)} residues")
            return residue_contact_data

        except Exception as e:
            logger.error(f"Residue contact numbers analysis failed: {e}")
            return {}

    def analyze_contact_distances(self, universe, trajectory, key_residues=None, cutoff=4.0, step=1):
        """
        Analyze distribution of contact distances between ligand and residues.

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str
            Path to trajectory file
        key_residues : list, optional
            List of residue names to analyze
        cutoff : float
            Contact cutoff distance in Angstroms (only pairs within this distance are included)
        step : int
            Frame step size

        Returns
        -------
        dict
            Dictionary with residue names as keys and distance arrays as values
            Format: {'ASP618': distances_array, 'ASN691': distances_array, ...}
        """
        try:
            # Load trajectory
            if trajectory:
                traj = md.load(trajectory, top=universe)
            else:
                traj = md.load(universe)

            if step > 1:
                traj = traj[::step]

            cutoff_nm = cutoff / 10.0

            # Default key residues
            if key_residues is None:
                key_residues = ['ASP618', 'ASP623', 'ASP760', 'ASN691', 'SER759',
                               'THR680', 'LYS551', 'ARG553', 'ARG555']
            else:
                # Normalize residue names to 3-letter format for MDTraj compatibility
                key_residues = [normalize_residue_to_3letter(res) for res in key_residues]

            # Get ligand atoms
            ligand_atoms = traj.topology.select("resname LIG")
            if len(ligand_atoms) == 0:
                logger.warning("No ligand atoms found")
                return {}

            logger.info(f"Analyzing contact distances for {len(key_residues)} residues")

            residue_distances = {}

            for residue in key_residues:
                try:
                    # Find residue atoms
                    residue_selections = [
                        f"resname {residue[:3]} and resid {residue[3:]}",
                        f"residue {residue[3:]}",
                        f"resSeq {residue[3:]} and resname {residue[:3]}"
                    ]

                    residue_atoms = []
                    for sel in residue_selections:
                        try:
                            atoms = traj.topology.select(sel)
                            if len(atoms) > 0:
                                residue_atoms = atoms
                                break
                        except:
                            continue

                    if len(residue_atoms) == 0:
                        continue

                    # Create atom pairs and calculate distances
                    atom_pairs = []
                    for lig_atom in ligand_atoms:
                        for res_atom in residue_atoms:
                            atom_pairs.append([lig_atom, res_atom])

                    # Calculate distances with parallel processing
                    with default_processor.configure_omp_threads():
                        distances = md.compute_distances(traj, atom_pairs)

                    # Get minimum distance per frame (closest ligand-residue approach)
                    # This gives one distance value per frame - no cutoff filtering
                    min_distances_per_frame = np.min(distances, axis=1) * 10.0  # Convert to Angstroms

                    if len(min_distances_per_frame) > 0:
                        residue_distances[residue] = min_distances_per_frame
                        logger.info(f"{residue}: {len(min_distances_per_frame)} frames, distance range {np.min(min_distances_per_frame):.1f}-{np.max(min_distances_per_frame):.1f} Å")

                except Exception as e:
                    logger.warning(f"Error analyzing distances for residue {residue}: {e}")
                    continue

            logger.info(f"Contact distances analysis complete: {len(residue_distances)} residues")
            return residue_distances

        except Exception as e:
            logger.error(f"Contact distances analysis failed: {e}")
            return {}

    def analyze_residue_distance_timeseries(self, universe, trajectory, key_residues=None, step=1, cache_name=None):
        """
        Analyze minimum distance timeseries between ligand and specific residues.

        Parameters
        ----------
        universe : str
            Path to topology file
        trajectory : str
            Path to trajectory file
        key_residues : list, optional
            List of residue names to analyze (e.g., ['ASP618', 'ARG555'])
        step : int
            Frame step size
        cache_name : str, optional
            Custom cache name for storing results

        Returns
        -------
        dict
            Dictionary with residue names as keys and distance timeseries as values
            Format: {'ASP618': distances_array[n_frames], 'ARG555': distances_array[n_frames], ...}
            Also includes 'times' key with time array in nanoseconds
        """
        try:
            # Create cache key
            if cache_name is None:
                traj_name = Path(trajectory).stem if trajectory else "topology_only"
                res_str = "_".join(key_residues[:3]) if key_residues else "default"
                cache_key = f"residue_distance_timeseries_{traj_name}_{res_str}_{step}step"
                cache_key = cache_key.replace(" ", "_").replace("(", "").replace(")", "")
            else:
                cache_key = cache_name

            cache_file = self._cache_dir / f"{cache_key}.pkl"

            # Check cache
            if cache_file.exists():
                logger.info(f"✓ Loading cached residue distance timeseries from {cache_file.name}")
                print(f"  ✓ Using cached distance timeseries: {cache_file.name}")
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)
            else:
                logger.info(f"✗ Cache miss, will compute and save to: {cache_file.name}")
                print(f"  ⚙ Computing distance timeseries (will cache to: {cache_file.name})")

            # Load trajectory
            if trajectory:
                traj = md.load(trajectory, top=universe)
            else:
                traj = md.load(universe)

            if step > 1:
                traj = traj[::step]

            # Default key residues
            if key_residues is None:
                key_residues = ['ASP618', 'ASP623', 'ASP760', 'ASN691', 'SER759',
                               'THR680', 'LYS551', 'ARG553', 'ARG555']
            else:
                # Normalize residue names to 3-letter format for MDTraj compatibility
                key_residues = [normalize_residue_to_3letter(res) for res in key_residues]

            # Get ligand atoms
            ligand_atoms = traj.topology.select("resname LIG")
            if len(ligand_atoms) == 0:
                logger.warning("No ligand atoms found")
                return {}

            logger.info(f"Analyzing distance timeseries for {len(key_residues)} residues over {traj.n_frames} frames")

            residue_timeseries = {}

            for residue in key_residues:
                try:
                    # Find residue atoms
                    residue_selections = [
                        f"resname {residue[:3]} and resid {residue[3:]}",
                        f"residue {residue[3:]}",
                        f"resSeq {residue[3:]} and resname {residue[:3]}"
                    ]

                    residue_atoms = []
                    for sel in residue_selections:
                        try:
                            atoms = traj.topology.select(sel)
                            if len(atoms) > 0:
                                residue_atoms = atoms
                                break
                        except:
                            continue

                    if len(residue_atoms) == 0:
                        logger.warning(f"Could not find residue {residue}")
                        continue

                    # Create atom pairs
                    atom_pairs = []
                    for lig_atom in ligand_atoms:
                        for res_atom in residue_atoms:
                            atom_pairs.append([lig_atom, res_atom])

                    # Calculate distances with parallel processing
                    with default_processor.configure_omp_threads():
                        distances = md.compute_distances(traj, atom_pairs)

                    # Get minimum distance per frame and convert to Angstroms
                    min_distances_per_frame = np.min(distances, axis=1) * 10.0

                    residue_timeseries[residue] = min_distances_per_frame
                    logger.info(f"{residue}: {len(min_distances_per_frame)} frames, range {np.min(min_distances_per_frame):.1f}-{np.max(min_distances_per_frame):.1f} Å")

                except Exception as e:
                    logger.warning(f"Error analyzing timeseries for residue {residue}: {e}")
                    continue

            # Add time array
            times = np.arange(traj.n_frames) * 0.5 * step
            residue_timeseries['times'] = times

            # Cache results
            with open(cache_file, 'wb') as f:
                pickle.dump(residue_timeseries, f)
            logger.info(f"Cached residue distance timeseries to {cache_file.name}")

            logger.info(f"Distance timeseries analysis complete: {len(residue_timeseries)-1} residues")
            return residue_timeseries

        except Exception as e:
            logger.error(f"Distance timeseries analysis failed: {e}")
            return {}
