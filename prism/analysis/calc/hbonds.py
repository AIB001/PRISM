import logging
import numpy as np
from typing import List, Tuple, Dict
from ..core.config import AnalysisConfig

try:
    import mdtraj as md

    MDTRAJ_AVAILABLE = True
except ImportError:
    MDTRAJ_AVAILABLE = False
    md = None

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

        self._topology_cache = {}
        self._hbond_pairs_cache = None
        self._hbond_triplets_cache = None

        self._residue_h_atoms = {}
        self._atom_residue_map = {}

        self._prescreen_cutoff = 8.0

    def set_trajectory_data(self, traj, ligand_residue, protein_residues=None):
        self.traj = traj
        self.ligand_residue = ligand_residue

        if protein_residues is not None:
            self.protein_residues = protein_residues
        else:
            self.protein_residues = []
            standard_aa = {
                "ALA",
                "ARG",
                "ASN",
                "ASP",
                "CYS",
                "GLN",
                "GLU",
                "GLY",
                "HIS",
                "ILE",
                "LEU",
                "LYS",
                "MET",
                "PHE",
                "PRO",
                "SER",
                "THR",
                "TRP",
                "TYR",
                "VAL",
            }
            for residue in traj.topology.residues:
                if residue.name in standard_aa:
                    self.protein_residues.append(residue)

        try:
            self._cache_topology_info()
            self._precompute_h_atom_mapping()
            self._prepare_hbond_analysis()
            logger.warning(f"Hydrogen bond analysis prepared with {len(self.protein_residues)} protein residues")
        except Exception as e:
            logger.warning(f"Failed to prepare hydrogen bond analysis: {e}")

    def _cache_topology_info(self):
        """Cache topology information to avoid repeated queries"""
        for atom in self.traj.topology.atoms:
            self._topology_cache[atom.index] = {
                "name": atom.name,
                "element": atom.element.symbol,
                "residue_name": atom.residue.name,
                "residue_seq": atom.residue.resSeq,
                "residue_index": atom.residue.index,
            }
            self._atom_residue_map[atom.index] = atom.residue.index

    def _precompute_h_atom_mapping(self):
        """Precompute hydrogen atom mapping for each residue"""
        for residue in self.traj.topology.residues:
            h_atoms = []
            for atom in residue.atoms:
                if atom.element.symbol == "H":
                    h_atoms.append(atom.index)
            self._residue_h_atoms[residue.index] = h_atoms

    def _prepare_hbond_analysis(self):
        try:
            if not hasattr(self, "traj") or self.traj is None:
                raise ValueError("No trajectory data available")

            frame0 = self.traj.xyz[0]

            if not self.protein_residues:
                logger.warning("No protein residues found for hydrogen bond analysis")
                return

            self._find_donors_acceptors_vectorized(frame0, self.protein_residues, is_protein=True)

            if self.ligand_residue is not None:
                # Validate that ligand_residue is a proper residue object, not a string
                if hasattr(self.ligand_residue, "index") and hasattr(self.ligand_residue, "atoms"):
                    self._find_donors_acceptors_vectorized(frame0, [self.ligand_residue], is_protein=False)
                else:
                    logger.warning(f"Invalid ligand residue object: {type(self.ligand_residue)}")
            else:
                logger.warning("No ligand residue found for hydrogen bond analysis")

            logger.warning(
                f"Found - Protein donors: {len(self.protein_donors)}, acceptors: {len(self.protein_acceptors)}"
            )
            logger.warning(f"Found - Ligand donors: {len(self.ligand_donors)}, acceptors: {len(self.ligand_acceptors)}")

            total_donors = len(self.protein_donors) + len(self.ligand_donors)
            total_acceptors = len(self.protein_acceptors) + len(self.ligand_acceptors)

            if total_donors == 0 or total_acceptors == 0:
                logger.warning(
                    f"Insufficient donors ({total_donors}) or acceptors ({total_acceptors}) for hydrogen bond analysis"
                )

        except Exception as e:
            logger.error(f"Failed to prepare hydrogen bond analysis: {e}")
            logger.warning("Hydrogen bond analysis will be skipped")

    def _find_donors_acceptors_vectorized(self, frame0: np.ndarray, residues: List, is_protein: bool):
        """Vectorized donor/acceptor finding"""
        donor_list = self.protein_donors if is_protein else self.ligand_donors
        acceptor_list = self.protein_acceptors if is_protein else self.ligand_acceptors

        for residue in residues:
            # Validate residue object
            if not hasattr(residue, "index") or not hasattr(residue, "atoms"):
                logger.warning(f"Invalid residue object: {type(residue)} - {residue}")
                continue

            residue_idx = residue.index
            h_indices = self._residue_h_atoms.get(residue_idx, [])

            no_indices = []
            for atom in residue.atoms:
                try:
                    if atom.element.symbol in ["N", "O"]:
                        no_indices.append(atom.index)
                except AttributeError:
                    continue

            if not no_indices:
                continue

            acceptor_list.extend(no_indices)

            if not h_indices:
                continue

            no_coords = frame0[no_indices]
            h_coords = frame0[h_indices]

            diff_vectors = no_coords[:, np.newaxis, :] - h_coords[np.newaxis, :, :]
            distances = np.linalg.norm(diff_vectors, axis=2)

            bonded_mask = distances < self.config.bond_length_threshold_nm
            donor_mask = np.any(bonded_mask, axis=1)

            donor_indices = np.array(no_indices)[donor_mask]
            donor_list.extend(donor_indices.tolist())

    def analyze_hydrogen_bonds(self, universe) -> List[Tuple[str, float, float, float]]:
        if universe is not None and hasattr(self, "traj"):
            logger.debug("Hydrogen bond analysis uses cached MDTraj trajectory; ignoring universe input.")

        if not hasattr(self, "traj") or self.traj is None:
            logger.warning("No MDTraj trajectory available for hydrogen bond analysis")
            return []

        hbond_pairs, hbond_triplets = self._prepare_hbond_pairs()

        if not hbond_pairs:
            logger.warning("No hydrogen bond pairs found")
            return []

        total_frames = self.traj.n_frames
        hbond_stats = self._analyze_hbonds_fully_vectorized(hbond_pairs, hbond_triplets)

        return self._calculate_hbond_frequencies(hbond_stats, total_frames)

    def extract_residue_hbond_timeseries(self, target_residues: List[str] = None) -> Dict[str, np.ndarray]:
        """
        Extract per-frame hydrogen bond time series for specific residues.

        Parameters
        ----------
        target_residues : list of str, optional
            List of residue names to extract (e.g., ['ASP623', 'ASN691']).
            If None, extracts all residues with H-bonds.

        Returns
        -------
        dict
            Dictionary mapping residue names to boolean arrays indicating H-bond presence per frame.
            Format: {residue_name: np.ndarray of shape (n_frames,)}
        """
        if not hasattr(self, "traj") or self.traj is None:
            logger.warning("No trajectory data available")
            return {}

        # Get H-bond pairs
        hbond_pairs, hbond_triplets = self._prepare_hbond_pairs()
        if not hbond_pairs:
            logger.warning("No hydrogen bond pairs found")
            return {}

        # Create keys for all H-bond pairs
        hbond_keys = [self._create_hbond_key_cached(pair) for pair in hbond_pairs]

        # Calculate H-bonds for all frames
        pairs_array = np.array(hbond_pairs)
        triplets_array = np.array(hbond_triplets)

        n_frames = self.traj.n_frames
        n_pairs = len(pairs_array)

        # Compute H-bond mask for all frames
        logger.info(f"Computing H-bond time series for {n_frames} frames...")

        try:
            # Load coordinates
            coords = self.traj.xyz

            # Calculate distances
            donor_coords = coords[:, pairs_array[:, 0], :]
            acceptor_coords = coords[:, pairs_array[:, 1], :]
            all_distances = np.linalg.norm(donor_coords - acceptor_coords, axis=2)

            # Calculate angles
            h_coords = coords[:, triplets_array[:, 0], :]
            donor_coords_triplet = coords[:, triplets_array[:, 1], :]
            acceptor_coords_triplet = coords[:, triplets_array[:, 2], :]

            vec_hd = donor_coords_triplet - h_coords
            vec_da = acceptor_coords_triplet - donor_coords_triplet

            dot_products = np.sum(vec_hd * vec_da, axis=2)
            norm_hd = np.linalg.norm(vec_hd, axis=2)
            norm_da = np.linalg.norm(vec_da, axis=2)

            denominator = norm_hd * norm_da
            valid_mask = denominator > 1e-10
            cos_angles = np.zeros((n_frames, n_pairs))
            cos_angles[valid_mask] = dot_products[valid_mask] / denominator[valid_mask]
            cos_angles = np.clip(cos_angles, -1.0, 1.0)
            all_angles_deg = np.degrees(np.arccos(cos_angles))

            # Apply H-bond criteria
            distance_mask = all_distances < self.config.hbond_distance_cutoff_nm
            angle_mask = all_angles_deg > self.config.hbond_angle_cutoff_deg
            hbond_mask = distance_mask & angle_mask  # Shape: (n_frames, n_pairs)

            # Group by residue
            residue_timeseries = {}

            for pair_idx, key in enumerate(hbond_keys):
                # Extract residue name from key
                # Key format: "ASP 623 (OD1) -> Ligand" or "Ligand -> ASP 623 (OD1)"
                residue_name = self._extract_residue_name_from_key(key)

                if residue_name is None:
                    continue

                # Filter by target residues if specified
                if target_residues is not None:
                    if not any(res in residue_name for res in target_residues):
                        continue

                # Store time series (OR operation if multiple bonds to same residue)
                if residue_name not in residue_timeseries:
                    residue_timeseries[residue_name] = hbond_mask[:, pair_idx].copy()
                else:
                    # Combine with existing (H-bond present if ANY bond to this residue exists)
                    residue_timeseries[residue_name] = np.logical_or(
                        residue_timeseries[residue_name], hbond_mask[:, pair_idx]
                    )

            logger.info(f"Extracted H-bond time series for {len(residue_timeseries)} residues")
            return residue_timeseries

        except Exception as e:
            logger.error(f"Error extracting H-bond time series: {e}")
            import traceback

            traceback.print_exc()
            return {}

    def _extract_residue_name_from_key(self, key: str) -> str:
        """
        Extract residue name from H-bond key.

        Parameters
        ----------
        key : str
            H-bond key like "ASP 623 (OD1) -> Ligand" or "Ligand -> ASN 691 (ND2)"

        Returns
        -------
        str or None
            Residue name like "ASP623" or None if not found
        """
        try:
            # Split by "->"
            parts = key.split(" -> ")

            # Find protein side (not "Ligand")
            protein_part = None
            for part in parts:
                if "Ligand" not in part:
                    protein_part = part.strip()
                    break

            if protein_part is None:
                return None

            # Extract residue name and number
            # Format: "ASP 623 (OD1)"
            tokens = protein_part.split()
            if len(tokens) >= 2:
                res_name = tokens[0]  # "ASP"
                res_num = tokens[1]  # "623"
                return f"{res_name}{res_num}"

            return None

        except Exception as e:
            logger.warning(f"Failed to extract residue name from key '{key}': {e}")
            return None

    def _prepare_hbond_pairs(self) -> Tuple[List, List]:
        """Prepare hydrogen bond pairs with caching"""
        if self._hbond_pairs_cache is not None:
            return self._hbond_pairs_cache, self._hbond_triplets_cache

        hbond_pairs = []
        hbond_triplets = []

        self._add_hbond_pairs_vectorized(
            self.protein_donors, self.ligand_acceptors, hbond_pairs, hbond_triplets, donor_is_protein=True
        )

        self._add_hbond_pairs_vectorized(
            self.ligand_donors, self.protein_acceptors, hbond_pairs, hbond_triplets, donor_is_protein=False
        )

        self._hbond_pairs_cache = hbond_pairs
        self._hbond_triplets_cache = hbond_triplets

        return hbond_pairs, hbond_triplets

    def _add_hbond_pairs_vectorized(
        self, donors: List[int], acceptors: List[int], hbond_pairs: List, hbond_triplets: List, donor_is_protein: bool
    ):
        """Add hydrogen bond pairs with distance prescreening"""
        donor_label = "protein" if donor_is_protein else "ligand"
        logger.debug("Adding hydrogen bond pairs for %s donors", donor_label)

        if not donors or not acceptors:
            return

        frame0 = self.traj.xyz[0]

        donor_coords = frame0[donors]
        acceptor_coords = frame0[acceptors]

        donor_acceptor_dists = np.linalg.norm(
            donor_coords[:, np.newaxis, :] - acceptor_coords[np.newaxis, :, :], axis=2
        )

        close_pairs = np.where(donor_acceptor_dists < self._prescreen_cutoff)
        close_donor_indices = close_pairs[0]
        close_acceptor_indices = close_pairs[1]

        donor_by_residue = {}
        for i, donor_idx in enumerate(donors):
            residue_idx = self._atom_residue_map[donor_idx]
            if residue_idx not in donor_by_residue:
                donor_by_residue[residue_idx] = []
            donor_by_residue[residue_idx].append((i, donor_idx))

        for residue_idx, residue_donor_info in donor_by_residue.items():
            h_atoms = self._residue_h_atoms.get(residue_idx, [])

            if not h_atoms:
                continue

            local_donor_indices = [info[0] for info in residue_donor_info]
            residue_donors = [info[1] for info in residue_donor_info]

            donor_coords_residue = frame0[residue_donors]
            h_coords = frame0[h_atoms]

            distances = np.linalg.norm(donor_coords_residue[:, np.newaxis, :] - h_coords[np.newaxis, :, :], axis=2)

            for local_idx, (donor_array_idx, donor_atom_idx) in enumerate(zip(local_donor_indices, residue_donors)):
                bonded_h_mask = distances[local_idx] < self.config.bond_length_threshold_nm
                bonded_h_indices = np.where(bonded_h_mask)[0]

                if len(bonded_h_indices) > 0:
                    closest_h_idx = bonded_h_indices[np.argmin(distances[local_idx][bonded_h_indices])]
                    hydrogen_atom_idx = h_atoms[closest_h_idx]

                    valid_acceptor_mask = close_donor_indices == donor_array_idx
                    valid_acceptor_local_indices = close_acceptor_indices[valid_acceptor_mask]

                    for acceptor_local_idx in valid_acceptor_local_indices:
                        acceptor_atom_idx = acceptors[acceptor_local_idx]
                        hbond_pairs.append([donor_atom_idx, acceptor_atom_idx])
                        hbond_triplets.append([hydrogen_atom_idx, donor_atom_idx, acceptor_atom_idx])

    def _analyze_hbonds_fully_vectorized(self, hbond_pairs: List, hbond_triplets: List) -> Dict:
        """Fully vectorized hydrogen bond analysis"""
        logger.warning(f"Analyzing hydrogen bonds for {self.traj.n_frames} frames...")

        if not hbond_pairs:
            return {}

        pairs_array = np.array(hbond_pairs)
        triplets_array = np.array(hbond_triplets)

        hbond_keys = [self._create_hbond_key_cached(pair) for pair in hbond_pairs]

        total_frames = self.traj.n_frames
        max_frames_for_analysis = 1000
        frame_step = max(1, total_frames // max_frames_for_analysis)
        frame_indices = np.arange(0, total_frames, frame_step)

        logger.warning(f"Analyzing {len(frame_indices)} frames (step={frame_step}) out of {total_frames} total frames")

        try:
            coords_subset = self.traj.xyz[frame_indices]
            n_frames, n_atoms, _ = coords_subset.shape
            n_pairs = len(pairs_array)

            donor_coords = coords_subset[:, pairs_array[:, 0], :]
            acceptor_coords = coords_subset[:, pairs_array[:, 1], :]

            all_distances = np.linalg.norm(donor_coords - acceptor_coords, axis=2)

            h_coords = coords_subset[:, triplets_array[:, 0], :]
            donor_coords_triplet = coords_subset[:, triplets_array[:, 1], :]
            acceptor_coords_triplet = coords_subset[:, triplets_array[:, 2], :]

            vec_hd = donor_coords_triplet - h_coords
            vec_da = acceptor_coords_triplet - donor_coords_triplet

            dot_products = np.sum(vec_hd * vec_da, axis=2)

            norm_hd = np.linalg.norm(vec_hd, axis=2)
            norm_da = np.linalg.norm(vec_da, axis=2)

            denominator = norm_hd * norm_da
            valid_mask = denominator > 1e-10
            cos_angles = np.zeros((n_frames, n_pairs))
            cos_angles[valid_mask] = dot_products[valid_mask] / denominator[valid_mask]

            cos_angles = np.clip(cos_angles, -1.0, 1.0)
            all_angles_deg = np.degrees(np.arccos(cos_angles))

            distance_mask = all_distances < self.config.hbond_distance_cutoff_nm
            angle_mask = all_angles_deg > self.config.hbond_angle_cutoff_deg
            hbond_mask = distance_mask & angle_mask

            hbond_stats = {}
            for pair_idx in range(n_pairs):
                valid_frames = np.where(hbond_mask[:, pair_idx])[0]

                if len(valid_frames) > 0:
                    key = hbond_keys[pair_idx]

                    if key not in hbond_stats:
                        hbond_stats[key] = {"count": 0, "distance": [], "angle": []}

                    hbond_stats[key]["count"] += len(valid_frames)
                    valid_distances = (all_distances[valid_frames, pair_idx] * 10.0).tolist()
                    valid_angles = all_angles_deg[valid_frames, pair_idx].tolist()

                    hbond_stats[key]["distance"].extend(valid_distances)
                    hbond_stats[key]["angle"].extend(valid_angles)

            logger.warning(f"Successfully analyzed {len(hbond_stats)} unique hydrogen bonds")
            return hbond_stats

        except Exception as e:
            logger.error(f"Error in fully vectorized analysis: {e}")
            return self._analyze_hbonds_custom_batch_fallback(hbond_pairs, hbond_triplets, hbond_keys)

    def _analyze_hbonds_custom_batch_fallback(
        self, hbond_pairs: List, hbond_triplets: List, hbond_keys: List[str]
    ) -> Dict:
        """Fallback batch processing method"""
        logger.warning("Using fallback batch processing method")

        batch_size = min(50, max(10, self.traj.n_frames // 20))
        frame_step = max(1, self.traj.n_frames // 500)
        frame_indices = list(range(0, self.traj.n_frames, frame_step))

        hbond_stats = {}
        pairs_array = np.array(hbond_pairs)
        triplets_array = np.array(hbond_triplets)

        for i in range(0, len(frame_indices), batch_size):
            batch_frames = frame_indices[i : i + batch_size]
            try:
                batch_stats = self._compute_hbonds_custom_batch(batch_frames, pairs_array, triplets_array, hbond_keys)

                for key, data in batch_stats.items():
                    if key not in hbond_stats:
                        hbond_stats[key] = {"count": 0, "distance": [], "angle": []}
                    hbond_stats[key]["count"] += data["count"]
                    hbond_stats[key]["distance"].extend(data["distance"])
                    hbond_stats[key]["angle"].extend(data["angle"])

            except Exception as e:
                logger.warning(f"Error processing batch {i//batch_size}: {e}")
                continue

        return hbond_stats

    def _compute_hbonds_custom_batch(
        self, frame_indices: List[int], pairs_array: np.ndarray, triplets_array: np.ndarray, hbond_keys: List[str]
    ) -> Dict:
        """Custom batch computation"""
        if len(pairs_array) == 0 or len(frame_indices) == 0:
            return {}

        batch_stats = {}

        try:
            batch_coords = self.traj.xyz[frame_indices]

            n_frames = len(frame_indices)
            n_pairs = len(pairs_array)

            all_distances = np.zeros((n_frames, n_pairs))
            all_angles = np.zeros((n_frames, n_pairs))

            for frame_idx in range(n_frames):
                coords = batch_coords[frame_idx]

                donor_coords = coords[pairs_array[:, 0]]
                acceptor_coords = coords[pairs_array[:, 1]]
                distances = np.linalg.norm(donor_coords - acceptor_coords, axis=1)
                all_distances[frame_idx] = distances

                h_coords = coords[triplets_array[:, 0]]
                donor_coords = coords[triplets_array[:, 1]]
                acceptor_coords = coords[triplets_array[:, 2]]

                vec_hd = donor_coords - h_coords
                vec_da = acceptor_coords - donor_coords

                dot_products = np.sum(vec_hd * vec_da, axis=1)
                norm_hd = np.linalg.norm(vec_hd, axis=1)
                norm_da = np.linalg.norm(vec_da, axis=1)

                denominator = norm_hd * norm_da
                valid_mask = denominator > 1e-10
                cos_angles = np.zeros(n_pairs)
                cos_angles[valid_mask] = dot_products[valid_mask] / denominator[valid_mask]

                cos_angles = np.clip(cos_angles, -1.0, 1.0)
                angles_rad = np.arccos(cos_angles)
                angles_deg = np.degrees(angles_rad)
                all_angles[frame_idx] = angles_deg

            distance_mask = all_distances < self.config.hbond_distance_cutoff_nm
            angle_mask = all_angles > self.config.hbond_angle_cutoff_deg
            hbond_mask = distance_mask & angle_mask

            for pair_idx in range(n_pairs):
                valid_frames = np.where(hbond_mask[:, pair_idx])[0]

                if len(valid_frames) > 0:
                    key = hbond_keys[pair_idx]

                    if key not in batch_stats:
                        batch_stats[key] = {"count": 0, "distance": [], "angle": []}

                    batch_stats[key]["count"] += len(valid_frames)
                    batch_stats[key]["distance"].extend((all_distances[valid_frames, pair_idx] * 10.0).tolist())
                    batch_stats[key]["angle"].extend(all_angles[valid_frames, pair_idx].tolist())

        except Exception as e:
            logger.warning(f"Error in custom batch computation: {e}")
            raise

        return batch_stats

    def _create_hbond_key_cached(self, pair: List[int]) -> str:
        """Create hydrogen bond key using cached topology info"""
        try:
            donor_idx, acceptor_idx = pair
            donor_info = self._topology_cache[donor_idx]
            acceptor_info = self._topology_cache[acceptor_idx]

            if donor_idx in self.protein_donors:
                return f"{donor_info['residue_name']} {donor_info['residue_seq']} " f"({donor_info['name']}) -> Ligand"
            else:
                return (
                    f"Ligand -> {acceptor_info['residue_name']} "
                    f"{acceptor_info['residue_seq']} ({acceptor_info['name']})"
                )
        except:
            return f"Unknown_HBond_{donor_idx}_{acceptor_idx}"

    def _calculate_hbond_frequencies(
        self, hbond_stats: Dict, total_frames: int
    ) -> List[Tuple[str, float, float, float]]:
        hbond_frequencies = []
        for key, data in hbond_stats.items():
            freq = data["count"] / total_frames
            avg_dist = np.mean(data["distance"]) if data["distance"] else 0
            avg_angle = np.mean(data["angle"]) if data["angle"] else 0
            hbond_frequencies.append((key, freq, avg_dist, avg_angle))

        hbond_frequencies.sort(key=lambda x: x[1], reverse=True)
        self.hbond_stats = hbond_stats
        return hbond_frequencies
