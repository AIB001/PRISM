import os
import json
import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

try:
    import mdtraj as md
    MDTRAJ_AVAILABLE = True
except ImportError:
    MDTRAJ_AVAILABLE = False
    md = None
logger = logging.getLogger(__name__)

class DataExporter:
    """Data export module"""
    
    @staticmethod
    def save_contact_data_csv(contact_proportions: Dict[str, float], output_path: str):
        try:
            df = pd.DataFrame([
                {
                    'residue_name': key.split()[0] if ' ' in key else key[:3],
                    'residue_number': int(''.join(filter(str.isdigit, key))) if any(c.isdigit() for c in key) else 0,
                    'residue_id': key,
                    'contact_proportion': value,
                    'contact_percentage': value * 100
                }
                for key, value in contact_proportions.items()
            ])
            df = df.sort_values('contact_proportion', ascending=False)
            df.to_csv(output_path, index=False, float_format='%.4f')
            logger.warning(f"Contact data saved: {output_path}")
            return df
        except Exception as e:
            logger.error(f"Failed to save contact data: {e}")
            return None
    
    @staticmethod
    def save_hbond_data_csv(hbond_frequencies: List[Tuple], hbond_stats: Dict, output_path: str):
        try:
            if not hbond_frequencies:
                logger.warning("No hydrogen bond data to save")
                return None
            
            data = []
            for hb in hbond_frequencies:
                key, freq, avg_dist, avg_angle = hb
                
                if ' -> ' in key:
                    donor_part, acceptor_part = key.split(' -> ')
                else:
                    donor_part, acceptor_part = key, "Unknown"
                
                detailed_stats = {}
                if key in hbond_stats:
                    stats = hbond_stats[key]
                    detailed_stats = {
                        'total_occurrences': stats['count'],
                        'distance_std': np.std(stats['distance']) if stats['distance'] else 0,
                        'angle_std': np.std(stats['angle']) if stats['angle'] else 0,
                        'distance_min': np.min(stats['distance']) if stats['distance'] else 0,
                        'distance_max': np.max(stats['distance']) if stats['distance'] else 0,
                        'angle_min': np.min(stats['angle']) if stats['angle'] else 0,
                        'angle_max': np.max(stats['angle']) if stats['angle'] else 0
                    }
                
                data.append({
                    'hbond_id': key,
                    'donor': donor_part,
                    'acceptor': acceptor_part,
                    'frequency': freq,
                    'frequency_percentage': freq * 100,
                    'avg_distance_angstrom': avg_dist,
                    'avg_angle_degrees': avg_angle,
                    **detailed_stats
                })
            
            df = pd.DataFrame(data)
            df = df.sort_values('frequency', ascending=False)
            df.to_csv(output_path, index=False, float_format='%.4f')
            logger.warning(f"Hydrogen bond data saved: {output_path}")
            return df
        except Exception as e:
            logger.error(f"Failed to save hydrogen bond data: {e}")
            return None
    
    @staticmethod
    def save_distance_data_csv(distance_stats: Dict[str, Dict[str, float]], output_path: str):
        try:
            data = []
            for residue_id, stats in distance_stats.items():
                data.append({
                    'residue_id': residue_id,
                    'mean_distance_angstrom': stats['mean'],
                    'min_distance_angstrom': stats['min'],
                    'max_distance_angstrom': stats['max'],
                    'std_distance_angstrom': stats['std']
                })
            
            df = pd.DataFrame(data)
            df.to_csv(output_path, index=False, float_format='%.4f')
            logger.warning(f"Distance data saved: {output_path}")
            return df
        except Exception as e:
            logger.error(f"Failed to save distance data: {e}")
            return None

    @staticmethod
    def save_multi_system_data(all_data: List[np.ndarray], all_times: List[np.ndarray], 
                              system_names: List[str], weights: List[float], 
                              analysis_type: str, parameters: Dict, output_dir: str):
        """Save multi-system analysis data"""
        try:
            os.makedirs(output_dir, exist_ok=True)
            
            multi_data = {
                'system_names': system_names,
                'molecular_weights': weights,
                'analysis_type': analysis_type,
                'parameters': parameters,
                'data': {}
            }
            
            for i, (name, data, times, weight) in enumerate(zip(system_names, all_data, all_times, weights)):
                multi_data['data'][name] = {
                    'contacts': data.tolist(),
                    'times': times.tolist(),
                    'molecular_weight': weight
                }
            
            data_path = os.path.join(output_dir, f"{analysis_type}_data.json")
            with open(data_path, 'w') as f:
                json.dump(multi_data, f, indent=2)
            logger.warning(f"Multi-system data saved: {data_path}")
            
            summary_data = []
            for name, data, times, weight in zip(system_names, all_data, all_times, weights):
                if len(data) > 0:
                    summary_data.append({
                        'system_name': name,
                        'molecular_weight': weight,
                        'final_contacts': data[-1] if len(data) > 0 else 0,
                        'max_contacts': np.max(data) if len(data) > 0 else 0,
                        'mean_contacts': np.mean(data) if len(data) > 0 else 0,
                        'simulation_length_ns': times[-1] if len(times) > 0 else 0
                    })
            
            if summary_data:
                summary_df = pd.DataFrame(summary_data)
                summary_df = summary_df.sort_values('molecular_weight')
                summary_path = os.path.join(output_dir, f"{analysis_type}_summary.csv")
                summary_df.to_csv(summary_path, index=False, float_format='%.2f')
                logger.warning(f"Multi-system summary saved: {summary_path}")
            
            return data_path, summary_path
        except Exception as e:
            logger.error(f"Failed to save multi-system data: {e}")
            return None, None

    @staticmethod
    def save_frame_by_frame_contacts(analyzer, output_path: str, selected_residues=None):
        try:
            if not analyzer.traj or not analyzer._new_method_results:
                logger.warning("No trajectory or contact results available")
                return None
                
            ligand_atoms = analyzer._new_method_results.get('ligand_atoms', [])
            contact_frequencies = analyzer._new_method_results.get('contact_frequencies', {})
            
            if not contact_frequencies:
                logger.warning("No contact data available")
                return None

            if selected_residues:
                filtered_contacts = {}
                for (lig_atom, res_id), freq in contact_frequencies.items():
                    formatted_res_id = f"{res_id[:3]} {res_id[3:]}"  
                    if any(selected_res in formatted_res_id for selected_res in selected_residues):
                        filtered_contacts[(lig_atom, res_id)] = freq
                contact_frequencies = filtered_contacts
            
            if not contact_frequencies:
                logger.warning("No contacts found for selected residues")
                return None

            atom_pairs = []
            pair_info = []

            protein_atom_map = {}  
            for atom in analyzer.traj.topology.atoms:
                if atom.residue.name in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                                       'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                                       'THR', 'TRP', 'TYR', 'VAL'] and atom.element.symbol != 'H':
                    res_id = f"{atom.residue.name}{atom.residue.resSeq}"
                    if res_id not in protein_atom_map:
                        protein_atom_map[res_id] = []
                    protein_atom_map[res_id].append(atom.index)

            for (lig_atom, res_id), freq in contact_frequencies.items():
                if res_id in protein_atom_map:
                    for prot_atom in protein_atom_map[res_id]:
                        atom_pairs.append([lig_atom, prot_atom])
                        pair_info.append((lig_atom, prot_atom, res_id, freq))
            
            if not atom_pairs:
                logger.warning("No atom pairs found for analysis")
                return None

            distances = md.compute_distances(analyzer.traj, atom_pairs, opt=True)
            contact_threshold = analyzer.config.contact_enter_threshold_nm
            
            frame_data = []
            for frame_idx in range(analyzer.traj.n_frames):
                frame_time = frame_idx * 0.5  
                
                for pair_idx, (lig_atom, prot_atom, res_id, expected_freq) in enumerate(pair_info):
                    distance_nm = distances[frame_idx, pair_idx]
                    distance_angstrom = distance_nm * 10.0
                    in_contact = distance_nm < contact_threshold

                    lig_atom_obj = analyzer.traj.topology.atom(lig_atom)
                    prot_atom_obj = analyzer.traj.topology.atom(prot_atom)
                    
                    frame_data.append({
                        'frame': frame_idx,
                        'time_ns': frame_time,
                        'ligand_atom_index': lig_atom,
                        'ligand_atom_name': lig_atom_obj.name,
                        'protein_atom_index': prot_atom,
                        'protein_atom_name': prot_atom_obj.name,
                        'protein_residue_name': prot_atom_obj.residue.name,
                        'protein_residue_number': prot_atom_obj.residue.resSeq,
                        'protein_residue_id': res_id,
                        'distance_nm': distance_nm,
                        'distance_angstrom': distance_angstrom,
                        'in_contact': in_contact,
                        'expected_contact_frequency': expected_freq
                    })

            df = pd.DataFrame(frame_data)
            df.to_csv(output_path, index=False, float_format='%.4f')
            logger.warning(f"Frame-by-frame contact data saved: {output_path} ({len(frame_data)} records)")
            
            return df
            
        except Exception as e:
            logger.error(f"Failed to save frame-by-frame contact data: {e}")
            return None
    
    @staticmethod
    def save_frame_by_frame_hbonds(analyzer, output_path: str, selected_hbonds=None):
        """Optimized frame-by-frame hydrogen bond data export"""
        try:
            if not hasattr(analyzer.hbond_analyzer, 'traj') or analyzer.hbond_analyzer.traj is None:
                logger.warning("No hydrogen bond analysis data available")
                return None
            
            hbond_analyzer = analyzer.hbond_analyzer
            traj = hbond_analyzer.traj

            hbond_pairs, hbond_triplets = hbond_analyzer._prepare_hbond_pairs()
        
            if not hbond_pairs:
                logger.warning("No hydrogen bond pairs found")
                return None

            if selected_hbonds:
                selected_keys = set(hb[0] for hb in selected_hbonds)
            
                filtered_pairs = []
                filtered_triplets = []
                filtered_keys = []

                for idx, pair in enumerate(hbond_pairs):
                    donor_idx, acceptor_idx = pair
                    donor_atom = traj.topology.atom(donor_idx)
                    acceptor_atom = traj.topology.atom(acceptor_idx)

                    if donor_idx in hbond_analyzer.protein_donors:
                        hbond_key = f"{donor_atom.residue.name} {donor_atom.residue.resSeq} ({donor_atom.name}) -> Ligand"
                    else:
                        hbond_key = f"Ligand -> {acceptor_atom.residue.name} {acceptor_atom.residue.resSeq} ({acceptor_atom.name})"

                    if hbond_key in selected_keys:
                        filtered_pairs.append(pair)
                        filtered_triplets.append(hbond_triplets[idx])
                        filtered_keys.append(hbond_key)

                hbond_pairs = filtered_pairs
                hbond_triplets = filtered_triplets
                hbond_keys = filtered_keys

                if not hbond_pairs:
                    logger.warning("No hydrogen bond pairs found for selected bonds")
                    return None
            else:
                hbond_keys = []
                for pair in hbond_pairs:
                    donor_idx, acceptor_idx = pair
                    donor_atom = traj.topology.atom(donor_idx)
                    acceptor_atom = traj.topology.atom(acceptor_idx)

                    if donor_idx in hbond_analyzer.protein_donors:
                        hbond_key = f"{donor_atom.residue.name} {donor_atom.residue.resSeq} ({donor_atom.name}) -> Ligand"
                    else:
                        hbond_key = f"Ligand -> {acceptor_atom.residue.name} {acceptor_atom.residue.resSeq} ({acceptor_atom.name})"

                    hbond_keys.append(hbond_key)

            logger.warning(f"Computing hydrogen bonds for {len(hbond_pairs)} pairs across {traj.n_frames} frames")

            try:
                all_distances = md.compute_distances(traj, hbond_pairs, periodic=False, opt=True)

                all_angles_rad = md.compute_angles(traj, hbond_triplets, periodic=False, opt=True)
                all_angles_deg = np.degrees(all_angles_rad)

                distance_mask = all_distances < hbond_analyzer.config.hbond_distance_cutoff_nm
                angle_mask = all_angles_deg > hbond_analyzer.config.hbond_angle_cutoff_deg
                is_hbond_matrix = distance_mask & angle_mask  

                logger.warning(f"Batch computation complete. Processing results...")

            except Exception as e:
                logger.error(f"Error in batch computation: {e}")
                return None

            frame_data = []
            n_frames, n_pairs = all_distances.shape

            pair_info = []
            for idx, (pair, key) in enumerate(zip(hbond_pairs, hbond_keys)):
                donor_idx, acceptor_idx = pair
                donor_atom = traj.topology.atom(donor_idx)
                acceptor_atom = traj.topology.atom(acceptor_idx)
            
                if donor_idx in hbond_analyzer.protein_donors:
                    donor_type = "protein"
                    acceptor_type = "ligand"
                else:
                    donor_type = "ligand"
                    acceptor_type = "protein"
            
                pair_info.append({
                'donor_idx': donor_idx,
                'acceptor_idx': acceptor_idx,
                'donor_name': donor_atom.name,
                'donor_residue_name': donor_atom.residue.name,
                'donor_residue_number': donor_atom.residue.resSeq,
                'donor_type': donor_type,
                'acceptor_name': acceptor_atom.name,
                'acceptor_residue_name': acceptor_atom.residue.name,
                'acceptor_residue_number': acceptor_atom.residue.resSeq,
                'acceptor_type': acceptor_type,
                'hbond_key': key
            })

            frame_indices = np.arange(n_frames)
            pair_indices = np.arange(n_pairs)

            frame_grid, pair_grid = np.meshgrid(frame_indices, pair_indices, indexing='ij')

            flat_frames = frame_grid.flatten()
            flat_pairs = pair_grid.flatten()
            flat_times = flat_frames * 0.5  
        
            flat_distances_nm = all_distances.flatten()
            flat_distances_angstrom = flat_distances_nm * 10.0
            flat_angles = all_angles_deg.flatten()
            flat_is_hbond = is_hbond_matrix.flatten()

            for i in range(len(flat_frames)):
                frame_idx = flat_frames[i]
                pair_idx = flat_pairs[i]
                info = pair_info[pair_idx]
            
                frame_data.append({
                'frame': frame_idx,
                'time_ns': flat_times[i],
                'donor_atom_index': info['donor_idx'],
                'donor_atom_name': info['donor_name'],
                'donor_residue_name': info['donor_residue_name'],
                'donor_residue_number': info['donor_residue_number'],
                'donor_type': info['donor_type'],
                'acceptor_atom_index': info['acceptor_idx'],
                'acceptor_atom_name': info['acceptor_name'],
                'acceptor_residue_name': info['acceptor_residue_name'],
                'acceptor_residue_number': info['acceptor_residue_number'],
                'acceptor_type': info['acceptor_type'],
                'hbond_key': info['hbond_key'],
                'distance_nm': flat_distances_nm[i],
                'distance_angstrom': flat_distances_angstrom[i],
                'angle_degrees': flat_angles[i],
                'is_hbond': flat_is_hbond[i]
            })
        
            if not frame_data:
                logger.warning("No hydrogen bond data generated")
                return None

            df = pd.DataFrame(frame_data)
            df.to_csv(output_path, index=False, float_format='%.4f')
        
            if selected_hbonds:
                logger.warning(f"Frame-by-frame hydrogen bond data saved: {output_path} ({len(frame_data)} records for {len(selected_hbonds)} selected H-bonds)")
            else:
                logger.warning(f"Frame-by-frame hydrogen bond data saved: {output_path} ({len(frame_data)} records)")
        
            return df
        
        except Exception as e:
            logger.error(f"Failed to save frame-by-frame hydrogen bond data: {e}")
            return None

    @staticmethod 
    def save_frame_by_frame_distances(analyzer, output_path: str, selected_residues=None, max_residues=5):
        """Save frame-by-frame distance data"""
        try:
            if not analyzer.traj:
                logger.warning("No trajectory data available")
                return None

            if hasattr(analyzer, '_new_method_results') and analyzer._new_method_results:
                contact_proportions = analyzer._new_method_results.get('residue_proportions', {})
                residue_avg_distances = analyzer._new_method_results.get('residue_avg_distances', {})
            else:
                logger.warning("No contact results available for distance analysis")
                return None

            if selected_residues:
                residues_to_analyze = []
                for sel_res in selected_residues:
                    found = False
                    for contact_res_id, prop in contact_proportions.items():
                        if ' ' in sel_res:
                            formatted_sel = sel_res.replace(' ', '')
                        else:
                            formatted_sel = sel_res
                        
                        if contact_res_id == formatted_sel or sel_res == contact_res_id:
                            residues_to_analyze.append((sel_res, prop))
                            found = True
                            break
                        
                    if not found:
                        logger.warning(f"Selected residue {sel_res} not found in contact results")
                    
            else:
                if residue_avg_distances:
                    sorted_by_distance = sorted(residue_avg_distances.items(), key=lambda x: x[1])
                    residues_to_analyze = []
                    for residue_id, avg_dist in sorted_by_distance[:max_residues]:
                        contact_prop = contact_proportions.get(residue_id, 0.0)
                        residues_to_analyze.append((residue_id, contact_prop))
                elif contact_proportions:
                    top_contacts = sorted(contact_proportions.items(), key=lambda x: x[1], reverse=True)
                    residues_to_analyze = top_contacts[:max_residues]
                else:
                    logger.warning("No contact proportions or distance data available")
                    return None

            if not residues_to_analyze:
                logger.warning("No residues selected for distance analysis")
                return None

            ligand_atoms = analyzer._new_method_results.get('ligand_atoms', [])
            if not ligand_atoms:
                logger.warning("No ligand atoms found")
                return None

            frame_data = []

            for frame_idx in range(analyzer.traj.n_frames):
                frame_time = frame_idx * 0.5
                frame_coords = analyzer.traj.xyz[frame_idx]  

                for residue_id, contact_prop in residues_to_analyze:
                    if isinstance(residue_id, str):
                        import re
                        match = re.match(r'([A-Z]{3})(\d+)', residue_id.replace(' ', ''))
                        if match:
                            res_name = match.group(1)
                            res_num = int(match.group(2))
                        else:
                            logger.warning(f"Could not parse residue ID: {residue_id}")
                            continue
                    else:
                        continue

                    residue_atoms = []
                    for atom in analyzer.traj.topology.atoms:
                        if (atom.residue.name == res_name and 
                            atom.residue.resSeq == res_num and 
                            atom.element.symbol != 'H'):
                            residue_atoms.append(atom.index)

                    if not residue_atoms:
                        continue

                    min_distances = []
                    for lig_atom in ligand_atoms:
                        lig_coord = frame_coords[lig_atom]
                        atom_distances = []
                        for res_atom in residue_atoms:
                            res_coord = frame_coords[res_atom]
                            dist = np.linalg.norm(lig_coord - res_coord)  
                            atom_distances.append(dist)
                        min_distances.append(min(atom_distances))

                    ligand_residue_min_dist = min(min_distances)
                    ligand_residue_min_dist_angstrom = ligand_residue_min_dist * 10.0

                    ligand_center = np.mean([frame_coords[atom] for atom in ligand_atoms], axis=0)
                    residue_center = np.mean([frame_coords[atom] for atom in residue_atoms], axis=0)
                    center_distance = np.linalg.norm(ligand_center - residue_center)
                    center_distance_angstrom = center_distance * 10.0

                    frame_data.append({
                        'frame': frame_idx,
                        'time_ns': frame_time,
                        'residue_id': f"{res_name} {res_num}",  
                        'residue_name': res_name,
                        'residue_number': res_num,
                        'contact_proportion': contact_prop,
                        'min_distance_nm': ligand_residue_min_dist,
                        'min_distance_angstrom': ligand_residue_min_dist_angstrom,
                        'center_distance_nm': center_distance,
                        'center_distance_angstrom': center_distance_angstrom,
                        'in_contact': ligand_residue_min_dist < analyzer.config.contact_enter_threshold_nm
                    })

            if not frame_data:
                logger.warning("No frame data generated")
                return None

            df = pd.DataFrame(frame_data)
            df.to_csv(output_path, index=False, float_format='%.4f')
            logger.warning(f"Frame-by-frame distance data saved: {output_path} ({len(frame_data)} records)")

            return df

        except Exception as e:
            logger.error(f"Failed to save frame-by-frame distance data: {e}")
            return None