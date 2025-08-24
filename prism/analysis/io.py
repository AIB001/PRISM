import os
import json
import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any

from .config import AnalysisConfig

logger = logging.getLogger(__name__)

def nm_to_angstrom(value):
    """Convert nanometers to angstroms"""
    return value * 10.0

class DataExporter:
    """Export analysis data to various formats"""
    
    @staticmethod
    def save_contact_data_csv(contact_proportions: Dict[str, float], output_path: str):
        """Save contact data to CSV"""
        try:
            df = pd.DataFrame([
                {
                    'residue_name': key.split()[0],
                    'residue_number': int(key.split()[1]),
                    'residue_id': key,
                    'contact_proportion': value,
                    'contact_percentage': value * 100
                }
                for key, value in contact_proportions.items()
            ])
            df = df.sort_values('contact_proportion', ascending=False)
            df.to_csv(output_path, index=False, float_format='%.4f')
            logger.info(f"Contact data saved: {output_path}")
            return df
        except Exception as e:
            logger.error(f"Failed to save contact data: {e}")
            return None
    
    @staticmethod
    def save_hbond_data_csv(hbond_frequencies: List[Tuple], output_path: str, hbond_stats: Dict):
        """Save hydrogen bond data to CSV"""
        try:
            if not hbond_frequencies:
                logger.warning("No hydrogen bond data to save")
                return None
            
            data = []
            for hb in hbond_frequencies:
                key, freq, avg_dist, avg_angle = hb
                
                if ' -> ' in key:
                    donor_part, acceptor_part = key.split(' -> ')
                elif ' → ' in key:
                    donor_part, acceptor_part = key.split(' → ')
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
            logger.info(f"Hydrogen bond data saved: {output_path}")
            return df
        except Exception as e:
            logger.error(f"Failed to save hydrogen bond data: {e}")
            return None
    
    @staticmethod
    def save_distance_data_csv(distance_stats: Dict[str, Dict[str, float]], output_path: str):
        """Save distance statistics to CSV"""
        try:
            data = []
            for residue_id, stats in distance_stats.items():
                data.append({
                    'residue_id': residue_id,
                    'mean_distance_angstrom': stats['mean'],
                    'min_distance_angstrom': stats['min'],
                    'max_distance_angstrom': stats['max'],
                    'std_distance_angstrom': stats['std'],
                    'contact_proportion': stats.get('contact_proportion', 0.0)
                })
            
            df = pd.DataFrame(data)
            df.to_csv(output_path, index=False, float_format='%.4f')
            logger.info(f"Distance data saved: {output_path}")
            return df
        except Exception as e:
            logger.error(f"Failed to save distance data: {e}")
            return None
    
    @staticmethod
    def save_detailed_contact_analysis(residue_id: str, distances: np.ndarray, 
                                     contact_states: np.ndarray, output_path: str,
                                     times: np.ndarray, config: AnalysisConfig):
        """Save detailed contact analysis for a specific residue"""
        try:
            distances_angstrom = distances * 10.0
            
            contact_fraction = np.mean(contact_states)
            
            contact_periods = []
            in_contact = False
            start_time = None
            
            for i, state in enumerate(contact_states):
                if state == 1 and not in_contact:
                    in_contact = True
                    start_time = times[i]
                elif state == 0 and in_contact:
                    in_contact = False
                    if start_time is not None:
                        contact_periods.append((start_time, times[i], times[i] - start_time))
            
            if in_contact and start_time is not None:
                contact_periods.append((start_time, times[-1], times[-1] - start_time))
            
            data = {
                'time_ns': times,
                'distance_angstrom': distances_angstrom,
                'in_contact': contact_states,
                'residue_id': residue_id
            }
            
            if len(distances_angstrom) > config.min_frames_for_smoothing:
                from .visualization import moving_average
                smoothed_distances = moving_average(distances_angstrom, config.smooth_window)
                smoothed_times = times[config.smooth_window-1:]
                
                padded_smooth = np.full(len(times), np.nan)
                padded_smooth[config.smooth_window-1:] = smoothed_distances
                data['distance_smoothed_angstrom'] = padded_smooth
            
            df = pd.DataFrame(data)
            df.to_csv(output_path, index=False, float_format='%.4f')
            
            summary_path = output_path.replace('.csv', '_summary.json')
            summary = {
                'residue_id': residue_id,
                'contact_fraction': float(contact_fraction),
                'distance_stats': {
                    'min_angstrom': float(np.min(distances_angstrom)),
                    'max_angstrom': float(np.max(distances_angstrom)),
                    'mean_angstrom': float(np.mean(distances_angstrom)),
                    'std_angstrom': float(np.std(distances_angstrom))
                },
                'contact_periods': [
                    {'start_ns': float(start), 'end_ns': float(end), 'duration_ns': float(duration)}
                    for start, end, duration in contact_periods
                ],
                'total_contact_time_ns': float(sum(duration for _, _, duration in contact_periods)),
                'number_of_contact_events': len(contact_periods),
                'simulation_length_ns': float(times[-1] - times[0])
            }
            
            with open(summary_path, 'w') as f:
                json.dump(summary, f, indent=2)
            
            logger.info(f"Detailed contact analysis saved: {output_path}")
            return df, summary
        except Exception as e:
            logger.error(f"Failed to save detailed contact analysis: {e}")
            return None, None
    
    @staticmethod
    def save_system_info(trajectory_manager, output_path: str):
        """Save system information"""
        try:
            system_info = trajectory_manager.get_system_info()
            system_info['analysis_parameters'] = trajectory_manager.config.__dict__
            
            with open(output_path, 'w') as f:
                json.dump(system_info, f, indent=2)
            
            logger.info(f"System information saved: {output_path}")
            return system_info
        except Exception as e:
            logger.error(f"Failed to save system information: {e}")
            return None
    
    @staticmethod
    def save_analysis_summary(contact_proportions: Dict[str, float], 
                            hbond_frequencies: List[Tuple], 
                            distance_stats: Dict[str, Dict[str, float]],  
                            output_path: str, config: AnalysisConfig, trajectory_manager):
        """Save analysis summary"""
        try:
            contact_stats = {
                'total_residues_analyzed': len(contact_proportions),
                'residues_with_contact': sum(1 for prop in contact_proportions.values() if prop > 0),
                'residues_with_significant_contact': sum(1 for prop in contact_proportions.values() if prop > 0.1),
                'max_contact_proportion': max(contact_proportions.values()) if contact_proportions else 0,
                'mean_contact_proportion': np.mean(list(contact_proportions.values())) if contact_proportions else 0,
                'top_10_contacts': dict(sorted(contact_proportions.items(), key=lambda x: x[1], reverse=True)[:10])
            }
            
            hbond_stats = {
                'total_hbonds_detected': len(hbond_frequencies),
                'significant_hbonds': sum(1 for hb in hbond_frequencies if hb[1] > 0.05),
                'max_hbond_frequency': max(hb[1] for hb in hbond_frequencies) if hbond_frequencies else 0,
                'mean_hbond_frequency': np.mean([hb[1] for hb in hbond_frequencies]) if hbond_frequencies else 0,
                'top_5_hbonds': [(hb[0], hb[1], hb[2], hb[3]) for hb in hbond_frequencies[:5]]
            }
            
            distance_stats_summary = {
                'residues_with_close_contact': sum(1 for stats in distance_stats.values() if stats['mean'] < nm_to_angstrom(config.distance_cutoff_nm)),
                'mean_min_distance': np.mean([stats['min'] for stats in distance_stats.values()]) if distance_stats else 0,
                'mean_mean_distance': np.mean([stats['mean'] for stats in distance_stats.values()]) if distance_stats else 0,
                'mean_max_distance': np.mean([stats['max'] for stats in distance_stats.values()]) if distance_stats else 0,
                'closest_residue': min(distance_stats.items(), key=lambda x: x[1]['min'])[0] if distance_stats else "None",
                'top_5_closest': dict(sorted(distance_stats.items(), key=lambda x: x[1]['min'])[:5])
            }
            
            summary = {
                'analysis_summary': {
                    'timestamp': pd.Timestamp.now().isoformat(),
                    'configuration': config.__dict__,
                    'contact_analysis': contact_stats,
                    'hydrogen_bond_analysis': hbond_stats,
                    'distance_analysis': distance_stats_summary
                }
            }
            
            with open(output_path, 'w') as f:
                json.dump(summary, f, indent=2)
            
            logger.info(f"Analysis summary saved: {output_path}")
            return summary
        except Exception as e:
            logger.error(f"Failed to save analysis summary: {e}")
            return None


class ReportGenerator:
    """Generate comprehensive analysis reports"""
    
    @staticmethod
    def generate_comprehensive_report(contact_proportions: Dict[str, float],
                                hbond_frequencies: List[Tuple],
                                distance_stats: Dict[str, Dict[str, float]],
                                output_dir: str, config: AnalysisConfig, 
                                trajectory_manager, hbond_stats: Dict):
        """Generate comprehensive report"""
        report_path = os.path.join(output_dir, "comprehensive_report.txt")
    
        try:
            with open(report_path, "w", encoding="utf-8") as f:
                f.write("=" * 80 + "\n")
                f.write("PROTEIN-LIGAND INTERACTION ANALYSIS RESULTS\n")
                f.write("=" * 80 + "\n\n")
            
                if trajectory_manager and trajectory_manager.traj is not None:
                    traj = trajectory_manager.traj
                    f.write("SYSTEM INFORMATION:\n")
                    f.write(f"  Trajectory frames: {traj.n_frames}\n")
                    total_time_ns = traj.time[-1] / 1000
                    f.write(f"  Total simulation time: {total_time_ns:.2f} ns\n")
                    f.write(f"  Ligand atoms: {len(trajectory_manager.ligand_atoms)}\n")
                    f.write(f"  Protein atoms: {len(trajectory_manager.protein_atoms)}\n")
                    f.write(f"  Protein residues: {len(trajectory_manager.protein_residues)}\n\n")
            
                f.write("CONTACT ANALYSIS RESULTS:\n")
                total_residues = len(contact_proportions)
                residues_with_contact = sum(1 for prop in contact_proportions.values() if prop > 0)
                residues_significant = sum(1 for prop in contact_proportions.values() if prop > 0.1)
                max_contact = max(contact_proportions.values()) if contact_proportions else 0
                mean_contact = np.mean(list(contact_proportions.values())) if contact_proportions else 0
            
                f.write(f"  Total residues analyzed: {total_residues}\n")
                f.write(f"  Residues with any contact: {residues_with_contact}\n")
                f.write(f"  Residues with >10% contact: {residues_significant}\n")
                f.write(f"  Maximum contact proportion: {max_contact * 100:.2f}%\n")
                f.write(f"  Mean contact proportion: {mean_contact * 100:.2f}%\n\n")
            
                f.write("  TOP 10 CONTACTING RESIDUES:\n")
                top_contacts = sorted(contact_proportions.items(), 
                                     key=lambda x: x[1], reverse=True)[:10]
                for i, (residue, proportion) in enumerate(top_contacts, 1):
                    f.write(f"    {i:2d}. {residue:<15} {proportion * 100:>8.2f}%\n")
                f.write("\n")
            
                if config.distance_analysis and distance_stats:
                    f.write("DISTANCE ANALYSIS RESULTS:\n")
                    cutoff_angstrom = nm_to_angstrom(config.distance_cutoff_nm)
                    residues_close = sum(1 for stats in distance_stats.values() 
                                        if stats['mean'] < cutoff_angstrom)
                
                    all_min_dists = [stats['min'] for stats in distance_stats.values()]
                    all_mean_dists = [stats['mean'] for stats in distance_stats.values()]
                    all_max_dists = [stats['max'] for stats in distance_stats.values()]
                
                    closest_residue, closest_stats = min(distance_stats.items(), 
                                                        key=lambda x: x[1]['min'], 
                                                        default=("None", {}))
                
                    f.write(f"  Residues with mean distance < {config.distance_cutoff_nm:.2f} nm: "
                           f"{residues_close}\n")
                    f.write(f"  Average minimum distance: {np.mean(all_min_dists):.2f} Å\n")
                    f.write(f"  Average mean distance: {np.mean(all_mean_dists):.2f} Å\n")
                    f.write(f"  Average maximum distance: {np.mean(all_max_dists):.2f} Å\n")
                    f.write(f"  Closest residue: {closest_residue} (min distance: {closest_stats.get('min', 0):.2f} Å)\n\n")
                
                    f.write("  TOP 5 CLOSEST RESIDUES:\n")
                    closest_residues = sorted(distance_stats.items(), 
                                             key=lambda x: x[1]['min'])[:5]
                    for i, (residue, stats) in enumerate(closest_residues, 1):
                        f.write(f"    {i:2d}. {residue:<15} Min: {stats['min']:.2f} Å  Mean: {stats['mean']:.2f} Å\n")
                    f.write("\n")
            
                f.write("HYDROGEN BOND ANALYSIS RESULTS:\n")
                total_hbonds = len(hbond_frequencies)
                significant_hbonds = sum(1 for hb in hbond_frequencies if hb[1] > 0.05)
                max_hb_freq = max(hb[1] for hb in hbond_frequencies) if hbond_frequencies else 0
                mean_hb_freq = np.mean([hb[1] for hb in hbond_frequencies]) if hbond_frequencies else 0
            
                f.write(f"  Total H-bonds detected: {total_hbonds}\n")
                f.write(f"  Significant H-bonds (>5%): {significant_hbonds}\n")
                f.write(f"  Maximum H-bond frequency: {max_hb_freq * 100:.2f}%\n")
                f.write(f"  Mean H-bond frequency: {mean_hb_freq * 100:.2f}%\n")
            
                if hbond_frequencies:
                    f.write("\n  TOP HYDROGEN BONDS:\n")
                    for i, hb in enumerate(hbond_frequencies[:10], 1):
                        key, freq, avg_dist, avg_angle = hb
                        f.write(f"    {i:2d}. {key}\n")
                        f.write(f"        Frequency: {freq * 100:>8.2f}%  Distance: {avg_dist:>6.2f} Å  Angle: {avg_angle:>6.1f}°\n")
                
                    f.write("\n  DETAILED HYDROGEN BOND STATISTICS:\n")
                    for hb in hbond_frequencies[:5]:
                        key, freq, avg_dist, avg_angle = hb
                        if key in hbond_stats:
                            stats = hbond_stats[key]
                            f.write(f"\n    {key}:\n")
                            f.write(f"      Total occurrences: {stats['count']}\n")
                            if stats['distance']:
                                dist_mean = np.mean(stats['distance'])
                                dist_std = np.std(stats['distance'])
                                dist_min = np.min(stats['distance'])
                                dist_max = np.max(stats['distance'])
                                f.write(f"      Distance: {dist_mean:.2f} ± {dist_std:.2f} Å\n")
                                f.write(f"      Distance range: {dist_min:.2f} - {dist_max:.2f} Å\n")
                            if stats['angle']:
                                angle_mean = np.mean(stats['angle'])
                                angle_std = np.std(stats['angle'])
                                angle_min = np.min(stats['angle'])
                                angle_max = np.max(stats['angle'])
                                f.write(f"      Angle: {angle_mean:.1f} ± {angle_std:.1f}°\n")
                                f.write(f"      Angle range: {angle_min:.1f} - {angle_max:.1f}°\n")

                f.write("\n" + "=" * 80 + "\n")
        
            logger.info(f"Comprehensive report generated: {os.path.basename(report_path)}")
            return report_path
        except Exception as e:
            logger.error(f"Failed to generate comprehensive report: {e}")
            return None # pyright: ignore[reportShadowedImports]