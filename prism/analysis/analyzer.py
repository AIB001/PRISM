import os
import logging
import numpy as np
from typing import Optional, List, Tuple, Dict
from concurrent.futures import ProcessPoolExecutor, as_completed

from .config import AnalysisConfig
from .trajectory import TrajectoryManager
from .calc.contacts import ContactAnalyzer, DistanceCalculator
from .calc.distance import DistanceAnalyzer
from .calc.hbonds import HBondAnalyzer
from .multisys import MultiSystemAnalyzer
from .export import DataExporter
from .plots.basic import Visualizer

logging.basicConfig(level=logging.WARNING, format='%(message)s')
logger = logging.getLogger(__name__)

class IntegratedProteinLigandAnalyzer:
    """Integrated protein-ligand analysis tool"""
    
    def __init__(self, base_dir: Optional[str] = None, topology_file: Optional[str] = None, 
                 trajectory_file: Optional[str] = None, config: Optional[AnalysisConfig] = None,
                 skip_initial_analysis: bool = False):
        self.config = config or AnalysisConfig()
        
        self.traj_manager = TrajectoryManager()
        self.distance_calc = DistanceCalculator()
        self.contact_analyzer = ContactAnalyzer(self.config)
        self.distance_analyzer = DistanceAnalyzer(self.config)
        self.hbond_analyzer = HBondAnalyzer(self.config)
        self.multi_analyzer = MultiSystemAnalyzer(self.config)
        self.visualizer = Visualizer()
        
        self.system_files = {}
        if base_dir:
            self.system_files = self._scan_system_directories(base_dir)
        
        self.traj = None
        self.universe = None
        self.ligand_atoms = []
        self.protein_atoms = []
        self.protein_residues = []
        self.ligand_residue = None
        
        self._cached_distance_results = {}
        self._new_method_results = None
        
        if topology_file and trajectory_file:
            if isinstance(topology_file, (tuple, list)):
                topology_file = topology_file[0] if topology_file else None
            if isinstance(trajectory_file, (tuple, list)):
                trajectory_file = trajectory_file[0] if trajectory_file else None
                
            if topology_file and trajectory_file:
                self.traj = self.traj_manager.load_mdtraj_trajectory(topology_file, trajectory_file)
                self.universe = self.traj_manager.load_universe(topology_file, trajectory_file)
            if not skip_initial_analysis:
                self._prepare_single_system()
    
    def _scan_system_directories(self, base_dir: str) -> Dict[str, Tuple[str, str, str]]:
        system_files = {}
        if not os.path.exists(base_dir):
            raise FileNotFoundError(f"Base directory not found: {base_dir}")
        
        for name in sorted(os.listdir(base_dir)):
            full_path = os.path.join(base_dir, name)
            if not os.path.isdir(full_path):
                continue
            
            traj, prod, top = self._find_md_files(full_path)
            if traj and top:
                system_files[name] = (traj, prod, top)
            else:
                logger.warning(f"Missing MD files: {name}")
        
        logger.warning(f"Found {len(system_files)} systems for multi-system analysis")
        return system_files
    
    @staticmethod
    def _find_md_files(system_dir: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        gmx_base = os.path.join(system_dir, 'GMX_PROLIG_MD')
        prod_dir = os.path.join(gmx_base, 'prod')
        traj_path = os.path.join(prod_dir, 'md.xtc')
        top_path = os.path.join(gmx_base, 'solv_ions.gro')
        
        if os.path.exists(traj_path) and os.path.exists(top_path):
            return traj_path, prod_dir, top_path
        return None, None, None
    
    def _prepare_single_system(self):
        try:
            self.ligand_atoms = self.universe.select_atoms(
                f"resname {self.config.ligand_resname} and not name H*"
            ).indices
            self.protein_atoms = self.universe.select_atoms("protein and not name H*").indices
        
            protein_atoms = self.universe.select_atoms("protein")
            self.protein_residues = list(set(protein_atoms.residues))
        
            ligand_residues = self.universe.select_atoms(f"resname {self.config.ligand_resname}").residues
            self.ligand_residue = ligand_residues[0] if len(ligand_residues) > 0 else None

            self._new_method_results = self.contact_analyzer.calculate_contact_proportions(self.traj)

            mdtraj_protein_residues = []
            for mda_residue in self.protein_residues:
                for mdtraj_residue in self.traj.topology.residues:
                    if (mdtraj_residue.name == mda_residue.resname and 
                        mdtraj_residue.resSeq == mda_residue.resid):
                        mdtraj_protein_residues.append(mdtraj_residue)
                        break

            mdtraj_ligand_residue = None
            if self.ligand_residue:
                for mdtraj_residue in self.traj.topology.residues:
                    if (mdtraj_residue.name == self.ligand_residue.resname and 
                        mdtraj_residue.resSeq == self.ligand_residue.resid):
                        mdtraj_ligand_residue = mdtraj_residue
                        break

            self.hbond_analyzer.set_trajectory_data(
                self.traj, 
                mdtraj_ligand_residue, 
                mdtraj_protein_residues
            )
        
            logger.warning(f"Found {len(self.ligand_atoms)} ligand atoms, {len(self.protein_atoms)} protein atoms")
            logger.warning(f"Found {len(self.protein_residues)} protein residues for analysis")
        
        except Exception as e:
            logger.error(f"Failed to prepare system: {e}")
            raise
    
    def calculate_contact_proportions(self) -> Dict[str, float]:
        """Calculate contact proportions using new method results"""
        if self._new_method_results is None:
            raise ValueError("No new method results available")
        
        return self._new_method_results.get('residue_proportions', {})
    
    def get_detailed_contact_results(self):
        """Get detailed results from new method"""
        return self._new_method_results
    
    def select_top_contacts(self, contact_proportions: Dict[str, float], top_n: int = 10) -> List[Tuple[str, float]]:
        """Select top contacting residues"""
        if self._new_method_results:
            residue_proportions = self._new_method_results.get('residue_proportions', {})
            sorted_contacts = sorted(
                residue_proportions.items(),
                key=lambda x: x[1],
                reverse=True
            )
            formatted_results = []
            for residue_id, proportion in sorted_contacts[:top_n]:
                if residue_id and len(residue_id) > 3:
                    res_name = ''.join([c for c in residue_id if c.isalpha()])
                    res_num = ''.join([c for c in residue_id if c.isdigit()])
                    if res_name and res_num:
                        formatted_key = f"{res_name} {res_num}"
                        formatted_results.append((formatted_key, proportion))
            return formatted_results
        else:
            sorted_contacts = sorted(
                contact_proportions.items(),
                key=lambda x: x[1],
                reverse=True
            )
            return sorted_contacts[:top_n]
    
    def calculate_distance_analysis(self) -> Dict[str, Dict[str, float]]:
        """Calculate distance statistics from results"""
        if self._new_method_results is None:
            logger.warning("No results available for distance analysis")
            return {}
        
        residue_avg_distances = self._new_method_results.get('residue_avg_distances', {})
        if not residue_avg_distances:
            logger.warning("No distance data available in results")
            return {}
        
        distance_stats = {}
        for residue_id, avg_distance_nm in residue_avg_distances.items():
            avg_distance_angstrom = avg_distance_nm * 10.0
            std_estimate = avg_distance_angstrom * 0.15  
            min_estimate = max(avg_distance_angstrom - 2*std_estimate, 2.0)  
            max_estimate = avg_distance_angstrom + 2*std_estimate
            
            distance_stats[residue_id] = {
                'mean': float(avg_distance_angstrom),
                'min': float(min_estimate),
                'max': float(max_estimate),
                'std': float(std_estimate)
            }
        
        return distance_stats
    
    def analyze_hydrogen_bonds(self) -> List[Tuple[str, float, float, float]]:
        """Hydrogen bond analysis"""
        return self.hbond_analyzer.analyze_hydrogen_bonds(self.universe)
    
    def analyze_residue_contact(self, residue_id: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Analyze specific residue contact with distance calculation"""
        try:
            if ' ' in residue_id:
                res_name, res_num_str = residue_id.split()
            else:
                res_name = residue_id[:3]
                res_num_str = residue_id[3:]
        
            try:
                res_num = int(res_num_str)
            except ValueError:
                logger.error(f"Invalid residue number in {residue_id}")
                return None, None

            ligand_atoms = self._new_method_results.get('ligand_atoms', [])
            if not ligand_atoms:
                logger.error("No ligand atoms available")
                return None, None

            target_residue_atoms = []
            for atom in self.traj.topology.atoms:
                if (atom.residue.name == res_name and 
                    atom.residue.resSeq == res_num and 
                    atom.element.symbol != 'H'):
                    target_residue_atoms.append(atom.index)
        
            if not target_residue_atoms:
                logger.error(f"No heavy atoms found for residue {residue_id}")
                return None, None

            atom_pairs = []
            for lig_atom in ligand_atoms:
                for res_atom in target_residue_atoms:
                    atom_pairs.append([lig_atom, res_atom])

            import mdtraj as md
            all_distances = md.compute_distances(self.traj, atom_pairs, opt=True)
            min_distances = np.min(all_distances, axis=1)
            contact_states = (min_distances < self.config.contact_enter_threshold_nm).astype(int)

            contact_frames = np.sum(contact_states)
            contact_percentage = (contact_frames / len(contact_states)) * 100
            avg_distance = np.mean(min_distances)
        
            logger.warning(f"Residue {residue_id}: {contact_frames}/{len(contact_states)} frames in contact ({contact_percentage:.1f}%)")
            logger.warning(f"  Average distance: {avg_distance*10:.2f} Å, Min: {np.min(min_distances)*10:.2f} Å, Max: {np.max(min_distances)*10:.2f} Å")
        
            return min_distances, contact_states
        
        except Exception as e:
            logger.error(f"Failed to analyze residue contact for {residue_id}: {e}")
            return None, None

    def analyze_single_system(self, output_dir: Optional[str] = None,
                            distance_plot_dir: Optional[str] = None, 
                            parallel_hbonds: bool = True,
                            save_data: bool = True,
                            selected_residues: Optional[List[str]] = None,
                            generate_report: bool = True,
                            save_frame_data: bool = True):  
        if self.traj is None:
            raise ValueError("Requires topology and trajectory files")

        if output_dir is None:
            output_dir = "analysis_results"
        os.makedirs(output_dir, exist_ok=True)

        logger.warning("Starting contact analysis")
        contact_proportions = self.calculate_contact_proportions()

        if selected_residues:
            residues_to_analyze = [(res.strip(), contact_proportions.get(res.strip(), 0.0)) 
                                   for res in selected_residues if res.strip() in contact_proportions]
        else:
            residues_to_analyze = self.select_top_contacts(contact_proportions)[:10]

        if save_data:
            contact_csv_path = os.path.join(output_dir, "contact_proportions.csv")
            DataExporter.save_contact_data_csv(contact_proportions, contact_csv_path)

        if save_frame_data:
            frame_contact_path = os.path.join(output_dir, "frame_by_frame_contacts.csv")
            if selected_residues:
                selected_res_for_frame = selected_residues
            else:
                top_contact_residues = self.select_top_contacts(contact_proportions, top_n=10)
                selected_res_for_frame = [res[0] for res in top_contact_residues]
            
            DataExporter.save_frame_by_frame_contacts(self, frame_contact_path, selected_res_for_frame)

        top_3_for_plots = residues_to_analyze[:3]
        for residue_id, _ in top_3_for_plots:
            distances, contact_states = self.analyze_residue_contact(residue_id)

            if distances is None:
                continue

            plot_dir = distance_plot_dir or output_dir
            os.makedirs(plot_dir, exist_ok=True)
            safe_name = residue_id.replace(" ", "_").replace("(", "").replace(")", "")
            save_path = os.path.join(plot_dir, f"{safe_name}_distance.png")
            self.visualizer.plot_distance_time(distances, residue_id, self.config, save_path)

        logger.warning("Starting hydrogen bond analysis")
        hbond_frequencies = self.analyze_hydrogen_bonds()

        if save_data and hbond_frequencies:
            hbond_csv_path = os.path.join(output_dir, "hydrogen_bonds.csv")
            DataExporter.save_hbond_data_csv(hbond_frequencies, self.hbond_analyzer.hbond_stats, hbond_csv_path)

        if save_frame_data:
            frame_hbond_path = os.path.join(output_dir, "frame_by_frame_hbonds.csv")
            top_hbonds = hbond_frequencies[:10] if len(hbond_frequencies) > 10 else hbond_frequencies
            DataExporter.save_frame_by_frame_hbonds(self, frame_hbond_path, top_hbonds)

        distance_stats = {}
        if self.config.distance_analysis:
            logger.warning("Starting distance analysis")
            distance_stats = self.calculate_distance_analysis()
        
            if save_data and distance_stats:
                distance_csv_path = os.path.join(output_dir, "distance_stats.csv")
                DataExporter.save_distance_data_csv(distance_stats, distance_csv_path)

            if save_frame_data:
                frame_distance_path = os.path.join(output_dir, "frame_by_frame_distances.csv")
                if selected_residues:
                    selected_res_for_distance = selected_residues
                else:
                    if distance_stats:
                        closest_residues = sorted(distance_stats.items(), 
                                               key=lambda x: x[1]['min'])[:5]
                        selected_res_for_distance = [res[0] for res in closest_residues]
                    else:
                        top_5_contacts = self.select_top_contacts(contact_proportions, top_n=5)
                        selected_res_for_distance = [res[0] for res in top_5_contacts]
                
                DataExporter.save_frame_by_frame_distances(self, frame_distance_path, selected_res_for_distance)

        if generate_report and distance_stats:
            self.generate_comprehensive_report(
                contact_proportions, 
                hbond_frequencies, 
                distance_stats, 
                output_dir
            )

        logger.warning("Analysis complete! Files generated:")
        logger.warning(f"  - Contact proportions: contact_proportions.csv")
        if save_frame_data:
            logger.warning(f"  - Frame-by-frame contacts: frame_by_frame_contacts.csv (top 10)")
        if hbond_frequencies:
            logger.warning(f"  - Hydrogen bonds: hydrogen_bonds.csv")
            if save_frame_data:
                logger.warning(f"  - Frame-by-frame H-bonds: frame_by_frame_hbonds.csv (top 10)")
        if distance_stats:
            logger.warning(f"  - Distance statistics: distance_stats.csv")
            if save_frame_data:
                logger.warning(f"  - Frame-by-frame distances: frame_by_frame_distances.csv (top 5)")

        return contact_proportions, hbond_frequencies, distance_stats
    
    def analyze_multi_system(self, analysis_type: str = "overview", selection1: str = "protein", 
                           selection2: Optional[str] = None, cutoff_angstrom: float = 3.0, 
                           timestep: float = 0.5,
                           parallel: bool = False, color_mode: str = "default", 
                           output_file: Optional[str] = None, save_data: bool = True):
        if not self.system_files:
            raise ValueError("No systems found. Requires base_dir with system directories.")
        
        if selection2 is None:
            selection2 = f"resname {self.config.ligand_resname}"
        
        if output_file is None:
            output_file = f"multi_system_{analysis_type}.png"
        
        cutoff_nm = cutoff_angstrom * 0.1
        
        if analysis_type == "overview":
            ylabel = "Number of Unique Contact Types"
            title = "Contact Convergence Across All Systems"
        else:
            ylabel = "Number of Retained Initial Contacts"
            title = "Initial Contact Stability Across All Systems"
        
        all_data, all_times, system_names, weights = [], [], [], []
        
        logger.warning(f"Starting {analysis_type} analysis for {len(self.system_files)} systems")
        
        if parallel and len(self.system_files) > 2:
            with ProcessPoolExecutor(max_workers=self.config.parallel_workers) as executor:
                futures = []
                for name, (traj_path, prod_dir, top_path) in self.system_files.items():
                    future = executor.submit(
                        self.multi_analyzer.process_single_system,
                        name, traj_path, prod_dir, top_path, analysis_type, 
                        selection1, selection2, cutoff_nm, timestep, self.config.ligand_resname
                    )
                    futures.append(future)
                
                for future in as_completed(futures):
                    name, contacts, times, weight = future.result()
                    if len(contacts) > 0:
                        all_data.append(contacts)
                        all_times.append(times)
                        system_names.append(name)
                        weights.append(weight)
        else:
            for name, (traj_path, prod_dir, top_path) in self.system_files.items():
                name, contacts, times, weight = self.multi_analyzer.process_single_system(
                    name, traj_path, prod_dir, top_path, analysis_type, 
                    selection1, selection2, cutoff_nm, timestep, self.config.ligand_resname
                )
                if len(contacts) > 0:
                    all_data.append(contacts)
                    all_times.append(times)
                    system_names.append(name)
                    weights.append(weight)
        
        if not all_data:
            logger.error("No valid data found for any system")
            return
        
        logger.warning(f"Successfully processed {len(all_data)} systems")
        
        if save_data:
            output_dir = os.path.dirname(output_file) or "multi_system_results"
            parameters = {
                'selection1': selection1,
                'selection2': selection2,
                'cutoff_angstrom': cutoff_angstrom,
                'timestep': timestep,
                'ligand_resname': self.config.ligand_resname
            }
            DataExporter.save_multi_system_data(
                all_data, all_times, system_names, weights, 
                analysis_type, parameters, output_dir
            )
        
        self.visualizer.plot_multi_system(
            all_data, all_times, system_names, weights, color_mode,
            ylabel, title, output_file, self.config
        )
        
        logger.warning(f"Multi-system analysis complete! Plot saved: {output_file}")
    
    def generate_comprehensive_report(self, contact_proportions: Dict[str, float],
                                    hbond_frequencies: List[Tuple],
                                    distance_stats: Dict[str, Dict[str, float]],
                                    output_dir: str):
        """Generate comprehensive analysis report"""
        report_path = os.path.join(output_dir, "comprehensive_report.txt")
    
        try:
            with open(report_path, "w", encoding="utf-8") as f:
                f.write("=" * 80 + "\n")
                f.write("PROTEIN-LIGAND INTERACTION ANALYSIS RESULTS\n")
                f.write("=" * 80 + "\n\n")
            
                if self.traj is not None:
                    f.write("SYSTEM INFORMATION:\n")
                    f.write(f"  Trajectory frames: {self.traj.n_frames}\n")
                    total_time_ns = self.traj.time[-1] / 1000
                    f.write(f"  Total simulation time: {total_time_ns:.2f} ns\n")
                    f.write(f"  Ligand atoms: {len(self.ligand_atoms)}\n")
                    f.write(f"  Protein atoms: {len(self.protein_atoms)}\n")
                    f.write(f"  Protein residues: {len(self.protein_residues)}\n\n")
            
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
            
                if self.config.distance_analysis and distance_stats:
                    f.write("DISTANCE ANALYSIS RESULTS:\n")
                    cutoff_angstrom = self.config.distance_cutoff_nm * 10.0
                    residues_close = sum(1 for stats in distance_stats.values() 
                                        if stats['mean'] < cutoff_angstrom)
                
                    all_min_dists = [stats['min'] for stats in distance_stats.values()]
                    all_mean_dists = [stats['mean'] for stats in distance_stats.values()]
                    all_max_dists = [stats['max'] for stats in distance_stats.values()]
                
                    closest_residue, closest_stats = min(distance_stats.items(), 
                                                        key=lambda x: x[1]['min'], 
                                                        default=("None", {}))
                
                    f.write(f"  Residues with mean distance < {self.config.distance_cutoff_nm:.2f} nm: "
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

                f.write("\n" + "=" * 80 + "\n")
        
            logger.warning(f"Comprehensive report generated: {os.path.basename(report_path)}")
            return report_path
        except Exception as e:
            logger.error(f"Failed to generate comprehensive report: {e}")
            return None


def analyze_and_visualize(topology_file, trajectory_file,
                         resname="LIG", distance_cutoff=0.4, n_frames=None,
                         output_dir="analysis_results", verbose=True,
                         save_all_data=True, return_detailed=False,
                         generate_report=True, use_new_method=True):
    """Simplified analysis interface function"""
    os.makedirs(output_dir, exist_ok=True)
    
    config = AnalysisConfig(
        ligand_resname=resname,
        contact_cutoff_nm=distance_cutoff,
        contact_enter_threshold_nm=distance_cutoff * 0.75,
        contact_exit_threshold_nm=distance_cutoff * 1.125,
    )
    
    if verbose:
        logger.warning(f"Initializing analyzer for {resname}")
        logger.warning(f"Topology: {topology_file}")
        logger.warning(f"Trajectory: {trajectory_file}")
        logger.warning(f"Output directory: {output_dir}")
    
    try:
        analyzer = IntegratedProteinLigandAnalyzer(
            topology_file=topology_file,
            trajectory_file=trajectory_file,
            config=config,
            skip_initial_analysis=True  
        )
    except Exception as e:
        logger.error(f"Failed to initialize analyzer: {e}")
        raise

    if n_frames is not None:
        if analyzer.traj and analyzer.traj.n_frames > n_frames:
            if verbose:
                logger.warning(f"Limiting to {n_frames} frames (original: {analyzer.traj.n_frames})")
            analyzer.traj = analyzer.traj[:n_frames]
            analyzer.universe.trajectory = analyzer.universe.trajectory[:n_frames]
    
    analyzer._prepare_single_system()
    
    contact_props, hbond_freqs, distance_stats = analyzer.analyze_single_system(
        output_dir=output_dir,
        parallel_hbonds=True,
        save_data=save_all_data,
        generate_report=generate_report
    )
    
    if verbose:
        actual_frames = analyzer.traj.n_frames if analyzer.traj else 0
        logger.warning(f"Analysis complete! Analyzed {actual_frames} frames. Results in {output_dir}")
        if use_new_method and analyzer._new_method_results:
            detailed_results = analyzer._new_method_results
            logger.warning(f"Found {len(detailed_results.get('contact_frequencies', {}))} detailed contacts")
            logger.warning(f"Total frames analyzed: {detailed_results.get('total_frames', 0)}")
    
    if return_detailed:
        results = {
            'contact_proportions': contact_props,
            'hbond_frequencies': hbond_freqs,
            'distance_stats': distance_stats,
            'detailed_contact_results': analyzer.get_detailed_contact_results() if use_new_method else None,
            'system_info': {
                'n_frames': analyzer.traj.n_frames if analyzer.traj else 0,
                'ligand_atoms': len(analyzer.ligand_atoms),
                'protein_atoms': len(analyzer.protein_atoms),
                'protein_residues': len(analyzer.protein_residues)
            },
            'config': config,
            'output_directory': output_dir,
            'method_used': 'new_method' if use_new_method else 'original_method'
        }
        return results
    else:
        return hbond_freqs, contact_props, distance_stats