#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Main trajectory analysis class for PRISM
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

from .config import AnalysisConfig
from .trajectory import TrajectoryManager
from .contact import ContactAnalyzer
from .hbond import HydrogenBondAnalyzer
from .distance import DistanceAnalyzer
# from .visualization import Visualizer
from .io import DataExporter, ReportGenerator
import numpy as np

logger = logging.getLogger(__name__)


class TrajAnalysis:
    """
    Main trajectory analysis class for PRISM
    
    This class provides a unified interface for analyzing MD trajectories
    of protein-ligand systems, including contact analysis, hydrogen bonds,
    and distance calculations.
    
    Examples
    --------
    >>> import prism as pm
    >>> traj = pm.TrajAnalysis("topology.gro", "trajectory.xtc")
    >>> contacts = traj.analyze_contacts()
    >>> hbonds = traj.analyze_hydrogen_bonds()
    >>> traj.generate_report("analysis_output")
    """
    
    def __init__(self, topology_file: str, trajectory_file: str, 
                 ligand_resname: str = "LIG", config: Optional[AnalysisConfig] = None):
        """
        Initialize trajectory analysis
        
        Parameters
        ----------
        topology_file : str
            Path to topology file (e.g., .gro, .pdb)
        trajectory_file : str
            Path to trajectory file (e.g., .xtc, .dcd)
        ligand_resname : str, optional
            Residue name of the ligand (default: "LIG")
        config : AnalysisConfig, optional
            Configuration object for analysis parameters
        """
        self.topology_file = os.path.abspath(topology_file)
        self.trajectory_file = os.path.abspath(trajectory_file)
        
        # Initialize configuration
        if config is None:
            config = AnalysisConfig(ligand_resname=ligand_resname)
        self.config = config
        
        # Initialize trajectory manager
        self.trajectory_manager = TrajectoryManager(
            self.topology_file, 
            self.trajectory_file, 
            self.config
        )
        
        if not self.trajectory_manager.validate_system():
            raise ValueError("System validation failed")
        
        # Initialize analyzers
        self.contact_analyzer = ContactAnalyzer(self.trajectory_manager, self.config)
        self.hbond_analyzer = HydrogenBondAnalyzer(self.trajectory_manager, self.config)
        self.distance_analyzer = DistanceAnalyzer(self.trajectory_manager, self.config)
        self.visualizer = Visualizer(self.trajectory_manager, self.config)
        
        # Results storage
        self.results = {
            'contacts': None,
            'hydrogen_bonds': None,
            'distances': None
        }
        
        logger.info(f"Initialized TrajAnalysis with {self.n_frames} frames")
    
    @property
    def n_frames(self) -> int:
        """Number of frames in trajectory"""
        return self.trajectory_manager.traj.n_frames
    
    @property
    def n_atoms(self) -> int:
        """Number of atoms in system"""
        return self.trajectory_manager.traj.n_atoms
    
    @property
    def ligand_atoms(self) -> List[int]:
        """Indices of ligand atoms"""
        return self.trajectory_manager.ligand_atoms
    
    @property
    def protein_atoms(self) -> List[int]:
        """Indices of protein atoms"""
        return self.trajectory_manager.protein_atoms
    
    @property
    def traj(self):
        """Direct access to mdtraj trajectory object"""
        return self.trajectory_manager.traj
    
    def analyze_contacts(self, save_data: bool = True, 
                        output_dir: Optional[str] = None) -> Dict[str, float]:
        """
        Analyze protein-ligand contacts
        
        Parameters
        ----------
        save_data : bool, optional
            Whether to save results to file
        output_dir : str, optional
            Directory to save results
            
        Returns
        -------
        Dict[str, float]
            Dictionary of residue contact proportions
        """
        logger.info("Starting contact analysis...")
        
        contact_proportions = self.contact_analyzer.calculate_contact_proportions()
        self.results['contacts'] = contact_proportions
        
        if save_data and output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, "contact_proportions.csv")
            exporter = DataExporter()
            exporter.save_contact_data_csv(contact_proportions, output_path)
        
        logger.info(f"Contact analysis complete. Found {len(contact_proportions)} contacts.")
        return contact_proportions
    
    def analyze_hydrogen_bonds(self, parallel: bool = True, save_data: bool = True,
                              output_dir: Optional[str] = None) -> List[Tuple]:
        """
        Analyze hydrogen bonds between protein and ligand
        
        Parameters
        ----------
        parallel : bool, optional
            Use parallel processing
        save_data : bool, optional
            Whether to save results to file
        output_dir : str, optional
            Directory to save results
            
        Returns
        -------
        List[Tuple]
            List of hydrogen bond data (name, frequency, distance, angle)
        """
        logger.info("Starting hydrogen bond analysis...")
        
        hbond_frequencies, hbond_stats = self.hbond_analyzer.analyze_hydrogen_bonds(parallel)
        self.results['hydrogen_bonds'] = (hbond_frequencies, hbond_stats)
        
        if save_data and output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, "hydrogen_bonds.csv")
            exporter = DataExporter()
            exporter.save_hbond_data_csv(hbond_frequencies, output_path, hbond_stats)
        
        logger.info(f"H-bond analysis complete. Found {len(hbond_frequencies)} H-bonds.")
        return hbond_frequencies
    
    def analyze_distances(self, save_data: bool = True,
                         output_dir: Optional[str] = None) -> Dict[str, Dict[str, float]]:
        """
        Analyze distances between ligand and protein residues
        
        Parameters
        ----------
        save_data : bool, optional
            Whether to save results to file
        output_dir : str, optional
            Directory to save results
            
        Returns
        -------
        Dict[str, Dict[str, float]]
            Distance statistics for each residue
        """
        logger.info("Starting distance analysis...")
        
        distance_stats = self.distance_analyzer.calculate_distance_analysis()
        self.results['distances'] = distance_stats
        
        if save_data and output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, "distance_stats.csv")
            exporter = DataExporter()
            exporter.save_distance_data_csv(distance_stats, output_path)
        
        logger.info(f"Distance analysis complete for {len(distance_stats)} residues.")
        return distance_stats
    
    def analyze_all(self, output_dir: str = "analysis_results", 
                   parallel: bool = True) -> Dict[str, Any]:
        """
        Run all analyses
        
        Parameters
        ----------
        output_dir : str, optional
            Output directory for results
        parallel : bool, optional
            Use parallel processing for H-bonds
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing all analysis results
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Run all analyses
        contacts = self.analyze_contacts(save_data=True, output_dir=output_dir)
        hbonds = self.analyze_hydrogen_bonds(parallel=parallel, save_data=True, 
                                            output_dir=output_dir)
        distances = self.analyze_distances(save_data=True, output_dir=output_dir)
        
        # Generate comprehensive report
        self.generate_report(output_dir)
        
        return {
            'contacts': contacts,
            'hydrogen_bonds': hbonds,
            'distances': distances
        }
    
    def analyze_residue_contact(self, residue_id: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Analyze contact timeline for a specific residue
        
        Parameters
        ----------
        residue_id : str
            Residue identifier (e.g., "GLU 123")
            
        Returns
        -------
        Tuple[Optional[np.ndarray], Optional[np.ndarray]]
            Distances and contact states arrays
        """
        return self.contact_analyzer.analyze_residue_contact(residue_id)
    
    def plot_contact_timeline(self, residue_id: str, output_path: Optional[str] = None):
        """
        Plot contact timeline for a specific residue
        
        Parameters
        ----------
        residue_id : str
            Residue identifier (e.g., "GLU 123")
        output_path : str, optional
            Path to save plot
        """
        distances, contact_states = self.contact_analyzer.analyze_residue_contact(residue_id)
        
        if distances is not None:
            self.visualizer.plot_distance_timeline(
                distances, contact_states, residue_id, output_path
            )
    
    def visualize_contacts(self, ligand_sdf: str, output_path: Optional[str] = None,
                          freq_threshold: float = 0.2, top_n: int = 10):
        """
        Visualize top protein-ligand contacts on 2D ligand structure
        
        Parameters
        ----------
        ligand_sdf : str
            Path to ligand SDF file
        output_path : str, optional
            Path to save visualization
        freq_threshold : float, optional
            Minimum contact frequency threshold
        top_n : int, optional
            Number of top contacts to show
        """
        self.visualizer.visualize_ligand_contacts(
            ligand_sdf, output_path, freq_threshold, top_n
        )
    
    def generate_report(self, output_dir: str = "analysis_results"):
        """
        Generate comprehensive analysis report
        
        Parameters
        ----------
        output_dir : str
            Directory to save report
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Ensure all analyses have been run
        if self.results['contacts'] is None:
            self.analyze_contacts(save_data=True, output_dir=output_dir)
        if self.results['hydrogen_bonds'] is None:
            self.analyze_hydrogen_bonds(save_data=True, output_dir=output_dir)
        if self.results['distances'] is None:
            self.analyze_distances(save_data=True, output_dir=output_dir)
        
        # Generate report
        report_generator = ReportGenerator()
        report_path = report_generator.generate_comprehensive_report(
            self.results['contacts'],
            self.results['hydrogen_bonds'][0] if self.results['hydrogen_bonds'] else [],
            self.results['distances'],
            output_dir,
            self.config,
            self.trajectory_manager,
            self.results['hydrogen_bonds'][1] if self.results['hydrogen_bonds'] else {}
        )
        
        logger.info(f"Comprehensive report generated: {report_path}")
        return report_path
    
    def get_system_info(self) -> Dict[str, Any]:
        """Get system information"""
        return self.trajectory_manager.get_system_info()
    
    def get_top_contacts(self, n: int = 10) -> List[Tuple[str, float]]:
        """Get top N contacting residues"""
        if self.results['contacts'] is None:
            self.analyze_contacts(save_data=False)
        
        return self.distance_analyzer.select_top_contacts(self.results['contacts'], n)
    
    def generate_html_visualization(self, ligand_file, output_file="contact_analysis.html"):
        """
        Generate interactive HTML visualization
        
        Parameters
        ----------
        ligand_file : str
            Path to ligand structure file (.sdf, .mol, .mol2)
        output_file : str
            Output HTML file path
            
        Returns
        -------
        str
            Path to generated HTML file
        """
        from .visualization import HTMLGenerator
        
        print("\nGenerating interactive HTML visualization...")
        
        # Use the loaded trajectory
        generator = HTMLGenerator(
            self.trajectory_file,
            self.topology_file,
            ligand_file
        )
        
        # Use existing trajectory if already loaded
        if self.traj is not None:
            generator.traj = self.traj
        
        output_path = generator.generate(output_file)
        
        return output_path