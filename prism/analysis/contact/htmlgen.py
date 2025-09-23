#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced Fixed Trajectory Analysis to Interactive HTML Generator
Integration with PRISM analysis module
"""

import argparse
import sys
from pathlib import Path

try:
    import mdtraj as md
    MDTRAJ_AVAILABLE = True
except ImportError:
    MDTRAJ_AVAILABLE = False
    md = None

# Import from local modules
from .contact_analyzer import FastContactAnalyzer
from .html_builder import HTMLBuilder
from .data_processor import DataProcessor
from .visualization_generator import VisualizationGenerator

class HTMLGenerator:
    """Generate interactive HTML visualization from trajectory analysis"""
    
    def __init__(self, trajectory_file, topology_file, ligand_file):
        """
        Initialize HTML generator
        
        Parameters
        ----------
        trajectory_file : str
            Path to trajectory file (.xtc, .dcd, etc.)
        topology_file : str
            Path to topology file (.pdb, .gro, etc.)
        ligand_file : str
            Path to ligand file (.sdf, .mol, .mol2)
        """
        self.trajectory_file = trajectory_file
        self.topology_file = topology_file
        self.ligand_file = ligand_file
        self.traj = None
        self.contact_results = None
        self.ligand_data = None
        self.contacts = None
        self.stats = None
        
    def analyze(self):
        """
        Run the complete analysis
        
        Returns
        -------
        dict
            Analysis results including contacts, ligand data, and statistics
        """
        print("=== Enhanced Contact Analysis ===")
        
        # Load trajectory
        print(f"Loading trajectory: {Path(self.trajectory_file).name}")
        self.traj = md.load(self.trajectory_file, top=self.topology_file)
        print(f"Loaded {self.traj.n_frames} frames, {self.traj.n_atoms} atoms")
        
        # Load and process ligand
        data_processor = DataProcessor()
        ligand_mol = data_processor.load_ligand_structure(self.ligand_file)
        
        # Analyze contacts
        print("\nAnalyzing contacts...")
        analyzer = FastContactAnalyzer(self.traj)
        self.contact_results = analyzer.calculate_contact_proportions()
        
        # Generate visualization data
        vis_generator = VisualizationGenerator()
        self.ligand_data = vis_generator.generate_ligand_data(self.contact_results, ligand_mol)
        self.contacts = vis_generator.generate_contact_data(self.contact_results, self.ligand_data)
        
        # Calculate statistics
        self.stats = data_processor.calculate_statistics(self.contacts)
        
        print(f"Found {self.stats['total_contacts']} significant contacts")
        
        return {
            'contact_results': self.contact_results,
            'ligand_data': self.ligand_data,
            'contacts': self.contacts,
            'stats': self.stats,
            'traj': self.traj
        }
    
    def generate(self, output_file="contact_analysis.html"):
        """
        Generate HTML file
        
        Parameters
        ----------
        output_file : str
            Path for output HTML file
            
        Returns
        -------
        str
            Path to generated HTML file
        """
        if self.contact_results is None:
            self.analyze()
        
        print("\nGenerating interactive HTML...")
        html_builder = HTMLBuilder()
        html_content = html_builder.generate_html(
            trajectory_file=Path(self.trajectory_file).name,
            topology_file=Path(self.topology_file).name,
            ligand_file=Path(self.ligand_file).name,
            ligand_name=self.contact_results['ligand_residue'].name if self.contact_results['ligand_residue'] else 'UNK',
            total_frames=self.contact_results['total_frames'],
            ligand_data=self.ligand_data,
            contacts=self.contacts,
            stats=self.stats
        )
        
        # Write output
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"Interactive HTML generated: {output_file}")
        print(f"Open the file in a web browser to view and adjust the layout")
        
        return output_file

def generate_html(trajectory, topology, ligand, output="contact_analysis.html"):
    """
    Convenience function to generate HTML visualization
    
    Parameters
    ----------
    trajectory : str
        Path to trajectory file
    topology : str
        Path to topology file
    ligand : str
        Path to ligand file
    output : str
        Output HTML file path
        
    Returns
    -------
    str
        Path to generated HTML file
    """
    generator = HTMLGenerator(trajectory, topology, ligand)
    return generator.generate(output)

def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(description="Enhanced trajectory analysis to interactive HTML")
    parser.add_argument("trajectory", help="Trajectory file (.xtc, .dcd, etc.)")
    parser.add_argument("topology", help="Topology file (.pdb, .gro, etc.)")
    parser.add_argument("ligand", help="Ligand file (.sdf, .mol, .mol2)")
    parser.add_argument("-o", "--output", default="contact_analysis.html", help="Output HTML file")
    
    args = parser.parse_args()
    
    try:
        generate_html(args.trajectory, args.topology, args.ligand, args.output)
        print(f"\n=== Complete! ===")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()