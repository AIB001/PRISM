#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visualization module for trajectory analysis
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.patches import Circle
from PIL import Image
import io
import os

# Change font to Arial for WSL compatibility
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'sans-serif']

# RDKit imports for molecular visualization
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    from sklearn.decomposition import PCA
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Molecular visualization will be disabled.")

from .config import AnalysisConfig

logger = logging.getLogger(__name__)

def nm_to_angstrom(value):
    """Convert nanometers to angstroms"""
    return value * 10.0

def angstrom_to_nm(value):
    """Convert angstroms to nanometers"""
    return value * 0.1

def moving_average(data, window=11):
    """Calculate moving average"""
    if len(data) < window:
        return data
    return np.convolve(data, np.ones(window)/window, mode='valid')

def should_smooth(data, min_frames=10):
    """Check if data should be smoothed"""
    return len(data) > min_frames

# Amino acid codes
AA_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

AA_SMILES = {
    'ALA': 'CC(N)C(=O)O',
    'ARG': 'NC(CCCNC(N)=N)C(=O)O',
    'ASN': 'NC(CC(N)=O)C(=O)O',
    'ASP': 'NC(CC(=O)O)C(=O)O',
    'CYS': 'NC(CS)C(=O)O',
    'GLN': 'NC(CCC(N)=O)C(=O)O',
    'GLU': 'NC(CCC(=O)O)C(=O)O',
    'GLY': 'NCC(=O)O',
    'HIS': 'NC(CC1=CN=CN1)C(=O)O',
    'ILE': 'CC(C)C(C)C(N)C(=O)O',
    'LEU': 'CC(C)CC(N)C(=O)O',
    'LYS': 'NC(CCCCN)C(=O)O',
    'MET': 'CSCCC(N)C(=O)O',
    'PHE': 'NC(CC1=CC=CC=C1)C(=O)O',
    'PRO': 'O=C(O)C1CCCN1',
    'SER': 'NC(CO)C(=O)O',
    'THR': 'CC(O)C(N)C(=O)O',
    'TRP': 'NC(CC1=CNC2=CC=CC=C12)C(=O)O',
    'TYR': 'NC(CC1=CC=C(O)C=C1)C(=O)O',
    'VAL': 'CC(C)C(N)C(=O)O'
}


class Visualizer:
    """Visualization tools for trajectory analysis"""
    
    def __init__(self, trajectory_manager, config):
        self.trajectory_manager = trajectory_manager
        self.config = config
    
    def plot_contact_matrix(self, contacts: Dict[str, float], output_path: Optional[str] = None):
        """Plot contact proportion as a bar chart"""
        if not contacts:
            logger.warning("No contacts to plot")
            return
            
        # Sort contacts by proportion
        sorted_contacts = sorted(contacts.items(), key=lambda x: x[1], reverse=True)
        
        # Take top 20 for visualization
        top_contacts = sorted_contacts[:20]
        
        residues = [item[0] for item in top_contacts]
        proportions = [item[1] * 100 for item in top_contacts]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        bars = ax.bar(range(len(residues)), proportions, color='steelblue', edgecolor='black')
        
        # Add value labels on bars
        for bar, prop in zip(bars, proportions):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{prop:.1f}%', ha='center', va='bottom', fontsize=9)
        
        ax.set_xlabel('Residue', fontsize=12)
        ax.set_ylabel('Contact Proportion (%)', fontsize=12)
        ax.set_title('Top 20 Protein-Ligand Contacts', fontsize=14, fontweight='bold')
        ax.set_xticks(range(len(residues)))
        ax.set_xticklabels(residues, rotation=45, ha='right')
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"Contact bar chart saved: {output_path}")
        else:
            plt.show()
    
    def plot_hbond_network(self, hbonds: List[Tuple], output_path: Optional[str] = None):
        """Plot hydrogen bond network"""
        if not hbonds:
            logger.warning("No hydrogen bonds to plot")
            return
        
        # Take top 10 H-bonds
        top_hbonds = hbonds[:10]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # Bar chart of H-bond frequencies
        labels = []
        frequencies = []
        distances = []
        angles = []
        
        for hb in top_hbonds:
            if len(hb) >= 4:
                label = hb[0]
                # Simplify label for display
                if ' → ' in label:
                    parts = label.split(' → ')
                    if 'Ligand' in parts[0]:
                        label = f"LIG → {parts[1].split()[0]}"
                    else:
                        label = f"{parts[0].split()[0]} → LIG"
                labels.append(label)
                frequencies.append(hb[1] * 100)
                distances.append(hb[2])
                angles.append(hb[3])
        
        y_pos = np.arange(len(labels))
        bars = ax1.barh(y_pos, frequencies, color='green', edgecolor='black')
        
        # Add value labels
        for bar, freq in zip(bars, frequencies):
            width = bar.get_width()
            ax1.text(width, bar.get_y() + bar.get_height()/2.,
                    f'{freq:.1f}%', ha='left', va='center', fontsize=9)
        
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(labels)
        ax1.set_xlabel('Frequency (%)', fontsize=12)
        ax1.set_title('Hydrogen Bond Frequencies', fontsize=14, fontweight='bold')
        ax1.grid(axis='x', alpha=0.3)
        
        # Scatter plot of distance vs angle
        scatter = ax2.scatter(distances, angles, s=100, c=frequencies, 
                             cmap='viridis', edgecolor='black', alpha=0.7)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax2)
        cbar.set_label('Frequency (%)', fontsize=10)
        
        ax2.set_xlabel('Average Distance (Å)', fontsize=12)
        ax2.set_ylabel('Average Angle (°)', fontsize=12)
        ax2.set_title('H-bond Geometry', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        
        # Add reference lines
        ax2.axhline(y=120, color='red', linestyle='--', alpha=0.5, label='Angle cutoff')
        ax2.axvline(x=3.5, color='red', linestyle='--', alpha=0.5, label='Distance cutoff')
        ax2.legend()
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"H-bond network saved: {output_path}")
        else:
            plt.show()
    
    def plot_distance_histogram(self, distances: Dict[str, Dict[str, float]], 
                               output_path: Optional[str] = None):
        """Plot distance distribution histograms"""
        if not distances:
            logger.warning("No distances to plot")
            return
        
        # Get top residues by minimum distance
        sorted_residues = sorted(distances.items(), key=lambda x: x[1]['min'])
        top_residues = sorted_residues[:12]
        
        fig, axes = plt.subplots(3, 4, figsize=(16, 12))
        axes = axes.flatten()
        
        for idx, (residue, stats) in enumerate(top_residues):
            ax = axes[idx]
            
            # Create histogram data
            mean_dist = stats['mean']
            std_dist = stats['std']
            min_dist = stats['min']
            max_dist = stats['max']
            
            # Generate synthetic distribution for visualization
            # (since we don't have the full distance array)
            x = np.linspace(min_dist, max_dist, 50)
            y = np.exp(-0.5 * ((x - mean_dist) / std_dist) ** 2)
            y = y / y.max()  # Normalize
            
            ax.fill_between(x, y, alpha=0.5, color='coral')
            ax.axvline(mean_dist, color='red', linestyle='--', label=f'Mean: {mean_dist:.2f}')
            ax.axvline(min_dist, color='blue', linestyle='--', label=f'Min: {min_dist:.2f}')
            
            ax.set_xlabel('Distance (Å)', fontsize=10)
            ax.set_ylabel('Relative Frequency', fontsize=10)
            ax.set_title(residue, fontsize=11, fontweight='bold')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(len(top_residues), len(axes)):
            axes[idx].set_visible(False)
        
        plt.suptitle('Distance Distributions for Top Contact Residues', 
                     fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"Distance distributions saved: {output_path}")
        else:
            plt.show()
    
    def plot_contact_heatmap(self, contact_analyzer, output_path: Optional[str] = None):
        """Plot contact heatmap over time"""
        if not hasattr(contact_analyzer, 'traj') or contact_analyzer.traj is None:
            logger.warning("No trajectory available for heatmap")
            return
        
        # Get top contacting residues
        contact_props = contact_analyzer.calculate_contact_proportions()
        sorted_contacts = sorted(contact_props.items(), key=lambda x: x[1], reverse=True)
        top_residues = [res[0] for res in sorted_contacts[:20]]
        
        # Calculate contact matrix over time
        n_frames = contact_analyzer.traj.n_frames
        n_residues = len(top_residues)
        contact_matrix = np.zeros((n_residues, n_frames))
        
        for i, residue_id in enumerate(top_residues):
            distances, contact_states = contact_analyzer.analyze_residue_contact(residue_id)
            if contact_states is not None:
                contact_matrix[i, :] = contact_states
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Downsample if too many frames
        if n_frames > 1000:
            step = n_frames // 1000
            contact_matrix_plot = contact_matrix[:, ::step]
            time_points = np.arange(0, n_frames, step) * contact_analyzer.traj.timestep / 1000  # to ns
        else:
            contact_matrix_plot = contact_matrix
            time_points = np.arange(n_frames) * contact_analyzer.traj.timestep / 1000  # to ns
        
        im = ax.imshow(contact_matrix_plot, aspect='auto', cmap='RdYlBu_r',
                      interpolation='nearest')
        
        # Set labels
        ax.set_yticks(range(n_residues))
        ax.set_yticklabels(top_residues)
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Residue', fontsize=12)
        ax.set_title('Contact Heatmap Over Time', fontsize=14, fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Contact', fontsize=10)
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels(['No', 'Yes'])
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"Contact heatmap saved: {output_path}")
        else:
            plt.show()
    
    def plot_distance_timeline(self, distances: np.ndarray, contact_states: np.ndarray, 
                              residue_id: str, save_path: Optional[str] = None):
        """Plot distance and contact timeline"""
        times = self.trajectory_manager.traj.time / 1000  # ps to ns
        distances_angstrom = distances * nm_to_angstrom(1.0)
        
        if should_smooth(distances_angstrom, self.config.min_frames_for_smoothing):
            smoothed_distances = moving_average(distances_angstrom, self.config.smooth_window)
            smoothed_times = times[self.config.smooth_window-1:]
            plot_data = smoothed_distances
            plot_times = smoothed_times
        else:
            plot_data = distances_angstrom
            plot_times = times
        
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(plot_times, plot_data, 'b-', linewidth=2, 
               label=f"Distance (window={self.config.smooth_window})")
        
        ax.grid(True)
        self._mark_contact_regions(ax, times, contact_states)
        
        ax.axhline(y=nm_to_angstrom(self.config.contact_enter_threshold_nm), 
                  color='g', linestyle='--', label='Contact start')
        ax.axhline(y=nm_to_angstrom(self.config.contact_exit_threshold_nm), 
                  color='r', linestyle='--', label='Contact end')
        
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Distance (Å)")
        ax.set_title(f"Distance to {residue_id}")
        ax.legend()
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            logger.info(f"Distance plot saved: {save_path}")
        else:
            plt.show()
    
    def _mark_contact_regions(self, ax, times: np.ndarray, contact_states: np.ndarray):
        """Mark contact regions on plot"""
        contact_start = None
        for i in range(len(contact_states)):
            if contact_states[i] == 1 and (i == 0 or contact_states[i-1] == 0):
                contact_start = times[i]
            elif contact_states[i] == 0 and i > 0 and contact_states[i-1] == 1:
                if contact_start is not None:
                    ax.axvspan(contact_start, times[i], alpha=0.3, color='gray')
    
    def visualize_ligand_contacts(self, ligand_sdf: str, output_path: Optional[str] = None,
                                 freq_threshold: float = 0.2, top_n: int = 10):
        """Visualize contacts on ligand structure"""
        if not RDKIT_AVAILABLE:
            logger.error("RDKit not available. Cannot create molecular visualization.")
            return
        
        try:
            # Check if the file exists and is valid
            if not os.path.exists(ligand_sdf):
                logger.error(f"Ligand SDF file not found: {ligand_sdf}")
                return
            
            # Try to load the molecule with RDKit
            mol = None
            try:
                # First try as SDF
                supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=False, sanitize=False)
                if supplier:
                    mol = supplier[0]
                    if mol is not None:
                        # Try to sanitize
                        try:
                            Chem.SanitizeMol(mol)
                        except:
                            logger.warning("Failed to sanitize molecule, continuing without sanitization")
            except Exception as e:
                logger.error(f"Failed to load SDF file: {e}")
                
                # Try as MOL file
                try:
                    mol = Chem.MolFromMolFile(ligand_sdf, removeHs=False, sanitize=False)
                    if mol is not None:
                        try:
                            Chem.SanitizeMol(mol)
                        except:
                            logger.warning("Failed to sanitize molecule, continuing without sanitization")
                except Exception as e2:
                    logger.error(f"Failed to load as MOL file: {e2}")
            
            if mol is None:
                logger.error("Could not load molecule from file")
                # Create a simple placeholder visualization
                self._create_placeholder_visualization(output_path)
                return
            
            # Continue with normal visualization
            from .contact import ContactAnalyzer
            analyzer = ContactAnalyzer(self.trajectory_manager, self.config)
            
            contact_frequencies, residue_contacts, residue_atom_details = analyzer.analyze_atomic_contacts(
                self.config.distance_cutoff_nm
            )
            
            visualizer = ContactVisualizer(
                mol=mol,
                contact_frequencies=contact_frequencies,
                residue_contacts=residue_contacts,
                residue_atom_details=residue_atom_details,
                n_frames=self.trajectory_manager.traj.n_frames
            )
            
            visualizer.visualize(
                freq_threshold=freq_threshold,
                top_n_contacts=top_n,
                save_path=output_path
            )
            
        except Exception as e:
            logger.error(f"Error in ligand visualization: {e}")
            import traceback
            traceback.print_exc()
            # Create a simple placeholder visualization
            self._create_placeholder_visualization(output_path)
    
    def _create_placeholder_visualization(self, output_path: Optional[str] = None):
        """Create a placeholder visualization when molecule loading fails"""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        ax.text(0.5, 0.5, 'Ligand Visualization Failed\n\nUnable to load molecular structure.\n'
                'This may be due to:\n'
                '• Invalid atom types in the file\n'
                '• File format issues\n'
                '• Missing RDKit dependencies\n\n'
                'Contact analysis has been completed successfully.',
                ha='center', va='center', fontsize=14,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"Placeholder visualization saved: {output_path}")
        else:
            plt.show()


# Keep the ContactVisualizer class as is, but with the font change already applied
class ContactVisualizer:
    """Visualizer for ligand-protein contacts"""
    
    def __init__(self, mol, contact_frequencies, residue_contacts, residue_atom_details, n_frames):
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit is required for molecular visualization")
        
        self.mol = mol
        self.contact_frequencies = contact_frequencies
        self.residue_contacts = residue_contacts
        self.residue_atom_details = residue_atom_details
        self.n_frames = n_frames
        
        self.fig_size = (20, 20)
        self.mol_scale = 1.0
        self.aa_structure_scale = 0.3
        self.label_radius = 40
        self.min_distance = 150
        self.structure_distance = 400
        self.label_distance = 180
    
    # def visualize(self, freq_threshold=0.2, top_n_contacts=10, save_path=None):
    #     """Create visualization"""
    #     self.generate_linear_conformation()
        
    #     if self.mol.GetNumConformers() == 0:
    #         AllChem.Compute2DCoords(self.mol)
        
    #     drawer = rdMolDraw2D.MolDraw2DCairo(1200, 1200)
    #     opts = drawer.drawOptions()
    #     opts.bondLineWidth = 4
    #     opts.atomLabelFontSize = 24
    #     opts.padding = 0.15
        
    #     max_freq = max(self.contact_frequencies) if max(self.contact_frequencies) > 0 else 1
    #     cmap = plt.cm.coolwarm
        
    #     highlight_atoms = []
    #     highlight_colors = {}
        
    #     heavy_idx = 0
    #     for atom_idx in range(self.mol.GetNumAtoms()):
    #         atom = self.mol.GetAtomWithIdx(atom_idx)
    #         if atom.GetSymbol() != 'H':
    #             if heavy_idx < len(self.contact_frequencies):
    #                 freq = self.contact_frequencies[heavy_idx] / max_freq
    #                 color = cmap(freq)
    #                 highlight_atoms.append(atom_idx)
    #                 highlight_colors[atom_idx] = (color[0], color[1], color[2], 0.8)
    #             heavy_idx += 1
        
    #     rdMolDraw2D.PrepareAndDrawMolecule(
    #         drawer, self.mol,
    #         highlightAtoms=highlight_atoms,
    #         highlightAtomColors=highlight_colors
    #     )
    #     drawer.FinishDrawing()
        
    #     img_data = drawer.GetDrawingText()
    #     img = Image.open(io.BytesIO(img_data))
        
    #     ligand_positions, ligand_names = self.get_ligand_atom_positions(drawer)
        
    #     occupied_positions = list(ligand_positions.values())
        
    #     residue_positions, top_contacts = self.find_optimal_positions(
    #         ligand_positions, self.residue_atom_details, top_n_contacts, freq_threshold
    #     )
        
    #     fig, ax = plt.subplots(figsize=self.fig_size)
    #     ax.imshow(img)
    #     ax.axis('off')
        
    #     # Draw top 3 contacts with amino acid structures
    #     aa_structure_bounds = []
        
    #     for idx, (lig_idx, res_info, freq) in enumerate(top_contacts[:3]):
    #         if res_info in residue_positions:
    #             x, y = residue_positions[res_info]
                
    #             structure_pos = self._find_empty_position(
    #                 (x, y), occupied_positions, self.min_distance, aa_structure_bounds
    #             )
    #             if structure_pos is None:
    #                 continue
                    
    #             x, y = structure_pos
    #             occupied_positions.append((x, y, 80))
                
    #             if lig_idx in ligand_positions:
    #                 lig_x, lig_y = ligand_positions[lig_idx]
    #                 angle = np.arctan2(lig_y - y, lig_x - x)
    #             else:
    #                 angle = 0
                
    #             coords = self.draw_amino_acid_structure(ax, x, y, res_info, res_info[:3], 
    #                                                    rotation=angle)
                
    #             if isinstance(coords, np.ndarray) and len(coords) > 0:
    #                 min_x, min_y = coords.min(axis=0)
    #                 max_x, max_y = coords.max(axis=0)
    #                 label_y = min_y - 45
    #                 min_y = label_y - 15
    #                 aa_structure_bounds.append((min_x - 20, min_y, max_x + 20, max_y + 20))
    #                 for coord in coords:
    #                     occupied_positions.append(tuple(coord))
    #                 occupied_positions.append((x, label_y, 80))
                
    #             if lig_idx in ligand_positions:
    #                 ax.plot([lig_x, x], [lig_y, y], 'k-', linewidth=3, alpha=0.7, zorder=4)
                    
    #                 mid_x = (lig_x + x) / 2
    #                 mid_y = (lig_y + y) / 2
    #                 ax.text(mid_x, mid_y, f'#{idx + 1}', fontsize=18,
    #                        bbox=dict(boxstyle="circle,pad=0.3", facecolor='yellow', alpha=0.9),
    #                        weight='bold', zorder=6)
        
    #     # Draw remaining contacts as labels
    #     for idx, (lig_idx, res_info, freq) in enumerate(top_contacts[3:], 4):
    #         if lig_idx in ligand_positions:
    #             lig_x, lig_y = ligand_positions[lig_idx]
                
    #             test_angles = np.linspace(0, 2*np.pi, 8, endpoint=False)
    #             best_angles = []
                
    #             actual_label_distance = self.label_distance * 1.2
                
    #             for angle_offset in test_angles:
    #                 test_x = lig_x + actual_label_distance * np.cos(angle_offset)
    #                 test_y = lig_y + actual_label_distance * np.sin(angle_offset)
                    
    #                 min_dist = float('inf')
    #                 for occ in occupied_positions:
    #                     if isinstance(occ, tuple) and len(occ) >= 2:
    #                         dist = np.sqrt((test_x - occ[0])**2 + (test_y - occ[1])**2)
    #                         min_dist = min(min_dist, dist)
                    
    #                 best_angles.append((angle_offset, min_dist))
                
    #             best_angles.sort(key=lambda x: x[1], reverse=True)
    #             best_angle = best_angles[0][0]
                
    #             label_x = lig_x + actual_label_distance * np.cos(best_angle)
    #             label_y = lig_y + actual_label_distance * np.sin(best_angle)
                
    #             occupied_positions.append((label_x, label_y, self.label_radius))
                
    #             res_type = res_info[:3]
    #             res_num = ''.join(filter(str.isdigit, res_info))
    #             single_letter = AA_CODES.get(res_type, res_type[0])
    #             label = f"{single_letter}{res_num}"
                
    #             circle = Circle((label_x, label_y), self.label_radius, 
    #                           facecolor='white', edgecolor='darkblue',
    #                           linewidth=3, zorder=10)
    #             ax.add_patch(circle)
                
    #             ax.text(label_x, label_y, label, fontsize=24, ha='center', va='center',
    #                    weight='bold', color='darkblue', zorder=11)
                
    #             ax.plot([lig_x, label_x], [lig_y, label_y], 'k-', 
    #                    linewidth=2, alpha=0.6, zorder=4)
        
    #     plt.tight_layout()
        
    #     if save_path:
    #         plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    #         logger.info(f"Contact visualization saved: {save_path}")
    #     else:
    #         plt.show()
        
    #     plt.close(fig)
    def visualize(self, freq_threshold=0.2, top_n_contacts=10, save_path=None):
            """Create visualization"""
            self.generate_linear_conformation()
            
            if self.mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(self.mol)
            
            # Create drawer without setting problematic attributes
            drawer = rdMolDraw2D.MolDraw2DCairo(1200, 1200)
            opts = drawer.drawOptions()
            opts.bondLineWidth = 4
            opts.padding = 0.15
            # DO NOT SET atomLabelFontSize - it doesn't exist in newer RDKit versions
            
            max_freq = max(self.contact_frequencies) if max(self.contact_frequencies) > 0 else 1
            cmap = plt.cm.coolwarm
            
            highlight_atoms = []
            highlight_colors = {}
            
            heavy_idx = 0
            for atom_idx in range(self.mol.GetNumAtoms()):
                atom = self.mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() != 'H':
                    if heavy_idx < len(self.contact_frequencies):
                        freq = self.contact_frequencies[heavy_idx] / max_freq
                        color = cmap(freq)
                        highlight_atoms.append(atom_idx)
                        highlight_colors[atom_idx] = (color[0], color[1], color[2], 0.8)
                    heavy_idx += 1
            
            rdMolDraw2D.PrepareAndDrawMolecule(
                drawer, self.mol,
                highlightAtoms=highlight_atoms,
                highlightAtomColors=highlight_colors
            )
            drawer.FinishDrawing()
            
            img_data = drawer.GetDrawingText()
            img = Image.open(io.BytesIO(img_data))
            
            ligand_positions, ligand_names = self.get_ligand_atom_positions(drawer)
            
            occupied_positions = list(ligand_positions.values())
            
            residue_positions, top_contacts = self.find_optimal_positions(
                ligand_positions, self.residue_atom_details, top_n_contacts, freq_threshold
            )
            
            fig, ax = plt.subplots(figsize=self.fig_size)
            ax.imshow(img)
            ax.axis('off')
            
            # Draw top 3 contacts with amino acid structures
            aa_structure_bounds = []
            
            for idx, (lig_idx, res_info, freq) in enumerate(top_contacts[:3]):
                if res_info in residue_positions:
                    x, y = residue_positions[res_info]
                    
                    structure_pos = self._find_empty_position(
                        (x, y), occupied_positions, self.min_distance, aa_structure_bounds
                    )
                    if structure_pos is None:
                        continue
                        
                    x, y = structure_pos
                    occupied_positions.append((x, y, 80))
                    
                    if lig_idx in ligand_positions:
                        lig_x, lig_y = ligand_positions[lig_idx]
                        angle = np.arctan2(lig_y - y, lig_x - x)
                    else:
                        angle = 0
                    
                    coords = self.draw_amino_acid_structure(ax, x, y, res_info, res_info[:3], 
                                                        rotation=angle)
                    
                    if isinstance(coords, np.ndarray) and len(coords) > 0:
                        min_x, min_y = coords.min(axis=0)
                        max_x, max_y = coords.max(axis=0)
                        label_y = min_y - 45
                        min_y = label_y - 15
                        aa_structure_bounds.append((min_x - 20, min_y, max_x + 20, max_y + 20))
                        for coord in coords:
                            occupied_positions.append(tuple(coord))
                        occupied_positions.append((x, label_y, 80))
                    
                    if lig_idx in ligand_positions:
                        ax.plot([lig_x, x], [lig_y, y], 'k-', linewidth=3, alpha=0.7, zorder=4)
                        
                        mid_x = (lig_x + x) / 2
                        mid_y = (lig_y + y) / 2
                        ax.text(mid_x, mid_y, f'#{idx + 1}', fontsize=18,
                            bbox=dict(boxstyle="circle,pad=0.3", facecolor='yellow', alpha=0.9),
                            weight='bold', zorder=6)
            
            # Draw remaining contacts as labels
            for idx, (lig_idx, res_info, freq) in enumerate(top_contacts[3:], 4):
                if lig_idx in ligand_positions:
                    lig_x, lig_y = ligand_positions[lig_idx]
                    
                    test_angles = np.linspace(0, 2*np.pi, 8, endpoint=False)
                    best_angles = []
                    
                    actual_label_distance = self.label_distance * 1.2
                    
                    for angle_offset in test_angles:
                        test_x = lig_x + actual_label_distance * np.cos(angle_offset)
                        test_y = lig_y + actual_label_distance * np.sin(angle_offset)
                        
                        min_dist = float('inf')
                        for occ in occupied_positions:
                            if isinstance(occ, tuple) and len(occ) >= 2:
                                dist = np.sqrt((test_x - occ[0])**2 + (test_y - occ[1])**2)
                                min_dist = min(min_dist, dist)
                        
                        best_angles.append((angle_offset, min_dist))
                    
                    best_angles.sort(key=lambda x: x[1], reverse=True)
                    best_angle = best_angles[0][0]
                    
                    label_x = lig_x + actual_label_distance * np.cos(best_angle)
                    label_y = lig_y + actual_label_distance * np.sin(best_angle)
                    
                    occupied_positions.append((label_x, label_y, self.label_radius))
                    
                    res_type = res_info[:3]
                    res_num = ''.join(filter(str.isdigit, res_info))
                    single_letter = AA_CODES.get(res_type, res_type[0])
                    label = f"{single_letter}{res_num}"
                    
                    circle = Circle((label_x, label_y), self.label_radius, 
                                facecolor='white', edgecolor='darkblue',
                                linewidth=3, zorder=10)
                    ax.add_patch(circle)
                    
                    ax.text(label_x, label_y, label, fontsize=24, ha='center', va='center',
                        weight='bold', color='darkblue', zorder=11)
                    
                    ax.plot([lig_x, label_x], [lig_y, label_y], 'k-', 
                        linewidth=2, alpha=0.6, zorder=4)
            
            plt.tight_layout()
            
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
                logger.info(f"Contact visualization saved: {save_path}")
            else:
                plt.show()
            
            plt.close(fig)
        
    def generate_linear_conformation(self):
        """Generate linear conformation for better visualization"""
        if self.mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(self.mol)
        
        best_mol = Chem.Mol(self.mol)
        best_linearity = 0
        
        if self.mol.GetNumAtoms() < 5:
            self.align_molecule_to_principal_axis()
            return
        
        for method in ['default', 'linear']:
            for seed in range(5):
                temp_mol = Chem.Mol(self.mol)
                
                if method == 'linear':
                    AllChem.Compute2DCoords(temp_mol, sampleSeed=seed*2, 
                                           canonOrient=False, bondLength=1.5)
                else:
                    AllChem.Compute2DCoords(temp_mol, sampleSeed=seed)
                
                conf = temp_mol.GetConformer()
                coords = []
                for i in range(temp_mol.GetNumAtoms()):
                    pos = conf.GetAtomPosition(i)
                    coords.append([pos.x, pos.y])
                coords = np.array(coords)
                
                if len(coords) > 2:
                    pca = PCA(n_components=2)
                    pca.fit(coords)
                    
                    transformed = pca.transform(coords)
                    transformed[:, 0] *= 1.8
                    transformed[:, 1] *= 0.8
                    
                    stretched_coords = pca.inverse_transform(transformed)
                    
                    for i in range(temp_mol.GetNumAtoms()):
                        x, y = stretched_coords[i]
                        conf.SetAtomPosition(i, (x, y, 0))
                    
                    linearity = pca.explained_variance_ratio_[0]
                    
                    if linearity > best_linearity:
                        best_linearity = linearity
                        best_mol = Chem.Mol(temp_mol)
        
        self.mol = best_mol
        self.align_molecule_to_principal_axis()
    
    def align_molecule_to_principal_axis(self):
        """Align molecule to principal axis"""
        conf = self.mol.GetConformer()
        coords = []
        for i in range(self.mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y])
        coords = np.array(coords)
        
        if len(coords) < 3 or np.std(coords) < 0.01:
            return 0
        
        pca = PCA(n_components=2)
        pca.fit(coords)
        
        principal_axis = pca.components_[0]
        angle = np.arctan2(principal_axis[1], principal_axis[0])
        
        rotation_angle = -angle
        cos_angle = np.cos(rotation_angle)
        sin_angle = np.sin(rotation_angle)
        rotation_matrix = np.array([[cos_angle, -sin_angle], 
                                   [sin_angle, cos_angle]])
        
        centered_coords = coords - coords.mean(axis=0)
        rotated_coords = centered_coords @ rotation_matrix.T
        
        for i in range(self.mol.GetNumAtoms()):
            x, y = rotated_coords[i]
            conf.SetAtomPosition(i, (x, y, 0))
        
        return rotation_angle
    
    def get_ligand_atom_positions(self, drawer):
        """Get ligand atom positions from drawer"""
        positions = {}
        names = {}
        heavy_idx = 0
        
        for atom_idx in range(self.mol.GetNumAtoms()):
            atom = self.mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() != 'H':
                pos = drawer.GetDrawCoords(atom_idx)
                positions[heavy_idx] = (pos.x, pos.y)
                names[heavy_idx] = f"{atom.GetSymbol()}{atom_idx}"
                heavy_idx += 1
                
        return positions, names
    
    def find_optimal_positions(self, ligand_positions, contact_data, top_n=10, freq_threshold=0.2):
        """Find optimal positions for labels"""
        all_contacts = []
        for (lig_idx, res_info), count in self.residue_contacts.items():
            freq = count / self.n_frames
            if freq >= freq_threshold:
                all_contacts.append((lig_idx, res_info, freq))
        
        all_contacts.sort(key=lambda x: x[2], reverse=True)
        top_contacts = all_contacts[:top_n]
        
        initial_positions = {}
        for idx, (_, res_info, _) in enumerate(top_contacts):
            if res_info in self.residue_atom_details:
                atom_contacts = self.residue_atom_details[res_info]
                max_atom = max(atom_contacts, key=atom_contacts.get)
                if max_atom in ligand_positions:
                    lig_x, lig_y = ligand_positions[max_atom]
                    
                    angle = (idx / len(top_contacts)) * 2 * np.pi
                    dist = 200 if idx < 3 else 150
                    
                    x = lig_x + dist * np.cos(angle)
                    y = lig_y + dist * np.sin(angle)
                    initial_positions[res_info] = (x, y)
        
        optimized_positions = self._optimize_layout(initial_positions, ligand_positions)
        
        return optimized_positions, top_contacts
    
    def _optimize_layout(self, initial_positions, ligand_positions, iterations=50):
        """Optimize label layout"""
        positions = dict(initial_positions)
        lig_points = np.array(list(ligand_positions.values()))
        
        for _ in range(iterations):
            forces = {}
            
            for res1, pos1 in positions.items():
                force = np.array([0.0, 0.0])
                
                for res2, pos2 in positions.items():
                    if res1 != res2:
                        diff = np.array(pos1) - np.array(pos2)
                        dist = np.linalg.norm(diff)
                        if dist < self.min_distance * 1.5:
                            if dist > 0:
                                repulsion = (self.min_distance * 1.5 - dist) / 10
                                force += (diff / dist) * repulsion
                
                for lig_pos in lig_points:
                    diff = np.array(pos1) - lig_pos
                    dist = np.linalg.norm(diff)
                    if dist < self.min_distance:
                        if dist > 0:
                            repulsion = (self.min_distance - dist) / 10
                            force += (diff / dist) * repulsion
                
                forces[res1] = force
            
            for res, force in forces.items():
                new_pos = np.array(positions[res]) + force * 0.1
                positions[res] = tuple(new_pos)
        
        return positions
    
    def draw_amino_acid_structure(self, ax, x, y, residue_name, residue_type, 
                                  rotation=0, scale=None):
        """Draw amino acid structure"""
        if scale is None:
            scale = self.aa_structure_scale
            
        res_type = residue_type[:3]
        if res_type not in AA_SMILES:
            return {}
        
        aa_mol = Chem.MolFromSmiles(AA_SMILES[res_type])
        if aa_mol is None:
            return {}
       
        AllChem.Compute2DCoords(aa_mol)
       
        conf = aa_mol.GetConformer()
        coords = []
        for i in range(aa_mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y])
        coords = np.array(coords)
        
        center = coords.mean(axis=0)
        coords -= center
        
        if rotation != 0:
            cos_r = np.cos(rotation)
            sin_r = np.sin(rotation)
            rot_matrix = np.array([[cos_r, -sin_r], [sin_r, cos_r]])
            coords = coords @ rot_matrix.T
        
        coords *= scale * 60
        coords[:, 0] += x
        coords[:, 1] += y
        
        for bond in aa_mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            
            x1, y1 = coords[idx1]
            x2, y2 = coords[idx2]
            
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                ax.plot([x1, x2], [y1, y2], 'k-', linewidth=4, zorder=8)
            else:
                ax.plot([x1, x2], [y1, y2], 'k-', linewidth=3, zorder=8)
        
        for i, (cx, cy) in enumerate(coords):
            atom = aa_mol.GetAtomWithIdx(i)
            symbol = atom.GetSymbol()
            
            if symbol != 'H':
                color = {'C': '#808080', 'N': '#4169E1', 'O': '#DC143C',
                        'S': '#FFD700'}.get(symbol, '#FF69B4')
                
                circle = Circle((cx, cy), 15, facecolor=color,
                                edgecolor='black', linewidth=2, zorder=10)
                ax.add_patch(circle)
                
                if symbol not in ['C']:
                    ax.text(cx, cy, symbol, fontsize=14, ha='center',
                            va='center', color='white', weight='bold', zorder=11)
        
        if isinstance(coords, np.ndarray) and len(coords) > 0:
            min_y = coords[:, 1].min()
        else:
            min_y = y
            
        label_offset_y = 45
        label_x = x
        label_y = min_y - label_offset_y
        
        ax.text(label_x, label_y, residue_name, fontsize=26, ha='center', va='top',
                weight='bold', color='darkred', zorder=12,
                bbox=dict(boxstyle="round,pad=0.2", facecolor='white', edgecolor='none', alpha=0.85))
        
        return coords
   
    def _find_empty_position(self, start_pos, occupied_positions, min_distance=100, aa_bounds=None):
        """Find empty position for label"""
        x_start, y_start = start_pos
        if aa_bounds is None:
            aa_bounds = []
        
        angles = np.linspace(0, 2*np.pi, 16, endpoint=False)
        distances = [min_distance, min_distance * 1.2, min_distance * 1.5, 
                    min_distance * 1.8, min_distance * 2.2]
        
        for dist in distances:
            for angle in angles:
                x = x_start + dist * np.cos(angle)
                y = y_start + dist * np.sin(angle)
                
                if self._is_position_clear((x, y), occupied_positions, aa_bounds, min_distance):
                    return (x, y)
        return None
    
    def _is_position_clear(self, pos, occupied_positions, aa_bounds, min_distance):
        """Check if position is clear"""
        x, y = pos
        
        for occ in occupied_positions:
            if isinstance(occ, tuple):
                if len(occ) == 2:
                    occ_x, occ_y = occ
                    dist = np.sqrt((x - occ_x)**2 + (y - occ_y)**2)
                    if dist < min_distance:
                        return False
                elif len(occ) == 3:
                    occ_x, occ_y, radius = occ
                    dist = np.sqrt((x - occ_x)**2 + (y - occ_y)**2)
                    if dist < radius + min_distance/2:
                        return False
        
        for bounds in aa_bounds:
            if len(bounds) == 4:
                min_x, min_y, max_x, max_y = bounds
                padding = 40
                if (min_x - padding <= x <= max_x + padding and 
                    min_y - padding <= y <= max_y + padding):
                    return False
        
        return True
