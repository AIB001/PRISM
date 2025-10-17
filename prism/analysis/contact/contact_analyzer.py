#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Contact Analysis Module
Fast contact analysis using vectorized operations
"""

import numpy as np
import mdtraj as md
from .utils import convert_numpy_types

class FastContactAnalyzer:
    """Fast contact analysis using vectorized operations"""
    
    def __init__(self, trajectory, config=None):
        self.traj = trajectory
        self.config = config or self._default_config()
        
    def _default_config(self):
        """Default configuration"""
        class Config:
            contact_enter_threshold_nm = 0.35
            contact_exit_threshold_nm = 0.4
            distance_cutoff_nm = 0.5
            min_frames_for_smoothing = 10
            smooth_window = 11
        return Config()
    
    def identify_ligand_residue(self):
        """Automatically identify ligand residue using PRISM utilities"""
        from ...utils.ligand import identify_ligand_residue
        return identify_ligand_residue(self.traj)
    
    def get_heavy_atoms(self, residue):
        """Get heavy atom indices from residue"""
        heavy_atoms = []
        for atom in residue.atoms:
            if atom.element.symbol != 'H':
                heavy_atoms.append(atom.index)
        return heavy_atoms
    
    def get_protein_heavy_atoms(self):
        """Get protein heavy atom indices"""
        protein_atoms = []
        standard_aa = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        for atom in self.traj.topology.atoms:
            if (atom.residue.name in standard_aa and 
                atom.element.symbol in ['C', 'N', 'O', 'S']):
                protein_atoms.append(atom.index)
        
        return protein_atoms
    
    def calculate_contact_proportions(self):
        """Calculate contact proportions using efficient vectorized operations"""
        ligand_residue = self.identify_ligand_residue()
        if not ligand_residue:
            raise ValueError("Could not identify ligand residue")
        
        print(f"Analyzing ligand: {ligand_residue.name}{ligand_residue.resSeq}")
        
        # Get atom indices
        ligand_atoms = self.get_heavy_atoms(ligand_residue)
        protein_atoms = self.get_protein_heavy_atoms()
        
        print(f"Ligand atoms: {len(ligand_atoms)}, Protein atoms: {len(protein_atoms)}")
        
        # Create atom pairs for distance calculation
        atom_pairs = []
        pair_info = []  # Store (ligand_atom_idx, protein_atom_idx, residue_id)
        
        for lig_atom in ligand_atoms:
            for prot_atom in protein_atoms:
                atom_pairs.append([lig_atom, prot_atom])
                residue = self.traj.topology.atom(prot_atom).residue
                residue_id = f"{residue.name}{residue.resSeq}"
                pair_info.append((lig_atom, prot_atom, residue_id))
        
        print(f"Calculating distances for {len(atom_pairs)} atom pairs...")
        
        # Sample frames for efficiency
        n_frames = self.traj.n_frames
        if n_frames > 1000:
            frame_indices = np.linspace(0, n_frames-1, 1000, dtype=int)
            traj_sample = self.traj[frame_indices]
        else:
            traj_sample = self.traj
            frame_indices = np.arange(n_frames)
        
        # Calculate all distances at once
        distances = md.compute_distances(traj_sample, atom_pairs, opt=True)
        print(f"Computed distances shape: {distances.shape}")
        
        # Find contacts
        contact_threshold = self.config.contact_enter_threshold_nm
        contacts = distances < contact_threshold
        
        # Count contacts per residue
        contact_counts = {}
        residue_contacts = {}
        ligand_atom_contacts = {}
        residue_ligand_atoms = {}  # Track which ligand atoms contact each residue
        contact_distances = {}
        
        for i, (lig_atom, prot_atom, residue_id) in enumerate(pair_info):
            contact_frames = int(np.sum(contacts[:, i]))
            
            if contact_frames > 0:
                # Calculate average distance when in contact
                contact_mask = contacts[:, i]
                avg_distance = float(np.mean(distances[contact_mask, i]))
                
                # Store detailed contact info
                key = (lig_atom, residue_id)
                if key not in contact_counts:
                    contact_counts[key] = 0
                    contact_distances[key] = []
                contact_counts[key] += contact_frames
                contact_distances[key].append(avg_distance)
                
                # Aggregate by residue
                if residue_id not in residue_contacts:
                    residue_contacts[residue_id] = 0
                    residue_ligand_atoms[residue_id] = {}
                residue_contacts[residue_id] += contact_frames
                
                # Track which ligand atoms contact this residue
                if lig_atom not in residue_ligand_atoms[residue_id]:
                    residue_ligand_atoms[residue_id][lig_atom] = 0
                residue_ligand_atoms[residue_id][lig_atom] += contact_frames
                
                # Aggregate by ligand atom
                if lig_atom not in ligand_atom_contacts:
                    ligand_atom_contacts[lig_atom] = 0
                ligand_atom_contacts[lig_atom] += contact_frames
        
        # Convert to frequencies
        n_frames_analyzed = len(frame_indices)
        contact_frequencies = {}
        avg_contact_distances = {}
        
        for key, count in contact_counts.items():
            freq = float(count / n_frames_analyzed)
            if freq >= 0.1:  # Only keep contacts with >10% frequency
                contact_frequencies[key] = freq
                avg_contact_distances[key] = float(np.mean(contact_distances[key]))
        
        # Calculate residue-level proportions
        residue_proportions = {}
        residue_avg_distances = {}
        residue_best_ligand_atoms = {}
        
        for residue_id, count in residue_contacts.items():
            # Find the ligand atom with most contacts to this residue
            if residue_id in residue_ligand_atoms:
                best_lig_atom = max(residue_ligand_atoms[residue_id].items(), key=lambda x: x[1])[0]
                residue_best_ligand_atoms[residue_id] = best_lig_atom
            
            # Count unique ligand atoms that contact this residue
            ligand_atoms_contacting = set()
            residue_distances = []
            for (lig_atom, res_id), freq in contact_frequencies.items():
                if res_id == residue_id:
                    ligand_atoms_contacting.add(lig_atom)
                    if (lig_atom, res_id) in avg_contact_distances:
                        residue_distances.append(avg_contact_distances[(lig_atom, res_id)])
            
            if ligand_atoms_contacting:
                # Frequency = contact frames / (total frames * ligand atoms involved)
                max_possible = n_frames_analyzed * len(ligand_atoms_contacting)
                proportion = float(count / max_possible)
                if proportion >= 0.02:  # Only keep residues with >2% contact
                    residue_proportions[residue_id] = proportion
                    if residue_distances:
                        residue_avg_distances[residue_id] = float(np.mean(residue_distances))
        
        print(f"Found {len(contact_frequencies)} significant atom-residue contacts")
        print(f"Found {len(residue_proportions)} contacting residues")
        
        # Convert all numpy types to Python types
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