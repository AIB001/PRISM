#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visualization Generator Module
Generate ligand and contact data for visualization
"""

import numpy as np
from .utils import convert_numpy_types  

# Try to import RDKit for molecular structure
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Molecular visualization will be simplified.")

class VisualizationGenerator:
    """Generate visualization data for HTML rendering"""
    
    def generate_ligand_data(self, contact_results, ligand_mol=None):
        """Generate ligand visualization data including hydrogens and both 2D/3D coordinates"""
        ligand_data = {'atoms': [], 'bonds': [], 'elements': []}
        
        if ligand_mol and RDKIT_AVAILABLE:
            # Generate both 2D and 3D coordinates
            AllChem.Compute2DCoords(ligand_mol)
            conf_2d = ligand_mol.GetConformer()
            
            # Store 2D coordinates
            coords_2d = []
            for atom_idx in range(ligand_mol.GetNumAtoms()):
                pos = conf_2d.GetAtomPosition(atom_idx)
                coords_2d.append([pos.x, pos.y])
            coords_2d = np.array(coords_2d)
            
            # Try to generate 3D coordinates
            mol_3d = Chem.Mol(ligand_mol)
            has_3d = False
            coords_3d = coords_2d.copy()
            
            try:
                AllChem.EmbedMolecule(mol_3d, randomSeed=42)
                if mol_3d.GetNumConformers() > 0:
                    AllChem.MMFFOptimizeMolecule(mol_3d)
                    conf_3d = mol_3d.GetConformer()
                    coords_3d = []
                    for atom_idx in range(mol_3d.GetNumAtoms()):
                        pos = conf_3d.GetAtomPosition(atom_idx)
                        coords_3d.append([pos.x, pos.y, pos.z])
                    coords_3d = np.array(coords_3d)
                    has_3d = True
            except:
                coords_3d = np.column_stack([coords_2d, np.zeros(len(coords_2d))])
            
            # Normalize 3D coordinates
            if has_3d:
                center_3d = coords_3d.mean(axis=0)
                coords_3d -= center_3d
                max_dist_3d = np.max(np.abs(coords_3d))
                if max_dist_3d > 0:
                    coords_3d = coords_3d / max_dist_3d * 50
            
            # Map ALL atoms (including hydrogens)
            atom_map = {}
            atom_idx = 0
            
            for mol_atom_idx in range(ligand_mol.GetNumAtoms()):
                atom = ligand_mol.GetAtomWithIdx(mol_atom_idx)
                
                # Count contacts for this atom
                contact_count = 0
                traj_atom_idx = -1
                if atom.GetSymbol() != 'H':
                    heavy_idx = len([a for a in range(mol_atom_idx) 
                                    if ligand_mol.GetAtomWithIdx(a).GetSymbol() != 'H'])
                    if heavy_idx < len(contact_results['ligand_atoms']):
                        traj_atom_idx = contact_results['ligand_atoms'][heavy_idx]
                        contact_count = contact_results['ligand_atom_contacts'].get(traj_atom_idx, 0)
                
                radius = 18 if atom.GetSymbol() != 'H' else 10
                
                ligand_data['atoms'].append({
                    'id': f'L{atom_idx}',
                    'x': float(coords_2d[mol_atom_idx][0] * 30),
                    'y': float(coords_2d[mol_atom_idx][1] * 30),
                    'x3d': float(coords_3d[mol_atom_idx][0]),
                    'y3d': float(coords_3d[mol_atom_idx][1]),
                    'z3d': float(coords_3d[mol_atom_idx][2]),
                    'element': str(atom.GetSymbol()),
                    'radius': radius,
                    'contacts': int(contact_count),
                    'trajIndex': traj_atom_idx
                })
                atom_map[mol_atom_idx] = atom_idx
                atom_idx += 1
            
            # Generate bonds
            for bond in ligand_mol.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                
                if begin_idx in atom_map and end_idx in atom_map:
                    ligand_data['bonds'].append([f'L{atom_map[begin_idx]}', f'L{atom_map[end_idx]}'])
        else:
            # Create simplified structure when RDKit not available
            n_atoms = min(len(contact_results['ligand_atoms']), 8)
            for i in range(n_atoms):
                angle = i * 2 * np.pi / n_atoms
                radius = 40
                z = 20 * np.sin(i * np.pi / 4)
                contact_count = contact_results['ligand_atom_contacts'].get(
                    contact_results['ligand_atoms'][i], 0) if i < len(contact_results['ligand_atoms']) else 0
                
                ligand_data['atoms'].append({
                    'id': f'L{i}',
                    'x': float(radius * np.cos(angle)),
                    'y': float(radius * np.sin(angle)),
                    'x3d': float(radius * np.cos(angle)),
                    'y3d': float(radius * np.sin(angle)),
                    'z3d': float(z),
                    'element': 'C',
                    'radius': 18,
                    'contacts': int(contact_count),
                    'trajIndex': contact_results['ligand_atoms'][i] if i < len(contact_results['ligand_atoms']) else -1
                })
            
            for i in range(n_atoms):
                ligand_data['bonds'].append([f'L{i}', f'L{(i+1)%n_atoms}'])
        
        return ligand_data
    
    def generate_contact_data(self, contact_results, ligand_data, max_contacts=20):
        """Generate contact data for HTML with proper alignment and TOP3 marking"""
        residue_prop = contact_results['residue_proportions']
        residue_distances = contact_results.get('residue_avg_distances', {})
        residue_best_ligand = contact_results.get('residue_best_ligand_atoms', {})
        
        # Sort residues by contact proportion
        sorted_residues = sorted(residue_prop.items(), key=lambda x: x[1], reverse=True)
        top_residues = sorted_residues[:max_contacts]
        
        # Identify global TOP3 residues based on frequency
        top3_residue_ids = set()
        if len(sorted_residues) >= 1:
            top3_residue_ids.add(sorted_residues[0][0])  # 1st
        if len(sorted_residues) >= 2:
            top3_residue_ids.add(sorted_residues[1][0])  # 2nd
        if len(sorted_residues) >= 3:
            top3_residue_ids.add(sorted_residues[2][0])  # 3rd
        
        contacts = []
        
        for idx, (residue_id, proportion) in enumerate(top_residues):
            # Get the best ligand atom for this residue
            best_ligand_atom_traj = residue_best_ligand.get(residue_id, contact_results['ligand_atoms'][0])
            
            # Find the corresponding ligand atom in our data
            best_ligand_atom_id = None
            for atom in ligand_data['atoms']:
                if atom.get('trajIndex') == best_ligand_atom_traj:
                    best_ligand_atom_id = atom['id']
                    break
            
            if not best_ligand_atom_id:
                best_ligand_atom_id = 'L0'
            
            # Extract residue type
            residue_type = residue_id[:3] if len(residue_id) >= 3 else residue_id
            
            # Get average distance
            avg_distance_nm = residue_distances.get(residue_id, 0.35)
            
            # Calculate pixel distance
            min_dist_nm = 0.2
            max_dist_nm = 0.8
            clamped_dist = min(max(avg_distance_nm, min_dist_nm), max_dist_nm)
            sqrt_dist = np.sqrt((clamped_dist - min_dist_nm) / (max_dist_nm - min_dist_nm))
            pixel_distance = 250 + sqrt_dist * 150
            
            # Check if this residue is in TOP3
            is_top3 = residue_id in top3_residue_ids
            
            contacts.append({
                'id': str(residue_id.replace(' ', '_')),
                'frequency': float(min(proportion, 1.0)),
                'residueType': str(residue_type),
                'ligandAtom': best_ligand_atom_id,
                'residue': str(residue_id),
                'avgDistance': float(avg_distance_nm),
                'pixelDistance': float(pixel_distance),
                'isTop3': is_top3  # Mark if this is a TOP3 contact
            })
        
        return contacts
