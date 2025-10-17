#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Data Processing Module
Handle ligand loading and statistics calculation
"""

# Try to import RDKit for molecular structure
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

class DataProcessor:
    """Process molecular data and calculate statistics"""
    
    def load_ligand_structure(self, ligand_path):
        """Load ligand structure for visualization"""
        ligand_mol = None
        if not RDKIT_AVAILABLE:
            return None
            
        try:
            if ligand_path.endswith('.sdf'):
                supplier = Chem.SDMolSupplier(ligand_path, removeHs=False)
                ligand_mol = supplier[0] if supplier else None
            elif ligand_path.endswith('.mol'):
                ligand_mol = Chem.MolFromMolFile(ligand_path, removeHs=False)
            elif ligand_path.endswith('.mol2'):
                ligand_mol = Chem.MolFromMol2File(ligand_path, removeHs=False)
                
            if ligand_mol:
                try:
                    Chem.SanitizeMol(ligand_mol)
                    print(f"Ligand loaded: {ligand_mol.GetNumAtoms()} atoms")
                except:
                    print("Warning: Could not sanitize ligand molecule")
        except Exception as e:
            print(f"Error loading ligand: {e}")
            
        return ligand_mol
    
    def calculate_statistics(self, contacts):
        """Calculate statistics from contact data"""
        total_contacts = len(contacts)
        high_freq = len([c for c in contacts if c['frequency'] > 0.6])
        max_freq = max([c['frequency'] for c in contacts]) if contacts else 0
        
        return {
            'total_contacts': total_contacts,
            'high_freq_contacts': high_freq,
            'max_freq_percent': int(max_freq * 100)
        }