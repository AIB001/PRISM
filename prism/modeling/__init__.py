#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Modeling Module - Advanced molecular modeling and system preparation
"""

from .base import ModelingBase
from .homology import HomologyModeler
from .docking import DockingEngine
from .conformer import ConformerGenerator
from .structure import StructureOptimizer

def create_model(target_sequence, template_pdb=None, ligand_smiles=None, **kwargs):
    """
    Create a molecular model for MD simulation.
    
    Parameters:
    -----------
    target_sequence : str
        Target protein sequence or path to FASTA file
    template_pdb : str, optional
        Path to template PDB file for homology modeling
    ligand_smiles : str, optional
        SMILES string for ligand generation
    **kwargs : optional
        Additional modeling parameters
        
    Returns:
    --------
    ModelingBase
        A modeling object with methods to generate structures
        
    Examples:
    ---------
    >>> import prism as pm
    >>> modeler = pm.modeling.create_model("MKLLIV...", template_pdb="template.pdb")
    >>> modeler.build_homology_model()
    """
    return HomologyModeler(target_sequence, template_pdb, ligand_smiles, **kwargs)

def dock_ligand(protein_pdb, ligand_smiles, binding_site=None, **kwargs):
    """
    Dock a ligand into a protein structure.
    
    Parameters:
    -----------
    protein_pdb : str
        Path to protein PDB file
    ligand_smiles : str
        SMILES string for ligand
    binding_site : dict, optional
        Binding site coordinates or residue selection
    **kwargs : optional
        Docking parameters
        
    Returns:
    --------
    DockingEngine
        A docking object with methods to perform docking
    """
    return DockingEngine(protein_pdb, ligand_smiles, binding_site, **kwargs)

def generate_conformers(smiles, num_conformers=100, **kwargs):
    """
    Generate conformers for a ligand.
    
    Parameters:
    -----------
    smiles : str
        SMILES string for ligand
    num_conformers : int, optional
        Number of conformers to generate (default: 100)
    **kwargs : optional
        Conformer generation parameters
        
    Returns:
    --------
    ConformerGenerator
        A conformer generator object
    """
    return ConformerGenerator(smiles, num_conformers, **kwargs)

def optimize_structure(input_file, force_field="amber99sb-ildn", **kwargs):
    """
    Optimize a protein or ligand structure.
    
    Parameters:
    -----------
    input_file : str
        Path to input structure file
    force_field : str, optional
        Force field for optimization (default: "amber99sb-ildn")
    **kwargs : optional
        Optimization parameters
        
    Returns:
    --------
    StructureOptimizer
        A structure optimizer object
    """
    return StructureOptimizer(input_file, force_field, **kwargs)

__all__ = [
    "ModelingBase",
    "HomologyModeler", 
    "DockingEngine",
    "ConformerGenerator",
    "StructureOptimizer",
    "create_model",
    "dock_ligand",
    "generate_conformers",
    "optimize_structure"
]