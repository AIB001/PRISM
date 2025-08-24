#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Base class for force field generators in PRISM
"""

from abc import ABC, abstractmethod
import os
from pathlib import Path


class ForceFieldGeneratorBase(ABC):
    """Abstract base class for force field generators"""
    
    def __init__(self, ligand_path, output_dir, overwrite=False):
        """
        Initialize the force field generator
        
        Parameters:
        -----------
        ligand_path : str
            Path to the ligand file
        output_dir : str
            Directory where output files will be stored
        overwrite : bool
            Whether to overwrite existing files
        """
        self.ligand_path = os.path.abspath(ligand_path)
        self.output_dir = os.path.abspath(output_dir)
        self.overwrite = overwrite
        self.ligand_name = Path(ligand_path).stem
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
    
    @abstractmethod
    def run(self):
        """
        Run the force field generation workflow
        
        Returns:
        --------
        str : Path to the output directory containing the generated files
        """
        pass
    
    @abstractmethod
    def get_output_dir_name(self):
        """
        Get the name of the output directory for this force field
        
        Returns:
        --------
        str : Directory name (e.g., 'LIG.amb2gmx' or 'LIG.openff2gmx')
        """
        pass
    
    def check_required_files(self, output_dir):
        """
        Check if all required files exist in the output directory
        
        Parameters:
        -----------
        output_dir : str
            Directory to check
        
        Returns:
        --------
        bool : True if all required files exist
        """
        required_files = [
            "LIG.gro",
            "LIG.itp",
            "LIG.top",
            "atomtypes_LIG.itp",
            "posre_LIG.itp"
        ]
        
        for filename in required_files:
            filepath = os.path.join(output_dir, filename)
            if not os.path.exists(filepath):
                return False
        
        return True