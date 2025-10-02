#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Base classes for PRISM modeling module
"""

import os
import logging
from abc import ABC, abstractmethod
from pathlib import Path

logger = logging.getLogger(__name__)


class ModelingBase(ABC):
    """
    Abstract base class for all modeling operations.
    
    Provides common functionality for structure preparation, validation,
    and output management.
    """
    
    def __init__(self, output_dir="modeling_output", **kwargs):
        """
        Initialize modeling base.
        
        Parameters:
        -----------
        output_dir : str
            Directory for output files
        **kwargs : optional
            Additional parameters
        """
        self.output_dir = os.path.abspath(output_dir)
        self.kwargs = kwargs
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Setup logging
        self._setup_logging()
        
        # Track generated files
        self.generated_files = {}
        
        logger.info(f"Initialized {self.__class__.__name__} with output_dir: {self.output_dir}")
    
    def _setup_logging(self):
        """Setup module-specific logging"""
        log_file = os.path.join(self.output_dir, "modeling.log")
        
        # Create file handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        file_handler.setFormatter(formatter)
        
        # Add handler to logger
        logger.addHandler(file_handler)
        logger.setLevel(logging.DEBUG)
    
    @abstractmethod
    def prepare_system(self):
        """
        Prepare the molecular system for modeling.
        
        Returns:
        --------
        dict
            Dictionary of prepared system components
        """
        pass
    
    @abstractmethod
    def validate_inputs(self):
        """
        Validate input files and parameters.
        
        Returns:
        --------
        bool
            True if all inputs are valid
        """
        pass
    
    def save_config(self, config_dict, filename="modeling_config.yaml"):
        """
        Save configuration to YAML file.
        
        Parameters:
        -----------
        config_dict : dict
            Configuration dictionary
        filename : str
            Output filename
        """
        import yaml
        
        config_path = os.path.join(self.output_dir, filename)
        
        with open(config_path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, indent=2)
        
        self.generated_files['config'] = config_path
        logger.info(f"Configuration saved to: {config_path}")
        
        return config_path
    
    def cleanup(self, keep_important=True):
        """
        Clean up temporary files.
        
        Parameters:
        -----------
        keep_important : bool
            Whether to keep important output files
        """
        if keep_important:
            # Only remove temporary files
            temp_patterns = ["*.tmp", "*.temp", "*~", "*.bak"]
            for pattern in temp_patterns:
                for temp_file in Path(self.output_dir).glob(pattern):
                    temp_file.unlink()
                    logger.info(f"Removed temporary file: {temp_file}")
        else:
            # Remove all files in output directory
            import shutil
            shutil.rmtree(self.output_dir)
            logger.info(f"Removed output directory: {self.output_dir}")
    
    def get_output_files(self):
        """
        Get dictionary of generated output files.
        
        Returns:
        --------
        dict
            Dictionary mapping file types to paths
        """
        # Update with files that exist
        existing_files = {}
        for file_type, file_path in self.generated_files.items():
            if os.path.exists(file_path):
                existing_files[file_type] = file_path
        
        return existing_files
    
    def summary(self):
        """Print a summary of the modeling operation"""
        print(f"\n{self.__class__.__name__} Summary:")
        print(f"Output directory: {self.output_dir}")
        
        files = self.get_output_files()
        if files:
            print("Generated files:")
            for file_type, file_path in files.items():
                print(f"  {file_type}: {os.path.basename(file_path)}")
        else:
            print("No files generated yet.")


class ModelingError(Exception):
    """Custom exception for modeling errors"""
    pass


class ValidationError(ModelingError):
    """Exception raised for input validation errors"""
    pass