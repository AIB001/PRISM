#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Core - High-level API for protein-ligand system building
"""

import os
import yaml
from pathlib import Path
from .builder import PRISMBuilder


class PRISMSystem:
    """
    High-level interface for building protein-ligand systems.
    
    This class provides a simplified API for the PRISMBuilder functionality.
    """
    
    def __init__(self, protein_path, ligand_path, output_dir="prism_output", 
                 config=None, ligand_forcefield="gaff", **kwargs):
        """
        Initialize a PRISM system.
        
        Parameters:
        -----------
        protein_path : str
            Path to protein PDB file
        ligand_path : str
            Path to ligand file (MOL2/SDF)
        output_dir : str, optional
            Output directory (default: "prism_output")
        config : str or dict, optional
            Configuration file path or dict
        ligand_forcefield : str, optional
            Ligand force field ("gaff" or "openff", default: "gaff")
        **kwargs : optional
            Additional parameters for PRISMBuilder
        """
        self.protein_path = os.path.abspath(protein_path)
        self.ligand_path = os.path.abspath(ligand_path)
        self.output_dir = os.path.abspath(output_dir)
        self.ligand_forcefield = ligand_forcefield
        
        # Validate input files
        self._validate_inputs()
        
        # Process configuration
        self.config_path = self._process_config(config)
        
        # Store additional parameters
        self.builder_kwargs = kwargs
        
        # Initialize builder (will be created when needed)
        self._builder = None
        
        # Track build status
        self._is_built = False
        self._build_output = None
        
        print(f"PRISM System initialized:")
        print(f"  Protein: {self.protein_path}")
        print(f"  Ligand: {self.ligand_path}")  
        print(f"  Output: {self.output_dir}")
        print(f"  Ligand FF: {self.ligand_forcefield}")
    
    def _validate_inputs(self):
        """Validate input files exist and have correct extensions"""
        if not os.path.exists(self.protein_path):
            raise FileNotFoundError(f"Protein file not found: {self.protein_path}")
        
        if not os.path.exists(self.ligand_path):
            raise FileNotFoundError(f"Ligand file not found: {self.ligand_path}")
        
        # Check file extensions
        protein_ext = Path(self.protein_path).suffix.lower()
        if protein_ext != '.pdb':
            print(f"Warning: Protein file extension is {protein_ext}, expected .pdb")
        
        ligand_ext = Path(self.ligand_path).suffix.lower()
        if ligand_ext not in ['.mol2', '.sdf', '.sd']:
            print(f"Warning: Ligand file extension is {ligand_ext}, expected .mol2 or .sdf")
    
    def _process_config(self, config):
        """Process configuration input"""
        if config is None:
            return None
        elif isinstance(config, str):
            # Config file path
            if os.path.exists(config):
                return os.path.abspath(config)
            else:
                print(f"Warning: Config file not found: {config}")
                return None
        elif isinstance(config, dict):
            # Config dictionary - save to temporary file
            config_file = os.path.join(self.output_dir, "temp_config.yaml")
            os.makedirs(self.output_dir, exist_ok=True)
            with open(config_file, 'w') as f:
                yaml.dump(config, f, default_flow_style=False)
            return config_file
        else:
            print(f"Warning: Invalid config type: {type(config)}")
            return None
    
    @property
    def builder(self):
        """Get or create the PRISMBuilder instance"""
        if self._builder is None:
            self._builder = PRISMBuilder(
                protein_path=self.protein_path,
                ligand_path=self.ligand_path,
                output_dir=self.output_dir,
                ligand_forcefield=self.ligand_forcefield,
                config_path=self.config_path,
                **self.builder_kwargs
            )
        return self._builder
    
    def build(self, cleanup=True):
        """
        Build the complete protein-ligand system.
        
        Parameters:
        -----------
        cleanup : bool, optional
            Whether to clean up temporary files (default: True)
            
        Returns:
        --------
        str
            Path to the output directory
        """
        print("\n" + "="*60)
        print("Building PRISM system...")
        print("="*60)
        
        try:
            # Run the complete workflow
            self._build_output = self.builder.run()
            self._is_built = True
            
            if not cleanup:
                print("Skipping cleanup (cleanup=False)")
            
            print(f"\nSystem build completed successfully!")
            print(f"Output directory: {self._build_output}")
            
            return self._build_output
            
        except Exception as e:
            print(f"Error building system: {e}")
            raise
    
    def build_step_by_step(self):
        """
        Build the system step by step with user control.
        
        Returns:
        --------
        dict
            Dictionary with paths to generated files
        """
        if self._is_built:
            print("System already built. Use rebuild() to rebuild.")
            return self.get_output_files()
        
        results = {}
        
        try:
            print("Step 1: Generating ligand force field...")
            lig_ff_dir = self.builder.generate_ligand_forcefield()
            results['ligand_forcefield'] = lig_ff_dir
            
            print("Step 2: Cleaning protein...")
            cleaned_protein = self.builder.clean_protein()
            results['cleaned_protein'] = cleaned_protein
            
            print("Step 3: Building GROMACS model...")
            model_dir = self.builder.build_model(cleaned_protein)
            results['model_directory'] = model_dir
            
            print("Step 4: Generating MDP files...")
            self.builder.generate_mdp_files()
            results['mdp_directory'] = self.builder.mdp_dir
            
            print("Step 5: Cleaning up...")
            self.builder.cleanup()
            
            self._is_built = True
            self._build_output = self.output_dir
            results['output_directory'] = self.output_dir
            
            return results
            
        except Exception as e:
            print(f"Error in step-by-step build: {e}")
            raise
    
    def rebuild(self, **kwargs):
        """
        Rebuild the system with optional parameter changes.
        
        Parameters:
        -----------
        **kwargs : optional
            Parameters to update before rebuilding
        """
        # Update parameters
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                self.builder_kwargs[key] = value
        
        # Reset builder to apply new parameters
        self._builder = None
        self._is_built = False
        
        # Rebuild
        return self.build()
    
    def get_output_files(self):
        """
        Get paths to important output files.
        
        Returns:
        --------
        dict
            Dictionary with paths to key files
        """
        if not self._is_built:
            print("System not built yet. Run build() first.")
            return {}
        
        files = {}
        
        # Main system files
        gromacs_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
        if os.path.exists(gromacs_dir):
            files['gromacs_directory'] = gromacs_dir
            files['system_gro'] = os.path.join(gromacs_dir, "solv_ions.gro")
            files['system_top'] = os.path.join(gromacs_dir, "topol.top")
        
        # MDP files
        mdp_dir = os.path.join(self.output_dir, "mdps")
        if os.path.exists(mdp_dir):
            files['mdp_directory'] = mdp_dir
            files['em_mdp'] = os.path.join(mdp_dir, "em.mdp")
            files['nvt_mdp'] = os.path.join(mdp_dir, "nvt.mdp")
            files['npt_mdp'] = os.path.join(mdp_dir, "npt.mdp")
            files['md_mdp'] = os.path.join(mdp_dir, "md.mdp")
        
        # Force field files
        if self.ligand_forcefield == "gaff":
            ff_dir = os.path.join(self.output_dir, "LIG.amb2gmx")
        else:
            ff_dir = os.path.join(self.output_dir, "LIG.openff2gmx")
        
        if os.path.exists(ff_dir):
            files['ligand_ff_directory'] = ff_dir
            files['ligand_itp'] = os.path.join(ff_dir, "LIG.itp")
            files['ligand_gro'] = os.path.join(ff_dir, "LIG.gro")
        
        # Configuration
        config_file = os.path.join(self.output_dir, "prism_config.yaml")
        if os.path.exists(config_file):
            files['config'] = config_file
        
        return files
    
    def generate_run_script(self, script_path=None):
        """
        Generate a shell script to run the MD simulation.
        
        Parameters:
        -----------
        script_path : str, optional
            Path for the script file (default: output_dir/run_md.sh)
            
        Returns:
        --------
        str
            Path to the generated script
        """
        if not self._is_built:
            raise RuntimeError("System not built yet. Run build() first.")
        
        if script_path is None:
            script_path = os.path.join(self.output_dir, "run_md.sh")
        
        script_content = f"""#!/bin/bash
# PRISM Generated MD Simulation Script
# Generated for system: {os.path.basename(self.protein_path)} + {os.path.basename(self.ligand_path)}

set -e  # Exit on any error

# Navigate to simulation directory
cd {os.path.join(self.output_dir, "GMX_PROLIG_MD")}

echo "Starting molecular dynamics simulation..."

# Energy minimization
echo "Step 1: Energy Minimization"
gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em

# NVT equilibration  
echo "Step 2: NVT Equilibration"
gmx grompp -f ../mdps/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt

# NPT equilibration
echo "Step 3: NPT Equilibration" 
gmx grompp -f ../mdps/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt

# Production MD
echo "Step 4: Production MD"
gmx grompp -f ../mdps/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1
gmx mdrun -deffnm md

echo "Simulation completed successfully!"
echo "Results are in: $(pwd)"
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make script executable
        os.chmod(script_path, 0o755)
        
        print(f"MD run script generated: {script_path}")
        return script_path
    
    def __repr__(self):
        """String representation of the system"""
        status = "Built" if self._is_built else "Not built"
        return (f"PRISMSystem(protein='{os.path.basename(self.protein_path)}', "
                f"ligand='{os.path.basename(self.ligand_path)}', "
                f"ff='{self.ligand_forcefield}', status='{status}')")
    
    def info(self):
        """Print detailed information about the system"""
        print(f"\nPRISM System Information:")
        print(f"  Protein file: {self.protein_path}")
        print(f"  Ligand file: {self.ligand_path}")
        print(f"  Output directory: {self.output_dir}")
        print(f"  Ligand force field: {self.ligand_forcefield}")
        print(f"  Build status: {'Complete' if self._is_built else 'Pending'}")
        
        if self._is_built:
            files = self.get_output_files()
            print(f"\nGenerated files:")
            for key, path in files.items():
                exists = "✓" if os.path.exists(path) else "✗"
                print(f"    {exists} {key}: {path}")
        
        if self.config_path:
            print(f"  Configuration: {self.config_path}")