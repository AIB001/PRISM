#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder - Main builder class for protein-ligand systems
"""

import os
import sys
import shutil
import argparse
from pathlib import Path

# Handle both package and direct script execution
if __name__ == "__main__" and __package__ is None:
    # Add parent directory to path for direct execution
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from prism.utils.environment import GromacsEnvironment
    from prism.utils.config import ConfigurationManager
    from prism.utils.mdp import MDPGenerator
    from prism.utils.system import SystemBuilder
    from prism.utils.protonation import ProteinProtonator
    from prism.forcefield.gaff import GAFFForceFieldGenerator
    from prism.forcefield.gaff2 import GAFF2ForceFieldGenerator
    from prism.forcefield.openff import OpenFFForceFieldGenerator
    from prism.forcefield.cgenff import CGenFFForceFieldGenerator
    from prism.forcefield.opls_aa import OPLSAAForceFieldGenerator
    from prism.forcefield.swissparam import MMFFForceFieldGenerator, MATCHForceFieldGenerator, HybridMMFFMATCHForceFieldGenerator
else:
    # Normal package imports
    from .utils.environment import GromacsEnvironment
    from .utils.config import ConfigurationManager
    from .utils.mdp import MDPGenerator
    from .utils.system import SystemBuilder
    from .utils.protonation import ProteinProtonator
    try:
        from .forcefield.gaff import GAFFForceFieldGenerator
        from .forcefield.gaff2 import GAFF2ForceFieldGenerator
        from .forcefield.openff import OpenFFForceFieldGenerator
        from .forcefield.cgenff import CGenFFForceFieldGenerator
        from .forcefield.opls_aa import OPLSAAForceFieldGenerator
        from .forcefield.swissparam import MMFFForceFieldGenerator, MATCHForceFieldGenerator, HybridMMFFMATCHForceFieldGenerator
    except ImportError:
        # For standalone usage
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        try:
            from prism.forcefield.gaff import GAFFForceFieldGenerator
            from prism.forcefield.gaff2 import GAFF2ForceFieldGenerator
            from prism.forcefield.openff import OpenFFForceFieldGenerator
            from prism.forcefield.cgenff import CGenFFForceFieldGenerator
            from prism.forcefield.opls_aa import OPLSAAForceFieldGenerator
            from prism.forcefield.swissparam import MMFFForceFieldGenerator, MATCHForceFieldGenerator, HybridMMFFMATCHForceFieldGenerator
        except ImportError:
            print("Error: Cannot import force field generators")
            print("Please check your PRISM installation")
            sys.exit(1)


class PRISMBuilder:
    """Complete system builder for protein-ligand MD simulations with multiple force field support"""

    def __init__(self, protein_path, ligand_path, output_dir, ligand_forcefield='gaff',
                 config_path=None, forcefield=None, water_model=None, overwrite=None, forcefield_path=None):
        """
        Initialize PRISM Builder with configuration support

        Parameters:
        -----------
        protein_path : str
            Path to the protein PDB file
        ligand_path : str
            Path to the ligand file (MOL2/SDF) - not used for cgenff
        output_dir : str
            Directory where output files will be stored
        ligand_forcefield : str
            Force field for ligand ('gaff', 'gaff2', 'openff', 'cgenff', 'opls', 'mmff', 'match', or 'hybrid')
        config_path : str, optional
            Path to configuration YAML file
        forcefield : str, optional
            Protein force field name (e.g., 'amber99sb', overrides config)
        water_model : str, optional
            Water model name (e.g., 'tip3p', overrides config)
        overwrite : bool, optional
            Whether to overwrite existing files (overrides config)
        forcefield_path : str, optional
            Path to CGenFF directory (required when ligand_forcefield='cgenff')
        """
        self.protein_path = os.path.abspath(protein_path)
        self.ligand_path = os.path.abspath(ligand_path)
        self.output_dir = os.path.abspath(output_dir)
        self.ligand_forcefield = ligand_forcefield.lower()
        self.forcefield_path = forcefield_path

        # Validate ligand force field
        if self.ligand_forcefield not in ['gaff', 'gaff2', 'openff', 'cgenff', 'opls', 'mmff', 'match', 'hybrid']:
            raise ValueError(f"Unsupported ligand force field: {self.ligand_forcefield}. Use 'gaff', 'gaff2', 'openff', 'cgenff', 'opls', 'mmff', 'match', or 'hybrid'")

        # Validate cgenff requires forcefield_path
        if self.ligand_forcefield == 'cgenff' and not forcefield_path:
            raise ValueError("CGenFF force field requires --forcefield-path to specify the directory containing downloaded CGenFF files")

        # Initialize GROMACS environment
        self.gromacs_env = GromacsEnvironment()

        # Process force field names before initializing config
        default_ff = forcefield if forcefield else "amber99sb"
        default_water = water_model if water_model else "tip3p"

        # Initialize configuration manager
        self.config_manager = ConfigurationManager(
            config_path,
            self.gromacs_env,
            forcefield_name=default_ff,
            water_model_name=default_water
        )
        self.config = self.config_manager.config

        # Add ligand force field to config
        self.config['ligand_forcefield'] = {
            'type': self.ligand_forcefield,
            'charge': self.config.get('simulation', {}).get('ligand_charge', 0)
        }

        # Override config with explicit parameters if provided
        if forcefield is not None:
            self.config_manager.update_forcefield_by_name(forcefield)
        if water_model is not None:
            self.config_manager.update_water_model_by_name(water_model)
        if overwrite is not None:
            self.config['general']['overwrite'] = overwrite

        # Extract configuration values
        self.overwrite = self.config['general']['overwrite']
        self.forcefield_idx = self.config['forcefield']['index']
        self.water_model_idx = self.config['water_model']['index']

        # Get force field and water model info
        self.forcefield = self._get_forcefield_info()
        self.water_model = self._get_water_model_info()

        # Extract names
        self.protein_name = Path(protein_path).stem
        self.ligand_name = Path(ligand_path).stem

        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

        # Subdirectories
        self.lig_ff_dir = None

        # Initialize sub-components
        self.mdp_generator = MDPGenerator(self.config, self.output_dir)
        self.system_builder = SystemBuilder(self.config, self.output_dir, self.overwrite)

        self._print_initialization_info()

    def _print_initialization_info(self):
        """Print initialization information"""
        print(f"\nInitialized PRISM Builder:")
        print(f"  GROMACS command: {self.gromacs_env.gmx_command}")
        print(f"  Protein: {self.protein_path}")
        print(f"  Ligand: {self.ligand_path}")
        print(f"  Ligand force field: {self.ligand_forcefield.upper()}")
        print(f"  Output directory: {self.output_dir}")
        print(f"  Protein force field: {self.forcefield['name']}")
        print(f"  Water model: {self.water_model['name']}")
        print(f"  Box distance: {self.config['box']['distance']} nm")
        print(f"  Temperature: {self.config['simulation']['temperature']} K")
        print(f"  pH: {self.config['simulation']['pH']}")
        print(f"  Production time: {self.config['simulation']['production_time_ns']} ns")

    def _get_forcefield_info(self):
        """Get force field information from config"""
        ff_idx = self.config['forcefield']['index']
        custom_ff = self.config['forcefield']['custom_forcefields']

        if ff_idx in custom_ff:
            return custom_ff[ff_idx]
        else:
            raise ValueError(f"Force field index {ff_idx} not found in configuration")

    def _get_water_model_info(self):
        """Get water model information from config"""
        wm_idx = self.config['water_model']['index']
        custom_wm = self.config['water_model']['custom_water_models']

        if wm_idx in custom_wm:
            return custom_wm[wm_idx]
        else:
            raise ValueError(f"Water model index {wm_idx} not found in configuration")

    def generate_ligand_forcefield(self):
        """Generate ligand force field using selected force field generator"""
        print(f"\n=== Generating Ligand Force Field ({self.ligand_forcefield.upper()}) ===")

        if self.ligand_forcefield == 'gaff':
            # Use GAFF force field generator
            generator = GAFFForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()

        elif self.ligand_forcefield == 'gaff2':
            # Use GAFF2 force field generator
            generator = GAFF2ForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()

        elif self.ligand_forcefield == 'openff':
            # Use OpenFF force field generator
            generator = OpenFFForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                charge=self.config['ligand_forcefield']['charge'],
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()

        elif self.ligand_forcefield == 'cgenff':
            # Use CGenFF force field generator
            generator = CGenFFForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                cgenff_dir=self.forcefield_path,
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()

        elif self.ligand_forcefield == 'opls':
            # Use OPLS-AA force field generator (via LigParGen)
            generator = OPLSAAForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                charge=self.config['ligand_forcefield']['charge'],
                charge_model='cm1a',  # or 'cm5' for more accuracy
                align_to_input=True,
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()

        elif self.ligand_forcefield == 'mmff':
            # Use MMFF force field generator (via SwissParam)
            generator = MMFFForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()

        elif self.ligand_forcefield == 'match':
            # Use MATCH force field generator (via SwissParam)
            generator = MATCHForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()

        elif self.ligand_forcefield == 'hybrid':
            # Use hybrid MMFF-based-MATCH force field generator (via SwissParam)
            generator = HybridMMFFMATCHForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()

        # Verify required files exist
        required_files = ["LIG.itp", "LIG.gro", "atomtypes_LIG.itp"]
        for filename in required_files:
            filepath = os.path.join(self.lig_ff_dir, filename)
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"Required file not found: {filepath}")

        print(f"Ligand force field files generated in: {self.lig_ff_dir}")
        return self.lig_ff_dir

    def clean_protein(self, ion_mode=None, distance_cutoff=None, keep_crystal_water=None, remove_artifacts=None):
        """
        Clean the protein PDB file with intelligent metal ion handling.

        Parameters
        ----------
        ion_mode : str, optional
            Ion handling mode. Options:
            - 'keep_all': Keep all metal ions (except water unless keep_crystal_water=True)
            - 'smart' (default): Keep structural metals (Zn, Mg, Ca, Fe, etc.), remove non-structural ions (Na, Cl, etc.)
            - 'remove_all': Remove all metal ions
            If not specified, reads from config or defaults to 'smart'
        distance_cutoff : float, optional
            Maximum distance (Angstroms) from protein for keeping metals.
            Metals farther than this will be removed even if they are structural.
            If not specified, reads from config or defaults to 5.0 Å
        keep_crystal_water : bool, optional
            Whether to keep crystal water molecules. Default: False
            If True, water molecules from the crystal structure are preserved.
            If not specified, reads from config or defaults to False
        remove_artifacts : bool, optional
            Whether to remove crystallization artifacts (GOL, EDO, PEG, NAG, etc.). Default: True
            If False, crystallization artifacts are kept in the output.
            If not specified, reads from config or defaults to True

        Returns
        -------
        str
            Path to cleaned (and optionally protonated) PDB file
        """
        from .utils.cleaner import ProteinCleaner

        print("\n=== Cleaning Protein ===")

        cleaned_pdb = os.path.join(self.output_dir, f"{self.protein_name}_clean.pdb")
        final_pdb = cleaned_pdb  # By default, return the cleaned PDB

        # Check if protonation optimization is requested
        optimize_protonation = self.config.get('protonation', {}).get('optimize', False)

        # If protonation is enabled, we'll create a separate protonated file
        if optimize_protonation:
            protonated_pdb = os.path.join(self.output_dir, f"{self.protein_name}_protonated.pdb")
            # Check if both cleaned and protonated files exist
            if os.path.exists(cleaned_pdb) and os.path.exists(protonated_pdb) and not self.overwrite:
                print(f"Using existing cleaned and protonated protein: {protonated_pdb}")
                return protonated_pdb
            final_pdb = protonated_pdb
        else:
            # No protonation - just check cleaned file
            if os.path.exists(cleaned_pdb) and not self.overwrite:
                print(f"Using existing cleaned protein: {cleaned_pdb}")
                return cleaned_pdb

        # Get parameters from config if not explicitly provided
        if ion_mode is None:
            ion_mode = self.config.get('protein_preparation', {}).get('ion_handling_mode', 'smart')
        if distance_cutoff is None:
            distance_cutoff = self.config.get('protein_preparation', {}).get('metal_distance_cutoff', 5.0)
        if keep_crystal_water is None:
            keep_crystal_water = self.config.get('protein_preparation', {}).get('keep_crystal_water', False)
        if remove_artifacts is None:
            remove_artifacts = self.config.get('protein_preparation', {}).get('remove_crystallization_artifacts', True)

        # Get custom metals from config
        keep_custom = self.config.get('protein_preparation', {}).get('keep_custom_metals', [])
        remove_custom = self.config.get('protein_preparation', {}).get('remove_custom_metals', [])

        # Initialize cleaner
        cleaner = ProteinCleaner(
            ion_mode=ion_mode,
            distance_cutoff=distance_cutoff,
            keep_crystal_water=keep_crystal_water,
            remove_artifacts=remove_artifacts,
            keep_custom_metals=keep_custom if keep_custom else None,
            remove_custom_metals=remove_custom if remove_custom else None,
            verbose=True
        )

        # Clean the protein
        cleaner.clean_pdb(self.protein_path, cleaned_pdb)

        # Post-process: Fix terminal hydrogen names for AMBER compatibility
        self._fix_hydrogen_names(cleaned_pdb)

        print(f"Protein cleaned and saved to: {cleaned_pdb}")

        # === NEW: Protonation optimization step ===
        if optimize_protonation:
            print("\n=== Optimizing Protein Protonation ===")

            # Get protonation config
            prot_config = self.config.get('protonation', {})
            ph = prot_config.get('ph', self.config.get('simulation', {}).get('pH', 7.0))
            his_state = prot_config.get('his_state', 'auto')
            preserve_h = prot_config.get('preserve_existing_h', False)

            try:
                # Initialize protonator
                protonator = ProteinProtonator(
                    ph=ph,
                    preserve_existing_h=preserve_h,
                    his_state=his_state,
                    verbose=True
                )

                # Run protonation optimization
                work_dir = os.path.join(self.output_dir, "protonation_work")
                os.makedirs(work_dir, exist_ok=True)

                results = protonator.optimize_hydrogens(
                    input_pdb=cleaned_pdb,
                    output_pdb=protonated_pdb,
                    work_dir=work_dir
                )

                print(f"Protonation optimization completed")
                print(f"Protonated protein saved to: {protonated_pdb}")

                # Validate the output
                validation = results.get('validation', {})
                if validation:
                    print(f"  Validation: {validation.get('residues_with_h', 0)} hydrogen atoms added")
                    if validation.get('his_residues'):
                        print(f"  Found {len(validation['his_residues'])} histidine residue(s)")

                return protonated_pdb

            except ImportError as e:
                print(f"\nWarning: Meeko not available for protonation optimization")
                print(f"  Error: {e}")
                print(f"  Install with: pip install meeko")
                print(f"  Continuing with cleaned protein without protonation optimization")
                return cleaned_pdb
            except Exception as e:
                print(f"\nWarning: Protonation optimization failed: {e}")
                print(f"  Continuing with cleaned protein without protonation optimization")
                return cleaned_pdb

        return final_pdb

    def _fix_hydrogen_names(self, pdb_file):
        """
        Fix terminal hydrogen names for AMBER compatibility.

        AMBER force fields expect specific naming for N-terminal hydrogens:
        H1, H2, H3 instead of HN1, HN2, HN3
        """
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        fixed_lines = []
        for line in lines:
            if line.startswith('ATOM'):
                # Fix terminal hydrogen names
                if 'HN1' in line:
                    line = line.replace('HN1', 'H1 ')
                elif 'HN2' in line:
                    line = line.replace('HN2', 'H2 ')
                elif 'HN3' in line:
                    line = line.replace('HN3', 'H3 ')
            fixed_lines.append(line)

        with open(pdb_file, 'w') as f:
            f.writelines(fixed_lines)

    def build_model(self, cleaned_protein):
        """Build the GROMACS model"""
        return self.system_builder.build(
            cleaned_protein,
            self.lig_ff_dir,
            self.forcefield_idx,
            self.water_model_idx,
            self.forcefield,  # Pass full force field info (name, dir, path)
            self.water_model  # Pass full water model info
        )

    def generate_mdp_files(self):
        """Generate MDP files for MD simulations"""
        self.mdp_generator.generate_all()

    def cleanup(self):
        """Clean up temporary files"""
        print("\n=== Cleaning up temporary files ===")

        # Cleanup directories for GAFF, GAFF2, and OpenFF
        cleanup_dirs = ["forcefield", "temp_openff"]

        for dir_name in cleanup_dirs:
            cleanup_dir = os.path.join(self.output_dir, dir_name)
            if os.path.exists(cleanup_dir):
                if self.ligand_forcefield in ['gaff', 'gaff2']:
                    temp_patterns = ["*.frcmod", "*.prep", "*.prmtop", "*.rst7", "*.log", "*.in",
                                     "ANTECHAMBER*", "ATOMTYPE*", "PREP*", "NEWPDB*", "sqm*", "leap*"]

                    for pattern in temp_patterns:
                        for file_path in Path(cleanup_dir).glob(pattern):
                            try:
                                os.remove(file_path)
                            except:
                                pass
                else:
                    # For OpenFF, we might want to keep the directory clean
                    # or remove it entirely if it's temporary
                    try:
                        shutil.rmtree(cleanup_dir)
                    except:
                        pass

        print("Cleanup completed")

    def save_config(self):
        """Save the current configuration to a file"""
        config_file = os.path.join(self.output_dir, "prism_config.yaml")
        self.config_manager.save_config(config_file)
        print(f"Configuration saved to: {config_file}")

    def generate_localrun_script(self):
        """Generate localrun.sh script for easy MD execution"""
        gmx_md_dir = os.path.join(self.output_dir, 'GMX_PROLIG_MD')
        if not os.path.exists(gmx_md_dir):
            print("Warning: GMX_PROLIG_MD directory not found, skipping localrun.sh generation")
            return None

        # Clean up Emacs backup files (#topol.top.1#, #topol.top.2#, etc.)
        import glob
        backup_pattern = os.path.join(gmx_md_dir, '#*#')
        backup_files = glob.glob(backup_pattern)
        if backup_files:
            print(f"Cleaning up {len(backup_files)} Emacs backup file(s)...")
            for backup_file in backup_files:
                try:
                    os.remove(backup_file)
                    print(f"  Removed: {os.path.basename(backup_file)}")
                except Exception as e:
                    print(f"  Warning: Could not remove {backup_file}: {e}")

        script_path = os.path.join(gmx_md_dir, 'localrun.sh')

        script_content = '''#!/bin/bash

######################################################
# SIMULATION PART
######################################################

# Energy Minimization (EM)
mkdir -p em
if [ -f ./em/em.gro ]; then
    echo "EM already completed, skipping..."
elif [ -f ./em/em.tpr ]; then
    echo "EM tpr file found, continuing from checkpoint..."
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v -cpi ./em/em.cpt
else
    echo "Starting EM from scratch..."
    gmx grompp -f ../mdps/em.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o ./em/em.tpr -maxwarn 999
    gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v
fi

# NVT Equilibration
mkdir -p nvt
if [ -f ./nvt/nvt.gro ]; then
    echo "NVT already completed, skipping..."
elif [ -f ./nvt/nvt.tpr ]; then
    echo "NVT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v -cpi ./nvt/nvt.cpt
else
    echo "Starting NVT from scratch..."
    gmx grompp -f ../mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro -p topol.top -o ./nvt/nvt.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v
fi

# NPT Equilibration
mkdir -p npt
if [ -f ./npt/npt.gro ]; then
    echo "NPT already completed, skipping..."
elif [ -f ./npt/npt.tpr ]; then
    echo "NPT tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v -cpi ./npt/npt.cpt
else
    echo "Starting NPT from scratch..."
    gmx grompp -f ../mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v
fi

# Production MD
mkdir -p prod
if [ -f ./prod/md.gro ]; then
    echo "Production MD already completed, skipping..."
elif [ -f ./prod/md.tpr ]; then
    echo "Production MD tpr file found, continuing from checkpoint..."
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v -cpi ./prod/md.cpt
else
    echo "Starting Production MD from scratch..."
    gmx grompp -f ../mdps/md.mdp -c ./npt/npt.gro -r ./npt/npt.gro -p topol.top -o ./prod/md.tpr -maxwarn 999
    gmx mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v
fi
'''

        with open(script_path, 'w') as f:
            f.write(script_content)

        # Make script executable
        os.chmod(script_path, 0o755)

        print(f"Local run script generated: {script_path}")
        return script_path

    def run(self):
        """Run the complete workflow"""
        print(f"\n{'='*60}")
        print("Starting PRISM Builder workflow")
        print(f"{'='*60}")

        try:
            # Save configuration for reference
            self.save_config()

            # Step 1: Generate ligand force field
            self.generate_ligand_forcefield()

            # Step 2: Clean protein
            cleaned_protein = self.clean_protein()

            # Step 3: Build model
            model_dir = self.build_model(cleaned_protein)
            if not model_dir:
                raise RuntimeError("Failed to build model")

            # Step 4: Generate MDP files
            self.generate_mdp_files()

            # Step 5: Cleanup
            self.cleanup()

            # Step 6: Generate local run script
            script_path = self.generate_localrun_script()

            print(f"\n{'='*60}")
            print("PRISM Builder workflow completed successfully!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {self.output_dir}")
            print(f"MD system files are in: {os.path.join(self.output_dir, 'GMX_PROLIG_MD')}")
            print(f"MDP files are in: {self.mdp_generator.mdp_dir}")
            print(f"Configuration saved in: {os.path.join(self.output_dir, 'prism_config.yaml')}")
            print(f"\nProtein force field used: {self.forcefield['name']}")
            print(f"Ligand force field used: {self.ligand_forcefield.upper()}")
            print(f"Water model used: {self.water_model['name']}")

            if script_path:
                gmx_md_dir = os.path.join(self.output_dir, 'GMX_PROLIG_MD')
                print(f"\n{'='*60}")
                print("Ready to run MD simulations!")
                print(f"{'='*60}")
                print(f"\nTo run the MD workflow:")
                print(f"  1. Navigate to the MD directory:")
                print(f"     cd {gmx_md_dir}")
                print(f"  2. Execute the simulation script:")
                print(f"     bash localrun.sh")
                print(f"\nThe script will run:")
                print(f"  - Energy Minimization (EM)")
                print(f"  - NVT Equilibration")
                print(f"  - NPT Equilibration")
                print(f"  - Production MD")
                print(f"\nNote: Adjust GPU and thread settings in localrun.sh as needed")

            return self.output_dir

        except Exception as e:
            print(f"\nError during PRISM Builder workflow: {e}")
            import traceback
            traceback.print_exc()
            raise


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="PRISM Builder - Build protein-ligand systems for GROMACS with multiple force field support",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Using GAFF force field with defaults (amber99sb, tip3p)
  prism protein.pdb ligand.mol2 -o output_dir

  # Using GAFF2 force field (improved version of GAFF)
  prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield gaff2

  # Using OpenFF force field with specific protein force field
  prism protein.pdb ligand.sdf -o output_dir --ligand-forcefield openff --forcefield amber14sb

  # Using OPLS-AA force field (via LigParGen server, requires internet)
  prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield opls

  # Using CGenFF force field (requires web-downloaded files)
  prism protein.pdb dummy.mol2 -o output_dir --ligand-forcefield cgenff --forcefield-path /path/to/cgenff_files

  # Using SwissParam force fields (via SwissParam server, requires internet)
  prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield mmff    # MMFF-based
  prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield match   # MATCH
  prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield hybrid  # Hybrid MMFF-based-MATCH

  # With custom configuration file
  prism protein.pdb ligand.mol2 -o output_dir --config config.yaml

  # Override specific parameters
  prism protein.pdb ligand.mol2 -o output_dir --forcefield amber99sb-ildn --water tip4p --temperature 300

  # Export default configuration template
  prism --export-config my_config.yaml

  # List available force fields
  prism --list-forcefields
        """
    )

    # Positional arguments (allow mixed order with options)
    parser.add_argument("protein", nargs='?', help="Path to protein PDB file")
    parser.add_argument("ligand", nargs='?', help="Path to ligand file (MOL2/SDF)")
    
    # Alternative way to specify input files
    parser.add_argument("--protein-file", "-pf", help="Path to protein PDB file (alternative to positional)")
    parser.add_argument("--ligand-file", "-lf", help="Path to ligand file (alternative to positional)")

    # Basic options
    basic = parser.add_argument_group('Basic options')
    basic.add_argument("--output", "-o", default="prism_output",
                       help="Output directory (default: prism_output)")
    basic.add_argument("--config", "-c",
                       help="Path to configuration YAML file")
    basic.add_argument("--overwrite", "-f", action="store_true",
                       help="Overwrite existing files")

    # Force field options
    ff_group = parser.add_argument_group('Force field options')
    ff_group.add_argument("--forcefield", "-ff", type=str, default="amber99sb",
                          help="Protein force field name (default: amber99sb)")
    ff_group.add_argument("--water", "-w", type=str, default="tip3p",
                          help="Water model name (default: tip3p)")
    ff_group.add_argument("--ligand-forcefield", "-lff", choices=['gaff', 'gaff2', 'openff', 'cgenff', 'opls', 'mmff', 'match', 'hybrid'], default='gaff',
                          help="Force field for ligand: gaff (default), gaff2 (improved), openff, cgenff, opls (OPLS-AA via LigParGen), mmff/match/hybrid (SwissParam)")
    ff_group.add_argument("--ligand-charge", type=int, default=0,
                          help="Net charge of ligand (default: 0)")
    ff_group.add_argument("--forcefield-path", "-ffp", type=str, default=None,
                          help="Path to force field directory (required for cgenff: directory containing web-downloaded CGenFF files)")

    # Box options
    box_group = parser.add_argument_group('Box options')
    box_group.add_argument("--box-distance", "-d", type=float, default=1.5,
                           help="Distance from protein to box edge in nm (default: 1.5)")
    box_group.add_argument("--box-shape", "-bs", choices=['cubic', 'dodecahedron', 'octahedron'], default='cubic',
                           help="Box shape (default: cubic)")
    box_group.add_argument("--no-center", action="store_true",
                           help="Don't center protein in box")

    # Simulation parameters
    sim_group = parser.add_argument_group('Simulation parameters')
    sim_group.add_argument("--temperature", "-t", type=float, default=310,
                           help="Temperature in K (default: 310)")
    sim_group.add_argument("--pressure", "-p", type=float, default=1.0,
                           help="Pressure in bar (default: 1.0)")
    sim_group.add_argument("--pH", type=float, default=7.0,
                           help="pH for protonation states (default: 7.0)")
    sim_group.add_argument("--production-ns", type=float, default=500,
                           help="Production time in ns (default: 500)")
    sim_group.add_argument("--dt", type=float, default=0.002,
                           help="Time step in ps (default: 0.002)")
    sim_group.add_argument("--nvt-ps", type=float, default=500,
                           help="NVT equilibration time in ps (default: 500)")
    sim_group.add_argument("--npt-ps", type=float, default=500,
                           help="NPT equilibration time in ps (default: 500)")

    # Ion options
    ion_group = parser.add_argument_group('Ion options')
    ion_group.add_argument("--no-neutralize", action="store_true",
                           help="Don't neutralize the system")
    ion_group.add_argument("--salt-concentration", "-sc", type=float, default=0.15,
                           help="Salt concentration in M (default: 0.15)")
    ion_group.add_argument("--positive-ion", "-pion", default="NA",
                           help="Positive ion type (default: NA)")
    ion_group.add_argument("--negative-ion", "-nion", default="CL",
                           help="Negative ion type (default: CL)")

    # Energy minimization
    em_group = parser.add_argument_group('Energy minimization')
    em_group.add_argument("--em-tolerance", type=float, default=200.0,
                          help="Energy minimization tolerance in kJ/mol/nm (default: 200.0)")
    em_group.add_argument("--em-steps", type=int, default=10000,
                          help="Maximum EM steps (default: 10000)")

    # Output options
    out_group = parser.add_argument_group('Output options')
    out_group.add_argument("--traj-interval", type=float, default=500,
                           help="Trajectory output interval in ps (default: 500)")
    out_group.add_argument("--energy-interval", type=float, default=10,
                           help="Energy output interval in ps (default: 10)")
    out_group.add_argument("--no-compressed", action="store_true",
                           help="Don't use compressed trajectory format")

    # Utility options
    util_group = parser.add_argument_group('Utility options')
    util_group.add_argument("--list-forcefields", action="store_true",
                            help="List available force fields and exit")
    util_group.add_argument("--install-forcefields", action="store_true",
                            help="Install additional force fields from PRISM configs to GROMACS")
    util_group.add_argument("--export-config", metavar="FILE",
                            help="Export default configuration to file and exit")
    util_group.add_argument("--gmx-command", default=None,
                            help="GROMACS command to use (auto-detected if not specified)")

    # MM/PBSA options
    mmpbsa_group = parser.add_argument_group('MM/PBSA options')
    mmpbsa_group.add_argument("--mmpbsa", "-pbsa", action="store_true",
                              help="Run MM/PBSA calculation instead of full MD workflow")
    mmpbsa_group.add_argument("--mode", "-m", choices=['single-frame', 'trajectory'], default='single-frame',
                              help="MM/PBSA mode: single-frame (docking pose) or trajectory (MD)")
    mmpbsa_group.add_argument("--structure", "-s", help="Path to structure file (GRO) for single-frame MM/PBSA")
    mmpbsa_group.add_argument("--system-dir", "-sd", help="Path to PRISM-generated GMX_PROLIG_MD directory")
    mmpbsa_group.add_argument("--mmpbsa-output", "-po", help="Output directory for MM/PBSA results (default: same as structure file directory)")

    args = parser.parse_args()

    # Handle MM/PBSA workflow (takes priority over normal workflow)
    if args.mmpbsa:
        if args.mode == 'single-frame':
            # Import MM/PBSA module
            try:
                if __name__ == "__main__" and __package__ is None:
                    from prism.mmpbsa import SingleFrameMMPBSA
                else:
                    from .mmpbsa import SingleFrameMMPBSA
            except ImportError as e:
                print(f"Error importing MM/PBSA module: {e}")
                print("Please ensure gmx_MMPBSA is installed: pip install gmx_MMPBSA")
                sys.exit(1)

            # Validate required arguments
            if not args.structure:
                parser.error("--structure is required for single-frame MM/PBSA")
            if not args.system_dir:
                parser.error("--system-dir is required for single-frame MM/PBSA")

            # Run single-frame MM/PBSA
            try:
                calculator = SingleFrameMMPBSA(
                    structure_file=args.structure,
                    system_dir=args.system_dir,
                    output_dir=args.mmpbsa_output,  # None defaults to structure file directory
                    overwrite=args.overwrite
                )
                calculator.run()
                sys.exit(0)
            except Exception as e:
                print(f"\nMM/PBSA calculation failed: {e}")
                import traceback
                traceback.print_exc()
                sys.exit(1)
        else:
            # Trajectory mode - to be implemented later
            print("Trajectory MM/PBSA mode will be implemented in the future")
            sys.exit(1)

    # Handle --export-config option
    if args.export_config:
        import pkg_resources
        try:
            # Try to get from package resources
            default_config = pkg_resources.resource_string('prism', 'configs/default_config.yaml').decode('utf-8')
        except:
            # Fallback to embedded string
            default_config = """# PRISM Default Configuration File
# This file serves as a template for creating custom configurations
# Copy this file and modify as needed for your simulations

general:
  overwrite: false
  # gmx_command is auto-detected

box:
  distance: 1.5
  shape: cubic
  center: true

simulation:
  temperature: 310
  pressure: 1.0
  pH: 7.0
  ligand_charge: 0
  production_time_ns: 500
  dt: 0.002
  equilibration_nvt_time_ps: 500
  equilibration_npt_time_ps: 500

ions:
  neutral: true
  concentration: 0.15
  positive_ion: NA
  negative_ion: CL

constraints:
  algorithm: lincs
  type: h-bonds
  lincs_iter: 1
  lincs_order: 4

energy_minimization:
  integrator: steep
  emtol: 200.0
  emstep: 0.01
  nsteps: 10000

output:
  trajectory_interval_ps: 500
  energy_interval_ps: 10
  log_interval_ps: 10
  compressed_trajectory: true

electrostatics:
  coulombtype: PME
  rcoulomb: 1.0
  pme_order: 4
  fourierspacing: 0.16

vdw:
  rvdw: 1.0
  dispcorr: EnerPres

temperature_coupling:
  tcoupl: V-rescale
  tc_grps:
    - Protein
    - Non-Protein
  tau_t:
    - 0.1
    - 0.1

pressure_coupling:
  pcoupl: C-rescale
  pcoupltype: isotropic
  tau_p: 1.0
  compressibility: 4.5e-05

# Note: Force field and water model are specified via command line
"""
        with open(args.export_config, 'w') as f:
            f.write(default_config)
        print(f"Default configuration exported to: {args.export_config}")
        sys.exit(0)

    # Handle --install-forcefields option
    if args.install_forcefields:
        try:
            # Import force field installer
            if __name__ == "__main__" and __package__ is None:
                from prism.utils.forcefield_installer import ForceFieldInstaller
            else:
                from .utils.forcefield_installer import ForceFieldInstaller

            installer = ForceFieldInstaller()
            installer.install_forcefields_interactive()
        except Exception as e:
            print(f"Error installing force fields: {e}")
            import traceback
            traceback.print_exc()
        sys.exit(0)

    # Handle --list-forcefields option
    if args.list_forcefields:
        try:
            # Import here to avoid circular imports
            if __name__ == "__main__" and __package__ is None:
                from prism.utils.environment import GromacsEnvironment
            else:
                from .utils.environment import GromacsEnvironment

            env = GromacsEnvironment()
            print("\nAvailable force fields:")
            for ff in env.list_force_fields():
                print(f"  - {ff}")

            # Show water models for default force field
            default_ff_idx = env.get_force_field_index("amber99sb")
            if default_ff_idx:
                print(f"\nWater models for amber99sb:")
                for wm in env.list_water_models(default_ff_idx):
                    print(f"  - {wm}")

            print("\nNote: Water models vary by force field. Specify a force field to see its water models.")
        except Exception as e:
            print(f"Error detecting force fields: {e}")
        sys.exit(0)

    # Check required arguments - handle alternative input methods
    protein_path = args.protein or args.protein_file
    ligand_path = args.ligand or args.ligand_file
    
    if not protein_path or not ligand_path:
        parser.error("Both protein and ligand files are required. Use positional arguments or --protein-file/--ligand-file flags")
    
    # Update args to use resolved paths
    args.protein = protein_path
    args.ligand = ligand_path

    # Build kwargs from command-line arguments
    kwargs = {}

    # Process command-line overrides for config
    config_overrides = {}

    # General overrides
    if args.gmx_command:
        config_overrides.setdefault('general', {})['gmx_command'] = args.gmx_command

    # Box overrides
    if args.box_distance != 1.5:
        config_overrides.setdefault('box', {})['distance'] = args.box_distance
    if args.box_shape != 'cubic':
        config_overrides.setdefault('box', {})['shape'] = args.box_shape
    if args.no_center:
        config_overrides.setdefault('box', {})['center'] = False

    # Simulation overrides
    if args.temperature != 310:
        config_overrides.setdefault('simulation', {})['temperature'] = args.temperature
    if args.pressure != 1.0:
        config_overrides.setdefault('simulation', {})['pressure'] = args.pressure
    if args.pH != 7.0:
        config_overrides.setdefault('simulation', {})['pH'] = args.pH
    if args.production_ns != 500:
        config_overrides.setdefault('simulation', {})['production_time_ns'] = args.production_ns
    if args.dt != 0.002:
        config_overrides.setdefault('simulation', {})['dt'] = args.dt
    if args.nvt_ps != 500:
        config_overrides.setdefault('simulation', {})['equilibration_nvt_time_ps'] = args.nvt_ps
    if args.npt_ps != 500:
        config_overrides.setdefault('simulation', {})['equilibration_npt_time_ps'] = args.npt_ps
    if args.ligand_charge != 0:
        config_overrides.setdefault('simulation', {})['ligand_charge'] = args.ligand_charge

    # Ion overrides
    if args.no_neutralize:
        config_overrides.setdefault('ions', {})['neutral'] = False
    if args.salt_concentration != 0.15:
        config_overrides.setdefault('ions', {})['concentration'] = args.salt_concentration
    if args.positive_ion != 'NA':
        config_overrides.setdefault('ions', {})['positive_ion'] = args.positive_ion
    if args.negative_ion != 'CL':
        config_overrides.setdefault('ions', {})['negative_ion'] = args.negative_ion

    # Energy minimization overrides
    if args.em_tolerance != 200.0:
        config_overrides.setdefault('energy_minimization', {})['emtol'] = args.em_tolerance
    if args.em_steps != 10000:
        config_overrides.setdefault('energy_minimization', {})['nsteps'] = args.em_steps

    # Output overrides
    if args.traj_interval != 500:
        config_overrides.setdefault('output', {})['trajectory_interval_ps'] = args.traj_interval
    if args.energy_interval != 10:
        config_overrides.setdefault('output', {})['energy_interval_ps'] = args.energy_interval
    if args.no_compressed:
        config_overrides.setdefault('output', {})['compressed_trajectory'] = False

    # Save config overrides if any
    if config_overrides:
        # Create temporary config file with overrides
        import tempfile
        import yaml

        # If user provided a config, load it first
        if args.config:
            with open(args.config, 'r') as f:
                base_config = yaml.safe_load(f)
        else:
            base_config = {}

        # Merge overrides
        for key, value in config_overrides.items():
            if key in base_config:
                base_config[key].update(value)
            else:
                base_config[key] = value

        # Write temporary config
        temp_config = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
        yaml.dump(base_config, temp_config, default_flow_style=False)
        temp_config.close()
        kwargs['config_path'] = temp_config.name
    elif args.config:
        kwargs['config_path'] = args.config

    # Create and run PRISM Builder
    builder = PRISMBuilder(
        args.protein,
        args.ligand,
        args.output,
        ligand_forcefield=args.ligand_forcefield,
        forcefield=args.forcefield,
        water_model=args.water,
        overwrite=args.overwrite,
        forcefield_path=args.forcefield_path,
        **kwargs
    )

    try:
        builder.run()
    finally:
        # Clean up temporary config if created
        if config_overrides and 'config_path' in kwargs:
            try:
                os.unlink(kwargs['config_path'])
            except:
                pass


if __name__ == "__main__":
    # Handle direct script execution
    if __package__ is None:
        # Fix imports for direct execution
        import sys
        import os
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    main()