#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Conformer generation for ligands and flexible molecules
"""

import os
import logging
from .base import ModelingBase, ModelingError, ValidationError

logger = logging.getLogger(__name__)


class ConformerGenerator(ModelingBase):
    """
    Generate and optimize molecular conformers for ligands.
    
    Supports multiple conformer generation methods and optimization protocols.
    """
    
    def __init__(self, smiles, num_conformers=100, output_dir="conformer_output", **kwargs):
        """
        Initialize conformer generator.
        
        Parameters:
        -----------
        smiles : str
            SMILES string for molecule
        num_conformers : int
            Number of conformers to generate
        output_dir : str
            Directory for output files
        **kwargs : optional
            Generation parameters (method, optimization, etc.)
        """
        super().__init__(output_dir, **kwargs)
        
        self.smiles = smiles
        self.num_conformers = num_conformers
        
        # Generation parameters
        self.method = kwargs.get('method', 'rdkit')
        self.optimization = kwargs.get('optimization', 'mmff')
        self.energy_window = kwargs.get('energy_window', 10.0)  # kcal/mol
        self.rmsd_threshold = kwargs.get('rmsd_threshold', 0.5)  # Angstrom
        
        # Validate inputs
        if not self.validate_inputs():
            raise ValidationError("Input validation failed")
        
        logger.info(f"ConformerGenerator initialized for SMILES: {self.smiles}")
    
    def validate_inputs(self):
        """
        Validate inputs for conformer generation.
        
        Returns:
        --------
        bool
            True if inputs are valid
        """
        # Check SMILES validity
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is None:
                logger.error(f"Invalid SMILES: {self.smiles}")
                return False
        except ImportError:
            logger.warning("RDKit not available - SMILES validation skipped")
        
        # Check parameter ranges
        if self.num_conformers <= 0:
            logger.error("Number of conformers must be positive")
            return False
        
        if self.energy_window <= 0:
            logger.error("Energy window must be positive")
            return False
        
        logger.info("Input validation passed")
        return True
    
    def prepare_system(self):
        """
        Prepare the conformer generation system.
        
        Returns:
        --------
        dict
            System preparation results
        """
        results = {}
        
        # Generate initial molecule
        mol_info = self._generate_molecule()
        results.update(mol_info)
        
        # Prepare for conformer generation
        prep_info = self._prepare_for_generation()
        results.update(prep_info)
        
        logger.info("Conformer system preparation completed")
        return results
    
    def _generate_molecule(self):
        """Generate RDKit molecule from SMILES"""
        logger.info("Generating molecule from SMILES...")
        
        results = {}
        
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            
            # Create molecule
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is None:
                raise ModelingError(f"Could not create molecule from SMILES: {self.smiles}")
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Calculate molecular properties
            properties = {
                'molecular_weight': Descriptors.MolWt(mol),
                'num_atoms': mol.GetNumAtoms(),
                'num_bonds': mol.GetNumBonds(),
                'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                'logp': Descriptors.MolLogP(mol)
            }
            
            # Save properties
            import yaml
            props_file = os.path.join(self.output_dir, "molecular_properties.yaml")
            with open(props_file, 'w') as f:
                yaml.dump(properties, f, default_flow_style=False)
            
            results['molecule_properties'] = props_file
            results['num_rotatable_bonds'] = properties['num_rotatable_bonds']
            
            self.generated_files['molecular_properties'] = props_file
            
            logger.info(f"Molecule prepared: MW={properties['molecular_weight']:.1f}, "
                       f"rotatable_bonds={properties['num_rotatable_bonds']}")
            
        except ImportError:
            logger.error("RDKit not available for molecule generation")
            raise ModelingError("RDKit required for conformer generation")
        except Exception as e:
            logger.error(f"Molecule generation failed: {e}")
            raise ModelingError(f"Molecule generation failed: {e}")
        
        return results
    
    def _prepare_for_generation(self):
        """Prepare parameters for conformer generation"""
        logger.info("Preparing conformer generation parameters...")
        
        results = {}
        
        # Adjust number of conformers based on flexibility
        rotatable_bonds = results.get('num_rotatable_bonds', 5)
        
        if rotatable_bonds > 10:
            # Highly flexible molecule - may need more conformers
            recommended_conformers = min(self.num_conformers * 2, 500)
            logger.info(f"Highly flexible molecule ({rotatable_bonds} rotatable bonds)")
            logger.info(f"Recommended conformers: {recommended_conformers}")
        elif rotatable_bonds < 3:
            # Rigid molecule - fewer conformers needed
            recommended_conformers = min(self.num_conformers, 50)
            logger.info(f"Rigid molecule ({rotatable_bonds} rotatable bonds)")
        else:
            recommended_conformers = self.num_conformers
        
        results['recommended_conformers'] = recommended_conformers
        
        return results
    
    def generate_conformers(self, prune_rmsd=True):
        """
        Generate conformers for the molecule.
        
        Parameters:
        -----------
        prune_rmsd : bool
            Whether to remove similar conformers based on RMSD
            
        Returns:
        --------
        str
            Path to generated conformers file
        """
        logger.info(f"Generating {self.num_conformers} conformers...")
        
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Create molecule
            mol = Chem.MolFromSmiles(self.smiles)
            mol = Chem.AddHs(mol)
            
            # Generate conformers
            conf_ids = AllChem.EmbedMultipleConfs(
                mol, 
                numConfs=self.num_conformers,
                randomSeed=42,
                pruneRmsThresh=self.rmsd_threshold if prune_rmsd else -1
            )
            
            logger.info(f"Generated {len(conf_ids)} initial conformers")
            
            # Optimize conformers
            if self.optimization == 'mmff':
                for conf_id in conf_ids:
                    AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
            elif self.optimization == 'uff':
                for conf_id in conf_ids:
                    AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
            
            # Calculate energies and filter
            conformer_energies = []
            for conf_id in conf_ids:
                if self.optimization == 'mmff':
                    props = AllChem.MMFFGetMoleculeProperties(mol)
                    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                    if ff is not None:
                        energy = ff.CalcEnergy()
                    else:
                        energy = float('inf')
                else:
                    energy = 0.0  # Placeholder for UFF
                
                conformer_energies.append((conf_id, energy))
            
            # Sort by energy and filter
            conformer_energies.sort(key=lambda x: x[1])
            
            # Keep conformers within energy window
            if conformer_energies:
                min_energy = conformer_energies[0][1]
                filtered_conformers = [
                    (conf_id, energy) for conf_id, energy in conformer_energies
                    if energy - min_energy <= self.energy_window
                ]
            else:
                filtered_conformers = conformer_energies
            
            logger.info(f"Filtered to {len(filtered_conformers)} conformers "
                       f"within {self.energy_window} kcal/mol")
            
            # Write conformers to file
            conformers_sdf = os.path.join(self.output_dir, "conformers.sdf")
            
            writer = Chem.SDWriter(conformers_sdf)
            for i, (conf_id, energy) in enumerate(filtered_conformers):
                mol.SetProp(f"Energy", f"{energy:.2f}")
                mol.SetProp(f"ConformerID", str(conf_id))
                mol.SetProp(f"Rank", str(i + 1))
                writer.write(mol, confId=conf_id)
            writer.close()
            
            # Save energy information
            energy_file = os.path.join(self.output_dir, "conformer_energies.txt")
            with open(energy_file, 'w') as f:
                f.write("ConformerID\tEnergy(kcal/mol)\tRelativeEnergy\tRank\n")
                for i, (conf_id, energy) in enumerate(filtered_conformers):
                    rel_energy = energy - min_energy if filtered_conformers else 0
                    f.write(f"{conf_id}\t{energy:.2f}\t{rel_energy:.2f}\t{i+1}\n")
            
            self.generated_files['conformers'] = conformers_sdf
            self.generated_files['energies'] = energy_file
            
            logger.info(f"Conformers generated: {conformers_sdf}")
            return conformers_sdf
            
        except ImportError:
            logger.error("RDKit not available for conformer generation")
            raise ModelingError("RDKit required for conformer generation")
        except Exception as e:
            logger.error(f"Conformer generation failed: {e}")
            raise ModelingError(f"Conformer generation failed: {e}")
    
    def select_best_conformer(self, criterion="energy"):
        """
        Select the best conformer based on specified criterion.
        
        Parameters:
        -----------
        criterion : str
            Selection criterion ("energy" or "diversity")
            
        Returns:
        --------
        str
            Path to best conformer file
        """
        if 'conformers' not in self.generated_files:
            raise ModelingError("No conformers available. Run generate_conformers() first.")
        
        logger.info(f"Selecting best conformer by {criterion}...")
        
        # Read conformers file and select best
        try:
            from rdkit import Chem
            
            supplier = Chem.SDMolSupplier(self.generated_files['conformers'])
            mols = [mol for mol in supplier if mol is not None]
            
            if not mols:
                raise ModelingError("No valid conformers found")
            
            if criterion == "energy":
                # Select lowest energy conformer (first one after sorting)
                best_mol = mols[0]
            else:
                # For diversity, select a representative conformer
                best_mol = mols[len(mols)//3]  # Middle conformer as representative
            
            # Write best conformer
            best_conformer_file = os.path.join(self.output_dir, "best_conformer.sdf")
            
            writer = Chem.SDWriter(best_conformer_file)
            writer.write(best_mol)
            writer.close()
            
            self.generated_files['best_conformer'] = best_conformer_file
            
            logger.info(f"Best conformer selected: {best_conformer_file}")
            return best_conformer_file
            
        except ImportError:
            logger.error("RDKit not available for conformer selection")
            raise ModelingError("RDKit required for conformer selection")
        except Exception as e:
            logger.error(f"Conformer selection failed: {e}")
            raise ModelingError(f"Conformer selection failed: {e}")