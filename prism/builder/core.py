#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder Core - Main PRISMBuilder class

This module contains the PRISMBuilder class that orchestrates all workflows.
Functionality is split across focused mixins:
- WorkflowMixin: Normal MD, REST2, MM/PBSA workflows
- FEPWorkflowMixin: FEP workflow for relative binding free energy
- LigandForceFieldMixin: Ligand force field generation
- SystemPreparationMixin: Protein cleaning, MDP generation, scripts
- PMFBuilderMixin: PMF workflow (steered MD + umbrella sampling)
"""

import os
from pathlib import Path

from ..utils.environment import GromacsEnvironment
from ..utils.config import ConfigurationManager
from ..utils.mdp import MDPGenerator
from ..utils.system import SystemBuilder
from ..utils.colors import (
    print_warning,
    success,
    path,
    number,
)

from .pmf import PMFBuilderMixin
from .workflow import WorkflowMixin
from .workflow_fep import FEPWorkflowMixin
from .ligand_ff import LigandForceFieldMixin
from .preparation import SystemPreparationMixin


class PRISMBuilder(
    WorkflowMixin,
    FEPWorkflowMixin,
    LigandForceFieldMixin,
    SystemPreparationMixin,
    PMFBuilderMixin,
):
    """Complete system builder for protein-ligand MD simulations with multiple force field support"""

    def __init__(
        self,
        protein_path,
        ligand_paths,
        output_dir,
        ligand_forcefield="gaff",
        config_path=None,
        forcefield=None,
        water_model=None,
        overwrite=None,
        forcefield_path=None,
        nter=None,
        cter=None,
        gaussian_method=None,
        do_optimization=False,
        resp_files=None,
        gaussian_nproc=16,
        gaussian_mem="4GB",
        pmf_mode=False,
        pullvec=None,
        pull_mode="pocket",
        box_extension=None,
        rest2_mode=False,
        t_ref=310.0,
        t_max=450.0,
        n_replicas=16,
        rest2_cutoff=0.5,
        mmpbsa_mode=False,
        mmpbsa_traj_ns=None,
        gmx2amber=False,
        fep_mode=False,
        mutant_ligand=None,
        fep_config=None,
        distance_cutoff=0.6,
        charge_cutoff=0.05,
        charge_reception="surround",
        charge_strategy="mean",
        lambda_windows=11,
        lambda_strategy="decoupled",
        lambda_distribution="nonlinear",
    ):
        """
        Initialize PRISM Builder with configuration support

        Parameters:
        -----------
        protein_path : str
            Path to the protein PDB file
        ligand_paths : str or list
            Path(s) to the ligand file(s) (MOL2/SDF) - str for single ligand, list for multiple
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
        nter : str, optional
            N-terminal type for pdb2gmx (e.g., 'NH3+', 'NH2'). Default: auto-select
        cter : str, optional
            C-terminal type for pdb2gmx (e.g., 'COO-', 'COOH'). Default: auto-select
        gaussian_method : str, optional
            Gaussian method for RESP charges: 'hf' or 'dft'. None disables Gaussian RESP.
        do_optimization : bool, optional
            Whether to perform geometry optimization before ESP (default: False)
        resp_files : str or list, optional
            Path(s) to existing RESP mol2 file(s) to use instead of AM1-BCC charges
        gaussian_nproc : int, optional
            Number of processors for Gaussian calculation (default: 16)
        gaussian_mem : str, optional
            Memory allocation for Gaussian calculation (default: '4GB')
        pmf_mode : bool, optional
            If True, build system for PMF calculations (default: False)
        pullvec : tuple of (int, int), optional
            User-defined pull vector as (protein_atom_index, ligand_atom_index).
            If None, uses pocket centroid -> ligand centroid (default).
        box_extension : tuple of (float, float, float), optional
            Box extension in (X, Y, Z) direction in nm for PMF mode.
            Default: (0.0, 0.0, 2.0) - extends Z by 2 nm for pulling space.
        rest2_mode : bool, optional
            If True, build standard MD system then generate REST2 setup (default: False)
        t_ref : float, optional
            REST2 reference temperature in K (default: 310.0)
        t_max : float, optional
            REST2 maximum effective temperature in K (default: 450.0)
        n_replicas : int, optional
            Number of REST2 replicas (default: 16)
        rest2_cutoff : float, optional
            Pocket detection cutoff in nm for REST2 (default: 0.5)
        mmpbsa_mode : bool, optional
            If True, build standard MD system then set up MM/PBSA calculation (default: False)
        mmpbsa_traj_ns : float, optional
            Production MD length in ns for trajectory-based MM/PBSA. If None (default),
            uses single-frame mode (no production MD). Only used when mmpbsa_mode=True.
        gmx2amber : bool, optional
            If True, use parmed + AMBER MMPBSA.py instead of gmx_MMPBSA (default: False).
            Requires AmberTools with parmed. Only used when mmpbsa_mode=True.
        fep_mode : bool, optional
            If True, enable FEP mode for relative binding free energy calculations (default: False).
            Requires mutant_ligand to be specified.
        mutant_ligand : str, optional
            Path to mutant ligand file for FEP calculations (required when fep_mode=True).
        fep_config : str, optional
            Path to FEP configuration YAML file.
        distance_cutoff : float, optional
            Distance cutoff for atom mapping in Angstroms (default: 0.6).
        charge_strategy : str, optional
            Charge strategy for common atoms: 'ref', 'mut', or 'mean' (default: 'mean').
        lambda_windows : int, optional
            Number of lambda windows for FEP calculations (default: 11).
        lambda_strategy : str, optional
            Lambda schedule strategy: 'coupled', 'decoupled', or 'custom' (default: 'decoupled').
        lambda_distribution : str, optional
            Lambda distribution type: 'linear', 'nonlinear', or 'quadratic' (default: 'nonlinear').
        """
        self.protein_path = os.path.abspath(protein_path)
        self.config_file = config_path

        # Handle both single ligand (string) and multiple ligands (list)
        if isinstance(ligand_paths, str):
            self.ligand_paths = [os.path.abspath(ligand_paths)]
        elif isinstance(ligand_paths, list):
            self.ligand_paths = [os.path.abspath(lp) for lp in ligand_paths]
        else:
            raise TypeError(f"ligand_paths must be str or list, got {type(ligand_paths)}")

        self.ligand_count = len(self.ligand_paths)
        self.output_dir = os.path.abspath(output_dir)
        self.ligand_forcefield = ligand_forcefield.lower()

        # Handle forcefield_path: can be None, single string, or list
        if forcefield_path is None:
            self.forcefield_paths = None
        elif isinstance(forcefield_path, str):
            self.forcefield_paths = [forcefield_path]
        elif isinstance(forcefield_path, list):
            self.forcefield_paths = forcefield_path
        else:
            raise TypeError(f"forcefield_path must be None, str, or list, got {type(forcefield_path)}")

        # Validate ligand force field
        if self.ligand_forcefield not in ["gaff", "gaff2", "openff", "cgenff", "opls", "mmff", "match", "hybrid"]:
            raise ValueError(
                f"Unsupported ligand force field: {self.ligand_forcefield}. Use 'gaff', 'gaff2', 'openff', 'cgenff', 'opls', 'mmff', 'match', or 'hybrid'"
            )

        # Validate cgenff requires forcefield_path
        if self.ligand_forcefield == "cgenff":
            if not self.forcefield_paths:
                raise ValueError(
                    "CGenFF force field requires --forcefield-path to specify the directory containing downloaded CGenFF files"
                )
            if len(self.forcefield_paths) != self.ligand_count:
                raise ValueError(
                    f"CGenFF requires one --forcefield-path per ligand. Expected {self.ligand_count}, got {len(self.forcefield_paths)}"
                )

        # Gaussian RESP charge options
        self.gaussian_method = gaussian_method  # None, 'hf', 'dft'
        self.do_optimization = do_optimization
        self.gaussian_nproc = gaussian_nproc
        self.gaussian_mem = gaussian_mem

        # PMF mode options
        self.pmf_mode = pmf_mode
        self.pullvec = pullvec  # (protein_atom_idx, ligand_atom_idx) or None
        self.pull_mode = pull_mode  # "pocket" or "whole_protein"
        self.box_extension = box_extension if box_extension else (0.0, 0.0, 2.0)

        # REST2 mode options
        self.rest2_mode = rest2_mode
        self.t_ref = t_ref
        self.t_max = t_max
        self.n_replicas = n_replicas
        self.rest2_cutoff = rest2_cutoff

        # MMPBSA mode options
        self.mmpbsa_mode = mmpbsa_mode
        self.mmpbsa_traj_ns = mmpbsa_traj_ns
        self.gmx2amber = gmx2amber

        # FEP mode options
        self.fep_mode = fep_mode
        self.mutant_ligand = mutant_ligand
        self.fep_config = fep_config
        self.distance_cutoff = distance_cutoff
        self.charge_cutoff = charge_cutoff
        self.charge_reception = charge_reception
        self.charge_strategy = charge_strategy
        self.lambda_windows = lambda_windows
        self.lambda_strategy = lambda_strategy
        self.lambda_distribution = lambda_distribution

        # Handle resp_files - normalize to list or None
        if resp_files is None:
            self.resp_files = None
        elif isinstance(resp_files, str):
            self.resp_files = [resp_files]
        else:
            self.resp_files = list(resp_files)

        # Validate resp_files count if provided
        if self.resp_files and len(self.resp_files) != self.ligand_count:
            print_warning(
                f"Number of RESP files ({len(self.resp_files)}) doesn't match number of ligands ({self.ligand_count})"
            )
            print_warning("RESP files will be matched to ligands in order, extras ignored")

        # Initialize GROMACS environment
        self.gromacs_env = GromacsEnvironment()

        # Process force field names before initializing config
        default_ff = forcefield if forcefield else "amber99sb"
        default_water = water_model if water_model else "tip3p"

        # Initialize configuration manager
        self.config_manager = ConfigurationManager(
            config_path, self.gromacs_env, forcefield_name=default_ff, water_model_name=default_water
        )
        self.config = self.config_manager.config

        # Add ligand force field to config
        self.config["ligand_forcefield"] = {
            "type": self.ligand_forcefield,
            "charge": self.config.get("simulation", {}).get("ligand_charge", 0),
        }

        # Override config with explicit parameters if provided
        if forcefield is not None:
            self.config_manager.update_forcefield_by_name(forcefield)
        if water_model is not None:
            self.config_manager.update_water_model_by_name(water_model)
        if overwrite is not None:
            self.config["general"]["overwrite"] = overwrite

        # Extract configuration values
        self.overwrite = self.config["general"]["overwrite"]

        # Load FEP-specific configuration file if provided
        if self.fep_config and os.path.exists(self.fep_config):
            import yaml

            print(f"Loading FEP configuration from: {self.fep_config}")
            with open(self.fep_config, "r") as f:
                fep_config_data = yaml.safe_load(f)

            # Merge FEP config into main config (under 'fep' key)
            if "fep" not in self.config:
                self.config["fep"] = {}
            self.config["fep"].update(fep_config_data)

            # Also merge selected FEP config sections to top level for script/MDP generation.
            for section in ["simulation", "electrostatics", "vdw", "output", "execution"]:
                if section in fep_config_data:
                    if section not in self.config:
                        self.config[section] = {}
                    self.config[section].update(fep_config_data[section])
                    if section == "simulation":
                        prod_time = self.config["simulation"].get("production_time_ns", "NOT SET")
                        print(f"  ✓ Loaded simulation config: production_time_ns={prod_time} ns")

            mapping_cfg = fep_config_data.get("mapping", {})
            lambda_cfg = fep_config_data.get("lambda", {})
            if "dist_cutoff" in mapping_cfg:
                self.distance_cutoff = mapping_cfg["dist_cutoff"]
            if "charge_cutoff" in mapping_cfg:
                self.charge_cutoff = mapping_cfg["charge_cutoff"]
            if "charge_reception" in mapping_cfg:
                self.charge_reception = mapping_cfg["charge_reception"]
            if "charge_common" in mapping_cfg:
                self.charge_strategy = mapping_cfg["charge_common"]
            if "windows" in lambda_cfg:
                self.lambda_windows = lambda_cfg["windows"]
            if "strategy" in lambda_cfg:
                self.lambda_strategy = lambda_cfg["strategy"]
            if "distribution" in lambda_cfg:
                self.lambda_distribution = lambda_cfg["distribution"]

            # Also merge execution config at top level if present
            if "execution" in fep_config_data:
                if "execution" not in self.config:
                    self.config["execution"] = {}
                self.config["execution"].update(fep_config_data["execution"])
                print(f"  ✓ Loaded execution config: mode={self.config['execution'].get('mode', 'standard')}")
                print(f"  ✓ GPU configuration: num_gpus={self.config['execution'].get('num_gpus', 1)}")
                print(f"  ✓ Parallel windows: {self.config['execution'].get('parallel_windows', 1)}")

        # Extract FEP configuration from config file (if available)
        fep_cfg = self.config.get("fep", {})
        if fep_cfg:
            # Override lambda parameters from config file if specified
            if "strategy" in fep_cfg and lambda_strategy == "decoupled":
                self.lambda_strategy = fep_cfg["strategy"]
            if "distribution" in fep_cfg and lambda_distribution == "nonlinear":
                self.lambda_distribution = fep_cfg["distribution"]
            if "lambda_windows" in fep_cfg and lambda_windows == 11:
                self.lambda_windows = fep_cfg["lambda_windows"]
        self.forcefield_idx = self.config["forcefield"]["index"]
        self.water_model_idx = self.config["water_model"]["index"]

        # Get force field and water model info
        self.forcefield = self._get_forcefield_info()
        self.water_model = self._get_water_model_info()

        # Extract names
        self.protein_name = Path(protein_path).stem
        self.ligand_names = [Path(lp).stem for lp in self.ligand_paths]

        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

        # Subdirectories - now a list for multiple ligands
        self.lig_ff_dirs = []

        # Terminal types for pdb2gmx
        self.nter = nter
        self.cter = cter

        # Initialize sub-components
        self.mdp_generator = MDPGenerator(self.config, self.output_dir)
        self.system_builder = SystemBuilder(
            self.config, self.output_dir, self.overwrite, pmf_mode=self.pmf_mode, box_extension=self.box_extension
        )

        self._print_initialization_info()

    def _print_initialization_info(self):
        """Print initialization information"""
        print(f"\nInitialized PRISM Builder:")
        print(f"  GROMACS command: {self.gromacs_env.gmx_command}")

        # Show PMF mode if enabled
        if self.pmf_mode:
            print(f"  Mode: {success('PMF (Steered MD)')}")
            print(
                f"  Box extension (X,Y,Z): ({self.box_extension[0]:.1f}, {self.box_extension[1]:.1f}, {self.box_extension[2]:.1f}) nm"
            )
            if self.pullvec:
                print(f"  Pull vector: Protein atom {self.pullvec[0]} -> Ligand atom {self.pullvec[1]}")
            else:
                mode_desc = (
                    "collision-based, whole protein" if self.pull_mode == "whole_protein" else "pocket clearance"
                )
                print(f"  Pull vector: Auto ({mode_desc})")

        # Show REST2 mode if enabled
        if self.rest2_mode:
            print(f"  Mode: {success('REST2 (Replica Exchange Solute Tempering)')}")
            print(f"  T_ref: {self.t_ref} K")
            print(f"  T_max: {self.t_max} K")
            print(f"  Replicas: {self.n_replicas}")
            print(f"  Cutoff: {self.rest2_cutoff} nm")

        # Show MMPBSA mode if enabled
        if self.mmpbsa_mode:
            if self.mmpbsa_traj_ns:
                print(f"  Mode: {success('MM/PBSA (Trajectory)')}")
                print(f"  Production MD: {self.mmpbsa_traj_ns} ns")
            else:
                print(f"  Mode: {success('MM/PBSA (Single-frame)')}")
            if self.gmx2amber:
                print(f"  Backend: AMBER MMPBSA.py (parmed conversion)")
            else:
                print(f"  Backend: gmx_MMPBSA")

        print(f"  Protein: {self.protein_path}")
        if self.ligand_count == 1:
            print(f"  Ligand: {self.ligand_paths[0]}")
        else:
            print(f"  Ligands ({number(self.ligand_count)}):")
            for i, lp in enumerate(self.ligand_paths, 1):
                print(f"    {number(i)}. {path(lp)}")
        print(f"  Ligand force field: {self.ligand_forcefield.upper()}")

        # Show Gaussian RESP options
        if self.gaussian_method or self.resp_files:
            if self.resp_files:
                print(f"  Charge method: RESP (from provided files)")
                for i, rf in enumerate(self.resp_files, 1):
                    print(f"    {number(i)}. {path(rf)}")
            else:
                print(f"  Charge method: Gaussian RESP ({self.gaussian_method.upper()}/6-31G*)")
                print(f"  Geometry optimization: {'Yes' if self.do_optimization else 'No'}")
                print(f"  Gaussian resources: {self.gaussian_nproc} CPUs, {self.gaussian_mem} memory")
        else:
            print(f"  Charge method: AM1-BCC (default)")

        print(f"  Output directory: {self.output_dir}")
        print(f"  Protein force field: {self.forcefield['name']}")
        print(f"  Water model: {self.water_model['name']}")
        print(f"  Box distance: {self.config['box']['distance']} nm")
        print(f"  Temperature: {self.config['simulation']['temperature']} K")
        print(f"  pH: {self.config['simulation']['pH']}")
        print(f"  Production time: {self.config['simulation']['production_time_ns']} ns")

    def _get_forcefield_info(self):
        """Get force field information from config"""
        ff_idx = self.config["forcefield"]["index"]
        custom_ff = self.config["forcefield"]["custom_forcefields"]

        if ff_idx in custom_ff:
            return custom_ff[ff_idx]
        else:
            raise ValueError(f"Force field index {ff_idx} not found in configuration")

    def _get_water_model_info(self):
        """Get water model information from config"""
        wm_idx = self.config["water_model"]["index"]
        custom_wm = self.config["water_model"]["custom_water_models"]

        if wm_idx in custom_wm:
            return custom_wm[wm_idx]
        else:
            raise ValueError(f"Water model index {wm_idx} not found in configuration")

    # All workflow orchestration methods are provided by mixins:
    # - WorkflowMixin: run(), run_normal(), run_rest2(), run_mmpbsa()
    # - FEPWorkflowMixin: run_fep() and all FEP helper methods
    # - LigandForceFieldMixin: generate_ligand_forcefield() and RESP workflow
    # - SystemPreparationMixin: clean_protein(), build_model(), generate_mdp_files(), etc.
    # - PMFBuilderMixin: run_pmf() and PMF-specific methods
