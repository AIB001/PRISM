#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder Core - Main PRISMBuilder class

This module contains the PRISMBuilder class with all core workflow methods.
PMF-related methods are provided by PMFBuilderMixin (builder/pmf.py).
"""

import os
import shutil
from pathlib import Path

from ..utils.environment import GromacsEnvironment
from ..utils.config import ConfigurationManager
from ..utils.mdp import MDPGenerator
from ..utils.system import SystemBuilder
from ..utils.colors import (
    print_header,
    print_subheader,
    print_step,
    print_success,
    print_error,
    print_warning,
    print_info,
    success,
    path,
    number,
)
from ..forcefield.gaff import GAFFForceFieldGenerator
from ..forcefield.gaff2 import GAFF2ForceFieldGenerator
from ..forcefield.openff import OpenFFForceFieldGenerator
from ..forcefield.cgenff import CGenFFForceFieldGenerator
from ..forcefield.opls_aa import OPLSAAForceFieldGenerator
from ..forcefield.swissparam import (
    MMFFForceFieldGenerator,
    MATCHForceFieldGenerator,
    HybridMMFFMATCHForceFieldGenerator,
)

from .pmf import PMFBuilderMixin


class PRISMBuilder(PMFBuilderMixin):
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

            # Also merge config sections to top level for write_fep_mdps
            # This allows write_fep_mdps to find production_time_ns, rcoulomb, etc.
            for section in ["simulation", "electrostatics", "vdw", "output"]:
                if section in fep_config_data:
                    if section not in self.config:
                        self.config[section] = {}
                    self.config[section].update(fep_config_data[section])
                    if section == "simulation":
                        prod_time = self.config["simulation"].get("production_time_ns", "NOT SET")
                        print(f"  ✓ Loaded simulation config: production_time_ns={prod_time} ns")

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

    def generate_ligand_forcefield(self):
        """Generate ligand force fields for all ligands using selected force field generator"""
        import concurrent.futures

        print_subheader(f"Generating Ligand Force Fields ({self.ligand_forcefield.upper()})")
        print(f"Processing {number(self.ligand_count)} ligand(s)...")

        # Determine force field suffix based on type
        if self.ligand_forcefield in ["gaff", "gaff2"]:
            ff_suffix = "amb2gmx"
        elif self.ligand_forcefield == "openff":
            ff_suffix = "openff2gmx"
        elif self.ligand_forcefield == "cgenff":
            ff_suffix = "cgenff2gmx"
        elif self.ligand_forcefield == "opls":
            ff_suffix = "opls2gmx"
        elif self.ligand_forcefield in ["mmff", "match", "hybrid"]:
            ff_suffix = "sp2gmx"  # SwissParam
        else:
            ff_suffix = "lig2gmx"

        # Create base directory for ligand force fields
        # Single ligand: place directly in output_dir (backward compatibility)
        # Multi-ligand: place in Ligand_Forcefield/ subdirectory (new feature)
        if self.ligand_count == 1:
            ligand_ff_base_dir = self.output_dir
        else:
            ligand_ff_base_dir = os.path.join(self.output_dir, "Ligand_Forcefield")
            os.makedirs(ligand_ff_base_dir, exist_ok=True)

        def generate_single_ligand_ff(ligand_idx, ligand_path):
            """Generate force field for a single ligand"""
            ligand_num = ligand_idx + 1  # 1-based numbering
            print(
                f"\n{success('→')} Processing ligand {number(ligand_num)}/{number(self.ligand_count)}: {path(os.path.basename(ligand_path))}"
            )

            # Create temporary output directory for this ligand
            temp_output_dir = os.path.join(self.output_dir, f"_temp_lig_{ligand_num}")
            os.makedirs(temp_output_dir, exist_ok=True)

            # Check for existing or provided RESP file
            resp_mol2 = self._find_resp_file(ligand_path, ligand_idx)

            # Determine charge mode for GAFF/GAFF2
            # Use 'gas' (fast) if RESP charges will be applied, 'bcc' otherwise
            if resp_mol2 or self.gaussian_method:
                charge_mode = "gas"
                print(f"  Using 'gas' charge mode (RESP charges will be applied)")
            else:
                charge_mode = "bcc"

            # Generate force field based on type
            if self.ligand_forcefield == "gaff":
                generator = GAFFForceFieldGenerator(
                    ligand_path, temp_output_dir, overwrite=self.overwrite, charge_mode=charge_mode
                )
            elif self.ligand_forcefield == "gaff2":
                generator = GAFF2ForceFieldGenerator(
                    ligand_path, temp_output_dir, overwrite=self.overwrite, charge_mode=charge_mode
                )
            elif self.ligand_forcefield == "openff":
                generator = OpenFFForceFieldGenerator(
                    ligand_path,
                    temp_output_dir,
                    charge=self.config["ligand_forcefield"]["charge"],
                    overwrite=self.overwrite,
                )
            elif self.ligand_forcefield == "cgenff":
                # For CGenFF, use the corresponding forcefield_path for this ligand
                cgenff_dir = self.forcefield_paths[ligand_idx] if self.forcefield_paths else None
                generator = CGenFFForceFieldGenerator(
                    ligand_path, temp_output_dir, cgenff_dir=cgenff_dir, overwrite=self.overwrite
                )
            elif self.ligand_forcefield == "opls":
                generator = OPLSAAForceFieldGenerator(
                    ligand_path,
                    temp_output_dir,
                    charge=self.config["ligand_forcefield"]["charge"],
                    charge_model="cm1a",
                    align_to_input=True,
                    overwrite=self.overwrite,
                )
            elif self.ligand_forcefield == "mmff":
                generator = MMFFForceFieldGenerator(ligand_path, temp_output_dir, overwrite=self.overwrite)
            elif self.ligand_forcefield == "match":
                generator = MATCHForceFieldGenerator(ligand_path, temp_output_dir, overwrite=self.overwrite)
            elif self.ligand_forcefield == "hybrid":
                generator = HybridMMFFMATCHForceFieldGenerator(ligand_path, temp_output_dir, overwrite=self.overwrite)

            temp_ff_dir = generator.run()

            # Create final directory with proper naming
            # Single ligand: LIG.xxx2gmx (no suffix, backward compatibility)
            # Multi-ligand: LIG.xxx2gmx_1, LIG.xxx2gmx_2, etc.
            if self.ligand_count == 1:
                final_ff_dir = os.path.join(ligand_ff_base_dir, f"LIG.{ff_suffix}")
            else:
                final_ff_dir = os.path.join(ligand_ff_base_dir, f"LIG.{ff_suffix}_{ligand_num}")

            # Move files from temp directory to final location
            if os.path.exists(final_ff_dir) and self.overwrite:
                shutil.rmtree(final_ff_dir)
            shutil.move(temp_ff_dir, final_ff_dir)

            # Clean up temp directory
            if os.path.exists(temp_output_dir):
                shutil.rmtree(temp_output_dir)

            # Only rename if we have multiple ligands (ligand_count > 1)
            # Single ligand keeps the name "LIG" for backward compatibility
            if self.ligand_count > 1:
                itp_file = os.path.join(final_ff_dir, "LIG.itp")
                gro_file = os.path.join(final_ff_dir, "LIG.gro")

                # Rename moleculetype in ITP file
                with open(itp_file, "r") as f:
                    itp_lines = f.readlines()

                with open(itp_file, "w") as f:
                    in_moleculetype = False
                    moleculetype_replaced = False
                    for line in itp_lines:
                        # Detect [ moleculetype ] section
                        if "[ moleculetype ]" in line:
                            in_moleculetype = True
                            f.write(line)
                            continue

                        # Replace the moleculetype name in the next non-comment, non-empty line
                        if in_moleculetype and not moleculetype_replaced:
                            if line.strip() and not line.strip().startswith(";"):
                                # Replace LIG with LIG_N (handle various spacing formats)
                                import re

                                line = re.sub(r"\bLIG\b", f"LIG_{ligand_num}", line, count=1)
                                moleculetype_replaced = True
                                in_moleculetype = False

                        # Replace residue names in [ atoms ] section
                        # Typical format: "     1 atomtype    1 LIG     O     ..."
                        if not line.strip().startswith(";") and " LIG " in line:
                            line = re.sub(r"(\s+)LIG(\s+)", rf"\1LIG_{ligand_num}\2", line)

                        f.write(line)

                # Rename residue in GRO file - more robust replacement
                with open(gro_file, "r") as f:
                    gro_lines = f.readlines()
                with open(gro_file, "w") as f:
                    for i, line in enumerate(gro_lines):
                        # Skip header (line 0) and atom count (line 1) and box vectors (last line)
                        if i > 1 and i < len(gro_lines) - 1:
                            # GRO format atom line: columns 0-4 (residue number), 5-9 (residue name), 10-14 (atom name)
                            if len(line) >= 20:  # Valid atom line
                                res_name = line[5:10].strip()
                                if res_name.startswith("LIG"):
                                    # Replace LIG with LIG_N, pad to 5 chars, right-aligned
                                    new_res_name = f"LIG_{ligand_num}"
                                    # GRO residue name is left-aligned in a 5-char field
                                    line = line[:5] + new_res_name.ljust(5)[:5] + line[10:]
                        f.write(line)

            # Verify required files exist
            required_files = ["LIG.itp", "LIG.gro"]
            for filename in required_files:
                filepath = os.path.join(final_ff_dir, filename)
                if not os.path.exists(filepath):
                    raise FileNotFoundError(f"Required file not found: {filepath}")

            # Handle RESP charges
            if resp_mol2:
                # Apply existing RESP charges
                self._apply_resp_charges(final_ff_dir, resp_mol2, ligand_num)
            elif self.gaussian_method and self.ligand_forcefield in ["gaff", "gaff2"]:
                # Generate Gaussian RESP workflow files
                resp_mol2 = self._handle_gaussian_workflow(ligand_path, final_ff_dir, ligand_idx)
                if resp_mol2:
                    self._apply_resp_charges(final_ff_dir, resp_mol2, ligand_num)

            print_success(f"  Ligand {number(ligand_num)} force field generated: {path(final_ff_dir)}")
            return final_ff_dir

        # Generate force fields in parallel using ThreadPoolExecutor
        print_info("Generating force fields in parallel...")
        self.lig_ff_dirs = []

        with concurrent.futures.ThreadPoolExecutor(max_workers=min(self.ligand_count, 4)) as executor:
            futures = {executor.submit(generate_single_ligand_ff, i, lp): i for i, lp in enumerate(self.ligand_paths)}

            for future in concurrent.futures.as_completed(futures):
                ligand_idx = futures[future]
                try:
                    ff_dir = future.result()
                    self.lig_ff_dirs.append((ligand_idx, ff_dir))
                except Exception as e:
                    print_error(f"Failed to generate force field for ligand {ligand_idx + 1}: {e}")
                    raise

        # Sort by ligand index to maintain order
        self.lig_ff_dirs.sort(key=lambda x: x[0])
        self.lig_ff_dirs = [ff_dir for _, ff_dir in self.lig_ff_dirs]

        print_success(f"All {number(self.ligand_count)} ligand force fields generated successfully")
        print(f"Force field directory: {path(ligand_ff_base_dir)}")

        return ligand_ff_base_dir

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
            If not specified, reads from config or defaults to 5.0 A
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
        from ..utils.cleaner import ProteinCleaner

        print_subheader("Cleaning Protein")

        cleaned_pdb = os.path.join(self.output_dir, f"{self.protein_name}_clean.pdb")

        # Check for existing cleaned file
        if os.path.exists(cleaned_pdb) and not self.overwrite:
            print(f"Using existing cleaned protein: {cleaned_pdb}")
            return cleaned_pdb

        # Get parameters from config if not explicitly provided
        if ion_mode is None:
            ion_mode = self.config.get("protein_preparation", {}).get("ion_handling_mode", "smart")
        if distance_cutoff is None:
            distance_cutoff = self.config.get("protein_preparation", {}).get("metal_distance_cutoff", 5.0)
        if keep_crystal_water is None:
            keep_crystal_water = self.config.get("protein_preparation", {}).get("keep_crystal_water", False)
        if remove_artifacts is None:
            remove_artifacts = self.config.get("protein_preparation", {}).get("remove_crystallization_artifacts", True)
        nterm_met = self.config.get("protein_preparation", {}).get("nterm_met", "keep")

        # Get custom metals from config
        keep_custom = self.config.get("protein_preparation", {}).get("keep_custom_metals", [])
        remove_custom = self.config.get("protein_preparation", {}).get("remove_custom_metals", [])

        # Initialize cleaner
        cleaner = ProteinCleaner(
            ion_mode=ion_mode,
            distance_cutoff=distance_cutoff,
            keep_crystal_water=keep_crystal_water,
            remove_artifacts=remove_artifacts,
            keep_custom_metals=keep_custom if keep_custom else None,
            remove_custom_metals=remove_custom if remove_custom else None,
            verbose=True,
            drop_nterm_met=nterm_met,
            forcefield_name=self.forcefield["name"] if self.forcefield else None,
        )

        # Clean the protein
        cleaner.clean_pdb(self.protein_path, cleaned_pdb)

        # Some prepared receptors keep non-canonical terminal cap atoms
        # (e.g. CAY/CY/OY or CAT/NT) on standard amino-acid residue names.
        # Strip them immediately so downstream pdbfixer/pdb2gmx sees a
        # canonical protein backbone.
        from prism.utils.cleaner import fix_terminal_atoms

        fix_terminal_atoms(
            cleaned_pdb,
            cleaned_pdb,
            force_field=self.forcefield["name"] if self.forcefield else None,
            verbose=True,
        )

        # Post-process: Fix terminal hydrogen names for AMBER compatibility
        self._fix_hydrogen_names(cleaned_pdb)

        print(f"Protein cleaned and saved to: {cleaned_pdb}")

        # Note: PROPKA HIS renaming (if --protonation propka) is applied later
        # in system.py's build() AFTER pdbfixer, to avoid pdbfixer reverting names.
        return cleaned_pdb

    def _fix_hydrogen_names(self, pdb_file):
        """
        Fix terminal hydrogen names for AMBER compatibility.

        AMBER force fields expect specific naming for N-terminal hydrogens:
        H1, H2, H3 instead of HN1, HN2, HN3
        """
        with open(pdb_file, "r") as f:
            lines = f.readlines()

        fixed_lines = []
        for line in lines:
            if line.startswith("ATOM"):
                # Fix terminal hydrogen names
                if "HN1" in line:
                    line = line.replace("HN1", "H1 ")
                elif "HN2" in line:
                    line = line.replace("HN2", "H2 ")
                elif "HN3" in line:
                    line = line.replace("HN3", "H3 ")
            fixed_lines.append(line)

        with open(pdb_file, "w") as f:
            f.writelines(fixed_lines)

    def _find_resp_file(self, ligand_path, ligand_idx):
        """
        Find RESP mol2 file for a ligand.

        Checks in order:
        1. Explicitly provided via --respfile
        2. Auto-detect ligand_resp.mol2 in same directory as ligand

        Parameters
        ----------
        ligand_path : str
            Path to ligand MOL2/SDF file
        ligand_idx : int
            Index of ligand (0-based)

        Returns
        -------
        str or None
            Path to RESP mol2 file, or None if not found
        """
        # Check --respfile parameter
        if self.resp_files and ligand_idx < len(self.resp_files):
            resp_file = self.resp_files[ligand_idx]
            if os.path.exists(resp_file):
                print_info(f"  Using provided RESP file: {path(resp_file)}")
                return resp_file
            else:
                print_warning(f"  RESP file not found: {resp_file}")

        # Auto-detect ligand_resp.mol2 in same directory
        ligand_dir = os.path.dirname(ligand_path)
        auto_resp = os.path.join(ligand_dir, "ligand_resp.mol2")
        if os.path.exists(auto_resp):
            print_info(f"  Auto-detected RESP file: {path(auto_resp)}")
            return auto_resp

        return None

    def _handle_gaussian_workflow(self, ligand_path, ff_dir, ligand_idx):
        """
        Handle Gaussian RESP workflow for a ligand.

        If Gaussian (g16) is available, runs the calculation automatically.
        Otherwise, generates input files and a run script for manual execution.

        Parameters
        ----------
        ligand_path : str
            Path to ligand MOL2/SDF file
        ff_dir : str
            Force field output directory
        ligand_idx : int
            Index of ligand (0-based)

        Returns
        -------
        str or None
            Path to ligand_resp.mol2 if generated, None otherwise
        """
        try:
            from ..gaussian import GaussianRESPWorkflow, check_gaussian_available
        except ImportError:
            print_error("Gaussian module not available")
            return None

        print_info(f"  Setting up Gaussian RESP workflow...")
        if ligand_idx is not None:
            print_info(f"  Ligand index: {ligand_idx}")

        # Create Gaussian working directory
        gaussian_dir = os.path.join(ff_dir, "gaussian_resp")
        os.makedirs(gaussian_dir, exist_ok=True)

        # Initialize workflow
        workflow = GaussianRESPWorkflow(
            ligand_path=ligand_path,
            output_dir=gaussian_dir,
            method=self.gaussian_method,
            do_optimization=self.do_optimization,
            nproc=self.gaussian_nproc,
            mem=self.gaussian_mem,
            verbose=True,
        )

        # Check if Gaussian is available
        if check_gaussian_available():
            print_info(f"  Gaussian (g16) detected, running RESP calculation...")
            resp_mol2 = workflow.run()
            if resp_mol2:
                print_success(f"  RESP charges calculated: {path(resp_mol2)}")
                return resp_mol2
            else:
                print_warning(f"  Gaussian calculation failed")
                return None
        else:
            # Generate files for manual execution
            files = workflow.prepare_files()
            print_warning(f"  Gaussian (g16) not found in PATH")
            print_warning(f"  Generated Gaussian input files in: {path(gaussian_dir)}")
            print_warning(f"  To calculate RESP charges:")
            print_warning(f"    1. Run: bash {files['script']}")
            print_warning(f"    2. Re-run PRISM with: --respfile {os.path.join(gaussian_dir, 'ligand_resp.mol2')}")
            return None

    def _apply_resp_charges(self, ff_dir, resp_mol2, ligand_num):
        """
        Apply RESP charges to ITP file.

        Parameters
        ----------
        ff_dir : str
            Force field directory containing LIG.itp
        resp_mol2 : str
            Path to RESP mol2 file with charges
        ligand_num : int
            Ligand number (1-based, for logging)
        """
        try:
            from ..gaussian import RESPChargeReplacer
        except ImportError:
            print_error("Gaussian module not available for charge replacement")
            return

        itp_path = os.path.join(ff_dir, "LIG.itp")
        if not os.path.exists(itp_path):
            print_error(f"  ITP file not found: {itp_path}")
            return

        print_info(f"  Applying RESP charges to ligand {number(ligand_num)}...")

        try:
            replacer = RESPChargeReplacer(resp_mol2, verbose=False)
            replacer.replace_itp_charges(itp_path, backup=True)
            print_success(f"  RESP charges applied from: {path(resp_mol2)}")
        except Exception as e:
            print_error(f"  Failed to apply RESP charges: {e}")

    def build_model(self, cleaned_protein):
        """Build the GROMACS model with support for multiple ligands"""
        return self.system_builder.build(
            cleaned_protein,
            self.lig_ff_dirs,  # Now a list of force field directories
            self.forcefield_idx,
            self.water_model_idx,
            self.forcefield,  # Pass full force field info (name, dir, path)
            self.water_model,  # Pass full water model info
            nter=self.nter,  # N-terminal type
            cter=self.cter,  # C-terminal type
        )

    def generate_mdp_files(self):
        """Generate MDP files for MD simulations"""
        self.mdp_generator.generate_all()

    def cleanup(self):
        """Clean up temporary files"""
        print_subheader("Cleaning up temporary files")

        # Cleanup directories for GAFF, GAFF2, and OpenFF
        cleanup_dirs = ["forcefield", "temp_openff"]

        for dir_name in cleanup_dirs:
            cleanup_dir = os.path.join(self.output_dir, dir_name)
            if os.path.exists(cleanup_dir):
                if self.ligand_forcefield in ["gaff", "gaff2"]:
                    temp_patterns = [
                        "*.frcmod",
                        "*.prep",
                        "*.prmtop",
                        "*.rst7",
                        "*.log",
                        "*.in",
                        "ANTECHAMBER*",
                        "ATOMTYPE*",
                        "PREP*",
                        "NEWPDB*",
                        "sqm*",
                        "leap*",
                    ]

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
        gmx_md_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
        if not os.path.exists(gmx_md_dir):
            print("Warning: GMX_PROLIG_MD directory not found, skipping localrun.sh generation")
            return None

        # Clean up Emacs backup files (#topol.top.1#, #topol.top.2#, etc.)
        import glob

        backup_pattern = os.path.join(gmx_md_dir, "#*#")
        backup_files = glob.glob(backup_pattern)
        if backup_files:
            print(f"Cleaning up {len(backup_files)} Emacs backup file(s)...")
            for backup_file in backup_files:
                try:
                    os.remove(backup_file)
                    print(f"  Removed: {os.path.basename(backup_file)}")
                except Exception as e:
                    print(f"  Warning: Could not remove {backup_file}: {e}")

        script_path = os.path.join(gmx_md_dir, "localrun.sh")

        script_content = """#!/bin/bash

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
"""

        with open(script_path, "w") as f:
            f.write(script_content)

        # Make script executable
        os.chmod(script_path, 0o755)

        print(f"Local run script generated: {script_path}")
        return script_path

    def run(self):
        """Run the complete workflow"""
        if self.pmf_mode:
            return self.run_pmf()
        elif self.rest2_mode:
            return self.run_rest2()
        elif self.mmpbsa_mode:
            return self.run_mmpbsa()
        elif self.fep_mode:
            return self.run_fep()
        else:
            return self.run_normal()

    def run_fep(self):
        """Run the FEP workflow: build standard MD systems, generate hybrid topology, create FEP scaffold"""
        from ..fep.modeling import FEPScaffoldBuilder
        from ..fep.gromacs.itp_builder import ITPBuilder

        print_header("PRISM FEP Builder Workflow")

        if not self.mutant_ligand:
            raise ValueError("FEP mode requires --mutant ligand file")

        try:
            # Create FEP output directory first
            # Avoid double-nesting: if output_dir already ends with GMX_PROLIG_FEP, use it directly
            if os.path.basename(os.path.abspath(self.output_dir)) == "GMX_PROLIG_FEP":
                fep_output = os.path.abspath(self.output_dir)
            else:
                fep_output = os.path.join(self.output_dir, "GMX_PROLIG_FEP")
            os.makedirs(fep_output, exist_ok=True)

            # Phase 1: Build bound system with reference ligand
            print_step(1, 5, "Building bound system with reference ligand")
            bound_output = os.path.join(fep_output, "_build/bound_md")

            self._build_standard_system(bound_output, use_protein=True)
            bound_system_dir = os.path.join(bound_output, "GMX_PROLIG_MD")

            # Save reference ligand FF directory from the freshly built bound system
            ref_ff_dir = self._resolve_generated_ligand_ff_dir(bound_output)
            print_success(f"Bound system built: {bound_system_dir}")
            print(f"  Reference ligand FF: {ref_ff_dir}")

            # Phase 2: Build unbound system (ligand only) with same box size as bound
            print_step(2, 5, "Building unbound system (ligand in water)")
            unbound_output = os.path.join(fep_output, "_build/unbound_md")

            # Read bound system box size
            bound_conf = os.path.join(bound_system_dir, "solv_ions.gro")
            box_size = self._read_box_size(bound_conf)
            print(f"  Using bound system box size: {box_size[0]:.3f} x {box_size[1]:.3f} x {box_size[2]:.3f} nm")

            self._build_standard_system(unbound_output, use_protein=False, box_size=box_size)
            unbound_system_dir = os.path.join(unbound_output, "GMX_PROLIG_MD")
            print_success(f"Unbound system built: {unbound_system_dir}")

            # Phase 3: Generate hybrid topology
            print_step(3, 5, "Generating hybrid topology via atom mapping")

            # ref_ff_dir was saved in Phase 1
            # Generate mutant ligand FF in _build directory
            mut_ff_output = os.path.join(fep_output, "_build/mutant_ligand_ff")
            self._generate_mutant_ligand_ff(self.mutant_ligand, mut_ff_output)
            mut_ff_dir = self._resolve_generated_ligand_ff_dir(mut_ff_output)

            # Build hybrid topology using ITPBuilder
            hybrid_output = os.path.join(fep_output, "common/hybrid")
            os.makedirs(hybrid_output, exist_ok=True)

            # Read ligand atoms from ITP and GRO files
            from ..fep.io import read_ligand_from_prism
            from ..fep.core.hybrid_topology import HybridTopologyBuilder, LigandTopologyInput

            # Debug: Check if files exist
            ref_itp = os.path.join(ref_ff_dir, "LIG.itp")
            ref_gro = os.path.join(ref_ff_dir, "LIG.gro")
            mut_itp = os.path.join(mut_ff_dir, "LIG.itp")
            mut_gro = os.path.join(mut_ff_dir, "LIG.gro")

            print(f"\n  Checking force field files:")
            print(f"    ref_itp exists: {os.path.exists(ref_itp)} - {ref_itp}")
            print(f"    ref_gro exists: {os.path.exists(ref_gro)} - {ref_gro}")
            print(f"    mut_itp exists: {os.path.exists(mut_itp)} - {mut_itp}")
            print(f"    mut_gro exists: {os.path.exists(mut_gro)} - {mut_gro}")

            ref_coord_source = self.ligand_paths[0] if self.ligand_paths else ref_gro
            if not os.path.exists(ref_coord_source):
                ref_coord_source = ref_gro

            mut_coord_source = self.mutant_ligand or mut_gro
            if not os.path.exists(mut_coord_source):
                mut_coord_source = mut_gro

            print("  Mapping coordinate sources:")
            print(f"    Reference coords: {ref_coord_source}")
            print(f"    Mutant coords:    {mut_coord_source}")

            ref_atoms = read_ligand_from_prism(itp_file=ref_itp, gro_file=ref_coord_source)

            mut_atoms = read_ligand_from_prism(itp_file=mut_itp, gro_file=mut_coord_source)

            # Perform atom mapping
            from ..fep.core.mapping import DistanceAtomMapper

            mapper = DistanceAtomMapper(dist_cutoff=self.distance_cutoff)
            mapping = mapper.map(ref_atoms, mut_atoms)

            # Build hybrid topology
            hybrid_builder = HybridTopologyBuilder(charge_strategy=self.charge_strategy)

            # Read ITP parameters for hybrid topology building
            from pathlib import Path

            ref_itp_data = ITPBuilder._parse_source_itp(Path(os.path.join(ref_ff_dir, "LIG.itp")).read_text())
            mut_itp_data = ITPBuilder._parse_source_itp(Path(os.path.join(mut_ff_dir, "LIG.itp")).read_text())

            params_a = LigandTopologyInput(
                masses={},
                bonds=ref_itp_data["sections"].get("bonds", []),
                pairs=ref_itp_data["sections"].get("pairs", []),
                angles=ref_itp_data["sections"].get("angles", []),
                dihedrals=ref_itp_data["sections"].get("dihedrals", []),
                impropers=ref_itp_data["sections"].get("impropers", []),
            )
            params_b = LigandTopologyInput(
                masses={},
                bonds=mut_itp_data["sections"].get("bonds", []),
                pairs=mut_itp_data["sections"].get("pairs", []),
                angles=mut_itp_data["sections"].get("angles", []),
                dihedrals=mut_itp_data["sections"].get("dihedrals", []),
                impropers=mut_itp_data["sections"].get("impropers", []),
            )

            hybrid_atoms = hybrid_builder.build(mapping, params_a, params_b, ref_atoms, mut_atoms)

            # Build hybrid bonded parameters from source ITPs
            # First write a temporary atoms-only ITP
            temp_itp = os.path.join(hybrid_output, "hybrid_atoms_temp.itp")
            ITPBuilder(hybrid_atoms, {}).write_itp(temp_itp, molecule_name="HYB")
            print(f"  Debug: temp_itp written: {os.path.exists(temp_itp)}")

            # Then build complete hybrid ITP with bonded terms
            hybrid_itp = os.path.join(hybrid_output, "hybrid.itp")
            hybrid_params = ITPBuilder.write_complete_hybrid_itp(
                output_path=hybrid_itp,
                hybrid_itp=temp_itp,
                ligand_a_itp=os.path.join(ref_ff_dir, "LIG.itp"),
                ligand_b_itp=os.path.join(mut_ff_dir, "LIG.itp"),
                molecule_name="HYB",
            )
            print(f"  Debug: hybrid.itp written: {os.path.exists(hybrid_itp)}")
            if os.path.exists(hybrid_itp):
                print(f"  Debug: hybrid.itp size: {os.path.getsize(hybrid_itp)} bytes")

            print_success(f"Hybrid topology: {hybrid_itp}")

            # Generate hybrid ligand structure file
            hybrid_gro = os.path.join(hybrid_output, "hybrid.gro")
            self._generate_hybrid_gro(hybrid_gro, hybrid_atoms, mapping, ref_ff_dir, mut_ff_dir)
            print_success(f"Hybrid structure: {hybrid_gro}")

            # Phase 4: Create FEP scaffold with complete systems
            print_step(4, 5, "Creating FEP scaffold with complete systems")

            fep_builder = FEPScaffoldBuilder(
                output_dir=fep_output,
                lambda_windows=self.lambda_windows,
                lambda_strategy=self.lambda_strategy,
                lambda_distribution=self.lambda_distribution,
                config=self.config,
                overwrite=False,
            )

            layout = fep_builder.build_from_components(
                receptor_pdb=self.protein_path,
                hybrid_itp=hybrid_itp,
                hybrid_gro=hybrid_gro,
                reference_ligand_dir=ref_ff_dir,
                mutant_ligand_dir=mut_ff_dir,
                bound_system_dir=bound_system_dir,
                unbound_system_dir=unbound_system_dir,
            )

            print_success(f"FEP scaffold created: {fep_output}")

            # Phase 5: Clean up intermediate build files
            print_step(5, 5, "Cleaning up intermediate build files")

            import shutil

            build_dir = os.path.join(fep_output, "_build")
            if os.path.exists(build_dir):
                shutil.rmtree(build_dir)
                print(f"  ✓ Removed build directory: {build_dir}")

            # Also clean up any files in output_dir root (outside GMX_PROLIG_FEP)
            cleanup_items = [
                os.path.join(self.output_dir, "mdps"),
            ]
            cleanup_items.extend(str(path) for path in Path(self.output_dir).glob("LIG.*") if path.is_dir())

            for item in cleanup_items:
                if os.path.exists(item):
                    if os.path.isdir(item):
                        shutil.rmtree(item)
                    else:
                        os.remove(item)
                    print(f"  ✓ Removed {os.path.basename(item)}")

            print_success("Cleanup complete")

            print_header("FEP Workflow Complete!")
            print(f"\n  FEP scaffold:       {path(fep_output)}")
            print(f"\n  To run FEP simulations:")
            print(f"  1. cd {fep_output}/bound && ./localrun.sh")
            print(f"  2. cd {fep_output}/unbound && ./localrun.sh")
            print(f"  3. Analyze with gmx bar or alchemical-analysis")

            return self.output_dir

        except Exception as e:
            print_error(f"FEP workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise

    def _generate_hybrid_gro(
        self, output_gro: str, hybrid_atoms: list, _mapping: "AtomMapping", ref_ff_dir: str, mut_ff_dir: str
    ) -> None:
        """
        Generate hybrid ligand GRO file with dummy atoms

        Parameters
        ----------
        output_gro : str
            Output path for hybrid GRO file
        hybrid_atoms : list
            List of HybridAtom objects from hybrid topology
        mapping : AtomMapping
            Atom mapping between reference and mutant ligands
        ref_ff_dir : str
            Reference ligand force field directory (contains LIG.gro)
        mut_ff_dir : str
            Mutant ligand force field directory (contains LIG.gro)
        """

        # Read reference and mutant ligand GRO files
        ref_gro_path = os.path.join(ref_ff_dir, "LIG.gro")
        mut_gro_path = os.path.join(mut_ff_dir, "LIG.gro")

        ref_coords = self._parse_gro(ref_gro_path)
        mut_coords = self._parse_gro(mut_gro_path)

        # Build coordinate lookup by atom name
        ref_coord_map = {atom["name"]: atom["coord"] for atom in ref_coords["atoms"]}
        mut_coord_map = {atom["name"]: atom["coord"] for atom in mut_coords["atoms"]}

        # Generate hybrid coordinates
        hybrid_atoms_data = []
        for hatom in hybrid_atoms:
            # Remove A/B suffixes from hybrid atom name to find original name
            base_name = hatom.name.rstrip("AB")

            if hatom.name.endswith("A"):
                # Transformed A atom (disappearing in state B)
                # Use coordinates from reference ligand
                coord = ref_coord_map.get(base_name, [0.0, 0.0, 0.0])
            elif hatom.name.endswith("B"):
                # Transformed B atom (appearing in state B)
                # Use coordinates from mutant ligand
                coord = mut_coord_map.get(base_name, [0.0, 0.0, 0.0])
            else:
                # Common atom
                # Use coordinates from reference ligand
                coord = ref_coord_map.get(hatom.name, [0.0, 0.0, 0.0])

            hybrid_atoms_data.append({"index": hatom.index, "name": hatom.name, "coord": coord})

        # Write hybrid GRO file
        with open(output_gro, "w") as f:
            f.write("Hybrid ligand structure generated by PRISM-FEP\n")
            f.write(f"{len(hybrid_atoms_data)}\n")

            # Get residue number from reference ligand
            if ref_coords["atoms"]:
                # Parse the first atom line to get residue number
                with open(ref_gro_path, "r") as ref_file:
                    ref_lines = ref_file.readlines()
                    first_atom_line = ref_lines[2]  # Skip title and atom count
                    res_num_str = first_atom_line[0:5].strip()
                    residue_number = int(res_num_str) if res_num_str else 1
            else:
                residue_number = 1

            for atom_data in hybrid_atoms_data:
                idx = atom_data["index"]
                name = atom_data["name"]
                x, y, z = atom_data["coord"]
                # GRO format: %5d%-5s%5s%5d%8.3f%8.3f%8.3f
                # residuenum (5 chars) + residuename (5 chars, left-justified)
                # + atomname (5 chars, right-justified) + atomnumber (5 chars) + x + y + z
                f.write(f"{residue_number:5d}LIG  {name:>5s}{idx:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")

            # Write box vectors (use default 1.0 nm, will be replaced by system box)
            f.write("   1.00000   1.00000   1.00000\n")

    def _parse_gro(self, gro_file: str) -> dict:
        """Parse GRO file and extract atoms and coordinates"""
        atoms = []
        with open(gro_file, "r") as f:
            lines = f.readlines()

        # Skip title line
        # Second line is atom count
        num_atoms = int(lines[1].strip())

        # Parse atom lines (format: %5d%-5s%5s%5d%8.3f%8.3f%8.3f)
        for i in range(num_atoms):
            line = lines[2 + i]
            # Skip box line (starts with space and has only numbers)
            if i == num_atoms - 1 and len(line.strip().split()) <= 3:
                break

            # More robust parsing using split
            parts = line.split()
            if len(parts) < 7:
                continue

            # GRO format: residue_number residue_name atom_name atom_number x y z
            # parts: [0]=residue_number, [1]=residue_name, [2]=atom_name, [3]=atom_number, [4]=x, [5]=y, [6]=z
            atom_name = parts[2].strip()
            try:
                x = float(parts[4])
                y = float(parts[5])
                z = float(parts[6])
                atoms.append({"name": atom_name, "coord": [x, y, z]})
            except (ValueError, IndexError):
                # Try fixed column format as fallback
                try:
                    atom_name = line[5:10].strip()
                    x = float(line[20:28].strip())
                    y = float(line[28:36].strip())
                    z = float(line[36:44].strip())
                    atoms.append({"name": atom_name, "coord": [x, y, z]})
                except (ValueError, IndexError):
                    continue

        return {"atoms": atoms}

    def _read_box_size(self, gro_file: str) -> tuple:
        """Read box size from GRO file last line"""
        with open(gro_file, "r") as f:
            lines = f.readlines()
            box_line = lines[-1].strip().split()
            return tuple(float(x) for x in box_line[:3])

    def _build_standard_system(self, output_dir: str, use_protein: bool, box_size: tuple = None):
        """Build standard MD system (bound or unbound)"""
        os.makedirs(output_dir, exist_ok=True)

        original_fep = self.fep_mode
        original_output = self.output_dir
        original_protein = self.protein_path
        original_ligand_paths = self.ligand_paths
        original_lig_ff_dirs = list(self.lig_ff_dirs)

        # Convert ligand paths to absolute paths BEFORE changing output_dir
        # to avoid path resolution issues
        if isinstance(self.ligand_paths, str):
            self.ligand_paths = os.path.abspath(self.ligand_paths)
        elif isinstance(self.ligand_paths, list):
            self.ligand_paths = [os.path.abspath(p) for p in self.ligand_paths]

        # Also convert protein path to absolute
        if self.protein_path:
            self.protein_path = os.path.abspath(self.protein_path)

        self.fep_mode = False
        self.output_dir = output_dir

        # Update sub-component output directories
        from pathlib import Path

        self.system_builder.output_dir = Path(output_dir)
        self.system_builder.model_dir = Path(output_dir) / "GMX_PROLIG_MD"
        self.system_builder.model_dir.mkdir(exist_ok=True)
        self.mdp_generator.output_dir = output_dir

        if use_protein:
            # Normal protein+ligand system
            self.run_normal()
        else:
            # Ligand-only system: skip protein processing
            self._build_ligand_only_system(output_dir, box_size=box_size)

        # Restore
        self.fep_mode = original_fep
        self.output_dir = original_output
        self.protein_path = original_protein
        self.ligand_paths = original_ligand_paths
        self.lig_ff_dirs = original_lig_ff_dirs
        self.system_builder.output_dir = Path(original_output)
        self.system_builder.model_dir = Path(original_output) / "GMX_PROLIG_MD"
        self.mdp_generator.output_dir = original_output

    def _build_ligand_only_system(self, output_dir: str, box_size: tuple = None):
        """Build ligand-only system (no protein) for unbound FEP leg.

        Parameters
        ----------
        output_dir : str
            Output directory for the system.
        box_size : tuple, optional
            Box size (x, y, z) in nm. If provided, uses this exact box size.
            If None, uses box_distance parameter to create box around ligand.
        """
        from pathlib import Path
        import shutil

        print_step(1, 5, "Generating ligand force field for unbound system")

        original_output = self.output_dir
        self.output_dir = output_dir
        self.generate_ligand_forcefield()
        self.output_dir = original_output

        model_dir = Path(output_dir) / "GMX_PROLIG_MD"
        model_dir.mkdir(exist_ok=True, parents=True)

        ligand_ff_dir = Path(self.lig_ff_dirs[0])
        ligand_gro = ligand_ff_dir / "LIG.gro"
        if not ligand_gro.exists():
            raise FileNotFoundError(f"Ligand structure not found: {ligand_gro}")

        print_step(2, 5, "Creating ligand-only topology")
        topol_path = model_dir / "topol.top"
        ff_name = self.forcefield["name"]
        water_model = self.water_model["name"]
        ligand_rel_dir = f"../{ligand_ff_dir.name}"

        with open(topol_path, "w") as f:
            f.write("; Ligand-only topology for FEP unbound leg\n")
            f.write(f'#include "{ff_name}.ff/forcefield.itp"\n')
            f.write(f'#include "{ligand_rel_dir}/atomtypes_LIG.itp"\n')
            f.write(f'#include "{ligand_rel_dir}/LIG.itp"\n')
            f.write(f'#include "{ff_name}.ff/{water_model}.itp"\n')
            f.write(f'#include "{ff_name}.ff/ions.itp"\n\n')
            f.write("[ system ]\n")
            f.write("Ligand in water\n\n")
            f.write("[ molecules ]\n")
            f.write("LIG    1\n")

        shutil.copy(ligand_gro, model_dir / "lig.gro")

        print_step(3, 5, "Creating simulation box")
        boxed_gro = model_dir / "lig_newbox.gro"
        if box_size:
            print(f"  Using specified box size: {box_size[0]:.3f} x {box_size[1]:.3f} x {box_size[2]:.3f} nm")
            self.system_builder._run_command(
                [
                    self.system_builder.gmx_command,
                    "editconf",
                    "-f",
                    str(model_dir / "lig.gro"),
                    "-o",
                    str(boxed_gro),
                    "-bt",
                    "cubic",
                    "-box",
                    str(box_size[0]),
                    str(box_size[1]),
                    str(box_size[2]),
                    "-c",
                ],
                work_dir=str(model_dir),
            )
        else:
            box_distance = self.config.get("system", {}).get("box_distance", 1.5)
            self.system_builder._run_command(
                [
                    self.system_builder.gmx_command,
                    "editconf",
                    "-f",
                    str(model_dir / "lig.gro"),
                    "-o",
                    str(boxed_gro),
                    "-bt",
                    "cubic",
                    "-d",
                    str(box_distance),
                    "-c",
                ],
                work_dir=str(model_dir),
            )

        print_step(4, 5, "Solvating system")
        solvated_gro = self.system_builder._solvate(str(boxed_gro), str(topol_path))

        print_step(5, 5, "Adding ions")
        self.system_builder._add_ions(solvated_gro, str(topol_path))

        print_success(f"Ligand-only system built in {model_dir}")

    def _resolve_generated_ligand_ff_dir(self, output_dir: str) -> str:
        """Return the generated ligand force-field directory for the current ligand FF."""
        ff_dirs = [path for path in Path(output_dir).glob("LIG.*") if path.is_dir()]
        if len(ff_dirs) != 1:
            raise FileNotFoundError(
                f"Expected exactly one generated ligand FF dir in {output_dir}, found {len(ff_dirs)}"
            )
        return str(ff_dirs[0])

    def _generate_mutant_ligand_ff(self, mutant_ligand: str, output_dir: str):
        """Generate force field for mutant ligand"""
        # Reuse existing force field generation logic

        # Temporarily save current state
        original_ligands = self.ligand_paths
        original_output = self.output_dir
        original_lig_ff_dirs = self.lig_ff_dirs

        # Set up for mutant ligand
        self.ligand_paths = [mutant_ligand]
        self.output_dir = output_dir

        # Generate FF
        self.generate_ligand_forcefield()

        # Restore state
        self.ligand_paths = original_ligands
        self.output_dir = original_output
        self.lig_ff_dirs = original_lig_ff_dirs

    def run_rest2(self):
        """Run the REST2 workflow: build standard MD first, then convert to REST2"""
        print_header("PRISM REST2 Builder Workflow")

        try:
            # Phase 1: Build standard MD system via normal workflow
            print_step(1, 2, "Building standard MD system (normal workflow)")
            self.run_normal()

            # Phase 2: Convert GMX_PROLIG_MD -> GMX_PROLIG_REST2
            print_step(2, 2, "Converting to REST2 replica exchange setup")

            from ..rest2 import REST2Workflow

            md_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
            rest2_dir = os.path.join(self.output_dir, "GMX_PROLIG_REST2")

            # Ligand residue name: PRISM always standardizes to 'LIG' (single)
            # or 'LIG_1','LIG_2',... (multi-ligand)
            if self.ligand_count == 1:
                lig_name = "LIG"
            else:
                lig_name = [f"LIG_{i}" for i in range(1, self.ligand_count + 1)]

            workflow = REST2Workflow(
                md_dir=md_dir,
                output_dir=rest2_dir,
                t_ref=self.t_ref,
                t_max=self.t_max,
                n_replicas=self.n_replicas,
                cutoff=self.rest2_cutoff,
                lig_name=lig_name,
                gmx=self.gromacs_env.gmx_command,
            )

            workflow.run()

            print_header("REST2 Setup Complete!")
            print(f"\n  Standard MD files: {path(md_dir)}")
            print(f"  REST2 files:       {path(rest2_dir)}")
            print(f"\n  To run REST2 simulations:")
            print(f"  1. cd {rest2_dir}")
            print(f"  2. bash rest2_run.sh")

            return self.output_dir

        except Exception as e:
            print_error(f"REST2 workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise

    def run_mmpbsa(self):
        """Run MMPBSA workflow: build standard MD first, then set up for MMPBSA.

        Two sub-modes:
        - Single-frame (default): EM->NVT->NPT->gmx_MMPBSA on equilibrated structure.
        - Trajectory-based (--mmpbsa-traj N): EM->NVT->NPT->Production MD->gmx_MMPBSA.
        """
        try:
            from ..mmpbsa.generator import MMPBSAGenerator
        except ImportError:
            from prism.mmpbsa.generator import MMPBSAGenerator

        is_trajectory = self.mmpbsa_traj_ns is not None
        mode_label = "Trajectory" if is_trajectory else "Single-frame"
        print_header(f"PRISM MM/PBSA Builder Workflow ({mode_label})")

        try:
            # For trajectory mode, override production MD length
            if is_trajectory:
                self.config["simulation"]["production_time_ns"] = self.mmpbsa_traj_ns

            # Phase 1: Build standard MD system via normal workflow
            print_step(1, 2, "Building standard MD system (normal workflow)")
            self.run_normal()

            # Phase 2: Convert GMX_PROLIG_MD -> GMX_PROLIG_MMPBSA
            print_step(2, 2, "Setting up MM/PBSA calculation files")
            md_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
            mmpbsa_dir = os.path.join(self.output_dir, "GMX_PROLIG_MMPBSA")

            # Handle existing MMPBSA directory from previous runs
            if os.path.exists(mmpbsa_dir):
                print_warning(f"Removing existing {os.path.basename(mmpbsa_dir)}/ from previous run")
                shutil.rmtree(mmpbsa_dir)

            os.rename(md_dir, mmpbsa_dir)

            # Clean up files that are irrelevant for MMPBSA workflow
            localrun_path = os.path.join(mmpbsa_dir, "localrun.sh")
            if os.path.exists(localrun_path):
                os.remove(localrun_path)

            # In single-frame mode, remove the unnecessary md.mdp
            if not is_trajectory:
                md_mdp = os.path.join(self.output_dir, "mdps", "md.mdp")
                if os.path.exists(md_mdp):
                    os.remove(md_mdp)

            # Note: -lm (ligand mol2) is NOT needed for gmx_MMPBSA when
            # -cp topol.top is provided, because the GROMACS topology already
            # contains all ligand parameters. gmx_MMPBSA converts internally
            # via parmed. Omitting -lm avoids parmchk2 failures caused by
            # atom type mismatches in the original mol2 file.

            # Generate mmpbsa.in and mmpbsa_run.sh via MMPBSAGenerator
            generator = MMPBSAGenerator(
                mmpbsa_dir=mmpbsa_dir,
                protein_ff_name=self.forcefield["name"],
                ligand_forcefield=self.ligand_forcefield,
                temperature=self.config["simulation"]["temperature"],
                traj_ns=self.mmpbsa_traj_ns,
                trajectory_interval_ps=self.config.get("output", {}).get("trajectory_interval_ps", 500),
                gmx2amber=self.gmx2amber,
            )
            generator.generate_input()
            generator.generate_script()

            mode_desc = "trajectory" if is_trajectory else "single-frame"
            print_success(f"Generated mmpbsa.in ({mode_desc} mode)")
            print_success(f"Generated mmpbsa_run.sh ({mode_desc} mode)")

            # Print completion info
            print_header("MM/PBSA Setup Complete!")
            print(f"\n  MMPBSA files: {path(mmpbsa_dir)}")
            if is_trajectory:
                print(f"  Sub-mode: Trajectory-based ({self.mmpbsa_traj_ns} ns production MD)")
            else:
                print(f"  Sub-mode: Single-frame (no production MD)")
            print(f"\n  To run MM/PBSA workflow:")
            print(f"  1. cd {mmpbsa_dir}")
            print(f"  2. bash mmpbsa_run.sh")

            return self.output_dir

        except Exception as e:
            print_error(f"MMPBSA workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise

    def run_normal(self):
        """Run the normal (non-PMF) workflow"""
        print_header("PRISM Builder Workflow")

        try:
            # Save configuration for reference
            self.save_config()

            # Step 1: Generate ligand force field
            print_step(1, 6, f"Generating ligand force field ({self.ligand_forcefield.upper()})")
            self.generate_ligand_forcefield()
            print_success(f"Ligand force field generated")

            # Step 2: Clean protein
            print_step(2, 6, "Cleaning protein structure")
            cleaned_protein = self.clean_protein()
            print_success("Protein structure cleaned")

            # Step 3: Build model
            print_step(3, 6, "Building GROMACS system")
            model_dir = self.build_model(cleaned_protein)
            if not model_dir:
                raise RuntimeError("Failed to build model")
            print_success("GROMACS system built")

            # Step 4: Generate MDP files
            print_step(4, 6, "Generating MD parameter files")
            self.generate_mdp_files()
            print_success("MDP files generated")

            # Step 5: Cleanup
            print_step(5, 6, "Cleaning up temporary files")
            self.cleanup()
            print_success("Cleanup completed")

            # Step 6: Generate local run script
            print_step(6, 6, "Generating simulation run script")
            script_path = self.generate_localrun_script()
            print_success("Run script generated")

            print_header("Workflow Complete!")
            print(f"\nOutput files are in: {path(self.output_dir)}")
            print(f"MD system files are in: {path(os.path.join(self.output_dir, 'GMX_PROLIG_MD'))}")
            print(f"MDP files are in: {path(self.mdp_generator.mdp_dir)}")
            print(f"Configuration saved in: {path(os.path.join(self.output_dir, 'prism_config.yaml'))}")
            print(f"\nProtein force field used: {number(self.forcefield['name'])}")
            print(f"Ligand force field used: {number(self.ligand_forcefield.upper())}")
            print(f"Water model used: {number(self.water_model['name'])}")

            if script_path:
                gmx_md_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
                print_header("Ready to Run MD Simulations!")
                print(f"\nTo run the MD workflow:")
                print(f"  1. Navigate to the MD directory:")
                print(f"     {path(f'cd {gmx_md_dir}')}")
                print(f"  2. Execute the simulation script:")
                print(f"     {number('bash localrun.sh')}")
                print(f"\nThe script will run:")
                print(f"  - Energy Minimization (EM)")
                print(f"  - NVT Equilibration")
                print(f"  - NPT Equilibration")
                print(f"  - Production MD")
                print(f"\nNote: Adjust GPU and thread settings in localrun.sh as needed")

            return self.output_dir

        except Exception as e:
            print_error(f"Workflow failed: {e}")
            import traceback

            traceback.print_exc()
            raise
