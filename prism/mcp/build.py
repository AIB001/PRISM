"""System building tools (MD, PMF, REST2, MM/PBSA)."""

import os
import json
import traceback
from typing import Optional

from ._common import _StdoutToStderr, _ensure_prism_importable, logger


def register(mcp):

    @mcp.tool()
    def build_system(
        protein_path: str,
        ligand_paths: str,
        output_dir: str = "prism_output",
        ligand_forcefield: str = "gaff2",
        forcefield: str = "amber14sb",
        water_model: str = "tip3p",
        protonation: str = "gromacs",
        temperature: float = 310.0,
        ph: float = 7.0,
        production_ns: float = 500.0,
        box_distance: float = 1.5,
        box_shape: str = "cubic",
        salt_concentration: float = 0.15,
        ligand_charge: int = 0,
        gaussian_method: Optional[str] = None,
        do_optimization: bool = False,
        overwrite: bool = False,
        pressure: float = 1.0,
        dt: float = 0.002,
        nvt_ps: float = 500.0,
        npt_ps: float = 500.0,
        nterm_met: str = "keep",
        em_tolerance: float = 200.0,
        em_steps: int = 10000,
    ) -> str:
        """Build a complete protein-ligand system for GROMACS molecular dynamics simulation.

        This is the main building tool. It takes a protein PDB and one or more ligand
        files, then performs the full workflow:
          1. Generate ligand force field parameters
          2. Clean and prepare the protein structure
          3. Build GROMACS system (topology, solvation, ions)
          4. Generate MDP parameter files (em, nvt, npt, md)
          5. Generate localrun.sh script for GPU-accelerated simulation

        Recommended defaults: amber14sb protein force field + gaff2 ligand force field.

        If the user has Gaussian (g16) installed, they can set gaussian_method to 'hf'
        or 'dft' to use high-precision RESP charges instead of the default AM1-BCC.
        If g16 is not available, PRISM will generate Gaussian input files and a script
        for the user to run manually on a machine with Gaussian.

        After building, the user runs:
          cd <output_dir>/GMX_PROLIG_MD && bash localrun.sh

        Args:
            protein_path: Absolute path to the protein PDB file.
            ligand_paths: Absolute path to ligand file(s). For multiple ligands,
                separate with commas: "/path/lig1.mol2,/path/lig2.mol2"
            output_dir: Directory for all output files. Default: "prism_output".
            ligand_forcefield: Ligand force field. Default: "gaff2" (recommended).
                Options: "gaff", "gaff2", "openff", "opls", "mmff", "match", "hybrid".
                Note: "cgenff" requires manual --forcefield-path setup via CLI.
            forcefield: Protein force field. Default: "amber14sb" (recommended).
                Common: "amber99sb", "amber14sb", "amber99sb-ildn", "charmm27".
                For amber19sb, use water_model="opc".
            water_model: Water model. Default: "tip3p".
                Options: "tip3p", "tip4p", "spc", "spce", "opc" (for amber19sb).
            protonation: Protonation method. Default: "gromacs".
                "gromacs": Let pdb2gmx handle protonation states (default HIE).
                "propka": Use PROPKA pKa prediction for intelligent per-residue
                ionizable residue states (HID/HIE/HIP, ASH/GLH, LYN/CYM/TYH). Requires propka package.
            temperature: Simulation temperature in Kelvin. Default: 310.
            ph: pH for protonation state assignment. Default: 7.0.
            production_ns: Production MD length in nanoseconds. Default: 500.
            box_distance: Distance from protein to box edge in nm. Default: 1.5.
            box_shape: Box shape. Default: "cubic". Options: "cubic", "dodecahedron", "octahedron".
            salt_concentration: NaCl concentration in mol/L. Default: 0.15.
            ligand_charge: Net formal charge of the ligand. Default: 0.
            gaussian_method: Enable Gaussian RESP charge calculation. Default: None (disabled).
                Set to "hf" for HF/6-31G* or "dft" for B3LYP/6-31G*.
                If Gaussian (g16) is installed, charges are calculated automatically.
                Otherwise, input files and a run script are generated.
            do_optimization: Perform geometry optimization before ESP. Default: false.
                Only used when gaussian_method is set.
            overwrite: Overwrite existing output files. Default: false.
            pressure: Simulation pressure in bar. Default: 1.0.
            dt: Integration time step in ps. Default: 0.002.
            nvt_ps: NVT equilibration time in ps. Default: 500.0.
            npt_ps: NPT equilibration time in ps. Default: 500.0.
            nterm_met: N-terminal methionine handling. Default: "keep".
                "keep": Keep N-terminal MET as-is.
                "drop": Remove N-terminal MET.
                "auto": Drop for CHARMM36 force fields, keep otherwise.
            em_tolerance: Energy minimization convergence tolerance in kJ/mol/nm. Default: 200.0.
            em_steps: Maximum number of energy minimization steps. Default: 10000.
        """
        logger.info(f"build_system: {protein_path} + {ligand_paths} -> {output_dir}")

        # --- Parse ligand paths (comma-separated string -> list) ---
        lig_list = [p.strip() for p in ligand_paths.split(",") if p.strip()]

        # --- Input validation ---
        errors = []
        if not os.path.isabs(protein_path):
            errors.append(f"protein_path must be an absolute path, got: {protein_path}")
        elif not os.path.exists(protein_path):
            errors.append(f"Protein file not found: {protein_path}")

        for lp in lig_list:
            if not os.path.isabs(lp):
                errors.append(f"ligand_path must be an absolute path, got: {lp}")
            elif not os.path.exists(lp):
                errors.append(f"Ligand file not found: {lp}")

        valid_lff = ["gaff", "gaff2", "openff", "opls", "mmff", "match", "hybrid"]
        if ligand_forcefield not in valid_lff:
            errors.append(f"Invalid ligand_forcefield '{ligand_forcefield}'. Options: {valid_lff}")

        if gaussian_method and gaussian_method not in ("hf", "dft"):
            errors.append(f"gaussian_method must be 'hf' or 'dft', got: {gaussian_method}")

        if nterm_met not in ("keep", "drop", "auto"):
            errors.append(f"nterm_met must be 'keep', 'drop', or 'auto', got: {nterm_met}")

        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        # --- Run PRISM builder ---
        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.builder import PRISMBuilder

                # Single ligand: pass string; multiple: pass list
                lig_input = lig_list[0] if len(lig_list) == 1 else lig_list

                builder = PRISMBuilder(
                    protein_path=protein_path,
                    ligand_paths=lig_input,
                    output_dir=output_dir,
                    ligand_forcefield=ligand_forcefield,
                    forcefield=forcefield,
                    water_model=water_model,
                    overwrite=overwrite,
                    gaussian_method=gaussian_method,
                    do_optimization=do_optimization,
                )

                # Apply simulation parameters
                builder.config["simulation"]["temperature"] = temperature
                builder.config["simulation"]["pH"] = ph
                builder.config["simulation"]["production_time_ns"] = production_ns
                builder.config["simulation"]["ligand_charge"] = ligand_charge
                builder.config["ligand_forcefield"]["charge"] = ligand_charge
                builder.config["box"]["distance"] = box_distance
                builder.config["box"]["shape"] = box_shape
                builder.config["ions"]["concentration"] = salt_concentration

                # Apply protonation method
                builder.config.setdefault("protonation", {})["method"] = protonation

                # Apply new parameters
                builder.config["simulation"]["pressure"] = pressure
                builder.config["simulation"]["dt"] = dt
                builder.config["simulation"]["equilibration_nvt_time_ps"] = nvt_ps
                builder.config["simulation"]["equilibration_npt_time_ps"] = npt_ps
                builder.config.setdefault("protein_preparation", {})["nterm_met"] = nterm_met
                builder.config.setdefault("energy_minimization", {})["emtol"] = em_tolerance
                builder.config["energy_minimization"]["nsteps"] = em_steps

                result_dir = builder.run()

            # --- Collect output info ---
            gmx_dir = os.path.join(result_dir, "GMX_PROLIG_MD")
            mdp_dir = os.path.join(result_dir, "mdps")
            files = {}
            if os.path.isdir(gmx_dir):
                files["system_coordinates"] = os.path.join(gmx_dir, "solv_ions.gro")
                files["topology"] = os.path.join(gmx_dir, "topol.top")
                files["run_script"] = os.path.join(gmx_dir, "localrun.sh")
            if os.path.isdir(mdp_dir):
                for name in ("em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"):
                    p = os.path.join(mdp_dir, name)
                    if os.path.exists(p):
                        files[name] = p

            return json.dumps(
                {
                    "success": True,
                    "output_dir": result_dir,
                    "gmx_dir": gmx_dir,
                    "message": (
                        f"System built successfully!\n"
                        f"To run the simulation:\n"
                        f"  cd {gmx_dir}\n"
                        f"  bash localrun.sh"
                    ),
                    "files": files,
                    "parameters": {
                        "protein_forcefield": forcefield,
                        "ligand_forcefield": ligand_forcefield,
                        "water_model": water_model,
                        "protonation": protonation,
                        "temperature_K": temperature,
                        "pH": ph,
                        "production_ns": production_ns,
                        "box_distance_nm": box_distance,
                        "salt_concentration_M": salt_concentration,
                        "gaussian_method": gaussian_method,
                        "pressure_bar": pressure,
                        "dt_ps": dt,
                        "nvt_ps": nvt_ps,
                        "npt_ps": npt_ps,
                        "nterm_met": nterm_met,
                        "em_tolerance": em_tolerance,
                        "em_steps": em_steps,
                    },
                },
                indent=2,
            )

        except Exception as e:
            logger.error(f"build_system failed: {e}\n{traceback.format_exc()}")
            return json.dumps(
                {
                    "success": False,
                    "error": str(e),
                    "traceback": traceback.format_exc(),
                },
                indent=2,
            )

    @mcp.tool()
    def build_pmf_system(
        protein_path: str,
        ligand_path: str,
        output_dir: str = "prism_pmf_output",
        ligand_forcefield: str = "gaff2",
        forcefield: str = "amber14sb",
        water_model: str = "tip3p",
        protonation: str = "gromacs",
        temperature: float = 310.0,
        ph: float = 7.0,
        salt_concentration: float = 0.15,
        ligand_charge: int = 0,
        gaussian_method: Optional[str] = None,
        do_optimization: bool = False,
        nterm_met: str = "keep",
        box_extension_z: float = 2.0,
        umbrella_time_ns: float = 10.0,
        umbrella_spacing: float = 0.12,
        overwrite: bool = False,
    ) -> str:
        """Build a system for PMF (Potential of Mean Force) calculation.

        PMF measures the binding free energy profile of a protein-ligand complex
        using steered molecular dynamics (SMD) and umbrella sampling with WHAM.

        Recommended: amber14sb + gaff2 (same as standard MD).

        The workflow:
        1. Generate ligand force field (gaff2 recommended)
        2. Clean protein structure
        3. Align complex so the unbinding pull vector is along the Z-axis
        4. Build GROMACS system with extended Z-box for pulling space
        5. Generate index file with pull/reference/freeze groups
        6. Generate SMD and umbrella sampling MDP files and run scripts

        After building, the user runs:
          cd <output_dir>/GMX_PROLIG_PMF
          bash smd_run.sh           # Step 1: Steered MD (EM -> NVT -> NPT -> pulling)
          bash umbrella_run.sh      # Step 2: Umbrella sampling + WHAM analysis

        Args:
            protein_path: Absolute path to protein PDB file.
            ligand_path: Absolute path to ligand file (MOL2 or SDF).
            output_dir: Output directory. Default: "prism_pmf_output".
            ligand_forcefield: Ligand force field. Default: "gaff2" (recommended).
            forcefield: Protein force field. Default: "amber14sb" (recommended).
            water_model: Water model. Default: "tip3p".
            protonation: Protonation method. Default: "gromacs".
                "gromacs": Let pdb2gmx handle protonation states (default HIE).
                "propka": Use PROPKA pKa prediction for intelligent per-residue
                ionizable residue states. Requires propka package.
            temperature: Simulation temperature in Kelvin. Default: 310.
            ph: pH for protonation state assignment. Default: 7.0.
            salt_concentration: NaCl concentration in mol/L. Default: 0.15.
            ligand_charge: Net formal charge of the ligand. Default: 0.
            gaussian_method: Enable Gaussian RESP charge calculation. Default: None.
                Set to "hf" for HF/6-31G* or "dft" for B3LYP/6-31G*.
            do_optimization: Perform geometry optimization before ESP. Default: false.
            nterm_met: N-terminal methionine handling. Default: "keep".
            box_extension_z: Extra Z-axis length (nm) for pulling space. Default: 2.0.
            umbrella_time_ns: Simulation time per umbrella window in ns. Default: 10.0.
            umbrella_spacing: Distance between umbrella windows in nm. Default: 0.12.
            overwrite: Overwrite existing files. Default: false.
        """
        logger.info(f"build_pmf_system: {protein_path} + {ligand_path}")

        errors = []
        if not os.path.isabs(protein_path):
            errors.append(f"protein_path must be an absolute path, got: {protein_path}")
        elif not os.path.exists(protein_path):
            errors.append(f"Protein file not found: {protein_path}")
        if not os.path.isabs(ligand_path):
            errors.append(f"ligand_path must be an absolute path, got: {ligand_path}")
        elif not os.path.exists(ligand_path):
            errors.append(f"Ligand file not found: {ligand_path}")
        if gaussian_method and gaussian_method not in ("hf", "dft"):
            errors.append(f"gaussian_method must be 'hf' or 'dft', got: {gaussian_method}")
        if nterm_met not in ("keep", "drop", "auto"):
            errors.append(f"nterm_met must be 'keep', 'drop', or 'auto', got: {nterm_met}")
        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.builder import PRISMBuilder

                builder = PRISMBuilder(
                    protein_path=protein_path,
                    ligand_paths=ligand_path,
                    output_dir=output_dir,
                    ligand_forcefield=ligand_forcefield,
                    forcefield=forcefield,
                    water_model=water_model,
                    overwrite=overwrite,
                    gaussian_method=gaussian_method,
                    do_optimization=do_optimization,
                    pmf_mode=True,
                    box_extension=(0.0, 0.0, box_extension_z),
                )

                # Apply shared simulation parameters
                builder.config["simulation"]["temperature"] = temperature
                builder.config["simulation"]["pH"] = ph
                builder.config["simulation"]["ligand_charge"] = ligand_charge
                builder.config["ligand_forcefield"]["charge"] = ligand_charge
                builder.config["ions"]["concentration"] = salt_concentration
                builder.config.setdefault("protonation", {})["method"] = protonation
                builder.config.setdefault("protein_preparation", {})["nterm_met"] = nterm_met

                # Apply PMF-specific parameters
                builder.config.setdefault("pmf", {})
                builder.config["pmf"]["umbrella_time_ns"] = umbrella_time_ns
                builder.config["pmf"]["umbrella_spacing"] = umbrella_spacing

                result_dir = builder.run()

            pmf_dir = os.path.join(result_dir, "GMX_PROLIG_PMF")
            return json.dumps(
                {
                    "success": True,
                    "output_dir": result_dir,
                    "pmf_dir": pmf_dir,
                    "message": (
                        f"PMF system built successfully!\n"
                        f"To run the PMF workflow:\n"
                        f"  cd {pmf_dir}\n"
                        f"  bash smd_run.sh        # Step 1: Steered MD\n"
                        f"  bash umbrella_run.sh   # Step 2: Umbrella sampling + WHAM"
                    ),
                    "parameters": {
                        "protein_forcefield": forcefield,
                        "ligand_forcefield": ligand_forcefield,
                        "water_model": water_model,
                        "protonation": protonation,
                        "temperature_K": temperature,
                        "salt_concentration_M": salt_concentration,
                        "gaussian_method": gaussian_method,
                        "box_extension_z_nm": box_extension_z,
                        "umbrella_time_ns": umbrella_time_ns,
                        "umbrella_spacing_nm": umbrella_spacing,
                    },
                },
                indent=2,
            )

        except Exception as e:
            logger.error(f"build_pmf_system failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e)}, indent=2)

    @mcp.tool()
    def build_rest2_system(
        protein_path: str,
        ligand_path: str,
        output_dir: str = "prism_rest2_output",
        ligand_forcefield: str = "gaff2",
        forcefield: str = "amber14sb",
        water_model: str = "tip3p",
        protonation: str = "gromacs",
        ph: float = 7.0,
        salt_concentration: float = 0.15,
        ligand_charge: int = 0,
        gaussian_method: Optional[str] = None,
        do_optimization: bool = False,
        nterm_met: str = "keep",
        t_ref: float = 310.0,
        t_max: float = 450.0,
        n_replicas: int = 16,
        rest2_cutoff: float = 0.5,
        overwrite: bool = False,
    ) -> str:
        """Build a system for REST2 (Replica Exchange with Solute Tempering v2).

        REST2 enhances conformational sampling of the protein-ligand binding pocket
        by running multiple replicas at different effective temperatures for the
        solute (protein + ligand), while keeping the solvent at reference temperature.
        This requires significantly more computational resources than standard MD.

        Recommended: amber14sb + gaff2.

        The workflow:
        1. Build a standard MD system (same as build_system)
        2. Convert to REST2 with scaled topologies for each replica

        After building, run:
          cd <output_dir>/GMX_PROLIG_REST2 && bash rest2_run.sh

        Args:
            protein_path: Absolute path to protein PDB file.
            ligand_path: Absolute path to ligand file (MOL2 or SDF).
            output_dir: Output directory. Default: "prism_rest2_output".
            ligand_forcefield: Ligand force field. Default: "gaff2" (recommended).
            forcefield: Protein force field. Default: "amber14sb" (recommended).
            water_model: Water model. Default: "tip3p".
            protonation: Protonation method. Default: "gromacs".
                "gromacs": Let pdb2gmx handle protonation states (default HIE).
                "propka": Use PROPKA pKa prediction for intelligent per-residue
                ionizable residue states. Requires propka package.
            ph: pH for protonation state assignment. Default: 7.0.
            salt_concentration: NaCl concentration in mol/L. Default: 0.15.
            ligand_charge: Net formal charge of the ligand. Default: 0.
            gaussian_method: Enable Gaussian RESP charge calculation. Default: None.
                Set to "hf" for HF/6-31G* or "dft" for B3LYP/6-31G*.
            do_optimization: Perform geometry optimization before ESP. Default: false.
            nterm_met: N-terminal methionine handling. Default: "keep".
            t_ref: Reference (physical) temperature in Kelvin. Default: 310.
            t_max: Maximum effective temperature in Kelvin. Default: 450.
            n_replicas: Number of REST2 replicas. Default: 16.
                More replicas = better exchange rate but more compute cost.
            rest2_cutoff: Distance cutoff (nm) for identifying binding pocket residues
                to include in the solute region. Default: 0.5.
            overwrite: Overwrite existing files. Default: false.
        """
        logger.info(f"build_rest2_system: {protein_path} + {ligand_path}")

        errors = []
        if not os.path.isabs(protein_path):
            errors.append(f"protein_path must be an absolute path, got: {protein_path}")
        elif not os.path.exists(protein_path):
            errors.append(f"Protein file not found: {protein_path}")
        if not os.path.isabs(ligand_path):
            errors.append(f"ligand_path must be an absolute path, got: {ligand_path}")
        elif not os.path.exists(ligand_path):
            errors.append(f"Ligand file not found: {ligand_path}")
        if gaussian_method and gaussian_method not in ("hf", "dft"):
            errors.append(f"gaussian_method must be 'hf' or 'dft', got: {gaussian_method}")
        if nterm_met not in ("keep", "drop", "auto"):
            errors.append(f"nterm_met must be 'keep', 'drop', or 'auto', got: {nterm_met}")
        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.builder import PRISMBuilder

                builder = PRISMBuilder(
                    protein_path=protein_path,
                    ligand_paths=ligand_path,
                    output_dir=output_dir,
                    ligand_forcefield=ligand_forcefield,
                    forcefield=forcefield,
                    water_model=water_model,
                    overwrite=overwrite,
                    gaussian_method=gaussian_method,
                    do_optimization=do_optimization,
                    rest2_mode=True,
                    t_ref=t_ref,
                    t_max=t_max,
                    n_replicas=n_replicas,
                    rest2_cutoff=rest2_cutoff,
                )

                # Apply shared simulation parameters
                builder.config["simulation"]["pH"] = ph
                builder.config["simulation"]["ligand_charge"] = ligand_charge
                builder.config["ligand_forcefield"]["charge"] = ligand_charge
                builder.config["ions"]["concentration"] = salt_concentration
                builder.config.setdefault("protonation", {})["method"] = protonation
                builder.config.setdefault("protein_preparation", {})["nterm_met"] = nterm_met

                result_dir = builder.run()

            rest2_dir = os.path.join(result_dir, "GMX_PROLIG_REST2")
            return json.dumps(
                {
                    "success": True,
                    "output_dir": result_dir,
                    "rest2_dir": rest2_dir,
                    "message": (
                        f"REST2 system built with {n_replicas} replicas "
                        f"(T: {t_ref}K - {t_max}K).\n"
                        f"To run:\n"
                        f"  cd {rest2_dir}\n"
                        f"  bash rest2_run.sh"
                    ),
                    "parameters": {
                        "protein_forcefield": forcefield,
                        "ligand_forcefield": ligand_forcefield,
                        "water_model": water_model,
                        "protonation": protonation,
                        "salt_concentration_M": salt_concentration,
                        "gaussian_method": gaussian_method,
                        "t_ref_K": t_ref,
                        "t_max_K": t_max,
                        "n_replicas": n_replicas,
                        "rest2_cutoff_nm": rest2_cutoff,
                    },
                },
                indent=2,
            )

        except Exception as e:
            logger.error(f"build_rest2_system failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e)}, indent=2)

    @mcp.tool()
    def build_mmpbsa_system(
        protein_path: str,
        ligand_path: str,
        output_dir: str = "prism_mmpbsa_output",
        ligand_forcefield: str = "gaff2",
        forcefield: str = "amber14sb",
        water_model: str = "tip3p",
        protonation: str = "gromacs",
        temperature: float = 310.0,
        ph: float = 7.0,
        salt_concentration: float = 0.15,
        ligand_charge: int = 0,
        gaussian_method: Optional[str] = None,
        do_optimization: bool = False,
        nterm_met: str = "keep",
        mmpbsa_traj_ns: Optional[float] = None,
        gmx2amber: bool = False,
        overwrite: bool = False,
    ) -> str:
        """Build a system for MM/PBSA binding free energy calculation.

        MM/PBSA estimates protein-ligand binding affinity from MD snapshots.
        Two sub-modes:
        1. Single-frame (default, mmpbsa_traj_ns=None):
           EM -> NVT -> NPT -> gmx_MMPBSA on equilibrated structure. Fast.
        2. Trajectory-based (set mmpbsa_traj_ns to a value):
           EM -> NVT -> NPT -> Production MD -> gmx_MMPBSA on trajectory frames.
           More accurate but requires longer simulation.

        Recommended: amber14sb + gaff2.

        After building, run:
          cd <output_dir>/GMX_PROLIG_MMPBSA && bash mmpbsa_run.sh

        Args:
            protein_path: Absolute path to protein PDB file.
            ligand_path: Absolute path to ligand file (MOL2 or SDF).
            output_dir: Output directory. Default: "prism_mmpbsa_output".
            ligand_forcefield: Ligand force field. Default: "gaff2" (recommended).
            forcefield: Protein force field. Default: "amber14sb" (recommended).
            water_model: Water model. Default: "tip3p".
            protonation: Protonation method. Default: "gromacs".
                "gromacs": Let pdb2gmx handle protonation states (default HIE).
                "propka": Use PROPKA pKa prediction for intelligent per-residue
                ionizable residue states. Requires propka package.
            temperature: Simulation temperature in Kelvin. Default: 310.
            ph: pH for protonation state assignment. Default: 7.0.
            salt_concentration: NaCl concentration in mol/L. Default: 0.15.
            ligand_charge: Net formal charge of the ligand. Default: 0.
            gaussian_method: Enable Gaussian RESP charge calculation. Default: None.
                Set to "hf" for HF/6-31G* or "dft" for B3LYP/6-31G*.
            do_optimization: Perform geometry optimization before ESP. Default: false.
            nterm_met: N-terminal methionine handling. Default: "keep".
            mmpbsa_traj_ns: Production MD length in ns for trajectory-based mode.
                If None (default), uses single-frame mode (no production MD).
            gmx2amber: Use AMBER MMPBSA.py via parmed instead of gmx_MMPBSA.
                Requires AmberTools with parmed. Default: false.
            overwrite: Overwrite existing files. Default: false.
        """
        logger.info(f"build_mmpbsa_system: {protein_path} + {ligand_path}")

        errors = []
        if not os.path.isabs(protein_path):
            errors.append(f"protein_path must be an absolute path, got: {protein_path}")
        elif not os.path.exists(protein_path):
            errors.append(f"Protein file not found: {protein_path}")
        if not os.path.isabs(ligand_path):
            errors.append(f"ligand_path must be an absolute path, got: {ligand_path}")
        elif not os.path.exists(ligand_path):
            errors.append(f"Ligand file not found: {ligand_path}")
        if gaussian_method and gaussian_method not in ("hf", "dft"):
            errors.append(f"gaussian_method must be 'hf' or 'dft', got: {gaussian_method}")
        if nterm_met not in ("keep", "drop", "auto"):
            errors.append(f"nterm_met must be 'keep', 'drop', or 'auto', got: {nterm_met}")
        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.builder import PRISMBuilder

                builder = PRISMBuilder(
                    protein_path=protein_path,
                    ligand_paths=ligand_path,
                    output_dir=output_dir,
                    ligand_forcefield=ligand_forcefield,
                    forcefield=forcefield,
                    water_model=water_model,
                    overwrite=overwrite,
                    gaussian_method=gaussian_method,
                    do_optimization=do_optimization,
                    mmpbsa_mode=True,
                    mmpbsa_traj_ns=mmpbsa_traj_ns,
                    gmx2amber=gmx2amber,
                )

                # Apply shared simulation parameters
                builder.config["simulation"]["temperature"] = temperature
                builder.config["simulation"]["pH"] = ph
                builder.config["simulation"]["ligand_charge"] = ligand_charge
                builder.config["ligand_forcefield"]["charge"] = ligand_charge
                builder.config["ions"]["concentration"] = salt_concentration
                builder.config.setdefault("protonation", {})["method"] = protonation
                builder.config.setdefault("protein_preparation", {})["nterm_met"] = nterm_met

                result_dir = builder.run()

            mode = "trajectory" if mmpbsa_traj_ns else "single-frame"
            mmpbsa_dir = os.path.join(result_dir, "GMX_PROLIG_MMPBSA")
            return json.dumps(
                {
                    "success": True,
                    "output_dir": result_dir,
                    "mmpbsa_dir": mmpbsa_dir,
                    "mode": mode,
                    "message": (
                        f"MM/PBSA system built ({mode} mode).\n"
                        f"To run:\n"
                        f"  cd {mmpbsa_dir}\n"
                        f"  bash mmpbsa_run.sh"
                    ),
                    "parameters": {
                        "protein_forcefield": forcefield,
                        "ligand_forcefield": ligand_forcefield,
                        "water_model": water_model,
                        "protonation": protonation,
                        "temperature_K": temperature,
                        "salt_concentration_M": salt_concentration,
                        "gaussian_method": gaussian_method,
                        "mode": mode,
                        "mmpbsa_traj_ns": mmpbsa_traj_ns,
                        "gmx2amber": gmx2amber,
                    },
                },
                indent=2,
            )

        except Exception as e:
            logger.error(f"build_mmpbsa_system failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e)}, indent=2)
