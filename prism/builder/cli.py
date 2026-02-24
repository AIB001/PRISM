#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder CLI - Command-line interface for PRISMBuilder

This module contains the main() function and argparse setup.
Extracted from the original monolithic builder.py for maintainability.
"""

import os
import sys
import argparse

from .core import PRISMBuilder


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

  # Multiple ligands (options must come AFTER all ligand files)
  prism protein.pdb ligand1.mol2 ligand2.mol2 -o output_dir -lff openff -ff amber14sb

  # Multiple ligands (alternative: use flags for more flexibility)
  prism -pf protein.pdb -lf ligand1.mol2 -lf ligand2.mol2 -o output_dir -lff gaff -ff amber14sb

  # Using OpenFF force field with specific protein force field
  prism protein.pdb ligand.sdf -o output_dir --ligand-forcefield openff --forcefield amber14sb

  # Using OPLS-AA force field (via LigParGen server, requires internet)
  prism protein.pdb ligand.mol2 -o output_dir --ligand-forcefield opls

  # Using CGenFF force field (requires web-downloaded files)
  # Single ligand:
  prism protein.pdb dummy.mol2 -o output_dir --ligand-forcefield cgenff --forcefield-path /path/to/cgenff_files
  # Multiple ligands (specify one --forcefield-path per ligand):
  prism -pf protein.pdb -lf ligand1.mol2 -lf ligand2.mol2 -o output_dir --ligand-forcefield cgenff -ffp /path/to/cgenff1 -ffp /path/to/cgenff2

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
        """,
    )

    # Positional arguments (must come before options, or use flags instead)
    parser.add_argument("protein", nargs="?", help="Path to protein PDB file")
    parser.add_argument(
        "ligands",
        nargs="*",
        help="Path(s) to ligand file(s) (MOL2/SDF). For multiple ligands, put all files before options or use --ligand-file flag",
    )

    # Alternative way to specify input files
    parser.add_argument("--protein-file", "-pf", help="Path to protein PDB file (alternative to positional)")
    parser.add_argument(
        "--ligand-file", "-lf", action="append", help="Path to ligand file (can be specified multiple times)"
    )

    # Basic options
    basic = parser.add_argument_group("Basic options")
    basic.add_argument("--output", "-o", default="prism_output", help="Output directory (default: prism_output)")
    basic.add_argument("--config", "-c", help="Path to configuration YAML file")
    basic.add_argument("--overwrite", "-f", action="store_true", help="Overwrite existing files")

    # Force field options
    ff_group = parser.add_argument_group("Force field options")
    ff_group.add_argument(
        "--forcefield", "-ff", type=str, default="amber99sb", help="Protein force field name (default: amber99sb)"
    )
    ff_group.add_argument(
        "--water", "-w", type=str, default=None, help="Water model name (default: auto-detect from force field)"
    )
    ff_group.add_argument(
        "--ligand-forcefield",
        "-lff",
        choices=["gaff", "gaff2", "openff", "cgenff", "opls", "mmff", "match", "hybrid"],
        default="gaff",
        help="Force field for ligand: gaff (default), gaff2 (improved), openff, cgenff, opls (OPLS-AA via LigParGen), mmff/match/hybrid (SwissParam)",
    )
    ff_group.add_argument(
        "--ligand-number",
        "-n",
        type=int,
        default=None,
        help="Number of ligands (default: auto-detect from input files)",
    )
    ff_group.add_argument(
        "--ligand-charge", type=int, default=0, help="Net charge of ligand (default: 0, applies to all ligands)"
    )
    ff_group.add_argument(
        "--forcefield-path",
        "-ffp",
        type=str,
        action="append",
        default=None,
        help="Path to force field directory (required for cgenff: directory containing web-downloaded CGenFF files). For multiple ligands, specify multiple times: -ffp path1 -ffp path2",
    )
    ff_group.add_argument(
        "--nter",
        type=str,
        default=None,
        help="N-terminal type for pdb2gmx (e.g., 'NH3+' for CHARMM36). Default: auto-select. Use when force field requires specific terminus.",
    )
    ff_group.add_argument(
        "--cter",
        type=str,
        default=None,
        help="C-terminal type for pdb2gmx (e.g., 'COO-' for CHARMM36). Default: auto-select. Use when force field requires specific terminus.",
    )

    # Gaussian RESP charge options
    gaussian_group = parser.add_argument_group("Gaussian RESP charge options")
    gaussian_group.add_argument(
        "--gaussian",
        "-g",
        choices=["hf", "dft"],
        default=None,
        help="Enable Gaussian RESP charge calculation: 'hf' (HF/6-31G*) or 'dft' (B3LYP/6-31G*). Default: disabled (uses AM1-BCC)",
    )
    gaussian_group.add_argument(
        "--isopt",
        choices=["true", "false"],
        default="false",
        help="Perform geometry optimization before ESP calculation (default: false)",
    )
    gaussian_group.add_argument(
        "--respfile",
        "-rf",
        action="append",
        default=None,
        help="Path to existing RESP mol2 file. For multiple ligands, specify once per ligand: -rf resp1.mol2 -rf resp2.mol2",
    )
    gaussian_group.add_argument(
        "--nproc", type=int, default=16, help="Number of processors for Gaussian calculation (default: 16)"
    )
    gaussian_group.add_argument(
        "--mem", type=str, default="4GB", help="Memory allocation for Gaussian calculation (default: 4GB)"
    )

    # Box options
    box_group = parser.add_argument_group("Box options")
    box_group.add_argument(
        "--box-distance", "-d", type=float, default=1.5, help="Distance from protein to box edge in nm (default: 1.5)"
    )
    box_group.add_argument(
        "--box-shape",
        "-bs",
        choices=["cubic", "dodecahedron", "octahedron"],
        default="cubic",
        help="Box shape (default: cubic)",
    )
    box_group.add_argument("--no-center", action="store_true", help="Don't center protein in box")

    # Simulation parameters
    sim_group = parser.add_argument_group("Simulation parameters")
    sim_group.add_argument("--temperature", "-t", type=float, default=310, help="Temperature in K (default: 310)")
    sim_group.add_argument("--pressure", "-p", type=float, default=1.0, help="Pressure in bar (default: 1.0)")
    sim_group.add_argument("--pH", type=float, default=7.0, help="pH for protonation states (default: 7.0)")
    sim_group.add_argument(
        "--protonation",
        "-proton",
        choices=["gromacs", "propka"],
        default="gromacs",
        help="Protonation method: 'gromacs' (default) or 'propka' (pKa-based ionizable residue protonation)",
    )
    sim_group.add_argument("--production-ns", type=float, default=500, help="Production time in ns (default: 500)")
    sim_group.add_argument("--dt", type=float, default=0.002, help="Time step in ps (default: 0.002)")
    sim_group.add_argument("--nvt-ps", type=float, default=500, help="NVT equilibration time in ps (default: 500)")
    sim_group.add_argument("--npt-ps", type=float, default=500, help="NPT equilibration time in ps (default: 500)")

    # Ion options
    ion_group = parser.add_argument_group("Ion options")
    ion_group.add_argument("--no-neutralize", action="store_true", help="Don't neutralize the system")
    ion_group.add_argument(
        "--salt-concentration", "-sc", type=float, default=0.15, help="Salt concentration in M (default: 0.15)"
    )
    ion_group.add_argument("--positive-ion", "-pion", default="NA", help="Positive ion type (default: NA)")
    ion_group.add_argument("--negative-ion", "-nion", default="CL", help="Negative ion type (default: CL)")

    # Protein preparation options
    prep_group = parser.add_argument_group("Protein preparation")
    prep_group.add_argument(
        "--nterm-met",
        choices=["keep", "drop", "auto"],
        default="keep",
        help="Handle N-terminal MET residues: keep (default), drop, auto (drop for CHARMM36*)",
    )

    # Energy minimization
    em_group = parser.add_argument_group("Energy minimization")
    em_group.add_argument(
        "--em-tolerance", type=float, default=200.0, help="Energy minimization tolerance in kJ/mol/nm (default: 200.0)"
    )
    em_group.add_argument("--em-steps", type=int, default=10000, help="Maximum EM steps (default: 10000)")

    # Output options
    out_group = parser.add_argument_group("Output options")
    out_group.add_argument(
        "--traj-interval", type=float, default=500, help="Trajectory output interval in ps (default: 500)"
    )
    out_group.add_argument(
        "--energy-interval", type=float, default=10, help="Energy output interval in ps (default: 10)"
    )
    out_group.add_argument("--no-compressed", action="store_true", help="Don't use compressed trajectory format")

    # Utility options
    util_group = parser.add_argument_group("Utility options")
    util_group.add_argument("--list-forcefields", action="store_true", help="List available force fields and exit")
    util_group.add_argument(
        "--install-forcefields",
        action="store_true",
        help="Install additional force fields from PRISM configs to GROMACS",
    )
    util_group.add_argument("--export-config", metavar="FILE", help="Export default configuration to file and exit")
    util_group.add_argument(
        "--gmx-command", default=None, help="GROMACS command to use (auto-detected if not specified)"
    )

    # MM/PBSA options
    mmpbsa_group = parser.add_argument_group("MM/PBSA options")
    mmpbsa_group.add_argument(
        "--mmpbsa",
        "-pbsa",
        action="store_true",
        help="Enable MM/PBSA mode. Default: single-frame (EM->NVT->NPT->gmx_MMPBSA, "
        "no production MD). Use --mmpbsa-traj for trajectory-based mode. "
        "Output goes to GMX_PROLIG_MMPBSA/ with mmpbsa_run.sh script.",
    )
    mmpbsa_group.add_argument(
        "--mmpbsa-traj",
        type=float,
        default=None,
        metavar="NS",
        help="Use trajectory-based MM/PBSA with specified production MD length (ns). "
        "Requires --mmpbsa. Default without this flag: single-frame MM/PBSA "
        "(no production MD). Analysis interval defaults to ~1 ns.",
    )
    mmpbsa_group.add_argument(
        "--gmx2amber",
        action="store_true",
        help="Use parmed + AMBER MMPBSA.py instead of gmx_MMPBSA. "
        "Requires AmberTools with parmed. No gmx_MMPBSA needed.",
    )

    # PMF (Steered MD) options
    pmf_group = parser.add_argument_group("PMF (Steered MD) options")
    pmf_group.add_argument(
        "--pmf",
        action="store_true",
        help="Enable PMF mode: build system for steered MD and PMF calculations. "
        "Aligns pull vector to Z-axis and extends box in Z direction.",
    )
    pmf_group.add_argument(
        "--pull-vector",
        "-pullvec",
        type=int,
        nargs=2,
        metavar=("PROT_ATOM", "LIG_ATOM"),
        help="Custom pull vector defined by atom indices: protein atom -> ligand atom. "
        "Default: pocket centroid -> ligand centroid (4A cutoff)",
    )
    pmf_group.add_argument(
        "--box-extension",
        "-boxext",
        type=float,
        nargs=3,
        metavar=("X", "Y", "Z"),
        default=[0.0, 0.0, 2.0],
        help="Box extension in X, Y, Z directions (nm) for PMF mode. "
        "Default: 0.0 0.0 2.0 (extends Z by 2 nm for pulling space)",
    )
    pmf_group.add_argument(
        "--umbrella-time", type=float, default=10.0, help="Simulation time per umbrella window in ns (default: 10.0)"
    )
    pmf_group.add_argument(
        "--umbrella-spacing", type=float, default=0.12, help="Spacing between umbrella windows in nm (default: 0.12)"
    )
    pmf_group.add_argument(
        "--wham-begin", type=int, default=1000, help="Time in ps to discard for WHAM equilibration (default: 1000)"
    )
    pmf_group.add_argument(
        "--wham-bootstrap",
        type=int,
        default=200,
        help="Number of bootstrap iterations for WHAM error estimation (default: 200)",
    )

    # REST2 (Replica Exchange Solute Tempering) options
    rest2_group = parser.add_argument_group("REST2 (Replica Exchange Solute Tempering) options")
    rest2_group.add_argument(
        "--rest2",
        action="store_true",
        help="Enable REST2 mode: build MD system then generate REST2 replica exchange setup",
    )
    rest2_group.add_argument("--t-ref", type=float, default=310.0, help="Reference temperature in K (default: 310)")
    rest2_group.add_argument(
        "--t-max", type=float, default=450.0, help="Maximum effective temperature in K (default: 450)"
    )
    rest2_group.add_argument("--replica-number", type=int, default=16, help="Number of REST2 replicas (default: 16)")
    rest2_group.add_argument(
        "--rest2-cutoff", type=float, default=0.5, help="Pocket detection cutoff in nm for REST2 (default: 0.5)"
    )

    args = parser.parse_args()

    # Handle --export-config option
    if args.export_config:
        import pkg_resources

        try:
            # Try to get from package resources
            default_config = pkg_resources.resource_string("prism", "configs/default_config.yaml").decode("utf-8")
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

protein_preparation:
  nterm_met: keep  # keep, drop, or auto (drop for CHARMM36*)

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
        with open(args.export_config, "w") as f:
            f.write(default_config)
        print(f"Default configuration exported to: {args.export_config}")
        sys.exit(0)

    # Handle --install-forcefields option
    if args.install_forcefields:
        try:
            from ..utils.forcefield_installer import ForceFieldInstaller

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
            from ..utils.environment import GromacsEnvironment

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

    # Collect ligand paths from both positional and flag arguments
    ligand_paths = []
    if args.ligands:  # Positional ligands
        ligand_paths.extend(args.ligands)
    if args.ligand_file:  # Flag-based ligands (can be specified multiple times)
        ligand_paths.extend(args.ligand_file)

    if not protein_path:
        parser.error("Protein file is required. Use positional argument or --protein-file flag")

    if not ligand_paths:
        parser.error(
            "At least one ligand file is required. Provide ligand path(s) as positional arguments or use --ligand-file flag"
        )

    # Validate ligand number
    actual_ligand_count = len(ligand_paths)
    if args.ligand_number is not None:
        if args.ligand_number != actual_ligand_count:
            parser.error(
                f"--ligand-number (-n) specified {args.ligand_number} but {actual_ligand_count} ligand file(s) provided"
            )

    # Validate --mmpbsa-traj requires --mmpbsa
    if args.mmpbsa_traj is not None and not args.mmpbsa:
        parser.error("--mmpbsa-traj requires --mmpbsa to be enabled")

    # Validate --gmx2amber requires --mmpbsa
    if args.gmx2amber and not args.mmpbsa:
        parser.error("--gmx2amber requires --mmpbsa to be enabled")

    # Update args to use resolved paths
    args.protein = protein_path
    args.ligands = ligand_paths
    args.ligand_number = actual_ligand_count

    # Build kwargs from command-line arguments
    kwargs = {}

    # Process command-line overrides for config
    config_overrides = {}

    # General overrides
    if args.gmx_command:
        config_overrides.setdefault("general", {})["gmx_command"] = args.gmx_command

    # Box overrides
    if args.box_distance != 1.5:
        config_overrides.setdefault("box", {})["distance"] = args.box_distance
    if args.box_shape != "cubic":
        config_overrides.setdefault("box", {})["shape"] = args.box_shape
    if args.no_center:
        config_overrides.setdefault("box", {})["center"] = False

    # Simulation overrides
    if args.temperature != 310:
        config_overrides.setdefault("simulation", {})["temperature"] = args.temperature
    if args.pressure != 1.0:
        config_overrides.setdefault("simulation", {})["pressure"] = args.pressure
    if args.pH != 7.0:
        config_overrides.setdefault("simulation", {})["pH"] = args.pH
    if args.protonation != "gromacs":
        config_overrides.setdefault("protonation", {})["method"] = args.protonation
    if args.production_ns != 500:
        config_overrides.setdefault("simulation", {})["production_time_ns"] = args.production_ns
    if args.dt != 0.002:
        config_overrides.setdefault("simulation", {})["dt"] = args.dt
    if args.nvt_ps != 500:
        config_overrides.setdefault("simulation", {})["equilibration_nvt_time_ps"] = args.nvt_ps
    if args.npt_ps != 500:
        config_overrides.setdefault("simulation", {})["equilibration_npt_time_ps"] = args.npt_ps
    if args.ligand_charge != 0:
        config_overrides.setdefault("simulation", {})["ligand_charge"] = args.ligand_charge

    # Ion overrides
    if args.no_neutralize:
        config_overrides.setdefault("ions", {})["neutral"] = False
    if args.salt_concentration != 0.15:
        config_overrides.setdefault("ions", {})["concentration"] = args.salt_concentration
    if args.positive_ion != "NA":
        config_overrides.setdefault("ions", {})["positive_ion"] = args.positive_ion
    if args.negative_ion != "CL":
        config_overrides.setdefault("ions", {})["negative_ion"] = args.negative_ion
    if args.nterm_met != "keep":
        config_overrides.setdefault("protein_preparation", {})["nterm_met"] = args.nterm_met

    # Energy minimization overrides
    if args.em_tolerance != 200.0:
        config_overrides.setdefault("energy_minimization", {})["emtol"] = args.em_tolerance
    if args.em_steps != 10000:
        config_overrides.setdefault("energy_minimization", {})["nsteps"] = args.em_steps

    # Output overrides
    if args.traj_interval != 500:
        config_overrides.setdefault("output", {})["trajectory_interval_ps"] = args.traj_interval
    if args.energy_interval != 10:
        config_overrides.setdefault("output", {})["energy_interval_ps"] = args.energy_interval
    if args.no_compressed:
        config_overrides.setdefault("output", {})["compressed_trajectory"] = False

    # PMF/Umbrella overrides
    if args.pmf:
        pmf_overrides = {}
        if args.umbrella_time != 10.0:
            pmf_overrides["umbrella_time_ns"] = args.umbrella_time
        if args.umbrella_spacing != 0.12:
            pmf_overrides["umbrella_spacing"] = args.umbrella_spacing
        if args.wham_begin != 1000:
            pmf_overrides["wham_begin"] = args.wham_begin
        if args.wham_bootstrap != 200:
            pmf_overrides["wham_bootstrap"] = args.wham_bootstrap
        if pmf_overrides:
            config_overrides.setdefault("pmf", {}).update(pmf_overrides)

    # Save config overrides if any
    if config_overrides:
        # Create temporary config file with overrides
        import tempfile
        import yaml

        # If user provided a config, load it first
        if args.config:
            with open(args.config, "r") as f:
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
        temp_config = tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False)
        yaml.dump(base_config, temp_config, default_flow_style=False)
        temp_config.close()
        kwargs["config_path"] = temp_config.name
    elif args.config:
        kwargs["config_path"] = args.config

    # Create and run PRISM Builder
    builder = PRISMBuilder(
        args.protein,
        args.ligands,  # Now a list of ligand paths
        args.output,
        ligand_forcefield=args.ligand_forcefield,
        forcefield=args.forcefield,
        water_model=args.water,
        overwrite=args.overwrite,
        forcefield_path=args.forcefield_path,
        nter=args.nter,  # N-terminal type
        cter=args.cter,  # C-terminal type
        gaussian_method=args.gaussian,  # Gaussian RESP method (hf/dft)
        do_optimization=(args.isopt == "true"),  # Geometry optimization
        resp_files=args.respfile,  # Existing RESP mol2 files
        gaussian_nproc=args.nproc,  # Gaussian processors
        gaussian_mem=args.mem,  # Gaussian memory
        pmf_mode=args.pmf,  # PMF mode flag
        pullvec=tuple(args.pull_vector) if args.pull_vector else None,  # Pull vector atom indices
        box_extension=tuple(args.box_extension) if args.box_extension else None,  # Box extension (X,Y,Z)
        rest2_mode=args.rest2,  # REST2 mode flag
        t_ref=args.t_ref,  # REST2 reference temperature
        t_max=args.t_max,  # REST2 max effective temperature
        n_replicas=args.replica_number,  # REST2 number of replicas
        rest2_cutoff=args.rest2_cutoff,  # REST2 pocket detection cutoff
        mmpbsa_mode=args.mmpbsa,  # MMPBSA mode flag
        mmpbsa_traj_ns=args.mmpbsa_traj if args.mmpbsa else None,  # Trajectory-based MMPBSA length (ns)
        gmx2amber=args.gmx2amber if args.mmpbsa else False,  # AMBER MMPBSA.py backend
        **kwargs,
    )

    try:
        builder.run()
    finally:
        # Clean up temporary config if created
        if config_overrides and "config_path" in kwargs:
            try:
                os.unlink(kwargs["config_path"])
            except:
                pass
