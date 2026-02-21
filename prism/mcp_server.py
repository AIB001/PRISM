#!/usr/bin/env python3
"""
PRISM MCP Server - AI Agent interface for protein-ligand system building.

This MCP server exposes PRISM's molecular dynamics system building capabilities
as tools that AI agents (Claude, ChatGPT, etc.) can call through natural
language conversation.

Usage:
    # Test with MCP Inspector (browser-based debug UI)
    mcp dev prism/mcp_server.py

    # Connect to Claude Code
    claude mcp add --transport stdio prism -- python prism/mcp_server.py
"""

import os
import sys
import json
import logging
import traceback
from typing import Optional

from mcp.server.fastmcp import FastMCP

# ==========================================================================
#  CRITICAL: Redirect all logging to stderr.
#
#  In stdio transport mode, stdout is exclusively used for MCP JSON-RPC
#  messages. ANY non-MCP output on stdout will corrupt the protocol.
#  PRISM internally uses print() extensively for progress output, so we
#  must redirect stdout to stderr before calling any PRISM code.
# ==========================================================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("prism-mcp")


# ==========================================================================
#  Initialize MCP Server
# ==========================================================================
mcp = FastMCP("prism")


# ==========================================================================
#  Helper: stdout suppressor
#
#  PRISM code uses print() everywhere for colored progress output.
#  When running as an MCP stdio server, those prints would go to stdout
#  and corrupt the JSON-RPC communication channel.
#
#  This context manager temporarily redirects sys.stdout -> sys.stderr
#  so all PRISM prints go to stderr (visible in logs, invisible to MCP).
# ==========================================================================
class _StdoutToStderr:
    """Context manager that redirects stdout to stderr."""

    def __enter__(self):
        self._original = sys.stdout
        sys.stdout = sys.stderr
        return self

    def __exit__(self, *args):
        sys.stdout = self._original


# ==========================================================================
#  Tool 1: Check Dependencies
# ==========================================================================
@mcp.tool()
def check_dependencies() -> str:
    """Check if all required dependencies for PRISM are installed and available.

    Checks for the following software/libraries:
    - GROMACS (molecular dynamics engine)
    - pdbfixer (protein structure preparation)
    - antechamber (AmberTools, for GAFF/GAFF2 force fields)
    - OpenFF toolkit (for OpenFF force field)
    - MDTraj (for trajectory analysis)
    - RDKit (for chemical informatics)

    Returns a JSON object with each dependency name mapped to true/false.
    """
    logger.info("Checking PRISM dependencies...")
    try:
        with _StdoutToStderr():
            # Add PRISM to path if needed
            prism_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if prism_root not in sys.path:
                sys.path.insert(0, prism_root)

            import prism as pm
            deps = pm.check_dependencies()

        logger.info(f"Dependencies: {deps}")
        return json.dumps(deps, indent=2)

    except Exception as e:
        logger.error(f"check_dependencies failed: {e}")
        return json.dumps({"error": str(e)})


# ==========================================================================
#  Tool 2: List Available Force Fields
# ==========================================================================
@mcp.tool()
def list_forcefields() -> str:
    """List all available force fields for protein-ligand system building.

    Returns two categories:
    1. Protein force fields: detected from the user's GROMACS installation
       (e.g., amber99sb, amber14sb, charmm27, etc.)
    2. Ligand force fields: supported by PRISM for small molecule parameterization
       (gaff, gaff2, openff, cgenff, opls, mmff, match, hybrid)

    Use this tool to help the user choose appropriate force fields before building.
    """
    logger.info("Listing available force fields...")
    try:
        with _StdoutToStderr():
            prism_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if prism_root not in sys.path:
                sys.path.insert(0, prism_root)

            from prism.forcefield import get_generator_info

            ligand_info = get_generator_info()

            # Try to get protein force fields from GROMACS
            protein_ffs = []
            try:
                from prism.utils.environment import GromacsEnvironment
                env = GromacsEnvironment()
                protein_ffs = env.list_force_fields()
            except Exception:
                protein_ffs = ["(GROMACS not detected - cannot list protein force fields)"]

        result = {
            "protein_forcefields": protein_ffs,
            "ligand_forcefields": {
                name: info["description"] for name, info in ligand_info.items()
            },
            "ligand_forcefield_cli_options": [
                "gaff", "gaff2", "openff", "cgenff",
                "opls", "mmff", "match", "hybrid",
            ],
        }
        return json.dumps(result, indent=2, ensure_ascii=False)

    except Exception as e:
        logger.error(f"list_forcefields failed: {e}")
        return json.dumps({"error": str(e)})


# ==========================================================================
#  Tool 3: Build Standard MD System (Normal Mode)
# ==========================================================================
@mcp.tool()
def build_system(
    protein_path: str,
    ligand_path: str,
    output_dir: str = "prism_output",
    ligand_forcefield: str = "gaff",
    forcefield: str = "amber99sb",
    water_model: str = "tip3p",
    temperature: float = 310.0,
    ph: float = 7.0,
    production_ns: float = 500.0,
    box_distance: float = 1.5,
    box_shape: str = "cubic",
    salt_concentration: float = 0.15,
    ligand_charge: int = 0,
    overwrite: bool = False,
) -> str:
    """Build a complete protein-ligand system for GROMACS molecular dynamics simulation.

    This is the main building tool. It takes a protein PDB file and a ligand file,
    then performs the full workflow:

    Step 1: Generate ligand force field parameters (GAFF/GAFF2/OpenFF/etc.)
    Step 2: Clean and prepare the protein structure (remove artifacts, handle metals)
    Step 3: Build GROMACS system:
            - pdb2gmx (protein topology)
            - Combine protein + ligand
            - Create simulation box
            - Solvate with water
            - Add ions to neutralize and set salt concentration
    Step 4: Generate MDP parameter files (em.mdp, nvt.mdp, npt.mdp, md.mdp)
    Step 5: Generate localrun.sh script for GPU-accelerated simulation

    After building, the user runs: cd GMX_PROLIG_MD && bash localrun.sh

    Args:
        protein_path: Absolute path to the protein PDB file.
        ligand_path: Absolute path to the ligand file (MOL2 or SDF format).
        output_dir: Directory for all output files. Default: "prism_output".
        ligand_forcefield: Force field for ligand parameterization.
            Options: "gaff" (default), "gaff2", "openff", "cgenff", "opls",
            "mmff", "match", "hybrid".
        forcefield: Protein force field. Default: "amber99sb".
            Common choices: "amber99sb", "amber14sb", "amber99sb-ildn", "charmm27".
        water_model: Water model. Default: "tip3p".
            Options: "tip3p", "tip4p", "spc", "spce".
        temperature: Simulation temperature in Kelvin. Default: 310 (body temperature).
        ph: pH for protonation state assignment. Default: 7.0.
        production_ns: Production MD length in nanoseconds. Default: 500.
        box_distance: Distance from protein to box edge in nm. Default: 1.5.
        box_shape: Box shape. Options: "cubic", "dodecahedron", "octahedron". Default: "cubic".
        salt_concentration: NaCl concentration in mol/L. Default: 0.15.
        ligand_charge: Net formal charge of the ligand. Default: 0.
        overwrite: Overwrite existing output files if they exist. Default: false.
    """
    logger.info(f"build_system: {protein_path} + {ligand_path} -> {output_dir}")

    # --- Input validation ---
    errors = []
    if not os.path.isabs(protein_path):
        errors.append(f"protein_path must be an absolute path, got: {protein_path}")
    elif not os.path.exists(protein_path):
        errors.append(f"Protein file not found: {protein_path}")

    if not os.path.isabs(ligand_path):
        errors.append(f"ligand_path must be an absolute path, got: {ligand_path}")
    elif not os.path.exists(ligand_path):
        errors.append(f"Ligand file not found: {ligand_path}")

    if errors:
        return json.dumps({"success": False, "errors": errors})

    # --- Run PRISM builder ---
    try:
        with _StdoutToStderr():
            prism_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if prism_root not in sys.path:
                sys.path.insert(0, prism_root)

            from prism.builder import PRISMBuilder

            builder = PRISMBuilder(
                protein_path=protein_path,
                ligand_paths=ligand_path,
                output_dir=output_dir,
                ligand_forcefield=ligand_forcefield,
                forcefield=forcefield,
                water_model=water_model,
                overwrite=overwrite,
            )

            # Update simulation parameters from user input
            builder.config["simulation"]["temperature"] = temperature
            builder.config["simulation"]["pH"] = ph
            builder.config["simulation"]["production_time_ns"] = production_ns
            builder.config["box"]["distance"] = box_distance
            builder.config["box"]["shape"] = box_shape
            builder.config["ions"]["concentration"] = salt_concentration
            builder.config["ligand_forcefield"]["charge"] = ligand_charge

            result_dir = builder.run()

        # --- Collect output file paths ---
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
                    "temperature_K": temperature,
                    "pH": ph,
                    "production_ns": production_ns,
                    "box_distance_nm": box_distance,
                    "salt_concentration_M": salt_concentration,
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


# ==========================================================================
#  Tool 4: Build PMF System (Steered MD + Umbrella Sampling)
# ==========================================================================
@mcp.tool()
def build_pmf_system(
    protein_path: str,
    ligand_path: str,
    output_dir: str = "prism_pmf_output",
    ligand_forcefield: str = "gaff",
    forcefield: str = "amber99sb",
    box_extension_z: float = 2.0,
    umbrella_time_ns: float = 10.0,
    overwrite: bool = False,
) -> str:
    """Build a system for PMF (Potential of Mean Force) calculation.

    Sets up steered molecular dynamics (SMD) and umbrella sampling for
    calculating the binding free energy profile of a protein-ligand complex.

    The workflow:
    1. Generate ligand force field
    2. Clean protein
    3. Align the complex so the unbinding pull vector is along the Z-axis
    4. Build GROMACS system with an extended box in Z for pulling space
    5. Generate index file with pull/reference/freeze atom groups
    6. Generate SMD and umbrella sampling MDP files and run scripts

    After building, the user runs two scripts sequentially:
      cd GMX_PROLIG_PMF && bash smd_run.sh && bash umbrella_run.sh

    Args:
        protein_path: Absolute path to protein PDB file.
        ligand_path: Absolute path to ligand file (MOL2/SDF).
        output_dir: Output directory. Default: "prism_pmf_output".
        ligand_forcefield: Ligand force field. Default: "gaff".
        forcefield: Protein force field. Default: "amber99sb".
        box_extension_z: Extra box length in Z direction (nm) for pulling space. Default: 2.0.
        umbrella_time_ns: Simulation time per umbrella window in ns. Default: 10.0.
        overwrite: Overwrite existing files. Default: false.
    """
    logger.info(f"build_pmf_system: {protein_path} + {ligand_path}")

    if not os.path.exists(protein_path):
        return json.dumps({"success": False, "error": f"Protein not found: {protein_path}"})
    if not os.path.exists(ligand_path):
        return json.dumps({"success": False, "error": f"Ligand not found: {ligand_path}"})

    try:
        with _StdoutToStderr():
            prism_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if prism_root not in sys.path:
                sys.path.insert(0, prism_root)

            from prism.builder import PRISMBuilder

            builder = PRISMBuilder(
                protein_path=protein_path,
                ligand_paths=ligand_path,
                output_dir=output_dir,
                ligand_forcefield=ligand_forcefield,
                forcefield=forcefield,
                overwrite=overwrite,
                pmf_mode=True,
                box_extension=(0.0, 0.0, box_extension_z),
            )
            builder.config.setdefault("pmf", {})
            builder.config["pmf"]["umbrella_time_ns"] = umbrella_time_ns

            result_dir = builder.run()

        pmf_dir = os.path.join(result_dir, "GMX_PROLIG_PMF")
        return json.dumps({
            "success": True,
            "output_dir": result_dir,
            "pmf_dir": pmf_dir,
            "message": (
                f"PMF system built successfully!\n"
                f"To run the PMF workflow:\n"
                f"  cd {pmf_dir}\n"
                f"  bash smd_run.sh      # Step 1: Steered MD\n"
                f"  bash umbrella_run.sh  # Step 2: Umbrella sampling + WHAM"
            ),
        }, indent=2)

    except Exception as e:
        logger.error(f"build_pmf_system failed: {e}\n{traceback.format_exc()}")
        return json.dumps({"success": False, "error": str(e)}, indent=2)


# ==========================================================================
#  Tool 5: Build REST2 System
# ==========================================================================
@mcp.tool()
def build_rest2_system(
    protein_path: str,
    ligand_path: str,
    output_dir: str = "prism_rest2_output",
    ligand_forcefield: str = "gaff",
    forcefield: str = "amber99sb",
    t_ref: float = 310.0,
    t_max: float = 450.0,
    n_replicas: int = 16,
    rest2_cutoff: float = 0.5,
    overwrite: bool = False,
) -> str:
    """Build a system for REST2 (Replica Exchange with Solute Tempering v2) enhanced sampling.

    REST2 enhances conformational sampling of the protein-ligand binding pocket
    by running multiple replicas at different effective temperatures for the
    solute region, while keeping the solvent at the reference temperature.

    The workflow:
    1. Build a standard MD system (same as build_system)
    2. Convert GMX_PROLIG_MD to GMX_PROLIG_REST2 with scaled topologies

    After building, run: cd GMX_PROLIG_REST2 && bash rest2_run.sh

    Args:
        protein_path: Absolute path to protein PDB file.
        ligand_path: Absolute path to ligand file (MOL2/SDF).
        output_dir: Output directory. Default: "prism_rest2_output".
        ligand_forcefield: Ligand force field. Default: "gaff".
        forcefield: Protein force field. Default: "amber99sb".
        t_ref: Reference (physical) temperature in Kelvin. Default: 310.
        t_max: Maximum effective temperature in Kelvin. Default: 450.
        n_replicas: Number of REST2 replicas. Default: 16.
        rest2_cutoff: Distance cutoff (nm) for identifying pocket residues. Default: 0.5.
        overwrite: Overwrite existing files. Default: false.
    """
    logger.info(f"build_rest2_system: {protein_path} + {ligand_path}")

    if not os.path.exists(protein_path):
        return json.dumps({"success": False, "error": f"Protein not found: {protein_path}"})
    if not os.path.exists(ligand_path):
        return json.dumps({"success": False, "error": f"Ligand not found: {ligand_path}"})

    try:
        with _StdoutToStderr():
            prism_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if prism_root not in sys.path:
                sys.path.insert(0, prism_root)

            from prism.builder import PRISMBuilder

            builder = PRISMBuilder(
                protein_path=protein_path,
                ligand_paths=ligand_path,
                output_dir=output_dir,
                ligand_forcefield=ligand_forcefield,
                forcefield=forcefield,
                overwrite=overwrite,
                rest2_mode=True,
                t_ref=t_ref,
                t_max=t_max,
                n_replicas=n_replicas,
                rest2_cutoff=rest2_cutoff,
            )
            result_dir = builder.run()

        rest2_dir = os.path.join(result_dir, "GMX_PROLIG_REST2")
        return json.dumps({
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
        }, indent=2)

    except Exception as e:
        logger.error(f"build_rest2_system failed: {e}\n{traceback.format_exc()}")
        return json.dumps({"success": False, "error": str(e)}, indent=2)


# ==========================================================================
#  Tool 6: Build MM/PBSA System
# ==========================================================================
@mcp.tool()
def build_mmpbsa_system(
    protein_path: str,
    ligand_path: str,
    output_dir: str = "prism_mmpbsa_output",
    ligand_forcefield: str = "gaff",
    forcefield: str = "amber99sb",
    mmpbsa_traj_ns: Optional[float] = None,
    gmx2amber: bool = False,
    overwrite: bool = False,
) -> str:
    """Build a system for MM/PBSA binding free energy calculation.

    MM/PBSA estimates protein-ligand binding affinity from MD snapshots.
    Two sub-modes are available:

    1. Single-frame mode (default): EM -> NVT -> NPT -> MM/PBSA on the
       equilibrated structure. Fast but less accurate.
    2. Trajectory mode (set mmpbsa_traj_ns): EM -> NVT -> NPT -> Production MD
       -> MM/PBSA on multiple trajectory frames. More accurate.

    After building, run: cd GMX_PROLIG_MMPBSA && bash mmpbsa_run.sh

    Args:
        protein_path: Absolute path to protein PDB file.
        ligand_path: Absolute path to ligand file (MOL2/SDF).
        output_dir: Output directory. Default: "prism_mmpbsa_output".
        ligand_forcefield: Ligand force field. Default: "gaff".
        forcefield: Protein force field. Default: "amber99sb".
        mmpbsa_traj_ns: Production MD length in ns for trajectory-based MM/PBSA.
            If not provided, uses single-frame mode (no production MD).
        gmx2amber: Use AMBER MMPBSA.py via parmed conversion instead of gmx_MMPBSA.
            Requires AmberTools. Default: false.
        overwrite: Overwrite existing files. Default: false.
    """
    logger.info(f"build_mmpbsa_system: {protein_path} + {ligand_path}")

    if not os.path.exists(protein_path):
        return json.dumps({"success": False, "error": f"Protein not found: {protein_path}"})
    if not os.path.exists(ligand_path):
        return json.dumps({"success": False, "error": f"Ligand not found: {ligand_path}"})

    try:
        with _StdoutToStderr():
            prism_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if prism_root not in sys.path:
                sys.path.insert(0, prism_root)

            from prism.builder import PRISMBuilder

            builder = PRISMBuilder(
                protein_path=protein_path,
                ligand_paths=ligand_path,
                output_dir=output_dir,
                ligand_forcefield=ligand_forcefield,
                forcefield=forcefield,
                overwrite=overwrite,
                mmpbsa_mode=True,
                mmpbsa_traj_ns=mmpbsa_traj_ns,
                gmx2amber=gmx2amber,
            )
            result_dir = builder.run()

        mode = "trajectory" if mmpbsa_traj_ns else "single-frame"
        mmpbsa_dir = os.path.join(result_dir, "GMX_PROLIG_MMPBSA")
        return json.dumps({
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
        }, indent=2)

    except Exception as e:
        logger.error(f"build_mmpbsa_system failed: {e}\n{traceback.format_exc()}")
        return json.dumps({"success": False, "error": str(e)}, indent=2)


# ==========================================================================
#  Resource: default configuration template
# ==========================================================================
@mcp.resource("prism://config/default")
def get_default_config() -> str:
    """Return the default PRISM configuration YAML as a reference template."""
    config_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "configs", "default_config.yaml"
    )
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            return f.read()
    return "# default_config.yaml not found"


# ==========================================================================
#  Entry point
# ==========================================================================
if __name__ == "__main__":
    logger.info("Starting PRISM MCP Server...")
    mcp.run(transport="stdio")
