#!/usr/bin/env python3
"""
PRISM MCP Server - AI Agent interface for protein-ligand system building.

This MCP server exposes PRISM's molecular dynamics system building capabilities
as tools that AI agents (Claude, etc.) can call through natural language.

Usage:
    # Test with MCP Inspector (browser-based debug UI)
    mcp dev prism/mcp_server.py

    # Connect to Claude Code
    claude mcp add --transport stdio prism -- python prism/mcp_server.py
"""

import os
import re
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


def _ensure_prism_importable():
    """Add PRISM root to sys.path if not already present."""
    prism_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if prism_root not in sys.path:
        sys.path.insert(0, prism_root)


# ==========================================================================
#  Tool 1: Check Dependencies
# ==========================================================================
@mcp.tool()
def check_dependencies() -> str:
    """Check if all required dependencies for PRISM are installed.

    Checks:
    - GROMACS (molecular dynamics engine, required)
    - pdbfixer (protein structure preparation, recommended)
    - AmberTools / antechamber (for GAFF/GAFF2 ligand force fields)
    - OpenFF toolkit (for OpenFF ligand force field)
    - PROPKA (for pKa-based histidine protonation prediction)
    - Gaussian g16 (for high-precision RESP charge calculation, optional)
    - MDTraj (for trajectory analysis)
    - RDKit (for chemical informatics)

    Call this tool first to verify the user's environment before building.

    Returns a JSON object with each dependency name mapped to true/false.
    """
    logger.info("Checking PRISM dependencies...")
    try:
        with _StdoutToStderr():
            _ensure_prism_importable()
            import subprocess, shutil

            deps = {
                "gromacs": False,
                "pdbfixer": False,
                "antechamber": False,
                "openff": False,
                "propka": False,
                "gaussian_g16": False,
                "mdtraj": False,
                "rdkit": False,
            }

            # GROMACS
            try:
                from prism.utils.environment import GromacsEnvironment
                GromacsEnvironment()
                deps["gromacs"] = True
            except Exception:
                pass

            # pdbfixer
            deps["pdbfixer"] = shutil.which("pdbfixer") is not None

            # AmberTools (antechamber)
            deps["antechamber"] = shutil.which("antechamber") is not None

            # OpenFF
            try:
                import openff.toolkit
                deps["openff"] = True
            except ImportError:
                pass

            # PROPKA
            try:
                import propka
                deps["propka"] = True
            except ImportError:
                pass

            # Gaussian g16
            deps["gaussian_g16"] = shutil.which("g16") is not None

            # MDTraj
            try:
                import mdtraj
                deps["mdtraj"] = True
            except ImportError:
                pass

            # RDKit
            try:
                import rdkit
                deps["rdkit"] = True
            except ImportError:
                pass

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
    1. Protein force fields: detected from the user's GROMACS installation.
       Recommended: amber14sb (well-validated, modern AMBER force field).
    2. Ligand force fields: supported by PRISM for small molecule parameterization.
       Recommended: gaff2 (improved GAFF with better torsion parameters).

    Common protein + ligand combinations:
    - amber14sb + gaff2  (recommended default)
    - amber99sb + gaff   (classic, widely used)
    - amber19sb + gaff2  (latest AMBER, requires OPC water: --water opc)
    - charmm36  + cgenff (CHARMM ecosystem)

    Use this tool to help the user choose appropriate force fields.
    """
    logger.info("Listing available force fields...")
    try:
        with _StdoutToStderr():
            _ensure_prism_importable()

            # Protein force fields from GROMACS
            protein_ffs = []
            try:
                from prism.utils.environment import GromacsEnvironment
                env = GromacsEnvironment()
                protein_ffs = env.list_force_fields()
            except Exception:
                protein_ffs = ["(GROMACS not detected - cannot list protein force fields)"]

            # Ligand force fields - hardcoded from PRISM's actual support
            ligand_ffs = {
                "gaff": "GAFF force field (AmberTools). Classic, widely used.",
                "gaff2": "GAFF2 force field (AmberTools). Improved torsion parameters. RECOMMENDED.",
                "openff": "Open Force Field (openff-toolkit). Modern, data-driven.",
                "cgenff": "CGenFF (CHARMM General FF). Requires web-downloaded files from cgenff.com. CLI only (use --forcefield-path).",
                "opls": "OPLS-AA (via LigParGen server). Requires internet.",
                "mmff": "MMFF-based (via SwissParam server). Requires internet.",
                "match": "MATCH (via SwissParam server). Requires internet.",
                "hybrid": "Hybrid MMFF-MATCH (via SwissParam server). Requires internet.",
            }

        result = {
            "protein_forcefields": protein_ffs,
            "ligand_forcefields": ligand_ffs,
            "recommended": {
                "protein": "amber14sb",
                "ligand": "gaff2",
                "water": "tip3p",
            },
        }
        return json.dumps(result, indent=2, ensure_ascii=False)

    except Exception as e:
        logger.error(f"list_forcefields failed: {e}")
        return json.dumps({"error": str(e)})


# ==========================================================================
#  Tool 3: Build Standard MD System
# ==========================================================================
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
            histidine states (HID/HIE/HIP). Requires propka package.
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
    # Note: cgenff requires --forcefield-path which MCP doesn't expose; use CLI instead

    if gaussian_method and gaussian_method not in ("hf", "dft"):
        errors.append(f"gaussian_method must be 'hf' or 'dft', got: {gaussian_method}")

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

        return json.dumps({
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
            },
        }, indent=2)

    except Exception as e:
        logger.error(f"build_system failed: {e}\n{traceback.format_exc()}")
        return json.dumps({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc(),
        }, indent=2)


# ==========================================================================
#  Tool 4: Build PMF System (Steered MD + Umbrella Sampling)
# ==========================================================================
@mcp.tool()
def build_pmf_system(
    protein_path: str,
    ligand_path: str,
    output_dir: str = "prism_pmf_output",
    ligand_forcefield: str = "gaff2",
    forcefield: str = "amber14sb",
    water_model: str = "tip3p",
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
        box_extension_z: Extra Z-axis length (nm) for pulling space. Default: 2.0.
        umbrella_time_ns: Simulation time per umbrella window in ns. Default: 10.0.
        umbrella_spacing: Distance between umbrella windows in nm. Default: 0.12.
        overwrite: Overwrite existing files. Default: false.
    """
    logger.info(f"build_pmf_system: {protein_path} + {ligand_path}")

    errors = []
    if not os.path.exists(protein_path):
        errors.append(f"Protein file not found: {protein_path}")
    if not os.path.exists(ligand_path):
        errors.append(f"Ligand file not found: {ligand_path}")
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
                pmf_mode=True,
                box_extension=(0.0, 0.0, box_extension_z),
            )

            # Apply PMF-specific parameters
            builder.config.setdefault("pmf", {})
            builder.config["pmf"]["umbrella_time_ns"] = umbrella_time_ns
            builder.config["pmf"]["umbrella_spacing"] = umbrella_spacing

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
                f"  bash smd_run.sh        # Step 1: Steered MD\n"
                f"  bash umbrella_run.sh   # Step 2: Umbrella sampling + WHAM"
            ),
            "parameters": {
                "protein_forcefield": forcefield,
                "ligand_forcefield": ligand_forcefield,
                "box_extension_z_nm": box_extension_z,
                "umbrella_time_ns": umbrella_time_ns,
                "umbrella_spacing_nm": umbrella_spacing,
            },
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
    ligand_forcefield: str = "gaff2",
    forcefield: str = "amber14sb",
    water_model: str = "tip3p",
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
    if not os.path.exists(protein_path):
        errors.append(f"Protein file not found: {protein_path}")
    if not os.path.exists(ligand_path):
        errors.append(f"Ligand file not found: {ligand_path}")
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
            "parameters": {
                "protein_forcefield": forcefield,
                "ligand_forcefield": ligand_forcefield,
                "t_ref_K": t_ref,
                "t_max_K": t_max,
                "n_replicas": n_replicas,
                "rest2_cutoff_nm": rest2_cutoff,
            },
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
    ligand_forcefield: str = "gaff2",
    forcefield: str = "amber14sb",
    water_model: str = "tip3p",
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
        mmpbsa_traj_ns: Production MD length in ns for trajectory-based mode.
            If None (default), uses single-frame mode (no production MD).
        gmx2amber: Use AMBER MMPBSA.py via parmed instead of gmx_MMPBSA.
            Requires AmberTools with parmed. Default: false.
        overwrite: Overwrite existing files. Default: false.
    """
    logger.info(f"build_mmpbsa_system: {protein_path} + {ligand_path}")

    errors = []
    if not os.path.exists(protein_path):
        errors.append(f"Protein file not found: {protein_path}")
    if not os.path.exists(ligand_path):
        errors.append(f"Ligand file not found: {ligand_path}")
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
            "parameters": {
                "protein_forcefield": forcefield,
                "ligand_forcefield": ligand_forcefield,
                "mode": mode,
                "mmpbsa_traj_ns": mmpbsa_traj_ns,
                "gmx2amber": gmx2amber,
            },
        }, indent=2)

    except Exception as e:
        logger.error(f"build_mmpbsa_system failed: {e}\n{traceback.format_exc()}")
        return json.dumps({"success": False, "error": str(e)}, indent=2)


# ==========================================================================
#  Resource: default configuration template
# ==========================================================================
@mcp.resource("prism://config/default")
def get_default_config() -> str:
    """Return the default PRISM configuration YAML as a reference template.

    This shows all configurable parameters with their default values.
    Users can export this to a file with: prism --export-config my_config.yaml
    """
    config_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "configs", "default_config.yaml"
    )
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            return f.read()
    return "# default_config.yaml not found"


# ==========================================================================
#  Tool 7: Validate Input Files
# ==========================================================================
@mcp.tool()
def validate_input_files(protein_path: str, ligand_path: str) -> str:
    """Pre-build validation of protein PDB and ligand files.

    Checks both files for common issues before building, saving time by
    catching problems early. Validates file existence, format, and content.

    For proteins: counts atoms, residues, chains, histidines, checks END record.
    For ligands: detects format (MOL2/SDF/PDB), counts atoms, checks structure.

    Call this BEFORE build_system() to verify inputs are correct.

    Args:
        protein_path: Absolute path to the protein PDB file.
        ligand_path: Absolute path to the ligand file (MOL2, SDF, or PDB).
    """
    logger.info(f"validate_input_files: {protein_path}, {ligand_path}")

    result = {"protein": {}, "ligand": {}}

    # --- Protein PDB validation ---
    prot = result["protein"]
    prot["path"] = protein_path
    prot["valid"] = False
    prot["warnings"] = []

    if not os.path.exists(protein_path):
        prot["warnings"].append(f"File not found: {protein_path}")
        result["ligand"] = _validate_ligand(ligand_path)
        return json.dumps(result, indent=2)

    if os.path.getsize(protein_path) == 0:
        prot["warnings"].append("File is empty")
        result["ligand"] = _validate_ligand(ligand_path)
        return json.dumps(result, indent=2)

    try:
        atom_count = 0
        chains = set()
        residue_ids = set()
        residue_names = set()
        his_count = 0
        has_end = False

        # Standard amino acids
        standard_aa = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL",
            # Common protonation variants
            "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
            "CYX", "ASH", "GLH", "LYN",
        }

        with open(protein_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    atom_count += 1
                    chain_id = line[21] if len(line) > 21 else " "
                    chains.add(chain_id.strip() or " ")
                    resname = line[17:20].strip()
                    resseq = line[22:26].strip()
                    residue_ids.add((chain_id, resseq, resname))
                    residue_names.add(resname)
                    if resname in ("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"):
                        his_count += 1
                elif line.startswith(("END", "ENDMDL")):
                    has_end = True

        if atom_count == 0:
            prot["warnings"].append("No ATOM/HETATM records found in PDB file")
        else:
            prot["valid"] = True

        # Count unique HIS residues (not atoms)
        his_residues = {(c, s) for c, s, r in residue_ids
                        if r in ("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP")}

        # Flag non-standard residues (exclude water and common ions)
        ignore_res = {"HOH", "WAT", "SOL", "NA", "CL", "K", "MG", "CA", "ZN"}
        non_standard = residue_names - standard_aa - ignore_res
        if non_standard:
            prot["warnings"].append(
                f"Non-standard residues found: {sorted(non_standard)}. "
                f"These may be ligands, cofactors, or modified residues."
            )

        if not has_end:
            prot["warnings"].append("No END/ENDMDL record found (may be truncated)")

        prot["atoms"] = atom_count
        prot["residues"] = len(residue_ids)
        prot["chains"] = sorted(chains)
        prot["his_count"] = len(his_residues)
        prot["has_end_record"] = has_end

    except Exception as e:
        prot["warnings"].append(f"Error reading PDB: {e}")

    # --- Ligand validation ---
    result["ligand"] = _validate_ligand(ligand_path)

    return json.dumps(result, indent=2)


def _validate_ligand(ligand_path: str) -> dict:
    """Validate a ligand file (MOL2, SDF, or PDB)."""
    lig = {
        "path": ligand_path,
        "valid": False,
        "warnings": [],
    }

    if not os.path.exists(ligand_path):
        lig["warnings"].append(f"File not found: {ligand_path}")
        return lig

    if os.path.getsize(ligand_path) == 0:
        lig["warnings"].append("File is empty")
        return lig

    ext = os.path.splitext(ligand_path)[1].lower()
    lig["format"] = ext.lstrip(".")

    try:
        with open(ligand_path, "r") as f:
            content = f.read()

        if ext == ".mol2":
            if "@<TRIPOS>ATOM" not in content:
                lig["warnings"].append("Missing @<TRIPOS>ATOM section in MOL2 file")
            else:
                # Count atoms between @<TRIPOS>ATOM and next section
                atom_section = content.split("@<TRIPOS>ATOM")[1]
                next_section = atom_section.split("@<TRIPOS>")[0] if "@<TRIPOS>" in atom_section else atom_section
                atom_lines = [l for l in next_section.strip().splitlines() if l.strip()]
                lig["atoms"] = len(atom_lines)
                lig["valid"] = len(atom_lines) > 0

        elif ext == ".sdf":
            lines = content.splitlines()
            if len(lines) >= 4:
                # SDF counts line is line 4 (index 3): "NNN MMM ..."
                counts_line = lines[3].strip()
                match = re.match(r"(\d+)\s+(\d+)", counts_line)
                if match:
                    lig["atoms"] = int(match.group(1))
                    lig["valid"] = True
                else:
                    lig["warnings"].append(f"Cannot parse SDF counts line: '{counts_line}'")
            else:
                lig["warnings"].append("SDF file too short (fewer than 4 lines)")
            if "$$$$" not in content:
                lig["warnings"].append("Missing $$$$ end marker in SDF file")

        elif ext == ".pdb":
            atom_count = sum(1 for line in content.splitlines()
                             if line.startswith(("ATOM", "HETATM")))
            lig["atoms"] = atom_count
            lig["valid"] = atom_count > 0
            if atom_count == 0:
                lig["warnings"].append("No ATOM/HETATM records in ligand PDB")

        else:
            lig["warnings"].append(
                f"Unrecognized ligand format: '{ext}'. "
                f"Supported: .mol2, .sdf, .pdb"
            )

    except Exception as e:
        lig["warnings"].append(f"Error reading ligand file: {e}")

    return lig


# ==========================================================================
#  Tool 8: Validate Build Output
# ==========================================================================
@mcp.tool()
def validate_build_output(output_dir: str, build_mode: str = "md") -> str:
    """Post-build validation — verify that PRISM produced all expected output files.

    Checks for required files, parses topology for molecule summary, and reads
    system atom count from the coordinate file.

    Call this AFTER a build completes to confirm everything is in order.

    Args:
        output_dir: The PRISM output directory path (same as passed to build_system).
        build_mode: The build type to validate. Default: "md".
            Options: "md", "pmf", "rest2", "mmpbsa".
    """
    logger.info(f"validate_build_output: {output_dir} (mode={build_mode})")

    result = {
        "valid": False,
        "build_mode": build_mode,
        "files_found": [],
        "files_missing": [],
        "topology_summary": {},
        "system_atoms": None,
        "warnings": [],
    }

    if not os.path.isdir(output_dir):
        result["warnings"].append(f"Output directory not found: {output_dir}")
        return json.dumps(result, indent=2)

    # --- Define expected files per build mode ---
    expected_files = {}

    if build_mode == "md":
        expected_files = {
            "solv_ions.gro": "GMX_PROLIG_MD/solv_ions.gro",
            "topol.top": "GMX_PROLIG_MD/topol.top",
            "localrun.sh": "GMX_PROLIG_MD/localrun.sh",
            "em.mdp": "mdps/em.mdp",
            "nvt.mdp": "mdps/nvt.mdp",
            "npt.mdp": "mdps/npt.mdp",
            "md.mdp": "mdps/md.mdp",
        }
    elif build_mode == "pmf":
        expected_files = {
            "solv_ions.gro": "GMX_PROLIG_PMF/solv_ions.gro",
            "topol.top": "GMX_PROLIG_PMF/topol.top",
            "smd_run.sh": "GMX_PROLIG_PMF/smd_run.sh",
            "umbrella_run.sh": "GMX_PROLIG_PMF/umbrella_run.sh",
            "smd.mdp": "mdps/smd.mdp",
            "umbrella.mdp": "mdps/umbrella.mdp",
        }
    elif build_mode == "rest2":
        expected_files = {
            "rest2_run.sh": "GMX_PROLIG_REST2/rest2_run.sh",
        }
    elif build_mode == "mmpbsa":
        expected_files = {
            "mmpbsa_run.sh": "GMX_PROLIG_MMPBSA/mmpbsa_run.sh",
            "mmpbsa.in": "GMX_PROLIG_MMPBSA/mmpbsa.in",
            "topol.top": "GMX_PROLIG_MMPBSA/topol.top",
        }
    else:
        result["warnings"].append(
            f"Unknown build_mode '{build_mode}'. Options: md, pmf, rest2, mmpbsa"
        )
        return json.dumps(result, indent=2)

    # --- Check each expected file ---
    for label, rel_path in expected_files.items():
        full_path = os.path.join(output_dir, rel_path)
        if os.path.exists(full_path) and os.path.getsize(full_path) > 0:
            result["files_found"].append(label)
        else:
            result["files_missing"].append(label)

    # --- REST2: also check for replica topologies ---
    if build_mode == "rest2":
        rest2_dir = os.path.join(output_dir, "GMX_PROLIG_REST2")
        if os.path.isdir(rest2_dir):
            replica_tops = sorted(
                f for f in os.listdir(rest2_dir)
                if re.match(r"topol_r\d+\.top", f)
            )
            if replica_tops:
                result["files_found"].append(f"replica_topologies ({len(replica_tops)} files)")
            else:
                result["files_missing"].append("topol_r*.top (replica topologies)")

    # --- Parse topology for molecule summary ---
    topol_path = None
    for rel in expected_files.values():
        if rel.endswith("topol.top"):
            topol_path = os.path.join(output_dir, rel)
            break

    if topol_path and os.path.exists(topol_path):
        molecules = _parse_topology_molecules(topol_path)
        if molecules:
            result["topology_summary"]["molecules"] = molecules

    # --- Read atom count from .gro file ---
    gro_path = None
    for rel in expected_files.values():
        if rel.endswith(".gro"):
            gro_path = os.path.join(output_dir, rel)
            break

    if gro_path and os.path.exists(gro_path):
        try:
            with open(gro_path, "r") as f:
                f.readline()  # title
                atom_line = f.readline().strip()
                result["system_atoms"] = int(atom_line)
        except (ValueError, StopIteration):
            result["warnings"].append("Could not parse atom count from .gro file")

    # --- Check prism_config.yaml ---
    config_path = os.path.join(output_dir, "prism_config.yaml")
    if os.path.exists(config_path):
        result["files_found"].append("prism_config.yaml")
    else:
        result["warnings"].append("prism_config.yaml not found (non-critical)")

    # --- Final validity ---
    result["valid"] = len(result["files_missing"]) == 0

    return json.dumps(result, indent=2)


def _parse_topology_molecules(topol_path: str) -> list:
    """Parse the [ molecules ] section from a GROMACS topology file."""
    molecules = []
    in_molecules = False
    try:
        with open(topol_path, "r") as f:
            for line in f:
                stripped = line.strip()
                # Skip comments and empty lines
                if not stripped or stripped.startswith(";"):
                    continue
                # Detect section headers
                if stripped.startswith("["):
                    section = stripped.strip("[] \t").lower()
                    in_molecules = (section == "molecules")
                    continue
                if in_molecules:
                    # Skip #include or preprocessor directives
                    if stripped.startswith("#"):
                        continue
                    parts = stripped.split()
                    if len(parts) >= 2:
                        try:
                            molecules.append({
                                "name": parts[0],
                                "count": int(parts[1]),
                            })
                        except ValueError:
                            pass
    except Exception:
        pass
    return molecules


# ==========================================================================
#  Tool 9: Check Topology
# ==========================================================================
@mcp.tool()
def check_topology(topology_path: str) -> str:
    """Deep inspection of a GROMACS topology file for troubleshooting.

    Parses the topology to extract force field, water model, molecule list,
    and all #include directives. Checks whether included files exist on disk.
    Flags common issues like missing solvent, ions, or ligand entries.

    Use this tool when a build produces unexpected results or when the user
    reports simulation errors related to topology.

    Args:
        topology_path: Absolute path to the GROMACS topology file (topol.top).
    """
    logger.info(f"check_topology: {topology_path}")

    result = {
        "path": topology_path,
        "forcefield": None,
        "water_model": None,
        "molecules": [],
        "includes": [],
        "total_molecule_entries": 0,
        "warnings": [],
    }

    if not os.path.exists(topology_path):
        result["warnings"].append(f"Topology file not found: {topology_path}")
        return json.dumps(result, indent=2)

    topol_dir = os.path.dirname(os.path.abspath(topology_path))

    # Resolve GMXLIB for include path checking
    gmxlib = os.environ.get("GMXLIB", "")
    gmxdata = os.environ.get("GMXDATA", "")
    gmx_share = ""
    if gmxdata:
        gmx_share = os.path.join(gmxdata, "top")

    def _include_exists(inc_file: str) -> bool:
        """Check if an #include file can be found."""
        # Check relative to topology directory
        if os.path.exists(os.path.join(topol_dir, inc_file)):
            return True
        # Check GMXLIB
        if gmxlib and os.path.exists(os.path.join(gmxlib, inc_file)):
            return True
        # Check GMXDATA/top
        if gmx_share and os.path.exists(os.path.join(gmx_share, inc_file)):
            return True
        return False

    try:
        with open(topology_path, "r") as f:
            content = f.read()

        # --- Extract #include directives ---
        include_pattern = re.compile(r'#include\s+"([^"]+)"')
        includes = include_pattern.findall(content)

        for inc in includes:
            exists = _include_exists(inc)
            result["includes"].append({"file": inc, "exists": exists})
            if not exists:
                result["warnings"].append(f"Include file not found: {inc}")

        # --- Detect force field from first include ---
        if includes:
            first_inc = includes[0]
            ff_match = re.match(r"([^/]+)\.ff/", first_inc)
            if ff_match:
                result["forcefield"] = ff_match.group(1)

        # --- Detect water model ---
        water_models = {"tip3p", "tip4p", "tip4pew", "tip5p", "spc", "spce", "opc"}
        for inc in includes:
            basename = os.path.splitext(os.path.basename(inc))[0].lower()
            if basename in water_models:
                result["water_model"] = basename
                break

        # --- Parse [ molecules ] section ---
        result["molecules"] = _parse_topology_molecules(topology_path)
        result["total_molecule_entries"] = len(result["molecules"])

        # --- Check for common issues ---
        mol_names = {m["name"] for m in result["molecules"]}

        if not any(n in mol_names for n in ("SOL", "HOH", "WAT", "TIP3", "SPC")):
            result["warnings"].append("System has no solvent (no SOL/HOH/WAT entry in [ molecules ])")

        if not any(n in mol_names for n in ("NA", "CL", "K", "Na+", "Cl-", "SOD", "CLA")):
            result["warnings"].append("System has no ions (no NA/CL entry in [ molecules ])")

        if not any(n for n in mol_names if n not in
                   ("Protein", "Protein_chain_A", "SOL", "HOH", "WAT",
                    "NA", "CL", "K", "Na+", "Cl-", "SOD", "CLA", "TIP3", "SPC")
                   and not n.startswith("Protein")):
            result["warnings"].append("No ligand found in topology (no non-protein/solvent/ion molecule)")

    except Exception as e:
        result["warnings"].append(f"Error parsing topology: {e}")

    return json.dumps(result, indent=2)


# ==========================================================================
#  Prompt: Agent System Prompt
# ==========================================================================
def _load_agent_prompt() -> str:
    """Load the agent prompt markdown from disk."""
    prompt_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "prompts", "agent_prompt.md"
    )
    with open(prompt_path, "r") as f:
        return f.read()


@mcp.resource("prism://prompts/agent")
def get_agent_prompt() -> str:
    """System prompt for AI agents using PRISM tools."""
    return _load_agent_prompt()


@mcp.prompt()
def prism_agent() -> str:
    """PRISM assistant system prompt.

    Injects the full PRISM workflow guide into the conversation, covering:
    environment setup (GROMACS 2026.0, conda, AmberTools), recommended
    workflow (validate → build → verify), force field choices, build modes
    (MD/PMF/REST2/MMPBSA), and troubleshooting guidance.

    Use this prompt to turn the AI agent into a PRISM-aware MD simulation
    assistant that knows exactly how to use every PRISM tool.
    """
    return _load_agent_prompt()


# ==========================================================================
#  Entry point
# ==========================================================================
if __name__ == "__main__":
    logger.info("Starting PRISM MCP Server...")
    mcp.run(transport="stdio")
