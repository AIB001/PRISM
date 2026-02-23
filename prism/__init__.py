#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM - Protein Receptor Interaction Simulation Modeler

A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
"""

import os
import multiprocessing

# Enable automatic OpenMP parallelization for trajectory analysis
# Set OMP_NUM_THREADS to use all available CPU cores by default
if "OMP_NUM_THREADS" not in os.environ:
    n_cores = multiprocessing.cpu_count()
    os.environ["OMP_NUM_THREADS"] = str(n_cores)

__version__ = "1.2.0"
__author__ = "PRISM Development Team"

from .builder import PRISMBuilder
from .core import PRISMSystem
from .sim import model
from .analysis import IntegratedProteinLigandAnalyzer as TrajAnalysis

try:
    from . import modeling
except ImportError:
    modeling = None

try:
    from . import pmf
except ImportError:
    pmf = None

try:
    from . import gaussian
except ImportError:
    gaussian = None

try:
    from . import rest2
except ImportError:
    rest2 = None


# High-level API functions
def system(protein_path, ligand_path, config=None, **kwargs):
    """
    Create a protein-ligand system for MD simulation.

    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)
    config : str or dict, optional
        Configuration file path or configuration dictionary
    **kwargs : optional
        Additional parameters (output_dir, ligand_forcefield, forcefield, etc.)

    Returns:
    --------
    PRISMSystem
        A system object with methods to build and run simulations

    Examples:
    ---------
    >>> import prism as pm
    >>> system = pm.system("protein.pdb", "ligand.mol2")
    >>> system.build()
    >>> system.generate_mdp_files()

    >>> # With custom configuration
    >>> system = pm.system("protein.pdb", "ligand.sdf",
    ...                    config="config.yaml",
    ...                    ligand_forcefield="openff")
    >>> system.build()
    """
    return PRISMSystem(protein_path, ligand_path, config=config, **kwargs)


def build_system(protein_path, ligand_path, output_dir="prism_output", **kwargs):
    """
    Build a complete protein-ligand system (one-step function).

    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)
    output_dir : str
        Output directory for generated files
    **kwargs : optional
        Additional parameters (forcefield, water_model, ligand_forcefield, etc.)

    Returns:
    --------
    str
        Path to the output directory

    Examples:
    ---------
    >>> import prism as pm
    >>> output_path = pm.build_system("protein.pdb", "ligand.mol2")

    >>> # With custom force fields
    >>> output_path = pm.build_system("protein.pdb", "ligand.mol2",
    ...                               forcefield="amber99sb-ildn",
    ...                               water_model="tip4p",
    ...                               ligand_forcefield="openff")
    """
    system_obj = PRISMSystem(protein_path, ligand_path, output_dir=output_dir, **kwargs)
    return system_obj.build()


def visualize_trajectory(trajectory, topology, ligand, output="contact_analysis.html", **kwargs):
    """
    Generate interactive HTML visualization of protein-ligand contacts from MD trajectory.

    This function creates an interactive HTML file that visualizes the contact
    patterns between protein and ligand throughout the trajectory.

    Parameters:
    -----------
    trajectory : str
        Path to trajectory file (.xtc, .dcd, .trr, etc.)
    topology : str
        Path to topology/structure file (.pdb, .gro, etc.)
    ligand : str
        Path to ligand structure file (.sdf, .mol, .mol2)
    output : str, optional
        Output HTML file path (default: "contact_analysis.html")
    **kwargs : optional
        Additional parameters for visualization

    Returns:
    --------
    str
        Path to the generated HTML file

    Examples:
    ---------
    >>> import prism as pm
    >>> # Basic visualization
    >>> pm.visualize_trajectory("md.xtc", "system.gro", "ligand.sdf")

    >>> # Custom output file
    >>> pm.visualize_trajectory("trajectory.xtc", "topology.pdb", "ligand.mol2",
    ...                        output="my_analysis.html")

    Notes:
    ------
    The generated HTML file is self-contained and can be opened in any modern
    web browser. It includes interactive features like:
    - Drag to reposition residues
    - Zoom in/out
    - Toggle between 2D and 3D views
    - Export high-resolution images
    """
    try:
        from .analysis.visualization import generate_html
    except ImportError as e:
        raise ImportError(
            "Visualization module not available. Missing dependencies: mdtraj.\n"
            "RECOMMENDED INSTALLATION:\n"
            "  mamba install -c conda-forge mdtraj\n"
            "  # OR: conda install -c conda-forge mdtraj\n"
            "  # OR: pip install mdtraj"
        ) from e

    return generate_html(trajectory, topology, ligand, output, **kwargs)


def analyze_trajectory(topology, trajectory, ligand_resname="LIG", output_dir="analysis_results"):
    """
    Perform comprehensive trajectory analysis with visualization.

    This is a convenience function that combines trajectory analysis with
    HTML visualization generation.

    Parameters:
    -----------
    topology : str
        Path to topology/structure file (.pdb, .gro, etc.)
    trajectory : str
        Path to trajectory file (.xtc, .dcd, .trr, etc.)
    ligand_resname : str, optional
        Residue name of ligand (default: "LIG")
    output_dir : str, optional
        Output directory for analysis results (default: "analysis_results")

    Returns:
    --------
    TrajAnalysis
        Analysis object with results

    Examples:
    ---------
    >>> import prism as pm
    >>> analysis = pm.analyze_trajectory("system.gro", "md.xtc")
    >>> # Results are saved in analysis_results/ directory

    >>> # Custom ligand residue name
    >>> analysis = pm.analyze_trajectory("complex.pdb", "trajectory.dcd",
    ...                                  ligand_resname="MOL")
    """
    traj = TrajAnalysis(topology, trajectory, ligand_resname=ligand_resname)
    traj.analyze_all(output_dir=output_dir)
    return traj


def check_dependencies():
    """
    Check if all required dependencies are available.

    Returns:
    --------
    dict
        Dictionary showing availability of each dependency

    Example:
    --------
    >>> import prism as pm
    >>> deps = pm.check_dependencies()
    >>> print(deps)
    {'gromacs': True, 'pdbfixer': True, 'antechamber': True, 'openff': False, 'mdtraj': True}
    """
    from .utils.environment import GromacsEnvironment
    import subprocess

    dependencies = {
        "gromacs": False,
        "pdbfixer": False,
        "antechamber": False,
        "openff": False,
        "mdtraj": False,
        "rdkit": False,
    }

    # Check GROMACS
    try:
        env = GromacsEnvironment()
        dependencies["gromacs"] = True
    except:
        pass

    # Check pdbfixer
    try:
        import pdbfixer

        dependencies["pdbfixer"] = True
    except:
        pass

    # Check for AmberTools (antechamber)
    try:
        subprocess.run(["antechamber", "-h"], capture_output=True, check=True)
        dependencies["antechamber"] = True
    except:
        pass

    # Check OpenFF
    try:
        import openff.toolkit

        dependencies["openff"] = True
    except:
        pass

    # Check MDTraj (for visualization)
    try:
        import mdtraj

        dependencies["mdtraj"] = True
    except:
        pass

    # Check RDKit (optional for enhanced visualization)
    try:
        import rdkit

        dependencies["rdkit"] = True
    except:
        pass

    return dependencies


def list_forcefields():
    """
    List available force fields from GROMACS installation.

    Returns:
    --------
    list
        List of available force field names

    Example:
    --------
    >>> import prism as pm
    >>> ffs = pm.list_forcefields()
    >>> for ff in ffs:
    ...     print(ff)
    amber99sb (index: 1)
    amber99sb-ildn (index: 2)
    amber14sb (index: 3)
    """
    from .utils.environment import GromacsEnvironment

    try:
        env = GromacsEnvironment()
        return env.list_force_fields()
    except Exception as e:
        print(f"Error detecting force fields: {e}")
        return []


def get_version():
    """Get PRISM version."""
    return __version__


def get_html_generator():
    """
    Get HTMLGenerator class for advanced visualization usage.

    Returns:
    --------
    HTMLGenerator
        The HTML generator class

    Example:
    --------
    >>> import prism as pm
    >>> HTMLGen = pm.get_html_generator()
    >>> generator = HTMLGen("trajectory.xtc", "topology.gro", "ligand.sdf")
    >>> generator.analyze()
    >>> generator.generate("output.html")
    """
    try:
        from .analysis.visualization import HTMLGenerator

        return HTMLGenerator
    except ImportError as e:
        raise ImportError(
            "HTMLGenerator not available. Missing dependencies: mdtraj.\n"
            "RECOMMENDED INSTALLATION:\n"
            "  mamba install -c conda-forge mdtraj\n"
            "  # OR: conda install -c conda-forge mdtraj\n"
            "  # OR: pip install mdtraj"
        ) from e


def process_trajectory(input_trajectory, output_trajectory, topology_file, **kwargs):
    """
    Process a trajectory to fix PBC artifacts and center protein-ligand complex.

    This function applies GROMACS trjconv processing to fix periodic boundary
    condition artifacts and center the system. Automatically handles DCD format
    conversion if needed.

    Parameters:
    -----------
    input_trajectory : str
        Path to input trajectory file (.dcd, .xtc, .trr, etc.)
    output_trajectory : str
        Path to output processed trajectory file (.xtc recommended)
    topology_file : str
        Path to topology file (.tpr, .gro, .pdb)
    center_selection : str, optional
        Selection for centering (default: "Protein")
    output_selection : str, optional
        Selection for output (default: "System")
    pbc_method : str, optional
        PBC method (default: "mol")
    unit_cell : str, optional
        Unit cell handling (default: "compact")
    overwrite : bool, optional
        Whether to overwrite existing files (default: False)

    Returns:
    --------
    str
        Path to the processed trajectory file

    Examples:
    ---------
    >>> import prism as pm
    >>> processed = pm.process_trajectory("repeat1.dcd", "repeat1_processed.xtc", "system.tpr")
    >>> print(f"Processed trajectory saved to: {processed}")

    >>> # With custom settings
    >>> processed = pm.process_trajectory("trajectory.xtc", "centered.xtc", "system.tpr",
    ...                                  center_selection="Protein",
    ...                                  pbc_method="atom",
    ...                                  overwrite=True)
    """
    from .analysis.trajectory_processor import TrajectoryProcessor

    processor = TrajectoryProcessor()
    return processor.process_trajectory(
        input_trajectory=input_trajectory, output_trajectory=output_trajectory, topology_file=topology_file, **kwargs
    )


def batch_process_trajectories(input_trajectories, output_dir, topology_file, **kwargs):
    """
    Process multiple trajectories in batch to fix PBC artifacts.

    This function processes multiple trajectory files using the same settings,
    applying PBC corrections and centering for each trajectory.

    Parameters:
    -----------
    input_trajectories : list of str
        List of input trajectory file paths
    output_dir : str
        Directory for output processed trajectories
    topology_file : str
        Path to topology file (.tpr, .gro, .pdb)
    prefix : str, optional
        Prefix for output filenames (default: "processed_")
    **kwargs : dict
        Additional arguments passed to process_trajectory

    Returns:
    --------
    list of str
        Paths to the processed trajectory files

    Examples:
    ---------
    >>> import prism as pm
    >>> inputs = ["repeat1.dcd", "repeat2.dcd", "repeat3.dcd"]
    >>> outputs = pm.batch_process_trajectories(inputs, "processed_trajs", "system.tpr")
    >>> print(f"Processed {len(outputs)} trajectories")
    """
    from .analysis.trajectory_processor import TrajectoryProcessor

    processor = TrajectoryProcessor()
    return processor.batch_process(
        input_trajectories=input_trajectories, output_dir=output_dir, topology_file=topology_file, **kwargs
    )


# Export main classes and functions
__all__ = [
    # Core classes
    "PRISMBuilder",
    "PRISMSystem",
    "TrajAnalysis",
    # High-level API functions
    "system",
    "build_system",
    "model",
    # Analysis functions
    "analyze_trajectory",
    "visualize_trajectory",
    "get_html_generator",
    "process_trajectory",
    "batch_process_trajectories",
    # Utility functions
    "check_dependencies",
    "list_forcefields",
    "get_version",
    # Module imports
    "modeling",
    "pmf",
    "gaussian",
    "rest2",
    # Version info
    "__version__",
]
