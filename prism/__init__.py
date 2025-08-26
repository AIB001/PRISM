# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# """
# PRISM - Protein Receptor Interaction Simulation Modeler

# A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
# """

# __version__ = "1.0.0"
# __author__ = "PRISM Development Team"

# from .builder import PRISMBuilder
# from .core import PRISMSystem
# from .sim import model
# from .analysis import TrajAnalysis

# # High-level API functions
# def system(protein_path, ligand_path, config=None, **kwargs):
#     """
#     Create a protein-ligand system for MD simulation.

#     Parameters:
#     -----------
#     protein_path : str
#         Path to protein PDB file
#     ligand_path : str
#         Path to ligand file (MOL2/SDF)
#     config : str or dict, optional
#         Configuration file path or configuration dictionary
#     **kwargs : optional
#         Additional parameters (output_dir, ligand_forcefield, forcefield, etc.)

#     Returns:
#     --------
#     PRISMSystem
#         A system object with methods to build and run simulations

#     Examples:
#     ---------
#     >>> import prism as pm
#     >>> system = pm.system("protein.pdb", "ligand.mol2")
#     >>> system.build()
#     >>> system.generate_mdp_files()

#     >>> # With custom configuration
#     >>> system = pm.system("protein.pdb", "ligand.sdf",
#     ...                    config="config.yaml",
#     ...                    ligand_forcefield="openff")
#     >>> system.build()
#     """
#     return PRISMSystem(protein_path, ligand_path, config=config, **kwargs)

# def build_system(protein_path, ligand_path, output_dir="prism_output", **kwargs):
#     """
#     Build a complete protein-ligand system (one-step function).

#     Parameters:
#     -----------
#     protein_path : str
#         Path to protein PDB file
#     ligand_path : str
#         Path to ligand file (MOL2/SDF)
#     output_dir : str
#         Output directory for generated files
#     **kwargs : optional
#         Additional parameters (forcefield, water_model, ligand_forcefield, etc.)

#     Returns:
#     --------
#     str
#         Path to the output directory

#     Examples:
#     ---------
#     >>> import prism as pm
#     >>> output_path = pm.build_system("protein.pdb", "ligand.mol2")

#     >>> # With custom force fields
#     >>> output_path = pm.build_system("protein.pdb", "ligand.mol2",
#     ...                               forcefield="amber99sb-ildn",
#     ...                               water_model="tip4p",
#     ...                               ligand_forcefield="openff")
#     """
#     system_obj = PRISMSystem(protein_path, ligand_path, output_dir=output_dir, **kwargs)
#     return system_obj.build()

# def check_dependencies():
#     """
#     Check if all required dependencies are available.

#     Returns:
#     --------
#     dict
#         Dictionary showing availability of each dependency

#     Example:
#     --------
#     >>> import prism as pm
#     >>> deps = pm.check_dependencies()
#     >>> print(deps)
#     {'gromacs': True, 'pdbfixer': True, 'antechamber': True, 'openff': False}
#     """
#     from .utils.environment import GromacsEnvironment
#     import subprocess

#     dependencies = {
#         'gromacs': False,
#         'pdbfixer': False,
#         'antechamber': False,
#         'openff': False
#     }

#     # Check GROMACS
#     try:
#         env = GromacsEnvironment()
#         dependencies['gromacs'] = True
#     except:
#         pass

#     # Check pdbfixer
#     try:
#         import pdbfixer
#         dependencies['pdbfixer'] = True
#     except:
#         pass

#     # Check for AmberTools (antechamber)
#     try:
#         subprocess.run(['antechamber', '-h'], capture_output=True, check=True)
#         dependencies['antechamber'] = True
#     except:
#         pass

#     # Check OpenFF
#     try:
#         import openff.toolkit
#         dependencies['openff'] = True
#     except:
#         pass

#     return dependencies

# def list_forcefields():
#     """
#     List available force fields from GROMACS installation.

#     Returns:
#     --------
#     list
#         List of available force field names

#     Example:
#     --------
#     >>> import prism as pm
#     >>> ffs = pm.list_forcefields()
#     >>> for ff in ffs:
#     ...     print(ff)
#     amber99sb (index: 1)
#     amber99sb-ildn (index: 2)
#     amber14sb (index: 3)
#     """
#     from .utils.environment import GromacsEnvironment

#     try:
#         env = GromacsEnvironment()
#         return env.list_force_fields()
#     except Exception as e:
#         print(f"Error detecting force fields: {e}")
#         return []

# def get_version():
#     """Get PRISM version."""
#     return __version__

# # Export main classes and functions
# __all__ = [
#     "PRISMBuilder",
#     "PRISMSystem",
#     "system",
#     "build_system",
#     "model",
#     "check_dependencies",
#     "list_forcefields",
#     "get_version",
#     "__version__",
#     "TrajAnalysis",
# ]

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM - Protein Receptor Interaction Simulation Modeler

A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
"""

__version__ = "1.0.0"
__author__ = "PRISM Development Team"

from .builder import PRISMBuilder
from .core import PRISMSystem
from .sim import model
from .analysis import TrajAnalysis

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
            "Visualization module not available. "
            "Please ensure mdtraj is installed: pip install mdtraj"
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
        'gromacs': False,
        'pdbfixer': False,
        'antechamber': False,
        'openff': False,
        'mdtraj': False,
        'rdkit': False
    }

    # Check GROMACS
    try:
        env = GromacsEnvironment()
        dependencies['gromacs'] = True
    except:
        pass

    # Check pdbfixer
    try:
        import pdbfixer
        dependencies['pdbfixer'] = True
    except:
        pass

    # Check for AmberTools (antechamber)
    try:
        subprocess.run(['antechamber', '-h'], capture_output=True, check=True)
        dependencies['antechamber'] = True
    except:
        pass

    # Check OpenFF
    try:
        import openff.toolkit
        dependencies['openff'] = True
    except:
        pass
    
    # Check MDTraj (for visualization)
    try:
        import mdtraj
        dependencies['mdtraj'] = True
    except:
        pass
    
    # Check RDKit (optional for enhanced visualization)
    try:
        import rdkit
        dependencies['rdkit'] = True
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

# Lazy import for visualization classes to avoid import errors if dependencies missing
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
            "HTMLGenerator not available. "
            "Please ensure mdtraj is installed: pip install mdtraj"
        ) from e

# Export main classes and functions
__all__ = [
    "PRISMBuilder",
    "PRISMSystem",
    "system",
    "build_system",
    "model",
    "check_dependencies",
    "list_forcefields",
    "get_version",
    "__version__",
    "TrajAnalysis",
    "visualize_trajectory",
    "analyze_trajectory",
    "get_html_generator",
]