#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Core trajectory processing functionality.
"""

import logging
from typing import Optional, List, Dict
from pathlib import Path

from .conversion import ConversionMixin
from .selection import SelectionMixin
from ...utils.environment import GromacsEnvironment

logger = logging.getLogger(__name__)


class TrajectoryProcessor(ConversionMixin, SelectionMixin):
    """
    Trajectory preprocessing utility for protein-ligand MD systems.

    This class provides methods to process trajectories using GROMACS trjconv
    to ensure proper periodic boundary condition handling and maintain
    protein-ligand contact integrity.
    """

    def __init__(self, topology_file: Optional[str] = None, gromacs_env: Optional[GromacsEnvironment] = None):
        """
        Initialize trajectory processor.

        Parameters
        ----------
        topology_file : str, optional
            Path to GROMACS topology file (.tpr)
        gromacs_env : GromacsEnvironment, optional
            Pre-initialized GROMACS environment
        """
        self.topology_file = topology_file

        # Check MDTraj availability for format conversion
        import importlib.util

        self.mdtraj_available = importlib.util.find_spec("mdtraj") is not None
        if not self.mdtraj_available:
            logger.warning("MDTraj not available. DCD conversion will not be supported.")

        # Initialize GROMACS environment
        if gromacs_env:
            self.gromacs_env = gromacs_env
        else:
            try:
                self.gromacs_env = GromacsEnvironment()
            except Exception as e:
                logger.error(f"Failed to initialize GROMACS environment: {e}")
                raise RuntimeError(
                    "GROMACS not available. Trajectory processing requires GROMACS installation.\n"
                    "RECOMMENDED INSTALLATION:\n"
                    "  mamba install -c conda-forge gromacs\n"
                    "  # OR: conda install -c conda-forge gromacs\n"
                    "  # OR: sudo apt-get install gromacs (Linux)\n"
                    "  # OR: brew install gromacs (macOS)"
                ) from e

        # Check trjconv availability
        self.trjconv_cmd = self._find_trjconv()

    def process_trajectory(
        self,
        input_trajectory: str,
        output_trajectory: str,
        topology_file: Optional[str] = None,
        center_selection: str = "auto",
        output_selection: str = "System",
        pbc_method: str = "atom",
        unit_cell: str = "compact",
        overwrite: bool = False,
    ) -> str:
        """
        Process a single trajectory to fix PBC artifacts and center protein-ligand complex.

        This method implements proper protein-ligand PBC processing following GROMACS best practices:
        1. First makes all molecules whole (-pbc whole)
        2. Then centers on ligand/protein complex (-pbc atom -center)

        This ensures the ligand stays in the binding pocket and maintains contact with receptor atoms.

        Parameters
        ----------
        input_trajectory : str
            Path to input trajectory file (.xtc, .dcd, .trr, etc.)
        output_trajectory : str
            Path to output trajectory file
        topology_file : str, optional
            Path to topology file (.tpr). Uses instance topology if not provided.
        center_selection : str
            Selection for centering. Options:
            - "auto": Automatically detect ligand, fallback to protein (default)
            - "ligand": Center on first non-standard residue (ligand)
            - "protein": Center on protein
            - Custom selection string
        output_selection : str
            Selection for output (default: "System")
        pbc_method : str
            PBC treatment method: "atom", "res", "mol" (default: "atom")
            - "atom": Most universal, keeps ligand in binding pocket
            - "res": For biological molecules with residues
            - "mol": Only if molecules are properly defined in topology
        unit_cell : str
            Unit cell representation: "compact", "tric", "rectangular" (default: "compact")
        overwrite : bool
            Whether to overwrite existing output file

        Returns
        -------
        str
            Path to processed trajectory file

        Examples
        --------
        >>> processor = TrajectoryProcessor("system.tpr")
        >>> processed = processor.process_trajectory("input.dcd", "output.xtc")
        >>> # Explicit ligand centering for protein-ligand complex
        >>> processed = processor.process_trajectory("input.dcd", "output.xtc",
        ...                                         center_selection="ligand")
        """
        # Validate inputs
        input_path = Path(input_trajectory)
        output_path = Path(output_trajectory)

        if not input_path.exists():
            raise FileNotFoundError(f"Input trajectory not found: {input_trajectory}")

        if output_path.exists() and not overwrite:
            logger.info(f"Output file exists, skipping: {output_trajectory}")
            return str(output_path)

        # Use provided topology or instance topology
        topo_file = topology_file or self.topology_file
        if not topo_file:
            raise ValueError("No topology file provided. Specify topology_file parameter or set instance topology.")

        if not Path(topo_file).exists():
            raise FileNotFoundError(f"Topology file not found: {topo_file}")

        # Handle DCD format conversion if needed
        trajectory_for_processing = input_trajectory
        temp_files_to_cleanup = []

        if input_path.suffix.lower() == ".dcd":
            logger.info("DCD format detected. Converting to XTC for GROMACS processing...")

            # Create temporary XTC file
            temp_xtc = output_path.parent / f"temp_{input_path.stem}.xtc"
            converted_file = self._convert_dcd_to_xtc(input_trajectory, topo_file, str(temp_xtc))
            trajectory_for_processing = converted_file
            temp_files_to_cleanup.append(temp_xtc)

        logger.info(f"Processing trajectory: {trajectory_for_processing}")
        logger.info(f"Output: {output_trajectory}")
        logger.info(f"Topology: {topo_file}")
        logger.info(f"PBC method: {pbc_method}, Unit cell: {unit_cell}")

        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Check if we need MDTraj for chain-specific selections
            needs_mdtraj = self._needs_mdtraj_processing(center_selection)

            if needs_mdtraj:
                # Use MDTraj for chain-specific PBC processing
                success = self._process_pbc_with_mdtraj(
                    topo_file,
                    trajectory_for_processing,
                    output_trajectory,
                    chain_selection=center_selection,
                    ligand_name=None,
                )
            else:
                # Use GROMACS two-step PBC processing
                success = self._process_pbc_two_step(
                    topo_file,
                    trajectory_for_processing,
                    output_trajectory,
                    center_selection,
                    output_selection,
                    pbc_method,
                    unit_cell,
                )

            if success and output_path.exists():
                logger.info(f"Successfully processed trajectory: {output_trajectory}")
                return str(output_path)
            else:
                raise RuntimeError(f"Failed to process trajectory: {input_trajectory}")

        finally:
            # Clean up temporary files
            for temp_file in temp_files_to_cleanup:
                if temp_file.exists():
                    temp_file.unlink()
                    logger.debug(f"Cleaned up temporary file: {temp_file}")

    def batch_process(
        self,
        input_trajectories: List[str],
        output_dir: str,
        topology_file: Optional[str] = None,
        output_suffix: str = "_processed",
        **kwargs,
    ) -> Dict[str, str]:
        """
        Process multiple trajectories in batch.

        Parameters
        ----------
        input_trajectories : List[str]
            List of input trajectory file paths
        output_dir : str
            Output directory for processed trajectories
        topology_file : str, optional
            Topology file path
        output_suffix : str
            Suffix to add to output filenames
        **kwargs
            Additional arguments passed to process_trajectory()

        Returns
        -------
        Dict[str, str]
            Mapping of input file -> output file paths

        Examples
        --------
        >>> processor = TrajectoryProcessor("system.tpr")
        >>> results = processor.batch_process(
        ...     ["repeat1.dcd", "repeat2.dcd", "repeat3.dcd"],
        ...     "processed_trajectories"
        ... )
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        results = {}

        for input_traj in input_trajectories:
            input_path = Path(input_traj)

            # Generate output filename
            base_name = input_path.stem
            extension = ".xtc"  # Always convert to XTC for consistency
            output_name = f"{base_name}{output_suffix}{extension}"
            output_file = output_path / output_name

            try:
                processed_file = self.process_trajectory(input_traj, str(output_file), topology_file, **kwargs)
                results[input_traj] = processed_file
                logger.info(f"✓ Processed: {input_traj} -> {processed_file}")

            except Exception as e:
                logger.error(f"✗ Failed to process {input_traj}: {e}")
                results[input_traj] = None

        logger.info(
            f"Batch processing complete: {len([v for v in results.values() if v])} / {len(input_trajectories)} successful"
        )
        return results

    def validate_processing(self, original_file: str, processed_file: str) -> Dict[str, any]:
        """
        Validate that trajectory processing was successful.

        Parameters
        ----------
        original_file : str
            Path to original trajectory
        processed_file : str
            Path to processed trajectory

        Returns
        -------
        Dict[str, any]
            Validation results including frame counts, file sizes, etc.
        """
        validation = {
            "original_exists": Path(original_file).exists(),
            "processed_exists": Path(processed_file).exists(),
            "original_size": 0,
            "processed_size": 0,
            "size_ratio": 0,
            "valid": False,
        }

        if validation["original_exists"]:
            validation["original_size"] = Path(original_file).stat().st_size

        if validation["processed_exists"]:
            validation["processed_size"] = Path(processed_file).stat().st_size

        if validation["original_size"] > 0 and validation["processed_size"] > 0:
            validation["size_ratio"] = validation["processed_size"] / validation["original_size"]
            # Consider valid if processed file is reasonable size (10-200% of original)
            validation["valid"] = 0.1 <= validation["size_ratio"] <= 2.0

        return validation

    def get_processing_info(self) -> Dict[str, any]:
        """Get information about the trajectory processing environment."""
        return {
            "gromacs_available": bool(self.gromacs_env),
            "trjconv_command": self.trjconv_cmd,
            "topology_file": self.topology_file,
            "gromacs_version": getattr(self.gromacs_env, "version", "unknown") if self.gromacs_env else None,
        }


def process_trajectory_simple(input_trajectory: str, output_trajectory: str, topology_file: str, **kwargs) -> str:
    """
    Simple function interface for trajectory processing.

    Parameters
    ----------
    input_trajectory : str
        Path to input trajectory file
    output_trajectory : str
        Path to output trajectory file
    topology_file : str
        Path to topology file (.tpr)
    **kwargs
        Additional processing options

    Returns
    -------
    str
        Path to processed trajectory

    Examples
    --------
    >>> processed = process_trajectory_simple(
    ...     "input.dcd", "output.xtc", "system.tpr"
    ... )
    """
    processor = TrajectoryProcessor(topology_file)
    return processor.process_trajectory(input_trajectory, output_trajectory, **kwargs)


def batch_process_trajectories(
    input_trajectories: List[str], output_dir: str, topology_file: str, **kwargs
) -> Dict[str, str]:
    """
    Simple function interface for batch trajectory processing.

    Parameters
    ----------
    input_trajectories : List[str]
        List of input trajectory files
    output_dir : str
        Output directory
    topology_file : str
        Path to topology file (.tpr)
    **kwargs
        Additional processing options

    Returns
    -------
    Dict[str, str]
        Mapping of input -> output files

    Examples
    --------
    >>> results = batch_process_trajectories(
    ...     ["repeat1.dcd", "repeat2.dcd"], "processed/", "system.tpr"
    ... )
    """
    processor = TrajectoryProcessor(topology_file)
    return processor.batch_process(input_trajectories, output_dir, **kwargs)
