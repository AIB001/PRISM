#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectory Processing Utilities for PRISM

This module provides utilities for preprocessing MD trajectories to ensure proper
periodic boundary condition handling and maintain protein-ligand contact integrity.
Integrates with GROMACS trjconv for robust trajectory processing.
"""

import os
import subprocess
import tempfile
import shutil
import logging
from pathlib import Path
from typing import List, Union, Optional, Dict, Tuple
import time

from ...utils.environment import GromacsEnvironment

logger = logging.getLogger(__name__)


class TrajectoryProcessor:
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
        try:
            import mdtraj
            self.mdtraj_available = True
        except ImportError:
            self.mdtraj_available = False
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

    def _find_trjconv(self) -> str:
        """Find trjconv executable."""
        possible_names = ['gmx', 'gmx_mpi', 'gromacs']

        for cmd_base in possible_names:
            if shutil.which(cmd_base):
                # Test if it's a full GROMACS installation with trjconv
                try:
                    result = subprocess.run(
                        [cmd_base, 'trjconv', '-h'],
                        capture_output=True, text=True, timeout=10
                    )
                    if result.returncode == 0:
                        return cmd_base
                except (subprocess.TimeoutExpired, FileNotFoundError):
                    continue

        # Try legacy individual commands
        if shutil.which('trjconv'):
            return 'trjconv'

        raise RuntimeError(
            "GROMACS trjconv not found. Please install GROMACS.\n"
            "RECOMMENDED INSTALLATION:\n"
            "  mamba install -c conda-forge gromacs\n"
            "  # OR: conda install -c conda-forge gromacs"
        )

    def _convert_dcd_to_xtc(self, dcd_file: str, topology_file: str, output_file: str) -> str:
        """
        Convert DCD file to XTC format using MDTraj.

        Parameters
        ----------
        dcd_file : str
            Path to input DCD file
        topology_file : str
            Path to topology file (PDB format for MDTraj)
        output_file : str
            Path to output XTC file

        Returns
        -------
        str
            Path to converted XTC file
        """
        if not self.mdtraj_available:
            raise RuntimeError(
                "MDTraj not available for DCD conversion.\n"
                "RECOMMENDED INSTALLATION:\n"
                "  mamba install -c conda-forge mdtraj\n"
                "  # OR: conda install -c conda-forge mdtraj\n"
                "  # OR: pip install mdtraj"
            )

        import mdtraj as md

        logger.info(f"Converting DCD to XTC: {dcd_file} -> {output_file}")

        # Find appropriate topology file for MDTraj
        topo_for_mdtraj = self._get_mdtraj_topology(topology_file)

        try:
            # Load DCD with MDTraj
            traj = md.load(dcd_file, top=topo_for_mdtraj)
            logger.info(f"Loaded DCD: {traj.n_frames} frames, {traj.n_atoms} atoms")

            # Save as XTC
            traj.save_xtc(output_file)
            logger.info(f"Saved XTC: {output_file}")

            return output_file

        except Exception as e:
            logger.error(f"Failed to convert DCD to XTC: {e}")
            raise RuntimeError(f"DCD conversion failed: {e}") from e

    def _get_mdtraj_topology(self, topology_file: str) -> str:
        """
        Get appropriate topology file for MDTraj.

        MDTraj needs PDB format, so if we have TPR, we need to find/create PDB.
        """
        topo_path = Path(topology_file)

        if topo_path.suffix.lower() == '.pdb':
            return topology_file

        # Look for PDB file in same directory
        pdb_files = list(topo_path.parent.glob('*.pdb'))
        if pdb_files:
            logger.info(f"Using PDB topology for MDTraj: {pdb_files[0]}")
            return str(pdb_files[0])

        raise FileNotFoundError(
            f"No PDB topology file found for MDTraj conversion. "
            f"Please provide a PDB file in the same directory as {topology_file}"
        )

    def process_trajectory(self,
                         input_trajectory: str,
                         output_trajectory: str,
                         topology_file: Optional[str] = None,
                         center_selection: str = "auto",
                         output_selection: str = "System",
                         pbc_method: str = "atom",
                         unit_cell: str = "compact",
                         overwrite: bool = False) -> str:
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

        if input_path.suffix.lower() == '.dcd':
            logger.info("DCD format detected. Converting to XTC for GROMACS processing...")

            # Create temporary XTC file
            temp_xtc = output_path.parent / f"temp_{input_path.stem}.xtc"
            converted_file = self._convert_dcd_to_xtc(
                input_trajectory, topo_file, str(temp_xtc)
            )
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
                    topo_file, trajectory_for_processing, output_trajectory,
                    chain_selection=center_selection, ligand_name=None
                )
            else:
                # Use GROMACS two-step PBC processing
                success = self._process_pbc_two_step(
                    topo_file, trajectory_for_processing, output_trajectory,
                    center_selection, output_selection, pbc_method, unit_cell
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

    def _process_pbc_two_step(self, topology: str, input_traj: str, output_traj: str,
                            center_selection: str, output_selection: str,
                            pbc_method: str, unit_cell: str) -> bool:
        """
        Two-step PBC processing for protein-ligand complexes.

        Step 1: Make all molecules whole (-pbc whole)
        Step 2: Center on ligand/protein and apply PBC (-pbc atom -center)
        """
        output_path = Path(output_traj)
        temp_whole_file = output_path.parent / f"temp_whole_{output_path.stem}.xtc"

        try:
            # Step 1: Make molecules whole
            logger.info("Step 1: Making molecules whole...")
            cmd_whole = self._build_trjconv_command(
                topology, input_traj, str(temp_whole_file), "whole", unit_cell, center=False
            )
            input_text_whole = "0\n"  # System selection for output
            success_whole = self._run_trjconv(cmd_whole, input_text_whole, str(temp_whole_file))

            if not success_whole or not temp_whole_file.exists():
                logger.error("Step 1 failed: Could not make molecules whole")
                return False

            # Step 2: Center and apply PBC
            logger.info("Step 2: Centering and applying PBC...")
            resolved_center = self._resolve_center_selection(center_selection, topology)
            cmd_center = self._build_trjconv_command(
                topology, str(temp_whole_file), output_traj, pbc_method, unit_cell, center=True
            )
            input_text_center = self._prepare_selections(resolved_center, output_selection)
            success_center = self._run_trjconv(cmd_center, input_text_center, output_traj)

            return success_center

        finally:
            # Clean up temporary file
            if temp_whole_file.exists():
                temp_whole_file.unlink()
                logger.debug(f"Cleaned up temporary file: {temp_whole_file}")

    def _process_pbc_with_mdtraj(self, topology: str, input_traj: str, output_traj: str,
                               chain_selection: Optional[str] = None, ligand_name: Optional[str] = None) -> bool:
        """
        Process PBC using MDTraj with chain-specific centering.

        Ensures functional chain and ligand stay together as a unit.
        This method is used when GROMACS cannot handle complex selections.

        Parameters
        ----------
        topology : str
            Path to topology file (PDB format)
        input_traj : str
            Path to input trajectory
        output_traj : str
            Path to output trajectory
        chain_selection : str, optional
            Chain selection (e.g., "P", "0", "A")
        ligand_name : str, optional
            Ligand residue name

        Returns
        -------
        bool
            True if processing succeeded
        """
        try:
            import mdtraj as md

            # Find PDB topology file for MDTraj
            topo_path = Path(topology)
            if topo_path.suffix.lower() != '.pdb':
                # Look for PDB file in same directory
                pdb_files = list(topo_path.parent.glob('*.pdb'))
                if not pdb_files:
                    logger.error("MDTraj requires PDB topology file, none found")
                    return False
                topology = str(pdb_files[0])

            logger.info("Using MDTraj for chain-specific PBC processing...")

            # Load trajectory
            logger.info(f"Loading trajectory with MDTraj: {input_traj}")
            traj = md.load(input_traj, top=topology)
            logger.info(f"Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")

            # Select functional chain
            chain_atoms = self._select_functional_chain(traj, chain_selection)
            if chain_atoms is None or len(chain_atoms) == 0:
                logger.warning("No functional chain found, using all protein atoms")
                chain_atoms = traj.topology.select("protein")

            # Select ligand atoms
            ligand_atoms = self._select_ligand_atoms(traj, ligand_name)

            # Create anchor molecules: chain + ligand together
            # Use numpy arrays directly - MDTraj expects atom index arrays
            anchor_molecules = []
            if len(chain_atoms) > 0:
                anchor_molecules.append(chain_atoms)
                logger.info(f"Anchoring on {len(chain_atoms)} chain atoms")

            if ligand_atoms is not None and len(ligand_atoms) > 0:
                anchor_molecules.append(ligand_atoms)
                logger.info(f"Anchoring on {len(ligand_atoms)} ligand atoms")

            if not anchor_molecules:
                logger.error("No anchor molecules found for PBC processing")
                return False

            # Apply PBC correction with two-step approach
            logger.info("Applying MDTraj PBC correction...")

            # Step 1: Make molecules whole
            logger.info("Step 1: Making molecules whole...")
            traj_whole = traj.make_molecules_whole()

            # Step 2: Center the entire trajectory (simple approach)
            logger.info("Step 2: Centering trajectory coordinates")
            traj_centered = traj_whole.center_coordinates(mass_weighted=False)

            # Save processed trajectory
            logger.info(f"Saving processed trajectory: {output_traj}")
            traj_centered.save_xtc(output_traj)

            logger.info("MDTraj PBC processing completed successfully")
            return True

        except ImportError:
            logger.error("MDTraj not available for chain-specific PBC processing")
            return False
        except Exception as e:
            logger.error(f"MDTraj PBC processing failed: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _select_functional_chain(self, traj, chain_selection: Optional[str] = None):
        """
        Select atoms from functional protein chain.

        Parameters
        ----------
        traj : mdtraj.Trajectory
            Trajectory object
        chain_selection : str, optional
            Chain specification (e.g., "P", "0", "A")

        Returns
        -------
        numpy.ndarray
            Atom indices for selected chain
        """
        try:
            if chain_selection is None:
                # Auto-detect: find chain with most protein atoms
                chain_atoms = self._detect_largest_protein_chain(traj)
                if chain_atoms is not None:
                    logger.info("Auto-detected largest protein chain")
                    return chain_atoms

            # Try to interpret chain selection
            if chain_selection.isdigit():
                # Numeric chain ID
                chainid = int(chain_selection)
                selection = f"protein and chainid {chainid}"
                logger.info(f"Selecting chain by ID: {selection}")
                return traj.topology.select(selection)

            else:
                # Chain letter - try to map to chainid
                chainid = self._map_chain_letter_to_id(traj, chain_selection)
                if chainid is not None:
                    selection = f"protein and chainid {chainid}"
                    logger.info(f"Mapped chain {chain_selection} to chainid {chainid}")
                    return traj.topology.select(selection)

            # Fallback: all protein atoms
            logger.warning(f"Could not resolve chain selection '{chain_selection}', using all protein")
            return traj.topology.select("protein")

        except Exception as e:
            logger.warning(f"Chain selection failed: {e}, using all protein")
            return traj.topology.select("protein")

    def _select_ligand_atoms(self, traj, ligand_name: Optional[str] = None):
        """
        Select ligand atoms using PRISM ligand detection.

        Parameters
        ----------
        traj : mdtraj.Trajectory
            Trajectory object
        ligand_name : str, optional
            Ligand residue name

        Returns
        -------
        numpy.ndarray or None
            Ligand atom indices
        """
        try:
            from ...utils.ligand import identify_ligand_residue

            # Use PRISM ligand detection
            ligand = identify_ligand_residue(traj, ligand_name)
            if ligand:
                # Select all atoms in ligand residue
                ligand_atoms = traj.topology.select(f"resname {ligand.name}")
                logger.info(f"Selected {len(ligand_atoms)} atoms from ligand {ligand.name}")
                return ligand_atoms
            else:
                logger.info("No ligand found")
                return None

        except Exception as e:
            logger.warning(f"Ligand selection failed: {e}")
            return None

    def _detect_largest_protein_chain(self, traj):
        """Find the chain with the most protein atoms."""
        try:
            chain_sizes = {}
            for chain in traj.topology.chains:
                protein_atoms = traj.topology.select(f"protein and chainid {chain.index}")
                if len(protein_atoms) > 0:
                    chain_sizes[chain.index] = len(protein_atoms)

            if chain_sizes:
                largest_chain = max(chain_sizes, key=chain_sizes.get)
                logger.info(f"Largest protein chain: {largest_chain} with {chain_sizes[largest_chain]} atoms")
                return traj.topology.select(f"protein and chainid {largest_chain}")

        except Exception as e:
            logger.warning(f"Could not detect largest chain: {e}")

        return None

    def _map_chain_letter_to_id(self, traj, chain_letter: str):
        """Map chain letter (A, B, P, etc.) to MDTraj chainid."""
        try:
            # Simple mapping: assume alphabetical order
            chain_map = {
                'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
                'P': None,  # Will try to detect
                'Q': None   # Will try to detect
            }

            if chain_letter in chain_map:
                chainid = chain_map[chain_letter]
                if chainid is not None:
                    # Verify chain exists
                    test_atoms = traj.topology.select(f"protein and chainid {chainid}")
                    if len(test_atoms) > 0:
                        return chainid

            # For special cases like P, Q, try to find by context
            if chain_letter in ['P', 'Q']:
                return self._find_special_chain(traj, chain_letter)

        except Exception as e:
            logger.warning(f"Chain letter mapping failed: {e}")

        return None

    def _find_special_chain(self, traj, chain_letter: str):
        """Find special chains like P (primary) or Q by context."""
        try:
            # Look for chains with ligand nearby
            from ...utils.ligand import identify_ligand_residue

            ligand = identify_ligand_residue(traj)
            if ligand:
                # Find chain closest to ligand
                for chain in traj.topology.chains:
                    chain_atoms = traj.topology.select(f"protein and chainid {chain.index}")
                    if len(chain_atoms) > 100:  # Substantial protein chain
                        logger.info(f"Found substantial protein chain {chain.index}")
                        return chain.index

        except Exception:
            pass

        return None

    def _needs_mdtraj_processing(self, center_selection: str) -> bool:
        """
        Determine if we need MDTraj for complex chain selections.

        Returns True if the selection requires MDTraj capabilities.
        """
        if not center_selection:
            return False

        # Chain-specific selections that GROMACS can't handle
        chain_indicators = [
            'chain',    # "chain P"
            'chainid',  # "chainid 0"
            'segid',    # "segid PR1"
        ]

        selection_lower = center_selection.lower()
        for indicator in chain_indicators:
            if indicator in selection_lower:
                logger.info(f"Detected chain selection '{center_selection}', using MDTraj")
                return True

        # Single letter chain selections
        if len(center_selection) == 1 and center_selection.isalpha():
            logger.info(f"Detected chain letter '{center_selection}', using MDTraj")
            return True

        # Numeric chain IDs
        if center_selection.isdigit():
            logger.info(f"Detected numeric chain ID '{center_selection}', using MDTraj")
            return True

        return False

    def _resolve_center_selection(self, center_selection: str, topology: str) -> str:
        """
        Resolve automatic center selection for protein-ligand systems.

        Returns GROMACS group name as string.
        """
        if center_selection == "auto":
            # Try to detect ligand first, fallback to protein
            ligand_name = self._detect_ligand_group(topology)
            if ligand_name:
                logger.info(f"Auto-detected ligand for centering: {ligand_name}")
                return ligand_name
            else:
                logger.info("No ligand detected, centering on protein")
                return "protein"
        elif center_selection == "ligand":
            ligand_name = self._detect_ligand_group(topology)
            if ligand_name:
                logger.info(f"Using ligand for centering: {ligand_name}")
                return ligand_name
            else:
                logger.warning("Ligand selection requested but no ligand detected, using protein")
                return "protein"
        else:
            return center_selection

    def _detect_ligand_group(self, topology: str) -> Optional[str]:
        """
        Detect ligand group name from GROMACS groups.

        Returns ligand group name (e.g., "LIG").
        """
        try:
            # Use PRISM ligand detection first
            from ...utils.ligand import identify_ligand_residue
            import mdtraj as md

            topo_path = Path(topology)
            pdb_files = list(topo_path.parent.glob('*.pdb'))
            if pdb_files:
                traj = md.load(str(pdb_files[0]))
                ligand = identify_ligand_residue(traj)
                if ligand:
                    logger.info(f"Detected ligand: {ligand.name}")
                    return ligand.name

        except Exception as e:
            logger.debug(f"Could not detect ligand: {e}")

        return None

    def _build_trjconv_command(self, topology: str, input_traj: str, output_traj: str,
                              pbc_method: str, unit_cell: str, center: bool = True) -> List[str]:
        """Build trjconv command with appropriate options."""
        cmd = []

        # Handle different GROMACS command formats
        if self.trjconv_cmd == 'trjconv':
            # Legacy individual command
            cmd = ['trjconv']
        else:
            # Modern gmx interface
            cmd = [self.trjconv_cmd, 'trjconv']

        # Add required parameters
        cmd.extend([
            '-s', topology,
            '-f', input_traj,
            '-o', output_traj,
            '-pbc', pbc_method,
            '-ur', unit_cell
        ])

        # Add centering only when requested (Step 2)
        if center:
            cmd.append('-center')

        return cmd

    def _prepare_selections(self, center_selection: str, output_selection: str) -> str:
        """
        Prepare input text for interactive selections.

        GROMACS will show available groups and we select by group name.
        """
        # For ligand names, use the name directly (GROMACS will find the group)
        if center_selection in ['LIG', 'MOL', 'UNL', 'DRG', 'INH', 'SUB', 'HET']:
            center_input = center_selection
            logger.info(f"Using ligand group '{center_selection}' for centering")
        elif center_selection == "protein":
            center_input = "Protein"  # Standard GROMACS protein group name
        else:
            center_input = center_selection

        # Output selection - always use System
        output_input = "System" if output_selection.lower() == "system" else "0"

        return f"{center_input}\n{output_input}\n"

    def _run_trjconv(self, cmd: List[str], input_text: str, output_file: str) -> bool:
        """Run trjconv command with error handling."""
        try:
            logger.debug(f"Running command: {' '.join(cmd)}")
            logger.debug(f"Input selections: {repr(input_text)}")

            # Run trjconv with input piped
            process = subprocess.run(
                cmd,
                input=input_text,
                text=True,
                capture_output=True,
                timeout=1800  # 30 minutes timeout for large files
            )

            if process.returncode == 0:
                logger.debug("trjconv completed successfully")
                logger.debug(f"stderr: {process.stderr}")
                return True
            else:
                logger.error(f"trjconv failed with return code {process.returncode}")
                logger.error(f"stdout: {process.stdout}")
                logger.error(f"stderr: {process.stderr}")
                return False

        except subprocess.TimeoutExpired:
            logger.error("trjconv command timed out (30 minutes)")
            return False
        except Exception as e:
            logger.error(f"Error running trjconv: {e}")
            return False

    def batch_process(self,
                     input_trajectories: List[str],
                     output_dir: str,
                     topology_file: Optional[str] = None,
                     output_suffix: str = "_processed",
                     **kwargs) -> Dict[str, str]:
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
                processed_file = self.process_trajectory(
                    input_traj, str(output_file), topology_file, **kwargs
                )
                results[input_traj] = processed_file
                logger.info(f"✓ Processed: {input_traj} -> {processed_file}")

            except Exception as e:
                logger.error(f"✗ Failed to process {input_traj}: {e}")
                results[input_traj] = None

        logger.info(f"Batch processing complete: {len([v for v in results.values() if v])} / {len(input_trajectories)} successful")
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
            'original_exists': Path(original_file).exists(),
            'processed_exists': Path(processed_file).exists(),
            'original_size': 0,
            'processed_size': 0,
            'size_ratio': 0,
            'valid': False
        }

        if validation['original_exists']:
            validation['original_size'] = Path(original_file).stat().st_size

        if validation['processed_exists']:
            validation['processed_size'] = Path(processed_file).stat().st_size

        if validation['original_size'] > 0 and validation['processed_size'] > 0:
            validation['size_ratio'] = validation['processed_size'] / validation['original_size']
            # Consider valid if processed file is reasonable size (10-200% of original)
            validation['valid'] = 0.1 <= validation['size_ratio'] <= 2.0

        return validation

    def get_processing_info(self) -> Dict[str, any]:
        """Get information about the trajectory processing environment."""
        return {
            'gromacs_available': bool(self.gromacs_env),
            'trjconv_command': self.trjconv_cmd,
            'topology_file': self.topology_file,
            'gromacs_version': getattr(self.gromacs_env, 'version', 'unknown') if self.gromacs_env else None,
        }


def process_trajectory_simple(input_trajectory: str,
                            output_trajectory: str,
                            topology_file: str,
                            **kwargs) -> str:
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


def batch_process_trajectories(input_trajectories: List[str],
                             output_dir: str,
                             topology_file: str,
                             **kwargs) -> Dict[str, str]:
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