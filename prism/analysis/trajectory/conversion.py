#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Trajectory format conversion utilities.
"""

import subprocess
import shutil
import logging
from pathlib import Path
from typing import List, Optional


logger = logging.getLogger(__name__)


class ConversionMixin:
    """Mixin for trajectory format conversion operations."""

    def _find_trjconv(self) -> str:
        """Find trjconv executable."""
        possible_names = ["gmx", "gmx_mpi", "gromacs"]

        for cmd_base in possible_names:
            if shutil.which(cmd_base):
                # Test if it's a full GROMACS installation with trjconv
                try:
                    result = subprocess.run([cmd_base, "trjconv", "-h"], capture_output=True, text=True, timeout=10)
                    if result.returncode == 0:
                        return cmd_base
                except (subprocess.TimeoutExpired, FileNotFoundError):
                    continue

        # Try legacy individual commands
        if shutil.which("trjconv"):
            return "trjconv"

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

        if topo_path.suffix.lower() == ".pdb":
            return topology_file

        # Look for PDB file in same directory
        pdb_files = list(topo_path.parent.glob("*.pdb"))
        if pdb_files:
            logger.info(f"Using PDB topology for MDTraj: {pdb_files[0]}")
            return str(pdb_files[0])

        raise FileNotFoundError(
            f"No PDB topology file found for MDTraj conversion. "
            f"Please provide a PDB file in the same directory as {topology_file}"
        )

    def _process_pbc_two_step(
        self,
        topology: str,
        input_traj: str,
        output_traj: str,
        center_selection: str,
        output_selection: str,
        pbc_method: str,
        unit_cell: str,
    ) -> bool:
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

    def _process_pbc_with_mdtraj(
        self,
        topology: str,
        input_traj: str,
        output_traj: str,
        chain_selection: Optional[str] = None,
        ligand_name: Optional[str] = None,
    ) -> bool:
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
            if topo_path.suffix.lower() != ".pdb":
                # Look for PDB file in same directory
                pdb_files = list(topo_path.parent.glob("*.pdb"))
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

    def _build_trjconv_command(
        self, topology: str, input_traj: str, output_traj: str, pbc_method: str, unit_cell: str, center: bool = True
    ) -> List[str]:
        """Build trjconv command with appropriate options."""
        cmd = []

        # Handle different GROMACS command formats
        if self.trjconv_cmd == "trjconv":
            # Legacy individual command
            cmd = ["trjconv"]
        else:
            # Modern gmx interface
            cmd = [self.trjconv_cmd, "trjconv"]

        # Add required parameters
        cmd.extend(["-s", topology, "-f", input_traj, "-o", output_traj, "-pbc", pbc_method, "-ur", unit_cell])

        # Add centering only when requested (Step 2)
        if center:
            cmd.append("-center")

        return cmd

    def _run_trjconv(self, cmd: List[str], input_text: str, output_file: str) -> bool:
        """Run trjconv command with error handling."""
        try:
            logger.debug(f"Running command: {' '.join(cmd)}")
            logger.debug(f"Input selections: {repr(input_text)}")
            logger.debug(f"Expected output file: {output_file}")

            # Run trjconv with input piped
            process = subprocess.run(
                cmd,
                input=input_text,
                text=True,
                capture_output=True,
                timeout=1800,  # 30 minutes timeout for large files
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
