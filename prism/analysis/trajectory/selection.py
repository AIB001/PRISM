#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Trajectory chain and center selection utilities.
"""

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class SelectionMixin:
    """Mixin for chain and center selection operations."""

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
                "A": 0,
                "B": 1,
                "C": 2,
                "D": 3,
                "E": 4,
                "F": 5,
                "P": None,  # Will try to detect
                "Q": None,  # Will try to detect
            }

            if chain_letter in chain_map:
                chainid = chain_map[chain_letter]
                if chainid is not None:
                    # Verify chain exists
                    test_atoms = traj.topology.select(f"protein and chainid {chainid}")
                    if len(test_atoms) > 0:
                        return chainid

            # For special cases like P, Q, try to find by context
            if chain_letter in ["P", "Q"]:
                return self._find_special_chain(traj, chain_letter)

        except Exception as e:
            logger.warning(f"Chain letter mapping failed: {e}")

        return None

    def _find_special_chain(self, traj, chain_letter: str):
        """Find special chains like P (primary) or Q by context."""
        try:
            logger.debug("Detecting special chain '%s' by ligand proximity", chain_letter)

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
            "chain",  # "chain P"
            "chainid",  # "chainid 0"
            "segid",  # "segid PR1"
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
            pdb_files = list(topo_path.parent.glob("*.pdb"))
            if pdb_files:
                traj = md.load(str(pdb_files[0]))
                ligand = identify_ligand_residue(traj)
                if ligand:
                    logger.info(f"Detected ligand: {ligand.name}")
                    return ligand.name

        except Exception as e:
            logger.debug(f"Could not detect ligand: {e}")

        return None

    def _prepare_selections(self, center_selection: str, output_selection: str) -> str:
        """
        Prepare input text for interactive selections.

        GROMACS will show available groups and we select by group name.
        """
        # For ligand names, use the name directly (GROMACS will find the group)
        if center_selection in ["LIG", "MOL", "UNL", "DRG", "INH", "SUB", "HET"]:
            center_input = center_selection
            logger.info(f"Using ligand group '{center_selection}' for centering")
        elif center_selection == "protein":
            center_input = "Protein"  # Standard GROMACS protein group name
        else:
            center_input = center_selection

        # Output selection - always use System
        output_input = "System" if output_selection.lower() == "system" else "0"

        return f"{center_input}\n{output_input}\n"
