#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MDTraj-based topology analysis utilities for molecular systems.

Provides chain detection and selection utilities using MDTraj for protein,
nucleic acid, and other molecular components. This is the MDTraj-compatible
version of the MDAnalysis topology utilities.
"""

import logging
from typing import Dict, Set, List, Optional, Union, Tuple
import numpy as np
import mdtraj as md

logger = logging.getLogger(__name__)


def detect_protein_chains_mdtraj(traj: md.Trajectory) -> Dict[str, str]:
    """
    Auto-detect protein chains in MDTraj trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        MDTraj trajectory object to analyze

    Returns
    -------
    Dict[str, str]
        Dictionary mapping chain names to CA selection strings
        e.g., {"Chain PR1": "protein and chainid PR1 and name CA"}
    """
    try:
        protein_chains = {}

        # Get all protein residues (standard amino acids)
        protein_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                          'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}

        # Find protein atoms
        protein_mask = np.array([residue.name in protein_residues for residue in traj.topology.residues])
        protein_residue_indices = np.where(protein_mask)[0]

        if len(protein_residue_indices) == 0:
            logger.warning("No protein residues found in trajectory")
            return protein_chains

        # Extract unique chain identifiers from protein residues
        unique_chainids = set()
        for res_idx in protein_residue_indices:
            residue = traj.topology.residue(res_idx)
            chain_id = residue.chain.index
            unique_chainids.add(chain_id)

        # Create CA selections for each detected protein chain
        for chain_id in sorted(unique_chainids):
            # Get CA atoms for this chain
            chain_ca_atoms = []
            for atom in traj.topology.atoms:
                if (atom.residue.name in protein_residues and
                    atom.residue.chain.index == chain_id and
                    atom.name == 'CA'):
                    chain_ca_atoms.append(atom.index)

            if len(chain_ca_atoms) > 0:
                chain_name = f"Chain {chain_id}"
                # Create selection string for this chain's CA atoms
                selection = f"chainid {chain_id} and name CA"
                protein_chains[chain_name] = selection
                logger.info(f"Detected protein {chain_name}: {len(chain_ca_atoms)} CA atoms")

        if len(protein_chains) == 0:
            logger.info("No distinct protein chains detected")

        return protein_chains

    except Exception as e:
        logger.warning(f"Error detecting protein chains: {e}")
        return {}


def detect_nucleic_chains_mdtraj(traj: md.Trajectory) -> Dict[str, str]:
    """
    Auto-detect nucleic acid chains in MDTraj trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        MDTraj trajectory object to analyze

    Returns
    -------
    Dict[str, str]
        Dictionary mapping chain names to P selection strings
        e.g., {"RNA Chain R": "nucleic and chainid R and name P"}
    """
    try:
        nucleic_chains = {}

        # Common nucleic acid residue names
        nucleic_residues = {'ADE', 'GUA', 'CYT', 'URA', 'THY', 'A', 'G', 'C', 'U', 'T',
                           'DA', 'DG', 'DC', 'DT', 'rA', 'rG', 'rC', 'rU'}

        # Find nucleic acid atoms
        nucleic_mask = np.array([residue.name in nucleic_residues for residue in traj.topology.residues])
        nucleic_residue_indices = np.where(nucleic_mask)[0]

        if len(nucleic_residue_indices) == 0:
            logger.info("No nucleic acid residues found in trajectory")
            return nucleic_chains

        # Extract unique chain identifiers from nucleic acid residues
        unique_chainids = set()
        for res_idx in nucleic_residue_indices:
            residue = traj.topology.residue(res_idx)
            chain_id = residue.chain.index
            unique_chainids.add(chain_id)

        # Create P selections for each detected nucleic chain
        for chain_id in sorted(unique_chainids):
            # Get P atoms for this chain
            chain_p_atoms = []
            for atom in traj.topology.atoms:
                if (atom.residue.name in nucleic_residues and
                    atom.residue.chain.index == chain_id and
                    atom.name == 'P'):
                    chain_p_atoms.append(atom.index)

            if len(chain_p_atoms) > 0:
                chain_name = f"Nucleic Chain {chain_id}"
                selection = f"chainid {chain_id} and name P"
                nucleic_chains[chain_name] = selection
                logger.info(f"Detected nucleic {chain_name}: {len(chain_p_atoms)} P atoms")
            else:
                # Try with "O5'" or "C5'" atoms if no P atoms found
                chain_o5_atoms = []
                for atom in traj.topology.atoms:
                    if (atom.residue.name in nucleic_residues and
                        atom.residue.chain.index == chain_id and
                        atom.name in ["O5'", "C5'"]):
                        chain_o5_atoms.append(atom.index)

                if len(chain_o5_atoms) > 0:
                    chain_name = f"Nucleic Chain {chain_id}"
                    selection = f"chainid {chain_id} and (name O5' or name C5')"
                    nucleic_chains[chain_name] = selection
                    logger.info(f"Detected nucleic {chain_name}: {len(chain_o5_atoms)} backbone atoms")

        return nucleic_chains

    except Exception as e:
        logger.warning(f"Error detecting nucleic acid chains: {e}")
        return {}


def detect_all_chains_mdtraj(traj: md.Trajectory) -> Dict[str, Dict[str, str]]:
    """
    Detect all molecular chains (protein and nucleic acid) in MDTraj trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        MDTraj trajectory object to analyze

    Returns
    -------
    Dict[str, Dict[str, str]]
        Nested dictionary with chain types and selections
        e.g., {
            "protein": {"Chain PR1": "chainid PR1 and name CA", ...},
            "nucleic": {"Nucleic Chain R": "chainid R and name P", ...}
        }
    """
    all_chains = {
        "protein": detect_protein_chains_mdtraj(traj),
        "nucleic": detect_nucleic_chains_mdtraj(traj)
    }

    total_protein = len(all_chains["protein"])
    total_nucleic = len(all_chains["nucleic"])
    logger.info(f"Total chains detected: {total_protein} protein, {total_nucleic} nucleic acid")

    return all_chains


def get_primary_protein_chain_mdtraj(traj: md.Trajectory) -> Optional[str]:
    """
    Get primary (largest) protein chain for alignment purposes.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        MDTraj trajectory object to analyze

    Returns
    -------
    Optional[str]
        Selection string for primary protein chain CA atoms, or None if no protein found
    """
    protein_chains = detect_protein_chains_mdtraj(traj)

    if not protein_chains:
        logger.warning("No protein chains found for primary chain selection")
        return None

    # Find chain with most CA atoms
    max_atoms = 0
    primary_chain = None
    primary_selection = None

    for chain_name, selection in protein_chains.items():
        # Count atoms in this selection
        atom_indices = traj.topology.select(selection)
        n_atoms = len(atom_indices)

        if n_atoms > max_atoms:
            max_atoms = n_atoms
            primary_chain = chain_name
            primary_selection = selection

    if primary_selection:
        logger.info(f"Primary protein chain: {primary_chain} ({max_atoms} CA atoms)")

    return primary_selection


def create_chain_selections_mdtraj(traj: md.Trajectory,
                                 include_protein: bool = True,
                                 include_nucleic: bool = True) -> Dict[str, str]:
    """
    Create selection strings for all detected chains.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        MDTraj trajectory object to analyze
    include_protein : bool
        Whether to include protein chains
    include_nucleic : bool
        Whether to include nucleic acid chains

    Returns
    -------
    Dict[str, str]
        Dictionary mapping chain names to selection strings
    """
    all_selections = {}

    if include_protein:
        protein_chains = detect_protein_chains_mdtraj(traj)
        all_selections.update(protein_chains)

    if include_nucleic:
        nucleic_chains = detect_nucleic_chains_mdtraj(traj)
        all_selections.update(nucleic_chains)

    logger.info(f"Created {len(all_selections)} chain selections")
    return all_selections


def get_chain_atom_indices_mdtraj(traj: md.Trajectory,
                                chain_selection: str) -> np.ndarray:
    """
    Get atom indices for a chain selection string.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        MDTraj trajectory object
    chain_selection : str
        MDTraj selection string (e.g., "chainid 4 and name CA")

    Returns
    -------
    np.ndarray
        Array of atom indices matching the selection
    """
    try:
        atom_indices = traj.topology.select(chain_selection)
        logger.debug(f"Selection '{chain_selection}' returned {len(atom_indices)} atoms")
        return atom_indices
    except Exception as e:
        logger.warning(f"Error selecting atoms with '{chain_selection}': {e}")
        return np.array([])


def validate_chain_selections_mdtraj(traj: md.Trajectory,
                                   chain_selections: Dict[str, str]) -> Dict[str, bool]:
    """
    Validate that chain selections return atoms.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        MDTraj trajectory object to validate against
    chain_selections : Dict[str, str]
        Dictionary of chain name to selection string

    Returns
    -------
    Dict[str, bool]
        Dictionary mapping chain names to validation results
    """
    validation_results = {}

    for chain_name, selection in chain_selections.items():
        atom_indices = get_chain_atom_indices_mdtraj(traj, selection)
        is_valid = len(atom_indices) > 0
        validation_results[chain_name] = is_valid

        if is_valid:
            logger.info(f"✓ Chain {chain_name}: {len(atom_indices)} atoms")
        else:
            logger.warning(f"✗ Chain {chain_name}: No atoms found with selection '{selection}'")

    return validation_results


def get_chain_statistics_mdtraj(traj: md.Trajectory) -> Dict[str, Dict[str, int]]:
    """
    Get statistics for all detected chains.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        MDTraj trajectory object to analyze

    Returns
    -------
    Dict[str, Dict[str, int]]
        Statistics for each chain type
        e.g., {"protein": {"n_chains": 2, "n_atoms": 1250}, ...}
    """
    stats = {}

    all_chains = detect_all_chains_mdtraj(traj)

    for chain_type, chains in all_chains.items():
        n_chains = len(chains)
        total_atoms = 0

        for chain_name, selection in chains.items():
            atom_indices = traj.topology.select(selection)
            total_atoms += len(atom_indices)

        stats[chain_type] = {
            "n_chains": n_chains,
            "n_atoms": total_atoms
        }

    return stats