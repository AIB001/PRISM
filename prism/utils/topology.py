#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Topology analysis utilities for molecular systems.

Provides chain detection and selection utilities for protein, nucleic acid,
and other molecular components.
"""

import logging
from typing import Dict, Set, List, Optional, Union, Tuple
import MDAnalysis as mda

logger = logging.getLogger(__name__)


def detect_protein_chains(universe: mda.Universe) -> Dict[str, str]:
    """
    Auto-detect protein chains in the universe.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe object to analyze

    Returns
    -------
    Dict[str, str]
        Dictionary mapping chain names to CA selection strings
        e.g., {"Chain PR1": "protein and chainid PR1 and name CA"}
    """
    try:
        protein_chains = {}

        # Get all unique chain/segment IDs for protein residues
        protein_atoms = universe.select_atoms("protein")
        if len(protein_atoms) == 0:
            logger.warning("No protein atoms found in universe")
            return protein_chains

        # Extract unique chain identifiers from protein atoms
        unique_chainids = set()
        for atom in protein_atoms:
            try:
                # Try different chain identifier attributes
                if hasattr(atom, 'chainid') and atom.chainid is not None:
                    unique_chainids.add(str(atom.chainid))
                elif hasattr(atom, 'segid') and atom.segid is not None and atom.segid != '':
                    unique_chainids.add(str(atom.segid))
                elif hasattr(atom.residue, 'chainid') and atom.residue.chainid is not None:
                    unique_chainids.add(str(atom.residue.chainid))
            except (AttributeError, TypeError):
                continue

        # Create CA selections for each detected protein chain
        for chainid in sorted(unique_chainids):
            if chainid and chainid != 'None':
                # Try chainid first, then segid as fallback
                selection_chainid = f"protein and chainid {chainid} and name CA"
                selection_segid = f"protein and segid {chainid} and name CA"

                test_atoms_chainid = universe.select_atoms(selection_chainid)
                test_atoms_segid = universe.select_atoms(selection_segid)

                if len(test_atoms_chainid) > 0:
                    chain_name = f"Chain {chainid}"
                    protein_chains[chain_name] = selection_chainid
                    logger.info(f"Detected protein {chain_name}: {len(test_atoms_chainid)} CA atoms")
                elif len(test_atoms_segid) > 0:
                    chain_name = f"Chain {chainid}"
                    protein_chains[chain_name] = selection_segid
                    logger.info(f"Detected protein {chain_name}: {len(test_atoms_segid)} CA atoms")

        if len(protein_chains) == 0:
            logger.info("No distinct protein chains detected")

        return protein_chains

    except Exception as e:
        logger.warning(f"Error detecting protein chains: {e}")
        return {}


def detect_nucleic_chains(universe: mda.Universe) -> Dict[str, str]:
    """
    Auto-detect nucleic acid chains in the universe.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe object to analyze

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

        # Create fallback selection string
        resname_selection = " or ".join([f"resname {res}" for res in nucleic_residues])

        # Try to select nucleic acids using nucleic keyword first
        try:
            nucleic_atoms = universe.select_atoms("nucleic")
        except Exception:
            # Fallback: select by residue name
            try:
                nucleic_atoms = universe.select_atoms(resname_selection)
            except Exception:
                logger.info("No nucleic acid atoms found")
                return nucleic_chains

        if len(nucleic_atoms) == 0:
            logger.info("No nucleic acid atoms found")
            return nucleic_chains

        # Extract unique chain identifiers from nucleic atoms
        unique_chainids = set()
        for atom in nucleic_atoms:
            try:
                if hasattr(atom, 'chainid') and atom.chainid is not None:
                    unique_chainids.add(str(atom.chainid))
                elif hasattr(atom, 'segid') and atom.segid is not None and atom.segid != '':
                    unique_chainids.add(str(atom.segid))
                elif hasattr(atom.residue, 'chainid') and atom.residue.chainid is not None:
                    unique_chainids.add(str(atom.residue.chainid))
            except (AttributeError, TypeError):
                continue

        # Create P selections for each detected nucleic chain
        for chainid in sorted(unique_chainids):
            if chainid and chainid != 'None':
                # Try different selection strategies
                selections_to_try = [
                    f"nucleic and chainid {chainid} and name P",
                    f"nucleic and segid {chainid} and name P",
                    f"({resname_selection}) and chainid {chainid} and name P",
                    f"({resname_selection}) and segid {chainid} and name P"
                ]

                for selection in selections_to_try:
                    try:
                        test_atoms = universe.select_atoms(selection)
                        if len(test_atoms) > 0:
                            chain_name = f"Nucleic Chain {chainid}"
                            nucleic_chains[chain_name] = selection
                            logger.info(f"Detected nucleic {chain_name}: {len(test_atoms)} P atoms")
                            break
                    except Exception:
                        continue

        return nucleic_chains

    except Exception as e:
        logger.warning(f"Error detecting nucleic acid chains: {e}")
        return {}


def detect_all_chains(universe: mda.Universe) -> Dict[str, Dict[str, str]]:
    """
    Detect all molecular chains (protein and nucleic acid) in the universe.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe object to analyze

    Returns
    -------
    Dict[str, Dict[str, str]]
        Nested dictionary with chain types and selections
        e.g., {
            "protein": {"Chain PR1": "protein and segid PR1 and name CA", ...},
            "nucleic": {"Nucleic Chain R": "nucleic and chainid R and name P", ...}
        }
    """
    all_chains = {
        "protein": detect_protein_chains(universe),
        "nucleic": detect_nucleic_chains(universe)
    }

    total_protein = len(all_chains["protein"])
    total_nucleic = len(all_chains["nucleic"])
    logger.info(f"Total chains detected: {total_protein} protein, {total_nucleic} nucleic acid")

    return all_chains


def get_primary_protein_chain(universe: mda.Universe) -> Optional[str]:
    """
    Get the primary (largest) protein chain for alignment purposes.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe object to analyze

    Returns
    -------
    Optional[str]
        Selection string for primary protein chain CA atoms, or None if no protein found
    """
    protein_chains = detect_protein_chains(universe)

    if not protein_chains:
        logger.warning("No protein chains found for primary chain selection")
        return None

    # Find the chain with most CA atoms
    max_atoms = 0
    primary_chain = None
    primary_selection = None

    for chain_name, selection in protein_chains.items():
        atoms = universe.select_atoms(selection)
        if len(atoms) > max_atoms:
            max_atoms = len(atoms)
            primary_chain = chain_name
            primary_selection = selection

    if primary_selection:
        logger.info(f"Primary protein chain: {primary_chain} ({max_atoms} CA atoms)")

    return primary_selection