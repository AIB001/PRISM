#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PDB file utility functions.

This module provides utility functions for manipulating PDB files,
including residue number renumbering and other common operations.
"""

import logging
from pathlib import Path
from typing import Optional, Dict

logger = logging.getLogger(__name__)


def renumber_residues(
    input_pdb: str,
    output_pdb: str,
    offset: int = 0,
    chain_map: Optional[Dict[str, int]] = None,
    start_at: Optional[int] = None,
    keep_ter: bool = True,
) -> Dict[str, any]:
    """
    Renumber residues in a PDB file.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output PDB file with renumbered residues
    offset : int
        Offset to add to residue numbers (default: 0)
    chain_map : dict, optional
        Dictionary mapping chain IDs to offsets, e.g., {'A': 0, 'B': 100}
        If provided, this takes precedence over the offset parameter.
    start_at : int, optional
        Renumber residues starting from this number (1, 2, 3, ...)
        Overrides offset if specified.
    keep_ter : bool
        Keep TER cards in output (default: True)

    Returns
    -------
    dict
        Dictionary with statistics about the renumbering

    Examples
    --------
    >>> # Add offset of 100 to all residues
    >>> renumber_residues("input.pdb", "output.pdb", offset=100)

    >>> # Renumber starting from 1
    >>> renumber_residues("input.pdb", "output.pdb", start_at=1)

    >>> # Different offsets per chain
    >>> renumber_residues("input.pdb", "output.pdb", chain_map={'A': 0, 'B': 500})
    """
    input_path = Path(input_pdb)
    output_path = Path(output_pdb)

    if not input_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

    if offset == 0 and start_at is None and chain_map is None:
        logger.warning("No renumbering specified (offset=0, no start_at, no chain_map)")

    stats = {"atoms_processed": 0, "residues_renumbered": set(), "chains_seen": set(), "num_ter": 0}

    # For start_at mode: track current residue number per chain
    current_nums = {}

    with open(input_path, "r") as f_in, open(output_path, "w") as f_out:
        for line in f_in:
            # Handle TER records
            if line.startswith("TER"):
                stats["num_ter"] += 1
                if keep_ter:
                    f_out.write(line)
                continue

            # Only process ATOM/HETATM lines
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                f_out.write(line)
                continue

            stats["atoms_processed"] += 1

            # Parse residue info
            # COLUMNS        DATA TYPE       FIELD       DEFINITION
            # ----------------------------------------------------
            # 1 - 6         Record name    "ATOM  "
            # 7 - 11        Integer       serial       Atom serial number
            # 13 - 16        Atom          name         Atom name
            # 17             Character      altLoc       Alternate location indicator
            # 18 - 20        Residue name  resName      Residue name
            # 22             Character      chainID      Chain identifier
            # 23 - 26        Integer       resSeq       Residue sequence number
            # 27             Character      iCode        Code for insertion of residues
            # ...
            res_name = line[17:20].strip()
            chain = line[21].strip()
            res_num_str = line[22:26].strip()
            ins_code = line[26] if len(line) > 26 else " "
            rest_of_line = line[27:] if len(line) > 27 else ""

            stats["chains_seen"].add(chain)

            try:
                old_num = int(res_num_str)
                res_key = f"{chain}:{old_num}:{ins_code}"

                # Determine new residue number
                if start_at is not None:
                    # Renumber sequentially
                    if chain not in current_nums:
                        current_nums[chain] = start_at
                    new_num = current_nums[chain]
                    # Check if we moved to a new residue
                    if res_key not in stats["residues_renumbered"]:
                        stats["residues_renumbered"].add(res_key)
                        current_nums[chain] += 1
                elif chain_map and chain in chain_map:
                    # Use chain-specific offset
                    new_num = old_num + chain_map[chain]
                    stats["residues_renumbered"].add(res_key)
                elif offset != 0:
                    # Use global offset
                    new_num = old_num + offset
                    stats["residues_renumbered"].add(res_key)
                else:
                    new_num = old_num

                # Write modified line
                new_line = (
                    line[:17]  # Through residue name
                    + f"{res_name:3s}"  # Residue name (right-aligned)
                    + chain  # Chain ID
                    + f"{new_num:>4}"  # New residue number (right-aligned)
                    + ins_code  # Insertion code
                    + rest_of_line
                )  # Rest of line
                f_out.write(new_line)

            except ValueError:
                # Couldn't parse residue number, keep original
                f_out.write(line)

    stats["residues_renumbered"] = len(stats["residues_renumbered"])
    stats["chains_seen"] = list(stats["chains_seen"])

    logger.info(f"Renumbered {stats['residues_renumbered']} residues " f"in {len(stats['chains_seen'])} chains")
    logger.info(f"Output: {output_path}")

    return stats


def apply_pka_predictions(
    input_pdb: str,
    output_pdb: str,
    predictions: Dict[str, Dict],
    ph: float = 7.0,
    ff_info: Optional[Dict[str, any]] = None,
) -> Dict[str, any]:
    """
    Apply PROPKA3 pKa predictions to a PDB file.

    This function renames residues based on their predicted protonation
    states at the given pH.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output PDB file with protonation states applied
    predictions : dict
        Dictionary from PROPKA3 with pKa predictions
        Format: {f"{chain}:{resnum}": {"resname": "HIS", "pka": 6.5, ...}}
    ph : float
        Target pH (default: 7.0)

    Returns
    -------
    dict
        Dictionary with renames applied and statistics
    """
    from prism.utils.protonation import PropkaProtonator, resolve_protonation_resname

    input_path = Path(input_pdb)
    output_path = Path(output_pdb)

    protonator = PropkaProtonator(ph=ph)

    # Determine protonation states
    renames = {}
    for key, pred in predictions.items():
        res_name = pred["resname"]
        pka = pred["pka"]

        if pka is None:
            continue

        desired, _ = protonator.determine_protonation_state(res_name, pka)
        state, _note = resolve_protonation_resname(res_name, desired, ff_info)

        if res_name != state and state != res_name + "H":
            renames[key] = {"from": res_name, "to": state}

    # Apply renames
    with open(input_path, "r") as f_in, open(output_path, "w") as f_out:
        for line in f_in:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                f_out.write(line)
                continue

            res_name = line[17:20].strip()
            chain = line[21].strip()
            res_num = line[22:26].strip()
            key = f"{chain}:{res_num}"

            if key in renames:
                new_res_name = renames[key]["to"]
                # Replace residue name (columns 17-19 in 0-indexed)
                new_line = line[:17] + f"{new_res_name:3s}" + line[20:]
                f_out.write(new_line)
            else:
                f_out.write(line)

    logger.info(f"Applied {len(renames)} protonation state renames")

    return {"renames": renames, "num_renamed": len(renames), "input": str(input_path), "output": str(output_path)}


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python pdb_utils.py <input.pdb> <output.pdb> [offset]")
        print("Examples:")
        print("  python pdb_utils.py input.pdb output.pdb 100     # Add offset of 100")
        print("  python pdb_utils.py input.pdb output.pdb 0        # No renumbering")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    offset = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    stats = renumber_residues(input_pdb, output_pdb, offset=offset)

    print(f"\nResidue renumbering complete:")
    print(f"  Atoms processed: {stats['atoms_processed']}")
    print(f"  Residues renumbered: {stats['residues_renumbered']}")
    print(f"  Chains seen: {stats['chains_seen']}")
    print(f"  TER cards: {stats['num_ter']}")
