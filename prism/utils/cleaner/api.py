#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Public API functions for protein cleaning.
"""

import os
from typing import Optional
from .core import ProteinCleaner


def clean_protein_pdb(
    input_pdb: str,
    output_pdb: str,
    ion_mode: str = "smart",
    distance_cutoff: float = 5.0,
    drop_nterm_met: Optional[str] = None,
    forcefield_name: Optional[str] = None,
    **kwargs,
) -> dict:
    """
    Convenience function to clean a protein PDB file.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output cleaned PDB file
    ion_mode : str
        Ion handling mode ('keep_all', 'smart', 'remove_all')
    distance_cutoff : float
        Maximum distance for keeping metals (Angstroms)
    **kwargs
        Additional arguments passed to ProteinCleaner

    Returns
    -------
    dict
        Cleaning statistics

    Examples
    --------
    >>> # Smart mode (default): keep structural metals, remove Na/Cl
    >>> stats = clean_protein_pdb('input.pdb', 'output.pdb')

    >>> # Keep all metals/ions
    >>> stats = clean_protein_pdb('input.pdb', 'output.pdb', ion_mode='keep_all')

    >>> # Remove all metals/ions
    >>> stats = clean_protein_pdb('input.pdb', 'output.pdb', ion_mode='remove_all')

    >>> # Custom distance cutoff
    >>> stats = clean_protein_pdb('input.pdb', 'output.pdb', distance_cutoff=8.0)
    """
    cleaner = ProteinCleaner(
        ion_mode=ion_mode,
        distance_cutoff=distance_cutoff,
        drop_nterm_met=drop_nterm_met,
        forcefield_name=forcefield_name,
        **kwargs,
    )
    return cleaner.clean_pdb(input_pdb, output_pdb)


def fix_terminal_atoms(pdb_file: str, output_file: str = None, force_field: str = None, verbose: bool = True) -> str:
    """
    Fix terminal atom naming and ordering for GROMACS compatibility.

    Different AMBER force fields use different C-terminal oxygen naming:
    - Standard AMBER (amber99sb, amber14sb, etc.): Single 'O' atom
    - Modified AMBER (amber99sb-star-ildn-*-mut, etc.): 'OC1' and 'OC2' atoms (carboxylate)

    This function handles both cases:
    1. Identifies C-terminal residues
    2. Removes duplicate backbone 'O'/'OXT' atoms
    3. Renames/creates proper terminal atoms based on force field:
       - Standard: Single 'O' atom after 'C'
       - Modified (mut): Two atoms 'OC1' and 'OC2' after 'C'

    Parameters
    ----------
    pdb_file : str
        Path to input PDB file
    output_file : str, optional
        Path to output file. If None, overwrites input file.
    force_field : str, optional
        Force field name to determine terminal atom naming.
        If contains 'mut' or 'star', uses OC1/OC2 naming.
        Otherwise uses single O atom (default).
    verbose : bool
        Print information about fixes

    Returns
    -------
    str
        Path to output file

    Examples
    --------
    >>> fix_terminal_atoms('protein.pdb', 'protein_fixed.pdb')
    >>> fix_terminal_atoms('protein.pdb', force_field='amber99sb-star-ildn-bsc1-mut')
    """
    if output_file is None:
        output_file = pdb_file

    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    # Determine terminal oxygen naming based on force field
    # Modified AMBER force fields (with 'mut' or 'star') use OC1/OC2
    # Standard AMBER uses single O atom
    use_oc1_oc2 = False
    if force_field:
        ff_lower = force_field.lower()
        if "mut" in ff_lower or "star" in ff_lower:
            use_oc1_oc2 = True
            if verbose:
                print(f"Detected modified AMBER force field ({force_field})")
                print("  Using OC1/OC2 terminal oxygen naming")

    # Backbone atom names that should appear before O in proper order
    BACKBONE_ATOMS = {"N", "CA", "C"}
    # Side chain atoms that indicate O is misplaced if O appears after them
    SIDECHAIN_INDICATORS = {
        "CB",
        "CG",
        "CG1",
        "CG2",
        "CD",
        "CD1",
        "CD2",
        "CE",
        "CE1",
        "CE2",
        "CE3",
        "CZ",
        "CZ2",
        "CZ3",
        "CH2",
        "NE",
        "NE1",
        "NE2",
        "ND1",
        "ND2",
        "NZ",
        "OD1",
        "OD2",
        "OE1",
        "OE2",
        "OG",
        "OG1",
        "OH",
        "SD",
        "SG",
    }

    all_lines = []
    with open(pdb_file, "r") as f:
        all_lines = f.readlines()

    # First pass: group atoms by residue
    residues = {}  # res_id -> list of (line, index, atom_name)

    for i, line in enumerate(all_lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            chain = line[21]
            res_id = f"{chain}:{res_name}:{res_num}"

            if res_id not in residues:
                residues[res_id] = []
            residues[res_id].append({"line": line, "index": i, "atom": atom_name, "res_name": res_name})

    # Second pass: identify C-terminal residues needing fixes
    c_terminal_residues = {}  # res_id -> {c_index, skip_indices, o_line}

    for res_id, atoms in residues.items():
        c_idx = None
        o_idx = None
        o_line = None
        has_sidechain = False
        oxt_idx = None

        # Find positions of key atoms
        for i, atom_info in enumerate(atoms):
            atom_name = atom_info["atom"]

            if atom_name == "C":
                c_idx = i
                c_line = atom_info["line"]  # Store C line for coordinate reference
            elif atom_name == "O":
                o_idx = i
                o_line = atom_info["line"]
            elif atom_name == "OXT":
                oxt_idx = i
                oxt_line = atom_info["line"]
            elif atom_name in SIDECHAIN_INDICATORS:
                has_sidechain = True

        # Determine if this is a C-terminal needing fix
        needs_fix = False
        skip_atom_indices = []  # Track all atoms to skip (may be multiple O/OXT)
        terminal_lines = []  # Lines to insert (may be 1 or 2 oxygen atoms)

        # Case 1: OXT exists (clear C-terminal)
        if oxt_idx is not None:
            needs_fix = True
            skip_atom_indices.append(oxt_idx)
            # CRITICAL FIX: If both O and OXT exist (PDBFixer added OXT), skip BOTH
            if o_idx is not None and o_idx != oxt_idx:
                skip_atom_indices.append(o_idx)

            # Generate terminal oxygen atom(s) based on force field
            source_line = atoms[oxt_idx]["line"]
            if use_oc1_oc2:
                # Modified AMBER: Create OC1 and OC2 (two terminal oxygens)
                oc1_line = source_line[:12] + " OC1" + source_line[16:]
                oc2_line = source_line[:12] + " OC2" + source_line[16:]
                terminal_lines = [oc1_line, oc2_line]
            else:
                # Standard AMBER: Single O atom
                o_line = source_line[:12] + " O  " + source_line[16:]
                terminal_lines = [o_line]

        # Case 2: O appears after side chain atoms (misplaced terminal O)
        elif o_idx is not None and c_idx is not None and has_sidechain:
            # Check if O appears after any side chain atom
            for i, atom_info in enumerate(atoms):
                if atom_info["atom"] in SIDECHAIN_INDICATORS and i < o_idx:
                    needs_fix = True
                    skip_atom_indices.append(o_idx)

                    # Generate terminal oxygen atom(s) based on force field
                    source_line = atoms[o_idx]["line"]
                    if use_oc1_oc2:
                        # Modified AMBER: Create OC1 and OC2
                        oc1_line = source_line[:12] + " OC1" + source_line[16:]
                        oc2_line = source_line[:12] + " OC2" + source_line[16:]
                        terminal_lines = [oc1_line, oc2_line]
                    else:
                        # Standard AMBER: Keep as O
                        terminal_lines = [source_line]
                    break

        if needs_fix and c_idx is not None and terminal_lines:
            c_terminal_residues[res_id] = {
                "c_index": atoms[c_idx]["index"],
                "skip_indices": [atoms[idx]["index"] for idx in skip_atom_indices],
                "terminal_lines": terminal_lines,  # Changed from 'o_line' to support multiple lines
                "res_name": atoms[0]["res_name"],
                "chain": res_id.split(":")[0],
                "res_num": res_id.split(":")[-1],
                "is_oxt": oxt_idx is not None,
                "use_oc1_oc2": use_oc1_oc2,
            }

    # Third pass: reconstruct file with corrected terminal residues
    skip_indices = set()
    insert_map = {}  # index -> lines_to_insert (may be multiple lines)
    fixed_count = 0

    for res_id, term_info in c_terminal_residues.items():
        # Skip ALL old O/OXT positions (may be multiple when both O and OXT exist)
        for skip_idx in term_info["skip_indices"]:
            skip_indices.add(skip_idx)
        # Insert corrected terminal oxygen(s) right after C
        insert_map[term_info["c_index"]] = term_info["terminal_lines"]
        fixed_count += 1

        if verbose:
            fix_type = "OXT → O" if term_info["is_oxt"] else "misplaced O"
            skip_count = len(term_info["skip_indices"])
            atoms_removed = f"{skip_count} atom(s)" if skip_count > 1 else "1 atom"

            if term_info["use_oc1_oc2"]:
                atoms_added = "OC1 and OC2"
            else:
                atoms_added = "O"

            print(
                f"  Fixed {fix_type} in {term_info['res_name']} {term_info['chain']}{term_info['res_num']} "
                f"(removed {atoms_removed}, added {atoms_added})"
            )

    # Write output with insertions
    with open(output_file, "w") as f:
        for i, line in enumerate(all_lines):
            if i not in skip_indices:
                f.write(line)
                # Insert terminal oxygen atom(s) after C atom if needed
                if i in insert_map:
                    for terminal_line in insert_map[i]:
                        f.write(terminal_line)

    if verbose:
        if fixed_count > 0:
            print(f"\nFixed {fixed_count} C-terminal residues with misplaced oxygen atoms")
        else:
            print("No terminal atom issues found")

    return output_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Clean protein PDB files with intelligent metal ion handling")
    parser.add_argument("input", help="Input PDB file")
    parser.add_argument("output", help="Output PDB file")
    parser.add_argument(
        "--mode",
        choices=["keep_all", "smart", "remove_all"],
        default="smart",
        help="Ion handling mode (default: smart)",
    )
    parser.add_argument(
        "--distance", type=float, default=5.0, help="Distance cutoff for keeping metals in Angstroms (default: 5.0)"
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress output messages")
    parser.add_argument("--fix-termini", action="store_true", help="Fix terminal atom names (OXT → O) for GROMACS")

    args = parser.parse_args()

    clean_protein_pdb(
        args.input, args.output, ion_mode=args.mode, distance_cutoff=args.distance, verbose=not args.quiet
    )

    if args.fix_termini:
        fix_terminal_atoms(args.output, verbose=not args.quiet)
