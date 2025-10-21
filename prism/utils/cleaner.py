#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Intelligent Protein Structure Cleaner for PRISM

This module provides smart cleaning of protein PDB files with intelligent
handling of metal ions and other heteroatoms. It can distinguish between
structural metals (e.g., Zn, Mg, Ca) that are essential for protein function
and non-structural ions (e.g., Na, Cl) that can be safely removed.

Features:
- Three ion handling modes: 'keep_all', 'smart', 'remove_all'
- Distance-based filtering to remove distant metals
- Preservation of important structural metals
- Removal of water molecules and non-essential ions
"""

import os
import numpy as np
from typing import List, Tuple, Set, Optional


class ProteinCleaner:
    """
    Intelligent protein structure cleaner with smart metal ion handling.

    This class processes PDB files to remove unwanted heteroatoms while
    intelligently preserving structural metals and other important cofactors.
    """

    # Structural metals that are typically important for protein function
    STRUCTURAL_METALS = {
        'ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'CO', 'NI', 'CD',
        'HG', 'ZN2', 'MG2', 'CA2', 'FE2', 'FE3', 'CU2', 'MN2'
    }

    # Non-structural ions that can usually be safely removed
    NON_STRUCTURAL_IONS = {
        'NA', 'CL', 'K', 'BR', 'I', 'F',
        'NA+', 'CL-', 'K+', 'BR-', 'I-',
        'SO4', 'PO4', 'NO3', 'CO3',
        'SUL', 'PHO', 'NIT', 'CAR'
    }

    # Crystallization artifacts and common additives (removed in smart mode)
    CRYSTALLIZATION_ARTIFACTS = {
        # Glycerol and polyols
        'GOL', 'EDO', 'MPD', 'MRD', 'PGO', 'PG4', 'BTB',
        # PEG oligomers
        'PEG', 'PGE', '1PE', 'P6G', 'P33', 'PE8', '2PE', 'XPE', '12P', '33O', 'P3E',
        # Sugars (unless covalently linked to protein)
        'NAG', 'NDG', 'BMA', 'MAN', 'GAL', 'FUC', 'GLC', 'SIA',
        # Detergents
        'DMS', 'BOG', 'LMT', 'C8E', 'DAO', 'SDS', 'LDA',
        # Other common additives
        'ACT', 'ACE', 'FMT', 'TRS', 'EPE', 'BCT', 'BME'
    }

    # Water residue names
    WATER_NAMES = {'HOH', 'WAT', 'H2O', 'TIP', 'SPC', 'SOL'}

    def __init__(self,
                 ion_mode: str = 'smart',
                 distance_cutoff: float = 5.0,
                 keep_crystal_water: bool = False,
                 remove_artifacts: bool = True,
                 keep_custom_metals: Optional[List[str]] = None,
                 remove_custom_metals: Optional[List[str]] = None,
                 verbose: bool = True):
        """
        Initialize the ProteinCleaner.

        Parameters
        ----------
        ion_mode : str
            Ion handling mode. Options:
            - 'keep_all': Keep all metal ions and heteroatoms (except water unless keep_crystal_water=True)
            - 'smart' (default): Keep structural metals, remove non-structural ions
            - 'remove_all': Remove all metal ions and heteroatoms
        distance_cutoff : float
            Maximum distance (in Angstroms) from protein for keeping metals.
            Metals farther than this from any protein heavy atom will be removed.
            Default: 5.0 Å
        keep_crystal_water : bool
            Whether to keep crystal water molecules. Default: False
            If True, water molecules are preserved in the output
        remove_artifacts : bool
            Whether to remove crystallization artifacts (GOL, EDO, PEG, NAG, etc.). Default: True
            If False, crystallization artifacts are kept in the output
        keep_custom_metals : list of str, optional
            Additional metal/ion names to always keep (e.g., ['MO', 'W'])
        remove_custom_metals : list of str, optional
            Additional metal/ion names to always remove
        verbose : bool
            Print detailed information during cleaning
        """
        if ion_mode not in ['keep_all', 'smart', 'remove_all']:
            raise ValueError(f"Invalid ion_mode: {ion_mode}. Must be 'keep_all', 'smart', or 'remove_all'")

        self.ion_mode = ion_mode
        self.distance_cutoff = distance_cutoff
        self.keep_crystal_water = keep_crystal_water
        self.remove_artifacts = remove_artifacts
        self.verbose = verbose

        # Build custom metal sets
        self.structural_metals = self.STRUCTURAL_METALS.copy()
        if keep_custom_metals:
            self.structural_metals.update(m.upper() for m in keep_custom_metals)

        self.non_structural_ions = self.NON_STRUCTURAL_IONS.copy()
        if remove_custom_metals:
            self.non_structural_ions.update(m.upper() for m in remove_custom_metals)

        # Statistics
        self.stats = {
            'total_atoms': 0,
            'protein_atoms': 0,
            'water_removed': 0,
            'water_kept': 0,
            'ions_removed': 0,
            'metals_kept': 0,
            'metals_removed_by_distance': 0,
            'artifacts_removed': 0
        }

    def clean_pdb(self, input_pdb: str, output_pdb: str) -> dict:
        """
        Clean a PDB file with intelligent ion handling.

        Parameters
        ----------
        input_pdb : str
            Path to input PDB file
        output_pdb : str
            Path to output cleaned PDB file

        Returns
        -------
        dict
            Statistics about the cleaning process
        """
        if not os.path.exists(input_pdb):
            raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")

        if self.verbose:
            print(f"\n=== Cleaning Protein Structure ===")
            print(f"Input: {input_pdb}")
            print(f"Output: {output_pdb}")
            print(f"Mode: {self.ion_mode}")
            print(f"Distance cutoff: {self.distance_cutoff} Å")
            print(f"Keep crystal water: {self.keep_crystal_water}")
            print(f"Remove crystallization artifacts: {self.remove_artifacts}")

        # Read and parse PDB file
        protein_lines, hetatm_lines = self._read_pdb(input_pdb)

        # Extract protein coordinates for distance calculations
        protein_coords = self._extract_protein_coords(protein_lines)

        # Process HETATM records
        kept_hetatm_lines = []
        if self.ion_mode != 'remove_all':
            kept_hetatm_lines = self._process_hetatm_records(hetatm_lines, protein_coords)

        # Write cleaned PDB
        self._write_pdb(output_pdb, protein_lines, kept_hetatm_lines, input_pdb)

        # Print statistics
        if self.verbose:
            self._print_statistics()

        return self.stats.copy()

    def _read_pdb(self, pdb_file: str) -> Tuple[List[str], List[str]]:
        """
        Read PDB file and separate protein and heteroatom records.

        Returns
        -------
        tuple
            (protein_lines, hetatm_lines)
        """
        protein_lines = []
        hetatm_lines = []

        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    protein_lines.append(line)
                    self.stats['total_atoms'] += 1
                    self.stats['protein_atoms'] += 1
                elif line.startswith('HETATM'):
                    hetatm_lines.append(line)
                    self.stats['total_atoms'] += 1

        return protein_lines, hetatm_lines

    def _extract_protein_coords(self, protein_lines: List[str]) -> np.ndarray:
        """
        Extract heavy atom coordinates from protein ATOM records.

        Parameters
        ----------
        protein_lines : list of str
            ATOM records from PDB file

        Returns
        -------
        np.ndarray
            Coordinates array of shape (n_atoms, 3)
        """
        coords = []
        for line in protein_lines:
            # Skip hydrogen atoms
            atom_name = line[12:16].strip()
            if atom_name.startswith('H'):
                continue

            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue

        return np.array(coords) if coords else np.array([]).reshape(0, 3)

    def _process_hetatm_records(self, hetatm_lines: List[str],
                                protein_coords: np.ndarray) -> List[str]:
        """
        Process HETATM records according to ion handling mode.

        Parameters
        ----------
        hetatm_lines : list of str
            HETATM records from PDB file
        protein_coords : np.ndarray
            Protein heavy atom coordinates

        Returns
        -------
        list of str
            HETATM records to keep
        """
        kept_lines = []

        for line in hetatm_lines:
            residue_name = line[17:20].strip().upper()
            atom_name = line[12:16].strip().upper()

            # Handle water molecules
            if residue_name in self.WATER_NAMES:
                if self.keep_crystal_water:
                    kept_lines.append(line)
                    self.stats['water_kept'] += 1
                else:
                    self.stats['water_removed'] += 1
                continue

            # Handle crystallization artifacts
            if residue_name in self.CRYSTALLIZATION_ARTIFACTS:
                if self.remove_artifacts:
                    # Remove if flag is True (default)
                    self.stats['artifacts_removed'] += 1
                else:
                    # Keep if flag is False
                    kept_lines.append(line)
                continue

            # Check if this is a metal/ion
            is_metal = self._is_metal_or_ion(residue_name, atom_name)

            if not is_metal:
                # Not a metal/ion or artifact - keep it (might be a cofactor or modified residue)
                kept_lines.append(line)
                continue

            # Handle metal/ion based on mode
            should_keep = False

            if self.ion_mode == 'keep_all':
                should_keep = True
            elif self.ion_mode == 'smart':
                should_keep = self._should_keep_metal_smart(residue_name, atom_name)
            # elif self.ion_mode == 'remove_all': should_keep remains False

            # Distance filtering for kept metals (only if cutoff is set)
            if should_keep and self.distance_cutoff is not None and len(protein_coords) > 0:
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    metal_coord = np.array([x, y, z])

                    min_dist = self._calculate_min_distance(metal_coord, protein_coords)

                    if min_dist > self.distance_cutoff:
                        should_keep = False
                        self.stats['metals_removed_by_distance'] += 1
                        if self.verbose:
                            print(f"  Removed {residue_name} at distance {min_dist:.2f} Å from protein")
                except (ValueError, IndexError):
                    # If we can't parse coordinates, skip this atom
                    should_keep = False

            if should_keep:
                kept_lines.append(line)
                self.stats['metals_kept'] += 1
            else:
                self.stats['ions_removed'] += 1

        return kept_lines

    def _is_metal_or_ion(self, residue_name: str, atom_name: str) -> bool:
        """
        Check if a residue/atom is a metal or ion.

        Parameters
        ----------
        residue_name : str
            Residue name from PDB
        atom_name : str
            Atom name from PDB

        Returns
        -------
        bool
            True if this is a metal or ion
        """
        # Check against known lists
        if residue_name in self.structural_metals or residue_name in self.non_structural_ions:
            return True

        # Also check atom name (sometimes metal is in atom name)
        if atom_name in self.structural_metals or atom_name in self.non_structural_ions:
            return True

        # Common metal ion patterns
        metal_patterns = ['ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'CO', 'NI', 'NA', 'K', 'CL']
        for pattern in metal_patterns:
            if pattern in residue_name or pattern in atom_name:
                return True

        return False

    def _should_keep_metal_smart(self, residue_name: str, atom_name: str) -> bool:
        """
        Determine if a metal should be kept in 'smart' mode.

        Parameters
        ----------
        residue_name : str
            Residue name from PDB
        atom_name : str
            Atom name from PDB

        Returns
        -------
        bool
            True if this metal should be kept
        """
        # Check if it's a structural metal
        if residue_name in self.structural_metals or atom_name in self.structural_metals:
            return True

        # Check if it's explicitly a non-structural ion
        if residue_name in self.non_structural_ions or atom_name in self.non_structural_ions:
            return False

        # Unknown metal - be conservative and keep it
        if self.verbose:
            print(f"  Warning: Unknown metal/ion {residue_name}/{atom_name}, keeping by default")
        return True

    def _calculate_min_distance(self, point: np.ndarray, coords: np.ndarray) -> float:
        """
        Calculate minimum distance from a point to a set of coordinates.

        Parameters
        ----------
        point : np.ndarray
            Single coordinate (3,)
        coords : np.ndarray
            Array of coordinates (n, 3)

        Returns
        -------
        float
            Minimum distance in Angstroms
        """
        if len(coords) == 0:
            return float('inf')

        distances = np.sqrt(np.sum((coords - point) ** 2, axis=1))
        return float(np.min(distances))

    def _convert_metals_to_atom(self, hetatm_lines: List[str]) -> Tuple[List[str], List[str]]:
        """
        Convert metal ion HETATM records to ATOM records for pdb2gmx compatibility.

        This allows pdb2gmx to directly recognize and include metal ions in the topology,
        avoiding the need for manual addition later.

        Parameters
        ----------
        hetatm_lines : list of str
            HETATM records from cleaned PDB

        Returns
        -------
        tuple
            (metal_atom_lines, other_hetatm_lines)
            - metal_atom_lines: Metal ions converted to ATOM format
            - other_hetatm_lines: Other HETATM records (cofactors, etc.)
        """
        metal_atom_lines = []
        other_hetatm_lines = []

        for line in hetatm_lines:
            residue_name = line[17:20].strip().upper()
            atom_name = line[12:16].strip().upper()

            # Check if this is a metal ion that should be converted
            is_metal = (residue_name in self.structural_metals or
                       atom_name in self.structural_metals or
                       residue_name in self.non_structural_ions or
                       atom_name in self.non_structural_ions)

            if is_metal:
                # Convert HETATM to ATOM
                # PDB format: columns 1-6 are record type
                atom_line = 'ATOM  ' + line[6:]
                metal_atom_lines.append(atom_line)

                if self.verbose:
                    chain = line[21] if len(line) > 21 else ' '
                    resnum = line[22:26].strip() if len(line) > 26 else ''
                    print(f"  Converting {residue_name} (chain {chain}, res {resnum}) from HETATM to ATOM")
            else:
                # Keep as HETATM (cofactors, modified residues, etc.)
                other_hetatm_lines.append(line)

        return metal_atom_lines, other_hetatm_lines

    def _group_by_chain(self, protein_lines: List[str], metal_lines: List[str]) -> dict:
        """
        Group protein and metal atoms by chain ID for proper PDB output.

        IMPORTANT: pdb2gmx requires each chain to have a consistent residue type.
        Since metals are type 'Ion' and proteins are type 'Protein', they cannot
        be in the same chain. We assign each metal to a NEW unique chain ID.

        Parameters
        ----------
        protein_lines : list of str
            Protein ATOM records
        metal_lines : list of str
            Metal ATOM records (converted from HETATM)

        Returns
        -------
        dict
            {chain_id: ([protein_atoms], [metal_atoms])}
        """
        import string

        chains = {}

        # Group protein atoms by chain
        for line in protein_lines:
            if len(line) > 21:
                chain_id = line[21]
                if chain_id not in chains:
                    chains[chain_id] = ([], [])
                chains[chain_id][0].append(line)

        # Assign metals to NEW unique chain IDs (not mixed with protein)
        # This is required because pdb2gmx doesn't allow mixed Protein/Ion types
        used_chains = set(chains.keys())
        available_chains = [c for c in string.ascii_uppercase + string.digits
                           if c not in used_chains]

        metal_chain_idx = 0
        for line in metal_lines:
            if len(line) > 21:
                original_chain = line[21]

                # Assign a new unique chain ID for this metal
                if metal_chain_idx < len(available_chains):
                    new_chain_id = available_chains[metal_chain_idx]
                else:
                    # Fallback: use numbers or special characters
                    new_chain_id = str(metal_chain_idx % 10)

                metal_chain_idx += 1

                # Modify the line to use the new chain ID
                # PDB format: chain ID is at position 21
                modified_line = line[:21] + new_chain_id + line[22:]

                # Create new chain entry for this metal
                if new_chain_id not in chains:
                    chains[new_chain_id] = ([], [])
                chains[new_chain_id][1].append(modified_line)

                if self.verbose:
                    res_name = line[17:20].strip()
                    resnum = line[22:26].strip()
                    print(f"  Reassigning {res_name} {resnum} from chain {original_chain} to chain {new_chain_id} (for pdb2gmx compatibility)")

        return chains

    def _write_pdb(self, output_pdb: str, protein_lines: List[str],
                   hetatm_lines: List[str], input_pdb: str):
        """
        Write cleaned PDB file with proper chain organization for pdb2gmx.

        Parameters
        ----------
        output_pdb : str
            Output file path
        protein_lines : list of str
            Protein ATOM records
        hetatm_lines : list of str
            Kept HETATM records
        input_pdb : str
            Original input file (for header information)
        """
        # Convert metal HETATM to ATOM for pdb2gmx compatibility
        metal_atom_lines, other_hetatm_lines = self._convert_metals_to_atom(hetatm_lines)

        # Group protein and metal atoms by chain for proper output
        chain_atoms = self._group_by_chain(protein_lines, metal_atom_lines)

        with open(output_pdb, 'w') as out:
            # Write header
            out.write(f"REMARK   Cleaned by PRISM ProteinCleaner\n")
            out.write(f"REMARK   Original file: {os.path.basename(input_pdb)}\n")
            out.write(f"REMARK   Ion mode: {self.ion_mode}\n")
            out.write(f"REMARK   Distance cutoff: {self.distance_cutoff} A\n")
            out.write(f"REMARK   Keep crystal water: {self.keep_crystal_water}\n")
            out.write(f"REMARK   Remove crystallization artifacts: {self.remove_artifacts}\n")
            if metal_atom_lines:
                out.write(f"REMARK   Converted {len(metal_atom_lines)} metal HETATM to ATOM for pdb2gmx\n")
                out.write(f"REMARK   Metals integrated into their respective protein chains\n")

            # Write chains with protein atoms followed by their metals
            for chain_id in sorted(chain_atoms.keys()):
                protein_atoms, metal_atoms = chain_atoms[chain_id]

                # Write protein atoms
                for line in protein_atoms:
                    out.write(line)

                # Write metal atoms immediately after protein atoms in the same chain
                # This keeps them together for pdb2gmx
                for line in metal_atoms:
                    out.write(line)

                # TER record after each chain
                out.write("TER\n")

            # Write other HETATM records (cofactors, modified residues, etc.)
            for line in other_hetatm_lines:
                out.write(line)

            # Final END record
            out.write("END\n")

    def _print_statistics(self):
        """Print cleaning statistics."""
        print("\n=== Cleaning Statistics ===")
        print(f"Total atoms read: {self.stats['total_atoms']}")
        print(f"Protein atoms: {self.stats['protein_atoms']}")
        print(f"Water molecules removed: {self.stats['water_removed']}")
        if self.stats['water_kept'] > 0:
            print(f"Water molecules kept: {self.stats['water_kept']}")
        if self.stats['artifacts_removed'] > 0:
            print(f"Crystallization artifacts removed: {self.stats['artifacts_removed']} atoms")
        print(f"Ions/metals removed: {self.stats['ions_removed']}")
        print(f"Metals kept: {self.stats['metals_kept']}")
        if self.stats['metals_removed_by_distance'] > 0:
            print(f"Metals removed by distance filter: {self.stats['metals_removed_by_distance']}")

        total_kept = self.stats['protein_atoms'] + self.stats['metals_kept'] + self.stats['water_kept']
        print(f"\nTotal atoms in output: {total_kept}")
        print("=" * 30)


def clean_protein_pdb(input_pdb: str, output_pdb: str,
                     ion_mode: str = 'smart',
                     distance_cutoff: float = 5.0,
                     **kwargs) -> dict:
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
        **kwargs
    )
    return cleaner.clean_pdb(input_pdb, output_pdb)


def fix_terminal_atoms(pdb_file: str, output_file: str = None, verbose: bool = True) -> str:
    """
    Fix terminal atom naming and ordering for GROMACS compatibility.

    GROMACS pdb2gmx expects C-terminal residues to have a single 'O' atom
    positioned right after the 'C' (carbonyl) atom. However, pdbfixer may:
    1. Add both backbone 'O' and terminal 'OXT'
    2. Place the terminal oxygen at the end of the residue
    3. Add O atom directly (not as OXT) at the end

    This function:
    1. Identifies C-terminal residues (by finding O/OXT atoms appearing after side chains)
    2. Removes duplicate backbone 'O' atoms
    3. Renames 'OXT' to 'O'
    4. Repositions the O atom to appear right after the C atom

    Parameters
    ----------
    pdb_file : str
        Path to input PDB file
    output_file : str, optional
        Path to output file. If None, overwrites input file.
    verbose : bool
        Print information about fixes

    Returns
    -------
    str
        Path to output file

    Examples
    --------
    >>> fix_terminal_atoms('protein.pdb', 'protein_fixed.pdb')
    >>> fix_terminal_atoms('protein.pdb')  # Overwrites input
    """
    if output_file is None:
        output_file = pdb_file

    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    # Backbone atom names that should appear before O in proper order
    BACKBONE_ATOMS = {'N', 'CA', 'C'}
    # Side chain atoms that indicate O is misplaced if O appears after them
    SIDECHAIN_INDICATORS = {'CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3',
                           'CZ', 'CZ2', 'CZ3', 'CH2', 'NE', 'NE1', 'NE2', 'ND1', 'ND2', 'NZ',
                           'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'SD', 'SG'}

    all_lines = []
    with open(pdb_file, 'r') as f:
        all_lines = f.readlines()

    # First pass: group atoms by residue
    residues = {}  # res_id -> list of (line, index, atom_name)

    for i, line in enumerate(all_lines):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            chain = line[21]
            res_id = f"{chain}:{res_name}:{res_num}"

            if res_id not in residues:
                residues[res_id] = []
            residues[res_id].append({'line': line, 'index': i, 'atom': atom_name, 'res_name': res_name})

    # Second pass: identify C-terminal residues needing fixes
    c_terminal_residues = {}  # res_id -> {c_index, o_line, o_index}

    for res_id, atoms in residues.items():
        c_idx = None
        o_idx = None
        o_line = None
        has_sidechain = False
        oxt_idx = None

        # Find positions of key atoms
        for i, atom_info in enumerate(atoms):
            atom_name = atom_info['atom']

            if atom_name == 'C':
                c_idx = i
            elif atom_name == 'O':
                o_idx = i
                o_line = atom_info['line']
            elif atom_name == 'OXT':
                oxt_idx = i
                # Rename OXT to O
                o_line = atom_info['line'][:12] + ' O  ' + atom_info['line'][16:]
            elif atom_name in SIDECHAIN_INDICATORS:
                has_sidechain = True

        # Determine if this is a C-terminal needing fix
        needs_fix = False

        # Case 1: OXT exists (clear C-terminal)
        if oxt_idx is not None:
            needs_fix = True
            o_idx = oxt_idx
        # Case 2: O appears after side chain atoms (misplaced terminal O)
        elif o_idx is not None and c_idx is not None and has_sidechain:
            # Check if O appears after any side chain atom
            for i, atom_info in enumerate(atoms):
                if atom_info['atom'] in SIDECHAIN_INDICATORS and i < o_idx:
                    needs_fix = True
                    break

        if needs_fix and c_idx is not None and o_line is not None:
            c_terminal_residues[res_id] = {
                'c_index': atoms[c_idx]['index'],
                'o_index': atoms[o_idx]['index'],
                'o_line': o_line,
                'res_name': atoms[0]['res_name'],
                'chain': res_id.split(':')[0],
                'res_num': res_id.split(':')[-1],
                'is_oxt': oxt_idx is not None
            }

    # Third pass: reconstruct file with corrected terminal residues
    skip_indices = set()
    insert_map = {}  # index -> line_to_insert
    fixed_count = 0

    for res_id, term_info in c_terminal_residues.items():
        # Skip the old O/OXT position
        skip_indices.add(term_info['o_index'])
        # Insert O right after C
        insert_map[term_info['c_index']] = term_info['o_line']
        fixed_count += 1

        if verbose:
            fix_type = "OXT → O" if term_info['is_oxt'] else "misplaced O"
            print(f"  Fixed {fix_type} in {term_info['res_name']} {term_info['chain']}{term_info['res_num']} "
                  f"(repositioned after C atom)")

    # Write output with insertions
    with open(output_file, 'w') as f:
        for i, line in enumerate(all_lines):
            if i not in skip_indices:
                f.write(line)
                # Insert O atom after C atom if needed
                if i in insert_map:
                    f.write(insert_map[i])

    if verbose:
        if fixed_count > 0:
            print(f"\nFixed {fixed_count} C-terminal residues with misplaced oxygen atoms")
        else:
            print("No terminal atom issues found")

    return output_file


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Clean protein PDB files with intelligent metal ion handling'
    )
    parser.add_argument('input', help='Input PDB file')
    parser.add_argument('output', help='Output PDB file')
    parser.add_argument('--mode', choices=['keep_all', 'smart', 'remove_all'],
                       default='smart', help='Ion handling mode (default: smart)')
    parser.add_argument('--distance', type=float, default=5.0,
                       help='Distance cutoff for keeping metals in Angstroms (default: 5.0)')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress output messages')
    parser.add_argument('--fix-termini', action='store_true',
                       help='Fix terminal atom names (OXT → O) for GROMACS')

    args = parser.parse_args()

    clean_protein_pdb(
        args.input,
        args.output,
        ion_mode=args.mode,
        distance_cutoff=args.distance,
        verbose=not args.quiet
    )

    if args.fix_termini:
        fix_terminal_atoms(args.output, verbose=not args.quiet)
