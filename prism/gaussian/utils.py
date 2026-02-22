#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utility functions for Gaussian RESP charge calculations.

This module provides helper functions for checking dependencies,
calculating molecular properties, and managing Gaussian calculations.
"""

import os
import subprocess
import shutil
from typing import Tuple, Optional

# Atomic symbols mapping (atomic number to element symbol)
ATOMIC_SYMBOLS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
    16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
    23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
    30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
    37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
    44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
    58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
    65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
    72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
    79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
    86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U'
}

# Reverse mapping: symbol to atomic number
SYMBOL_TO_NUMBER = {v: k for k, v in ATOMIC_SYMBOLS.items()}


def normalize_element_symbol(element: str) -> str:
    """
    Normalize element symbol to standard format (first letter uppercase, rest lowercase).

    Examples: "CL" -> "Cl", "BR" -> "Br", "c" -> "C", "na" -> "Na"

    Parameters
    ----------
    element : str
        Raw element symbol

    Returns
    -------
    str
        Normalized element symbol
    """
    element = element.strip()
    if len(element) == 0:
        return element
    elif len(element) == 1:
        return element.upper()
    else:
        return element[0].upper() + element[1:].lower()


# Number of valence electrons for common elements
VALENCE_ELECTRONS = {
    'H': 1, 'He': 2,
    'Li': 1, 'Be': 2, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8,
    'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8,
    'K': 1, 'Ca': 2, 'Br': 7, 'I': 7
}


def check_gaussian_available() -> bool:
    """
    Check if Gaussian (g16) is available in the system PATH.

    Returns
    -------
    bool
        True if g16 is available, False otherwise
    """
    return shutil.which('g16') is not None


def check_antechamber_available() -> bool:
    """
    Check if antechamber is available in the system PATH.

    Returns
    -------
    bool
        True if antechamber is available, False otherwise
    """
    return shutil.which('antechamber') is not None


def calculate_charge_multiplicity(mol2_file: str, verbose: bool = False) -> Tuple[int, int]:
    """
    Calculate the net charge and spin multiplicity of a molecule from MOL2 file.

    The net charge is calculated from the sum of partial charges in the MOL2 file.
    The multiplicity is determined by the number of electrons:
    - Even electrons: singlet (multiplicity = 1)
    - Odd electrons: doublet (multiplicity = 2)

    Parameters
    ----------
    mol2_file : str
        Path to the MOL2 file
    verbose : bool
        Whether to print debug information

    Returns
    -------
    Tuple[int, int]
        (charge, multiplicity) tuple
    """
    total_charge = 0.0
    total_electrons = 0
    in_atom_section = False
    atom_count = 0

    with open(mol2_file, 'r') as f:
        for line in f:
            line_stripped = line.strip()

            if line_stripped.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            elif line_stripped.startswith("@<TRIPOS>"):
                in_atom_section = False
                continue

            if in_atom_section and line_stripped:
                parts = line.split()
                if len(parts) >= 6:
                    atom_count += 1
                    # Extract element from atom type
                    atom_type = parts[5]
                    if '.' in atom_type:
                        element = normalize_element_symbol(atom_type.split('.')[0])
                    else:
                        element = normalize_element_symbol(atom_type)

                    # Count electrons from element
                    if element in SYMBOL_TO_NUMBER:
                        total_electrons += SYMBOL_TO_NUMBER[element]
                    else:
                        # Try to extract element from atom name
                        atom_name = parts[1]
                        found = False
                        for length in [2, 1]:
                            if len(atom_name) >= length:
                                alpha_chars = ''.join(c for c in atom_name[:length+1] if c.isalpha())
                                if len(alpha_chars) >= length:
                                    candidate = normalize_element_symbol(alpha_chars[:length])
                                    if candidate in SYMBOL_TO_NUMBER:
                                        total_electrons += SYMBOL_TO_NUMBER[candidate]
                                        found = True
                                        break
                        if not found and verbose:
                            print(f"  [WARNING] Could not determine element for atom '{atom_name}' type '{atom_type}'")

                    # Extract charge if present
                    if len(parts) >= 9:
                        try:
                            charge = float(parts[8])
                            total_charge += charge
                        except ValueError:
                            pass

    # Round charge to nearest integer
    net_charge = int(round(total_charge))

    # Calculate actual electrons (total electrons - charge)
    # If charge is +1, we have 1 less electron
    # If charge is -1, we have 1 more electron
    actual_electrons = total_electrons - net_charge

    # Multiplicity: 1 for even electrons (singlet), 2 for odd (doublet)
    if actual_electrons % 2 == 0:
        multiplicity = 1  # Singlet (closed-shell)
    else:
        multiplicity = 2  # Doublet (open-shell)

    if verbose:
        print(f"  [DEBUG] Atoms: {atom_count}, Total electrons: {total_electrons}")
        print(f"  [DEBUG] Sum of partial charges: {total_charge:.4f} -> Net charge: {net_charge}")
        print(f"  [DEBUG] Actual electrons: {actual_electrons} ({'even' if actual_electrons % 2 == 0 else 'odd'})")
        print(f"  [DEBUG] Multiplicity: {multiplicity}")

    return net_charge, multiplicity


def calculate_charge_multiplicity_rdkit(mol2_file: str, verbose: bool = False) -> Tuple[int, int]:
    """
    Calculate charge and multiplicity using RDKit (more accurate).

    Parameters
    ----------
    mol2_file : str
        Path to the MOL2 file
    verbose : bool
        Whether to print debug information

    Returns
    -------
    Tuple[int, int]
        (charge, multiplicity) tuple

    Raises
    ------
    ImportError
        If RDKit is not available
    """
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("RDKit is required for this function. Install with: pip install rdkit")

    # Read MOL2 file
    mol = Chem.MolFromMol2File(mol2_file, removeHs=False)
    if mol is None:
        if verbose:
            print("  [DEBUG] RDKit failed to read MOL2, falling back to simple method")
        # Fallback to simple method
        return calculate_charge_multiplicity(mol2_file, verbose=verbose)

    # Get formal charge
    net_charge = Chem.GetFormalCharge(mol)

    # Count electrons
    total_electrons = 0
    for atom in mol.GetAtoms():
        total_electrons += atom.GetAtomicNum()

    # Adjust for charge
    actual_electrons = total_electrons - net_charge

    # Determine multiplicity
    # Check for radical electrons
    num_radical_electrons = 0
    for atom in mol.GetAtoms():
        num_radical_electrons += atom.GetNumRadicalElectrons()

    if num_radical_electrons > 0:
        multiplicity = num_radical_electrons + 1
    elif actual_electrons % 2 == 0:
        multiplicity = 1  # Singlet
    else:
        multiplicity = 2  # Doublet

    if verbose:
        print(f"  [DEBUG] RDKit: Atoms: {mol.GetNumAtoms()}, Total electrons: {total_electrons}")
        print(f"  [DEBUG] RDKit: Formal charge: {net_charge}, Radical electrons: {num_radical_electrons}")
        print(f"  [DEBUG] RDKit: Actual electrons: {actual_electrons} ({'even' if actual_electrons % 2 == 0 else 'odd'})")
        print(f"  [DEBUG] RDKit: Multiplicity: {multiplicity}")

    return net_charge, multiplicity


def get_charge_multiplicity(mol2_file: str, use_rdkit: bool = False, verbose: bool = False) -> Tuple[int, int]:
    """
    Get charge and multiplicity from MOL2 file.

    Uses simple MOL2 parsing by default (more reliable for electron counting).
    RDKit can be enabled but may give incorrect electron counts for some molecules.

    The multiplicity is automatically determined based on electron count:
    - Even electrons: singlet (multiplicity = 1)
    - Odd electrons: doublet (multiplicity = 2)

    Parameters
    ----------
    mol2_file : str
        Path to the MOL2 file
    use_rdkit : bool
        Whether to use RDKit (default: False, simple parsing is more reliable)
    verbose : bool
        Whether to print debug information (default: False)

    Returns
    -------
    Tuple[int, int]
        (charge, multiplicity) tuple
    """
    # Always calculate with simple method first (more reliable for electron counting)
    simple_charge, simple_mult = calculate_charge_multiplicity(mol2_file, verbose=verbose)

    if use_rdkit:
        try:
            rdkit_charge, rdkit_mult = calculate_charge_multiplicity_rdkit(mol2_file, verbose=verbose)
            # If results differ, warn and prefer simple method
            if simple_mult != rdkit_mult:
                if verbose:
                    print(f"  [WARNING] RDKit and simple method disagree on multiplicity!")
                    print(f"  [WARNING] Simple: {simple_mult}, RDKit: {rdkit_mult}")
                    print(f"  [WARNING] Using simple method result (more reliable for electron counting)")
        except ImportError:
            if verbose:
                print("  [DEBUG] RDKit not available")

    return simple_charge, simple_mult


def run_gaussian(gjf_file: str, timeout: Optional[int] = None) -> Tuple[bool, str]:
    """
    Run Gaussian calculation.

    Parameters
    ----------
    gjf_file : str
        Path to Gaussian input file (.gjf)
    timeout : int, optional
        Timeout in seconds (default: None = no timeout)

    Returns
    -------
    Tuple[bool, str]
        (success, log_file_path) tuple
    """
    if not check_gaussian_available():
        raise RuntimeError("Gaussian (g16) not found in PATH")

    log_file = os.path.splitext(gjf_file)[0] + ".log"

    try:
        subprocess.run(
            ['g16', gjf_file],
            cwd=os.path.dirname(gjf_file) or '.',
            timeout=timeout,
            capture_output=True,
            text=True
        )

        # Check if calculation completed normally
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                content = f.read()
                if "Normal termination" in content:
                    return True, log_file

        return False, log_file

    except subprocess.TimeoutExpired:
        return False, log_file
    except Exception as e:
        return False, str(e)


def extract_optimized_coords_from_log(log_file: str) -> list:
    """
    Extract optimized coordinates from Gaussian log file.

    Parameters
    ----------
    log_file : str
        Path to Gaussian log file

    Returns
    -------
    list
        List of (element, x, y, z) tuples
    """
    atoms = []

    with open(log_file, 'r') as f:
        content = f.read()

    # Find the last "Standard orientation" section
    import re
    pattern = r"Standard orientation:.*?-+\n.*?-+\n(.*?)-+\n"
    matches = re.findall(pattern, content, re.DOTALL)

    if not matches:
        # Try "Input orientation" as fallback
        pattern = r"Input orientation:.*?-+\n.*?-+\n(.*?)-+\n"
        matches = re.findall(pattern, content, re.DOTALL)

    if matches:
        # Use the last orientation found
        coord_block = matches[-1]
        for line in coord_block.strip().split('\n'):
            parts = line.split()
            if len(parts) >= 6:
                try:
                    atomic_num = int(parts[1])
                    x = float(parts[3])
                    y = float(parts[4])
                    z = float(parts[5])

                    element = ATOMIC_SYMBOLS.get(atomic_num, 'X')
                    atoms.append((element, x, y, z))
                except (ValueError, IndexError):
                    continue

    return atoms


def verify_gaussian_termination(log_file: str) -> bool:
    """
    Verify that Gaussian calculation terminated normally.

    Parameters
    ----------
    log_file : str
        Path to Gaussian log file

    Returns
    -------
    bool
        True if calculation terminated normally
    """
    if not os.path.exists(log_file):
        return False

    with open(log_file, 'r') as f:
        content = f.read()

    return "Normal termination" in content
