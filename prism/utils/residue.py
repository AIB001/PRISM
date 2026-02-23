"""
Amino acid residue name conversion utilities.

This module provides functions to convert between 3-letter and 1-letter amino acid codes,
enabling flexible residue name formatting throughout PRISM.
"""

import re
from typing import Union, List

# Three-letter to one-letter amino acid code mapping
AA_3TO1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

# One-letter to three-letter amino acid code mapping (reverse)
AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


def convert_residue_to_1letter(residue_name: str) -> str:
    """
    Convert residue name from 3-letter to 1-letter code.

    Parameters
    ----------
    residue_name : str
        Residue name in 3-letter format (e.g., "ASP618", "SER759")

    Returns
    -------
    str
        Residue name in 1-letter format (e.g., "D618", "S759")
        If already in 1-letter format or unrecognized, returns unchanged

    Examples
    --------
    >>> convert_residue_to_1letter("ASP618")
    'D618'
    >>> convert_residue_to_1letter("SER759")
    'S759'
    >>> convert_residue_to_1letter("D618")  # Already 1-letter
    'D618'
    """
    # Pattern: 3-letter code followed by digits
    match = re.match(r"^([A-Z]{3})(\d+)$", residue_name)

    if match:
        three_letter = match.group(1)
        residue_num = match.group(2)

        if three_letter in AA_3TO1:
            return f"{AA_3TO1[three_letter]}{residue_num}"

    # Return unchanged if not matching pattern or not recognized
    return residue_name


def convert_residue_to_3letter(residue_name: str) -> str:
    """
    Convert residue name from 1-letter to 3-letter code.

    Parameters
    ----------
    residue_name : str
        Residue name in 1-letter format (e.g., "D618", "S759")

    Returns
    -------
    str
        Residue name in 3-letter format (e.g., "ASP618", "SER759")
        If already in 3-letter format or unrecognized, returns unchanged

    Examples
    --------
    >>> convert_residue_to_3letter("D618")
    'ASP618'
    >>> convert_residue_to_3letter("S759")
    'SER759'
    >>> convert_residue_to_3letter("ASP618")  # Already 3-letter
    'ASP618'
    """
    # Pattern: 1-letter code followed by digits
    match = re.match(r"^([A-Z])(\d+)$", residue_name)

    if match:
        one_letter = match.group(1)
        residue_num = match.group(2)

        if one_letter in AA_1TO3:
            return f"{AA_1TO3[one_letter]}{residue_num}"

    # Return unchanged if not matching pattern or not recognized
    return residue_name


def format_residue(residue_name: str, format: str = "1letter") -> str:
    """
    Format residue name according to specified format.

    Parameters
    ----------
    residue_name : str
        Residue name in any format
    format : str
        Target format: "1letter" or "3letter"

    Returns
    -------
    str
        Formatted residue name

    Examples
    --------
    >>> format_residue("ASP618", "1letter")
    'D618'
    >>> format_residue("D618", "3letter")
    'ASP618'
    >>> format_residue("ASP618", "3letter")
    'ASP618'
    """
    if format == "1letter":
        return convert_residue_to_1letter(residue_name)
    elif format == "3letter":
        return convert_residue_to_3letter(residue_name)
    else:
        raise ValueError(f"Unknown residue format: {format}. Must be '1letter' or '3letter'")


def format_residue_list(residue_names: Union[List[str], tuple], format: str = "1letter") -> List[str]:
    """
    Format a list of residue names according to specified format.

    Parameters
    ----------
    residue_names : list or tuple
        List of residue names
    format : str
        Target format: "1letter" or "3letter"

    Returns
    -------
    list
        List of formatted residue names

    Examples
    --------
    >>> format_residue_list(["ASP618", "SER759", "ASN691"], "1letter")
    ['D618', 'S759', 'N691']
    >>> format_residue_list(["D618", "S759"], "3letter")
    ['ASP618', 'SER759']
    """
    return [format_residue(res, format) for res in residue_names]


def normalize_residue_to_3letter(residue_name: str) -> str:
    """
    Normalize residue name to 3-letter format for MDTraj/MDAnalysis selections.

    This function ensures that residue names are in 3-letter format required
    for trajectory analysis selections, regardless of user input format.

    Parameters
    ----------
    residue_name : str
        Residue name in any format

    Returns
    -------
    str
        Residue name in 3-letter format

    Examples
    --------
    >>> normalize_residue_to_3letter("D618")
    'ASP618'
    >>> normalize_residue_to_3letter("ASP618")
    'ASP618'
    """
    return convert_residue_to_3letter(residue_name)
