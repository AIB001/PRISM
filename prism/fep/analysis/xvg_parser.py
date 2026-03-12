#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GROMACS XVG File Parser

Parses dhdl.xvg files output by GROMACS FEP calculations.
"""

from typing import Dict, List
import numpy as np


class XVGParser:
    """
    GROMACS XVG file parser

    Parses dhdl.xvg files containing dH/dλ values from FEP simulations.

    Notes
    -----
    XVG file format:
    - Comment lines start with @, #, or ;
    - Data columns separated by whitespace
    - First column is usually time, subsequent columns are dH/dλ for different λ windows
    """

    def __init__(self):
        self.data = None
        self.metadata = {}

    def parse(self, xvg_file: str) -> Dict:
        """
        Parse GROMACS XVG file

        Parameters
        ----------
        xvg_file : str
            Path to dhdl.xvg file

        Returns
        -------
        Dict
            Parsed data containing:
            - 'lambdas': List of λ values
            - 'dhdl': 2D array of dH/dλ values [time, lambda]
            - 'metadata': File metadata (title, legend, etc.)

        Raises
        ------
        NotImplementedError
            To be implemented in next phase
        """
        raise NotImplementedError("XVG file parsing to be implemented")

    def _extract_metadata(self, lines: List[str]) -> Dict:
        """
        Extract metadata from XVG file header

        Parameters
        ----------
        lines : List[str]
            File lines (excluding data section)

        Returns
        -------
        Dict
            Metadata dictionary with keys like 'title', 'xaxis', 'yaxis', 'legend'
        """
        raise NotImplementedError("_extract_metadata to be implemented")

    def _parse_data_section(self, lines: List[str]) -> np.ndarray:
        """
        Parse numerical data section

        Parameters
        ----------
        lines : List[str]
            Data lines (numerical values only)

        Returns
        -------
        np.ndarray
            2D array of shape (n_timesteps, n_observables)
        """
        raise NotImplementedError("_parse_data_section to be implemented")
