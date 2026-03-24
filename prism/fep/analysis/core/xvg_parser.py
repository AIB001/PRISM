#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GROMACS XVG File Parser

Parses dhdl.xvg files output by GROMACS FEP calculations.
"""

from typing import Dict, List, Tuple, Union, Any
from pathlib import Path
import logging
import numpy as np
import pandas as pd


def find_xvg_file(window_dir: Path) -> Path:
    """Find the XVG file for a lambda window."""
    for candidate in ("dhdl.xvg", "prod.xvg", "prod/dhdl.xvg"):
        candidate_path = window_dir / candidate
        if candidate_path.exists():
            return candidate_path
    raise FileNotFoundError(f"No dhdl.xvg/prod.xvg file found in {window_dir}")


def parse_leg_data(
    leg_dir: Path,
    leg_name: str,
    estimator_name: str,
    temperature: float,
    gmx_module: Any,
    logger: logging.Logger,
) -> List[pd.DataFrame]:
    """Parse all lambda-window XVG files for one leg."""
    window_dirs = sorted(leg_dir.glob("window_*"))
    if not window_dirs:
        raise FileNotFoundError(
            f"No lambda window directories found in {leg_dir}. Expected format: window_00, window_01, ..."
        )

    logger.info(f"Found {len(window_dirs)} lambda windows in {leg_name} leg")

    datasets: List[pd.DataFrame] = []
    for window_dir in window_dirs:
        try:
            xvg_file = find_xvg_file(window_dir)
        except FileNotFoundError:
            logger.warning(f"dhdl.xvg/prod.xvg not found in {window_dir}, skipping")
            continue

        try:
            if estimator_name == "TI":
                df = gmx_module.extract_dHdl(xvg_file, T=temperature)
            else:
                df = gmx_module.extract_u_nk(xvg_file, T=temperature)
            datasets.append(df)
            logger.debug(f"Parsed {xvg_file}: {len(df)} frames")
        except Exception as exc:
            logger.error(f"Error parsing {xvg_file}: {exc}")
            raise

    if not datasets:
        raise FileNotFoundError(f"No valid dhdl.xvg files found in {leg_dir}")

    return datasets


class XVGParser:
    """
    GROMACS XVG file parser

    Parses dhdl.xvg files containing dH/dλ values from FEP simulations.

    Notes
    -----
    XVG file format:
    - Comment lines start with @, #, or ;
    - Data columns separated by whitespace
    - First column is usually time (ps), subsequent columns are dH/dλ for different λ windows
    - Legend information in @ s{} legend comments
    """

    def __init__(self):
        self.data = None
        self.metadata = {}
        self._raw_data = None

    def parse(self, xvg_file: Union[str, Path]) -> Dict:
        """
        Parse GROMACS XVG file

        Parameters
        ----------
        xvg_file : Union[str, Path]
            Path to dhdl.xvg file

        Returns
        -------
        Dict
            Parsed data containing:
            - 'time': Time values (ps)
            - 'dhdl': Dictionary of dH/dλ values for each λ
            - 'lambdas': List of λ values
            - 'metadata': File metadata (title, legend, etc.)

        Raises
        ------
        FileNotFoundError
            If xvg_file doesn't exist
        ValueError
            If file format is invalid
        """
        xvg_file = Path(xvg_file)
        if not xvg_file.exists():
            raise FileNotFoundError(f"XVG file not found: {xvg_file}")

        # Read all lines
        with open(xvg_file, "r") as f:
            lines = f.readlines()

        # Extract metadata and data
        metadata_lines, data_lines = self._split_metadata_data(lines)
        self.metadata = self._extract_metadata(metadata_lines)

        # Parse data section
        self._raw_data = self._parse_data_section(data_lines)

        # Organize by lambda
        self.data = self._organize_by_lambda(self._raw_data)

        return self.data

    def _split_metadata_data(self, lines: List[str]) -> Tuple[List[str], List[str]]:
        """
        Split file into metadata and data sections

        Parameters
        ----------
        lines : List[str]
            All file lines

        Returns
        -------
        Tuple[List[str], List[str]]
            (metadata_lines, data_lines)
        """
        metadata_lines = []
        data_lines = []

        in_data = False
        for line in lines:
            stripped = line.strip()
            # Skip empty lines
            if not stripped:
                continue

            # Check if we're in data section
            # Data lines start with a number (possibly negative)
            if stripped[0].isdigit() or stripped[0] == "-":
                in_data = True

            if in_data:
                data_lines.append(line)
            else:
                metadata_lines.append(line)

        return metadata_lines, data_lines

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
        metadata = {"title": "", "xaxis": "", "yaxis": "", "legend": []}

        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith(";"):
                continue

            if stripped.startswith("@"):
                parts = stripped[1:].strip().split(maxsplit=1)
                if len(parts) < 2:
                    continue

                key, value = parts
                key = key.strip().lower()

                if key == "title":
                    metadata["title"] = value
                elif key == "xaxis":
                    metadata["xaxis"] = value
                elif key == "yaxis":
                    metadata["yaxis"] = value
                elif key.startswith("s"):
                    # Legend entry
                    metadata["legend"].append(value)

        return metadata

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

        Raises
        ------
        ValueError
            If data format is invalid
        """
        if not lines:
            raise ValueError("No data lines found in XVG file")

        data = []
        for line in lines:
            # Split by whitespace and convert to float
            try:
                values = [float(x) for x in line.strip().split()]
                if values:  # Skip empty lines
                    data.append(values)
            except ValueError as e:
                raise ValueError(f"Invalid data line: {line}") from e

        if not data:
            raise ValueError("No valid data found in XVG file")

        # Ensure all rows have same number of columns
        n_cols = len(data[0])
        for i, row in enumerate(data):
            if len(row) != n_cols:
                raise ValueError(f"Inconsistent column count at line {i}: {len(row)} vs {n_cols}")

        return np.array(data)

    def _organize_by_lambda(self, data: np.ndarray) -> Dict:
        """
        Organize parsed data by lambda values

        Parameters
        ----------
        data : np.ndarray
            Raw data array [time, col1, col2, ...]

        Returns
        -------
        Dict
            Organized data with time, dhdl by lambda, lambdas, metadata
        """
        # First column is time
        time = data[:, 0]

        # Remaining columns are dH/dλ values
        # Column 2 is typically dH/dλ, columns 3+ might be components
        dhdl_data = data[:, 1] if data.shape[1] > 1 else np.zeros(len(time))

        # Extract lambda values from legend or use defaults
        lambdas = self._extract_lambdas_from_metadata()

        result = {
            "time": time,
            "dhdl": {0.0: dhdl_data},  # Default to single lambda
            "lambdas": lambdas,
            "metadata": self.metadata,
        }

        return result

    def _extract_lambdas_from_metadata(self) -> List[float]:
        """
        Extract lambda values from metadata

        Returns
        -------
        List[float]
            List of lambda values
        """
        # Try to extract from legend
        # Legend format is typically "dH/dλ (λ=0.00)"
        lambdas = []

        for legend_entry in self.metadata.get("legend", []):
            if "lambda" in legend_entry.lower() or "λ" in legend_entry:
                # Try to extract lambda value
                import re

                match = re.search(r"λ\s*=\s*([\d.]+)", legend_entry)
                if match:
                    try:
                        lambdas.append(float(match.group(1)))
                    except ValueError:
                        pass

        return lambdas if lambdas else [0.0]  # Default to lambda=0

    def to_alchemlyb_dataframe(self, _temperature: float = 310.0) -> pd.DataFrame:
        """
        Convert parsed data to alchemlyb DataFrame format

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin (default: 310.0)

        Returns
        -------
        pd.DataFrame
            DataFrame with columns for time, dH/dλ, etc.
        """
        if self.data is None:
            raise ValueError("No data loaded. Call parse() first.")

        # Create DataFrame
        data_dict = {"time": self.data["time"], "dH/dλ": self.data["dhdl"].get(0.0, np.zeros(len(self.data["time"])))}

        df = pd.DataFrame(data_dict)

        # Set index to time (alchemlyb convention)
        df = df.set_index("time")

        return df
