#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP Analyzer - Main analyzer class for FEP calculations

This module provides the main FEPAnalyzer class that orchestrates
the entire FEP analysis workflow, from parsing GROMACS output files
to generating comprehensive HTML reports.

Example Usage:
-------------
    from prism.fep.analysis import FEPAnalyzer

    analyzer = FEPAnalyzer(
        bound_dir='fep_project/bound',
        unbound_dir='fep_project/unbound',
        temperature=310.0,
        estimator='MBAR'
    )

    results = analyzer.analyze()
    report_path = analyzer.generate_html_report('analysis_report.html')

Author: PRISM Team
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
import numpy as np

try:
    from alchemlyb.estimators import MBAR, BAR, TI
    from alchemlyb.parsing import gmx

    ALCHEMLYB_AVAILABLE = True
except ImportError:
    ALCHEMLYB_AVAILABLE = False
    MBAR = None
    BAR = None
    TI = None
    gmx = None
    logging.warning("alchemlyb not available. Install with: pip install alchemlyb")


@dataclass
class FEResults:
    """
    Container for FEP analysis results

    Attributes
    ----------
    delta_g : float
        Free energy difference (kcal/mol)
    delta_g_error : float
        Standard error of delta_g (kcal/mol)
    delta_g_bound : float
        Free energy difference for bound leg (kcal/mol)
    delta_g_unbound : float
        Free energy difference for unbound leg (kcal/mol)
    delta_g_components : Dict[str, float]
        Free energy decomposition by component (elec, vdw, etc.)
    convergence : Dict
        Convergence diagnostics
    metadata : Dict
        Analysis metadata (temperature, lambda windows, etc.)
    """

    delta_g: float = 0.0
    delta_g_error: float = 0.0
    delta_g_bound: float = 0.0
    delta_g_unbound: float = 0.0
    delta_g_components: Dict[str, float] = field(default_factory=dict)
    convergence: Dict = field(default_factory=dict)
    metadata: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """Convert results to dictionary"""
        return {
            "delta_g": self.delta_g,
            "delta_g_error": self.delta_g_error,
            "delta_g_bound": self.delta_g_bound,
            "delta_g_unbound": self.delta_g_unbound,
            "delta_g_components": self.delta_g_components,
            "convergence": self.convergence,
            "metadata": self.metadata,
        }


class FEPAnalyzer:
    """
    Main analyzer for FEP calculations

    This class orchestrates the complete FEP analysis workflow:

    1. Parse GROMACS output files (dhdl.xvg) from bound and unbound legs
    2. Compute free energy differences using BAR, MBAR, or TI
    3. Perform convergence analysis
    4. Generate comprehensive HTML report

    Parameters
    ----------
    bound_dir : Union[str, Path]
        Path to bound leg directory containing lambda window subdirectories
    unbound_dir : Union[str, Path]
        Path to unbound leg directory containing lambda window subdirectories
    temperature : float, optional
        Simulation temperature in Kelvin (default: 310.0)
    estimator : str, optional
        Free energy estimator: 'MBAR', 'BAR', or 'TI' (default: 'MBAR')
    energy_components : List[str], optional
        Energy components to analyze (default: ['elec', 'vdw'])

    Examples
    --------
    >>> analyzer = FEPAnalyzer(
    ...     bound_dir='fep_project/bound',
    ...     unbound_dir='fep_project/unbound',
    ...     temperature=310.0,
    ...     estimator='MBAR'
    ... )
    >>> results = analyzer.analyze()
    >>> print(f"ΔG = {results.delta_g:.2f} ± {results.delta_g_error:.2f} kcal/mol")
    """

    # Estimator mapping
    ESTIMATORS = {
        "MBAR": MBAR if ALCHEMLYB_AVAILABLE else None,
        "BAR": BAR if ALCHEMLYB_AVAILABLE else None,
        "TI": TI if ALCHEMLYB_AVAILABLE else None,
    }

    def __init__(
        self,
        bound_dir: Union[str, Path],
        unbound_dir: Union[str, Path],
        temperature: float = 310.0,
        estimator: str = "MBAR",
        backend: str = "alchemlyb",
        energy_components: Optional[List[str]] = None,
        allow_mock: bool = False,  # For testing without alchemlyb
    ):
        """
        Initialize FEP analyzer

        Parameters
        ----------
        bound_dir : Union[str, Path]
            Path to bound leg directory
        unbound_dir : Union[str, Path]
            Path to unbound leg directory
        temperature : float, optional
            Temperature in Kelvin (default: 310.0)
        estimator : str, optional
            Estimator method: 'MBAR', 'BAR', or 'TI' (default: 'MBAR')
        backend : str, optional
            Backend for analysis: 'alchemlyb' or 'gmx_bar' (default: 'alchemlyb')
            - 'alchemlyb': Use alchemlyb library (supports TI, BAR, MBAR)
            - 'gmx_bar': Use GROMACS gmx bar command (only supports BAR)
        energy_components : List[str], optional
            Energy components to analyze (default: ['elec', 'vdw'])
        allow_mock : bool, optional
            Allow mock mode for testing without alchemlyb (default: False)

        Raises
        ------
        ValueError
            If estimator is not available or required dependencies are missing
        """
        self.bound_dir = Path(bound_dir)
        self.unbound_dir = Path(unbound_dir)
        self.temperature = temperature
        self.estimator_name = estimator.upper()
        self.backend = backend.lower()
        self.energy_components = energy_components or ["elec", "vdw"]
        self.allow_mock = allow_mock

        # Validate backend
        if self.backend not in ["alchemlyb", "gmx_bar"]:
            raise ValueError(f"Unknown backend: {backend}. Available: ['alchemlyb', 'gmx_bar']")

        # Validate backend+estimator compatibility
        if self.backend == "gmx_bar" and self.estimator_name != "BAR":
            raise ValueError(f"gmx_bar backend only supports BAR estimator, got {estimator}")

        # Validate estimator for alchemlyb backend
        if self.backend == "alchemlyb":
            if not ALCHEMLYB_AVAILABLE and not allow_mock:
                raise ValueError("alchemlyb is required for FEP analysis. " "Install with: pip install alchemlyb")

            if self.estimator_name not in self.ESTIMATORS:
                raise ValueError(f"Unknown estimator: {estimator}. " f"Available: {list(self.ESTIMATORS.keys())}")

            self.estimator = self.ESTIMATORS[self.estimator_name]
        else:
            self.estimator = None  # gmx_bar doesn't use alchemlyb estimators

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # Analysis results storage
        self.results: Optional[FEResults] = None
        self.raw_data: Dict = {}

    def analyze(self) -> FEResults:
        """
        Run complete FEP analysis workflow

        This method:
        1. Parses GROMACS output files from bound and unbound legs
        2. Computes free energy differences using the selected estimator and backend
        3. Performs convergence analysis
        4. Decomposes free energy by components

        Returns
        -------
        FEResults
            Analysis results containing free energy differences and errors

        Raises
        ------
        FileNotFoundError
            If required input files are missing
        ValueError
            If data format is incorrect
        """
        # Use mock mode if alchemlyb not available
        if self.backend == "alchemlyb" and not ALCHEMLYB_AVAILABLE and self.allow_mock:
            return self._analyze_mock()

        self.logger.info(f"Starting FEP analysis with {self.backend} backend ({self.estimator_name})")

        # Step 1: Parse input files (for alchemlyb) or prepare file lists (for gmx_bar)
        if self.backend == "alchemlyb":
            self.logger.info("Parsing GROMACS output files...")
            bound_data = self._parse_leg_data(self.bound_dir, "bound")
            unbound_data = self._parse_leg_data(self.unbound_dir, "unbound")

            # Step 2: Compute free energy differences
            self.logger.info(f"Computing free energy differences using {self.estimator_name}...")
            delta_g_bound, error_bound, fitted_bound = self._compute_free_energy_alchemlyb(bound_data)
            delta_g_unbound, error_unbound, fitted_unbound = self._compute_free_energy_alchemlyb(unbound_data)

            # Store fitted models for plotting
            self.raw_data["fitted_bound"] = fitted_bound
            self.raw_data["fitted_unbound"] = fitted_unbound
            self.raw_data["lambda_profiles"] = self._build_lambda_profiles(
                fitted_bound,
                fitted_unbound,
                self.estimator_name,
            )
        else:  # gmx_bar backend
            self.logger.info("Computing free energy differences using gmx bar...")
            delta_g_bound, error_bound = self._compute_free_energy_gmx_bar(self.bound_dir, "bound")
            delta_g_unbound, error_unbound = self._compute_free_energy_gmx_bar(self.unbound_dir, "unbound")

            # For gmx_bar, we don't parse data into memory
            bound_data = []
            unbound_data = []
            self.raw_data["lambda_profiles"] = None

        # Binding free energy = unbound - bound
        delta_g = delta_g_unbound - delta_g_bound
        delta_g_error = np.sqrt(error_bound**2 + error_unbound**2)

        # Step 3: Convergence analysis (if enough data)
        convergence = self._analyze_convergence(bound_data, unbound_data)

        # Step 4: Energy decomposition (if components available)
        components = self._decompose_energies(bound_data, unbound_data)

        # Store results
        self.results = FEResults(
            delta_g=delta_g,
            delta_g_error=delta_g_error,
            delta_g_bound=delta_g_bound,
            delta_g_unbound=delta_g_unbound,
            delta_g_components=components,
            convergence=convergence,
            metadata={
                "temperature": self.temperature,
                "estimator": self.estimator_name,
                "backend": self.backend,
                "n_lambda_windows": len(bound_data) if bound_data else self._count_windows(self.bound_dir),
            },
        )

        self.logger.info(f"Analysis complete: ΔG = {delta_g:.2f} ± {delta_g_error:.2f} kcal/mol")
        return self.results

    def _analyze_mock(self) -> FEResults:
        """
        Mock analysis for testing without alchemlyb

        Returns
        -------
        FEResults
            Mock results with synthetic data
        """
        self.logger.info("Running mock analysis (alchemlyb not available)")

        # Count lambda windows
        bound_windows = len(list(self.bound_dir.glob("window_*")))
        unbound_windows = len(list(self.unbound_dir.glob("window_*")))

        # Generate synthetic results
        delta_g_bound = -2.5  # kcal/mol
        delta_g_unbound = -4.8  # kcal/mol
        delta_g = delta_g_unbound - delta_g_bound
        delta_g_error = 0.15  # kcal/mol

        self.results = FEResults(
            delta_g=delta_g,
            delta_g_error=delta_g_error,
            delta_g_bound=delta_g_bound,
            delta_g_unbound=delta_g_unbound,
            delta_g_components={"elec": -1.2, "vdw": -1.1},
            convergence={"n_samples": 5000, "time_convergence": {"status": "mock_data"}},
            metadata={
                "temperature": self.temperature,
                "estimator": self.estimator_name,
                "n_lambda_windows": bound_windows,
                "mock": True,
            },
        )

        self.logger.info(f"Mock analysis complete: ΔG = {delta_g:.2f} ± {delta_g_error:.2f} kcal/mol")
        return self.results

    def _parse_leg_data(self, leg_dir: Path, leg_name: str) -> List:
        """
        Parse GROMACS output files for a leg (bound or unbound)

        Parameters
        ----------
        leg_dir : Path
            Path to leg directory
        leg_name : str
            Name of leg ('bound' or 'unbound')

        Returns
        -------
        List
            List of alchemlyb DataFrame objects, one per lambda window

        Raises
        ------
        FileNotFoundError
            If dhdl.xvg files are not found
        """
        # Find lambda window directories
        window_dirs = sorted(leg_dir.glob("window_*"))
        if not window_dirs:
            raise FileNotFoundError(
                f"No lambda window directories found in {leg_dir}. " f"Expected format: window_00, window_01, ..."
            )

        self.logger.info(f"Found {len(window_dirs)} lambda windows in {leg_name} leg")

        # Parse dhdl.xvg from each window
        datasets = []
        for window_dir in window_dirs:
            # Check multiple possible filenames (PRISM uses prod.xvg, legacy uses dhdl.xvg)
            xvg_file = None
            for candidate in ["dhdl.xvg", "prod.xvg", "prod/dhdl.xvg"]:
                candidate_path = window_dir / candidate
                if candidate_path.exists():
                    xvg_file = candidate_path
                    break

            if xvg_file is None:
                self.logger.warning(f"dhdl.xvg/prod.xvg not found in {window_dir}, skipping")
                continue

            try:
                # Parse using alchemlyb
                # Choose extraction method based on estimator
                if self.estimator_name == "TI":
                    df = gmx.extract_dHdl(xvg_file, T=self.temperature)
                else:  # BAR or MBAR
                    df = gmx.extract_u_nk(xvg_file, T=self.temperature)
                datasets.append(df)
                self.logger.debug(f"Parsed {xvg_file}: {len(df)} frames")
            except Exception as e:
                self.logger.error(f"Error parsing {xvg_file}: {e}")
                raise

        if not datasets:
            raise FileNotFoundError(f"No valid dhdl.xvg files found in {leg_dir}")

        return datasets

    def _compute_free_energy_alchemlyb(self, datasets: List) -> Tuple[float, float]:
        """
        Compute free energy difference using selected estimator

        Parameters
        ----------
        datasets : List
            List of alchemlyb DataFrames

        Returns
        -------
        Tuple[float, float, object]
            (delta_g, error) in kcal/mol, and fitted estimator object
        """
        import pandas as pd

        # alchemlyb estimators require a single concatenated DataFrame
        combined = pd.concat(datasets)
        fitted = self.estimator().fit(combined)

        # Get delta_g and error
        # Note: alchemlyb returns results in kJ/mol, need to convert to kcal/mol
        # delta_f_ is a DataFrame with states as both rows and columns
        # iloc[0, -1] gives the free energy from first state to last state
        delta_g_kj = float(fitted.delta_f_.iloc[0, -1])  # First to last state
        delta_g_kcal = delta_g_kj / 4.184  # Convert kJ to kcal

        # Error estimate (if available)
        if hasattr(fitted, "d_delta_f_"):
            error_kj = float(fitted.d_delta_f_.iloc[0, -1])
            error_kcal = error_kj / 4.184
        else:
            error_kcal = 0.0
            self.logger.warning("Error estimate not available for this estimator")

        return delta_g_kcal, error_kcal, fitted

    def _compute_free_energy_gmx_bar(self, leg_dir: Path, leg_name: str) -> Tuple[float, float]:
        """
        Compute free energy difference using GROMACS gmx bar command

        Parameters
        ----------
        leg_dir : Path
            Path to leg directory
        leg_name : str
            Name of leg ('bound' or 'unbound')

        Returns
        -------
        Tuple[float, float]
            (delta_g, error) in kcal/mol
        """
        import subprocess
        import re

        # Find all dhdl.xvg files
        xvg_files = []
        window_dirs = sorted(leg_dir.glob("window_*"))

        for window_dir in window_dirs:
            # Check multiple possible filenames
            for candidate in ["dhdl.xvg", "prod.xvg", "prod/dhdl.xvg"]:
                xvg_file = window_dir / candidate
                if xvg_file.exists():
                    xvg_files.append(str(xvg_file))
                    break

        if not xvg_files:
            raise FileNotFoundError(f"No dhdl.xvg files found in {leg_dir}")

        self.logger.info(f"Found {len(xvg_files)} XVG files for {leg_name} leg")

        # Run gmx bar
        # Output to a temporary file
        output_file = leg_dir / f"bar_{leg_name}.xvg"

        cmd = [
            "gmx",
            "bar",
            "-f",
            *xvg_files,
            "-o",
            str(output_file),
            "-temp",
            str(self.temperature),
            "-prec",
            "6",
        ]

        self.logger.info(f"Running: gmx bar -f ... -o {output_file}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Parse output from stdout (gmx bar prints results to stdout)
            # Format: "Total free energy difference: X.XXX +/- Y.YYY kT"
            # Or: "Delta G = X.XXX +/- Y.YYY kT"
            stdout = result.stdout

            # Try to find the final result line
            # gmx bar output format varies, look for patterns like:
            # "Total Delta G:  XXXX +/- YYYY kT"
            # "Final estimate: XXXX +/- YYYY"

            delta_g_kt = None
            error_kt = None

            # Pattern 1: Look for "total" line (most reliable)
            # Format: "total      0 -     31,   DG -25.137419 +/- 15.177305"
            for line in stdout.split("\n"):
                if line.strip().startswith("total"):
                    # Extract numbers using regex
                    match = re.search(r"DG\s+([+-]?\d+\.?\d*)\s*\+/-\s*([+-]?\d+\.?\d*)", line)
                    if match:
                        delta_g_kt = float(match.group(1))
                        error_kt = float(match.group(2))
                        self.logger.info(f"Found total result: {delta_g_kt:.2f} ± {error_kt:.2f} kT")
                        break

            # Pattern 2: Fallback - look for "Total Delta G" or "Final estimate"
            if delta_g_kt is None:
                for line in stdout.split("\n"):
                    if "Delta G" in line or "estimate" in line:
                        match = re.search(r"([+-]?\d+\.?\d*)\s*\+/-\s*([+-]?\d+\.?\d*)", line)
                        if match:
                            delta_g_kt = float(match.group(1))
                            error_kt = float(match.group(2))
                            break

            if delta_g_kt is None:
                # Fallback: parse the output file
                if output_file.exists():
                    with open(output_file) as f:
                        for line in f:
                            if not line.startswith("#") and not line.startswith("@"):
                                parts = line.split()
                                if len(parts) >= 2:
                                    try:
                                        delta_g_kt = float(parts[0])
                                        error_kt = float(parts[1])
                                        break
                                    except ValueError:
                                        continue

            if delta_g_kt is None:
                raise ValueError("Could not parse gmx bar output")

            # Convert from kT to kJ/mol
            # kT at temperature T: kT = R * T, where R = 8.314 J/(mol·K)
            R = 8.314 / 1000  # kJ/(mol·K)
            kT_to_kj = R * self.temperature

            delta_g_kj = delta_g_kt * kT_to_kj
            error_kj = error_kt * kT_to_kj if error_kt else 0.0

            # Convert to kcal/mol
            delta_g_kcal = delta_g_kj / 4.184
            error_kcal = error_kj / 4.184

            self.logger.info(f"gmx bar result: {delta_g_kcal:.2f} ± {error_kcal:.2f} kcal/mol")

            return delta_g_kcal, error_kcal

        except subprocess.CalledProcessError as e:
            self.logger.error(f"gmx bar failed: {e}")
            self.logger.error(f"stderr: {e.stderr}")
            raise RuntimeError(f"gmx bar command failed: {e}")

    def _build_lambda_profiles(
        self,
        fitted_bound: Any,
        fitted_unbound: Any,
        estimator_name: str,
    ) -> Dict[str, Optional[Dict[str, List[float]]]]:
        """
        Build lambda-dependent profiles for bound and unbound legs.

        Parameters
        ----------
        fitted_bound : Any
            Fitted estimator for bound leg.
        fitted_unbound : Any
            Fitted estimator for unbound leg.
        estimator_name : str
            Estimator name (TI, BAR, MBAR).

        Returns
        -------
        Dict[str, Optional[Dict[str, List[float]]]]
            Dictionary with bound/unbound lambda profiles.
        """
        estimator = estimator_name.upper()
        return {
            "bound": self._extract_lambda_data(fitted_bound, estimator),
            "unbound": self._extract_lambda_data(fitted_unbound, estimator),
        }

    def _extract_lambda_data(self, fitted_model: Any, estimator_name: str) -> Optional[Dict[str, List[float]]]:
        """
        Extract lambda-dependent data from a fitted estimator.

        Parameters
        ----------
        fitted_model : Any
            Fitted estimator instance.
        estimator_name : str
            Estimator name (TI, BAR, MBAR).

        Returns
        -------
        Optional[Dict[str, List[float]]]
            Lambda profile with lambdas, dhdl, cumulative_dg, state_dg.
        """
        if fitted_model is None:
            return None

        lambda_data: Dict[str, List[float]] = {
            "lambdas": [],
            "dhdl": [],
            "cumulative_dg": [],
            "state_dg": [],
        }

        if estimator_name == "TI" and hasattr(fitted_model, "dhdl"):
            dhdl_df = fitted_model.dhdl
            n_states = len(dhdl_df)
            lambda_data["lambdas"] = list(range(n_states))
            lambda_data["dhdl"] = [float(v) for v in dhdl_df.sum(axis=1).values]
            if hasattr(fitted_model, "delta_f_"):
                lambda_data["cumulative_dg"] = [float(v) for v in fitted_model.delta_f_.iloc[0].values]
        elif hasattr(fitted_model, "delta_f_"):
            delta_f = fitted_model.delta_f_
            n_states = delta_f.shape[0]
            if n_states > 1:
                lambda_data["lambdas"] = [float(i) / (n_states - 1) for i in range(n_states)]
            else:
                lambda_data["lambdas"] = [0.0]

            for i in range(n_states):
                if i == 0:
                    lambda_data["state_dg"].append(0.0)
                else:
                    lambda_data["state_dg"].append(float(delta_f.iloc[0, i]))

            lambda_data["cumulative_dg"] = lambda_data["state_dg"].copy()
            lambda_data["dhdl"] = [0.0] * n_states

        if not lambda_data["lambdas"]:
            return None

        return lambda_data

    def _analyze_convergence(self, bound_data: List, unbound_data: List) -> Dict:  # noqa: ARG002
        """
        Analyze convergence of free energy calculations

        Parameters
        ----------
        bound_data : List
            Bound leg datasets
        unbound_data : List
            Unbound leg datasets (not used in current implementation)

        Returns
        -------
        Dict
            Convergence metrics
        """
        convergence = {
            "hysteresis": None,
            "time_convergence": {},
            "n_samples": sum(len(df) for df in bound_data) if bound_data else 0,
        }

        # Time convergence analysis (subsample by time)
        # TODO: Implement actual convergence analysis
        convergence["time_convergence"]["status"] = "not_implemented"

        return convergence

    def _decompose_energies(self, bound_data: List, unbound_data: List) -> Dict[str, float]:  # noqa: ARG002
        """
        Decompose free energy by components (elec, vdw, etc.)

        Parameters
        ----------
        bound_data : List
            Bound leg datasets
        unbound_data : List
            Unbound leg datasets (not used in current implementation)

        Returns
        -------
        Dict[str, float]
            Free energy decomposition by component
        """
        components = {}

        # Check if datasets have component columns
        if bound_data and "dH/dλ" in bound_data[0].columns:
            # Simple decomposition - just total dH/dλ
            # Full decomposition would require parsing separate component files
            components["total"] = self.results.delta_g if self.results else 0.0

        return components

    def generate_html_report(self, output_path: Union[str, Path]) -> Path:
        """
        Generate comprehensive HTML report

        Parameters
        ----------
        output_path : Union[str, Path]
            Path to output HTML file

        Returns
        -------
        Path
            Path to generated HTML file

        Raises
        ------
        RuntimeError
            If analyze() has not been called first
        """
        if self.results is None:
            raise RuntimeError("Must call analyze() before generating report")

        from prism.fep.analysis.report import HTMLReportGenerator

        generator = HTMLReportGenerator(self.results, self.raw_data)
        html_path = generator.generate(output_path)

        self.logger.info(f"HTML report generated: {html_path}")
        return html_path

    def save_results(self, output_path: Union[str, Path]) -> Path:
        """
        Save analysis results to JSON file

        Parameters
        ----------
        output_path : Union[str, Path]
            Path to output JSON file

        Returns
        -------
        Path
            Path to saved JSON file
        """
        if self.results is None:
            raise RuntimeError("Must call analyze() before saving results")

        output_path = Path(output_path)
        with open(output_path, "w") as f:
            json.dump(self.results.to_dict(), f, indent=2)

        self.logger.info(f"Results saved to {output_path}")
        return output_path
