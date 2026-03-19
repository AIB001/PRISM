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
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
import numpy as np

try:
    import alchemlyb
    from alchemlyb.estimators import MBAR, BAR, TI
    from alchemlyb.preprocessing import slicing
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
        energy_components : List[str], optional
            Energy components to analyze (default: ['elec', 'vdw'])
        allow_mock : bool, optional
            Allow mock mode for testing without alchemlyb (default: False)

        Raises
        ------
        ValueError
            If estimator is not available or alchemlyb is not installed
        """
        self.bound_dir = Path(bound_dir)
        self.unbound_dir = Path(unbound_dir)
        self.temperature = temperature
        self.estimator_name = estimator.upper()
        self.energy_components = energy_components or ["elec", "vdw"]
        self.allow_mock = allow_mock

        # Validate estimator
        if not ALCHEMLYB_AVAILABLE and not allow_mock:
            raise ValueError("alchemlyb is required for FEP analysis. " "Install with: pip install alchemlyb")

        if self.estimator_name not in self.ESTIMATORS:
            raise ValueError(f"Unknown estimator: {estimator}. " f"Available: {list(self.ESTIMATORS.keys())}")

        self.estimator = self.ESTIMATORS[self.estimator_name]

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
        2. Computes free energy differences using the selected estimator
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
        if not ALCHEMLYB_AVAILABLE and self.allow_mock:
            return self._analyze_mock()

        self.logger.info(f"Starting FEP analysis with {self.estimator_name}")

        # Step 1: Parse input files
        self.logger.info("Parsing GROMACS output files...")
        bound_data = self._parse_leg_data(self.bound_dir, "bound")
        unbound_data = self._parse_leg_data(self.unbound_dir, "unbound")

        # Step 2: Compute free energy differences
        self.logger.info(f"Computing free energy differences using {self.estimator_name}...")
        delta_g_bound, error_bound = self._compute_free_energy(bound_data)
        delta_g_unbound, error_unbound = self._compute_free_energy(unbound_data)

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
                "n_lambda_windows": len(bound_data),
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
            xvg_file = window_dir / "dhdl.xvg"
            if not xvg_file.exists():
                # Check for production run output
                xvg_file = window_dir / "prod" / "dhdl.xvg"

            if not xvg_file.exists():
                self.logger.warning(f"dhdl.xvg not found in {window_dir}, skipping")
                continue

            try:
                # Parse using alchemlyb
                df = gmx.extract_dHdl(xvg_file, T=self.temperature)
                datasets.append(df)
                self.logger.debug(f"Parsed {xvg_file}: {len(df)} frames")
            except Exception as e:
                self.logger.error(f"Error parsing {xvg_file}: {e}")
                raise

        if not datasets:
            raise FileNotFoundError(f"No valid dhdl.xvg files found in {leg_dir}")

        return datasets

    def _compute_free_energy(self, datasets: List) -> Tuple[float, float]:
        """
        Compute free energy difference using selected estimator

        Parameters
        ----------
        datasets : List
            List of alchemlyb DataFrames

        Returns
        -------
        Tuple[float, float]
            (delta_g, error) in kcal/mol
        """
        # Subsample if necessary (alchemlyb handles this automatically)
        # Fit estimator
        fitted = self.estimator().fit(datasets)

        # Get delta_g and error
        # Note: alchemlyb returns results in kJ/mol, need to convert to kcal/mol
        delta_g_kj = fitted.delta_f_.iloc[-1]  # Last lambda point
        delta_g_kcal = delta_g_kj / 4.184  # Convert kJ to kcal

        # Error estimate (if available)
        if hasattr(fitted, "d_delta_f_"):
            error_kj = fitted.d_delta_f_.iloc[-1]
            error_kcal = error_kj / 4.184
        else:
            error_kcal = 0.0
            self.logger.warning("Error estimate not available for this estimator")

        return delta_g_kcal, error_kcal

    def _analyze_convergence(self, bound_data: List, unbound_data: List) -> Dict:
        """
        Analyze convergence of free energy calculations

        Parameters
        ----------
        bound_data : List
            Bound leg datasets
        unbound_data : List
            Unbound leg datasets

        Returns
        -------
        Dict
            Convergence metrics
        """
        convergence = {
            "hysteresis": None,
            "time_convergence": {},
            "n_samples": sum(len(df) for df in bound_data),
        }

        # Time convergence analysis (subsample by time)
        try:
            from alchemlyb.convergence import forward_backward_convergence

            # This is a placeholder - actual implementation would analyze
            # convergence as a function of simulation time
            convergence["time_convergence"]["status"] = "not_implemented"
        except ImportError:
            convergence["time_convergence"]["status"] = "alchemlyb_convergence_not_available"

        return convergence

    def _decompose_energies(self, bound_data: List, unbound_data: List) -> Dict[str, float]:
        """
        Decompose free energy by components (elec, vdw, etc.)

        Parameters
        ----------
        bound_data : List
            Bound leg datasets
        unbound_data : List
            Unbound leg datasets

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

        generator = HTMLReportGenerator(self.results)
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
