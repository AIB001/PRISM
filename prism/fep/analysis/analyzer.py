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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from .convergence import compute_bootstrap, compute_time_convergence
from .estimators import compute_free_energy_alchemlyb, compute_free_energy_gmx_bar, summarize_repeat_results
from .profiles import build_lambda_profiles
from .xvg_parser import parse_leg_data

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
    repeat_results : List[Dict[str, float]]
        Individual results for each repeat (if multiple repeats analyzed)
        Format: [{"repeat": 1, "bound": 54.07, "unbound": 49.84, "ddG": 4.23}, ...]
    n_repeats : int
        Number of repeats analyzed
    lambda_profiles : Optional[Dict]
        Lambda-dependent profiles for plotting (dg vs lambda, dhdl vs lambda)
    """

    delta_g: float = 0.0
    delta_g_error: float = 0.0
    delta_g_bound: float = 0.0
    delta_g_unbound: float = 0.0
    delta_g_components: Dict[str, float] = field(default_factory=dict)
    convergence: Dict = field(default_factory=dict)
    metadata: Dict = field(default_factory=dict)
    repeat_results: List[Dict[str, float]] = field(default_factory=list)
    n_repeats: int = 1
    lambda_profiles: Optional[Dict[str, Any]] = None

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
            "repeat_results": self.repeat_results,
            "n_repeats": self.n_repeats,
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
        bound_dir: Union[str, Path, List[Union[str, Path]]],
        unbound_dir: Union[str, Path, List[Union[str, Path]]],
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
        # Support single directory or list of directories (multiple repeats)
        if isinstance(bound_dir, (list, tuple)):
            self.bound_dirs = [Path(d) for d in bound_dir]
        else:
            self.bound_dirs = [Path(bound_dir)]

        if isinstance(unbound_dir, (list, tuple)):
            self.unbound_dirs = [Path(d) for d in unbound_dir]
        else:
            self.unbound_dirs = [Path(unbound_dir)]

        # For backward compatibility, keep the old attributes
        self.bound_dir = self.bound_dirs[0]
        self.unbound_dir = self.unbound_dirs[0]

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
                raise ValueError("alchemlyb is required for FEP analysis. Install with: pip install alchemlyb")

            if self.estimator_name not in self.ESTIMATORS:
                raise ValueError(f"Unknown estimator: {estimator}. Available: {list(self.ESTIMATORS.keys())}")

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

        If multiple repeats are provided (via bound_dir/unbound_dir as lists),
        all repeats are analyzed and results are averaged.

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
        self.logger.info(
            f"Processing {len(self.bound_dirs)} bound repeats and {len(self.unbound_dirs)} unbound repeats"
        )

        # Check if we have multiple repeats
        has_multiple_repeats = len(self.bound_dirs) > 1 or len(self.unbound_dirs) > 1

        # Step 1: Parse input files (for alchemlyb) or prepare file lists (for gmx_bar)
        if self.backend == "alchemlyb":
            self.logger.info("Parsing GROMACS output files...")

            # Process all repeats
            bound_results = []
            unbound_results = []
            fitted_bounds = []
            fitted_unbounds = []
            bound_data_list = []
            unbound_data_list = []

            for i, bound_dir in enumerate(self.bound_dirs):
                repeat_label = f"repeat{i + 1}" if has_multiple_repeats else "bound"
                self.logger.debug(f"Processing {repeat_label} from {bound_dir}")
                bdata = parse_leg_data(
                    bound_dir,
                    "bound",
                    self.estimator_name,
                    self.temperature,
                    gmx,
                    self.logger,
                )
                bound_data_list.append(bdata)
                dg, err, fitted, overlap_m = compute_free_energy_alchemlyb(self.estimator, bdata, self.logger)
                bound_results.append((dg, err))
                fitted_bounds.append(fitted)
                if overlap_m is not None and self.raw_data.get("overlap_matrix_bound") is None:
                    self.raw_data["overlap_matrix_bound"] = overlap_m

            for i, unbound_dir in enumerate(self.unbound_dirs):
                repeat_label = f"repeat{i + 1}" if has_multiple_repeats else "unbound"
                self.logger.debug(f"Processing {repeat_label} from {unbound_dir}")
                udata = parse_leg_data(
                    unbound_dir,
                    "unbound",
                    self.estimator_name,
                    self.temperature,
                    gmx,
                    self.logger,
                )
                unbound_data_list.append(udata)
                dg, err, fitted, overlap_m = compute_free_energy_alchemlyb(self.estimator, udata, self.logger)
                unbound_results.append((dg, err))
                fitted_unbounds.append(fitted)
                if overlap_m is not None and self.raw_data.get("overlap_matrix_unbound") is None:
                    self.raw_data["overlap_matrix_unbound"] = overlap_m

            summary = summarize_repeat_results(bound_results, unbound_results)
            delta_g_bound = summary["delta_g_bound"]
            error_bound = summary["error_bound"]
            delta_g_unbound = summary["delta_g_unbound"]
            error_unbound = summary["error_unbound"]
            repeat_results = summary["repeat_results"]
            n_repeats = summary["n_repeats"]
            repeat_statistics = summary["repeat_statistics"]
            delta_g_bound_values = summary["delta_g_bound_values"]
            delta_g_unbound_values = summary["delta_g_unbound_values"]

            # Store fitted models for plotting
            self.raw_data["bound_data_list"] = bound_data_list
            self.raw_data["unbound_data_list"] = unbound_data_list
            self.raw_data["fitted_bound"] = fitted_bounds[0] if len(fitted_bounds) == 1 else fitted_bounds
            self.raw_data["fitted_unbound"] = fitted_unbounds[0] if len(fitted_unbounds) == 1 else fitted_unbounds
            self.raw_data["lambda_profiles"] = build_lambda_profiles(
                fitted_bounds if has_multiple_repeats else fitted_bounds[0],
                fitted_unbounds if has_multiple_repeats else fitted_unbounds[0],
                self.estimator_name,
            )

            self.raw_data["time_convergence"] = compute_time_convergence(
                bound_data_list, unbound_data_list, self.estimator, n_points=10
            )
            self.raw_data["bootstrap"] = compute_bootstrap(
                bound_data_list, unbound_data_list, self.estimator, n_bootstrap=50, fraction=0.8
            )

            if has_multiple_repeats:
                self.logger.info(f"Bound repeats: {delta_g_bound_values}")
                self.logger.info(f"Unbound repeats: {delta_g_unbound_values}")
                self.logger.info(f"Repeat statistics: {repeat_statistics}")
        else:  # gmx_bar backend
            self.logger.info("Computing free energy differences using gmx bar...")
            delta_g_bound, error_bound = compute_free_energy_gmx_bar(
                self.bound_dir, "bound", self.temperature, self.logger
            )
            delta_g_unbound, error_unbound = compute_free_energy_gmx_bar(
                self.unbound_dir, "unbound", self.temperature, self.logger
            )

            # For gmx_bar, we don't parse data into memory
            bound_data_list = []
            unbound_data_list = []
            repeat_results = []
            n_repeats = 1
            repeat_statistics = {}
            self.raw_data["bound_data_list"] = []
            self.raw_data["unbound_data_list"] = []
            self.raw_data["lambda_profiles"] = None
            self.raw_data["time_convergence"] = None
            self.raw_data["bootstrap"] = None

        # Binding free energy = bound - unbound
        delta_g = delta_g_bound - delta_g_unbound
        delta_g_error = np.sqrt(error_bound**2 + error_unbound**2)

        # Step 3: Convergence analysis (if enough data)
        convergence = self._build_convergence_summary()

        # Step 4: Energy decomposition (if components available)
        components = self._build_energy_decomposition()

        n_windows = self._count_lambda_windows(bound_data_list if self.backend == "alchemlyb" else None)

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
                "n_lambda_windows": n_windows,
                "repeat_statistics": repeat_statistics,
            },
            repeat_results=repeat_results,
            n_repeats=n_repeats,
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

        # Generate synthetic results
        delta_g_bound = -2.5  # kcal/mol
        delta_g_unbound = -4.8  # kcal/mol
        delta_g = delta_g_bound - delta_g_unbound
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

    def _build_convergence_summary(self) -> Dict[str, Any]:
        """Build convergence metadata from analysis outputs already stored in raw_data."""
        time_convergence = self.raw_data.get("time_convergence")
        n_samples = self._count_total_samples()
        return {
            "hysteresis": None,
            "time_convergence": time_convergence or {},
            "n_samples": n_samples,
        }

    def _build_energy_decomposition(self) -> Dict[str, float]:
        """Return available energy-component data for the report."""
        return {}

    def _count_total_samples(self) -> int:
        """Count all parsed samples across repeats and lambda windows."""
        bound_data_list = self.raw_data.get("bound_data_list")
        if not isinstance(bound_data_list, list):
            return 0

        sample_count = 0
        for repeat_data in bound_data_list:
            if not isinstance(repeat_data, list):
                continue
            sample_count += sum(len(dataset) for dataset in repeat_data)
        return sample_count

    def _count_lambda_windows(self, bound_data_list: Optional[List[List[Any]]]) -> int:
        """Count lambda windows from parsed data, falling back to directory inspection."""
        if bound_data_list:
            first_repeat = bound_data_list[0]
            if isinstance(first_repeat, list):
                return len(first_repeat)
        return len(list(self.bound_dir.glob("window_*")))

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


@dataclass
class MultiEstimatorResults:
    """
    Container for multi-estimator FEP analysis results

    Stores results from multiple estimators (TI/BAR/MBAR) along with
    comparison metrics to help users understand method agreement and
    choose the best results.

    Attributes
    ----------
    methods : Dict[str, FEResults]
        Results from each estimator, keyed by estimator name.
        Example: {'TI': FEResults, 'BAR': FEResults, 'MBAR': FEResults}
    comparison : Dict[str, Any]
        Comparison metrics across estimators:
        - delta_g_range: Range of ΔG values (max - min)
        - delta_g_std: Standard deviation of ΔG values
        - delta_g_mean: Mean ΔG across methods
        - agreement: 'good' (<0.5), 'moderate' (0.5-1.0), 'poor' (>1.0)
        - best_method: Recommended estimator (usually MBAR)
        - diverged: True if std > 1.0 kcal/mol
    metadata : Dict[str, Any]
        Analysis metadata including:
        - temperature: Simulation temperature
        - estimators_used: List of estimator names that succeeded
        - n_estimators: Number of estimators
        - parallel_execution: Whether estimators were run in parallel
    """

    methods: Dict[str, FEResults] = field(default_factory=dict)
    comparison: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """
        Convert results to dictionary for JSON serialization

        Returns
        -------
        Dict
            Dictionary with 'methods', 'comparison', and 'metadata' keys
        """
        return {
            "methods": {name: results.to_dict() for name, results in self.methods.items()},
            "comparison": self.comparison,
            "metadata": self.metadata,
        }


class FEPMultiEstimatorAnalyzer:
    """
    Run FEP analysis with multiple estimators

    This class enables running multiple FEP estimators (TI/BAR/MBAR) on the
    same data, comparing results, and generating a comprehensive report with
    tab-based switching between methods.

    Key features:
    - Parse XVG files once (shared across all estimators) - major performance optimization
    - Automatic comparison metrics and divergence detection
    - Unified JSON output with all results
    - Optional parallel execution (disabled by default due to pandas pickling limitations)

    Parameters
    ----------
    bound_dir : Union[str, Path, List[Union[str, Path]]]
        Bound leg directory/directories (supports multiple repeats)
    unbound_dir : Union[str, Path, List[Union[str, Path]]]
        Unbound leg directory/directories (supports multiple repeats)
    estimators : list of str, optional
        List of estimator names to run (default: ['TI', 'BAR', 'MBAR'])
    temperature : float, optional
        Simulation temperature in Kelvin (default: 310.0)
    parallel : bool, optional
        Whether to run estimators in parallel (default: False).
        Note: Parallel execution is experimental and may not work on all systems
        due to pandas DataFrame serialization limitations.
    backend : str, optional
        Analysis backend - only 'alchemlyb' supported for multi-estimator mode

    Examples
    --------
    >>> analyzer = FEPMultiEstimatorAnalyzer(
    ...     bound_dir='fep/bound',
    ...     unbound_dir='fep/unbound',
    ...     estimators=['TI', 'BAR', 'MBAR']
    ... )
    >>> multi_results = analyzer.analyze()
    >>> print(f"TI:  {multi_results.methods['TI'].delta_g:.2f}")
    >>> print(f"BAR: {multi_results.methods['BAR'].delta_g:.2f}")
    >>> print(f"MBAR: {multi_results.methods['MBAR'].delta_g:.2f}")

    For CLI usage with --all-estimators flag, see prism.fep.analysis.cli
    """

    ESTIMATORS = {
        "MBAR": MBAR if ALCHEMLYB_AVAILABLE else None,
        "BAR": BAR if ALCHEMLYB_AVAILABLE else None,
        "TI": TI if ALCHEMLYB_AVAILABLE else None,
    }

    def __init__(
        self,
        bound_dir: Union[str, Path, List[Union[str, Path]]],
        unbound_dir: Union[str, Path, List[Union[str, Path]]],
        estimators: Optional[list] = None,
        temperature: float = 310.0,
        parallel: bool = False,
        backend: str = "alchemlyb",
    ):
        """
        Initialize multi-estimator analyzer

        Parameters
        ----------
        bound_dir : str, Path, or list
            Bound leg directory/directories
        unbound_dir : str, Path, or list
            Unbound leg directory/directories
        estimators : list of str, optional
            Estimator names to run (default: ['TI', 'BAR', 'MBAR'])
        temperature : float, optional
            Temperature in Kelvin (default: 310.0)
        parallel : bool, optional
            Run estimators in parallel (default: False).
            Note: Parallel execution is experimental due to pandas serialization limitations.
        backend : str, optional
            Analysis backend (default: 'alchemlyb', only option for multi-estimator)

        Raises
        ------
        ValueError
            If backend is not 'alchemlyb' or estimators are unavailable
        """
        # Support single directory or list of directories (multiple repeats)
        if isinstance(bound_dir, (list, tuple)):
            self.bound_dirs = [Path(d) for d in bound_dir]
        else:
            self.bound_dirs = [Path(bound_dir)]

        if isinstance(unbound_dir, (list, tuple)):
            self.unbound_dirs = [Path(d) for d in unbound_dir]
        else:
            self.unbound_dirs = [Path(unbound_dir)]

        # Set default estimators if not specified
        if estimators is None:
            estimators = ["TI", "BAR", "MBAR"]

        # Validate backend
        if backend != "alchemlyb":
            raise ValueError(f"Multi-estimator mode only supports 'alchemlyb' backend, got {backend}")

        self.estimators = estimators
        self.temperature = temperature
        self.parallel = parallel
        self.backend = backend

        # Validate estimators
        for estimator_name in estimators:
            if estimator_name not in self.ESTIMATORS:
                raise ValueError(f"Unknown estimator: {estimator_name}. Available: {list(self.ESTIMATORS.keys())}")
            if self.ESTIMATORS[estimator_name] is None:
                raise ValueError(f"Estimator {estimator_name} not available. Install alchemlyb.")

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # Results storage
        self.results: Optional[MultiEstimatorResults] = None

    def analyze(self) -> MultiEstimatorResults:
        """
        Run analysis with all estimators

        This method:
        1. Parses XVG files once (shared across all estimators)
        2. Runs each estimator in parallel (or sequentially)
        3. Builds comparison metrics
        4. Returns MultiEstimatorResults

        Returns
        -------
        MultiEstimatorResults
            Results from all estimators + comparison metrics

        Raises
        ------
        FileNotFoundError
            If required input files are missing
        ValueError
            If estimator fails or data format is incorrect
        """
        self.logger.info(f"Starting multi-estimator FEP analysis with {len(self.estimators)} estimators")
        self.logger.info(
            f"Processing {len(self.bound_dirs)} bound repeats and {len(self.unbound_dirs)} unbound repeats"
        )

        # Step 1: Parse XVG files for each estimator (different formats required)
        # Note: TI uses dH/dλ format, while BAR/MBAR use u_nk format
        # We can't truly "parse once, use many" because estimators need different data formats
        self.logger.info("Parsing XVG files for each estimator...")
        results = {}

        for estimator_name in self.estimators:
            self.logger.info(f"Parsing data for {estimator_name} estimator...")

            # Parse data for this specific estimator
            bound_data_list = []
            unbound_data_list = []

            for i, bound_dir in enumerate(self.bound_dirs):
                repeat_label = f"repeat{i + 1}" if len(self.bound_dirs) > 1 else "bound"
                self.logger.debug(f"Parsing {repeat_label} from {bound_dir}")
                bdata = parse_leg_data(bound_dir, "bound", estimator_name, self.temperature, gmx, self.logger)
                bound_data_list.append(bdata)

            for i, unbound_dir in enumerate(self.unbound_dirs):
                repeat_label = f"repeat{i + 1}" if len(self.unbound_dirs) > 1 else "unbound"
                self.logger.debug(f"Parsing {repeat_label} from {unbound_dir}")
                udata = parse_leg_data(unbound_dir, "unbound", estimator_name, self.temperature, gmx, self.logger)
                unbound_data_list.append(udata)

            # Run this estimator
            try:
                estimator_result = self._run_single_estimator(estimator_name, bound_data_list, unbound_data_list)
                results[estimator_name] = estimator_result
                self.logger.info(f"✓ {estimator_name} completed: ΔG = {estimator_result.delta_g:.3f} kcal/mol")
            except Exception as exc:
                self.logger.error(f"✗ {estimator_name} failed: {exc}")
                # Don't re-raise - continue with other estimators

        # Check if any estimators succeeded
        if not results:
            raise RuntimeError(
                f"All {len(self.estimators)} estimators failed. "
                f"Check logs for details. This may indicate incompatible data format."
            )

        # Step 3: Build comparison metrics
        comparison = self._build_comparison_metrics(results)

        # Step 4: Build metadata
        metadata = {
            "temperature": self.temperature,
            "estimators_used": list(results.keys()),
            "n_estimators": len(results),
            "parallel_execution": self.parallel,
            "backend": self.backend,
        }

        # Step 5: Store results
        self.results = MultiEstimatorResults(methods=results, comparison=comparison, metadata=metadata)

        self.logger.info(f"Multi-estimator analysis complete:")
        self.logger.info(f"  Estimators: {', '.join(results.keys())}")
        self.logger.info(f"  ΔG range: {comparison['delta_g_range']:.3f} kcal/mol")
        self.logger.info(f"  Agreement: {comparison['agreement']}")

        return self.results

    def _run_estimators_parallel(self, bound_data_list: list, unbound_data_list: list) -> Dict[str, FEResults]:
        """Run estimators in parallel using ProcessPoolExecutor"""
        from concurrent.futures import ProcessPoolExecutor, as_completed

        results = {}
        self.logger.info(f"Running {len(self.estimators)} estimators in parallel...")

        with ProcessPoolExecutor(max_workers=len(self.estimators)) as executor:
            futures = {}
            for estimator_name in self.estimators:
                future = executor.submit(self._run_single_estimator, estimator_name, bound_data_list, unbound_data_list)
                futures[future] = estimator_name

            for future in as_completed(futures):
                estimator_name = futures[future]
                try:
                    results[estimator_name] = future.result()
                    self.logger.info(
                        f"✓ {estimator_name} completed: ΔG = {results[estimator_name].delta_g:.3f} kcal/mol"
                    )
                except Exception as exc:
                    self.logger.error(f"✗ {estimator_name} failed: {exc}")
                    # Don't re-raise - continue with other estimators

        if not results:
            raise RuntimeError("All estimators failed. Check logs for details.")

        return results

    def _run_estimators_sequential(self, bound_data_list: list, unbound_data_list: list) -> Dict[str, FEResults]:
        """Run estimators sequentially (for debugging or platforms without multiprocessing)"""
        results = {}
        self.logger.info(f"Running {len(self.estimators)} estimators sequentially...")

        for estimator_name in self.estimators:
            try:
                results[estimator_name] = self._run_single_estimator(estimator_name, bound_data_list, unbound_data_list)
                self.logger.info(f"✓ {estimator_name} completed: ΔG = {results[estimator_name].delta_g:.3f} kcal/mol")
            except Exception as exc:
                self.logger.error(f"✗ {estimator_name} failed: {exc}")

        if not results:
            raise RuntimeError("All estimators failed. Check logs for details.")

        return results

    def _run_single_estimator(self, estimator_name: str, bound_data_list: list, unbound_data_list: list) -> FEResults:
        """
        Run a single estimator (called in parallel or sequential)

        Parameters
        ----------
        estimator_name : str
            Name of estimator ('TI', 'BAR', or 'MBAR')
        bound_data_list : list
            List of bound leg datasets (one per repeat)
        unbound_data_list : list
            List of unbound leg datasets (one per repeat)

        Returns
        -------
        FEResults
            Analysis results from this estimator
        """
        from .estimators import compute_free_energy_alchemlyb, summarize_repeat_results

        estimator_cls = self.ESTIMATORS[estimator_name]

        # Process all repeats for bound leg
        bound_results = []
        unbound_results = []
        fitted_bounds = []
        fitted_unbounds = []

        for bdata in bound_data_list:
            dg, err, fitted, overlap_m = compute_free_energy_alchemlyb(estimator_cls, bdata, self.logger)
            bound_results.append((dg, err))
            fitted_bounds.append(fitted)

        for udata in unbound_data_list:
            dg, err, fitted, overlap_m = compute_free_energy_alchemlyb(estimator_cls, udata, self.logger)
            unbound_results.append((dg, err))
            fitted_unbounds.append(fitted)

        # Summarize repeat results
        summary = summarize_repeat_results(bound_results, unbound_results)

        # Build lambda profiles for plotting
        from .profiles import build_lambda_profiles

        lambda_profiles = build_lambda_profiles(fitted_bounds, fitted_unbounds, estimator_name)

        # Calculate binding free energy
        delta_g = summary["delta_g_bound"] - summary["delta_g_unbound"]
        delta_g_error = np.sqrt(summary["error_bound"] ** 2 + summary["error_unbound"] ** 2)

        # Count lambda windows
        n_windows = len(bound_data_list[0]) if bound_data_list else 0

        # Build results
        results = FEResults(
            delta_g=delta_g,
            delta_g_error=delta_g_error,
            delta_g_bound=summary["delta_g_bound"],
            delta_g_unbound=summary["delta_g_unbound"],
            delta_g_components={},  # Not implemented yet
            convergence={},  # Not implemented yet
            metadata={
                "temperature": self.temperature,
                "estimator": estimator_name,
                "backend": self.backend,
                "n_lambda_windows": n_windows,
                "repeat_statistics": summary["repeat_statistics"],
            },
            repeat_results=summary["repeat_results"],
            n_repeats=summary["n_repeats"],
            lambda_profiles=lambda_profiles,  # Add lambda profiles for plotting
        )

        return results

    def _build_comparison_metrics(self, results: Dict[str, FEResults]) -> Dict[str, Any]:
        """
        Build comparison metrics across estimators

        Parameters
        ----------
        results : Dict[str, FEResults]
            Results from each estimator

        Returns
        -------
        Dict[str, Any]
            Comparison metrics including range, std, agreement, etc.
        """
        if not results:
            return {}

        delta_g_values = [r.delta_g for r in results.values()]

        import numpy as np

        delta_g_range = float(max(delta_g_values) - min(delta_g_values))
        delta_g_std = float(np.std(delta_g_values))
        delta_g_mean = float(np.mean(delta_g_values))

        # Determine agreement level
        if delta_g_std < 0.5:
            agreement = "good"
        elif delta_g_std < 1.0:
            agreement = "moderate"
        else:
            agreement = "poor"

        # Check if results diverged significantly
        diverged = delta_g_std > 1.0

        return {
            "delta_g_range": delta_g_range,
            "delta_g_std": delta_g_std,
            "delta_g_mean": delta_g_mean,
            "agreement": agreement,
            "best_method": "MBAR",  # Could be made dynamic based on overlap matrix
            "diverged": diverged,
        }

    def generate_html_report(self, output_path: Union[str, Path]) -> Path:
        """
        Generate multi-estimator comparison HTML report

        Parameters
        ----------
        output_path : str or Path
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

        from prism.fep.analysis.report import MultiEstimatorReportGenerator

        generator = MultiEstimatorReportGenerator(self.results)
        html_path = generator.generate(output_path)

        self.logger.info(f"Multi-estimator HTML report generated: {html_path}")
        return html_path

    def save_results(self, output_path: Union[str, Path]) -> Path:
        """
        Save multi-estimator results to JSON file

        Parameters
        ----------
        output_path : str or Path
            Path to output JSON file

        Returns
        -------
        Path
            Path to saved JSON file

        Raises
        ------
        RuntimeError
            If analyze() has not been called first
        """
        if self.results is None:
            raise RuntimeError("Must call analyze() before saving results")

        output_path = Path(output_path)
        with open(output_path, "w") as f:
            json.dump(self.results.to_dict(), f, indent=2)

        self.logger.info(f"Multi-estimator results saved to {output_path}")
        return output_path
