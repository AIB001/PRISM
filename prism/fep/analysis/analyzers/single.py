#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single-Estimator FEP Analyzer

Analyzes FEP calculations using a single estimator (TI, BAR, or MBAR).

Example Usage:
-------------
    from prism.fep.analysis.analyzers import FEPAnalyzer

    analyzer = FEPAnalyzer(
        bound_dir='fep_project/bound',
        unbound_dir='fep_project/unbound',
        temperature=310.0,
        estimator='MBAR'
    )

    results = analyzer.analyze()
    print(f"ΔG = {results.delta_g:.2f} ± {results.delta_g_error:.2f} kcal/mol")
"""

import json
import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Union


from ..core.models import FEResults
from ..core.convergence import compute_bootstrap, compute_time_convergence
from ..core.estimators import compute_free_energy_alchemlyb, compute_free_energy_gmx_bar, summarize_repeat_results
from ..core.profiles import build_lambda_profiles
from ..core.xvg_parser import parse_leg_data

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
        skip_bootstrap: bool = False,
        skip_time_convergence: bool = False,
    ):
        """
        Initialize FEP analyzer

        Parameters
        ----------
        bound_dir : str or Path or list
            Directory path(s) for bound leg (supports multiple repeats)
        unbound_dir : str or Path or list
            Directory path(s) for unbound leg (supports multiple repeats)
        temperature : float
            Simulation temperature in Kelvin
        estimator : str
            Estimator type ('TI', 'BAR', or 'MBAR')
        backend : str
            Analysis backend ('alchemlyb' or 'gmx_bar')
        energy_components : list, optional
            Energy components to analyze
        allow_mock : bool
            Allow mock results for testing (default: False)
        skip_bootstrap : bool
            Skip bootstrap error estimation (saves ~10-30 min, default: False)
        skip_time_convergence : bool
            Skip time convergence analysis (saves ~5-10 min, default: False)
        """
        self.bound_dirs = [Path(bound_dir)] if not isinstance(bound_dir, list) else [Path(d) for d in bound_dir]
        self.unbound_dirs = [Path(unbound_dir)] if not isinstance(unbound_dir, list) else [Path(d) for d in unbound_dir]
        self.temperature = temperature
        self.estimator = estimator.upper()
        self.backend = backend.lower()
        self.energy_components = energy_components or ["elec", "vdw"]
        self.allow_mock = allow_mock
        self.skip_bootstrap = skip_bootstrap
        self.skip_time_convergence = skip_time_convergence

        if len(self.bound_dirs) != len(self.unbound_dirs):
            raise ValueError("Number of bound and unbound directories must match")

        # Validate estimator
        if self.estimator not in self.ESTIMATORS:
            raise ValueError(f"Unknown estimator: {self.estimator}. Choose from {list(self.ESTIMATORS.keys())}")

        if self.ESTIMATORS[self.estimator] is None and not allow_mock:
            raise ImportError(f"Estimator {self.estimator} requires alchemlyb. Install with: pip install alchemlyb")

        # Validate backend
        if self.backend == "gmx_bar" and self.estimator != "BAR":
            raise ValueError("gmx_bar backend only supports BAR estimator")

        # Setup logger
        self.logger = logging.getLogger(__name__)

        # Storage for raw data and results
        self.raw_data: Dict[str, Any] = {}
        self.results: Optional[FEResults] = None

    def analyze(self) -> FEResults:
        """
        Run complete FEP analysis workflow

        Returns
        -------
        FEResults
            Analysis results including free energies, errors, and convergence metrics
        """
        # Check for mock mode
        if self.allow_mock and not ALCHEMLYB_AVAILABLE:
            self.logger.warning("Using mock results (alchemlyb not available)")
            return self._analyze_mock()

        try:
            estimator_cls = self.ESTIMATORS[self.estimator]

            # Parse data files
            bound_data_list = []
            unbound_data_list = []

            for i, bound_dir in enumerate(self.bound_dirs):
                repeat_label = f"repeat{i + 1}" if len(self.bound_dirs) > 1 else "bound"
                self.logger.debug(f"Parsing {repeat_label} from {bound_dir}")
                bdata = parse_leg_data(bound_dir, "bound", self.estimator, self.temperature, gmx, self.logger)
                bound_data_list.append(bdata)

            for i, unbound_dir in enumerate(self.unbound_dirs):
                repeat_label = f"repeat{i + 1}" if len(self.unbound_dirs) > 1 else "unbound"
                self.logger.debug(f"Parsing {repeat_label} from {unbound_dir}")
                udata = parse_leg_data(unbound_dir, "unbound", self.estimator, self.temperature, gmx, self.logger)
                unbound_data_list.append(udata)

            # Compute free energies per repeat
            bound_results = []
            unbound_results = []
            fitted_bounds = []
            fitted_unbounds = []

            if self.backend == "gmx_bar":
                for bound_dir in self.bound_dirs:
                    bound_results.append(compute_free_energy_gmx_bar(bound_dir, "bound", self.temperature, self.logger))
                for unbound_dir in self.unbound_dirs:
                    unbound_results.append(
                        compute_free_energy_gmx_bar(unbound_dir, "unbound", self.temperature, self.logger)
                    )
            else:
                for bdata in bound_data_list:
                    dg, err, fitted, _overlap = compute_free_energy_alchemlyb(estimator_cls, bdata, self.logger)
                    bound_results.append((dg, err))
                    fitted_bounds.append(fitted)
                for udata in unbound_data_list:
                    dg, err, fitted, _overlap = compute_free_energy_alchemlyb(estimator_cls, udata, self.logger)
                    unbound_results.append((dg, err))
                    fitted_unbounds.append(fitted)

            # Build repeat statistics and final binding free energy
            summary = summarize_repeat_results(bound_results, unbound_results)
            delta_g_bound = summary["delta_g_bound"]
            delta_g_unbound = summary["delta_g_unbound"]
            delta_g_error = math.sqrt(summary["error_bound"] ** 2 + summary["error_unbound"] ** 2)
            delta_g = delta_g_bound - delta_g_unbound

            self.logger.info(f"FEP Analysis Results ({self.estimator}):")
            self.logger.info(f"  ΔG_bound   = {delta_g_bound:.3f} kcal/mol")
            self.logger.info(f"  ΔG_unbound = {delta_g_unbound:.3f} kcal/mol")
            self.logger.info(f"  ΔΔG_bind   = {delta_g:.3f} ± {delta_g_error:.3f} kcal/mol")

            # Build lambda profiles for plotting
            lambda_profiles = (
                None
                if self.backend == "gmx_bar"
                else build_lambda_profiles(fitted_bounds, fitted_unbounds, self.estimator)
            )

            # Compute convergence metrics (optional)
            if self.backend == "gmx_bar" or self.skip_time_convergence:
                self.logger.info("Skipping time convergence analysis")
                self.raw_data["time_convergence"] = None
            else:
                self.raw_data["time_convergence"] = compute_time_convergence(
                    bound_data_list, unbound_data_list, estimator_cls, n_points=10
                )

            if self.backend == "gmx_bar" or self.skip_bootstrap:
                self.logger.info("Skipping bootstrap analysis")
                self.raw_data["bootstrap"] = None
            else:
                self.logger.info("Running bootstrap analysis...")
                self.raw_data["bootstrap"] = compute_bootstrap(
                    bound_data_list,
                    unbound_data_list,
                    estimator_cls,
                    n_bootstrap=50,
                    fraction=0.8,
                )

            self.raw_data["lambda_profiles"] = lambda_profiles

            # Store results
            self.results = FEResults(
                delta_g=delta_g,
                delta_g_error=delta_g_error,
                delta_g_bound=summary["delta_g_bound"],
                delta_g_unbound=summary["delta_g_unbound"],
                delta_g_components=self._build_energy_decomposition(),
                convergence=self._build_convergence_summary(),
                metadata={
                    "temperature": self.temperature,
                    "estimator": self.estimator,
                    "backend": self.backend,
                    "n_lambda_windows": self._count_lambda_windows(bound_data_list),
                    "repeat_statistics": summary.get("repeat_statistics", {}),
                },
                repeat_results=summary["repeat_results"],
                n_repeats=summary["n_repeats"],
                lambda_profiles=lambda_profiles,
                time_convergence=self.raw_data["time_convergence"],
                bootstrap=self.raw_data["bootstrap"],
            )

            return self.results
        except Exception as exc:
            if self.allow_mock:
                self.logger.warning(f"Analysis failed; falling back to mock results: {exc}")
                return self._analyze_mock()
            raise

    def _analyze_mock(self) -> FEResults:
        """Generate mock results for testing"""
        self.raw_data = {
            "lambda_profiles": None,
            "time_convergence": None,
            "bootstrap": None,
        }
        self.results = FEResults(
            delta_g=4.23,
            delta_g_error=0.15,
            delta_g_bound=54.07,
            delta_g_unbound=49.84,
            delta_g_components={"elec": 2.1, "vdw": 2.13},
            convergence={"status": "mock"},
            metadata={
                "temperature": self.temperature,
                "estimator": self.estimator,
                "backend": "mock",
                "n_lambda_windows": 11,
            },
            repeat_results=[],
            n_repeats=1,
            lambda_profiles=None,
            time_convergence=None,
            bootstrap=None,
        )
        return self.results

    def _build_convergence_summary(self) -> Dict[str, Any]:
        """Build convergence summary dictionary"""
        summary = {
            "time_convergence_performed": self.raw_data.get("time_convergence") is not None,
            "bootstrap_performed": self.raw_data.get("bootstrap") is not None,
        }

        if self.raw_data.get("time_convergence"):
            summary["time_convergence"] = {
                "n_points": self.raw_data["time_convergence"].get("n_points", 0),
                "converged": self.raw_data["time_convergence"].get("converged", False),
            }

        if self.raw_data.get("bootstrap"):
            summary["bootstrap"] = {
                "n_bootstrap": self.raw_data["bootstrap"].get("n_bootstrap", 0),
                "mean_delta_g": self.raw_data["bootstrap"].get("mean_delta_g", 0.0),
                "stderr_delta_g": self.raw_data["bootstrap"].get("stderr_delta_g", 0.0),
            }

        return summary

    def _build_energy_decomposition(self) -> Dict[str, float]:
        """Build energy decomposition by component"""
        # This would require additional data from XVG files
        # For now, return empty dict
        return {}

    def _count_total_samples(self) -> int:
        """Count total number of samples across all lambda windows"""
        if not self.results:
            return 0

        total = 0
        if self.results.lambda_profiles:
            for profile in self.results.lambda_profiles.get("bound", []):
                total += len(profile.get("energies", []))
            for profile in self.results.lambda_profiles.get("unbound", []):
                total += len(profile.get("energies", []))

        return total

    def _count_lambda_windows(self, bound_data_list: Optional[List[List[Any]]]) -> int:
        """Count number of lambda windows"""
        if not bound_data_list or not bound_data_list[0]:
            return 0

        # Count from first repeat's data
        return len(bound_data_list[0])

    def generate_html_report(self, output_path: Union[str, Path]) -> Path:
        """
        Generate HTML report

        Parameters
        ----------
        output_path : str or Path
            Output HTML file path

        Returns
        -------
        Path
            Path to generated report
        """
        from ..reports import HTMLReportGenerator

        if self.results is None:
            raise RuntimeError("No results to generate report. Run analyze() first.")

        generator = HTMLReportGenerator(self.results, raw_data=self.raw_data)
        generator.generate(output_path)

        self.logger.info(f"HTML report saved to: {output_path}")
        return Path(output_path)

    def save_results(self, output_path: Union[str, Path]) -> Path:
        """
        Save analysis results to JSON file

        Parameters
        ----------
        output_path : str or Path
            Output JSON file path

        Returns
        -------
        Path
            Path to saved results
        """
        if self.results is None:
            raise RuntimeError("No results to save. Run analyze() first.")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            json.dump(self.results.to_dict(), f, indent=2)

        self.logger.info(f"Results saved to: {output_path}")
        return output_path
