#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Multi-Estimator FEP Analyzer

Analyzes FEP calculations using multiple estimators (TI, BAR, MBAR)
and generates comparison metrics.

Example Usage:
-------------
    from prism.fep.analysis.analyzers import FEPMultiEstimatorAnalyzer

    multi = FEPMultiEstimatorAnalyzer(
        bound_dirs=['bound1', 'bound2'],
        unbound_dirs=['unbound1', 'unbound2']
    )

    results = multi.analyze()
    print(f"ΔG range: {results.comparison['delta_g_range']:.3f} kcal/mol")
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from ..core.models import FEResults, MultiEstimatorResults
from ..core.convergence import compute_bootstrap, compute_time_convergence
from ..core.xvg_parser import parse_leg_data

try:
    from alchemlyb.estimators import MBAR, BAR, TI

    ALCHEMLYB_AVAILABLE = True
except ImportError:
    ALCHEMLYB_AVAILABLE = False
    MBAR = None
    BAR = None
    TI = None
    logging.warning("alchemlyb not available. Install with: pip install alchemlyb")


class FEPMultiEstimatorAnalyzer:
    """
    Multi-Estimator FEP Analyzer

    Runs FEP analysis with multiple estimators (TI, BAR, MBAR) in parallel
    and compares results to assess convergence and method agreement.

    Features:
    ---------
    - Parse XVG files once (shared across all estimators) - major performance optimization
    - Automatic comparison metrics and divergence detection
    - Unified JSON output with all results
    - Progress bars for long-running operations

    Parameters
    ----------
    bound_dirs : list of str or Path
        List of bound leg directories (one per repeat)
    unbound_dirs : list of str or Path
        List of unbound leg directories (one per repeat)
    estimators : list of str, optional
        List of estimator names to run (default: ['TI', 'BAR', 'MBAR'])
    temperature : float, optional
        Simulation temperature in Kelvin (default: 310.0)
    backend : str, optional
        Analysis backend - only 'alchemlyb' supported for multi-estimator mode
    show_progress : bool, optional
        Show progress bars during analysis (default: True)

    Examples
    --------
    >>> multi = FEPMultiEstimatorAnalyzer(
    ...     bound_dirs=['fep/bound_rep1', 'fep/bound_rep2'],
    ...     unbound_dirs=['fep/unbound_rep1', 'fep/unbound_rep2'],
    ...     estimators=['TI', 'BAR', 'MBAR']
    ... )
    >>> results = multi.analyze()
    >>> print(f"Best estimator: {results.comparison['best_method']}")
    """

    # Estimator mapping
    ESTIMATORS = {
        "TI": TI if ALCHEMLYB_AVAILABLE else None,
        "BAR": BAR if ALCHEMLYB_AVAILABLE else None,
        "MBAR": MBAR if ALCHEMLYB_AVAILABLE else None,
    }

    def __init__(
        self,
        bound_dirs: Union[str, Path, List[Union[str, Path]]],
        unbound_dirs: Union[str, Path, List[Union[str, Path]]],
        estimators: Optional[list] = None,
        temperature: float = 310.0,
        backend: str = "alchemlyb",
        show_progress: bool = True,
        skip_bootstrap: bool = False,
        skip_time_convergence: bool = False,
        cache_file: Optional[Union[str, Path]] = None,
    ):
        """
        Initialize multi-estimator analyzer

        Parameters
        ----------
        bound_dirs : str or Path or list
            Directory path(s) for bound leg (supports multiple repeats)
        unbound_dirs : str or Path or list
            Directory path(s) for unbound leg (supports multiple repeats)
        estimators : list, optional
            Estimator names to run (default: ['TI', 'BAR', 'MBAR'])
        temperature : float, optional
            Temperature in Kelvin (default: 310.0)
        backend : str, optional
            Analysis backend (default: 'alchemlyb', only option for multi-estimator)
        show_progress : bool, optional
            Show progress bars (default: True)
        skip_bootstrap : bool, optional
            Skip bootstrap error estimation (saves ~30-60 min, default: False)
        skip_time_convergence : bool, optional
            Skip time convergence analysis (saves ~10-15 min, default: False)
        cache_file : str or Path, optional
            Path to cache file for saving/loading analysis results.
            If provided and file exists, results will be loaded from cache
            instead of re-running analysis. After analysis, results will be
            saved to this file for future use.

        Raises
        ------
        ValueError
            If invalid estimators specified or backend not supported
        ImportError
            If alchemlyb not available
        """
        # Normalize directory paths
        self.bound_dirs = [Path(bound_dirs)] if not isinstance(bound_dirs, list) else [Path(d) for d in bound_dirs]
        self.unbound_dirs = (
            [Path(unbound_dirs)] if not isinstance(unbound_dirs, list) else [Path(d) for d in unbound_dirs]
        )

        # Validate inputs
        if len(self.bound_dirs) != len(self.unbound_dirs):
            raise ValueError("Number of bound_dirs must match number of unbound_dirs")

        # Set defaults
        self.estimators = estimators or ["TI", "BAR", "MBAR"]
        self.temperature = temperature
        self.backend = backend
        self.show_progress = show_progress
        self.skip_bootstrap = skip_bootstrap
        self.skip_time_convergence = skip_time_convergence
        self.cache_file = Path(cache_file) if cache_file else None

        # Validate estimators
        for estimator_name in self.estimators:
            if estimator_name not in self.ESTIMATORS:
                raise ValueError(f"Unknown estimator: {estimator_name}. Choose from {list(self.ESTIMATORS.keys())}")
            if self.ESTIMATORS[estimator_name] is None:
                raise ImportError(f"Estimator {estimator_name} requires alchemlyb. Install with: pip install alchemlyb")

        # Validate backend
        if self.backend != "alchemlyb":
            raise ValueError(f"Backend '{self.backend}' not supported. Use 'alchemlyb' for multi-estimator mode.")

        # Setup logger
        self.logger = logging.getLogger(__name__)

        # Storage
        self.results: Optional[MultiEstimatorResults] = None

        # Import gmx for XVG parsing
        global gmx
        from alchemlyb.parsing import gmx

    def analyze(self) -> MultiEstimatorResults:
        """
        Run analysis with all estimators

        This method:
        1. Checks for cached results and loads if available
        2. Parses XVG files for each estimator with progress bars
        3. Runs each estimator sequentially
        4. Builds comparison metrics
        5. Saves results to cache file if specified
        6. Returns MultiEstimatorResults

        Returns
        -------
        MultiEstimatorResults
            Results from all estimators with comparison metrics
        """
        # Check cache first
        if self.cache_file and self.cache_file.exists():
            self.logger.info(f"Loading cached results from {self.cache_file}")
            try:
                import json

                with open(self.cache_file, "r") as f:
                    cache_data = json.load(f)

                # Verify cache matches current configuration
                cache_metadata = cache_data.get("metadata", {})
                cache_estimators = cache_metadata.get("estimators_used", [])
                cache_temp = cache_metadata.get("temperature")

                if set(cache_estimators) == set(self.estimators) and cache_temp == self.temperature:
                    self.logger.info("Cache matches current configuration - loading results")
                    self.results = MultiEstimatorResults.from_dict(cache_data)
                    return self.results
                else:
                    self.logger.warning(
                        f"Cache mismatch: cache has {cache_estimators} @ {cache_temp}K, "
                        f"current has {self.estimators} @ {self.temperature}K. Re-running analysis."
                    )
            except Exception as e:
                self.logger.warning(f"Failed to load cache: {e}. Re-running analysis.")

        # Step 1: Parse XVG files...
        self.logger.info("Parsing XVG files...")

        if self.show_progress:
            try:
                from tqdm import tqdm

                estimator_iter = tqdm(self.estimators, desc="Parsing XVG files", unit="estimator")
            except ImportError:
                estimator_iter = self.estimators
        else:
            estimator_iter = self.estimators

        # Store parsed data for each estimator
        parsed_data = {}

        for estimator_name in estimator_iter:
            if self.show_progress and hasattr(estimator_iter, "set_description"):
                estimator_iter.set_description(f"  {estimator_name}")

            # Parse data for this specific estimator
            bound_data_list = []
            unbound_data_list = []

            # Parse bound legs
            for i, bound_dir in enumerate(self.bound_dirs):
                repeat_label = f"repeat{i + 1}" if len(self.bound_dirs) > 1 else "bound"
                if self.show_progress and hasattr(estimator_iter, "write"):
                    estimator_iter.write(f"  Parsing {repeat_label} from {bound_dir}")
                else:
                    self.logger.debug(f"Parsing {repeat_label} from {bound_dir}")
                bdata = parse_leg_data(bound_dir, "bound", estimator_name, self.temperature, gmx, self.logger)
                bound_data_list.append(bdata)

            # Parse unbound legs
            for i, unbound_dir in enumerate(self.unbound_dirs):
                repeat_label = f"repeat{i + 1}" if len(self.unbound_dirs) > 1 else "unbound"
                if self.show_progress and hasattr(estimator_iter, "write"):
                    estimator_iter.write(f"  Parsing {repeat_label} from {unbound_dir}")
                else:
                    self.logger.debug(f"Parsing {repeat_label} from {unbound_dir}")
                udata = parse_leg_data(unbound_dir, "unbound", estimator_name, self.temperature, gmx, self.logger)
                unbound_data_list.append(udata)

            parsed_data[estimator_name] = (bound_data_list, unbound_data_list)

        # Step 2: Run all estimators
        results = {}
        for estimator_name, (bound_data_list, unbound_data_list) in parsed_data.items():
            try:
                if self.show_progress:
                    self.logger.info(f"Running {estimator_name} estimator...")

                estimator_result = self._run_single_estimator(estimator_name, bound_data_list, unbound_data_list)
                results[estimator_name] = estimator_result

                if self.show_progress:
                    self.logger.info(f"✓ {estimator_name} completed: ΔG = {estimator_result.delta_g:.3f} kcal/mol")
            except Exception as exc:
                self.logger.error(f"✗ {estimator_name} failed: {exc}")
                # Don't re-raise - continue with other estimators

        if not results:
            raise RuntimeError(
                "All estimators failed. Check your input data and ensure XVG files are correctly formatted. "
                f"Check logs for details. This may indicate incompatible data format."
            )

        # Step 3: Compute convergence metrics for each estimator
        if self.show_progress:
            self.logger.info("Computing convergence metrics...")

        for estimator_name, estimator_result in results.items():
            bound_data_list, unbound_data_list = parsed_data[estimator_name]
            estimator_class = self.ESTIMATORS[estimator_name]

            # Time convergence analysis
            if self.skip_time_convergence:
                if self.show_progress:
                    self.logger.info(f"  Skipping time convergence analysis for {estimator_name}")
                estimator_result.time_convergence = None
            else:
                try:
                    estimator_result.time_convergence = compute_time_convergence(
                        bound_data_list,
                        unbound_data_list,
                        estimator_class,
                        n_points=10,
                    )
                except Exception as exc:
                    self.logger.warning(f"Time convergence analysis failed for {estimator_name}: {exc}")
                    estimator_result.time_convergence = None

            # Bootstrap analysis
            if self.skip_bootstrap:
                if self.show_progress:
                    self.logger.info(f"  Skipping bootstrap analysis for {estimator_name}")
                estimator_result.bootstrap = None
            else:
                if self.show_progress:
                    self.logger.info(f"  Running bootstrap analysis for {estimator_name}...")
                try:
                    estimator_result.bootstrap = compute_bootstrap(
                        bound_data_list,
                        unbound_data_list,
                        estimator_class,
                        n_bootstrap=50,
                        fraction=0.8,
                    )
                except Exception as exc:
                    self.logger.warning(f"Bootstrap analysis failed for {estimator_name}: {exc}")
                    estimator_result.bootstrap = None

        # Step 4: Build comparison metrics
        comparison = self._build_comparison_metrics(results)

        # Step 5: Build metadata
        metadata = {
            "temperature": self.temperature,
            "estimators_used": list(results.keys()),
            "n_estimators": len(results),
            "backend": self.backend,
        }

        # Step 6: Store results
        self.results = MultiEstimatorResults(methods=results, comparison=comparison, metadata=metadata)

        self.logger.info(f"Multi-estimator analysis complete:")
        self.logger.info(f"  Estimators: {', '.join(results.keys())}")
        self.logger.info(f"  ΔG range: {comparison['delta_g_range']:.3f} kcal/mol")
        self.logger.info(f"  Agreement: {comparison['agreement']}")

        # Save to cache if specified
        if self.cache_file:
            self.logger.info(f"Saving results to cache: {self.cache_file}")
            try:
                import json

                self.cache_file.parent.mkdir(parents=True, exist_ok=True)
                with open(self.cache_file, "w") as f:
                    json.dump(self.results.to_dict(), f, indent=2)
                self.logger.info("Cache saved successfully")
            except Exception as e:
                self.logger.warning(f"Failed to save cache: {e}")

        return self.results

    def _run_single_estimator(self, estimator_name: str, bound_data_list: list, unbound_data_list: list) -> FEResults:
        """
        Run a single estimator

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
        from ..core.estimators import compute_free_energy_alchemlyb, summarize_repeat_results
        from ..core.profiles import build_lambda_profiles

        estimator_cls = self.ESTIMATORS[estimator_name]

        # Process all repeats for bound leg
        bound_results = []
        unbound_results = []
        fitted_bounds = []
        fitted_unbounds = []

        # Fit bound leg data
        for bdata in bound_data_list:
            dg, err, fitted, overlap_m = compute_free_energy_alchemlyb(estimator_cls, bdata, self.logger)
            bound_results.append((dg, err))
            fitted_bounds.append(fitted)

        # Fit unbound leg data
        for udata in unbound_data_list:
            dg, err, fitted, overlap_m = compute_free_energy_alchemlyb(estimator_cls, udata, self.logger)
            unbound_results.append((dg, err))
            fitted_unbounds.append(fitted)

        # Summarize repeat results
        summary = summarize_repeat_results(bound_results, unbound_results)

        # Build lambda profiles for plotting
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
            convergence={},
            metadata={
                "temperature": self.temperature,
                "estimator": estimator_name,
                "backend": self.backend,
                "n_lambda_windows": n_windows,
                "repeat_statistics": summary["repeat_statistics"],
            },
            repeat_results=summary["repeat_results"],
            n_repeats=summary["n_repeats"],
            lambda_profiles=lambda_profiles,
        )

        return results

    def _build_comparison_metrics(self, results: Dict[str, FEResults]) -> Dict[str, Any]:
        """
        Build comparison metrics between estimators

        Parameters
        ----------
        results : dict
            Dictionary of estimator name to FEResults

        Returns
        -------
        dict
            Comparison metrics including range, std, mean, agreement, etc.
        """
        # Extract delta G values
        delta_g_values = [r.delta_g for r in results.values()]

        # Compute statistics
        delta_g_range = max(delta_g_values) - min(delta_g_values)
        delta_g_std = np.std(delta_g_values)
        delta_g_mean = np.mean(delta_g_values)

        # Determine agreement level
        if delta_g_std < 0.5:
            agreement = "good"
        elif delta_g_std < 1.0:
            agreement = "moderate"
        else:
            agreement = "poor"

        # Check for divergence
        diverged = bool(delta_g_std > 1.0)

        # Select best method (prefer MBAR if available)
        if "MBAR" in results:
            best_method = "MBAR"
        elif "BAR" in results:
            best_method = "BAR"
        else:
            best_method = list(results.keys())[0]

        return {
            "delta_g_range": delta_g_range,
            "delta_g_std": delta_g_std,
            "delta_g_mean": delta_g_mean,
            "agreement": agreement,
            "diverged": diverged,
            "best_method": best_method,
            "n_estimators": len(results),
        }
