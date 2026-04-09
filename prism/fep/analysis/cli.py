#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Command-line interface for FEP analysis

Provides CLI for analyzing FEP calculations from GROMACS output files.

Usage:
------
    prism --fep-analyze \\
        --bound-dir fep_project/bound \\
        --unbound-dir fep_project/unbound \\
        --output report.html \\
        --estimator MBAR \\
        --temperature 310

Author: PRISM Team
"""

import argparse
import logging
import sys
from pathlib import Path

from prism.fep.analysis import FEPAnalyzer, FEPMultiEstimatorAnalyzer, MultiEstimatorReportGenerator


def setup_logging(verbose: bool = False) -> None:
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def parse_arguments(argv=None) -> argparse.Namespace:
    """
    Parse command-line arguments

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog="prism --fep-analyze",
        description="Analyze FEP calculations from GROMACS output files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
--------
  # Basic analysis with MBAR estimator
  prism --fep-analyze \\
      --bound-dir fep_project/bound \\
      --unbound-dir fep_project/unbound \\
      --output report.html

  # Use BAR estimator with custom temperature
  prism --fep-analyze \\
      --bound-dir fep_project/bound \\
      --unbound-dir fep_project/unbound \\
      --estimator BAR \\
      --temperature 300 \\
      --output report.html

  # Run multiple estimators and generate comparison report
  prism --fep-analyze \\
      --bound-dir fep_project/bound \\
      --unbound-dir fep_project/unbound \\
      --estimator TI BAR MBAR \\
      --output report.html

  # Run all estimators (convenience flag)
  prism --fep-analyze \\
      --bound-dir fep_project/bound \\
      --unbound-dir fep_project/unbound \\
      --all-estimators \\
      --output report.html

  # Save results to JSON as well
  prism --fep-analyze \\
      --bound-dir fep_project/bound \\
      --unbound-dir fep_project/unbound \\
      --output report.html \\
      --json results.json

For more information, see: https://github.com/your-repo/prism
        """,
    )

    # Required arguments
    required = parser.add_argument_group("Required Arguments")
    required.add_argument(
        "--bound-dir",
        type=str,
        required=True,
        metavar="PATH",
        help="Path to bound leg directory containing lambda window subdirectories (window_00, window_01, ...)",
    )
    required.add_argument(
        "--unbound-dir",
        type=str,
        required=True,
        metavar="PATH",
        help="Path to unbound leg directory containing lambda window subdirectories (window_00, window_01, ...)",
    )

    # Output options
    output = parser.add_argument_group("Output Options")
    output.add_argument(
        "--output",
        "-o",
        type=str,
        default="fep_analysis_report.html",
        metavar="FILE",
        help="Output HTML report file (default: fep_analysis_report.html)",
    )
    output.add_argument(
        "--json", type=str, default=None, metavar="FILE", help="Also save results to JSON file (optional)"
    )

    # Analysis options
    analysis = parser.add_argument_group("Analysis Options")
    analysis.add_argument(
        "--estimator",
        "-e",
        nargs="+",
        choices=["MBAR", "BAR", "TI"],
        default=["MBAR"],
        help="Free energy estimator method(s). Can specify multiple: --estimator TI BAR MBAR (default: MBAR)",
    )
    analysis.add_argument(
        "--all-estimators",
        action="store_true",
        help="Run all estimators (TI, BAR, MBAR) and generate comparison report",
    )
    analysis.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars (useful for logging to file)",
    )
    analysis.add_argument(
        "--backend",
        type=str,
        choices=["alchemlyb", "gmx_bar"],
        default="alchemlyb",
        help="Backend for analysis (default: alchemlyb). 'alchemlyb' supports TI/BAR/MBAR; 'gmx_bar' only supports BAR",
    )
    analysis.add_argument(
        "--temperature",
        type=float,
        default=310.0,
        metavar="TEMP",
        help="Simulation temperature in Kelvin (default: 310.0)",
    )
    analysis.add_argument(
        "--components",
        type=str,
        nargs="+",
        default=["elec", "vdw"],
        metavar="COMP",
        help="Energy components to analyze (default: elec vdw)",
    )
    analysis.add_argument(
        "--bootstrap-n-jobs",
        type=int,
        default=1,
        metavar="N",
        help="Number of parallel workers for bootstrap analysis (default: 1 = serial)",
    )

    # Utility options
    utility = parser.add_argument_group("Utility Options")
    utility.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    utility.add_argument("--version", action="version", version="PRISM-FEP Analysis 0.1.0")

    return parser.parse_args(argv)


def validate_arguments(args: argparse.Namespace) -> None:
    """
    Validate command-line arguments

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments

    Raises
    ------
    FileNotFoundError
        If input directories don't exist
    ValueError
        If arguments are invalid
    """
    # Check bound directory
    bound_dir = Path(args.bound_dir)
    if not bound_dir.exists():
        raise FileNotFoundError(f"Bound directory not found: {args.bound_dir}")

    # Check for lambda windows
    window_dirs = list(bound_dir.glob("window_*"))
    if not window_dirs:
        raise ValueError(
            f"No lambda window directories found in {args.bound_dir}. Expected format: window_00, window_01, ..."
        )

    # Check unbound directory
    unbound_dir = Path(args.unbound_dir)
    if not unbound_dir.exists():
        raise FileNotFoundError(f"Unbound directory not found: {args.unbound_dir}")

    # Check for lambda windows
    window_dirs = list(unbound_dir.glob("window_*"))
    if not window_dirs:
        raise ValueError(
            f"No lambda window directories found in {args.unbound_dir}. Expected format: window_00, window_01, ..."
        )

    # Validate temperature
    if args.temperature <= 0:
        raise ValueError(f"Temperature must be positive: {args.temperature}")

    # Validate multi-estimator mode backend
    if args.all_estimators or len(args.estimator) > 1:
        if args.backend != "alchemlyb":
            raise ValueError(
                f"Multi-estimator mode only supports 'alchemlyb' backend. "
                f"Got: {args.backend}. Please use --backend alchemlyb or omit this flag."
            )


def main(argv=None) -> int:
    """
    Main entry point for CLI

    Returns
    -------
    int
        Exit code (0 for success, 1 for error)
    """
    try:
        # Parse arguments
        args = parse_arguments(argv)

        # Setup logging
        setup_logging(args.verbose)
        logger = logging.getLogger(__name__)

        # Validate arguments
        validate_arguments(args)

        # Detect multi-estimator mode
        if args.all_estimators or len(args.estimator) > 1:
            # Multi-estimator mode
            estimators = ["TI", "BAR", "MBAR"] if args.all_estimators else args.estimator
            multi_mode = True
        else:
            # Single estimator mode
            estimators = [args.estimator[0]] if args.estimator else ["MBAR"]
            multi_mode = False

        # Print banner
        logger.info("=" * 70)
        logger.info("PRISM-FEP Analysis")
        logger.info("=" * 70)
        logger.info(f"Bound directory: {args.bound_dir}")
        logger.info(f"Unbound directory: {args.unbound_dir}")
        if multi_mode:
            logger.info(f"Estimators: {', '.join(estimators)}")
        else:
            logger.info(f"Estimator: {estimators[0]}")
        logger.info(f"Temperature: {args.temperature} K")
        logger.info("=" * 70)

        if multi_mode:
            # Multi-estimator mode
            # Validate backend (only alchemlyb supported for multi-estimator)
            if args.backend != "alchemlyb":
                raise ValueError(
                    f"Multi-estimator mode only supports 'alchemlyb' backend. "
                    f"Got: {args.backend}. Please use --backend alchemlyb or omit this flag."
                )

            # Create multi-estimator analyzer
            analyzer = FEPMultiEstimatorAnalyzer(
                bound_dirs=args.bound_dir,
                unbound_dirs=args.unbound_dir,
                estimators=estimators,
                temperature=args.temperature,
                backend=args.backend,
                show_progress=not args.no_progress,
                bootstrap_n_jobs=args.bootstrap_n_jobs,
            )

            # Run analysis
            logger.info("Starting FEP analysis with multiple estimators...")
            multi_results = analyzer.analyze()

            # Print summary for all estimators
            logger.info("")
            logger.info("=" * 70)
            logger.info("RESULTS SUMMARY - Multi-Estimator Comparison")
            logger.info("=" * 70)
            for name, results in multi_results.methods.items():
                logger.info(
                    f"{name}: ΔG = {results.delta_g:.2f} ± {results.delta_g_error:.2f} kcal/mol "
                    f"(bound: {results.delta_g_bound:.2f}, unbound: {results.delta_g_unbound:.2f})"
                )
            logger.info("-" * 70)
            logger.info(
                f"Agreement: {multi_results.comparison['agreement']} "
                f"(std: {multi_results.comparison['delta_g_std']:.2f} kcal/mol)"
            )
            logger.info("=" * 70)

            # Generate HTML report
            logger.info(f"Generating HTML report: {args.output}")
            generator = MultiEstimatorReportGenerator(multi_results)
            report_path = generator.generate(args.output)
            logger.info(f"✓ HTML report saved: {report_path}")

            # Save JSON if requested
            if args.json:
                logger.info(f"Saving JSON results: {args.json}")
                import json

                json_path = Path(args.json)
                with open(json_path, "w") as f:
                    json.dump(multi_results.to_dict(), f, indent=2)
                logger.info(f"✓ JSON results saved: {json_path}")

        else:
            # Single estimator mode (existing behavior)
            # Create analyzer
            analyzer = FEPAnalyzer(
                bound_dir=args.bound_dir,
                unbound_dir=args.unbound_dir,
                temperature=args.temperature,
                estimator=estimators[0],
                backend=args.backend,
                energy_components=args.components,
                bootstrap_n_jobs=args.bootstrap_n_jobs,
            )

            # Run analysis
            logger.info("Starting FEP analysis...")
            results = analyzer.analyze()

            # Print summary
            logger.info("")
            logger.info("=" * 70)
            logger.info("RESULTS SUMMARY")
            logger.info("=" * 70)
            logger.info(f"Binding Free Energy (ΔG): {results.delta_g:.2f} ± {results.delta_g_error:.2f} kcal/mol")
            logger.info(f"  Bound leg:   {results.delta_g_bound:.2f} kcal/mol")
            logger.info(f"  Unbound leg: {results.delta_g_unbound:.2f} kcal/mol")
            logger.info("=" * 70)

            # Generate HTML report
            logger.info(f"Generating HTML report: {args.output}")
            report_path = analyzer.generate_html_report(args.output)
            logger.info(f"✓ HTML report saved: {report_path}")

            # Save JSON if requested
            if args.json:
                logger.info(f"Saving JSON results: {args.json}")
                json_path = analyzer.save_results(args.json)
                logger.info(f"✓ JSON results saved: {json_path}")

        logger.info("")
        logger.info("Analysis complete!")

        return 0

    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        return 1
    except ValueError as e:
        logging.error(f"Invalid argument: {e}")
        return 1
    except Exception as e:
        logging.exception(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
