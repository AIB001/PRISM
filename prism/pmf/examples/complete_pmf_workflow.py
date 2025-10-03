#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Complete PMF Workflow Example

This script demonstrates the recommended two-API approach for PMF calculations:
1. PMFBuilder - System remodeling and preparation
2. PMFRunner - Complete PMF workflow execution

Usage:
    python complete_pmf_workflow.py --input ./gromacssim --output ./pmf_results
"""

import sys
import argparse
from pathlib import Path

# Add PRISM to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

try:
    # Use only the two recommended APIs
    from prism.pmf import PMFBuilder, PMFRunner
except ImportError as e:
    print(f"Error importing PRISM PMF module: {e}")
    print("Please ensure PRISM is properly installed.")
    sys.exit(1)


def main():
    """Complete PMF workflow using two-API design"""
    parser = argparse.ArgumentParser(
        description='Complete PMF Workflow - Two API Design',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete workflow with default settings
  python complete_pmf_workflow.py --input ./gromacssim --output ./pmf_results

  # With custom configuration
  python complete_pmf_workflow.py --input ./gromacssim --output ./pmf_results \\
      --config pmf_config.yaml

  # Script-only mode (generate scripts for server execution)
  python complete_pmf_workflow.py --input ./gromacssim --output ./pmf_results \\
      --script-only
        """
    )
    parser.add_argument('--input', '-i', required=True,
                       help='Input MD results directory')
    parser.add_argument('--output', '-o', required=True,
                       help='Output directory for PMF results')
    parser.add_argument('--config', '-c', type=str, default=None,
                       help='PMF configuration YAML file')
    parser.add_argument('--auto-equilibrate', action='store_true', default=True,
                       help='Auto-run equilibration after remodeling (default)')
    parser.add_argument('--script-only', dest='auto_equilibrate', action='store_false',
                       help='Generate equilibration scripts only (for server execution)')
    parser.add_argument('--steps', nargs='+', default=None,
                       help='Specific workflow steps to run (default: all)')

    args = parser.parse_args()

    print("="*70)
    print("PRISM PMF COMPLETE WORKFLOW")
    print("="*70)
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Config: {args.config or 'default'}")
    print(f"Equilibration Mode: {'auto-run' if args.auto_equilibrate else 'script-only'}")
    print()

    # ========================================================================
    # STEP 1: System Remodeling with PMFBuilder
    # ========================================================================
    print("="*70)
    print("STEP 1: SYSTEM REMODELING")
    print("="*70)

    remodel_output = Path(args.output) / "pmf_system"

    builder = PMFBuilder(
        md_results_dir=args.input,
        output_dir=str(remodel_output)
    )

    print("Building PMF-optimized system...")
    remodel_results = builder.build(equilibrate=args.auto_equilibrate)

    print(f"\nRemodeling completed:")
    print(f"  System directory: {remodel_results['system_dir']}")
    print(f"  Equilibrated: {remodel_results.get('equilibrated', False)}")

    if not args.auto_equilibrate:
        print("\nEquilibration scripts generated. Please run:")
        print(f"  cd {remodel_results['system_dir']}")
        print(f"  bash run_equilibration_local.sh")
        print("\nThen run PMF workflow with the equilibrated system.")
        return

    # ========================================================================
    # STEP 2: PMF Workflow with PMFRunner
    # ========================================================================
    print("\n" + "="*70)
    print("STEP 2: PMF WORKFLOW")
    print("="*70)

    pmf_system_dir = Path(remodel_results['system_dir'])
    pmf_output = Path(args.output) / "pmf_calculations"

    # Load configuration if provided
    config = args.config

    runner = PMFRunner(config=config)

    print("Running complete PMF workflow...")
    print("This includes: SMD -> Umbrella Sampling -> WHAM Analysis")
    print()

    try:
        workflow_results = runner.run_complete_workflow(
            md_system_dir=str(pmf_system_dir),
            output_dir=str(pmf_output),
            steps=args.steps
        )

        # ========================================================================
        # STEP 3: Display Results
        # ========================================================================
        print("\n" + "="*70)
        print("PMF WORKFLOW COMPLETED")
        print("="*70)

        print(f"Status: {workflow_results.get('status', 'unknown')}")
        print(f"Steps executed: {', '.join(workflow_results.get('steps_executed', []))}")

        if 'binding_energy' in workflow_results:
            be = workflow_results['binding_energy']
            unit = workflow_results.get('analysis', {}).get('energy_unit', 'kcal/mol')
            print(f"\nBinding Energy: {be['value']:.2f} +/- {be.get('error', 0):.2f} {unit}")
            print(f"  Minimum at: {be.get('min_distance', 'N/A')} nm")
            print(f"  Maximum at: {be.get('max_distance', 'N/A')} nm")

        print(f"\nResults directory: {workflow_results['output_directory']}")

        if 'analysis' in workflow_results:
            analysis = workflow_results['analysis']
            if 'plots' in analysis:
                print("\nGenerated plots:")
                for plot_name, plot_path in analysis['plots'].items():
                    print(f"  - {plot_name}: {plot_path}")

        print("\nComplete workflow finished successfully!")

    except Exception as e:
        print(f"\nError during PMF workflow: {e}")
        print("Check log files for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
