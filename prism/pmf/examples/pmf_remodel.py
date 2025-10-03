#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF System Remodeling Script

Config-driven approach for PMF system preparation.
All parameters are defined in configuration files.

Usage:
    # Use default local configuration
    python pmf_remodel.py -i ./gromacssim -o ./pmf_output

    # Use cluster configuration
    python pmf_remodel.py -i ./gromacssim -o ./pmf_output -c default_cluster_config.yaml

    # Use custom configuration
    python pmf_remodel.py -i ./gromacssim -o ./pmf_output -c my_custom.yaml
"""

import sys
import yaml
import argparse
from pathlib import Path

# Add PRISM to path - try multiple possible locations
script_path = Path(__file__).resolve()

# Strategy 1: Check if already installed (can import prism)
try:
    import prism
    # Already in path, no need to add
except ImportError:
    # Strategy 2: Find prism package relative to script location
    possible_prism_roots = [
        script_path.parent.parent.parent,  # If script in prism/pmf/examples/
        Path.cwd(),                         # Current working directory
        Path.cwd().parent,                  # Parent of current directory
    ]

    prism_found = False
    for prism_root in possible_prism_roots:
        if (prism_root / 'prism' / '__init__.py').exists():
            sys.path.insert(0, str(prism_root))
            prism_found = True
            break

    if not prism_found:
        print("Error: Cannot locate PRISM package")
        print("Please run this script from the PRISM root directory or install PRISM")
        print(f"\nSearched in:")
        for root in possible_prism_roots:
            print(f"  - {root}")
        sys.exit(1)

try:
    from prism.pmf import PMFBuilder
except ImportError as e:
    print(f"Error importing PRISM components: {e}")
    print("Please ensure PRISM is properly installed and accessible.")
    print(f"\nTried to find prism module in:")
    for root in possible_prism_roots:
        print(f"  - {root}")
    sys.exit(1)


def load_config(config_path: str) -> dict:
    """
    Load configuration with intelligent path search

    Search order:
    1. Absolute/relative path as provided
    2. Relative to current working directory
    3. In prism/pmf/configs/remodel_config/ directory (package installation)
    4. With .yaml extension added
    """
    # 1. If absolute path or relative path from cwd exists as-is
    config_file = Path(config_path)
    if config_file.is_absolute() and config_file.exists():
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)

    # 2. Try relative to current working directory
    config_file = Path.cwd() / config_path
    if config_file.exists():
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)

    # 3. Search in package's remodel_config directory
    # Find prism.pmf module location dynamically
    try:
        import prism.pmf
        pmf_module_path = Path(prism.pmf.__file__).parent
        config_dir = pmf_module_path / 'configs' / 'remodel_config'
    except (ImportError, AttributeError):
        # Fallback to script location if import fails
        config_dir = Path(__file__).parent.parent / 'configs' / 'remodel_config'

    config_file = config_dir / config_path
    if config_file.exists():
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)

    # 4. Try with .yaml extension in package directory
    if not config_path.endswith('.yaml'):
        config_file = config_dir / f"{config_path}.yaml"
        if config_file.exists():
            with open(config_file, 'r') as f:
                return yaml.safe_load(f)

    # Not found - show all searched locations
    searched_paths = [
        str(Path(config_path).resolve()),
        str(Path.cwd() / config_path),
        str(config_dir / config_path),
    ]
    if not config_path.endswith('.yaml'):
        searched_paths.append(str(config_dir / f"{config_path}.yaml"))

    raise FileNotFoundError(
        f"Configuration file not found: {config_path}\n"
        f"Searched in:\n" + "\n".join(f"  - {p}" for p in searched_paths)
    )


def display_config_info(config: dict):
    """Display configuration information"""
    print("\nConfiguration Summary:")
    print("-" * 50)

    # Execution mode
    exec_mode = config.get('execution', {}).get('mode', 'auto-equilibrate')
    print(f"Execution Mode: {exec_mode}")

    # Force field (show if specified, otherwise auto-detect)
    general = config.get('general', {})
    if 'forcefield' in general:
        print(f"Force Field: {general['forcefield']} (specified)")
    else:
        print(f"Force Field: auto-detect from MD system")

    if 'water_model' in general:
        print(f"Water Model: {general['water_model']} (specified)")
    else:
        print(f"Water Model: auto-detect from MD system")

    # PMF parameters
    box_config = config.get('box', {})
    print(f"Pulling Distance: {box_config.get('pulling_distance', 2.0)} nm")

    pmf_system = config.get('pmf_system', {})
    print(f"Reference Group: {pmf_system.get('reference_group', 'Protein')}")
    print(f"Moving Group: {pmf_system.get('moving_group', 'LIG')}")

    # Equilibration info
    if exec_mode == 'auto-equilibrate':
        print("\nEquilibration: Will run automatically")
    else:
        print("\nEquilibration: Scripts will be generated")
        if 'slurm' in config:
            slurm = config['slurm']
            print(f"  SLURM Partition: {slurm.get('partition', 'N/A')}")
            print(f"  CPUs per Task: {slurm.get('cpus_per_task', 'N/A')}")
            print(f"  Time Limit: {slurm.get('time', 'N/A')}")

    print("-" * 50)


def display_results(results: dict, execution_mode: str):
    """Display results based on execution mode"""
    print("\n" + "=" * 70)
    print("PMF SYSTEM REMODELING COMPLETED")
    print("=" * 70)

    print(f"\nStatus: {results.get('status', 'completed')}")
    print(f"System Directory: {results['system_dir']}")

    # Alignment info
    if 'alignment' in results:
        alignment = results['alignment']
        print(f"\nAlignment:")
        print(f"  Protein Centroid: {alignment.get('protein_centroid', 'N/A')}")
        print(f"  Ligand Centroid: {alignment.get('ligand_centroid', 'N/A')}")
        print(f"  Distance: {alignment.get('distance', 'N/A')} nm")

    # Final system info
    if 'final_system' in results:
        final_sys = results['final_system']
        print(f"\nFinal System:")
        print(f"  Structure: {final_sys.get('structure', 'N/A')}")
        print(f"  Topology: {final_sys.get('topology', 'N/A')}")
        print(f"  Equilibrated: {final_sys.get('equilibrated', False)}")

    # Mode-specific output
    if execution_mode == 'auto-equilibrate':
        print("\n" + "=" * 70)
        print("AUTO-EQUILIBRATE MODE - COMPLETED")
        print("=" * 70)
        print("\nSystem has been rebuilt and equilibrated automatically.")
        print("\nNext Steps:")
        print("  1. Verify equilibration results")
        print("  2. Proceed with PMF workflow using the equilibrated system")

    else:  # script-only
        print("\n" + "=" * 70)
        print("SCRIPT-ONLY MODE - SCRIPTS GENERATED")
        print("=" * 70)
        print("\nEquilibration scripts have been generated.")
        print("\nGenerated Scripts:")
        print(f"  - run_equilibration_local.sh")
        print(f"  - run_equilibration_slurm.sh")
        print("\nNext Steps:")
        print("  1. Upload system to server:")
        print(f"     scp -r {results['system_dir']} user@server:/work/")
        print("  2. Run equilibration on server:")
        print("     cd /work/GMX_PMF_SYSTEM")
        print("     bash run_equilibration_local.sh")
        print("     # OR")
        print("     sbatch run_equilibration_slurm.sh")
        print("  3. Proceed with PMF workflow after equilibration")

    print(f"\nAll files available in: {results['system_dir']}")


def main():
    """Main entry point for PMF system remodeling"""
    parser = argparse.ArgumentParser(
        description='PRISM PMF System Remodeling - Config-Driven Approach',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default local configuration (auto-equilibrate)
  python pmf_remodel.py -i ./md_results -o ./pmf_system

  # Use default cluster configuration (script-only)
  python pmf_remodel.py -i ./md_results -o ./pmf_system -c default_cluster_config.yaml

  # Use custom configuration
  python pmf_remodel.py -i ./md_results -o ./pmf_system -c my_custom.yaml

Available default configurations:
  - default_local_config.yaml    : Local development (auto-equilibrate)
  - default_cluster_config.yaml  : Cluster deployment (script-only)

Configuration files location: prism/pmf/configs/remodel_config/
        """
    )

    parser.add_argument('--input', '-i', required=True,
                       help='Input MD results directory')
    parser.add_argument('--output', '-o', required=True,
                       help='Output directory for PMF system')
    parser.add_argument('--config', '-c', type=str,
                       default='default_local_config.yaml',
                       help='Configuration file (default: default_local_config.yaml)')

    args = parser.parse_args()

    try:
        print("=" * 70)
        print("PRISM PMF SYSTEM REMODELING")
        print("=" * 70)
        print(f"\nInput MD Directory: {args.input}")
        print(f"Output Directory: {args.output}")
        print(f"Configuration File: {args.config}")

        # Load configuration
        config = load_config(args.config)
        print(f"\nConfiguration loaded successfully")

        # Display configuration summary
        display_config_info(config)

        # Get execution mode from config
        execution_mode = config.get('execution', {}).get('mode', 'auto-equilibrate')
        auto_equilibrate = (execution_mode == 'auto-equilibrate')

        # Get frame from config
        frame = config.get('extraction', {}).get('frame_number', -1)

        # Initialize PMF Builder
        print("\n" + "=" * 70)
        print("BUILDING PMF SYSTEM")
        print("=" * 70)

        pmf_builder = PMFBuilder(
            md_results_dir=args.input,
            output_dir=args.output,
            config=config
        )

        # Build system
        results = pmf_builder.build(
            frame=frame,
            equilibrate=auto_equilibrate
        )

        # Display results
        display_results(results, execution_mode)

    except FileNotFoundError as e:
        print(f"\nError: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\nPMF system build failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
