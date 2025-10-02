#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF System Remodeling Script

This script uses PRISM's built-in PMFBuilder for proper PMF system preparation.
It follows PRISM's established architecture and workflows.

Usage:
    python pmf_remodel.py --input ./gromacssim --output ./pmf_output --config config.yaml
"""

import sys
import argparse
from pathlib import Path

# Add PRISM to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

try:
    from prism.pmf.builders.pmf_builder import PMFBuilder
    from prism.utils.config import ConfigurationManager
except ImportError as e:
    print(f"❌ Error importing PRISM components: {e}")
    print("Please ensure PRISM is properly installed and accessible.")
    sys.exit(1)

def main():
    """Main entry point for PRISM PMF system remodeling"""
    parser = argparse.ArgumentParser(description='PRISM PMF System Remodeling Script')
    parser.add_argument('--input', '-i', required=True,
                       help='Input MD results directory (e.g., gromacssim)')
    parser.add_argument('--output', '-o', required=True,
                       help='Output directory for PMF system')
    parser.add_argument('--config', '-c', type=str, default=None,
                       help='Path to PMF configuration YAML file')
    parser.add_argument('--frame', '-f', type=int, default=-1,
                       help='Frame to extract from trajectory (-1 for last frame)')
    parser.add_argument('--equilibrate', action='store_true', default=True,
                       help='Run equilibration after system rebuild')
    parser.add_argument('--no-equilibrate', dest='equilibrate', action='store_false',
                       help='Skip equilibration, only generate scripts')

    args = parser.parse_args()

    try:
        print("="*70)
        print("PRISM PMF SYSTEM REMODELING")
        print("="*70)
        print(f"Input MD Directory: {args.input}")
        print(f"Output Directory: {args.output}")
        print(f"Configuration: {args.config or 'default'}")
        print(f"Frame: {args.frame}")
        print(f"Equilibration: {'enabled' if args.equilibrate else 'disabled'}")
        print()

        # Load configuration
        config = {}
        if args.config:
            config_manager = ConfigurationManager()
            config = config_manager.load_config(args.config)
            print(f"✓ Loaded configuration from: {args.config}")
        else:
            print("ℹ  Using default configuration")

        # Initialize PMF Builder
        pmf_builder = PMFBuilder(
            md_results_dir=args.input,
            output_dir=args.output,
            config=config
        )

        print(f"✓ Initialized PMF Builder")
        print(f"  MD Results Directory: {pmf_builder.md_results_dir}")
        print(f"  Output Directory: {pmf_builder.output_dir}")

        # Build PMF system
        print("\n" + "="*50)
        print("BUILDING PMF SYSTEM")
        print("="*50)

        results = pmf_builder.build(
            frame=args.frame,
            equilibrate=args.equilibrate
        )

        # Display results
        print("\n" + "="*70)
        print("PMF SYSTEM BUILD COMPLETED")
        print("="*70)

        print(f"✓ Status: {results.get('status', 'completed')}")
        print(f"✓ System Directory: {results['system_dir']}")

        if 'final_system' in results:
            final_sys = results['final_system']
            print(f"✓ Final Structure: {final_sys.get('structure', 'N/A')}")
            print(f"✓ Final Topology: {final_sys.get('topology', 'N/A')}")
            print(f"✓ Equilibrated: {final_sys.get('equilibrated', False)}")

        if 'alignment' in results:
            alignment = results['alignment']
            print(f"✓ Protein Centroid: {alignment.get('protein_centroid', 'N/A')}")
            print(f"✓ Ligand Centroid: {alignment.get('ligand_centroid', 'N/A')}")
            print(f"✓ Distance: {alignment.get('distance', 'N/A')} nm")

        if 'equilibration' in results and args.equilibrate:
            eq_results = results['equilibration']
            if 'script_files' in eq_results:
                print("\n📋 Equilibration Scripts Generated:")
                for script_type, script_path in eq_results['script_files'].items():
                    print(f"  {script_type}: {script_path}")
            elif 'final_system' in eq_results:
                print(f"✓ Equilibration completed: {eq_results['final_system']}")

        print("\n🎉 PMF system is ready for steered molecular dynamics!")

        if not args.equilibrate:
            print("\n📝 Next Steps:")
            print("1. Review the generated equilibration scripts")
            print("2. Run equilibration (EM → NVT → NPT)")
            print("3. Proceed with SMD simulations")
        else:
            print("\n📝 Next Steps:")
            print("1. Verify equilibration results")
            print("2. Proceed with SMD simulations using the final system")

        print(f"\n📂 All files available in: {results['system_dir']}")

    except Exception as e:
        print(f"❌ PMF system build failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()