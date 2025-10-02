#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF API Usage Examples

This example demonstrates the clean, unified API for PMF calculations,
replacing shell scripts with elegant Python calls.
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import prism.pmf as pmf


def example_complete_workflow():
    """Example: Complete PMF workflow with unified API"""
    print("=== Complete PMF Workflow Example ===\n")
    
    try:
        # One-line complete PMF workflow
        results = pmf.run_pmf_workflow(
            md_system_dir="./test/4xb4_md_results",
            output_dir="./pmf_calculation"
        )
        
        print("✅ PMF workflow completed!")
        print(f"Duration: {results['duration']}")
        print(f"Status: {results['status']}")
        
        if 'binding_energy' in results:
            be = results['binding_energy']
            print(f"Binding Energy: {be['value']:.2f} kcal/mol")
        
        print(f"Results directory: {results['output_directory']}")
        
    except Exception as e:
        print(f"Note: {e} (Expected without actual MD results)")


def example_custom_config():
    """Example: PMF workflow with custom configuration"""
    print("=== Custom Configuration Example ===\n")
    
    # Create custom configuration
    custom_config = {
        'builder': {
            'rebuild_system': True,
            'protein_forcefield': 'charmm27',
            'ligand_forcefield': 'gaff',
            'z_extension': 3.0  # Larger pulling space
        },
        'smd': {
            'pull_rate': 0.01,      # Faster pulling
            'pull_k': 1500.0,       # Stronger restraint
            'nsteps': 1250000       # Shorter simulation (2.5 ns)
        },
        'umbrella': {
            'sample_interval_near': 0.05,    # Dense sampling
            'production_time_ps': 15000      # 15 ns per window
        }
    }
    
    try:
        results = pmf.run_pmf_workflow(
            md_system_dir="./md_system",
            output_dir="./custom_pmf",
            config=custom_config
        )
        
        print("✅ Custom PMF workflow completed!")
        
    except Exception as e:
        print(f"Note: {e} (Expected without actual MD results)")


def example_step_by_step():
    """Example: Step-by-step PMF execution"""
    print("=== Step-by-Step Execution Example ===\n")
    
    try:
        # Step 1: System builder only
        builder_results = pmf.run_pmf_step(
            md_system_dir="./md_system",
            output_dir="./step_by_step_pmf", 
            step="builder"
        )
        print("✅ Step 1 (Builder) completed")
        
        # Step 2: SMD preparation
        smd_results = pmf.run_pmf_step(
            md_system_dir="./md_system",
            output_dir="./step_by_step_pmf",
            step="smd"
        )
        print("✅ Step 2 (SMD) prepared")
        print(f"Manual command: {smd_results['smd']['manual_command']}")
        
        # Step 3: Umbrella sampling (after SMD completion)
        # umbrella_results = pmf.run_pmf_step(..., step="umbrella")
        
        # Step 4: Analysis (after umbrella completion)
        # analysis_results = pmf.run_pmf_step(..., step="analysis")
        
    except Exception as e:
        print(f"Note: {e} (Expected without actual MD results)")


def example_configuration_templates():
    """Example: Using configuration templates"""
    print("=== Configuration Templates Example ===\n")
    
    # Create different configuration templates
    templates = ['default', 'fast', 'accurate']
    
    for template in templates:
        config_file = f"pmf_config_{template}.yaml"
        
        try:
            pmf.create_pmf_config(config_file, template=template)
            print(f"✅ Created {template} configuration: {config_file}")
            
            # Show how to use the template
            print(f"Usage: pmf.run_pmf_workflow(md_dir, output_dir, config='{config_file}')")
            
        except Exception as e:
            print(f"Template creation failed: {e}")


def example_pmf_runner_class():
    """Example: Using PMFRunner class for advanced control"""
    print("=== PMFRunner Class Example ===\n")
    
    try:
        # Initialize runner with configuration
        runner = pmf.PMFRunner(config="pmf_config_default.yaml")
        
        # Run specific steps
        results = runner.run_complete_workflow(
            md_system_dir="./md_system",
            output_dir="./advanced_pmf",
            steps=['builder', 'smd']  # Only run builder and SMD
        )
        
        print("✅ Selective workflow completed")
        print(f"Steps executed: {results['steps_executed']}")
        
    except Exception as e:
        print(f"Note: {e} (Expected without actual MD results)")


def show_api_comparison():
    """Show the difference between old and new approaches"""
    print("=== API Comparison ===\n")
    
    print("❌ Old approach (shell script):")
    print("./scripts/run_pmf.sh -s ./md_system -o ./pmf_calc -c config.yaml")
    print("  - Requires shell scripting")
    print("  - Hard to integrate with Python workflows")
    print("  - Limited error handling")
    print("  - Platform dependent")
    print()
    
    print("✅ New approach (unified Python API):")
    print("import prism.pmf as pmf")
    print("results = pmf.run_pmf_workflow('./md_system', './pmf_calc')")
    print("  - Pure Python integration")
    print("  - Comprehensive error handling")
    print("  - Programmatic access to results")
    print("  - Cross-platform compatibility")
    print("  - Easy to embed in larger workflows")
    print()


def main():
    """Main function demonstrating all PMF API examples"""
    print("PRISM PMF Unified API Examples")
    print("=" * 50)
    print()
    
    # Show API comparison
    show_api_comparison()
    
    # Example 1: Complete workflow
    example_complete_workflow()
    print()
    
    # Example 2: Custom configuration
    example_custom_config()
    print()
    
    # Example 3: Step-by-step execution
    example_step_by_step()
    print()
    
    # Example 4: Configuration templates
    example_configuration_templates()
    print()
    
    # Example 5: Advanced PMFRunner usage
    example_pmf_runner_class()
    print()
    
    print("=== Recommended Usage Patterns ===")
    print()
    print("1. Quick PMF calculation:")
    print("   results = pmf.run_pmf_workflow('./md_system', './pmf_results')")
    print()
    print("2. Custom parameters:")
    print("   config = {'smd': {'pull_rate': 0.01}, 'umbrella': {'production_time_ps': 15000}}")
    print("   results = pmf.run_pmf_workflow('./md_system', './pmf_results', config=config)")
    print()
    print("3. Step-by-step control:")
    print("   pmf.run_pmf_step('./md_system', './pmf_results', 'builder')")
    print("   pmf.run_pmf_step('./md_system', './pmf_results', 'smd')")
    print("   # ... manual execution ...")
    print("   pmf.run_pmf_step('./md_system', './pmf_results', 'analysis')")
    print()
    print("4. Integration with existing workflows:")
    print("   # Standard PRISM workflow")
    print("   system = prism.system('protein.pdb', 'ligand.mol2')")
    print("   build_results = system.build()")
    print("   # Seamless PMF calculation")
    print("   pmf_results = pmf.run_pmf_workflow(build_results['output_dir'], './pmf')")
    print()
    
    print("✅ PMF API examples completed!")
    print("The unified Python API provides a clean, flexible interface for PMF calculations.")


if __name__ == "__main__":
    main()