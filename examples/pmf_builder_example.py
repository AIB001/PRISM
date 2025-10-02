#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Builder Example

This example demonstrates how to use the PMF Builder to reconstruct 
protein-ligand systems from MD results with PMF-optimized geometry.
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import prism.pmf as pmf


def example_pmf_builder_usage():
    """Example of using PMF Builder to reconstruct systems"""
    
    print("=== PRISM PMF Builder Example ===\n")
    
    # Example 1: Direct PMF Builder usage
    print("1. Direct PMF Builder Usage:")
    print("-" * 30)
    
    try:
        # Create PMF builder instance
        builder = pmf.pmf_builder(
            md_results_dir="./test/4xb4_md_results",  # MD results directory
            output_dir="./pmf_system",                # PMF system output
            z_extension=2.5,                          # Extended Z-axis length
            reference_group='Protein',
            moving_group='LIG'
        )
        
        print(f"✓ PMF Builder created")
        print(f"  Source: ./test/4xb4_md_results")
        print(f"  Output: ./pmf_system")
        print(f"  Z extension: 2.5 nm")
        
        # Build PMF-optimized system with full equilibration
        # build_results = builder.build(frame=-1, equilibrate=True)
        print("  → Ready to build and equilibrate PMF-optimized system")
        print("    • Extract structures from MD trajectory") 
        print("    • Calculate centroids and Z-axis alignment")
        print("    • Rebuild system with extended Z-box")
        print("    • Run EM → NVT → NPT equilibration")
        
    except Exception as e:
        print(f"  Note: {e}")
        print("  (This is expected without actual MD results)")
    
    print()
    
    # Example 2: PMF System with integrated builder
    print("2. PMF System with Integrated Builder:")
    print("-" * 40)
    
    try:
        # Create PMF system with rebuild option enabled
        pmf_system = pmf.pmf_system(
            system_dir="./test/4xb4_md_results",
            output_dir="./pmf_calculations",
            rebuild_system=True,  # Enable system reconstruction
            config={
                'box': {'z_extension': 3.0},
                'pmf_system': {
                    'reference_group': 'Protein',
                    'moving_group': 'LIG'
                }
            }
        )
        
        print(f"✓ PMF System created with builder enabled")
        print(f"  Source: ./test/4xb4_md_results") 
        print(f"  Output: ./pmf_calculations")
        print(f"  Rebuild: Enabled")
        
        # Rebuild system would extract, align, extend box, and equilibrate
        # rebuild_results = pmf_system.rebuild_system()
        print("  → Ready to rebuild, equilibrate and optimize system")
        print("    • System reconstruction with Z-axis alignment")
        print("    • Automatic EM/NVT/NPT equilibration")
        
        # Then run PMF calculations  
        # pmf_results = pmf_system.run(mode='auto')
        print("  → Ready to run PMF calculations on equilibrated system")
        
    except Exception as e:
        print(f"  Note: {e}")
        print("  (This is expected without actual MD results)")
    
    print()
    
    # Example 3: Complete workflow
    print("3. Complete PMF Workflow with Builder:")
    print("-" * 40)
    
    workflow_example = '''
# Complete workflow example:
import prism

# Step 1: Build initial system (standard PRISM)
system = prism.system("protein.pdb", "ligand.mol2")
md_results = system.build()

# Step 2: Create PMF system with optimization
pmf = prism.pmf.pmf_system(
    md_results['output_dir'], 
    "./pmf_analysis",
    rebuild_system=True,
    z_extension=2.0
)

# Step 3: Rebuild system for PMF (extract, align, extend)
rebuild_results = pmf.rebuild_system()
print(f"Protein centroid: {rebuild_results['alignment']['protein_centroid']}")
print(f"Ligand centroid: {rebuild_results['alignment']['ligand_centroid']}")
print(f"Initial distance: {rebuild_results['alignment']['initial_distance']:.3f} nm")

# Step 4: Run PMF calculation
results = pmf.run(mode='auto')
print(f"Binding energy: {results['binding_energy']['value']:.2f} kcal/mol")
'''
    
    print(workflow_example)
    
    print("=== PMF Builder Features ===\n")
    
    features = [
        "✓ Structure Extraction: Separates Protein and LIG from MD trajectory",
        "✓ Centroid Calculation: Computes geometric centers of both components", 
        "✓ Z-axis Alignment: Aligns protein-ligand vector to Z-axis for pulling",
        "✓ Box Extension: Extends Z-axis to accommodate steered MD simulation",
        "✓ System Rebuild: Creates complete GROMACS system optimized for PMF",
        "✓ Complete Equilibration: Automatic EM → NVT → NPT equilibration",
        "✓ Flexible Control: Full or step-by-step equilibration control",
        "✓ Status Monitoring: Real-time equilibration status and validation", 
        "✓ PRISM Integration: Seamlessly integrated with PRISM workflow patterns"
    ]
    
    for feature in features:
        print(feature)
    
    print("\n=== Configuration Options ===\n")
    
    config_options = {
        'reference_group': 'Group name for reference (usually "Protein")',
        'moving_group': 'Group name for moving part (usually "LIG")',  
        'box.z_extension': 'Additional Z-axis length for pulling (nm)',
        'extraction.frame_number': 'Frame to extract (-1 for last frame)',
        'pmf_system.align_axis': 'Axis for alignment (default: "z")',
        'pmf_system.center_system': 'Whether to center system (default: True)',
        'equilibration.em.nsteps': 'Energy minimization steps (default: 50000)',
        'equilibration.nvt.nsteps': 'NVT equilibration steps (default: 50000)',
        'equilibration.npt.nsteps': 'NPT equilibration steps (default: 100000)',
        'equilibration.*.temperature': 'Simulation temperature (default: 310.0 K)'
    }
    
    for key, desc in config_options.items():
        print(f"  {key}: {desc}")
    
    print("\n✅ PMF Builder ready for protein-ligand PMF calculations!")


if __name__ == "__main__":
    example_pmf_builder_usage()