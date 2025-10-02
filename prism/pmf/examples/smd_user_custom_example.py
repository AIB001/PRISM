#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SMD Module - Full User Customization Examples

This example demonstrates how to use the fully customizable SMD module
where all parameters are user-defined without automatic calculations.
"""

import sys
from pathlib import Path

# Add PRISM to path
sys.path.insert(0, str(Path(__file__).parents[3]))

from prism.pmf.methods.smd import create_smd_module

def example_basic_usage():
    """Example 1: Basic usage with required parameters"""
    print("=== Example 1: Basic SMD Setup ===")

    smd = create_smd_module(
        pull_rate=0.005,      # 0.005 nm/ps pulling rate
        nsteps=5000000,       # 5 million steps = 10 ns simulation
        dt=0.002,             # 2 fs time step
        pull_k=1000.0         # 1000 kJ/mol/nm² force constant
    )

    print(f"Created SMD module with:")
    print(f"  - Pull rate: {smd.pull_rate} nm/ps")
    print(f"  - Total steps: {smd.nsteps:,}")
    print(f"  - Time step: {smd.dt} ps")
    print(f"  - Force constant: {smd.pull_k} kJ/mol/nm²")
    print(f"  - Total simulation time: {smd.nsteps * smd.dt / 1000:.1f} ns")
    print()

def example_custom_parameters():
    """Example 2: Custom parameters for specific needs"""
    print("=== Example 2: Custom Parameters ===")

    # For faster testing
    smd_fast = create_smd_module(
        pull_rate=0.02,       # Faster pulling
        nsteps=1000000,       # 1 million steps = 2 ns
        dt=0.002,
        pull_k=500.0,         # Lower force constant
        temperature=300.0,    # Room temperature
        pressure=1.5          # Higher pressure
    )

    print(f"Fast SMD setup:")
    print(f"  - Pull rate: {smd_fast.pull_rate} nm/ps (fast)")
    print(f"  - Simulation time: {smd_fast.nsteps * smd_fast.dt / 1000:.1f} ns (short)")
    print(f"  - Temperature: {smd_fast.temperature} K")
    print(f"  - Pressure: {smd_fast.pressure} bar")
    print()

def example_high_precision():
    """Example 3: High precision setup"""
    print("=== Example 3: High Precision Setup ===")

    smd_precise = create_smd_module(
        pull_rate=0.001,      # Very slow pulling
        nsteps=25000000,      # 25 million steps = 50 ns
        dt=0.001,             # 1 fs time step (more precise)
        pull_k=2000.0,        # High force constant
        temperature=310.0,    # Body temperature
        pressure=1.0          # Standard pressure
    )

    print(f"High precision SMD setup:")
    print(f"  - Pull rate: {smd_precise.pull_rate} nm/ps (very slow)")
    print(f"  - Simulation time: {smd_precise.nsteps * smd_precise.dt / 1000:.1f} ns (long)")
    print(f"  - Time step: {smd_precise.dt} ps (precise)")
    print(f"  - Force constant: {smd_precise.pull_k} kJ/mol/nm² (strong)")
    print()

def example_config_file():
    """Example 4: Using configuration file + user overrides"""
    print("=== Example 4: Config File + User Overrides ===")

    # Configuration from file
    config = {
        'smd': {
            'pull_rate': 0.003,
            'nsteps': 8000000,
            'dt': 0.002,
            'pull_k': 1500.0
        }
    }

    # User overrides specific parameters
    smd_mixed = create_smd_module(
        config=config,
        pull_rate=0.008,      # Override config pull_rate
        temperature=320.0     # Add temperature not in config
    )

    print(f"Mixed configuration:")
    print(f"  - Pull rate: {smd_mixed.pull_rate} nm/ps (user override)")
    print(f"  - Steps: {smd_mixed.nsteps:,} (from config)")
    print(f"  - dt: {smd_mixed.dt} ps (from config)")
    print(f"  - Pull k: {smd_mixed.pull_k} kJ/mol/nm² (from config)")
    print(f"  - Temperature: {smd_mixed.temperature} K (user added)")
    print()

def example_parameter_relationships():
    """Example 5: Understanding parameter relationships"""
    print("=== Example 5: Parameter Relationships ===")

    scenarios = [
        ("Quick test", 0.01, 2000000, 0.002),
        ("Standard", 0.005, 5000000, 0.002),
        ("Precise", 0.002, 15000000, 0.001),
    ]

    for name, rate, steps, dt in scenarios:
        simulation_time_ns = steps * dt / 1000
        expected_distance = rate * simulation_time_ns * 1000  # nm

        print(f"{name} scenario:")
        print(f"  - Pull rate: {rate} nm/ps")
        print(f"  - Steps: {steps:,}")
        print(f"  - Total time: {simulation_time_ns:.1f} ns")
        print(f"  - Expected distance: {expected_distance:.1f} nm")
        print(f"  - Distance per ns: {rate * 1000:.1f} nm/ns")
        print()

if __name__ == "__main__":
    print("SMD Module - Full User Customization Examples")
    print("=" * 50)
    print()

    example_basic_usage()
    example_custom_parameters()
    example_high_precision()
    example_config_file()
    example_parameter_relationships()

    print("Key Points:")
    print("- User must specify: pull_rate, nsteps, dt, pull_k")
    print("- No automatic calculations or defaults for core parameters")
    print("- Full control over simulation parameters")
    print("- Temperature and pressure have sensible defaults (310K, 1bar)")
    print("- User parameters always override config file values")