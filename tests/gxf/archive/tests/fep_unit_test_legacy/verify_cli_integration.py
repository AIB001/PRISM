#!/usr/bin/env python3
"""
Archived verification script for historical FEP CLI integration checks.
"""

import subprocess
import sys


def test_cli_help() -> bool:
    """Verify that historical FEP-related options are present in CLI help."""
    print("Testing CLI help text...")
    result = subprocess.run(["prism", "--help"], capture_output=True, text=True)

    required_options = [
        "--fep",
        "--mutant",
        "--fep-config",
        "--distance-cutoff",
        "--charge-strategy",
        "--lambda-windows",
    ]

    for opt in required_options:
        if opt in result.stdout:
            print(f"  ✓ {opt} found in help")
        else:
            print(f"  ✗ {opt} NOT found in help")
            return False
    return True


def test_prism_builder() -> bool:
    """Verify PRISMBuilder accepts the historical FEP parameters."""
    print("\nTesting PRISMBuilder FEP parameters...")

    try:
        from prism.builder.core import PRISMBuilder
    except ImportError as exc:
        print(f"  ✗ Import error: {exc}")
        return False

    try:
        builder = PRISMBuilder(
            protein_path="/tmp/protein.pdb",
            ligand_paths=["/tmp/ligand.mol2"],
            output_dir="/tmp/output",
            fep_mode=True,
            mutant_ligand="/tmp/mutant.mol2",
            distance_cutoff=0.5,
            charge_strategy="mean",
            lambda_windows=15,
        )
    except Exception as exc:
        print(f"  ✗ Error: {exc}")
        return False

    checks = [
        builder.fep_mode is True,
        builder.mutant_ligand == "/tmp/mutant.mol2",
        builder.distance_cutoff == 0.5,
        builder.charge_strategy == "mean",
        builder.lambda_windows == 15,
    ]
    passed = all(checks)
    print("  ✓ PRISMBuilder accepts all FEP parameters" if passed else "  ✗ PRISMBuilder parameter check failed")
    return passed


def test_fep_imports() -> bool:
    """Verify FEP modules can be imported."""
    print("\nTesting FEP module imports...")
    try:
        from prism.fep import Atom, AtomMapping, FEstimator, HybridAtom, ITPBuilder, XVGParser
    except ImportError as exc:
        print(f"  ✗ Import error: {exc}")
        return False

    _ = (Atom, AtomMapping, HybridAtom, ITPBuilder, XVGParser, FEstimator)
    print("  ✓ All FEP modules imported successfully")
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("FEP CLI Integration Verification")
    print("=" * 60)

    results = [
        ("CLI Help", test_cli_help()),
        ("PRISMBuilder", test_prism_builder()),
        ("Module Imports", test_fep_imports()),
    ]

    print("\n" + "=" * 60)
    print("Summary:")
    print("=" * 60)

    all_passed = True
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {name}")
        if not passed:
            all_passed = False

    print("=" * 60)
    if all_passed:
        print("\n✓ All integration tests passed!")
        sys.exit(0)

    print("\n✗ Some integration tests failed")
    sys.exit(1)
