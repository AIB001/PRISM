#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Verify case.yaml and config.conf consistency across all test cases

This script checks that basic FEP parameters are consistent between
case.yaml and config.conf files in all test directories.
"""

import sys
from pathlib import Path
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from prism.fep.config import read_fep_config


def check_config_consistency(test_dir):
    """
    Check if case.yaml and config.conf have consistent FEP parameters

    Parameters
    ----------
    test_dir : str
        Path to test directory containing case.yaml and config.conf

    Returns
    -------
    dict
        Comparison results
    """
    case_yaml = Path(test_dir) / "case.yaml"
    config_conf = Path(test_dir) / "config.conf"

    if not case_yaml.exists():
        return {"status": "skip", "reason": f"case.yaml not found: {case_yaml}"}
    if not config_conf.exists():
        return {"status": "skip", "reason": f"config.conf not found: {config_conf}"}

    # Read config.conf using FEbuilder parser
    try:
        conf_params = read_fep_config(str(config_conf))
    except Exception as e:
        return {"status": "error", "reason": f"Failed to read config.conf: {e}"}

    # Read case.yaml
    try:
        with open(case_yaml) as f:
            yaml_data = yaml.safe_load(f)
    except Exception as e:
        return {"status": "error", "reason": f"Failed to read case.yaml: {e}"}

    # Extract model parameters from case.yaml
    if "model" not in yaml_data:
        return {"status": "error", "reason": "No [model] section in case.yaml"}

    model = yaml_data["model"]

    # Compare key parameters
    comparisons = []

    # charge_common
    yaml_charge_common = model.get("charge_common", "mean")
    conf_charge_common = conf_params.get("charge_common", "mean")
    charge_common_match = yaml_charge_common == conf_charge_common
    comparisons.append(
        {
            "param": "charge_common",
            "yaml": yaml_charge_common,
            "config": conf_charge_common,
            "match": charge_common_match,
        }
    )

    # charge_reception
    yaml_charge_reception = model.get("charge_reception", "surround")
    conf_charge_reception = conf_params.get("charge_reception", "surround")
    charge_reception_match = yaml_charge_reception == conf_charge_reception
    comparisons.append(
        {
            "param": "charge_reception",
            "yaml": yaml_charge_reception,
            "config": conf_charge_reception,
            "match": charge_reception_match,
        }
    )

    # distance
    yaml_distance = model.get("distance", 10)
    conf_distance = conf_params.get("distance", 10)
    distance_match = abs(yaml_distance - conf_distance) < 0.1
    comparisons.append({"param": "distance", "yaml": yaml_distance, "config": conf_distance, "match": distance_match})

    # Check if all match
    all_match = all(c["match"] for c in comparisons)

    return {"status": "ok" if all_match else "mismatch", "comparisons": comparisons, "all_match": all_match}


def main():
    """Check all test directories"""

    print("=" * 70)
    print("Verifying case.yaml and config.conf consistency")
    print("=" * 70)

    # Find all test directories
    test_base = Path("tests/gxf/FEP/test")

    if not test_base.exists():
        print(f"Test base directory not found: {test_base}")
        return 1

    # Find all directories containing case.yaml
    test_dirs = []
    for case_file in test_base.rglob("case.yaml"):
        test_dir = case_file.parent
        if test_dir.name != "unit_test":  # Skip unit_test for now
            test_dirs.append(test_dir)

    print(f"\nFound {len(test_dirs)} test directories\n")

    # Check each directory
    results = {}
    mismatch_count = 0
    error_count = 0

    for test_dir in sorted(test_dirs):
        rel_path = test_dir.relative_to(test_base)
        print(f"[{rel_path}]")

        result = check_config_consistency(test_dir)
        results[str(rel_path)] = result

        if result["status"] == "ok":
            print(f"  ✓ All parameters match")
            for comp in result["comparisons"]:
                print(f"    {comp['param']}: {comp['yaml']}")
        elif result["status"] == "mismatch":
            print(f"  ✗ Parameters don't match:")
            mismatch_count += 1
            for comp in result["comparisons"]:
                status = "✓" if comp["match"] else "✗"
                print(f"    {status} {comp['param']}: yaml={comp['yaml']}, config={comp['config']}")
        elif result["status"] == "error":
            print(f"  ✗ Error: {result['reason']}")
            error_count += 1
        else:
            print(f"  ⊘ Skipped: {result['reason']}")

        print()

    # Summary
    print("=" * 70)
    print("Summary:")
    print(f"  Total checked: {len(test_dirs)}")
    print(f"  ✓ All match: {len(test_dirs) - mismatch_count - error_count}")
    print(f"  ✗ Mismatch: {mismatch_count}")
    print(f"  ✗ Error: {error_count}")
    print("=" * 70)

    return 0 if mismatch_count == 0 and error_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
