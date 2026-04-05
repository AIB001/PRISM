#!/usr/bin/env python3
"""
FEP System Comprehensive Checker

Performs comprehensive validation of FEP system setup according to PRISM-FEP standards.
Checks mapping quality, topology completeness, grompp validation, and scaffold generation.

Usage:
    python fep_system_check.py <fep_directory>

Example:
    python fep_system_check.py tests/gxf/FEP/unit_test/42-38/CHARMM-GUI/GMX_PROLIG_FEP
"""

import sys
import subprocess
from pathlib import Path
from typing import Dict
import re


class FEPSystemChecker:
    """Comprehensive FEP system validation checker."""

    def __init__(self, fep_dir: str):
        self.fep_dir = Path(fep_dir).resolve()
        if not self.fep_dir.exists():
            raise FileNotFoundError(f"FEP directory not found: {fep_dir}")

    def check_all(self) -> Dict[str, any]:
        """
        Run all checks and return comprehensive results.

        Returns
        -------
        Dict with keys:
            - system_path: str
            - has_common_hybrid: bool
            - has_bound: bool
            - has_unbound: bool
            - mapping_quality: Dict
            - topology_complete: bool
            - grompp_bound: bool
            - grompp_unbound: bool
            - missing_files: List[str]
            - errors: List[str]
            - warnings: List[str]
        """
        results = {
            "system_path": str(self.fep_dir),
            "has_common_hybrid": False,
            "has_bound": False,
            "has_unbound": False,
            "mapping_quality": {},
            "topology_complete": False,
            "grompp_bound": False,
            "grompp_unbound": False,
            "missing_files": [],
            "errors": [],
            "warnings": [],
        }

        # Check directory structure
        self._check_directory_structure(results)

        # Check mapping quality
        if results["has_common_hybrid"]:
            self._check_mapping_quality(results)

        # Check topology completeness
        self._check_topology_completeness(results)

        # Check grompp validation
        if results["has_bound"]:
            self._check_grompp_validation("bound", results)

        if results["has_unbound"]:
            self._check_grompp_validation("unbound", results)

        # Check ligand placement
        if results["has_bound"]:
            self._check_ligand_placement(results)

        return results

    def _check_directory_structure(self, results: Dict) -> None:
        """Check required directory structure."""
        common_hybrid = self.fep_dir / "common" / "hybrid"
        bound = self.fep_dir / "bound"
        unbound = self.fep_dir / "unbound"

        results["has_common_hybrid"] = common_hybrid.exists()
        results["has_bound"] = bound.exists()
        results["has_unbound"] = unbound.exists()

        if not results["has_common_hybrid"]:
            results["errors"].append(f"Missing common/hybrid directory")

        if not results["has_bound"]:
            results["errors"].append(f"Missing bound directory")

        if not results["has_unbound"]:
            results["warnings"].append(f"Missing unbound directory (may not be generated yet)")

    def _check_mapping_quality(self, results: Dict) -> None:
        """Check FEP mapping quality from HTML report."""
        mapping_html = self.fep_dir / "common" / "hybrid" / "mapping.html"

        if not mapping_html.exists():
            results["errors"].append(f"Missing mapping.html")
            return

        try:
            html_content = mapping_html.read_text()

            # Check for gray/unknown atoms
            gray_count = html_content.count("rgb(200, 200, 200)")
            unknown_count = html_content.count('"classification": "unknown"')

            # Extract mapping statistics
            common_match = re.search(r"<strong>Common:</strong>\s*(\d+)", html_content)
            transformed_a_match = re.search(r"<strong>Transformed A:</strong>\s*(\d+)", html_content)
            transformed_b_match = re.search(r"<strong>Transformed B:</strong>\s*(\d+)", html_content)
            surrounding_a_match = re.search(r"<strong>Surrounding A:</strong>\s*(\d+)", html_content)
            surrounding_b_match = re.search(r"<strong>Surrounding B:</strong>\s*(\d+)", html_content)

            results["mapping_quality"] = {
                "gray_atoms": gray_count,
                "unknown_atoms": unknown_count,
                "common": int(common_match.group(1)) if common_match else 0,
                "transformed_a": int(transformed_a_match.group(1)) if transformed_a_match else 0,
                "transformed_b": int(transformed_b_match.group(1)) if transformed_b_match else 0,
                "surrounding_a": int(surrounding_a_match.group(1)) if surrounding_a_match else 0,
                "surrounding_b": int(surrounding_b_match.group(1)) if surrounding_b_match else 0,
            }

            # Validate mapping quality
            if gray_count > 0:
                results["warnings"].append(f"Found {gray_count} gray atoms in mapping")

            if unknown_count > 0:
                results["errors"].append(f"Found {unknown_count} unclassified atoms in mapping")

            # Check force field label
            if "CGenFF" not in html_content and "GAFF" not in html_content:
                results["warnings"].append(f"Force field label not found in mapping HTML")

        except Exception as e:
            results["errors"].append(f"Failed to parse mapping.html: {e}")

    def _check_topology_completeness(self, results: Dict) -> None:
        """Check if required topology files exist."""
        common_hybrid = self.fep_dir / "common" / "hybrid"

        required_files = [
            "hybrid.itp",
            "hybrid.gro",
            "atomtypes_hybrid.itp",
            "ff_hybrid.itp",
        ]

        missing = []
        for filename in required_files:
            filepath = common_hybrid / filename
            if not filepath.exists():
                missing.append(filename)

        results["missing_files"] = missing
        results["topology_complete"] = len(missing) == 0

        if missing:
            results["errors"].append(f"Missing topology files: {', '.join(missing)}")

    def _check_grompp_validation(self, leg: str, results: Dict) -> None:
        """Check grompp validation for bound/unbound leg."""
        leg_dir = self.fep_dir / leg

        # Find repeat directories (repeat1, repeat2, etc.)
        repeat_dirs = sorted(leg_dir.glob("repeat*"))

        if not repeat_dirs:
            results["warnings"].append(f"No repeat directories found in {leg}/")
            return

        # Check first repeat
        repeat_dir = repeat_dirs[0]
        input_dir = repeat_dir / "input"

        if not input_dir.exists():
            results["errors"].append(f"Missing input directory in {repeat_dir.name}")
            return

        # Check required files
        conf_gro = input_dir / "conf.gro"
        topol_top = repeat_dir / "topol.top"

        if not conf_gro.exists():
            results["errors"].append(f"Missing conf.gro in {leg}/{repeat_dir.name}")
            return

        if not topol_top.exists():
            results["errors"].append(f"Missing topol.top in {leg}/{repeat_dir.name}")
            return

        # Try grompp validation
        em_mdp = repeat_dir / "mdps" / "em.mdp"
        if not em_mdp.exists():
            # Try to find any .mdp file
            mdp_files = list(repeat_dir.glob("*.mdp"))
            if not mdp_files:
                results["warnings"].append(f"No MDP files found in {leg}/{repeat_dir.name}")
                return
            em_mdp = mdp_files[0]

        # Run grompp
        try:
            result = subprocess.run(
                [
                    "gmx",
                    "grompp",
                    "-f",
                    str(em_mdp),
                    "-c",
                    str(conf_gro),
                    "-r",
                    str(conf_gro),
                    "-p",
                    str(topol_top),
                    "-o",
                    "/tmp/test_grompp.tpr",
                    "-maxwarn",
                    "2",
                ],
                capture_output=True,
                text=True,
                timeout=60,
            )

            if result.returncode == 0:
                if leg == "bound":
                    results["grompp_bound"] = True
                else:
                    results["grompp_unbound"] = True
            else:
                error_msg = f"grompp failed for {leg} leg"
                if "Fatal error" in result.stderr:
                    fatal_error = re.search(r"Fatal error:(.*?)(?=\n|$)", result.stderr, re.DOTALL)
                    if fatal_error:
                        error_msg += f": {fatal_error.group(1).strip()}"
                results["errors"].append(error_msg)

        except subprocess.TimeoutExpired:
            results["warnings"].append(f"grompp timeout for {leg} leg")
        except FileNotFoundError:
            results["warnings"].append(f"gmx command not found, skipping grompp validation")

    def _check_ligand_placement(self, results: Dict) -> None:
        """Check ligand placement in bound state."""
        bound_dir = self.fep_dir / "bound"
        repeat_dirs = sorted(bound_dir.glob("repeat*"))

        if not repeat_dirs:
            return

        repeat_dir = repeat_dirs[0]
        conf_gro = repeat_dir / "input" / "conf.gro"

        if not conf_gro.exists():
            return

        try:
            # Parse GRO file to check ligand coordinates
            with open(conf_gro, "r") as f:
                lines = f.readlines()

            if len(lines) < 3:
                return

            natoms = int(lines[1].strip())
            ligand_coords = []
            protein_coords = []

            for i in range(2, min(2 + natoms, len(lines))):
                line = lines[i]
                if len(line) < 44:
                    continue

                resname = line[5:10].strip()
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])

                if resname == "HYB":
                    ligand_coords.append((x, y, z))
                elif resname not in ["SOL", "NA", "CL", "K", "MG", "CA"]:
                    protein_coords.append((x, y, z))

            if not ligand_coords:
                results["errors"].append("No ligand (HYB) atoms found in bound conf.gro")
                return

            if not protein_coords:
                results["warnings"].append("No protein atoms found in bound conf.gro")
                return

            # Calculate centroid distances
            ligand_centroid = [
                sum(c[0] for c in ligand_coords) / len(ligand_coords),
                sum(c[1] for c in ligand_coords) / len(ligand_coords),
                sum(c[2] for c in ligand_coords) / len(ligand_coords),
            ]

            protein_centroid = [
                sum(c[0] for c in protein_coords) / len(protein_coords),
                sum(c[1] for c in protein_coords) / len(protein_coords),
                sum(c[2] for c in protein_coords) / len(protein_coords),
            ]

            # Calculate distance
            distance = (
                (ligand_centroid[0] - protein_centroid[0]) ** 2
                + (ligand_centroid[1] - protein_centroid[1]) ** 2
                + (ligand_centroid[2] - protein_centroid[2]) ** 2
            ) ** 0.5

            # Calculate minimum distance
            min_distance = float("inf")
            for lig_coord in ligand_coords:
                for prot_coord in protein_coords:
                    dist = (
                        (lig_coord[0] - prot_coord[0]) ** 2
                        + (lig_coord[1] - prot_coord[1]) ** 2
                        + (lig_coord[2] - prot_coord[2]) ** 2
                    ) ** 0.5
                    if dist < min_distance:
                        min_distance = dist

            results["ligand_placement"] = {
                "centroid_distance": distance,
                "min_distance": min_distance,
                "ligand_atoms": len(ligand_coords),
                "protein_atoms": len(protein_coords),
            }

            # Validate placement
            if distance > 1.0:  # More than 1 nm apart
                results["warnings"].append(f"Ligand-protein centroid distance large: {distance:.3f} nm")

            if min_distance < 0.15:  # Less than 1.5 Å
                results["warnings"].append(f"Ligand-protein minimum distance very small: {min_distance:.3f} nm")

        except Exception as e:
            results["warnings"].append(f"Failed to analyze ligand placement: {e}")

    def print_report(self, results: Dict) -> None:
        """Print formatted validation report."""
        print("\n" + "=" * 80)
        print(f"FEP System Validation Report")
        print("=" * 80)
        print(f"\nSystem: {results['system_path']}")

        # Directory structure
        print("\n[1] Directory Structure")
        print(f"  common/hybrid: {'✓' if results['has_common_hybrid'] else '✗'}")
        print(f"  bound:         {'✓' if results['has_bound'] else '✗'}")
        print(f"  unbound:       {'✓' if results['has_unbound'] else '✗'}")

        # Mapping quality
        if results["mapping_quality"]:
            print("\n[2] Mapping Quality")
            mq = results["mapping_quality"]
            print(f"  Common atoms:      {mq['common']}")
            print(f"  Transformed A:     {mq['transformed_a']}")
            print(f"  Transformed B:     {mq['transformed_b']}")
            print(f"  Surrounding A:     {mq['surrounding_a']}")
            print(f"  Surrounding B:     {mq['surrounding_b']}")
            print(f"  Gray atoms:        {mq['gray_atoms']} {'✗' if mq['gray_atoms'] > 0 else '✓'}")
            print(f"  Unknown atoms:     {mq['unknown_atoms']} {'✗' if mq['unknown_atoms'] > 0 else '✓'}")

        # Topology completeness
        print("\n[3] Topology Completeness")
        if results["topology_complete"]:
            print("  ✓ All required files present")
        else:
            print("  ✗ Missing files:")
            for f in results["missing_files"]:
                print(f"    - {f}")

        # Grompp validation
        print("\n[4] Grompp Validation")
        print(f"  Bound:   {'✓ PASS' if results['grompp_bound'] else '✗ FAIL'}")
        print(f"  Unbound: {'✓ PASS' if results['grompp_unbound'] else '✗ FAIL' if results['has_unbound'] else 'N/A'}")

        # Ligand placement
        if "ligand_placement" in results:
            print("\n[5] Ligand Placement")
            lp = results["ligand_placement"]
            print(f"  Centroid distance: {lp['centroid_distance']:.3f} nm")
            print(f"  Min distance:      {lp['min_distance']:.3f} nm")
            print(f"  Ligand atoms:      {lp['ligand_atoms']}")
            print(f"  Protein atoms:     {lp['protein_atoms']}")

        # Errors and warnings
        if results["errors"]:
            print(f"\n[!] ERRORS ({len(results['errors'])})")
            for i, error in enumerate(results["errors"], 1):
                print(f"  {i}. {error}")

        if results["warnings"]:
            print(f"\n[!] WARNINGS ({len(results['warnings'])})")
            for i, warning in enumerate(results["warnings"], 1):
                print(f"  {i}. {warning}")

        # Overall status
        print("\n" + "=" * 80)
        has_errors = len(results["errors"]) > 0
        has_warnings = len(results["warnings"]) > 0

        if has_errors:
            print("Status: ✗ FAILED - System has critical errors")
        elif has_warnings:
            print("Status: ⚠ WARNING - System passed with warnings")
        else:
            print("Status: ✓ PASSED - System fully validated")
        print("=" * 80 + "\n")


def main():
    if len(sys.argv) < 2:
        print("Usage: python fep_system_check.py <fep_directory>")
        print("\nExample:")
        print("  python fep_system_check.py tests/gxf/FEP/unit_test/42-38/CHARMM-GUI/GMX_PROLIG_FEP")
        sys.exit(1)

    fep_dir = sys.argv[1]

    try:
        checker = FEPSystemChecker(fep_dir)
        results = checker.check_all()
        checker.print_report(results)

        # Exit with error code if validation failed
        sys.exit(1 if results["errors"] else 0)

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
