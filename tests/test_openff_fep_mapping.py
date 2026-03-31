#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenFF ForceFieldGenerator test with end-to-end FEP mapping
Tests OpenFF force field generation and FEP mapping visualization
"""

import pytest
import sys
from pathlib import Path

# Add PRISM to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from prism.forcefield.openff import OpenFFForceFieldGenerator
from prism.fep.visualize.mapping import visualize_mapping_png
from prism.fep.visualize.html import visualize_mapping_html
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.io import read_ligand_from_prism

try:
    from openff.toolkit import Molecule  # noqa: F401

    OPENFF_TOOLKIT_AVAILABLE = True
except ImportError:
    OPENFF_TOOLKIT_AVAILABLE = False


class TestOpenFFWithFEPMapping:
    """Test OpenFF force field generation with FEP mapping"""

    @pytest.fixture(autouse=True)
    def _require_openff_toolkit(self):
        if not OPENFF_TOOLKIT_AVAILABLE:
            pytest.skip("OpenFF toolkit is not available in this environment")

    @pytest.fixture
    def test_data_dir(self):
        """Get test data directory"""
        return Path("tests/gxf/FEP/test/FEP+T4/oMeEtPh-EtPh")

    def test_openff_generation_oMeEtPh(self, test_data_dir, tmp_path):
        """
        Test OpenFF force field generation for oMeEtPh molecule
        """
        mol2_file = test_data_dir / "oMeEtPh.mol2"

        if not mol2_file.exists():
            pytest.skip(f"Test data not found: {mol2_file}")

        print(f"\n{'='*70}")
        print("Testing OpenFF Force Field Generation: oMeEtPh")
        print(f"{'='*70}")

        # Create output directory
        output_dir = tmp_path / "openff_oMeEtPh"
        output_dir.mkdir(exist_ok=True)

        try:
            # Generate OpenFF force field
            generator = OpenFFForceFieldGenerator(
                ligand_path=str(mol2_file),
                output_dir=str(output_dir),
                charge=0,
                forcefield="openff-2.1.0",
                overwrite=True,
            )

            result_dir = generator.run()

            # Verify output files
            required_files = ["LIG.gro", "LIG.itp", "LIG.top", "atomtypes_LIG.itp", "posre_LIG.itp"]

            print(f"\n✓ OpenFF force field generated: {result_dir}")
            print("  Checking output files:")

            for filename in required_files:
                filepath = Path(result_dir) / filename
                assert filepath.exists(), f"Missing required file: {filename}"
                print(f"    ✓ {filename}")

            print(f"\n{'='*70}")
            print("✓ OpenFF generation test passed for oMeEtPh")
            print(f"{'='*70}\n")

        except ImportError as e:
            pytest.skip(f"OpenFF dependencies not available: {e}")

    def test_openff_generation_EtPh(self, test_data_dir, tmp_path):
        """
        Test OpenFF force field generation for EtPh molecule
        """
        mol2_file = test_data_dir / "EtPh.mol2"

        if not mol2_file.exists():
            pytest.skip(f"Test data not found: {mol2_file}")

        print(f"\n{'='*70}")
        print("Testing OpenFF Force Field Generation: EtPh")
        print(f"{'='*70}")

        # Create output directory
        output_dir = tmp_path / "openff_EtPh"
        output_dir.mkdir(exist_ok=True)

        try:
            # Generate OpenFF force field
            generator = OpenFFForceFieldGenerator(
                ligand_path=str(mol2_file),
                output_dir=str(output_dir),
                charge=0,
                forcefield="openff-2.1.0",
                overwrite=True,
            )

            result_dir = generator.run()

            # Verify output files
            required_files = ["LIG.gro", "LIG.itp", "LIG.top", "atomtypes_LIG.itp", "posre_LIG.itp"]

            print(f"\n✓ OpenFF force field generated: {result_dir}")
            print("  Checking output files:")

            for filename in required_files:
                filepath = Path(result_dir) / filename
                assert filepath.exists(), f"Missing required file: {filename}"
                print(f"    ✓ {filename}")

            print(f"\n{'='*70}")
            print("✓ OpenFF generation test passed for EtPh")
            print(f"{'='*70}\n")

        except ImportError as e:
            pytest.skip(f"OpenFF dependencies not available: {e}")

    @pytest.mark.slow
    def test_end_to_end_fep_mapping(self, test_data_dir, tmp_path):
        """
        End-to-end test: OpenFF generation + FEP mapping visualization

        This tests the complete workflow:
        1. Generate OpenFF force fields for both molecules
        2. Visualize FEP mapping between them
        """
        mol2_a = test_data_dir / "oMeEtPh.mol2"
        mol2_b = test_data_dir / "EtPh.mol2"

        if not mol2_a.exists() or not mol2_b.exists():
            pytest.skip(f"Test data not found")

        print(f"\n{'='*70}")
        print("End-to-End Test: OpenFF + FEP Mapping")
        print(f"{'='*70}")

        # Create output directories
        output_dir_a = tmp_path / "openff_oMeEtPh"
        output_dir_b = tmp_path / "openff_EtPh"
        output_dir_a.mkdir(exist_ok=True)
        output_dir_b.mkdir(exist_ok=True)

        try:
            # Step 1: Generate OpenFF force fields
            print("\n[Step 1] Generating OpenFF force fields...")

            print("  → Generating oMeEtPh force field...")
            generator_a = OpenFFForceFieldGenerator(
                ligand_path=str(mol2_a),
                output_dir=str(output_dir_a),
                charge=0,
                forcefield="openff-2.1.0",
                overwrite=True,
            )
            result_dir_a = generator_a.run()
            print(f"    ✓ Generated: {result_dir_a}")

            print("  → Generating EtPh force field...")
            generator_b = OpenFFForceFieldGenerator(
                ligand_path=str(mol2_b),
                output_dir=str(output_dir_b),
                charge=0,
                forcefield="openff-2.1.0",
                overwrite=True,
            )
            result_dir_b = generator_b.run()
            print(f"    ✓ Generated: {result_dir_b}")

            # Step 2: Visualize FEP mapping
            print("\n[Step 2] Visualizing FEP mapping...")

            output_png = tmp_path / "oMeEtPh-EtPh_openff_mapping.png"

            # Read ligand atoms from PRISM-generated files
            itp_a = Path(result_dir_a) / "LIG.itp"
            gro_a = Path(result_dir_a) / "LIG.gro"
            itp_b = Path(result_dir_b) / "LIG.itp"
            gro_b = Path(result_dir_b) / "LIG.gro"

            atoms_a = read_ligand_from_prism(str(itp_a), str(gro_a))
            atoms_b = read_ligand_from_prism(str(itp_b), str(gro_b))

            print(f"    → Loaded {len(atoms_a)} atoms from oMeEtPh")
            print(f"    → Loaded {len(atoms_b)} atoms from EtPh")

            # Create atom mapping
            mapper = DistanceAtomMapper()
            mapping = mapper.map(atoms_a, atoms_b)

            print(f"    → Mapped atoms:")
            print(f"      - Common: {len(mapping.common)}")
            print(f"      - Transformed A: {len(mapping.transformed_a)}")
            print(f"      - Transformed B: {len(mapping.transformed_b)}")
            print(f"      - Surrounding A: {len(mapping.surrounding_a)}")
            print(f"      - Surrounding B: {len(mapping.surrounding_b)}")

            # Visualize mapping (PNG)
            visualize_mapping_png(
                mapping=mapping,
                pdb_a=str(test_data_dir / "oMeEtPh.pdb"),
                pdb_b=str(test_data_dir / "EtPh.pdb"),
                mol2_a=str(mol2_a),
                mol2_b=str(mol2_b),
                output_path=str(output_png),
            )

            assert output_png.exists(), "FEP mapping visualization not generated"
            print(f"    ✓ Mapping visualization saved: {output_png}")

            # Verify file size (should be a valid PNG)
            file_size = output_png.stat().st_size
            assert file_size > 1000, f"PNG file too small ({file_size} bytes)"
            print(f"    ✓ PNG file size: {file_size} bytes")

            # Step 3: Generate HTML visualization
            print("\n[Step 3] Generating HTML visualization...")

            output_html = tmp_path / "oMeEtPh-EtPh_openff_mapping.html"

            visualize_mapping_html(
                mapping=mapping,
                pdb_a=str(test_data_dir / "oMeEtPh.pdb"),
                pdb_b=str(test_data_dir / "EtPh.pdb"),
                mol2_a=str(mol2_a),
                mol2_b=str(mol2_b),
                atoms_a=atoms_a,
                atoms_b=atoms_b,
                output_path=str(output_html),
                title="OpenFF FEP Mapping: oMeEtPh → EtPh",
                ligand_a_name="oMeEtPh",
                ligand_b_name="EtPh",
            )

            assert output_html.exists(), "HTML visualization not generated"
            html_size = output_html.stat().st_size
            print(f"    ✓ HTML visualization saved: {output_html}")
            print(f"    ✓ HTML file size: {html_size} bytes")

            print(f"\n{'='*70}")
            print("✓ End-to-end test passed!")
            print(f"  - OpenFF force fields generated for both molecules")
            print(f"  - FEP mapping visualization created (PNG + HTML)")
            print(f"  - PNG output: {output_png}")
            print(f"  - HTML output: {output_html}")
            print(f"{'='*70}\n")

        except ImportError as e:
            pytest.skip(f"Dependencies not available: {e}")


class TestOpenFFCompatibility:
    """Test OpenFF compatibility with PRISM standards"""

    @pytest.fixture(autouse=True)
    def _require_openff_toolkit(self):
        if not OPENFF_TOOLKIT_AVAILABLE:
            pytest.skip("OpenFF toolkit is not available in this environment")

    def test_openff_output_structure(self, tmp_path):
        """
        Test that OpenFF generates PRISM-compatible output structure
        """
        # Create a minimal but valid MOL2 file for testing (ethane)
        mol2_content = """@<TRIPOS>MOLECULE
ethane
 8 7 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3     1  LIG1       -0.1200
      2 C2          1.5400    0.0000    0.0000 C.3     1  LIG1       -0.1200
      3 H1         -0.3633    1.0272    0.0000 H       1  LIG1        0.0600
      4 H2         -0.3633   -0.5136   -0.8892 H       1  LIG1        0.0600
      5 H3         -0.3633   -0.5136    0.8892 H       1  LIG1        0.0600
      6 H4          1.9033    0.5136   -0.8892 H       1  LIG1        0.0600
      7 H5          1.9033    0.5136    0.8892 H       1  LIG1        0.0600
      8 H6          1.9033   -1.0272    0.0000 H       1  LIG1        0.0600
@<TRIPOS>BOND
     1     1     2    1
     2     1     3    1
     3     1     4    1
     4     1     5    1
     5     2     6    1
     6     2     7    1
     7     2     8    1
"""
        mol2_file = tmp_path / "test.mol2"
        mol2_file.write_text(mol2_content)

        output_dir = tmp_path / "openff_test"
        output_dir.mkdir(exist_ok=True)

        try:
            generator = OpenFFForceFieldGenerator(
                ligand_path=str(mol2_file), output_dir=str(output_dir), charge=0, overwrite=True
            )

            result_dir = generator.run()

            # Check output directory name
            assert "LIG.openff2gmx" in result_dir, "Output directory should follow PRISM naming convention"

            # Check required files
            required_files = ["LIG.gro", "LIG.itp", "LIG.top", "atomtypes_LIG.itp", "posre_LIG.itp"]

            for filename in required_files:
                filepath = Path(result_dir) / filename
                assert filepath.exists(), f"Missing PRISM-standard file: {filename}"

            print("✓ OpenFF output structure is PRISM-compatible")

        except ImportError as e:
            pytest.skip(f"OpenFF dependencies not available: {e}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
