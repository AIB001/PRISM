"""
Generate HTML visualization for any FEP test case using FEPConfig

Usage:
    python generate_fep_html.py <test_dir>

Example:
    python generate_fep_html.py tests/gxf/FEP/unit_test/39-8
"""

import sys
from pathlib import Path
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.io import read_ligand_from_prism
from prism.fep.visualize.html import visualize_mapping_html
from prism.fep.config import FEPConfig


def generate_html(test_dir: str):
    """Generate HTML visualization for FEP test case"""

    test_dir = Path(test_dir)
    if not test_dir.exists():
        print(f"Error: Test directory not found: {test_dir}")
        return

    # Load configuration using FEPConfig
    fep_config = FEPConfig(str(test_dir))
    mapping_params = fep_config.get_mapping_params()
    ff_type = fep_config.get_forcefield_type()

    print(f"Generating HTML for: {test_dir.name}")
    print(f"  Force Field: {ff_type}")
    print(f"  Mapping params: {mapping_params}")

    # Determine output directory based on force field type
    ff_type_lower = ff_type.lower()
    if "gaff" in ff_type_lower:
        output_dir = test_dir / "gaff_test_output"
    elif "openff" in ff_type_lower:
        output_dir = test_dir / "openff_test_output"
    elif "opls" in ff_type_lower:
        output_dir = test_dir / "opls_test_output"
    elif "cgenff" in ff_type_lower:
        output_dir = test_dir / "opls_test_output"  # CGenFF uses OPLS output
    else:
        output_dir = test_dir / f"{ff_type_lower}_test_output"

    # Try to find ITP/GRO files
    itp_a = None
    gro_a = None
    itp_b = None
    gro_b = None

    # Check for PRISM output
    if output_dir.exists():
        # Look for ligand-specific subdirectories (e.g., opls_39, opls_8, gaff_ligand_a, etc.)
        for subdir in sorted(output_dir.iterdir()):
            if subdir.is_dir():
                # Try different GROMACS directory names
                possible_names = [
                    "LIG.opls2gmx",
                    "LIG.gaff2gmx",
                    "LIG.amb2gmx",
                    "LIG.openff2gmx",
                    "LIG.cgenff2gmx",
                    "LIG.gromacs",
                ]
                gromacs_dir = None
                for name in possible_names:
                    test_path = subdir / name
                    if test_path.exists():
                        gromacs_dir = test_path
                        break

                if gromacs_dir and gromacs_dir.exists():
                    itp_file = gromacs_dir / "LIG.itp"
                    gro_file = gromacs_dir / "LIG.gro"
                    if itp_file.exists() and gro_file.exists():
                        if itp_a is None:
                            itp_a = itp_file
                            gro_a = gro_file
                        elif itp_b is None:
                            itp_b = itp_file
                            gro_b = gro_file

    # Fallback: check for manually placed ITP files
    if itp_a is None or not itp_a.exists():
        itp_a = output_dir / "ligand_a.itp"
        gro_a = output_dir / "ligand_a.gro"
    if itp_b is None or not itp_b.exists():
        itp_b = output_dir / "ligand_b.itp"
        gro_b = output_dir / "ligand_b.gro"

    if not itp_a or not itp_a.exists():
        print(f"Error: Cannot find ITP file for ligand A")
        print(f"  Searched in: {output_dir}")
        return

    if not itp_b or not itp_b.exists():
        print(f"Error: Cannot find ITP file for ligand B")
        print(f"  Searched in: {output_dir}")
        return

    # Read ligands
    atoms_a = read_ligand_from_prism(str(itp_a), str(gro_a))
    atoms_b = read_ligand_from_prism(str(itp_b), str(gro_b))

    print(f"  Ligand A: {len(atoms_a)} atoms")
    print(f"  Ligand B: {len(atoms_b)} atoms")

    # Create mapping
    mapper = DistanceAtomMapper(**mapping_params)
    mapping = mapper.map(atoms_a, atoms_b)

    print(
        f"  Mapping: Common={len(mapping.common)}, "
        f"Transformed=A{len(mapping.transformed_a)}+B{len(mapping.transformed_b)}, "
        f"Surrounding=A{len(mapping.surrounding_a)}"
    )

    # Get PDB files for visualization
    pdb_files = sorted(test_dir.glob("*.pdb"))
    if len(pdb_files) < 2:
        print(f"\nError: Need at least 2 PDB files for visualization")
        print(f"  Found {len(pdb_files)} PDB file(s) in {test_dir}")
        print(f"  Please place ligand PDB files in the test directory")
        return

    pdb_a = str(pdb_files[0])
    pdb_b = str(pdb_files[1])
    print(f"  PDB A: {Path(pdb_a).name}")
    print(f"  PDB B: {Path(pdb_b).name}")

    # Generate HTML
    output_html = output_dir / f"{test_dir.name}_fep_mapping.html"

    # Get HTML config from FEPConfig
    html_config = fep_config.get_html_config()

    # Add case path for HTML display
    test_dir_abs = test_dir.resolve()
    cwd_abs = Path.cwd().resolve()
    try:
        rel_path = test_dir_abs.relative_to(cwd_abs)
    except ValueError:
        rel_path = test_dir_abs  # Use absolute path if not relative
    html_config["case"] = {"path": str(rel_path)}

    visualize_mapping_html(
        mapping,
        pdb_a=pdb_a,
        pdb_b=pdb_b,
        atoms_a=atoms_a,
        atoms_b=atoms_b,
        output_path=str(output_html),
        title=f"FEP Mapping: {test_dir.name}",
        ligand_a_name="Ligand A",
        ligand_b_name="Ligand B",
        config=html_config,
    )

    print(f"\n✓ HTML saved to: {output_html}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_fep_html.py <test_dir>")
        print("\nExample:")
        print("  python generate_fep_html.py tests/gxf/FEP/unit_test/39-8")
        print("\nAvailable test cases:")
        test_dirs = [
            "tests/gxf/FEP/unit_test/39-8",
            "tests/gxf/FEP/unit_test/42-38",
            "tests/gxf/FEP/unit_test/oMeEtPh-EtPh",
        ]
        for d in test_dirs:
            p = Path(d)
            if p.exists():
                print(f"  - {d}")
        sys.exit(1)

    test_dir = sys.argv[1]
    generate_html(test_dir)
