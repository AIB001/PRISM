"""
Generate OPLS force field for 39-8 test case
"""

from pathlib import Path
from prism.forcefield.opls_aa import OPLSAAForceFieldGenerator
from prism.fep.io import read_ligand_from_prism
from prism.fep.core.mapping import DistanceAtomMapper


def generate_39_8_opls():
    """Generate OPLS data for 39-8 system"""

    test_dir = Path(__file__).parent / "39-8"
    output_dir = test_dir / "opls_test_output"
    output_dir.mkdir(exist_ok=True)

    # Input files (use PDB instead of mol2 for correct coordinates)
    mol2_a = test_dir / "39.pdb"
    mol2_b = test_dir / "8.pdb"

    if not mol2_a.exists():
        print(f"Input file not found: {mol2_a}")
        return

    print("=" * 60)
    print("OPLS-AA FEP Mapping Test: 39 → 8")
    print("=" * 60)

    # Generate OPLS for ligand A (39)
    print("\n[Step 1] Generating OPLS-AA force fields...")
    output_dir_a = output_dir / "opls_39"

    if not output_dir_a.exists():
        print(f"  → Generating 39 force field...")
        opls_gen_a = OPLSAAForceFieldGenerator(
            ligand_path=str(mol2_a), output_dir=str(output_dir_a), charge_model="cm1a", align_to_input=True
        )
        opls_gen_a.run()
    else:
        print(f"  ✓ 39 force field already exists")

    # Generate OPLS for ligand B (8)
    output_dir_b = output_dir / "opls_8"

    if not output_dir_b.exists():
        print(f"  → Generating 8 force field...")
        opls_gen_b = OPLSAAForceFieldGenerator(
            ligand_path=str(mol2_b), output_dir=str(output_dir_b), charge_model="cm1a", align_to_input=True
        )
        opls_gen_b.run()
    else:
        print(f"  ✓ 8 force field already exists")

    # Read ligands
    print("\n[Step 2] Reading ligands...")
    itp_a = output_dir_a / "LIG.opls2gmx" / "LIG.itp"
    gro_a = output_dir_a / "LIG.opls2gmx" / "LIG.gro"
    itp_b = output_dir_b / "LIG.opls2gmx" / "LIG.itp"
    gro_b = output_dir_b / "LIG.opls2gmx" / "LIG.gro"

    atoms_a = read_ligand_from_prism(str(itp_a), str(gro_a))
    atoms_b = read_ligand_from_prism(str(itp_b), str(gro_b))

    print(f"  → Ligand A (39): {len(atoms_a)} atoms")
    print(f"  → Ligand B (8): {len(atoms_b)} atoms")

    # Create mapping
    print("\n[Step 3] Creating mapping...")
    mapper = DistanceAtomMapper(dist_cutoff=0.6, charge_cutoff=0.05)
    mapping = mapper.map(atoms_a, atoms_b)

    print(f"  Common: {len(mapping.common)}")
    print(f"  Transformed A: {len(mapping.transformed_a)}")
    print(f"  Transformed B: {len(mapping.transformed_b)}")
    print(f"  Surrounding A: {len(mapping.surrounding_a)}")

    print("\n" + "=" * 60)
    print("✅ OPLS data generation complete!")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print(f"\nNow you can run test_charge_common_39_8.py")


if __name__ == "__main__":
    generate_39_8_opls()
