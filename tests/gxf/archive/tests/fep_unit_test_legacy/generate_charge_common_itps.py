"""
Generate ITP files for all charge_common strategies with 39-8 system
"""

from pathlib import Path
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.core.hybrid_topology import HybridTopologyBuilder
from prism.fep.gromacs.itp_builder import ITPBuilder
from prism.fep.io import read_ligand_from_prism


def generate_charge_common_itps():
    """Generate ITP files for all three charge_common strategies"""

    test_dir = Path(__file__).parent / "39-8"
    opls_dir = test_dir / "opls_test_output"

    itp_a = opls_dir / "opls_39" / "LIG.opls2gmx" / "LIG.itp"
    gro_a = opls_dir / "opls_39" / "LIG.opls2gmx" / "LIG.gro"
    itp_b = opls_dir / "opls_8" / "LIG.opls2gmx" / "LIG.itp"
    gro_b = opls_dir / "opls_8" / "LIG.opls2gmx" / "LIG.gro"

    if not itp_a.exists():
        print("OPLS test data not found")
        return

    # Read ligands
    atoms_a = read_ligand_from_prism(str(itp_a), str(gro_a))
    atoms_b = read_ligand_from_prism(str(itp_b), str(gro_b))

    print(f"System: 39 -> 8")
    print(f"Ligand A: {len(atoms_a)} atoms")
    print(f"Ligand B: {len(atoms_b)} atoms")

    # Create mapping
    mapper = DistanceAtomMapper(dist_cutoff=0.6, charge_cutoff=0.05)
    mapping = mapper.map(atoms_a, atoms_b)

    print(f"\nMapping:")
    print(f"  Common: {len(mapping.common)}")
    print(f"  Transformed A: {len(mapping.transformed_a)}")
    print(f"  Transformed B: {len(mapping.transformed_b)}")
    print(f"  Surrounding A: {len(mapping.surrounding_a)}")

    # Create params
    params_a = {"masses": {atom.atom_type: 12.011 if atom.element == "C" else 1.008 for atom in atoms_a}}
    params_b = {"masses": {atom.atom_type: 12.011 if atom.element == "C" else 1.008 for atom in atoms_b}}

    # Generate ITP for each strategy
    strategies = ["mean", "ref", "mut"]

    for strategy in strategies:
        print(f"\n{'='*60}")
        print(f"Generating ITP for charge_common = '{strategy}'")
        print(f"{'='*60}")

        builder = HybridTopologyBuilder(charge_strategy=strategy, charge_reception="surround", recharge_hydrogen=False)
        hybrid_atoms = builder.build(mapping, params_a, params_b, atoms_a, atoms_b)

        # Build ITP file
        itp_builder = ITPBuilder(hybrid_atoms, hybrid_params={})
        output_path = opls_dir / f"hybrid_charge_{strategy}.itp"
        itp_builder.write_itp(str(output_path), molecule_name="HYB")

        print(f"✓ ITP file generated: {output_path}")

        # Verify
        original_charge_a = sum(atom.charge for atom in atoms_a)
        original_charge_b = sum(atom.charge for atom in atoms_b)
        final_charge_a = sum(atom.state_a_charge for atom in hybrid_atoms)
        final_charge_b = sum(
            atom.state_b_charge if atom.state_b_charge is not None else atom.state_a_charge for atom in hybrid_atoms
        )

        print(f"Charge conservation:")
        print(f"  ΔA = {abs(final_charge_a - original_charge_a):.2e}")
        print(f"  ΔB = {abs(final_charge_b - original_charge_b):.2e}")

    print(f"\n{'='*60}")
    print("✓ All ITP files generated!")
    print(f"{'='*60}")


if __name__ == "__main__":
    generate_charge_common_itps()
