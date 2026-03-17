"""
Generate HTML visualization for charge_common='ref' strategy with 39-8 system
"""

from pathlib import Path
from prism.fep.core.mapping import DistanceAtomMapper
from prism.fep.core.hybrid_topology import HybridTopologyBuilder
from prism.fep.io import read_ligand_from_prism
from prism.fep.visualize.html import visualize_mapping_html
from prism.fep.config import FEPConfig


def generate_ref_html():
    """Generate HTML visualization for ref strategy"""

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

    # Build hybrid topology with ref strategy
    params_a = {"masses": {atom.atom_type: 12.011 if atom.element == "C" else 1.008 for atom in atoms_a}}
    params_b = {"masses": {atom.atom_type: 12.011 if atom.element == "C" else 1.008 for atom in atoms_b}}

    # Load configuration using FEPConfig
    fep_config = FEPConfig(str(test_dir))
    mapping_params = fep_config.get_mapping_params()

    builder = HybridTopologyBuilder(
        charge_strategy=mapping_params.get("charge_common", "mean"),
        charge_reception=mapping_params.get("charge_reception", "surround"),
        recharge_hydrogen=False,
    )
    hybrid_atoms = builder.build(mapping, params_a, params_b, atoms_a, atoms_b)

    # Generate HTML visualization
    output_html = opls_dir / "39-8_ref_mapping.html"

    # Get HTML config from FEPConfig
    html_config = fep_config.get_html_config()

    # Add case path for HTML display
    html_config["case"] = {"path": str(test_dir.relative_to(Path.cwd()))}

    visualize_mapping_html(
        mapping,
        pdb_a=str(test_dir / "39.pdb"),
        pdb_b=str(test_dir / "8.pdb"),
        atoms_a=atoms_a,
        atoms_b=atoms_b,
        output_path=str(output_html),
        title="FEP Mapping: 39 → 8 (charge_common='ref')",
        ligand_a_name="39",
        ligand_b_name="8",
        config=html_config,
    )

    print(f"\n✓ HTML visualization generated: {output_html}")


if __name__ == "__main__":
    generate_ref_html()
