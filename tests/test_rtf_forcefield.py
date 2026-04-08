from pathlib import Path

from prism.forcefield.rtf import RTFForceFieldGenerator


def test_rtf_generator_uses_prm_nonbonded_parameters(tmp_path):
    fixture_dir = Path("tests/gxf/FEP/unit_test/p38-19-24/input/24")
    generator = RTFForceFieldGenerator(
        rtf_file=fixture_dir / "24.rtf",
        prm_file=fixture_dir / "24.prm",
        pdb_file=fixture_dir / "24.pdb",
        output_dir=tmp_path,
        overwrite=True,
    )

    lig_dir = Path(generator.run())
    atomtypes = (lig_dir / "atomtypes_LIG.itp").read_text()
    lines = {
        line.split()[0]: line.split() for line in atomtypes.splitlines() if line and not line.startswith(("[", ";"))
    }

    assert lines["FGR1"][1] == "FGR1"
    assert lines["FGR1"][2] == "18.99800"
    assert lines["FGR1"][5] == "3.02906e-01"
    assert lines["FGR1"][6] == "5.02080e-01"
    assert lines["HGP1"][2] == "1.00800"
