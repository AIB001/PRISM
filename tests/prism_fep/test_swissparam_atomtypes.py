from pathlib import Path

from prism.forcefield.swissparam import SwissParamForceFieldGenerator


def test_generate_atomtypes_file_uses_reference_ligand_itp_parameters(tmp_path: Path):
    ligand = tmp_path / "ligand.mol2"
    ligand.write_text("@<TRIPOS>MOLECULE\nLIG\n")

    generator = SwissParamForceFieldGenerator(str(ligand), str(tmp_path), "both", overwrite=True)
    Path(generator.lig_ff_dir).mkdir(parents=True)
    (Path(generator.lig_ff_dir) / "ligand.itp").write_text(
        "\n".join(
            [
                "[ atomtypes ]",
                "; name at.num mass charge ptype sigma epsilon",
                "CB      6   12.0110  0.0  A  0.355005  0.292880",
                "HCMM    1    1.0079  0.0  A  0.235197  0.092048",
                "",
                "[ atoms ]",
                "; nr type resnr resid atom cgnr charge mass",
                "1 CB   1 BBE1 C1   1  0.0 12.0110",
                "2 HCMM 1 BBE1 H1   2  0.0  1.0079",
                "",
            ]
        )
        + "\n"
    )

    rtf_data = {
        "atoms": {
            "C1": {"type": "C321", "charge": 0.0},
            "H1": {"type": "HGA3", "charge": 0.0},
        }
    }

    generator._generate_atomtypes_file(rtf_data)
    content = (Path(generator.lig_ff_dir) / "atomtypes_LIG.itp").read_text()

    assert "C321    C321  12.01100" in content
    assert "HGA3    HGA3   1.00790" in content
    assert "2.351970e-01" in content
    assert "9.204800e-02" in content
