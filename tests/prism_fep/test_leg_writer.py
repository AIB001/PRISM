from pathlib import Path

from prism.fep.modeling.leg_writer import LegWriter
from prism.utils.system.topology import _should_vendor_forcefield


def test_copy_and_modify_topology_only_rewrites_molecules_section(tmp_path: Path):
    output_dir = tmp_path / "fep"
    hybrid_dir = output_dir / "common" / "hybrid"
    hybrid_dir.mkdir(parents=True)
    (hybrid_dir / "hybrid.itp").write_text("; hybrid topology\n")

    source_top = tmp_path / "topol.top"
    source_top.write_text(
        "\n".join(
            [
                '#include "amber99sb.ff/forcefield.itp"',
                '#include "LIG.itp"',
                "",
                "[ moleculetype ]",
                "; Name            nrexcl",
                "Protein_chain_A   3",
                "",
                "[ atoms ]",
                "1  opls_001  1  PRO  N  1  -0.3  14.0",
                "",
                "[ system ]",
                "Test system",
                "",
                "[ molecules ]",
                "; Compound        #mols",
                "Protein_chain_A   1",
                "LIG               1",
                "SOL               10",
                "",
            ]
        )
        + "\n"
    )

    target_top = tmp_path / "bound" / "repeat1" / "topol.top"
    target_top.parent.mkdir(parents=True)

    writer = LegWriter(output_dir=output_dir)
    writer.copy_and_modify_topology(source_top, target_top)

    content = target_top.read_text()
    assert content.count("hybrid.itp") == 1
    assert '#include "../../fep/common/hybrid/hybrid.itp"' in content

    before_atoms = content.split("[ atoms ]", 1)[0]
    assert "HYB      1" not in before_atoms

    molecules = writer.extract_molecules_from_topology(target_top)
    assert molecules == [("Protein_chain_A", 1), ("HYB", 1), ("SOL", 10)]


def test_copy_and_modify_topology_preserves_ligand_position_in_molecules_section(tmp_path: Path):
    output_dir = tmp_path / "fep"
    hybrid_dir = output_dir / "common" / "hybrid"
    hybrid_dir.mkdir(parents=True)
    (hybrid_dir / "hybrid.itp").write_text("; hybrid topology\n")

    source_top = tmp_path / "topol.top"
    source_top.write_text(
        "\n".join(
            [
                '#include "charmm36m-mut.ff/forcefield.itp"',
                '#include "LIG.itp"',
                "",
                "[ system ]",
                "Protein in water",
                "",
                "[ molecules ]",
                "; Compound        #mols",
                "LIG               1",
                "Protein_chain_P   1",
                "SOL               10",
                "NA                1",
                "CL                2",
                "",
            ]
        )
        + "\n"
    )

    target_top = tmp_path / "bound" / "repeat1" / "topol.top"
    target_top.parent.mkdir(parents=True)

    writer = LegWriter(output_dir=output_dir)
    writer.copy_and_modify_topology(source_top, target_top)

    molecules = writer.extract_molecules_from_topology(target_top)
    assert molecules == [("HYB", 1), ("Protein_chain_P", 1), ("SOL", 10), ("NA", 1), ("CL", 2)]


def test_copy_topology_itp_files_skips_forcefield_directories(tmp_path: Path):
    source_dir = tmp_path / "source"
    target_dir = tmp_path / "target"
    source_dir.mkdir()
    target_dir.mkdir()

    (source_dir / "protein.itp").write_text("; protein\n")
    ff_dir = source_dir / "charmm36m-mut.ff"
    ff_dir.mkdir()
    (ff_dir / "forcefield.itp").write_text("; ff\n")

    writer = LegWriter(output_dir=tmp_path / "fep")
    writer.copy_topology_itp_files(source_dir, target_dir)

    assert (target_dir / "protein.itp").exists()
    assert not (target_dir / "charmm36m-mut.ff").exists()


def test_should_vendor_forcefield_only_for_non_installed_paths(tmp_path: Path):
    ff_dir = tmp_path / "custom" / "amber99sb.ff"
    ff_dir.mkdir(parents=True)

    assert _should_vendor_forcefield({"path": str(ff_dir), "source": "GMXLIB"}, tmp_path) is False
    assert _should_vendor_forcefield({"path": str(ff_dir), "source": "GROMACS installation"}, tmp_path) is False
    assert _should_vendor_forcefield({"path": str(ff_dir), "source": "custom"}, tmp_path) is True
    assert (
        _should_vendor_forcefield({"path": str(tmp_path / "amber99sb.ff"), "source": "current directory"}, tmp_path)
        is False
    )
