from pathlib import Path

import numpy as np

from prism.forcefield.common.io import (
    detect_coordinate_file_format,
    read_gro_coordinates,
    read_mol2_coordinates,
    read_pdb_coordinates,
    write_gro_coordinates,
)


def test_read_gro_coordinates_parses_atom_name_index_and_coords(tmp_path: Path):
    gro = tmp_path / "ligand.gro"
    gro.write_text(
        "\n".join(
            [
                "Test GRO",
                "    2",
                "    1LIG     C1    1   0.123   0.456   0.789",
                "    1LIG     H2    2   0.100   0.200   0.300",
                "   1.00000   1.00000   1.00000",
            ]
        )
        + "\n"
    )

    coords = read_gro_coordinates(gro)

    assert len(coords) == 2
    assert coords[0][0] == 1
    assert coords[0][1] == "C1"
    assert np.allclose(coords[0][2], np.array([1.23, 4.56, 7.89]))
    assert coords[1][0] == 2
    assert coords[1][1] == "H2"
    assert np.allclose(coords[1][2], np.array([1.0, 2.0, 3.0]))


def test_write_gro_coordinates_roundtrip_preserves_values(tmp_path: Path):
    gro = tmp_path / "written.gro"
    source = [
        (1, "CA", np.array([1.0, 2.0, 3.0])),
        (2, "HB2", np.array([4.5, 5.5, 6.5])),
    ]

    write_gro_coordinates(gro, source, residue_name="LIG", title="Roundtrip")
    roundtrip = read_gro_coordinates(gro)

    assert [item[0] for item in roundtrip] == [1, 2]
    assert [item[1] for item in roundtrip] == ["CA", "HB2"]
    assert np.allclose(roundtrip[0][2], source[0][2])
    assert np.allclose(roundtrip[1][2], source[1][2])


def test_read_pdb_coordinates_preserves_fixed_column_alignment(tmp_path: Path):
    pdb = tmp_path / "ligand.pdb"
    pdb.write_text(
        "\n".join(
            [
                "HETATM    1  C1  LIG A   1      12.345  23.456  34.567  1.00  0.00           C",
                "HETATM   12  O2  LIG A   1      -1.000   0.500   9.999  1.00  0.00           O",
                "END",
            ]
        )
        + "\n"
    )

    coords = read_pdb_coordinates(pdb)

    assert coords[0][0] == 1
    assert coords[0][1] == "C1"
    assert np.allclose(coords[0][2], np.array([12.345, 23.456, 34.567]))
    assert coords[1][0] == 12
    assert coords[1][1] == "O2"
    assert np.allclose(coords[1][2], np.array([-1.0, 0.5, 9.999]))


def test_read_mol2_coordinates_reads_atom_block(tmp_path: Path):
    mol2 = tmp_path / "ligand.mol2"
    mol2.write_text(
        "\n".join(
            [
                "@<TRIPOS>MOLECULE",
                "LIG",
                " 2 0 0 0 0",
                "SMALL",
                "NO_CHARGES",
                "@<TRIPOS>ATOM",
                "      1 C1         1.0000    2.0000    3.0000 C.3    1 LIG       0.0000",
                "      2 O2         4.0000    5.0000    6.0000 O.2    1 LIG       0.0000",
                "@<TRIPOS>BOND",
            ]
        )
        + "\n"
    )

    coords = read_mol2_coordinates(mol2)

    assert len(coords) == 2
    assert coords[0][0] == 1
    assert coords[0][1] == "C1"
    assert np.allclose(coords[0][2], np.array([1.0, 2.0, 3.0]))
    assert coords[1][0] == 2
    assert coords[1][1] == "O2"
    assert np.allclose(coords[1][2], np.array([4.0, 5.0, 6.0]))


def test_detect_coordinate_file_format_by_extension_and_content(tmp_path: Path):
    gro = tmp_path / "coords.dat"
    gro.write_text("Title\n    1\n    1LIG     C1    1   0.123   0.456   0.789\n   1.00000   1.00000   1.00000\n")
    pdb = tmp_path / "coords.pdb"
    pdb.write_text("HETATM    1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00           C\n")
    mol2 = tmp_path / "coords.mol2"
    mol2.write_text("@<TRIPOS>MOLECULE\nLIG\n")

    assert detect_coordinate_file_format(gro) == "gro"
    assert detect_coordinate_file_format(pdb) == "pdb"
    assert detect_coordinate_file_format(mol2) == "mol2"
