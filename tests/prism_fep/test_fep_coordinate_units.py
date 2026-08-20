import numpy as np

from prism.fep.coordinates import angstrom_to_nm, nm_to_angstrom
from prism.fep.core.mapping import Atom
from prism.fep.io import (
    _load_structure_coordinates,
    read_mol2_atoms,
    read_rtf_for_fep,
    write_ligand_to_pdb,
)


def test_load_structure_coordinates_normalizes_pdb_to_nm(tmp_path):
    pdb = tmp_path / "lig.pdb"
    pdb.write_text("ATOM      1  C1  LIG A   1      12.340  23.450  34.560  1.00  0.00           C\nEND\n")

    coords_by_id, coords_by_name, elements_by_id = _load_structure_coordinates(str(pdb))

    expected = np.array([1.234, 2.345, 3.456])
    assert np.allclose(coords_by_id[1], expected)
    assert np.allclose(coords_by_name["C1"], expected)
    assert elements_by_id[1] == "C"


def test_coordinate_unit_helpers_round_trip_scalars_and_vectors():
    assert angstrom_to_nm(12.34) == 1.234
    assert nm_to_angstrom(1.234) == 12.34
    assert np.allclose(angstrom_to_nm((12.34, 23.45, 34.56)), (1.234, 2.345, 3.456))
    assert np.allclose(nm_to_angstrom(np.array([1.234, 2.345, 3.456])), np.array([12.34, 23.45, 34.56]))


def test_load_structure_coordinates_preserves_gro_nm(tmp_path):
    gro = tmp_path / "lig.gro"
    gro.write_text(
        "Test ligand\n" "1\n" "    1LIG     C1    1   1.234   2.345   3.456\n" "   4.00000   4.00000   4.00000\n"
    )

    coords_by_id, coords_by_name, _ = _load_structure_coordinates(str(gro))
    expected = np.array([1.234, 2.345, 3.456])
    assert np.allclose(coords_by_id[1], expected)
    assert np.allclose(coords_by_name["C1"], expected)


def test_read_mol2_atoms_converts_angstrom_to_nm(tmp_path):
    mol2 = tmp_path / "lig.mol2"
    mol2.write_text(
        "@<TRIPOS>MOLECULE\n"
        "LIG\n"
        "1 0 0 0 0\n"
        "SMALL\n"
        "NO_CHARGES\n"
        "@<TRIPOS>ATOM\n"
        "      1 C1         12.3400   23.4500   34.5600 C.ar       1 LIG       -0.1200\n"
        "@<TRIPOS>BOND\n"
    )

    atoms = read_mol2_atoms(str(mol2))

    assert len(atoms) == 1
    assert np.allclose(atoms[0].coord, np.array([1.234, 2.345, 3.456]))


def test_read_rtf_for_fep_converts_pdb_angstrom_to_nm(tmp_path):
    rtf = tmp_path / "lig.rtf"
    pdb = tmp_path / "lig.pdb"
    rtf.write_text("ATOM C1 CG2R61 -0.1200\n")
    pdb.write_text("HETATM    1  C1  LIG A   1      12.340  23.450  34.560  1.00  0.00           C\nEND\n")

    atoms = read_rtf_for_fep(str(rtf), str(pdb))

    assert len(atoms) == 1
    assert np.allclose(atoms[0].coord, np.array([1.234, 2.345, 3.456]))


def test_write_ligand_to_pdb_converts_nm_to_angstrom(tmp_path):
    output = tmp_path / "lig.pdb"
    atoms = [
        Atom(
            name="C1",
            element="C",
            coord=np.array([1.234, 2.345, 3.456]),
            charge=-0.12,
            atom_type="ca",
            index=1,
        )
    ]

    write_ligand_to_pdb(atoms, str(output))
    pdb_line = output.read_text().splitlines()[1]

    assert float(pdb_line[30:38]) == 12.340
    assert float(pdb_line[38:46]) == 23.450
    assert float(pdb_line[46:54]) == 34.560
