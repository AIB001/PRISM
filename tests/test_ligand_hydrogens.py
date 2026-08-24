"""Ligand hydrogen validation before and after parameterization.

Structure-based generative models, docking outputs and PDB-extracted ligands
all commonly carry heavy atoms only. Handing such a file to antechamber does
not reliably fail: AM1-BCC aborts on the odd electron count for some molecules,
but for others the pipeline runs to completion and writes a topology with no
hydrogens at all. In a measured run, 19 of 41 heavy-atom-only ligands "built
successfully" and every one of the resulting topologies had zero hydrogens --
including an experimentally determined crystal ligand.
"""

import tempfile
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem
from rdkit.Chem import AllChem

from prism.forcefield.gaff import GAFFForceFieldGenerator


def _write_ligand(directory, smiles, with_hydrogens):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(molecule, randomSeed=20260819)
    if not with_hydrogens:
        molecule = Chem.RemoveHs(molecule)
    path = Path(directory) / ("with_h.sdf" if with_hydrogens else "no_h.sdf")
    Chem.MolToMolFile(molecule, str(path))
    return path


def _generator(ligand_path):
    return GAFFForceFieldGenerator(
        ligand_path=str(ligand_path), output_dir=tempfile.mkdtemp()
    )


# N-methylbenzamide: 9 hydrogens once complete.
DRUG_LIKE = "CNC(=O)c1ccccc1"
# Hexafluorobenzene genuinely has none, and must not be rejected.
NO_HYDROGENS_BY_NATURE = "Fc1c(F)c(F)c(F)c(F)c1F"


def test_heavy_atom_only_ligand_is_rejected_before_parameterization(tmp_path):
    ligand = _write_ligand(tmp_path, DRUG_LIKE, with_hydrogens=False)
    generator = _generator(ligand)

    assert generator.count_missing_hydrogens() == 9
    with pytest.raises(ValueError) as error:
        generator.require_explicit_hydrogens()
    message = str(error.value)
    assert "missing 9 hydrogen" in message
    # The message has to tell the user what to do about it.
    assert "AddHs" in message and "addCoords" in message


def test_complete_ligand_passes(tmp_path):
    ligand = _write_ligand(tmp_path, DRUG_LIKE, with_hydrogens=True)
    generator = _generator(ligand)

    assert generator.count_missing_hydrogens() == 0
    generator.require_explicit_hydrogens()


def test_a_molecule_with_no_hydrogens_of_its_own_is_not_rejected(tmp_path):
    """Hexafluorobenzene is complete without any; the check must not fire."""
    ligand = _write_ligand(tmp_path, NO_HYDROGENS_BY_NATURE, with_hydrogens=False)
    generator = _generator(ligand)

    assert generator.count_missing_hydrogens() == 0
    generator.require_explicit_hydrogens()


def test_check_is_skipped_rather_than_guessed_without_rdkit(tmp_path, monkeypatch):
    """RDKit is optional; without it the question cannot be answered."""
    ligand = _write_ligand(tmp_path, DRUG_LIKE, with_hydrogens=False)
    generator = _generator(ligand)

    import builtins

    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "rdkit" or name.startswith("rdkit."):
            raise ImportError("rdkit disabled for this test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    assert generator.count_missing_hydrogens() is None
    # Unknown must not mean rejected.
    generator.require_explicit_hydrogens()


def test_hydrogen_free_topology_is_rejected(tmp_path):
    """The other half of the failure: a produced topology with no hydrogens."""
    ligand = _write_ligand(tmp_path, DRUG_LIKE, with_hydrogens=True)
    generator = _generator(ligand)

    lig_dir = tmp_path / "LIG.amb2gmx"
    lig_dir.mkdir()
    (lig_dir / "LIG.itp").write_text(
        "[ moleculetype ]\n"
        " LIG  3\n\n"
        "[ atoms ]\n"
        ";   nr  type  resi  res  atom  cgnr charge mass\n"
        "     1    ca     1  LIG    C1     1  -0.10 12.01\n"
        "     2    ca     1  LIG    C2     2  -0.10 12.01\n"
        "     3     n     1  LIG    N1     3  -0.40 14.01\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="no hydrogens"):
        generator.verify_topology_hydrogens(str(lig_dir))


def test_topology_with_hydrogens_is_accepted(tmp_path):
    ligand = _write_ligand(tmp_path, DRUG_LIKE, with_hydrogens=True)
    generator = _generator(ligand)

    lig_dir = tmp_path / "LIG.amb2gmx"
    lig_dir.mkdir()
    (lig_dir / "LIG.itp").write_text(
        "[ atoms ]\n"
        "     1    ca     1  LIG    C1     1  -0.10 12.01\n"
        "     2    ha     1  LIG    H1     2   0.10  1.008\n",
        encoding="utf-8",
    )
    generator.verify_topology_hydrogens(str(lig_dir))
