"""Chemical counterexamples for the structural peptide exclusion rule."""

import pytest

Chem = pytest.importorskip("rdkit.Chem")
from prism.design.curation.peptides import classify_peptide


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCC(=O)NCC(=O)O", "peptide"),  # Gly-Gly
        ("N[C@@H](C)C(=O)NCC(=O)N[C@@H](C)C(=O)O", "peptide"),
        ("O=C1NCC(=O)NC1", "peptide"),  # cyclic Gly-Gly
        ("NCC(=O)N1CCC[C@H]1C(=O)O", "peptide"),  # proline
        ("N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O", "peptide"),  # gamma-Glu
        ("CNCC(=O)N(C)CC(=O)O", "peptide"),  # peptoid
        ("CC(C)(C)OC(=O)NCC(=O)NCC(=O)O", "peptide"),  # protected Gly-Gly
        ("O=C(OCc1ccccc1)NCC(=O)NCC(=O)O", "peptide"),  # Cbz protection
        ("CC(=O)NCC(=O)O", "nonpeptide"),  # N-acetylglycine: ONE residue
        ("NCC(=O)O", "nonpeptide"),  # glycine alone
        ("CC(=O)Nc1ccc(NC(C)=O)cc1", "nonpeptide"),  # two ordinary amides
        ("NC(=O)N", "nonpeptide"),  # urea
        ("O=C1NCCC1", "peptide_like_review"),  # lactam -> one gamma-AA, not proven peptide
    ],
)
def test_virtual_hydrolysis_distinguishes_peptides(smiles, expected):
    assert classify_peptide(Chem.MolFromSmiles(smiles))["status"] == expected


def test_fragment_evidence_and_missing_graph():
    result = classify_peptide(Chem.MolFromSmiles("NCC(=O)NCC(=O)O"))
    assert len(result["peptide_bonds"]) == 1
    assert result["largest_connected_amino_acid_group"] == 2
    assert len(result["hydrolysis_fragments"]) == 2
    assert all(f["alpha_amino_acid_backbones"] == 1 for f in result["hydrolysis_fragments"])
    assert classify_peptide(None)["status"] == "unknown"
