"""Scientifically consequential invariants of M1; no pretrained weights needed."""

import copy

import pytest

torch = pytest.importorskip("torch")
Chem = pytest.importorskip("rdkit.Chem")
from torch import nn

from prism.design.curation.prepare import functional_groups
from prism.design.m1.conditioning import NativeTextConditioner, TextAdapter, install_lora
from prism.design.m1.data import chemical_keys, grouped_split, validate_cut, validate_split
from prism.design.m1.losses import alignment_set_loss, x1_losses
from prism.design.m1.frame import core_frames
from prism.design.m1.native import core_center_prior, verify_preprocessed_molecule
from prism.design.m1.train import EpochStream, lr_factor, restore_rng, rng_state, seed_all


def test_phosphate_matches_count_one_chemical_centre():
    g = functional_groups(Chem.MolFromSmiles("COP(=O)(O)O"))
    assert len(g["phosphate_or_phosphonate"]) == 1
    assert len(g["phosphate_or_phosphonate"][0]) == 5
    assert len(functional_groups(Chem.MolFromSmiles("O=P(O)(O)OP(=O)(O)O"))["phosphate_or_phosphonate"]) == 2


def test_adapter_starts_exactly_at_baseline_then_gets_text_gradients():
    torch.manual_seed(5)
    a = TextAdapter(16, 8)
    inv, text = torch.randn(2, 5, 16), torch.randn(2, 3, 8)
    mask, valid = torch.ones(2, 3, dtype=torch.bool), torch.ones(2, 5, dtype=torch.bool)
    before, _ = a(inv, text, mask, valid)
    assert torch.equal(before, inv)
    (before.square().mean()).backward()
    assert a.film.weight.grad.abs().sum() > 0
    assert a.attention.out_proj.weight.grad.abs().sum() > 0
    optimizer = torch.optim.SGD(a.parameters(), lr=0.1)
    optimizer.step()
    optimizer.zero_grad()
    result, _ = a(inv, text, mask, valid)
    other, _ = a(inv, -text, mask, valid)
    assert not torch.allclose(result, other)
    result.square().mean().backward()
    assert a.text_projection.weight.grad.abs().sum() > 0


def test_null_dropout_is_independent_of_dropped_text_and_padding():
    a = TextAdapter(16, 8)
    nn.init.normal_(a.film.weight)
    nn.init.normal_(a.attention.out_proj.weight)
    x, t = torch.randn(1, 4, 16), torch.randn(1, 3, 8)
    valid = torch.tensor([[True, True, False, False]])
    mask = torch.tensor([[True, True, False]])
    out, _ = a(x, t, mask, valid, torch.tensor([True]))
    other, _ = a(x, t * 900, mask, valid, torch.tensor([True]))
    assert torch.equal(out, other) and torch.isfinite(out).all()
    assert torch.equal(out[:, 2:], x[:, 2:])
    padded = torch.cat((t, torch.full((1, 2, 8), 1e4)), 1)
    expected, _ = a(x, t, mask, valid)
    actual, _ = a(x, padded, torch.cat((mask, torch.zeros(1, 2, dtype=torch.bool)), 1), valid)
    assert torch.allclose(actual, expected, atol=1e-6)


def test_explicit_attachment_survives_text_dropout_and_rejects_padding():
    a = TextAdapter(16, 8)
    a.attachment_vector.data.fill_(0.5)
    x, t = torch.randn(1, 4, 16), torch.randn(1, 3, 8)
    valid = torch.tensor([[True, True, True, False]])
    anchor = torch.tensor([[False, True, False, False]])
    y, _ = a(x, t, torch.ones(1, 3, dtype=torch.bool), valid, torch.tensor([True]), anchor)
    assert torch.allclose(y[:, 1], x[:, 1] + 0.5)
    assert torch.equal(y[:, [0, 2, 3]], x[:, [0, 2, 3]])
    with pytest.raises(ValueError):
        a(x, t, torch.ones(1, 3, dtype=torch.bool), valid, attachment_mask=torch.tensor([[False, False, False, True]]))


def toy_loss_inputs():
    n, classes = 4, 3
    target = {"coords": torch.zeros(1, n, 3), "mask": torch.tensor([[1, 1, 1, 0]], dtype=torch.bool)}
    target.update(
        {
            k: torch.nn.functional.one_hot(
                torch.zeros((1, n, n) if k == "bonds" else (1, n), dtype=torch.long), classes
            ).float()
            for k in ("atomics", "charges", "bonds")
        }
    )
    pred = {k: torch.zeros_like(v, requires_grad=True) for k, v in target.items() if k != "mask"}
    fixed = torch.tensor([[1, 1, 0, 0]], dtype=torch.bool)
    return pred, target, fixed


def test_loss_includes_boundary_bonds_excludes_fixed_diagonal_and_padding():
    p, t, f = toy_loss_inputs()
    loss, _ = x1_losses(p, t, f)
    loss.backward()
    assert p["bonds"].grad[0, 0, 2].abs().sum() > 0  # F-U
    assert p["bonds"].grad[0, 0, 1].abs().sum() == 0  # F-F
    assert p["bonds"].grad[0, 2, 2].abs().sum() == 0  # self
    assert p["bonds"].grad[0, 2, 3].abs().sum() == 0  # pad
    assert p["atomics"].grad[0, :2].abs().sum() == 0
    assert p["atomics"].grad[0, 2].abs().sum() > 0
    p2 = {k: v.detach().clone() for k, v in p.items()}
    p2["coords"][0, :2] = 10000
    p2["coords"][0, 3] = -10000
    assert torch.equal(x1_losses(p2, t, f)[0], loss.detach())
    p2["bonds"][0, 0, 2, 1] = 10
    assert x1_losses(p2, t, f)[0] > loss.detach()


def test_no_variable_or_padding_in_fixed_mask_fails():
    p, t, _ = toy_loss_inputs()
    with pytest.raises(ValueError):
        x1_losses(p, t, torch.tensor([[1, 1, 1, 0]]))
    with pytest.raises(ValueError):
        x1_losses(p, t, torch.tensor([[1, 0, 0, 1]]))


def test_alignment_unknown_not_negative_symmetry_positive_set():
    a = torch.ones(1, 3, 2, requires_grad=True)
    positive = torch.tensor([[[1, 0], [1, 0], [0, 0]]], dtype=torch.bool)
    valid = torch.tensor([[True, True, False]])
    observed = torch.tensor([[True, False]])
    assert abs(float(alignment_set_loss(a, positive, observed, valid))) < 1e-6
    empty = alignment_set_loss(a, positive, torch.zeros_like(observed), valid)
    empty.backward()
    assert a.grad.abs().sum() == 0
    with pytest.raises(ValueError):
        alignment_set_loss(a, positive, torch.ones_like(observed), valid)


def record(ident, chemical, scaffold, protein):
    return {
        "record_id": ident,
        "complex_group": ident,
        "split_keys": {"chemical_identity": chemical, "scaffold": scaffold, "protein_groups": [protein]},
    }


def test_connected_split_transitive_leakage_and_independent_acyclic_groups():
    rows = [record("a", "A", "S", "P"), record("b", "B", "S", "Q"), record("c", "C", "T", "Q")]
    rows += [record(str(i), str(i), str(i), str(i)) for i in range(20)]
    split = grouped_split(rows)
    validate_split(rows, split)
    assert len({split["assignments"][key]["split"] for key in "abc"}) == 1
    bad = copy.deepcopy(split)
    bad["assignments"]["b"]["split"] = "test" if bad["assignments"]["a"]["split"] != "test" else "train"
    with pytest.raises(ValueError):
        validate_split(rows, bad)
    with pytest.raises(ValueError):
        grouped_split(rows[:3])
    assert (
        chemical_keys(Chem.MolFromSmiles("CCCC"))["scaffold"] != chemical_keys(Chem.MolFromSmiles("CCCO"))["scaffold"]
    )


def test_split_group_normalization_does_not_modify_training_molecule():
    a = Chem.MolFromSmiles("N[C@@H](C)C(=O)O")
    before = Chem.MolToSmiles(a)
    b = Chem.MolFromSmiles("[NH3+][C@H](C)C(=O)[O-]")
    assert chemical_keys(a) == chemical_keys(b)
    assert Chem.MolToSmiles(a) == before


def test_prior_is_independent_of_withheld_target_coordinates():
    from types import SimpleNamespace

    class Mol(SimpleNamespace):
        def _copy_with(self, **changes):
            return Mol(**dict(vars(self), **changes))

    f = Mol(coords=torch.randn(5, 3))
    target = Mol(coords=torch.randn(5, 3))
    mask = torch.tensor([True, True, False, False, False])
    first = core_center_prior(None, f, target, mask, mask)
    target.coords[~mask] = 1e9
    second = core_center_prior(None, f, target, mask, mask)
    assert torch.equal(first.coords, second.coords)
    assert torch.allclose(first.coords[~mask].mean(0), target.coords[mask].mean(0), atol=1e-6)


def test_core_frame_is_proper_rotation_covariant_and_has_no_target_leak():
    xyz = torch.tensor([[[1.0, 2.0, 0.0], [2.0, 3.0, 1.0], [1.0, 4.0, 1.0], [6.0, 6.0, 7.0]]])
    fixed = torch.tensor([[True, True, True, False]])
    b = core_frames(xyz, fixed)
    rot = torch.tensor([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    rotated = xyz @ rot
    assert torch.allclose(xyz @ b, rotated @ core_frames(rotated, fixed), atol=1e-6)
    assert torch.allclose(torch.linalg.det(b), torch.ones(1), atol=1e-6)
    xyz[:, -1] = 1e9
    assert torch.equal(b, core_frames(xyz, fixed))
    with pytest.raises(ValueError):
        core_frames(torch.zeros_like(xyz), fixed)


def test_rng_and_sample_order_resume_exactly_without_workers():
    seed_all(44)
    stream = EpochStream(7)
    stream.take(5)
    state, order = rng_state(), copy.deepcopy(stream.state)
    expected = stream.take(13)
    noise = torch.randn(5)
    restore_rng(state)
    resumed = EpochStream(7, order)
    assert expected == resumed.take(13)
    assert torch.equal(noise, torch.randn(5))
    assert lr_factor(0, 10, 100) == 0.1 and lr_factor(100, 10, 100) == 0


def test_native_stereo_is_perceived_from_verified_coordinates_without_repair():
    from rdkit.Chem import AllChem

    source = Chem.AddHs(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    assert AllChem.EmbedMolecule(source, randomSeed=27) == 0
    actual = Chem.Mol(source)
    Chem.RemoveStereochemistry(actual)
    verify_preprocessed_molecule(actual, source)
    assert all(a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for a in actual.GetAtoms())
    wrong = Chem.Mol(source)
    for atom in wrong.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            atom.InvertChirality()
    with pytest.raises(ValueError):
        verify_preprocessed_molecule(actual, wrong)


@pytest.mark.parametrize(
    "smiles,fixed,variable,anchor,new_anchor",
    [
        ("CCC/C=C/CCCCCC", [1, 2, 3], list(range(4, 12)), 3, 4),
        ("CCC[C@H](N)CCCCC", [1, 2, 3, 4, 6, 7, 8, 9, 10], [5], 4, 5),
        ("CCNC(=O)CCCC", [1, 2, 3], [4, 5, 6, 7, 8, 9], 3, 4),
    ],
)
def test_unsafe_training_cuts_are_rejected(smiles, fixed, variable, anchor, new_anchor):
    mol = Chem.MolFromSmiles(smiles)
    with pytest.raises(ValueError):
        validate_cut(
            mol,
            list(range(1, mol.GetNumAtoms() + 1)),
            dict(fixed=fixed, variable=variable, attachment_fixed=anchor, attachment_variable=new_anchor),
        )


def test_lora_exact_targets_start_at_base_and_only_adapter_parameters_train():
    m = nn.Module()
    m.gen = nn.Module()
    m.gen.ligand_dec = nn.Module()
    m.gen.ligand_dec.inv_proj = nn.Linear(5, 7)
    m.requires_grad_(False)
    x = torch.randn(3, 5)
    expected = m.gen.ligand_dec.inv_proj(x).detach()
    adapters = install_lora(m, ["gen.ligand_dec.inv_proj"], rank=2)
    actual = m.gen.ligand_dec.inv_proj(x)
    assert torch.equal(actual, expected)
    actual.square().sum().backward()
    lora = adapters["gen.ligand_dec.inv_proj"]
    assert lora.b.grad.abs().sum() > 0 and lora.base.weight.grad is None
    with pytest.raises(ValueError):
        install_lora(m, ["gen.ligand_dec.missing_inv"])


def test_text_hook_requires_scoped_context():
    class Embed(nn.Module):
        def forward(self, x, *, atom_mask):
            return x, "edges"

    m = nn.Module()
    m.gen = nn.Module()
    m.gen.ligand_dec = nn.Module()
    m.gen.ligand_dec.inv_emb = Embed()
    adapter = NativeTextConditioner(m, 16, 8)
    x, valid = torch.randn(1, 3, 16), torch.ones(1, 3, dtype=torch.bool)
    with pytest.raises(RuntimeError):
        m.gen.ligand_dec.inv_emb(x, atom_mask=valid)
    with adapter.condition(torch.randn(1, 2, 8), torch.ones(1, 2, dtype=torch.bool)):
        assert torch.equal(m.gen.ligand_dec.inv_emb(x, atom_mask=valid)[0], x)
    adapter.close()
    assert torch.equal(m.gen.ligand_dec.inv_emb(x, atom_mask=valid)[0], x)
