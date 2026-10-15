"""Training bridge for the pinned FLOWR.root implementation, without patching it."""

from collections import OrderedDict
import copy
import hashlib
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem
import torch

from ..backends import check_backend, load_backend
from ..chem import canonical_smiles, heavy_copy, load_reference
from ..io import verify_ref


def core_center_prior(_self, from_mol, to_mol, to_mol_mask, from_mol_mask):
    """Noise centre depends exclusively on the supplied fixed core.

    Upstream's local training prior uses the *withheld fragment* centroid. That
    would leak the target pose in a reconstruction experiment.
    """
    if not to_mol_mask.any() or not (~from_mol_mask).any():
        raise ValueError("M1 local reconstruction requires fixed and variable atoms")
    coords = from_mol.coords.clone()
    coords[~from_mol_mask] += to_mol.coords[to_mol_mask].mean(0) - coords[~from_mol_mask].mean(0)
    return from_mol._copy_with(coords=coords)


def make_training_interpolant(eval_interpolant, vocab, charges, hybridization):
    from flowr.data.interpolate import ComplexInterpolant

    class ReviewedInterpolant(ComplexInterpolant):
        masks = None

        # Preserve the upstream method's keyword arguments for API compatibility.
        def _determine_modes_and_extract_masks(self, to_mols, rdkit_mols, interaction_mask):  # noqa: ARG002
            if self.masks is None or len(self.masks) != len(to_mols):
                raise ValueError("Explicit per-example masks are required")
            for mol, mask in zip(to_mols, self.masks):
                if mask.shape != (mol.seq_length,) or not mask.any() or mask.all():
                    raise ValueError("Native atom/mask mismatch")
            return [("fragment_inpainting", mask, True) for mask in self.masks]

        _align_prior_to_fragment = core_center_prior

    return ReviewedInterpolant(
        eval_interpolant.prior_sampler,
        ligand_coord_interpolation="linear",
        ligand_type_interpolation=eval_interpolant.type_interpolation,
        ligand_bond_interpolation=eval_interpolant.bond_interpolation,
        pocket_noise="fix",
        fragment_inpainting=True,
        inference=False,
        ligand_time_alpha=1.0,
        ligand_time_beta=1.0,
        mixed_uniform_beta_time=False,
        split_continuous_discrete_time=False,
        batch_ot=False,
        rotation_alignment=False,
        permutation_alignment=False,
        anisotropic_prior=False,
        ref_ligand_com_prior=False,
        flow_interactions=False,
        interaction_conditional=False,
        vocab=vocab,
        vocab_charges=charges,
        vocab_hybridization=hybridization,
    )


def to_device(value, device):
    if isinstance(value, torch.Tensor):
        return value.to(device)
    if isinstance(value, dict):
        return {k: to_device(v, device) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return type(value)(to_device(v, device) for v in value)
    return value


def rotate_coordinates(value, rotation):
    if isinstance(value, dict):
        return {
            k: v @ rotation if isinstance(v, torch.Tensor) and k.endswith("coords") else rotate_coordinates(v, rotation)
            for k, v in value.items()
        }
    if isinstance(value, (tuple, list)):
        return type(value)(rotate_coordinates(v, rotation) for v in value)
    return value


def verify_preprocessed_molecule(actual, expected):
    if (
        actual is None
        or [a.GetAtomicNum() for a in actual.GetAtoms()] != [a.GetAtomicNum() for a in expected.GetAtoms()]
        or not np.allclose(actual.GetConformer().GetPositions(), expected.GetConformer().GetPositions(), atol=1e-3)
    ):
        raise ValueError("Native preprocessing changed atom order/elements/coordinates")
    # GeometricMol carries stereo in 3D coordinates, not RDKit chiral tags.
    # Perceive those tags on the verification copy; no graph/coordinate repair.
    actual = Chem.Mol(actual)
    Chem.AssignStereochemistryFrom3D(actual, replaceExistingTags=True)
    if canonical_smiles(actual) != canonical_smiles(expected):
        raise ValueError(
            "Native preprocessing chemistry/3D stereo differs: "
            + canonical_smiles(actual)
            + " != "
            + canonical_smiles(expected)
        )


class NativeBridge:
    def __init__(self, backend_path, manifest_root, workdir, first_record, device="cuda"):
        torch.set_float32_matmul_precision("highest")
        self.backend = load_backend(backend_path)
        self.binding = check_backend(self.backend)
        source = self.backend["source"]
        sys.path.insert(0, source)
        import flowr
        import flowr.gen.generate_from_pdb as api
        from flowr.data.datamodules import GeometricInterpolantDM

        if Path(source).resolve() not in Path(flowr.__file__).resolve().parents:
            raise ValueError("FLOWR import does not come from the pinned checkout")
        self.api, self.root, self.workdir, self.device = api, Path(manifest_root), Path(workdir), device
        self.workdir.mkdir(parents=True, exist_ok=True)
        ref, pocket = self.stage(first_record)
        argv = [
            "flowr-root",
            "--pdb_file",
            str(pocket),
            "--ligand_file",
            str(ref),
            "--arch",
            "pocket",
            "--pocket_type",
            "holo",
            "--pocket_noise",
            "fix",
            "--ckpt_path",
            self.backend["checkpoint"],
            "--save_dir",
            str(workdir),
            "--gpus",
            "1" if device == "cuda" else "0",
            "--num_workers",
            "0",
            "--batch_cost",
            "1",
            "--no_ligand_valence_repair",
            "--fragment_inpainting",
        ]
        old = sys.argv
        try:
            sys.argv = argv
            self.args = api.get_args()
        finally:
            sys.argv = old
        unsafe = torch.serialization.get_unsafe_globals_in_checkpoint(self.backend["checkpoint"])
        if unsafe:
            raise ValueError("Unreviewed pickle globals in base checkpoint")
        loaded = api.load_model(self.args)
        self.model, self.hparams, self.vocab, self.charges, self.hybridization, self.aromatic, _, _ = loaded
        h = self.hparams
        if (
            not h.get("fragment_inpainting")
            or not h.get("remove_hs")
            or h.get("coord_scale") != 1.0
            or h.get("flow_interactions")
            or h.get("interaction_conditional")
        ):
            raise ValueError("Unsupported checkpoint training semantics")
        self.model.to(device).eval().requires_grad_(False)
        from .frame import CoreFrame

        self.core_frame = CoreFrame(self.model)
        self.transform, evaluation = api.util.load_util(
            self.args, h, self.vocab, self.charges, self.hybridization, self.aromatic
        )
        self.interpolant = make_training_interpolant(evaluation, self.vocab, self.charges, self.hybridization)
        self.dm = GeometricInterpolantDM(
            None, None, None, 1, 1, train_interpolant=self.interpolant, num_workers=0, pad_to_bucket=False
        )
        self.cache = OrderedDict()

    def stage(self, record):
        folder = self.workdir / hashlib.sha256(record["record_id"].encode()).hexdigest()[:16]
        folder.mkdir(exist_ok=True)
        source = verify_ref(self.root, record["files"]["reference"])
        mol = heavy_copy(load_reference(source))
        before = canonical_smiles(mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        ref = folder / "reference.sdf"
        with Chem.SDWriter(str(ref)) as writer:
            writer.write(mol)
        reread = Chem.MolFromMolFile(str(ref), removeHs=False)
        if canonical_smiles(reread) != before or not np.allclose(
            reread.GetConformer().GetPositions(), mol.GetConformer().GetPositions(), atol=1e-5
        ):
            raise ValueError("Staging changed reference chemistry or frame")
        original = verify_ref(self.root, record["files"]["pocket"])
        pocket = folder / "pocket.pdb"
        pocket.write_text(
            "".join(line for line in original.read_text().splitlines(keepends=True) if not line.startswith("CONECT"))
        )
        return ref, pocket

    def system(self, record):
        ident = record["record_id"]
        if ident in self.cache:
            self.cache.move_to_end(ident)
            return self.cache[ident]
        ref, pocket = self.stage(record)
        args = copy.copy(self.args)
        args.ligand_file, args.pdb_file = str(ref), str(pocket)
        system = self.api.util.load_data_from_pdb(
            args,
            remove_hs=self.hparams["remove_hs"],
            remove_aromaticity=self.hparams["remove_aromaticity"],
            ligand_idx=0,
            chain_id=None,
            canonicalize_conformer=False,
        )
        actual = system.ligand.to_rdkit(self.vocab, self.charges, self.hybridization, sanitise=True)
        expected = Chem.MolFromMolFile(str(ref), removeHs=False)
        verify_preprocessed_molecule(actual, expected)
        system = self.transform(system)
        self.cache[ident] = system
        if len(self.cache) > 64:
            self.cache.popitem(last=False)
        return system

    def batch(self, records, mask_indices=None, rotate=False):
        systems = [self.system(r) for r in records]
        masks, attachments = [], []
        for i, record in enumerate(records):
            idx = int(torch.randint(len(record["mask_candidates"]), (1,))) if mask_indices is None else mask_indices[i]
            fixed = set(record["mask_candidates"][idx]["fixed"])
            masks.append(torch.tensor([v in fixed for v in record["atom_maps"]], dtype=torch.bool))
            attachments.append(
                torch.tensor(
                    [v == record["mask_candidates"][idx]["attachment_fixed"] for v in record["atom_maps"]],
                    dtype=torch.bool,
                )
            )
        self.interpolant.masks = masks
        batch = self.dm._collate(systems, "train")
        batch[1]["attachment_mask"] = torch.nn.utils.rnn.pad_sequence(attachments, batch_first=True)
        # A proper SO(3) rotation, never a chirality-changing reflection.
        if rotate:
            q, _ = torch.linalg.qr(torch.randn(3, 3))
            q[:, 0] *= torch.linalg.det(q)
            batch = rotate_coordinates(batch, q)
        return to_device(batch, self.device)

    def predict(self, batch, sc_probability=0.5):
        prior, data, interp, times = batch
        builder, model = self.model.builder, self.model
        pocket = builder.extract_pocket_from_complex(data)
        target = builder.extract_ligand_from_complex(data)
        ligand = builder.extract_ligand_from_complex(interp)
        base = builder.extract_ligand_from_complex(prior)
        fixed = data["fragment_mask"].bool()
        for value in (target, ligand, base):
            value["fragment_mask"], value["fragment_mode"] = fixed, data["fragment_mode"]
        # No hidden target/interactions in self conditioning. Native inference
        # accepts this same zero-filled initial state.
        condition = None
        if model.self_condition:
            condition = {
                "coords": base["coords"],
                "atomics": torch.zeros_like(ligand["atomics"]),
                "bonds": torch.zeros_like(ligand["bonds"]),
            }
            if model.sc_charges or model.self_condition_mode == "residual":
                condition["charges"] = torch.zeros_like(ligand["charges"])
            if model.add_feats:
                condition["hybridization"] = torch.zeros_like(ligand["hybridization"])
            if sc_probability and torch.rand((), device=self.device) < sc_probability:
                with torch.no_grad():
                    previous = model(ligand, pocket, times.T, cond_batch=condition, training=True)
                    condition = {k: previous[k].softmax(-1) if k != "coords" else previous[k] for k in condition}
                    if model.self_condition_mode == "residual":
                        condition.update(mask=ligand["mask"], _sc_fired=True)
        prediction = model(ligand, pocket, times.T, cond_batch=condition, training=True)
        return prediction, target, fixed
