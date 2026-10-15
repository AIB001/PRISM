"""Install trained text weights in the actual native sampling path."""

from contextlib import contextmanager
from types import MethodType

import torch

from .conditioning import NativeTextConditioner, install_lora
from .encoder import FrozenTextEncoder
from .frame import CoreFrame
from .native import core_center_prior
from .train import restore_trainable


class TextConditionedSampling:
    def __init__(
        self, model, checkpoint, instruction, device, base_sha256, source_commit, interpolant, attachment_index, n_atoms
    ):
        torch.set_float32_matmul_precision("highest")
        state = torch.load(checkpoint, map_location="cpu", weights_only=True)
        if state.get("format") != "prism-design-m1-checkpoint/1":
            raise ValueError("Unsupported M1 checkpoint format")
        binding = state["binding"]
        if binding["base_sha256"] != base_sha256 or binding["source_commit"] != source_commit:
            raise ValueError("M1 checkpoint requires its original base/source")
        runtime = state["runtime"]
        if runtime["task"] != "fragment_reconstruction" or runtime["prior"] != "fixed_core_centroid":
            raise ValueError("Unsupported M1 training semantics")
        if runtime.get("coordinate_frame") != "right_handed_fixed_core/1":
            raise ValueError("M1 weights predate the validated coordinate-frame protocol")
        if runtime.get("attachment_conditioning") != "fixed_attachment_node/1":
            raise ValueError("M1 weights require explicit attachment conditioning")
        self.core_frame = CoreFrame(model)
        if model.gen.ligand_dec.d_inv != runtime["invariant_dim"]:
            raise ValueError("Decoder dimensions do not match M1 weights")
        self.encoder = FrozenTextEncoder(runtime["text_encoder"], "cpu")
        self.conditioner = NativeTextConditioner(
            model, runtime["invariant_dim"], runtime["text_dim"], runtime["heads"]
        ).to(device)
        lora = install_lora(model, runtime["lora_targets"], runtime["lora_rank"], runtime["lora_alpha"])
        restore_trainable(state["trainable"], self.conditioner, lora)
        self.conditioner.eval()
        model.eval()
        tokens, mask = self.encoder.encode([instruction])
        self.tokens, self.mask = tokens.to(device), mask.to(device)
        self.attachment_mask = torch.zeros(1, n_atoms, dtype=torch.bool, device=device)
        self.attachment_mask[0, attachment_index] = True
        if not instruction.strip():
            self.mask.zero_()  # learned null branch, not an arbitrary empty BERT sentence
        # Train and inference share the core-only prior; no hidden edited pose.
        interpolant._align_prior_to_fragment = MethodType(core_center_prior, interpolant)
        interpolant.rotation_alignment = False
        interpolant.permutation_alignment = False
        interpolant.anisotropic_prior = False
        interpolant.ref_ligand_com_prior = False
        self.metadata = {
            "native_text_conditioning": True,
            "optimizer_step": state["step"],
            "training_fingerprint": state["fingerprint"],
            "text_tokens": int(mask.sum()),
            "prior": "fixed_core_centroid",
            "guidance_scale": 1.0,
            "efficacy_validated": False,
        }
        self.metadata["coordinate_frame"] = runtime["coordinate_frame"]
        self.metadata["attachment_conditioning"] = runtime["attachment_conditioning"]
        self.metadata["attachment_index_0based"] = attachment_index

    @contextmanager
    def condition(self):
        with self.conditioner.condition(self.tokens, self.mask, attachment_mask=self.attachment_mask):
            yield
