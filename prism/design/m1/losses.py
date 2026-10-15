"""Masked x1 objectives; fixed/variable boundary bonds remain supervised."""

import torch
from torch.nn import functional as F


def _mean_per_example(values, mask):
    if values.shape != mask.shape:
        raise ValueError("Loss mask shape mismatch")
    counts = mask.flatten(1).sum(1)
    if not (counts > 0).all():
        raise ValueError("Every training example needs supervised entries")
    return ((values * mask).flatten(1).sum(1) / counts).mean()


def x1_losses(prediction, target, fixed_mask, weights=None):
    """Uniform time, direct clean-coordinate prediction; NOT a velocity target.

    Coordinate values use the native checkpoint coordinate scale. Categorical
    targets are one-hot native vocabularies. No padding, diagonal, or F-F bond
    contributes; both absent and present U-U/F-U edges contribute once.
    """
    weights = weights or {"coords": 1.0, "atomics": 1.0, "charges": 1.0, "bonds": 1.0}
    valid = target["mask"].bool()
    if fixed_mask.shape != valid.shape or (fixed_mask.bool() & ~valid).any():
        raise ValueError("Fixed mask includes padding or has incorrect shape")
    variable = valid & ~fixed_mask.bool()
    n = valid.shape[1]
    pair = valid[:, :, None] & valid[:, None, :] & (variable[:, :, None] | variable[:, None, :])
    pair &= torch.ones(n, n, device=valid.device, dtype=torch.bool).triu(1)
    losses = {"coords": _mean_per_example((prediction["coords"] - target["coords"]).square().sum(-1), variable)}
    for name in ("atomics", "charges", "bonds"):
        logits, labels = prediction[name], target[name].argmax(-1)
        value = F.cross_entropy(logits.reshape(-1, logits.shape[-1]), labels.flatten(), reduction="none").reshape(
            labels.shape
        )
        losses[name] = _mean_per_example(value, pair if name == "bonds" else variable)
    total = sum(losses[name] * weights[name] for name in losses)
    if not torch.isfinite(total):
        raise FloatingPointError("Nonfinite M1 loss")
    return total, losses


def alignment_set_loss(attention, positive_atoms, observed, atom_mask):
    """Optional set-valued token-to-atom grounding. Unknown is never negative.

    attention [B,N,T] is normalized across valid atoms for each annotated token.
    positive_atoms [B,N,T] can include symmetry-equivalent atoms. Zero observed
    tokens returns a differentiable zero; this dataset supplies no such labels.
    """
    p = attention * atom_mask[..., None]
    p = p / p.sum(1, keepdim=True).clamp_min(1e-12)
    positive = positive_atoms.bool() & atom_mask[..., None].bool()
    known = observed.bool()
    if (known & ~positive.any(1)).any():
        raise ValueError("Observed alignment needs a nonempty valid positive set")
    mass = (p * positive).sum(1).clamp_min(1e-12)
    return -(mass.log() * known).sum() / known.sum().clamp_min(1)
