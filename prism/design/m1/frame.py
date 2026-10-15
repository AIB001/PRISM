"""SO(3)-covariant core frame around the orientation-sensitive native network.

The pinned decoder broadcasts scalar time embeddings onto xyz channels
(flowr/models/pocket.py, LigandDecoder.forward). An invariant-only text adapter
does not remove that native orientation bias. All coordinates are rotated into
a deterministic frame from the KNOWN mapped core, then outputs rotate back.
"""

from contextvars import ContextVar

import torch


def core_frames(coords, fixed, min_length=1e-4, min_sine=1e-3):
    """Columns are a right-handed orthonormal basis in the supplied frame.

    Core atom order is part of the bound input contract. Choose the first
    distinct pair and first noncollinear third atom, with detached decisions.
    Withheld coordinates never participate. Degeneracy fails explicitly.
    """
    if coords.ndim != 3 or fixed.shape != coords.shape[:2]:
        raise ValueError("Core frame coordinate/mask shape mismatch")
    bases = []
    for xyz, mask in zip(coords.detach(), fixed.bool()):
        points = xyz[mask]
        if len(points) < 3:
            raise ValueError("A core coordinate frame needs three fixed atoms")
        vectors = points[1:] - points[0]
        lengths = torch.linalg.vector_norm(vectors, dim=-1)
        distinct = (lengths > min_length).nonzero().flatten()
        if not len(distinct):
            raise ValueError("Degenerate coincident fixed core")
        first = distinct[0]
        u = vectors[first] / lengths[first]
        orthogonal = vectors - (vectors @ u)[:, None] * u
        transverse = torch.linalg.vector_norm(orthogonal, dim=-1)
        candidates = (
            ((lengths > min_length) & (transverse / lengths.clamp_min(min_length) > min_sine)).nonzero().flatten()
        )
        if not len(candidates):
            raise ValueError("Collinear fixed core has no supported coordinate frame")
        second = candidates[0]
        v = orthogonal[second] / transverse[second]
        w = torch.linalg.cross(u, v)
        bases.append(torch.stack((u, v, w), dim=-1))
    return torch.stack(bases)


class CoreFrame:
    """Hooks preserve the pretrained parameter names and weights exactly.

    Native generation can supply a precomputed pocket representation in the
    original frame. Discard that cache and recompute it after rotation. Do not
    reuse orientation-dependent pocket features across frames.
    """

    def __init__(self, native_model):
        self._basis = ContextVar(f"prism_core_frame_{id(self)}", default=None)
        self._pre = native_model.gen.register_forward_pre_hook(self._enter, with_kwargs=True)
        self._post = native_model.gen.register_forward_hook(self._leave, with_kwargs=True, always_call=True)

    def _enter(self, _module, args, kwargs):
        coords = args[0]
        fixed = kwargs.get("inpaint_mask")
        if fixed is None:
            raise ValueError("M1 canonical core frame requires an explicit fixed mask")
        # Gram-Schmidt and frame transforms are small and numerically sensitive.
        # Keep these operations in float64 even under BF16 autocast; the native
        # network still receives its original dtype. This does not add a target.
        basis = core_frames(coords.double(), fixed)
        self._basis.set(basis)
        values = dict(kwargs)
        for key in ("cond_coords", "pocket_coords"):
            if values.get(key) is not None:
                values[key] = torch.bmm(values[key].double(), basis).to(values[key].dtype)
        if values.get("pocket_equis") is not None or values.get("pocket_invs") is not None:
            values["pocket_equis"] = values["pocket_invs"] = None
        return (torch.bmm(coords.double(), basis).to(coords.dtype), *args[1:]), values

    def _leave(self, _module, _args, _kwargs, output):
        basis = self._basis.get()
        self._basis.set(None)
        if output is None:
            return output
        if basis is None:
            raise RuntimeError("Unbalanced M1 core-frame hooks")
        result = dict(output)
        result["coords"] = torch.bmm(output["coords"].double(), basis.transpose(1, 2)).to(output["coords"].dtype)
        return result

    def close(self):
        self._pre.remove()
        self._post.remove()
