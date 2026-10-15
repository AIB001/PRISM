"""Text changes invariant node features, never Cartesian coordinate channels."""

from contextlib import contextmanager
from contextvars import ContextVar
import math

import torch
from torch import nn


class TextAdapter(nn.Module):
    def __init__(self, invariant_dim, text_dim, heads=4):
        super().__init__()
        if invariant_dim % heads:
            raise ValueError("invariant_dim must be divisible by heads")
        self.norm = nn.LayerNorm(invariant_dim)
        self.text_projection = nn.Linear(text_dim, invariant_dim)
        self.film = nn.Linear(text_dim, invariant_dim * 2)
        self.attention = nn.MultiheadAttention(invariant_dim, heads, batch_first=True, dropout=0.0)
        self.null_token = nn.Parameter(torch.zeros(1, 1, text_dim))
        self.attachment_vector = nn.Parameter(torch.zeros(1, 1, invariant_dim))
        # One zero projection per branch; no additional zero gate, which would
        # kill all first-step gradients. Initial outputs match the frozen
        # backbone with the same core-frame wrapper.
        nn.init.zeros_(self.film.weight)
        nn.init.zeros_(self.film.bias)
        nn.init.zeros_(self.attention.out_proj.weight)
        nn.init.zeros_(self.attention.out_proj.bias)

    def forward(self, invariant, tokens, token_mask, atom_mask, drop_text=None, attachment_mask=None):
        if tokens.ndim != 3 or token_mask.shape != tokens.shape[:2]:
            raise ValueError("Text must be [batch,tokens,hidden] with a matching mask")
        if invariant.shape[:2] != atom_mask.shape or tokens.shape[0] != invariant.shape[0]:
            raise ValueError("Atom/text batch mismatch")
        token_mask = token_mask.bool()
        # Empty instructions get one learned null token, never an all-masked
        # softmax. Classifier-free dropout removes every original token.
        empty = ~token_mask.any(-1)
        if drop_text is not None:
            empty = empty | drop_text.bool()
        tokens = torch.cat((self.null_token.expand(tokens.shape[0], -1, -1), tokens), dim=1)
        token_mask = torch.cat((empty[:, None], token_mask & ~empty[:, None]), dim=1)
        pooled = (tokens * token_mask[..., None]).sum(1) / token_mask.sum(1, keepdim=True)
        scale, shift = self.film(pooled).chunk(2, dim=-1)
        values = self.text_projection(tokens)
        update, weights = self.attention(
            self.norm(invariant),
            values,
            values,
            key_padding_mask=~token_mask,
            need_weights=True,
            average_attn_weights=True,
        )
        residual = scale[:, None] * self.norm(invariant) + shift[:, None] + update
        if attachment_mask is not None:
            if (
                attachment_mask.shape != atom_mask.shape
                or not (attachment_mask.sum(-1) == 1).all()
                or (attachment_mask.bool() & ~atom_mask.bool()).any()
            ):
                raise ValueError("Exactly one valid fixed attachment node is required")
            residual = residual + attachment_mask[..., None] * self.attachment_vector
        return invariant + residual * atom_mask[..., None], weights[..., 1:]


class NativeTextConditioner(nn.Module):
    """Scoped hook into the pinned native decoder's invariant embedding.

    Kept separate from the backbone state dict. Do not use activation
    checkpointing across the context boundary. No no_grad around the frozen
    backbone: gradients must propagate back through it to these adapters.
    """

    def __init__(self, native_model, invariant_dim, text_dim, heads=4):
        super().__init__()
        self.adapter = TextAdapter(invariant_dim, text_dim, heads)
        self._context = ContextVar(f"prism_text_{id(self)}", default=None)
        self.last_attention = None
        self._hook = native_model.gen.ligand_dec.inv_emb.register_forward_hook(self._inject, with_kwargs=True)

    def _inject(self, _module, _args, kwargs, output):
        context = self._context.get()
        if context is None:
            raise RuntimeError("M1 decoder used without an explicit text/null context")
        invariant, edges = output
        tokens, token_mask, drop_text, attachment_mask = context
        if tokens.shape[0] == 1 and invariant.shape[0] > 1:
            tokens = tokens.expand(invariant.shape[0], -1, -1)
            token_mask = token_mask.expand(invariant.shape[0], -1)
            if drop_text is not None:
                drop_text = drop_text.expand(invariant.shape[0])
            if attachment_mask is not None:
                attachment_mask = attachment_mask.expand(invariant.shape[0], -1)
        atom_mask = kwargs.get("atom_mask")
        if atom_mask is None:
            raise RuntimeError("Pinned invariant-embedding signature changed")
        updated, self.last_attention = self.adapter(
            invariant, tokens, token_mask, atom_mask, drop_text, attachment_mask
        )
        return updated, edges

    @contextmanager
    def condition(self, tokens, mask, drop_text=None, attachment_mask=None):
        token = self._context.set((tokens, mask, drop_text, attachment_mask))
        try:
            yield
        finally:
            self._context.reset(token)

    def close(self):
        self._hook.remove()


class LoRALinear(nn.Module):
    def __init__(self, linear, rank=8, alpha=16.0, dropout=0.0):
        super().__init__()
        if not isinstance(linear, nn.Linear) or rank < 1:
            raise ValueError("LoRA requires an explicit Linear module and positive rank")
        self.base = linear.requires_grad_(False)
        self.a = nn.Parameter(torch.empty(rank, linear.in_features, device=linear.weight.device))
        self.b = nn.Parameter(torch.zeros(linear.out_features, rank, device=linear.weight.device))
        nn.init.kaiming_uniform_(self.a, a=math.sqrt(5))
        self.scale = alpha / rank
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        return self.base(x) + ((self.dropout(x) @ self.a.T) @ self.b.T) * self.scale


def install_lora(native_model, exact_module_names, rank=8, alpha=16.0):
    if len(set(exact_module_names)) != len(exact_module_names):
        raise ValueError("Duplicate LoRA targets")
    modules = dict(native_model.named_modules())
    # Explicit, reviewed invariant-only target names. No fuzzy match that can
    # silently train zero parameters or accidentally alter vector operations.
    for name in exact_module_names:
        if not name.startswith("gen.ligand_dec.") or "inv" not in name or not isinstance(modules.get(name), nn.Linear):
            raise ValueError(f"Not an explicit decoder invariant Linear target: {name}")
    installed = {}
    for name in exact_module_names:
        parent_name, child_name = name.rsplit(".", 1)
        layer = LoRALinear(modules[name], rank, alpha)
        setattr(native_model.get_submodule(parent_name), child_name, layer)
        installed[name] = layer
    return installed
