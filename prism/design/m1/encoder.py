"""Frozen local text encoder with explicit model, tokenizer and length binding."""

from pathlib import Path

import torch
from transformers import AutoModel, AutoTokenizer

from ..io import sha256


class FrozenTextEncoder:
    def __init__(self, config, device="cpu"):
        self.config = dict(config)
        directory = Path(config["local_path"]).expanduser().resolve()
        if not directory.is_dir() or not config.get("revision"):
            raise ValueError("A downloaded, revision-pinned encoder is required")
        self.max_length = config.get("max_length", 128)
        if not 4 <= self.max_length <= 128:
            raise ValueError("This multilingual encoder protocol allows 4..128 tokens")
        self.tokenizer = AutoTokenizer.from_pretrained(directory, local_files_only=True, trust_remote_code=False)
        self.model = (
            AutoModel.from_pretrained(directory, local_files_only=True, trust_remote_code=False, use_safetensors=True)
            .to(device)
            .eval()
            .requires_grad_(False)
        )
        self.device = device
        self.dim = self.model.config.hidden_size
        self.binding = {
            "model_id": config["model_id"],
            "revision": config["revision"],
            "max_length": self.max_length,
            "hidden_size": self.dim,
            "files": {p.name: sha256(p) for p in sorted(directory.iterdir()) if p.is_file()},
        }
        if config.get("binding") is not None and config["binding"] != self.binding:
            raise ValueError("Text model/tokenizer files changed")

    @torch.no_grad()
    def encode(self, texts):
        # Never silently discard the end of a constraint/caption.
        encoded = self.tokenizer(texts, padding=True, truncation=False, return_tensors="pt")
        if encoded.input_ids.shape[1] > self.max_length:
            raise ValueError(f"Text exceeds {self.max_length} tokens; shorten it explicitly")
        encoded = {k: v.to(self.device) for k, v in encoded.items()}
        states = self.model(**encoded).last_hidden_state
        return states.detach(), encoded["attention_mask"].bool()
