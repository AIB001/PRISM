"""Single-device, resumable adapter/LoRA training on native FLOWR x1 outputs."""

import argparse
from contextlib import nullcontext
import hashlib
import json
import math
import os
from pathlib import Path
import random
import time

import numpy as np
import torch

from ..io import read_json, sha256, write_json
from .conditioning import NativeTextConditioner, install_lora
from .data import read_records, validate_split
from .encoder import FrozenTextEncoder
from .losses import x1_losses
from .native import NativeBridge

DEFAULTS = {
    "seed": 20261013,
    "device": "cuda",
    "max_steps": 1000,
    "batch_size": 1,
    "gradient_accumulation": 8,
    "adapter_lr": 0.0001,
    "lora_lr": 0.00001,
    "weight_decay": 0.01,
    "warmup_steps": 50,
    "clip_grad_norm": 1.0,
    "text_dropout": 0.1,
    "self_condition_probability": 0.5,
    "precision": "float32",
    "validation_interval": 50,
    "validation_limit": 0,
    "early_stopping_patience": 10,
    "heads": 4,
    "lora_targets": [],
    "lora_rank": 8,
    "lora_alpha": 16.0,
    "rotate": True,
    "loss_weights": {"coords": 1.0, "atomics": 1.0, "charges": 1.0, "bonds": 1.0},
}
REQUIRED = {"backend_config", "manifest", "split", "text_encoder", "output"}


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def rng_state():
    state = np.random.get_state()
    return {
        "python": random.getstate(),
        "numpy": [state[0], state[1].tolist(), int(state[2]), int(state[3]), float(state[4])],
        "torch": torch.get_rng_state(),
        "cuda": torch.cuda.get_rng_state_all() if torch.cuda.is_available() else [],
    }


def restore_rng(state):
    random.setstate(state["python"])
    s = state["numpy"]
    np.random.set_state((s[0], np.asarray(s[1], dtype=np.uint32), s[2], s[3], s[4]))
    torch.set_rng_state(state["torch"].cpu())
    if state["cuda"]:
        torch.cuda.set_rng_state_all([v.cpu() for v in state["cuda"]])


class EpochStream:
    """Permutation + cursor are saved; there is no unsaved worker/prefetch RNG."""

    def __init__(self, size, state=None):
        self.size = size
        self.state = state or {"order": torch.randperm(size).tolist(), "cursor": 0, "epoch": 0}

    def take(self, count):
        result = []
        while len(result) < count:
            if self.state["cursor"] == self.size:
                self.state = {
                    "order": torch.randperm(self.size).tolist(),
                    "cursor": 0,
                    "epoch": self.state["epoch"] + 1,
                }
            idx = self.state["cursor"]
            take = min(count - len(result), self.size - idx)
            result += self.state["order"][idx : idx + take]
            self.state["cursor"] += take
        return result


def lr_factor(step, warmup, total):
    if step < warmup:
        return (step + 1) / max(warmup, 1)
    progress = min(1.0, (step - warmup) / max(total - warmup, 1))
    return 0.5 * (1 + math.cos(math.pi * progress))


def trainable_state(conditioner, lora):
    return {
        "conditioner": {k: v.detach().cpu().clone() for k, v in conditioner.state_dict().items()},
        "lora": {
            name: {k: getattr(layer, k).detach().cpu().clone() for k in ("a", "b")} for name, layer in lora.items()
        },
    }


def restore_trainable(state, conditioner, lora):
    conditioner.load_state_dict(state["conditioner"], strict=True)
    if set(state["lora"]) != set(lora):
        raise ValueError("LoRA target set differs from checkpoint")
    with torch.no_grad():
        for name, layer in lora.items():
            for key in ("a", "b"):
                getattr(layer, key).copy_(state["lora"][name][key])


def save_checkpoint(path, state):
    temporary = Path(str(path) + ".tmp")
    torch.save(state, temporary)
    os.replace(temporary, path)


def load_config(path):
    original = read_json(path)
    if not REQUIRED <= original.keys() or set(original) - REQUIRED - DEFAULTS.keys():
        raise ValueError("Unknown/incomplete training config")
    config = dict(DEFAULTS, **original)
    for key in ("backend_config", "manifest", "split", "output"):
        p = Path(config[key]).expanduser()
        config[key] = str((p if p.is_absolute() else Path(path).resolve().parent / p).resolve())
    for key in ("max_steps", "batch_size", "gradient_accumulation", "validation_interval", "early_stopping_patience"):
        if type(config[key]) is not int or config[key] <= 0:
            raise ValueError(key + " must be a positive integer")
    if config["precision"] not in {"float32", "bfloat16"} or config["device"] not in {"cpu", "cuda"}:
        raise ValueError("Unsupported device or precision; FP16 without loss scaling is not allowed")
    if any(not 0 <= config[k] <= 1 for k in ("text_dropout", "self_condition_probability")):
        raise ValueError("Probabilities must be in [0,1]")
    if (
        config["adapter_lr"] <= 0
        or config["lora_lr"] <= 0
        or config["clip_grad_norm"] <= 0
        or not 0 <= config["warmup_steps"] < config["max_steps"]
    ):
        raise ValueError("Invalid LR/gradient/warmup configuration")
    if (
        set(config["loss_weights"]) != {"coords", "atomics", "charges", "bonds"}
        or min(config["loss_weights"].values()) <= 0
    ):
        raise ValueError("Every supported x1 objective needs positive weight")
    return config


def run_training(config_path, resume=None, stop_after_step=None):
    config = load_config(config_path)
    output = Path(config["output"])
    if output.exists() and any(output.iterdir()) and resume is None:
        raise ValueError("Nonempty training output; explicitly resume or choose a new directory")
    output.mkdir(parents=True, exist_ok=True)
    if config["device"] == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA requested but unavailable")
    if config["precision"] == "bfloat16" and (config["device"] != "cuda" or not torch.cuda.is_bf16_supported()):
        raise ValueError("BF16 requires a supported CUDA device")
    seed_all(config["seed"])
    records = read_records(config["manifest"])
    split = read_json(config["split"])
    if split.get("manifest_sha256") != sha256(config["manifest"]):
        raise ValueError("Split is bound to a different manifest")
    validate_split(records, split)
    train = [r for r in records if split["assignments"][r["record_id"]]["split"] == "train"]
    val = [r for r in records if split["assignments"][r["record_id"]]["split"] == "val"]
    if config["validation_limit"]:
        val = val[: config["validation_limit"]]
    # Test records are deliberately never evaluated during fitting/model choice.
    bridge = NativeBridge(
        config["backend_config"], Path(config["manifest"]).parent, output / "native_inputs", train[0], config["device"]
    )
    encoder = FrozenTextEncoder(config["text_encoder"], device="cpu")
    width = bridge.model.gen.ligand_dec.d_inv
    conditioner = NativeTextConditioner(bridge.model, width, encoder.dim, config["heads"]).to(config["device"])
    lora = install_lora(bridge.model, config["lora_targets"], config["lora_rank"], config["lora_alpha"])
    groups = [{"params": list(conditioner.parameters()), "lr": config["adapter_lr"]}]
    if lora:
        groups.append({"params": [p for layer in lora.values() for p in (layer.a, layer.b)], "lr": config["lora_lr"]})
    parameters = [p for group in groups for p in group["params"]]
    optimizer = torch.optim.AdamW(groups, weight_decay=config["weight_decay"])
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lambda s: lr_factor(s, config["warmup_steps"], config["max_steps"])
    )
    config_binding = {k: v for k, v in config.items() if k != "output"}
    binding = {
        "config": config_binding,
        "manifest_sha256": sha256(config["manifest"]),
        "split_sha256": sha256(config["split"]),
        "base_sha256": bridge.backend["checkpoint_sha256"],
        "source_commit": bridge.backend["source_commit"],
        "text_encoder": encoder.binding,
        "torch_version": str(torch.__version__),
        "record_schema_sha256": sha256(Path(__file__).parents[1] / "schemas" / "m1-record-v1.json"),
        "implementation": {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
    }
    fingerprint = hashlib.sha256(json.dumps(binding, sort_keys=True).encode()).hexdigest()
    stream, step, best, bad_validations = EpochStream(len(train)), 0, float("inf"), 0
    if resume:
        state = torch.load(resume, map_location="cpu", weights_only=True)
        if state["format"] != "prism-design-m1-checkpoint/1" or state["fingerprint"] != fingerprint:
            raise ValueError("Resume binding/config/data/base/encoder differs")
        restore_trainable(state["trainable"], conditioner, lora)
        optimizer.load_state_dict(state["optimizer"])
        scheduler.load_state_dict(state["scheduler"])
        stream = EpochStream(len(train), state["stream"])
        step, best, bad_validations = state["step"], state["best_validation"], state["bad_validations"]
        restore_rng(state["rng"])
    initial = trainable_state(conditioner, lora)
    text_cache = {}

    def texts(batch):
        values = []
        for r in batch:
            if r["record_id"] not in text_cache:
                text_cache[r["record_id"]] = encoder.encode([r["text"]])
            values.append(text_cache[r["record_id"]])
        max_length = max(t.shape[1] for t, _ in values)
        t = torch.zeros(len(batch), max_length, encoder.dim)
        m = torch.zeros(len(batch), max_length, dtype=torch.bool)
        for i, (v, mask) in enumerate(values):
            t[i, : v.shape[1]], m[i, : v.shape[1]] = v[0], mask[0]
        return t.to(config["device"]), m.to(config["device"])

    def precision():
        return torch.autocast("cuda", dtype=torch.bfloat16) if config["precision"] == "bfloat16" else nullcontext()

    def compute(batch_records, training):
        native = bridge.batch(
            batch_records,
            rotate=config["rotate"] and training,
            mask_indices=None if training else [0] * len(batch_records),
        )
        tokens, mask = texts(batch_records)
        drop = torch.rand(len(batch_records), device=config["device"]) < config["text_dropout"] if training else None
        with conditioner.condition(tokens, mask, drop, attachment_mask=native[1]["attachment_mask"]), precision():
            prediction, target, fixed = bridge.predict(
                native, config["self_condition_probability"] if training else 0.0
            )
            return x1_losses(prediction, target, fixed, config["loss_weights"])

    def validate():
        saved = rng_state()
        try:
            seed_all(config["seed"] + 900001)
            conditioner.eval()
            metrics = []
            with torch.no_grad():
                for r in val:
                    total, loss = compute([r], False)
                    metrics.append({"total": float(total), **{k: float(v) for k, v in loss.items()}})
            return {k: sum(row[k] for row in metrics) / len(metrics) for k in metrics[0]}
        finally:
            restore_rng(saved)
            conditioner.train()

    metadata = {
        "format": "prism-design-m1-checkpoint/1",
        "fingerprint": fingerprint,
        "binding": binding,
        "runtime": {
            "invariant_dim": width,
            "text_dim": encoder.dim,
            "heads": config["heads"],
            "text_encoder": dict(config["text_encoder"], binding=encoder.binding),
            "lora_targets": config["lora_targets"],
            "lora_rank": config["lora_rank"],
            "lora_alpha": config["lora_alpha"],
            "task": "fragment_reconstruction",
            "prior": "fixed_core_centroid",
            "coordinate_frame": "right_handed_fixed_core/1",
            "attachment_conditioning": "fixed_attachment_node/1",
        },
    }
    write_json(output / "run_config.json", config)
    write_json(
        output / "provenance.json",
        {
            "fingerprint": fingerprint,
            "binding": binding,
            "train_records": len(train),
            "validation_records_used": len(val),
            "test_records_used": 0,
            "trainable_parameters": sum(p.numel() for p in parameters),
            "base_total_parameters": sum(p.numel() for p in bridge.model.parameters()),
            "torch_version": str(torch.__version__),
            "training_kind": "adapter_lora" if lora else "adapter_only",
        },
    )
    start = time.monotonic()
    last_validation = None
    last_grad = None
    stop = min(config["max_steps"], stop_after_step) if stop_after_step is not None else config["max_steps"]
    with (output / "metrics.jsonl").open("a", encoding="utf-8") as log:
        while step < stop:
            conditioner.train()
            optimizer.zero_grad(set_to_none=True)
            running = 0.0
            # Fixed full accumulation windows, including the last optimizer
            # step. There is no underweighted partial last batch.
            for _ in range(config["gradient_accumulation"]):
                batch_records = [train[i] for i in stream.take(config["batch_size"])]
                total, losses = compute(batch_records, True)
                (total / config["gradient_accumulation"]).backward()
                running += float(total.detach()) / config["gradient_accumulation"]
            last_grad = float(
                torch.nn.utils.clip_grad_norm_(parameters, config["clip_grad_norm"], error_if_nonfinite=True)
            )
            if not last_grad > 0:
                raise ValueError("Text/LoRA branch has no gradient")
            if any(p.grad is not None for name, p in bridge.model.named_parameters() if not p.requires_grad):
                raise ValueError("Frozen backbone accumulated gradients")
            optimizer.step()
            scheduler.step()
            step += 1
            event = {
                "step": step,
                "train_loss": running,
                "gradient_norm_before_clip": last_grad,
                "lr_next_step": [group["lr"] for group in optimizer.param_groups],
                "epoch": stream.state["epoch"],
            }
            improved = False
            if step % config["validation_interval"] == 0 or step == config["max_steps"]:
                last_validation = validate()
                event["validation"] = last_validation
                improved = last_validation["total"] < best
                best = min(best, last_validation["total"])
                bad_validations = 0 if improved else bad_validations + 1
            state = dict(
                metadata,
                trainable=trainable_state(conditioner, lora),
                optimizer=optimizer.state_dict(),
                scheduler=scheduler.state_dict(),
                step=step,
                best_validation=best,
                bad_validations=bad_validations,
                stream=stream.state,
                rng=rng_state(),
            )
            save_checkpoint(output / "last.pt", state)
            if improved:
                save_checkpoint(output / "best.pt", state)
            log.write(json.dumps(event, allow_nan=False) + "\n")
            log.flush()
            print(json.dumps(event), flush=True)
            if bad_validations >= config["early_stopping_patience"]:
                break
    final_state = trainable_state(conditioner, lora)
    changed = sum(not torch.equal(v, initial["conditioner"][k]) for k, v in final_state["conditioner"].items())
    report = {
        "status": "COMPLETED" if step >= config["max_steps"] else "CHECKPOINTED",
        "optimizer_steps": step,
        "changed_adapter_tensors_this_invocation": changed,
        "last_gradient_norm": last_grad,
        "best_validation": best if math.isfinite(best) else None,
        "validation": last_validation,
        "wall_seconds": time.monotonic() - start,
        "scientific_efficacy_claim": False,
        "peak_cuda_allocated_bytes": torch.cuda.max_memory_allocated() if config["device"] == "cuda" else None,
    }
    write_json(output / "training_summary.json", report)
    conditioner.close()
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--resume")
    parser.add_argument("--stop-after-step", type=int)
    args = parser.parse_args(argv)
    print(json.dumps(run_training(args.config, args.resume, args.stop_after_step), indent=2))


if __name__ == "__main__":
    main()
