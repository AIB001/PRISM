"""Check every record with the real pinned tokenizer and FLOWR preprocessing."""

import argparse
from collections import Counter
import json
from pathlib import Path

from ..io import read_json, sha256, write_json
from .data import read_records, validate_split
from .encoder import FrozenTextEncoder
from .native import NativeBridge
from .train import load_config


def preflight(config_path, output):
    config = load_config(config_path)
    manifest = Path(config["manifest"])
    directory = Path(output).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    records = read_records(manifest)
    bridge = NativeBridge(config["backend_config"], manifest.parent, directory / "inputs_cache", records[0], "cpu")
    encoder = FrozenTextEncoder(config["text_encoder"])
    from .frame import core_frames
    import torch

    accepted, failures, audit = [], [], []
    for r in records:
        try:
            system = bridge.system(r)
            for cut in r["mask_candidates"]:
                fixed_ids = set(cut["fixed"])
                mask_for_frame = torch.tensor([[a in fixed_ids for a in r["atom_maps"]]])
                core_frames(system.ligand.coords[None], mask_for_frame)
            tokens, mask = encoder.encode([r["text"]])
            # Real collation catches native vocabularies/bonds/mask polarity.
            native = bridge.batch([r], [0], rotate=False)
            prior, target, interp, times = native
            fixed = target["fragment_mask"].bool()
            p = bridge.model.builder.extract_ligand_from_complex(prior)
            t = bridge.model.builder.extract_ligand_from_complex(target)
            if not torch.allclose(p["coords"][fixed], t["coords"][fixed], atol=1e-6):
                raise ValueError("Fixed-core coordinates changed in prior")
            pair = fixed[:, :, None] & fixed[:, None, :]
            # Native targets use all-zero rows for absent/self edges; priors
            # use the categorical one-hot "no bond" row. Compare chemistry.
            if not torch.equal(p["bonds"].argmax(-1)[pair], t["bonds"].argmax(-1)[pair]):
                raise ValueError("Fixed-core bonds changed in prior")
            accepted.append(r)
            audit.append(
                {
                    "id": r["record_id"],
                    "status": "PASS",
                    "text_tokens": int(mask.sum()),
                    "fixed_atoms": int(fixed.sum()),
                    "heavy_atoms": int(t["mask"].sum()),
                }
            )
        except Exception as exc:
            failures.append({"id": r["record_id"], "status": "EXCLUDED", "error": type(exc).__name__ + ": " + str(exc)})
            if len(failures) <= 3:
                print(json.dumps(failures[-1]), flush=True)
        if (len(audit) + len(failures)) % 25 == 0:
            print(
                json.dumps({"checked": len(audit) + len(failures), "accepted": len(audit), "excluded": len(failures)}),
                flush=True,
            )
    write_json(directory / "record_checks.json", audit)
    write_json(directory / "excluded.json", failures)
    write_json(
        directory / "invariant_linear_modules.json",
        {
            name: [module.in_features, module.out_features]
            for name, module in bridge.model.named_modules()
            if name.startswith("gen.ligand_dec.") and "inv" in name and isinstance(module, torch.nn.Linear)
        },
    )
    if not accepted:
        raise ValueError("No record passed native preflight; inspect excluded.json")
    # Preserve relative input paths and original group assignments. A subset
    # never grants permission to loosen the original leakage groups.
    destination = manifest.parent / "native_records.jsonl"
    with destination.open("w", encoding="utf-8") as handle:
        for r in accepted:
            handle.write(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n")
    split = read_json(config["split"])
    ids = {r["record_id"] for r in accepted}
    split["assignments"] = {k: v for k, v in split["assignments"].items() if k in ids}
    split["counts"] = dict(Counter(v["split"] for v in split["assignments"].values()))
    split["manifest_sha256"] = sha256(destination)
    split["parent_split_sha256"] = sha256(config["split"])
    validate_split(accepted, split)
    write_json(manifest.parent / "native_split.json", split)
    write_json(directory / "record_checks.json", audit)
    write_json(directory / "excluded.json", failures)
    summary = {
        "checked": len(records),
        "accepted": len(accepted),
        "excluded": len(failures),
        "split_counts": split["counts"],
        "base_sha256": bridge.backend["checkpoint_sha256"],
        "source_commit": bridge.backend["source_commit"],
        "encoder": encoder.binding,
        "source_manifest_sha256": sha256(manifest),
        "native_manifest_sha256": sha256(destination),
    }
    write_json(directory / "summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    preflight(args.config, args.output)
