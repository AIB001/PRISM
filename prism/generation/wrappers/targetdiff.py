"""Stable PRISM-to-TargetDiff command contract.

This module is intended to run inside an environment containing the official
TargetDiff checkout. It delegates sampling to the upstream script while
providing explicit paths and a reproducible per-request sampling config.
"""

import argparse
import collections
import os
import pickle
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import yaml

# These wrappers are invoked two ways: as a package module during testing, and
# as a standalone script by every generation config (the model environments do
# not have PRISM installed). A plain relative import works only in the first
# case, so fall back to a path-based import for the second.
try:
    from . import _device as _device_support
except ImportError:  # standalone script invocation
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _device as _device_support  # type: ignore[no-redef]


_SAFE_PICKLE_GLOBALS: Dict[Tuple[str, str], Any] = {}


class _RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str):
        try:
            return _SAFE_PICKLE_GLOBALS[(module, name)]
        except KeyError as exc:
            raise pickle.UnpicklingError(
                f"Checkpoint contains forbidden pickle global: {module}.{name}"
            ) from exc


class _RestrictedPickleModule:
    __name__ = "prism_restricted_pickle"
    Unpickler = _RestrictedUnpickler
    load = staticmethod(pickle.load)
    loads = staticmethod(pickle.loads)
    dump = staticmethod(pickle.dump)
    dumps = staticmethod(pickle.dumps)


def _sanitize_checkpoint(source: Path, output: Path) -> Path:
    """Load a checkpoint through a strict pickle allowlist and rewrite it safely."""
    import torch
    from easydict import EasyDict

    _SAFE_PICKLE_GLOBALS.clear()
    _SAFE_PICKLE_GLOBALS.update(
        {
            ("collections", "OrderedDict"): collections.OrderedDict,
            ("easydict", "EasyDict"): EasyDict,
            ("torch", "FloatStorage"): torch.FloatStorage,
            ("torch._utils", "_rebuild_tensor_v2"): torch._utils._rebuild_tensor_v2,
        }
    )
    checkpoint = torch.load(
        str(source),
        map_location="cpu",
        pickle_module=_RestrictedPickleModule,
    )
    if not isinstance(checkpoint, dict) or not {"config", "model"}.issubset(checkpoint):
        raise SystemExit("TargetDiff checkpoint must contain 'config' and 'model' mappings")
    if not isinstance(checkpoint["model"], (dict, collections.OrderedDict)):
        raise SystemExit("TargetDiff checkpoint 'model' must be a state dictionary")
    torch.save(checkpoint, str(output))
    return output


def _device(value: str) -> str:
    """Resolve --device against the host; ``auto`` falls back to CPU."""
    return _device_support.resolve(value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run official TargetDiff pocket sampling")
    parser.add_argument("--targetdiff-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--pocket", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--num-samples", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", type=_device, default="auto")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--num-steps", type=int, default=1000)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.targetdiff_root.expanduser().resolve()
    script = root / "scripts" / "sample_for_pocket.py"
    if not script.is_file():
        raise SystemExit(f"TargetDiff sampling script not found: {script}")
    checkpoint = args.checkpoint.expanduser().resolve()
    pocket = args.pocket.expanduser().resolve()
    if not checkpoint.is_file():
        raise SystemExit(f"TargetDiff checkpoint not found: {checkpoint}")
    if not pocket.is_file() or pocket.suffix.lower() != ".pdb":
        raise SystemExit(f"TargetDiff requires a pocket PDB: {pocket}")
    if args.num_samples <= 0 or args.batch_size <= 0 or args.num_steps <= 0:
        raise SystemExit("num-samples, batch-size, and num-steps must be positive")

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    safe_checkpoint = _sanitize_checkpoint(checkpoint, output / "trusted_checkpoint.pt")
    config = {
        "model": {"checkpoint": str(safe_checkpoint)},
        "sample": {
            "seed": args.seed,
            "num_samples": args.num_samples,
            "num_steps": args.num_steps,
            "pos_only": False,
            "center_pos_mode": "protein",
            "sample_num_atoms": "prior",
        },
    }
    config_path = output / "prism_targetdiff_sampling.yml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    command = [
        sys.executable,
        str(script),
        str(config_path),
        "--pdb_path",
        str(pocket),
        "--device",
        args.device,
        "--batch_size",
        str(min(args.batch_size, args.num_samples)),
        "--result_path",
        str(output),
        "--num_samples",
        str(args.num_samples),
    ]
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(root) if not existing_pythonpath else str(root) + os.pathsep + existing_pythonpath
    )
    completed = subprocess.run(command, cwd=str(root), env=environment, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
