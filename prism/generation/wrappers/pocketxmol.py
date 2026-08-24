"""Stable PRISM-to-PocketXMol SBDD command contract."""

import argparse
import collections
import csv
import math
import os
import pickle
import shutil
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

_UPSTREAM_LAUNCHER = r'''"""Run PocketXMol inference without importing its training-only Lightning stack."""
import os
import runpy
import sys
import types
from pathlib import Path

root = Path(os.environ["PRISM_POCKETXMOL_ROOT"]).resolve()
sys.path.insert(0, str(root))

from utils.transforms import FeaturizeMol, FeaturizePocket


class DataModule:
    """Inference-only subset of PocketXMol's official training DataModule."""

    def __init__(self, config):
        self.config = config

    def get_featurizers(self):
        featurizer = FeaturizeMol(self.config.transforms.featurizer)
        if "featurizer_pocket" in self.config.transforms:
            return [
                FeaturizePocket(self.config.transforms.featurizer_pocket),
                featurizer,
            ]
        return [featurizer]

    def get_in_dims(self, featurizers=None):
        if featurizers is None:
            featurizers = self.get_featurizers()
        dimensions = {
            "num_node_types": featurizers[-1].num_node_types,
            "num_edge_types": featurizers[-1].num_edge_types,
        }
        if len(featurizers) == 2:
            dimensions["pocket_in_dim"] = featurizers[0].feature_dim
        return dimensions


training_stub = types.ModuleType("scripts.train_pl")
training_stub.DataModule = DataModule
sys.modules["scripts.train_pl"] = training_stub

script = root / "scripts" / "sample_use.py"
sys.argv[0] = str(script)
runpy.run_path(str(script), run_name="__main__")
'''


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
    """Load the official checkpoint with a strict pickle allowlist and rewrite it."""
    import torch
    from easydict import EasyDict

    _SAFE_PICKLE_GLOBALS.clear()
    _SAFE_PICKLE_GLOBALS.update(
        {
            ("argparse", "Namespace"): argparse.Namespace,
            ("collections", "OrderedDict"): collections.OrderedDict,
            ("easydict", "EasyDict"): EasyDict,
            ("torch", "DoubleStorage"): torch.DoubleStorage,
            ("torch", "FloatStorage"): torch.FloatStorage,
            ("torch._utils", "_rebuild_tensor_v2"): torch._utils._rebuild_tensor_v2,
        }
    )
    checkpoint = torch.load(
        str(source), map_location="cpu", pickle_module=_RestrictedPickleModule
    )
    if not isinstance(checkpoint, dict) or "state_dict" not in checkpoint:
        raise SystemExit("PocketXMol checkpoint must contain a 'state_dict' mapping")
    if not isinstance(checkpoint["state_dict"], (dict, collections.OrderedDict)):
        raise SystemExit("PocketXMol checkpoint 'state_dict' must be a mapping")
    torch.save(checkpoint, str(output))
    return output


def _device(value: str) -> str:
    """Resolve --device against the host; ``auto`` falls back to CPU."""
    return _device_support.resolve(value)


def _positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("value must be a positive finite number")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run official PocketXMol SBDD sampling")
    parser.add_argument("--pocketxmol-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--train-config", type=Path, required=True)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--reference-ligand", default="")
    parser.add_argument("--center", type=float, nargs=3, required=True)
    parser.add_argument("--radius", type=_positive_float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--num-samples", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", type=_device, default="auto")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--num-steps", type=int, default=100)
    parser.add_argument("--mean-atoms", type=_positive_float, default=28.0)
    parser.add_argument("--std-atoms", type=_positive_float, default=2.0)
    parser.add_argument("--mode", choices=("simple", "refine"), default="simple")
    return parser


def _task_config(args: argparse.Namespace, protein: Path, reference: Optional[Path]) -> dict:
    pocket_args: Dict[str, Any] = {"radius": args.radius}
    if reference is not None and reference.suffix.lower() == ".sdf":
        pocket_args["ref_ligand_path"] = str(reference)
    else:
        pocket_args["pocket_coord"] = list(args.center)

    config: Dict[str, Any] = {
        "sample": {
            "seed": args.seed,
            "batch_size": min(args.batch_size, args.num_samples),
            "num_mols": args.num_samples,
            "save_traj_prob": 0,
        },
        "data": {
            "protein_path": str(protein),
            "is_pep": False,
            "pocket_args": pocket_args,
            "pocmol_args": {"data_id": "prism_sbdd", "pdbid": protein.stem},
        },
        "transforms": {
            "featurizer_pocket": {"center": list(args.center)},
            "variable_mol_size": {
                "name": "variable_mol_size",
                "num_atoms_distri": {
                    "strategy": "mol_atoms_based",
                    "mean": {"coef": 0, "bias": args.mean_atoms},
                    "std": {"coef": 0, "bias": args.std_atoms},
                    "min": 5,
                },
            },
        },
    }
    if args.mode == "simple":
        config.update(
            {
                "task": {"name": "sbdd", "transform": {"name": "sbdd"}},
                "noise": {
                    "name": "sbdd",
                    "num_steps": args.num_steps,
                    "prior": "from_train",
                    "level": {
                        "name": "advance",
                        "min": 0.0,
                        "max": 1.0,
                        "step2level": {
                            "scale_start": 0.99999,
                            "scale_end": 0.00001,
                            "width": 3,
                        },
                    },
                },
            }
        )
    else:
        config.update(
            {
                "task": {
                    "name": "sbdd",
                    "transform": {"name": "ar", "part1_pert": "small"},
                },
                "noise": {
                    "name": "maskfill",
                    "num_steps": args.num_steps,
                    "ar_config": {
                        "strategy": "refine",
                        "r": 3,
                        "threshold_node": 0.98,
                        "threshold_pos": 0.91,
                        "threshold_bond": 0.98,
                        "max_ar_step": 10,
                        "change_init_step": 1,
                    },
                    "prior": {"part1": "from_train", "part2": "from_train"},
                    "level": {
                        "part1": {"name": "uniform", "min": 0.6, "max": 1.0},
                        "part2": {
                            "name": "advance",
                            "min": 0.0,
                            "max": 1.0,
                            "step2level": {
                                "scale_start": 0.99999,
                                "scale_end": 0.00001,
                                "width": 3,
                            },
                        },
                    },
                },
            }
        )
    return config


def _append_properties(source: Path, destination: Path, row: Dict[str, str]) -> None:
    body = source.read_text(encoding="utf-8", errors="replace")
    body = body.rsplit("$$$$", 1)[0].rstrip()
    properties = []
    for column in ("cfd_traj", "cfd_pos", "cfd_node", "cfd_edge"):
        value = row.get(column, "")
        if value:
            properties.append(f">  <POCKETXMOL_{column.upper()}>\n{value}\n")
    destination.write_text(
        body + "\n" + "\n".join(properties) + "\n$$$$\n", encoding="utf-8"
    )


def _collect_outputs(upstream: Path, destination: Path, limit: int) -> int:
    destination.mkdir(exist_ok=True)
    count = 0
    for info_path in sorted(upstream.glob("**/gen_info.csv")):
        log_dir = info_path.parent
        molecule_dir = log_dir / f"{log_dir.name}_SDF"
        with info_path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("tag", "").strip():
                    continue
                source = molecule_dir / row["filename"]
                if not source.is_file() or source.suffix.lower() != ".sdf":
                    continue
                count += 1
                _append_properties(
                    source, destination / f"candidate_{count:06d}.sdf", row
                )
                if count >= limit:
                    return count
    return count


def _write_upstream_launcher(path: Path) -> Path:
    path.write_text(_UPSTREAM_LAUNCHER, encoding="utf-8")
    return path


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.pocketxmol_root.expanduser().resolve()
    script = root / "scripts" / "sample_use.py"
    checkpoint = args.checkpoint.expanduser().resolve()
    train_config = args.train_config.expanduser().resolve()
    protein = args.protein.expanduser().resolve()
    if not script.is_file():
        raise SystemExit(f"PocketXMol sampling script not found: {script}")
    if not checkpoint.is_file():
        raise SystemExit(f"PocketXMol checkpoint not found: {checkpoint}")
    if not train_config.is_file():
        raise SystemExit(f"PocketXMol train config not found: {train_config}")
    if not protein.is_file() or protein.suffix.lower() != ".pdb":
        raise SystemExit(f"PocketXMol requires a protein PDB: {protein}")
    if args.num_samples <= 0 or args.batch_size <= 0 or args.num_steps <= 0:
        raise SystemExit("num-samples, batch-size, and num-steps must be positive")

    reference_value = args.reference_ligand.strip()
    reference = Path(reference_value).expanduser().resolve() if reference_value else None
    if reference is not None and not reference.is_file():
        raise SystemExit(f"PocketXMol reference ligand not found: {reference}")

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    runtime_model = output / "runtime_model"
    checkpoint_dir = runtime_model / "checkpoints"
    train_config_dir = runtime_model / "train_config"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    train_config_dir.mkdir(parents=True, exist_ok=True)
    safe_checkpoint = _sanitize_checkpoint(
        checkpoint, checkpoint_dir / "pocketxmol.ckpt"
    )
    shutil.copy2(train_config, train_config_dir / "train.yml")
    model_config = output / "prism_pocketxmol_model.yml"
    model_config.write_text(
        yaml.safe_dump({"model": {"checkpoint": str(safe_checkpoint)}}, sort_keys=False),
        encoding="utf-8",
    )
    task_config = output / "prism_pocketxmol_task.yml"
    task_config.write_text(
        yaml.safe_dump(_task_config(args, protein, reference), sort_keys=False),
        encoding="utf-8",
    )

    upstream = output / "upstream"
    upstream.mkdir(exist_ok=True)
    launcher = _write_upstream_launcher(output / "run_pocketxmol_upstream.py")
    command = [
        sys.executable,
        str(launcher),
        "--config_task",
        str(task_config),
        "--config_model",
        str(model_config),
        "--outdir",
        str(upstream),
        "--device",
        args.device,
        "--batch_size",
        str(min(args.batch_size, args.num_samples)),
        "--num_workers",
        "0",
    ]
    environment = os.environ.copy()
    environment["PRISM_POCKETXMOL_ROOT"] = str(root)
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(root)
        if not existing_pythonpath
        else str(root) + os.pathsep + existing_pythonpath
    )
    completed = subprocess.run(command, cwd=str(root), env=environment, check=False)
    if completed.returncode != 0:
        return completed.returncode
    count = _collect_outputs(upstream, output / "sdf", args.num_samples)
    if count == 0:
        print("PocketXMol produced no successful SDF molecules", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
