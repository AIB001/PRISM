"""Stable PRISM-to-MolCRAFT custom-pocket inference contract."""

import argparse
import collections
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
    from . import _checkpoint
    from . import _device as _device_support
except ImportError:  # standalone script invocation
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _checkpoint  # type: ignore[no-redef]
    import _device as _device_support  # type: ignore[no-redef]


_SAFE_PICKLE_GLOBALS: Dict[Tuple[str, str], Any] = {}

_UPSTREAM_LAUNCHER = r'''"""Run the unmodified MolCRAFT sampler in an isolated work directory."""
import os
import sys
import types
from pathlib import Path

root = Path(os.environ["PRISM_MOLCRAFT_ROOT"]).resolve()
runtime = Path(os.environ["PRISM_MOLCRAFT_RUNTIME"]).resolve()
output = Path(os.environ["PRISM_MOLCRAFT_OUTPUT"]).resolve()
sys.path.insert(0, str(root))
os.chdir(runtime)

# The upstream sampling file imports optional docking/PoseCheck helpers for a
# post-processing class that PRISM never calls. Keep those evaluation packages
# out of the generation runtime while preserving the model and sampler itself.
class DisabledEvaluation:
    def __init__(self, *args, **kwargs):
        raise RuntimeError("Optional MolCRAFT evaluation is disabled by PRISM")


posecheck_stub = types.ModuleType("posecheck")
posecheck_stub.PoseCheck = DisabledEvaluation
sys.modules["posecheck"] = posecheck_stub
docking_stub = types.ModuleType("core.evaluation.docking_vina")
docking_stub.VinaDockingTask = DisabledEvaluation
sys.modules["core.evaluation.docking_vina"] = docking_stub
metrics_stub = types.ModuleType("core.evaluation.metrics")
metrics_stub.CondMolGenMetric = type(
    "UnusedCondMolGenMetric", (), {"__init__": lambda self, *args, **kwargs: None}
)
sys.modules["core.evaluation.metrics"] = metrics_stub
wandb_stub = types.ModuleType("wandb")
sys.modules["wandb"] = wandb_stub
visualization_stub = types.ModuleType("core.evaluation.visualization")
visualization_stub.visualize = lambda *args, **kwargs: None
visualization_stub.visualize_chain = lambda *args, **kwargs: None
sys.modules["core.evaluation.visualization"] = visualization_stub

import core.callbacks.validation_callback_for_sample as sample_callback
import sample_for_pocket as upstream

sample_callback.OUT_DIR = str(output)
original_config = upstream.Config


def runtime_config(path):
    cfg = original_config(path)
    cfg.seed = int(os.environ["PRISM_MOLCRAFT_SEED"])
    cfg.evaluation.batch_size = int(os.environ["PRISM_MOLCRAFT_BATCH_SIZE"])
    cfg.accounting.logdir = str(runtime / "logs")
    cfg.accounting.test_outputs_dir = str(runtime / "test_outputs")
    return cfg


upstream.Config = runtime_config
upstream.call(
    os.environ["PRISM_MOLCRAFT_PROTEIN"],
    os.environ["PRISM_MOLCRAFT_REFERENCE"],
    ckpt_path=os.environ["PRISM_MOLCRAFT_CHECKPOINT"],
    num_samples=int(os.environ["PRISM_MOLCRAFT_NUM_SAMPLES"]),
    sample_steps=int(os.environ["PRISM_MOLCRAFT_SAMPLE_STEPS"]),
    sample_num_atoms=os.environ["PRISM_MOLCRAFT_SAMPLE_NUM_ATOMS"],
    seed=int(os.environ["PRISM_MOLCRAFT_SEED"]),
)
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


def _sanitize_checkpoint(source: Path, output: Path, config_path: Path) -> Path:
    """Restricted-load the official checkpoint and materialize its embedded config."""
    # torch < 2.0 silently discards ``pickle_module``, so the allowlist is
    # enforced here by walking the pickle opcodes before torch touches the
    # file. Nothing in the checkpoint executes until this passes.
    _checkpoint.verify(Path(source), extra_allowed=set())

    import torch

    allowed = {
        ("collections", "OrderedDict"): collections.OrderedDict,
        ("torch._utils", "_rebuild_tensor_v2"): torch._utils._rebuild_tensor_v2,
    }
    for storage_name in (
        "BFloat16Storage",
        "BoolStorage",
        "ByteStorage",
        "CharStorage",
        "DoubleStorage",
        "FloatStorage",
        "HalfStorage",
        "IntStorage",
        "LongStorage",
        "ShortStorage",
    ):
        storage = getattr(torch, storage_name, None)
        if storage is not None:
            allowed[("torch", storage_name)] = storage
    _SAFE_PICKLE_GLOBALS.clear()
    _SAFE_PICKLE_GLOBALS.update(allowed)
    checkpoint = torch.load(
        str(source), map_location="cpu", pickle_module=_RestrictedPickleModule
    )
    if not isinstance(checkpoint, dict):
        raise SystemExit("MolCRAFT checkpoint must be a mapping")
    state_dict = checkpoint.get("state_dict")
    hyper_parameters = checkpoint.get("hyper_parameters")
    if not isinstance(state_dict, (dict, collections.OrderedDict)):
        raise SystemExit("MolCRAFT checkpoint must contain a state_dict mapping")
    if not isinstance(hyper_parameters, dict):
        raise SystemExit("MolCRAFT checkpoint must contain hyper_parameters")

    runtime_config = dict(hyper_parameters)
    dynamics = dict(runtime_config.get("dynamics", {}))
    runtime_config["dynamics"] = dynamics
    ligand_embedding = state_dict.get("dynamics.ligand_atom_emb.weight")
    try:
        checkpoint_input_dim = int(ligand_embedding.shape[1])
        ligand_feature_dim = int(dynamics["ligand_atom_feature_dim"])
        time_embedding_dim = int(dynamics.get("time_emb_dim", 0))
    except (AttributeError, KeyError, TypeError, ValueError, IndexError) as exc:
        raise SystemExit("MolCRAFT checkpoint has an invalid ligand embedding contract") from exc
    configured_input_dim = ligand_feature_dim + time_embedding_dim
    if checkpoint_input_dim == ligand_feature_dim and time_embedding_dim > 0:
        # The published checkpoint predates the current sampler's concatenated
        # time feature even though its embedded config says time_emb_dim=1.
        # Shape-derived correction is required to instantiate the exact model
        # represented by the published weights.
        dynamics["time_emb_dim"] = 0
        print(
            "MolCRAFT compatibility: disabled the config-only time embedding "
            f"to match checkpoint input width {checkpoint_input_dim}",
            file=sys.stderr,
        )
    elif checkpoint_input_dim != configured_input_dim:
        raise SystemExit(
            "MolCRAFT ligand embedding shape is incompatible with its embedded config: "
            f"checkpoint={checkpoint_input_dim}, config={configured_input_dim}"
        )
    accounting = dict(runtime_config.get("accounting", {}))
    runtime_config["accounting"] = accounting
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        yaml.safe_dump(runtime_config, sort_keys=False), encoding="utf-8"
    )
    torch.save(checkpoint, str(output))
    return output


def _device(value: str) -> str:
    """Resolve --device against the host. This model has no working CPU path."""
    return _device_support.resolve(
        value,
        allow_cpu=False,
        cpu_unsupported_message=(
            "MolCRAFT's official sampler builds a Lightning Trainer without an accelerator setting, so its device resolves to CUDA regardless"
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run official MolCRAFT pocket sampling")
    parser.add_argument("--molcraft-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--reference-ligand", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--num-samples", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", type=_device, default="auto")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--sample-steps", type=int, default=100)
    parser.add_argument("--sample-num-atoms", choices=("prior", "ref"), default="prior")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.molcraft_root.expanduser().resolve()
    source_script = root / "sample_for_pocket.py"
    checkpoint = args.checkpoint.expanduser().resolve()
    protein = args.protein.expanduser().resolve()
    reference = args.reference_ligand.expanduser().resolve()
    if not source_script.is_file():
        raise SystemExit(f"MolCRAFT sampling script not found: {source_script}")
    if not checkpoint.is_file():
        raise SystemExit(f"MolCRAFT checkpoint not found: {checkpoint}")
    if not protein.is_file() or protein.suffix.lower() != ".pdb":
        raise SystemExit(f"MolCRAFT requires a protein/pocket PDB: {protein}")
    if not reference.is_file() or reference.suffix.lower() != ".sdf":
        raise SystemExit(f"MolCRAFT requires a reference SDF: {reference}")
    if args.num_samples <= 0 or args.batch_size <= 0 or args.sample_steps <= 0:
        raise SystemExit("num-samples, batch-size, and sample-steps must be positive")

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    runtime = output / "runtime"
    checkpoint_dir = runtime / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    safe_checkpoint = _sanitize_checkpoint(
        checkpoint,
        checkpoint_dir / "last.ckpt",
        checkpoint_dir / "config.yaml",
    )
    upstream_output = output / "upstream"
    upstream_output.mkdir(exist_ok=True)
    launcher = output / "run_molcraft_upstream.py"
    launcher.write_text(_UPSTREAM_LAUNCHER, encoding="utf-8")

    environment = os.environ.copy()
    environment_lib = str(Path(sys.prefix) / "lib")
    inherited_library_path = environment.get("LD_LIBRARY_PATH")
    environment["LD_LIBRARY_PATH"] = (
        environment_lib
        if not inherited_library_path
        else environment_lib + os.pathsep + inherited_library_path
    )
    environment.update(
        {
            "PRISM_MOLCRAFT_ROOT": str(root),
            "PRISM_MOLCRAFT_RUNTIME": str(runtime),
            "PRISM_MOLCRAFT_OUTPUT": str(upstream_output),
            "PRISM_MOLCRAFT_CHECKPOINT": str(safe_checkpoint),
            "PRISM_MOLCRAFT_PROTEIN": str(protein),
            "PRISM_MOLCRAFT_REFERENCE": str(reference),
            "PRISM_MOLCRAFT_NUM_SAMPLES": str(args.num_samples),
            "PRISM_MOLCRAFT_SEED": str(args.seed),
            "PRISM_MOLCRAFT_BATCH_SIZE": str(min(args.batch_size, args.num_samples)),
            "PRISM_MOLCRAFT_SAMPLE_STEPS": str(args.sample_steps),
            "PRISM_MOLCRAFT_SAMPLE_NUM_ATOMS": args.sample_num_atoms,
        }
    )
    if args.device == "cpu":
        environment["CUDA_VISIBLE_DEVICES"] = ""
    elif args.device.startswith("cuda:"):
        environment["CUDA_VISIBLE_DEVICES"] = args.device.split(":", 1)[1]

    completed = subprocess.run(
        [sys.executable, str(launcher)],
        cwd=str(runtime),
        env=environment,
        check=False,
    )
    if completed.returncode != 0:
        return completed.returncode

    sdf_output = output / "sdf"
    sdf_output.mkdir(exist_ok=True)
    sources = sorted(upstream_output.glob("*.sdf"), key=lambda path: path.name)
    for index, source in enumerate(sources[: args.num_samples], start=1):
        shutil.copy2(source, sdf_output / f"candidate_{index:06d}.sdf")
    if not sources:
        print("MolCRAFT produced no SDF molecules", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
