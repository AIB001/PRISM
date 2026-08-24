"""Stable PRISM-to-DiffSBDD custom-pocket inference contract.

DiffSBDD is an equivariant diffusion model trained on CrossDocked, the same
corpus as TargetDiff and Pocket2Mol.  Its official entry point,
``generate_ligands.py``, accepts a reference ligand either as ``<chain>:<resi>``
inside the receptor PDB or as a standalone SDF path, so PRISM's normalized
protein/reference layout maps onto it directly.

Two upstream capabilities are deliberately left off by default:
``--sanitize`` discards molecules that fail RDKit sanitization, and ``--relax``
runs a force-field relaxation.  Enabling either would silently change what the
model is measured on, so PRISM exposes them as explicit opt-ins and reports
raw model output otherwise.
"""

import argparse
import collections
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

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

_UPSTREAM_LAUNCHER = r'''"""Run the unmodified DiffSBDD entry point with a fixed seed."""
import json
import os
import random
import runpy
import sys
import types
from pathlib import Path

root = Path(os.environ["PRISM_DIFFSBDD_ROOT"]).resolve()
sys.path.insert(0, str(root))
os.chdir(root)

# lightning_modules pulls in two module-scope dependencies that only the
# training and visualisation paths use: wandb for experiment tracking, and
# imageio through analysis.visualization for writing trajectory GIFs. Every
# generation-path call site passes wandb=None and never renders a frame, so
# stubbing both keeps them out of the inference runtime without touching the
# model, the sampler, or the molecule builder.
for _optional in ("wandb", "imageio"):
    sys.modules.setdefault(_optional, types.ModuleType(_optional))

import torch

seed = int(os.environ["PRISM_DIFFSBDD_SEED"])
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
import numpy as np
np.random.seed(seed)
random.seed(seed)

script = root / "generate_ligands.py"
sys.argv = [str(script)] + json.loads(os.environ["PRISM_DIFFSBDD_ARGS"])
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
    """Load the published checkpoint through an allowlist and rewrite it."""
    # torch < 2.0 silently discards ``pickle_module``, so the allowlist is
    # enforced here by walking the pickle opcodes before torch touches the
    # file. Nothing in the checkpoint executes until this passes.
    _checkpoint.verify(Path(source), extra_allowed={("argparse", "Namespace")})

    import torch

    _SAFE_PICKLE_GLOBALS.clear()
    _SAFE_PICKLE_GLOBALS.update(
        {
            ("collections", "OrderedDict"): collections.OrderedDict,
            ("torch._utils", "_rebuild_tensor_v2"): torch._utils._rebuild_tensor_v2,
        }
    )
    for storage in (
        "BFloat16Storage", "BoolStorage", "ByteStorage", "CharStorage",
        "DoubleStorage", "FloatStorage", "HalfStorage", "IntStorage",
        "LongStorage", "ShortStorage",
    ):
        if hasattr(torch, storage):
            _SAFE_PICKLE_GLOBALS[("torch", storage)] = getattr(torch, storage)

    checkpoint = torch.load(
        str(source), map_location="cpu", pickle_module=_RestrictedPickleModule
    )
    if not isinstance(checkpoint, dict) or "state_dict" not in checkpoint:
        raise SystemExit("DiffSBDD checkpoint must contain a 'state_dict' mapping")
    if not isinstance(checkpoint.get("hyper_parameters", {}), dict):
        raise SystemExit("DiffSBDD checkpoint 'hyper_parameters' must be a mapping")
    torch.save(checkpoint, str(output))
    return output


def _device(value: str) -> str:
    """Resolve --device against the host; ``auto`` falls back to CPU."""
    return _device_support.resolve(value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run official DiffSBDD ligand generation")
    parser.add_argument("--diffsbdd-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument(
        "--reference-ligand",
        type=Path,
        required=True,
        help="SDF defining the pocket; upstream reads it directly",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--num-samples", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", type=_device, default="auto")
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument(
        "--timesteps",
        type=int,
        default=None,
        help="Denoising steps; upstream default is the model's trained schedule",
    )
    parser.add_argument(
        "--num-nodes-lig",
        type=int,
        default=None,
        help=(
            "Fix the heavy-atom count. Left unset, the model samples the size "
            "from the distribution it learned for a pocket of this volume"
        ),
    )
    parser.add_argument(
        "--upstream-sanitize",
        action="store_true",
        help="Let upstream discard molecules that fail RDKit sanitization (off by default)",
    )
    parser.add_argument(
        "--upstream-relax",
        action="store_true",
        help="Let upstream force-field relax each molecule (off by default)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.diffsbdd_root.expanduser().resolve()
    script = root / "generate_ligands.py"
    if not script.is_file():
        raise SystemExit(f"DiffSBDD entry point not found: {script}")
    checkpoint = args.checkpoint.expanduser().resolve()
    protein = args.protein.expanduser().resolve()
    reference = args.reference_ligand.expanduser().resolve()
    if not checkpoint.is_file():
        raise SystemExit(f"DiffSBDD checkpoint not found: {checkpoint}")
    if not protein.is_file() or protein.suffix.lower() != ".pdb":
        raise SystemExit(f"DiffSBDD requires a protein PDB: {protein}")
    if not reference.is_file() or reference.suffix.lower() != ".sdf":
        raise SystemExit(f"DiffSBDD requires a reference ligand SDF: {reference}")
    if args.num_samples <= 0 or args.batch_size <= 0:
        raise SystemExit("num-samples and batch-size must be positive")
    # Upstream asserts n_samples % batch_size == 0.
    batch_size = args.batch_size
    if args.num_samples % batch_size:
        batch_size = args.num_samples
        print(
            f"DiffSBDD: batch size {args.batch_size} does not divide "
            f"{args.num_samples}; using {batch_size}",
            file=sys.stderr,
        )

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    trusted = _sanitize_checkpoint(checkpoint, output / "trusted_checkpoint.ckpt")

    sdf_output = output / "sdf"
    sdf_output.mkdir(exist_ok=True)
    outfile = sdf_output / "candidates.sdf"

    upstream_args = [
        str(trusted),
        "--pdbfile", str(protein),
        "--ref_ligand", str(reference),
        "--outfile", str(outfile),
        "--n_samples", str(args.num_samples),
        "--batch_size", str(batch_size),
    ]
    if args.timesteps is not None:
        upstream_args += ["--timesteps", str(args.timesteps)]
    if args.num_nodes_lig is not None:
        upstream_args += ["--num_nodes_lig", str(args.num_nodes_lig)]
    if args.upstream_sanitize:
        upstream_args.append("--sanitize")
    if args.upstream_relax:
        upstream_args.append("--relax")

    launcher = output / "run_diffsbdd_upstream.py"
    launcher.write_text(_UPSTREAM_LAUNCHER, encoding="utf-8")

    environment = os.environ.copy()
    environment["PRISM_DIFFSBDD_ROOT"] = str(root)
    environment["PRISM_DIFFSBDD_SEED"] = str(args.seed)
    environment["PRISM_DIFFSBDD_ARGS"] = json.dumps(upstream_args)
    existing = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = str(root) if not existing else str(root) + os.pathsep + existing
    if args.device.startswith("cuda:"):
        environment["CUDA_VISIBLE_DEVICES"] = args.device.split(":", 1)[1]
    elif args.device == "cpu":
        environment["CUDA_VISIBLE_DEVICES"] = ""

    completed = subprocess.run(
        [sys.executable, str(launcher)], cwd=str(root), env=environment, check=False
    )
    if completed.returncode != 0:
        return completed.returncode

    if not outfile.is_file() or not outfile.read_text(encoding="utf-8").strip():
        print("DiffSBDD produced no SDF molecules", file=sys.stderr)
        return 2

    records = outfile.read_text(encoding="utf-8").count("$$$$")
    (output / "run_summary.json").write_text(
        json.dumps(
            {
                "requested": args.num_samples,
                "records_written": records,
                "batch_size": batch_size,
                "timesteps": args.timesteps,
                "num_nodes_lig": args.num_nodes_lig,
                "upstream_sanitize": args.upstream_sanitize,
                "upstream_relax": args.upstream_relax,
                "seed": args.seed,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        f"DiffSBDD: requested {args.num_samples}, wrote {records} record(s)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
