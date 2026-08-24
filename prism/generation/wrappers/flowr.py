"""Stable PRISM-to-FLOWR PDB generation contract."""

import argparse
import collections
import os
import pickle
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path
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
from typing import Any, Dict, Optional, Sequence, Tuple


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
    """Restricted-load the official checkpoint before invoking Lightning."""
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
        raise SystemExit("FLOWR checkpoint must be a mapping")
    if not isinstance(checkpoint.get("state_dict"), dict):
        raise SystemExit("FLOWR checkpoint must contain a state_dict mapping")
    if not isinstance(checkpoint.get("hyper_parameters"), dict):
        raise SystemExit("FLOWR checkpoint must contain hyper_parameters")
    torch.save(checkpoint, str(output))
    return output


def _device(value: str) -> str:
    """Resolve --device against the host. FLOWR's generator is CUDA-only."""
    return _device_support.resolve(
        value,
        allow_cpu=False,
        cpu_unsupported_message="FLOWR's official generator requires CUDA",
    )


def _write_launcher(path: Path) -> Path:
    """Write an isolated launcher that omits unused optional evaluation stacks."""
    source = textwrap.dedent(
        """\
        import runpy
        import sys
        import types

        import torch
        from torchmetrics import Metric


        class _NullMetric(Metric):
            full_state_update = False

            def __init__(self, *args, **kwargs):
                super().__init__()
                self.add_state("count", default=torch.tensor(0), dist_reduce_fx="sum")

            def update(self, *args, **kwargs):
                self.count += 1

            def compute(self):
                return torch.tensor(float("nan"), device=self.count.device)


        metrics = types.ModuleType("flowr.util.metrics")
        for name in (
            "AtomStability",
            "AverageEnergy",
            "AverageOptRmsd",
            "AverageStrainEnergy",
            "DistributionDistance",
            "EnergyValidity",
            "GenBench3DSB",
            "GenBench3DValidity",
            "MolecularAccuracy",
            "MolecularPairRMSD",
            "MoleculeStability",
            "Novelty",
            "PoseBustersValidity",
            "Uniqueness",
            "Validity",
        ):
            setattr(metrics, name, _NullMetric)
        metrics.ALLOWED_VALENCIES = {}
        metrics._is_valid_valence = lambda *args, **kwargs: False

        def _unused_interaction_recovery(*args, **kwargs):
            raise RuntimeError(
                "PRISM FLOWR generation does not enable optional interaction recovery"
            )

        metrics.interaction_recovery_per_complex = _unused_interaction_recovery
        sys.modules["flowr.util.metrics"] = metrics

        pymol = types.ModuleType("pymol")
        pymol.cmd = None
        sys.modules["pymol"] = pymol

        class _UnavailableUnit:
            def __mul__(self, other):
                return self

            __rmul__ = __mul__

            def __truediv__(self, other):
                return self

            __rtruediv__ = __truediv__

            def __pow__(self, other):
                return self


        class _UnavailablePDBFile:
            def __init__(self, *args, **kwargs):
                raise RuntimeError(
                    "PRISM FLOWR generation does not enable optional OpenMM optimization"
                )


        openmm = types.ModuleType("openmm")
        openmm_app = types.ModuleType("openmm.app")
        openmm_app.PDBFile = _UnavailablePDBFile
        openmm.app = openmm_app
        openmm.unit = types.SimpleNamespace(
            angstroms=_UnavailableUnit(),
            kilocalories_per_mole=_UnavailableUnit(),
            kilojoule_per_mole=_UnavailableUnit(),
            nanometer=_UnavailableUnit(),
            picoseconds=_UnavailableUnit(),
        )
        sys.modules["openmm"] = openmm
        sys.modules["openmm.app"] = openmm_app

        pdbinf = types.ModuleType("pdbinf")
        pdbinf.STANDARD_AA_DOC = {}
        pdbinf._pdbinf = types.SimpleNamespace()
        sys.modules["pdbinf"] = pdbinf

        runpy.run_module("flowr.gen.generate_from_pdb", run_name="__main__")
        """
    )
    path.write_text(source, encoding="utf-8")
    return path


def _prepare_protein(source: Path, runtime: Path) -> Path:
    """Drop PDB CONECT records that can reference atoms removed upstream."""
    if source.suffix.lower() != ".pdb":
        return source
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines(True)
    if not any(line.startswith("CONECT") for line in lines):
        return source
    output = runtime / source.name
    output.write_text(
        "".join(line for line in lines if not line.startswith("CONECT")),
        encoding="utf-8",
    )
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run official FLOWR PDB sampling")
    parser.add_argument("--flowr-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--reference-ligand", type=Path, required=True)
    parser.add_argument("--pocket-mode", choices=("reference", "prepared"), required=True)
    parser.add_argument("--pocket-cutoff", type=float, default=6.0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--num-samples", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", type=_device, default="auto")
    parser.add_argument("--batch-cost", type=int, default=100)
    parser.add_argument("--max-sample-iter", type=int, default=20)
    parser.add_argument("--num-workers", type=int, default=4)
    parser.add_argument("--integration-steps", type=int, default=100)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.flowr_root.expanduser().resolve()
    module = root / "flowr" / "gen" / "generate_from_pdb.py"
    checkpoint = args.checkpoint.expanduser().resolve()
    protein = args.protein.expanduser().resolve()
    reference = args.reference_ligand.expanduser().resolve()
    if not module.is_file():
        raise SystemExit(f"FLOWR generation module not found: {module}")
    if not checkpoint.is_file():
        raise SystemExit(f"FLOWR checkpoint not found: {checkpoint}")
    if not protein.is_file() or protein.suffix.lower() not in {".pdb", ".cif"}:
        raise SystemExit(f"FLOWR requires a PDB/CIF protein or pocket: {protein}")
    if not reference.is_file() or reference.suffix.lower() not in {".sdf", ".pdb"}:
        raise SystemExit(f"FLOWR requires a reference SDF or PDB: {reference}")
    if (
        args.num_samples <= 0
        or args.batch_cost <= 0
        or args.max_sample_iter <= 0
        or args.num_workers <= 0
        or args.integration_steps <= 0
        or args.pocket_cutoff <= 0
    ):
        raise SystemExit("FLOWR numeric sampling parameters must be positive")

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    runtime = output / "runtime"
    runtime.mkdir(exist_ok=True)
    safe_checkpoint = _sanitize_checkpoint(checkpoint, runtime / "flowr.ckpt")
    prepared_protein = _prepare_protein(protein, runtime)
    launcher = _write_launcher(runtime / "launch_flowr.py")
    upstream = output / "upstream"
    upstream.mkdir(exist_ok=True)
    command = [
        sys.executable,
        str(launcher),
        "--pdb_file",
        str(prepared_protein),
        "--ligand_file",
        str(reference),
        "--pocket_cutoff",
        str(args.pocket_cutoff),
        "--gpus",
        "1",
        "--num_workers",
        str(args.num_workers),
        "--batch_cost",
        str(args.batch_cost),
        "--arch",
        "pocket",
        "--pocket_type",
        "holo",
        "--ckpt_path",
        str(safe_checkpoint),
        "--save_dir",
        str(upstream),
        "--max_sample_iter",
        str(args.max_sample_iter),
        "--sample_n_molecules_per_target",
        str(args.num_samples),
        "--seed",
        str(args.seed),
        "--categorical_strategy",
        "uniform-sample",
        "--integration_steps",
        str(args.integration_steps),
        "--filter_valid_unique",
        "--sample_mol_sizes",
    ]
    if args.pocket_mode == "reference":
        command.append("--cut_pocket")

    environment = os.environ.copy()
    environment_lib = str(Path(sys.prefix) / "lib")
    inherited_library_path = environment.get("LD_LIBRARY_PATH")
    environment["LD_LIBRARY_PATH"] = (
        environment_lib
        if not inherited_library_path
        else environment_lib + os.pathsep + inherited_library_path
    )
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(root)
        if not existing_pythonpath
        else str(root) + os.pathsep + existing_pythonpath
    )
    environment["CUDA_VISIBLE_DEVICES"] = args.device.split(":", 1)[1]
    completed = subprocess.run(
        command, cwd=str(root), env=environment, check=False
    )
    if completed.returncode != 0:
        return completed.returncode

    source = upstream / f"samples_{protein.stem}.sdf"
    if not source.is_file():
        print(f"FLOWR produced no expected SDF file: {source}", file=sys.stderr)
        return 2
    sdf_output = output / "sdf"
    sdf_output.mkdir(exist_ok=True)
    shutil.copy2(source, sdf_output / "candidates.sdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
