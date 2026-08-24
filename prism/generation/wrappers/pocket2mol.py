"""Stable PRISM-to-Pocket2Mol command contract."""

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


def _natural_key(path: Path) -> Tuple[Any, ...]:
    """Sort ``2.sdf`` before ``10.sdf`` by comparing digit runs numerically."""
    return tuple(
        int(part) if part.isdigit() else part
        for part in re.split(r"(\d+)", path.name)
    )


# Pocket2Mol's sampling loop is ``while len(pool.finished) < num_samples`` and
# its step counter is the number of atoms added, so a molecule that finishes at
# step k has roughly k heavy atoms.  Stopping as soon as enough molecules have
# accumulated therefore stops at small sizes: on 1IEP the run reached step 22 of
# 50 and the largest molecule had exactly 22 heavy atoms, against 37 for the
# crystal ligand.  PRISM removes ``num_samples`` from the stopping condition by
# passing this sentinel, lets ``max_steps`` bound the size, and selects the
# requested count afterwards.
_UNBOUNDED_SAMPLES = 10**9

SELECTION_MODES = ("stratified", "largest", "target", "completion")


def _heavy_atom_count(path: Path) -> Optional[int]:
    """Atom count from the V2000 counts line; Pocket2Mol emits heavy atoms only."""
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return None
    if len(lines) < 4:
        return None
    try:
        return int(lines[3][:3])
    except ValueError:
        return None


def _select_indices(
    counts: Sequence[Optional[int]], wanted: int, mode: str, target: Optional[int]
) -> list:
    """Choose which finished molecules to deliver, by an explicit criterion."""
    total = len(counts)
    if wanted >= total:
        return list(range(total))
    if mode == "completion":
        # Reproduces the upstream stopping behavior, including its size bias.
        return list(range(wanted))

    def size_key(index):
        return (counts[index] is None, counts[index] or 0, index)

    if mode == "largest":
        return sorted(sorted(range(total), key=size_key)[-wanted:])
    if mode == "target":
        if target is None:
            raise SystemExit("--select target requires --target-heavy-atoms")
        ranked = sorted(
            range(total),
            key=lambda i: (counts[i] is None, abs((counts[i] or 0) - target), i),
        )
        return sorted(ranked[:wanted])

    # Stratified: walk the size-sorted list at even intervals so the delivered
    # set spans the model's full size distribution instead of its lower tail.
    ordered = sorted(range(total), key=size_key)
    step = (total - 1) / (wanted - 1) if wanted > 1 else 0.0
    chosen = []
    seen = set()
    for position in range(wanted):
        index = ordered[min(total - 1, round(position * step))]
        if index not in seen:
            seen.add(index)
            chosen.append(index)
    # Rounding can land on the same molecule twice; top up from what is left.
    for index in ordered:
        if len(chosen) >= wanted:
            break
        if index not in seen:
            seen.add(index)
            chosen.append(index)
    return sorted(chosen)


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
    """Load a checkpoint through a strict pickle allowlist and rewrite it."""
    # torch < 2.0 silently discards ``pickle_module``, so the allowlist is
    # enforced here by walking the pickle opcodes before torch touches the
    # file. Nothing in the checkpoint executes until this passes.
    _checkpoint.verify(Path(source), extra_allowed={("easydict", "EasyDict")})

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
        str(source), map_location="cpu", pickle_module=_RestrictedPickleModule
    )
    if not isinstance(checkpoint, dict) or not {"config", "model"}.issubset(checkpoint):
        raise SystemExit("Pocket2Mol checkpoint must contain 'config' and 'model' mappings")
    if not isinstance(checkpoint["model"], (dict, collections.OrderedDict)):
        raise SystemExit("Pocket2Mol checkpoint 'model' must be a state dictionary")
    torch.save(checkpoint, str(output))
    return output


def _recover_finished_from_pool(upstream_output: Path, pocket2mol_root: Path) -> int:
    """Write SDFs from the last pool checkpoint after an upstream crash.

    ``sample_for_pdb.py`` raises ``ValueError: probabilities do not sum to 1``
    inside ``np.random.choice`` once its beam queue drains to a few candidates,
    because ``logp_to_rank_prob`` stops normalizing exactly.  Running to
    ``max_steps`` instead of stopping at a sample count makes that state
    reachable, and upstream only guards the loop against ``KeyboardInterrupt``
    -- so the process dies before its "Save sdf mols" block and every finished
    molecule is lost.  The loop does checkpoint the pool each step, so recover
    from that instead of losing the run.

    Returns the number of molecules written, or 0 when nothing is recoverable.
    """
    candidates = [
        path
        for path in upstream_output.glob("**/samples_*.pt")
        if path.stem != "samples_init"
    ]
    if not candidates:
        return 0
    try:
        import torch
        from rdkit import Chem
    except ImportError as exc:
        # Recovery is best-effort; it must never mask the original failure.
        print(f"Pocket2Mol: cannot recover without {exc.name}", file=sys.stderr)
        return 0

    # The pool pickles Pocket2Mol's own ``utils.*`` data classes. The subprocess
    # gets the source tree through PYTHONPATH; this process needs it too, or the
    # unpickler fails with "No module named 'utils'".
    root_path = str(pocket2mol_root)
    if root_path not in sys.path:
        sys.path.insert(0, root_path)

    def checkpoint_rank(path: Path):
        # ``samples_all.pt`` is the complete final dump; else take the last step.
        if path.stem == "samples_all":
            return (1, 0)
        suffix = path.stem.rsplit("_", 1)[-1]
        return (0, int(suffix)) if suffix.isdigit() else (-1, 0)

    latest = max(candidates, key=checkpoint_rank)
    # This file is PRISM's own output from this run, written moments ago inside
    # the run directory. The restricted allowlist guards untrusted *model
    # checkpoints*; it does not apply to our own intermediate state.
    try:
        pool = torch.load(str(latest), map_location="cpu", weights_only=False)
    except Exception as exc:  # noqa: BLE001 - recovery must not mask the original failure
        print(f"Pocket2Mol: could not read {latest.name}: {exc}", file=sys.stderr)
        return 0

    finished = pool.get("finished") if isinstance(pool, dict) else getattr(pool, "finished", None)
    if not finished:
        return 0

    sdf_dir = latest.parent / "SDF"
    sdf_dir.mkdir(exist_ok=True)
    written = 0
    for entry in finished:
        molecule = getattr(entry, "rdmol", None)
        if molecule is None:
            continue
        try:
            Chem.MolToMolFile(molecule, str(sdf_dir / f"{written}.sdf"))
        except Exception:  # noqa: BLE001 - skip a single unwritable molecule
            continue
        written += 1
    print(
        f"Pocket2Mol: upstream sampling crashed; recovered {written} finished "
        f"molecule(s) from {latest.name}",
        file=sys.stderr,
    )
    return written


def _device(value: str) -> str:
    """Resolve --device against the host. This model has no working CPU path."""
    return _device_support.resolve(
        value,
        allow_cpu=False,
        cpu_unsupported_message=(
            "Pocket2Mol's official sampler hardcodes .to('cuda') in models/maskfill.py, so it cannot run on CPU"
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run official Pocket2Mol PDB sampling")
    parser.add_argument("--pocket2mol-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--center", type=float, nargs=3, required=True)
    parser.add_argument("--bbox-size", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--num-samples", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", type=_device, default="auto")
    parser.add_argument("--beam-size", type=int, default=300)
    parser.add_argument(
        "--max-steps",
        type=int,
        default=50,
        help=(
            "Atom-growth steps. This, not --num-samples, bounds molecule size: "
            "a molecule finishing at step k has roughly k heavy atoms. Runtime "
            "is about beam-size * 1.5 s per step"
        ),
    )
    parser.add_argument(
        "--select",
        choices=SELECTION_MODES,
        default="stratified",
        help=(
            "How to choose --num-samples molecules from the full run. "
            "stratified (default): span the whole size distribution. "
            "largest: the biggest molecules. "
            "target: closest to --target-heavy-atoms. "
            "completion: upstream order, which is biased towards small molecules"
        ),
    )
    parser.add_argument(
        "--target-heavy-atoms",
        type=int,
        default=None,
        help="Desired heavy-atom count for --select target (e.g. the reference ligand's)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.pocket2mol_root.expanduser().resolve()
    script = root / "sample_for_pdb.py"
    if not script.is_file():
        raise SystemExit(f"Pocket2Mol sampling script not found: {script}")
    checkpoint = args.checkpoint.expanduser().resolve()
    protein = args.protein.expanduser().resolve()
    if not checkpoint.is_file():
        raise SystemExit(f"Pocket2Mol checkpoint not found: {checkpoint}")
    if not protein.is_file() or protein.suffix.lower() != ".pdb":
        raise SystemExit(f"Pocket2Mol requires a protein PDB: {protein}")
    if (
        args.num_samples <= 0
        or args.bbox_size <= 0
        or args.beam_size <= 0
        or args.max_steps <= 0
    ):
        raise SystemExit(
            "num-samples, bbox-size, beam-size, and max-steps must be positive"
        )

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    safe_checkpoint = _sanitize_checkpoint(
        checkpoint, output / "trusted_checkpoint.pt"
    )
    config = {
        "model": {"checkpoint": str(safe_checkpoint)},
        "sample": {
            "seed": args.seed,
            # Deliberately not args.num_samples -- see _UNBOUNDED_SAMPLES.
            "num_samples": (
                args.num_samples
                if args.select == "completion"
                else _UNBOUNDED_SAMPLES
            ),
            "beam_size": args.beam_size,
            "max_steps": args.max_steps,
            "threshold": {
                "focal_threshold": 0.5,
                "pos_threshold": 0.25,
                "element_threshold": 0.3,
                "hasatom_threshold": 0.6,
                "bond_threshold": 0.4,
            },
        },
    }
    config_path = output / "prism_pocket2mol_sampling.yml"
    config_path.write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    upstream_output = output / "upstream"
    upstream_output.mkdir(exist_ok=True)
    center = ",".join(str(value) for value in args.center)
    command = [
        sys.executable,
        str(script),
        "--pdb_path",
        str(protein),
        # Pocket2Mol parses the center as one comma-separated string.  Use the
        # ``--option=value`` form so a negative first coordinate is not
        # mistaken for another command-line option by its argparse parser.
        f"--center={center}",
        "--bbox_size",
        str(args.bbox_size),
        "--config",
        str(config_path),
        "--device",
        args.device,
        "--outdir",
        str(upstream_output),
    ]
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(root)
        if not existing_pythonpath
        else str(root) + os.pathsep + existing_pythonpath
    )
    completed = subprocess.run(command, cwd=str(root), env=environment, check=False)
    recovered = 0
    if completed.returncode != 0:
        recovered = _recover_finished_from_pool(upstream_output, root)
        if not recovered:
            return completed.returncode

    # Pocket2Mol names its output ``0.sdf`` .. ``N.sdf`` without zero padding,
    # so a plain lexicographic sort orders 10 before 2.  Natural order is also
    # completion order, which is what makes the size correlation visible.
    sources = sorted(upstream_output.glob("**/SDF/*.sdf"), key=_natural_key)
    if not sources:
        print("Pocket2Mol produced no SDF molecules", file=sys.stderr)
        return 2

    # Every finished molecule is kept, in completion order, for provenance.
    all_output = output / "sdf_all"
    all_output.mkdir(exist_ok=True)
    for index, source in enumerate(sources, start=1):
        shutil.copy2(source, all_output / f"mol_{index:06d}.sdf")

    counts = [_heavy_atom_count(source) for source in sources]
    chosen = _select_indices(counts, args.num_samples, args.select, args.target_heavy_atoms)

    sdf_output = output / "sdf"
    sdf_output.mkdir(exist_ok=True)
    for index, source_index in enumerate(chosen, start=1):
        shutil.copy2(sources[source_index], sdf_output / f"candidate_{index:06d}.sdf")

    known = [count for count in counts if count is not None]
    selected = [counts[i] for i in chosen if counts[i] is not None]
    summary = {
        "requested": args.num_samples,
        "finished_total": len(sources),
        "selected": len(chosen),
        "selection_mode": args.select,
        "target_heavy_atoms": args.target_heavy_atoms,
        "max_steps": args.max_steps,
        "beam_size": args.beam_size,
        "upstream_num_samples": config["sample"]["num_samples"],
        "upstream_return_code": completed.returncode,
        "recovered_from_checkpoint": recovered or None,
        "heavy_atoms_all": {
            "min": min(known) if known else None,
            "median": sorted(known)[len(known) // 2] if known else None,
            "max": max(known) if known else None,
        },
        "heavy_atoms_selected": {
            "min": min(selected) if selected else None,
            "median": sorted(selected)[len(selected) // 2] if selected else None,
            "max": max(selected) if selected else None,
        },
        "completion_order_heavy_atoms": counts,
        "selected_completion_indices": chosen,
        "note": (
            "Completion order correlates with molecule size because each sampling "
            "step adds one atom. --select controls which part of that distribution "
            "is delivered; sdf_all/ holds every finished molecule."
        ),
    }
    (output / "selection.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Pocket2Mol: {len(sources)} finished (heavy atoms "
        f"{summary['heavy_atoms_all']['min']}-{summary['heavy_atoms_all']['max']}), "
        f"delivered {len(chosen)} via '{args.select}' "
        f"(heavy atoms {summary['heavy_atoms_selected']['min']}-"
        f"{summary['heavy_atoms_selected']['max']})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
