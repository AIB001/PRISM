"""Command-line interface for ligand generation orchestration."""

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from . import weights
from .errors import GenerationError
from .handoff import export_md_inputs
from .registry import expand_models, model_capabilities, model_ids
from .runner import GenerationRunner, load_generation_config
from .types import GenerationRequest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="prism generate",
        description="Generate candidate 3D ligands using isolated protein-pocket model environments.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  prism generate --model targetdiff --protein protein.pdb --pocket pocket.pdb -o generated_ligand
  prism --generation true --model pocket2mol --protein protein.pdb --pocket reference.sdf -o generated_ligand
  prism generate --model flowr --protein protein.pdb --pocket pocket.pdb --reference-ligand reference.sdf -o generated_ligand
  prism generate --model all --protein protein.pdb --pocket pocket.yaml -o generated_ligand --dry-run
""",
    )
    parser.add_argument(
        "--model",
        action="append",
        default=[],
        metavar="MODEL",
        help="Model ID; repeat for multiple models or use 'all'",
    )
    parser.add_argument("--protein", type=Path, help="Protein structure (PDB/CIF/mmCIF)")
    parser.add_argument("--pocket", type=Path, help="Pocket input (PDB/SDF/MOL/MOL2/YAML/TXT)")
    parser.add_argument(
        "--pocket-kind",
        choices=["auto", "structure", "reference_ligand", "center", "residues"],
        default="auto",
        help="Scientific meaning of --pocket (default: infer from extension/schema)",
    )
    parser.add_argument(
        "--reference-ligand",
        type=Path,
        help=(
            "Conditioning ligand required by reference-guided models. Direct models "
            "only use a ligand when it is supplied as --pocket to locate the pocket"
        ),
    )
    parser.add_argument("--output", "-o", type=Path, default=Path("generated_ligand"))
    parser.add_argument("--num-samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--device",
        default="auto",
        help="Device passed to wrappers: auto, cpu, cuda, or cuda:N",
    )
    parser.add_argument("--generation-config", type=Path, help="Docker/Conda model backend YAML")
    parser.add_argument(
        "--qc",
        dest="qc_level",
        choices=["off", "basic", "standard", "full"],
        default="full",
        help=(
            "Quality-control depth. basic: integrity, sanitization, hydrogens. "
            "standard: adds geometry and chemical plausibility. full (default): "
            "adds pose checks against the receptor. Failing molecules are always "
            "quarantined with a repair proposal, never silently repaired"
        ),
    )
    parser.add_argument(
        "--no-hydrogens",
        dest="add_hydrogens",
        action="store_false",
        help=(
            "Skip writing candidates_h.sdf. Hydrogen addition preserves every "
            "heavy-atom position but commits to one valence-based protonation state"
        ),
    )
    parser.add_argument("--resume", action="store_true", help="Reuse successful matching model runs")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace only the selected models' existing run directories",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and stage the request without starting model environments",
    )
    parser.add_argument("--list-models", action="store_true", help="List registered model IDs and exit")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.list_models:
        capabilities = model_capabilities()
        for model in model_ids():
            mode = capabilities[model]["generation_mode"].replace("_", "-")
            if capabilities[model]["requires_reference_ligand"]:
                formats = ", ".join(
                    capabilities[model]["accepted_reference_ligand_suffixes"]
                )
                reference = f"required ({formats})"
            else:
                reference = "not required; ligand pocket supported"
            print(f"{model}\t{mode}\treference ligand: {reference}")
        return 0
    if not args.model:
        parser.error("at least one --model is required")
    if args.protein is None:
        parser.error("--protein is required")
    if args.pocket is None:
        parser.error("--pocket is required")
    try:
        models = tuple(expand_models(value.lower() for value in args.model))
        config = load_generation_config(args.generation_config)
        request = GenerationRequest(
            protein=args.protein.expanduser().resolve(),
            pocket=args.pocket.expanduser().resolve(),
            output=args.output.expanduser().resolve(),
            models=models,
            num_samples=args.num_samples,
            seed=args.seed,
            device=args.device,
            pocket_kind=args.pocket_kind,
            reference_ligand=args.reference_ligand.expanduser().resolve() if args.reference_ligand else None,
            resume=args.resume,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            qc_level=args.qc_level,
            add_hydrogens=args.add_hydrogens,
        )
        manifest = GenerationRunner(request, config).run()
    except GenerationError as exc:
        parser.error(f"{exc.code}: {exc}")
    print(
        f"Generation {manifest['status'].lower()}: "
        f"{manifest['candidate_count']} candidate(s) in {request.output}"
    )
    _print_quality_summary(manifest)
    return 1 if manifest["status"] == "FAILED" else 0


def build_prepare_md_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="prism prepare-md",
        description=(
            "Export MD-ready ligands from a finished generation run. Raw model "
            "output carries no hydrogens and builds into a hydrogen-free "
            "topology without any error, so only verified hydrogenated "
            "candidates are exported."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  prism prepare-md generated_ligand -o md_inputs
  prism prepare-md generated_ligand -o md_inputs --only-pass --limit 10
  prism prepare-md generated_ligand -o md_inputs --model flowr --model pocketxmol
""",
    )
    parser.add_argument(
        "run_dir",
        type=Path,
        help="Output directory of a finished 'prism generate' run",
    )
    parser.add_argument("--output", "-o", type=Path, default=Path("md_inputs"))
    parser.add_argument(
        "--model",
        action="append",
        default=[],
        metavar="MODEL",
        help="Export only this model's candidates; repeat for several",
    )
    parser.add_argument(
        "--only-pass",
        action="store_true",
        help=(
            "Export only candidates graded PASS. Warnings are advisory and "
            "build successfully in most measured cases, so they are included "
            "by default"
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Export at most this many candidates, best QC grade first",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace an existing export directory"
    )
    return parser


def prepare_md_main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_prepare_md_parser()
    args = parser.parse_args(argv)
    try:
        manifest = export_md_inputs(
            args.run_dir,
            args.output,
            models=args.model or None,
            only_pass=args.only_pass,
            limit=args.limit,
            overwrite=args.overwrite,
        )
    except GenerationError as exc:
        parser.error(f"{exc.code}: {exc}")

    print(
        f"Exported {manifest['exported_count']} MD-ready ligand(s) to "
        f"{args.output}"
    )
    if manifest["status_counts"]:
        print(
            "  QC grades: "
            + ", ".join(
                f"{status}={count}"
                for status, count in sorted(manifest["status_counts"].items())
            )
        )
    if manifest["skipped_count"]:
        print(
            f"  {manifest['skipped_count']} skipped: "
            + ", ".join(
                f"{reason}={count}"
                for reason, count in sorted(manifest["skipped_reasons"].items())
            )
            + f" -- see {args.output / 'skipped.csv'}"
        )
    if manifest["exported_count"]:
        print(f"\nBuild the systems with:\n  {manifest['build_command']}")
    return 0 if manifest["exported_count"] else 1


def build_weights_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="prism weights",
        description=(
            "Inspect, download and verify the generative models' checkpoints. "
            "The checkpoints are about a gigabyte in total and are not part of a "
            "PRISM install; they are fetched into $PRISM_MODELS_DIR on demand."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  prism weights list
  prism weights list --verify
  prism weights path --model pocketxmol
  prism weights download --model pocketxmol --model molcraft
  prism weights verify
""",
    )
    subparsers = parser.add_subparsers(dest="command")
    for name, help_text in (
        ("list", "Show every artifact with its size, licence and local status"),
        ("path", "Print the local path of each artifact"),
        ("download", "Fetch the artifacts PRISM is licensed to mirror"),
        ("verify", "Re-hash every present artifact and report corruption"),
    ):
        sub = subparsers.add_parser(name, help=help_text, description=help_text)
        sub.add_argument(
            "--model",
            action="append",
            default=[],
            metavar="MODEL",
            help="Restrict to this model; repeat for several (default: all)",
        )
        if name == "list":
            sub.add_argument(
                "--verify",
                action="store_true",
                help="Re-hash present artifacts instead of only checking their size",
            )
        if name == "download":
            sub.add_argument(
                "--force", action="store_true", help="Download again even if the artifact is already valid"
            )
    return parser


def _selected_models(parser: argparse.ArgumentParser, values: Sequence[str]) -> Optional[list]:
    if not values:
        return None
    try:
        return expand_models(value.lower() for value in values)
    except GenerationError as exc:
        parser.error(f"{exc.code}: {exc}")


def _print_weights_header() -> None:
    print(f"Weights directory: {weights.models_dir()}  (override with {weights.MODELS_DIR_ENV})")
    base = weights.mirror_base_url()
    print(f"Mirror: {base}" if base else f"Mirror: not configured (set {weights.MIRROR_ENV})")
    print()


def _download_progress(label: str, done: int, total: int) -> None:
    if not sys.stderr.isatty():
        return
    fraction = min(1.0, done / total) if total else 0.0
    sys.stderr.write(
        f"\r  {label}  {fraction * 100:5.1f}%  "
        f"{weights.format_size(done)} / {weights.format_size(total)}   "
    )
    sys.stderr.flush()


def weights_main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_weights_parser()
    args = parser.parse_args(argv)
    command = args.command or "list"
    models = _selected_models(parser, args.model) if hasattr(args, "model") else None

    try:
        if command in {"list", "verify"}:
            report = weights.artifact_report(models, verify=getattr(args, "verify", command == "verify"))
            _print_weights_header()
            specs = {spec.model: spec for spec in weights.model_specs(models)}
            print(f"{'MODEL':<12}{'ARTIFACT':<14}{'SIZE':>10}  {'STATUS':<16}{'LICENCE':<20}SOURCE")
            for artifact, status in report:
                spec = specs[artifact.model]
                source = "PRISM mirror" if artifact.mirror_asset else "upstream only"
                print(
                    f"{artifact.model:<12}{artifact.name:<14}"
                    f"{weights.format_size(artifact.size_bytes):>10}  {status:<16}"
                    f"{spec.weights_license:<20}{source}"
                )
            outstanding = [(artifact, status) for artifact, status in report if status != "present"]
            print()
            if outstanding:
                total = sum(artifact.size_bytes for artifact, _ in outstanding)
                print(
                    f"{len(outstanding)} of {len(report)} artifact(s) unavailable "
                    f"({weights.format_size(total)}). Run 'prism weights download' for the "
                    "mirrored ones; the rest are listed with upstream instructions when a "
                    "run needs them."
                )
                return 1 if command == "verify" else 0
            print(f"All {len(report)} artifact(s) present.")
            return 0

        if command == "path":
            for spec in weights.model_specs(models):
                for artifact in spec.artifacts.values():
                    print(f"{artifact.model}/{artifact.name}\t{artifact.path}")
            return 0

        if command == "download":
            specs = weights.model_specs(models)
            mirrorable = [spec for spec in specs if spec.redistributable]
            upstream_only = [spec for spec in specs if not spec.redistributable]
            if not mirrorable:
                print(weights.describe_missing(weights.missing_artifacts(models)) or "Nothing to download.")
                return 1
            _print_weights_header()
            for spec in mirrorable:
                for artifact in spec.artifacts.values():
                    status = weights.artifact_status(artifact, verify=True)
                    if status == "present" and not args.force:
                        print(f"  {artifact.model}/{artifact.name}  already present")
                        continue
                    target = weights.download_artifact(
                        artifact, force=args.force, progress=_download_progress
                    )
                    if sys.stderr.isatty():
                        sys.stderr.write("\r" + " " * 78 + "\r")
                        sys.stderr.flush()
                    print(f"  {artifact.model}/{artifact.name}  -> {target}")
                if spec.attribution:
                    print(f"      {spec.attribution}")
            for spec in upstream_only:
                print()
                print(
                    f"{spec.title} is not redistributed by PRISM (licence: {spec.weights_license}); "
                    "fetch it from upstream:"
                )
                if spec.upstream_url:
                    print(f"  {spec.upstream_url}")
                for artifact in spec.artifacts.values():
                    print(f"  -> {artifact.path}")
            return 0
    except GenerationError as exc:
        parser.error(f"{exc.code}: {exc}")

    parser.error(f"unknown weights command: {command}")


def _print_quality_summary(manifest: dict) -> None:
    """Surface QC outcomes, especially samples lost before PRISM saw them."""
    quality = manifest.get("quality") or {}
    if not quality or quality.get("level") == "off":
        return
    counts = quality.get("status_counts") or {}
    if counts:
        print(
            "  QC ("
            + quality["level"]
            + "): "
            + ", ".join(f"{status}={count}" for status, count in sorted(counts.items()))
        )
    quarantine_count = quality.get("quarantined_count", 0)
    if quarantine_count:
        proposals = quality.get("repair_proposals") or {}
        print(
            f"  {quarantine_count} quarantined; repair proposals: "
            f"{proposals.get('unique_minimal', 0)} unique, "
            f"{proposals.get('ambiguous', 0)} ambiguous, "
            f"{proposals.get('none', 0)} none. "
            "None applied -- review models/<model>/quarantine/"
        )
    for model, attrition in sorted((quality.get("attrition") or {}).items()):
        dropped = attrition.get("dropped_before_sdf", 0)
        if dropped:
            print(
                f"  {model}: {dropped} of {attrition.get('requested')} sample(s) were "
                "discarded by the model before any SDF was written; PRISM cannot "
                "inspect or repair those"
            )


if __name__ == "__main__":
    raise SystemExit(main())
