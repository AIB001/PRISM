"""Command-line interface for ligand generation orchestration."""

import argparse
from pathlib import Path
from typing import Optional, Sequence

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
