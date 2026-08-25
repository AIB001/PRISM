"""Structure-based ligand generation tools.

The generation models run in their own isolated Conda, Docker, or Slurm
environments; this module only exposes PRISM's orchestration and quality
control over them, so every tool needs a generation YAML that points at a
configured backend.
"""

import json
import traceback

from ._common import _StdoutToStderr, _ensure_prism_importable, logger


def register(mcp):
    @mcp.tool()
    def list_generation_models() -> str:
        """List the structure-based ligand generation models PRISM can orchestrate.

        Reports each model's generation mode, so an agent can tell which ones
        need a conditioning ligand:

        - ``direct`` (Pocket2Mol, TargetDiff, PocketXMol): generate from pocket
          geometry alone. A ligand may still be supplied as the pocket in order
          to locate the site.
        - ``reference_guided`` (MolCRAFT, FLOWR, DiffSBDD): a reference ligand
          is a required conditioning input.

        No model runs here; this is a capability query and needs no config.

        Returns:
            JSON with one entry per model ID.
        """
        try:
            _ensure_prism_importable()
            with _StdoutToStderr():
                from prism.generation.registry import model_capabilities, model_ids

                capabilities = model_capabilities()
            return json.dumps(
                {"models": {model: capabilities[model] for model in model_ids()}},
                indent=2,
            )
        except Exception as exc:
            logger.error("list_generation_models failed: %s", exc)
            return json.dumps({"error": str(exc), "traceback": traceback.format_exc()})

    @mcp.tool()
    def generate_ligands(
        protein_path: str,
        pocket_path: str,
        generation_config: str,
        models: str = "diffsbdd",
        output_dir: str = "generated_ligand",
        num_samples: int = 100,
        seed: int = 0,
        device: str = "auto",
        pocket_kind: str = "auto",
        reference_ligand: str = "",
        qc_level: str = "full",
        add_hydrogens: bool = True,
        dry_run: bool = False,
        resume: bool = False,
        overwrite: bool = False,
    ) -> str:
        """Generate candidate 3D ligands for a binding pocket, with quality control.

        Requires a generation YAML (``generation_config``) that enables at least
        one model and gives its backend, environment, pinned checkpoint, and
        SHA256. Start from the examples in ``prism/configs/``. Without it no
        model can run, because the model environments live outside PRISM.

        Pocket input is interpreted by scientific meaning: PDB as a cropped
        pocket structure, SDF/MOL/MOL2 as a reference ligand, TXT as a residue
        list, YAML as an explicit center/residue/reference specification. Set
        ``pocket_kind`` to make that explicit.

        Every candidate is quality-checked. Molecules that fail chemically are
        quarantined with a written repair proposal rather than deleted or
        silently repaired: a bond order that sanitizes is not evidence it is the
        one the model intended. Warnings (geometry, pose, plausibility) are
        reported but do not withhold a molecule.

        Use ``dry_run=True`` first to validate inputs and capability contracts
        without starting a model.

        Args:
            protein_path: Absolute path to the receptor PDB/CIF.
            pocket_path: Absolute path to the pocket definition.
            generation_config: Absolute path to the generation YAML.
            models: Comma-separated model IDs, or "all".
            output_dir: Directory for candidates, reports, and the manifest.
            num_samples: Candidates requested per model.
            seed: Random seed passed to each wrapper.
            device: "auto", "cpu", "cuda", or "cuda:N".
            pocket_kind: auto, structure, reference_ligand, center, or residues.
            reference_ligand: Conditioning ligand for reference-guided models.
            qc_level: off, basic, standard, or full (full adds receptor pose checks).
            add_hydrogens: Also write candidates_h.sdf beside the raw output.
            dry_run: Validate and stage only.
            resume: Reuse matching successful model runs.
            overwrite: Replace the selected models' run directories.

        Returns:
            JSON with per-model status, candidate counts, quality summary,
            attrition, and the paths of the produced files.
        """
        try:
            _ensure_prism_importable()
            with _StdoutToStderr():
                from pathlib import Path

                from prism.generation.cli import main as _unused  # noqa: F401
                from prism.generation.registry import expand_models
                from prism.generation.runner import GenerationRunner, load_generation_config
                from prism.generation.types import GenerationRequest

                selected = tuple(
                    expand_models(
                        value.strip().lower()
                        for value in models.split(",")
                        if value.strip()
                    )
                )
                config = load_generation_config(Path(generation_config).expanduser())
                request = GenerationRequest(
                    protein=Path(protein_path).expanduser().resolve(),
                    pocket=Path(pocket_path).expanduser().resolve(),
                    output=Path(output_dir).expanduser().resolve(),
                    models=selected,
                    num_samples=num_samples,
                    seed=seed,
                    device=device,
                    pocket_kind=pocket_kind,
                    reference_ligand=(
                        Path(reference_ligand).expanduser().resolve()
                        if reference_ligand
                        else None
                    ),
                    qc_level=qc_level,
                    add_hydrogens=add_hydrogens,
                    dry_run=dry_run,
                    resume=resume,
                    overwrite=overwrite,
                )
                manifest = GenerationRunner(request, config).run()
                output = request.output

            return json.dumps(
                {
                    "status": manifest["status"],
                    "candidate_count": manifest["candidate_count"],
                    "models": {
                        model: {
                            "status": result["status"],
                            "requested": result["requested"],
                            "generated": result["generated"],
                            "quarantined": len(result.get("quarantined", [])),
                            "attrition": result.get("attrition", {}),
                            "error_code": result.get("error_code"),
                            "error_message": result.get("error_message"),
                        }
                        for model, result in manifest["models"].items()
                    },
                    "quality": manifest.get("quality", {}),
                    "outputs": {
                        "manifest": str(output / "manifest.json"),
                        "candidates_sdf": str(output / "candidates.sdf"),
                        "candidates_with_hydrogens_sdf": str(output / "candidates_h.sdf"),
                        "summary_csv": str(output / "summary.csv"),
                        "qc_report_csv": str(output / "qc_report.csv"),
                        "qc_issues_jsonl": str(output / "qc_issues.jsonl"),
                        "quarantine_note": (
                            "Failing molecules and their repair proposals are under "
                            "models/<model>/quarantine/; PRISM never applies one."
                        ),
                    },
                },
                indent=2,
            )
        except Exception as exc:
            logger.error("generate_ligands failed: %s", exc)
            return json.dumps(
                {
                    "error": str(exc),
                    "error_code": getattr(exc, "code", None),
                    "traceback": traceback.format_exc(),
                }
            )

    @mcp.tool()
    def prepare_generated_for_md(
        run_dir: str,
        output_dir: str = "md_inputs",
        models: str = "",
        only_pass: bool = False,
        limit: int = 0,
        overwrite: bool = False,
    ) -> str:
        """Export MD-ready ligands from a finished generation run.

        Generative models emit heavy atoms only. Handing that output straight to
        ``build_system`` does not fail: across 82 measured gaff2 builds it
        produced a topology with **zero hydrogens in 41 of 41 cases** -- correct
        file layout, wrong chemistry, no warning at any stage. Hydrogenated
        input gave the right hydrogen count in 35 of 41.

        This tool selects the candidates a run cleared, verifies each
        hydrogenated file actually carries hydrogens on disk, copies them into
        ``output_dir/ligands/``, and returns the exact ``prism`` command that
        builds them. It does not run the build itself.

        Molecules that quality control withheld are never exported. Every
        exclusion is written to ``skipped.csv`` with its reason.

        Args:
            run_dir: Output directory of a finished generate_ligands run.
            output_dir: Directory for the exported ligands and manifest.
            models: Comma-separated model IDs to export; empty means all.
            only_pass: Export only candidates graded PASS. Warnings are
                advisory and build successfully in most measured cases, so
                they are included by default.
            limit: Export at most this many, best QC grade first. 0 means all.
            overwrite: Replace an existing export directory.

        Returns:
            JSON with the exported ligand paths, QC grade counts, skip reasons,
            and the build command.
        """
        try:
            _ensure_prism_importable()
            with _StdoutToStderr():
                from pathlib import Path

                from prism.generation.handoff import export_md_inputs

                selected_models = [
                    value.strip().lower()
                    for value in models.split(",")
                    if value.strip()
                ]
                manifest = export_md_inputs(
                    Path(run_dir).expanduser(),
                    Path(output_dir).expanduser(),
                    models=selected_models or None,
                    only_pass=only_pass,
                    limit=limit if limit > 0 else None,
                    overwrite=overwrite,
                )
            return json.dumps(manifest, indent=2)
        except Exception as exc:
            logger.error("prepare_generated_for_md failed: %s", exc)
            return json.dumps(
                {
                    "error": str(exc),
                    "error_code": getattr(exc, "code", None),
                    "traceback": traceback.format_exc(),
                }
            )

    @mcp.tool()
    def check_generated_ligands(
        sdf_path: str,
        protein_path: str = "",
        reference_ligand: str = "",
        qc_level: str = "full",
        max_records: int = 500,
    ) -> str:
        """Quality-check an existing SDF of generated ligands without regenerating.

        Runs the same checks the generation pipeline applies: structural
        integrity, sanitization with a decoded failure reason, hydrogen
        addability, geometry (bond lengths, angles, aromatic planarity),
        chemical plausibility (property panel, PAINS/BRENK alerts, reactive
        groups, stereochemistry), and — when a receptor is supplied — steric
        agreement with the protein.

        Thresholds are calibrated so that experimentally determined ligands pass
        with no flags. Nothing is modified on disk.

        Args:
            sdf_path: Absolute path to a single- or multi-record SDF.
            protein_path: Receptor PDB; enables pose checks. Optional.
            reference_ligand: SDF whose centroid defines the expected pocket
                position, for the pocket-offset check. Optional.
            qc_level: off, basic, standard, or full.
            max_records: Stop after this many records.

        Returns:
            JSON with per-status counts, the flag histogram, and the worst cases.
        """
        try:
            _ensure_prism_importable()
            with _StdoutToStderr():
                from collections import Counter
                from pathlib import Path

                from prism.generation.postprocess import _sdf_records
                from prism.generation.quality import (
                    QualityContext,
                    assess,
                    load_protein_coordinates,
                    rdkit_available,
                )

                if not rdkit_available():
                    return json.dumps(
                        {"error": "RDKit is required for quality control and is not installed"}
                    )

                protein = (
                    load_protein_coordinates(Path(protein_path).expanduser())
                    if protein_path
                    else None
                )
                center = None
                if reference_ligand:
                    from rdkit import Chem

                    molecule = Chem.MolFromMolFile(
                        str(Path(reference_ligand).expanduser()), sanitize=False
                    )
                    if molecule is not None and molecule.GetNumConformers():
                        conformer = molecule.GetConformer()
                        points = [
                            (
                                conformer.GetAtomPosition(i).x,
                                conformer.GetAtomPosition(i).y,
                                conformer.GetAtomPosition(i).z,
                            )
                            for i in range(molecule.GetNumAtoms())
                            if molecule.GetAtomWithIdx(i).GetAtomicNum() > 1
                        ]
                        if points:
                            center = tuple(
                                sum(axis) / len(points) for axis in zip(*points)
                            )
                context = QualityContext(
                    level=qc_level, protein_coordinates=protein, pocket_center=center
                )

                statuses = Counter()
                flags = Counter()
                worst = []
                total = 0
                text = Path(sdf_path).expanduser().read_text(
                    encoding="utf-8", errors="replace"
                )
                for index, record in enumerate(_sdf_records(text), start=1):
                    if total >= max_records:
                        break
                    total += 1
                    report = assess(record, context)
                    statuses[report.status] += 1
                    for flag in report.flags:
                        flags[flag.split(":")[0]] += 1
                    if report.status.startswith("FAIL") or report.flags:
                        worst.append(
                            {
                                "record": index,
                                "status": report.status,
                                "flags": report.flags,
                                "smiles": report.canonical_smiles,
                                "repair": (
                                    report.repair.to_dict() if report.repair else None
                                ),
                            }
                        )

            # Severity first: the listing is truncated, and a chemistry failure
            # pushed past the cut would read as "nothing failed" even though
            # status_counts says otherwise.
            from prism.generation.quality import STATUS_ORDER

            def severity(entry):
                status = entry["status"]
                rank = STATUS_ORDER.index(status) if status in STATUS_ORDER else -1
                return (-rank, -len(entry["flags"]))

            worst.sort(key=severity)
            shown = 25

            return json.dumps(
                {
                    "records_checked": total,
                    "qc_level": context.level,
                    "pose_checks_enabled": context.wants_pose,
                    "status_counts": dict(statuses),
                    "flag_counts": dict(flags.most_common()),
                    "problem_record_count": len(worst),
                    "problem_records_shown": min(shown, len(worst)),
                    "problem_records": worst[:shown],
                    "note": (
                        "Warnings are advisory. A FAIL_CHEMISTRY record carries a "
                        "repair proposal that PRISM does not apply."
                    ),
                },
                indent=2,
            )
        except Exception as exc:
            logger.error("check_generated_ligands failed: %s", exc)
            return json.dumps({"error": str(exc), "traceback": traceback.format_exc()})
