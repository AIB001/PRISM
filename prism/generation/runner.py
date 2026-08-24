"""Generation workflow orchestration and reproducibility metadata."""

import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from importlib import resources
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

import yaml

from .adapters.base import ModelAdapter
from .errors import CapabilityError, ConfigurationError, GenerationError, InputError
from .pocket import normalize_pocket
from .postprocess import QC_REPORT_FIELDS, collect_candidates, write_qc_report
from .quality import QualityContext, load_protein_coordinates, rdkit_available
from .registry import build_adapters
from .types import (
    CandidateRecord,
    GenerationRequest,
    ModelRunResult,
    PocketSpec,
    QuarantineRecord,
)


def load_generation_config(path: Optional[Path] = None) -> Dict[str, Any]:
    try:
        if path is None:
            with resources.open_text("prism.configs", "generation.yaml", encoding="utf-8") as handle:
                value = yaml.safe_load(handle)
        else:
            value = yaml.safe_load(path.expanduser().read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        location = path if path is not None else "bundled generation.yaml"
        raise ConfigurationError(f"Unable to read generation config '{location}': {exc}") from exc
    if not isinstance(value, dict):
        raise ConfigurationError("Generation config must contain a top-level YAML mapping")
    return value


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _request_hash(request: GenerationRequest, pocket: PocketSpec, config: Mapping[str, Any]) -> str:
    payload = {
        "request": request.to_dict(),
        "pocket": pocket.to_dict(),
        "protein_sha256": _file_hash(request.protein),
        "pocket_sha256": _file_hash(request.pocket),
        "reference_ligand_sha256": (
            _file_hash(pocket.reference_ligand) if pocket.reference_ligand is not None else None
        ),
        "config": config,
    }
    payload["request"].pop("resume", None)
    payload["request"].pop("overwrite", None)
    payload["request"].pop("dry_run", None)
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()


class GenerationRunner:
    def __init__(
        self,
        request: GenerationRequest,
        config: Mapping[str, Any],
        adapters: Optional[Mapping[str, ModelAdapter]] = None,
    ):
        self.request = request
        self.config = dict(config)
        self.adapters = dict(adapters) if adapters is not None else build_adapters(config, request.models)

    def _validate_request(self) -> None:
        if not self.request.protein.is_file():
            raise InputError(f"Protein file does not exist: {self.request.protein}")
        if self.request.protein.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
            raise InputError("Protein input must be PDB, CIF, or mmCIF")
        if self.request.num_samples <= 0:
            raise InputError("--num-samples must be a positive integer")
        if not self.request.models:
            raise InputError("At least one generation model is required")
        missing = sorted(set(self.request.models) - set(self.adapters))
        if missing:
            raise ConfigurationError("No adapter supplied for model(s): " + ", ".join(missing))

    def _validate_reference_ligand_contract(self, pocket: PocketSpec) -> None:
        """Validate the shared conditioning input before creating output files."""
        required = sorted(
            model
            for model in self.request.models
            if self.adapters[model].requires_reference_ligand
        )
        if required and pocket.reference_ligand is None:
            raise CapabilityError(
                "Reference ligand required by selected model(s): "
                + ", ".join(required)
                + ". Provide --reference-ligand LIGAND.sdf, or use the ligand as "
                "--pocket with --pocket-kind reference_ligand.",
                code="REFERENCE_LIGAND_REQUIRED",
            )
        for model in required:
            self.adapters[model].check_reference_ligand(pocket)

        # A ligand used as the pocket definition is meaningful for direct
        # models because it supplies the center/bounding geometry. A separate
        # conditioning ligand is not, unless at least one selected model is
        # reference-guided.
        if (
            pocket.reference_ligand is not None
            and pocket.kind != "reference_ligand"
            and not required
        ):
            direct = ", ".join(sorted(self.request.models))
            raise CapabilityError(
                "--reference-ligand is not needed by the selected direct generation "
                f"model(s): {direct}. Remove it, or pass that ligand as --pocket "
                "when it is intended to define the pocket geometry.",
                code="REFERENCE_LIGAND_NOT_NEEDED",
            )

    def _prepare_output(self, request_hash: str) -> None:
        output = self.request.output
        if output.exists() and not output.is_dir():
            raise InputError(f"Output path is not a directory: {output}")
        manifest_path = output / "manifest.json"
        if manifest_path.exists() and not (self.request.resume or self.request.overwrite):
            raise InputError(
                f"Generation output already exists: {output}. Use --resume or --overwrite."
            )
        if self.request.resume and manifest_path.exists():
            try:
                previous = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise InputError(f"Unable to read existing manifest for resume: {exc}") from exc
            if previous.get("request_hash") != request_hash:
                raise InputError("Existing output does not match this request; resume was refused")
        output.mkdir(parents=True, exist_ok=True)

    def _stage_inputs(self, pocket: PocketSpec) -> Dict[str, Path]:
        inputs = self.request.output / "inputs"
        inputs.mkdir(parents=True, exist_ok=True)
        staged = {
            "protein": inputs / f"protein{self.request.protein.suffix.lower()}",
            "pocket": inputs / f"pocket{self.request.pocket.suffix.lower()}",
        }
        shutil.copy2(self.request.protein, staged["protein"])
        shutil.copy2(self.request.pocket, staged["pocket"])
        if pocket.reference_ligand is not None:
            staged["reference_ligand"] = inputs / f"reference{pocket.reference_ligand.suffix.lower()}"
            if pocket.reference_ligand.resolve() != staged["pocket"].resolve():
                shutil.copy2(pocket.reference_ligand, staged["reference_ligand"])
        _write_json(inputs / "pocket.json", pocket.to_dict())
        return staged

    @staticmethod
    def _status_result(model: str, run_dir: Path) -> Optional[ModelRunResult]:
        status_path = run_dir / "status.json"
        if not status_path.is_file():
            return None
        try:
            data = json.loads(status_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        result_data = dict(data.get("result", {}))
        if result_data.get("model") != model:
            return None
        try:
            result_data["candidates"] = [
                CandidateRecord(**item) for item in result_data.get("candidates", [])
            ]
            result_data["quarantined"] = [
                QuarantineRecord(**item) for item in result_data.get("quarantined", [])
            ]
            return ModelRunResult(**result_data)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _quality_summary(
        context: QualityContext,
        candidates,
        quarantined,
        results: Mapping[str, ModelRunResult],
    ) -> Dict[str, Any]:
        """Aggregate QC counters, including where samples were lost."""
        status_counts: Dict[str, int] = {}
        for candidate in candidates:
            status_counts[candidate.audit_status] = status_counts.get(candidate.audit_status, 0) + 1
        for item in quarantined:
            status_counts[item.audit_status] = status_counts.get(item.audit_status, 0) + 1
        return {
            "level": context.level,
            "rdkit_available": rdkit_available(),
            "pose_checks_enabled": context.wants_pose,
            "thresholds": dict(context.thresholds),
            "status_counts": status_counts,
            "quarantined_count": len(quarantined),
            "repair_proposals": {
                "unique_minimal": sum(
                    1 for item in quarantined if item.repair_confidence == "unique_minimal"
                ),
                "ambiguous": sum(
                    1 for item in quarantined if item.repair_confidence == "ambiguous"
                ),
                "none": sum(1 for item in quarantined if item.repair_confidence == "none"),
                "applied": 0,
                "policy": "detect-only; PRISM never applies a speculative repair",
            },
            "attrition": {model: dict(result.attrition) for model, result in results.items()},
        }

    def _quality_context(self, pocket: PocketSpec, staged_paths: Mapping[str, Path]) -> QualityContext:
        """Build the QC settings shared by every model in this request.

        Pose checks need the receptor.  When the protein cannot be read as
        coordinates the level degrades to ``standard`` rather than failing the
        run, because a missing pose check must not cost the user their molecules.
        """
        level = self.request.qc_level
        if level == "off" or not rdkit_available():
            return QualityContext(level="off")
        protein_coordinates = None
        if level == "full":
            protein = staged_paths.get("protein", self.request.protein)
            if Path(protein).suffix.lower() == ".pdb":
                protein_coordinates = load_protein_coordinates(Path(protein))
            if protein_coordinates is None:
                level = "standard"
        return QualityContext(
            level=level,
            protein_coordinates=protein_coordinates,
            pocket_center=pocket.center,
        )

    def _resume_result(self, model: str, run_dir: Path) -> Optional[ModelRunResult]:
        if not self.request.resume:
            return None
        result = self._status_result(model, run_dir)
        if result is None or result.status not in {"SUCCEEDED", "PARTIAL"}:
            return None
        if not (run_dir / "candidates.sdf").is_file():
            return None
        return result

    def run(self) -> Dict[str, Any]:
        self._validate_request()
        pocket = normalize_pocket(
            self.request.pocket,
            kind=self.request.pocket_kind,
            reference_ligand=self.request.reference_ligand,
        )
        self._validate_reference_ligand_contract(pocket)
        request_hash = _request_hash(self.request, pocket, self.config)
        self._prepare_output(request_hash)
        staged_paths = self._stage_inputs(pocket)
        output = self.request.output
        (output / "request.yaml").write_text(
            yaml.safe_dump(self.request.to_dict(), sort_keys=True), encoding="utf-8"
        )
        (output / "generation_config.yaml").write_text(
            yaml.safe_dump(self.config, sort_keys=True), encoding="utf-8"
        )

        started_at = _utc_now()
        quality_context = self._quality_context(pocket, staged_paths)
        results: Dict[str, ModelRunResult] = {}
        validation_failures = []
        for model in self.request.models:
            run_dir = output / "models" / model
            resumed = self._resume_result(model, run_dir)
            if resumed is not None:
                results[model] = resumed
                continue
            if run_dir.exists() and self.request.overwrite:
                shutil.rmtree(run_dir)
            run_dir.mkdir(parents=True, exist_ok=True)
            history = [{"state": "PENDING", "at": _utc_now()}]
            if self.request.dry_run:
                try:
                    self.adapters[model].check_capabilities(pocket)
                    result = ModelRunResult(model=model, status="PLANNED", requested=self.request.num_samples)
                except GenerationError as exc:
                    result = ModelRunResult(
                        model=model,
                        status="FAILED",
                        requested=self.request.num_samples,
                        error_code=exc.code,
                        error_message=str(exc),
                    )
                history.append({"state": result.status, "at": _utc_now()})
                _write_json(run_dir / "status.json", {"status": result.status, "history": history, "result": result.to_dict()})
                results[model] = result
                continue
            try:
                history.append({"state": "PREPARING", "at": _utc_now()})
                adapter = self.adapters[model]
                adapter.check_capabilities(pocket)
                history.append({"state": "RUNNING", "at": _utc_now()})
                execution = adapter.run(self.request, pocket, output, run_dir, staged_paths)
                history.append({"state": "COLLECTING", "at": _utc_now()})
                candidates, quarantined, candidate_failures, attrition = collect_candidates(
                    model,
                    run_dir,
                    adapter.output_patterns(),
                    self.request.num_samples,
                    context=quality_context,
                    add_hydrogens=self.request.add_hydrogens,
                )
                validation_failures.extend(
                    {"model": model, "code": "INVALID_CANDIDATE", "message": message}
                    for message in candidate_failures
                )
                status = "SUCCEEDED" if len(candidates) >= self.request.num_samples else "PARTIAL"
                history.append({"state": "VALIDATING", "at": _utc_now()})
                result = ModelRunResult(
                    model=model,
                    status=status,
                    requested=self.request.num_samples,
                    generated=len(candidates),
                    candidates=candidates,
                    command=execution.command,
                    return_code=execution.return_code,
                    provenance=adapter.provenance(execution),
                    quarantined=quarantined,
                    attrition=attrition,
                )
            except GenerationError as exc:
                result = ModelRunResult(
                    model=model,
                    status="FAILED",
                    requested=self.request.num_samples,
                    error_code=exc.code,
                    error_message=str(exc),
                )
            except Exception as exc:  # keep independent models running
                result = ModelRunResult(
                    model=model,
                    status="FAILED",
                    requested=self.request.num_samples,
                    error_code="UNEXPECTED_ERROR",
                    error_message=str(exc),
                )
            history.append({"state": result.status, "at": _utc_now()})
            _write_json(run_dir / "status.json", {"status": result.status, "history": history, "result": result.to_dict()})
            results[model] = result

        # ``--overwrite`` is scoped to the explicitly selected model
        # directories.  Preserve valid status records for every unselected
        # model so the aggregate manifest, SDF, and CSV remain a complete view
        # of the output directory after a single-model retry.
        if self.request.overwrite:
            models_dir = output / "models"
            for run_dir in sorted(models_dir.iterdir()) if models_dir.is_dir() else []:
                model = run_dir.name
                if model in results or not run_dir.is_dir():
                    continue
                preserved = self._status_result(model, run_dir)
                if preserved is not None:
                    results[model] = preserved

        candidates = [candidate for result in results.values() for candidate in result.candidates]
        quarantined = [item for result in results.values() for item in result.quarantined]
        (output / "candidates.sdf").write_text(
            "".join(Path(candidate.path).read_text(encoding="utf-8") for candidate in candidates),
            encoding="utf-8",
        )
        hydrogenated = [
            candidate.hydrogenated_path
            for candidate in candidates
            if candidate.hydrogenated_path and Path(candidate.hydrogenated_path).is_file()
        ]
        if hydrogenated:
            # Written alongside, never in place of, candidates.sdf: adding
            # hydrogens preserves every heavy-atom position but still commits to
            # one valence-based protonation state.
            (output / "candidates_h.sdf").write_text(
                "".join(Path(path).read_text(encoding="utf-8") for path in hydrogenated),
                encoding="utf-8",
            )
        with (output / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "candidate_id",
                    "model",
                    "path",
                    "valid_3d",
                    "audit_status",
                    "sanitize_status",
                    "quality_flags",
                    "canonical_smiles",
                    "hydrogenated_path",
                    "native_score_name",
                    "native_score",
                ],
            )
            writer.writeheader()
            for candidate in candidates:
                data = candidate.to_dict()
                data["quality_flags"] = ";".join(candidate.quality_flags)
                writer.writerow({key: data[key] for key in writer.fieldnames})

        # One aggregate QC table across every model in the output directory.
        aggregate_rows = []
        for model in sorted(results):
            model_report = output / "models" / model / "qc_report.csv"
            if not model_report.is_file():
                continue
            with model_report.open(encoding="utf-8", newline="") as handle:
                aggregate_rows.extend(dict(row) for row in csv.DictReader(handle))
        write_qc_report(output, aggregate_rows)

        failures = [
            {
                "model": result.model,
                "code": result.error_code,
                "message": result.error_message,
            }
            for result in results.values()
            if result.status == "FAILED"
        ] + validation_failures
        (output / "failures.jsonl").write_text(
            "".join(json.dumps(item, sort_keys=True) + "\n" for item in failures), encoding="utf-8"
        )
        statuses = {result.status for result in results.values()}
        if statuses == {"PLANNED"}:
            overall_status = "PLANNED"
        elif statuses <= {"SUCCEEDED"}:
            overall_status = "SUCCEEDED"
        elif statuses <= {"SUCCEEDED", "PARTIAL"}:
            overall_status = "PARTIAL" if "PARTIAL" in statuses else "SUCCEEDED"
        elif statuses == {"FAILED"}:
            overall_status = "FAILED"
        else:
            overall_status = "PARTIAL"
        manifest = {
            "schema_version": 2,
            "status": overall_status,
            "request_hash": request_hash,
            "started_at": started_at,
            "finished_at": _utc_now(),
            "request": self.request.to_dict(),
            "pocket": pocket.to_dict(),
            "models": {model: result.to_dict() for model, result in results.items()},
            "candidate_count": len(candidates),
            "quality": self._quality_summary(quality_context, candidates, quarantined, results),
        }
        _write_json(output / "manifest.json", manifest)
        return manifest
