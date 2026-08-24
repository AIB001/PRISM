"""Base protocol shared by all external model adapters."""

from pathlib import Path
from typing import Any, Dict, Mapping, Sequence, Set

from ..errors import CapabilityError, ConfigurationError
from ..execution import ExecutionResult, execute
from ..types import GenerationRequest, PocketSpec


class ModelAdapter:
    model_id = ""
    # ``direct`` models generate from pocket geometry alone. ``reference_guided``
    # models use a ligand as a required conditioning input, not merely as a
    # convenient way to locate the pocket.
    generation_mode = "direct"
    accepted_pocket_kinds: Set[str] = set()
    requires_reference_ligand = False
    accepted_reference_ligand_suffixes: Set[str] = set()

    def __init__(self, config: Mapping[str, Any]):
        self.config = dict(config)

    def check_capabilities(self, pocket: PocketSpec) -> None:
        if pocket.kind not in self.accepted_pocket_kinds:
            supported = ", ".join(sorted(self.accepted_pocket_kinds))
            raise CapabilityError(
                f"{self.model_id} does not accept pocket kind '{pocket.kind}'; supported: {supported}"
            )
        self.check_reference_ligand(pocket)

    def check_reference_ligand(self, pocket: PocketSpec) -> None:
        if self.requires_reference_ligand and pocket.reference_ligand is None:
            raise CapabilityError(
                f"{self.model_id} requires a reference ligand",
                code="REFERENCE_LIGAND_REQUIRED",
            )
        if (
            self.requires_reference_ligand
            and pocket.reference_ligand is not None
            and self.accepted_reference_ligand_suffixes
            and pocket.reference_ligand.suffix.lower()
            not in self.accepted_reference_ligand_suffixes
        ):
            supported = ", ".join(sorted(self.accepted_reference_ligand_suffixes))
            raise CapabilityError(
                f"{self.model_id} reference ligand must use one of: {supported}",
                code="REFERENCE_LIGAND_FORMAT_UNSUPPORTED",
            )

    def output_patterns(self) -> Sequence[str]:
        patterns = self.config.get("output_patterns", ["raw/**/*.sdf"])
        if not isinstance(patterns, list) or not patterns or not all(isinstance(item, str) for item in patterns):
            raise ConfigurationError("'output_patterns' must be a non-empty list of glob strings")
        return patterns

    def prepare_inputs(
        self,
        _request: GenerationRequest,
        _pocket: PocketSpec,
        _run_dir: Path,
        staged_paths: Mapping[str, Path],
    ) -> Mapping[str, Path]:
        return staged_paths

    def command_context(
        self,
        _request: GenerationRequest,
        _pocket: PocketSpec,
        _prepared_paths: Mapping[str, Path],
    ) -> Mapping[str, Any]:
        """Return model-specific values for configured command placeholders."""
        return {}

    def provenance(self, execution: ExecutionResult) -> Dict[str, Any]:
        return {
            "backend": self.config.get("backend"),
            "image": self.config.get("image"),
            "environment": self.config.get("environment"),
            "slurm": self.config.get("slurm"),
            "model_commit": self.config.get("model_commit"),
            "artifacts": execution.artifact_sha256,
        }

    def run(
        self,
        request: GenerationRequest,
        pocket: PocketSpec,
        output_root: Path,
        run_dir: Path,
        staged_paths: Mapping[str, Path],
    ) -> ExecutionResult:
        self.check_capabilities(pocket)
        raw_dir = run_dir / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)
        prepared_paths = dict(
            self.prepare_inputs(request, pocket, run_dir, staged_paths)
        )
        if not self.requires_reference_ligand and pocket.kind != "reference_ligand":
            # A conditioning ligand may be shared by another selected model,
            # but it must not silently change a direct model's request.
            prepared_paths.pop("reference_ligand", None)
        context: Dict[str, Any] = {
            "protein": prepared_paths["protein"],
            "pocket": prepared_paths["pocket"],
            "output": raw_dir,
            "num_samples": request.num_samples,
            "seed": request.seed,
            "device": request.device,
            "pocket_kind": pocket.kind,
            "center_x": pocket.center[0] if pocket.center else "",
            "center_y": pocket.center[1] if pocket.center else "",
            "center_z": pocket.center[2] if pocket.center else "",
            "radius": pocket.radius if pocket.radius is not None else "",
            "bbox_x": pocket.bbox_size[0] if pocket.bbox_size else "",
            "bbox_y": pocket.bbox_size[1] if pocket.bbox_size else "",
            "bbox_z": pocket.bbox_size[2] if pocket.bbox_size else "",
            "reference_ligand": prepared_paths.get("reference_ligand", ""),
        }
        context.update({key: value for key, value in prepared_paths.items() if key not in context})
        extra_context = self.command_context(request, pocket, prepared_paths)
        overlap = sorted(set(context) & set(extra_context))
        if overlap:
            raise ConfigurationError(
                "Model-specific command placeholders conflict with built-ins: "
                + ", ".join(overlap)
            )
        context.update(extra_context)
        return execute(self.config, context, output_root, run_dir, request.device)
