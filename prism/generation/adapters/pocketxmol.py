import math
from pathlib import Path
from typing import Any, Mapping

from ..errors import ConfigurationError, InputError
from ..types import GenerationRequest, PocketSpec
from .base import ModelAdapter


class PocketXMolAdapter(ModelAdapter):
    model_id = "pocketxmol"
    accepted_pocket_kinds = {"center", "reference_ligand", "structure"}

    def prepare_inputs(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        _run_dir: Path,
        staged_paths: Mapping[str, Path],
    ) -> Mapping[str, Path]:
        if staged_paths["protein"].suffix.lower() != ".pdb":
            raise InputError("PocketXMol custom SBDD inference requires a PDB protein")
        if pocket.center is None:
            raise InputError("PocketXMol requires a pocket center")
        return staged_paths

    def command_context(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        _prepared_paths: Mapping[str, Path],
    ) -> Mapping[str, Any]:
        if pocket.radius is not None:
            radius = pocket.radius
        else:
            key = (
                "default_reference_radius"
                if pocket.kind == "reference_ligand"
                else "default_center_radius"
            )
            fallback = 10.0 if pocket.kind == "reference_ligand" else 15.0
            configured = self.config.get(key, fallback)
            if isinstance(configured, bool):
                raise ConfigurationError(f"PocketXMol '{key}' must be positive")
            try:
                radius = float(configured)
            except (TypeError, ValueError) as exc:
                raise ConfigurationError(f"PocketXMol '{key}' must be positive") from exc
        if not math.isfinite(radius) or radius <= 0:
            raise ConfigurationError("PocketXMol pocket radius must be positive")
        return {"pocketxmol_radius": radius}
