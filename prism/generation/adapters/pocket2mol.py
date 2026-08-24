import math
from pathlib import Path
from typing import Any, Mapping

from ..errors import ConfigurationError, InputError
from ..types import GenerationRequest, PocketSpec
from .base import ModelAdapter


def _reference_heavy_atoms(path: Path) -> int:
    """Heavy-atom count from a V2000 counts line, or 0 when unavailable.

    Pocket2Mol's molecule size tracks its sampling step count rather than the
    pocket, so the reference ligand's size is the natural target for
    ``--select target``. Only MOL/SDF expose the count without a chemistry
    toolkit; anything else leaves the placeholder empty and the config can fall
    back to a fixed value.
    """
    if path is None or path.suffix.lower() not in {".sdf", ".mol"}:
        return 0
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return 0
    if len(lines) < 4 or "V3000" in lines[3]:
        return 0
    try:
        return int(lines[3][:3])
    except ValueError:
        return 0


class Pocket2MolAdapter(ModelAdapter):
    model_id = "pocket2mol"
    accepted_pocket_kinds = {"center", "reference_ligand", "structure"}

    def prepare_inputs(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        _run_dir: Path,
        staged_paths: Mapping[str, Path],
    ) -> Mapping[str, Path]:
        if staged_paths["protein"].suffix.lower() != ".pdb":
            raise InputError("Pocket2Mol custom-pocket inference requires a PDB protein")
        if pocket.center is None:
            raise InputError("Pocket2Mol requires a pocket center")
        return staged_paths

    def command_context(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        _prepared_paths: Mapping[str, Path],
    ) -> Mapping[str, Any]:
        configured = self.config.get("default_bbox_size", 23.0)
        if isinstance(configured, bool):
            raise ConfigurationError("Pocket2Mol 'default_bbox_size' must be positive")
        try:
            default_bbox = float(configured)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError("Pocket2Mol 'default_bbox_size' must be positive") from exc
        if not math.isfinite(default_bbox) or default_bbox <= 0:
            raise ConfigurationError("Pocket2Mol 'default_bbox_size' must be positive")

        # The official entry point accepts one cubic side length. A requested
        # radius takes precedence; otherwise never crop tighter than the
        # official 23 Å default or the normalized pocket's largest extent.
        if pocket.radius is not None:
            bbox_size = 2.0 * pocket.radius
        elif pocket.bbox_size is not None:
            bbox_size = max(default_bbox, *pocket.bbox_size)
        else:
            bbox_size = default_bbox
        return {
            "bbox_size": bbox_size,
            "pocket2mol_target_heavy_atoms": _reference_heavy_atoms(pocket.reference_ligand),
        }
