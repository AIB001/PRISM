from pathlib import Path
from typing import Mapping

from ..protein import extract_pocket_pdb
from ..types import GenerationRequest, PocketSpec
from .base import ModelAdapter


class TargetDiffAdapter(ModelAdapter):
    model_id = "targetdiff"
    accepted_pocket_kinds = {"center", "reference_ligand", "residues", "structure"}

    def prepare_inputs(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        run_dir: Path,
        staged_paths: Mapping[str, Path],
    ) -> Mapping[str, Path]:
        if pocket.kind == "structure":
            return staged_paths
        prepared = dict(staged_paths)
        prepared["pocket"] = extract_pocket_pdb(
            staged_paths["protein"], pocket, run_dir / "prepared" / "pocket.pdb"
        )
        return prepared
