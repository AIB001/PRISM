import math
from pathlib import Path
from typing import Any, Mapping

from ..errors import ConfigurationError, InputError
from ..protein import extract_pocket_pdb
from ..types import GenerationRequest, PocketSpec
from .base import ModelAdapter


class FLOWRAdapter(ModelAdapter):
    model_id = "flowr"
    generation_mode = "reference_guided"
    accepted_pocket_kinds = {"center", "reference_ligand", "residues", "structure"}
    requires_reference_ligand = True
    accepted_reference_ligand_suffixes = {".pdb", ".sdf"}

    def prepare_inputs(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        run_dir: Path,
        staged_paths: Mapping[str, Path],
    ) -> Mapping[str, Path]:
        reference = staged_paths.get("reference_ligand")
        if reference is None or reference.suffix.lower() not in {".sdf", ".pdb"}:
            raise InputError("FLOWR requires a reference ligand in SDF or PDB format")

        prepared = dict(staged_paths)
        if pocket.kind == "reference_ligand":
            return prepared
        if pocket.kind == "structure":
            if staged_paths["pocket"].suffix.lower() != ".pdb":
                raise InputError("FLOWR structure pockets must be PDB files")
            prepared["protein"] = staged_paths["pocket"]
            return prepared

        if staged_paths["protein"].suffix.lower() != ".pdb":
            raise InputError("FLOWR center/residue preparation requires a PDB protein")
        prepared["protein"] = extract_pocket_pdb(
            staged_paths["protein"],
            pocket,
            run_dir / "prepared" / "pocket.pdb",
        )
        return prepared

    def command_context(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        _prepared_paths: Mapping[str, Path],
    ) -> Mapping[str, Any]:
        configured = self.config.get("pocket_cutoff", 6.0)
        if isinstance(configured, bool):
            raise ConfigurationError("FLOWR 'pocket_cutoff' must be positive")
        try:
            cutoff = float(configured)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError("FLOWR 'pocket_cutoff' must be positive") from exc
        if not math.isfinite(cutoff) or cutoff <= 0:
            raise ConfigurationError("FLOWR 'pocket_cutoff' must be positive")
        return {
            "flowr_pocket_cutoff": cutoff,
            "flowr_pocket_mode": "reference" if pocket.kind == "reference_ligand" else "prepared",
        }
