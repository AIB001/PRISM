from pathlib import Path
from typing import Any, Mapping

from ..errors import InputError
from ..types import GenerationRequest, PocketSpec
from .base import ModelAdapter


class DiffSBDDAdapter(ModelAdapter):
    """Equivariant diffusion SBDD model (Schneuing et al.).

    Upstream defines the pocket from a reference ligand -- either a
    ``<chain>:<resi>`` selector inside the receptor PDB or a standalone SDF.
    PRISM always passes the SDF, so the ligand is a genuine conditioning input
    rather than a way of locating the pocket, which makes this a
    reference-guided model in the same sense as MolCRAFT and FLOWR.
    """

    model_id = "diffsbdd"
    generation_mode = "reference_guided"
    accepted_pocket_kinds = {"reference_ligand", "structure", "center"}
    requires_reference_ligand = True
    accepted_reference_ligand_suffixes = {".sdf"}

    def prepare_inputs(
        self,
        _request: GenerationRequest,
        _pocket: PocketSpec,
        _run_dir: Path,
        staged_paths: Mapping[str, Path],
    ) -> Mapping[str, Path]:
        if staged_paths["protein"].suffix.lower() != ".pdb":
            raise InputError("DiffSBDD pocket extraction requires a PDB protein")
        return staged_paths

    def command_context(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        _prepared_paths: Mapping[str, Path],
    ) -> Mapping[str, Any]:
        # Exposed so a configuration may pin the ligand size to the reference
        # rather than sampling it from the learned size distribution.
        return {"diffsbdd_reference_heavy_atoms": _reference_heavy_atoms(pocket.reference_ligand)}


def _reference_heavy_atoms(path: Path) -> int:
    """Heavy-atom count from a V2000 counts line, or 0 when unavailable."""
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
