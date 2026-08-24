from pathlib import Path
from typing import Mapping

from ..errors import InputError
from ..protein import extract_pocket_pdb
from ..types import GenerationRequest, PocketSpec
from .base import ModelAdapter


class MolCRAFTAdapter(ModelAdapter):
    model_id = "molcraft"
    generation_mode = "reference_guided"
    accepted_pocket_kinds = {"center", "reference_ligand", "residues", "structure"}
    requires_reference_ligand = True
    accepted_reference_ligand_suffixes = {".sdf"}

    def prepare_inputs(
        self,
        _request: GenerationRequest,
        pocket: PocketSpec,
        run_dir: Path,
        staged_paths: Mapping[str, Path],
    ) -> Mapping[str, Path]:
        reference = staged_paths.get("reference_ligand")
        if reference is None or reference.suffix.lower() != ".sdf":
            raise InputError("MolCRAFT custom-pocket inference requires a reference SDF")

        prepared = dict(staged_paths)
        if pocket.kind == "reference_ligand":
            if staged_paths["protein"].suffix.lower() != ".pdb":
                raise InputError("MolCRAFT custom-pocket inference requires a PDB protein")
            return prepared
        if pocket.kind == "structure":
            if staged_paths["pocket"].suffix.lower() != ".pdb":
                raise InputError("MolCRAFT structure pockets must be PDB files")
            prepared["protein"] = staged_paths["pocket"]
            return prepared

        if staged_paths["protein"].suffix.lower() != ".pdb":
            raise InputError("MolCRAFT center/residue preparation requires a PDB protein")
        prepared["protein"] = extract_pocket_pdb(
            staged_paths["protein"],
            pocket,
            run_dir / "prepared" / "pocket.pdb",
        )
        return prepared
