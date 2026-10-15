"""Compile a hash-bound annotation into an explicit backend-independent edit."""

from pathlib import Path
import shutil

import numpy as np
from rdkit import Chem, rdBase

from .chem import canonical_smiles, chemistry_copy, heavy_copy, map_to_index
from .errors import require
from .io import file_ref, write_json


def growth_prior(grounded):
    """Propose a direction from an actual reference C-H vector, never a hard claim."""
    spec = grounded.spec
    amap = spec["editing"]["attachment_rules"][0]["retained_atom_map"]
    idx = map_to_index(grounded.molecule)[amap]
    mol = Chem.AddHs(chemistry_copy(grounded.molecule), addCoords=True)
    positions = mol.GetConformer().GetPositions()
    atom = mol.GetAtomWithIdx(idx)
    candidates = []
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() != 1:
            continue
        vector = positions[neighbor.GetIdx()] - positions[idx]
        length = np.linalg.norm(vector)
        if length < 0.1:
            continue
        direction = vector / length
        center = positions[idx] + 1.6 * direction
        distance = np.linalg.norm(grounded.protein_xyz - center, axis=1).min()
        candidates.append((float(distance), neighbor.GetIdx(), center, direction))
    require(candidates, "No valid C-H direction for growth prior", "GROWTH_DIRECTION_UNKNOWN")
    _, hydrogen_idx, center, direction = max(candidates, key=lambda v: (v[0], -v[1]))
    return {
        "attachment_atom_map": amap,
        "center_A": center.tolist(),
        "unit_direction": direction.tolist(),
        "hydrogen_index_in_prepared_reference": hydrogen_idx,
        "method": "reference_C-H_direction_with_largest_minimum_protein_distance_at_1.6A",
        "enforcement": "soft_prior_hint; exact_attachment_checked_after_generation",
    }


def compile_design(grounded, directory):
    require(
        grounded.execution_checked,
        "Run input validation with for_execution=True before compilation",
        "EXECUTION_NOT_VALIDATED",
    )
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    require(not (directory / "compiled_design.json").exists(), "Compiled request already exists", "OUTPUT_EXISTS")
    inputs = directory / "inputs"
    inputs.mkdir(exist_ok=True)
    molecule = heavy_copy(grounded.molecule)
    mapping = map_to_index(molecule)
    selectors = grounded.spec["selectors"]
    before = canonical_smiles(molecule)
    source_positions = molecule.GetConformer().GetPositions().copy()
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    reference = inputs / "reference_backend.sdf"
    with Chem.SDWriter(str(reference)) as writer:
        writer.write(molecule)
    reread = Chem.MolFromMolFile(str(reference), removeHs=False)
    require(
        canonical_smiles(reread) == before
        and np.allclose(reread.GetConformer().GetPositions(), source_positions, atol=1e-5),
        "Backend reference changed chemistry or pose",
        "STAGING_MISMATCH",
    )
    staged = {"reference_ligand": file_ref(reference)}
    for role in ("receptor", "pocket"):
        target = inputs / f"{role}.pdb"
        shutil.copyfile(grounded.files[role], target)
        staged[role] = file_ref(target)
    # Upstream reconstructs standard protein bonds. A stale CONECT entry can
    # refer to atoms absent in a supplied pocket. Preserve the source separately.
    backend_pocket = inputs / "pocket_backend.pdb"
    text = (inputs / "pocket.pdb").read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)
    backend_pocket.write_text("".join(line for line in lines if not line.startswith("CONECT")), encoding="utf-8")
    staged["backend_pocket"] = file_ref(backend_pocket)
    spec = grounded.spec
    rules = spec["editing"]["attachment_rules"][0]
    compiled = {
        "schema": "prism-design-compiled/1",
        "design_id": spec["design_id"] + "-m0",
        "parent_design_id": spec["design_id"],
        "task": spec["task"],
        "annotation": file_ref(grounded.document.path),
        "annotation_root": str(grounded.document.root),
        "inputs": staged,
        "reference_atom_map_to_backend_index_0based": {str(k): v for k, v in mapping.items()},
        "fixed_indices_0based": [mapping[m] for m in selectors["core"]["atom_maps"]],
        "editable_indices_0based": [mapping[m] for m in selectors["editable"]["atom_maps"]],
        "attachment_index_0based": mapping[rules["retained_atom_map"]],
        "new_fragment_heavy_atoms_range": spec["editing"]["new_fragment_heavy_atoms_range"],
        "constraints": spec["constraints"],
        "coordinate_frame": "original_receptor_frame_no_ligand_alignment",
        "rdkit_version": rdBase.rdkitVersion,
        "growth_prior": None,
        "transformations": {
            "ligand": "remove_H_and_atom_map_labels; preserve_heavy_order_graph_and_coordinates",
            "pocket": "remove_CONECT_only",
            "removed_CONECT_count": sum(line.startswith("CONECT") for line in lines),
        },
        "enforcement": {
            "fixed_subgraph": "native_conditional_request_plus_independent_verification",
            "attachment_charge_size_pose_stereo": "independent_postgeneration_verification",
            "descriptor_objectives": "independent_postgeneration_scoring",
        },
        "input_validation": grounded.report(),
    }
    if spec["task"] == "fragment_growing":
        prior = growth_prior(grounded)
        path = inputs / "growth_prior.xyz"
        path.write_text(
            "1\nPRISM reference C-H growth hint\nC " + " ".join(f"{x:.6f}" for x in prior["center_A"]) + "\n",
            encoding="utf-8",
        )
        prior["file"] = file_ref(path)
        compiled["growth_prior"] = prior
    write_json(directory / "compiled_design.json", compiled)
    return compiled
