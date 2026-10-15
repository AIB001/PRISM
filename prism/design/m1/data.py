"""Portable, hash-bound reconstruction records and leakage-aware group splits.

No PDBbind-specific paths or labels are required by the trainer. An observed
caption is supervision for reconstruction, never an improvement preference.
"""

from collections import defaultdict
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

from ..chem import chemistry_copy, connected_subset, heavy_copy, load_reference
from ..geometry import pocket_review
from ..io import _nonfinite, _pairs, read_json, sha256, verify_ref, write_json

SCHEMA = "prism-design-m1-record/1"
FIELDS = {
    "schema",
    "record_id",
    "complex_group",
    "task",
    "text",
    "label_kind",
    "review",
    "files",
    "atom_maps",
    "mask_candidates",
    "split_keys",
    "provenance",
}


@lru_cache(maxsize=1)
def record_validator():
    from jsonschema import Draft202012Validator

    return Draft202012Validator(read_json(Path(__file__).parents[1] / "schemas" / "m1-record-v1.json"))


def chemical_keys(mol):
    # Used ONLY for grouping. Never replace the training graph/charge/pose.
    clean = Chem.Mol(mol)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    clean = rdMolStandardize.Cleanup(clean)
    clean = rdMolStandardize.FragmentParent(clean)
    clean = rdMolStandardize.Uncharger().uncharge(clean)
    clean = rdMolStandardize.TautomerEnumerator().Canonicalize(clean)
    Chem.RemoveStereochemistry(clean)  # Enantiomers/tautomers never straddle splits.
    identity = Chem.MolToSmiles(clean)
    core = MurckoScaffold.GetScaffoldForMol(clean)
    # Empty acyclic scaffolds must NOT connect every acyclic molecule.
    scaffold = Chem.MolToSmiles(core) if core.GetNumAtoms() else "acyclic:" + identity
    return {"chemical_identity": identity, "scaffold": scaffold}


def validate_cut(mol, maps, cut):
    index = {v: i for i, v in enumerate(maps)}
    clean = chemistry_copy(mol)
    fixed = {index[i] for i in cut["fixed"]}
    variable = {index[i] for i in cut["variable"]}
    edges = [b for b in mol.GetBonds() if (b.GetBeginAtomIdx() in fixed) != (b.GetEndAtomIdx() in fixed)]
    if (
        not connected_subset(mol, fixed)
        or not connected_subset(mol, variable)
        or len(edges) != 1
        or edges[0].GetBondType() != Chem.BondType.SINGLE
        or edges[0].IsInRing()
        or {edges[0].GetBeginAtomIdx(), edges[0].GetEndAtomIdx()}
        != {index[cut["attachment_fixed"]], index[cut["attachment_variable"]]}
    ):
        raise ValueError("Mask does not describe a connected single-boundary fragment")
    boundary = {index[cut["attachment_fixed"]], index[cut["attachment_variable"]]}
    if clean.GetAtomWithIdx(index[cut["attachment_fixed"]]).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
        raise ValueError("A fixed boundary stereocentre requires a dedicated protocol")
    for smarts in ("[CX3](=[OX1])-[NX3]", "[SX4](=[OX1])(=[OX1])-[NX3]"):
        if any({m[0], m[-1]} == boundary for m in clean.GetSubstructMatches(Chem.MolFromSmarts(smarts))):
            raise ValueError("Protected amide/sulfonamide boundary")
    for bond in clean.GetBonds():
        if bond.GetStereo() not in {Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY}:
            defining = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), *bond.GetStereoAtoms()}
            if defining & fixed and not defining <= fixed:
                raise ValueError("Boundary crosses a double-bond stereo definition")


def validate_record(record, root, check_files=True):
    errors = list(record_validator().iter_errors(record))
    if errors:
        raise ValueError("M1 JSON schema: " + errors[0].message)
    if set(record) != FIELDS or record["schema"] != SCHEMA:
        raise ValueError("Unknown or incomplete M1 record schema")
    if record["task"] != "fragment_reconstruction" or record["label_kind"] != "observed_structure_caption":
        raise ValueError("Only observed fragment reconstruction is implemented; no synthetic improvement labels")
    if not all(isinstance(record[k], str) and record[k].strip() for k in ("record_id", "complex_group", "text")):
        raise ValueError("Record identity and caption must be nonempty strings")
    if len(record["text"]) > 4000:
        raise ValueError("Caption exceeds the record length limit")
    review = record["review"]
    if (
        review.get("status") != "completed"
        or review.get("decision") != "nonpeptide"
        or review.get("quality_pass") is not True
    ):
        raise ValueError("Training requires completed nonpeptide review and a passed source-quality gate")
    if set(record["files"]) != {"reference", "pocket", "protein"}:
        raise ValueError("Reference, supplied pocket and protein are required")
    if not record["split_keys"].get("protein_groups"):
        raise ValueError("Protein grouping evidence is required")
    maps = record["atom_maps"]
    if not maps or any(type(v) is not int or v <= 0 for v in maps) or len(set(maps)) != len(maps):
        raise ValueError("Invalid stable source atom IDs")
    candidates = record["mask_candidates"]
    if not candidates:
        raise ValueError("No reviewed fragment-reconstruction mask")
    for cut in candidates:
        if set(cut) != {"fixed", "variable", "attachment_fixed", "attachment_variable"}:
            raise ValueError("Unknown mask fields")
        f, u = cut["fixed"], cut["variable"]
        if (
            len(f) != len(set(f))
            or len(u) != len(set(u))
            or set(f) & set(u)
            or set(f) | set(u) != set(maps)
            or len(f) < 5
            or not 2 <= len(u) <= 12
            or len(u) / len(maps) > 0.4
            or cut["attachment_fixed"] not in f
            or cut["attachment_variable"] not in u
        ):
            raise ValueError("Invalid fixed/variable partition or attachment")
    if not check_files:
        return record
    files = {key: verify_ref(root, value) for key, value in record["files"].items()}
    mol = heavy_copy(load_reference(files["reference"]))
    if [a.GetAtomMapNum() for a in mol.GetAtoms()] != maps:
        raise ValueError("Heavy-atom order/source maps changed")
    if chemical_keys(mol) != {k: record["split_keys"][k] for k in ("chemical_identity", "scaffold")}:
        raise ValueError("Chemical grouping keys disagree with the reference")
    for cut in candidates:
        validate_cut(mol, maps, cut)
    if not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise ValueError("Nonfinite target coordinates")
    geometry, _ = pocket_review(files["protein"], files["pocket"], mol.GetConformer().GetPositions())
    if geometry["nearest_protein_heavy_distance_A"] < 2.0 or geometry["metals_within_3A"]:
        raise ValueError("Source complex needs a clash/covalent/metal-specific protocol")
    return record


def read_records(manifest, check_files=True):
    manifest = Path(manifest).resolve()
    records = []
    seen = set()
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        r = json.loads(line, object_pairs_hook=_pairs, parse_constant=_nonfinite)
        validate_record(r, manifest.parent, check_files)
        if r["record_id"] in seen:
            raise ValueError("Duplicate record_id: " + r["record_id"])
        seen.add(r["record_id"])
        records.append(r)
    if not records:
        raise ValueError("Empty M1 dataset")
    return records


def grouped_split(records, seed=20261013, ratios=(0.8, 0.1, 0.1)):
    """Connected components of identity OR scaffold OR protein OR complex.

    Never split a large component to achieve an attractive fraction. This is
    exact-protein/UniProt protection unless supplied groups are homolog clusters.
    """
    if len(ratios) != 3 or min(ratios) <= 0 or abs(sum(ratios) - 1) > 1e-9:
        raise ValueError("Three positive train/val/test fractions must sum to one")
    parent = list(range(len(records)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owners = {}
    for i, r in enumerate(records):
        keys = r["split_keys"]
        groups = [
            ("complex", r["complex_group"]),
            ("identity", keys["chemical_identity"]),
            ("scaffold", keys["scaffold"]),
        ]
        groups += [("protein", p) for p in keys["protein_groups"]]
        for key in groups:
            if key in owners:
                parent[find(i)] = find(owners[key])
            owners[key] = i
    components = defaultdict(list)
    for i, r in enumerate(records):
        components[find(i)].append(r["record_id"])
    if len(components) < 3:
        raise ValueError("Fewer than three independent components; no honest three-way split")
    stable = lambda ids: hashlib.sha256((str(seed) + "|" + "|".join(sorted(ids))).encode()).hexdigest()
    ordered = sorted(components.values(), key=lambda ids: (-len(ids), stable(ids)))
    names, counts, result = ("train", "val", "test"), [0, 0, 0], {}
    targets = [len(records) * fraction for fraction in ratios]
    for component in ordered:
        idx = max(range(3), key=lambda j: (targets[j] - counts[j]) / targets[j])
        for ident in component:
            result[ident] = {"split": names[idx], "component": stable(component)}
        counts[idx] += len(component)
    if min(counts) == 0:
        raise ValueError("Grouping produced an empty split; change the dataset, not the grouping guarantee")
    return {
        "seed": seed,
        "ratios_requested": list(ratios),
        "counts": dict(zip(names, counts)),
        "n_components": len(ordered),
        "largest_component": max(map(len, ordered)),
        "assignments": result,
        "protocol": "connected chemical identity/scaffold + supplied protein groups + complex",
        "homology_disjoint_claim": False,
        "pretrained_backbone_contamination": "unknown; not a decontaminated generalization benchmark",
    }


def validate_split(records, split):
    assignment = split["assignments"]
    if set(assignment) != {r["record_id"] for r in records}:
        raise ValueError("Split does not cover exactly this dataset")
    groups = {}
    for r in records:
        part = assignment[r["record_id"]]["split"]
        if part not in {"train", "val", "test"}:
            raise ValueError("Unknown split")
        k = r["split_keys"]
        for key in [
            ("complex", r["complex_group"]),
            ("chemical", k["chemical_identity"]),
            ("scaffold", k["scaffold"]),
        ] + [("protein", p) for p in k["protein_groups"]]:
            if key in groups and groups[key] != part:
                raise ValueError(f"Cross-split leakage: {key}")
            groups[key] = part
    if {v["split"] for v in assignment.values()} != {"train", "val", "test"}:
        raise ValueError("Empty split")


def prepare_manifest(manifest, split_path, seed=20261013):
    records = read_records(manifest)
    split = grouped_split(records, seed)
    split["manifest_sha256"] = sha256(manifest)
    validate_split(records, split)
    write_json(split_path, split)
    return split
