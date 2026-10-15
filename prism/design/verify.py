"""Evaluate raw generated graphs/poses without repairing or aligning them."""

import numpy as np
from rdkit import Chem, rdBase

from .chem import canonical_smiles, chemistry_copy, connected_subset, descriptors, heavy_copy, map_to_index


def _compare(measured, operator, requested):
    if operator == "eq":
        return measured == requested
    if operator == "le":
        return measured <= requested
    if operator == "between":
        return requested[0] <= measured <= requested[1]
    return False


def _check(ident, priority, status, measured=None, requested=None, detail=None):
    return {
        "id": ident,
        "priority": priority,
        "status": status,
        "measured": measured,
        "requested": requested,
        "detail": detail,
    }


def _mapping_candidates(grounded, candidate, max_matches):
    ref = heavy_copy(grounded.molecule)
    ref_maps = map_to_index(ref)
    core_maps = grounded.spec["selectors"]["core"]["atom_maps"]
    query = Chem.MolFromSmarts(
        Chem.MolFragmentToSmarts(ref, atomsToUse=[ref_maps[m] for m in core_maps], isomericSmarts=True)
    )
    query_maps = [a.GetAtomMapNum() for a in query.GetAtoms()]
    matches = candidate.GetSubstructMatches(query, uniquify=False, useChirality=True, maxMatches=max_matches)
    source_xyz = ref.GetConformer().GetPositions()
    target_xyz = candidate.GetConformer().GetPositions()
    retained = grounded.spec["editing"]["attachment_rules"][0]["retained_atom_map"]
    results = []
    for match in matches:
        mapping = dict(zip(query_maps, match))
        core = set(match)
        unchanged = True
        for amap, idx in mapping.items():
            a, b = ref.GetAtomWithIdx(ref_maps[amap]), candidate.GetAtomWithIdx(idx)
            if (a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsAromatic(), a.GetIsotope()) != (
                b.GetAtomicNum(),
                b.GetFormalCharge(),
                b.GetIsAromatic(),
                b.GetIsotope(),
            ):
                unchanged = False
        # A substructure match alone need not be an induced subgraph: disallow
        # extra bonds between fixed atoms as well as changed bond orders.
        for i, amap in enumerate(core_maps):
            for bmap in core_maps[i + 1 :]:
                a = ref.GetBondBetweenAtoms(ref_maps[amap], ref_maps[bmap])
                b = candidate.GetBondBetweenAtoms(mapping[amap], mapping[bmap])
                if (a is None) != (b is None) or (a is not None and a.GetBondType() != b.GetBondType()):
                    unchanged = False
        if not unchanged:
            continue
        outside = set(range(candidate.GetNumAtoms())) - core
        boundary = []
        for bond in candidate.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if (a in core) != (b in core):
                boundary.append((a if a in core else b, b if a in core else a, bond.GetBondTypeAsDouble()))
        correct_boundary = (
            len(boundary) == 1
            and boundary[0][0] == mapping[retained]
            and boundary[0][2] == 1
            and connected_subset(candidate, outside)
        )
        delta = np.array([target_xyz[mapping[m]] - source_xyz[ref_maps[m]] for m in core_maps])
        rmsd = float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))
        results.append(
            {
                "mapping": mapping,
                "core_rmsd_A": rmsd,
                "boundary_correct": correct_boundary,
                "boundary": boundary,
                "new_fragment_heavy_atoms": len(outside),
            }
        )
    # Legal boundary mappings take precedence over a low RMSD at a wrong site.
    results.sort(key=lambda r: (not r["boundary_correct"], r["core_rmsd_A"]))
    return results, len(matches) >= max_matches


def verify_candidate(grounded, molecule, max_matches=2048):
    spec = grounded.spec
    result = {
        "schema": "prism-design-validation/1",
        "design_id": spec["design_id"],
        "evaluation_context": "raw_generated_pose",
        "ligand_superposition": False,
        "chemical_repair": False,
        "rdkit_version": rdBase.rdkitVersion,
        "checks": [],
        "core_mapping": None,
        "design_hard_pass": False,
        "all_soft_objectives_pass": False,
        "selection_pass": False,
        "experimental_activity_improvement": None,
    }
    try:
        if molecule is None:
            raise ValueError("No decoded molecule")
        mol = chemistry_copy(molecule)
        if len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError("Generated molecule has multiple components")
        mol = heavy_copy(mol, keep_maps=False)
        if mol.GetNumConformers() != 1 or not np.isfinite(mol.GetConformer().GetPositions()).all():
            raise ValueError("Missing or non-finite raw 3D coordinates")
    except Exception as exc:
        for c in spec["constraints"]:
            result["checks"].append(
                _check(
                    c["id"], c["priority"], "fail" if c["kind"] == "chemical_validity" else "unknown", detail=str(exc)
                )
            )
        result["status"] = "REJECTED"
        return result

    values = descriptors(mol)
    baseline = descriptors(grounded.molecule)
    result.update({"canonical_smiles": canonical_smiles(mol), "descriptors": values, "reference_descriptors": baseline})
    maps, truncated = _mapping_candidates(grounded, mol, max_matches)
    best = maps[0] if maps else None
    result["mapping_search_truncated"] = truncated
    result["core_mapping"] = best
    unspecified = [
        idx
        for idx, code in Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
        if code == "?"
    ]
    core_indices = set(best["mapping"].values()) if best else set()
    # A formerly achiral growth atom can become a new center while remaining
    # in the fixed atom set. Count it too; excluding all fixed atoms is wrong.
    reference = chemistry_copy(grounded.molecule)
    reference_unspecified_maps = {
        grounded.molecule.GetAtomWithIdx(idx).GetAtomMapNum()
        for idx, code in Chem.FindMolChiralCenters(reference, includeUnassigned=True, useLegacyImplementation=False)
        if code == "?"
    }
    preexisting_unspecified = {best["mapping"][m] for m in reference_unspecified_maps if best and m in best["mapping"]}
    new_unspecified = len(set(unspecified) - preexisting_unspecified)
    result["fixed_candidate_indices"] = sorted(core_indices)
    for c in spec["constraints"]:
        kind = c["kind"]
        unknown = False
        if kind == "chemical_validity":
            measured = True
        elif kind == "fixed_subgraph":
            measured = best is not None
            unknown = best is None and truncated
        elif kind == "single_attachment_boundary":
            measured = 1 if best and best["boundary_correct"] else 0
            unknown = best is None and truncated
        elif kind == "formal_charge":
            measured = values["formal_charge"]
        elif kind == "new_unspecified_stereocenter_count":
            measured = new_unspecified
        elif kind == "core_rmsd":
            measured = best["core_rmsd_A"] if best else None
            unknown = best is None
        elif kind == "new_fragment_heavy_atoms":
            measured = best["new_fragment_heavy_atoms"] if best else None
            unknown = best is None
        elif kind == "descriptor_delta":
            measured = values[c["descriptor"]] - baseline[c["descriptor"]]
        else:
            measured, unknown = None, True
        status = "unknown" if unknown else "pass" if _compare(measured, c["operator"], c["value"]) else "fail"
        result["checks"].append(_check(c["id"], c["priority"], status, measured, c["value"]))

    xyz = mol.GetConformer().GetPositions()
    closest = float(np.linalg.norm(xyz[:, None, :] - grounded.protein_xyz[None, :, :], axis=2).min())
    periodic = Chem.GetPeriodicTable()
    bad_bonds = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        radii = periodic.GetRcovalent(mol.GetAtomWithIdx(a).GetAtomicNum()) + periodic.GetRcovalent(
            mol.GetAtomWithIdx(b).GetAtomicNum()
        )
        ratio = float(np.linalg.norm(xyz[a] - xyz[b]) / radii)
        if not 0.6 <= ratio <= 1.45:
            bad_bonds.append({"indices_0based": [a, b], "covalent_radius_ratio": ratio})
    changed = result["canonical_smiles"] != canonical_smiles(grounded.molecule)
    pose_ok = closest >= 2.0 and not bad_bonds
    result["selection_quality"] = {
        "nearest_protein_heavy_distance_A": closest,
        "bad_bond_lengths": bad_bonds,
        "pose_screen_pass": pose_ok,
        "changed_from_reference": changed,
        "protocol": "pilot-extra-screen/1; min protein distance 2A; covalent radius ratio 0.60..1.45",
        "not_a_binding_energy_or_full_pose_validation": True,
    }
    result["design_hard_pass"] = all(c["status"] == "pass" for c in result["checks"] if c["priority"] == "hard")
    soft = [c for c in result["checks"] if c["priority"] == "soft"]
    result["all_soft_objectives_pass"] = all(c["status"] == "pass" for c in soft)
    result["soft_objective_count"] = len(soft)
    result["selection_pass"] = result["design_hard_pass"] and result["all_soft_objectives_pass"] and changed and pose_ok
    result["status"] = "PASS" if result["selection_pass"] else "REJECTED"
    return result
