"""Bind a reviewed annotation to real atoms before it becomes a run request."""

import hashlib
from dataclasses import dataclass

import numpy as np
from rdkit import Chem, rdBase

from .chem import canonical_smiles, chemistry_copy, connected_subset, load_reference, map_to_index
from .errors import DesignError, require
from .geometry import pocket_review
from .io import read_json, verify_ref
from .schema import load_design

PROTOCOLS = {
    "chemical_validity": ("boolean", "rdkit-sanitize-and-single-component/1"),
    "fixed_subgraph": ("boolean", "mapped-core-graph-and-stereo/1"),
    "single_attachment_boundary": ("count", "mapped-boundary/1"),
    "formal_charge": ("elementary_charge", "rdkit-formal-charge/1"),
    "new_unspecified_stereocenter_count": ("count", "mapped-stereo-difference/1"),
    "core_rmsd": ("angstrom", "receptor-frame-core-rmsd-no-ligand-alignment/1"),
    "new_fragment_heavy_atoms": ("count", "mapped-new-fragment-size/1"),
}


@dataclass
class GroundedDesign:
    document: object
    molecule: object
    files: dict
    geometry: dict
    protein_xyz: object
    execution_checked: bool = False

    @property
    def spec(self):
        return self.document.data

    def report(self):
        return {
            "status": "VALIDATED",
            "design_id": self.spec["design_id"],
            "task": self.spec["task"],
            "spec_sha256": self.document.sha256,
            "review_tier": self.spec["review"]["review_tier"],
            "eligible_for_compilation": self.execution_checked,
            "atom_map_to_source_index_0based": map_to_index(self.molecule),
            "geometry": self.geometry,
            "rdkit_version": rdBase.rdkitVersion,
            "affinity_or_design_improvement_verified": False,
        }


def validate_selectors(spec, mol):
    maps = map_to_index(mol)
    core_list = spec["selectors"]["core"]["atom_maps"]
    edit_list = spec["selectors"]["editable"]["atom_maps"]
    core, edit = set(core_list), set(edit_list)
    require(core | edit == set(maps), "Masks must cover exactly the reference heavy atoms", "MASK_COVERAGE")
    require(not core & edit, "Fixed and editable masks overlap", "MASK_OVERLAP")
    require(connected_subset(mol, [maps[m] for m in core]), "Fixed core is disconnected", "CORE_DISCONNECTED")
    require(len(core) >= 5, "Pilot requires at least five fixed heavy atoms", "CORE_TOO_SMALL")
    rules = spec["editing"]["attachment_rules"][0]
    retained, replaced = rules["retained_atom_map"], rules["replaced_atom_map"]
    require(retained in core, "Retained attachment must belong to the fixed core", "ATTACHMENT_INVALID")
    if spec["task"] == "fragment_replacement":
        require(
            connected_subset(mol, [maps[m] for m in edit]), "Editable fragment is disconnected", "EDIT_DISCONNECTED"
        )
        require(2 <= len(edit) <= 12 and len(edit) <= 0.4 * len(maps), "Edit exceeds pilot size limits", "EDIT_SIZE")
        require(replaced in edit, "Replaced attachment must belong to the edited fragment", "ATTACHMENT_INVALID")
        boundaries = [
            b
            for b in mol.GetBonds()
            if b.GetBeginAtom().GetAtomicNum() > 1
            and b.GetEndAtom().GetAtomicNum() > 1
            and ((b.GetBeginAtom().GetAtomMapNum() in core) != (b.GetEndAtom().GetAtomMapNum() in core))
        ]
        require(len(boundaries) == 1, "Expected exactly one edit boundary", "ATTACHMENT_COUNT")
        b = boundaries[0]
        require(
            {b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum()} == {retained, replaced},
            "Attachment does not identify the actual boundary",
            "ATTACHMENT_INVALID",
        )
        require(
            not b.IsInRing() and b.GetBondTypeAsDouble() == 1, "Only non-ring single cuts supported", "CUT_UNSUPPORTED"
        )
        clean = chemistry_copy(mol)
        require(
            clean.GetAtomWithIdx(maps[retained]).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED,
            "A cut at a specified stereocenter needs a boundary stereo protocol",
            "CUT_STEREO_UNSUPPORTED",
        )
        fixed_indices = {maps[m] for m in core}
        for stereo_bond in clean.GetBonds():
            if stereo_bond.GetStereo() not in {Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY}:
                defining = {stereo_bond.GetBeginAtomIdx(), stereo_bond.GetEndAtomIdx(), *stereo_bond.GetStereoAtoms()}
                require(
                    not defining & fixed_indices or defining <= fixed_indices,
                    "Cut crosses a double-bond stereo definition",
                    "CUT_STEREO_UNSUPPORTED",
                )
        for smarts in ("[CX3](=[OX1])-[NX3]", "[SX4](=[OX1])(=[OX1])-[NX3]"):
            for match in mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
                require(
                    {match[0], match[-1]} != {maps[retained], maps[replaced]},
                    "Amide/sulfonamide bond cannot be the cut",
                    "PROTECTED_BOND",
                )
        coords = mol.GetConformer().GetPositions()
        require(
            np.allclose(
                coords[maps[replaced]] - coords[maps[retained]], rules["reference_boundary_vector_A"], atol=1e-5
            ),
            "Boundary vector does not match reference coordinates",
            "BOUNDARY_VECTOR",
        )
    else:
        require(not edit and replaced is None, "Growth must retain every source heavy atom", "GROWTH_MASK")
        require(spec["selectors"]["growth_site"]["atom_maps"] == [retained], "Growth selector mismatch", "GROWTH_MASK")
        clean = chemistry_copy(mol)
        atom = clean.GetAtomWithIdx(maps[retained])
        require(
            atom.GetAtomicNum() == 6 and atom.GetTotalNumHs(includeNeighbors=True) > 0,
            "Growth site must be a carbon bearing hydrogen",
            "GROWTH_VALENCE",
        )
        require(
            atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED,
            "Growth at a preassigned tetrahedral center is unsupported",
            "GROWTH_STEREO",
        )


def validate_constraints(spec, mol):
    ids = [c["id"] for c in spec["constraints"]]
    require(len(set(ids)) == len(ids), "Constraint IDs must be unique", "DUPLICATE_CONSTRAINT")
    for kind, (unit, protocol) in PROTOCOLS.items():
        found = [c for c in spec["constraints"] if c["kind"] == kind]
        require(
            len(found) == 1 and found[0]["priority"] == "hard",
            f"Required hard constraint missing/duplicated: {kind}",
            "CONSTRAINT_REQUIRED",
        )
        require(
            found[0]["unit"] == unit and found[0]["verifier_protocol_ref"] == protocol,
            f"Unsupported unit/protocol for {kind}",
            "PROTOCOL_UNSUPPORTED",
        )
        selector = "core" if kind in {"fixed_subgraph", "core_rmsd"} else None
        context = "raw_generated_pose" if kind == "core_rmsd" else "chemical_graph"
        require(
            found[0]["selector_ref"] == selector and found[0]["evaluation_context"] == context,
            f"Unsupported selector/evaluation context for {kind}",
            "PROTOCOL_UNSUPPORTED",
        )
    for c in spec["constraints"]:
        kind = c["kind"]
        if kind == "descriptor_delta":
            valid = {
                "crippen_logp": ("dimensionless", "rdkit-crippen-delta/1"),
                "rotatable_bonds_strict": ("count", "rdkit-strict-rotors-delta/1"),
            }
            require(
                c["descriptor"] in valid and (c["unit"], c["verifier_protocol_ref"]) == valid[c["descriptor"]],
                "Unsupported descriptor protocol",
                "PROTOCOL_UNSUPPORTED",
            )
            require(
                c["evaluation_context"] == "chemical_graph"
                and c["selector_ref"] is None
                and c["operator"] in {"le", "eq"}
                and type(c["value"]) in {float, int},
                "Descriptor deltas require a numeric whole-molecule graph target",
                "CONSTRAINT_CONFLICT",
            )
        elif kind == "formal_charge":
            require(
                c["value"] == Chem.GetFormalCharge(mol) and c["operator"] == "eq",
                "Charge differs from reference state",
                "CHARGE_MISMATCH",
            )
        elif kind in {"fixed_subgraph", "chemical_validity"}:
            require(
                c["value"] is True and c["operator"] == "eq",
                "Validity and core preservation must be required",
                "CONSTRAINT_CONFLICT",
            )
        elif kind == "single_attachment_boundary":
            require(
                c["value"] == 1 and c["operator"] == "eq", "Only one new boundary bond supported", "CONSTRAINT_CONFLICT"
            )
        elif kind == "new_unspecified_stereocenter_count":
            require(
                c["value"] == 0 and c["operator"] == "eq",
                "Unspecified new stereo is unsupported",
                "CONSTRAINT_CONFLICT",
            )
        elif kind == "core_rmsd":
            require(
                c["operator"] == "le"
                and 0 < c["value"] <= 1
                and c["selector_ref"] == "core"
                and c["evaluation_context"] == "raw_generated_pose",
                "Invalid core RMSD policy",
                "CONSTRAINT_CONFLICT",
            )
        elif kind == "new_fragment_heavy_atoms":
            require(
                c["operator"] == "between"
                and c["value"] == spec["editing"]["new_fragment_heavy_atoms_range"]
                and 0 < c["value"][0] <= c["value"][1] <= 12,
                "Fragment size range conflicts",
                "CONSTRAINT_CONFLICT",
            )


def ground_design(document, for_execution=False):
    spec = document.data
    require(spec["task"] != "input_review", "This record requests input review, not generation", "INPUT_REVIEW_ONLY")
    files = {}
    for field, ref in spec["inputs"].items():
        if isinstance(ref, dict) and "sha256" in ref:
            files[field] = verify_ref(document.root, ref)
    require(
        {"receptor", "pocket", "reference_ligand", "atom_mapping", "evidence"} <= files.keys(),
        "Generation requires receptor, pocket, mapped reference, atom mapping and evidence",
        "MISSING_INPUT",
    )
    for ref in [
        spec["review"]["review_record"],
        spec["provenance"]["selection_manifest"],
        spec["provenance"]["manual_decisions"],
    ]:
        verify_ref(document.root, ref)
    review = read_json(verify_ref(document.root, spec["review"]["review_record"]))
    require(
        review["pdb_id"] == spec["complex_id"] and review["review_tier"] == spec["review"]["review_tier"],
        "Review identity/tier does not match DesignSpec",
        "REVIEW_BINDING",
    )
    require(
        verify_ref(document.root, review["evidence"]) == files["evidence"],
        "Review points to different evidence",
        "REVIEW_BINDING",
    )
    evidence = read_json(files["evidence"])
    require(evidence["pdb_id"] == spec["complex_id"], "Evidence belongs to another complex", "EVIDENCE_BINDING")
    require(
        set(evidence["quality_flags"]) == set(spec["review"]["quality_flags"]),
        "Source quality flags were changed without updating evidence",
        "EVIDENCE_BINDING",
    )
    require(
        review["objective_code"] == spec["objective"]["code"],
        "Objective differs from the reviewed decision",
        "REVIEW_BINDING",
    )
    choice = review["selected_candidate_evidence"]
    rule = spec["editing"]["attachment_rules"][0]
    if spec["task"] == "fragment_replacement":
        require(
            set(choice["fixed_atom_maps"]) == set(spec["selectors"]["core"]["atom_maps"])
            and set(choice["editable_atom_maps"]) == set(spec["selectors"]["editable"]["atom_maps"])
            and choice["retained_attachment_atom_map"] == rule["retained_atom_map"]
            and choice["replaced_attachment_atom_map"] == rule["replaced_atom_map"],
            "Edit masks/attachment differ from the reviewed decision",
            "REVIEW_BINDING",
        )
    else:
        require(
            choice["attachment_atom_map"] == rule["retained_atom_map"],
            "Growth site differs from the reviewed decision",
            "REVIEW_BINDING",
        )
    mol = load_reference(files["reference_ligand"])
    require(
        hashlib.sha256(canonical_smiles(mol).encode()).hexdigest() == spec["inputs"]["chemical_state_id"],
        "Reference chemistry/stereo does not match chemical_state_id",
        "CHEMICAL_STATE_MISMATCH",
    )
    require(
        spec["inputs"]["receptor_state_id"] == spec["inputs"]["receptor"]["sha256"]
        and spec["inputs"]["pose_id"] == spec["inputs"]["reference_ligand"]["sha256"],
        "Receptor/pose state identity mismatch",
        "STATE_ID_MISMATCH",
    )
    mapping = read_json(files["atom_mapping"])
    maps = map_to_index(mol)
    require(mapping["pdb_id"] == spec["complex_id"], "Atom map belongs to another complex", "ATOM_MAP_INVALID")
    require(
        len(mapping["atoms"]) == len(maps) and {a["atom_map"] for a in mapping["atoms"]} == set(maps),
        "Atom map table coverage mismatch or duplicates",
        "ATOM_MAP_INVALID",
    )
    coords = mol.GetConformer().GetPositions()
    for atom in mapping["atoms"]:
        idx = maps[atom["atom_map"]]
        require(
            mol.GetAtomWithIdx(idx).GetSymbol() == atom["element"]
            and np.allclose(coords[idx], atom["xyz_angstrom"], atol=1e-5),
            "Atom map element/coordinates disagree",
            "ATOM_MAP_INVALID",
        )
    validate_selectors(spec, mol)
    validate_constraints(spec, mol)
    geom, protein_xyz = pocket_review(files["receptor"], files["pocket"], coords[list(maps.values())])
    if for_execution:
        require(
            spec["review"]["review_tier"] == "candidate_for_m0_pilot",
            "Only reviewed M0 pilot candidates may execute; curate a new reviewed annotation first",
            "REVIEW_REQUIRED",
        )
        require(
            geom["nearest_protein_heavy_distance_A"] >= 2 and not geom["metals_within_3A"],
            "Short protein contact or metal context requires a separate protocol",
            "SPECIAL_CHEMISTRY",
        )
        require(
            not (set(spec["review"]["quality_flags"]) - {"MOL2_PARSE_FAILED"}),
            "Unresolved source quality flags prevent pilot execution",
            "SOURCE_QUALITY_REVIEW",
        )
    return GroundedDesign(document, mol, files, geom, protein_xyz, for_execution)


def validate_design(path, input_root=None, for_execution=False):
    try:
        grounded = ground_design(load_design(path, input_root), for_execution=for_execution)
        return grounded.report()
    except (DesignError, OSError, ValueError, KeyError, TypeError, IndexError) as exc:
        return {
            "status": "REJECTED",
            "error_code": getattr(exc, "code", "INVALID_INPUT"),
            "message": str(exc),
            "eligible_for_compilation": False,
        }
