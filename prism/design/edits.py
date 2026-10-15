"""Execute bounded LLM graph edits. No generated coordinates or graph repair."""

from rdkit import Chem, rdBase

from .chem import canonical_smiles, chemistry_copy, connected_subset, descriptors, heavy_copy, map_to_index
from .errors import DesignError, require
from .verify import _compare


def _stereo_atoms(mol):
    """Protect every explicitly specified center and its defining neighbours."""
    result = set()
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            result.update([atom.GetAtomMapNum()] + [a.GetAtomMapNum() for a in atom.GetNeighbors()])
    for bond in mol.GetBonds():
        if bond.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
            result.update(
                mol.GetAtomWithIdx(i).GetAtomMapNum()
                for i in [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), *bond.GetStereoAtoms()]
            )
    return result


def _unspecified_stereo(mol):
    clean = chemistry_copy(mol, keep_maps=True)
    # Atom-map labels are not chemical substituents for stereo perception.
    maps = [a.GetAtomMapNum() for a in clean.GetAtoms()]
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    result = set()
    for stereo in Chem.FindPotentialStereo(clean):
        if stereo.specified != Chem.StereoSpecified.Unspecified:
            continue
        if stereo.type == Chem.StereoType.Atom_Tetrahedral:
            result.add(("atom", maps[stereo.centeredOn]))
        elif stereo.type == Chem.StereoType.Bond_Double:
            b = clean.GetBondWithIdx(stereo.centeredOn)
            result.add(("bond", *sorted([maps[b.GetBeginAtomIdx()], maps[b.GetEndAtomIdx()]])))
        else:
            raise DesignError("Unsupported stereochemistry in proposed graph", "EDIT_STEREO_UNSUPPORTED")
    return result


def execute_edits(grounded, operations):
    """Only modify the reviewed editable region or add at its fixed attachment.

    Atom maps remain stable across deletion; new identifiers cannot reuse any
    source identifier. V1 deliberately excludes aromatic and stereo rewrites.
    """
    require(
        isinstance(operations, list) and 1 <= len(operations) <= 48, "Expected 1..48 edit operations", "EDIT_OPERATIONS"
    )
    ref = heavy_copy(grounded.molecule)
    source_maps = map_to_index(ref)
    core = set(grounded.spec["selectors"]["core"]["atom_maps"])
    editable = set(grounded.spec["selectors"]["editable"]["atom_maps"])
    retained = grounded.spec["editing"]["attachment_rules"][0]["retained_atom_map"]
    protected_stereo = _stereo_atoms(ref)
    rw = Chem.RWMol(ref)
    used = set(source_maps)
    mutable = set(editable)
    # Explicit H counts inherited from an SDF must not leave e.g. a mutated
    # carbon carrying three hydrogens after becoming oxygen.
    for operation in operations:
        name = operation["op"]
        maps = map_to_index(rw)
        ids = [operation["atom_map"]] if "atom_map" in operation else operation["atom_maps"]
        require(all(type(m) is int for m in ids), "Atom maps must be integers", "ATOM_MAP_INVALID")
        require(
            not set(ids) & protected_stereo,
            "Editing specified stereo or its defining neighbours is unsupported in v1",
            "EDIT_STEREO_UNSUPPORTED",
        )
        if name == "add_atom":
            amap = ids[0]
            require(
                amap not in used and amap > max(source_maps),
                "New atom maps must exceed all source maps and cannot be reused",
                "EDIT_NEW_MAP",
            )
            require(type(operation["formal_charge"]) is int, "Charge must be an integer", "EDIT_CHARGE")
            atom = Chem.Atom(operation["element"])
            atom.SetFormalCharge(operation["formal_charge"])
            atom.SetAtomMapNum(amap)
            rw.AddAtom(atom)
            used.add(amap)
            mutable.add(amap)
            continue
        require(all(m in maps for m in ids), "Edit references an absent atom map", "EDIT_UNKNOWN_ATOM")
        if name in {"delete_atom", "mutate_atom"}:
            amap = ids[0]
            require(amap in mutable and amap not in core, "Cannot change a fixed atom", "EDIT_FIXED_CORE")
            atom = rw.GetAtomWithIdx(maps[amap])
            require(
                not atom.GetIsAromatic(), "Aromatic atom rewrites are unsupported in v1", "EDIT_AROMATIC_UNSUPPORTED"
            )
            if name == "delete_atom":
                rw.RemoveAtom(maps[amap])
                mutable.remove(amap)
            else:
                require(type(operation["formal_charge"]) is int, "Charge must be an integer", "EDIT_CHARGE")
                atom.SetAtomicNum(Chem.GetPeriodicTable().GetAtomicNumber(operation["element"]))
                atom.SetFormalCharge(operation["formal_charge"])
                atom.SetNumExplicitHs(0)
                atom.SetNoImplicit(False)
        elif name in {"add_bond", "delete_bond", "change_bond"}:
            require(len(set(ids)) == 2, "Self bonds are forbidden", "EDIT_BOND")
            fixed = set(ids) & core
            require(
                not fixed or (fixed == {retained} and len(set(ids) & mutable) == 1),
                "Only the declared fixed attachment may connect to edited atoms",
                "EDIT_FIXED_CORE",
            )
            require(set(ids) <= core | mutable, "Bond outside the editable region", "EDIT_SCOPE")
            require(
                all(not rw.GetAtomWithIdx(maps[m]).GetIsAromatic() for m in set(ids) - {retained}),
                "Aromatic bond rewrites are unsupported in v1",
                "EDIT_AROMATIC_UNSUPPORTED",
            )
            bond = rw.GetBondBetweenAtoms(*[maps[m] for m in ids])
            if name == "add_bond":
                require(bond is None, "Bond already exists", "EDIT_BOND_EXISTS")
            else:
                require(bond is not None, "Bond does not exist", "EDIT_BOND_MISSING")
                require(not bond.GetIsAromatic(), "Cannot rewrite aromatic bonds", "EDIT_AROMATIC_UNSUPPORTED")
            if fixed and name != "delete_bond":
                require(operation["order"] == 1, "Attachment must be a single bond", "EDIT_ATTACHMENT")
            if name != "add_bond":
                rw.RemoveBond(*[maps[m] for m in ids])
            if name != "delete_bond":
                order = operation["order"]
                require(type(order) is int and order in (1, 2, 3), "Unsupported bond order", "EDIT_BOND")
                rw.AddBond(
                    *[maps[m] for m in ids],
                    {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[order],
                )
        else:
            raise DesignError("Unknown edit operation", "EDIT_OPERATIONS")
        require(rw.GetNumAtoms() <= len(core) + 48, "Intermediate graph exceeds edit budget", "EDIT_SIZE")
    mol = rw.GetMol()
    mol.RemoveAllConformers()  # Inherited coordinates are NOT an edited pose.
    pending_boundary = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        if (a in core) != (b in core):
            pending_boundary.append((a if a in core else b, bond.GetBondTypeAsDouble()))
    require(
        pending_boundary == [(retained, 1.0)],
        "Product must retain exactly one declared single-bond attachment",
        "EDIT_ATTACHMENT",
    )
    # Normal implicit-H valence completion, not heavy-atom graph repair.
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() in mutable or atom.GetAtomMapNum() == retained:
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)
    try:
        mol = chemistry_copy(mol, keep_maps=True)
    except Exception as exc:
        raise DesignError("Edited graph fails chemical sanitization", "EDIT_CHEMISTRY") from exc
    require(len(Chem.GetMolFrags(mol)) == 1, "Edited molecule is disconnected", "EDIT_DISCONNECTED")
    require(
        not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()), "Radical products are unsupported", "EDIT_RADICAL"
    )
    mapping = map_to_index(mol)
    require(core <= set(mapping), "Fixed atoms were removed", "EDIT_FIXED_CORE")
    for amap in core:
        a, b = ref.GetAtomWithIdx(source_maps[amap]), mol.GetAtomWithIdx(mapping[amap])
        signature = lambda x: (
            x.GetAtomicNum(),
            x.GetFormalCharge(),
            x.GetIsotope(),
            x.GetIsAromatic(),
            x.GetChiralTag(),
        )
        require(signature(a) == signature(b), "Fixed atom identity/stereo changed", "EDIT_FIXED_CORE")
    for i, amap in enumerate(sorted(core)):
        for bmap in sorted(core)[i + 1 :]:
            a = ref.GetBondBetweenAtoms(source_maps[amap], source_maps[bmap])
            b = mol.GetBondBetweenAtoms(mapping[amap], mapping[bmap])
            signature = lambda x: None if x is None else (x.GetBondType(), x.GetStereo())
            require(signature(a) == signature(b), "Fixed bond identity/stereo changed", "EDIT_FIXED_CORE")
    outside = set(mapping) - core
    boundary = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        if (a in core) != (b in core):
            boundary.append((a if a in core else b, bond.GetBondTypeAsDouble()))
    require(
        boundary == [(retained, 1.0)] and connected_subset(mol, [mapping[m] for m in outside]),
        "Product must have one connected fragment at the declared single-bond attachment",
        "EDIT_ATTACHMENT",
    )
    old_unspecified = {s for s in _unspecified_stereo(ref) if set(s[1:]) <= core}
    require(
        not (_unspecified_stereo(mol) - old_unspecified),
        "Edit introduces unspecified tetrahedral or double-bond stereo",
        "EDIT_UNSPECIFIED_STEREO",
    )
    require(
        canonical_smiles(mol) != canonical_smiles(ref), "Proposal leaves chemical identity unchanged", "EDIT_NO_CHANGE"
    )
    baseline, measured = descriptors(ref), descriptors(mol)
    checks = []
    for constraint in grounded.spec["constraints"]:
        kind = constraint["kind"]
        values = {
            "chemical_validity": True,
            "fixed_subgraph": True,
            "single_attachment_boundary": 1,
            "formal_charge": measured["formal_charge"],
            "new_unspecified_stereocenter_count": 0,
            "new_fragment_heavy_atoms": len(outside),
        }
        value = values.get(kind)
        if kind == "descriptor_delta":
            value = measured[constraint["descriptor"]] - baseline[constraint["descriptor"]]
        if constraint["evaluation_context"] == "raw_generated_pose":
            status = "not_evaluated_requires_3d"
        else:
            require(kind in values or kind == "descriptor_delta", "Unsupported graph constraint", "EDIT_CONSTRAINT")
            status = "pass" if _compare(value, constraint["operator"], constraint["value"]) else "fail"
        checks.append(
            {
                "id": constraint["id"],
                "priority": constraint["priority"],
                "status": status,
                "measured": value,
                "requested": constraint["value"],
            }
        )
    require(
        not any(c["priority"] == "hard" and c["status"] == "fail" for c in checks),
        "Edited graph violates a hard DesignSpec constraint",
        "EDIT_HARD_CONSTRAINT",
    )
    report = {
        "status": "GRAPH_VALIDATED",
        "canonical_smiles": canonical_smiles(mol),
        "mapped_smiles": Chem.MolToSmiles(mol),
        "descriptors": measured,
        "reference_descriptors": baseline,
        "checks": checks,
        "graph_hard_constraints_pass": True,
        "all_measurable_soft_objectives_pass": all(c["status"] == "pass" for c in checks if c["priority"] == "soft"),
        "pose_validated": False,
        "synthesizability_verified": False,
        "activity_improvement_verified": False,
        "text_goal_satisfaction_verified": False,
        "rdkit_version": rdBase.rdkitVersion,
        "hydrogen_policy": "recompute_implicit_H_on_editable_atoms_and_attachment; no_heavy_graph_repair",
    }
    return mol, report
