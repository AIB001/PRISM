"""Synthetic chemistry only: no PDBbind structures or model weights in tests."""

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
pytest.importorskip("jsonschema")
from rdkit.Chem import AllChem

from prism.design.chem import canonical_smiles, heavy_copy, map_to_index
from prism.design.compile import compile_design
from prism.design.errors import DesignError
from prism.design.io import file_ref, read_json, write_json
from prism.design.schema import load_design
from prism.design.validate import PROTOCOLS, ground_design, validate_constraints, validate_design, validate_selectors
from prism.design.verify import verify_candidate


def make_case(root, growth=False):
    root.mkdir(parents=True, exist_ok=True)
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCc1ccccc1"))
    assert AllChem.EmbedMolecule(mol, randomSeed=71) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(10 + 2 * i)  # Identity is deliberately not row+1.
    ligand = root / "reference.sdf"
    with Chem.SDWriter(str(ligand)) as writer:
        writer.write(mol)
    # Bind exactly the serialized, rounded reference coordinates.
    from prism.design.chem import load_reference

    mol = load_reference(ligand)
    maps = map_to_index(mol)
    xyz = mol.GetConformer().GetPositions()
    plane_z = max(xyz[i, 2] for i in maps.values()) + 3.2
    pdb = (
        "".join(
            f"ATOM  {j:5d}  CA  ALA A{j:4d}    {xyz[i, 0]:8.3f}{xyz[i, 1]:8.3f}{plane_z:8.3f}  1.00 20.00           C  \n"
            for j, i in enumerate(maps.values(), 1)
        )
        + "END\n"
    )
    (root / "protein.pdb").write_text(pdb)
    (root / "pocket.pdb").write_text(pdb)
    (root / "original.mol2").write_text("Synthetic test placeholder; never parsed as MOL2\n")
    write_json(root / "selection_manifest.json", {"synthetic_unit_test": True})
    (root / "decisions.tsv").write_text("synthetic unit test\n")
    write_json(root / "evidence.json", {"pdb_id": "0tst", "quality_flags": []})
    write_json(
        root / "atom_map.json",
        {
            "pdb_id": "0tst",
            "atoms": [
                {"atom_map": m, "element": mol.GetAtomWithIdx(i).GetSymbol(), "xyz_angstrom": xyz[i].tolist()}
                for m, i in maps.items()
            ],
        },
    )
    edit = [] if growth else [10, 12, 14]
    core = [m for m in maps if m not in edit]
    retained, replaced = (18, None) if growth else (16, 14)
    objective = "small_growth" if growth else "local_shape"
    choice = (
        {"attachment_atom_map": retained}
        if growth
        else {
            "fixed_atom_maps": core,
            "editable_atom_maps": edit,
            "retained_attachment_atom_map": retained,
            "replaced_attachment_atom_map": replaced,
        }
    )
    write_json(
        root / "review.json",
        {
            "pdb_id": "0tst",
            "review_tier": "candidate_for_m0_pilot",
            "objective_code": objective,
            "evidence": file_ref(root / "evidence.json"),
            "selected_candidate_evidence": choice,
        },
    )
    constraints = []
    values = {
        "chemical_validity": True,
        "fixed_subgraph": True,
        "single_attachment_boundary": 1,
        "formal_charge": 0,
        "new_unspecified_stereocenter_count": 0,
        "core_rmsd": 0.5,
        "new_fragment_heavy_atoms": [2, 6],
    }
    for kind, (unit, protocol) in PROTOCOLS.items():
        constraints.append(
            {
                "id": kind,
                "kind": kind,
                "priority": "hard",
                "operator": "le" if kind == "core_rmsd" else "between" if kind == "new_fragment_heavy_atoms" else "eq",
                "value": values[kind],
                "unit": unit,
                "verifier_protocol_ref": protocol,
                "selector_ref": "core" if kind in {"core_rmsd", "fixed_subgraph"} else None,
                "descriptor": None,
                "evaluation_context": "raw_generated_pose" if kind == "core_rmsd" else "chemical_graph",
                "missing_policy": "reject_or_mark_unknown_never_pass",
                "origin": "assistant_proposed_pilot_requirement",
                "calibration_status": "pilot_policy_not_experimentally_calibrated",
                "source_evidence_refs": ["inputs.evidence"],
            }
        )
    selectors = {
        "core": {"kind": "reference_heavy_atom_maps", "atom_maps": core},
        "editable": {"kind": "reference_heavy_atom_maps", "atom_maps": edit},
    }
    if growth:
        selectors["growth_site"] = {"kind": "reference_heavy_atom_maps", "atom_maps": [retained]}
    receptor = file_ref(root / "protein.pdb")
    reference = file_ref(ligand)
    spec = {
        "schema": "prism-design-annotation/0.1",
        "design_id": "synthetic-test-only",
        "complex_id": "0tst",
        "selection_rank": 1,
        "parent_design_id": None,
        "language": "zh-CN",
        "illustrative_fixture": False,
        "stage": "grounded_draft",
        "task": "fragment_growing" if growth else "fragment_replacement",
        "prompt": "Synthetic test only",
        "objective": {
            "code": objective,
            "description_zh": "测试",
            "claim_scope": "unverified_design_requirement",
            "biological_activity_claim": False,
            "improvement_observed": None,
        },
        "inputs": {
            "receptor": receptor,
            "pocket": file_ref(root / "pocket.pdb"),
            "reference_ligand": reference,
            "atom_mapping": file_ref(root / "atom_map.json"),
            "evidence": file_ref(root / "evidence.json"),
            "original_ligand_sdf": reference,
            "original_ligand_mol2": file_ref(root / "original.mol2"),
            "coordinate_unit": "angstrom",
            "coordinate_frame_id": "synthetic_frame",
            "chemical_state_id": hashlib.sha256(canonical_smiles(mol).encode()).hexdigest(),
            "receptor_state_id": receptor["sha256"],
            "pose_id": reference["sha256"],
        },
        "selectors": selectors,
        "editing": {
            "fixed_selector": "core",
            "editable_selector": "editable",
            "growth_site_selector": "growth_site" if growth else None,
            "fixed_graph": True,
            "fixed_coordinates": False,
            "preserve_specified_stereochemistry": True,
            "heavy_atoms_only": True,
            "hydrogen_policy": "test",
            "new_fragment_heavy_atoms_range": [2, 6],
            "attachment_rules": [
                {
                    "retained_atom_map": retained,
                    "replaced_atom_map": replaced,
                    "allowed_bond_orders": [1],
                    "maximum_new_connections": 1,
                    "reference_boundary_vector_A": None
                    if growth
                    else (xyz[maps[replaced]] - xyz[maps[retained]]).tolist(),
                    "vector_is_measured_replacement_hint_not_new_pose": not growth,
                }
            ],
        },
        "constraints": constraints,
        "backend": {
            "id": None,
            "source_commit": None,
            "checkpoint_sha256": None,
            "capability_check": "not_run",
            "required_capabilities": ["fixed_subgraph"],
            "atom_indexing": "adapter_must_resolve_atom_maps_to_backend_indices",
        },
        "sampling": {
            "requested_slots": 0,
            "status": "not_run",
            "proposed_pilot_slots_after_integration": 16,
            "seed_base": 20261013,
        },
        "postprocessing": {
            "preserve_raw_output": True,
            "allow_silent_graph_repair": False,
            "allow_untracked_coordinate_change": False,
        },
        "review": {
            "review_record": file_ref(root / "review.json"),
            "review_tier": "candidate_for_m0_pilot",
            "rationale_zh": "测试",
            "case_specific_checks_zh": "测试",
            "quality_flags": [],
            "semantic_review_status": "current_assistant_reviewed_extracted_evidence",
            "expert_review_status": "not_performed",
            "assay_metadata_status": "local_companion_csv_unverified",
            "design_feasibility_status": "not_tested",
        },
        "provenance": {
            "selection_manifest": dict(file_ref(root / "selection_manifest.json"), path="selection_manifest.json"),
            "manual_decisions": file_ref(root / "decisions.tsv"),
            "label_origin": "assistant_proposed_design_plus_program_measured_evidence",
            "requested_model": "gpt-6-astra",
            "runtime_model_id": "not_exposed_to_artifact_builder",
            "separate_model_api_invoked": False,
            "rdkit_version": "test",
        },
        "use_policy": {
            "executable_in_current_prism": False,
            "optimization_pair_available": False,
            "affinity_training_eligible": False,
            "reconstruction_mask_candidate": True,
            "train_validation_test_split": None,
            "leakage_review_required": True,
        },
        "execution_blockers": ["annotation remains a draft; compiler creates a separate execution artifact"],
    }
    path = root / "complexes" / "0tst" / "designspec.json"
    write_json(path, spec)
    return path


@pytest.fixture
def case(tmp_path):
    return make_case(tmp_path / "annotation")


@pytest.fixture
def grounded(case):
    return ground_design(load_design(case), for_execution=True)


def changed_candidate(grounded):
    mol = heavy_copy(grounded.molecule, keep_maps=False)
    atom = mol.GetAtomWithIdx(0)
    atom.SetAtomicNum(7)
    atom.SetNumExplicitHs(0)
    atom.SetNoImplicit(False)
    Chem.SanitizeMol(mol)
    return mol


def checks(result):
    return {c["id"]: c for c in result["checks"]}


def test_binding_masks_compile_preserves_frame_and_not_annotation(case, grounded, tmp_path):
    before = case.read_bytes()
    compiled = compile_design(grounded, tmp_path / "compiled")
    assert load_design(case).root == case.parents[2]
    assert compiled["fixed_indices_0based"] == [3, 4, 5, 6, 7, 8]
    assert compiled["editable_indices_0based"] == [0, 1, 2]
    staged = Chem.MolFromMolFile(compiled["inputs"]["reference_ligand"]["path"], removeHs=False)
    assert np.allclose(
        staged.GetConformer().GetPositions(), heavy_copy(grounded.molecule).GetConformer().GetPositions()
    )
    assert all(a.GetAtomMapNum() == 0 for a in staged.GetAtoms())
    assert canonical_smiles(staged) == canonical_smiles(grounded.molecule)
    assert not Chem.FindMolChiralCenters(staged, includeUnassigned=True)
    assert case.read_bytes() == before


@pytest.mark.parametrize("field", ["pocket", "reference_ligand", "atom_mapping", "evidence", "receptor"])
def test_changed_bound_input_rejected(case, field):
    path = Path(read_json(case)["inputs"][field]["path"])
    path.write_bytes(path.read_bytes() + b"\n")
    assert validate_design(case, for_execution=True)["error_code"] == "INPUT_HASH_MISMATCH"


def test_changed_review_decision_rejected(case):
    data = read_json(case)
    data["editing"]["attachment_rules"][0]["retained_atom_map"] = 18
    write_json(case, data)
    assert validate_design(case, for_execution=True)["error_code"] == "REVIEW_BINDING"


@pytest.mark.parametrize("value, code", [(' {"a": 1, "a": 2}', "DUPLICATE_JSON_KEY"), ('{"x":NaN}', "NONFINITE_JSON")])
def test_strict_json(tmp_path, value, code):
    path = tmp_path / "bad.json"
    path.write_text(value)
    with pytest.raises(DesignError) as error:
        read_json(path)
    assert error.value.code == code


def test_missing_reference_is_rejected_without_keyerror(case):
    data = read_json(case)
    data["inputs"]["reference_ligand"] = None
    write_json(case, data)
    assert validate_design(case)["status"] == "REJECTED"


def test_schema_unknown_field_and_no_boolean_number_confusion(case):
    data = read_json(case)
    data["unrecognized"] = 1
    write_json(case, data)
    assert validate_design(case)["error_code"] == "SCHEMA_VALIDATION_FAILED"


def test_compile_requires_execution_review(case, tmp_path):
    with pytest.raises(DesignError, match="for_execution"):
        compile_design(ground_design(load_design(case)), tmp_path / "out")
    data = read_json(case)
    review = read_json(data["review"]["review_record"]["path"])
    review["review_tier"] = data["review"]["review_tier"] = "needs_review"
    write_json(data["review"]["review_record"]["path"], review)
    data["review"]["review_record"] = file_ref(data["review"]["review_record"]["path"])
    write_json(case, data)
    assert validate_design(case, for_execution=True)["error_code"] == "REVIEW_REQUIRED"


def test_source_flags_cannot_be_cleared_from_spec(case):
    data = read_json(case)
    data["review"]["quality_flags"] = ["STEREO_UNRESOLVED"]
    write_json(case, data)
    assert validate_design(case, for_execution=True)["error_code"] == "EVIDENCE_BINDING"


def test_masks_and_constraint_units_fail_closed(grounded):
    spec = copy.deepcopy(grounded.spec)
    spec["selectors"]["editable"]["atom_maps"].append(16)
    with pytest.raises(DesignError, match="overlap"):
        validate_selectors(spec, grounded.molecule)
    spec = copy.deepcopy(grounded.spec)
    spec["constraints"][5]["unit"] = "dimensionless"
    with pytest.raises(DesignError, match="unit/protocol"):
        validate_constraints(spec, grounded.molecule)
    spec["constraints"] = spec["constraints"][:-1]
    with pytest.raises(DesignError):
        validate_constraints(spec, grounded.molecule)


def test_candidate_atom_reordering_and_reference_is_not_a_new_design(grounded):
    mol = changed_candidate(grounded)
    result = verify_candidate(grounded, mol)
    assert result["selection_pass"]
    reordered = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    r = verify_candidate(grounded, reordered)
    assert r["selection_pass"] and r["core_mapping"]["core_rmsd_A"] == pytest.approx(0)
    unchanged = verify_candidate(grounded, grounded.molecule)
    assert unchanged["design_hard_pass"] and not unchanged["selection_pass"]


def test_raw_pose_displacement_cannot_be_hidden_by_alignment(grounded):
    mol = changed_candidate(grounded)
    conf = mol.GetConformer()
    for i, xyz in enumerate(conf.GetPositions()):
        conf.SetAtomPosition(i, xyz + np.array([0, 0, 6.0]))
    result = verify_candidate(grounded, mol)
    assert checks(result)["core_rmsd"]["status"] == "fail"
    assert result["core_mapping"]["core_rmsd_A"] == pytest.approx(6.0)
    assert not result["ligand_superposition"]


def test_charge_and_core_graph_rejected(grounded):
    mol = changed_candidate(grounded)
    mol.GetAtomWithIdx(0).SetFormalCharge(1)
    assert checks(verify_candidate(grounded, mol))["formal_charge"]["status"] == "fail"
    mol = changed_candidate(grounded)
    mol.GetAtomWithIdx(6).SetAtomicNum(7)
    assert checks(verify_candidate(grounded, mol))["fixed_subgraph"]["status"] == "fail"


def test_multiple_components_and_missing_pose_are_not_passes(grounded):
    multi = Chem.CombineMols(changed_candidate(grounded), Chem.MolFromSmiles("O"))
    assert not verify_candidate(grounded, multi)["design_hard_pass"]
    mol = changed_candidate(grounded)
    mol.RemoveAllConformers()
    result = verify_candidate(grounded, mol)
    assert checks(result)["core_rmsd"]["status"] == "unknown"
    assert not result["selection_pass"]


def test_soft_objective_measured_separately(grounded):
    constraint = copy.deepcopy(grounded.spec["constraints"][0])
    constraint.update(
        id="lower_logp",
        kind="descriptor_delta",
        descriptor="crippen_logp",
        priority="soft",
        operator="le",
        value=-0.5,
        unit="dimensionless",
        verifier_protocol_ref="rdkit-crippen-delta/1",
    )
    grounded.spec["constraints"].append(constraint)
    r = verify_candidate(grounded, changed_candidate(grounded))
    assert r["design_hard_pass"] and r["all_soft_objectives_pass"]
    assert checks(r)["lower_logp"]["measured"] < -0.5
    constraint["value"] = -100
    r = verify_candidate(grounded, changed_candidate(grounded))
    assert r["design_hard_pass"] and not r["selection_pass"]


def grow_candidate(grounded, attachment=4):
    source = heavy_copy(grounded.molecule, keep_maps=False)
    rw = Chem.RWMol(source)
    a, b = rw.AddAtom(Chem.Atom(6)), rw.AddAtom(Chem.Atom(6))
    rw.AddBond(attachment, a, Chem.BondType.SINGLE)
    rw.AddBond(a, b, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    conf = mol.GetConformer()
    # New coordinates affect only the pose screen, not the graph/attachment assertions.
    origin = conf.GetAtomPosition(attachment)
    conf.SetAtomPosition(a, (origin.x, origin.y, origin.z - 1.5))
    conf.SetAtomPosition(b, (origin.x + 1.4, origin.y, origin.z - 2.0))
    return mol


def test_growth_compiles_and_wrong_site_is_rejected(tmp_path):
    g = ground_design(load_design(make_case(tmp_path / "annotation", growth=True)), for_execution=True)
    compiled = compile_design(g, tmp_path / "compiled")
    assert len(compiled["fixed_indices_0based"]) == 9 and not compiled["editable_indices_0based"]
    assert compiled["growth_prior"]["attachment_atom_map"] == 18
    valid = verify_candidate(g, grow_candidate(g))
    assert valid["design_hard_pass"]
    invalid = verify_candidate(g, grow_candidate(g, attachment=0))
    assert not invalid["design_hard_pass"]
    assert checks(invalid)["single_attachment_boundary"]["status"] == "fail"


def test_new_unspecified_center_at_preexisting_growth_atom_is_counted(tmp_path):
    g = ground_design(load_design(make_case(tmp_path / "annotation", growth=True)), for_execution=True)
    # Add an ethyl group to the existing propyl CH2, creating a new stereocenter.
    g.spec["editing"]["attachment_rules"][0]["retained_atom_map"] = 12
    mol = grow_candidate(g, attachment=1)
    result = verify_candidate(g, mol)
    assert checks(result)["new_unspecified_stereocenter_count"]["measured"] == 1
    assert not result["selection_pass"]


def test_slot_accounting_failed_missing_duplicate_and_md_handoff(grounded, tmp_path):
    from prism.design.runner import evaluate_slots, summarize_slots
    from prism.generation.handoff import select_for_md

    output = tmp_path / "out"
    backend = output / "backend"
    backend.mkdir(parents=True)
    molecule = changed_candidate(grounded)
    events = []
    for slot in [1, 2]:
        path = backend / f"raw{slot}.sdf"
        with Chem.SDWriter(str(path)) as writer:
            writer.write(molecule)
        events.append({"slot": slot, "status": "DECODED", "sdf": str(path)})
    events.append({"slot": 3, "status": "DECODE_FAILED", "sdf": None})
    (backend / "attempts.jsonl").write_text("".join(json.dumps(e) + "\n" for e in events))
    records, candidates = evaluate_slots(grounded, backend, output, 4)
    counts = summarize_slots(records)
    assert counts["raw_slots"] == 3 and counts["decoded"] == 2 and counts["decode_failed"] == 1
    assert counts["not_attempted"] == 1 and counts["unique_selected"] == 1
    assert counts["duplicate_selected_chemistry"] == 1
    assert records[-1]["validation"] is None
    assert candidates[0]["hydrogenated_path"]
    write_json(output / "manifest.json", {"models": {"flowr_root": {"candidates": candidates}}})
    selected, skipped = select_for_md(output)
    assert len(selected) == 1 and not skipped
    assert selected[0].hydrogen_count > 0


def test_budget_and_invalid_ledger_are_not_silently_truncated(tmp_path):
    from prism.design.runner import _attempts

    (tmp_path / "attempts.jsonl").write_text('{"slot":2,"status":"DECODE_FAILED"}\n')
    with pytest.raises(DesignError, match="ledger"):
        _attempts(tmp_path, 1)


def test_cli_routing_is_lazy_and_validate_works(case):
    result = subprocess.run(
        [sys.executable, "-m", "prism.design", "validate", "--spec", str(case), "--for-execution"],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["eligible_for_compilation"]
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from prism.cli import main; import sys; main(['design','--help']); assert 'prism.builder' not in sys.modules",
        ],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0 and "prism design" in result.stdout


def make_backend_config(tmp_path):
    from prism.design.backends import SUPPORTED_COMMIT

    source = tmp_path / "flowr"
    source.mkdir(exist_ok=True)
    python = tmp_path / "venv" / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.symlink_to(sys.executable)
    checkpoint = tmp_path / "weights.ckpt"
    checkpoint.write_bytes(b"synthetic-not-a-model")
    data = {
        "backend": "flowr_root",
        "python": str(python),
        "source": str(source),
        "source_commit": SUPPORTED_COMMIT,
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": file_ref(checkpoint)["sha256"],
    }
    config = tmp_path / "backend.json"
    write_json(config, data)
    return config


def test_backend_pin_changed_weights_and_venv_symlink(tmp_path, monkeypatch):
    from prism.design import backends

    config_path = make_backend_config(tmp_path)
    config = backends.load_backend(config_path)
    assert Path(config["python"]).is_symlink()
    monkeypatch.setattr(
        backends.subprocess,
        "check_output",
        lambda command, **kwargs: backends.SUPPORTED_COMMIT if "rev-parse" in command else "",
    )
    assert backends.check_backend(config)["status"] == "VERIFIED_FILES"
    Path(config["checkpoint"]).write_bytes(b"changed")
    with pytest.raises(DesignError) as error:
        backends.check_backend(config)
    assert error.value.code == "CHECKPOINT_HASH_MISMATCH"
    data = read_json(config_path)
    data["source_commit"] = "0" * 40
    write_json(config_path, data)
    with pytest.raises(DesignError) as error:
        backends.load_backend(config_path)
    assert error.value.code == "BACKEND_REVISION"


def test_backend_does_not_fallback_to_unconditional_generation(tmp_path):
    from prism.design.backends import load_backend

    path = make_backend_config(tmp_path)
    data = read_json(path)
    data["backend"] = "flowr"
    write_json(path, data)
    with pytest.raises(DesignError) as error:
        load_backend(path)
    assert error.value.code == "BACKEND_UNSUPPORTED"


def test_backend_timeout_keeps_process_record(tmp_path, monkeypatch):
    from prism.design import backends

    config = backends.load_backend(make_backend_config(tmp_path))
    monkeypatch.setattr(backends, "check_backend", lambda config: {"test_only": True})

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], 1)

    monkeypatch.setattr(backends.subprocess, "run", timeout)
    out = tmp_path / "backend_run"
    result = backends.execute_backend(config, tmp_path / "compiled.json", out, 2, 1)
    assert result["status"] == "BACKEND_TIMEOUT"
    assert read_json(out / "process.json")["return_code"] is None


def test_failed_backend_cannot_be_reported_as_completed(case, tmp_path, monkeypatch):
    from prism.design import runner

    config_path = make_backend_config(tmp_path)
    monkeypatch.setattr(runner, "execute_backend", lambda *args: {"status": "BACKEND_TIMEOUT", "return_code": None})
    result = runner.run_design(case, config_path, tmp_path / "run", samples=3)
    assert result["status"] == "FAILED"
    assert result["counts"]["requested"] == 3 and result["counts"]["raw_slots"] == 0
    assert result["counts"]["not_attempted"] == 3
    manifest = read_json(tmp_path / "run" / "manifest.json")
    assert manifest["models"]["flowr_root"]["candidates"] == []
    with pytest.raises(DesignError) as error:
        runner.run_design(case, config_path, tmp_path / "run", samples=3)
    assert error.value.code == "OUTPUT_EXISTS"


def test_fixed_stereochemistry_is_verified_in_independent_graph_match(grounded):
    # Keep a complete tetrahedral definition within the core, away from the cut.
    reference = Chem.MolFromSmiles("CCC[C@H](F)COc1ccccc1")
    AllChem.EmbedMolecule(reference, randomSeed=12)
    for atom in reference.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    g = copy.deepcopy(grounded)
    g.molecule = reference
    g.spec["selectors"]["core"]["atom_maps"] = list(range(3, reference.GetNumAtoms() + 1))
    g.spec["selectors"]["editable"]["atom_maps"] = [1, 2]
    g.spec["editing"]["attachment_rules"][0].update(retained_atom_map=3, replaced_atom_map=2)
    good = verify_candidate(g, reference)
    assert checks(good)["fixed_subgraph"]["status"] == "pass"
    flipped = Chem.Mol(reference)
    flipped.GetAtomWithIdx(3).InvertChirality()
    bad = verify_candidate(g, flipped)
    assert checks(bad)["fixed_subgraph"]["status"] == "fail"


def test_pocket_from_another_site_in_same_receptor_is_rejected(tmp_path):
    from prism.design.geometry import pocket_review

    def line(i, xyz):
        return f"ATOM  {i:5d}  CA  ALA A{i:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00           C  \n"

    nearby, distant = line(1, (0, 0, 3.2)), line(2, (50, 50, 50))
    protein = tmp_path / "protein.pdb"
    protein.write_text(nearby + distant)
    pocket = tmp_path / "pocket.pdb"
    pocket.write_text(distant)
    with pytest.raises(DesignError) as error:
        pocket_review(protein, pocket, np.array([[0.0, 0.0, 0.0]]))
    assert error.value.code == "POCKET_OVERLAP"


def test_wrong_constraint_context_is_rejected(grounded):
    spec = copy.deepcopy(grounded.spec)
    spec["constraints"][1]["selector_ref"] = "editable"
    with pytest.raises(DesignError) as error:
        validate_constraints(spec, grounded.molecule)
    assert error.value.code == "PROTOCOL_UNSUPPORTED"


def test_m0_qc_warning_is_withheld_without_losing_design_result(grounded, tmp_path, monkeypatch):
    from prism.design.runner import evaluate_slots, summarize_slots
    from prism.generation import quality

    backend = tmp_path / "backend"
    backend.mkdir()
    path = backend / "raw.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(changed_candidate(grounded))
    (backend / "attempts.jsonl").write_text(json.dumps({"slot": 1, "status": "DECODED", "sdf": str(path)}) + "\n")
    monkeypatch.setattr(quality, "assess", lambda *args: quality.QualityReport(status="WARN_GEOMETRY"))
    records, candidates = evaluate_slots(grounded, backend, tmp_path, 1)
    assert records[0]["validation"]["selection_pass"]
    assert not candidates and summarize_slots(records)["unique_selected"] == 0
