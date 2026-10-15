"""Planning boundaries, chemical execution and a real local HTTP transport.

Provider responses here are fixtures, never evidence of LLM design quality.
"""

import copy
import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

Chem = pytest.importorskip("rdkit.Chem")
pytest.importorskip("jsonschema")
from test_design import make_case
from prism.design.chem import canonical_smiles, map_to_index
from prism.design.cli import main
from prism.design.errors import DesignError
from prism.design.edits import execute_edits
from prism.design.io import file_ref, read_json, sha256, write_json
from prism.design.llm import load_llm_config
from prism.design.planner import (
    check_plan,
    load_planning_request,
    parse_proposal,
    plan_design,
    prepare_planning,
    validate_proposal,
)


@pytest.fixture
def case(tmp_path):
    spec = make_case(tmp_path / "source")
    request = prepare_planning(
        spec, tmp_path / "request", "Explore a polar terminal group; preserve the aromatic core."
    )
    return Path(request["request"])


def proposal_for(request):
    data = read_json(request)
    return {
        "schema": "prism-design-edit-proposal/1",
        "context_sha256": data["context_sha256"],
        "decision": "propose",
        "summary": "Test a terminal oxygen substitution as a polarity hypothesis.",
        "evidence_ids": ["atom:10", "bond:10:12"],
        "operations": [{"op": "mutate_atom", "atom_map": 10, "element": "O", "formal_charge": 0}],
        "uncertainties": ["No synthesis route or binding improvement was measured."],
        "unfulfilled_goals": [],
    }


def test_prepare_is_grounded_and_replayable(case):
    grounded, context, messages = load_planning_request(case)
    assert context["selectors"]["editable"]["atom_maps"] == [10, 12, 14]
    assert "[CH3:10]" in context["reference_mapped_smiles"]
    assert {e["kind"] for e in context["evidence"]} >= {"reference_distance_contact", "measured_reference_atom"}
    assert len(messages) == 2
    assert str(grounded.document.root) not in messages[1]["content"]
    before = sha256(grounded.document.path)
    molecule, report = validate_proposal(grounded, context, proposal_for(case))
    assert canonical_smiles(molecule) == "OCCc1ccccc1"
    assert molecule.GetNumConformers() == 0
    assert report["pose_validated"] is False
    assert any(c["status"] == "not_evaluated_requires_3d" for c in report["checks"])
    assert before == sha256(grounded.document.path)


@pytest.mark.parametrize(
    "mutation,code",
    [
        (lambda p: p.update(context_sha256="0" * 64), "PLANNER_BINDING"),
        (lambda p: p.update(evidence_ids=["invented_hydrogen_bond"]), "PLANNER_EVIDENCE"),
        (lambda p: p.update(evidence_ids=["ligand"]), "PLANNER_EVIDENCE"),
        (lambda p: p["operations"][0].update(atom_map=999), "PLANNER_EVIDENCE"),
        (lambda p: p["operations"][0].update(atom_map=16), "EDIT_FIXED_CORE"),
        (lambda p: p["operations"][0].update(element="F", atom_map=12), "EDIT_CHEMISTRY"),
        (lambda p: p["operations"][0].update(formal_charge=-1), "EDIT_HARD_CONSTRAINT"),
        (lambda p: p["operations"][0].update(element="C"), "EDIT_NO_CHANGE"),
        (lambda p: p["operations"][0].update(atom_map=True), "PLANNER_SCHEMA"),
        (lambda p: p["operations"][0].update(atom_map=10.0), "ATOM_MAP_INVALID"),
        (lambda p: p.update(shell="anything"), "PLANNER_SCHEMA"),
    ],
)
def test_rejects_invalid_proposals(case, mutation, code):
    grounded, context, _ = load_planning_request(case)
    proposal = proposal_for(case)
    mutation(proposal)
    # Preserve an actual evidence citation for the core-mutation test.
    if proposal["operations"][0]["atom_map"] == 16:
        proposal["evidence_ids"] = ["atom:16"]
    with pytest.raises(DesignError) as exc:
        validate_proposal(grounded, context, proposal)
    assert exc.value.code == code


def test_add_delete_uses_atom_maps_not_shifted_rows(case):
    grounded, context, _ = load_planning_request(case)
    proposal = proposal_for(case)
    proposal["operations"] = [
        {"op": "delete_atom", "atom_map": 10},
        {"op": "add_atom", "atom_map": 100, "element": "O", "formal_charge": 0},
        {"op": "add_bond", "atom_maps": [100, 12], "order": 1},
    ]
    mol, report = validate_proposal(grounded, context, proposal)
    assert canonical_smiles(mol) == "OCCc1ccccc1"
    assert 10 not in map_to_index(mol) and 100 in map_to_index(mol)
    assert report["graph_hard_constraints_pass"]


@pytest.mark.parametrize(
    "operations,code",
    [
        ([{"op": "delete_bond", "atom_maps": [10, 12]}], "EDIT_DISCONNECTED"),
        ([{"op": "add_bond", "atom_maps": [10, 18], "order": 1}], "EDIT_FIXED_CORE"),
        ([{"op": "add_bond", "atom_maps": [10, 16], "order": 1}], "EDIT_ATTACHMENT"),
        ([{"op": "delete_atom", "atom_map": 10}, {"op": "delete_atom", "atom_map": 12}], "EDIT_HARD_CONSTRAINT"),
        ([{"op": "add_bond", "atom_maps": [10, 12], "order": 1}], "EDIT_BOND_EXISTS"),
        ([{"op": "delete_bond", "atom_maps": [10, 14]}], "EDIT_BOND_MISSING"),
        ([{"op": "change_bond", "atom_maps": [12, 14], "order": 2}], "EDIT_UNSPECIFIED_STEREO"),
        (
            [
                {"op": "add_atom", "atom_map": 10, "element": "O", "formal_charge": 0},
                {"op": "add_bond", "atom_maps": [10, 12], "order": 1},
            ],
            "EDIT_NEW_MAP",
        ),
    ],
)
def test_graph_failures(case, operations, code):
    grounded, context, _ = load_planning_request(case)
    proposal = proposal_for(case)
    proposal["operations"] = operations
    with pytest.raises(DesignError) as exc:
        validate_proposal(grounded, context, proposal)
    assert exc.value.code == code


def test_valid_bond_order_change(case):
    grounded, context, _ = load_planning_request(case)
    proposal = proposal_for(case)
    proposal["operations"] = [{"op": "change_bond", "atom_maps": [10, 12], "order": 2}]
    mol, _ = validate_proposal(grounded, context, proposal)
    assert canonical_smiles(mol) == "C=CCc1ccccc1"


def test_new_tetrahedral_stereo_is_rejected(case):
    grounded, context, _ = load_planning_request(case)
    proposal = proposal_for(case)
    proposal["operations"] = [
        {"op": "add_atom", "atom_map": 100, "element": "F", "formal_charge": 0},
        {"op": "add_bond", "atom_maps": [12, 100], "order": 1},
    ]
    with pytest.raises(DesignError) as exc:
        validate_proposal(grounded, context, proposal)
    assert exc.value.code == "EDIT_UNSPECIFIED_STEREO"


def test_remote_fixed_stereo_preserved_and_defining_neighbours_protected(case):
    grounded, _, _ = load_planning_request(case)
    mol = Chem.MolFromSmiles("CCCc1ccc([C@H](F)Cl)cc1")
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(10 + 2 * i)
    grounded.molecule = mol
    grounded.spec["selectors"]["core"]["atom_maps"] = [a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetIdx() >= 3]
    edited, _ = execute_edits(grounded, [{"op": "mutate_atom", "atom_map": 10, "element": "O", "formal_charge": 0}])
    assert canonical_smiles(edited) == canonical_smiles(Chem.MolFromSmiles("OCCc1ccc([C@H](F)Cl)cc1"))
    fluorine_map = next(a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetSymbol() == "F")
    with pytest.raises(DesignError) as exc:
        execute_edits(grounded, [{"op": "mutate_atom", "atom_map": fluorine_map, "element": "Cl", "formal_charge": 0}])
    assert exc.value.code == "EDIT_STEREO_UNSUPPORTED"


def test_soft_failure_is_reported_without_claiming_goal_success(case):
    grounded, context, _ = load_planning_request(case)
    constraint = copy.deepcopy(grounded.spec["constraints"][0])
    constraint.update(
        id="impossible_logp_goal",
        kind="descriptor_delta",
        descriptor="crippen_logp",
        priority="soft",
        operator="le",
        value=-100,
    )
    grounded.spec["constraints"].append(constraint)
    _, report = validate_proposal(grounded, context, proposal_for(case))
    assert report["graph_hard_constraints_pass"]
    assert not report["all_measurable_soft_objectives_pass"]
    assert not report["text_goal_satisfaction_verified"]


def test_aromatic_editable_rewrite_is_explicitly_unsupported(case):
    grounded, _, _ = load_planning_request(case)
    # Unit-level alternate scope only; this is not a reviewed execution request.
    grounded.spec["selectors"]["core"]["atom_maps"].remove(18)
    grounded.spec["selectors"]["editable"]["atom_maps"].append(18)
    with pytest.raises(DesignError) as exc:
        execute_edits(grounded, [{"op": "mutate_atom", "atom_map": 18, "element": "N", "formal_charge": 0}])
    assert exc.value.code == "EDIT_AROMATIC_UNSUPPORTED"


def test_growth_at_reviewed_attachment(tmp_path):
    spec = make_case(tmp_path / "source", growth=True)
    prepared = prepare_planning(spec, tmp_path / "request")
    request = Path(prepared["request"])
    grounded, context, _ = load_planning_request(request)
    proposal = proposal_for(request)
    proposal["evidence_ids"] = ["atom:18"]
    proposal["operations"] = [
        {"op": "add_atom", "atom_map": 100, "element": "C", "formal_charge": 0},
        {"op": "add_atom", "atom_map": 101, "element": "O", "formal_charge": 0},
        {"op": "add_bond", "atom_maps": [18, 100], "order": 1},
        {"op": "add_bond", "atom_maps": [100, 101], "order": 1},
    ]
    mol, _ = validate_proposal(grounded, context, proposal)
    assert mol.GetNumHeavyAtoms() == 11


def test_imported_response_and_bound_artifacts(case, tmp_path):
    response = tmp_path / "response.json"
    write_json(response, proposal_for(case))
    result = plan_design(case, tmp_path / "plan", response_path=response)
    assert result["status"] == "GRAPH_VALIDATED"
    assert check_plan(result["plan"])["graph_validation"]["pose_validated"] is False
    plan = read_json(result["plan"])
    assert plan["origin"]["mode"] == "imported_response"
    assert read_json(tmp_path / "plan" / "attempts.json")[0]["metadata"]["separate_model_api_invoked"] is False
    (tmp_path / "plan" / "edited_graph.smi").write_text("CC\n")
    with pytest.raises(DesignError, match="SHA-256"):
        check_plan(result["plan"])


def test_abstention_is_not_a_successful_edit(case, tmp_path):
    proposal = proposal_for(case)
    proposal.update(decision="abstain", operations=[], evidence_ids=[])
    response = tmp_path / "response.json"
    write_json(response, proposal)
    result = plan_design(case, tmp_path / "plan", response_path=response)
    assert result["status"] == "ABSTAINED"
    assert read_json(result["plan"])["graph_preview"] is None
    assert check_plan(result["plan"])["status"] == "ABSTAINED"


def test_failure_is_saved_and_not_repaired_silently(case, tmp_path):
    response = tmp_path / "response.txt"
    response.write_text('{"operations": [], "operations": []}')
    result = plan_design(case, tmp_path / "plan", response_path=response)
    assert result["status"] == "REJECTED"
    assert result["last_error"]["error_code"] == "DUPLICATE_JSON_KEY"
    assert not (tmp_path / "plan" / "edited_graph.smi").exists()


def test_tampered_context_rejected_even_with_updated_file_hash(case):
    manifest = read_json(case)
    context_path = Path(manifest["context"]["path"])
    context = read_json(context_path)
    context["reference_smiles"] = "CC"
    write_json(context_path, context)
    manifest["context"] = file_ref(context_path)
    write_json(case, manifest)
    with pytest.raises(DesignError) as exc:
        load_planning_request(case)
    assert exc.value.code == "PLANNER_BINDING"


def test_reject_changed_source_before_model_call(case, tmp_path):
    spec = read_json(case)["base_design"]["path"]
    with open(spec, "a") as handle:
        handle.write("\n")
    with pytest.raises(DesignError) as exc:
        plan_design(case, tmp_path / "plan", response_path=tmp_path / "missing.json")
    assert exc.value.code == "INPUT_HASH_MISMATCH"


def test_candidate_must_match_exact_planned_graph(case, tmp_path):
    grounded, _, _ = load_planning_request(case)
    response = tmp_path / "response.json"
    write_json(response, proposal_for(case))
    result = plan_design(case, tmp_path / "plan", response_path=response)
    with pytest.raises(DesignError) as exc:
        check_plan(result["plan"], grounded.files["reference_ligand"])
    assert exc.value.code == "PLAN_GRAPH_MISMATCH"


def test_matching_candidate_uses_existing_pose_verifier(case, tmp_path):
    from test_design import changed_candidate

    grounded, _, _ = load_planning_request(case)
    response = tmp_path / "response.json"
    proposal = proposal_for(case)
    proposal["operations"][0]["element"] = "N"
    write_json(response, proposal)
    result = plan_design(case, tmp_path / "plan", response_path=response)
    candidate = changed_candidate(grounded)
    path = tmp_path / "posed.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(candidate)
    good = check_plan(result["plan"], path)
    assert good["status"] == "PASS"
    assert good["pose_validation"]["design_hard_pass"]
    for i, pos in enumerate(candidate.GetConformer().GetPositions()):
        candidate.GetConformer().SetAtomPosition(i, pos + [0, 0, 6])
    with Chem.SDWriter(str(path)) as writer:
        writer.write(candidate)
    bad = check_plan(result["plan"], path)
    assert bad["status"] == "REJECTED"
    assert not bad["pose_validation"]["design_hard_pass"]


@pytest.mark.parametrize(
    "text,code",
    [
        ('{"number": NaN}', "NONFINITE_JSON"),
        ('{"number": 1e999}', "NONFINITE_JSON"),
        ('{"x": 1, "x": 2}', "DUPLICATE_JSON_KEY"),
        ("prefix {} suffix", "INVALID_JSON"),
    ],
)
def test_strict_response_parser(text, code):
    with pytest.raises(DesignError) as exc:
        parse_proposal(text)
    assert exc.value.code == code


def test_fenced_json(case):
    assert parse_proposal("```json\n" + json.dumps(proposal_for(case)) + "\n```")["decision"] == "propose"


@pytest.fixture
def http_server():
    calls, responses = [], []

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            calls.append(
                {
                    "path": self.path,
                    "body": json.loads(self.rfile.read(int(self.headers["Content-Length"]))),
                    "authorization": self.headers.get("Authorization"),
                }
            )
            code, data = responses.pop(0)
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps(data).encode())

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_port}/v1", calls, responses
    server.shutdown()
    server.server_close()
    thread.join()


def _config(tmp_path, url, **kwargs):
    path = tmp_path / "llm.json"
    write_json(
        path,
        dict(
            {
                "schema": "prism-design-planner-config/1",
                "base_url": url,
                "model": "test-fixture-not-a-real-model",
                "api_key_env": None,
            },
            **kwargs,
        ),
    )
    return path


def _completion(proposal, finish="stop"):
    return {
        "id": "fixture",
        "model": "fixture",
        "choices": [{"finish_reason": finish, "message": {"role": "assistant", "content": json.dumps(proposal)}}],
        "usage": {"total_tokens": 20},
    }


def test_live_http_protocol_and_bounded_validation_retry(case, tmp_path, http_server, monkeypatch):
    url, calls, responses = http_server
    monkeypatch.setenv("PRISM_TEST_LLM_KEY", "test-secret-do-not-log")
    config = _config(tmp_path, url, max_attempts=2, api_key_env="PRISM_TEST_LLM_KEY")
    bad = proposal_for(case)
    bad["operations"][0]["element"] = "C"
    responses.extend([(200, _completion(bad)), (200, _completion(proposal_for(case)))])
    result = plan_design(case, tmp_path / "plan", llm_config=config)
    assert result["status"] == "GRAPH_VALIDATED" and result["attempts"] == 2
    assert len(calls) == 2 and calls[0]["path"] == "/v1/chat/completions"
    assert calls[0]["authorization"] == "Bearer test-secret-do-not-log"
    assert calls[0]["body"]["response_format"] == {"type": "json_object"}
    assert "EDIT_NO_CHANGE" in calls[1]["body"]["messages"][-1]["content"]
    for path in (tmp_path / "plan").iterdir():
        assert "test-secret-do-not-log" not in path.read_text()
    assert check_plan(result["plan"])["status"] == "GRAPH_VALIDATED"


def test_validation_retry_budget_is_exhausted_and_preserved(case, tmp_path, http_server):
    url, calls, responses = http_server
    config = _config(tmp_path, url, max_attempts=2)
    bad = proposal_for(case)
    bad["operations"][0]["element"] = "C"
    responses.extend([(200, _completion(bad)), (200, _completion(bad))])
    result = plan_design(case, tmp_path / "plan", llm_config=config)
    assert result["status"] == "REJECTED" and len(calls) == 2
    ledger = read_json(tmp_path / "plan" / "attempts.json")
    assert all(e["error_code"] == "EDIT_NO_CHANGE" for e in ledger)
    assert all(e["response"] for e in ledger)


def test_redirects_do_not_forward_authorization(tmp_path, monkeypatch):
    from prism.design.llm import _NoRedirect
    from urllib.request import Request

    with pytest.raises(DesignError) as exc:
        _NoRedirect().redirect_request(
            Request("https://example.com", headers={"Authorization": "Bearer test"}),
            None,
            307,
            "redirect",
            {},
            "https://another.example.com",
        )
    assert exc.value.code == "LLM_REDIRECT"


def test_source_change_during_model_call_rejects_result(case, tmp_path, monkeypatch):
    spec = read_json(case)["base_design"]["path"]
    content = json.dumps(proposal_for(case))

    def fake_call(config, messages):
        with open(spec, "a") as handle:
            handle.write("\n")
        return content, {"separate_model_api_invoked": False}

    monkeypatch.setattr("prism.design.planner.chat_completion", fake_call)
    config = _config(tmp_path, "http://127.0.0.1:8000/v1", max_attempts=3)
    result = plan_design(case, tmp_path / "plan", llm_config=config)
    assert result["status"] == "REJECTED"
    assert result["attempts"] == 1
    assert result["last_error"]["error_code"] == "INPUT_HASH_MISMATCH"


def test_http_failure_no_retry_and_no_provider_body_leak(case, tmp_path, http_server):
    url, calls, responses = http_server
    config = _config(tmp_path, url, max_attempts=3)
    responses.append((401, {"error": "sensitive-provider-error-body"}))
    result = plan_design(case, tmp_path / "plan", llm_config=config)
    assert result["status"] == "REJECTED" and len(calls) == 1
    assert "sensitive-provider-error-body" not in (tmp_path / "plan" / "attempts.json").read_text()


@pytest.mark.parametrize("finish", ["length", "tool_calls", "content_filter"])
def test_incomplete_response_rejected(case, tmp_path, http_server, finish):
    url, calls, responses = http_server
    config = _config(tmp_path, url)
    responses.append((200, _completion(proposal_for(case), finish)))
    result = plan_design(case, tmp_path / "plan", llm_config=config)
    assert result["last_error"]["error_code"] == "LLM_INCOMPLETE"
    assert len(calls) == 1


def test_missing_key_does_not_send_http(case, tmp_path, http_server, monkeypatch):
    url, calls, _ = http_server
    monkeypatch.delenv("PRISM_TEST_LLM_KEY", raising=False)
    config = _config(tmp_path, url, api_key_env="PRISM_TEST_LLM_KEY")
    result = plan_design(case, tmp_path / "plan", llm_config=config)
    assert result["last_error"]["error_code"] == "LLM_CREDENTIALS_MISSING"
    assert calls == []


@pytest.mark.parametrize(
    "url",
    [
        "http://example.com/v1",
        "https://key@example.com/v1",
        "https://example.com/v1?api_key=secret",
        "https://example.com/chat/completions",
    ],
)
def test_bad_endpoint_config(tmp_path, url):
    with pytest.raises(DesignError):
        load_llm_config(_config(tmp_path, url))


def test_cli_plan_and_replay(case, tmp_path, capsys):
    response = tmp_path / "response.json"
    write_json(response, proposal_for(case))
    assert main(["plan", "--request", str(case), "--response", str(response), "--output", str(tmp_path / "plan")]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "GRAPH_VALIDATED"
    assert main(["plan-check", "--plan", result["plan"]]) == 0


def test_existing_cli_validate_still_dispatches(case, capsys):
    spec = read_json(case)["base_design"]["path"]
    assert main(["validate", "--spec", spec]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "VALIDATED"


def test_top_level_cli_prepare_routes_to_planner(tmp_path, capsys):
    from prism.cli import main as prism_main

    spec = make_case(tmp_path / "source")
    goal = tmp_path / "goal.txt"
    goal.write_text("Test a polar terminal group", encoding="utf-8")
    assert (
        prism_main(
            [
                "design",
                "plan-prepare",
                "--spec",
                str(spec),
                "--goal-file",
                str(goal),
                "--output",
                str(tmp_path / "request"),
            ]
        )
        == 0
    )
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "PREPARED"


def test_existing_output_is_never_overwritten(case, tmp_path):
    response = tmp_path / "response.json"
    write_json(response, proposal_for(case))
    plan_design(case, tmp_path / "plan", response_path=response)
    with pytest.raises(DesignError) as exc:
        plan_design(case, tmp_path / "plan", response_path=response)
    assert exc.value.code == "OUTPUT_EXISTS"
