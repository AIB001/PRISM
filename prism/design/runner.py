"""Auditable edit runs with a denominator of raw slots, including failures."""

from collections import Counter
import csv
import json
from pathlib import Path
import shutil

import numpy as np
from rdkit import Chem

from .backends import execute_backend, load_backend
from .chem import heavy_copy
from .compile import compile_design
from .errors import DesignError, require
from .io import file_ref, read_json, verify_ref, write_json
from .schema import load_design
from .validate import ground_design
from .verify import verify_candidate


def _csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _attempts(directory, samples):
    """Missing attempts are NOT_ATTEMPTED, never invented decode failures."""
    path = directory / "attempts.jsonl"
    events = []
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            try:
                event = json.loads(line)
            except json.JSONDecodeError as exc:
                raise DesignError("Invalid raw slot ledger", "LEDGER_INVALID") from exc
            require(
                isinstance(event, dict)
                and event.get("slot") == len(events) + 1
                and event.get("status") in {"DECODED", "DECODE_FAILED", "SERIALIZE_FAILED"},
                "Nonsequential/invalid raw slot ledger",
                "LEDGER_INVALID",
            )
            events.append(event)
    require(len(events) <= samples, "Backend exceeded the raw sampling budget", "LEDGER_INVALID")
    events.extend(
        {"slot": slot, "status": "NOT_ATTEMPTED", "sdf": None, "error": "Backend did not report this slot"}
        for slot in range(len(events) + 1, samples + 1)
    )
    return events


def evaluate_slots(grounded, backend_dir, output, samples):
    """Evaluate all slots and emit only unique, passing candidates for MD handoff."""
    from prism.generation.quality import QualityContext, assess, hydrogenate_record
    from prism.generation.types import CandidateRecord

    backend_dir, output = Path(backend_dir), Path(output)
    events = _attempts(backend_dir, samples)
    context = QualityContext(
        level="full",
        protein_coordinates=grounded.protein_xyz,
        pocket_center=tuple(np.mean(heavy_copy(grounded.molecule).GetConformer().GetPositions(), axis=0)),
    )
    records, selected, seen = [], [], {}
    (output / "validation").mkdir(exist_ok=True)
    (output / "candidates").mkdir(exist_ok=True)
    for event in events:
        slot = event["slot"]
        record = dict(
            event,
            candidate_id=f"flowr_root_{slot:04d}",
            validation=None,
            qc=None,
            selected=False,
            duplicate_of=None,
            hydrogenation=None,
        )
        if event["status"] == "DECODED":
            source = Path(event["sdf"]).resolve()
            require(
                backend_dir.resolve() in source.parents and source.is_file(),
                "Decoded slot points outside backend outputs or is missing",
                "LEDGER_INVALID",
            )
            record["raw_file"] = file_ref(source)
            text = source.read_text(encoding="utf-8")
            supplier = Chem.SDMolSupplier(str(source), sanitize=False, removeHs=False)
            molecule = supplier[0] if len(supplier) == 1 else None
            validation = verify_candidate(grounded, molecule)
            qc = assess(text, context)
            record.update(validation=validation, qc=qc.to_dict())
            accepted = validation["selection_pass"] and qc.status == "PASS"
            smiles = validation.get("canonical_smiles")
            if accepted and smiles in seen:
                record["duplicate_of"] = seen[smiles]
            elif accepted:
                seen[smiles] = record["candidate_id"]
                record["selected"] = True
                target = output / "candidates" / f"{record['candidate_id']}.sdf"
                shutil.copyfile(source, target)
                hydrogens, info = hydrogenate_record(text)
                record["hydrogenation"] = info
                hydrogen_path = None
                if hydrogens is not None:
                    hpath = target.with_name(target.stem + "_with_H.sdf")
                    hpath.write_text(hydrogens, encoding="utf-8")
                    hydrogen_path = str(hpath)
                candidate = CandidateRecord(
                    candidate_id=record["candidate_id"],
                    model="flowr_root",
                    path=str(target),
                    source_path=str(source),
                    valid_3d=True,
                    sanitize_status=qc.sanitize_status,
                    audit_status=qc.status,
                    quality_flags=qc.flags,
                    canonical_smiles=smiles,
                    hydrogenated_path=hydrogen_path,
                    quality=qc.to_dict(),
                ).to_dict()
                candidate["design_validation"] = {"path": str(output / "validation" / f"slot_{slot:04d}.json")}
                selected.append(candidate)
        validation_path = output / "validation" / f"slot_{slot:04d}.json"
        write_json(validation_path, record)
        records.append(record)
    # Bind completed validation records, including MD hydrogenation provenance.
    for candidate in selected:
        candidate["design_validation"] = file_ref(candidate["design_validation"]["path"])
    return records, selected


def summarize_slots(records):
    validations = [r["validation"] for r in records if r["validation"] is not None]
    reasons = Counter(c["id"] + ":" + c["status"] for v in validations for c in v["checks"] if c["status"] != "pass")
    reasons.update(
        "unchanged_reference"
        for v in validations
        if v.get("selection_quality", {}).get("changed_from_reference") is False
    )
    reasons.update(
        "pose_screen_failed" for v in validations if v.get("selection_quality", {}).get("pose_screen_pass") is False
    )
    reasons.update("prism_qc:" + r["qc"]["status"] for r in records if r["qc"] and r["qc"]["status"] != "PASS")
    return {
        "requested": len(records),
        "raw_slots": sum(r["status"] != "NOT_ATTEMPTED" for r in records),
        "decoded": sum(r["status"] == "DECODED" for r in records),
        "decode_failed": sum(r["status"] == "DECODE_FAILED" for r in records),
        "serialize_failed": sum(r["status"] == "SERIALIZE_FAILED" for r in records),
        "not_attempted": sum(r["status"] == "NOT_ATTEMPTED" for r in records),
        "hard_constraints_passed": sum(v["design_hard_pass"] for v in validations),
        "design_selection_passed": sum(v["selection_pass"] for v in validations),
        "design_and_qc_passed": sum(
            bool(r["validation"] and r["validation"]["selection_pass"] and r["qc"]["status"] == "PASS") for r in records
        ),
        "unique_selected": sum(r["selected"] for r in records),
        "duplicate_selected_chemistry": sum(r["duplicate_of"] is not None for r in records),
        "md_hydrogen_files": sum(
            bool(
                r["selected"] and r["hydrogenation"] and r["hydrogenation"]["status"] in {"added", "explicit_complete"}
            )
            for r in records
        ),
        "failure_counts_nonexclusive": dict(reasons),
    }


def run_design(spec, backend_config, output, samples=16, seed=20261013, input_root=None):
    require(type(samples) is int and samples > 0, "Sample count must be positive", "INVALID_BUDGET")
    require(type(seed) is int and 0 <= seed < 2**32, "Seed must be a uint32", "INVALID_SEED")
    grounded = ground_design(load_design(spec, input_root), for_execution=True)
    config = load_backend(backend_config)
    if grounded.spec["task"] == "fragment_growing":
        lo, hi = grounded.spec["editing"]["new_fragment_heavy_atoms_range"]
        require(lo <= config["grow_size"] <= hi, "Configured growth size is outside the DesignSpec", "GROWTH_SIZE")
    output = Path(output).expanduser().resolve()
    require(not output.exists(), f"Output already exists; choose a new run directory: {output}", "OUTPUT_EXISTS")
    output.mkdir(parents=True)
    write_json(output / "status.json", {"status": "RUNNING", "requested": samples, "seed": seed})
    try:
        compiled = compile_design(grounded, output)
        try:
            process = execute_backend(config, output / "compiled_design.json", output / "backend", samples, seed)
        except (DesignError, OSError) as exc:
            process = {
                "status": "BACKEND_PREFLIGHT_FAILED",
                "error_code": getattr(exc, "code", "IO_ERROR"),
                "message": str(exc),
            }
        records, candidates = evaluate_slots(grounded, output / "backend", output, samples)
        counts = summarize_slots(records)
        complete = process["status"] == "COMPLETED" and counts["raw_slots"] == samples
        status = "COMPLETED" if complete else "FAILED"
        summary = {
            "schema": "prism-design-run/1",
            "status": status,
            "complex_id": grounded.spec["complex_id"],
            "task": grounded.spec["task"],
            "seed": seed,
            "compiled_design": file_ref(output / "compiled_design.json"),
            "counts": counts,
            "backend_process": process,
            "input_review": grounded.report(),
            "selection_policy": "all hard constraints + all measurable soft objectives + changed molecule + pose screen + PRISM PASS QC; deduplicate accepted chemical identities",
            "deduplication_scope": "within_case_canonical_isomeric_smiles; all_raw_poses_retained",
            "experimental_activity_improvement": None,
            "independent_generalization_benchmark": False,
            "training_performed": False,
            "raw_outputs_repaired": False,
        }
        write_json(output / "results.json", {"summary": summary, "slots": records})
        _csv(
            output / "slots.csv",
            [
                {
                    "slot": r["slot"],
                    "raw_status": r["status"],
                    "hard_pass": (r["validation"] or {}).get("design_hard_pass", False),
                    "design_selection_pass": (r["validation"] or {}).get("selection_pass", False),
                    "qc_status": (r["qc"] or {}).get("status", "NOT_CHECKED"),
                    "selected": r["selected"],
                    "duplicate_of": r["duplicate_of"],
                    "smiles": (r["validation"] or {}).get("canonical_smiles", ""),
                    "failed_checks": ";".join(
                        c["id"] + ":" + c["status"]
                        for c in (r["validation"] or {}).get("checks", [])
                        if c["status"] != "pass"
                    ),
                }
                for r in records
            ],
        )
        write_json(
            output / "manifest.json",
            {
                "schema": "prism-design-md-handoff/1",
                "status": status,
                "request": {
                    "protein": compiled["inputs"]["receptor"]["path"],
                    "pocket": compiled["inputs"]["pocket"]["path"],
                    "output": str(output),
                    "num_samples": samples,
                    "seed": seed,
                    "models": ["flowr_root"],
                },
                "design_results": file_ref(output / "results.json"),
                "models": {
                    "flowr_root": {
                        "model": "flowr_root",
                        "status": status,
                        "requested": samples,
                        "generated": len(candidates),
                        "candidates": candidates,
                        "attrition": counts,
                        "provenance": process.get("provenance", {}),
                    }
                },
            },
        )
        write_json(output / "status.json", summary)
        return summary
    except Exception as exc:
        write_json(
            output / "status.json",
            {
                "status": "FAILED",
                "requested": samples,
                "error_code": getattr(exc, "code", "RUN_ERROR"),
                "message": str(exc),
                "note": "Raw backend ledger and files, if present, remain available for audit",
            },
        )
        raise


def run_pilot(queue, backend_config, output, seed=20261013):
    queue = Path(queue).expanduser().resolve()
    entries = read_json(queue)["first_smoke_cases"]
    require(
        entries and len({e["pdb_id"] for e in entries}) == len(entries), "Empty/duplicate pilot cases", "QUEUE_INVALID"
    )
    jobs = []
    for entry in entries:
        spec = verify_ref(queue.parent, entry["designspec"])
        document = load_design(spec)
        require(
            document.data["complex_id"] == entry["pdb_id"] and document.data["task"] == entry["task"],
            "Queue identity/task mismatch",
            "QUEUE_INVALID",
        )
        require(
            type(entry["proposed_candidate_slots"]) is int and entry["proposed_candidate_slots"] > 0,
            "Invalid queue budget",
            "QUEUE_INVALID",
        )
        jobs.append((entry, spec))
    output = Path(output).expanduser().resolve()
    require(not output.exists(), "Pilot output already exists; choose a new directory", "OUTPUT_EXISTS")
    output.mkdir(parents=True)
    write_json(
        output / "request.json", {"queue": file_ref(queue), "backend_config": file_ref(backend_config), "seed": seed}
    )
    summaries = []
    for index, (entry, spec) in enumerate(jobs):
        case = entry["pdb_id"]
        print(
            f"PRISM-Design {index + 1}/{len(jobs)}: {case}, {entry['proposed_candidate_slots']} raw slots", flush=True
        )
        try:
            result = run_design(spec, backend_config, output / case, entry["proposed_candidate_slots"], seed + index)
        except Exception as exc:
            result = {
                "complex_id": case,
                "status": "FAILED",
                "task": entry["task"],
                "error_code": getattr(exc, "code", "RUN_ERROR"),
                "message": str(exc),
            }
        summaries.append(result)
        write_json(output / "pilot_summary.json", {"status": "RUNNING", "cases": summaries})
        print(f"  {result['status']}: {result.get('counts', result.get('message'))}", flush=True)
    totals = Counter()
    for result in summaries:
        totals.update({key: value for key, value in result.get("counts", {}).items() if isinstance(value, int)})
    report = {
        "schema": "prism-design-pilot/1",
        "status": "COMPLETED" if all(r["status"] == "COMPLETED" for r in summaries) else "FAILED",
        "queue": file_ref(queue),
        "planned_slots": sum(e["proposed_candidate_slots"] for e in entries),
        "totals": dict(totals),
        "cases": summaries,
        "training_performed": False,
        "independent_generalization_benchmark": False,
        "experimental_activity_improvement": None,
    }
    write_json(output / "pilot_summary.json", report)
    _csv(
        output / "pilot_summary.csv",
        [
            {
                "complex_id": r["complex_id"],
                "task": r["task"],
                "status": r["status"],
                **{key: r.get("counts", {}).get(key) for key in totals},
            }
            for r in summaries
        ],
    )
    return report
