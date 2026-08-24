"""Collect model-native SDF files into PRISM's normalized output layout.

Collection is also where quality control happens.  Records that pass, or that
only raise warnings, become candidates.  Records that fail are *quarantined*
rather than deleted: the original bytes are preserved and, where a minimal
valence correction exists, a repair proposal is written next to them.  PRISM
never applies such a proposal itself -- see ``prism.generation.quality``.
"""

import csv
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from . import quality
from .errors import ExecutionFailure
from .quality import QualityContext, QualityReport
from .types import CandidateRecord, QuarantineRecord


def _natural_key(path: Path) -> Tuple[Any, ...]:
    """Order ``2.sdf`` before ``10.sdf``.

    Several models number their output without zero padding. Lexicographic
    order would then interleave the indices, and truncating that order at
    ``--num-samples`` would deliver an arbitrary subset rather than the first
    molecules the model produced.
    """
    return tuple(
        (1, int(part), "") if part.isdigit() else (0, 0, part)
        for part in re.split(r"(\d+)", str(path))
    )


def _sdf_records(text: str) -> Iterable[str]:
    for index, fragment in enumerate(text.split("$$$$")):
        if fragment.strip():
            # Preserve the three-line molfile header: an empty molecule title is
            # valid and is emitted by TargetDiff/RDKit. Only discard the newline
            # belonging to a preceding record delimiter.
            if index and fragment.startswith("\r\n"):
                fragment = fragment[2:]
            elif index and fragment.startswith("\n"):
                fragment = fragment[1:]
            # Keep the body verbatim apart from guaranteeing one final newline.
            # Stripping trailing newlines would delete the blank line that
            # terminates a trailing SD data field -- which is part of the
            # format, not padding -- and would make well-formed model output
            # look defective to the quality checks downstream.
            if not fragment.endswith("\n"):
                fragment += "\n"
            yield fragment + "$$$$\n"


def _coordinates(record: str) -> Tuple[List[Tuple[float, float, float]], bool]:
    lines = record.splitlines()
    if len(lines) < 5:
        return [], False
    try:
        atom_count = int(lines[3][:3])
    except ValueError:
        return [], False
    coordinates = []
    for line in lines[4 : 4 + atom_count]:
        try:
            coordinate = (float(line[0:10]), float(line[10:20]), float(line[20:30]))
        except ValueError:
            return [], False
        if not all(math.isfinite(value) for value in coordinate):
            return [], False
        coordinates.append(coordinate)
    dimension_marked = len(lines) > 1 and "3D" in lines[1].upper()
    non_planar = any(abs(coordinate[2]) > 1e-6 for coordinate in coordinates)
    return coordinates, bool(coordinates) and (dimension_marked or non_planar)


def _with_properties(record: str, candidate_id: str, model: str) -> str:
    body = record.rsplit("$$$$", 1)[0].rstrip()
    # ``rstrip`` removes the blank line that terminates a trailing SD data
    # field, so it has to be restored before appending more fields -- otherwise
    # PRISM reintroduces exactly the layout defect it normalizes FLOWR to fix.
    separator = "\n\n" if any(line.startswith(">") for line in body.splitlines()) else "\n"
    return (
        f"{body}{separator}"
        f">  <PRISM_CANDIDATE_ID>\n{candidate_id}\n\n"
        f">  <PRISM_MODEL>\n{model}\n\n$$$$\n"
    )


def _sdf_property(record: str, name: str) -> Optional[str]:
    match = re.search(
        rf"^>\s*<\s*{re.escape(name)}\s*>\s*\r?\n([^\r\n]*)",
        record,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    return match.group(1).strip() if match else None


def _native_score(record: str, model: str) -> Tuple[Optional[str], Optional[float]]:
    properties = {"pocketxmol": ("cfd_traj", "POCKETXMOL_CFD_TRAJ")}
    specification = properties.get(model)
    if specification is None:
        return None, None
    score_name, property_name = specification
    value = _sdf_property(record, property_name)
    if value is None:
        return score_name, None
    try:
        score = float(value)
    except ValueError:
        return score_name, None
    return (score_name, score) if math.isfinite(score) else (score_name, None)


# --------------------------------------------------------------------------
# Quality-control reporting
# --------------------------------------------------------------------------

QC_REPORT_FIELDS = [
    "candidate_id",
    "model",
    "disposition",
    "audit_status",
    "flags",
    "sanitize_status",
    "parse_mode",
    "canonical_smiles",
    "source_file",
    "source_record",
    "path",
    "atom_count",
    "bond_count",
    "fragments",
    "explicit_h",
    "hydrogens_added",
    "heavy_atom_max_displacement_a",
    "formal_charge",
    "molecular_weight",
    "clogp",
    "tpsa",
    "rotatable_bonds",
    "qed",
    "synthetic_accessibility",
    "lipinski_violations",
    "stereocenters",
    "unspecified_stereocenters",
    "structural_alerts",
    "reactive_groups",
    "bad_bond_lengths",
    "internal_clashes",
    "angle_outliers",
    "worst_angle_deviation_deg",
    "aromatic_planarity_rms_a",
    "min_protein_distance_a",
    "protein_hard_clashes",
    "protein_close_contacts",
    "pocket_center_offset_a",
    "repair_confidence",
    "repair_variant_count",
]


def _report_row(
    candidate_id: str,
    model: str,
    disposition: str,
    report: QualityReport,
    source_file: Path,
    source_record: int,
    path: Optional[Path],
) -> Dict[str, Any]:
    """Flatten a QualityReport into one CSV row."""

    def get(section: Dict[str, Any], key: str) -> Any:
        value = section.get(key)
        if isinstance(value, list):
            return ",".join(str(item) for item in value)
        return "" if value is None else value

    row = {name: "" for name in QC_REPORT_FIELDS}
    row.update(
        {
            "candidate_id": candidate_id,
            "model": model,
            "disposition": disposition,
            "audit_status": report.status,
            "flags": ";".join(report.flags),
            "sanitize_status": report.sanitize_status,
            "parse_mode": report.parse_mode,
            "canonical_smiles": report.canonical_smiles,
            "source_file": str(source_file),
            "source_record": source_record,
            "path": str(path) if path is not None else "",
        }
    )
    for section, keys in (
        (report.integrity, ("atom_count", "bond_count", "fragments", "explicit_h")),
        (report.hydrogen, ("hydrogens_added", "heavy_atom_max_displacement_a")),
        (
            report.plausibility,
            (
                "formal_charge",
                "molecular_weight",
                "clogp",
                "tpsa",
                "rotatable_bonds",
                "qed",
                "synthetic_accessibility",
                "lipinski_violations",
                "stereocenters",
                "unspecified_stereocenters",
                "structural_alerts",
                "reactive_groups",
            ),
        ),
        (
            report.geometry,
            (
                "bad_bond_lengths",
                "internal_clashes",
                "angle_outliers",
                "worst_angle_deviation_deg",
                "aromatic_planarity_rms_a",
            ),
        ),
        (
            report.pose,
            (
                "min_protein_distance_a",
                "protein_hard_clashes",
                "protein_close_contacts",
                "pocket_center_offset_a",
            ),
        ),
    ):
        for key in keys:
            row[key] = get(section, key)
    if report.repair is not None:
        row["repair_confidence"] = report.repair.confidence
        row["repair_variant_count"] = len(report.repair.variants)
    return row


def _write_repair_proposal(
    quarantine_dir: Path, candidate_id: str, report: QualityReport
) -> Optional[Path]:
    """Write the proposal and every variant, clearly separated from candidates.

    Variants live under ``quarantine/`` and are never concatenated into
    ``candidates.sdf``.  Applying one is a decision for the user.
    """
    if report.repair is None:
        return None
    proposal_path = quarantine_dir / f"{candidate_id}.repair_proposal.json"
    payload = dict(report.repair.to_dict())
    payload["candidate_id"] = candidate_id
    payload["sanitize_status"] = report.sanitize_status
    payload["policy"] = (
        "PRISM does not apply speculative repairs. Review the variants, then "
        "copy the chosen one into your working set yourself."
    )
    variant_paths = []
    for index, variant in enumerate(report.repair.variants, start=1):
        variant_path = quarantine_dir / f"{candidate_id}.variant_{index:02d}.sdf"
        variant_path.write_text(
            variant["molblock"].rstrip("\n")
            + f"\n>  <PRISM_REPAIR_VARIANT>\n{index}\n\n"
            + f">  <PRISM_REPAIR_EDITS>\n{variant['edits']}\n\n"
            + f">  <PRISM_REPAIR_CONFIDENCE>\n{report.repair.confidence}\n\n"
            + ">  <PRISM_REPAIR_APPLIED>\nfalse\n\n$$$$\n",
            encoding="utf-8",
        )
        variant_paths.append(str(variant_path))
    payload["variant_paths"] = variant_paths
    proposal_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return proposal_path


def collect_candidates(
    model: str,
    run_dir: Path,
    patterns: Sequence[str],
    limit: int,
    context: Optional[QualityContext] = None,
    add_hydrogens: bool = True,
) -> Tuple[List[CandidateRecord], List[QuarantineRecord], List[str], Dict[str, Any]]:
    """Validate, quality-check, and split every discovered SDF record.

    Returns the accepted candidates, the quarantined records, human-readable
    failure messages, and an attrition summary describing where samples were
    lost between the request and the usable output.
    """
    context = context or QualityContext()
    source_files = []
    for pattern in patterns:
        source_files.extend(path for path in run_dir.glob(pattern) if path.is_file())
    source_files = sorted(set(source_files), key=_natural_key)
    if not source_files:
        raise ExecutionFailure(
            f"No SDF output matched {list(patterns)} in {run_dir}",
            code="NO_MODEL_OUTPUT",
        )

    molecules_dir = run_dir / "molecules"
    molecules_dir.mkdir(parents=True, exist_ok=True)
    hydrogen_dir = run_dir / "molecules_h"
    quarantine_dir = run_dir / "quarantine"

    candidates: List[CandidateRecord] = []
    quarantined: List[QuarantineRecord] = []
    failures: List[str] = []
    report_rows: List[Dict[str, Any]] = []
    hydrogen_blocks: List[str] = []
    total_records = 0
    surplus = 0

    for source_file in source_files:
        try:
            records = list(_sdf_records(source_file.read_text(encoding="utf-8", errors="replace")))
        except OSError as exc:
            failures.append(f"{source_file}: unable to read: {exc}")
            continue
        for record_index, record in enumerate(records, start=1):
            total_records += 1
            coordinates, valid_3d = _coordinates(record)
            report = quality.assess(record, context)

            if report.status == "NOT_CHECKED":
                # RDKit unavailable or QC disabled: fall back to the header and
                # coordinate parsing that PRISM has always performed.
                usable = bool(coordinates) and valid_3d
                if not usable:
                    report.flags = ["invalid_3d_sdf_record"]
            else:
                usable = report.usable

            if len(candidates) >= limit and usable:
                # The model oversampled. Keep counting so the attrition summary
                # stays honest, but do not emit more candidates than requested.
                surplus += 1
                report_rows.append(
                    _report_row(
                        f"{model}_surplus_{surplus:06d}",
                        model,
                        "surplus",
                        report,
                        source_file,
                        record_index,
                        None,
                    )
                )
                continue

            if not usable:
                quarantine_dir.mkdir(parents=True, exist_ok=True)
                quarantine_id = f"{model}_q{len(quarantined) + 1:06d}"
                quarantine_path = quarantine_dir / f"{quarantine_id}.sdf"
                # Preserve the original bytes; a quarantined record is evidence.
                quarantine_path.write_text(record, encoding="utf-8")
                proposal_path = _write_repair_proposal(quarantine_dir, quarantine_id, report)
                reason = "; ".join(report.flags) or "invalid 3D SDF record"
                quarantined.append(
                    QuarantineRecord(
                        candidate_id=quarantine_id,
                        model=model,
                        path=str(quarantine_path),
                        source_path=str(source_file),
                        audit_status=report.status,
                        reason=reason,
                        quality_flags=list(report.flags),
                        repair_confidence=(
                            report.repair.confidence if report.repair is not None else "none"
                        ),
                        repair_variant_count=(
                            len(report.repair.variants) if report.repair is not None else 0
                        ),
                        repair_proposal_path=str(proposal_path) if proposal_path else None,
                        quality=report.to_dict(),
                    )
                )
                failures.append(
                    f"{source_file} record {record_index}: {report.status}: {reason}"
                )
                report_rows.append(
                    _report_row(
                        quarantine_id,
                        model,
                        "quarantined",
                        report,
                        source_file,
                        record_index,
                        quarantine_path,
                    )
                )
                continue

            candidate_id = f"{model}_{len(candidates) + 1:06d}"
            candidate_path = molecules_dir / f"candidate_{len(candidates) + 1:06d}.sdf"

            emitted = record
            if report.status == "WARN_SDF_FORMAT":
                # Layout-only repair: the molecular graph is untouched.
                normalized = quality.normalize_record(record)
                if normalized is not None:
                    emitted = normalized
            candidate_path.write_text(
                _with_properties(emitted, candidate_id, model), encoding="utf-8"
            )

            hydrogenated_path = None
            if add_hydrogens and report.status != "NOT_CHECKED":
                block, hydrogen_info = quality.hydrogenate_record(emitted)
                report.hydrogen = hydrogen_info
                if block is not None:
                    hydrogen_dir.mkdir(parents=True, exist_ok=True)
                    target = hydrogen_dir / f"candidate_{len(candidates) + 1:06d}_h.sdf"
                    target.write_text(
                        _with_properties(block, candidate_id, model), encoding="utf-8"
                    )
                    hydrogen_blocks.append(target.read_text(encoding="utf-8"))
                    hydrogenated_path = str(target)

            native_score_name, native_score = _native_score(record, model)
            candidates.append(
                CandidateRecord(
                    candidate_id=candidate_id,
                    model=model,
                    path=str(candidate_path),
                    source_path=str(source_file),
                    valid_3d=valid_3d or report.integrity.get("is_3d", False),
                    sanitize_status=(
                        "passed" if report.sanitize_status == "passed" else report.sanitize_status
                    ),
                    native_score_name=native_score_name,
                    native_score=native_score,
                    audit_status=report.status,
                    quality_flags=list(report.flags),
                    canonical_smiles=report.canonical_smiles,
                    hydrogenated_path=hydrogenated_path,
                    quality=report.to_dict(),
                )
            )
            report_rows.append(
                _report_row(
                    candidate_id,
                    model,
                    "candidate",
                    report,
                    source_file,
                    record_index,
                    candidate_path,
                )
            )

    if not candidates:
        raise ExecutionFailure(
            "The model produced no valid 3D SDF molecules"
            + (
                f" ({len(quarantined)} quarantined for review in {quarantine_dir})"
                if quarantined
                else ""
            ),
            code="NO_VALID_MOLECULES",
        )

    combined_path = run_dir / "candidates.sdf"
    combined_path.write_text(
        "".join(Path(candidate.path).read_text(encoding="utf-8") for candidate in candidates),
        encoding="utf-8",
    )
    if hydrogen_blocks:
        (run_dir / "candidates_h.sdf").write_text("".join(hydrogen_blocks), encoding="utf-8")

    write_qc_report(run_dir, report_rows)

    attrition = {
        "requested": limit,
        "records_written_by_model": total_records,
        "records_accepted": len(candidates),
        "records_quarantined": len(quarantined),
        "records_surplus": surplus,
        # Models that perceive bond orders from coordinates drop their own
        # reconstruction failures before writing any SDF. PRISM cannot inspect
        # or repair those; it can only make the loss visible.
        "dropped_before_sdf": max(0, limit - total_records),
    }
    return candidates, quarantined, failures, attrition


def write_qc_report(directory: Path, rows: Sequence[Dict[str, Any]]) -> None:
    """Write ``qc_report.csv`` and the non-passing subset as ``qc_issues.jsonl``."""
    if not rows:
        return
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "qc_report.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=QC_REPORT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    issues = [row for row in rows if row.get("audit_status") not in {"PASS", "NOT_CHECKED"}]
    (directory / "qc_issues.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True, default=str) + "\n" for row in issues),
        encoding="utf-8",
    )
