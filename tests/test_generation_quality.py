"""Quality control for generated ligands.

The central invariant under test: PRISM detects and proposes, but never applies
a speculative chemical repair. A molecule whose bond orders do not sanitize is
quarantined with its original bytes intact and a written proposal beside it; it
never reaches ``candidates.sdf``.
"""

import json
from pathlib import Path

import pytest

from prism.generation import quality
from prism.generation.postprocess import collect_candidates
from prism.generation.quality import QualityContext, assess


# A clean, MMFF-optimized N-methylbenzamide. Used as the calibration control:
# a check that flags this molecule is too strict to be useful.
CLEAN = """
     RDKit          3D

 10 10  0  0  0  0  0  0  0  0999 V2000
   -3.0559    0.2419   -0.0575 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4092    1.3905    0.3974 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0145    1.4465    0.4105 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2608    0.3485   -0.0228 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9152   -0.7974   -0.4942 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3113   -0.8489   -0.5068 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2249    0.4520    0.0091 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7692    1.5485   -0.0549 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9174   -0.7342    0.1527 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.3490   -0.7266    0.3198 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  4  7  1  0
  7  8  2  0
  7  9  1  0
  9 10  1  0
  6  1  1  0
M  END
$$$$
"""

# An ether oxygen carrying one double and one single bond: explicit valence 3.
# Exactly one minimal edit restores a valid model, mirroring pocketxmol_006.
UNIQUE_REPAIR = """valence-error
  PRISM    3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4000    0.0000    0.1000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1000    1.2000    0.2000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
M  END
$$$$
"""

# A nitrogen double-bonded to two different carbons. Reducing either bond is an
# equally minimal, equally neutral edit yielding a different molecule, so the
# intended structure is not recoverable from the 3D record. Mirrors the
# pocketxmol_010 triazole, which had three such variants.
AMBIGUOUS_REPAIR = """ambiguous
  PRISM    3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3000    0.0000    0.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7000    1.2000    0.2000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    2.6000    0.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  2  0
  3  4  1  0
M  END
$$$$
"""

# A carbon with two double bonds and two single bonds: valence 6. No single
# edit can fix it. MolCRAFT discards this class of failure before writing SDF.
UNREPAIRABLE = """pentavalent
  PRISM    3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7000    1.3000    0.2000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7000   -1.3000    0.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2000    0.1000    1.5000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
  1  5  2  0
M  END
$$$$
"""

rdkit_required = pytest.mark.skipif(
    not quality.rdkit_available(), reason="RDKit is an optional PRISM dependency"
)


def _unterminated_property(record: str) -> str:
    """Reproduce FLOWR's SD data fields, which omit the required blank line."""
    return record.rstrip().replace("$$$$", ">  <SCORE>\n0.42\n$$$$") + "\n"


def _write_run(tmp_path, *records):
    run_dir = tmp_path / "model-run"
    raw = run_dir / "raw"
    raw.mkdir(parents=True)
    for index, record in enumerate(records, start=1):
        (raw / f"result_{index:03d}.sdf").write_text(record, encoding="utf-8")
    return run_dir


# --------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------


@rdkit_required
def test_clean_molecule_passes_every_check():
    """Thresholds must not fire on a well-formed, relaxed structure."""
    report = assess(CLEAN, QualityContext(level="standard"))

    assert report.status == "PASS"
    assert report.flags == []
    assert report.sanitize_status == "passed"
    assert report.canonical_smiles == "CNC(=O)c1ccccc1"
    assert report.geometry["bad_bond_lengths"] == 0
    assert report.geometry["internal_clashes"] == 0


# --------------------------------------------------------------------------
# Repair proposals are detect-only
# --------------------------------------------------------------------------


@rdkit_required
def test_unique_minimal_repair_is_proposed_not_applied(tmp_path):
    run_dir = _write_run(tmp_path, UNIQUE_REPAIR, CLEAN)

    candidates, quarantined, failures, _ = collect_candidates(
        "pocketxmol", run_dir, ["raw/*.sdf"], 5
    )

    # The broken molecule never enters the candidate stream.
    assert [candidate.canonical_smiles for candidate in candidates] == ["CNC(=O)c1ccccc1"]
    assert "COC" not in (run_dir / "candidates.sdf").read_text(encoding="utf-8")
    assert failures

    assert len(quarantined) == 1
    record = quarantined[0]
    assert record.audit_status == "FAIL_CHEMISTRY"
    assert record.repair_confidence == "unique_minimal"
    assert record.repair_variant_count == 1

    # The original bytes are preserved as evidence, unmodified.
    assert (run_dir / "quarantine" / f"{record.candidate_id}.sdf").read_text(
        encoding="utf-8"
    ) == UNIQUE_REPAIR

    proposal = json.loads(open(record.repair_proposal_path, encoding="utf-8").read())
    assert proposal["applied"] is False
    assert proposal["confidence"] == "unique_minimal"
    assert proposal["variants"][0]["smiles"] == "COC"
    assert "bond 0-1: DOUBLE -> SINGLE" in proposal["variants"][0]["edits"]

    # The proposed structure exists only under quarantine/, for human review.
    variant = tmp_path / proposal["variant_paths"][0]
    assert variant.exists()
    assert "PRISM_REPAIR_APPLIED" in variant.read_text(encoding="utf-8")
    assert "false" in variant.read_text(encoding="utf-8")


@rdkit_required
def test_ambiguous_repair_lists_every_variant_and_picks_none(tmp_path):
    run_dir = _write_run(tmp_path, AMBIGUOUS_REPAIR, CLEAN)

    _, quarantined, _, _ = collect_candidates("pocketxmol", run_dir, ["raw/*.sdf"], 5)

    record = quarantined[0]
    assert record.repair_confidence == "ambiguous"
    assert record.repair_variant_count == 2

    proposal = json.loads(open(record.repair_proposal_path, encoding="utf-8").read())
    assert {variant["smiles"] for variant in proposal["variants"]} == {"CC=NC", "C=NCC"}
    assert proposal["applied"] is False
    assert "cannot choose" in proposal["reason"]
    # Every variant is written out; none is promoted over the others.
    assert len(proposal["variant_paths"]) == 2


@rdkit_required
def test_unrepairable_molecule_records_no_variants(tmp_path):
    run_dir = _write_run(tmp_path, UNREPAIRABLE, CLEAN)

    _, quarantined, _, _ = collect_candidates("molcraft", run_dir, ["raw/*.sdf"], 5)

    record = quarantined[0]
    assert record.repair_confidence == "none"
    assert record.repair_variant_count == 0
    proposal = json.loads(open(record.repair_proposal_path, encoding="utf-8").read())
    assert proposal["variant_paths"] == []
    assert "regenerate instead of repairing" in proposal["reason"]


# --------------------------------------------------------------------------
# Safe transformations
# --------------------------------------------------------------------------


@rdkit_required
def test_record_splitting_preserves_a_well_formed_trailing_data_field(tmp_path):
    """The splitter must not manufacture the defect the checker looks for.

    A trailing SD data field is terminated by a blank line. Stripping trailing
    newlines while splitting records deletes it, so correct model output would
    be reported as malformed and then "repaired" for no reason.
    """
    from prism.generation.postprocess import _sdf_records

    well_formed = CLEAN.rstrip().replace("$$$$", ">  <SCORE>\n0.42\n\n$$$$") + "\n"
    assert quality.sdf_data_fields_well_formed(well_formed)

    for count in (1, 3):
        records = list(_sdf_records(well_formed * count))
        assert len(records) == count
        for record in records:
            assert quality.sdf_data_fields_well_formed(record), record

    run_dir = _write_run(tmp_path, well_formed)
    candidates, quarantined, _, _ = collect_candidates("flowr", run_dir, ["raw/*.sdf"], 5)
    assert quarantined == []
    assert candidates[0].audit_status == "PASS"
    assert "sdf_data_fields_not_terminated" not in candidates[0].quality_flags


@rdkit_required
def test_sdf_format_repair_does_not_change_the_molecular_graph(tmp_path):
    """FLOWR's missing blank line is a layout defect, so normalizing is safe."""
    record = _unterminated_property(CLEAN)
    report = assess(record, QualityContext(level="standard"))
    assert report.status == "WARN_SDF_FORMAT"
    assert "sdf_data_fields_not_terminated" in report.flags

    run_dir = _write_run(tmp_path, record)
    candidates, quarantined, _, _ = collect_candidates("flowr", run_dir, ["raw/*.sdf"], 5)

    # A format warning is not a reason to withhold a chemically valid molecule.
    assert quarantined == []
    assert candidates[0].audit_status == "WARN_SDF_FORMAT"
    assert candidates[0].canonical_smiles == "CNC(=O)c1ccccc1"

    normalized = (run_dir / "candidates.sdf").read_text(encoding="utf-8")
    assert quality.sdf_data_fields_well_formed(normalized)
    assert assess(normalized, QualityContext(level="basic")).canonical_smiles == (
        "CNC(=O)c1ccccc1"
    )


@rdkit_required
def test_hydrogens_are_written_beside_candidates_without_moving_heavy_atoms(tmp_path):
    run_dir = _write_run(tmp_path, CLEAN)

    candidates, _, _, _ = collect_candidates("flowr", run_dir, ["raw/*.sdf"], 5)

    candidate = candidates[0]
    assert candidate.hydrogenated_path is not None
    # candidates.sdf keeps the model's own protonation state.
    assert (run_dir / "candidates.sdf").read_text(encoding="utf-8").count(" H ") == 0
    assert (run_dir / "candidates_h.sdf").is_file()

    hydrogen = candidate.quality["hydrogen"]
    assert hydrogen["hydrogens_added"] == 9
    assert hydrogen["heavy_atom_max_displacement_a"] == 0.0
    assert "pH-dependent protomers" in hydrogen["note"]


@rdkit_required
def test_hydrogen_output_can_be_disabled(tmp_path):
    run_dir = _write_run(tmp_path, CLEAN)

    candidates, _, _, _ = collect_candidates(
        "flowr", run_dir, ["raw/*.sdf"], 5, add_hydrogens=False
    )

    assert candidates[0].hydrogenated_path is None
    assert not (run_dir / "candidates_h.sdf").exists()


# --------------------------------------------------------------------------
# Pose plausibility
# --------------------------------------------------------------------------


@rdkit_required
def test_protein_clash_is_flagged_but_does_not_withhold_the_molecule(tmp_path):
    """The molecule is valid; only its placement is doubtful, so it stays."""
    import numpy as np

    clashing = np.array([[3.35, -0.73, 0.32]])  # on top of the terminal methyl
    context = QualityContext(
        level="full", protein_coordinates=clashing, pocket_center=(0.0, 0.0, 0.0)
    )

    report = assess(CLEAN, context)

    assert report.status == "WARN_POSE"
    assert report.usable
    assert report.pose["protein_hard_clashes"] >= 1
    assert report.pose["min_protein_distance_a"] < 2.0
    assert any(flag.startswith("protein_clash:") for flag in report.flags)


@rdkit_required
def test_pocket_center_offset_is_reported(tmp_path):
    context = QualityContext(level="full", pocket_center=(50.0, 0.0, 0.0))
    context.protein_coordinates = None

    report = assess(CLEAN, context)

    # Without a receptor the level degrades, so no pose section is produced.
    assert report.pose == {}

    import numpy as np

    context = QualityContext(
        level="full",
        protein_coordinates=np.array([[100.0, 100.0, 100.0]]),
        pocket_center=(50.0, 0.0, 0.0),
    )
    report = assess(CLEAN, context)
    assert report.pose["pocket_center_offset_a"] > 2.5
    assert any(flag.startswith("pocket_center_offset:") for flag in report.flags)


def test_protein_coordinates_skip_hydrogens_and_bad_records(tmp_path):
    pdb = tmp_path / "protein.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1      10.000  11.000  12.000  1.00 20.00           N\n"
        "ATOM      2  H   ALA A   1      10.500  11.500  12.500  1.00 20.00           H\n"
        "ATOM      3  CA  ALA A   1         bad     bad     bad  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    coordinates = quality.load_protein_coordinates(pdb)

    assert coordinates is not None
    assert len(coordinates) == 1


# --------------------------------------------------------------------------
# Sample accounting
# --------------------------------------------------------------------------


@rdkit_required
def test_quarantined_records_do_not_consume_the_sample_budget(tmp_path):
    run_dir = _write_run(tmp_path, UNIQUE_REPAIR, CLEAN, CLEAN)

    candidates, quarantined, _, attrition = collect_candidates(
        "pocketxmol", run_dir, ["raw/*.sdf"], 2
    )

    # Asking for two usable molecules yields two, not one plus a reject.
    assert len(candidates) == 2
    assert len(quarantined) == 1
    assert attrition["records_accepted"] == 2
    assert attrition["records_quarantined"] == 1
    assert attrition["records_written_by_model"] == 3


@rdkit_required
def test_attrition_reports_samples_dropped_before_any_sdf(tmp_path):
    """Geometry-perception models discard failures upstream, invisibly."""
    run_dir = _write_run(tmp_path, CLEAN, CLEAN)

    _, _, _, attrition = collect_candidates("targetdiff", run_dir, ["raw/*.sdf"], 5)

    assert attrition["requested"] == 5
    assert attrition["records_written_by_model"] == 2
    assert attrition["dropped_before_sdf"] == 3


@rdkit_required
def test_oversampling_is_counted_as_surplus(tmp_path):
    run_dir = _write_run(tmp_path, CLEAN, CLEAN, CLEAN)

    candidates, _, _, attrition = collect_candidates("pocket2mol", run_dir, ["raw/*.sdf"], 1)

    assert len(candidates) == 1
    assert attrition["records_surplus"] == 2
    assert attrition["dropped_before_sdf"] == 0


@rdkit_required
def test_candidates_follow_model_numbering_not_lexicographic_order(tmp_path):
    """Truncating a lexicographic order at --num-samples delivers a subset.

    Pocket2Mol emits molecules in completion order, and completion order is
    correlated with size, so the wrong ordering silently changes which
    molecules the user receives.
    """
    run_dir = tmp_path / "model-run"
    raw = run_dir / "raw"
    raw.mkdir(parents=True)
    for index in range(12):
        (raw / f"{index}.sdf").write_text(
            CLEAN.replace("\n     RDKit", f"mol{index}\n     RDKit", 1), encoding="utf-8"
        )

    candidates, _, _, attrition = collect_candidates(
        "pocket2mol", run_dir, ["raw/*.sdf"], 5
    )

    assert [Path(c.source_path).name for c in candidates] == [f"{i}.sdf" for i in range(5)]
    assert attrition["records_surplus"] == 7


@rdkit_required
def test_qc_report_covers_every_record(tmp_path):
    import csv

    run_dir = _write_run(tmp_path, UNIQUE_REPAIR, CLEAN)

    collect_candidates("pocketxmol", run_dir, ["raw/*.sdf"], 5)

    with (run_dir / "qc_report.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["disposition"] for row in rows} == {"candidate", "quarantined"}
    assert {row["audit_status"] for row in rows} == {"PASS", "FAIL_CHEMISTRY"}

    issues = (run_dir / "qc_issues.jsonl").read_text(encoding="utf-8").strip().splitlines()
    assert len(issues) == 1
    assert json.loads(issues[0])["audit_status"] == "FAIL_CHEMISTRY"


# --------------------------------------------------------------------------
# Optional dependency
# --------------------------------------------------------------------------


def test_quality_degrades_when_rdkit_is_unavailable(monkeypatch, tmp_path):
    """RDKit is optional; without it PRISM keeps its pre-QC behavior."""
    monkeypatch.setattr(quality, "rdkit_available", lambda: False)

    report = assess(UNIQUE_REPAIR, QualityContext(level="full"))
    assert report.status == "NOT_CHECKED"

    run_dir = _write_run(tmp_path, UNIQUE_REPAIR)
    candidates, quarantined, _, _ = collect_candidates("pocketxmol", run_dir, ["raw/*.sdf"], 5)

    # No chemistry judgment is made, so nothing is withheld on chemical grounds.
    assert len(candidates) == 1
    assert quarantined == []
    assert candidates[0].audit_status == "NOT_CHECKED"


def test_qc_can_be_switched_off(tmp_path):
    report = assess(UNIQUE_REPAIR, QualityContext(level="off"))

    assert report.status == "NOT_CHECKED"
    assert report.repair is None
