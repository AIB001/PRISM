"""Hand-off from a finished generation run to MD system building.

The invariant under test: a ligand only leaves this layer if the file about to
be handed over carries explicit hydrogens. Generative models emit heavy atoms
only, and parameterization does not reject that input -- it writes a topology
with no hydrogens at all, correct in layout and wrong in chemistry. Measured
over 82 gaff2 builds: 0 of 41 raw inputs produced a correct topology, against
35 of 41 hydrogenated ones.

The second invariant: nothing is dropped silently. Every candidate the hand-off
withholds is reported with a reason.
"""

import json
from pathlib import Path

import pytest

from prism.generation.errors import GenerationError
from prism.generation.handoff import (
    count_explicit_hydrogens,
    export_md_inputs,
    select_for_md,
)


# Methanol with explicit hydrogens. Column positions matter: the no-RDKit
# fallback reads the element symbol out of the fixed-width V2000 atom block.
WITH_HYDROGENS = """methanol-h
  PRISM    3D

  6  5  0  0  0  0  0  0  0  0999 V2000
   -0.3480    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0150    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7180    0.9420   -0.4190 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7180   -0.0640    1.0280 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7180   -0.8780   -0.6090 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.3250    0.0000   -0.9060 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
M  END
$$$$
"""

# The same molecule as every model actually emits it: heavy atoms only.
HEAVY_ATOMS_ONLY = """methanol
  PRISM    3D

  2  1  0  0  0  0  0  0  0  0999 V2000
   -0.3480    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0150    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
"""


def _write_run(tmp_path, candidates, *, protein="/data/receptor.pdb"):
    """Build a generation output directory the way the runner leaves one.

    ``candidates`` is a list of dicts with keys: model, status, flags (opt),
    hydrogens (bool or "missing_file" or "absent_record").
    """
    run_dir = tmp_path / "generated_ligand"
    models = {}
    counters = {}

    for spec in candidates:
        model = spec["model"]
        counters[model] = counters.get(model, 0) + 1
        index = counters[model]
        candidate_id = f"{model}_{index:06d}"

        molecules = run_dir / "models" / model / "molecules"
        molecules.mkdir(parents=True, exist_ok=True)
        heavy = molecules / f"candidate_{index:06d}.sdf"
        heavy.write_text(HEAVY_ATOMS_ONLY, encoding="utf-8")

        hydrogens = spec.get("hydrogens", True)
        hydrogenated_path = None
        if hydrogens != "absent_record":
            hydrogen_dir = run_dir / "models" / model / "hydrogens"
            hydrogen_dir.mkdir(parents=True, exist_ok=True)
            target = hydrogen_dir / f"candidate_{index:06d}_h.sdf"
            hydrogenated_path = str(target)
            if hydrogens is True:
                target.write_text(WITH_HYDROGENS, encoding="utf-8")
            elif hydrogens == "no_hydrogens":
                # The manifest claims a hydrogenated copy; the file has none.
                target.write_text(HEAVY_ATOMS_ONLY, encoding="utf-8")
            # "missing_file": recorded in the manifest, never written

        models.setdefault(model, {"status": "SUCCEEDED", "candidates": []})
        models[model]["candidates"].append(
            {
                "candidate_id": candidate_id,
                "model": model,
                "path": str(heavy),
                "source_path": str(heavy),
                "valid_3d": True,
                "audit_status": spec["status"],
                "quality_flags": spec.get("flags", []),
                "canonical_smiles": "CO",
                "hydrogenated_path": hydrogenated_path,
            }
        )

    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "status": "SUCCEEDED",
                "request": {"protein": protein},
                "models": models,
                "candidate_count": len(candidates),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return run_dir


# --------------------------------------------------------------------------
# The hydrogen invariant
# --------------------------------------------------------------------------


def test_a_ligand_without_hydrogens_is_never_exported(tmp_path):
    """The measured silent failure: parameterization accepts this and emits a
    hydrogen-free topology. The manifest says hydrogenation happened; the file
    disagrees, and the file is what gets built."""
    run_dir = _write_run(
        tmp_path,
        [
            {"model": "flowr", "status": "PASS", "hydrogens": True},
            {"model": "flowr", "status": "PASS", "hydrogens": "no_hydrogens"},
        ],
    )

    manifest = export_md_inputs(run_dir, tmp_path / "md_inputs")

    assert manifest["exported_count"] == 1
    assert manifest["skipped_reasons"] == {"NO_EXPLICIT_HYDROGENS": 1}
    exported = list((tmp_path / "md_inputs" / "ligands").iterdir())
    assert [path.name for path in exported] == ["flowr_000001.sdf"]
    assert "H " in exported[0].read_text(encoding="utf-8")


def test_every_exported_ligand_carries_hydrogens(tmp_path):
    run_dir = _write_run(
        tmp_path, [{"model": "molcraft", "status": "PASS"} for _ in range(3)]
    )

    manifest = export_md_inputs(run_dir, tmp_path / "md_inputs")

    assert manifest["exported_count"] == 3
    for ligand in manifest["ligands"]:
        assert ligand["hydrogen_count"] == 4
        assert count_explicit_hydrogens(Path(ligand["ligand_path"])) == 4


def test_a_run_generated_without_hydrogens_exports_nothing_and_says_why(tmp_path):
    """``prism generate --no-hydrogens`` leaves nothing this layer can hand on."""
    run_dir = _write_run(
        tmp_path,
        [{"model": "targetdiff", "status": "PASS", "hydrogens": "absent_record"}],
    )

    manifest = export_md_inputs(run_dir, tmp_path / "md_inputs")

    assert manifest["exported_count"] == 0
    assert manifest["skipped_reasons"] == {"NO_HYDROGEN_FILE": 1}
    detail = (tmp_path / "md_inputs" / "skipped.csv").read_text(encoding="utf-8")
    assert "--no-hydrogens" in detail


def test_a_hydrogen_file_recorded_but_missing_is_reported_not_ignored(tmp_path):
    run_dir = _write_run(
        tmp_path,
        [{"model": "diffsbdd", "status": "PASS", "hydrogens": "missing_file"}],
    )

    manifest = export_md_inputs(run_dir, tmp_path / "md_inputs")

    assert manifest["exported_count"] == 0
    assert manifest["skipped_reasons"] == {"NO_HYDROGEN_FILE": 1}


def test_hydrogen_counting_falls_back_to_the_molfile_when_rdkit_is_absent(
    tmp_path, monkeypatch
):
    """RDKit is optional. Refusing to prepare MD inputs without it would be
    worse than reading the columns the format already guarantees."""
    import builtins

    real_import = builtins.__import__

    def no_rdkit(name, *args, **kwargs):
        if name == "rdkit" or name.startswith("rdkit."):
            raise ImportError("rdkit disabled for this test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_rdkit)

    with_h = tmp_path / "with_h.sdf"
    with_h.write_text(WITH_HYDROGENS, encoding="utf-8")
    without_h = tmp_path / "without_h.sdf"
    without_h.write_text(HEAVY_ATOMS_ONLY, encoding="utf-8")

    assert count_explicit_hydrogens(with_h) == 4
    assert count_explicit_hydrogens(without_h) == 0


# --------------------------------------------------------------------------
# Quality-based selection
# --------------------------------------------------------------------------


def test_quality_control_failures_are_never_handed_to_md(tmp_path):
    """Quarantined originals built successfully in 0 of the measured cases."""
    run_dir = _write_run(
        tmp_path,
        [
            {"model": "diffsbdd", "status": "PASS"},
            {"model": "diffsbdd", "status": "FAIL_CHEMISTRY"},
            {"model": "diffsbdd", "status": "FAIL_STRUCTURE"},
        ],
    )

    selected, skipped = select_for_md(run_dir)

    assert [candidate.audit_status for candidate in selected] == ["PASS"]
    assert {item.reason for item in skipped} == {"FAILED_QC"}


def test_warnings_are_advisory_and_still_reach_md(tmp_path):
    """A warning describes geometry or plausibility, not invalid chemistry;
    withholding those molecules would discard most of every run."""
    run_dir = _write_run(
        tmp_path,
        [
            {"model": "pocketxmol", "status": "PASS"},
            {"model": "pocketxmol", "status": "WARN_GEOMETRY"},
            {"model": "pocketxmol", "status": "WARN_PLAUSIBILITY"},
            {"model": "pocketxmol", "status": "WARN_POSE"},
        ],
    )

    selected, skipped = select_for_md(run_dir)

    assert len(selected) == 4
    assert skipped == []


def test_only_pass_narrows_the_selection_and_records_the_rest(tmp_path):
    run_dir = _write_run(
        tmp_path,
        [
            {"model": "pocketxmol", "status": "PASS"},
            {"model": "pocketxmol", "status": "WARN_GEOMETRY"},
        ],
    )

    selected, skipped = select_for_md(run_dir, only_pass=True)

    assert [candidate.audit_status for candidate in selected] == ["PASS"]
    assert [item.reason for item in skipped] == ["NOT_PASS"]


def test_limit_keeps_the_best_grades_not_an_arbitrary_slice(tmp_path):
    """Truncating collection order would deliver whichever model was collected
    first. QC grade predicts build success, so the ranking has to lead."""
    run_dir = _write_run(
        tmp_path,
        [
            {"model": "aaa", "status": "WARN_POSE"},
            {"model": "aaa", "status": "WARN_GEOMETRY"},
            {"model": "zzz", "status": "PASS"},
        ],
    )

    selected, skipped = select_for_md(run_dir, limit=1)

    assert [candidate.audit_status for candidate in selected] == ["PASS"]
    assert [candidate.model for candidate in selected] == ["zzz"]
    assert {item.reason for item in skipped} == {"OVER_LIMIT"}


def test_flag_count_breaks_ties_within_a_grade(tmp_path):
    run_dir = _write_run(
        tmp_path,
        [
            {"model": "m", "status": "WARN_GEOMETRY", "flags": ["a", "b", "c"]},
            {"model": "m", "status": "WARN_GEOMETRY", "flags": ["a"]},
        ],
    )

    selected, _ = select_for_md(run_dir, limit=1)

    assert selected[0].quality_flags == ["a"]


def test_model_filter_selects_one_model_and_reports_the_others(tmp_path):
    run_dir = _write_run(
        tmp_path,
        [
            {"model": "flowr", "status": "PASS"},
            {"model": "targetdiff", "status": "PASS"},
        ],
    )

    selected, skipped = select_for_md(run_dir, models=["flowr"])

    assert [candidate.model for candidate in selected] == ["flowr"]
    assert [item.reason for item in skipped] == ["MODEL_NOT_SELECTED"]


# --------------------------------------------------------------------------
# Output contract
# --------------------------------------------------------------------------


def test_the_generation_run_is_left_untouched(tmp_path):
    """The hand-off copies; it never edits a molecule in place."""
    run_dir = _write_run(tmp_path, [{"model": "flowr", "status": "PASS"}])
    before = {
        path: path.read_bytes()
        for path in sorted(run_dir.rglob("*"))
        if path.is_file()
    }

    export_md_inputs(run_dir, tmp_path / "md_inputs")

    after = {
        path: path.read_bytes()
        for path in sorted(run_dir.rglob("*"))
        if path.is_file()
    }
    assert before == after


def test_the_build_command_names_the_protein_and_every_exported_ligand(tmp_path):
    run_dir = _write_run(
        tmp_path,
        [
            {"model": "flowr", "status": "PASS"},
            {"model": "flowr", "status": "PASS"},
        ],
        protein="/data/1iep.pdb",
    )

    manifest = export_md_inputs(run_dir, tmp_path / "md_inputs")

    command = manifest["build_command"]
    assert "-pf /data/1iep.pdb" in command
    assert command.count("-lf ") == 2
    for ligand in manifest["ligands"]:
        assert ligand["ligand_path"] in command


def test_md_manifest_records_the_selection_criteria(tmp_path):
    run_dir = _write_run(tmp_path, [{"model": "flowr", "status": "PASS"}])

    export_md_inputs(run_dir, tmp_path / "md_inputs", only_pass=True, limit=5)

    manifest = json.loads(
        (tmp_path / "md_inputs" / "md_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["selection"] == {"models": "all", "only_pass": True, "limit": 5}
    assert manifest["generation_run"] == str(run_dir.resolve())


def test_an_existing_export_is_not_overwritten_without_saying_so(tmp_path):
    run_dir = _write_run(tmp_path, [{"model": "flowr", "status": "PASS"}])
    output = tmp_path / "md_inputs"
    export_md_inputs(run_dir, output)

    with pytest.raises(GenerationError) as error:
        export_md_inputs(run_dir, output)
    assert error.value.code == "MD_EXPORT_EXISTS"

    manifest = export_md_inputs(run_dir, output, overwrite=True)
    assert manifest["exported_count"] == 1


def test_a_directory_that_is_not_a_generation_run_is_refused_clearly(tmp_path):
    (tmp_path / "not-a-run").mkdir()

    with pytest.raises(GenerationError) as error:
        select_for_md(tmp_path / "not-a-run")
    assert error.value.code == "NO_GENERATION_MANIFEST"


def test_a_corrupt_manifest_is_refused_clearly(tmp_path):
    run_dir = tmp_path / "broken"
    run_dir.mkdir()
    (run_dir / "manifest.json").write_text("{not json", encoding="utf-8")

    with pytest.raises(GenerationError) as error:
        select_for_md(run_dir)
    assert error.value.code == "INVALID_GENERATION_MANIFEST"


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def test_prepare_md_cli_exports_and_prints_the_build_command(tmp_path, capsys):
    from prism.cli import main

    run_dir = _write_run(
        tmp_path,
        [
            {"model": "flowr", "status": "PASS"},
            {"model": "flowr", "status": "FAIL_CHEMISTRY"},
        ],
    )
    output = tmp_path / "md_inputs"

    exit_code = main(["prepare-md", str(run_dir), "-o", str(output)])

    assert exit_code == 0
    captured = capsys.readouterr().out
    assert "Exported 1 MD-ready ligand(s)" in captured
    assert "FAILED_QC=1" in captured
    assert "prism -pf" in captured
    assert (output / "ligands" / "flowr_000001.sdf").is_file()
    assert (output / "md_manifest.json").is_file()


def test_prepare_md_cli_reports_failure_when_nothing_is_exportable(tmp_path, capsys):
    from prism.cli import main

    run_dir = _write_run(
        tmp_path, [{"model": "flowr", "status": "PASS", "hydrogens": "no_hydrogens"}]
    )

    exit_code = main(["prepare-md", str(run_dir), "-o", str(tmp_path / "md_inputs")])

    assert exit_code == 1
    assert "NO_EXPLICIT_HYDROGENS=1" in capsys.readouterr().out
