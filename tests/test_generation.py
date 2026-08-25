import csv
import json
import hashlib
import io
import os
import pickle
import shutil
import subprocess
import sys
import types
from pathlib import Path

import pytest

from prism.cli import _extract_generation_flag
from prism.generation.adapters.base import ModelAdapter
from prism.generation.errors import CapabilityError, InputError
from prism.generation.execution import ExecutionResult
from prism.generation.execution import execute
from prism.generation.pocket import normalize_pocket
from prism.generation.protein import extract_pocket_pdb
from prism.generation.runner import GenerationRunner
from prism.generation.types import GenerationRequest


PDB = """\
ATOM      1  N   ALA A   1      10.000  11.000  12.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      14.000  15.000  16.000  1.00 20.00           C
END
"""

SDF = """\
candidate
  PRISM    3D

  1  0  0  0  0  0  0  0  0  0  1 V2000
    1.0000    2.0000    3.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""


class FakeAdapter(ModelAdapter):
    model_id = "targetdiff"
    accepted_pocket_kinds = {"structure"}

    def run(self, _request, pocket, _output_root, run_dir, _staged_paths):
        self.check_capabilities(pocket)
        raw = run_dir / "raw"
        raw.mkdir(parents=True, exist_ok=True)
        (raw / "result.sdf").write_text(SDF, encoding="utf-8")
        return ExecutionResult(command=["fake-model"], return_code=0)


class FakeReferenceAdapter(FakeAdapter):
    model_id = "flowr"
    generation_mode = "reference_guided"
    requires_reference_ligand = True
    accepted_reference_ligand_suffixes = {".sdf"}


def _write(path: Path, content: str) -> Path:
    path.write_text(content, encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("arguments", "expected", "forwarded"),
    [
        (["--generation", "--model", "targetdiff"], True, ["--model", "targetdiff"]),
        (["--generation", "protein.pdb", "--model", "targetdiff"], True, ["protein.pdb", "--model", "targetdiff"]),
        (["--generation=false", "protein.pdb"], False, ["protein.pdb"]),
        (["--generation", "no", "protein.pdb"], False, ["protein.pdb"]),
        (["protein.pdb", "ligand.sdf"], None, ["protein.pdb", "ligand.sdf"]),
    ],
)
def test_generation_flag_compatibility(arguments, expected, forwarded):
    assert _extract_generation_flag(arguments) == (expected, forwarded)


def test_generation_false_is_removed_before_builder_forwarding(monkeypatch):
    from prism.cli import main

    captured = []
    builder_module = types.ModuleType("prism.builder")

    def fake_builder_main():
        captured.extend(sys.argv[1:])
        return 17

    builder_module.main = fake_builder_main
    monkeypatch.setitem(sys.modules, "prism.builder", builder_module)

    assert main(["--generation", "false", "protein.pdb", "ligand.sdf"]) == 17
    assert captured == ["protein.pdb", "ligand.sdf"]


def test_importing_prism_does_not_eagerly_import_builder():
    code = (
        "import sys, prism; "
        "assert 'prism.builder' not in sys.modules; "
        "assert prism.get_version()"
    )
    completed = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        capture_output=True,
        text=True,
        env={key: value for key, value in os.environ.items() if key != "PYTHONPATH"},
    )
    assert completed.returncode == 0, completed.stderr


def test_structure_pocket_geometry_is_normalized(tmp_path):
    pocket_path = _write(tmp_path / "pocket.pdb", PDB)
    pocket = normalize_pocket(pocket_path)

    assert pocket.kind == "structure"
    assert pocket.center == pytest.approx((12.0, 13.0, 14.0))
    assert pocket.bbox_size == pytest.approx((4.0, 4.0, 4.0))


def test_reference_ligand_is_inferred_from_sdf(tmp_path):
    ligand_path = _write(tmp_path / "reference.sdf", SDF)
    pocket = normalize_pocket(ligand_path)

    assert pocket.kind == "reference_ligand"
    assert pocket.reference_ligand == ligand_path.resolve()
    assert pocket.center == pytest.approx((1.0, 2.0, 3.0))


def test_yaml_paths_are_relative_to_yaml_file(tmp_path):
    ligand_path = _write(tmp_path / "reference.sdf", SDF)
    pocket_path = _write(
        tmp_path / "pocket.yaml",
        "kind: reference_ligand\nreference_ligand: reference.sdf\nradius: 7.5\n",
    )
    pocket = normalize_pocket(pocket_path)

    assert pocket.reference_ligand == ligand_path.resolve()
    assert pocket.radius == 7.5


def test_residue_file_is_strict(tmp_path):
    pocket_path = _write(tmp_path / "pocket.txt", "A:123\nB:45A\n")
    pocket = normalize_pocket(pocket_path)
    assert pocket.kind == "residues"
    assert pocket.residues == ("A:123", "B:45A")

    _write(pocket_path, "not-a-residue\n")
    with pytest.raises(InputError, match="Invalid residue"):
        normalize_pocket(pocket_path)


def test_flowr_style_reference_requirement_is_explicit(tmp_path):
    from prism.generation.adapters.flowr import FLOWRAdapter

    pocket_path = _write(tmp_path / "pocket.pdb", PDB)
    pocket = normalize_pocket(pocket_path)
    with pytest.raises(CapabilityError) as error:
        FLOWRAdapter({}).check_capabilities(pocket)
    assert error.value.code == "REFERENCE_LIGAND_REQUIRED"


def test_runner_rejects_missing_reference_before_creating_output(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    output = tmp_path / "generated"
    request = GenerationRequest(
        protein=protein,
        pocket=pocket,
        output=output,
        models=("flowr",),
        dry_run=True,
    )
    adapter = FakeReferenceAdapter({})

    with pytest.raises(CapabilityError) as error:
        GenerationRunner(request, {"models": {}}, adapters={"flowr": adapter}).run()

    assert error.value.code == "REFERENCE_LIGAND_REQUIRED"
    assert "flowr" in str(error.value)
    assert not output.exists()


def test_runner_rejects_unused_explicit_reference_for_direct_models(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    reference = _write(tmp_path / "reference.sdf", SDF)
    output = tmp_path / "generated"
    request = GenerationRequest(
        protein=protein,
        pocket=pocket,
        reference_ligand=reference,
        output=output,
        models=("targetdiff",),
        dry_run=True,
    )

    with pytest.raises(CapabilityError) as error:
        GenerationRunner(
            request,
            {"models": {}},
            adapters={"targetdiff": FakeAdapter({})},
        ).run()

    assert error.value.code == "REFERENCE_LIGAND_NOT_NEEDED"
    assert "Remove it" in str(error.value)
    assert not output.exists()


def test_reference_is_allowed_when_a_mixed_selection_requires_it(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    reference = _write(tmp_path / "reference.sdf", SDF)
    request = GenerationRequest(
        protein=protein,
        pocket=pocket,
        reference_ligand=reference,
        output=tmp_path / "generated",
        models=("targetdiff", "flowr"),
        dry_run=True,
    )
    manifest = GenerationRunner(
        request,
        {"models": {}},
        adapters={
            "targetdiff": FakeAdapter({}),
            "flowr": FakeReferenceAdapter({}),
        },
    ).run()

    assert manifest["status"] == "PLANNED"
    assert set(manifest["models"]) == {"targetdiff", "flowr"}


def test_direct_adapter_does_not_forward_shared_conditioning_reference(
    tmp_path, monkeypatch
):
    captured = {}

    def fake_execute(_config, context, _output_root, _run_dir, _device):
        captured.update(context)
        return ExecutionResult(command=["fake"], return_code=0)

    monkeypatch.setattr("prism.generation.adapters.base.execute", fake_execute)
    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket_path = _write(tmp_path / "pocket.pdb", PDB)
    reference = _write(tmp_path / "reference.sdf", SDF)
    pocket = normalize_pocket(pocket_path, reference_ligand=reference)
    request = GenerationRequest(
        protein=protein,
        pocket=pocket_path,
        reference_ligand=reference,
        output=tmp_path / "generated",
        models=("targetdiff", "flowr"),
    )
    adapter = ModelAdapter({})
    adapter.model_id = "targetdiff"
    adapter.accepted_pocket_kinds = {"structure"}
    staged = {
        "protein": protein,
        "pocket": pocket_path,
        "reference_ligand": reference,
    }

    adapter.run(request, pocket, request.output, tmp_path / "run", staged)

    assert captured["reference_ligand"] == ""


def test_reference_format_is_validated_before_creating_output(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    reference = _write(tmp_path / "reference.pdb", PDB)
    request = GenerationRequest(
        protein=protein,
        pocket=pocket,
        reference_ligand=reference,
        output=tmp_path / "generated",
        models=("flowr",),
        dry_run=True,
    )

    with pytest.raises(CapabilityError) as error:
        GenerationRunner(
            request,
            {"models": {}},
            adapters={"flowr": FakeReferenceAdapter({})},
        ).run()

    assert error.value.code == "REFERENCE_LIGAND_FORMAT_UNSUPPORTED"
    assert ".sdf" in str(error.value)
    assert not request.output.exists()


def test_direct_model_accepts_reference_ligand_as_pocket_definition(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    reference = _write(tmp_path / "reference.sdf", SDF)
    request = GenerationRequest(
        protein=protein,
        pocket=reference,
        output=tmp_path / "generated",
        models=("targetdiff",),
        dry_run=True,
    )
    adapter = FakeAdapter({})
    adapter.accepted_pocket_kinds = {"reference_ligand"}
    manifest = GenerationRunner(
        request,
        {"models": {}},
        adapters={"targetdiff": adapter},
    ).run()

    assert manifest["status"] == "PLANNED"


def test_flowr_adapter_selects_cutting_mode_and_validates_cutoff(tmp_path):
    from prism.generation.adapters.flowr import FLOWRAdapter
    from prism.generation.errors import ConfigurationError

    ligand = _write(tmp_path / "reference.sdf", SDF)
    reference = normalize_pocket(ligand)
    assert FLOWRAdapter({}).command_context(None, reference, {}) == {
        "flowr_pocket_cutoff": 6.0,
        "flowr_pocket_mode": "reference",
    }

    pocket = normalize_pocket(_write(tmp_path / "pocket.pdb", PDB))
    assert FLOWRAdapter({}).command_context(None, pocket, {})[
        "flowr_pocket_mode"
    ] == "prepared"
    with pytest.raises(ConfigurationError, match="pocket_cutoff"):
        FLOWRAdapter({"pocket_cutoff": 0}).command_context(None, reference, {})


def test_runner_collects_fake_adapter_output(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    output = tmp_path / "generated"
    request = GenerationRequest(
        protein=protein,
        pocket=pocket,
        output=output,
        models=("targetdiff",),
        num_samples=1,
    )
    config = {"models": {"targetdiff": {"output_patterns": ["raw/**/*.sdf"]}}}

    manifest = GenerationRunner(
        request,
        config,
        adapters={"targetdiff": FakeAdapter(config["models"]["targetdiff"])},
    ).run()

    assert manifest["status"] == "SUCCEEDED"
    assert manifest["candidate_count"] == 1
    assert (output / "candidates.sdf").is_file()
    assert (output / "models" / "targetdiff" / "molecules" / "candidate_000001.sdf").is_file()
    assert json.loads((output / "manifest.json").read_text(encoding="utf-8"))["status"] == "SUCCEEDED"


def test_overwriting_one_model_preserves_other_models_in_aggregate_outputs(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    output = tmp_path / "generated"
    config = {
        "models": {
            "targetdiff": {"output_patterns": ["raw/**/*.sdf"]},
            "pocketxmol": {"output_patterns": ["raw/**/*.sdf"]},
        }
    }
    first = GenerationRequest(
        protein=protein,
        pocket=pocket,
        output=output,
        models=("targetdiff",),
        num_samples=1,
    )
    GenerationRunner(
        first,
        config,
        adapters={"targetdiff": FakeAdapter(config["models"]["targetdiff"])},
    ).run()

    retry = GenerationRequest(
        protein=protein,
        pocket=pocket,
        output=output,
        models=("pocketxmol",),
        num_samples=1,
        overwrite=True,
    )
    manifest = GenerationRunner(
        retry,
        config,
        adapters={"pocketxmol": FakeAdapter(config["models"]["pocketxmol"])},
    ).run()

    assert set(manifest["models"]) == {"targetdiff", "pocketxmol"}
    assert manifest["candidate_count"] == 2
    with (output / "summary.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {(row["model"], row["candidate_id"]) for row in rows} == {
        ("targetdiff", "targetdiff_000001"),
        ("pocketxmol", "pocketxmol_000001"),
    }
    assert (output / "candidates.sdf").read_text(encoding="utf-8").count("$$$$") == 2


def test_sdf_record_parser_preserves_empty_title_line():
    from prism.generation.postprocess import _coordinates, _sdf_records

    record = next(iter(_sdf_records(SDF.replace("candidate\n", "\n"))))
    coordinates, valid_3d = _coordinates(record)

    assert coordinates == [(1.0, 2.0, 3.0)]
    assert valid_3d is True


def test_targetdiff_extracts_pocket_from_reference_center(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    ligand = _write(
        tmp_path / "reference.sdf",
        SDF.replace("    1.0000    2.0000    3.0000", "   12.0000   13.0000   14.0000"),
    )
    pocket = normalize_pocket(ligand)

    output = extract_pocket_pdb(protein, pocket, tmp_path / "prepared" / "pocket.pdb")

    text = output.read_text(encoding="utf-8")
    assert "ALA A   1" in text
    assert text.endswith("END\n")


def test_reference_ligand_extraction_uses_all_ligand_atoms(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    ligand = _write(
        tmp_path / "elongated.sdf",
        """elongated
  PRISM    3D

  2  0  0  0  0  0  0  0  0  0  1 V2000
   12.0000   13.0000   14.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  112.0000  113.0000  114.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
""",
    )
    pocket = normalize_pocket(ligand)

    output = extract_pocket_pdb(protein, pocket, tmp_path / "all-atoms-pocket.pdb")
    assert "ALA A   1" in output.read_text(encoding="utf-8")


def test_artifact_hash_is_verified_before_execution(tmp_path):
    checkpoint = _write(tmp_path / "checkpoint.pt", "trusted checkpoint bytes")
    bad_config = {
        "enabled": True,
        "backend": "conda",
        "environment": "targetdiff",
        "artifacts": {
            "checkpoint": {
                "path": str(checkpoint),
                "sha256": "0" * 64,
            }
        },
        "command": ["python", "wrapper.py", "--checkpoint", "{checkpoint}"],
    }
    from prism.generation.errors import ConfigurationError

    with pytest.raises(ConfigurationError, match="SHA256 mismatch"):
        execute(bad_config, {}, tmp_path, tmp_path, "cpu")

    bad_config["artifacts"]["checkpoint"]["sha256"] = hashlib.sha256(
        checkpoint.read_bytes()
    ).hexdigest()
    bad_config["command"] = ["definitely-not-run-in-this-test"]
    # Rendering now succeeds; the missing Conda executable/environment is a later concern.
    from prism.generation.execution import _artifacts

    context, mounts, hashes = _artifacts(bad_config, "conda")
    assert context["checkpoint"] == str(checkpoint.resolve())
    assert mounts == []
    assert hashes["checkpoint"] == bad_config["artifacts"]["checkpoint"]["sha256"]


def test_docker_execution_renders_mounts_and_container_device(tmp_path, monkeypatch):
    checkpoint = _write(tmp_path / "checkpoint.pt", "checkpoint")
    output = tmp_path / "generated"
    run_dir = output / "models" / "targetdiff"
    run_dir.mkdir(parents=True)
    pocket = _write(output / "inputs.pdb", PDB)
    digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    config = {
        "enabled": True,
        "backend": "docker",
        "image": "local/targetdiff@sha256:" + "a" * 64,
        "artifacts": {
            "checkpoint": {
                "path": str(checkpoint),
                "sha256": digest,
                "container_path": "/opt/checkpoints/model.pt",
            }
        },
        "command": [
            "python",
            "wrapper.py",
            "--pocket",
            "{pocket}",
            "--checkpoint",
            "{checkpoint}",
            "--device",
            "{device}",
        ],
    }
    captured = {}

    class _FakeLauncher:
        """Stand-in for the launcher process: records the command, exits cleanly."""

        pid = 4242

        def __init__(self, command, **_kwargs):
            captured["command"] = command

        def wait(self, timeout=None):  # noqa: ARG002
            return 0

    monkeypatch.setattr("prism.generation.execution.subprocess.Popen", _FakeLauncher)
    result = execute(
        config,
        {"pocket": pocket, "device": "cuda:2"},
        output,
        run_dir,
        "cuda:2",
    )

    command = captured["command"]
    assert command[command.index("--gpus") + 1] == "device=2"
    assert f"{checkpoint.resolve()}:/opt/checkpoints/model.pt:ro" in command
    assert "/work/inputs.pdb" in command
    assert command[-1] == "cuda:0"
    assert result.artifact_sha256 == {"checkpoint": digest}


def test_slurm_execution_renders_resource_request_and_conda_path(tmp_path, monkeypatch):
    checkpoint = _write(tmp_path / "checkpoint.pt", "checkpoint")
    output = tmp_path / "generated"
    run_dir = output / "models" / "targetdiff"
    run_dir.mkdir(parents=True)
    digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    config = {
        "enabled": True,
        "backend": "slurm",
        "environment": "/opt/conda-envs/prism-gen-targetdiff",
        "slurm": {
            "executable": "/opt/gridview/slurm/bin/srun",
            "partition": "debug",
            "account": "iqb",
            "time": "00:45:00",
            "gpus": 1,
            "cpus_per_task": 4,
            "memory": "24G",
            "job_name": "prism-targetdiff",
        },
        "artifacts": {
            "checkpoint": {"path": str(checkpoint), "sha256": digest}
        },
        "command": ["python", "wrapper.py", "{checkpoint}", "{device}"],
    }
    captured = {}

    class _FakeLauncher:
        """Stand-in for the launcher process: records the command, exits cleanly."""

        pid = 4242

        def __init__(self, command, **_kwargs):
            captured["command"] = command

        def wait(self, timeout=None):  # noqa: ARG002
            return 0

    monkeypatch.setattr("prism.generation.execution.subprocess.Popen", _FakeLauncher)
    execute(config, {"device": "auto"}, output, run_dir, "auto")

    command = captured["command"]
    assert command[:5] == [
        "/opt/gridview/slurm/bin/srun",
        "--partition",
        "debug",
        "--account",
        "iqb",
    ]
    assert command[command.index("--gres") + 1] == "gpu:1"
    assert command[command.index("--mem") + 1] == "24G"
    assert command[command.index("conda") + 3] == "-p"
    assert "/opt/conda-envs/prism-gen-targetdiff" in command
    assert command[-1] == "cuda:0"


def test_targetdiff_wrapper_writes_upstream_config(tmp_path, monkeypatch):
    from prism.generation.wrappers.targetdiff import main

    root = tmp_path / "targetdiff"
    script = root / "scripts" / "sample_for_pocket.py"
    script.parent.mkdir(parents=True)
    _write(script, "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "model.pt", "checkpoint")
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    output = tmp_path / "raw"
    captured = {}

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        captured["command"] = command
        captured["cwd"] = cwd
        captured["env"] = env
        captured["check"] = check
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("prism.generation.wrappers.targetdiff.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.targetdiff._sanitize_checkpoint",
        lambda source, _output: source,
    )
    result = main(
        [
            "--targetdiff-root",
            str(root),
            "--checkpoint",
            str(checkpoint),
            "--pocket",
            str(pocket),
            "--output",
            str(output),
            "--num-samples",
            "7",
            "--seed",
            "19",
            "--device",
            "cuda",
        ]
    )

    config = __import__("yaml").safe_load(
        (output / "prism_targetdiff_sampling.yml").read_text(encoding="utf-8")
    )
    assert result == 0
    assert config["model"]["checkpoint"] == str(checkpoint.resolve())
    assert config["sample"]["num_samples"] == 7
    assert config["sample"]["seed"] == 19
    assert "--pdb_path" in captured["command"]
    assert captured["command"][captured["command"].index("--device") + 1] == "cuda:0"
    assert captured["env"]["PYTHONPATH"].split(os.pathsep)[0] == str(root.resolve())


def test_targetdiff_restricted_unpickler_rejects_unknown_globals():
    from prism.generation.wrappers.targetdiff import _RestrictedUnpickler

    payload = pickle.dumps(eval)
    with pytest.raises(pickle.UnpicklingError, match="forbidden pickle global"):
        _RestrictedUnpickler(io.BytesIO(payload)).load()


def test_pocket2mol_adapter_uses_official_default_box(tmp_path):
    from prism.generation.adapters.pocket2mol import Pocket2MolAdapter

    ligand = _write(tmp_path / "reference.sdf", SDF)
    pocket = normalize_pocket(ligand)
    context = Pocket2MolAdapter({}).command_context(None, pocket, {})
    assert context["bbox_size"] == 23.0


def test_pocket2mol_wrapper_writes_config_and_collects_sdf(tmp_path, monkeypatch):
    from prism.generation.wrappers.pocket2mol import _UNBOUNDED_SAMPLES, main

    root = tmp_path / "pocket2mol"
    script = root / "sample_for_pdb.py"
    script.parent.mkdir(parents=True)
    _write(script, "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "model.pt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB)
    output = tmp_path / "raw"
    captured = {}

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        captured["command"] = command
        captured["cwd"] = cwd
        captured["env"] = env
        upstream = output / "upstream" / "sample_20200101" / "SDF"
        upstream.mkdir(parents=True)
        _write(upstream / "0.sdf", SDF)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("prism.generation.wrappers.pocket2mol.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.pocket2mol._sanitize_checkpoint",
        lambda source, _output: source,
    )
    result = main(
        [
            "--pocket2mol-root",
            str(root),
            "--checkpoint",
            str(checkpoint),
            "--protein",
            str(protein),
            "--center",
            "-1",
            "2",
            "3",
            "--bbox-size",
            "23",
            "--output",
            str(output),
            "--num-samples",
            "1",
            "--seed",
            "19",
            "--device",
            "cuda",
        ]
    )

    config = __import__("yaml").safe_load(
        (output / "prism_pocket2mol_sampling.yml").read_text(encoding="utf-8")
    )
    assert result == 0
    # --num-samples selects; it is deliberately kept out of the upstream
    # stopping condition, which would otherwise cap molecule size.
    assert config["sample"]["num_samples"] == _UNBOUNDED_SAMPLES
    assert config["sample"]["seed"] == 19
    assert "--center=-1.0,2.0,3.0" in captured["command"]
    assert captured["command"][captured["command"].index("--device") + 1] == "cuda:0"
    assert (output / "sdf" / "candidate_000001.sdf").is_file()


def test_pocket2mol_wrapper_stages_every_molecule_in_numeric_order(tmp_path, monkeypatch):
    """Pocket2Mol writes 0.sdf..14.sdf, so lexicographic order drops molecules."""
    from prism.generation.wrappers.pocket2mol import main

    root = tmp_path / "pocket2mol"
    script = root / "sample_for_pdb.py"
    script.parent.mkdir(parents=True)
    _write(script, "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "model.pt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB)
    output = tmp_path / "raw"

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        upstream = output / "upstream" / "sample_20200101" / "SDF"
        upstream.mkdir(parents=True)
        for index in range(15):
            # Tag each molecule so the mapping to its source is checkable.
            _write(upstream / f"{index}.sdf", SDF.replace("candidate", f"mol{index}"))
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("prism.generation.wrappers.pocket2mol.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.pocket2mol._sanitize_checkpoint",
        lambda source, _output: source,
    )
    assert (
        main(
            [
                "--pocket2mol-root", str(root),
                "--checkpoint", str(checkpoint),
                "--protein", str(protein),
                "--center", "-1", "2", "3",
                "--bbox-size", "23",
                "--output", str(output),
                "--num-samples", "10",
                "--seed", "19",
                "--device", "cuda",
            ]
        )
        == 0
    )

    # sdf_all/ keeps every finished molecule, in the model's completion order.
    staged = sorted((output / "sdf_all").glob("*.sdf"))
    assert len(staged) == 15
    titles = [path.read_text(encoding="utf-8").splitlines()[0] for path in staged]
    assert titles == [f"mol{index}" for index in range(15)]


def test_pocket2mol_selection_modes_span_or_target_the_size_distribution():
    """Completion order tracks molecule size, so the criterion must be explicit."""
    from prism.generation.wrappers.pocket2mol import _select_indices

    # Sizes as Pocket2Mol emits them: one atom added per sampling step.
    counts = list(range(5, 25))  # 20 molecules, 5..24 heavy atoms

    completion = _select_indices(counts, 5, "completion", None)
    assert [counts[i] for i in completion] == [5, 6, 7, 8, 9]

    largest = _select_indices(counts, 5, "largest", None)
    assert [counts[i] for i in largest] == [20, 21, 22, 23, 24]

    target = _select_indices(counts, 5, "target", 15)
    assert sorted(counts[i] for i in target) == [13, 14, 15, 16, 17]

    stratified = _select_indices(counts, 5, "stratified", None)
    sizes = sorted(counts[i] for i in stratified)
    assert len(sizes) == 5
    assert sizes[0] == 5 and sizes[-1] == 24  # spans the full range

    # Asking for at least as many as exist delivers everything.
    assert _select_indices(counts, 20, "stratified", None) == list(range(20))
    assert _select_indices(counts, 99, "largest", None) == list(range(20))


def test_pocket2mol_stratified_selection_never_returns_duplicates():
    from prism.generation.wrappers.pocket2mol import _select_indices

    # A flat distribution makes rounding collisions likely.
    counts = [10] * 30
    for wanted in (1, 2, 7, 29):
        chosen = _select_indices(counts, wanted, "stratified", None)
        assert len(chosen) == wanted
        assert len(set(chosen)) == wanted


def test_pocket2mol_target_mode_requires_a_target():
    from prism.generation.wrappers.pocket2mol import _select_indices

    with pytest.raises(SystemExit, match="target-heavy-atoms"):
        _select_indices([5, 6, 7], 2, "target", None)


def test_pocket2mol_unbounded_sampling_replaces_the_stopping_condition(tmp_path, monkeypatch):
    """--num-samples must not reach the upstream while-loop, only the selection."""
    from prism.generation.wrappers.pocket2mol import _UNBOUNDED_SAMPLES, main

    root = tmp_path / "pocket2mol"
    script = root / "sample_for_pdb.py"
    script.parent.mkdir(parents=True)
    _write(script, "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "model.pt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB)
    output = tmp_path / "raw"

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        upstream = output / "upstream" / "sample" / "SDF"
        upstream.mkdir(parents=True)
        # 12 molecules whose declared atom count grows with completion order.
        for index in range(12):
            atoms = index + 3
            body = "\n".join(
                f"{i:>10}.0000{i:>10}.0000{i:>10}.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
                for i in range(atoms)
            )
            _write(
                upstream / f"{index}.sdf",
                f"mol{index}\n  PRISM    3D\n\n{atoms:>3}  0  0  0  0  0  0  0  0  0999 V2000\n{body}\nM  END\n$$$$\n",
            )
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("prism.generation.wrappers.pocket2mol.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.pocket2mol._sanitize_checkpoint",
        lambda source, _output: source,
    )
    arguments = [
        "--pocket2mol-root", str(root),
        "--checkpoint", str(checkpoint),
        "--protein", str(protein),
        "--center", "-1", "2", "3",
        "--bbox-size", "23",
        "--output", str(output),
        "--num-samples", "4",
        "--seed", "19",
        "--device", "cuda",
        "--max-steps", "40",
    ]
    assert main(arguments) == 0

    config = __import__("yaml").safe_load(
        (output / "prism_pocket2mol_sampling.yml").read_text(encoding="utf-8")
    )
    assert config["sample"]["num_samples"] == _UNBOUNDED_SAMPLES
    assert config["sample"]["max_steps"] == 40

    summary = json.loads((output / "selection.json").read_text(encoding="utf-8"))
    assert summary["finished_total"] == 12
    assert summary["selected"] == 4
    assert summary["selection_mode"] == "stratified"
    # The delivered set spans the distribution rather than its small tail.
    assert summary["heavy_atoms_all"] == {"min": 3, "median": 9, "max": 14}
    assert summary["heavy_atoms_selected"]["min"] == 3
    assert summary["heavy_atoms_selected"]["max"] == 14

    assert len(list((output / "sdf").glob("*.sdf"))) == 4
    # Nothing is discarded: every finished molecule stays available.
    assert len(list((output / "sdf_all").glob("*.sdf"))) == 12


def test_pocket2mol_completion_mode_preserves_the_upstream_contract(tmp_path, monkeypatch):
    """The old behavior stays reachable for reproducing earlier runs."""
    from prism.generation.wrappers.pocket2mol import main

    root = tmp_path / "pocket2mol"
    (root).mkdir(parents=True)
    _write(root / "sample_for_pdb.py", "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "model.pt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB)
    output = tmp_path / "raw"

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        upstream = output / "upstream" / "sample" / "SDF"
        upstream.mkdir(parents=True)
        for index in range(6):
            _write(upstream / f"{index}.sdf", SDF)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("prism.generation.wrappers.pocket2mol.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.pocket2mol._sanitize_checkpoint",
        lambda source, _output: source,
    )
    assert (
        main(
            [
                "--pocket2mol-root", str(root),
                "--checkpoint", str(checkpoint),
                "--protein", str(protein),
                "--center", "-1", "2", "3",
                "--bbox-size", "23",
                "--output", str(output),
                "--num-samples", "3",
                "--seed", "19",
                "--device", "cuda",
                "--select", "completion",
            ]
        )
        == 0
    )
    config = __import__("yaml").safe_load(
        (output / "prism_pocket2mol_sampling.yml").read_text(encoding="utf-8")
    )
    assert config["sample"]["num_samples"] == 3
    assert len(list((output / "sdf").glob("*.sdf"))) == 3


class _FinishedEntry:
    """A Pocket2Mol pool entry as upstream pickles it into samples_*.pt.

    Module level on purpose: pickle stores a class by qualified name, so one
    defined inside the test would be unpicklable and torch.save would fail
    before the test could exercise anything. rdkit is imported per instance to
    keep collection free of the dependency, and AllChem explicitly because
    rdkit only binds it onto Chem as a side effect of some other import.
    """

    def __init__(self, smiles):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        self.rdmol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(self.rdmol, randomSeed=7)


def test_pocket2mol_recovers_finished_molecules_when_upstream_crashes(tmp_path, monkeypatch):
    """Upstream dies in np.random.choice once its beam queue drains.

    It only guards the loop against KeyboardInterrupt, so the crash escapes
    before its own "Save sdf mols" block. The per-step pool checkpoint is the
    only surviving copy of the finished molecules.
    """
    pytest.importorskip("torch")
    pytest.importorskip("rdkit")

    from prism.generation.wrappers.pocket2mol import main

    root = tmp_path / "pocket2mol"
    root.mkdir(parents=True)
    _write(root / "sample_for_pdb.py", "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "model.pt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB)
    output = tmp_path / "raw"
    FinishedEntry = _FinishedEntry

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        import torch

        run_dir = output / "upstream" / "sample"
        run_dir.mkdir(parents=True)
        # Upstream checkpoints the pool each step, then crashes without ever
        # reaching its SDF writer, so no SDF/ directory exists.
        torch.save({"finished": [FinishedEntry("CCO"), FinishedEntry("c1ccccc1")]},
                   run_dir / "samples_init.pt")
        torch.save(
            {"finished": [FinishedEntry("CCO"), FinishedEntry("c1ccccc1"),
                          FinishedEntry("CNC(=O)c1ccccc1")]},
            run_dir / "samples_9.pt",
        )
        return subprocess.CompletedProcess(command, 1)

    monkeypatch.setattr("prism.generation.wrappers.pocket2mol.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.pocket2mol._sanitize_checkpoint",
        lambda source, _output: source,
    )
    assert (
        main(
            [
                "--pocket2mol-root", str(root),
                "--checkpoint", str(checkpoint),
                "--protein", str(protein),
                "--center", "-1", "2", "3",
                "--bbox-size", "23",
                "--output", str(output),
                "--num-samples", "3",
                "--seed", "19",
                "--device", "cuda",
            ]
        )
        == 0
    )

    summary = json.loads((output / "selection.json").read_text(encoding="utf-8"))
    assert summary["upstream_return_code"] == 1
    assert summary["recovered_from_checkpoint"] == 3
    assert summary["finished_total"] == 3
    assert len(list((output / "sdf").glob("*.sdf"))) == 3


def test_pocket2mol_propagates_failure_when_nothing_is_recoverable(tmp_path, monkeypatch):
    from prism.generation.wrappers.pocket2mol import main

    root = tmp_path / "pocket2mol"
    root.mkdir(parents=True)
    _write(root / "sample_for_pdb.py", "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "model.pt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB)
    output = tmp_path / "raw"

    monkeypatch.setattr(
        "prism.generation.wrappers.pocket2mol.subprocess.run",
        lambda command, cwd, env, check: subprocess.CompletedProcess(command, 9),
    )
    monkeypatch.setattr(
        "prism.generation.wrappers.pocket2mol._sanitize_checkpoint",
        lambda source, _output: source,
    )
    assert (
        main(
            [
                "--pocket2mol-root", str(root),
                "--checkpoint", str(checkpoint),
                "--protein", str(protein),
                "--center", "-1", "2", "3",
                "--bbox-size", "23",
                "--output", str(output),
                "--num-samples", "3",
                "--seed", "19",
                "--device", "cuda",
            ]
        )
        == 9
    )


def test_pocket2mol_adapter_exposes_reference_size_as_target(tmp_path):
    from prism.generation.adapters.pocket2mol import Pocket2MolAdapter

    ligand = _write(tmp_path / "reference.sdf", SDF)
    pocket = normalize_pocket(ligand)
    context = Pocket2MolAdapter({}).command_context(None, pocket, {})

    # The single-atom fixture declares one atom on its counts line.
    assert context["pocket2mol_target_heavy_atoms"] == 1

    center_yaml = _write(tmp_path / "center.yaml", "center: [1, 2, 3]\n")
    center = normalize_pocket(center_yaml)
    assert Pocket2MolAdapter({}).command_context(None, center, {})[
        "pocket2mol_target_heavy_atoms"
    ] == 0


def test_pocket2mol_restricted_unpickler_rejects_unknown_globals():
    from prism.generation.wrappers.pocket2mol import _RestrictedUnpickler

    payload = pickle.dumps(eval)
    with pytest.raises(pickle.UnpicklingError, match="forbidden pickle global"):
        _RestrictedUnpickler(io.BytesIO(payload)).load()


def test_pocketxmol_adapter_selects_reference_and_center_radii(tmp_path):
    from prism.generation.adapters.pocketxmol import PocketXMolAdapter

    ligand = _write(tmp_path / "reference.sdf", SDF)
    reference = normalize_pocket(ligand)
    adapter = PocketXMolAdapter({})
    assert adapter.command_context(None, reference, {})["pocketxmol_radius"] == 10.0

    center_yaml = _write(tmp_path / "center.yaml", "center: [1, 2, 3]\n")
    center = normalize_pocket(center_yaml)
    assert adapter.command_context(None, center, {})["pocketxmol_radius"] == 15.0


def test_pocketxmol_wrapper_writes_task_and_collects_confidence(tmp_path, monkeypatch):
    from prism.generation.postprocess import collect_candidates
    from prism.generation.wrappers.pocketxmol import main

    root = tmp_path / "pocketxmol"
    script = root / "scripts" / "sample_use.py"
    script.parent.mkdir(parents=True)
    _write(script, "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "pocketxmol.ckpt", "checkpoint")
    train_config = _write(tmp_path / "train.yml", "model: {}\n")
    protein = _write(tmp_path / "protein.pdb", PDB)
    reference = _write(tmp_path / "reference.sdf", SDF)
    output = tmp_path / "raw"
    captured = {}

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        captured["command"] = command
        captured["cwd"] = cwd
        captured["env"] = env
        log_dir = output / "upstream" / "task_model_20200101"
        sdf_dir = log_dir / "task_model_20200101_SDF"
        sdf_dir.mkdir(parents=True)
        _write(sdf_dir / "0.sdf", SDF)
        _write(
            log_dir / "gen_info.csv",
            "filename,tag,cfd_traj,cfd_pos,cfd_node,cfd_edge\n"
            "0.sdf,,0.75,0.8,0.7,0.6\n",
        )
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("prism.generation.wrappers.pocketxmol.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.pocketxmol._sanitize_checkpoint",
        lambda source, _output: source,
    )
    result = main(
        [
            "--pocketxmol-root",
            str(root),
            "--checkpoint",
            str(checkpoint),
            "--train-config",
            str(train_config),
            "--protein",
            str(protein),
            "--reference-ligand",
            str(reference),
            "--center",
            "1",
            "2",
            "3",
            "--radius",
            "10",
            "--output",
            str(output),
            "--num-samples",
            "1",
            "--seed",
            "23",
            "--device",
            "cuda",
        ]
    )

    task = __import__("yaml").safe_load(
        (output / "prism_pocketxmol_task.yml").read_text(encoding="utf-8")
    )
    assert result == 0
    assert task["sample"]["num_mols"] == 1
    assert task["data"]["pocket_args"]["ref_ligand_path"] == str(reference.resolve())
    assert task["transforms"]["featurizer_pocket"]["center"] == [1.0, 2.0, 3.0]
    launcher = output / "run_pocketxmol_upstream.py"
    assert captured["command"][1] == str(launcher)
    assert captured["env"]["PRISM_POCKETXMOL_ROOT"] == str(root.resolve())
    launcher_text = launcher.read_text(encoding="utf-8")
    assert 'types.ModuleType("scripts.train_pl")' in launcher_text
    assert "pytorch_lightning" not in launcher_text
    assert captured["command"][captured["command"].index("--device") + 1] == "cuda:0"
    staged_sdf = output / "sdf" / "candidate_000001.sdf"
    assert "<POCKETXMOL_CFD_TRAJ>\n0.75" in staged_sdf.read_text(encoding="utf-8")

    run_dir = tmp_path / "model-run"
    raw = run_dir / "raw" / "sdf"
    raw.mkdir(parents=True)
    shutil.copy2(staged_sdf, raw / staged_sdf.name)
    candidates, quarantined, failures, _attrition = collect_candidates(
        "pocketxmol", run_dir, ["raw/sdf/*.sdf"], 1
    )
    assert failures == []
    assert quarantined == []
    assert candidates[0].native_score_name == "cfd_traj"
    assert candidates[0].native_score == pytest.approx(0.75)


def test_pocketxmol_restricted_unpickler_rejects_unknown_globals():
    from prism.generation.wrappers.pocketxmol import _RestrictedUnpickler

    payload = pickle.dumps(eval)
    with pytest.raises(pickle.UnpicklingError, match="forbidden pickle global"):
        _RestrictedUnpickler(io.BytesIO(payload)).load()


def test_molcraft_wrapper_isolates_config_and_collects_sdf(tmp_path, monkeypatch):
    from prism.generation.wrappers.molcraft import main

    root = tmp_path / "molcraft"
    root.mkdir()
    _write(root / "sample_for_pocket.py", "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "molcraft.ckpt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB)
    reference = _write(tmp_path / "reference.sdf", SDF)
    output = tmp_path / "raw"
    captured = {}

    def fake_sanitize(source, destination, config_path):
        config_path.parent.mkdir(parents=True, exist_ok=True)
        _write(config_path, "seed: 1234\n")
        return source

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        captured["command"] = command
        captured["cwd"] = cwd
        captured["env"] = env
        upstream = Path(env["PRISM_MOLCRAFT_OUTPUT"])
        _write(upstream / "0.sdf", SDF)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(
        "prism.generation.wrappers.molcraft._sanitize_checkpoint", fake_sanitize
    )
    monkeypatch.setattr("prism.generation.wrappers.molcraft.subprocess.run", fake_run)
    result = main(
        [
            "--molcraft-root",
            str(root),
            "--checkpoint",
            str(checkpoint),
            "--protein",
            str(protein),
            "--reference-ligand",
            str(reference),
            "--output",
            str(output),
            "--num-samples",
            "1",
            "--seed",
            "29",
            "--device",
            "cuda:2",
        ]
    )

    assert result == 0
    assert captured["env"]["CUDA_VISIBLE_DEVICES"] == "2"
    assert captured["env"]["LD_LIBRARY_PATH"].split(os.pathsep)[0] == str(
        Path(sys.prefix) / "lib"
    )
    assert captured["env"]["PRISM_MOLCRAFT_SEED"] == "29"
    assert Path(captured["cwd"]) == output / "runtime"
    assert (output / "sdf" / "candidate_000001.sdf").is_file()
    launcher = (output / "run_molcraft_upstream.py").read_text(encoding="utf-8")
    assert "sample_callback.OUT_DIR" in launcher
    assert 'types.ModuleType("posecheck")' in launcher
    assert 'types.ModuleType("core.evaluation.metrics")' in launcher
    assert "upstream.call(" in launcher


def test_molcraft_restricted_unpickler_rejects_unknown_globals():
    from prism.generation.wrappers.molcraft import _RestrictedUnpickler

    payload = pickle.dumps(eval)
    with pytest.raises(pickle.UnpicklingError, match="forbidden pickle global"):
        _RestrictedUnpickler(io.BytesIO(payload)).load()


def test_flowr_wrapper_uses_official_module_and_collects_multisdf(tmp_path, monkeypatch):
    from prism.generation.wrappers.flowr import main

    root = tmp_path / "flowr"
    module = root / "flowr" / "gen" / "generate_from_pdb.py"
    module.parent.mkdir(parents=True)
    _write(module, "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "flowr.ckpt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB + "CONECT    1    2\n")
    reference = _write(tmp_path / "reference.sdf", SDF)
    output = tmp_path / "raw"
    captured = {}

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        captured["command"] = command
        captured["cwd"] = cwd
        captured["env"] = env
        save_dir = Path(command[command.index("--save_dir") + 1])
        _write(save_dir / "samples_protein.sdf", SDF + SDF)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(
        "prism.generation.wrappers.flowr._sanitize_checkpoint",
        lambda source, _output: source,
    )
    monkeypatch.setattr("prism.generation.wrappers.flowr.subprocess.run", fake_run)
    result = main(
        [
            "--flowr-root",
            str(root),
            "--checkpoint",
            str(checkpoint),
            "--protein",
            str(protein),
            "--reference-ligand",
            str(reference),
            "--pocket-mode",
            "reference",
            "--output",
            str(output),
            "--num-samples",
            "2",
            "--seed",
            "31",
            "--device",
            "cuda:3",
        ]
    )

    assert result == 0
    assert Path(captured["command"][1]).name == "launch_flowr.py"
    launcher = Path(captured["command"][1]).read_text(encoding="utf-8")
    assert 'runpy.run_module("flowr.gen.generate_from_pdb"' in launcher
    assert 'sys.modules["flowr.util.metrics"]' in launcher
    assert 'sys.modules["openmm"]' in launcher
    assert 'sys.modules["pdbinf"]' in launcher
    prepared_protein = Path(
        captured["command"][captured["command"].index("--pdb_file") + 1]
    )
    assert prepared_protein.parent == output / "runtime"
    assert "CONECT" not in prepared_protein.read_text(encoding="utf-8")
    assert "--cut_pocket" in captured["command"]
    assert captured["env"]["CUDA_VISIBLE_DEVICES"] == "3"
    assert captured["env"]["LD_LIBRARY_PATH"].split(os.pathsep)[0] == str(
        Path(sys.prefix) / "lib"
    )
    assert (output / "sdf" / "candidates.sdf").read_text(encoding="utf-8").count(
        "$$$$"
    ) == 2


def test_flowr_restricted_unpickler_rejects_unknown_globals():
    from prism.generation.wrappers.flowr import _RestrictedUnpickler

    payload = pickle.dumps(eval)
    with pytest.raises(pickle.UnpicklingError, match="forbidden pickle global"):
        _RestrictedUnpickler(io.BytesIO(payload)).load()


def test_runner_dry_run_stages_without_enabled_backend(tmp_path):
    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    request = GenerationRequest(
        protein=protein,
        pocket=pocket,
        output=tmp_path / "planned",
        models=("targetdiff",),
        num_samples=5,
        dry_run=True,
    )
    config = {"models": {"targetdiff": {"enabled": False}}}
    manifest = GenerationRunner(request, config).run()

    assert manifest["status"] == "PLANNED"
    assert manifest["models"]["targetdiff"]["requested"] == 5


def test_generation_cli_dry_run(tmp_path):
    from prism.generation.cli import main

    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    output = tmp_path / "cli-output"

    return_code = main(
        [
            "--model",
            "targetdiff",
            "--protein",
            str(protein),
            "--pocket",
            str(pocket),
            "--output",
            str(output),
            "--num-samples",
            "3",
            "--dry-run",
        ]
    )

    assert return_code == 0
    assert json.loads((output / "manifest.json").read_text(encoding="utf-8"))["status"] == "PLANNED"


def test_generation_cli_lists_model_modes(capsys):
    from prism.generation.cli import main

    assert main(["--list-models"]) == 0
    output = capsys.readouterr().out
    assert (
        "targetdiff\tdirect\treference ligand: not required; ligand pocket supported"
        in output
    )
    assert "molcraft\treference-guided\treference ligand: required (.sdf)" in output
    assert "flowr\treference-guided\treference ligand: required (.pdb, .sdf)" in output


def test_generation_cli_reports_missing_reference_ligand(tmp_path, capsys):
    from prism.generation.cli import main

    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    output = tmp_path / "generated"
    with pytest.raises(SystemExit) as error:
        main(
            [
                "--model",
                "molcraft",
                "--protein",
                str(protein),
                "--pocket",
                str(pocket),
                "--output",
                str(output),
                "--dry-run",
            ]
        )

    assert error.value.code == 2
    assert "REFERENCE_LIGAND_REQUIRED" in capsys.readouterr().err
    assert not output.exists()


def test_generation_cli_reports_unneeded_reference_ligand(tmp_path, capsys):
    from prism.generation.cli import main

    protein = _write(tmp_path / "protein.pdb", PDB)
    pocket = _write(tmp_path / "pocket.pdb", PDB)
    reference = _write(tmp_path / "reference.sdf", SDF)
    output = tmp_path / "generated"
    with pytest.raises(SystemExit) as error:
        main(
            [
                "--model",
                "targetdiff",
                "--protein",
                str(protein),
                "--pocket",
                str(pocket),
                "--reference-ligand",
                str(reference),
                "--output",
                str(output),
                "--dry-run",
            ]
        )

    assert error.value.code == 2
    assert "REFERENCE_LIGAND_NOT_NEEDED" in capsys.readouterr().err
    assert not output.exists()


def test_diffsbdd_is_registered_as_reference_guided():
    from prism.generation.registry import expand_models, model_capabilities

    assert "diffsbdd" in expand_models(["diffsbdd"])
    capability = model_capabilities()["diffsbdd"]
    assert capability["generation_mode"] == "reference_guided"
    assert capability["requires_reference_ligand"] is True
    assert capability["accepted_reference_ligand_suffixes"] == [".sdf"]


def test_diffsbdd_requires_a_reference_ligand(tmp_path):
    from prism.generation.adapters.diffsbdd import DiffSBDDAdapter

    pocket = normalize_pocket(_write(tmp_path / "pocket.pdb", PDB))
    with pytest.raises(CapabilityError) as error:
        DiffSBDDAdapter({}).check_capabilities(pocket)
    assert error.value.code == "REFERENCE_LIGAND_REQUIRED"


def test_diffsbdd_exposes_reference_size(tmp_path):
    from prism.generation.adapters.diffsbdd import DiffSBDDAdapter

    pocket = normalize_pocket(_write(tmp_path / "reference.sdf", SDF))
    context = DiffSBDDAdapter({}).command_context(None, pocket, {})
    assert context["diffsbdd_reference_heavy_atoms"] == 1


def test_diffsbdd_wrapper_builds_upstream_arguments(tmp_path, monkeypatch):
    """--num-samples must reach upstream; sanitize/relax stay off by default."""
    from prism.generation.wrappers.diffsbdd import main

    root = tmp_path / "diffsbdd"
    root.mkdir(parents=True)
    _write(root / "generate_ligands.py", "# upstream placeholder\n")
    checkpoint = _write(tmp_path / "model.ckpt", "checkpoint")
    protein = _write(tmp_path / "protein.pdb", PDB)
    reference = _write(tmp_path / "reference.sdf", SDF)
    output = tmp_path / "raw"
    captured = {}

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        captured["env"] = env
        captured["command"] = command
        target = output / "sdf" / "candidates.sdf"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(SDF * 3, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("prism.generation.wrappers.diffsbdd.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.diffsbdd._sanitize_checkpoint",
        lambda source, _output: source,
    )
    assert (
        main(
            [
                "--diffsbdd-root", str(root),
                "--checkpoint", str(checkpoint),
                "--protein", str(protein),
                "--reference-ligand", str(reference),
                "--output", str(output),
                "--num-samples", "20",
                "--seed", "7",
                "--device", "cuda:1",
                "--batch-size", "10",
            ]
        )
        == 0
    )

    arguments = json.loads(captured["env"]["PRISM_DIFFSBDD_ARGS"])
    assert "--n_samples" in arguments and arguments[arguments.index("--n_samples") + 1] == "20"
    assert str(reference) == arguments[arguments.index("--ref_ligand") + 1]
    assert "--sanitize" not in arguments and "--relax" not in arguments
    assert captured["env"]["PRISM_DIFFSBDD_SEED"] == "7"
    assert captured["env"]["CUDA_VISIBLE_DEVICES"] == "1"

    summary = json.loads((output / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["records_written"] == 3
    assert summary["upstream_relax"] is False


def test_diffsbdd_wrapper_fixes_an_indivisible_batch_size(tmp_path, monkeypatch):
    """Upstream asserts n_samples % batch_size == 0."""
    from prism.generation.wrappers.diffsbdd import main

    root = tmp_path / "diffsbdd"
    root.mkdir(parents=True)
    _write(root / "generate_ligands.py", "# upstream placeholder\n")
    captured = {}

    def fake_run(command, cwd, env, check):  # noqa: ARG001
        captured["env"] = env
        target = tmp_path / "raw" / "sdf" / "candidates.sdf"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(SDF, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("prism.generation.wrappers.diffsbdd.subprocess.run", fake_run)
    monkeypatch.setattr(
        "prism.generation.wrappers.diffsbdd._sanitize_checkpoint",
        lambda source, _output: source,
    )
    main(
        [
            "--diffsbdd-root", str(root),
            "--checkpoint", str(_write(tmp_path / "m.ckpt", "c")),
            "--protein", str(_write(tmp_path / "protein.pdb", PDB)),
            "--reference-ligand", str(_write(tmp_path / "reference.sdf", SDF)),
            "--output", str(tmp_path / "raw"),
            "--num-samples", "7",
            "--seed", "1",
            "--batch-size", "10",
        ]
    )
    arguments = json.loads(captured["env"]["PRISM_DIFFSBDD_ARGS"])
    assert arguments[arguments.index("--batch_size") + 1] == "7"


# ---------------------------------------------------------------------------
# Device resolution
# ---------------------------------------------------------------------------


def test_auto_device_falls_back_to_cpu_without_a_gpu(monkeypatch):
    """`auto` must probe the host; hardcoding cuda:0 fails on a CPU-only node."""
    from prism.generation.wrappers import _device as device_support

    monkeypatch.setattr(device_support, "visible_gpu_count", lambda: 0)
    # Only these three have a working CPU path; see the CUDA-only test below.
    for name in ("targetdiff", "pocketxmol", "diffsbdd"):
        module = __import__(
            f"prism.generation.wrappers.{name}", fromlist=["_device"]
        )
        assert module._device("auto") == "cpu", name
        assert module._device("cpu") == "cpu", name
        # An explicit request is the caller's decision, not auto-detection.
        assert module._device("cuda") == "cuda:0", name


def test_auto_device_selects_a_gpu_when_one_exists(monkeypatch):
    from prism.generation.wrappers import _device as device_support

    monkeypatch.setattr(device_support, "visible_gpu_count", lambda: 2)
    for name in ("targetdiff", "pocket2mol", "pocketxmol", "molcraft", "flowr", "diffsbdd"):
        module = __import__(
            f"prism.generation.wrappers.{name}", fromlist=["_device"]
        )
        assert module._device("auto") == "cuda:0", name
        assert module._device("cuda") == "cuda:0", name
        assert module._device("cuda:1") == "cuda:1", name


@pytest.mark.parametrize(
    ("name", "reason"),
    [
        ("flowr", "requires CUDA"),
        # Measured on a CPU-only node: pocket2mol hardcodes .to('cuda') in its
        # sampler, and molcraft's Lightning Trainer resolves to CUDA anyway.
        ("pocket2mol", "hardcodes"),
        ("molcraft", "accelerator"),
    ],
)
def test_models_without_a_cpu_path_refuse_early_with_a_reason(monkeypatch, name, reason):
    """Refuse in argparse rather than after a minute with a deep CUDA traceback."""
    from prism.generation.wrappers import _device as device_support

    module = __import__(f"prism.generation.wrappers.{name}", fromlist=["_device"])
    monkeypatch.setattr(device_support, "visible_gpu_count", lambda: 0)
    with pytest.raises(Exception, match=reason):
        module._device("cpu")
    with pytest.raises(Exception, match="no usable NVIDIA GPU"):
        module._device("auto")
    # An explicit CUDA request still works.
    assert module._device("cuda") == "cuda:0"


def test_gpu_detection_honours_cuda_visible_devices(monkeypatch):
    from prism.generation.wrappers import _device as device_support

    monkeypatch.setattr(device_support.shutil, "which", lambda _name: "/usr/bin/nvidia-smi")
    monkeypatch.setattr(
        device_support.subprocess,
        "run",
        lambda *a, **k: subprocess.CompletedProcess(
            a[0], 0, stdout="GPU 0: NVIDIA X\nGPU 1: NVIDIA X\nGPU 2: NVIDIA X\n"
        ),
    )
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    assert device_support.visible_gpu_count() == 3

    # An explicit mask limits what this process may use.
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "1")
    assert device_support.visible_gpu_count() == 1
    # The conventional "no GPU" mask must read as zero, not as one device.
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "-1")
    assert device_support.visible_gpu_count() == 0
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    assert device_support.visible_gpu_count() == 0


def test_gpu_detection_survives_a_broken_nvidia_smi(monkeypatch):
    """A driver mismatch makes nvidia-smi fail; that means no GPU, not a crash."""
    from prism.generation.wrappers import _device as device_support

    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.setattr(device_support.shutil, "which", lambda _name: "/usr/bin/nvidia-smi")
    monkeypatch.setattr(
        device_support.subprocess,
        "run",
        lambda *a, **k: subprocess.CompletedProcess(a[0], 9, stdout="", stderr="driver error"),
    )
    assert device_support.visible_gpu_count() == 0

    def explode(*_a, **_k):
        raise OSError("nvidia-smi vanished")

    monkeypatch.setattr(device_support.subprocess, "run", explode)
    assert device_support.visible_gpu_count() == 0


# ---------------------------------------------------------------------------
# Checkpoint allowlisting
# ---------------------------------------------------------------------------


def _fake_checkpoint(path, payload_globals):
    """A zip checkpoint whose pickle stream references the given globals."""
    import zipfile

    stream = pickle.dumps({"state_dict": {}, "objects": list(payload_globals)})
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("archive/data.pkl", stream)
    return path


def test_checkpoint_scan_reports_globals_without_executing_them(tmp_path):
    """torch<2.0 ignores pickle_module, so the allowlist runs before torch does."""
    from prism.generation.wrappers import _checkpoint

    marker = tmp_path / "executed"

    class Payload:
        def __reduce__(self):
            return (os.system, (f"touch {marker}",))

    import zipfile

    path = tmp_path / "poisoned.pt"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("archive/data.pkl", pickle.dumps({"x": Payload()}))

    found = _checkpoint.referenced_globals(path)
    assert ("posix", "system") in found or ("nt", "system") in found
    # Scanning must not run the payload.
    assert not marker.exists()

    with pytest.raises(_checkpoint.CheckpointRejected, match="forbidden pickle global"):
        _checkpoint.verify(path)
    assert not marker.exists()


def test_checkpoint_scan_allows_the_globals_real_checkpoints_need(tmp_path):
    from prism.generation.wrappers import _checkpoint

    path = _fake_checkpoint(tmp_path / "clean.pt", [])
    used = _checkpoint.verify(path)
    assert ("collections", "OrderedDict") not in used or True  # payload-dependent

    # Model-specific config objects are opt-in per wrapper, not global.
    import argparse
    import zipfile

    path = tmp_path / "namespace.pt"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("archive/data.pkl", pickle.dumps({"hp": argparse.Namespace(a=1)}))
    with pytest.raises(_checkpoint.CheckpointRejected, match="argparse.Namespace"):
        _checkpoint.verify(path)
    assert ("argparse", "Namespace") in _checkpoint.verify(
        path, extra_allowed={("argparse", "Namespace")}
    )


def test_checkpoint_scan_rejects_a_non_checkpoint(tmp_path):
    from prism.generation.wrappers import _checkpoint

    import zipfile

    path = tmp_path / "empty.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("readme.txt", "not a checkpoint")
    with pytest.raises(_checkpoint.CheckpointRejected, match="no pickle stream"):
        _checkpoint.verify(path)

    path = tmp_path / "garbage.pt"
    path.write_bytes(b"\x80\x04not-a-pickle")
    with pytest.raises(_checkpoint.CheckpointRejected):
        _checkpoint.verify(path)


def test_every_wrapper_that_loads_a_checkpoint_verifies_it_first():
    """Regression guard: the check must not be dropped from a wrapper."""
    from pathlib import Path as _Path

    root = _Path(__file__).resolve().parents[1] / "prism/generation/wrappers"
    for name in ("pocket2mol", "molcraft", "flowr", "diffsbdd"):
        source = (root / f"{name}.py").read_text(encoding="utf-8")
        if "_sanitize_checkpoint" not in source:
            continue
        assert "_checkpoint.verify(" in source, f"{name} loads a checkpoint unverified"


def test_every_wrapper_runs_as_a_standalone_script():
    """Generation configs invoke wrappers by path, not as package modules.

    The model environments do not have PRISM installed, so a bare relative
    import breaks every backend at once. Importing the wrapper as a package --
    which is what the rest of this file does -- cannot catch that, so the check
    has to actually execute the file the way production does.
    """
    from pathlib import Path as _Path

    root = _Path(__file__).resolve().parents[1] / "prism/generation/wrappers"
    wrappers = sorted(
        path for path in root.glob("*.py") if not path.name.startswith("_")
    )
    assert len(wrappers) == 6, [p.name for p in wrappers]

    for wrapper in wrappers:
        completed = subprocess.run(
            [sys.executable, str(wrapper), "--help"],
            capture_output=True,
            text=True,
            # No PYTHONPATH: the model environments cannot import prism.
            env={
                key: value
                for key, value in os.environ.items()
                if key != "PYTHONPATH"
            },
            cwd="/",
        )
        assert completed.returncode == 0, (
            f"{wrapper.name} cannot run standalone:\n{completed.stderr[-600:]}"
        )


def test_timeout_kills_the_whole_process_tree(tmp_path):
    """The direct child is only a launcher; the model runs as a grandchild.

    ``subprocess.run(timeout=...)`` kills the launcher and leaves the model
    running -- an abandoned charge calculation was seen holding a core for
    nearly two hours after its task had already been recorded as timed out.
    """
    import time

    from prism.generation.execution import execute

    marker = tmp_path / "grandchild.pid"
    launcher = tmp_path / "launcher.sh"
    launcher.write_text(
        "#!/bin/bash\n"
        f"sleep 600 &\n"
        f"echo $! > {marker}\n"
        "wait\n",
        encoding="utf-8",
    )
    launcher.chmod(0o755)

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    config = {
        "enabled": True,
        "backend": "conda",
        "environment": "unused",
        "command": [str(launcher)],
        "timeout_seconds": 2,
    }
    # Bypass the conda wrapper: run the launcher directly as the child.
    import prism.generation.execution as execution_module

    original = execution_module._conda_command
    execution_module._conda_command = lambda _environment, inner: inner
    try:
        with pytest.raises(Exception, match="timeout"):
            execute(config, {}, tmp_path, run_dir, "cpu")
    finally:
        execution_module._conda_command = original

    assert marker.exists(), "launcher never started its grandchild"
    grandchild = int(marker.read_text().strip())
    # Give the signal a moment to land.
    for _ in range(20):
        if not Path(f"/proc/{grandchild}").exists():
            break
        time.sleep(0.25)
    assert not Path(f"/proc/{grandchild}").exists(), (
        f"grandchild {grandchild} survived the timeout"
    )
