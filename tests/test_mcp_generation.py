"""MCP generation tools.

These wrap the same orchestration and quality-control code the CLI uses, so
the tests here cover the MCP layer's own contract: tools are registered,
results are JSON, failures are reported as data rather than raised, and the
truncated problem listing shows the most severe records.
"""

import json

import pytest

pytest.importorskip("mcp")

from prism.mcp import generation as generation_tools


class _Collector:
    """Stand-in for FastMCP that captures the registered tool callables."""

    def __init__(self):
        self.tools = {}

    def tool(self, *_args, **_kwargs):
        def decorator(function):
            self.tools[function.__name__] = function
            return function

        return decorator


@pytest.fixture(scope="module")
def tools():
    collector = _Collector()
    generation_tools.register(collector)
    return collector.tools


# Methanol, with and without the explicit hydrogens an MD build requires.
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

HEAVY_ONLY = """methanol
  PRISM    3D

  2  1  0  0  0  0  0  0  0  0999 V2000
   -0.3480    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0150    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
"""

CLEAN = """
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2000    1.3000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2000   -1.2000    0.1000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
M  END
$$$$
"""

# Oxygen carrying one single and one double bond: explicit valence 3.
BROKEN = """valence-error
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


def test_every_generation_tool_is_registered(tools):
    assert set(tools) == {
        "list_generation_models",
        "generate_ligands",
        "check_generated_ligands",
        "prepare_generated_for_md",
    }


def test_list_generation_models_reports_conditioning_requirements(tools):
    payload = json.loads(tools["list_generation_models"]())

    models = payload["models"]
    assert set(models) >= {
        "pocket2mol", "targetdiff", "pocketxmol", "molcraft", "flowr", "diffsbdd",
    }
    # An agent has to be able to tell which models need a conditioning ligand.
    assert models["diffsbdd"]["generation_mode"] == "reference_guided"
    assert models["diffsbdd"]["requires_reference_ligand"] is True
    assert models["targetdiff"]["generation_mode"] == "direct"
    assert models["targetdiff"]["requires_reference_ligand"] is False


def test_check_generated_ligands_summarises_a_real_sdf(tmp_path):
    pytest.importorskip("rdkit")
    collector = _Collector()
    generation_tools.register(collector)
    sdf = tmp_path / "candidates.sdf"
    sdf.write_text(CLEAN * 2 + BROKEN, encoding="utf-8")

    payload = json.loads(
        collector.tools["check_generated_ligands"](sdf_path=str(sdf), qc_level="standard")
    )

    assert payload["records_checked"] == 3
    assert payload["status_counts"].get("FAIL_CHEMISTRY") == 1
    assert payload["flag_counts"].get("sanitize_failed") == 1
    # The failing record carries its proposal, and PRISM has not applied it.
    failing = [r for r in payload["problem_records"] if r["status"] == "FAIL_CHEMISTRY"]
    assert len(failing) == 1
    assert failing[0]["repair"]["confidence"] == "unique_minimal"
    assert failing[0]["repair"]["applied"] is False


def test_problem_records_put_failures_before_warnings(tmp_path):
    """The listing is truncated; a failure pushed past the cut reads as clean."""
    pytest.importorskip("rdkit")
    collector = _Collector()
    generation_tools.register(collector)
    sdf = tmp_path / "many.sdf"
    # 40 warning-only records, then the single chemistry failure last.
    warning = CLEAN.rstrip().replace("$$$$", ">  <SCORE>\n1\n$$$$") + "\n"
    sdf.write_text(warning * 40 + BROKEN, encoding="utf-8")

    payload = json.loads(
        collector.tools["check_generated_ligands"](sdf_path=str(sdf), qc_level="standard")
    )

    assert payload["records_checked"] == 41
    assert payload["problem_record_count"] > payload["problem_records_shown"]
    assert payload["problem_records"][0]["status"] == "FAIL_CHEMISTRY"


def _generation_run(tmp_path, specs):
    """Minimal finished-run layout: manifest plus per-candidate files."""
    run_dir = tmp_path / "generated_ligand"
    candidates = []
    for index, (status, body) in enumerate(specs, start=1):
        molecules = run_dir / "models" / "flowr" / "molecules"
        hydrogens = run_dir / "models" / "flowr" / "molecules_h"
        molecules.mkdir(parents=True, exist_ok=True)
        hydrogens.mkdir(parents=True, exist_ok=True)
        heavy = molecules / f"candidate_{index:06d}.sdf"
        heavy.write_text(HEAVY_ONLY, encoding="utf-8")
        hydrogenated = hydrogens / f"candidate_{index:06d}_h.sdf"
        hydrogenated.write_text(body, encoding="utf-8")
        candidates.append(
            {
                "candidate_id": f"flowr_{index:06d}",
                "model": "flowr",
                "path": str(heavy),
                "source_path": str(heavy),
                "valid_3d": True,
                "audit_status": status,
                "quality_flags": [],
                "canonical_smiles": "CO",
                "hydrogenated_path": str(hydrogenated),
            }
        )
    (run_dir / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "status": "SUCCEEDED",
                "request": {"protein": "/data/receptor.pdb"},
                "models": {"flowr": {"status": "SUCCEEDED", "candidates": candidates}},
            }
        ),
        encoding="utf-8",
    )
    return run_dir


def test_prepare_generated_for_md_exports_only_hydrogenated_ligands(tmp_path):
    """The tool's whole purpose: what it returns must be buildable. Raw model
    output builds into a hydrogen-free topology with no error at any stage."""
    collector = _Collector()
    generation_tools.register(collector)
    run_dir = _generation_run(
        tmp_path,
        [
            ("PASS", WITH_HYDROGENS),
            ("PASS", HEAVY_ONLY),
            ("FAIL_CHEMISTRY", WITH_HYDROGENS),
        ],
    )

    payload = json.loads(
        collector.tools["prepare_generated_for_md"](
            run_dir=str(run_dir), output_dir=str(tmp_path / "md_inputs")
        )
    )

    assert payload["exported_count"] == 1
    assert payload["skipped_reasons"] == {
        "NO_EXPLICIT_HYDROGENS": 1,
        "FAILED_QC": 1,
    }
    assert payload["ligands"][0]["hydrogen_count"] == 4
    assert "prism -pf" in payload["build_command"]


def test_tools_report_failures_as_data_not_exceptions(tools, tmp_path):
    """MCP tools must return an error payload; a raise would break the session."""
    payload = json.loads(tools["check_generated_ligands"](sdf_path=str(tmp_path / "missing.sdf")))
    assert "error" in payload

    payload = json.loads(
        tools["generate_ligands"](
            protein_path=str(tmp_path / "missing.pdb"),
            pocket_path=str(tmp_path / "missing.sdf"),
            generation_config=str(tmp_path / "missing.yaml"),
        )
    )
    assert "error" in payload


def test_generate_ligands_rejects_an_unknown_model(tools, tmp_path):
    payload = json.loads(
        tools["generate_ligands"](
            protein_path=str(tmp_path / "p.pdb"),
            pocket_path=str(tmp_path / "l.sdf"),
            generation_config=str(tmp_path / "c.yaml"),
            models="not-a-model",
        )
    )
    assert "error" in payload
    assert "not-a-model" in payload["error"]


def test_server_registers_generation_alongside_the_existing_tools():
    import asyncio

    from prism.mcp import mcp

    names = {tool.name for tool in asyncio.run(mcp.list_tools())}
    assert {
        "list_generation_models",
        "generate_ligands",
        "check_generated_ligands",
        "prepare_generated_for_md",
    } <= names
    # The generation tools are additions, not replacements.
    assert {"build_system", "analyze_rmsd"} <= names
