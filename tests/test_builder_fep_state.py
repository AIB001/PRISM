"""Regression tests for temporary builder state restoration in FEP helpers."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from prism.builder.workflow_fep import FEPWorkflowMixin
from prism.fep.modeling import e2e
from prism.fep.visualize.reporting import MappingReportService


class _DummySystemBuilder:
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.model_dir = self.output_dir / "GMX_PROLIG_MD"


class _DummyMDPGenerator:
    def __init__(self, output_dir):
        self.output_dir = output_dir


class _DummyBuilder(FEPWorkflowMixin):
    def __init__(self, tmp_path):
        self.fep_mode = True
        self.output_dir = str(tmp_path / "original_out")
        self.protein_path = str(tmp_path / "protein.pdb")
        self.ligand_paths = [str(tmp_path / "ligand.mol2")]
        self.lig_ff_dirs = [str(tmp_path / "LIG.amb2gmx")]
        self.ligand_forcefield = "gaff2"
        self.forcefield_paths = ["ff_ref", "ff_mut"]
        self.system_builder = _DummySystemBuilder(self.output_dir)
        self.mdp_generator = _DummyMDPGenerator(self.output_dir)

    def run_normal(self):
        raise RuntimeError("normal-build-failed")

    def _build_ligand_only_system(self, _output_dir, box_size=None):
        _ = box_size
        raise RuntimeError("ligand-only-build-failed")

    def generate_ligand_forcefield(self):
        raise RuntimeError("ff-generation-failed")


def test_build_standard_system_restores_state_after_run_normal_failure(tmp_path):
    builder = _DummyBuilder(tmp_path)
    original_output = builder.output_dir
    original_protein = builder.protein_path
    original_ligands = list(builder.ligand_paths)
    original_ff_dirs = list(builder.lig_ff_dirs)
    original_model_dir = builder.system_builder.model_dir
    original_mdp_output = builder.mdp_generator.output_dir

    with pytest.raises(RuntimeError, match="normal-build-failed"):
        builder._build_standard_system(str(tmp_path / "fep_bound"), use_protein=True)

    assert builder.fep_mode is True
    assert builder.output_dir == original_output
    assert builder.protein_path == original_protein
    assert builder.ligand_paths == original_ligands
    assert builder.lig_ff_dirs == original_ff_dirs
    assert builder.system_builder.output_dir == Path(original_output)
    assert builder.system_builder.model_dir == original_model_dir
    assert builder.mdp_generator.output_dir == original_mdp_output


def test_build_standard_system_restores_state_after_ligand_only_failure(tmp_path):
    builder = _DummyBuilder(tmp_path)
    original_output = builder.output_dir
    original_protein = builder.protein_path
    original_ligands = list(builder.ligand_paths)
    original_ff_dirs = list(builder.lig_ff_dirs)

    with pytest.raises(RuntimeError, match="ligand-only-build-failed"):
        builder._build_standard_system(str(tmp_path / "fep_unbound"), use_protein=False)

    assert builder.fep_mode is True
    assert builder.output_dir == original_output
    assert builder.protein_path == original_protein
    assert builder.ligand_paths == original_ligands
    assert builder.lig_ff_dirs == original_ff_dirs
    assert builder.system_builder.output_dir == Path(original_output)
    assert builder.system_builder.model_dir == Path(original_output) / "GMX_PROLIG_MD"
    assert builder.mdp_generator.output_dir == original_output


def test_generate_mutant_ligand_ff_restores_state_after_failure(tmp_path):
    builder = _DummyBuilder(tmp_path)
    original_output = builder.output_dir
    original_ligands = list(builder.ligand_paths)
    original_ff_dirs = list(builder.lig_ff_dirs)

    with pytest.raises(RuntimeError, match="ff-generation-failed"):
        builder._generate_mutant_ligand_ff(str(tmp_path / "mutant.mol2"), str(tmp_path / "mutant_out"))

    assert builder.output_dir == original_output
    assert builder.ligand_paths == original_ligands
    assert builder.lig_ff_dirs == original_ff_dirs


def test_cleanup_fep_build_artifacts_keeps_initial_ligand_ff_dirs(tmp_path):
    builder = _DummyBuilder(tmp_path)
    fep_root = tmp_path / "GMX_PROLIG_FEP"
    build_dir = fep_root / "_build"
    ref_ff = build_dir / "bound_md" / "LIG.openff2gmx"
    mut_ff = build_dir / "mutant_ligand_ff" / "LIG.openff2gmx"
    bound_md = build_dir / "bound_md" / "GMX_PROLIG_MD"
    unbound_md = build_dir / "unbound_md" / "GMX_PROLIG_MD"
    weird_nested = fep_root / "GMX_PROLIG_FEP"

    for directory in (ref_ff, mut_ff, bound_md, unbound_md, weird_nested):
        directory.mkdir(parents=True, exist_ok=True)

    builder._cleanup_fep_build_artifacts(str(fep_root))

    assert ref_ff.exists()
    assert mut_ff.exists()
    assert not bound_md.exists()
    assert not unbound_md.exists()
    assert weird_nested.exists()


def test_resolve_generated_ligand_ff_dir_prefers_registry_primary_dir(tmp_path):
    builder = _DummyBuilder(tmp_path)
    output_dir = tmp_path / "case"
    (output_dir / "LIG.amb2gmx").mkdir(parents=True)
    (output_dir / "LIG.opls2gmx").mkdir()

    assert builder._resolve_generated_ligand_ff_dir(str(output_dir)).endswith("LIG.amb2gmx")


def test_resolve_generated_ligand_ff_dir_prefers_suffix_one_for_multi_ligand(tmp_path):
    builder = _DummyBuilder(tmp_path)
    builder.ligand_forcefield = "mmff"
    output_dir = tmp_path / "case"
    lig_root = output_dir / "Ligand_Forcefield"
    (lig_root / "LIG.sp2gmx_2").mkdir(parents=True)
    (lig_root / "LIG.sp2gmx_1").mkdir()

    assert builder._resolve_generated_ligand_ff_dir(str(output_dir)).endswith("LIG.sp2gmx_1")


def test_resolve_ligand_ff_artifact_uses_nested_registry_dirs(tmp_path):
    builder = _DummyBuilder(tmp_path)
    ff_dir = tmp_path / "ff"
    nested = ff_dir / "LIG.openff2gmx"
    nested.mkdir(parents=True)
    (nested / "LIG.gro").write_text("gro")

    resolved = builder._resolve_ligand_ff_artifact(str(ff_dir), "LIG.gro")
    assert resolved.endswith("LIG.openff2gmx/LIG.gro")


def test_mapping_report_service_preserves_cgenff_origin_labels(tmp_path):
    service = MappingReportService()
    gui_dir = tmp_path / "gui"
    web_dir = tmp_path / "web"
    (gui_dir / "gromacs").mkdir(parents=True)
    (gui_dir / "gromacs" / "LIG.itp").write_text("; gui")
    (web_dir / "charmm36.ff").mkdir(parents=True)
    (web_dir / "charmm36.ff" / "charmm36.itp").write_text("; ff")

    assert service._describe_mapping_forcefield("cgenff", [gui_dir], "cgenff") == "CGENFF (CHARMM-GUI)"
    assert service._describe_mapping_forcefield("cgenff", [], "cgenff") == "CGENFF (WEBSITE)"


def test_e2e_build_uses_hybrid_service_instead_of_template_fallback(monkeypatch, tmp_path):
    calls = {}

    class DummyHybridService:
        def __init__(self, **kwargs):
            calls["kwargs"] = kwargs

        def build_from_forcefield_dirs(self, **kwargs):
            calls["build"] = kwargs
            out = Path(kwargs["hybrid_output_dir"])
            out.mkdir(parents=True, exist_ok=True)
            hybrid_itp = out / "hybrid.itp"
            hybrid_itp.write_text("; hybrid")
            return SimpleNamespace(
                hybrid_itp=hybrid_itp,
                mapping=SimpleNamespace(common=[1, 2], transformed_a=[3], transformed_b=[4]),
            )

    class DummyLayout:
        root = Path("root")
        bound_dir = Path("bound")
        unbound_dir = Path("unbound")

    class DummyBuilder:
        @staticmethod
        def _normalize_output_dir(output_dir):
            return Path(output_dir)

        def __init__(self, **kwargs):
            calls["builder_init"] = kwargs

        def build_from_components(self, **kwargs):
            calls["scaffold"] = kwargs
            return DummyLayout()

    monkeypatch.setattr(e2e, "HybridBuildService", DummyHybridService)
    monkeypatch.setattr(e2e, "FEPScaffoldBuilder", DummyBuilder)

    result = e2e.build_fep_system_from_prism_ligands(
        receptor_pdb="receptor.pdb",
        reference_ligand_dir="ref_ff",
        mutant_ligand_dir="mut_ff",
        output_dir=str(tmp_path / "fep_out"),
    )

    assert calls["build"]["reference_ligand_dir"] == "ref_ff"
    assert calls["build"]["mutant_ligand_dir"] == "mut_ff"
    assert calls["scaffold"]["hybrid_itp"].endswith("temp_hybrid/hybrid.itp")
    assert result == DummyLayout.root
