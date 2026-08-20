from prism.forcefield.registry import (
    discover_ligand_directory,
    get_forcefield_output_subdirs,
    get_primary_forcefield_output_subdir,
    get_forcefield_spec,
    iter_existing_ligand_output_dirs,
    iter_known_ligand_output_subdirs,
    normalize_forcefield_name,
    resolve_ligand_artifact,
    resolve_ligand_artifact_dir,
)
from prism.forcefield import get_generator_info
from prism.builder.cli import _default_fep_output_dir
from prism.sim.utils import find_forcefield_dir


def test_registry_normalizes_aliases():
    assert normalize_forcefield_name("charmm_gui") == "charmm-gui"
    assert normalize_forcefield_name("opls-aa") == "opls"
    assert normalize_forcefield_name("both") == "hybrid"


def test_registry_exposes_legacy_and_canonical_output_dirs():
    assert get_forcefield_output_subdirs("gaff2") == ("LIG.amb2gmx",)
    assert get_forcefield_output_subdirs("mmff") == ("LIG.sp2gmx", "LIG.mmff2gmx")
    assert "LIG.match2gmx" in get_forcefield_output_subdirs("match")
    assert get_primary_forcefield_output_subdir("mmff") == "LIG.sp2gmx"


def test_registry_marks_generic_atom_name_backends():
    assert get_forcefield_spec("openff").generic_atom_names is True
    assert get_forcefield_spec("gaff2").generic_atom_names is False


def test_iter_known_ligand_output_subdirs_is_deduplicated():
    subdirs = iter_known_ligand_output_subdirs()
    assert len(subdirs) == len(set(subdirs))
    assert "LIG.amb2gmx" in subdirs
    assert "LIG.openff2gmx" in subdirs
    assert "LIG.sp2gmx" in subdirs


def test_resolve_ligand_artifact_prefers_direct_file(tmp_path):
    direct = tmp_path / "LIG.itp"
    direct.write_text("; direct")
    assert resolve_ligand_artifact(tmp_path, "LIG.itp") == direct


def test_resolve_ligand_artifact_checks_forcefield_subdirs(tmp_path):
    ff_dir = tmp_path / "LIG.openff2gmx"
    ff_dir.mkdir()
    itp = ff_dir / "LIG.itp"
    itp.write_text("; openff")

    assert resolve_ligand_artifact(tmp_path, "LIG.itp", forcefield="openff") == itp


def test_resolve_ligand_artifact_accepts_legacy_swissparam_layout(tmp_path):
    ff_dir = tmp_path / "LIG.sp2gmx"
    ff_dir.mkdir()
    gro = ff_dir / "LIG.gro"
    gro.write_text("gro")

    assert resolve_ligand_artifact(tmp_path, "LIG.gro", forcefield="mmff") == gro


def test_resolve_ligand_artifact_dir_returns_parent_directory(tmp_path):
    ff_dir = tmp_path / "LIG.charmm2gmx"
    ff_dir.mkdir()
    (ff_dir / "LIG.itp").write_text("; itp")

    assert resolve_ligand_artifact_dir(tmp_path, "LIG.itp") == ff_dir


def test_iter_existing_ligand_output_dirs_respects_registry_order(tmp_path):
    (tmp_path / "LIG.openff2gmx").mkdir()
    (tmp_path / "LIG.amb2gmx").mkdir()

    dirs = iter_existing_ligand_output_dirs(tmp_path)
    assert [path.name for path in dirs][:2] == ["LIG.amb2gmx", "LIG.openff2gmx"]


def test_discover_ligand_directory_prefers_registry_layout(tmp_path):
    ff_dir = tmp_path / "LIG.opls2gmx"
    ff_dir.mkdir()
    (ff_dir / "LIG.itp").write_text("; itp")

    assert discover_ligand_directory(tmp_path) == ff_dir


def test_cli_default_fep_output_dir_uses_canonical_ligand_name():
    assert _default_fep_output_dir("amber14sb_OL15", "charmm_gui") == "amber14sb_ol15-mut_charmm_gui"
    assert _default_fep_output_dir("amber14sb_OL15", "opls-aa") == "amber14sb_ol15-mut_opls"


def test_get_generator_info_uses_registry_output_dirs():
    info = get_generator_info()
    assert info["GAFF2"]["output_dir"] == "LIG.amb2gmx"
    assert info["MMFF"]["output_dir"] == "LIG.mmff2gmx"
    assert info["OpenFF"]["output_dir"] == "LIG.openff2gmx"


def test_find_forcefield_dir_accepts_opls_and_swissparam_layouts(tmp_path):
    opls_dir = tmp_path / "LIG.opls2gmx"
    opls_dir.mkdir()
    assert find_forcefield_dir(tmp_path) == str(opls_dir)
