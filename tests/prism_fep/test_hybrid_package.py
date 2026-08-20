from prism.fep.modeling.hybrid_package import HybridPackageBuilder
from prism.fep.modeling.core import FEPScaffoldBuilder


def test_dedupe_moleculetype_sections_removes_duplicate_block():
    builder = HybridPackageBuilder()
    content = """; header

[ moleculetype ]
; Name            nrexcl
HYB  3

[ moleculetype ]
; Name            nrexcl
HYB                  3

[ atoms ]
1  C  1 LIG C1 1 0.0 12.0
"""

    cleaned = builder._dedupe_moleculetype_sections(content)
    assert cleaned.count("[ moleculetype ]") == 1
    assert cleaned.count("HYB") == 1
    assert "[ atoms ]" in cleaned


def test_parse_defaults_block_stops_before_includes():
    builder = HybridPackageBuilder()
    top = """[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333

#include "atomtypes_LIG.itp"
#include "LIG.itp"
"""

    block = builder._parse_defaults_block(top)
    assert "#include" not in block
    assert "[ defaults ]" in block
    assert "0.8333" in block


def test_collect_cgenff_hybrid_bonded_itp_uses_registry_subdirs(tmp_path):
    builder = HybridPackageBuilder()
    ref_dir = tmp_path / "ref"
    mut_dir = tmp_path / "mut"
    (ref_dir / "LIG.cgenff2gmx").mkdir(parents=True)
    (mut_dir / "LIG.amb2gmx").mkdir(parents=True)
    (ref_dir / "LIG.cgenff2gmx" / "cgenff_bonded_LIG.itp").write_text("[ bondtypes ]\nA B 1 0.1 100\n")
    (mut_dir / "LIG.amb2gmx" / "cgenff_bonded_LIG.itp").write_text("[ angletypes ]\nA B C 1 100 10\n")

    content = builder._build_cgenff_hybrid_bonded_itp(ref_dir, mut_dir)
    assert "[ bondtypes ]" in content
    assert "[ angletypes ]" in content


def test_discover_ligand_dir_uses_registry_known_layouts(tmp_path):
    scaffold = FEPScaffoldBuilder(output_dir=str(tmp_path / "out"))
    lig_dir = tmp_path / "lig"
    (lig_dir / "LIG.openff2gmx").mkdir(parents=True)
    (lig_dir / "LIG.openff2gmx" / "LIG.itp").write_text("; itp")

    assert scaffold._discover_ligand_dir(lig_dir) == lig_dir / "LIG.openff2gmx"
