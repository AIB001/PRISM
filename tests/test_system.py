"""Tests for prism.utils.system — SystemBuilder pure-logic methods.

Only tests methods that do NOT require GROMACS (file parsing, text processing).
"""

import pytest
from prism.utils.system import SystemBuilder


# ---------------------------------------------------------------------------
# Helper: build a minimal PDB ATOM line
# ---------------------------------------------------------------------------
def _atom_line(serial, name, resname, chain, resseq, x=0.0, y=0.0, z=0.0):
    """Build a correctly formatted PDB ATOM record."""
    return (
        f"ATOM  {serial:5d} {name:<4s} {resname:<3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]:>2s}  \n"
    )


def _hetatm_line(serial, name, resname, chain, resseq, x=0.0, y=0.0, z=0.0):
    """Build a correctly formatted PDB HETATM record."""
    return (
        f"HETATM{serial:5d} {name:<4s} {resname:<3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]:>2s}  \n"
    )


# ---------------------------------------------------------------------------
# Fixture: minimal SystemBuilder (no GROMACS needed)
# ---------------------------------------------------------------------------
@pytest.fixture
def builder(tmp_path):
    """SystemBuilder with minimal config, model_dir = tmp_path/GMX_PROLIG_MD."""
    config = {
        "general": {"gmx_command": "gmx"},
        "simulation": {"pH": 7.0},
        "box": {"distance": 1.5, "shape": "cubic", "center": True},
        "ions": {"neutral": True, "concentration": 0.15, "positive_ion": "NA", "negative_ion": "CL"},
    }
    sb = SystemBuilder(config, str(tmp_path), overwrite=True)
    return sb


# ===========================================================================
# TestCountHistidines
# ===========================================================================
class TestCountHistidines:
    def test_counts_his_residues(self, builder, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "N", "HIS", "A", 10)
            + _atom_line(2, "CA", "HIS", "A", 10)  # same residue
            + _atom_line(3, "N", "HIS", "A", 20)
            + "END\n"
        )
        assert builder._count_histidines(str(pdb)) == 2

    def test_ignores_hie_hid_hip(self, builder, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "N", "HIE", "A", 10)
            + _atom_line(2, "N", "HID", "A", 20)
            + _atom_line(3, "N", "HIP", "A", 30)
            + "END\n"
        )
        assert builder._count_histidines(str(pdb)) == 0

    def test_unique_residues_only(self, builder, tmp_path):
        """Multiple atoms in same HIS residue count as 1."""
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "N", "HIS", "A", 10)
            + _atom_line(2, "CA", "HIS", "A", 10)
            + _atom_line(3, "C", "HIS", "A", 10)
            + _atom_line(4, "O", "HIS", "A", 10)
            + "END\n"
        )
        assert builder._count_histidines(str(pdb)) == 1

    def test_no_histidines(self, builder, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(_atom_line(1, "N", "ALA", "A", 1) + _atom_line(2, "N", "GLY", "A", 2) + "END\n")
        assert builder._count_histidines(str(pdb)) == 0

    def test_ignores_hetatm(self, builder, tmp_path):
        """_count_histidines only checks ATOM records, not HETATM."""
        pdb = tmp_path / "test.pdb"
        pdb.write_text(_hetatm_line(1, "N", "HIS", "A", 10) + "END\n")
        assert builder._count_histidines(str(pdb)) == 0

    def test_multi_chain(self, builder, tmp_path):
        """HIS in different chains count separately."""
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "N", "HIS", "A", 10)
            + _atom_line(2, "N", "HIS", "B", 10)  # same resnum, diff chain
            + "END\n"
        )
        assert builder._count_histidines(str(pdb)) == 2


# ===========================================================================
# TestRenameHisForCmap
# ===========================================================================
class TestRenameHisForCmap:
    def test_renames_when_cmap_exists(self, builder, tmp_path):
        """HIS → HIE when force field has cmap.itp."""
        ff_dir = tmp_path / "amber19sb.ff"
        ff_dir.mkdir()
        (ff_dir / "cmap.itp").write_text("; cmap data\n")

        pdb = tmp_path / "test.pdb"
        pdb.write_text(_atom_line(1, "N", "HIS", "A", 10) + _atom_line(2, "CA", "HIS", "A", 10) + "END\n")

        ff_info = {"path": str(ff_dir)}
        result = builder._rename_his_for_cmap(str(pdb), ff_info)

        assert result is True
        content = pdb.read_text()
        for line in content.splitlines():
            if line.startswith("ATOM") and line[22:26].strip() == "10":
                assert line[17:20] == "HIE"

    def test_skips_when_no_cmap(self, builder, tmp_path):
        """No renaming if force field has no cmap.itp (e.g. amber14sb)."""
        ff_dir = tmp_path / "amber14sb.ff"
        ff_dir.mkdir()
        # No cmap.itp

        pdb = tmp_path / "test.pdb"
        original = _atom_line(1, "N", "HIS", "A", 10) + "END\n"
        pdb.write_text(original)

        ff_info = {"path": str(ff_dir)}
        result = builder._rename_his_for_cmap(str(pdb), ff_info)

        assert result is False
        assert pdb.read_text() == original  # unchanged

    def test_skips_when_no_ff_info(self, builder, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(_atom_line(1, "N", "HIS", "A", 10) + "END\n")

        assert builder._rename_his_for_cmap(str(pdb), None) is False
        assert builder._rename_his_for_cmap(str(pdb), {}) is False

    def test_skips_already_renamed(self, builder, tmp_path):
        """HIE/HID/HIP residues are not matched (line[17:20] != 'HIS')."""
        ff_dir = tmp_path / "amber19sb.ff"
        ff_dir.mkdir()
        (ff_dir / "cmap.itp").write_text("; cmap\n")

        pdb = tmp_path / "test.pdb"
        pdb.write_text(_atom_line(1, "N", "HIE", "A", 10) + _atom_line(2, "N", "HIP", "A", 20) + "END\n")

        ff_info = {"path": str(ff_dir)}
        result = builder._rename_his_for_cmap(str(pdb), ff_info)

        assert result is False  # nothing to rename

    def test_only_renames_remaining_his(self, builder, tmp_path):
        """When PROPKA already renamed some, _rename_his_for_cmap only touches leftover HIS."""
        ff_dir = tmp_path / "amber19sb.ff"
        ff_dir.mkdir()
        (ff_dir / "cmap.itp").write_text("; cmap\n")

        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "N", "HIP", "A", 10)  # PROPKA-renamed, skip
            + _atom_line(2, "N", "HIS", "A", 20)  # remaining, should rename
            + "END\n"
        )

        ff_info = {"path": str(ff_dir)}
        result = builder._rename_his_for_cmap(str(pdb), ff_info)

        assert result is True
        content = pdb.read_text()
        lines = [l for l in content.splitlines() if l.startswith("ATOM")]
        assert lines[0][17:20] == "HIP"  # untouched
        assert lines[1][17:20] == "HIE"  # renamed


# ===========================================================================
# TestFixImproperFromGromppErrors
# ===========================================================================
class TestFixImproperFromGromppErrors:
    def test_parses_single_error(self, builder):
        """Single error line → 1 line commented out."""
        topol = builder.model_dir / "topol.top"
        topol.write_text(
            "; header\n"  # line 1
            "[ dihedrals ]\n"  # line 2
            "   1    2    3    4     4\n"  # line 3 — this should be removed
            "[ bonds ]\n"  # line 4
        )

        stderr = "ERROR 1 [file topol.top, line 3]:\n" "  No default Per. Imp. Dih. types\n"
        fixed = builder._fix_improper_from_grompp_errors(stderr)

        assert fixed == 1
        content = topol.read_text()
        assert "; REMOVED by PRISM" in content
        # Line 3 should be commented, line 2 and 4 untouched
        lines = content.splitlines()
        assert lines[1] == "[ dihedrals ]"
        assert lines[2].startswith("; REMOVED by PRISM")
        assert lines[3] == "[ bonds ]"

    def test_parses_multiple_errors_same_file(self, builder):
        topol = builder.model_dir / "topol.top"
        topol.write_text(
            "line1\n"
            "line2\n"
            "   1    2    3    4     4\n"  # line 3
            "   5    6    7    8     4\n"  # line 4
            "   9   10   11   12     4\n"  # line 5
            "line6\n"
        )

        stderr = (
            "ERROR 1 [file topol.top, line 3]:\n"
            "  No default Per. Imp. Dih. types\n"
            "ERROR 2 [file topol.top, line 5]:\n"
            "  No default Per. Imp. Dih. types\n"
        )
        fixed = builder._fix_improper_from_grompp_errors(stderr)

        assert fixed == 2
        lines = topol.read_text().splitlines()
        assert lines[2].startswith("; REMOVED")  # line 3
        assert not lines[3].startswith("; REMOVED")  # line 4 untouched
        assert lines[4].startswith("; REMOVED")  # line 5

    def test_parses_itp_file_errors(self, builder):
        """Errors referencing an .itp file (not topol.top)."""
        itp = builder.model_dir / "topol_Protein_chain_A.itp"
        itp.write_text("line1\nline2\nline3\n")

        stderr = "ERROR 1 [file topol_Protein_chain_A.itp, line 2]:\n" "  No default Per. Imp. Dih. types\n"
        fixed = builder._fix_improper_from_grompp_errors(stderr)

        assert fixed == 1
        lines = itp.read_text().splitlines()
        assert lines[1].startswith("; REMOVED")

    def test_returns_zero_when_no_match(self, builder):
        stderr = "Some other GROMACS error message\n"
        assert builder._fix_improper_from_grompp_errors(stderr) == 0

    def test_returns_zero_for_empty_stderr(self, builder):
        assert builder._fix_improper_from_grompp_errors("") == 0

    def test_creates_backup(self, builder):
        topol = builder.model_dir / "topol.top"
        topol.write_text("line1\nline2\nline3\n")

        stderr = "ERROR 1 [file topol.top, line 2]:\n" "  No default Per. Imp. Dih. types\n"
        builder._fix_improper_from_grompp_errors(stderr)

        backup = builder.model_dir / "topol.top.backup"
        assert backup.exists()
        assert backup.read_text() == "line1\nline2\nline3\n"

    def test_backup_not_overwritten(self, builder):
        """Second fix pass should not overwrite existing backup."""
        topol = builder.model_dir / "topol.top"
        topol.write_text("original_line1\noriginal_line2\noriginal_line3\n")

        stderr1 = "ERROR 1 [file topol.top, line 2]:\n" "  No default Per. Imp. Dih. types\n"
        builder._fix_improper_from_grompp_errors(stderr1)

        # Now topol.top is modified; run a second fix pass
        stderr2 = "ERROR 1 [file topol.top, line 3]:\n" "  No default Per. Imp. Dih. types\n"
        builder._fix_improper_from_grompp_errors(stderr2)

        # Backup should still contain the ORIGINAL content
        backup = builder.model_dir / "topol.top.backup"
        assert "original_line1" in backup.read_text()
        assert "; REMOVED" not in backup.read_text()

    def test_skips_missing_file(self, builder):
        """Error referencing a non-existent file → skip, return 0."""
        stderr = "ERROR 1 [file nonexistent.itp, line 5]:\n" "  No default Per. Imp. Dih. types\n"
        assert builder._fix_improper_from_grompp_errors(stderr) == 0


# ===========================================================================
# TestExtractMetalsFromPdb
# ===========================================================================
class TestExtractMetalsFromPdb:
    def test_extracts_zinc(self, builder, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "N", "ALA", "A", 1) + _hetatm_line(100, "ZN", "ZN", "A", 50, x=10.0, y=20.0, z=30.0) + "END\n"
        )
        metals = builder._extract_metals_from_pdb(str(pdb))

        assert len(metals) == 1
        assert metals[0]["residue_name"] == "ZN"
        assert metals[0]["coords"] == (10.0, 20.0, 30.0)

    def test_extracts_atom_record_metal(self, builder, tmp_path):
        """Metals converted from HETATM to ATOM by cleaner should still be found."""
        pdb = tmp_path / "test.pdb"
        pdb.write_text(_atom_line(1, "ZN", "ZN", "M", 1, x=5.0, y=6.0, z=7.0) + "END\n")
        metals = builder._extract_metals_from_pdb(str(pdb))
        assert len(metals) == 1

    def test_skips_protein_residues(self, builder, tmp_path):
        """Standard amino acids must not be mistaken for metals."""
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "CA", "ALA", "A", 1)  # CA atom in ALA, not calcium
            + _atom_line(2, "CA", "GLY", "A", 2)
            + "END\n"
        )
        metals = builder._extract_metals_from_pdb(str(pdb))
        assert len(metals) == 0

    def test_skips_his_variants(self, builder, tmp_path):
        """HIE, HID, HIP are protein residues, not metals."""
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "N", "HIE", "A", 10)
            + _atom_line(2, "N", "HID", "A", 20)
            + _atom_line(3, "N", "HIP", "A", 30)
            + "END\n"
        )
        metals = builder._extract_metals_from_pdb(str(pdb))
        assert len(metals) == 0

    def test_multiple_metals(self, builder, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _hetatm_line(1, "ZN", "ZN", "A", 50, x=1.0, y=2.0, z=3.0)
            + _hetatm_line(2, "MG", "MG", "A", 51, x=4.0, y=5.0, z=6.0)
            + _hetatm_line(3, "FE", "FE", "B", 52, x=7.0, y=8.0, z=9.0)
            + "END\n"
        )
        metals = builder._extract_metals_from_pdb(str(pdb))
        assert len(metals) == 3
        names = {m["residue_name"] for m in metals}
        assert names == {"ZN", "MG", "FE"}

    def test_no_metals(self, builder, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(_atom_line(1, "N", "ALA", "A", 1) + _atom_line(2, "CA", "ALA", "A", 1) + "END\n")
        assert builder._extract_metals_from_pdb(str(pdb)) == []
