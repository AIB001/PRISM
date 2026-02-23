"""Tests for prism.utils.protonation — PropkaProtonator class."""

import pytest
from unittest.mock import patch, MagicMock


# ---------------------------------------------------------------------------
# Helper: build a minimal PDB ATOM line (80 chars, fixed-width)
# ---------------------------------------------------------------------------
def _atom_line(serial, name, resname, chain, resseq, x=0.0, y=0.0, z=0.0):
    """Build a correctly formatted PDB ATOM record."""
    return (
        f"ATOM  {serial:5d} {name:<4s} {resname:<3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]:>2s}  \n"
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def pdb_with_his(tmp_path):
    """PDB file containing two HIS residues on chain A (res 10 and 20)."""
    pdb = tmp_path / "protein.pdb"
    lines = [
        _atom_line(1, "N", "ALA", "A", 1),
        _atom_line(2, "CA", "ALA", "A", 1),
        _atom_line(3, "N", "HIS", "A", 10),
        _atom_line(4, "CA", "HIS", "A", 10),
        _atom_line(5, "N", "HIS", "A", 20),
        _atom_line(6, "CA", "HIS", "A", 20),
        _atom_line(7, "N", "GLY", "A", 30),
        "END\n",
    ]
    pdb.write_text("".join(lines))
    return pdb


@pytest.fixture
def pdb_no_his(tmp_path):
    """PDB file with zero histidines."""
    pdb = tmp_path / "nohis.pdb"
    lines = [
        _atom_line(1, "N", "ALA", "A", 1),
        _atom_line(2, "CA", "GLY", "A", 2),
        "END\n",
    ]
    pdb.write_text("".join(lines))
    return pdb


@pytest.fixture
def pdb_mixed_his(tmp_path):
    """PDB with one HIS and one already-renamed HIE."""
    pdb = tmp_path / "mixed.pdb"
    lines = [
        _atom_line(1, "N", "HIS", "A", 10),  # should be renamed
        _atom_line(2, "N", "HIE", "A", 20),  # already renamed, skip
        "END\n",
    ]
    pdb.write_text("".join(lines))
    return pdb


# ---------------------------------------------------------------------------
# TestPropkaProtonatorInit
# ---------------------------------------------------------------------------
class TestPropkaProtonatorInit:
    """Constructor and dependency check."""

    def test_default_ph(self):
        from prism.utils.protonation import PropkaProtonator

        p = PropkaProtonator()
        assert p.ph == 7.0
        assert p.verbose is False

    def test_custom_ph(self):
        from prism.utils.protonation import PropkaProtonator

        p = PropkaProtonator(ph=4.5, verbose=True)
        assert p.ph == 4.5
        assert p.verbose is True

    def test_import_error_when_propka_missing(self):
        """If propka is not installed, _check_propka_available should raise."""
        from prism.utils.protonation import PropkaProtonator

        obj = object.__new__(PropkaProtonator)
        with patch("builtins.__import__", side_effect=ImportError("No module")):
            with pytest.raises(ImportError, match="PROPKA is not available"):
                obj._check_propka_available()


# ---------------------------------------------------------------------------
# TestRenameHistidines
# ---------------------------------------------------------------------------
class TestRenameHistidines:
    """PDB renaming logic in rename_histidines (mock predict_his_states)."""

    def _make_protonator(self, his_states):
        """Create protonator with mocked predict_his_states."""
        from prism.utils.protonation import PropkaProtonator

        p = PropkaProtonator(ph=7.0)
        p.predict_his_states = MagicMock(return_value=his_states)
        return p

    def test_renames_his_to_hie(self, pdb_with_his, tmp_path):
        states = {("A", "10"): "HIE", ("A", "20"): "HIE"}
        p = self._make_protonator(states)
        out = tmp_path / "out.pdb"
        stats = p.rename_histidines(str(pdb_with_his), str(out))

        content = out.read_text()
        assert "HIE" in content
        assert stats["total_his"] == 2
        assert len(stats["renamed"]) == 2

    def test_renames_his_to_hip(self, pdb_with_his, tmp_path):
        states = {("A", "10"): "HIP", ("A", "20"): "HIP"}
        p = self._make_protonator(states)
        out = tmp_path / "out.pdb"
        stats = p.rename_histidines(str(pdb_with_his), str(out))

        content = out.read_text()
        assert "HIP" in content
        # ALA and GLY lines should not contain HIP
        for line in content.splitlines():
            if "ALA" in line or "GLY" in line:
                assert "HIP" not in line

    def test_skips_non_his_residues(self, pdb_with_his, tmp_path):
        states = {("A", "10"): "HIE", ("A", "20"): "HIP"}
        p = self._make_protonator(states)
        out = tmp_path / "out.pdb"
        p.rename_histidines(str(pdb_with_his), str(out))

        # ALA at residue 1 and GLY at residue 30 should be untouched
        for line in out.read_text().splitlines():
            if line.startswith("ATOM") and "   1 " in line[22:27]:
                assert line[17:20] == "ALA"

    def test_skips_already_renamed_hie(self, pdb_mixed_his, tmp_path):
        """HIE residues are not matched by line[17:20]=='HIS', so skip."""
        states = {("A", "10"): "HIP"}  # only res 10 (HIS) should be renamed
        p = self._make_protonator(states)
        out = tmp_path / "out.pdb"
        stats = p.rename_histidines(str(pdb_mixed_his), str(out))

        assert stats["renamed"] == {("A", "10"): "HIP"}
        # Res 20 was already HIE, should remain HIE (not HIP)
        for line in out.read_text().splitlines():
            if line.startswith("ATOM") and line[22:26].strip() == "20":
                assert line[17:20] == "HIE"

    def test_no_histidines_returns_empty(self, pdb_no_his, tmp_path):
        p = self._make_protonator({})
        out = tmp_path / "out.pdb"
        stats = p.rename_histidines(str(pdb_no_his), str(out))

        assert stats["total_his"] == 0
        assert stats["renamed"] == {}
        assert out.read_text() == pdb_no_his.read_text()

    def test_in_place_rename(self, pdb_with_his):
        """When input == output, file is modified in place."""
        states = {("A", "10"): "HIP", ("A", "20"): "HIE"}
        p = self._make_protonator(states)

        original_content = pdb_with_his.read_text()
        stats = p.rename_histidines(str(pdb_with_his), str(pdb_with_his))

        new_content = pdb_with_his.read_text()
        assert new_content != original_content
        assert "HIP" in new_content
        assert "HIE" in new_content
        assert len(stats["renamed"]) == 2

    def test_stats_structure(self, pdb_with_his, tmp_path):
        states = {("A", "10"): "HIE", ("A", "20"): "HIP"}
        p = self._make_protonator(states)
        out = tmp_path / "out.pdb"
        stats = p.rename_histidines(str(pdb_with_his), str(out))

        assert "total_his" in stats
        assert "renamed" in stats
        assert "states" in stats
        assert stats["states"] is states
