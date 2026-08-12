"""Focused regression tests for standard PTM staging."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
PTM_DIR = ROOT / "prism" / "ptm"
CLEANER_DIR = ROOT / "prism" / "utils" / "cleaner"


def _load_ptm_package():
    """Load only ``prism/ptm`` without importing PRISM's optional heavy deps."""
    package_name = "_prism_ptm_under_test"
    if package_name in sys.modules:
        return sys.modules[package_name]
    spec = importlib.util.spec_from_file_location(
        package_name,
        PTM_DIR / "__init__.py",
        submodule_search_locations=[str(PTM_DIR)],
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[package_name] = module
    spec.loader.exec_module(module)
    return module


def _load_cleaner_package():
    package_name = "_prism_cleaner_under_test"
    if package_name in sys.modules:
        return sys.modules[package_name]
    spec = importlib.util.spec_from_file_location(
        package_name,
        CLEANER_DIR / "__init__.py",
        submodule_search_locations=[str(CLEANER_DIR)],
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[package_name] = module
    spec.loader.exec_module(module)
    return module




ptm = _load_ptm_package()


def _pdb_atom(
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resid: int,
    x: float,
    y: float = 0.0,
    z: float = 0.0,
    *,
    record: str = "ATOM",
    altloc: str = "",
    icode: str = "",
    occupancy: float = 1.0,
) -> str:
    element = next((c for c in atom_name if c.isalpha()), "C")
    return (
        f"{record:<6}{serial:5d} {atom_name:>4s}{altloc:1s}{resname:>3s} {chain:1s}"
        f"{resid:4d}{icode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}{20.00:6.2f}          {element:>2s}\n"
    )


def _write_pdb(path: Path, lines: list[str]) -> Path:
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def _rtp_residue_charge(resname: str) -> float:
    rtp = ROOT / "prism" / "configs" / "forcefield" / "charmm36-jul2022.ff" / "aminoacids.rtp"
    lines = rtp.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == f"[ {resname} ]")
    atoms_start = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "[ atoms ]")
    charge = 0.0
    for line in lines[atoms_start + 1 :]:
        stripped = line.strip()
        if stripped.startswith("["):
            break
        if not stripped or stripped.startswith(";"):
            continue
        charge += float(stripped.split()[2])
    return charge


def _rtp_residue_atoms(resname: str) -> set[str]:
    rtp = ROOT / "prism" / "configs" / "forcefield" / "charmm36-jul2022.ff" / "aminoacids.rtp"
    lines = rtp.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == f"[ {resname} ]")
    atoms_start = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "[ atoms ]")
    atoms = set()
    for line in lines[atoms_start + 1 :]:
        stripped = line.strip()
        if stripped.startswith("["):
            break
        if stripped and not stripped.startswith(";"):
            atoms.add(stripped.split()[0])
    return atoms


@pytest.mark.parametrize(
    ("code", "charge", "expected_resname"),
    [
        ("SEP", -2, "SP2"),
        ("SEP", -1, "SEP"),
        ("TPO", -2, "THP2"),
        ("TPO", -1, "TPO"),
        ("PTR", -2, "TP2"),
        ("PTR", -1, "PTR"),
    ],
)
def test_charmm_phosphorylation_mapping_matches_bundled_rtp_charge(code, charge, expected_resname):
    definition = ptm.get_ptm(code)
    assert definition.charmm_name_for_charge(charge) == expected_resname
    assert _rtp_residue_charge(expected_resname) == pytest.approx(charge)


@pytest.mark.parametrize(
    ("code", "target"),
    [
        ("SEP", "SP2"),
        ("SEP", "SEP"),
        ("TPO", "THP2"),
        ("TPO", "TPO"),
        ("PTR", "TP2"),
        ("PTR", "PTR"),
        ("MLZ", "MLZ"),
        ("MLY", "MLY"),
        ("M3L", "M3L"),
        ("ALY", "ALY"),
        ("2MR", "2MR"),
    ],
)
def test_required_ptm_heavy_atoms_match_bundled_rtp_delta(code, target):
    definition = ptm.get_ptm(code)
    target_atoms = _rtp_residue_atoms(target)
    parent_atoms = _rtp_residue_atoms(definition.parent)
    added_heavy_atoms = {
        atom for atom in target_atoms - parent_atoms if not atom.startswith("H")
    }

    assert set(definition.required_heavy_atoms_for(target)) == added_heavy_atoms


def test_unsupported_explicit_charge_does_not_fall_back_to_default_residue():
    assert ptm.get_ptm("SEP").charmm_name_for_charge(-7) is None
    assert ptm.get_ptm("SEP").amber_name_for_charge(-7) is None


def test_2mr_catalog_name_matches_symmetric_dimethylarginine_rtp():
    assert ptm.get_ptm("2MR").name == "Symmetric dimethyl-arginine"


def test_default_phosphoserine_staging_uses_dianionic_rtp_and_atom_alias(tmp_path):
    source = _write_pdb(
        tmp_path / "input.pdb",
        [
            _pdb_atom(1, "o3p", "SEP", "A", 7, 0.0),
            _pdb_atom(2, "p", "SEP", "A", 7, 1.0),
            _pdb_atom(3, "o1p", "SEP", "A", 7, 2.0),
            _pdb_atom(4, "o2p", "SEP", "A", 7, 3.0),
        ],
    )
    output = tmp_path / "output.pdb"
    config = ptm.PTMConfig(
        requests=[ptm.PTMRequest(chain="A", resid=7, code="SEP")],
        disulfides="none",
    )

    ptm.PTMStager(config, "charmm36-jul2022", verbose=False).apply(
        str(source), str(output), work_dir=str(tmp_path)
    )

    atoms = [
        line
        for line in output.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    ]
    assert {atom[17:20].strip() for atom in atoms} == {"SP2"}
    assert {atom[12:16].strip() for atom in atoms} == {"OT", "P", "O1P", "O2P"}


def test_yaml_string_charge_is_normalized_to_int():
    config = ptm.parse_ptm_yaml(
        {
            "residues": [{"chain": "A", "resid": "7", "code": "sep", "charge": "-1"}],
            "disulfides": "none",
        }
    )
    assert config.requests == [ptm.PTMRequest(chain="A", resid=7, code="SEP", charge=-1)]
    assert config.requests[0].describe() == "A7 -> SEP (charge -1)"


def test_ptm_yaml_section_must_be_a_mapping():
    with pytest.raises(ValueError, match="must be a mapping"):
        ptm.parse_ptm_yaml(["A:7:SEP"])


@pytest.mark.parametrize("residues", [7, False, 0, {}, "A:7:SEP"])
def test_ptm_yaml_residues_must_be_a_sequence_of_mappings(residues):
    with pytest.raises(ValueError, match=r"ptm\.residues.*sequence"):
        ptm.parse_ptm_yaml({"residues": residues})


@pytest.mark.parametrize(
    "entry",
    [
        {"chain": "A", "resid": 7.5, "code": "SEP"},
        {"chain": "A", "resid": True, "code": "SEP"},
        {"chain": "A", "resid": 7, "code": "SEP", "charge": -1.5},
        {"chain": "A", "resid": 7, "code": "SEP", "charge": True},
    ],
)
def test_yaml_rejects_non_integer_residue_ids_and_charges(entry):
    with pytest.raises(ValueError, match="Invalid"):
        ptm.parse_ptm_yaml({"residues": [entry], "disulfides": "none"})




@pytest.mark.parametrize(("resname", "resid", "message"), [("GLY", 7, "incompatible"), ("GLY", 8, "not found")])
def test_ptm_target_must_exist_and_match_parent_or_modified_residue(tmp_path, resname, resid, message):
    source = _write_pdb(tmp_path / f"{message}.pdb", [_pdb_atom(1, "CA", resname, "A", 7, 0.0)])
    output = tmp_path / f"{message}_out.pdb"
    config = ptm.PTMConfig(
        requests=[ptm.PTMRequest(chain="A", resid=resid, code="SEP")],
        disulfides="none",
    )

    with pytest.raises(ValueError, match=message):
        ptm.PTMStager(config, "charmm36-jul2022", verbose=False).apply(str(source), str(output))


def test_ptm_target_requires_modification_heavy_atoms_before_writing_output(tmp_path):
    source = _write_pdb(
        tmp_path / "plain_serine.pdb",
        [_pdb_atom(1, "CA", "SER", "A", 7, 0.0)],
    )
    output = tmp_path / "plain_serine_out.pdb"
    config = ptm.PTMConfig(
        requests=[ptm.PTMRequest(chain="A", resid=7, code="SEP")],
        disulfides="none",
    )

    with pytest.raises(ValueError, match=r"missing required SP2 heavy atom"):
        ptm.PTMStager(config, "charmm36-jul2022", verbose=False).apply(
            str(source), str(output)
        )

    assert not output.exists()


@pytest.mark.parametrize(
    ("ptm_request", "message"),
    [
        (ptm.PTMRequest(chain="A", resid=7, code="BOG"), "Unknown PTM"),
        (ptm.PTMRequest(chain="A", resid=7, code="HYP"), "validated parameters"),
        (ptm.PTMRequest(chain="A", resid=7, code="SEP", charge=-7), "charge -7"),
    ],
)
def test_invalid_ptm_requests_fail_before_writing_output(tmp_path, ptm_request, message):
    source = _write_pdb(tmp_path / "invalid_request.pdb", [_pdb_atom(1, "CA", "SER", "A", 7, 0.0)])
    output = tmp_path / "invalid_request_out.pdb"

    with pytest.raises(ValueError, match=message):
        ptm.PTMStager(
            ptm.PTMConfig(requests=[ptm_request], disulfides="none"),
            "charmm36-jul2022",
            verbose=False,
        ).apply(str(source), str(output))

    assert not output.exists()


@pytest.mark.parametrize(
    ("code", "parent"),
    [("SEP", "SER"), ("TPO", "THR"), ("PTR", "TYR")],
)
def test_amber_monoanionic_phosphorylation_variants_fail_route_preflight(tmp_path, code, parent):
    source = _write_pdb(tmp_path / f"{code}.pdb", [_pdb_atom(1, "CA", parent, "A", 7, 0.0)])
    output = tmp_path / f"{code}_amber_out.pdb"

    with pytest.raises(RuntimeError, match=r"tleap\+"):
        ptm.PTMStager(
            ptm.PTMConfig(
                requests=[ptm.PTMRequest(chain="A", resid=7, code=code, charge=-1)],
                disulfides="none",
            ),
            "amber14sb",
            verbose=False,
        ).apply(str(source), str(output))

    assert not output.exists()


def test_requested_hetatm_is_promoted_before_cleaning_without_changing_order(tmp_path):
    source = _write_pdb(
        tmp_path / "raw.pdb",
        [
            _pdb_atom(1, "CA", "ALA", "A", 1, 0.0),
            _pdb_atom(2, "P", "SEP", "A", 2, 1.0, record="HETATM"),
            _pdb_atom(3, "CA", "GLY", "A", 3, 2.0),
        ],
    )
    output = tmp_path / "preclean.pdb"
    config = ptm.PTMConfig(
        requests=[ptm.PTMRequest(chain="A", resid=2, code="SEP")],
        disulfides="none",
    )

    promoted = ptm.promote_requested_ptm_records(str(source), str(output), config)

    coordinate_lines = [
        line for line in output.read_text(encoding="utf-8").splitlines() if line.startswith(("ATOM", "HETATM"))
    ]
    assert promoted == 1
    assert [line[:6].strip() for line in coordinate_lines] == ["ATOM", "ATOM", "ATOM"]
    assert [int(line[22:26]) for line in coordinate_lines] == [1, 2, 3]

    cleaned = tmp_path / "cleaned.pdb"
    cleaner = _load_cleaner_package().ProteinCleaner(verbose=False)
    cleaner.clean_pdb(str(output), str(cleaned))
    cleaned_lines = cleaned.read_text(encoding="utf-8").splitlines()
    first_ter = cleaned_lines.index("TER")
    protein_lines = [line for line in cleaned_lines[:first_ter] if line.startswith("ATOM")]
    assert [int(line[22:26]) for line in protein_lines] == [1, 2, 3]


def _metal_pdb_atom(rec, serial, atom, resname, chain, resid, x, y, z, elem):
    """Emit a fixed-column PDB coordinate line.

    The residue name is written into the 4-character field (columns 18-21) so
    both 3- and 4-letter names (e.g. THP2) land in the same columns.
    """
    return (
        f"{rec:<6}{serial:>5} {atom:<4} {resname:<4}{chain:1}{resid:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}{'':>10}{elem:>2}\n"
    )


def test_ptm_residue_atoms_are_not_extracted_as_metals(tmp_path):
    """PTM residues are built into the protein chain by pdb2gmx, so they must
    never be treated as metal ions.

    Their standard atom names Calpha ``CA`` and Cdelta ``CD`` collide with the
    calcium/cadmium element symbols. Without an explicit skip the metal
    processor emits e.g. ``M3L`` / ``TP2`` into the topology ``[ molecules ]``
    as a nonexistent moleculetype, which breaks ``gmx grompp`` at the ion step
    for every PTM system. This is a regression guard for that fix.
    """
    from prism.utils.system.metals import MetalProcessorMixin

    ptm_pdb = tmp_path / "ptm.pdb"
    ptm_pdb.write_text(
        _metal_pdb_atom("ATOM", 1, "N", "M3L", "P", 4, 36.201, 16.048, 59.925, "N")
        + _metal_pdb_atom("ATOM", 2, "CA", "M3L", "P", 4, 37.428, 16.832, 59.852, "C")
        + _metal_pdb_atom("ATOM", 3, "CD", "M3L", "P", 4, 39.127, 15.847, 56.462, "C")
        + _metal_pdb_atom("ATOM", 4, "CM1", "M3L", "P", 4, 41.476, 15.566, 54.559, "C")
        # phosphotyrosine (dianion) also carries a Calpha "CA"
        + _metal_pdb_atom("ATOM", 5, "CA", "TP2", "A", 1158, 10.0, 10.0, 10.0, "C")
        # phosphothreonine dianion THP2 is a 4-character name; its Calpha/Cdelta
        # must still be skipped when the resname is read from the 4-char field.
        + _metal_pdb_atom("ATOM", 6, "CA", "THP2", "A", 160, 20.0, 20.0, 20.0, "C")
        + _metal_pdb_atom("ATOM", 7, "CD", "THP2", "A", 160, 21.0, 20.0, 20.0, "C")
        + "END\n",
        encoding="utf-8",
    )
    assert MetalProcessorMixin()._extract_metals_from_pdb(str(ptm_pdb)) == []

    # A genuine metal ion must STILL be detected (guard against over-skipping).
    zn_pdb = tmp_path / "zn.pdb"
    zn_pdb.write_text(
        _metal_pdb_atom("HETATM", 6, "ZN", "ZN", "A", 301, 20.0, 20.0, 20.0, "ZN")
        + "END\n",
        encoding="utf-8",
    )
    assert [m["residue_name"] for m in MetalProcessorMixin()._extract_metals_from_pdb(str(zn_pdb))] == ["ZN"]


def test_cli_ptm_build_allows_protein_only(monkeypatch, tmp_path):
    """A protein-only PTM build must be accepted by the CLI without a ligand,
    mirroring the membrane protein-only path. Regression guard for the
    ``--ptm``/``--ssbond`` relaxation of the "ligand is required" check.
    """
    import sys as _sys

    from prism.builder import cli as builder_cli

    captured = {}

    class FakeBuilder:
        def __init__(self, protein_path, ligand_paths, output_dir, **kwargs):
            captured.update(
                protein_path=protein_path,
                ligand_paths=ligand_paths,
                output_dir=output_dir,
                **kwargs,
            )

        def run(self):
            return "ok"

    monkeypatch.setattr(builder_cli, "PRISMBuilder", FakeBuilder)
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM\n", encoding="utf-8")
    monkeypatch.setattr(
        _sys,
        "argv",
        [
            "prism",
            str(protein),
            "--forcefield",
            "charmm36-jul2022",
            "--ptm",
            "A:145:SEP",
            "-o",
            str(tmp_path / "out"),
        ],
    )

    builder_cli.main()

    assert captured["ligand_paths"] == []
    assert captured["ptm_config"] is not None
    assert [r.code for r in captured["ptm_config"].requests] == ["SEP"]
