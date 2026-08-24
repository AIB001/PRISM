"""Calibration guards for the quality-control thresholds.

Every threshold in ``prism.generation.quality`` is a claim about where correct
chemistry ends and a defect begins. The only way to keep such a claim honest is
to test it against molecules that are known to be right, so this file holds two
answer keys:

* experimentally determined ligands from the PDB, which must pass every
  geometric check with no flags at all. The set deliberately includes the
  upper tail of what real structures show -- 2XNB reaches a 26.1 deg bond
  angle deviation and 1D3P a 0.035 A aromatic ripple -- because a set of only
  well-behaved ligands would not notice a threshold tightened into the range
  real crystallography occupies;
* approved drugs, which must not be dismissed as reactive.

Both were added after a measured failure: the reactive-group SMARTS originally
flagged 5 of 20 marketed drugs (amidines and guanidines read as Schiff bases,
hydroxamic acids as labile N-O, hydrazides as free hydrazines). Nothing had
caught it because the checks were only ever run against generated molecules,
where a false positive is indistinguishable from a real defect.
"""

from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem, RDLogger

from prism.generation.quality import QualityContext, REACTIVE_SMARTS, assess

RDLogger.DisableLog("rdApp.*")

CRYSTAL_LIGANDS = Path(__file__).parent / "data" / "crystal_ligands"

# Marketed drugs whose motifs are easy to mistake for reactive groups.
APPROVED_DRUGS = {
    "metformin": "CN(C)C(=N)NC(=N)N",                       # guanidine
    "benzamidine": "NC(=N)c1ccccc1",                        # amidine
    "vorinostat": "ONC(=O)CCCCCCC(=O)Nc1ccccc1",            # hydroxamic acid
    "isoniazid": "NNC(=O)c1ccncc1",                         # hydrazide
    "celecoxib": "Cc1ccc(cc1)-c1cc(nn1-c1ccc(cc1)S(N)(=O)=O)C(F)(F)F",  # pyrazole
    "fluconazole": "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F",    # triazoles
    "warfarin": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",    # 4-hydroxycoumarin
    "ranitidine": "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1",
    "sildenafil": "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1",
    "omeprazole": "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "ciprofloxacin": "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
    "losartan": "CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)-c1ccccc1-c1nnn[nH]1",
    "allopurinol": "O=c1[nH]cnc2[nH]ncc12",
    "tamoxifen": "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
}

# Genuinely labile motifs that must still be caught.
TRUE_POSITIVES = {
    "schiff_base": ("C(=Nc1ccccc1)c1ccccc1", "acyclic_imine"),
    "ketimine": ("CC(C)=NC", "acyclic_imine"),
    "free_hydrazine": ("NNCc1ccccc1", "hydrazine"),
    "o_methylhydroxylamine": ("CON", "hydroxylamine"),
    "allene": ("C=C=Cc1ccccc1", "allene"),
    "acyl_chloride": ("O=C(Cl)c1ccccc1", "acyl_halide"),
    "epoxide": ("C1OC1c1ccccc1", "epoxide_aziridine"),
    "disulfide": ("CSSCc1ccccc1", "disulfide"),
}


def _patterns():
    compiled = {}
    for name, smarts in REACTIVE_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        assert pattern is not None, f"{name} has invalid SMARTS: {smarts}"
        compiled[name] = pattern
    return compiled


def _crystal_ligands():
    paths = sorted(CRYSTAL_LIGANDS.glob("*.sdf"))
    assert paths, f"no crystal ligands in {CRYSTAL_LIGANDS}"
    return paths


def test_every_reactive_smarts_compiles():
    _patterns()


@pytest.mark.parametrize("path", _crystal_ligands(), ids=lambda p: p.stem)
def test_crystal_ligands_pass_every_check(path):
    """An experimentally determined pose is the answer key: it must be clean."""
    report = assess(path.read_text(encoding="utf-8"), QualityContext(level="standard"))

    assert report.status == "PASS", f"{path.stem}: {report.flags}"
    assert report.flags == []
    assert report.sanitize_status == "passed"
    assert report.integrity["fragments"] == 1
    assert report.integrity["dummy_atoms"] == 0
    assert report.geometry["bad_bond_lengths"] == 0
    assert report.geometry["internal_clashes"] == 0
    # The angle threshold has to sit above what real structures actually show.
    assert report.geometry["worst_angle_deviation_deg"] < 35.0
    assert report.geometry["aromatic_planarity_rms_a"] < 0.10
    assert report.plausibility["unspecified_stereocenters"] == 0


@pytest.mark.parametrize("path", _crystal_ligands(), ids=lambda p: p.stem)
def test_crystal_ligands_accept_hydrogens_without_moving(path):
    from prism.generation.quality import hydrogenate_record

    block, info = hydrogenate_record(path.read_text(encoding="utf-8"))

    assert block is not None
    assert info["status"] == "added"
    assert info["hydrogens_added"] > 0
    # Adding hydrogens must never perturb the experimental pose.
    assert info["heavy_atom_max_displacement_a"] == 0.0


@pytest.mark.parametrize("name,smiles", sorted(APPROVED_DRUGS.items()))
def test_approved_drugs_are_not_called_reactive(name, smiles):
    """A pattern that dismisses a marketed drug is too broad to be useful."""
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None, name

    hits = sorted(
        group for group, pattern in _patterns().items()
        if molecule.HasSubstructMatch(pattern)
    )
    # Hydralazine is a genuine free arylhydrazine and is deliberately absent
    # from this set; nothing else here should match.
    assert hits == [], f"{name} wrongly flagged as {hits}"


@pytest.mark.parametrize("name,case", sorted(TRUE_POSITIVES.items()))
def test_genuinely_labile_motifs_are_still_caught(name, case):
    """Tightening the patterns must not blind them to real hazards."""
    smiles, expected = case
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None, name

    pattern = _patterns()[expected]
    assert molecule.HasSubstructMatch(pattern), f"{name} no longer matches {expected}"
