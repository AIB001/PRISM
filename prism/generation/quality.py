"""Quality control for generated ligands.

Policy
------
Generation models hand PRISM molecules of very different provenance.  Some
emit only heavy-atom coordinates and let an upstream perception step guess the
bonds; others emit bond orders and formal charges directly.  The defects differ
accordingly, so this module separates what may be *applied* from what may only
be *reported*:

``apply``
    Transformations that cannot change the molecular graph or the generated
    pose.  Only SDF data-field normalization and hydrogen addition qualify, and
    hydrogen addition is always written to a separate file.

``detect``
    Everything else, **including speculative valence repair**.  When a molecule
    fails sanitization this module searches for minimal edits that would make it
    valid and records them as a proposal with every equally-ranked variant.  The
    proposal is never substituted into the candidate stream: a bond order that
    sanitizes is not evidence that it is the bond order the model intended.

RDKit is an optional PRISM dependency.  When it is missing every RDKit-backed
check degrades to ``NOT_CHECKED`` instead of failing the run.

Thresholds
----------
Default thresholds are calibrated against experimentally determined ligands
rather than textbook ideals.  For the 8c7y crystal ligand the observed values
are: worst bond angle deviation 15.9 deg, aromatic ring planarity RMS 0.019 A,
closest protein heavy-atom contact 2.84 A.  A check that flags a crystal pose is
too strict to be useful, so the defaults sit well outside those values.
"""

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


# Ordered from healthiest to worst; ``worst_status`` relies on this order.
STATUS_ORDER = (
    "PASS",
    "WARN_SDF_FORMAT",
    "WARN_PLAUSIBILITY",
    "WARN_GEOMETRY",
    "WARN_POSE",
    "FAIL_STRUCTURE",
    "FAIL_CHEMISTRY",
    "FAIL_PARSE",
)
PASSING_STATUSES = frozenset(status for status in STATUS_ORDER if not status.startswith("FAIL"))

QC_LEVELS = ("off", "basic", "standard", "full")

DEFAULT_THRESHOLDS: Dict[str, float] = {
    # Bond length as a multiple of the summed covalent radii.
    "min_bond_ratio": 0.60,
    "max_bond_ratio": 1.45,
    # Non-bonded intramolecular overlap, also as a covalent-radius multiple.
    "clash_ratio": 0.55,
    # Bond angle deviation from the hybridization ideal, in degrees.  The
    # crystal reference reaches 15.9 deg, so 20 deg counts an outlier and 35 deg
    # promotes the molecule to a geometry warning.
    "angle_outlier_deg": 20.0,
    "angle_warn_deg": 35.0,
    # RMS out-of-plane deviation of aromatic ring atoms, in Angstrom.
    "aromatic_planarity_a": 0.10,
    # Ligand heavy atom to protein heavy atom distance, in Angstrom.
    "protein_hard_clash_a": 2.0,
    "protein_close_contact_a": 2.8,
    # Ligand centroid displacement from the requested pocket center.
    "pocket_offset_warn_a": 2.5,
}

# Hydrolytically or oxidatively unstable groups, plus perception artifacts that
# are common in geometry-derived bond orders (allene, acyclic imine, enol).
# Reported only: whether any of these disqualifies a candidate is a project
# policy decision, not a correctness question.
REACTIVE_SMARTS: Dict[str, str] = {
    "acyl_halide": "[CX3](=O)[F,Cl,Br,I]",
    "anhydride": "[CX3](=O)[OX2][CX3]=O",
    "peroxide": "[OX2][OX2]",
    "epoxide_aziridine": "[C;r3][O,N;r3]",
    "azide": "[N-]=[N+]=N",
    "isocyanate": "[NX2]=[CX2]=[OX1]",
    "n_nitroso": "[NX3][NX2]=O",
    # Free hydrazine only. An acyl-stabilized N-N is a hydrazide, which is a
    # marketed-drug motif (isoniazid), and an aromatic N-N is a pyrazole or
    # triazole ring rather than a reactive group.
    "hydrazine": "[NX3;!$(N[CX3]=[OX1]);!a][NX3;!$(N[CX3]=[OX1]);!a]",
    "gem_diol": "[CX4]([OX2H])[OX2H]",
    "hemiketal": "[CX4;!$(C=O)]([OX2H1])[OX2,NX3]",
    "aminal": "[CX4]([NX3])[NX3]",
    "enol": "[CX3]=[CX3][OX2H]",
    "ketene": "[CX3]=[CX2]=[OX1]",
    # Schiff base only. Requiring the imine carbon to carry no N/O/S excludes
    # amidines, guanidines and imidates -- benzamidine is the classic thrombin
    # warhead and metformin is a guanidine, neither of which is labile.
    "acyclic_imine": "[CX3;!R;!$(C~[NX3]);!$(C~[OX2]);!$(C~[SX2])]=[NX2;!R;!$(N~[!#6;!#1])]",
    "thiol": "[SX2H]",
    "disulfide": "[SX2][SX2]",
    "alpha_halo_carbonyl": "[CX3](=O)[CX4][F,Cl,Br,I]",
    "allene": "[CX3]=[CX2]=[CX3]",
    "alkyne_in_ring": "[CX2;R]#[CX2;R]",
    "heteroatom_in_4_ring": "[O,N,S;r4]",
    # Free hydroxylamine only. A bare N-O match also flags hydroxamic acids,
    # which are an approved-drug pharmacophore (vorinostat), so the acyl-
    # stabilized form is excluded.
    "hydroxylamine": "[NX3;!$(N[CX3]=[OX1,SX1]);!$(N=*);!R][OX2H1,OX2H0;!$(O[CX3]=O)]",
}


def rdkit_available() -> bool:
    """Return True when the optional RDKit dependency can be imported."""
    try:
        import rdkit  # noqa: F401
        from rdkit import Chem  # noqa: F401
    except (ImportError, ModuleNotFoundError):
        return False
    return True


def worst_status(statuses: Sequence[str]) -> str:
    known = [status for status in statuses if status in STATUS_ORDER]
    if not known:
        return "NOT_CHECKED"
    return max(known, key=STATUS_ORDER.index)


@dataclass
class RepairProposal:
    """A minimal-edit repair that PRISM deliberately does not apply."""

    confidence: str
    reason: str
    variants: List[Dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "confidence": self.confidence,
            "reason": self.reason,
            # Carried on the proposal itself so every consumer -- manifest,
            # CSV, MCP response -- sees that PRISM did not apply it.
            "applied": False,
            "variant_count": len(self.variants),
            "variants": [
                {key: value for key, value in variant.items() if key != "molblock"}
                for variant in self.variants
            ],
        }


@dataclass
class QualityContext:
    """Inputs and settings shared by every candidate in one model run."""

    level: str = "full"
    protein_coordinates: Optional[Any] = None
    pocket_center: Optional[Tuple[float, float, float]] = None
    thresholds: Dict[str, float] = field(default_factory=lambda: dict(DEFAULT_THRESHOLDS))

    def __post_init__(self):
        if self.level not in QC_LEVELS:
            raise ValueError(f"Unknown QC level '{self.level}'; expected one of {QC_LEVELS}")
        merged = dict(DEFAULT_THRESHOLDS)
        merged.update(self.thresholds or {})
        self.thresholds = merged

    @property
    def wants_geometry(self) -> bool:
        return self.level in {"standard", "full"}

    @property
    def wants_plausibility(self) -> bool:
        return self.level in {"standard", "full"}

    @property
    def wants_pose(self) -> bool:
        return self.level == "full" and self.protein_coordinates is not None


@dataclass
class QualityReport:
    status: str = "NOT_CHECKED"
    flags: List[str] = field(default_factory=list)
    parse_mode: str = "not_checked"
    sanitize_status: str = "not_checked"
    integrity: Dict[str, Any] = field(default_factory=dict)
    hydrogen: Dict[str, Any] = field(default_factory=dict)
    geometry: Dict[str, Any] = field(default_factory=dict)
    plausibility: Dict[str, Any] = field(default_factory=dict)
    pose: Dict[str, Any] = field(default_factory=dict)
    repair: Optional[RepairProposal] = None
    canonical_smiles: str = ""

    @property
    def usable(self) -> bool:
        return self.status in PASSING_STATUSES or self.status == "NOT_CHECKED"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "status": self.status,
            "flags": list(self.flags),
            "parse_mode": self.parse_mode,
            "sanitize_status": self.sanitize_status,
            "canonical_smiles": self.canonical_smiles,
            "integrity": dict(self.integrity),
            "hydrogen": dict(self.hydrogen),
            "geometry": dict(self.geometry),
            "plausibility": dict(self.plausibility),
            "pose": dict(self.pose),
            "repair": self.repair.to_dict() if self.repair is not None else None,
        }


# --------------------------------------------------------------------------
# Text-level checks that do not require RDKit
# --------------------------------------------------------------------------


def sdf_data_fields_well_formed(record: str) -> bool:
    """Every ``> <TAG>`` block must be terminated by a blank line.

    FLOWR emits property values without the terminating blank line the CTfile
    specification requires.  A single record still reads, but concatenating the
    records breaks sequential parsing, so PRISM normalizes the file.
    """
    lines = record.splitlines()
    for index, line in enumerate(lines):
        if not line.startswith(">"):
            continue
        terminated = False
        for following in lines[index + 1 :]:
            if following == "":
                terminated = True
                break
            if following.startswith(">") or following == "$$$$":
                break
        if not terminated:
            return False
    return True


def declared_counts(record: str) -> Tuple[Optional[int], Optional[int]]:
    """Read the atom/bond counts from the V2000 counts line."""
    lines = record.splitlines()
    if len(lines) < 4 or "V3000" in lines[3]:
        return None, None
    try:
        return int(lines[3][0:3]), int(lines[3][3:6])
    except ValueError:
        return None, None


# --------------------------------------------------------------------------
# RDKit-backed parsing and sanitization
# --------------------------------------------------------------------------


def _mol_block(record: str) -> str:
    return record.rsplit("$$$$", 1)[0]


def parse_unsanitized(record: str):
    """Parse without sanitization so broken molecules stay inspectable."""
    from rdkit import Chem, rdBase

    block = _mol_block(record)
    with rdBase.BlockLogs():
        mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False, strictParsing=True)
    if mol is not None:
        return mol, "strict"
    with rdBase.BlockLogs():
        mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False, strictParsing=False)
    return mol, ("lenient" if mol is not None else "failed")


def sanitize(mol):
    """Sanitize a copy and return the decoded failure reason on error.

    ``Chem.MolFromMolBlock(sanitize=True)`` only reports ``None``.  Decoding the
    ``SANITIZE_*`` flags tells the user whether the molecule failed on valence,
    kekulization, or aromaticity, which is what makes a repair proposal
    reviewable.
    """
    from rdkit import Chem, rdBase

    candidate = Chem.Mol(mol)
    try:
        with rdBase.BlockLogs():
            status = Chem.SanitizeMol(candidate, catchErrors=True)
    except Exception as exc:  # RDKit raises different types across releases
        return None, f"exception:{type(exc).__name__}"
    if status == Chem.SanitizeFlags.SANITIZE_NONE:
        return candidate, "passed"
    names = []
    value = int(status)
    for name in dir(Chem.SanitizeFlags):
        if not name.startswith("SANITIZE_") or name in {"SANITIZE_ALL", "SANITIZE_NONE"}:
            continue
        if value & int(getattr(Chem.SanitizeFlags, name)):
            names.append(name)
    return None, "+".join(names) if names else str(value)


def _positions(mol) -> List[Tuple[float, float, float]]:
    if mol.GetNumConformers() == 0:
        return []
    conformer = mol.GetConformer()
    points = []
    for index in range(mol.GetNumAtoms()):
        position = conformer.GetAtomPosition(index)
        points.append((position.x, position.y, position.z))
    return points


# --------------------------------------------------------------------------
# Structural integrity
# --------------------------------------------------------------------------


def integrity_checks(mol, record: str) -> Dict[str, Any]:
    from rdkit import Chem

    points = _positions(mol)
    finite = bool(points) and all(math.isfinite(value) for point in points for value in point)
    lines = record.splitlines()
    # A single atom, or a planar fragment, has no z spread. The molfile
    # dimensional code on line 2 is the authoritative marker in that case.
    dimension_marked = len(lines) > 1 and "3D" in lines[1].upper()
    duplicates = 0
    if finite:
        for first in range(len(points)):
            for second in range(first + 1, len(points)):
                if math.dist(points[first], points[second]) < 0.10:
                    duplicates += 1
    atom_count = mol.GetNumAtoms()
    bond_count = mol.GetNumBonds()
    declared_atoms, declared_bonds = declared_counts(record)
    return {
        "atom_count": atom_count,
        "bond_count": bond_count,
        "declared_atoms": declared_atoms,
        "declared_bonds": declared_bonds,
        "counts_match": declared_atoms in {None, atom_count}
        and declared_bonds in {None, bond_count},
        "fragments": len(Chem.GetMolFrags(mol, sanitizeFrags=False)),
        "dummy_atoms": sum(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()),
        "radical_electrons": sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()),
        "explicit_h": sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()),
        "finite_coordinates": finite,
        "duplicate_coordinate_pairs": duplicates,
        "is_3d": finite
        and (
            dimension_marked
            or (max(p[2] for p in points) - min(p[2] for p in points)) > 1e-3
        ),
    }


def _structurally_complete(integrity: Dict[str, Any]) -> bool:
    # Bond count is not asserted directly: a lone atom legitimately has none,
    # and any other missing bond shows up as more than one fragment.
    return bool(
        integrity["counts_match"]
        and integrity["atom_count"] > 0
        and integrity["finite_coordinates"]
        and integrity["is_3d"]
        and integrity["fragments"] == 1
        and integrity["duplicate_coordinate_pairs"] == 0
        and integrity["dummy_atoms"] == 0
    )


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------


def geometry_checks(mol, thresholds: Dict[str, float]) -> Dict[str, Any]:
    """Bond lengths, non-bonded overlap, bond angles, aromatic planarity.

    Bond lengths alone are not sufficient: geometry-perception models can place
    every bond at a plausible distance while leaving the angles badly distorted.
    """
    from rdkit import Chem

    points = _positions(mol)
    result: Dict[str, Any] = {
        "bad_bond_lengths": 0,
        "min_bond_ratio": None,
        "max_bond_ratio": None,
        "internal_clashes": 0,
        "angle_outliers": 0,
        "worst_angle_deviation_deg": None,
        "aromatic_planarity_rms_a": None,
    }
    if not points or not all(math.isfinite(v) for p in points for v in p):
        return result

    table = Chem.GetPeriodicTable()
    bonded = set()
    ratios = []
    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonded.add((min(begin, end), max(begin, end)))
        radius = table.GetRcovalent(mol.GetAtomWithIdx(begin).GetAtomicNum()) + table.GetRcovalent(
            mol.GetAtomWithIdx(end).GetAtomicNum()
        )
        if radius <= 0:
            continue
        ratio = math.dist(points[begin], points[end]) / radius
        ratios.append(ratio)
        if ratio < thresholds["min_bond_ratio"] or ratio > thresholds["max_bond_ratio"]:
            result["bad_bond_lengths"] += 1
    if ratios:
        result["min_bond_ratio"] = round(min(ratios), 4)
        result["max_bond_ratio"] = round(max(ratios), 4)

    clashes = 0
    for first in range(mol.GetNumAtoms()):
        first_radius = table.GetRcovalent(mol.GetAtomWithIdx(first).GetAtomicNum())
        for second in range(first + 1, mol.GetNumAtoms()):
            if (first, second) in bonded:
                continue
            second_radius = table.GetRcovalent(mol.GetAtomWithIdx(second).GetAtomicNum())
            limit = thresholds["clash_ratio"] * (first_radius + second_radius)
            if limit > 0 and math.dist(points[first], points[second]) < limit:
                clashes += 1
    result["internal_clashes"] = clashes

    ideal = {
        Chem.HybridizationType.SP: 180.0,
        Chem.HybridizationType.SP2: 120.0,
        Chem.HybridizationType.SP3: 109.5,
    }
    rings = mol.GetRingInfo().AtomRings()
    outliers = 0
    worst = 0.0
    for atom in mol.GetAtoms():
        neighbors = [item.GetIdx() for item in atom.GetNeighbors()]
        if len(neighbors) < 2:
            continue
        target = ideal.get(atom.GetHybridization())
        if target is None:
            continue
        member_of = [len(ring) for ring in rings if atom.GetIdx() in ring]
        if member_of and min(member_of) < 5:
            # Three- and four-membered rings are legitimately strained.
            continue
        origin = points[atom.GetIdx()]
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                first = [points[neighbors[i]][k] - origin[k] for k in range(3)]
                second = [points[neighbors[j]][k] - origin[k] for k in range(3)]
                norm = math.sqrt(sum(v * v for v in first)) * math.sqrt(sum(v * v for v in second))
                if norm <= 0:
                    continue
                cosine = sum(a * b for a, b in zip(first, second)) / norm
                angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
                deviation = abs(angle - target)
                worst = max(worst, deviation)
                if deviation > thresholds["angle_outlier_deg"]:
                    outliers += 1
    result["angle_outliers"] = outliers
    result["worst_angle_deviation_deg"] = round(worst, 2)

    planarity = 0.0
    for ring in rings:
        if not all(mol.GetAtomWithIdx(index).GetIsAromatic() for index in ring):
            continue
        planarity = max(planarity, _ring_planarity_rms([points[index] for index in ring]))
    result["aromatic_planarity_rms_a"] = round(planarity, 4)
    return result


def _ring_planarity_rms(ring_points: Sequence[Tuple[float, float, float]]) -> float:
    """RMS distance of ring atoms from their best-fit plane."""
    import numpy as np

    matrix = np.asarray(ring_points, dtype=float)
    matrix = matrix - matrix.mean(axis=0)
    try:
        normal = np.linalg.svd(matrix)[2][-1]
    except np.linalg.LinAlgError:
        return 0.0
    return float(np.sqrt(np.mean((matrix @ normal) ** 2)))


# --------------------------------------------------------------------------
# Chemical plausibility
# --------------------------------------------------------------------------

_FILTER_CATALOG = None
_REACTIVE_PATTERNS = None


def _filter_catalog():
    global _FILTER_CATALOG
    if _FILTER_CATALOG is None:
        from rdkit.Chem import FilterCatalog
        from rdkit.Chem.FilterCatalog import FilterCatalogParams

        params = FilterCatalogParams()
        for name in ("PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH", "ZINC"):
            params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, name))
        _FILTER_CATALOG = FilterCatalog.FilterCatalog(params)
    return _FILTER_CATALOG


def _reactive_patterns():
    global _REACTIVE_PATTERNS
    if _REACTIVE_PATTERNS is None:
        from rdkit import Chem

        compiled = {}
        for name, smarts in REACTIVE_SMARTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None:
                compiled[name] = pattern
        _REACTIVE_PATTERNS = compiled
    return _REACTIVE_PATTERNS


def _synthetic_accessibility(mol) -> Optional[float]:
    """RDKit Contrib SA score; absent in some minimal RDKit builds."""
    try:
        import os
        import sys

        from rdkit.Chem import RDConfig

        contrib = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if contrib not in sys.path:
            sys.path.append(contrib)
        import sascorer

        return round(float(sascorer.calculateScore(mol)), 3)
    except Exception:
        return None


def plausibility_checks(mol) -> Dict[str, Any]:
    """Property panel, structural alerts, reactive groups, and stereochemistry.

    Everything here is advisory.  The 8c7y crystal ligand itself violates one
    Lipinski rule, so property ranges describe a molecule rather than judge it.
    """
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

    molecular_weight = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    donors = Lipinski.NumHDonors(mol)
    acceptors = Lipinski.NumHAcceptors(mol)
    rotatable = Lipinski.NumRotatableBonds(mol)

    alerts = []
    try:
        alerts = sorted({entry.GetDescription() for entry in _filter_catalog().GetMatches(mol)})
    except Exception:
        alerts = []

    reactive = sorted(
        name for name, pattern in _reactive_patterns().items() if mol.HasSubstructMatch(pattern)
    )

    stereo_centers = 0
    unspecified_stereo = 0
    try:
        stereo = Chem.Mol(mol)
        Chem.AssignStereochemistryFrom3D(stereo)
        centers = Chem.FindMolChiralCenters(
            stereo, includeUnassigned=True, useLegacyImplementation=False
        )
        stereo_centers = len(centers)
        unspecified_stereo = sum(1 for _, label in centers if label == "?")
    except Exception:
        stereo_centers = -1
        unspecified_stereo = -1

    try:
        drug_likeness = round(QED.qed(mol), 4)
    except Exception:
        drug_likeness = None

    return {
        "molecular_weight": round(molecular_weight, 3),
        "clogp": round(logp, 3),
        "tpsa": round(tpsa, 2),
        "h_bond_donors": donors,
        "h_bond_acceptors": acceptors,
        "rotatable_bonds": rotatable,
        "rings": mol.GetRingInfo().NumRings(),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 4),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "qed": drug_likeness,
        "synthetic_accessibility": _synthetic_accessibility(mol),
        "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        "lipinski_violations": sum(
            [molecular_weight > 500, logp > 5, donors > 5, acceptors > 10]
        ),
        "veber_violations": sum([rotatable > 10, tpsa > 140]),
        "structural_alerts": alerts,
        "reactive_groups": reactive,
        "stereocenters": stereo_centers,
        "unspecified_stereocenters": unspecified_stereo,
    }


# --------------------------------------------------------------------------
# Pose plausibility against the receptor
# --------------------------------------------------------------------------


def load_protein_coordinates(path: Path):
    """Heavy-atom coordinates from a PDB file, or None when unreadable."""
    try:
        import numpy as np
    except ImportError:
        return None
    points = []
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
        if element == "H":
            continue
        try:
            points.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            continue
    return np.asarray(points, dtype=float) if points else None


def pose_checks(mol, context: QualityContext) -> Dict[str, Any]:
    """Steric agreement between the generated pose and the receptor.

    A molecule can be perfectly valid in isolation and still overlap the
    protein.  Two of the fifty-one audited molecules did, at 1.73 A and 1.88 A,
    against a 2.84 A closest contact for the crystal ligand.  Fraction-buried is
    deliberately not reported: every audited molecule scored 1.00, so it
    separates nothing.
    """
    import numpy as np

    result: Dict[str, Any] = {
        "min_protein_distance_a": None,
        "protein_hard_clashes": 0,
        "protein_close_contacts": 0,
        "pocket_center_offset_a": None,
    }
    points = [
        position
        for index, position in enumerate(_positions(mol))
        if mol.GetAtomWithIdx(index).GetAtomicNum() > 1
    ]
    if not points:
        return result
    ligand = np.asarray(points, dtype=float)

    if context.pocket_center is not None:
        offset = float(np.linalg.norm(ligand.mean(axis=0) - np.asarray(context.pocket_center)))
        result["pocket_center_offset_a"] = round(offset, 3)

    protein = context.protein_coordinates
    if protein is None or len(protein) == 0:
        return result

    # Restrict the receptor to a shell around the ligand so the pairwise
    # distance matrix stays small for large structures.
    centroid = ligand.mean(axis=0)
    reach = float(np.linalg.norm(ligand - centroid, axis=1).max()) + 6.0
    nearby = protein[np.linalg.norm(protein - centroid, axis=1) <= reach]
    if len(nearby) == 0:
        result["min_protein_distance_a"] = None
        return result

    distances = np.linalg.norm(ligand[:, None, :] - nearby[None, :, :], axis=2)
    result["min_protein_distance_a"] = round(float(distances.min()), 3)
    result["protein_hard_clashes"] = int((distances < context.thresholds["protein_hard_clash_a"]).sum())
    result["protein_close_contacts"] = int(
        (distances < context.thresholds["protein_close_contact_a"]).sum()
    )
    return result


# --------------------------------------------------------------------------
# Repair proposals -- generated, recorded, never applied
# --------------------------------------------------------------------------


def propose_repairs(mol) -> RepairProposal:
    """Search minimal local edits that would make a molecule sanitizable.

    The search is intentionally narrow: one bond-order reduction, optionally
    paired with one heteroatom charge or ``[nH]`` assignment.  Wholesale
    reconstruction from coordinates is not attempted, because a proposal the
    user cannot read is a proposal the user cannot judge.

    All equally-ranked variants are returned.  When more than one survives the
    molecule is ambiguous by construction -- a 3D SDF does not record which
    tautomer the model meant -- and the caller must not pick one.
    """
    from rdkit import Chem

    bond_edits = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            bond_edits.append((bond.GetIdx(), Chem.BondType.SINGLE, 1))
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            bond_edits.append((bond.GetIdx(), Chem.BondType.DOUBLE, 1))
            bond_edits.append((bond.GetIdx(), Chem.BondType.SINGLE, 2))

    atom_edits = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and any(
            bond.GetBondType() == Chem.BondType.AROMATIC for bond in atom.GetBonds()
        ):
            atom_edits.append((atom.GetIdx(), "explicit_h", 1))
        if atom.GetAtomicNum() in {7, 8, 15, 16}:
            atom_edits.append((atom.GetIdx(), "formal_charge", 1))
            atom_edits.append((atom.GetIdx(), "formal_charge", -1))

    solutions = []
    for bond_edit in [None, *bond_edits]:
        for atom_edit in [None, *atom_edits]:
            if bond_edit is None and atom_edit is None:
                continue
            candidate = Chem.Mol(mol)
            edits = []
            order_delta = 0
            if bond_edit is not None:
                index, new_type, delta = bond_edit
                bond = candidate.GetBondWithIdx(index)
                previous = bond.GetBondType()
                bond.SetBondType(new_type)
                bond.SetIsAromatic(False)
                edits.append(
                    f"bond {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()}: {previous} -> {new_type}"
                )
                order_delta = delta
            if atom_edit is not None:
                index, operation, value = atom_edit
                atom = candidate.GetAtomWithIdx(index)
                if operation == "explicit_h":
                    atom.SetNumExplicitHs(value)
                    atom.SetNoImplicit(True)
                    edits.append(f"atom {index}: assign [nH]")
                else:
                    atom.SetFormalCharge(atom.GetFormalCharge() + value)
                    edits.append(f"atom {index}: formal charge {value:+d}")
            repaired, _ = sanitize(candidate)
            if repaired is None:
                continue
            charge = sum(a.GetFormalCharge() for a in repaired.GetAtoms())
            radicals = sum(a.GetNumRadicalElectrons() for a in repaired.GetAtoms())
            rank = (
                int(bond_edit is not None) + int(atom_edit is not None),
                abs(charge),
                radicals,
                order_delta,
            )
            solutions.append((rank, repaired, edits))

    if not solutions:
        return RepairProposal(
            confidence="none",
            reason="No single-edit correction restores a valid valence model. "
            "The molecular graph itself is likely wrong; regenerate instead of repairing.",
        )

    solutions.sort(key=lambda item: item[0])
    best_rank = solutions[0][0]
    unique: Dict[str, Dict[str, str]] = {}
    for rank, repaired, edits in solutions:
        if rank != best_rank:
            continue
        smiles = Chem.MolToSmiles(repaired, canonical=True)
        unique.setdefault(
            smiles,
            {
                "smiles": smiles,
                "edits": "; ".join(edits),
                "molblock": Chem.MolToMolBlock(repaired),
            },
        )
    variants = list(unique.values())
    if len(variants) == 1:
        reason = (
            "Exactly one minimal edit restores a valid valence model. It remains a "
            "proposal: sanitizing is not evidence that this is the bond order the "
            "model intended. Review before use."
        )
        confidence = "unique_minimal"
    else:
        reason = (
            f"{len(variants)} equally minimal edits restore a valid valence model. "
            "A 3D SDF does not record which tautomer or bond assignment was intended, "
            "so PRISM cannot choose. Pick one manually or discard the candidate."
        )
        confidence = "ambiguous"
    return RepairProposal(confidence=confidence, reason=reason, variants=variants)


# --------------------------------------------------------------------------
# Safe transformations
# --------------------------------------------------------------------------


def normalize_record(record: str) -> Optional[str]:
    """Rewrite a record through RDKit so its SD data fields are well formed.

    This preserves the molecular graph, the coordinates, and every property; it
    only repairs the file layout.
    """
    from rdkit import Chem

    mol, _ = parse_unsanitized(record)
    if mol is None:
        return None
    sanitized, status = sanitize(mol)
    if sanitized is None:
        return None
    for key, value in _record_properties(record).items():
        sanitized.SetProp(key, value)
    return Chem.MolToMolBlock(sanitized) + _property_block(_record_properties(record)) + "$$$$\n"


def _record_properties(record: str) -> Dict[str, str]:
    properties: Dict[str, str] = {}
    for match in re.finditer(
        r"^>\s*<([^>]+)>[^\n]*\r?\n([^\r\n]*)", record, flags=re.MULTILINE
    ):
        properties[match.group(1).strip()] = match.group(2).strip()
    return properties


def _property_block(properties: Dict[str, str]) -> str:
    return "".join(f">  <{key}>\n{value}\n\n" for key, value in properties.items())


def hydrogenate_record(record: str) -> Tuple[Optional[str], Dict[str, Any]]:
    """Add explicit hydrogens without moving a single heavy atom.

    ``AddHs(addCoords=True)`` builds hydrogen positions geometrically from the
    existing heavy-atom frame, so the generated pose is preserved exactly. This
    is valence-based protonation only: it does not enumerate pH-dependent
    protomers and does not choose a tautomer.
    """
    from rdkit import Chem

    info: Dict[str, Any] = {
        "status": "unavailable",
        "hydrogens_added": None,
        "heavy_atom_max_displacement_a": None,
    }
    mol, _ = parse_unsanitized(record)
    if mol is None:
        return None, info
    sanitized, _ = sanitize(mol)
    if sanitized is None:
        info["status"] = "not_attempted_invalid_chemistry"
        return None, info
    try:
        hydrogenated = Chem.AddHs(Chem.Mol(sanitized), addCoords=True)
        Chem.SanitizeMol(hydrogenated)
    except Exception:
        info["status"] = "add_h_failed"
        return None, info

    before = _positions(sanitized)
    after = _positions(hydrogenated)
    displacement = 0.0
    if before and after:
        displacement = max(math.dist(before[i], after[i]) for i in range(len(before)))
    if not all(math.isfinite(v) for point in after for v in point):
        info["status"] = "add_h_failed"
        return None, info

    added = hydrogenated.GetNumAtoms() - sanitized.GetNumAtoms()
    info.update(
        {
            "status": "explicit_complete" if added == 0 else "added",
            "hydrogens_added": added,
            "heavy_atom_max_displacement_a": round(displacement, 6),
            "note": "RDKit valence-based; pH-dependent protomers and tautomers not enumerated",
        }
    )
    for key, value in _record_properties(record).items():
        hydrogenated.SetProp(key, value)
    block = Chem.MolToMolBlock(hydrogenated) + _property_block(_record_properties(record)) + "$$$$\n"
    return block, info


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def assess(record: str, context: Optional[QualityContext] = None) -> QualityReport:
    """Run every enabled check on one SDF record and classify the result."""
    context = context or QualityContext()
    if context.level == "off" or not rdkit_available():
        return QualityReport(status="NOT_CHECKED")

    from rdkit import Chem

    report = QualityReport()
    well_formed = sdf_data_fields_well_formed(record)
    mol, parse_mode = parse_unsanitized(record)
    report.parse_mode = parse_mode
    if mol is None:
        report.status = "FAIL_PARSE"
        report.flags = ["molfile_unparsable"]
        report.sanitize_status = "not_attempted"
        return report

    report.integrity = integrity_checks(mol, record)
    sanitized, sanitize_status = sanitize(mol)
    report.sanitize_status = sanitize_status

    flags: List[str] = []
    statuses: List[str] = ["PASS"]

    if not well_formed:
        flags.append("sdf_data_fields_not_terminated")
        statuses.append("WARN_SDF_FORMAT")

    if sanitized is None:
        # Speculative repair is proposed, never applied.
        report.repair = propose_repairs(mol)
        flags.append(f"sanitize_failed:{sanitize_status}")
        report.status = "FAIL_CHEMISTRY"
        report.flags = flags
        return report

    if not _structurally_complete(report.integrity):
        for key, message in (
            ("counts_match", "counts_line_mismatch"),
            ("finite_coordinates", "non_finite_coordinates"),
            ("is_3d", "not_three_dimensional"),
        ):
            if not report.integrity[key]:
                flags.append(message)
        if report.integrity["fragments"] != 1:
            flags.append(f"disconnected_fragments:{report.integrity['fragments']}")
        if report.integrity["dummy_atoms"]:
            flags.append("dummy_atoms_present")
        if report.integrity["duplicate_coordinate_pairs"]:
            flags.append("overlapping_atoms")
        report.status = "FAIL_STRUCTURE"
        report.flags = flags
        return report

    report.canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(sanitized), canonical=True)

    if context.wants_geometry:
        report.geometry = geometry_checks(sanitized, context.thresholds)
        if report.geometry["bad_bond_lengths"]:
            flags.append(f"bond_length_outliers:{report.geometry['bad_bond_lengths']}")
            statuses.append("WARN_GEOMETRY")
        if report.geometry["internal_clashes"]:
            flags.append(f"internal_clashes:{report.geometry['internal_clashes']}")
            statuses.append("WARN_GEOMETRY")
        worst_angle = report.geometry["worst_angle_deviation_deg"]
        if worst_angle is not None and worst_angle > context.thresholds["angle_warn_deg"]:
            flags.append(f"bond_angle_distortion:{worst_angle}deg")
            statuses.append("WARN_GEOMETRY")
        planarity = report.geometry["aromatic_planarity_rms_a"]
        if planarity is not None and planarity > context.thresholds["aromatic_planarity_a"]:
            flags.append(f"aromatic_ring_nonplanar:{planarity}A")
            statuses.append("WARN_GEOMETRY")

    if context.wants_plausibility:
        report.plausibility = plausibility_checks(sanitized)
        if report.plausibility["reactive_groups"]:
            flags.append("reactive_groups:" + ",".join(report.plausibility["reactive_groups"]))
            statuses.append("WARN_PLAUSIBILITY")
        if report.plausibility["unspecified_stereocenters"] > 0:
            flags.append(
                f"unspecified_stereocenters:{report.plausibility['unspecified_stereocenters']}"
            )
            statuses.append("WARN_PLAUSIBILITY")

    if context.wants_pose:
        report.pose = pose_checks(sanitized, context)
        if report.pose["protein_hard_clashes"]:
            flags.append(f"protein_clash:{report.pose['protein_hard_clashes']}")
            statuses.append("WARN_POSE")
        offset = report.pose["pocket_center_offset_a"]
        if offset is not None and offset > context.thresholds["pocket_offset_warn_a"]:
            flags.append(f"pocket_center_offset:{offset}A")
            statuses.append("WARN_POSE")

    _, hydrogen_info = hydrogenate_record(record)
    report.hydrogen = hydrogen_info

    report.status = worst_status(statuses)
    report.flags = flags
    return report
