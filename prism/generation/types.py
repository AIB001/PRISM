"""Serializable data models used by the generation orchestration layer."""

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


@dataclass(frozen=True)
class PocketSpec:
    kind: str
    source: Path
    center: Optional[Tuple[float, float, float]] = None
    radius: Optional[float] = None
    bbox_size: Optional[Tuple[float, float, float]] = None
    residues: Tuple[str, ...] = ()
    reference_ligand: Optional[Path] = None

    def to_dict(self) -> Dict[str, Any]:
        result = asdict(self)
        result["source"] = str(self.source)
        if self.reference_ligand is not None:
            result["reference_ligand"] = str(self.reference_ligand)
        return result


@dataclass(frozen=True)
class GenerationRequest:
    protein: Path
    pocket: Path
    output: Path
    models: Tuple[str, ...]
    num_samples: int = 100
    seed: int = 0
    device: str = "auto"
    pocket_kind: str = "auto"
    reference_ligand: Optional[Path] = None
    resume: bool = False
    overwrite: bool = False
    dry_run: bool = False
    # "full" also checks the generated pose against the receptor; it degrades to
    # "standard" automatically when the protein cannot be read.
    qc_level: str = "full"
    add_hydrogens: bool = True

    def to_dict(self) -> Dict[str, Any]:
        result = asdict(self)
        for key in ("protein", "pocket", "output", "reference_ligand"):
            if result[key] is not None:
                result[key] = str(result[key])
        return result


@dataclass
class CandidateRecord:
    candidate_id: str
    model: str
    path: str
    source_path: str
    valid_3d: bool
    sanitize_status: str = "not_checked"
    native_score_name: Optional[str] = None
    native_score: Optional[float] = None
    # Quality-control fields default to a neutral value so that status.json
    # files written by earlier PRISM versions still resume.
    audit_status: str = "NOT_CHECKED"
    quality_flags: List[str] = field(default_factory=list)
    canonical_smiles: str = ""
    hydrogenated_path: Optional[str] = None
    quality: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class QuarantineRecord:
    """A record withheld from the candidate stream, preserved for review."""

    candidate_id: str
    model: str
    path: str
    source_path: str
    audit_status: str
    reason: str
    quality_flags: List[str] = field(default_factory=list)
    repair_confidence: str = "none"
    repair_variant_count: int = 0
    repair_proposal_path: Optional[str] = None
    quality: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ModelRunResult:
    model: str
    status: str
    requested: int
    generated: int = 0
    candidates: List[CandidateRecord] = field(default_factory=list)
    error_code: Optional[str] = None
    error_message: Optional[str] = None
    command: List[str] = field(default_factory=list)
    return_code: Optional[int] = None
    provenance: Dict[str, Any] = field(default_factory=dict)
    quarantined: List[QuarantineRecord] = field(default_factory=list)
    # Models that perceive bonds from geometry discard their own failures before
    # writing any SDF.  Reporting the gap between what was requested, what the
    # model wrote, and what survived QC is the only way that loss is visible.
    attrition: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
