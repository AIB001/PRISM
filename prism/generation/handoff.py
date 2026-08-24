"""Prepare finished generation runs for MD system building.

The generation pipeline writes every accepted candidate twice: the model's own
heavy-atom output under ``models/<model>/molecules/`` and, when hydrogen
addition succeeds, a hydrogenated copy under ``models/<model>/molecules_h/``.
Only the second is usable as an MD input.

Measured over 82 gaff2 builds, handing raw model output to system building
produced a topology with **zero hydrogens in 41 of 41 cases** -- correct file
layout, wrong chemistry, and no warning at any stage.  Hydrogenated input
produced the correct hydrogen count in 35 of 41.  The failure is silent, so the
check belongs on the path between the two modules rather than in either one.

This module is deliberately one-directional: it reads a finished run and writes
MD-ready inputs plus the exact command that consumes them.  It does not import
the builder, so generation stays independent of how -- or whether -- a system is
built afterwards.
"""

import csv
import json
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .errors import InputError
from .quality import STATUS_ORDER

#: Written next to the exported ligands so a consumer can tell what was chosen.
MD_MANIFEST_NAME = "md_manifest.json"
SKIPPED_REPORT_NAME = "skipped.csv"

SKIPPED_FIELDS = ("candidate_id", "model", "audit_status", "reason", "detail")


@dataclass
class MDCandidate:
    """One candidate cleared for MD system building."""

    candidate_id: str
    model: str
    ligand_path: str
    source_path: str
    audit_status: str
    hydrogen_count: int
    quality_flags: List[str] = field(default_factory=list)
    canonical_smiles: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class SkippedCandidate:
    """A candidate withheld from the MD hand-off, with the reason preserved.

    Every exclusion is recorded.  A hand-off that quietly delivered fewer
    ligands than the run produced would read as "these are the candidates"
    rather than "these are the candidates that survived".
    """

    candidate_id: str
    model: str
    audit_status: str
    reason: str
    detail: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _load_manifest(run_dir: Path) -> Dict[str, Any]:
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file():
        raise InputError(
            f"No manifest.json in '{run_dir}'. Point --run-dir at the output "
            "directory of a finished 'prism generate' run.",
            code="NO_GENERATION_MANIFEST",
        )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise InputError(
            f"manifest.json in '{run_dir}' is not valid JSON: {exc}",
            code="INVALID_GENERATION_MANIFEST",
        ) from exc
    if not isinstance(manifest, dict) or "models" not in manifest:
        raise InputError(
            f"manifest.json in '{run_dir}' has no 'models' section; it was not "
            "written by 'prism generate'.",
            code="INVALID_GENERATION_MANIFEST",
        )
    return manifest


def _candidate_records(manifest: Dict[str, Any]) -> Iterable[Tuple[str, Dict[str, Any]]]:
    for model in sorted(manifest.get("models") or {}):
        result = manifest["models"][model] or {}
        for record in result.get("candidates") or []:
            yield model, record


def count_explicit_hydrogens(path: Path) -> Optional[int]:
    """Count hydrogens actually present in an SDF file.

    The manifest records that hydrogenation succeeded; this reads the file that
    is about to be handed over.  The two can disagree -- a run directory may be
    copied, pruned, or regenerated between the two steps -- and the whole point
    of this layer is that a missing hydrogen is invisible downstream.

    Returns ``None`` when the file cannot be parsed at all.
    """
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None

    try:
        from rdkit import Chem, rdBase
    except (ImportError, ModuleNotFoundError):
        return _count_hydrogens_from_molfile(text)

    try:
        with rdBase.BlockLogs():
            molecule = Chem.MolFromMolBlock(str(text), sanitize=False, removeHs=False)
        if molecule is None:
            return _count_hydrogens_from_molfile(text)
        return sum(1 for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 1)
    except Exception:
        return _count_hydrogens_from_molfile(text)


def _count_hydrogens_from_molfile(text: str) -> Optional[int]:
    """Count hydrogens from a V2000 atom block without RDKit.

    RDKit is optional for the hand-off: refusing to prepare MD inputs because a
    chemistry toolkit is missing would be worse than reading the fixed-width
    columns the format already guarantees.
    """
    lines = text.splitlines()
    if len(lines) < 4:
        return None
    try:
        atom_count = int(lines[3][0:3])
    except (ValueError, IndexError):
        return None
    if atom_count <= 0 or len(lines) < 4 + atom_count:
        return None
    hydrogens = 0
    for line in lines[4 : 4 + atom_count]:
        # V2000 atom line: x(10) y(10) z(10) then the element symbol.
        symbol = line[31:34].strip() if len(line) >= 34 else line[30:].split()[:1]
        if isinstance(symbol, list):
            symbol = symbol[0] if symbol else ""
        if symbol == "H":
            hydrogens += 1
    return hydrogens


def _rank(candidate: MDCandidate) -> Tuple[int, int, str]:
    """Order candidates best-first.

    ``--limit`` truncates this order, so it has to put the most trustworthy
    molecules first rather than whichever model happened to be collected first.
    QC grade predicts build success directly: PASS built successfully in 100% of
    the measured cases, a geometry warning beyond 35 deg in 73%.
    """
    status = candidate.audit_status
    severity = STATUS_ORDER.index(status) if status in STATUS_ORDER else len(STATUS_ORDER)
    return (severity, len(candidate.quality_flags), candidate.candidate_id)


def select_for_md(
    run_dir: Path,
    *,
    models: Optional[Sequence[str]] = None,
    only_pass: bool = False,
    limit: Optional[int] = None,
) -> Tuple[List[MDCandidate], List[SkippedCandidate]]:
    """Choose the candidates of a finished run that can be built into a system.

    Returns the selected candidates alongside every exclusion and its reason.
    Nothing is written and nothing is copied: each returned ``ligand_path``
    still points into the run directory.  ``export_md_inputs`` copies them out.
    """
    run_dir = Path(run_dir).expanduser().resolve()
    manifest = _load_manifest(run_dir)
    wanted = {model.strip().lower() for model in models} if models else None

    eligible: List[MDCandidate] = []
    skipped: List[SkippedCandidate] = []

    for model, record in _candidate_records(manifest):
        candidate_id = record.get("candidate_id", "")
        status = record.get("audit_status", "NOT_CHECKED")

        if wanted is not None and model.lower() not in wanted:
            skipped.append(
                SkippedCandidate(candidate_id, model, status, "MODEL_NOT_SELECTED")
            )
            continue
        if status.startswith("FAIL"):
            skipped.append(
                SkippedCandidate(
                    candidate_id,
                    model,
                    status,
                    "FAILED_QC",
                    "quality control withheld this molecule; see quarantine/",
                )
            )
            continue
        if only_pass and status != "PASS":
            skipped.append(
                SkippedCandidate(candidate_id, model, status, "NOT_PASS")
            )
            continue

        hydrogenated = record.get("hydrogenated_path")
        if not hydrogenated:
            skipped.append(
                SkippedCandidate(
                    candidate_id,
                    model,
                    status,
                    "NO_HYDROGEN_FILE",
                    "the run recorded no hydrogenated copy; regenerate without --no-hydrogens",
                )
            )
            continue
        source = Path(hydrogenated)
        if not source.is_file():
            skipped.append(
                SkippedCandidate(
                    candidate_id, model, status, "NO_HYDROGEN_FILE", f"missing: {source}"
                )
            )
            continue

        hydrogens = count_explicit_hydrogens(source)
        if hydrogens is None:
            skipped.append(
                SkippedCandidate(
                    candidate_id, model, status, "UNREADABLE_LIGAND", f"cannot parse {source}"
                )
            )
            continue
        if hydrogens == 0:
            # The measured failure mode: parameterization accepts this file and
            # emits a hydrogen-free topology without complaining.
            skipped.append(
                SkippedCandidate(
                    candidate_id,
                    model,
                    status,
                    "NO_EXPLICIT_HYDROGENS",
                    f"{source} carries no hydrogens; building it would yield a "
                    "hydrogen-free topology",
                )
            )
            continue

        eligible.append(
            MDCandidate(
                candidate_id=candidate_id or f"{model}_unnamed",
                model=model,
                ligand_path=str(source),
                source_path=str(source),
                audit_status=status,
                hydrogen_count=hydrogens,
                quality_flags=list(record.get("quality_flags") or []),
                canonical_smiles=record.get("canonical_smiles", ""),
            )
        )

    eligible.sort(key=_rank)
    if limit is not None and limit >= 0 and len(eligible) > limit:
        for candidate in eligible[limit:]:
            skipped.append(
                SkippedCandidate(
                    candidate.candidate_id,
                    candidate.model,
                    candidate.audit_status,
                    "OVER_LIMIT",
                    f"ranked beyond --limit {limit}",
                )
            )
        eligible = eligible[:limit]
    return eligible, skipped


def format_build_command(
    protein: Optional[str],
    ligand_paths: Sequence[str],
    output_dir: str = "MD",
    ligand_forcefield: str = "gaff2",
    forcefield: str = "amber14sb",
) -> str:
    """Render the ``prism`` invocation that turns the export into MD systems.

    Emitted rather than executed.  Building is long-running and has its own
    failure modes; coupling it to the hand-off would make a parameterization
    error look like a generation error.
    """
    protein_argument = protein or "<protein.pdb>"
    ligands = " ".join(f"-lf {path}" for path in ligand_paths) or "-lf <ligand.sdf>"
    return (
        f"prism -pf {protein_argument} {ligands} -o {output_dir} "
        f"-lff {ligand_forcefield} -ff {forcefield}"
    )


def export_md_inputs(
    run_dir: Path,
    output_dir: Path,
    *,
    models: Optional[Sequence[str]] = None,
    only_pass: bool = False,
    limit: Optional[int] = None,
    overwrite: bool = False,
) -> Dict[str, Any]:
    """Copy MD-ready ligands out of a generation run and describe the result.

    Each exported file is the hydrogenated candidate, verified to carry
    hydrogens on disk.  The originals in the run directory are never modified.
    """
    run_dir = Path(run_dir).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    manifest = _load_manifest(run_dir)
    selected, skipped = select_for_md(
        run_dir, models=models, only_pass=only_pass, limit=limit
    )

    ligands_dir = output_dir / "ligands"
    if ligands_dir.exists() and any(ligands_dir.iterdir()):
        if not overwrite:
            raise InputError(
                f"'{ligands_dir}' already contains files. Pass --overwrite to replace them.",
                code="MD_EXPORT_EXISTS",
            )
        shutil.rmtree(ligands_dir)
    ligands_dir.mkdir(parents=True, exist_ok=True)

    exported: List[MDCandidate] = []
    for candidate in selected:
        target = ligands_dir / f"{candidate.candidate_id}.sdf"
        shutil.copyfile(candidate.source_path, target)
        exported.append(
            MDCandidate(
                candidate_id=candidate.candidate_id,
                model=candidate.model,
                ligand_path=str(target),
                source_path=candidate.source_path,
                audit_status=candidate.audit_status,
                hydrogen_count=candidate.hydrogen_count,
                quality_flags=list(candidate.quality_flags),
                canonical_smiles=candidate.canonical_smiles,
            )
        )

    with (output_dir / SKIPPED_REPORT_NAME).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SKIPPED_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(item.to_dict() for item in skipped)

    protein = (manifest.get("request") or {}).get("protein")
    ligand_paths = [candidate.ligand_path for candidate in exported]
    status_counts: Dict[str, int] = {}
    for candidate in exported:
        status_counts[candidate.audit_status] = status_counts.get(candidate.audit_status, 0) + 1
    skip_counts: Dict[str, int] = {}
    for item in skipped:
        skip_counts[item.reason] = skip_counts.get(item.reason, 0) + 1

    md_manifest = {
        "schema_version": 1,
        "generation_run": str(run_dir),
        "generation_status": manifest.get("status"),
        "protein": protein,
        "selection": {
            "models": list(models) if models else "all",
            "only_pass": only_pass,
            "limit": limit,
        },
        "exported_count": len(exported),
        "skipped_count": len(skipped),
        "status_counts": status_counts,
        "skipped_reasons": skip_counts,
        "ligands": [candidate.to_dict() for candidate in exported],
        "build_command": format_build_command(protein, ligand_paths),
        "note": (
            "Ligands are the hydrogenated copies, verified to carry explicit "
            "hydrogens. Raw model output has none and builds into a "
            "hydrogen-free topology without any error."
        ),
    }
    (output_dir / MD_MANIFEST_NAME).write_text(
        json.dumps(md_manifest, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    return md_manifest
