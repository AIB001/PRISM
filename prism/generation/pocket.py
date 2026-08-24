"""Normalize heterogeneous pocket inputs into a single :class:`PocketSpec`."""

import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

from .errors import InputError
from .types import PocketSpec


POCKET_KINDS = {"auto", "structure", "reference_ligand", "center", "residues"}
_REFERENCE_SUFFIXES = {".sdf", ".mol", ".mol2"}
_ALLOWED_YAML_KEYS = {
    "kind",
    "center",
    "radius",
    "bbox_size",
    "residues",
    "reference_ligand",
}
_RESIDUE_PATTERN = re.compile(r"^[A-Za-z0-9]:-?\d+[A-Za-z]?$")


def _triple(value: Any, name: str, allow_scalar: bool = False) -> Optional[Tuple[float, float, float]]:
    if value is None:
        return None
    if allow_scalar and isinstance(value, (int, float)) and not isinstance(value, bool):
        values: Sequence[Any] = (value, value, value)
    elif isinstance(value, (list, tuple)) and len(value) == 3:
        values = value
    else:
        suffix = " or a positive number" if allow_scalar else ""
        raise InputError(f"'{name}' must be a list of three numbers{suffix}")
    try:
        result = tuple(float(item) for item in values)
    except (TypeError, ValueError) as exc:
        raise InputError(f"'{name}' must contain only numbers") from exc
    if not all(math.isfinite(item) for item in result):
        raise InputError(f"'{name}' must contain only finite numbers")
    if name == "bbox_size" and any(item <= 0 for item in result):
        raise InputError("'bbox_size' values must be positive")
    return result  # type: ignore[return-value]


def _coordinates_from_pdb(path: Path) -> List[Tuple[float, float, float]]:
    coordinates = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    coordinates.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError as exc:
                    raise InputError(f"Invalid coordinate record in '{path}'") from exc
    return coordinates


def _coordinates_from_mol(path: Path) -> List[Tuple[float, float, float]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 4:
        raise InputError(f"MOL/SDF file '{path}' is too short")
    try:
        atom_count = int(lines[3][:3])
    except ValueError as exc:
        raise InputError(f"Unable to read the V2000 counts line in '{path}'") from exc
    coordinates = []
    for line in lines[4 : 4 + atom_count]:
        try:
            coordinates.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
        except ValueError as exc:
            raise InputError(f"Invalid atom coordinate in '{path}'") from exc
    if len(coordinates) != atom_count:
        raise InputError(f"MOL/SDF file '{path}' contains fewer atoms than declared")
    return coordinates


def _coordinates_from_mol2(path: Path) -> List[Tuple[float, float, float]]:
    coordinates = []
    in_atoms = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>"):
                if in_atoms:
                    break
                continue
            if in_atoms and line.strip():
                fields = line.split()
                if len(fields) < 5:
                    raise InputError(f"Invalid MOL2 atom record in '{path}'")
                try:
                    coordinates.append(tuple(float(item) for item in fields[2:5]))
                except ValueError as exc:
                    raise InputError(f"Invalid MOL2 atom coordinate in '{path}'") from exc
    return coordinates  # type: ignore[return-value]


def coordinates_from_structure(path: Path) -> List[Tuple[float, float, float]]:
    """Read finite atom coordinates from supported pocket structure formats."""
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        coordinates = _coordinates_from_pdb(path)
    elif suffix in {".sdf", ".mol"}:
        coordinates = _coordinates_from_mol(path)
    elif suffix == ".mol2":
        coordinates = _coordinates_from_mol2(path)
    else:
        raise InputError(f"Coordinate extraction is not supported for '{suffix}' files")
    if not coordinates:
        raise InputError(f"No atom coordinates were found in '{path}'")
    if not all(math.isfinite(value) for coordinate in coordinates for value in coordinate):
        raise InputError(f"Non-finite atom coordinate found in '{path}'")
    return coordinates


def _geometry(coordinates: Iterable[Tuple[float, float, float]]) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    points = list(coordinates)
    center = tuple(sum(point[axis] for point in points) / len(points) for axis in range(3))
    bbox = tuple(max(point[axis] for point in points) - min(point[axis] for point in points) for axis in range(3))
    return center, bbox  # type: ignore[return-value]


def _residues(value: Any) -> Tuple[str, ...]:
    if not isinstance(value, list) or not value:
        raise InputError("'residues' must be a non-empty list")
    result = tuple(str(item).strip() for item in value)
    invalid = [item for item in result if not _RESIDUE_PATTERN.fullmatch(item)]
    if invalid:
        raise InputError("Invalid residue identifier(s): " + ", ".join(invalid))
    return result


def _load_residue_file(path: Path) -> Tuple[str, ...]:
    entries = []
    for line in path.read_text(encoding="utf-8").splitlines():
        value = line.split("#", 1)[0].strip()
        if value:
            entries.append(value)
    return _residues(entries)


def _load_yaml(path: Path) -> Dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise InputError(f"Unable to read pocket YAML '{path}': {exc}") from exc
    if not isinstance(value, dict):
        raise InputError("Pocket YAML must contain a top-level mapping")
    unknown = sorted(set(value) - _ALLOWED_YAML_KEYS)
    if unknown:
        raise InputError("Unknown pocket YAML field(s): " + ", ".join(unknown))
    return value


def normalize_pocket(
    source: Path,
    kind: str = "auto",
    reference_ligand: Optional[Path] = None,
) -> PocketSpec:
    """Validate a pocket input and return its normalized scientific meaning."""
    source = source.expanduser().resolve()
    if not source.is_file():
        raise InputError(f"Pocket file does not exist: {source}")
    if kind not in POCKET_KINDS:
        raise InputError(f"Unknown pocket kind: {kind}")

    explicit_reference = reference_ligand.expanduser().resolve() if reference_ligand else None
    if explicit_reference is not None and not explicit_reference.is_file():
        raise InputError(f"Reference ligand does not exist: {explicit_reference}")

    suffix = source.suffix.lower()
    data: Dict[str, Any] = {}
    if suffix in {".yaml", ".yml"}:
        data = _load_yaml(source)
        yaml_kind = data.get("kind")
        if yaml_kind is not None and yaml_kind not in POCKET_KINDS - {"auto"}:
            raise InputError(f"Unknown pocket kind in YAML: {yaml_kind}")
        inferred_kind = yaml_kind or (
            "reference_ligand"
            if "reference_ligand" in data
            else "center"
            if "center" in data
            else "residues"
            if "residues" in data
            else None
        )
        if inferred_kind is None:
            raise InputError("Pocket YAML must define kind, center, residues, or reference_ligand")
    elif suffix in _REFERENCE_SUFFIXES:
        inferred_kind = "reference_ligand"
    elif suffix == ".pdb":
        inferred_kind = "structure"
    elif suffix == ".txt":
        inferred_kind = "residues"
    else:
        raise InputError(
            f"Unsupported pocket extension '{suffix}'. Expected PDB, SDF, MOL, MOL2, YAML, or TXT"
        )

    effective_kind = inferred_kind if kind == "auto" else kind
    if suffix not in {".yaml", ".yml"} and kind != "auto":
        valid_pairs = {
            "structure": {".pdb"},
            "reference_ligand": _REFERENCE_SUFFIXES,
            "residues": {".txt"},
        }
        allowed = valid_pairs.get(kind)
        if allowed is not None and suffix not in allowed:
            raise InputError(f"Pocket file '{source.name}' is incompatible with kind '{kind}'")

    center = _triple(data.get("center"), "center")
    bbox_size = _triple(data.get("bbox_size"), "bbox_size", allow_scalar=True)
    radius_value = data.get("radius")
    if radius_value is not None:
        if isinstance(radius_value, bool):
            raise InputError("'radius' must be a positive number")
        try:
            radius = float(radius_value)
        except (TypeError, ValueError) as exc:
            raise InputError("'radius' must be a positive number") from exc
        if not math.isfinite(radius) or radius <= 0:
            raise InputError("'radius' must be a positive finite number")
    else:
        radius = None

    residues: Tuple[str, ...] = ()
    if "residues" in data:
        residues = _residues(data["residues"])
    elif suffix == ".txt":
        residues = _load_residue_file(source)

    yaml_reference = data.get("reference_ligand")
    if yaml_reference is not None:
        if not isinstance(yaml_reference, str) or not yaml_reference.strip():
            raise InputError("'reference_ligand' must be a non-empty path")
        yaml_reference_path = (source.parent / yaml_reference).resolve()
        if not yaml_reference_path.is_file():
            raise InputError(f"Reference ligand does not exist: {yaml_reference_path}")
    else:
        yaml_reference_path = None

    normalized_reference = explicit_reference or yaml_reference_path
    if effective_kind == "reference_ligand" and normalized_reference is None:
        normalized_reference = source

    geometry_source = normalized_reference if effective_kind == "reference_ligand" else None
    if effective_kind == "structure" and suffix == ".pdb":
        geometry_source = source
    if geometry_source is not None and (center is None or bbox_size is None):
        calculated_center, calculated_bbox = _geometry(coordinates_from_structure(geometry_source))
        center = center or calculated_center
        bbox_size = bbox_size or calculated_bbox

    if effective_kind == "center" and center is None:
        raise InputError("A center pocket requires a three-coordinate 'center'")
    if effective_kind == "residues" and not residues:
        raise InputError("A residues pocket requires at least one residue")
    if effective_kind == "reference_ligand" and normalized_reference is None:
        raise InputError("A reference_ligand pocket requires a ligand file")

    return PocketSpec(
        kind=effective_kind,
        source=source,
        center=center,
        radius=radius,
        bbox_size=bbox_size,
        residues=residues,
        reference_ligand=normalized_reference,
    )
