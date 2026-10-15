"""Validate the annotation contract without importing chemistry or GPU code."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

from .errors import DesignError, require
from .io import read_json, resolve_path, sha256

ANNOTATION_SCHEMA = "prism-design-annotation/0.1"


@dataclass(frozen=True)
class DesignDocument:
    path: Path
    root: Path
    data: Dict[str, Any]
    sha256: str


def validate_schema(data):
    try:
        from jsonschema import Draft202012Validator
    except ImportError as exc:
        raise DesignError(
            "DesignSpec validation requires jsonschema; install prism-md[design]", "MISSING_DEPENDENCY"
        ) from exc
    schema = read_json(Path(__file__).parent / "schemas" / "annotation-v0.1.json")
    errors = sorted(Draft202012Validator(schema).iter_errors(data), key=lambda e: str(list(e.path)))
    if errors:
        message = "; ".join(f"{'.'.join(map(str, e.path))}: {e.message}" for e in errors[:5])
        raise DesignError(message, "SCHEMA_VALIDATION_FAILED")


def load_design(path, input_root=None):
    path = Path(path).expanduser().resolve()
    data = read_json(path)
    validate_schema(data)
    anchor = data["provenance"]["selection_manifest"]
    roots = [Path(input_root).resolve()] if input_root is not None else list(path.parents)
    root = None
    for candidate in roots:
        target = resolve_path(candidate, anchor["path"])
        if target.is_file() and sha256(target) == anchor["sha256"]:
            root = candidate
            break
    require(
        root is not None, "Cannot locate the hash-bound annotation root; use --input-root", "ANNOTATION_ROOT_NOT_FOUND"
    )
    return DesignDocument(path, root, data, sha256(path))
