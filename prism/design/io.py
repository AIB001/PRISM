"""Strict JSON and immutable file bindings shared by compiler and verifier."""

import hashlib
import json
import os
import re
import tempfile
from pathlib import Path

from .errors import DesignError, require


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON key: " + key, "DUPLICATE_JSON_KEY")
        result[key] = value
    return result


def _nonfinite(value):
    raise DesignError("Non-finite JSON number: " + value, "NONFINITE_JSON")


def parse_json(text):
    """Parse a JSON value from a model response with the same rules as files."""
    try:
        data = json.loads(text, object_pairs_hook=_pairs, parse_constant=_nonfinite)
    except json.JSONDecodeError as exc:
        raise DesignError("Invalid JSON response", "INVALID_JSON") from exc
    try:
        json.dumps(data, allow_nan=False)  # Also catches overflow such as 1e999.
    except ValueError as exc:
        raise DesignError("Non-finite JSON number", "NONFINITE_JSON") from exc
    return data


def read_json(path):
    try:
        return json.loads(
            Path(path).read_text(encoding="utf-8-sig"), object_pairs_hook=_pairs, parse_constant=_nonfinite
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise DesignError(f"Cannot read JSON {path}: {exc}", "INVALID_JSON") from exc


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    fd, name = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(payload)
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def resolve_path(root, value):
    if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", value):
        value = "/mnt/" + value[0].lower() + "/" + value[3:].replace("\\", "/")
    path = Path(value).expanduser()
    return (path if path.is_absolute() else Path(root) / path).resolve()


def verify_ref(root, ref):
    path = resolve_path(root, ref["path"])
    require(path.is_file(), f"Missing bound file: {path}", "MISSING_INPUT")
    require(sha256(path) == ref["sha256"], f"SHA-256 mismatch: {path}", "INPUT_HASH_MISMATCH")
    return path


def file_ref(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": sha256(path)}
