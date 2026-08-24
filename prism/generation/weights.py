"""Discovery, verification and download of the generative models' checkpoints.

The six models orchestrated by :mod:`prism.generation` need roughly a gigabyte
of published checkpoints. Those bytes are deliberately **not** part of a PRISM
install: the package ships the orchestration code plus the manifest that
describes every artifact, and the checkpoints are fetched on demand into a
cache directory the user controls.

``prism/data/model_weights/manifest.json`` is the single source of truth for
what an artifact is called, where it belongs on disk, how large it is and what
its SHA256 must be. It also records the upstream licence, because PRISM may
only mirror artifacts whose licence permits redistribution -- everything else
is marked ``redistributable: false`` and has to be fetched from upstream by the
user, with the manifest telling them exactly where it goes.

The on-disk layout intentionally mirrors both the upstream release layout and
the cluster deployment layout (``<models_dir>/<model>/...``), so an existing
``prism-models`` directory works as-is by pointing ``PRISM_MODELS_DIR`` at it.
"""

import hashlib
import json
import os
import shutil
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .errors import ConfigurationError, GenerationError


#: Environment variable pointing at the directory that holds the checkpoints.
MODELS_DIR_ENV = "PRISM_MODELS_DIR"

#: Environment variable overriding the manifest's mirror base URL.
MIRROR_ENV = "PRISM_MODELS_MIRROR"

_MANIFEST_NAME = "manifest.json"
_MANIFEST_SUBDIR = os.path.join("data", "model_weights")

_BLOCK_SIZE = 4 * 1024 * 1024

#: Schemes a mirror may use. ``file`` is allowed on purpose: the integration
#: cluster has no outbound internet, so an offline mirror is a local directory.
_ALLOWED_SCHEMES = {"https", "file"}

#: Hosts for which plain HTTP is tolerated (local test servers only).
_LOCAL_HOSTS = {"localhost", "127.0.0.1", "::1", "[::1]"}


class WeightsNotAvailableError(ConfigurationError):
    """A required checkpoint is absent, truncated or corrupt on disk."""

    code = "MODEL_WEIGHTS_MISSING"


class WeightsDownloadError(GenerationError):
    """A mirror download could not be completed."""

    code = "MODEL_WEIGHTS_DOWNLOAD_FAILED"


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ArtifactSpec:
    """One file a model needs, as described by the manifest."""

    model: str
    name: str
    relpath: str
    sha256: str
    size_bytes: int
    mirror_asset: Optional[str]

    @property
    def path(self) -> Path:
        return models_dir() / self.relpath


@dataclass(frozen=True)
class ModelSpec:
    """Every artifact of one model, plus the provenance that gates mirroring."""

    model: str
    title: str
    weights_license: str
    redistributable: bool
    upstream_url: str
    upstream_instructions: str
    source_repo: str
    source_commit: Optional[str]
    attribution: str
    artifacts: Dict[str, ArtifactSpec]

    @property
    def size_bytes(self) -> int:
        return sum(artifact.size_bytes for artifact in self.artifacts.values())


def manifest_path() -> Path:
    """Return the location of the bundled weights manifest."""
    return Path(__file__).resolve().parent.parent / _MANIFEST_SUBDIR / _MANIFEST_NAME


def models_dir() -> Path:
    """Return the directory holding downloaded checkpoints.

    ``PRISM_MODELS_DIR`` wins so a site can share one read-only copy between
    users; otherwise the XDG cache location is used. The directory is not
    created here -- only a download needs it to exist.
    """
    override = os.environ.get(MODELS_DIR_ENV)
    if override and override.strip():
        return Path(override.strip()).expanduser()
    cache_home = os.environ.get("XDG_CACHE_HOME")
    if cache_home and cache_home.strip():
        return Path(cache_home.strip()).expanduser() / "prism" / "models"
    return Path.home() / ".cache" / "prism" / "models"


def _relpath(model: str, name: str, value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"Weights manifest: {model}.{name} requires a non-empty 'relpath'")
    candidate = value.strip()
    if candidate.startswith("/") or candidate.startswith("\\") or ":" in candidate.split("/", 1)[0]:
        raise ConfigurationError(f"Weights manifest: {model}.{name} 'relpath' must be relative: {candidate}")
    parts = candidate.split("/")
    if any(part in {"", ".", ".."} for part in parts):
        raise ConfigurationError(f"Weights manifest: {model}.{name} 'relpath' must be normalized: {candidate}")
    return candidate


def _artifact_spec(model: str, name: str, raw: Any, redistributable: bool) -> ArtifactSpec:
    if not isinstance(raw, dict):
        raise ConfigurationError(f"Weights manifest: {model}.{name} must be a JSON object")
    digest = raw.get("sha256")
    if not isinstance(digest, str) or len(digest) != 64 or not all(c in "0123456789abcdefABCDEF" for c in digest):
        raise ConfigurationError(f"Weights manifest: {model}.{name} requires a 64-character hex 'sha256'")
    size = raw.get("size_bytes")
    if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        raise ConfigurationError(f"Weights manifest: {model}.{name} requires a positive integer 'size_bytes'")
    asset = raw.get("mirror_asset")
    if asset is not None and (not isinstance(asset, str) or not asset.strip() or "/" in asset):
        raise ConfigurationError(f"Weights manifest: {model}.{name} 'mirror_asset' must be a bare file name")
    if asset is not None and not redistributable:
        # The licence gate is a manifest invariant, not a runtime policy: a
        # non-redistributable artifact must never acquire a mirror asset by an
        # editing accident.
        raise ConfigurationError(
            f"Weights manifest: {model}.{name} declares a mirror asset but {model} is not redistributable"
        )
    return ArtifactSpec(
        model=model,
        name=name,
        relpath=_relpath(model, name, raw.get("relpath")),
        sha256=digest.lower(),
        size_bytes=size,
        mirror_asset=asset,
    )


def load_manifest(path: Optional[Path] = None) -> Dict[str, Any]:
    """Read and validate the weights manifest."""
    location = manifest_path() if path is None else Path(path).expanduser()
    if not location.is_file():
        raise WeightsNotAvailableError(
            "No model weights manifest found at {0}.\n"
            "The manifest ships inside the PRISM package; a missing file usually means "
            "PRISM was installed as a wheel built without package data. Reinstall from "
            "source with 'pip install -e .'.".format(location)
        )
    try:
        raw = json.loads(location.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"Unable to read weights manifest '{location}': {exc}") from exc
    if not isinstance(raw, dict) or not isinstance(raw.get("models"), dict):
        raise ConfigurationError(f"Weights manifest '{location}' must contain a 'models' object")

    models: Dict[str, ModelSpec] = {}
    for model, entry in raw["models"].items():
        if not isinstance(entry, dict):
            raise ConfigurationError(f"Weights manifest: '{model}' must be a JSON object")
        artifacts_raw = entry.get("artifacts")
        if not isinstance(artifacts_raw, dict) or not artifacts_raw:
            raise ConfigurationError(f"Weights manifest: '{model}' requires a non-empty 'artifacts' object")
        redistributable = entry.get("redistributable")
        if not isinstance(redistributable, bool):
            raise ConfigurationError(f"Weights manifest: '{model}' requires a boolean 'redistributable'")
        source = entry.get("source") or {}
        models[model] = ModelSpec(
            model=model,
            title=str(entry.get("title") or model),
            weights_license=str(entry.get("weights_license") or "unspecified"),
            redistributable=redistributable,
            upstream_url=str(entry.get("upstream_url") or ""),
            upstream_instructions=str(entry.get("upstream_instructions") or ""),
            source_repo=str(source.get("repo") or ""),
            source_commit=source.get("commit") or None,
            attribution=str(entry.get("attribution") or ""),
            artifacts={
                name: _artifact_spec(model, name, value, redistributable)
                for name, value in artifacts_raw.items()
            },
        )
    mirror = raw.get("mirror") or {}
    return {"models": models, "mirror_base_url": mirror.get("base_url"), "raw": raw, "path": location}


def model_specs(models: Optional[Iterable[str]] = None, path: Optional[Path] = None) -> List[ModelSpec]:
    """Return the manifest entries for ``models`` (all of them by default)."""
    manifest = load_manifest(path)
    known: Dict[str, ModelSpec] = manifest["models"]
    if models is None:
        return list(known.values())
    selected = []
    for model in models:
        spec = known.get(model)
        if spec is None:
            raise ConfigurationError(
                f"No weights manifest entry for model '{model}'. Known: " + ", ".join(sorted(known))
            )
        selected.append(spec)
    return selected


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(_BLOCK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_status(artifact: ArtifactSpec, verify: bool = False) -> str:
    """Return ``present``, ``missing``, ``size_mismatch`` or ``corrupt``.

    The size check is free and catches the common failure -- an interrupted
    copy. Hashing 600 MB is not free, so it only happens when asked for.
    """
    path = artifact.path
    if not path.is_file():
        return "missing"
    try:
        if path.stat().st_size != artifact.size_bytes:
            return "size_mismatch"
    except OSError:
        return "missing"
    if not verify:
        return "present"
    try:
        return "present" if sha256_of(path) == artifact.sha256 else "corrupt"
    except OSError:
        return "missing"


def artifact_report(
    models: Optional[Iterable[str]] = None, verify: bool = False, path: Optional[Path] = None
) -> List[Tuple[ArtifactSpec, str]]:
    """Return ``(artifact, status)`` for every artifact of the given models."""
    report = []
    for spec in model_specs(models, path):
        for artifact in spec.artifacts.values():
            report.append((artifact, artifact_status(artifact, verify=verify)))
    return report


def missing_artifacts(
    models: Optional[Iterable[str]] = None, verify: bool = False, path: Optional[Path] = None
) -> List[Tuple[ArtifactSpec, str]]:
    """Return only the artifacts that are not usable as they stand."""
    return [(artifact, status) for artifact, status in artifact_report(models, verify, path) if status != "present"]


# ---------------------------------------------------------------------------
# User-facing reporting
# ---------------------------------------------------------------------------


def format_size(size_bytes: int) -> str:
    value = float(size_bytes)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024.0 or unit == "GiB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{value:.1f} GiB"


_STATUS_TEXT = {
    "missing": "not downloaded",
    "size_mismatch": "wrong size (incomplete download?)",
    "corrupt": "SHA256 mismatch",
}


def describe_missing(entries: Sequence[Tuple[ArtifactSpec, str]], path: Optional[Path] = None) -> str:
    """Compose the message shown when generation cannot find its weights."""
    if not entries:
        return ""
    specs = {spec.model: spec for spec in model_specs(None, path)}
    affected = list(dict.fromkeys(artifact.model for artifact, _ in entries))

    lines = ["Model weights are not available locally.", ""]
    for artifact, status in entries:
        lines.append(
            f"  {artifact.model}/{artifact.name}  {format_size(artifact.size_bytes)}  "
            f"-- {_STATUS_TEXT.get(status, status)}"
        )
        lines.append(f"      {artifact.path}")
    lines.append("")

    mirrorable = [model for model in affected if specs[model].redistributable]
    upstream_only = [model for model in affected if not specs[model].redistributable]

    if mirrorable:
        lines.append("Download the weights PRISM is licensed to mirror:")
        lines.append("  prism weights download " + " ".join(f"--model {model}" for model in mirrorable))
        lines.append("")
    for model in upstream_only:
        spec = specs[model]
        lines.append(
            f"{spec.title} weights are not redistributed by PRISM "
            f"(licence: {spec.weights_license}). Fetch them from upstream and save them at "
            f"the path above:"
        )
        if spec.upstream_url:
            lines.append(f"  {spec.upstream_url}")
        if spec.upstream_instructions:
            lines.append(f"  {spec.upstream_instructions}")
        lines.append("")

    lines.append(f"Weights directory: {models_dir()}  (override with {MODELS_DIR_ENV})")
    lines.append("Inspect the full inventory with: prism weights list")
    return "\n".join(lines)


def ensure_available(models: Iterable[str], verify: bool = False, path: Optional[Path] = None) -> None:
    """Raise :class:`WeightsNotAvailableError` unless every artifact is usable."""
    entries = missing_artifacts(models, verify, path)
    if entries:
        raise WeightsNotAvailableError(describe_missing(entries, path))


# ---------------------------------------------------------------------------
# Bridging a generation config back to the manifest
# ---------------------------------------------------------------------------


def artifact_for_path(candidate: Path, path: Optional[Path] = None) -> Optional[ArtifactSpec]:
    """Return the manifest artifact a configured path refers to, if any.

    A generation config points at absolute paths. When one of them is absent,
    the difference between "you configured the wrong path" and "you have not
    downloaded the weights yet" is exactly whether the path is the one this
    manifest would have used, so that is what is matched here.
    """
    try:
        resolved = Path(candidate).expanduser().resolve(strict=False)
        root = models_dir().expanduser().resolve(strict=False)
    except OSError:
        return None
    for spec in model_specs(None, path):
        for artifact in spec.artifacts.values():
            try:
                if (root / artifact.relpath).resolve(strict=False) == resolved:
                    return artifact
            except OSError:
                continue
    return None


def configured_gaps(model: str, model_config: Mapping[str, Any]) -> List[Tuple[str, Path]]:
    """Return ``(artifact_name, path)`` for configured artifacts that are absent."""
    artifacts = model_config.get("artifacts")
    if not isinstance(artifacts, dict):
        return []
    gaps = []
    for name, specification in artifacts.items():
        if not isinstance(specification, dict):
            continue
        raw_path = specification.get("path")
        if not isinstance(raw_path, str) or not raw_path.strip():
            continue
        resolved = Path(raw_path).expanduser()
        if not resolved.is_file():
            gaps.append((str(name), resolved))
    return gaps


def describe_configured_gaps(gaps: Mapping[str, Sequence[Tuple[str, Path]]]) -> str:
    """Compose the preflight message for configured-but-absent artifacts.

    Paths that the manifest recognises are reported with the download hint;
    anything else is reported as a plain configuration problem, because PRISM
    has no way to obtain a file the user pointed somewhere of their own.
    """
    known: List[Tuple[ArtifactSpec, str]] = []
    unknown: List[Tuple[str, str, Path]] = []
    for model, entries in gaps.items():
        for name, candidate in entries:
            artifact = artifact_for_path(candidate)
            if artifact is None:
                unknown.append((model, name, candidate))
            else:
                known.append((artifact, "missing"))

    sections = []
    if known:
        sections.append(describe_missing(known))
    if unknown:
        lines = ["The generation config points at artifacts that do not exist:", ""]
        for model, name, candidate in unknown:
            lines.append(f"  {model}/{name}")
            lines.append(f"      {candidate}")
        lines.append("")
        lines.append(
            "These paths are outside the PRISM weights directory, so 'prism weights download' "
            "cannot supply them. Correct the config, or point its artifact paths at "
            f"${{{MODELS_DIR_ENV}}} so PRISM can manage them."
        )
        sections.append("\n".join(lines))
    return "\n\n".join(sections)


def preflight(config: Mapping[str, Any], models: Iterable[str]) -> None:
    """Check every enabled model's artifacts before any of them is launched.

    Without this the first model fails after the inputs have been staged and
    the remaining five are never even examined, so a user with three missing
    checkpoints discovers them one run at a time.
    """
    model_configs = config.get("models")
    if not isinstance(model_configs, dict):
        return
    gaps: Dict[str, List[Tuple[str, Path]]] = {}
    for model in models:
        model_config = model_configs.get(model)
        if not isinstance(model_config, dict) or model_config.get("enabled") is not True:
            # A disabled backend has its own, more specific error.
            continue
        found = configured_gaps(model, model_config)
        if found:
            gaps[model] = found
    if gaps:
        raise WeightsNotAvailableError(describe_configured_gaps(gaps))


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def mirror_base_url(path: Optional[Path] = None) -> Optional[str]:
    """Return the configured mirror base URL, honouring the env override."""
    override = os.environ.get(MIRROR_ENV)
    if override and override.strip():
        return override.strip()
    base = load_manifest(path)["mirror_base_url"]
    return base.strip() if isinstance(base, str) and base.strip() else None


def _validate_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme in _ALLOWED_SCHEMES:
        return url
    if parsed.scheme == "http" and parsed.hostname in _LOCAL_HOSTS:
        return url
    raise WeightsDownloadError(
        f"Refusing to download over '{parsed.scheme or 'an unspecified scheme'}': "
        f"mirrors must use https (or file:// for an offline mirror). URL: {url}"
    )


def download_artifact(
    artifact: ArtifactSpec,
    base_url: Optional[str] = None,
    force: bool = False,
    progress: Optional[Callable[[str, int, int], None]] = None,
    path: Optional[Path] = None,
) -> Path:
    """Fetch one artifact into the weights directory and verify it.

    The download lands in a ``.part`` file and is only moved into place after
    its SHA256 matches, so an interrupted download can never masquerade as a
    valid checkpoint.
    """
    target = artifact.path
    if not force and artifact_status(artifact, verify=True) == "present":
        return target

    spec = {model.model: model for model in model_specs(None, path)}[artifact.model]
    if artifact.mirror_asset is None:
        raise WeightsNotAvailableError(
            f"{spec.title} weights are not redistributed by PRISM (licence: {spec.weights_license}).\n"
            f"Fetch them from upstream and save them as:\n  {target}\n"
            + (f"  {spec.upstream_url}\n" if spec.upstream_url else "")
            + (f"  {spec.upstream_instructions}" if spec.upstream_instructions else "")
        )

    resolved_base = base_url if base_url is not None else mirror_base_url(path)
    if not resolved_base:
        raise WeightsDownloadError(
            "No weights mirror is configured yet, so PRISM cannot download "
            f"{artifact.model}/{artifact.name} for you.\n"
            f"Set {MIRROR_ENV} to a mirror base URL, or fetch the file from upstream and save it as:\n"
            f"  {target}\n"
            + (f"  {spec.upstream_url}" if spec.upstream_url else "")
        )
    url = _validate_url(resolved_base.rstrip("/") + "/" + artifact.mirror_asset)

    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + ".part")
    digest = hashlib.sha256()
    downloaded = 0
    try:
        with urllib.request.urlopen(url) as response, partial.open("wb") as handle:  # noqa: S310 - scheme checked
            while True:
                block = response.read(_BLOCK_SIZE)
                if not block:
                    break
                handle.write(block)
                digest.update(block)
                downloaded += len(block)
                if progress is not None:
                    progress(f"{artifact.model}/{artifact.name}", downloaded, artifact.size_bytes)
    except (urllib.error.URLError, OSError) as exc:
        partial.unlink(missing_ok=True)
        raise WeightsDownloadError(f"Unable to download {artifact.model}/{artifact.name} from {url}: {exc}") from exc

    actual = digest.hexdigest()
    if actual != artifact.sha256:
        partial.unlink(missing_ok=True)
        raise WeightsDownloadError(
            f"{artifact.model}/{artifact.name} SHA256 mismatch: expected {artifact.sha256}, got {actual}. "
            "The download was discarded."
        )
    if downloaded != artifact.size_bytes:
        partial.unlink(missing_ok=True)
        raise WeightsDownloadError(
            f"{artifact.model}/{artifact.name} size mismatch: expected {artifact.size_bytes} bytes, "
            f"got {downloaded}. The download was discarded."
        )
    shutil.move(str(partial), str(target))
    return target


def download_models(
    models: Iterable[str],
    force: bool = False,
    progress: Optional[Callable[[str, int, int], None]] = None,
    path: Optional[Path] = None,
) -> List[Path]:
    """Download every artifact of the given models."""
    written = []
    for spec in model_specs(models, path):
        for artifact in spec.artifacts.values():
            written.append(download_artifact(artifact, force=force, progress=progress, path=path))
    return written
