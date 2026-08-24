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

An artifact declares at most one acquisition route, and which one it gets is a
question about the publisher rather than about the licence:

``mirror_asset``
    Fetched from the PRISM mirror. Reserved for checkpoints published through
    Google Drive, where scripted download is defeated by virus-scan
    interstitials and per-file quotas.
``download_url``
    Fetched straight from the publisher, for checkpoints on hosts that serve
    stable direct URLs. Nothing is mirrored that can be had this way.
``archive`` + ``archive_member``
    The publisher ships a bundle. PRISM downloads it once, extracts only the
    members it needs, verifies each, and discards the bundle.

No route at all means the user places the file by hand; ``upstream_url`` and
``upstream_instructions`` then say how.

The on-disk layout intentionally mirrors both the upstream release layout and
the cluster deployment layout (``<models_dir>/<model>/...``), so an existing
``prism-models`` directory works as-is by pointing ``PRISM_MODELS_DIR`` at it.
"""

import hashlib
import http.client
import json
import os
import shutil
import tarfile
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .errors import ConfigurationError, GenerationError


#: Environment variable pointing at the directory that holds the checkpoints.
MODELS_DIR_ENV = "PRISM_MODELS_DIR"

#: Environment variable overriding the manifest's mirror base URL.
MIRROR_ENV = "PRISM_MODELS_MIRROR"

_MANIFEST_NAME = "manifest.json"
_MANIFEST_SUBDIR = os.path.join("data", "model_weights")

#: Block size for reading a local file, where a large block is simply faster.
_BLOCK_SIZE = 4 * 1024 * 1024

#: Block size for reading from a socket. Deliberately much smaller: ``read(n)``
#: blocks until it has ``n`` bytes, so a 4 MiB block on a 20 KB/s link means
#: three minutes between writes -- the progress bar sits still, the ``.part``
#: file stays at zero, and a connection dropped mid-block loses all of it. This
#: is roughly ten seconds' worth of a slow link.
_STREAM_BLOCK = 256 * 1024

#: How many times a transfer is attempted before giving up. Retries resume from
#: the bytes already on disk, so this is a budget for dropped connections rather
#: than for re-downloading the file five times.
_DOWNLOAD_ATTEMPTS = 5

#: Seconds a socket may go without delivering data before the attempt is
#: abandoned. This is a per-read timeout, not a deadline for the transfer: a
#: slow link keeps working, a stalled one does not hang forever. Without it a
#: connection that opens and then delivers nothing blocks indefinitely, which
#: is exactly how a download that will never finish presents itself.
_SOCKET_TIMEOUT = 60.0

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
class ArchiveSpec:
    """A publisher's bundle that some artifacts are extracted from."""

    model: str
    name: str
    url: str
    size_bytes: int
    md5: Optional[str]
    note: str


@dataclass(frozen=True)
class ArtifactSpec:
    """One file a model needs, as described by the manifest."""

    model: str
    name: str
    relpath: str
    sha256: str
    size_bytes: int
    mirror_asset: Optional[str]
    download_url: Optional[str] = None
    archive: Optional[str] = None
    archive_member: Optional[str] = None

    @property
    def path(self) -> Path:
        return models_dir() / self.relpath

    @property
    def is_fetchable(self) -> bool:
        """Whether PRISM can obtain this file without the user's help."""
        return bool(self.mirror_asset or self.download_url or self.archive)

    @property
    def source_label(self) -> str:
        """One word for where the bytes come from, for the inventory table."""
        if self.mirror_asset:
            return "PRISM mirror"
        if self.download_url:
            return "publisher"
        if self.archive:
            return "publisher bundle"
        return "manual"


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
    archives: Dict[str, ArchiveSpec] = field(default_factory=dict)

    @property
    def size_bytes(self) -> int:
        return sum(artifact.size_bytes for artifact in self.artifacts.values())

    @property
    def is_fetchable(self) -> bool:
        return all(artifact.is_fetchable for artifact in self.artifacts.values())


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


def _manifest_url(label: str, url: str) -> str:
    """Apply the download scheme policy while reading the manifest.

    Same rule as a runtime download, but a bad URL here is a broken manifest
    rather than a failed transfer, so it is reported as one.
    """
    try:
        return _validate_url(url)
    except WeightsDownloadError as exc:
        raise ConfigurationError(f"Weights manifest: {label} {exc}") from exc


def _archive_spec(model: str, name: str, raw: Any) -> ArchiveSpec:
    if not isinstance(raw, dict):
        raise ConfigurationError(f"Weights manifest: {model} archive '{name}' must be a JSON object")
    url = raw.get("url")
    if not isinstance(url, str) or not url.strip():
        raise ConfigurationError(f"Weights manifest: {model} archive '{name}' requires a 'url'")
    size = raw.get("size_bytes")
    if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        raise ConfigurationError(
            f"Weights manifest: {model} archive '{name}' requires a positive integer 'size_bytes'"
        )
    md5 = raw.get("md5")
    if md5 is not None and (not isinstance(md5, str) or len(md5) != 32):
        raise ConfigurationError(f"Weights manifest: {model} archive '{name}' 'md5' must be 32 hex characters")
    return ArchiveSpec(
        model=model,
        name=name,
        url=_manifest_url(f"{model} archive '{name}'", url.strip()),
        size_bytes=size,
        md5=md5.lower() if isinstance(md5, str) else None,
        note=str(raw.get("note") or ""),
    )


def _artifact_spec(
    model: str,
    name: str,
    raw: Any,
    redistributable: bool,
    archives: Mapping[str, ArchiveSpec],
) -> ArtifactSpec:
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

    url = raw.get("download_url")
    if url is not None and (not isinstance(url, str) or not url.strip()):
        raise ConfigurationError(f"Weights manifest: {model}.{name} 'download_url' must be a non-empty string")
    url = _manifest_url(f"{model}.{name}", url.strip()) if isinstance(url, str) else None

    archive = raw.get("archive")
    member = raw.get("archive_member")
    if archive is not None:
        if not isinstance(archive, str) or archive not in archives:
            known = ", ".join(sorted(archives)) or "none declared"
            raise ConfigurationError(
                f"Weights manifest: {model}.{name} names archive '{archive}', "
                f"which {model} does not declare (archives: {known})"
            )
        if not isinstance(member, str) or not member.strip():
            raise ConfigurationError(
                f"Weights manifest: {model}.{name} uses an archive and so requires an 'archive_member'"
            )
        member = member.strip()
    elif member is not None:
        raise ConfigurationError(
            f"Weights manifest: {model}.{name} declares an 'archive_member' without naming an 'archive'"
        )

    routes = [route for route in (asset, url, archive) if route]
    if len(routes) > 1:
        # Two routes to the same bytes is not a richer manifest, it is an
        # unanswered question about which one is authoritative.
        raise ConfigurationError(
            f"Weights manifest: {model}.{name} declares more than one acquisition route; "
            "give it exactly one of mirror_asset, download_url or archive"
        )

    return ArtifactSpec(
        model=model,
        name=name,
        relpath=_relpath(model, name, raw.get("relpath")),
        sha256=digest.lower(),
        size_bytes=size,
        mirror_asset=asset,
        download_url=url,
        archive=archive,
        archive_member=member,
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
        archives_raw = entry.get("archives") or {}
        if not isinstance(archives_raw, dict):
            raise ConfigurationError(f"Weights manifest: '{model}' 'archives' must be a JSON object")
        archives = {name: _archive_spec(model, name, value) for name, value in archives_raw.items()}
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
                name: _artifact_spec(model, name, value, redistributable, archives)
                for name, value in artifacts_raw.items()
            },
            archives=archives,
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


def _digest_of(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(_BLOCK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_of(path: Path) -> str:
    return _digest_of(path, "sha256")


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

    fetchable = [model for model in affected if specs[model].is_fetchable]
    manual = [model for model in affected if not specs[model].is_fetchable]

    if fetchable:
        lines.append("PRISM can fetch these for you:")
        lines.append("  prism weights download " + " ".join(f"--model {model}" for model in fetchable))
        lines.append("")
    for model in manual:
        spec = specs[model]
        lines.append(manual_reason(spec) + ". Download them and save them at the path above:")
        if spec.upstream_url:
            lines.append(f"  {spec.upstream_url}")
        if spec.upstream_instructions:
            lines.append(f"  {spec.upstream_instructions}")
        lines.append("")

    lines.append(f"Weights directory: {models_dir()}  (override with {MODELS_DIR_ENV})")
    lines.append("Inspect the full inventory with: prism weights list")
    return "\n".join(lines)


def manual_reason(spec: ModelSpec) -> str:
    """Say why PRISM will not fetch a model's weights itself.

    Two quite different situations look identical on disk, and the difference
    decides what the user should do next. A licence that forbids redistribution
    is permanent; an artifact that merely has no route in the manifest is a gap
    PRISM could close. Naming the licence when the licence is not the reason
    would send someone looking for permission they already have.
    """
    if not spec.redistributable:
        return f"{spec.title} weights are not redistributed by PRISM (licence: {spec.weights_license})"
    return f"{spec.title} weights have no automatic download source"


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


def _resumes_at(response: Any, offset: int) -> bool:
    """Whether ``response`` continues the file from ``offset``.

    A partial response has to be both partial *and* the part that was asked
    for: ``Content-Range`` is what says so, and a server that answers 206 with
    some other range would otherwise be appended to blindly.
    """
    if getattr(response, "status", None) != 206:
        return False
    header = response.headers.get("Content-Range", "")
    _, _, span = header.partition("bytes ")
    start, _, _ = span.partition("-")
    try:
        return int(start.strip()) == offset
    except ValueError:
        return False


def _stream(
    url: str,
    destination: Path,
    label: str,
    expected_size: int,
    progress: Optional[Callable[[str, int, int], None]],
    algorithm: str = "sha256",
) -> Tuple[str, int]:
    """Stream ``url`` into ``destination``, returning ``(hexdigest, written)``.

    Servers hosting hundred-megabyte checkpoints drop connections, and a
    checkpoint that has to restart from zero on a slow link never finishes at
    all. Each attempt therefore resumes with a ``Range`` request from whatever
    is already on disk, and only a server that ignores the range forces a
    restart.

    The digest is computed by reading the finished file rather than by hashing
    blocks as they arrive: across a resumed transfer an incremental digest and
    the bytes actually on disk can disagree, and the digest is the one thing
    here that has to be trustworthy.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.unlink(missing_ok=True)
    written = 0
    last_error = ""

    for attempt in range(1, _DOWNLOAD_ATTEMPTS + 1):
        request = urllib.request.Request(url)
        if written:
            request.add_header("Range", f"bytes={written}-")
        try:
            with urllib.request.urlopen(request, timeout=_SOCKET_TIMEOUT) as response:  # noqa: S310 - scheme checked
                if written and not _resumes_at(response, written):
                    # The server ignored the range, or answered with a
                    # different one; appending either would splice two copies
                    # together. The SHA256 would catch that at the end, but
                    # only after paying for the whole transfer.
                    written = 0
                with destination.open("ab" if written else "wb") as handle:
                    while True:
                        block = response.read(_STREAM_BLOCK)
                        if not block:
                            break
                        handle.write(block)
                        written += len(block)
                        if progress is not None:
                            progress(label, written, expected_size)
        except (urllib.error.URLError, http.client.HTTPException, OSError) as exc:
            # A server that announces a Content-Length and then closes early
            # raises IncompleteRead, which is an HTTPException and not an
            # OSError -- and that is precisely the failure this loop exists to
            # survive, so it has to be named explicitly.
            last_error = f"{type(exc).__name__}: {exc}"
        else:
            if not expected_size or written >= expected_size:
                return _digest_of(destination, algorithm), written
            last_error = f"the connection closed after {written} of {expected_size} bytes"

        # Trust the filesystem over the loop counter about how much survived.
        written = destination.stat().st_size if destination.is_file() else 0
        if attempt < _DOWNLOAD_ATTEMPTS and progress is not None:
            progress(f"{label} (retry {attempt}/{_DOWNLOAD_ATTEMPTS - 1})", written, expected_size)

    destination.unlink(missing_ok=True)
    raise WeightsDownloadError(
        f"Unable to download {label} from {url} after {_DOWNLOAD_ATTEMPTS} attempts: {last_error}"
    )


def _accept_or_discard(artifact: ArtifactSpec, partial: Path, actual: str, written: int) -> Path:
    """Move a staged file into place, or delete it and explain why not."""
    if actual != artifact.sha256:
        partial.unlink(missing_ok=True)
        raise WeightsDownloadError(
            f"{artifact.model}/{artifact.name} SHA256 mismatch: expected {artifact.sha256}, got {actual}. "
            "The download was discarded."
        )
    if written != artifact.size_bytes:
        partial.unlink(missing_ok=True)
        raise WeightsDownloadError(
            f"{artifact.model}/{artifact.name} size mismatch: expected {artifact.size_bytes} bytes, "
            f"got {written}. The download was discarded."
        )
    target = artifact.path
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(partial), str(target))
    return target


def _fetch_file(artifact: ArtifactSpec, url: str, progress) -> Path:
    """Download one artifact straight to its final path, via a ``.part`` file."""
    partial = artifact.path.with_name(artifact.path.name + ".part")
    label = f"{artifact.model}/{artifact.name}"
    actual, written = _stream(url, partial, label, artifact.size_bytes, progress)
    return _accept_or_discard(artifact, partial, actual, written)


# ---------------------------------------------------------------------------
# Publisher bundles
# ---------------------------------------------------------------------------

#: Where a bundle is cached while its members are extracted. It lives under the
#: weights directory so it inherits the same disk budget the user chose.
_ARCHIVE_SUBDIR = ".archives"


def archive_cache_path(archive: ArchiveSpec) -> Path:
    name = Path(urllib.parse.urlparse(archive.url).path).name or f"{archive.name}.tar.gz"
    return models_dir() / _ARCHIVE_SUBDIR / f"{archive.model}-{name}"


def _ensure_archive(archive: ArchiveSpec, progress) -> Path:
    """Return a local copy of the publisher's bundle, downloading if needed."""
    cached = archive_cache_path(archive)
    if cached.is_file() and cached.stat().st_size == archive.size_bytes:
        if archive.md5 is None or _digest_of(cached, "md5") == archive.md5:
            return cached
        cached.unlink(missing_ok=True)

    partial = cached.with_name(cached.name + ".part")
    label = f"{archive.model} bundle"
    algorithm = "md5" if archive.md5 else "sha256"
    actual, written = _stream(archive.url, partial, label, archive.size_bytes, progress, algorithm)
    if written != archive.size_bytes:
        partial.unlink(missing_ok=True)
        raise WeightsDownloadError(
            f"{label} size mismatch: expected {archive.size_bytes} bytes, got {written}. "
            "The download was discarded."
        )
    if archive.md5 and actual != archive.md5:
        partial.unlink(missing_ok=True)
        raise WeightsDownloadError(
            f"{label} MD5 mismatch: expected {archive.md5}, got {actual}. The download was discarded."
        )
    shutil.move(str(partial), str(cached))
    return cached


def _member_key(name: str, wanted: Mapping[str, ArtifactSpec]) -> Optional[str]:
    """Match a tar member against the members we are looking for.

    Matching on a suffix rather than the exact string is deliberate: publishers
    re-roll bundles with or without a leading ``./`` or a top-level directory,
    and the SHA256 check downstream makes a loose match safe -- the wrong file
    cannot survive it.
    """
    candidate = name.lstrip("./")
    for key in wanted:
        if candidate == key or candidate.endswith("/" + key):
            return key
    return None


def _extract_members(
    archive_path: Path, wanted: Mapping[str, ArtifactSpec], progress
) -> Dict[str, Path]:
    """Extract the wanted members in a single pass over the bundle.

    Members are streamed out by name and verified individually; the archive is
    never unpacked wholesale, so a hostile path inside it has nothing to write
    to.
    """
    remaining = dict(wanted)
    written: Dict[str, Path] = {}
    try:
        with tarfile.open(archive_path, "r:*") as tar:
            for member in tar:
                if not remaining:
                    break
                if not member.isfile():
                    continue
                key = _member_key(member.name, remaining)
                if key is None:
                    continue
                artifact = remaining.pop(key)
                source = tar.extractfile(member)
                if source is None:  # pragma: no cover - isfile() already excludes these
                    continue
                partial = artifact.path.with_name(artifact.path.name + ".part")
                partial.parent.mkdir(parents=True, exist_ok=True)
                digest = hashlib.sha256()
                count = 0
                with partial.open("wb") as handle:
                    for block in iter(lambda: source.read(_BLOCK_SIZE), b""):
                        handle.write(block)
                        digest.update(block)
                        count += len(block)
                        if progress is not None:
                            progress(f"{artifact.model}/{artifact.name}", count, artifact.size_bytes)
                written[key] = _accept_or_discard(artifact, partial, digest.hexdigest(), count)
    except (tarfile.TarError, OSError) as exc:
        raise WeightsDownloadError(f"Unable to read the bundle at {archive_path}: {exc}") from exc

    if remaining:
        missing = ", ".join(sorted(remaining))
        raise WeightsDownloadError(
            f"The bundle at {archive_path} does not contain: {missing}.\n"
            "The publisher may have re-rolled it with a different layout; the manifest's "
            "'archive_member' entries need updating."
        )
    return written


def _fetch_from_archive(artifact: ArtifactSpec, spec: ModelSpec, force: bool, progress) -> Path:
    """Extract ``artifact`` -- and any sibling still absent -- from one bundle.

    Downloading a 600 MB bundle once per artifact would be absurd when two of
    them come out of the same file, so every outstanding member is taken in the
    same pass and the bundle is dropped once nothing else needs it.
    """
    archive = spec.archives[artifact.archive]
    wanted = {
        sibling.archive_member: sibling
        for sibling in spec.artifacts.values()
        if sibling.archive == artifact.archive
        and (sibling.name == artifact.name or force or artifact_status(sibling, verify=True) != "present")
    }
    archive_path = _ensure_archive(archive, progress)
    extracted = _extract_members(archive_path, wanted, progress)

    from_archive = [a for a in spec.artifacts.values() if a.archive == artifact.archive]
    if all(artifact_status(a) == "present" for a in from_archive):
        archive_path.unlink(missing_ok=True)
        parent = archive_path.parent
        try:
            next(parent.iterdir())
        except StopIteration:
            parent.rmdir()
        except OSError:
            pass
    return extracted[artifact.archive_member]


def download_artifact(
    artifact: ArtifactSpec,
    base_url: Optional[str] = None,
    force: bool = False,
    progress: Optional[Callable[[str, int, int], None]] = None,
    path: Optional[Path] = None,
) -> Path:
    """Fetch one artifact into the weights directory and verify it.

    The bytes land in a ``.part`` file and are only moved into place once the
    SHA256 matches, so an interrupted download can never masquerade as a valid
    checkpoint. An artifact that comes out of a publisher bundle may bring its
    siblings with it, since they share one download.
    """
    target = artifact.path
    if not force and artifact_status(artifact, verify=True) == "present":
        return target

    spec = {model.model: model for model in model_specs(None, path)}[artifact.model]

    if artifact.archive is not None:
        return _fetch_from_archive(artifact, spec, force, progress)

    if artifact.download_url is not None:
        return _fetch_file(artifact, artifact.download_url, progress)

    if artifact.mirror_asset is None:
        raise WeightsNotAvailableError(
            manual_reason(spec) + ".\n"
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
    return _fetch_file(artifact, _validate_url(resolved_base.rstrip("/") + "/" + artifact.mirror_asset), progress)


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
