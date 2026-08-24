"""Weights manifest, discovery, preflight and download."""

import hashlib
import http.server
import io
import json
import tarfile
import threading
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

from prism.generation import weights
from prism.generation.cli import weights_main
from prism.generation.errors import ConfigurationError
from prism.generation.execution import _artifacts
from prism.generation.registry import model_ids
from prism.generation.runner import load_generation_config
from prism.generation.weights import WeightsDownloadError, WeightsNotAvailableError


_LOOPBACK = {"127.0.0.1", "localhost", "::1"}


@pytest.fixture(autouse=True)
def no_internet(monkeypatch):
    """Fail loudly rather than fetch from a real publisher.

    A fixture that stripped one acquisition route but left ``download_url``
    behind once turned this suite into a 600 MiB download from Zenodo and a
    2.5-hour run that still reported the failure it was written to catch.
    Every transfer here belongs to the local mirror fixture.
    """
    real_urlopen = urllib.request.urlopen

    def guarded(request, *args, **kwargs):
        url = getattr(request, "full_url", request)
        host = urllib.parse.urlsplit(url).hostname
        if host not in _LOOPBACK:
            raise AssertionError(
                f"the test suite tried to download from {host!r}; point the "
                "manifest under test at the local mirror fixture instead"
            )
        return real_urlopen(request, *args, **kwargs)

    monkeypatch.setattr(urllib.request, "urlopen", guarded)


def test_the_suite_cannot_reach_the_internet():
    """A guard that silently passes is worse than no guard at all."""
    with pytest.raises(AssertionError, match="tried to download"):
        urllib.request.urlopen("https://zenodo.org/records/8183747")


# ---------------------------------------------------------------------------
# The shipped manifest
# ---------------------------------------------------------------------------


def _manifest_with(tmp_path, mutate):
    """Write a copy of the shipped manifest with one edit applied.

    Licence facts change; the rules about them must not. Tests that need a
    non-redistributable model build one here instead of pinning whichever model
    happens to be unmirrored today.
    """
    raw = json.loads(weights.manifest_path().read_text(encoding="utf-8"))
    mutate(raw)
    target = tmp_path / "manifest.json"
    target.write_text(json.dumps(raw), encoding="utf-8")
    return target


def _make_upstream_only(raw, model):
    """Strip every acquisition route from a model, as a licence change would.

    All three routes have to go: a model that keeps its ``download_url`` is
    still fetchable, and the tests below are about what PRISM says when it can
    obtain nothing at all.
    """
    entry = raw["models"][model]
    entry["redistributable"] = False
    entry.pop("archives", None)
    for artifact in entry["artifacts"].values():
        artifact["mirror_asset"] = None
        artifact.pop("download_url", None)
        artifact.pop("archive", None)
        artifact.pop("archive_member", None)


def test_manifest_ships_with_the_package():
    assert weights.manifest_path().is_file(), (
        "prism/data/model_weights/manifest.json must ship with the package; "
        "without it no model can report where its checkpoint belongs"
    )


def test_manifest_covers_every_registered_model():
    described = {spec.model for spec in weights.model_specs()}
    assert described == set(model_ids())


def test_every_artifact_is_fully_specified():
    for spec in weights.model_specs():
        assert spec.artifacts, f"{spec.model} declares no artifacts"
        for artifact in spec.artifacts.values():
            assert len(artifact.sha256) == 64
            assert artifact.sha256 == artifact.sha256.lower()
            assert artifact.size_bytes > 0
            assert not artifact.relpath.startswith("/")
            assert ".." not in artifact.relpath.split("/")
            assert artifact.relpath.startswith(f"{spec.model}/"), (
                "an artifact must live under its own model directory so that "
                "PRISM_MODELS_DIR can be pointed at an existing prism-models tree"
            )


def test_only_redistributable_models_carry_a_mirror_asset():
    """The licence gate is a manifest invariant, not a runtime decision."""
    for spec in weights.model_specs():
        for artifact in spec.artifacts.values():
            if artifact.mirror_asset is not None:
                assert spec.redistributable, (
                    f"{spec.model}/{artifact.name} would be mirrored by PRISM although "
                    f"{spec.model} is marked non-redistributable"
                )


def test_every_model_documents_its_upstream_source():
    """Mirroring a checkpoint does not remove the obligation to say where it came from."""
    for spec in weights.model_specs():
        assert spec.upstream_url or spec.upstream_instructions, (
            f"{spec.model} must document an upstream source"
        )


def test_every_mirrored_model_carries_an_attribution():
    for spec in weights.model_specs():
        if any(artifact.mirror_asset for artifact in spec.artifacts.values()):
            assert spec.attribution, (
                f"{spec.model} is redistributed by PRISM, so it must name whom to credit"
            )


def test_every_artifact_declares_exactly_one_acquisition_route():
    """Two routes to one file leaves open which of them is authoritative."""
    for spec in weights.model_specs():
        for artifact in spec.artifacts.values():
            routes = [
                name
                for name, value in (
                    ("mirror_asset", artifact.mirror_asset),
                    ("download_url", artifact.download_url),
                    ("archive", artifact.archive),
                )
                if value
            ]
            assert len(routes) == 1, (
                f"{spec.model}/{artifact.name} declares {len(routes)} routes ({routes or 'none'}); "
                "every artifact needs exactly one"
            )


def test_the_shipped_manifest_can_fetch_every_artifact():
    """No user of the six shipped models should have to fetch a file by hand."""
    for spec in weights.model_specs():
        for artifact in spec.artifacts.values():
            assert artifact.is_fetchable, (
                f"{spec.model}/{artifact.name} has no acquisition route, so "
                "'prism weights download' would leave it to the user"
            )


def test_publisher_downloads_go_over_https():
    """A checkpoint fetched in the clear can be swapped in transit."""
    for spec in weights.model_specs():
        for archive in spec.archives.values():
            assert archive.url.startswith("https://"), f"{spec.model} bundle is not https"
        for artifact in spec.artifacts.values():
            if artifact.download_url:
                assert artifact.download_url.startswith("https://"), (
                    f"{spec.model}/{artifact.name} is not fetched over https"
                )


def test_mirroring_is_reserved_for_what_cannot_be_fetched_from_the_publisher():
    """PRISM hosts only the checkpoints whose publisher defeats scripted download.

    Every mirrored byte is one PRISM has to store, serve and keep in step with
    upstream, so the mirror is the fallback and not the default. A new
    ``mirror_asset`` on a model that already had a working publisher URL is the
    regression this guards against.
    """
    mirrored = {
        spec.model
        for spec in weights.model_specs()
        for artifact in spec.artifacts.values()
        if artifact.mirror_asset
    }
    assert mirrored == {"targetdiff", "pocket2mol", "molcraft"}, (
        "these three publish through Google Drive, which cannot be scripted; "
        "the other three serve stable Zenodo URLs and must not be mirrored"
    )


def test_the_shipped_manifest_points_at_a_mirror(monkeypatch):
    monkeypatch.delenv(weights.MIRROR_ENV, raising=False)
    base = weights.mirror_base_url()
    assert base and base.startswith("https://"), "the weights mirror must be configured and https"


def test_mirror_asset_names_are_unique():
    """Every asset lands in one flat release, so two models cannot claim one name."""
    assets = [
        artifact.mirror_asset
        for spec in weights.model_specs()
        for artifact in spec.artifacts.values()
        if artifact.mirror_asset
    ]
    assert len(assets) == len(set(assets))


def test_manifest_rejects_a_mirror_asset_on_a_non_redistributable_model(tmp_path):
    """The licence gate fires on the combination, whatever the shipped manifest says."""

    def revoke(raw):
        model = next(
            name
            for name, entry in raw["models"].items()
            if any(artifact.get("mirror_asset") for artifact in entry["artifacts"].values())
        )
        raw["models"][model]["redistributable"] = False

    with pytest.raises(ConfigurationError, match="not redistributable"):
        weights.load_manifest(_manifest_with(tmp_path, revoke))


def test_manifest_rejects_two_acquisition_routes(tmp_path):
    def double_up(raw):
        artifact = raw["models"]["flowr"]["artifacts"]["checkpoint"]
        artifact["mirror_asset"] = "flowr-flowr_noHs.ckpt"

    with pytest.raises(ConfigurationError, match="more than one acquisition route"):
        weights.load_manifest(_manifest_with(tmp_path, double_up))


def test_manifest_rejects_an_archive_that_is_not_declared(tmp_path):
    def dangle(raw):
        raw["models"]["flowr"]["artifacts"]["checkpoint"].pop("download_url")
        raw["models"]["flowr"]["artifacts"]["checkpoint"]["archive"] = "weights"
        raw["models"]["flowr"]["artifacts"]["checkpoint"]["archive_member"] = "flowr.ckpt"

    with pytest.raises(ConfigurationError, match="which flowr does not declare"):
        weights.load_manifest(_manifest_with(tmp_path, dangle))


def test_manifest_rejects_an_archive_member_with_no_archive(tmp_path):
    def orphan(raw):
        raw["models"]["flowr"]["artifacts"]["checkpoint"]["archive_member"] = "flowr.ckpt"

    with pytest.raises(ConfigurationError, match="without naming an 'archive'"):
        weights.load_manifest(_manifest_with(tmp_path, orphan))


def test_manifest_rejects_an_insecure_publisher_url(tmp_path):
    """The scheme rule applies to publisher URLs, not only to the mirror."""

    def downgrade(raw):
        raw["models"]["flowr"]["artifacts"]["checkpoint"]["download_url"] = (
            "http://zenodo.org/records/15737419/files/flowr_noHs.ckpt"
        )

    with pytest.raises(ConfigurationError, match="Refusing to download"):
        weights.load_manifest(_manifest_with(tmp_path, downgrade))


def test_manifest_rejects_an_insecure_bundle_url(tmp_path):
    def downgrade(raw):
        raw["models"]["pocketxmol"]["archives"]["weights"]["url"] = "ftp://example.invalid/bundle.tar.gz"

    with pytest.raises(ConfigurationError, match="Refusing to download"):
        weights.load_manifest(_manifest_with(tmp_path, downgrade))


@pytest.mark.parametrize("relpath", ["/abs/x.ckpt", "../escape.ckpt", "a/../../b.ckpt"])
def test_manifest_rejects_paths_that_escape_the_weights_directory(tmp_path, relpath):
    raw = json.loads(weights.manifest_path().read_text(encoding="utf-8"))
    raw["models"]["flowr"]["artifacts"]["checkpoint"]["relpath"] = relpath
    poisoned = tmp_path / "manifest.json"
    poisoned.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ConfigurationError):
        weights.load_manifest(poisoned)


# ---------------------------------------------------------------------------
# Where the weights live
# ---------------------------------------------------------------------------


def test_models_dir_prefers_the_explicit_override(tmp_path, monkeypatch):
    monkeypatch.setenv(weights.MODELS_DIR_ENV, str(tmp_path / "shared"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    assert weights.models_dir() == tmp_path / "shared"


def test_models_dir_falls_back_to_the_xdg_cache(tmp_path, monkeypatch):
    monkeypatch.delenv(weights.MODELS_DIR_ENV, raising=False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    assert weights.models_dir() == tmp_path / "cache" / "prism" / "models"


def test_models_dir_defaults_under_the_home_cache(tmp_path, monkeypatch):
    monkeypatch.delenv(weights.MODELS_DIR_ENV, raising=False)
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    assert weights.models_dir() == tmp_path / ".cache" / "prism" / "models"


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------


@pytest.fixture()
def weights_dir(tmp_path, monkeypatch):
    root = tmp_path / "models"
    root.mkdir()
    monkeypatch.setenv(weights.MODELS_DIR_ENV, str(root))
    return root


def _artifact(model="pocketxmol", name="train_config"):
    return weights.model_specs([model])[0].artifacts[name]


def _materialize(artifact, payload=None):
    """Write a file that satisfies the manifest, or a corrupted one."""
    target = artifact.path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(payload if payload is not None else b"x" * artifact.size_bytes)
    return target


def test_absent_artifact_reports_missing(weights_dir):
    assert weights.artifact_status(_artifact()) == "missing"


def test_truncated_artifact_is_caught_without_hashing(weights_dir):
    artifact = _artifact()
    _materialize(artifact, b"short")
    assert weights.artifact_status(artifact) == "size_mismatch"


def test_right_size_wrong_bytes_is_only_caught_by_verify(weights_dir):
    artifact = _artifact()
    _materialize(artifact)
    assert weights.artifact_status(artifact) == "present"
    assert weights.artifact_status(artifact, verify=True) == "corrupt"


def test_ensure_available_lists_every_missing_artifact(weights_dir):
    with pytest.raises(WeightsNotAvailableError) as excinfo:
        weights.ensure_available(["pocketxmol", "flowr"])
    message = str(excinfo.value)
    assert "pocketxmol/checkpoint" in message
    assert "flowr/checkpoint" in message
    assert "prism weights download --model pocketxmol --model flowr" in message
    assert str(weights_dir) in message


def test_ensure_available_names_the_upstream_source_when_a_model_is_not_mirrored(weights_dir, tmp_path):
    manifest = _manifest_with(tmp_path, lambda raw: _make_upstream_only(raw, "flowr"))
    with pytest.raises(WeightsNotAvailableError) as excinfo:
        weights.ensure_available(["flowr"], path=manifest)
    message = str(excinfo.value)
    assert "not redistributed by PRISM" in message
    assert "zenodo.org/records/15737419" in message
    assert "prism weights download" not in message.split("not redistributed by PRISM")[0]


def test_ensure_available_is_silent_when_everything_is_there(weights_dir):
    for artifact in weights.model_specs(["pocketxmol"])[0].artifacts.values():
        _materialize(artifact)
    weights.ensure_available(["pocketxmol"])


# ---------------------------------------------------------------------------
# Config placeholders
# ---------------------------------------------------------------------------


def test_models_dir_placeholder_is_expanded_when_a_config_is_loaded(weights_dir, tmp_path):
    config_file = tmp_path / "generation.yaml"
    config_file.write_text(
        "models:\n"
        "  flowr:\n"
        "    enabled: true\n"
        "    artifacts:\n"
        "      checkpoint:\n"
        '        path: "${PRISM_MODELS_DIR}/flowr/checkpoints/flowr_noHs.ckpt"\n',
        encoding="utf-8",
    )
    config = load_generation_config(config_file)
    path = config["models"]["flowr"]["artifacts"]["checkpoint"]["path"]
    assert path == str(weights_dir / "flowr/checkpoints/flowr_noHs.ckpt")


def test_wrappers_placeholder_resolves_to_the_installed_wrappers(tmp_path):
    config_file = tmp_path / "generation.yaml"
    config_file.write_text('models:\n  flowr:\n    command: ["${PRISM_WRAPPERS_DIR}/flowr.py"]\n', encoding="utf-8")
    command = load_generation_config(config_file)["models"]["flowr"]["command"]
    assert Path(command[0]).is_file()


def test_unrelated_environment_variables_are_not_expanded(tmp_path, monkeypatch):
    """A config string can become a command argument; it gets two names, not the environment."""
    monkeypatch.setenv("SECRET_TOKEN", "hunter2")
    config_file = tmp_path / "generation.yaml"
    config_file.write_text('models:\n  flowr:\n    command: ["${SECRET_TOKEN}"]\n', encoding="utf-8")
    assert load_generation_config(config_file)["models"]["flowr"]["command"] == ["${SECRET_TOKEN}"]


def test_the_shipped_portable_config_resolves_to_real_weight_paths(weights_dir):
    config = load_generation_config(
        Path(__file__).resolve().parent.parent / "prism" / "configs" / "generation.local.example.yaml"
    )
    manifest = {
        artifact.path
        for spec in weights.model_specs()
        for artifact in spec.artifacts.values()
    }
    configured = {
        Path(specification["path"])
        for model_config in config["models"].values()
        for specification in model_config["artifacts"].values()
    }
    assert configured == manifest, "the portable config and the manifest must agree on every path"


def test_the_shipped_portable_config_matches_the_manifest_hashes():
    config = load_generation_config(
        Path(__file__).resolve().parent.parent / "prism" / "configs" / "generation.local.example.yaml"
    )
    for spec in weights.model_specs():
        for name, artifact in spec.artifacts.items():
            assert config["models"][spec.model]["artifacts"][name]["sha256"] == artifact.sha256


# ---------------------------------------------------------------------------
# Preflight and the execution choke point
# ---------------------------------------------------------------------------


def _config_for(model, artifact):
    return {
        "models": {
            model: {
                "enabled": True,
                "backend": "conda",
                "environment": "prism-gen-test",
                "command": ["python", "-c", "pass"],
                "artifacts": {artifact.name: {"path": str(artifact.path), "sha256": artifact.sha256}},
            }
        }
    }


def test_preflight_names_the_download_command(weights_dir):
    artifact = _artifact("pocketxmol", "checkpoint")
    with pytest.raises(WeightsNotAvailableError) as excinfo:
        weights.preflight(_config_for("pocketxmol", artifact), ["pocketxmol"])
    assert "prism weights download --model pocketxmol" in str(excinfo.value)


def test_preflight_ignores_disabled_models(weights_dir):
    artifact = _artifact("pocketxmol", "checkpoint")
    config = _config_for("pocketxmol", artifact)
    config["models"]["pocketxmol"]["enabled"] = False
    weights.preflight(config, ["pocketxmol"])


def test_preflight_reports_all_models_at_once(weights_dir):
    config = {"models": {}}
    for model in ("pocketxmol", "flowr", "diffsbdd"):
        artifact = _artifact(model, "checkpoint")
        config["models"].update(_config_for(model, artifact)["models"])
    with pytest.raises(WeightsNotAvailableError) as excinfo:
        weights.preflight(config, ["pocketxmol", "flowr", "diffsbdd"])
    message = str(excinfo.value)
    for model in ("pocketxmol", "flowr", "diffsbdd"):
        assert f"{model}/checkpoint" in message


def test_a_path_outside_the_weights_directory_is_a_plain_config_error(weights_dir, tmp_path):
    config = {
        "models": {
            "flowr": {
                "enabled": True,
                "artifacts": {"checkpoint": {"path": str(tmp_path / "elsewhere.ckpt"), "sha256": "0" * 64}},
            }
        }
    }
    with pytest.raises(WeightsNotAvailableError) as excinfo:
        weights.preflight(config, ["flowr"])
    message = str(excinfo.value)
    assert "outside the PRISM weights directory" in message
    assert "prism weights download" not in message.split("outside the PRISM weights directory")[0]


def test_execution_maps_a_managed_path_to_the_download_hint(weights_dir):
    artifact = _artifact("pocketxmol", "checkpoint")
    config = _config_for("pocketxmol", artifact)["models"]["pocketxmol"]
    with pytest.raises(WeightsNotAvailableError) as excinfo:
        _artifacts(config, "conda")
    assert excinfo.value.code == "MODEL_WEIGHTS_MISSING"
    assert "prism weights download --model pocketxmol" in str(excinfo.value)


def test_execution_keeps_the_old_error_for_an_unmanaged_path(weights_dir, tmp_path):
    config = {"artifacts": {"checkpoint": {"path": str(tmp_path / "nope.pt"), "sha256": "0" * 64}}}
    with pytest.raises(ConfigurationError) as excinfo:
        _artifacts(config, "conda")
    assert excinfo.value.code == "MODEL_NOT_CONFIGURED"
    assert "does not exist" in str(excinfo.value)


def test_weights_error_is_still_a_configuration_error(weights_dir):
    """Existing handlers catch ConfigurationError; the new error must not slip past them."""
    assert issubclass(WeightsNotAvailableError, ConfigurationError)


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


class _QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):  # noqa: ARG002 - silence the test server
        pass


@pytest.fixture()
def mirror(tmp_path):
    """Serve a directory over http://127.0.0.1 as a stand-in for the real mirror."""
    root = tmp_path / "mirror"
    root.mkdir()
    handler = type("Handler", (_QuietHandler,), {"directory": str(root)})

    def factory(*args, **kwargs):
        return handler(*args, directory=str(root), **kwargs)

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), factory)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield root, f"http://127.0.0.1:{server.server_address[1]}"
    finally:
        server.shutdown()
        server.server_close()


def _mirrored(payload, sha256=None, size_bytes=None):
    """A mirror-routed artifact of a size a test can serve.

    The three mirrored checkpoints are 33-45 MiB, so transport tests describe
    their own artifact rather than moving that much data through a test server.
    It still belongs to a real model, because ``download_artifact`` looks the
    model up in the manifest to explain itself when a fetch is refused.
    """
    return weights.ArtifactSpec(
        model="targetdiff",
        name="checkpoint",
        relpath="targetdiff/checkpoints/pretrained_diffusion.pt",
        sha256=sha256 or hashlib.sha256(payload).hexdigest(),
        size_bytes=len(payload) if size_bytes is None else size_bytes,
        mirror_asset="targetdiff-pretrained_diffusion.pt",
    )


def _publish(mirror_root, artifact, payload):
    (mirror_root / artifact.mirror_asset).write_bytes(payload)


def test_download_verifies_and_installs(weights_dir, mirror):
    root, base = mirror
    payload = b"pretrained diffusion weights" * 64
    artifact = _mirrored(payload)
    _publish(root, artifact, payload)
    written = weights.download_artifact(artifact, base_url=base)
    assert written.read_bytes() == payload
    assert weights.artifact_status(artifact, verify=True) == "present"


def test_a_corrupt_download_is_discarded(weights_dir, mirror):
    root, base = mirror
    payload = b"y" * 4096
    artifact = _mirrored(payload, sha256="0" * 64)
    _publish(root, artifact, payload)
    with pytest.raises(WeightsDownloadError, match="SHA256 mismatch"):
        weights.download_artifact(artifact, base_url=base)
    assert not artifact.path.exists(), "a failed download must not be left in place"
    assert not artifact.path.with_name(artifact.path.name + ".part").exists()


def test_a_publisher_download_needs_no_mirror(weights_dir, mirror, monkeypatch):
    """The Zenodo-hosted checkpoints are fetched with the mirror unconfigured."""
    root, base = mirror
    monkeypatch.delenv(weights.MIRROR_ENV, raising=False)
    payload = b"flowr checkpoint" * 32
    (root / "flowr_noHs.ckpt").write_bytes(payload)
    artifact = weights.ArtifactSpec(
        model="flowr",
        name="checkpoint",
        relpath="flowr/checkpoints/flowr_noHs.ckpt",
        sha256=hashlib.sha256(payload).hexdigest(),
        size_bytes=len(payload),
        mirror_asset=None,
        download_url=f"{base}/flowr_noHs.ckpt",
    )
    assert weights.download_artifact(artifact).read_bytes() == payload


def test_download_refuses_a_non_redistributable_artifact(weights_dir, mirror, tmp_path):
    _, base = mirror
    manifest = _manifest_with(tmp_path, lambda raw: _make_upstream_only(raw, "flowr"))
    artifact = weights.model_specs(["flowr"], manifest)[0].artifacts["checkpoint"]
    with pytest.raises(WeightsNotAvailableError) as excinfo:
        weights.download_artifact(artifact, base_url=base, path=manifest)
    message = str(excinfo.value)
    assert "not redistributed by PRISM" in message
    assert "zenodo.org/records/15737419" in message
    assert not artifact.path.exists()


def test_download_refuses_an_insecure_mirror(weights_dir):
    artifact = _mirrored(b"unused")
    with pytest.raises(WeightsDownloadError, match="Refusing to download"):
        weights.download_artifact(artifact, base_url="http://example.com/weights")


def test_download_explains_itself_when_no_mirror_is_configured(weights_dir, tmp_path, monkeypatch):
    monkeypatch.delenv(weights.MIRROR_ENV, raising=False)
    manifest = _manifest_with(tmp_path, lambda raw: raw["mirror"].update({"base_url": None}))
    artifact = weights.model_specs(["molcraft"], manifest)[0].artifacts["checkpoint"]
    with pytest.raises(WeightsDownloadError) as excinfo:
        weights.download_artifact(artifact, path=manifest)
    message = str(excinfo.value)
    assert "No weights mirror is configured" in message
    assert str(artifact.path) in message


def test_mirror_env_overrides_the_manifest(weights_dir, monkeypatch):
    monkeypatch.setenv(weights.MIRROR_ENV, "https://example.invalid/prism-weights")
    assert weights.mirror_base_url() == "https://example.invalid/prism-weights"


# ---------------------------------------------------------------------------
# Interrupted transfers
#
# A 600 MB checkpoint on a slow link is dropped mid-transfer often enough that
# restarting from zero means never finishing. These tests hold the server end
# of that story.
# ---------------------------------------------------------------------------


@pytest.fixture()
def flaky_mirror():
    """Serve one payload, dropping the first connection halfway through."""
    state = {
        "payload": b"",
        "drop_first": True,
        "drop_every": False,
        "honour_range": True,
        "lie_about_range": False,
        "starts": [],
    }

    class Handler(http.server.BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *args):  # noqa: ARG002 - silence the test server
            pass

        def do_GET(self):  # noqa: N802 - BaseHTTPRequestHandler's spelling
            payload = state["payload"]
            header = self.headers.get("Range")
            start = int(header.split("=", 1)[1].split("-", 1)[0]) if header else 0
            if header and state["lie_about_range"]:
                # 206, but for a range nobody asked for.
                start = 0
                ranged = True
            else:
                ranged = bool(header) and state["honour_range"]
            state["starts"].append(start if ranged else 0)
            body = payload[start:] if ranged else payload
            if ranged:
                self.send_response(206)
                self.send_header(
                    "Content-Range", f"bytes {start}-{len(payload) - 1}/{len(payload)}"
                )
            else:
                self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            if state["drop_every"] or (state["drop_first"] and len(state["starts"]) == 1):
                # Announce the full length, then hang up halfway: the client
                # sees IncompleteRead, which is what a real dropped transfer
                # looks like.
                self.wfile.write(body[: len(body) // 2])
                self.close_connection = True
                return
            self.wfile.write(body)

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield state, f"http://127.0.0.1:{server.server_address[1]}"
    finally:
        server.shutdown()
        server.server_close()


def test_a_dropped_transfer_resumes_from_what_arrived(weights_dir, flaky_mirror):
    state, base = flaky_mirror
    payload = bytes(range(256)) * 256
    state["payload"] = payload
    artifact = _mirrored(payload)

    written = weights.download_artifact(artifact, base_url=base)

    assert written.read_bytes() == payload
    assert state["starts"] == [0, len(payload) // 2], (
        "the second attempt must ask for the remaining bytes, not the whole file"
    )


def test_a_server_that_ignores_range_restarts_cleanly(weights_dir, flaky_mirror):
    """A restart must replace the partial file, never be appended to it."""
    state, base = flaky_mirror
    payload = bytes(range(256)) * 256
    state["payload"] = payload
    state["honour_range"] = False
    artifact = _mirrored(payload)

    written = weights.download_artifact(artifact, base_url=base)

    assert written.read_bytes() == payload, "a spliced file would be longer and hash wrong"
    assert state["starts"] == [0, 0]


def test_a_partial_response_for_the_wrong_range_is_not_appended(weights_dir, flaky_mirror):
    """206 is not enough; it has to be the range that was asked for."""
    state, base = flaky_mirror
    payload = bytes(range(256)) * 256
    state["payload"] = payload
    state["lie_about_range"] = True
    artifact = _mirrored(payload)

    written = weights.download_artifact(artifact, base_url=base)

    assert written.read_bytes() == payload


def test_a_transfer_that_never_completes_installs_nothing(weights_dir, flaky_mirror, monkeypatch):
    """Retries are a budget, not a loop: giving up leaves no half a checkpoint."""
    state, base = flaky_mirror
    payload = b"a" * 4096
    state["payload"] = payload
    state["drop_every"] = True
    monkeypatch.setattr(weights, "_DOWNLOAD_ATTEMPTS", 2)
    artifact = _mirrored(payload)

    with pytest.raises(WeightsDownloadError, match="after 2 attempts"):
        weights.download_artifact(artifact, base_url=base)

    assert not artifact.path.exists()
    assert not artifact.path.with_name(artifact.path.name + ".part").exists()


# ---------------------------------------------------------------------------
# Publisher bundles
#
# PocketXMol publishes every trained model in one 611 MiB tarball; PRISM wants
# two files out of it. These tests cover the single pass that takes them.
# ---------------------------------------------------------------------------

_CKPT_MEMBER = "data/trained_models/pxm/checkpoints/pocketxmol.ckpt"
_CONFIG_MEMBER = "data/trained_models/pxm/train_config/train.yml"


def _write_bundle(path, members, prefix="./model_weights/"):
    """Build a tarball whose members sit under a wrapper directory.

    Publishers re-roll bundles with and without a leading ``./`` or a top-level
    directory, so the fixture uses the awkward layout on purpose -- matching an
    exact string would pass here and fail on the real file.
    """
    with tarfile.open(path, "w:gz") as tar:
        for name, payload in members.items():
            info = tarfile.TarInfo(prefix + name)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    return path.read_bytes()


def _bundle_manifest(tmp_path, served, members, **overrides):
    """A manifest whose PocketXMol bundle is the local tarball."""

    def rewrite(raw):
        entry = raw["models"]["pocketxmol"]
        entry["archives"]["weights"] = {
            "url": served["url"],
            "size_bytes": len(served["bytes"]),
            "md5": hashlib.md5(served["bytes"]).hexdigest(),
            "note": "test bundle",
            **overrides,
        }
        for name, member in (("checkpoint", _CKPT_MEMBER), ("train_config", _CONFIG_MEMBER)):
            payload = members[member]
            entry["artifacts"][name]["sha256"] = hashlib.sha256(payload).hexdigest()
            entry["artifacts"][name]["size_bytes"] = len(payload)

    return _manifest_with(tmp_path, rewrite)


@pytest.fixture()
def bundle(tmp_path, mirror):
    """Serve a two-member tarball and return the manifest that describes it."""
    root, base = mirror
    members = {_CKPT_MEMBER: b"checkpoint bytes" * 16, _CONFIG_MEMBER: b"train:\n  seed: 0\n"}
    payload = _write_bundle(root / "model_weights.tar.gz", members)
    served = {"url": f"{base}/model_weights.tar.gz", "bytes": payload}
    return members, _bundle_manifest(tmp_path, served, members), served


def test_a_bundle_yields_every_member_it_was_downloaded_for(weights_dir, bundle):
    members, manifest, _ = bundle
    spec = weights.model_specs(["pocketxmol"], manifest)[0]

    weights.download_artifact(spec.artifacts["checkpoint"], path=manifest)

    assert spec.artifacts["checkpoint"].path.read_bytes() == members[_CKPT_MEMBER]
    assert spec.artifacts["train_config"].path.read_bytes() == members[_CONFIG_MEMBER], (
        "both members come out of one download; leaving the sibling behind "
        "would cost a second 611 MiB transfer"
    )


def test_the_bundle_is_deleted_once_nothing_needs_it(weights_dir, bundle):
    _, manifest, _ = bundle
    spec = weights.model_specs(["pocketxmol"], manifest)[0]

    weights.download_artifact(spec.artifacts["checkpoint"], path=manifest)

    cached = weights.archive_cache_path(spec.archives["weights"])
    assert not cached.exists(), "611 MiB must not stay behind to keep 232 MiB"
    assert not cached.parent.exists()


def test_a_bundle_missing_its_member_is_reported_against_the_manifest(weights_dir, tmp_path, mirror):
    """A re-rolled bundle is a manifest problem, and the message says so."""
    root, base = mirror
    members = {_CKPT_MEMBER: b"checkpoint bytes" * 16, _CONFIG_MEMBER: b"train:\n  seed: 0\n"}
    payload = _write_bundle(root / "model_weights.tar.gz", {_CONFIG_MEMBER: members[_CONFIG_MEMBER]})
    manifest = _bundle_manifest(
        tmp_path, {"url": f"{base}/model_weights.tar.gz", "bytes": payload}, members
    )
    spec = weights.model_specs(["pocketxmol"], manifest)[0]

    with pytest.raises(WeightsDownloadError) as excinfo:
        weights.download_artifact(spec.artifacts["checkpoint"], path=manifest)

    message = str(excinfo.value)
    assert "does not contain" in message
    assert _CKPT_MEMBER in message
    assert "'archive_member'" in message
    assert not spec.artifacts["checkpoint"].path.exists()
    assert spec.artifacts["train_config"].path.exists(), (
        "the member that was found verified against its own SHA256, so it is "
        "kept -- discarding it would only buy another download of the bundle"
    )


def test_a_tampered_bundle_is_refused_before_anything_is_extracted(weights_dir, tmp_path, mirror):
    root, base = mirror
    members = {_CKPT_MEMBER: b"checkpoint bytes" * 16, _CONFIG_MEMBER: b"train:\n  seed: 0\n"}
    payload = _write_bundle(root / "model_weights.tar.gz", members)
    manifest = _bundle_manifest(
        tmp_path,
        {"url": f"{base}/model_weights.tar.gz", "bytes": payload},
        members,
        md5="0" * 32,
    )
    spec = weights.model_specs(["pocketxmol"], manifest)[0]

    with pytest.raises(WeightsDownloadError, match="MD5 mismatch"):
        weights.download_artifact(spec.artifacts["checkpoint"], path=manifest)
    assert not spec.artifacts["checkpoint"].path.exists()


def test_a_cached_bundle_is_not_downloaded_twice(weights_dir, bundle):
    """An interrupted extraction must not pay for the bundle again."""
    members, manifest, served = bundle
    spec = weights.model_specs(["pocketxmol"], manifest)[0]
    cached = weights.archive_cache_path(spec.archives["weights"])
    cached.parent.mkdir(parents=True, exist_ok=True)
    cached.write_bytes(served["bytes"])

    calls = []

    def _explode(*args, **kwargs):
        calls.append(args)
        raise AssertionError("the cached bundle should have been reused")

    original = weights._stream
    weights._stream = _explode
    try:
        weights.download_artifact(spec.artifacts["checkpoint"], path=manifest)
    finally:
        weights._stream = original

    assert not calls
    assert spec.artifacts["checkpoint"].path.read_bytes() == members[_CKPT_MEMBER]


def _patched_specs(artifact, sha256, size_bytes):
    """Return a model_specs replacement with one artifact's hash/size adjusted."""
    original = weights.model_specs

    def patched(models=None, path=None):
        specs = original(models, path)
        adjusted = []
        for spec in specs:
            if spec.model != artifact.model:
                adjusted.append(spec)
                continue
            artifacts = dict(spec.artifacts)
            current = artifacts[artifact.name]
            artifacts[artifact.name] = weights.ArtifactSpec(
                model=current.model,
                name=current.name,
                relpath=current.relpath,
                sha256=sha256,
                size_bytes=size_bytes,
                mirror_asset=current.mirror_asset,
            )
            adjusted.append(
                weights.ModelSpec(
                    model=spec.model,
                    title=spec.title,
                    weights_license=spec.weights_license,
                    redistributable=spec.redistributable,
                    upstream_url=spec.upstream_url,
                    upstream_instructions=spec.upstream_instructions,
                    source_repo=spec.source_repo,
                    source_commit=spec.source_commit,
                    attribution=spec.attribution,
                    artifacts=artifacts,
                )
            )
        return adjusted

    return patched


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_weights_list_reports_every_artifact(weights_dir, capsys):
    assert weights_main(["list"]) == 0
    out = capsys.readouterr().out
    for spec in weights.model_specs():
        assert spec.model in out
    assert str(weights_dir) in out
    assert "unavailable" in out


def test_weights_list_defaults_to_list(weights_dir, capsys):
    assert weights_main([]) == 0
    assert "MODEL" in capsys.readouterr().out


def test_weights_path_prints_resolved_paths(weights_dir, capsys):
    assert weights_main(["path", "--model", "molcraft"]) == 0
    out = capsys.readouterr().out.strip().splitlines()
    assert len(out) == 1
    assert out[0].endswith(str(weights_dir / "molcraft/checkpoints/molcraft.ckpt"))


def test_weights_verify_fails_when_something_is_missing(weights_dir):
    assert weights_main(["verify", "--model", "diffsbdd"]) == 1


def test_weights_verify_passes_on_a_complete_model(weights_dir, monkeypatch):
    artifact = _artifact("diffsbdd", "checkpoint")
    payload = b"z" * 64
    monkeypatch.setattr(
        weights, "model_specs", _patched_specs(artifact, hashlib.sha256(payload).hexdigest(), len(payload))
    )
    _materialize(weights.model_specs(["diffsbdd"])[0].artifacts["checkpoint"], payload)
    assert weights_main(["verify", "--model", "diffsbdd"]) == 0


def test_weights_download_of_an_upstream_only_model_explains_itself(
    weights_dir, tmp_path, monkeypatch, capsys
):
    manifest = _manifest_with(tmp_path, lambda raw: _make_upstream_only(raw, "flowr"))
    monkeypatch.setattr(weights, "manifest_path", lambda: manifest)
    assert weights_main(["download", "--model", "flowr"]) == 1
    assert "not redistributed by PRISM" in capsys.readouterr().out


def test_weights_download_fetches_and_then_reports_what_it_kept(
    weights_dir, bundle, monkeypatch, capsys
):
    members, manifest, _ = bundle
    monkeypatch.setattr(weights, "manifest_path", lambda: manifest)

    assert weights_main(["download", "--model", "pocketxmol"]) == 0
    first = capsys.readouterr().out
    assert "the publisher ships one bundle of" in first
    spec = weights.model_specs(["pocketxmol"], manifest)[0]
    assert spec.artifacts["checkpoint"].path.read_bytes() == members[_CKPT_MEMBER]

    assert weights_main(["download", "--model", "pocketxmol"]) == 0
    second = capsys.readouterr().out
    assert second.count("already present") == 2
    assert "the publisher ships one bundle of" not in second, (
        "announcing a 611 MiB download that will not happen reads as a threat "
        "to re-fetch what is already on disk"
    )


def test_weights_rejects_an_unknown_model(weights_dir):
    with pytest.raises(SystemExit):
        weights_main(["list", "--model", "not-a-model"])
