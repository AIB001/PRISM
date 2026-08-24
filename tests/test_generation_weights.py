"""Weights manifest, discovery, preflight and download."""

import hashlib
import http.server
import json
import threading
from pathlib import Path

import pytest

from prism.generation import weights
from prism.generation.cli import weights_main
from prism.generation.errors import ConfigurationError
from prism.generation.execution import _artifacts
from prism.generation.registry import model_ids
from prism.generation.runner import load_generation_config
from prism.generation.weights import WeightsDownloadError, WeightsNotAvailableError


# ---------------------------------------------------------------------------
# The shipped manifest
# ---------------------------------------------------------------------------


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


def test_non_redistributable_models_tell_the_user_where_to_get_the_file():
    for spec in weights.model_specs():
        if not spec.redistributable:
            assert spec.upstream_url or spec.upstream_instructions, (
                f"{spec.model} cannot be mirrored, so it must document an upstream source"
            )


def test_manifest_rejects_a_mirror_asset_on_a_non_redistributable_model(tmp_path):
    raw = json.loads(weights.manifest_path().read_text(encoding="utf-8"))
    raw["models"]["flowr"]["artifacts"]["checkpoint"]["mirror_asset"] = "flowr.ckpt"
    poisoned = tmp_path / "manifest.json"
    poisoned.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="not redistributable"):
        weights.load_manifest(poisoned)


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
    # pocketxmol may be mirrored, flowr may not: the message has to say both.
    assert "prism weights download --model pocketxmol" in message
    assert "not redistributed by PRISM" in message
    assert str(weights_dir) in message


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


def _publish(mirror_root, artifact, payload):
    (mirror_root / artifact.mirror_asset).write_bytes(payload)


def test_download_verifies_and_installs(weights_dir, mirror, monkeypatch):
    root, base = mirror
    artifact = _artifact("pocketxmol", "train_config")
    payload = artifact.path.name.encode() * 3
    digest = hashlib.sha256(payload).hexdigest()
    monkeypatch.setattr(
        weights,
        "model_specs",
        _patched_specs(artifact, sha256=digest, size_bytes=len(payload)),
    )
    artifact = weights.model_specs(["pocketxmol"])[0].artifacts["train_config"]
    _publish(root, artifact, payload)
    written = weights.download_artifact(artifact, base_url=base)
    assert written.read_bytes() == payload
    assert weights.artifact_status(artifact, verify=True) == "present"


def test_a_corrupt_download_is_discarded(weights_dir, mirror):
    root, base = mirror
    artifact = _artifact("pocketxmol", "train_config")
    _publish(root, artifact, b"y" * artifact.size_bytes)
    with pytest.raises(WeightsDownloadError, match="SHA256 mismatch"):
        weights.download_artifact(artifact, base_url=base)
    assert not artifact.path.exists(), "a failed download must not be left in place"
    assert not artifact.path.with_name(artifact.path.name + ".part").exists()


def test_download_refuses_a_non_redistributable_artifact(weights_dir, mirror):
    _, base = mirror
    artifact = _artifact("flowr", "checkpoint")
    with pytest.raises(WeightsNotAvailableError) as excinfo:
        weights.download_artifact(artifact, base_url=base)
    message = str(excinfo.value)
    assert "not redistributed by PRISM" in message
    assert "zenodo.org/records/15737419" in message


def test_download_refuses_an_insecure_mirror(weights_dir):
    artifact = _artifact("pocketxmol", "train_config")
    with pytest.raises(WeightsDownloadError, match="Refusing to download"):
        weights.download_artifact(artifact, base_url="http://example.com/weights")


def test_download_explains_itself_when_no_mirror_is_configured(weights_dir, monkeypatch):
    monkeypatch.delenv(weights.MIRROR_ENV, raising=False)
    artifact = _artifact("pocketxmol", "train_config")
    with pytest.raises(WeightsDownloadError) as excinfo:
        weights.download_artifact(artifact)
    message = str(excinfo.value)
    assert "No weights mirror is configured" in message
    assert str(artifact.path) in message


def test_mirror_env_overrides_the_manifest(weights_dir, monkeypatch):
    monkeypatch.setenv(weights.MIRROR_ENV, "https://example.invalid/prism-weights")
    assert weights.mirror_base_url() == "https://example.invalid/prism-weights"


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


def test_weights_download_of_an_upstream_only_model_explains_itself(weights_dir, capsys):
    assert weights_main(["download", "--model", "flowr"]) == 1
    assert "not redistributed by PRISM" in capsys.readouterr().out


def test_weights_rejects_an_unknown_model(weights_dir):
    with pytest.raises(SystemExit):
        weights_main(["list", "--model", "not-a-model"])
