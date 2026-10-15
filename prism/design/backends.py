"""Pinned, isolated FLOWR.root execution; generic generation is not an edit API."""

import os
from pathlib import Path
import subprocess
import time

from .errors import DesignError, require
from .io import file_ref, read_json, resolve_path, sha256, write_json

SUPPORTED_COMMIT = "739a33c2aab74d074dd24605b9a2e34dc809ae9c"
DEFAULTS = {"steps": 100, "batch_cost": 1, "grow_size": 2, "device": "cuda", "timeout_seconds": 1800, "m1": None}


def load_backend(path):
    path = Path(path).expanduser().resolve()
    data = read_json(path)
    required = {"backend", "python", "source", "source_commit", "checkpoint", "checkpoint_sha256"}
    require(isinstance(data, dict) and required <= data.keys(), "Incomplete backend config", "BACKEND_CONFIG")
    require(not (data.keys() - required - DEFAULTS.keys()), "Unknown backend config options", "BACKEND_CONFIG")
    data = dict(DEFAULTS, **data)
    require(
        data["backend"] == "flowr_root", "Only the verified flowr_root edit adapter is supported", "BACKEND_UNSUPPORTED"
    )
    require(
        data["source_commit"] == SUPPORTED_COMMIT,
        "This adapter requires its tested FLOWR.root source revision",
        "BACKEND_REVISION",
    )
    for key in ("steps", "batch_cost", "grow_size", "timeout_seconds"):
        require(type(data[key]) is int and data[key] > 0, f"{key} must be a positive integer", "BACKEND_CONFIG")
    require(data["device"] in {"cuda", "cpu"}, "device must be cuda or cpu", "BACKEND_CONFIG")
    require(
        isinstance(data["checkpoint_sha256"], str) and len(data["checkpoint_sha256"]) == 64,
        "Checkpoint SHA-256 is required",
        "BACKEND_CONFIG",
    )
    for key in ("source", "checkpoint"):
        data[key] = str(resolve_path(path.parent, data[key]))
    # Resolving the venv Python symlink would execute the base interpreter and
    # lose its environment. Keep the entry path, while resolving data files.
    python = Path(data["python"]).expanduser()
    data["python"] = str((python if python.is_absolute() else path.parent / python).absolute())
    if data["m1"] is not None:
        m1 = data["m1"]
        require(
            isinstance(m1, dict)
            and set(m1) == {"checkpoint", "checkpoint_sha256", "instruction_file", "instruction_sha256"},
            "M1 requires a bound adapter checkpoint and instruction file",
            "BACKEND_CONFIG",
        )
        for key in ("checkpoint", "instruction_file"):
            m1[key] = str(resolve_path(path.parent, m1[key]))
    data["config_file"] = file_ref(path)
    return data


def check_backend(config):
    for role in ("python", "checkpoint"):
        require(Path(config[role]).is_file(), f"Missing {role}: {config[role]}", "BACKEND_MISSING")
    source = Path(config["source"])
    require(source.is_dir(), "FLOWR.root source checkout is missing", "BACKEND_MISSING")
    try:
        commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
        dirty = subprocess.check_output(
            ["git", "-C", str(source), "status", "--porcelain", "--untracked-files=no"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError) as exc:
        raise DesignError("Cannot verify FLOWR.root checkout", "BACKEND_REVISION") from exc
    require(
        commit == config["source_commit"] == SUPPORTED_COMMIT and not dirty,
        "FLOWR.root revision differs or tracked source has modifications",
        "BACKEND_REVISION",
    )
    require(
        sha256(config["checkpoint"]) == config["checkpoint_sha256"],
        "Checkpoint hash mismatch",
        "CHECKPOINT_HASH_MISMATCH",
    )
    if config.get("m1"):
        for path_key, digest_key in (("checkpoint", "checkpoint_sha256"), ("instruction_file", "instruction_sha256")):
            require(
                sha256(config["m1"][path_key]) == config["m1"][digest_key],
                "M1 input hash mismatch",
                "INPUT_HASH_MISMATCH",
            )
    return {
        "backend": "flowr_root",
        "status": "VERIFIED_FILES",
        "source_commit": commit,
        "source": str(source),
        "source_tracked_clean": True,
        "checkpoint": {"path": config["checkpoint"], "sha256": config["checkpoint_sha256"]},
        "config_file": config["config_file"],
        "runtime_bridge": file_ref(Path(__file__).with_name("flowr_runtime.py")),
        "note": "Native mask, checkpoint capability and preprocessing are probed inside each actual run",
    }


def backend_command(config, compiled, output, samples, seed, probe_only=False):
    command = [
        config["python"],
        str(Path(__file__).with_name("flowr_runtime.py")),
        "--compiled",
        str(compiled),
        "--checkpoint",
        config["checkpoint"],
        "--source-root",
        config["source"],
        "--output",
        str(output),
        "--samples",
        str(samples),
        "--seed",
        str(seed),
        "--steps",
        str(config["steps"]),
        "--batch-cost",
        str(config["batch_cost"]),
        "--grow-size",
        str(config["grow_size"]),
        "--device",
        config["device"],
    ]
    if probe_only:
        command.append("--probe-only")
    if config.get("m1"):
        m1 = config["m1"]
        command += [
            "--m1-checkpoint",
            m1["checkpoint"],
            "--m1-instruction-file",
            m1["instruction_file"],
            "--m1-checkpoint-sha256",
            m1["checkpoint_sha256"],
            "--m1-instruction-sha256",
            m1["instruction_sha256"],
            "--base-sha256",
            config["checkpoint_sha256"],
            "--source-commit",
            config["source_commit"],
        ]
    return command


def execute_backend(config, compiled, output, samples, seed):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    provenance = check_backend(config)
    command = backend_command(config, compiled, output, samples, seed)
    write_json(output / "launch.json", {"command": command, "provenance": provenance})
    environment = os.environ.copy()
    environment.pop("PYTHONHOME", None)
    environment["PYTHONPATH"] = config["source"]
    environment["PYTHONUNBUFFERED"] = "1"
    start = time.monotonic()
    status = "COMPLETED"
    with (output / "backend.log").open("w", encoding="utf-8") as log:
        try:
            process = subprocess.run(
                command,
                cwd=config["source"],
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=config["timeout_seconds"],
                check=False,
            )
            returncode = process.returncode
            if returncode:
                status = "BACKEND_FAILED"
        except subprocess.TimeoutExpired:
            returncode, status = None, "BACKEND_TIMEOUT"
        except OSError as exc:
            log.write(str(exc) + "\n")
            returncode, status = None, "BACKEND_LAUNCH_FAILED"
    result = {
        "status": status,
        "return_code": returncode,
        "wall_seconds": time.monotonic() - start,
        "command": command,
        "provenance": provenance,
    }
    write_json(output / "process.json", result)
    return result
