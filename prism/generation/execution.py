"""Subprocess execution isolated behind Docker, Conda, or Slurm."""

import os
import subprocess
import hashlib
import re
import signal
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Tuple

from .errors import ConfigurationError, ExecutionFailure
from .weights import WeightsNotAvailableError, artifact_for_path, describe_missing


@dataclass
class ExecutionResult:
    command: List[str]
    return_code: int
    artifact_sha256: Dict[str, str] = field(default_factory=dict)


class _StrictFormat(dict):
    def __missing__(self, key: str) -> str:
        raise ConfigurationError(f"Unknown command placeholder: {{{key}}}")


def _render_command(template: Any, context: Mapping[str, Any]) -> List[str]:
    if not isinstance(template, list) or not template:
        raise ConfigurationError("Configured model command must be a non-empty YAML list")
    rendered = []
    values = _StrictFormat({key: str(value) for key, value in context.items()})
    for item in template:
        if not isinstance(item, (str, int, float)):
            raise ConfigurationError("Every configured command item must be a string or number")
        try:
            rendered.append(str(item).format_map(values))
        except ValueError as exc:
            raise ConfigurationError(f"Invalid command template item '{item}': {exc}") from exc
    return rendered


def _container_device_args(device: str) -> List[str]:
    normalized = device.strip().lower()
    if normalized == "cpu":
        return []
    if normalized in {"auto", "cuda"}:
        return ["--gpus", "all"]
    if normalized.startswith("cuda:"):
        identifier = normalized.split(":", 1)[1]
        if not identifier.isdigit():
            raise ConfigurationError(f"Invalid CUDA device: {device}")
        return ["--gpus", f"device={identifier}"]
    raise ConfigurationError(f"Invalid device: {device}")


def _conda_command(environment: Any, inner_command: List[str]) -> List[str]:
    if not isinstance(environment, str) or not environment.strip():
        raise ConfigurationError("A Conda environment must be a non-empty name or path")
    selector = "-p" if "/" in environment else "-n"
    return ["conda", "run", "--no-capture-output", selector, environment] + inner_command


def _slurm_token(value: Any, name: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise ConfigurationError(f"Slurm '{name}' must contain only letters, numbers, '.', '_', or '-'")
    return value


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, bool):
        raise ConfigurationError(f"Slurm '{name}' must be a positive integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"Slurm '{name}' must be a positive integer") from exc
    if result <= 0:
        raise ConfigurationError(f"Slurm '{name}' must be a positive integer")
    return result


def _slurm_command(
    model_config: Mapping[str, Any], command_context: Mapping[str, Any], run_dir: Path
) -> List[str]:
    slurm = model_config.get("slurm")
    if not isinstance(slurm, dict):
        raise ConfigurationError("A Slurm backend requires a 'slurm' YAML mapping")
    unknown = sorted(
        set(slurm)
        - {
            "executable",
            "partition",
            "account",
            "qos",
            "time",
            "gpus",
            "cpus_per_task",
            "memory",
            "job_name",
        }
    )
    if unknown:
        raise ConfigurationError("Unknown Slurm configuration field(s): " + ", ".join(unknown))

    executable = slurm.get("executable", "srun")
    if not isinstance(executable, str) or not executable.strip():
        raise ConfigurationError("Slurm 'executable' must be a non-empty path or command")
    partition = _slurm_token(slurm.get("partition"), "partition")
    account = _slurm_token(slurm.get("account"), "account")
    time_limit = slurm.get("time")
    if not isinstance(time_limit, str) or not re.fullmatch(r"(?:\d+-)?\d{1,2}:\d{2}:\d{2}", time_limit):
        raise ConfigurationError("Slurm 'time' must use [D-]HH:MM:SS format")
    memory = slurm.get("memory")
    if not isinstance(memory, str) or not re.fullmatch(r"\d+(?:[KMGTP])?", memory, re.IGNORECASE):
        raise ConfigurationError("Slurm 'memory' must be a number with an optional K/M/G/T/P suffix")
    gpus = _positive_integer(slurm.get("gpus", 1), "gpus")
    cpus = _positive_integer(slurm.get("cpus_per_task", 4), "cpus_per_task")
    job_name = _slurm_token(slurm.get("job_name", "prism-generation"), "job_name")

    command = [
        executable,
        "--partition",
        partition,
        "--account",
        account,
        "--time",
        time_limit,
        "--gres",
        f"gpu:{gpus}",
        "--cpus-per-task",
        str(cpus),
        "--mem",
        memory,
        "--job-name",
        job_name,
        "--chdir",
        str(run_dir.resolve()),
        "--kill-on-bad-exit=1",
    ]
    qos = slurm.get("qos")
    if qos is not None:
        command.extend(["--qos", _slurm_token(qos, "qos")])

    slurm_context = dict(command_context)
    if "device" in slurm_context:
        slurm_context["device"] = "cuda:0"
    inner_command = _render_command(model_config.get("command"), slurm_context)
    command.extend(_conda_command(model_config.get("environment"), inner_command))
    return command


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _missing_artifact_error(name: str, source: Path) -> ConfigurationError:
    """Distinguish an unconfigured path from an undownloaded checkpoint.

    Both are "the file is not there", but only one of them has a remedy PRISM
    can offer, so a path the weights manifest recognises gets the download
    instructions instead of a bare "does not exist".
    """
    artifact = artifact_for_path(source)
    if artifact is None:
        return ConfigurationError(f"Artifact '{name}' does not exist: {source}")
    return WeightsNotAvailableError(describe_missing([(artifact, "missing")]))


def _artifacts(
    model_config: Mapping[str, Any], backend: Any
) -> Tuple[Dict[str, str], List[Tuple[Path, str]], Dict[str, str]]:
    configured = model_config.get("artifacts", {})
    if not isinstance(configured, dict):
        raise ConfigurationError("Model 'artifacts' must be a YAML mapping")
    context: Dict[str, str] = {}
    mounts: List[Tuple[Path, str]] = []
    hashes: Dict[str, str] = {}
    for name, specification in configured.items():
        if not isinstance(name, str) or not name.isidentifier():
            raise ConfigurationError(f"Invalid artifact name: {name!r}")
        if not isinstance(specification, dict):
            raise ConfigurationError(f"Artifact '{name}' must be a YAML mapping")
        source_value = specification.get("path")
        expected_hash = specification.get("sha256")
        if not isinstance(source_value, str) or not source_value.strip():
            raise ConfigurationError(f"Artifact '{name}' requires a non-empty 'path'")
        if not isinstance(expected_hash, str) or len(expected_hash) != 64:
            raise ConfigurationError(f"Artifact '{name}' requires a 64-character 'sha256'")
        source = Path(source_value).expanduser().resolve()
        if not source.is_file():
            raise _missing_artifact_error(name, source)
        actual_hash = _sha256(source)
        if actual_hash.lower() != expected_hash.lower():
            message = f"Artifact '{name}' SHA256 mismatch: expected {expected_hash.lower()}, got {actual_hash}"
            managed = artifact_for_path(source)
            if managed is not None:
                message += (
                    "\nThis file is managed by PRISM; re-fetch it with "
                    f"'prism weights download --model {managed.model} --force' or check the whole "
                    "inventory with 'prism weights verify'."
                )
            raise ConfigurationError(message)
        hashes[name] = actual_hash
        if backend == "docker":
            target = specification.get("container_path")
            if not isinstance(target, str) or not target.startswith("/"):
                raise ConfigurationError(
                    f"Docker artifact '{name}' requires an absolute 'container_path'"
                )
            context[name] = target
            mounts.append((source, target))
        else:
            context[name] = str(source)
    return context, mounts, hashes



def _terminate_process_group(process: "subprocess.Popen") -> None:
    """Signal the launcher's whole process group, then the launcher itself.

    ``start_new_session=True`` makes the launcher a process-group leader, so a
    single ``killpg`` reaches the model and anything it spawned. SIGTERM first
    gives a wrapper the chance to clean up its own children; SIGKILL follows if
    the group is still alive.
    """
    for signal_number, grace in ((signal.SIGTERM, 5.0), (signal.SIGKILL, 5.0)):
        try:
            os.killpg(os.getpgid(process.pid), signal_number)
        except (ProcessLookupError, PermissionError, OSError):
            # Already gone, or the group could not be resolved; fall back to
            # the direct child so the launcher is never left running.
            try:
                process.kill()
            except OSError:
                pass
            break
        try:
            process.wait(timeout=grace)
            break
        except subprocess.TimeoutExpired:
            continue
    try:
        process.wait(timeout=5.0)
    except subprocess.TimeoutExpired:
        pass


def execute(
    model_config: Mapping[str, Any],
    context: Mapping[str, Any],
    output_root: Path,
    run_dir: Path,
    device: str,
) -> ExecutionResult:
    """Render and execute one configured model wrapper without invoking a shell."""
    if model_config.get("enabled") is not True:
        raise ConfigurationError(
            "Model backend is disabled. Configure and enable it in a generation YAML file."
        )
    backend = model_config.get("backend")
    command_context = dict(context)
    artifact_context, artifact_mounts, artifact_hashes = _artifacts(model_config, backend)
    overlap = sorted(set(command_context) & set(artifact_context))
    if overlap:
        raise ConfigurationError("Artifact placeholder conflicts with built-in context: " + ", ".join(overlap))
    command_context.update(artifact_context)
    environment_values = model_config.get("env", {})
    if not isinstance(environment_values, dict):
        raise ConfigurationError("Model 'env' must be a YAML mapping")

    if backend == "conda":
        inner_command = _render_command(model_config.get("command"), command_context)
        command = _conda_command(model_config.get("environment"), inner_command)
    elif backend == "slurm":
        command = _slurm_command(model_config, command_context, run_dir)
    elif backend == "docker":
        image = model_config.get("image")
        if not isinstance(image, str) or not image.strip():
            raise ConfigurationError("A Docker backend requires a pinned 'image'")
        if model_config.get("allow_unpinned_image") is not True and not re.search(
            r"@sha256:[0-9a-fA-F]{64}$", image
        ):
            raise ConfigurationError(
                "Docker image must be pinned by digest (name@sha256:...) unless "
                "'allow_unpinned_image: true' is explicitly configured"
            )
        normalized_device = device.strip().lower()
        if "device" in command_context:
            command_context["device"] = "cpu" if normalized_device == "cpu" else "cuda:0"
        container_context = {}
        for key, value in command_context.items():
            value_path = Path(str(value))
            try:
                relative = value_path.resolve().relative_to(output_root.resolve())
                container_context[key] = str(Path("/work") / relative)
            except (OSError, ValueError):
                container_context[key] = value
        inner_command = _render_command(model_config.get("command"), container_context)
        command = ["docker", "run", "--rm"]
        command.extend(_container_device_args(device))
        for key, value in sorted(environment_values.items()):
            command.extend(["--env", f"{key}={value}"])
        for source, target in artifact_mounts:
            command.extend(["--volume", f"{source}:{target}:ro"])
        command.extend(
            [
                "--volume",
                f"{output_root.resolve()}:/work",
                "--workdir",
                str(Path("/work") / run_dir.resolve().relative_to(output_root.resolve())),
                image,
            ]
        )
        command.extend(inner_command)
    else:
        raise ConfigurationError("Model backend must be 'docker', 'conda', or 'slurm'")

    environment = os.environ.copy()
    environment.update({str(key): str(value) for key, value in environment_values.items()})

    logs_dir = run_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    timeout = model_config.get("timeout_seconds")
    if timeout is not None and (not isinstance(timeout, int) or timeout <= 0):
        raise ConfigurationError("'timeout_seconds' must be a positive integer")

    try:
        with (logs_dir / "stdout.log").open("w", encoding="utf-8") as stdout_handle, (
            logs_dir / "stderr.log"
        ).open("w", encoding="utf-8") as stderr_handle:
            # The direct child is only a launcher -- ``conda run``, ``srun`` or
            # ``docker run`` -- and the model itself is a grandchild. On timeout
            # ``subprocess.run`` kills the direct child only, which leaves the
            # real computation running: an abandoned AM1-BCC job was observed
            # burning a core for nearly two hours after its task was recorded as
            # timed out. Running the launcher in its own session means the whole
            # tree can be signalled at once.
            process = subprocess.Popen(
                command,
                cwd=str(run_dir),
                env=environment,
                stdout=stdout_handle,
                stderr=stderr_handle,
                start_new_session=True,
            )
            try:
                return_code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                _terminate_process_group(process)
                raise ExecutionFailure(
                    f"Model execution exceeded its {timeout}-second timeout"
                ) from None
            except BaseException:
                # A cancelled run must not leave the model behind either.
                _terminate_process_group(process)
                raise
    except FileNotFoundError as exc:
        raise ExecutionFailure(f"Execution backend was not found: {command[0]}") from exc
    except OSError as exc:
        raise ExecutionFailure(f"Unable to start model process: {exc}") from exc

    if return_code != 0:
        raise ExecutionFailure(
            f"Model process exited with code {return_code}; see {logs_dir / 'stderr.log'}",
            code="MODEL_PROCESS_FAILED",
        )
    return ExecutionResult(
        command=command,
        return_code=return_code,
        artifact_sha256=artifact_hashes,
    )
