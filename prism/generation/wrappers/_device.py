"""Runtime device selection shared by every model wrapper.

``auto`` must mean "use a GPU if this machine has one", not "assume cuda:0".
The wrappers run on whatever node the scheduler gave them, and a hardcoded
device id fails outright on a CPU-only host -- which is also what the project's
compute-resource rule forbids.

Detection deliberately prefers ``nvidia-smi`` over ``torch.cuda.is_available()``:
it is the same probe the generated run scripts use, it costs milliseconds
instead of a multi-second CUDA context import, and it works identically in
every model environment regardless of which torch build is installed.
"""

import argparse
import os
import shutil
import subprocess
from typing import Optional


def visible_gpu_count() -> int:
    """Number of usable NVIDIA GPUs, honoring CUDA_VISIBLE_DEVICES."""
    masked = os.environ.get("CUDA_VISIBLE_DEVICES")
    if masked is not None:
        entries = [item.strip() for item in masked.split(",") if item.strip()]
        if not entries or entries == ["-1"]:
            return 0
    if shutil.which("nvidia-smi") is None:
        return 0
    try:
        completed = subprocess.run(
            ["nvidia-smi", "-L"],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return 0
    if completed.returncode != 0:
        return 0
    detected = sum(1 for line in completed.stdout.splitlines() if line.startswith("GPU "))
    if masked is not None:
        # The mask renumbers devices, so the process sees at most len(entries).
        return min(detected, len([item for item in masked.split(",") if item.strip()]))
    return detected


def resolve(
    value: str,
    *,
    allow_cpu: bool = True,
    cpu_unsupported_message: Optional[str] = None,
) -> str:
    """Normalize a --device value, resolving ``auto`` against the real host.

    Args:
        value: auto, cpu, cuda, or cuda:N.
        allow_cpu: False for models whose official generator is CUDA-only.
        cpu_unsupported_message: Error text used when CPU is rejected.
    """
    normalized = value.strip().lower()
    message = cpu_unsupported_message or "this model's official generator requires CUDA"

    if normalized == "auto":
        if visible_gpu_count() > 0:
            return "cuda:0"
        if allow_cpu:
            return "cpu"
        raise argparse.ArgumentTypeError(
            f"--device auto found no usable NVIDIA GPU, and {message}"
        )
    # An explicit cuda request is the caller's decision and is passed through
    # unchecked: some CUDA containers ship without nvidia-smi, and refusing a
    # deliberate request on a failed probe would be worse than letting the model
    # raise its own error. Only ``auto`` needs the host to be interrogated.
    if normalized == "cuda":
        return "cuda:0"
    if normalized.startswith("cuda:") and normalized[5:].isdigit():
        return normalized
    if normalized == "cpu":
        if allow_cpu:
            return "cpu"
        raise argparse.ArgumentTypeError(f"--device cpu is not supported: {message}")

    options = "auto, cpu, cuda, or cuda:N" if allow_cpu else "auto, cuda, or cuda:N"
    raise argparse.ArgumentTypeError(f"device must be {options}")
