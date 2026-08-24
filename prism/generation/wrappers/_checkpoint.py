"""Checkpoint allowlisting that actually holds on every supported torch.

Passing ``pickle_module=`` to ``torch.load`` looks like it installs a
restricted unpickler, but torch < 2.0 discards it:

    if weights_only:
        if pickle_module is not None:
            raise RuntimeError(...)
    else:
        pickle_module = pickle          # torch 1.13.1, serialization.py

``weights_only`` defaults to False, so on torch 1.13 -- which five of the six
model environments use -- the custom ``find_class`` is never called and the
allowlist silently does nothing.  ``weights_only=True`` is not an option
either: every published checkpoint here carries a non-tensor config object
(EasyDict, argparse.Namespace), so it refuses to load them.

So the allowlist is enforced *before* torch is involved, by walking the pickle
opcode stream with ``pickletools.genops``.  That reports every global the
stream references without executing a single opcode, which is a stronger
guarantee than any unpickler hook: nothing runs until the checkpoint has been
cleared.
"""

import pickletools
import zipfile
from pathlib import Path
from typing import Iterable, Set, Tuple

Global = Tuple[str, str]

# Shared by every model: tensor storage rebuilding and ordered state dicts.
BASE_ALLOWED: Set[Global] = {
    ("collections", "OrderedDict"),
    ("torch._utils", "_rebuild_tensor_v2"),
    ("torch._utils", "_rebuild_parameter"),
} | {
    ("torch", f"{kind}Storage")
    for kind in (
        "BFloat16", "Bool", "Byte", "Char", "ComplexDouble", "ComplexFloat",
        "Double", "Float", "Half", "Int", "Long", "Short",
    )
}

# Opcodes that push a string usable as a STACK_GLOBAL operand.
_STRING_OPS = {
    "SHORT_BINUNICODE", "BINUNICODE", "BINUNICODE8", "UNICODE",
    "SHORT_BINSTRING", "BINSTRING", "STRING",
}
# Opcodes that leave the operand stack untouched for our purposes.
_TRANSPARENT_OPS = {"MEMOIZE", "FRAME", "PROTO"}


class CheckpointRejected(Exception):
    """The checkpoint references a pickle global outside the allowlist."""


def referenced_globals(path: Path) -> Set[Global]:
    """Every ``(module, name)`` the checkpoint's pickle streams reference.

    Executes nothing: the opcode stream is only walked.
    """
    found: Set[Global] = set()
    streams = []
    try:
        with zipfile.ZipFile(path) as archive:
            names = [name for name in archive.namelist() if name.endswith(".pkl")]
            if not names:
                raise CheckpointRejected(
                    f"{path.name} contains no pickle stream; it is not a torch checkpoint"
                )
            streams = [archive.read(name) for name in names]
    except zipfile.BadZipFile:
        # A legacy (pre-zip) checkpoint is one bare pickle stream.
        streams = [path.read_bytes()]

    for data in streams:
        stack = []
        try:
            operations = list(pickletools.genops(data))
        except Exception as exc:  # noqa: BLE001 - malformed stream is a rejection
            raise CheckpointRejected(
                f"{path.name} has an unreadable pickle stream: {exc}"
            ) from exc
        for operation, argument, _position in operations:
            name = operation.name
            if name == "GLOBAL":
                module, _, attribute = str(argument).partition(" ")
                found.add((module, attribute))
            elif name == "STACK_GLOBAL":
                if len(stack) >= 2:
                    found.add((str(stack[-2]), str(stack[-1])))
                stack.clear()
            elif name in _STRING_OPS:
                stack.append(argument)
            elif name not in _TRANSPARENT_OPS:
                stack.clear()
    return found


def verify(path: Path, extra_allowed: Iterable[Global] = ()) -> Set[Global]:
    """Reject the checkpoint unless every referenced global is allowlisted.

    Returns the set of globals it actually uses, so a caller can log what a
    trusted checkpoint legitimately needs.
    """
    allowed = set(BASE_ALLOWED) | set(extra_allowed)
    used = referenced_globals(Path(path))
    forbidden = sorted(used - allowed)
    if forbidden:
        listed = ", ".join(f"{module}.{name}" for module, name in forbidden)
        raise CheckpointRejected(
            f"Checkpoint contains forbidden pickle global(s): {listed}. "
            "Loading it would execute that code. Verify the source before use."
        )
    return used
