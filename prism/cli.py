"""Top-level command router for PRISM."""

import sys
from contextlib import contextmanager
from typing import Iterator, List, Optional, Sequence, Tuple


_TRUE_VALUES = {"1", "true", "yes", "on"}
_FALSE_VALUES = {"0", "false", "no", "off"}


def _parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in _TRUE_VALUES:
        return True
    if normalized in _FALSE_VALUES:
        return False
    raise ValueError(
        "--generation expects true or false "
        "(accepted aliases: 1/0, yes/no, on/off)"
    )


def _extract_generation_flag(argv: Sequence[str]) -> Tuple[Optional[bool], List[str]]:
    """Remove and parse the legacy ``--generation`` flag."""
    result: List[str] = []
    generation: Optional[bool] = None
    index = 0
    while index < len(argv):
        argument = argv[index]
        if argument.startswith("--generation="):
            if generation is not None:
                raise ValueError("--generation may only be specified once")
            generation = _parse_bool(argument.split("=", 1)[1])
        elif argument == "--generation":
            if generation is not None:
                raise ValueError("--generation may only be specified once")
            next_value = argv[index + 1].strip().lower() if index + 1 < len(argv) else ""
            if next_value in _TRUE_VALUES | _FALSE_VALUES:
                generation = _parse_bool(argv[index + 1])
                index += 1
            else:
                generation = True
        else:
            result.append(argument)
        index += 1
    return generation, result


@contextmanager
def _forwarded_argv(argv: Sequence[str]) -> Iterator[None]:
    original = sys.argv
    sys.argv = [original[0]] + list(argv)
    try:
        yield
    finally:
        sys.argv = original


def main(argv: Optional[Sequence[str]] = None):
    """Route generation commands and preserve the existing Builder CLI."""
    arguments = list(sys.argv[1:] if argv is None else argv)

    if arguments and arguments[0] == "generate":
        from .generation.cli import main as generation_main

        return generation_main(arguments[1:])

    if arguments and arguments[0] == "prepare-md":
        from .generation.cli import prepare_md_main

        return prepare_md_main(arguments[1:])

    try:
        generation, forwarded = _extract_generation_flag(arguments)
    except ValueError as exc:
        raise SystemExit(f"prism: error: {exc}") from exc

    if generation:
        from .generation.cli import main as generation_main

        return generation_main(forwarded)

    from .builder import main as builder_main

    if generation is False or argv is not None:
        with _forwarded_argv(forwarded):
            return builder_main()
    return builder_main()


if __name__ == "__main__":
    main()
