#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PTM request specification parsing (CLI and YAML).

CLI form (repeatable)::

    --ptm A:145:SEP            # chain A, residue 145 -> phosphoserine (default charge)
    --ptm A:32:TPO:-1          # explicit protonation charge
    --ssbond auto              # auto-detect disulfides (default) | none

YAML form::

    ptm:
      disulfides: auto                 # auto | none | [[A, 12], [A, 40]]
      amber_phospho_ff: phosaa19SB     # phosaa19SB | phosaa14SB | phosaa10
      residues:
        - {chain: A, resid: 145, code: SEP}
        - {chain: A, resid: 32,  code: TPO, charge: -1}
"""

from __future__ import annotations

from dataclasses import MISSING, dataclass, field, fields
from numbers import Integral
import re
from typing import List, Optional, Tuple, Union

from .catalog import get_ptm, is_known_ptm


#: Chain id assumed when a request (or a PDB record) carries no chain.
DEFAULT_CHAIN = "A"

#: Disulfide handling: detect SG-SG pairs geometrically, or leave cysteines
#: alone.  An explicit list of ``[chain, resid]`` pairs restricts detection to
#: those positions.
DEFAULT_DISULFIDE_MODE = "auto"
DISULFIDE_MODES = (DEFAULT_DISULFIDE_MODE, "none")

#: AMBER phosphorylation parameter sets accepted for the ``tleap`` route.
#: Declared once here so the YAML parser, the CLI ``--phospho-ff`` choices and
#: :meth:`PTMConfig.validation_errors` cannot drift apart.
DEFAULT_PHOSPHO_FF = "phosaa19SB"
PHOSPHO_PARAMETER_SETS = (DEFAULT_PHOSPHO_FF, "phosaa14SB", "phosaa10")
_PHOSPHO_BY_LOWER = {name.lower(): name for name in PHOSPHO_PARAMETER_SETS}


def _normalise_integer(value, *, context: str) -> int:
    """Accept genuine integers and integer strings without truncation."""
    if isinstance(value, bool):
        raise ValueError(f"Invalid {context} {value!r}; expected an integer.")
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, str) and re.fullmatch(r"[+-]?\d+", value.strip()):
        return int(value.strip())
    raise ValueError(f"Invalid {context} {value!r}; expected an integer.")


def _normalise_charge(value, *, context: str) -> Optional[int]:
    if value is None:
        return None
    try:
        return _normalise_integer(value, context=f"charge for {context}")
    except ValueError as exc:
        raise ValueError(f"Invalid charge {value!r} for {context}; expected an integer.") from exc


def _require_known_code(code, *, context: str) -> str:
    """Reject a PTM code the catalog does not define, in both parsers alike."""
    normalized = str(code).strip().upper()
    if not is_known_ptm(normalized):
        raise ValueError(
            f"Unknown PTM code {normalized!r} in {context}. "
            "Known codes: see `prism --list-ptms`."
        )
    return normalized


def _require_phospho_ff(value) -> str:
    """Return the canonical spelling of an AMBER phospho parameter set."""
    if isinstance(value, str):
        canonical = _PHOSPHO_BY_LOWER.get(value.strip().lower())
        if canonical is not None:
            return canonical
    raise ValueError(
        f"Unknown PTM phosphorylation parameter set {value!r}; "
        f"choose one of: {', '.join(PHOSPHO_PARAMETER_SETS)}."
    )


def _require_disulfide_mode(value) -> Union[str, List[Tuple[str, int]]]:
    """Normalise a disulfide request to a mode name or explicit residue pairs.

    An unrecognised mode name used to fall through to geometric auto-detection,
    so ``ssbond="off"`` silently formed every bond it could find.
    """
    if value is None:
        return "none"
    if isinstance(value, bool):
        # ``disulfides: false``/``true`` are the YAML-native spellings the
        # stager already honours; normalise them to the mode names.
        return DEFAULT_DISULFIDE_MODE if value else "none"
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized not in DISULFIDE_MODES:
            raise ValueError(
                f"Unknown PTM disulfide mode {value!r}; choose one of: "
                f"{', '.join(DISULFIDE_MODES)}, or list explicit [chain, resid] pairs."
            )
        return normalized
    if isinstance(value, (list, tuple)):
        pairs: List[Tuple[str, int]] = []
        for index, entry in enumerate(value):
            if (
                isinstance(entry, (str, bytes, dict))
                or not isinstance(entry, (list, tuple))
                or len(entry) != 2
            ):
                raise ValueError(
                    f"Invalid PTM disulfide entry at index {index}; expected a [chain, resid] pair."
                )
            chain, resid = entry
            try:
                # Never ``int(...)`` a residue id: 12.9 would silently bond
                # residue 12 instead of reporting the malformed request.
                resid = _normalise_integer(resid, context="disulfide residue id")
            except ValueError as exc:
                raise ValueError(
                    f"Invalid residue id {resid!r} in PTM disulfide entry {index}."
                ) from exc
            pairs.append((str(chain).strip() or DEFAULT_CHAIN, resid))
        return pairs
    raise ValueError(
        f"Invalid PTM disulfides value {value!r}; expected "
        f"{' or '.join(DISULFIDE_MODES)}, or a list of [chain, resid] pairs."
    )


@dataclass(frozen=True)
class PTMRequest:
    """A single requested modification at one residue position."""

    chain: str
    resid: int
    code: str
    charge: Optional[int] = None

    def __post_init__(self):
        object.__setattr__(self, "chain", str(self.chain or "").strip() or DEFAULT_CHAIN)
        object.__setattr__(
            self,
            "resid",
            _normalise_integer(self.resid, context="PTM residue id"),
        )
        object.__setattr__(self, "code", str(self.code).strip().upper())
        object.__setattr__(
            self,
            "charge",
            _normalise_charge(self.charge, context=f"PTM {self.code or '<unknown>'}"),
        )

    def describe(self) -> str:
        c = f" (charge {self.charge:+d})" if self.charge is not None else ""
        return f"{self.chain}{self.resid} -> {self.code}{c}"


@dataclass
class PTMConfig:
    """Resolved PTM configuration for a build."""

    requests: List[PTMRequest] = field(default_factory=list)
    #: ``"auto"`` | ``"none"`` | explicit list of ``(chain, resid)`` pairs to bond.
    disulfides: Union[str, List[Tuple[str, int]]] = DEFAULT_DISULFIDE_MODE
    amber_phospho_ff: str = DEFAULT_PHOSPHO_FF

    @property
    def enabled(self) -> bool:
        return bool(self.requests) or self.disulfides not in ("none", None, False)

    def validation_errors(self) -> List[str]:
        """Return fatal request errors that must be resolved before staging."""
        problems: List[str] = []
        for r in self.requests:
            d = get_ptm(r.code)
            if d is None:
                problems.append(f"Unknown PTM code '{r.code}' at {r.chain}{r.resid}")
            elif not d.validated:
                problems.append(
                    f"PTM '{r.code}' ({d.name}) has no canonical validated parameters: {d.note}"
                )
            elif r.charge is not None:
                supported_charges = {
                    d.default_charge,
                    *d.charmm_charge_variants.keys(),
                    *d.amber_charge_variants.keys(),
                }
                if r.charge not in supported_charges:
                    problems.append(
                        f"PTM '{r.code}' does not define charge {r.charge:+d}; "
                        f"available charges: {sorted(supported_charges)}"
                    )
        # Directly constructed configs never pass through a parser, so this is
        # where their disulfide/phospho selections are checked.
        for value, requirement in (
            (self.disulfides, _require_disulfide_mode),
            (self.amber_phospho_ff, _require_phospho_ff),
        ):
            try:
                requirement(value)
            except ValueError as exc:
                problems.append(str(exc))
        return problems

    def validate(self) -> List[str]:
        """Alias for :meth:`validation_errors`.

        New code should use ``validation_errors``/``raise_for_errors``, the
        names :class:`~prism.membrane.config.MembraneConfig` uses.
        """
        return self.validation_errors()

    def raise_for_errors(self) -> None:
        errors = self.validation_errors()
        if errors:
            raise ValueError("Invalid PTM configuration: " + "; ".join(errors))


def ptm_defaults() -> dict:
    """Return the default value of every :class:`PTMConfig` field.

    Both parsers read their fallbacks from here (and a generated YAML template
    can too), so each default is written once -- in the dataclass -- instead of
    being hand-copied into every parser, where the copies drift.  Mirrors
    :func:`prism.membrane.config.membrane_defaults`.
    """
    defaults = {}
    for spec in fields(PTMConfig):
        # Containers are rebuilt per call: the caller may hand one to a config.
        if spec.default_factory is not MISSING:
            defaults[spec.name] = spec.default_factory()
        else:
            defaults[spec.name] = spec.default
    return defaults


def parse_ptm_cli(ptm_args: Optional[List[str]], ssbond: Optional[str], phospho_ff: Optional[str]) -> PTMConfig:
    """Build a :class:`PTMConfig` from CLI arguments."""
    defaults = ptm_defaults()
    requests: List[PTMRequest] = []
    for spec in ptm_args or []:
        parts = [p.strip() for p in str(spec).split(":")]
        if len(parts) not in (3, 4):
            raise ValueError(
                f"Invalid --ptm '{spec}'. Expected CHAIN:RESID:CODE[:CHARGE], e.g. A:145:SEP or A:32:TPO:-1"
            )
        chain, resid_s, code = parts[0], parts[1], parts[2].upper()
        try:
            resid = _normalise_integer(resid_s, context="residue id")
        except ValueError as exc:
            raise ValueError(f"Invalid residue id '{resid_s}' in --ptm '{spec}'") from exc
        charge = None
        if len(parts) >= 4 and parts[3] != "":
            try:
                charge = _normalise_charge(parts[3], context=f"--ptm '{spec}'")
            except ValueError as exc:
                raise ValueError(f"Invalid charge '{parts[3]}' in --ptm '{spec}'") from exc
        code = _require_known_code(code, context=f"--ptm '{spec}'")
        requests.append(PTMRequest(chain=chain or DEFAULT_CHAIN, resid=resid, code=code, charge=charge))

    return PTMConfig(
        requests=requests,
        disulfides=_require_disulfide_mode(ssbond or defaults["disulfides"]),
        amber_phospho_ff=_require_phospho_ff(phospho_ff or defaults["amber_phospho_ff"]),
    )


def parse_ptm_yaml(cfg: Optional[dict]) -> PTMConfig:
    """Build a :class:`PTMConfig` from a parsed YAML ``ptm:`` mapping."""
    if cfg is None:
        return PTMConfig()
    if not isinstance(cfg, dict):
        raise ValueError("The YAML 'ptm' section must be a mapping.")

    defaults = ptm_defaults()

    residue_entries = cfg.get("residues", [])
    if residue_entries is None:
        residue_entries = []
    elif not isinstance(residue_entries, (list, tuple)):
        raise ValueError("The YAML 'ptm.residues' value must be a sequence of mappings.")

    requests: List[PTMRequest] = []
    for index, r in enumerate(residue_entries):
        if not isinstance(r, dict):
            raise ValueError(f"Invalid PTM residue entry at index {index}; expected a mapping.")
        # Unknown codes are rejected here, exactly as ``parse_ptm_cli`` rejects
        # them, so a YAML typo fails at parse time instead of mid-build.
        code = _require_known_code(r.get("code", ""), context=f"PTM YAML entry {index}")
        if "resid" not in r:
            raise ValueError(f"PTM residue entry at index {index} is missing 'resid'.")
        try:
            resid = _normalise_integer(r["resid"], context="residue id")
        except ValueError as exc:
            raise ValueError(f"Invalid residue id {r['resid']!r} in PTM YAML entry {index}.") from exc
        requests.append(
            PTMRequest(
                chain=str(r.get("chain", DEFAULT_CHAIN)).strip() or DEFAULT_CHAIN,
                resid=resid,
                code=code,
                charge=_normalise_charge(r.get("charge"), context=f"PTM YAML entry {index}"),
            )
        )

    return PTMConfig(
        requests=requests,
        disulfides=_require_disulfide_mode(cfg.get("disulfides", defaults["disulfides"])),
        amber_phospho_ff=_require_phospho_ff(
            cfg.get("amber_phospho_ff", defaults["amber_phospho_ff"])
        ),
    )
