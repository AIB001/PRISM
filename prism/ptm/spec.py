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
      amber_phospho_ff: phosaa19SB     # phosaa19SB | phosaa14SB
      residues:
        - {chain: A, resid: 145, code: SEP}
        - {chain: A, resid: 32,  code: TPO, charge: -1}
"""

from __future__ import annotations

from dataclasses import dataclass, field
from numbers import Integral
import re
from typing import List, Optional, Tuple, Union

from .catalog import get_ptm, is_known_ptm


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


@dataclass(frozen=True)
class PTMRequest:
    """A single requested modification at one residue position."""

    chain: str
    resid: int
    code: str
    charge: Optional[int] = None

    def __post_init__(self):
        object.__setattr__(self, "chain", str(self.chain or "").strip() or "A")
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
    # disulfides: "auto" | "none" | explicit list of (chain, resid) pairs to bond
    disulfides: Union[str, List[Tuple[str, int]]] = "auto"
    amber_phospho_ff: str = "phosaa19SB"

    @property
    def enabled(self) -> bool:
        return bool(self.requests) or self.disulfides not in ("none", None, False)

    def validate(self) -> List[str]:
        """Return user-request errors that must be resolved before staging."""
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
        if self.amber_phospho_ff not in ("phosaa19SB", "phosaa14SB", "phosaa10"):
            problems.append(f"Unknown amber_phospho_ff '{self.amber_phospho_ff}'")
        return problems


def parse_ptm_cli(ptm_args: Optional[List[str]], ssbond: Optional[str], phospho_ff: Optional[str]) -> PTMConfig:
    """Build a :class:`PTMConfig` from CLI arguments."""
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
        if not is_known_ptm(code):
            raise ValueError(
                f"Unknown PTM code '{code}'. Known codes: see `prism --list-ptms`."
            )
        requests.append(PTMRequest(chain=chain or "A", resid=resid, code=code, charge=charge))

    disulfides: Union[str, List[Tuple[str, int]]] = (ssbond or "auto").lower()
    return PTMConfig(
        requests=requests,
        disulfides=disulfides,
        amber_phospho_ff=(phospho_ff or "phosaa19SB"),
    )


def parse_ptm_yaml(cfg: Optional[dict]) -> PTMConfig:
    """Build a :class:`PTMConfig` from a parsed YAML ``ptm:`` mapping."""
    if cfg is None:
        return PTMConfig(requests=[], disulfides="auto")
    if not isinstance(cfg, dict):
        raise ValueError("The YAML 'ptm' section must be a mapping.")

    residue_entries = cfg.get("residues", [])
    if residue_entries is None:
        residue_entries = []
    elif not isinstance(residue_entries, (list, tuple)):
        raise ValueError("The YAML 'ptm.residues' value must be a sequence of mappings.")

    requests: List[PTMRequest] = []
    for index, r in enumerate(residue_entries):
        if not isinstance(r, dict):
            raise ValueError(f"Invalid PTM residue entry at index {index}; expected a mapping.")
        code = str(r.get("code", "")).strip().upper()
        if "resid" not in r:
            raise ValueError(f"PTM residue entry at index {index} is missing 'resid'.")
        try:
            resid = _normalise_integer(r["resid"], context="residue id")
        except ValueError as exc:
            raise ValueError(f"Invalid residue id {r['resid']!r} in PTM YAML entry {index}.") from exc
        requests.append(
            PTMRequest(
                chain=str(r.get("chain", "A")).strip() or "A",
                resid=resid,
                code=code,
                charge=_normalise_charge(r.get("charge"), context=f"PTM YAML entry {index}"),
            )
        )

    disulfides = cfg.get("disulfides", "auto")
    if isinstance(disulfides, str):
        disulfides = disulfides.lower()
    elif isinstance(disulfides, list):
        disulfides = [(str(c), int(i)) for c, i in disulfides]

    return PTMConfig(
        requests=requests,
        disulfides=disulfides,
        amber_phospho_ff=str(cfg.get("amber_phospho_ff", "phosaa19SB")),
    )
