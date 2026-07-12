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
from typing import List, Optional, Tuple, Union

from .catalog import get_ptm, is_known_ptm


@dataclass(frozen=True)
class PTMRequest:
    """A single requested modification at one residue position."""

    chain: str
    resid: int
    code: str
    charge: Optional[int] = None

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
        """Return a list of human-readable warnings/errors (empty if clean)."""
        problems: List[str] = []
        for r in self.requests:
            d = get_ptm(r.code)
            if d is None:
                problems.append(f"Unknown PTM code '{r.code}' at {r.chain}{r.resid}")
            elif not d.validated:
                problems.append(
                    f"PTM '{r.code}' ({d.name}) has no canonical validated parameters: {d.note}"
                )
        if self.amber_phospho_ff not in ("phosaa19SB", "phosaa14SB", "phosaa10"):
            problems.append(f"Unknown amber_phospho_ff '{self.amber_phospho_ff}'")
        return problems


def parse_ptm_cli(ptm_args: Optional[List[str]], ssbond: Optional[str], phospho_ff: Optional[str]) -> PTMConfig:
    """Build a :class:`PTMConfig` from CLI arguments."""
    requests: List[PTMRequest] = []
    for spec in ptm_args or []:
        parts = [p.strip() for p in str(spec).split(":")]
        if len(parts) < 3:
            raise ValueError(
                f"Invalid --ptm '{spec}'. Expected CHAIN:RESID:CODE[:CHARGE], e.g. A:145:SEP or A:32:TPO:-1"
            )
        chain, resid_s, code = parts[0], parts[1], parts[2].upper()
        try:
            resid = int(resid_s)
        except ValueError:
            raise ValueError(f"Invalid residue id '{resid_s}' in --ptm '{spec}'")
        charge = None
        if len(parts) >= 4 and parts[3] != "":
            try:
                charge = int(parts[3])
            except ValueError:
                raise ValueError(f"Invalid charge '{parts[3]}' in --ptm '{spec}'")
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
    if not cfg:
        return PTMConfig(requests=[], disulfides="auto")

    requests: List[PTMRequest] = []
    for r in cfg.get("residues", []) or []:
        code = str(r.get("code", "")).upper()
        requests.append(
            PTMRequest(
                chain=str(r.get("chain", "A")),
                resid=int(r["resid"]),
                code=code,
                charge=r.get("charge"),
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
