#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Membrane-build configuration.

Holds the lipid composition, orientation source, and bilayer/solvent box
settings for an all-atom protein-in-membrane build. Parsed from CLI flags or a
YAML ``membrane:`` section.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from numbers import Integral, Real
from typing import List, Optional


# Common lipids understood by PACKMOL-Memgen (AMBER Lipid21 names). This is a
# convenience whitelist for validation/help; PACKMOL-Memgen supports many more.
COMMON_LIPIDS = {
    "POPC", "POPE", "POPS", "POPG", "POPA", "POPI",
    "DOPC", "DOPE", "DOPS", "DPPC", "DMPC", "DLPC",
    "CHL1", "CHOL", "PSM",  # cholesterol / sphingomyelin
}

# Lipid force fields and the protein-FF families they pair with.
LIPID_FF_FAMILY = {
    "lipid17": "amber",     # AMBER Lipid17 (legacy PACKMOL-Memgen option)
    "lipid21": "amber",     # AMBER Lipid21 (PRISM default; selected explicitly)
    "slipids": "amber",     # Slipids-2020 (AMBER-compatible)
    "charmm36": "charmm",   # CHARMM36 lipids
}

# ``--parametrize`` in PACKMOL-Memgen produces an AMBER topology.  PRISM can
# therefore select the AMBER lipid force fields directly, but the advertised
# Slipids/CHARMM entries need a separate parameterisation backend and must not
# be silently passed off as Lipid21.
PACKMOL_PARAMETRIZED_LIPID_FFS = {"lipid17", "lipid21"}

# PACKMOL-Memgen does not consume GROMACS force-field directory names.  Keep
# the translation in one place so a membrane request cannot silently report a
# GROMACS force field while PACKMOL-Memgen parameterizes with its defaults.
PACKMOL_PROTEIN_FF_ALIASES = {
    "amber14sb": "ff14SB",
    "ff14sb": "ff14SB",
    "amber19sb": "ff19SB",
    "ff19sb": "ff19SB",
}
PACKMOL_WATER_MODEL_ALIASES = {
    "tip3p": "tip3p",
    "opc": "opc",
}
PACKMOL_FORCEFIELD_PAIRS = {
    ("ff14SB", "tip3p"),
    ("ff19SB", "opc"),
}

ORIENTATION_METHODS = {"opm", "ppm", "preoriented", "memembed"}
DEFAULT_PRODUCTION_NS = 500.0
MEMBRANE_TIMESTEP_PS = 0.002


def _require_bool(value, field_name: str) -> bool:
    """Return a real boolean without accepting truthy strings or integers."""
    if type(value) is not bool:
        raise ValueError(f"Membrane {field_name} must be a boolean (true or false).")
    return value


def _require_finite_number(
    value,
    field_name: str,
    *,
    allow_zero: bool = False,
) -> float:
    """Parse a YAML/CLI number without silently accepting strings, NaN, or Inf."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"Membrane {field_name} must be a finite number.")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Membrane {field_name} must be finite (NaN/Inf are not allowed).")
    if number < 0 or (number == 0 and not allow_zero):
        qualifier = "non-negative" if allow_zero else "positive"
        raise ValueError(f"Membrane {field_name} must be {qualifier}.")
    return number


def membrane_production_steps(value) -> int:
    """Return production steps, rejecting durations shorter than one step."""
    duration_ns = _require_finite_number(value, "production duration")
    steps = int(duration_ns * 1000 / MEMBRANE_TIMESTEP_PS)
    if steps < 1:
        raise ValueError(
            "Membrane production duration must cover at least one "
            f"{MEMBRANE_TIMESTEP_PS:g} ps integration step."
        )
    return steps


def _require_ratio(value) -> List[int]:
    """Parse an integer lipid-ratio sequence without lossy ``int(...)`` casts."""
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise ValueError("Membrane lipid ratio must be a sequence of positive integers.")
    ratio = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, Integral):
            raise ValueError("Membrane lipid ratios must be positive integers; floats and strings are not allowed.")
        item = int(item)
        if item <= 0:
            raise ValueError("Membrane lipid ratios must be positive integers.")
        ratio.append(item)
    return ratio


def _require_lipids(value) -> List[str]:
    if isinstance(value, str):
        value = [value]
    if isinstance(value, (bytes, dict)) or not isinstance(value, (list, tuple)):
        raise ValueError("Membrane lipids must be a lipid name or a sequence of lipid names.")
    lipids = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise ValueError("Membrane lipid names must be non-empty strings.")
        lipid = item.strip().upper()
        if ":" in lipid or any(character.isspace() for character in lipid):
            raise ValueError(f"Invalid membrane lipid name {item!r}.")
        lipids.append(lipid)
    if not lipids:
        raise ValueError("At least one membrane lipid must be specified.")
    return lipids


def _require_choice(value, field_name: str, choices) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Membrane {field_name} must be one of: {', '.join(sorted(choices))}.")
    normalized = value.strip().lower()
    if normalized not in choices:
        raise ValueError(
            f"Unknown membrane {field_name} {value!r}; choose one of: {', '.join(sorted(choices))}."
        )
    return normalized


def _require_nonempty_string(value, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Membrane {field_name} must be a non-empty string.")
    return value.strip()


@dataclass
class MembraneConfig:
    """Configuration for building a protein-membrane system."""

    enabled: bool = False
    #: Lipid composition (one or more lipid names).
    lipids: List[str] = field(default_factory=lambda: ["POPC"])
    #: Molar ratio per lipid (parallel to ``lipids``); defaults to equal parts.
    ratio: List[int] = field(default_factory=list)
    #: Orientation source: "opm" (fetch by PDB id), "ppm" (local/server),
    #: "preoriented" (input already membrane-aligned), or "memembed".
    orient: str = "preoriented"
    #: PDB id to fetch a pre-oriented structure from the OPM database.
    pdb_id: Optional[str] = None
    #: Lipid force field (selects the protein-FF family pairing).
    lipid_ff: str = "lipid21"
    #: Protein and water force fields used by PACKMOL-Memgen parameterization.
    #: These are synchronized with PRISMBuilder's resolved global settings.
    protein_forcefield: str = "ff14SB"
    water_model: str = "tip3p"
    #: Salt concentration (mol/L) for the aqueous phase.
    salt_concentration: float = 0.15
    positive_ion: str = "Na+"
    negative_ion: str = "Cl-"
    #: Water-layer thickness above/below the bilayer (Angstrom; --dist_wat).
    water_thickness: float = 17.5
    #: Minimum protein-to-box distance in the membrane plane (Angstrom; --dist).
    xy_distance: float = 17.0
    #: Run a short AMBER (sander/pmemd) minimization on the packed system
    #: before conversion (works around GROMACS-crash regressions).
    minimize: bool = True
    temperature: float = 303.15
    production_ns: float = DEFAULT_PRODUCTION_NS

    def family(self) -> str:
        """Protein force-field family implied by the lipid force field."""
        lipid_ff = self.lipid_ff.lower() if isinstance(self.lipid_ff, str) else ""
        try:
            return LIPID_FF_FAMILY[lipid_ff]
        except KeyError as exc:
            raise ValueError(f"Unknown lipid force field {self.lipid_ff!r}.") from exc

    def resolved_ratio(self) -> List[int]:
        if not isinstance(self.lipids, (list, tuple)) or not self.lipids:
            raise ValueError("At least one membrane lipid must be specified.")
        ratio = _require_ratio(self.ratio)
        if not ratio:
            return [1] * len(self.lipids)
        if len(ratio) != len(self.lipids):
            raise ValueError(
                f"Lipid ratio count ({len(ratio)}) must match lipid count ({len(self.lipids)})."
            )
        return ratio

    def packmol_forcefield_pair(self):
        """Return validated PACKMOL-Memgen protein/water force-field names."""
        protein_key = str(self.protein_forcefield or "").strip().lower()
        water_key = str(self.water_model or "").strip().lower()
        try:
            protein_ff = PACKMOL_PROTEIN_FF_ALIASES[protein_key]
        except KeyError as exc:
            supported = ", ".join(sorted({"amber14sb", "amber19sb"}))
            raise ValueError(
                f"Unsupported membrane protein force field {self.protein_forcefield!r}; "
                f"supported PRISM names: {supported}."
            ) from exc
        try:
            water_ff = PACKMOL_WATER_MODEL_ALIASES[water_key]
        except KeyError as exc:
            supported = ", ".join(sorted(PACKMOL_WATER_MODEL_ALIASES))
            raise ValueError(
                f"Unsupported membrane water model {self.water_model!r}; supported models: {supported}."
            ) from exc
        if (protein_ff, water_ff) not in PACKMOL_FORCEFIELD_PAIRS:
            raise ValueError(
                "Unsupported membrane protein/water force-field pairing "
                f"{self.protein_forcefield!r}/{self.water_model!r}; use "
                "amber14sb/tip3p or amber19sb/opc."
            )
        return protein_ff, water_ff

    def validation_errors(self) -> List[str]:
        """Return fatal configuration errors that must block external tools."""
        problems = []

        if type(self.enabled) is not bool:
            problems.append("Membrane enabled must be a boolean (true or false).")
        if type(self.minimize) is not bool:
            problems.append("Membrane minimize must be a boolean (true or false).")

        try:
            lipids = _require_lipids(self.lipids)
        except ValueError as exc:
            problems.append(str(exc))
            lipids = []

        try:
            orient = _require_choice(self.orient, "orientation method", ORIENTATION_METHODS)
        except ValueError as exc:
            problems.append(str(exc))
            orient = ""

        lipid_ff = self.lipid_ff.strip().lower() if isinstance(self.lipid_ff, str) else ""
        if lipid_ff not in LIPID_FF_FAMILY:
            problems.append(f"Unknown lipid force field {self.lipid_ff!r}.")
        elif lipid_ff not in PACKMOL_PARAMETRIZED_LIPID_FFS:
            problems.append(
                "The automated membrane backend does not support "
                f"lipid_ff={self.lipid_ff!r}; use lipid17 or lipid21."
            )

        if orient == "opm" and (not isinstance(self.pdb_id, str) or not self.pdb_id.strip()):
            problems.append("orient=opm requires a non-empty membrane PDB id.")
        elif self.pdb_id is not None and (not isinstance(self.pdb_id, str) or not self.pdb_id.strip()):
            problems.append("Membrane PDB id must be a non-empty string when provided.")

        try:
            ratio = _require_ratio(self.ratio)
            if ratio and lipids and len(ratio) != len(lipids):
                problems.append(
                    f"Lipid ratio count ({len(ratio)}) must match lipid count ({len(lipids)})."
                )
        except ValueError as exc:
            problems.append(str(exc))

        for value, field_name, allow_zero in (
            (self.salt_concentration, "salt concentration", True),
            (self.water_thickness, "water thickness", False),
            (self.xy_distance, "protein-to-box distance", False),
            (self.temperature, "temperature", False),
        ):
            try:
                _require_finite_number(value, field_name, allow_zero=allow_zero)
            except ValueError as exc:
                problems.append(str(exc))

        try:
            membrane_production_steps(self.production_ns)
        except ValueError as exc:
            problems.append(str(exc))

        for value, field_name in (
            (self.positive_ion, "positive ion"),
            (self.negative_ion, "negative ion"),
        ):
            try:
                _require_nonempty_string(value, field_name)
            except ValueError as exc:
                problems.append(str(exc))

        if lipid_ff in PACKMOL_PARAMETRIZED_LIPID_FFS:
            try:
                self.packmol_forcefield_pair()
            except ValueError as exc:
                problems.append(str(exc))
        return problems

    def validation_warnings(self) -> List[str]:
        """Return advisory messages which do not make a build scientifically invalid."""
        warnings = []
        if isinstance(self.lipids, (list, tuple)):
            for lipid in self.lipids:
                if isinstance(lipid, str) and lipid.strip().upper() not in COMMON_LIPIDS:
                    warnings.append(
                        f"Lipid '{lipid}' is not in the common-lipid whitelist "
                        "(it may still be valid for PACKMOL-Memgen)."
                    )
        return warnings

    def validate(self) -> List[str]:
        """Return all fatal errors followed by advisory warnings (legacy API)."""
        return self.validation_errors() + self.validation_warnings()

    def raise_for_errors(self) -> None:
        errors = self.validation_errors()
        if errors:
            raise ValueError("Invalid membrane configuration: " + "; ".join(errors))


def parse_membrane_cli(
    membrane: bool,
    lipid: Optional[List[str]],
    lipid_ratio: Optional[List[int]],
    orient: Optional[str],
    pdb_id: Optional[str],
    lipid_ff: Optional[str],
    salt_concentration: float = 0.15,
    temperature: float = 303.15,
    production_ns: float = DEFAULT_PRODUCTION_NS,
    positive_ion: str = "Na+",
    negative_ion: str = "Cl-",
    protein_forcefield: str = "amber14sb",
    water_model: str = "tip3p",
) -> MembraneConfig:
    """Build a :class:`MembraneConfig` from CLI arguments."""
    config = MembraneConfig(
        enabled=_require_bool(membrane, "enabled"),
        lipids=_require_lipids(lipid or ["POPC"]),
        ratio=_require_ratio(lipid_ratio),
        orient=_require_choice(orient or "preoriented", "orientation method", ORIENTATION_METHODS),
        pdb_id=pdb_id.strip() if isinstance(pdb_id, str) else pdb_id,
        lipid_ff=_require_choice(lipid_ff or "lipid21", "lipid force field", LIPID_FF_FAMILY),
        salt_concentration=_require_finite_number(
            salt_concentration, "salt concentration", allow_zero=True
        ),
        positive_ion=_require_nonempty_string(positive_ion, "positive ion"),
        negative_ion=_require_nonempty_string(negative_ion, "negative ion"),
        protein_forcefield=_require_nonempty_string(protein_forcefield, "protein force field"),
        water_model=_require_nonempty_string(water_model, "water model"),
        temperature=_require_finite_number(temperature, "temperature"),
        production_ns=_require_finite_number(production_ns, "production duration"),
    )
    config.raise_for_errors()
    return config


def parse_membrane_yaml(
    cfg: Optional[dict],
    temperature: float = 303.15,
    production_ns: float = DEFAULT_PRODUCTION_NS,
) -> MembraneConfig:
    """Build a :class:`MembraneConfig` from a parsed YAML ``membrane:`` mapping."""
    if cfg is None:
        return MembraneConfig(enabled=False)
    if not isinstance(cfg, dict):
        raise ValueError("The YAML 'membrane' section must be a mapping.")

    pdb_id = cfg.get("pdb_id")
    if isinstance(pdb_id, str):
        pdb_id = pdb_id.strip()

    if "lipids" in cfg:
        lipid_value = cfg["lipids"]
    elif "lipid" in cfg:
        lipid_value = cfg["lipid"]
    else:
        lipid_value = ["POPC"]

    config = MembraneConfig(
        enabled=_require_bool(cfg.get("enabled", True), "enabled"),
        lipids=_require_lipids(lipid_value),
        ratio=_require_ratio(cfg.get("ratio", [])),
        orient=_require_choice(
            cfg.get("orient", "preoriented"), "orientation method", ORIENTATION_METHODS
        ),
        pdb_id=pdb_id,
        lipid_ff=_require_choice(
            cfg.get("lipid_ff", "lipid21"), "lipid force field", LIPID_FF_FAMILY
        ),
        protein_forcefield=_require_nonempty_string(
            cfg.get("protein_forcefield", cfg.get("forcefield", "amber14sb")),
            "protein force field",
        ),
        water_model=_require_nonempty_string(
            cfg.get("water_model", cfg.get("water", "tip3p")), "water model"
        ),
        salt_concentration=_require_finite_number(
            cfg.get("salt_concentration", 0.15), "salt concentration", allow_zero=True
        ),
        positive_ion=_require_nonempty_string(cfg.get("positive_ion", "Na+"), "positive ion"),
        negative_ion=_require_nonempty_string(cfg.get("negative_ion", "Cl-"), "negative ion"),
        water_thickness=_require_finite_number(
            cfg.get("water_thickness", 17.5), "water thickness"
        ),
        xy_distance=_require_finite_number(
            cfg.get("xy_distance", 17.0), "protein-to-box distance"
        ),
        minimize=_require_bool(cfg.get("minimize", True), "minimize"),
        temperature=_require_finite_number(cfg.get("temperature", temperature), "temperature"),
        production_ns=_require_finite_number(
            cfg.get("production_ns", cfg.get("production_time_ns", production_ns)),
            "production duration",
        ),
    )
    config.raise_for_errors()
    return config
