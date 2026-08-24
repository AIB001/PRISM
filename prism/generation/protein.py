"""Dependency-free PDB pocket extraction for external generation models."""

import math
from pathlib import Path
from typing import List, Sequence, Set, Tuple

from .errors import InputError
from .pocket import coordinates_from_structure
from .types import PocketSpec


Coordinate = Tuple[float, float, float]
ResidueKey = Tuple[str, str, str]


def _distance(first: Coordinate, second: Coordinate) -> float:
    return math.sqrt(sum((left - right) ** 2 for left, right in zip(first, second)))


def _pdb_atoms(path: Path) -> List[Tuple[str, ResidueKey, Coordinate]]:
    atoms = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise InputError(f"Unable to read protein PDB '{path}': {exc}") from exc
    for line in lines:
        if not line.startswith("ATOM  "):
            continue
        alternate = line[16:17]
        if alternate not in {" ", "A"}:
            continue
        try:
            coordinate = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError as exc:
            raise InputError(f"Invalid ATOM coordinate in protein PDB '{path}'") from exc
        if not all(math.isfinite(value) for value in coordinate):
            raise InputError(f"Non-finite ATOM coordinate in protein PDB '{path}'")
        key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip())
        atoms.append((line, key, coordinate))
    if not atoms:
        raise InputError(f"No protein ATOM records found in '{path}'")
    return atoms


def _residue_keys(entries: Sequence[str]) -> Set[ResidueKey]:
    keys = set()
    for entry in entries:
        chain, residue = entry.split(":", 1)
        insertion = residue[-1] if residue[-1].isalpha() else ""
        number = residue[:-1] if insertion else residue
        keys.add((chain, number, insertion))
    return keys


def extract_pocket_pdb(protein: Path, pocket: PocketSpec, destination: Path) -> Path:
    """Select whole residues and write a TargetDiff-compatible pocket PDB."""
    if protein.suffix.lower() != ".pdb":
        raise InputError(
            "Automatic pocket extraction currently requires a PDB protein; "
            "provide a cropped pocket PDB for CIF/mmCIF proteins"
        )
    atoms = _pdb_atoms(protein)
    selected: Set[ResidueKey]
    if pocket.kind == "residues":
        requested = _residue_keys(pocket.residues)
        available = {key for _line, key, _coordinate in atoms}
        missing = sorted(requested - available)
        if missing:
            formatted = [f"{chain}:{number}{insertion}" for chain, number, insertion in missing]
            raise InputError("Pocket residues not found in protein: " + ", ".join(formatted))
        selected = requested
    elif pocket.kind == "reference_ligand" and pocket.reference_ligand is not None:
        radius = pocket.radius if pocket.radius is not None else 10.0
        ligand_coordinates = coordinates_from_structure(pocket.reference_ligand)
        selected = {
            key
            for _line, key, coordinate in atoms
            if any(_distance(coordinate, ligand_coordinate) <= radius for ligand_coordinate in ligand_coordinates)
        }
    elif pocket.center is not None:
        radius = pocket.radius if pocket.radius is not None else 10.0
        selected = {
            key for _line, key, coordinate in atoms if _distance(coordinate, pocket.center) <= radius
        }
    else:
        raise InputError(f"Cannot extract a pocket PDB from pocket kind '{pocket.kind}'")
    if not selected:
        raise InputError("Pocket extraction selected no protein residues")

    destination.parent.mkdir(parents=True, exist_ok=True)
    output_lines = ["HEADER    PRISM GENERATED POCKET"]
    output_lines.extend(line for line, key, _coordinate in atoms if key in selected)
    output_lines.append("END")
    destination.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return destination
