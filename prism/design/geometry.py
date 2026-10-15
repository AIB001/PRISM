"""PDB atom identity and distances, independent of the generation backend."""

import math
from pathlib import Path

import numpy as np

from .errors import require

AMINO = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP ASH GLH CYX CYM LYN MSE".split()
)
METALS = {"ZN", "MG", "MN", "CA", "FE", "CO", "NI", "CU", "CD", "HG", "NA", "K"}


def pdb_atoms(path):
    atoms = {}
    model = "1"
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("MODEL"):
            model = line[10:14].strip() or "1"
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        resname = line[17:20].strip()
        name = line[12:16].strip()
        element = line[76:78].strip().upper()
        require(bool(element), f"PDB atom has no element: {name}", "PDB_ELEMENT_MISSING")
        if element in {"H", "D"}:
            continue
        try:
            xyz = [float(line[i : i + 8]) for i in (30, 38, 46)]
            occupancy = float(line[54:60].strip() or "1")
        except ValueError as exc:
            raise ValueError(f"Invalid PDB coordinate/occupancy in {path}") from exc
        require(all(math.isfinite(x) for x in xyz + [occupancy]), "Non-finite PDB value", "PDB_NONFINITE")
        key = (model, line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), resname, name)
        atom = {
            "key": key,
            "xyz": xyz,
            "element": element,
            "occupancy": occupancy,
            "altloc": line[16:17].strip(),
            "protein": line.startswith("ATOM") or resname in AMINO,
        }
        rank = lambda a: (a["occupancy"], a["altloc"] == "", a["altloc"] == "A")
        if key not in atoms or rank(atom) > rank(atoms[key]):
            atoms[key] = atom
    require(bool(atoms), f"No PDB heavy atoms: {path}", "PDB_EMPTY")
    return list(atoms.values())


def pocket_review(protein_path, pocket_path, ligand_xyz):
    protein_all = pdb_atoms(protein_path)
    protein = [a for a in protein_all if a["protein"]]
    pocket = [a for a in pdb_atoms(pocket_path) if a["protein"]]
    require(protein and pocket, "Empty protein/pocket polymer atoms", "PDB_EMPTY")
    by_key = {a["key"]: a for a in protein}
    matched = [(a, by_key[a["key"]]) for a in pocket if a["key"] in by_key]
    fraction = len(matched) / len(pocket)
    displacement = max((np.linalg.norm(np.array(a["xyz"]) - b["xyz"]) for a, b in matched), default=float("inf"))
    require(
        fraction >= 0.95 and displacement <= 0.05,
        "Pocket identity/coordinate frame disagrees with receptor",
        "POCKET_FRAME_MISMATCH",
    )
    xyz = np.array([a["xyz"] for a in protein])
    nearest = np.linalg.norm(np.asarray(ligand_xyz)[:, None, :] - xyz[None, :, :], axis=2).min(axis=1)
    pocket_xyz = np.array([a["xyz"] for a in pocket])
    pocket_nearest = np.linalg.norm(np.asarray(ligand_xyz)[:, None, :] - pocket_xyz[None, :, :], axis=2).min(axis=1)
    require(float(np.mean(pocket_nearest <= 4)) >= 0.2, "Ligand does not overlap the supplied pocket", "POCKET_OVERLAP")
    metals = [a for a in protein_all if not a["protein"] and a["element"] in METALS]
    near_metals = [
        ":".join(a["key"]) for a in metals if np.linalg.norm(np.asarray(ligand_xyz) - a["xyz"], axis=1).min() <= 3
    ]
    return {
        "pocket_key_match_fraction": fraction,
        "pocket_max_displacement_A": float(displacement),
        "nearest_protein_heavy_distance_A": float(nearest.min()),
        "ligand_fraction_within_4A": float(np.mean(nearest <= 4)),
        "ligand_fraction_within_pocket_4A": float(np.mean(pocket_nearest <= 4)),
        "metals_within_3A": near_metals,
    }, xyz
