#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Membrane-normal orientation of a protein (z-axis = bilayer normal).

Orientation is a prerequisite for building a bilayer around a protein. Three
routes are supported, in increasing order of cost:

* **OPM** — if the structure's PDB id is in the Orientations of Proteins in
  Membranes database, fetch the pre-oriented coordinates over HTTP (zero compute).
* **PPM** — de-novo orientation via the PPM 3.0 server / standalone (gated).
* **preoriented** — trust that the input is already membrane-aligned.

Only the OPM fetch is implemented as a network call here; PPM/MEMEMBED are gated
hooks with actionable guidance, because they require an external server or a
license-restricted binary (see PRISM_Capability_Expansion_Roadmap.md).
"""

from __future__ import annotations

import math
import os
from typing import Optional


OPM_PDB_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb_id}.pdb"
OPM_ALT_URL = "https://opm.phar.umich.edu/pdb/{pdb_id}.pdb"
MAX_PEPTIDE_CN_DISTANCE_ANGSTROM = 2.0


def fetch_opm(pdb_id: str, out_path: str, timeout: int = 30) -> Optional[str]:
    """Fetch a membrane-oriented structure from the OPM database by PDB id.

    Returns the path on success, or None if the structure is not in OPM / the
    download fails. The downloaded PDB carries DUMMY boundary atoms marking the
    membrane planes — strip them before topology generation.
    """
    import urllib.request
    import urllib.error

    pdb_id = pdb_id.strip().lower()
    for tmpl in (OPM_PDB_URL, OPM_ALT_URL):
        url = tmpl.format(pdb_id=pdb_id)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PRISM/membrane"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            # OPM files can carry a long metadata header.  Limiting the check to
            # the first 5 kB rejects otherwise valid large structures.
            has_coordinate_record = any(
                line.startswith((b"ATOM  ", b"HETATM"))
                for line in (data or b"").splitlines()
            )
            if has_coordinate_record:
                with open(out_path, "wb") as fh:
                    fh.write(data)
                return out_path
        except (urllib.error.URLError, urllib.error.HTTPError, OSError):
            continue
    return None


def strip_opm_dummy_atoms(in_pdb: str, out_pdb: str) -> int:
    """Remove OPM/PPM DUMMY boundary atoms (membrane-plane markers).

    Returns the number of DUMMY atoms removed.
    """
    removed = 0
    with open(in_pdb, "r") as fh:
        lines = fh.readlines()
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            atom = line[12:16].strip()
            if resname in ("DUM", "DUMMY") or atom in ("N", "O") and resname == "DUM":
                removed += 1
                continue
        out.append(line)
    with open(out_pdb, "w") as fh:
        fh.writelines(out)
    return removed


def insert_ter_at_chain_breaks(
    in_pdb: str,
    out_pdb: str,
    max_cn_distance: float = MAX_PEPTIDE_CN_DISTANCE_ANGSTROM,
) -> int:
    """Insert ``TER`` before discontinuous peptide-coordinate segments.

    Some OPM entries retain one chain identifier across unresolved residue
    ranges but omit a ``TER`` record. PACKMOL-Memgen/tleap then renumbers the
    observed residues consecutively and creates a peptide bond across the gap.
    The resulting topology can contain excluded atom pairs more than a
    nanometre apart and is not safe to simulate.

    A break is inserted between adjacent peptide-like residues when their chain
    changes or the preceding C to following N distance exceeds
    ``max_cn_distance``. Discontinuous residue numbering is used only when that
    geometry cannot be measured. Existing ``TER`` and MODEL boundaries are
    preserved. The input may equal the output path.
    """
    if not math.isfinite(max_cn_distance) or max_cn_distance <= 0:
        raise ValueError("max_cn_distance must be a positive finite value")

    with open(in_pdb, "r") as fh:
        lines = fh.readlines()

    residues = []
    residue_by_key = {}
    model = 0
    for index, line in enumerate(lines):
        if line.startswith("MODEL"):
            model += 1
            continue
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 27:
            continue
        try:
            resid = int(line[22:26])
        except (ValueError, IndexError):
            continue
        key = (model, line[21], resid, line[26] if len(line) > 26 else " ")
        residue = residue_by_key.get(key)
        if residue is None:
            residue = {
                "key": key,
                "first_index": index,
                "last_index": index,
                "atoms": {},
                "has_atom_record": False,
            }
            residue_by_key[key] = residue
            residues.append(residue)
        residue["last_index"] = index
        if line.startswith("ATOM"):
            residue["has_atom_record"] = True
        atom_name = line[12:16].strip().upper()
        try:
            coordinate = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except (ValueError, IndexError):
            coordinate = None
        if coordinate is not None and atom_name not in residue["atoms"]:
            residue["atoms"][atom_name] = coordinate

    insertion_indices = set()
    for previous, current in zip(residues, residues[1:]):
        previous_model, previous_chain, previous_resid, previous_icode = previous["key"]
        current_model, current_chain, current_resid, current_icode = current["key"]
        if current_model != previous_model:
            continue
        between = lines[previous["last_index"] + 1 : current["first_index"]]
        if any(line.startswith(("TER", "ENDMDL", "MODEL")) for line in between):
            continue
        # ATOM records identify polymer residues.  Retain support for modified
        # amino acids written as HETATM when they carry a complete peptide
        # backbone, but do not infer breaks across ordinary ligands/solvent.
        previous_is_polymer = previous["has_atom_record"] or {
            "N", "CA", "C"
        }.issubset(previous["atoms"])
        current_is_polymer = current["has_atom_record"] or {
            "N", "CA", "C"
        }.issubset(current["atoms"])
        if not (previous_is_polymer and current_is_polymer):
            continue

        chain_changed = current_chain != previous_chain
        numbering_contiguous = (
            current_chain == previous_chain
            and (
                current_resid == previous_resid + 1
                or (
                    current_resid == previous_resid
                    and current_icode.strip()
                    and current_icode != previous_icode
                )
            )
        )
        c_coord = previous["atoms"].get("C")
        n_coord = current["atoms"].get("N")
        geometry_available = c_coord is not None and n_coord is not None
        geometry_broken = (
            geometry_available and math.dist(c_coord, n_coord) > max_cn_distance
        )

        # Residue numbers are identifiers, not connectivity: a valid polymer
        # can be numbered with gaps.  Prefer measured C-N geometry whenever it
        # is available, and use numbering only as a conservative fallback.
        if (
            chain_changed
            or geometry_broken
            or (not geometry_available and not numbering_contiguous)
        ):
            insertion_indices.add(current["first_index"])

    output = []
    inserted = 0
    for index, line in enumerate(lines):
        if index in insertion_indices and not (
            output and output[-1].startswith(("TER", "ENDMDL", "MODEL"))
        ):
            output.append("TER\n")
            inserted += 1
        output.append(line)
    with open(out_pdb, "w") as fh:
        fh.writelines(output)
    return inserted


def orient_protein(input_pdb: str, output_pdb: str, method: str, pdb_id: Optional[str] = None) -> dict:
    """Orient a protein along the membrane normal.

    Returns ``{"oriented_pdb", "method", "note"}``. Automated orientation
    requests fail closed: an un-oriented input must never be relabelled as
    preoriented merely because a network service or backend is unavailable.
    """
    import shutil

    method = (method or "preoriented").lower()
    note = ""

    if method not in {"opm", "ppm", "preoriented", "memembed"}:
        raise ValueError(f"Unsupported membrane orientation method: {method!r}")
    if method == "opm" and not pdb_id:
        raise ValueError("OPM orientation requires a PDB id.")

    if method == "opm":
        fetched = fetch_opm(pdb_id, output_pdb)
        if fetched:
            n = strip_opm_dummy_atoms(output_pdb, output_pdb)
            breaks = insert_ter_at_chain_breaks(output_pdb, output_pdb)
            return {"oriented_pdb": output_pdb, "method": "opm",
                    "note": (
                        f"Fetched oriented structure {pdb_id} from OPM "
                        f"(removed {n} DUMMY atoms; inserted {breaks} TER chain-break records)."
                    )}
        raise RuntimeError(
            f"OPM fetch failed for '{pdb_id}'. Check network/database availability, "
            "or orient the structure externally and explicitly use 'preoriented'."
        )

    if method == "ppm":
        raise NotImplementedError(
            "PPM 3.0 orientation is not automated by PRISM. Use the PPM web server or "
            "standalone 'immers' binary, then pass its output with 'preoriented'."
        )

    if method == "memembed":
        note = "Orientation delegated to PACKMOL-Memgen's MEMEMBED stage."

    # Explicitly preoriented input, or MEMEMBED delegated to PACKMOL-Memgen.
    if os.path.abspath(input_pdb) != os.path.abspath(output_pdb):
        shutil.copy2(input_pdb, output_pdb)
    breaks = insert_ter_at_chain_breaks(output_pdb, output_pdb)
    if breaks:
        break_note = f" Inserted {breaks} TER chain-break record(s)."
    else:
        break_note = ""
    return {
        "oriented_pdb": output_pdb,
        "method": method if method == "memembed" else "preoriented",
        "note": (note or "Input treated as preoriented.") + break_note,
    }
