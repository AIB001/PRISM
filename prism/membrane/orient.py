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

import os
from typing import Optional


OPM_PDB_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb_id}.pdb"
OPM_ALT_URL = "https://opm.phar.umich.edu/pdb/{pdb_id}.pdb"


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
            if data and b"ATOM" in data[:5000]:
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


def orient_protein(input_pdb: str, output_pdb: str, method: str, pdb_id: Optional[str] = None) -> dict:
    """Orient a protein along the membrane normal.

    Returns ``{"oriented_pdb", "method", "note"}``. Falls back to copying the
    input (treated as preoriented) when an automated route is unavailable.
    """
    import shutil

    method = (method or "preoriented").lower()
    note = ""

    if method == "opm" and pdb_id:
        fetched = fetch_opm(pdb_id, output_pdb)
        if fetched:
            n = strip_opm_dummy_atoms(output_pdb, output_pdb)
            return {"oriented_pdb": output_pdb, "method": "opm",
                    "note": f"Fetched oriented structure {pdb_id} from OPM (removed {n} DUMMY atoms)."}
        note = f"OPM fetch failed for '{pdb_id}'; treating input as preoriented."

    if method == "ppm":
        note = ("PPM 3.0 orientation requires the web server (https://opm.phar.umich.edu/ppm_server3) "
                "or the standalone 'immers' binary; run it externally and pass the oriented PDB with "
                "--membrane-orient preoriented.")

    # preoriented / fallback
    if os.path.abspath(input_pdb) != os.path.abspath(output_pdb):
        shutil.copy2(input_pdb, output_pdb)
    return {"oriented_pdb": output_pdb, "method": "preoriented", "note": note or "Input treated as preoriented."}
