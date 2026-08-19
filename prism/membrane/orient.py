#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Membrane-normal orientation of a protein (z-axis = bilayer normal).

Orientation is a prerequisite for building a bilayer around a protein. The
routes are:

* **auto** — measure the input against the bilayer frame and dispatch: keep
  coordinates that are already framed, run MEMEMBED on everything else. This is
  the route a bare ``--membrane`` should take, because it is the only one that
  is correct for both kinds of input.
* **MEMEMBED** — fit the membrane slab de novo with the ``memembed`` binary
  (~5 s; ships with AmberTools alongside packmol-memgen).
* **OPM** — if the structure's PDB id is in the Orientations of Proteins in
  Membranes database, fetch the pre-oriented coordinates over HTTP (zero compute).
* **PPM** — de-novo orientation via the PPM 3.0 server / standalone (gated,
  because it needs an external server).
* **preoriented** — trust that the input is already membrane-aligned.

Whichever route runs, "trust" is checked rather than assumed: PACKMOL-Memgen
always builds the bilayer centred on z = 0 and places the solute at its *input*
coordinates (``fixed 0. 0. 0 0. 0. 0.``). A structure taken straight from the
RCSB is essentially never in that frame, and packing it produces a protein
floating in bulk water — a scientifically worthless system that otherwise builds
"successfully" with no warning. :func:`membrane_alignment_report` is that
measurement; it is what ``auto`` dispatches on, and what every route is checked
against afterwards. The check is not redundant after an orientation run:
MEMEMBED fits a hydrophobic belt, so a protein that has none (a soluble protein
handed to ``--membrane``) still comes back unframed, and a fit that failed is
exactly as wrong as an un-oriented input.

Everything this module reports is returned as notes carrying an explicit
severity, never as a message that prefixes its own: only the caller knows how a
message is rendered, and "the protein is nowhere near the bilayer" must not
print at the weight of "a TER record was inserted".
"""

from __future__ import annotations

import math
import os
import statistics
from typing import Dict, List, Optional, Tuple

from . import config as config_vocab
from .config import DEFAULT_CENTER_TOLERANCE_A


#: MEMEMBED binary.  It ships with AmberTools alongside packmol-memgen, which is
#: how PRISM used to reach it -- but only as a side effect of launching a full
#: bilayer pack.  Running it directly costs ~2 s and makes the orientation
#: available to *every* route, including the ones that never invoke
#: packmol-memgen at all.
MEMEMBED_EXECUTABLE = "memembed"

#: Residue names of the pseudo-atoms that mark a fitted membrane's boundary
#: planes.  OPM and PPM write them, and MEMEMBED writes its own slab at
#: z = +-15 A the same way.  None of them may reach a force field, so every
#: orientation route strips them; :func:`strip_opm_dummy_atoms` is where that
#: happens for all three.
MEMBRANE_MARKER_RESNAMES = ("DUM", "DUMMY")

OPM_PDB_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb_id}.pdb"
OPM_ALT_URL = "https://opm.phar.umich.edu/pdb/{pdb_id}.pdb"
MAX_PEPTIDE_CN_DISTANCE_ANGSTROM = 2.0

#: Half-thickness of the lipid region PACKMOL-Memgen builds around z = 0.  A
#: default POPC-like patch occupies roughly z in [-23, +23] A including head
#: groups; only used to describe the target frame in diagnostics.
BILAYER_HALF_THICKNESS_A = 23.0
#: Half-width of the z window used to isolate the membrane-spanning part of a
#: structure.  Alignment is judged on this slab, never on the whole molecule:
#: real membrane proteins carry extramembrane domains (a GPCR ectodomain can be
#: 100 A tall) that would drag a whole-molecule centre far off the midplane and
#: raise a false alarm on a perfectly framed structure.  The value tracks the
#: hydrophobic core rather than the full lipid region, because soluble domains
#: pack against the head groups but cannot enter the core — widening it would
#: readmit exactly the mass this slab exists to exclude.
MEMBRANE_SLAB_HALF_WIDTH_A = 18.0
#: A solute whose total z extent is smaller than one bilayer cannot straddle the
#: midplane, so the two-sided density test below is meaningless for it.  Such
#: inputs (short peptides, single-residue test stubs) are judged on position
#: alone.
MIN_SPANNING_EXTENT_A = 30.0
#: Share of the slab atoms each leaflet side must carry before the structure
#: counts as straddling the midplane.  Loose enough to tolerate the strongly
#: asymmetric TM regions of real transporters, tight enough that a protein
#: merely grazing the bilayer with a loop is still rejected.
MIN_LEAFLET_ATOM_FRACTION = 0.15

#: Severity ranking of the notes :func:`orient_protein` returns, worst last.
#: Severity is carried as data rather than as a literal "WARNING:" inside the
#: message: the caller prefixes and colours notes itself, so a message that
#: prefixed itself was both duplicated ("WARNING: ... WARNING: ...") and
#: indistinguishable from a cosmetic remark once printed.
SEVERITY_ORDER = ("info", "warning", "error")
#: Chain breaks are listed individually so the user can find them in the input,
#: but a badly fragmented model must not bury the rest of the report.
MAX_LISTED_CHAIN_BREAKS = 5

#: Residue names counted as protein for the purpose of *measuring where the
#: solute is*.  That is this module's question, and it is not the same question
#: the builder's thermostat-group classifier asks, so the two sets are composed
#: separately from the shared vocabulary in :mod:`prism.membrane.config` rather
#: than being one flat list maintained twice.
#:
#: Wider than the classifier's protein set by exactly ``UNK``: an unidentified
#: polymer residue is still part of the thing being aligned, and leaving it out
#: biases the measurement.  The classifier must not make the same assumption --
#: it has to put the atom in a thermostat group and fails closed instead.
#:
#: Includes the modified residues :mod:`prism.ptm` produces (``SEP``/``TPO``/
#: ``PTR`` and friends).  Without them PRISM measured the alignment of its own
#: PTM output on an atom set that omitted the modified residues.
PROTEIN_RESIDUE_NAMES = config_vocab.SOLUTE_MEASUREMENT_RESIDUE_NAMES
#: Explicitly non-solute residues, used only when falling back to raw ATOM
#: records for a structure that uses none of the names above.  Water, bath ions,
#: structural metals and lipids all belong here: none of them should move the
#: measured centre of the solute.
NON_PROTEIN_RESIDUE_NAMES = config_vocab.NON_SOLUTE_RESIDUE_NAMES


def _read_protein_z_values(pdb_path: str) -> Tuple[List[float], str]:
    """Return the z coordinates of the protein atoms in ``pdb_path``.

    Selection is by residue name rather than by record type: modified residues
    are routinely written as HETATM, while packed systems write lipids and ions
    that must not dilute the measurement.  If a structure uses none of the known
    residue names (exotic or renamed force-field libraries), fall back to ATOM
    records that are demonstrably not solvent/ion/lipid — an unrecognised
    protein must still be measured rather than waved through unchecked.
    """
    named: List[float] = []
    fallback: List[float] = []
    try:
        with open(pdb_path, "r", errors="replace") as fh:
            for line in fh:
                # A short line cannot carry a z field; skipping it keeps a
                # truncated or hand-written file measurable instead of fatal.
                if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
                    continue
                try:
                    z = float(line[46:54])
                except ValueError:
                    continue
                if not math.isfinite(z):
                    continue
                resname = line[17:20].strip().upper()
                if resname in PROTEIN_RESIDUE_NAMES:
                    named.append(z)
                elif line.startswith("ATOM") and resname not in NON_PROTEIN_RESIDUE_NAMES:
                    fallback.append(z)
    except OSError:
        return [], "unreadable"

    if named:
        return named, "protein-residue-names"
    if fallback:
        return fallback, "unnamed-atom-records"
    return [], "none"


def membrane_alignment_report(
    pdb_path: str,
    *,
    tolerance: float = DEFAULT_CENTER_TOLERANCE_A,
    slab_half_width: float = MEMBRANE_SLAB_HALF_WIDTH_A,
) -> Dict[str, object]:
    """Measure how well ``pdb_path`` sits in the bilayer frame (z = 0 midplane).

    The verdict deliberately ignores extramembrane mass.  Atoms within
    ``slab_half_width`` of the midplane approximate the membrane-spanning
    region, and a structure is considered aligned when that region carries
    substantial density in *both* leaflets and its median z lies within
    ``tolerance`` of the midplane.  A whole-molecule centre-of-mass test would
    reject correctly framed proteins with large one-sided domains, which is the
    one false positive this function must not produce.

    Returns a dict with ``n_atoms``, ``min_z``/``max_z``/``median_z``,
    ``straddles_midplane`` and ``aligned``.  Coordinate-derived entries are
    ``None`` when nothing measurable was found (``measurable`` is then False);
    an unmeasurable file is reported, never silently declared aligned.
    """
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be a non-negative finite value")
    if not math.isfinite(slab_half_width) or slab_half_width <= 0:
        raise ValueError("slab_half_width must be a positive finite value")

    z_values, selection = _read_protein_z_values(pdb_path)
    report: Dict[str, object] = {
        "path": pdb_path,
        "n_atoms": len(z_values),
        "selection": selection,
        "tolerance": float(tolerance),
        "slab_half_width": float(slab_half_width),
        "bilayer_half_thickness": BILAYER_HALF_THICKNESS_A,
        "min_z": None,
        "max_z": None,
        "median_z": None,
        "span_z": None,
        "slab_median_z": None,
        "n_slab_atoms": 0,
        "n_below": 0,
        "n_above": 0,
        "straddles_midplane": False,
        "measurable": False,
        "aligned": False,
        "reason": "",
    }
    if not z_values:
        report["reason"] = (
            "no protein atoms could be read from the file"
            if selection != "unreadable"
            else "the file could not be read"
        )
        return report

    min_z, max_z = min(z_values), max(z_values)
    median_z = statistics.median(z_values)
    slab = [z for z in z_values if abs(z) <= slab_half_width]
    n_below = sum(1 for z in slab if z < 0.0)
    n_above = sum(1 for z in slab if z > 0.0)
    report.update(
        {
            "min_z": min_z,
            "max_z": max_z,
            "median_z": median_z,
            "span_z": max_z - min_z,
            "n_slab_atoms": len(slab),
            "n_below": n_below,
            "n_above": n_above,
            "measurable": True,
        }
    )
    if slab:
        report["slab_median_z"] = statistics.median(slab)

    if (max_z - min_z) < MIN_SPANNING_EXTENT_A:
        # Too thin to cross a bilayer: "straddles" is not a meaningful question,
        # so fall back to asking whether it sits on the midplane at all.  This
        # keeps short peptides and minimal test structures usable while still
        # catching one parked hundreds of angstroms away.
        report["straddles_midplane"] = min_z <= 0.0 <= max_z
        report["aligned"] = abs(median_z) <= tolerance
        span = max_z - min_z
        if report["aligned"]:
            report["reason"] = (
                f"structure spans only {span:.1f} A in z (too thin to cross a bilayer) but "
                f"sits on the midplane, median z = {median_z:+.1f} A"
            )
        else:
            report["reason"] = (
                f"structure sits {median_z:+.1f} A from the midplane, beyond the "
                f"{tolerance:.1f} A tolerance (and spans only {span:.1f} A in z, too thin "
                "to cross a bilayer)"
            )
        return report

    if not slab:
        report["reason"] = (
            f"no protein atom lies within {slab_half_width:.0f} A of the midplane, so the "
            "structure is entirely outside the bilayer PACKMOL-Memgen will build"
        )
        return report

    leaflet_fraction = min(n_below, n_above) / float(len(slab))
    slab_median = float(report["slab_median_z"])
    report["straddles_midplane"] = leaflet_fraction >= MIN_LEAFLET_ATOM_FRACTION
    if not report["straddles_midplane"]:
        side = "above" if n_above >= n_below else "below"
        report["reason"] = (
            f"the membrane-spanning region is one-sided: only "
            f"{leaflet_fraction * 100:.1f}% of the {len(slab)} atoms near the midplane lie "
            f"on the thin side, i.e. the structure sits {side} the bilayer rather than "
            "through it"
        )
        return report
    if abs(slab_median) > tolerance:
        report["reason"] = (
            f"the membrane-spanning region is offset by {slab_median:+.1f} A from the "
            f"midplane, beyond the {tolerance:.1f} A tolerance"
        )
        return report

    report["aligned"] = True
    report["reason"] = (
        f"membrane-spanning region centred at z = {slab_median:+.1f} A with "
        f"{n_below}/{n_above} atoms below/above the midplane"
    )
    return report


def _describe_offset(report: Dict[str, object]) -> str:
    """Phrase how far the structure would end up from the lipid region."""
    min_z, max_z = report.get("min_z"), report.get("max_z")
    if min_z is None or max_z is None:
        return "outside the lipid region"
    if min_z > BILAYER_HALF_THICKNESS_A:
        return f"{min_z - BILAYER_HALF_THICKNESS_A:.1f} A above the lipid region, in bulk water"
    if max_z < -BILAYER_HALF_THICKNESS_A:
        return f"{-BILAYER_HALF_THICKNESS_A - max_z:.1f} A below the lipid region, in bulk water"
    return "partly outside the lipid region"


def _measurement_line(report: Dict[str, object]) -> str:
    """State what was measured, so the verdict can be checked rather than believed."""
    return (
        f"  Measured over {report['n_atoms']} protein atoms: z from {report['min_z']:.1f} "
        f"to {report['max_z']:.1f} A, median {report['median_z']:.1f} A "
        f"({report['n_slab_atoms']} atom(s) within "
        f"{float(report['slab_half_width']):.0f} A of the midplane)."
    )


def _packing_frame_line(report: Dict[str, object]) -> str:
    """Explain why being outside the frame is fatal rather than cosmetic."""
    half = BILAYER_HALF_THICKNESS_A
    return (
        f"  PACKMOL-Memgen always builds the bilayer centred on z = 0 (lipid region "
        f"approximately z in [-{half:.0f}, +{half:.0f}] A) and keeps the solute at its "
        f"input coordinates, so this structure would be packed "
        f"{_describe_offset(report)} instead of embedded in it."
    )


def _orientation_alternatives(
    *, memembed_ready: bool = True, include_memembed: bool = True
) -> List[str]:
    """The ways a user can get a membrane-framed structure, as message lines.

    ``memembed_ready`` False replaces the "let PRISM do it" bullet with the
    install step it depends on: offering a route that is guaranteed to fail on
    this machine spends the user's next run on a second error message.
    ``include_memembed`` False drops that bullet entirely, for the one caller
    that is reporting a MEMEMBED run which already happened and did not work.
    """
    if not include_memembed:
        automatic: List[str] = []
    elif memembed_ready:
        automatic = [
            "    * let PRISM orient it automatically with MEMEMBED -- a few seconds, and it "
            "works with every bilayer route:",
            "        --membrane-orient memembed",
        ]
    else:
        automatic = [
            "    * install AmberTools so PRISM can orient it automatically with MEMEMBED "
            "(conda install -c conda-forge ambertools), then rerun -- a few seconds per "
            "structure, and it works with every bilayer route",
        ]
    return automatic + [
        "    * fetch a pre-oriented copy from OPM, if this structure is deposited there:",
        "        --membrane-orient opm --membrane-pdbid <pdb id>",
        "    * orient the structure externally with the PPM 3.0 server "
        "(https://opm.umich.edu/ppm_server) or download the matching OPM entry "
        "(https://opm.umich.edu/), then pass that file unchanged with "
        "--membrane-orient preoriented",
    ]


def _alignment_failure_message(
    report: Dict[str, object], *, memembed_ready: bool = True
) -> str:
    """Build an error a user can act on without reading the source."""
    return "\n".join(
        [
            f"Structure is not aligned to the membrane frame: {report['reason']}.",
            _measurement_line(report),
            _packing_frame_line(report),
            "  Fix it in one of these ways:",
        ]
        + _orientation_alternatives(memembed_ready=memembed_ready)
    )


def validate_membrane_alignment(
    pdb_path: str,
    tolerance: float = DEFAULT_CENTER_TOLERANCE_A,
) -> Dict[str, object]:
    """Raise :class:`ValueError` unless ``pdb_path`` is in the bilayer frame.

    Returns the :func:`membrane_alignment_report` dict on success.  A file with
    no readable protein atoms is *not* misalignment — the measurement simply
    could not be made, and raising here would mislabel unrelated input errors —
    so it is returned with ``measurable`` False for the caller to surface as a
    warning; the build fails later on its own, with a message about the real
    problem.
    """
    report = membrane_alignment_report(pdb_path, tolerance=tolerance)
    if report["measurable"] and not report["aligned"]:
        raise ValueError(_alignment_failure_message(report))
    return report


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


def memembed_available() -> bool:
    """True when the MEMEMBED binary can be run directly."""
    import shutil as _shutil

    return _shutil.which(MEMEMBED_EXECUTABLE) is not None


def run_memembed(
    input_pdb: str,
    output_pdb: str,
    *,
    threads: Optional[int] = None,
    timeout: int = 1800,
) -> int:
    """Orient *input_pdb* along the membrane normal with MEMEMBED.

    Writes the oriented solute to *output_pdb* with MEMEMBED's DUM slab markers
    removed, and returns how many of those markers were stripped.

    PRISM used to obtain this orientation only as a side effect of running
    packmol-memgen, which meant asking for ``orient=memembed`` forced the whole
    build down the packmol bilayer route.  That is backwards: the structures
    that need automatic orientation are exactly the ones not in OPM -- new
    depositions -- and they were therefore barred from the fast route.  MEMEMBED
    is a standalone binary and takes a couple of seconds, so PRISM runs it
    itself and hands the oriented coordinates to whichever route was asked for.
    """
    import shutil as _shutil
    import subprocess

    executable = _shutil.which(MEMEMBED_EXECUTABLE)
    if executable is None:
        raise RuntimeError(
            f"{MEMEMBED_EXECUTABLE} not found on PATH (it ships with AmberTools alongside "
            "packmol-memgen). Install AmberTools, or orient the structure externally and pass "
            "it with orient=preoriented."
        )

    input_pdb = os.path.abspath(input_pdb)
    output_pdb = os.path.abspath(output_pdb)
    # MEMEMBED writes beside its input unless told otherwise, and PRISM's input
    # may live in a directory the caller does not own. Keep the raw output next
    # to the destination so the DUM strip below reads a file we placed.
    raw_output = output_pdb + ".memembed"

    command = [executable, "-o", raw_output]
    if threads and int(threads) > 1:
        command += ["-a", str(int(threads))]
    command.append(input_pdb)

    try:
        proc = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"{MEMEMBED_EXECUTABLE} could not be run: {exc}") from exc

    if proc.returncode != 0 or not os.path.isfile(raw_output):
        detail = (proc.stderr or proc.stdout or "").strip().splitlines()
        tail = detail[-1] if detail else f"exit code {proc.returncode}"
        raise RuntimeError(f"{MEMEMBED_EXECUTABLE} failed to orient {os.path.basename(input_pdb)}: {tail}")

    try:
        removed = strip_opm_dummy_atoms(raw_output, output_pdb)
    finally:
        try:
            os.remove(raw_output)
        except OSError:
            pass
    return removed


def strip_opm_dummy_atoms(in_pdb: str, out_pdb: str) -> int:
    """Remove OPM/PPM/MEMEMBED DUMMY boundary atoms (membrane-plane markers).

    Returns the number of DUMMY atoms removed.
    """
    removed = 0
    with open(in_pdb, "r") as fh:
        lines = fh.readlines()
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            # OPM/PPM write the two membrane planes as pseudo-atoms named N and
            # O inside residue DUM/DUMMY.  The residue name alone identifies
            # them; an extra atom-name clause could only ever match a subset of
            # the same records, so it never changed the outcome.
            resname = line[17:20].strip().upper()
            if resname in MEMBRANE_MARKER_RESNAMES:
                removed += 1
                continue
        out.append(line)
    with open(out_pdb, "w") as fh:
        fh.writelines(out)
    return removed


def _residue_records(lines: List[str]) -> List[Dict[str, object]]:
    """Group the coordinate lines of a PDB into residues, in file order."""
    residues: List[Dict[str, object]] = []
    residue_by_key: Dict[tuple, Dict[str, object]] = {}
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
                "resname": line[17:20].strip().upper(),
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
    return residues


def _scan_chain_breaks(
    lines: List[str], max_cn_distance: float
) -> List[Dict[str, object]]:
    """Return one record per peptide junction that must be split by a ``TER``."""
    residues = _residue_records(lines)
    breaks: List[Dict[str, object]] = []
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
            # The two new termini end up where the backbone was cut, so record
            # that z: whether the charged caps tleap adds are harmless depends
            # entirely on it (see :func:`_buried_chain_break_note`).
            anchors = [
                coordinate for coordinate in (c_coord, n_coord) if coordinate is not None
            ] or list(previous["atoms"].values()) + list(current["atoms"].values())
            breaks.append(
                {
                    "line_index": current["first_index"],
                    "chain": previous_chain,
                    "prev_resname": previous["resname"],
                    "prev_resid": previous_resid,
                    "next_chain": current_chain,
                    "next_resname": current["resname"],
                    "next_resid": current_resid,
                    "gap_distance": (
                        math.dist(c_coord, n_coord) if geometry_available else None
                    ),
                    "midpoint_z": (
                        statistics.fmean(a[2] for a in anchors) if anchors else None
                    ),
                }
            )
    return breaks


def find_chain_breaks(
    in_pdb: str,
    max_cn_distance: float = MAX_PEPTIDE_CN_DISTANCE_ANGSTROM,
) -> List[Dict[str, object]]:
    """Locate the junctions :func:`insert_ter_at_chain_breaks` would split.

    Each record carries ``chain``, the flanking ``prev_resname``/``prev_resid``
    and ``next_resname``/``next_resid``, the measured ``gap_distance`` (C to N,
    ``None`` when those backbone atoms are absent) and ``midpoint_z``.  Callers
    need the identities and the z position, not just a count: a break is a
    decision about where two charged chain termini will be created, and only
    the caller knows whether that place is inside the bilayer.

    Run this *before* inserting the ``TER`` records — afterwards the same scan
    honours those records and reports nothing.
    """
    if not math.isfinite(max_cn_distance) or max_cn_distance <= 0:
        raise ValueError("max_cn_distance must be a positive finite value")
    with open(in_pdb, "r") as fh:
        lines = fh.readlines()
    return _scan_chain_breaks(lines, max_cn_distance)


def _describe_chain_break(record: Dict[str, object]) -> str:
    """Name one break so the user can find it in their own input file."""
    chain = str(record["chain"] or "").strip() or "?"
    gap = record["gap_distance"]
    gap_text = f"C-N gap {gap:.1f} A" if gap is not None else "backbone gap unmeasurable"
    z = record["midpoint_z"]
    z_text = f"z = {z:+.1f} A" if z is not None else "z unknown"
    return (
        f"{chain}:{record['prev_resname']}{record['prev_resid']}"
        f"-{record['next_resname']}{record['next_resid']} at {z_text} ({gap_text})"
    )


def _buried_chain_break_note(
    breaks: List[Dict[str, object]],
    *,
    slab_half_width: float = MEMBRANE_SLAB_HALF_WIDTH_A,
) -> str:
    """Report the breaks that land inside the hydrophobic core, or "".

    tleap builds every ``TER``-separated fragment as its own chain and caps it
    with the charged ff14SB terminal templates (NH3+ / COO-).  In an
    extramembrane loop that is a well-understood approximation; in the middle of
    a discontinuous TM helix it buries a charge pair in the lipid core, which
    distorts the local structure far more than the stretched bond the ``TER``
    exists to prevent.  The count alone cannot distinguish the two cases, so the
    z position of each break is what gets reported.

    Only breaks *within* a chain qualify: a chain boundary is a real terminus of
    the model, charged with or without our ``TER``, so reporting it as an
    introduced artifact would cry wolf on every TM complex.
    """
    buried = [
        record
        for record in breaks
        if record["chain"] == record["next_chain"]
        and record["midpoint_z"] is not None
        and abs(float(record["midpoint_z"])) < slab_half_width
    ]
    if not buried:
        return ""
    listed = [_describe_chain_break(record) for record in buried[:MAX_LISTED_CHAIN_BREAKS]]
    if len(buried) > MAX_LISTED_CHAIN_BREAKS:
        listed.append(f"and {len(buried) - MAX_LISTED_CHAIN_BREAKS} more")
    return (
        f"{len(buried)} chain break(s) lie inside the membrane slab "
        f"(|z| < {slab_half_width:.0f} A): " + "; ".join(listed) + ". "
        "tleap caps each chain-break fragment with charged termini (NH3+/COO-), so these "
        "bury a charge pair in the hydrophobic core. Model the missing loop(s) and rerun, "
        "or continue knowing the transmembrane region carries an artificial charge pair."
    )


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

    Returns the number of records inserted; use :func:`find_chain_breaks` for
    the identity and position of each break.
    """
    if not math.isfinite(max_cn_distance) or max_cn_distance <= 0:
        raise ValueError("max_cn_distance must be a positive finite value")

    with open(in_pdb, "r") as fh:
        lines = fh.readlines()

    insertion_indices = {
        record["line_index"] for record in _scan_chain_breaks(lines, max_cn_distance)
    }
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


def _note_payload(entries: List[Tuple[str, str]]) -> Dict[str, object]:
    """Collapse ``(text, severity)`` pairs into the reported note fields.

    ``note`` keeps the single flat string older callers already print;
    ``notes``/``severity`` let a caller that cares render "this structure is not
    in the bilayer at all" differently from "a TER was inserted".
    """
    notes = [
        {"text": text, "severity": severity} for text, severity in entries if text
    ]
    severity = "info"
    for item in notes:
        if SEVERITY_ORDER.index(item["severity"]) > SEVERITY_ORDER.index(severity):
            severity = item["severity"]
    return {
        "note": " ".join(str(item["text"]) for item in notes),
        "notes": notes,
        "severity": severity,
    }


def _auto_missing_memembed_message(report: Dict[str, object]) -> str:
    """Explain that ``auto`` diagnosed the problem but cannot fix it here."""
    return "\n".join(
        [
            f"Structure is not aligned to the membrane frame: {report['reason']}.",
            f"  orient=auto measured this and would have oriented it with MEMEMBED, but the "
            f"'{MEMEMBED_EXECUTABLE}' binary is not on PATH (it ships with AmberTools "
            f"alongside packmol-memgen).",
            _measurement_line(report),
            _packing_frame_line(report),
            "  Fix it in one of these ways:",
        ]
        + _orientation_alternatives(memembed_ready=False)
    )


def _auto_unorientable_message(
    before: Dict[str, object], after: Dict[str, object]
) -> str:
    """Report that MEMEMBED ran and the result still is not membrane-framed.

    The generic :func:`_alignment_failure_message` is wrong here, and wrong in
    the expensive direction: its first suggestion is to run MEMEMBED, which is
    what just failed.  MEMEMBED orients a protein by fitting a slab to its
    hydrophobic belt, so the single most likely reason it produced nothing
    usable is that there is no belt to fit — i.e. the structure has no
    transmembrane region and does not belong in a bilayer at all.  That
    possibility is named first, because every fix below it is wasted effort for
    the user who pointed ``--membrane`` at a soluble protein.
    """
    return "\n".join(
        [
            "MEMEMBED could not put this structure in the membrane frame, and orient=auto "
            "has no route left to try.",
            f"  Input was measured as un-framed ({before['reason']}), so MEMEMBED was run; "
            f"its output is still un-framed: {after['reason']}.",
            f"  After orientation: z from {after['min_z']:.1f} to {after['max_z']:.1f} A, "
            f"{after['n_below']}/{after['n_above']} atoms below/above the midplane among the "
            f"{after['n_slab_atoms']} within {float(after['slab_half_width']):.0f} A of it.",
            "  The likeliest explanation is that this protein has no transmembrane region: "
            "MEMEMBED orients by fitting a membrane slab to a hydrophobic belt, and a soluble "
            "protein has none, so the fit it returns is meaningless. Check that this is the "
            "membrane protein you meant to build -- if it is soluble, drop --membrane and "
            "build it in a water box instead.",
            "  If it really is a membrane protein (a monotopic or peripheral one, which "
            "MEMEMBED fits poorly), orient it externally and pass the result unchanged:",
        ]
        + _orientation_alternatives(include_memembed=False)
    )


def _resolve_auto_orientation(
    input_pdb: str, *, tolerance: float
) -> Tuple[str, List[Tuple[str, str]], Dict[str, object]]:
    """Decide which route ``auto`` takes for *input_pdb*.

    Returns ``(method, notes, report)``: the concrete method to run, the
    ``(text, severity)`` entries explaining the choice, and the measurement it
    was made from.

    The decision is a measurement, not a guess, and it is made *before* any
    coordinate is written: a structure whose membrane-spanning region already
    straddles z = 0 is kept exactly as it is (running MEMEMBED on it would cost
    ~5 s and perturb a frame the user may have prepared deliberately), and
    anything else is handed to MEMEMBED.  An input that cannot be measured is
    passed through rather than oriented: MEMEMBED's verdict on a file PRISM
    could not read would be a membrane fit nobody has checked, and the real
    problem — an unreadable or protein-free file — surfaces downstream with a
    message about itself.

    This runs regardless of ``validate``: it chooses the route rather than
    approving the output, and turning validation off must not turn ``auto``
    into a coin flip.
    """
    report = membrane_alignment_report(input_pdb, tolerance=tolerance)
    if not report["measurable"]:
        return (
            "preoriented",
            [(
                f"orient=auto could not measure the input against the membrane frame "
                f"({report['reason']}), so it was passed through unchanged and MEMEMBED was "
                f"not run; orienting a file PRISM cannot read would only produce an "
                f"unverifiable membrane fit.",
                "warning",
            )],
            report,
        )
    if report["aligned"]:
        # The closing caveat is not padding.  This branch is the only report a
        # user gets when they point --membrane at a soluble protein that happens
        # to straddle z = 0, and "already in the membrane frame" read as "PRISM
        # confirmed this is a membrane protein" is precisely the misreading that
        # ends in a bilayer packed around a globular protein.  Position is what
        # was measured, so position is what is claimed.
        return (
            "preoriented",
            [(
                f"orient=auto kept the input coordinates unchanged: the structure already "
                f"sits in the bilayer frame ({report['reason']}; measured over "
                f"{report['n_atoms']} protein atoms, z from {report['min_z']:.1f} to "
                f"{report['max_z']:.1f} A), so MEMEMBED was not run. This checks where the "
                f"structure sits, not whether it is a membrane protein -- a soluble protein "
                f"centred on z = 0 passes it too.",
                "info",
            )],
            report,
        )
    if not memembed_available():
        # Fail closed, and say what was measured: 'auto' must never fall back to
        # treating an un-framed structure as preoriented just because the tool
        # that would have fixed it is missing.
        raise RuntimeError(_auto_missing_memembed_message(report))
    return (
        "memembed",
        [(
            f"orient=auto moved the input into the membrane frame with MEMEMBED: as given it "
            f"was not membrane-framed ({report['reason']}; measured over {report['n_atoms']} "
            f"protein atoms, z from {report['min_z']:.1f} to {report['max_z']:.1f} A, median "
            f"{report['median_z']:.1f} A). Every coordinate PRISM writes from here on is "
            f"MEMEMBED's, not the input's.",
            "info",
        )],
        report,
    )


def orient_protein(
    input_pdb: str,
    output_pdb: str,
    method: str,
    pdb_id: Optional[str] = None,
    *,
    validate: bool = True,
    tolerance: float = DEFAULT_CENTER_TOLERANCE_A,
) -> dict:
    """Orient a protein along the membrane normal.

    ``method`` is ``"auto"``, ``"preoriented"``, ``"memembed"``, ``"opm"`` or
    ``"ppm"``.  ``auto`` is not a backend of its own: it measures the input
    against the bilayer frame and dispatches to ``preoriented`` (already framed
    — coordinates untouched) or ``memembed`` (anything else), which is what
    makes a bare ``--membrane`` correct for an RCSB download and for an
    OPM/PPM-prepared file alike.  See :func:`_resolve_auto_orientation`.

    Returns ``{"oriented_pdb", "method", "requested_method", "note", "notes",
    "severity"}``, where ``method`` is the route that actually ran (never
    ``"auto"``), ``notes`` itemises the flat ``note`` string as
    ``{"text", "severity"}`` and ``severity`` is the worst of them
    (``"info"``/``"warning"``/``"error"``).  Automated orientation requests fail
    closed: an un-oriented input must never be relabelled as preoriented merely
    because a network service or backend is unavailable.

    With ``validate`` (default), the prepared structure is checked against the
    bilayer frame and rejected if it does not belong there — neither trusting
    the input nor running an orientation tool may be a silent promise.  ``opm``
    output is framed by construction and only annotated.  Set ``validate=False``
    to downgrade the check to a note.
    """
    import shutil

    method = (method or "preoriented").lower()
    note = ""

    if method not in {"auto", "opm", "ppm", "preoriented", "memembed"}:
        raise ValueError(f"Unsupported membrane orientation method: {method!r}")
    if method == "opm" and not pdb_id:
        raise ValueError("OPM orientation requires a PDB id.")
    # Reject a nonsensical tolerance before any file is written, so a caller
    # cannot end up with output produced under an unchecked criterion.
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be a non-negative finite value")

    requested_method = method
    auto_entries: List[Tuple[str, str]] = []
    auto_input_report: Optional[Dict[str, object]] = None
    if method == "auto":
        method, auto_entries, auto_input_report = _resolve_auto_orientation(
            input_pdb, tolerance=tolerance
        )

    if method == "opm":
        fetched = fetch_opm(pdb_id, output_pdb)
        if fetched:
            n = strip_opm_dummy_atoms(output_pdb, output_pdb)
            # Measure the junctions before splitting them; the same scan reports
            # nothing once the TER records are in place.
            chain_breaks = find_chain_breaks(output_pdb)
            breaks = insert_ter_at_chain_breaks(output_pdb, output_pdb)
            entries = [
                (
                    f"Fetched oriented structure {pdb_id} from OPM "
                    f"(removed {n} DUMMY atoms; inserted {breaks} TER chain-break records).",
                    "info",
                ),
                # OPM coordinates are in the bilayer frame by construction, so
                # the break z values can be read against the membrane directly.
                (_buried_chain_break_note(chain_breaks), "warning"),
            ]
            # A disagreement with our own check means the measurement is off,
            # not the structure.  Record it instead of blocking a good input.
            if validate:
                try:
                    report = membrane_alignment_report(output_pdb, tolerance=tolerance)
                except (OSError, ValueError):
                    report = None
                if report is not None and report["measurable"] and not report["aligned"]:
                    entries.append((
                        "Note: the OPM coordinates did not pass PRISM's membrane-frame check "
                        f"({report['reason']}); proceeding because OPM entries are oriented at source.",
                        "warning",
                    ))
            return {
                "oriented_pdb": output_pdb,
                "method": "opm",
                "requested_method": requested_method,
                **_note_payload(entries),
            }
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
        # Run MEMEMBED here rather than letting packmol-memgen do it internally.
        # Delegating tied the orientation to the packmol bilayer route, which
        # locked the fast route out of exactly the structures that need
        # automatic orientation (anything not in OPM).
        removed = run_memembed(input_pdb, output_pdb)
        note = f"Oriented with MEMEMBED (removed {removed} DUM slab marker(s))."
        copied = True
    else:
        # Explicitly preoriented input: trusted as-is, and checked below.
        copied = os.path.abspath(input_pdb) != os.path.abspath(output_pdb)
        if copied:
            shutil.copy2(input_pdb, output_pdb)
    # Measure the junctions before splitting them; the same scan reports nothing
    # once the TER records are in place.
    chain_breaks = find_chain_breaks(output_pdb)
    breaks = insert_ter_at_chain_breaks(output_pdb, output_pdb)
    # The dispatch note comes first: it is the reason the rest of this report
    # describes the coordinates it does.  It also replaces the generic
    # "treated as preoriented" line, which would otherwise say the same thing
    # with less information.
    entries: List[Tuple[str, str]] = list(auto_entries)
    if note:
        entries.append((note, "info"))
    elif not auto_entries:
        entries.append(("Input treated as preoriented.", "info"))
    if breaks:
        entries.append((f"Inserted {breaks} TER chain-break record(s).", "info"))

    # MEMEMBED produces genuinely oriented coordinates, so its output is checked
    # like any other: a fit that failed to find a membrane is exactly what this
    # catches, and it used to be invisible because the orientation happened out
    # of sight inside packmol-memgen.
    auto_ran_memembed = auto_input_report is not None and method == "memembed"
    if validate:
        try:
            report = validate_membrane_alignment(output_pdb, tolerance)
        except ValueError:
            # Measure before deleting: after MEMEMBED has already run, the
            # generic failure message opens by suggesting MEMEMBED, so 'auto'
            # needs the numbers to say what really happened instead.
            replacement = (
                _auto_unorientable_message(
                    auto_input_report,
                    membrane_alignment_report(output_pdb, tolerance=tolerance),
                )
                if auto_ran_memembed
                else None
            )
            # Fail closed: do not leave a copy that looks like a prepared,
            # membrane-ready structure.  An in-place run touches the
            # caller's own file, so that one is left alone.
            if copied:
                try:
                    os.remove(output_pdb)
                except OSError:
                    pass
            if replacement is not None:
                # 'from None': the suppressed message's first suggestion is to
                # run MEMEMBED, which is the step that just failed here.
                raise ValueError(replacement) from None
            raise
    else:
        report = membrane_alignment_report(output_pdb, tolerance=tolerance)
    if not report["measurable"]:
        entries.append(
            (f"Membrane alignment NOT verified: {report['reason']}.", "warning")
        )
    elif not report["aligned"]:
        # Only reachable with validate=False.  This is the worst state the
        # module can report — the bilayer will be built around empty space —
        # so it is flagged as an error rather than left to look like the
        # cosmetic notes it is printed next to.
        entries.append((
            f"Membrane alignment NOT verified (validation disabled) and the structure "
            f"is measurably outside the bilayer frame "
            f"({report['reason']}); z {report['min_z']:.1f} to {report['max_z']:.1f} A "
            f"against a bilayer at z in "
            f"[-{BILAYER_HALF_THICKNESS_A:.0f}, +{BILAYER_HALF_THICKNESS_A:.0f}] A. "
            f"PACKMOL-Memgen keeps the solute at these coordinates, so it would be "
            f"packed {_describe_offset(report)} instead of embedded in the bilayer, and "
            f"the resulting system answers no membrane question. "
            + (
                # Suggesting MEMEMBED here would be advice to repeat the step
                # that produced these coordinates.
                "MEMEMBED already ran on this structure (orient=auto) and did not find a "
                "membrane in it, so it may have no transmembrane region at all: check that "
                "--membrane is the right mode for this protein, and if it is, orient it "
                "externally (--membrane-orient opm/preoriented)."
                if auto_ran_memembed
                else "Re-enable validation, or orient the structure "
                "(--membrane-orient opm/memembed)."
            ),
            "error",
        ))
    else:
        entries.append(
            (f"Membrane alignment verified ({report['reason']}).", "info")
        )
        # Only meaningful once the structure is known to be in the bilayer
        # frame; against an arbitrary frame the break z values say nothing.
        entries.append((_buried_chain_break_note(chain_breaks), "warning"))
    return {
        "oriented_pdb": output_pdb,
        # The route that ran, so a caller can tell whether the coordinates it
        # got back are the ones it handed in; ``requested_method`` keeps the
        # "auto" ask visible next to that answer.
        "method": method if method == "memembed" else "preoriented",
        "requested_method": requested_method,
        **_note_payload(entries),
    }
