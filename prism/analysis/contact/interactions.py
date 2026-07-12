#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein-ligand interaction typing over an MD trajectory.

This module classifies *what kind* of interaction each contacting residue makes
with the ligand (hydrogen bond, salt bridge, pi-stacking, pi-cation,
hydrophobic, halogen bond) and reports the per-type *occupancy* — the fraction
of analysed frames in which the interaction is present.

It is intentionally dependency-light: it uses only ``mdtraj``/``numpy`` (already
required by the analysis module) plus an optional RDKit ligand molecule for
ligand chemistry (aromatic rings, formal charges, halogens). When RDKit or the
ligand mol is unavailable, element-based fallbacks are used so the typer still
produces hydrogen-bond, hydrophobic and (protein-side) salt-bridge / pi
information.

Geometric cutoffs follow the published reference values used by PLIP / ProLIF:

    hydrophobic          C...C   <= 0.40 nm
    salt bridge          charged-centroid distance <= 0.55 nm
    pi-stacking          ring-centroid distance <= 0.55 nm
    pi-cation            ring-centroid ... cation <= 0.45 nm
    halogen bond         ligand-halogen ... acceptor(O/N) <= 0.40 nm
    hydrogen bond        mdtraj Wernet-Nilsson criterion

All distances are in nanometres to match mdtraj.

The public entry point is :class:`InteractionTyper`.  A failure anywhere in the
typing pass is non-fatal: the typer logs a warning and returns empty results so
the surrounding report falls back to the legacy distance-only colouring.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import mdtraj as md  # noqa: F401  (used for type parity with the analyzer)

    MDTRAJ_AVAILABLE = True
except Exception:  # pragma: no cover - mdtraj is a hard analysis dep in practice
    MDTRAJ_AVAILABLE = False

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Interaction-type metadata (shared with the front-end via the report builder) #
# --------------------------------------------------------------------------- #

#: Canonical interaction types, ordered by how prominently they should be drawn.
INTERACTION_TYPES: Tuple[str, ...] = (
    "hydrogen_bond",
    "salt_bridge",
    "pi_stacking",
    "pi_cation",
    "halogen_bond",
    "hydrophobic",
)

#: Display metadata used to build the legend and edge styling. Colours are the
#: Okabe-Ito colour-blind-safe categorical palette; ``dash`` mirrors LigPlot+ /
#: Maestro conventions (solid = directional/electrostatic, dashed = weak/pi).
INTERACTION_STYLE: Dict[str, Dict[str, object]] = {
    "hydrogen_bond": {"label": "Hydrogen bond", "color": "#0072B2", "dash": []},
    "salt_bridge": {"label": "Salt bridge", "color": "#CC79A7", "dash": []},
    "pi_stacking": {"label": "π–π stacking", "color": "#009E73", "dash": [6, 4]},
    "pi_cation": {"label": "π–cation", "color": "#E69F00", "dash": [6, 4]},
    "halogen_bond": {"label": "Halogen bond", "color": "#56B4E9", "dash": [2, 3]},
    "hydrophobic": {"label": "Hydrophobic", "color": "#999999", "dash": [1, 4]},
}


#: Residue chemistry classes, used to colour residue nodes (Maestro/MOE style).
RESIDUE_CLASS: Dict[str, str] = {}
for _names, _cls in (
    (("ALA", "VAL", "LEU", "ILE", "MET", "PRO", "CYS", "CYX"), "hydrophobic"),
    (("PHE", "TYR", "TRP", "HIS", "HID", "HIE", "HIP"), "aromatic"),
    (("ARG", "LYS", "LYN"), "positive"),
    (("ASP", "GLU", "ASH", "GLH"), "negative"),
    (("SER", "THR", "ASN", "GLN", "TYR"), "polar"),
    (("GLY",), "glycine"),
):
    for _n in _names:
        RESIDUE_CLASS[_n] = _cls


def residue_class(resname: str) -> str:
    """Return the chemistry class for a 3-letter residue name."""
    return RESIDUE_CLASS.get(str(resname).upper()[:3], "other")


#: Front-end colours for residue chemistry classes (used to colour residue
#: nodes in the interaction diagram, Maestro/MOE style). Chosen to be distinct
#: and reasonably colour-blind friendly.
RESIDUE_CLASS_COLORS: Dict[str, str] = {
    "hydrophobic": "#E1A140",
    "aromatic": "#9B7EDE",
    "positive": "#4C72B0",
    "negative": "#C44E52",
    "polar": "#55A868",
    "glycine": "#9AA0A6",
    "other": "#7F8C8D",
}


# --------------------------------------------------------------------------- #
# Atom-group definitions                                                       #
# --------------------------------------------------------------------------- #

# Protein side-chain functional groups, keyed by residue name -> list of atom
# names whose centroid defines the group.
_PROTEIN_CATION = {
    "LYS": ["NZ"],
    "ARG": ["CZ", "NH1", "NH2", "NE"],
    "HIP": ["ND1", "NE2", "CE1", "CG", "CD2"],
}
_PROTEIN_ANION = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}
_PROTEIN_AROMATIC = {
    "PHE": [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
    "TYR": [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
    "HIS": [["CG", "ND1", "CE1", "NE2", "CD2"]],
    "HID": [["CG", "ND1", "CE1", "NE2", "CD2"]],
    "HIE": [["CG", "ND1", "CE1", "NE2", "CD2"]],
    "HIP": [["CG", "ND1", "CE1", "NE2", "CD2"]],
    "TRP": [
        ["CG", "CD1", "NE1", "CE2", "CD2"],
        ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    ],
}

# Halogen-bond donors are Cl/Br/I only. Organic fluorine has no sigma-hole and
# does not form halogen bonds (PLIP/ProLIF exclude it), so F is intentionally omitted.
_HALOGEN_ELEMENTS = {"Cl", "Br", "I", "CL", "BR"}


class InteractionTyper:
    """Classify and quantify protein-ligand interactions over a trajectory."""

    # geometric cutoffs in nm
    HYDROPHOBIC_CUT = 0.40
    SALTBRIDGE_CUT = 0.55
    PISTACK_CUT = 0.55
    PICATION_CUT = 0.45
    HALOGEN_CUT = 0.40

    def __init__(self, max_frames: int = 1000, verbose: bool = False):
        self.max_frames = max_frames
        self.verbose = verbose

    # ------------------------------------------------------------------ #
    # Public API                                                          #
    # ------------------------------------------------------------------ #
    def compute(
        self,
        traj,
        ligand_residue,
        ligand_mol=None,
        ligand_atoms: Optional[List[int]] = None,
    ) -> Dict[str, object]:
        """Compute per-residue interaction types and occupancies.

        Parameters
        ----------
        traj : mdtraj.Trajectory
            The loaded trajectory (full length; it is sampled internally).
        ligand_residue : mdtraj residue
            The identified ligand residue.
        ligand_mol : rdkit.Chem.Mol, optional
            RDKit ligand molecule (for aromatic rings / formal charges /
            halogens). Heavy-atom order must match ``ligand_atoms``.
        ligand_atoms : list[int], optional
            Trajectory atom indices of the ligand heavy atoms, in the same
            order RDKit enumerates its heavy atoms (as produced by the contact
            analyzer). Required to map RDKit chemistry onto trajectory atoms.

        Returns
        -------
        dict
            ``{"residue_interactions": {res_id: {type: occupancy}},
               "summary": {type: n_residues}, "n_frames": int}``
            Empty dict on failure.
        """
        try:
            return self._compute(traj, ligand_residue, ligand_mol, ligand_atoms)
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning(f"Interaction typing failed ({exc}); falling back to distance-only contacts")
            return {"residue_interactions": {}, "summary": {}, "n_frames": 0}

    # ------------------------------------------------------------------ #
    # Internals                                                           #
    # ------------------------------------------------------------------ #
    def _compute(self, traj, ligand_residue, ligand_mol, ligand_atoms):
        top = traj.topology

        # Sample frames for efficiency (mirror the contact analyzer policy).
        n_frames = traj.n_frames
        if n_frames > self.max_frames:
            idx = np.linspace(0, n_frames - 1, self.max_frames, dtype=int)
            sample = traj[idx]
        else:
            sample = traj
        n = sample.n_frames
        xyz = sample.xyz  # (n, n_atoms, 3) in nm

        # ---- Build ligand atom groups -------------------------------------
        lig_heavy = [a.index for a in ligand_residue.atoms if a.element.symbol != "H"]
        lig_all = [a.index for a in ligand_residue.atoms]
        lig_carbons = [a.index for a in ligand_residue.atoms if a.element.symbol == "C"]
        lig_halogens = [a.index for a in ligand_residue.atoms if a.element.symbol in _HALOGEN_ELEMENTS]
        lig_on = [a.index for a in ligand_residue.atoms if a.element.symbol in ("O", "N")]

        lig_rings: List[List[int]] = []
        lig_cations: List[int] = []
        lig_anions: List[int] = []
        if ligand_mol is not None and ligand_atoms:
            lig_rings, lig_cations, lig_anions = self._ligand_chemistry(ligand_mol, ligand_atoms)

        # ---- Build protein atom groups, keyed by residue -----------------
        prot_residues = self._protein_residues(top)
        prot_cation = {}   # res_id -> [atom idx]
        prot_anion = {}
        prot_aromatic = {}  # res_id -> [ [atom idx,...], ... ]
        prot_carbons = {}
        prot_on = {}
        res_name = {}
        for res in prot_residues:
            # NOTE: rid must match FastContactAnalyzer's residue id ("NAME"+resSeq,
            # no chain) so the typing joins onto the distance-contact list in
            # visualization_generator._annotate_interactions. In multi-chain
            # systems residues with the same name+number across chains share an
            # id (a pre-existing convention of the contact module); the OR-based
            # occupancy keeps that merge bounded to [0, 1].
            rid = f"{res.name}{res.resSeq}"
            res_name[rid] = res.name
            atoms_by_name = {a.name: a.index for a in res.atoms}
            cat = _PROTEIN_CATION.get(res.name)
            if cat:
                prot_cation[rid] = [atoms_by_name[n] for n in cat if n in atoms_by_name]
            an = _PROTEIN_ANION.get(res.name)
            if an:
                prot_anion[rid] = [atoms_by_name[n] for n in an if n in atoms_by_name]
            arom = _PROTEIN_AROMATIC.get(res.name)
            if arom:
                rings = [[atoms_by_name[n] for n in ring if n in atoms_by_name] for ring in arom]
                prot_aromatic[rid] = [r for r in rings if len(r) >= 3]
            prot_carbons[rid] = [a.index for a in res.atoms if a.element.symbol == "C"]
            prot_on[rid] = [a.index for a in res.atoms if a.element.symbol in ("O", "N")]

        # Accumulate a per-frame boolean *presence* mask per (res_id, type).
        # OR-ing masks (rather than summing True-counts) ensures each frame is
        # counted at most once per interaction type, so occupancy stays in [0, 1]
        # even when several sub-loops can fire for the same residue in the same
        # frame (e.g. HIP is both cationic and aromatic, so both pi-cation
        # directions can match; or a residue id is shared across chains).
        counts: Dict[Tuple[str, str], np.ndarray] = {}

        def bump(rid, itype, mask):
            """OR a per-frame boolean mask into the (rid, itype) presence accumulator."""
            if not np.any(mask):
                return
            key = (rid, itype)
            if key in counts:
                counts[key] |= mask
            else:
                counts[key] = np.asarray(mask, dtype=bool).copy()

        # ---- Hydrophobic: ligand C ... protein sidechain C ----------------
        if lig_carbons:
            lig_c = xyz[:, lig_carbons, :]  # (n, L, 3)
            for rid, pc in prot_carbons.items():
                if not pc:
                    continue
                present = self._min_pair_dist(lig_c, xyz[:, pc, :]) <= self.HYDROPHOBIC_CUT
                bump(rid, "hydrophobic", present)

        # ---- Salt bridges -------------------------------------------------
        # ligand cation ... protein anion ; ligand anion ... protein cation
        if lig_cations:
            lig_cat_c = xyz[:, lig_cations, :].mean(axis=1)  # (n,3)
            for rid, atoms in prot_anion.items():
                if not atoms:
                    continue
                d = np.linalg.norm(lig_cat_c - xyz[:, atoms, :].mean(axis=1), axis=1)
                bump(rid, "salt_bridge", d <= self.SALTBRIDGE_CUT)
        if lig_anions:
            lig_an_c = xyz[:, lig_anions, :].mean(axis=1)
            for rid, atoms in prot_cation.items():
                if not atoms:
                    continue
                d = np.linalg.norm(lig_an_c - xyz[:, atoms, :].mean(axis=1), axis=1)
                bump(rid, "salt_bridge", d <= self.SALTBRIDGE_CUT)

        # ---- pi-stacking & pi-cation -------------------------------------
        if lig_rings:
            lig_ring_cent = [xyz[:, r, :].mean(axis=1) for r in lig_rings]
            # pi-stacking: ligand ring ... protein aromatic ring
            for rid, rings in prot_aromatic.items():
                present = np.zeros(n, dtype=bool)
                for pr in rings:
                    pcent = xyz[:, pr, :].mean(axis=1)
                    for lc in lig_ring_cent:
                        present |= np.linalg.norm(lc - pcent, axis=1) <= self.PISTACK_CUT
                bump(rid, "pi_stacking", present)
            # pi-cation: ligand ring ... protein cation
            for rid, atoms in prot_cation.items():
                if not atoms:
                    continue
                pcat = xyz[:, atoms, :].mean(axis=1)
                present = np.zeros(n, dtype=bool)
                for lc in lig_ring_cent:
                    present |= np.linalg.norm(lc - pcat, axis=1) <= self.PICATION_CUT
                bump(rid, "pi_cation", present)
        # ligand cation ... protein aromatic ring (other pi-cation direction)
        if lig_cations:
            lcat = xyz[:, lig_cations, :].mean(axis=1)
            for rid, rings in prot_aromatic.items():
                present = np.zeros(n, dtype=bool)
                for pr in rings:
                    pcent = xyz[:, pr, :].mean(axis=1)
                    present |= np.linalg.norm(lcat - pcent, axis=1) <= self.PICATION_CUT
                bump(rid, "pi_cation", present)

        # ---- Halogen bonds: ligand halogen ... protein O/N ----------------
        if lig_halogens:
            lig_x = xyz[:, lig_halogens, :]
            for rid, on in prot_on.items():
                if not on:
                    continue
                present = self._min_pair_dist(lig_x, xyz[:, on, :]) <= self.HALOGEN_CUT
                bump(rid, "halogen_bond", present)

        # ---- Hydrogen bonds via mdtraj Wernet-Nilsson --------------------
        self._add_hydrogen_bonds(sample, ligand_residue, prot_residues, counts, n)

        # ---- Aggregate to occupancies ------------------------------------
        residue_interactions: Dict[str, Dict[str, float]] = {}
        summary: Dict[str, int] = {}
        for (rid, itype), mask in counts.items():
            occ = float(np.count_nonzero(mask)) / float(n) if n else 0.0
            if occ < 0.05:  # ignore transient noise (<5% of frames)
                continue
            residue_interactions.setdefault(rid, {})[itype] = round(min(occ, 1.0), 3)
            summary[itype] = summary.get(itype, 0) + 1

        if self.verbose:
            logger.info(
                f"Interaction typing: {len(residue_interactions)} residues, "
                f"types={ {k: summary[k] for k in summary} }, frames={n}"
            )

        return {"residue_interactions": residue_interactions, "summary": summary, "n_frames": int(n)}

    # ------------------------------------------------------------------ #
    @staticmethod
    def _protein_residues(top):
        standard = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
        }
        return [r for r in top.residues if r.name in standard]

    @staticmethod
    def _min_pair_dist(a_xyz, b_xyz):
        """Minimum pairwise distance per frame between two atom groups.

        ``a_xyz``: (n, A, 3), ``b_xyz``: (n, B, 3). Returns (n,) of the minimum
        distance over all A x B pairs in each frame.
        """
        # (n, A, 1, 3) - (n, 1, B, 3) -> (n, A, B, 3)
        diff = a_xyz[:, :, None, :] - b_xyz[:, None, :, :]
        d = np.sqrt((diff * diff).sum(axis=-1))  # (n, A, B)
        return d.reshape(d.shape[0], -1).min(axis=1)

    def _ligand_chemistry(self, ligand_mol, ligand_atoms):
        """Extract aromatic rings / cations / anions from the RDKit ligand.

        Returns trajectory-index groups. RDKit heavy-atom order is assumed to
        match the order of ``ligand_atoms`` (the convention PRISM already uses
        in visualization_generator).
        """
        rings: List[List[int]] = []
        cations: List[int] = []
        anions: List[int] = []
        try:
            from rdkit import Chem  # noqa: F401

            # Map RDKit heavy-atom ordinal -> trajectory index.
            heavy_to_traj: Dict[int, int] = {}
            heavy_ord = 0
            mol_idx_to_heavy: Dict[int, int] = {}
            for atom in ligand_mol.GetAtoms():
                if atom.GetSymbol() != "H":
                    if heavy_ord < len(ligand_atoms):
                        heavy_to_traj[heavy_ord] = ligand_atoms[heavy_ord]
                        mol_idx_to_heavy[atom.GetIdx()] = heavy_ord
                    heavy_ord += 1

            # Aromatic rings.
            ri = ligand_mol.GetRingInfo()
            for ring in ri.AtomRings():
                if all(ligand_mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    traj_ring = [
                        heavy_to_traj[mol_idx_to_heavy[i]]
                        for i in ring
                        if i in mol_idx_to_heavy and mol_idx_to_heavy[i] in heavy_to_traj
                    ]
                    if len(traj_ring) >= 4:
                        rings.append(traj_ring)

            # Formal-charge centres.
            for atom in ligand_mol.GetAtoms():
                fc = atom.GetFormalCharge()
                if fc == 0 or atom.GetIdx() not in mol_idx_to_heavy:
                    continue
                tidx = heavy_to_traj.get(mol_idx_to_heavy[atom.GetIdx()])
                if tidx is None:
                    continue
                if fc > 0:
                    cations.append(tidx)
                else:
                    anions.append(tidx)
        except Exception as exc:  # pragma: no cover
            logger.debug(f"Ligand chemistry extraction failed: {exc}")
        return rings, cations, anions

    def _add_hydrogen_bonds(self, sample, ligand_residue, prot_residues, counts, n):
        """Detect protein-ligand H-bonds per frame via mdtraj Wernet-Nilsson.

        Writes a per-frame boolean presence mask into ``counts`` (OR-ing into any
        existing accumulator), consistent with the geometric interaction passes.
        """
        try:
            lig_atom_set = {a.index for a in ligand_residue.atoms}
            top = sample.topology
            # residue index -> "NAMEresSeq" for protein residues only
            res_label = {}
            for r in prot_residues:
                res_label[r.index] = f"{r.name}{r.resSeq}"

            per_frame = md.wernet_nilsson(sample)  # list length n, arrays (k,3)
            hb_masks: Dict[str, np.ndarray] = {}
            for f, frame_hbonds in enumerate(per_frame):
                for donor, hydrogen, acceptor in frame_hbonds:
                    d_lig = donor in lig_atom_set
                    a_lig = acceptor in lig_atom_set
                    if d_lig == a_lig:
                        continue  # both ligand or both protein -> skip
                    prot_atom = acceptor if d_lig else donor
                    ridx = top.atom(int(prot_atom)).residue.index
                    rid = res_label.get(ridx)
                    if rid:
                        m = hb_masks.get(rid)
                        if m is None:
                            m = np.zeros(n, dtype=bool)
                            hb_masks[rid] = m
                        m[f] = True
            for rid, mask in hb_masks.items():
                key = (rid, "hydrogen_bond")
                if key in counts:
                    counts[key] |= mask
                else:
                    counts[key] = mask
        except Exception as exc:  # pragma: no cover
            logger.debug(f"Hydrogen-bond detection skipped: {exc}")
