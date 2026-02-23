#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein protonation state prediction using PROPKA.

This module uses PROPKA to predict per-residue pKa values and determine
protonation states of ionizable residues (HIS, ASP, GLU, LYS, CYS, TYR)
based on the protein's 3D environment and a target pH.

The renamed residues are then passed to GROMACS pdb2gmx, which regenerates
standardized hydrogens.

For histidine protonation states, the module predicts HID/HIE/HIP based on
pKa values. For other residues, it can predict alternative protonation states
(e.g., ASH, GLH, LYN, CYM, TYH) when the pH differs significantly from the
residue's pKa.

References:
    - Olsson et al. (2011) PROPKA3: Consistent Treatment of Internal
      and Surface Residues in Empirical pKa predictions.
      J. Chem. Theory Comput. 7, 525-537.
    - Søndergaard et al. (2011) Improved pKa predictions through Poisson-
      Boltzmann coupled PROPKA. J. Chem. Theory Comput. 7, 3968-3979.
"""

import os
import re
import subprocess
import logging
import tempfile
import shutil
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional, Tuple, Any, Set

logger = logging.getLogger(__name__)


# Mapping from PROPKA residue names to GROMACS/standard names
PROPKA_TO_GROMACS = {
    "HIS": "HIS",  # Unspecified, will be determined by hydrogens
    "HID": "HID",  # ND1 protonated
    "HIE": "HIE",  # NE2 protonated
    "HIP": "HIP",  # Both protonated
    "ASP": "ASP",  # Deprotonated (COO-)
    "ASH": "ASH",  # Protonated (COOH)
    "GLU": "GLU",  # Deprotonated (COO-)
    "GLH": "GLH",  # Protonated (COOH)
    "LYS": "LYS",  # Protonated (NH3+)
    "LYN": "LYN",  # Deprotonated (NH2)
    "CYS": "CYS",  # Protonated (SH)
    "CYM": "CYM",  # Deprotonated (S-)
    "TYR": "TYR",  # Deprotonated (O-)
    "TYH": "TYH",  # Protonated (OH)
    "ARG": "ARG",  # Always protonated at pH < 12
    "N+ ": "N+",  # N-terminus
    "C- ": "C-",  # C-terminus
}


# Preferred residue name aliases per force field (checked against aminoacids.rtp)
PROTONATION_ALIASES = {
    # Histidine: Amber uses HID/HIE/HIP, CHARMM uses HSD/HSE/HSP
    "HID": ["HID", "HSD"],
    "HIE": ["HIE", "HSE"],
    "HIP": ["HIP", "HSP"],
    # Protonated acidic residues
    "ASH": ["ASH", "ASPP", "ASPH"],
    "GLH": ["GLH", "GLUP", "GLUH"],
    # Deprotonated cysteine
    "CYM": ["CYM", "CYD"],
    # Deprotonated tyrosine (rare)
    "TYH": ["TYH", "TYM"],
}

PROTONATION_FALLBACKS = {
    # If neutral lysine is unavailable, fall back to LYS (charged)
    "LYN": "LYS",
    # If deprotonated tyrosine is unavailable, fall back to TYR
    "TYH": "TYR",
}


def _ff_rtp_path(ff_info: Optional[Dict[str, Any]]) -> Optional[Path]:
    """Get aminoacids.rtp path for a force field info dict."""
    if not ff_info:
        return None
    ff_path = ff_info.get("path") or ff_info.get("dir")
    if not ff_path:
        return None
    ff_dir = Path(ff_path)
    if not ff_dir.is_dir():
        return None
    rtp = ff_dir / "aminoacids.rtp"
    return rtp if rtp.exists() else None


@lru_cache(maxsize=16)
def _load_rtp_residue_names(rtp_path_str: str) -> Set[str]:
    """Parse residue names from aminoacids.rtp."""
    names: Set[str] = set()
    try:
        with open(rtp_path_str, "r") as f:
            for line in f:
                match = re.match(r"^\s*\[\s*([A-Za-z0-9+\-]+)\s*\]", line)
                if match:
                    names.add(match.group(1))
    except Exception as exc:
        logger.warning(f"Failed to parse rtp file {rtp_path_str}: {exc}")
    return names


def get_ff_residue_names(ff_info: Optional[Dict[str, Any]]) -> Optional[Set[str]]:
    """Return residue names defined in aminoacids.rtp for a force field."""
    rtp_path = _ff_rtp_path(ff_info)
    if not rtp_path:
        return None
    return _load_rtp_residue_names(str(rtp_path))


def resolve_protonation_resname(
    original_resname: str, desired_state: str, ff_info: Optional[Dict[str, Any]]
) -> Tuple[str, Optional[str]]:
    """
    Resolve a protonation state name to a force-field-compatible residue name.

    Returns (resolved_name, note). If note is not None, a mapping/fallback occurred.
    """
    ff_residues = get_ff_residue_names(ff_info)
    if not ff_residues:
        return desired_state, None

    ff_name = ""
    if ff_info:
        ff_name = str(ff_info.get("name") or ff_info.get("dir") or "").lower()
    is_charmm = "charmm" in ff_name

    if desired_state in ff_residues:
        # CHARMM: HIP is phosphohistidine; prefer HSP for protonated His
        if is_charmm and desired_state == "HIP" and "HSP" in ff_residues:
            return "HSP", "HIP->HSP (CHARMM protonated His)"
        return desired_state, None

    aliases = list(PROTONATION_ALIASES.get(desired_state, []))
    if is_charmm:
        if desired_state == "HIP":
            aliases = ["HSP"]
        elif desired_state == "HID":
            aliases = ["HSD", "HID"]
        elif desired_state == "HIE":
            aliases = ["HSE", "HIE"]

    for alt in aliases:
        if alt in ff_residues:
            return alt, f"{desired_state}->{alt}"

    if original_resname and original_resname in ff_residues:
        return original_resname, f"{desired_state}->{original_resname} (fallback)"

    fallback = PROTONATION_FALLBACKS.get(desired_state)
    if fallback and fallback in ff_residues:
        return fallback, f"{desired_state}->{fallback} (fallback)"

    return desired_state, f"{desired_state} (unmapped)"


def default_histidine_name(ff_info: Optional[Dict[str, Any]]) -> str:
    """
    Pick a neutral histidine name supported by the selected force field.
    Preference: HIE (Amber) -> HSE (CHARMM) -> HIS -> HSD.
    """
    ff_residues = get_ff_residue_names(ff_info)
    if not ff_residues:
        return "HIE"
    for candidate in ("HIE", "HSE", "HIS", "HSD"):
        if candidate in ff_residues:
            return candidate
    return "HIS"


class PropkaProtonator:
    """
    Predict protonation states using PROPKA pKa calculations.

    PROPKA is an empirical pKa predictor that estimates per-residue pKa
    values from the protein's 3D structure.

    This class provides both simple histidine-focused prediction (for backward
    compatibility) and full residue protonation prediction (all ionizable
    residues).
    """

    def __init__(self, ph: float = 7.0, verbose: bool = False):
        """
        Initialize the PROPKA protonator.

        Parameters
        ----------
        ph : float
            Target pH for protonation state prediction (default: 7.0)
        verbose : bool
            Enable verbose logging (default: False)
        """
        self.ph = ph
        self.verbose = verbose
        self.window_size = 1.0
        self.propka_path = self._find_propka()
        self._check_propka_available()

    def _find_propka(self) -> str:
        """Find PROPKA executable in PATH or common locations."""
        import shutil

        candidates = ["propka3", "propka31", "propka"]

        for name in candidates:
            path = shutil.which(name)
            if path:
                logger.info(f"Found PROPKA at: {path}")
                return path

        # Check conda environment
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            for name in candidates:
                path = Path(conda_prefix) / "bin" / name
                if path.exists():
                    logger.info(f"Found PROPKA at: {path}")
                    return str(path)

        raise RuntimeError(
            "PROPKA not found. Install with:\n" "  conda install -c conda-forge propka\n" "Or pip install propka"
        )

    def _check_propka_available(self):
        """Check if PROPKA is installed and working."""
        try:
            result = subprocess.run([self.propka_path, "--version"], capture_output=True, text=True, timeout=10)
            logger.info(f"PROPKA version: {result.stdout.strip() or 'unknown'}")
        except Exception as e:
            raise RuntimeError(f"PROPKA check failed: {e}\n" f"Executable: {self.propka_path}") from e

    def predict_his_states(self, pdb_file: str) -> Dict[tuple, str]:
        """
        Predict histidine protonation states from pKa values.

        This is a simplified interface for backward compatibility that returns
        only histidine residues in the format expected by the Builder.

        Parameters
        ----------
        pdb_file : str
            Path to input PDB file

        Returns
        -------
        dict
            Mapping of (chain, resnum) -> protonation state ("HIE", "HID", or "HIP")
        """
        import propka.run

        pdb_path = Path(pdb_file)
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        logger.info(f"Running PROPKA pKa prediction on {pdb_path.name} at pH {self.ph}")

        mol = propka.run.single(str(pdb_path), optargs=["--quiet"], write_pka=False)
        conf = mol.conformations[mol.conformation_names[0]]

        his_states = {}
        for group in conf.get_titratable_groups():
            if group.residue_type == "HIS":
                chain = group.atom.chain_id.strip()
                resnum = str(group.atom.res_num)

                if group.pka_value > self.ph:
                    state = "HIP"  # Doubly protonated (pKa above pH)
                else:
                    state = "HIE"  # Neutral, NE2 protonated (default for neutral HIS)

                his_states[(chain, resnum)] = state

                if self.verbose:
                    logger.info(f"  HIS {chain}:{resnum} pKa={group.pka_value:.2f} -> {state}")

        return his_states

    def predict_pka(self, input_pdb: str, work_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        Run PROPKA to predict pKa values for all ionizable residues.

        This is the full-featured method that handles all ionizable residues:
        HIS, ASP, GLU, LYS, CYS, TYR, ARG.

        Parameters
        ----------
        input_pdb : str
            Path to input PDB file
        work_dir : str, optional
            Working directory for PROPKA output (default: temp directory)

        Returns
        -------
        dict
            Dictionary with pKa predictions for each ionizable residue.
            Format: { "chain:resnum": {"resname": str, "chain": str, "resnum": str, "pka": float} }
        """
        input_path = Path(input_pdb).resolve()

        if not input_path.exists():
            raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

        # Create working directory
        if work_dir is None:
            temp_dir = tempfile.mkdtemp(prefix="propka_pka_")
            work_path = Path(temp_dir)
            cleanup = True
        else:
            work_path = Path(work_dir)
            work_path.mkdir(parents=True, exist_ok=True)
            cleanup = False

        logger.info(f"Running PROPKA for {input_path.name}")
        logger.info(f"Working directory: {work_path}")

        try:
            # Run PROPKA
            cmd = [self.propka_path, str(input_path), "-o", str(self.ph), "--quiet"]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=str(work_path))

            if result.returncode != 0:
                raise RuntimeError(f"PROPKA failed: {result.stderr}")

            # Parse output files
            pka_file = work_path / f"{input_path.stem}.pka"
            propka_input = work_path / f"{input_path.stem}.propka_input"

            predictions = {}

            if pka_file.exists():
                predictions = self._parse_pka_file(pka_file)
            else:
                logger.warning(f"PROPKA .pka file not found, trying .propka_input")

            if propka_input.exists():
                predictions.update(self._parse_propka_input(propka_input))

            if not predictions:
                logger.warning("No pKa predictions found. PROPKA may have failed silently.")

            return predictions

        finally:
            if cleanup:
                shutil.rmtree(work_path, ignore_errors=True)

    def _parse_pka_file(self, pka_file: Path) -> Dict[str, Any]:
        """Parse PROPKA .pka output file.

        The SUMMARY section has format:
            RESIDUE CHAIN NUM pKa    model-pKa
            HIS     A     69   5.13    6.50
        """
        predictions = {}
        in_summary = False

        with open(pka_file, "r") as f:
            for line in f:
                line = line.strip()

                # Start of summary section
                if "SUMMARY OF THIS PREDICTION" in line:
                    in_summary = True
                    continue

                # Skip header line in summary
                if in_summary and "Group" in line and "pKa" in line:
                    continue

                # Parse summary lines
                if in_summary:
                    parts = line.split()
                    if len(parts) >= 4:
                        res_name = parts[0]
                        res_num = parts[1]
                        chain = parts[2]
                        pka_str = parts[3]

                        try:
                            pka = float(pka_str)
                            key = f"{chain}:{res_num}"
                            predictions[key] = {"resname": res_name, "chain": chain, "resnum": res_num, "pka": pka}
                        except ValueError:
                            continue

        return predictions

    def _parse_propka_input(self, input_file: Path) -> Dict[str, Any]:
        """Parse PROPKA .propka_input file for detailed predictions."""
        predictions = {}

        with open(input_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) >= 4:
                    try:
                        res_name = parts[0]
                        chain = parts[1]
                        res_num = parts[2]
                        pka = float(parts[3])

                        key = f"{chain}:{res_num}"
                        predictions[key] = {"resname": res_name, "chain": chain, "resnum": res_num, "pka": pka}
                    except (ValueError, IndexError):
                        continue

        return predictions

    def determine_protonation_state(self, res_name: str, pka: float) -> Tuple[str, bool]:
        """
        Determine protonation state based on pKa and target pH.

        Parameters
        ----------
        res_name : str
            Residue name (e.g., "HIS", "ASP", "GLU", "LYS")
        pka : float
            Predicted pKa value

        Returns
        -------
        tuple (state_name, is_certain)
            state_name: GROMACS-compatible residue name
            is_certain: True if |pH - pKa| >= window_size
        """
        ph = self.ph
        delta = ph - pka
        is_certain = abs(delta) >= self.window_size

        if res_name == "HIS":
            if ph > pka + 1.0:
                return "HIP", is_certain
            elif ph < pka - 1.0:
                return "HIE", is_certain
            else:
                return "HIE", False

        elif res_name in ("ASP", "GLU"):
            if ph > pka:
                return res_name, is_certain  # Deprotonated (COO-)
            else:
                return ("ASH" if res_name == "ASP" else "GLH"), is_certain  # Protonated (COOH)

        elif res_name == "LYS":
            if ph < pka:
                return "LYS", is_certain  # Protonated (NH3+)
            else:
                return "LYN", is_certain  # Deprotonated (NH2)

        elif res_name == "CYS":
            if ph < pka:
                return "CYS", is_certain  # Protonated (SH)
            else:
                return "CYM", is_certain  # Deprotonated (S-)

        elif res_name == "TYR":
            if ph < pka:
                return "TYR", is_certain  # Protonated (OH)
            else:
                return "TYH", is_certain  # Deprotonated (O-)

        elif res_name == "ARG":
            return "ARG", is_certain

        else:
            return res_name, is_certain

    def rename_histidines(self, pdb_file: str, output_file: str, ff_info: Optional[Dict[str, Any]] = None) -> Dict:
        """
        Predict protonation states and rename HIS residues in PDB file.

        This method only modifies HIS residues (HID/HIE/HIP) for backward
        compatibility. For full residue protonation, use optimize_protein_protonation().

        Parameters
        ----------
        pdb_file : str
            Path to input PDB file
        output_file : str
            Path to output PDB file with renamed residues
        ff_info : dict, optional
            Force field info dict (used to map residue names to available rtp entries)

        Returns
        -------
        dict
            Statistics with keys: 'total_his', 'renamed', 'states'
        """
        his_states = self.predict_his_states(pdb_file)

        stats = {
            "total_his": len(his_states),
            "renamed": {},
            "states": his_states,
        }

        if not his_states:
            if pdb_file != output_file:
                shutil.copy2(pdb_file, output_file)
            return stats

        # Read PDB and rename HIS residues
        with open(pdb_file, "r") as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20] == "HIS" and line[20] == " ":
                chain = line[21].strip()
                resnum = line[22:26].strip()
                key = (chain, resnum)

                if key in his_states:
                    desired = his_states[key]
                    new_name, note = resolve_protonation_resname("HIS", desired, ff_info)
                    if note and self.verbose:
                        logger.info(f"  HIS {chain}:{resnum} {note}")
                    line = line[:17] + f"{new_name:3s}" + line[20:]
                    stats["renamed"][key] = new_name

            new_lines.append(line)

        with open(output_file, "w") as f:
            f.writelines(new_lines)

        return stats

    def optimize_protein_protonation(
        self,
        input_pdb: str,
        output_pdb: str,
        work_dir: Optional[str] = None,
        ff_info: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Predict and apply protonation states for all ionizable residues.

        This method handles all ionizable residues (HIS, ASP, GLU, LYS, CYS, TYR)
        and renames them according to their predicted protonation states at the
        target pH.

        Parameters
        ----------
        input_pdb : str
            Path to input PDB file
        output_pdb : str
            Path to output PDB file with protonation states applied
        work_dir : str, optional
            Working directory for PROPKA output
        ff_info : dict, optional
            Force field info dict (used to map residue names to available rtp entries)

        Returns
        -------
        dict
            Results with keys: 'predictions', 'renamed', 'summary'
        """
        # Get pKa predictions
        predictions = self.predict_pka(input_pdb, work_dir)

        # Determine protonation states
        renames = {}
        summary = {"normal": [], "altered": []}

        for key, pred in predictions.items():
            resname = pred["resname"]
            pka = pred["pka"]
            chain = pred["chain"]
            resnum = pred["resnum"]

            desired_state, certain = self.determine_protonation_state(resname, pka)
            mapped_state, note = resolve_protonation_resname(resname, desired_state, ff_info)
            if note:
                msg = f"Residue {chain}:{resnum} {resname}: {note}"
                if "fallback" in note or "unmapped" in note:
                    logger.warning(msg)
                elif self.verbose:
                    logger.info(msg)
            new_state = mapped_state

            if new_state != resname:
                renames[(chain, resnum, resname)] = new_state
                summary["altered"].append(f"{chain}:{resnum} {resname}->{new_state} (pKa={pka:.2f})")
            else:
                summary["normal"].append(f"{chain}:{resnum} {resname} (pKa={pka:.2f})")

        # Apply renames to PDB
        if renames:
            with open(input_pdb, "r") as f:
                lines = f.readlines()

            new_lines = []
            for line in lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    line_resname = line[17:20].strip()
                    chain = line[21].strip()
                    resnum = line[22:26].strip()

                    key = (chain, resnum, line_resname)
                    if key in renames:
                        new_name = renames[key]
                        line = line[:17] + f"{new_name:3s}" + line[20:]

                new_lines.append(line)

            with open(output_pdb, "w") as f:
                f.writelines(new_lines)
        else:
            # If input and output are the same path, no action needed
            if Path(input_pdb).resolve() != Path(output_pdb).resolve():
                shutil.copy2(input_pdb, output_pdb)

        return {"predictions": predictions, "renamed": renames, "summary": summary}


def optimize_protein_protonation_propka(
    input_pdb: str, output_pdb: str, ph: float = 7.0, ff_info: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Convenience function to optimize protein protonation using PROPKA.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output PDB file
    ph : float
        Target pH (default: 7.0)
    ff_info : dict, optional
        Force field info dict (used to map residue names to available rtp entries)

    Returns
    -------
    dict
        Results from protonation prediction
    """
    protonator = PropkaProtonator(ph=ph, verbose=True)
    return protonator.optimize_protein_protonation(input_pdb, output_pdb, ff_info=ff_info)
