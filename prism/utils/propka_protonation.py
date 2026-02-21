#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein protonation state prediction using PROPKA3.

This module provides a wrapper around PROPKA3 to predict pKa values and
protonation states of ionizable residues in protein structures at a given pH.

Unlike Meeko (which only matches existing hydrogens), PROPKA3 actually
predicts pKa values based on the 3D structural environment and determines
the most probable protonation state at the specified pH.

References:
    - Olsson et al. (2011) PROPKA3: Consistent Treatment of Internal
      and Surface Residues in Empirical pKa predictions.
      J. Chem. Theory Comput. 7, 525-537.
    - Søndergaard et al. (2011) Improved pKa predictions through Poisson-
      Boltzmann coupled PROPKA. J. Chem. Theory Comput. 7, 3968-3979.
"""

import os
import subprocess
import logging
import tempfile
import shutil
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Any
import re

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
    "N+ ": "N+",   # N-terminus
    "C- ": "C-",   # C-terminus
}

# Ionizable residues handled by PROPKA
IONIZABLE_RESIDUES = {
    "ARG": {"pka_ref": 12.0, "always_positive": True},
    "ASP": {"pka_ref": 3.9, "acidic": True},
    "CYS": {"pka_ref": 8.3, "acidic": True},
    "GLU": {"pka_ref": 4.1, "acidic": True},
    "HIS": {"pka_ref": 6.5, "can_be_positive": True},
    "LYS": {"pka_ref": 10.5, "always_positive": True},
    "TYR": {"pka_ref": 10.1, "acidic": True},
    "N+ ": {"pka_ref": 8.0, "always_positive": True},  # N-terminus
    "C- ": {"pka_ref": 3.5, "acidic": True},        # C-terminus
}


class PropkaProtonator:
    """
    Protein protonation state predictor using PROPKA3.

    This class uses PROPKA3 to predict pKa values of ionizable residues
    and determines the protonation state at a specified pH.

    Parameters
    ----------
    ph : float
        Target pH for protonation state prediction (default: 7.0)
    propka_path : str, optional
        Path to propka3 executable (default: searches PATH)
    window_size : float
        pH window size around pKa for uncertain predictions (default: 1.0)
        If |pH - pKa| < window_size, the prediction is considered uncertain.
    verbose : bool
        Enable verbose logging (default: False)
    """

    def __init__(self,
                 ph: float = 7.0,
                 propka_path: Optional[str] = None,
                 window_size: float = 1.0,
                 verbose: bool = False):
        self.ph = ph
        self.window_size = window_size
        self.verbose = verbose

        # Find PROPKA3 executable
        self.propka_path = propka_path or self._find_propka()
        self._check_propka_available()

    def _find_propka(self) -> str:
        """Find PROPKA3 executable in PATH or common locations."""
        # Check common names
        candidates = ["propka3", "propka31", "propka"]

        for name in candidates:
            path = shutil.which(name)
            if path:
                logger.info(f"Found PROPKA3 at: {path}")
                return path

        # Check conda environment
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            for name in candidates:
                path = Path(conda_prefix) / "bin" / name
                if path.exists():
                    logger.info(f"Found PROPKA3 at: {path}")
                    return str(path)

        raise RuntimeError(
            "PROPKA3 not found. Install with:\n"
            "  conda install -c conda-forge propka\n"
            "Or specify path with propka_path parameter."
        )

    def _check_propka_available(self):
        """Check if PROPKA3 is available and working."""
        try:
            result = subprocess.run(
                [self.propka_path, "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            logger.info(f"PROPKA version: {result.stdout.strip() or 'unknown'}")
        except Exception as e:
            raise RuntimeError(
                f"PROPKA3 check failed: {e}\n"
                f"Executable: {self.propka_path}"
            ) from e

    def predict_pka(self,
                   input_pdb: str,
                   work_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        Run PROPKA3 to predict pKa values.

        Parameters
        ----------
        input_pdb : str
            Path to input PDB file
        work_dir : str, optional
            Working directory for PROPKA3 output (default: temp directory)

        Returns
        -------
        dict
            Dictionary with pKa predictions for each ionizable residue
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

        logger.info(f"Running PROPKA3 for {input_path.name}")
        logger.info(f"Working directory: {work_path}")

        try:
            # Run PROPKA3
            cmd = [
                self.propka_path,
                str(input_path),
                "-o", str(self.ph),
                "--quiet"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                cwd=str(work_path)
            )

            if result.returncode != 0:
                raise RuntimeError(f"PROPKA3 failed: {result.stderr}")

            # Parse output files
            pka_file = work_path / f"{input_path.stem}.pka"
            propka_input = work_path / f"{input_path.stem}.propka_input"

            predictions = {}

            if pka_file.exists():
                predictions = self._parse_pka_file(pka_file)
            else:
                logger.warning(f"PROPKA3 .pka file not found, trying .propka_input")

            if propka_input.exists():
                predictions.update(self._parse_propka_input(propka_input))

            if not predictions:
                logger.warning("No pKa predictions found. PROPKA3 may have failed silently.")

            return predictions

        finally:
            if cleanup:
                shutil.rmtree(work_path, ignore_errors=True)

    def _parse_pka_file(self, pka_file: Path) -> Dict[str, Any]:
        """Parse PROPKA3 .pka output file.

        The SUMMARY section has format:
            RESIDUE CHAIN NUM pKa    model-pKa
            HIS     A     69   5.13    6.50
        """
        predictions = {}
        in_summary = False

        with open(pka_file, 'r') as f:
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
                    # Format: "HIS   69 A     5.13       6.50"
                    # or: "ASP   7 A     4.10       3.80"
                    parts = line.split()
                    if len(parts) >= 4:
                        res_name = parts[0]
                        res_num = parts[1]
                        chain = parts[2]
                        pka_str = parts[3]

                        # Validate pKa is numeric
                        try:
                            pka = float(pka_str)
                            key = f"{chain}:{res_num}"
                            predictions[key] = {
                                "resname": res_name,
                                "chain": chain,
                                "resnum": res_num,
                                "pka": pka
                            }
                        except ValueError:
                            continue

        return predictions

    def _parse_propka_input(self, input_file: Path) -> Dict[str, Any]:
        """Parse PROPKA3 .propka_input file for detailed predictions."""
        predictions = {}

        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # Format: RESNAME CHAIN RESNUM PKA
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        res_name = parts[0]
                        chain = parts[1]
                        res_num = parts[2]
                        pka = float(parts[3])

                        key = f"{chain}:{res_num}"
                        predictions[key] = {
                            "resname": res_name,
                            "chain": chain,
                            "resnum": res_num,
                            "pka": pka
                        }
                    except (ValueError, IndexError):
                        continue

        return predictions

    def determine_protonation_state(self,
                                   res_name: str,
                                   pka: float) -> Tuple[str, bool]:
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
            # For histidine, we need to determine HID vs HIE vs HIP
            # This is simplified; actual state depends on local environment
            if ph > pka + 1.0:
                # Fully protonated (positive)
                return "HIP", is_certain
            elif ph < pka - 1.0:
                # Fully deprotonated (neutral)
                # Default to HIE (most common at pH 7)
                return "HIE", is_certain
            else:
                # Transition zone - check hydrogen positions later
                # Default to HIE for pH ~7
                return "HIE", False

        elif res_name in ("ASP", "GLU"):
            # Acidic: pH > pKa → deprotonated (COO-, normal)
            #        pH < pKa → protonated (COOH, unusual at pH 7)
            if ph > pka:
                return res_name, is_certain  # Deprotonated (COO-)
            else:
                return f"{res_name[0]}SH", is_certain  # Protonated (COOH)

        elif res_name == "LYS":
            # Basic: pH < pKa → protonated (NH3+, normal)
            #        pH > pKa → deprotonated (NH2, unusual at pH 7)
            if ph < pka:
                return "LYS", is_certain  # Protonated (NH3+)
            else:
                return "LYN", is_certain  # Deprotonated (NH2)

        elif res_name == "CYS":
            # pH < pKa → protonated (SH, normal)
            # pH > pKa → deprotonated (S-, unusual at pH 7)
            if ph < pka:
                return "CYS", is_certain  # Protonated (SH)
            else:
                return "CYM", is_certain  # Deprotonated (S-)

        elif res_name == "TYR":
            # pH < pKa → protonated (OH, normal)
            # pH > pKa → deprotonated (O-, unusual at pH 7)
            if ph < pka:
                return "TYR", is_certain  # Protonated (OH)
            else:
                return "TYH", is_certain  # Deprotonated (O-)

        elif res_name == "ARG":
            # Arg always protonated below pH ~12
            return "ARG", is_certain

        else:
            return res_name, is_certain

    def optimize_protein_protonation(self,
                                    input_pdb: str,
                                    output_pdb: str,
                                    work_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        Predict and apply protonation states to a protein structure.

        Parameters
        ----------
        input_pdb : str
            Path to input PDB file
        output_pdb : str
            Path to output PDB file with protonation states applied
        work_dir : str, optional
            Working directory for PROPKA3 output

        Returns
        -------
        dict
            Results dictionary with:
            - pka_predictions: pKa values for each residue
            - protonation_states: predicted states at target pH
            - uncertain: list of residues with uncertain predictions
            - renames: dictionary of residue renames applied
        """
        # Get pKa predictions
        pka_predictions = self.predict_pka(input_pdb, work_dir)

        if not pka_predictions:
            logger.warning("No pKa predictions obtained from PROPKA3")
            return {
                "pka_predictions": {},
                "protonation_states": {},
                "uncertain": [],
                "renames": {}
            }

        # Determine protonation states
        protonation_states = {}
        uncertain = []
        renames = {}

        for key, pred in pka_predictions.items():
            res_name = pred["resname"]
            pka = pred["pka"]

            if pka is None:
                logger.warning(f"No pKa predicted for {key}")
                continue

            state, is_certain = self.determine_protonation_state(res_name, pka)

            protonation_states[key] = {
                "original_resname": res_name,
                "protonation_state": state,
                "pka": pka,
                "ph": self.ph,
                "certain": is_certain
            }

            if not is_certain:
                uncertain.append(key)

            # Track if renaming needed (e.g., HIS -> HID)
            if res_name != state and state != res_name + "H":
                renames[key] = {"from": res_name, "to": state}

        # Apply renames to PDB file
        self._apply_protonation_renames(input_pdb, output_pdb, renames)

        return {
            "pka_predictions": pka_predictions,
            "protonation_states": protonation_states,
            "uncertain": uncertain,
            "renames": renames,
            "input": str(input_pdb),
            "output": str(output_pdb)
        }

    def _apply_protonation_renames(self,
                                   input_pdb: str,
                                   output_pdb: str,
                                   renames: Dict[str, Dict[str, str]]):
        """Apply protonation state renames to PDB file."""
        input_path = Path(input_pdb)
        output_path = Path(output_pdb)

        with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
            for line in f_in:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    f_out.write(line)
                    continue

                res_name = line[17:20].strip()
                chain = line[21].strip()
                res_num = line[22:26].strip()
                key = f"{chain}:{res_num}"

                if key in renames:
                    new_res_name = renames[key]["to"]
                    # Replace residue name (columns 18-20 in 1-indexed, 17-19 in 0-indexed)
                    new_line = line[:17] + f"{new_res_name:3s}" + line[20:]
                    f_out.write(new_line)
                else:
                    f_out.write(line)

        logger.info(f"Applied {len(renames)} protonation state renames")

    def print_summary(self, results: Dict[str, Any]):
        """Print a summary of protonation predictions."""
        print(f"\n{'='*60}")
        print(f"PROPKA3 Protonation Prediction Summary (pH {self.ph})")
        print(f"{'='*60}")

        if not results["protonation_states"]:
            print("No predictions available.")
            return

        # Group by residue type
        by_type = {}
        for key, state in results["protonation_states"].items():
            orig = state["original_resname"]
            if orig not in by_type:
                by_type[orig] = []
            by_type[orig].append((key, state))

        # Print ionizable residues
        for res_type in ["HIS", "ASP", "GLU", "CYS", "LYS", "TYR"]:
            if res_type not in by_type:
                continue

            print(f"\n{res_type} residues:")
            for key, state in by_type[res_type]:
                cert = "✓" if state["certain"] else "?"
                print(f"  {key}: {state['original_resname']} -> {state['protonation_state']} "
                      f"(pKa={state['pka']:.2f}, pH={state['ph']:.1f}) {cert}")

        if results["uncertain"]:
            print(f"\n⚠ Uncertain predictions (|pH-pKa| < {self.window_size}):")
            for key in results["uncertain"]:
                state = results["protonation_states"][key]
                print(f"  {key}: pKa={state['pka']:.2f}, state={state['protonation_state']}")

        print(f"\n{'='*60}\n")


def optimize_protein_protonation_propka(input_pdb: str,
                                      output_pdb: str,
                                      ph: float = 7.0,
                                      propka_path: Optional[str] = None,
                                      window_size: float = 1.0,
                                      work_dir: Optional[str] = None,
                                      verbose: bool = True) -> Dict[str, Any]:
    """
    Convenience function to optimize protein protonation states using PROPKA3.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output PDB file with optimized protonation states
    ph : float
        Target pH (default: 7.0)
    propka_path : str, optional
        Path to PROPKA3 executable
    window_size : float
        pH window for uncertain predictions (default: 1.0)
    work_dir : str, optional
        Working directory for PROPKA3 output
    verbose : bool
        Print summary (default: True)

    Returns
    -------
    dict
        Results dictionary with predictions and renames

    Examples
    --------
    >>> results = optimize_protein_protonation_propka("protein.pdb", "protein_prot.pdb", ph=7.4)
    >>> print(results["protonation_states"]["A:69"])
    """
    protonator = PropkaProtonator(
        ph=ph,
        propka_path=propka_path,
        window_size=window_size,
        verbose=verbose
    )

    results = protonator.optimize_protein_protonation(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        work_dir=work_dir
    )

    if verbose:
        protonator.print_summary(results)

    return results


if __name__ == "__main__":
    import sys

    # Example usage
    if len(sys.argv) < 3:
        print("Usage: python propka_protonation.py <input.pdb> <output.pdb> [pH]")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    ph = float(sys.argv[3]) if len(sys.argv) > 3 else 7.0

    results = optimize_protein_protonation_propka(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        ph=ph
    )
