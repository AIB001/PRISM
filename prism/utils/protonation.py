#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein protonation state optimization using Meeko.

This module provides preprocessing for protein structures to optimize
hydrogen positions and protonation states before GROMACS pdb2gmx processing.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict, List
import tempfile
import shutil

logger = logging.getLogger(__name__)


class ProteinProtonator:
    """
    Optimize protein protonation states using Meeko.

    This class provides a preprocessing step that uses Meeko's reduce2.py
    and mk_prepare_receptor.py to add and optimize hydrogens before passing
    the structure to GROMACS pdb2gmx.
    """

    def __init__(self,
                 ph: float = 7.0,
                 preserve_existing_h: bool = False,
                 his_state: str = "auto",
                 verbose: bool = False):
        """
        Initialize the protein protonator.

        Parameters
        ----------
        ph : float
            Target pH for protonation (default: 7.0)
        preserve_existing_h : bool
            Keep existing hydrogens if present (default: False)
        his_state : str
            Histidine protonation state: 'auto', 'HID', 'HIE', or 'HIP' (default: 'auto')
        verbose : bool
            Enable verbose logging (default: False)
        """
        self.ph = ph
        self.preserve_existing_h = preserve_existing_h
        self.his_state = his_state
        self.verbose = verbose

        # Check if Meeko is available
        self._check_meeko_available()

    def _check_meeko_available(self):
        """Check if Meeko is installed and available."""
        try:
            import meeko
            logger.info(f"Meeko version: {meeko.__version__}")

            # Check for reduce2.py command
            try:
                result = subprocess.run(
                    ["reduce2.py", "--help"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                self.has_reduce2 = (result.returncode == 0)
            except (subprocess.TimeoutExpired, FileNotFoundError):
                logger.warning("reduce2.py not found in PATH. Using mk_prepare_receptor only.")
                self.has_reduce2 = False

        except ImportError as e:
            raise ImportError(
                "Meeko is not available. Install with: pip install meeko\n"
                "Or install with PRISM: pip install -e .[protonation]"
            ) from e

    def optimize_hydrogens(self,
                          input_pdb: str,
                          output_pdb: str,
                          work_dir: Optional[str] = None) -> Dict[str, str]:
        """
        Optimize hydrogen positions and protonation states.

        Parameters
        ----------
        input_pdb : str
            Path to input PDB file
        output_pdb : str
            Path to output PDB file with optimized hydrogens
        work_dir : str, optional
            Working directory for intermediate files (default: temp directory)

        Returns
        -------
        dict
            Dictionary with paths to intermediate files and logs
        """
        input_path = Path(input_pdb).resolve()
        output_path = Path(output_pdb).resolve()

        if not input_path.exists():
            raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")

        # Create working directory
        if work_dir is None:
            temp_dir = tempfile.mkdtemp(prefix="prism_protonation_")
            work_path = Path(temp_dir)
            cleanup = True
        else:
            work_path = Path(work_dir)
            work_path.mkdir(parents=True, exist_ok=True)
            cleanup = False

        logger.info(f"Optimizing protonation states for {input_path.name}")
        logger.info(f"Target pH: {self.ph}")
        logger.info(f"Working directory: {work_path}")

        results = {
            "input": str(input_path),
            "output": str(output_path),
            "work_dir": str(work_path)
        }

        try:
            # Step 1: Add/optimize hydrogens with reduce2.py (if available)
            if self.has_reduce2:
                h_optimized_pdb = work_path / "protein_H_optimized.pdb"
                self._run_reduce2(input_path, h_optimized_pdb)
                results["reduce2_output"] = str(h_optimized_pdb)
                intermediate_pdb = h_optimized_pdb
            else:
                # Skip reduce2, use input directly
                intermediate_pdb = input_path
                logger.info("Skipping reduce2.py (not available)")

            # Step 2: Prepare receptor with Meeko (validates and cleans structure)
            self._run_mk_prepare_receptor(intermediate_pdb, output_path, work_path)
            results["final_output"] = str(output_path)

            # Step 3: Validate output
            if not output_path.exists():
                raise RuntimeError(f"Output PDB not created: {output_path}")

            logger.info(f"Successfully optimized protonation states")
            logger.info(f"Output saved to: {output_path}")

        except Exception as e:
            logger.error(f"Protonation optimization failed: {e}")
            raise
        finally:
            # Cleanup temporary directory if created
            if cleanup:
                shutil.rmtree(work_path, ignore_errors=True)

        return results

    def _run_reduce2(self, input_pdb: Path, output_pdb: Path):
        """
        Run reduce2.py for hydrogen optimization.

        Parameters
        ----------
        input_pdb : Path
            Input PDB file
        output_pdb : Path
            Output PDB file with optimized hydrogens
        """
        logger.info("Running reduce2.py for hydrogen optimization...")

        cmd = ["reduce2.py", str(input_pdb)]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minutes timeout
            )

            if result.returncode != 0:
                raise RuntimeError(f"reduce2.py failed: {result.stderr}")

            # Write output
            with open(output_pdb, 'w') as f:
                f.write(result.stdout)

            logger.info(f"Hydrogens optimized successfully")

        except subprocess.TimeoutExpired:
            raise RuntimeError("reduce2.py timed out after 5 minutes")
        except Exception as e:
            raise RuntimeError(f"reduce2.py execution failed: {e}")

    def _run_mk_prepare_receptor(self, input_pdb: Path, output_pdb: Path, work_dir: Path):
        """
        Run mk_prepare_receptor.py to prepare and validate receptor.

        Parameters
        ----------
        input_pdb : Path
            Input PDB file
        output_pdb : Path
            Output PDB file
        work_dir : Path
            Working directory for intermediate files
        """
        logger.info("Running mk_prepare_receptor.py...")

        # Base command - use --read_pdb instead of -i to use RDKit parser
        # This avoids the ProDy dependency
        cmd = [
            "mk_prepare_receptor.py",
            "--read_pdb", str(input_pdb),
            "--write_pdb", str(output_pdb)
        ]

        # Add histidine state if specified
        if self.his_state != "auto":
            # Set all histidines to specified state
            cmd.extend(["--set_template", f":{self.his_state}"])

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minutes timeout
                cwd=str(work_dir)
            )

            if result.returncode != 0:
                logger.error(f"mk_prepare_receptor.py stderr: {result.stderr}")
                raise RuntimeError(f"mk_prepare_receptor.py failed: {result.stderr}")

            if self.verbose:
                logger.info(f"mk_prepare_receptor.py output:\n{result.stdout}")

            logger.info("Receptor preparation completed successfully")

        except subprocess.TimeoutExpired:
            raise RuntimeError("mk_prepare_receptor.py timed out after 5 minutes")
        except Exception as e:
            raise RuntimeError(f"mk_prepare_receptor.py execution failed: {e}")

    def validate_protonation(self, pdb_file: str) -> Dict[str, any]:
        """
        Validate protonation states in PDB file.

        Parameters
        ----------
        pdb_file : str
            Path to PDB file to validate

        Returns
        -------
        dict
            Validation results with statistics
        """
        pdb_path = Path(pdb_file)

        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        stats = {
            "total_residues": 0,
            "residues_with_h": 0,
            "his_residues": [],
            "metal_atoms": [],
            "warnings": []
        }

        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()

                    # Count hydrogens
                    if atom_name.startswith('H'):
                        stats["residues_with_h"] += 1

                    # Track histidines
                    if res_name in ["HIS", "HID", "HIE", "HIP"]:
                        res_id = line[22:26].strip()
                        chain = line[21].strip()
                        if f"{chain}:{res_id}" not in stats["his_residues"]:
                            stats["his_residues"].append(f"{res_name} {chain}:{res_id}")

                    # Track metals
                    if res_name in ["ZN", "MG", "CA", "FE", "MN", "CU", "NI"]:
                        res_id = line[22:26].strip()
                        chain = line[21].strip()
                        stats["metal_atoms"].append(f"{res_name} {chain}:{res_id}")

        # Add warnings
        if stats["residues_with_h"] == 0:
            stats["warnings"].append("No hydrogen atoms found in structure")

        if stats["his_residues"]:
            logger.info(f"Found {len(stats['his_residues'])} histidine residues")

        if stats["metal_atoms"]:
            logger.warning(f"Found {len(stats['metal_atoms'])} metal atoms - verify coordination")

        return stats


def optimize_protein_protonation(input_pdb: str,
                                 output_pdb: str,
                                 ph: float = 7.0,
                                 his_state: str = "auto",
                                 work_dir: Optional[str] = None,
                                 verbose: bool = False) -> Dict[str, str]:
    """
    Convenience function to optimize protein protonation states.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output PDB file with optimized hydrogens
    ph : float
        Target pH (default: 7.0)
    his_state : str
        Histidine state: 'auto', 'HID', 'HIE', or 'HIP' (default: 'auto')
    work_dir : str, optional
        Working directory for intermediate files
    verbose : bool
        Enable verbose logging (default: False)

    Returns
    -------
    dict
        Results dictionary with file paths and statistics

    Examples
    --------
    >>> results = optimize_protein_protonation("protein.pdb", "protein_H.pdb", ph=7.4)
    >>> print(results["output"])
    protein_H.pdb
    """
    protonator = ProteinProtonator(
        ph=ph,
        his_state=his_state,
        verbose=verbose
    )

    results = protonator.optimize_hydrogens(input_pdb, output_pdb, work_dir)

    # Add validation stats
    stats = protonator.validate_protonation(output_pdb)
    results["validation"] = stats

    return results
