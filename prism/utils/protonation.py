#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein protonation state prediction using PROPKA.

This module uses PROPKA to predict per-residue pKa values and determine
histidine protonation states (HID/HIE/HIP) based on the protein's 3D
environment and a target pH. The renamed residues are then passed to
GROMACS pdb2gmx, which regenerates standardized hydrogens.
"""

import logging
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)


class PropkaProtonator:
    """
    Predict histidine protonation states using PROPKA pKa calculations.

    PROPKA is an empirical pKa predictor that estimates per-residue pKa
    values from the protein's 3D structure. This class uses those predictions
    to rename HIS residues to HID/HIE/HIP before GROMACS pdb2gmx processing.
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
        self._check_propka_available()

    def _check_propka_available(self):
        """Check if PROPKA is installed and importable."""
        try:
            import propka
            logger.info(f"PROPKA version: {propka.__version__}")
        except ImportError:
            raise ImportError(
                "PROPKA is not available.\n"
                "Install with: pip install propka\n"
                "Or install with PRISM: pip install -e .[protonation]"
            )

    def predict_his_states(self, pdb_file: str) -> Dict[tuple, str]:
        """
        Predict histidine protonation states from pKa values.

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

    def rename_histidines(self, pdb_file: str, output_file: str) -> Dict:
        """
        Predict protonation states and rename HIS residues in PDB file.

        Parameters
        ----------
        pdb_file : str
            Path to input PDB file
        output_file : str
            Path to output PDB file with renamed residues

        Returns
        -------
        dict
            Statistics with keys: 'total_his', 'renamed', 'states'
        """
        his_states = self.predict_his_states(pdb_file)

        stats = {
            'total_his': len(his_states),
            'renamed': {},
            'states': his_states,
        }

        if not his_states:
            # No histidines found; just copy the file
            if pdb_file != output_file:
                import shutil
                shutil.copy2(pdb_file, output_file)
            return stats

        # Read PDB and rename HIS residues
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            if (line.startswith('ATOM') or line.startswith('HETATM')) \
                    and line[17:20] == 'HIS' and line[20] == ' ':
                chain = line[21].strip()
                resnum = line[22:26].strip()
                key = (chain, resnum)

                if key in his_states:
                    new_name = his_states[key]
                    line = line[:17] + f"{new_name:3s}" + line[20:]
                    stats['renamed'][key] = new_name

            new_lines.append(line)

        with open(output_file, 'w') as f:
            f.writelines(new_lines)

        return stats
