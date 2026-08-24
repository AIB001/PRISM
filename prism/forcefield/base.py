#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Base class for force field generators in PRISM
"""

from abc import ABC, abstractmethod
import os
from pathlib import Path


class ForceFieldGeneratorBase(ABC):
    """Abstract base class for force field generators"""

    def __init__(self, ligand_path, output_dir, overwrite=False):
        """
        Initialize the force field generator

        Parameters:
        -----------
        ligand_path : str
            Path to the ligand file
        output_dir : str
            Directory where output files will be stored
        overwrite : bool
            Whether to overwrite existing files
        """
        self.ligand_path = os.path.abspath(ligand_path)
        self.output_dir = os.path.abspath(output_dir)
        self.overwrite = overwrite
        self.ligand_name = Path(ligand_path).stem

        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

    @abstractmethod
    def run(self):
        """
        Run the force field generation workflow

        Returns:
        --------
        str : Path to the output directory containing the generated files
        """
        pass

    @abstractmethod
    def get_output_dir_name(self):
        """
        Get the name of the output directory for this force field

        Returns:
        --------
        str : Directory name (e.g., 'LIG.amb2gmx' or 'LIG.openff2gmx')
        """
        pass

    @property
    def capabilities(self):
        """
        Declare force field capabilities

        Returns:
        --------
        dict : Capability flags
            - supports_fep (bool): Whether this force field supports FEP calculations
            - supports_pmf (bool): Whether this force field supports PMF calculations
            - requires_external (bool): Whether external service/download is required
            - has_protein_ff (bool): Whether this includes protein force field
            - charge_method (str): Charge calculation method (e.g., 'AM1-BCC', 'RESP', 'EEM')
        """
        return {
            "supports_fep": True,
            "supports_pmf": True,
            "requires_external": False,
            "has_protein_ff": False,
            "charge_method": "unknown",
        }

    def count_missing_hydrogens(self):
        """How many hydrogens the ligand file is missing, or None if unknown.

        Returns 0 for a molecule that genuinely has none (hexafluorobenzene, a
        halide ion) and None when RDKit -- an optional PRISM dependency -- is
        unavailable to answer the question.
        """
        try:
            from rdkit import Chem, rdBase
        except (ImportError, ModuleNotFoundError):
            return None

        suffix = Path(self.ligand_path).suffix.lower()
        readers = {
            ".sdf": lambda p: Chem.MolFromMolFile(p, sanitize=False, removeHs=False),
            ".mol": lambda p: Chem.MolFromMolFile(p, sanitize=False, removeHs=False),
            ".mol2": lambda p: Chem.MolFromMol2File(p, sanitize=False, removeHs=False),
            ".pdb": lambda p: Chem.MolFromPDBFile(p, sanitize=False, removeHs=False),
        }
        reader = readers.get(suffix)
        if reader is None:
            return None
        try:
            with rdBase.BlockLogs():
                molecule = reader(self.ligand_path)
                if molecule is None:
                    return None
                molecule.UpdatePropertyCache(strict=False)
                Chem.FastFindRings(molecule)
                # AddHs adds exactly the hydrogens the valence model says are
                # missing, so the difference is the shortfall.
                return Chem.AddHs(molecule).GetNumAtoms() - molecule.GetNumAtoms()
        except Exception:
            return None

    def expected_hydrogen_count(self):
        """Hydrogens the finished topology should contain, or None if unknown.

        Distinct from :meth:`count_missing_hydrogens`, which reports only what
        the *input file* lacks: a correctly hydrogenated ligand is missing none
        but still needs its hydrogens to survive into the topology.
        """
        try:
            from rdkit import Chem, rdBase
        except (ImportError, ModuleNotFoundError):
            return None

        suffix = Path(self.ligand_path).suffix.lower()
        readers = {
            ".sdf": lambda p: Chem.MolFromMolFile(p, sanitize=False, removeHs=False),
            ".mol": lambda p: Chem.MolFromMolFile(p, sanitize=False, removeHs=False),
            ".mol2": lambda p: Chem.MolFromMol2File(p, sanitize=False, removeHs=False),
            ".pdb": lambda p: Chem.MolFromPDBFile(p, sanitize=False, removeHs=False),
        }
        reader = readers.get(suffix)
        if reader is None:
            return None
        try:
            with rdBase.BlockLogs():
                molecule = reader(self.ligand_path)
                if molecule is None:
                    return None
                molecule.UpdatePropertyCache(strict=False)
                Chem.FastFindRings(molecule)
                complete = Chem.AddHs(molecule)
                return sum(1 for a in complete.GetAtoms() if a.GetAtomicNum() == 1)
        except Exception:
            return None

    def require_explicit_hydrogens(self):
        """Refuse a ligand whose hydrogens are absent from the file.

        Structure-based generative models emit heavy atoms only, and so do many
        docking outputs and PDB-extracted ligands. Handing such a file to
        antechamber does not reliably fail: AM1-BCC aborts on the resulting odd
        electron count for some molecules, but for others the pipeline runs to
        completion and writes a topology with no hydrogens at all -- correct
        file layout, wrong chemistry, no warning. Catching it here costs
        milliseconds instead of a wrong production run.
        """
        missing = self.count_missing_hydrogens()
        if not missing:
            return
        raise ValueError(
            f"Ligand '{Path(self.ligand_path).name}' is missing {missing} hydrogen(s). "
            "Parameterization needs explicit hydrogens: without them the charge "
            "calculation either aborts on an odd electron count or silently "
            "produces a hydrogen-free topology.\n"
            "  Add them first, keeping the heavy-atom coordinates:\n"
            "    from rdkit import Chem\n"
            "    mol = Chem.AddHs(Chem.MolFromMolFile('ligand.sdf'), addCoords=True)\n"
            "    Chem.MolToMolFile(mol, 'ligand_h.sdf')\n"
            "  PRISM's generation module already writes this as candidates_h.sdf.\n"
            "  For a physiological protonation state, use a pH-aware tool "
            "(Dimorphite-DL, OpenBabel -p 7.4) instead."
        )

    def verify_topology_hydrogens(self, lig_dir):
        """Warn when a produced topology has no hydrogens but should.

        The second half of the same failure: even with a correct input, a
        toolchain hiccup can yield a heavy-atom-only topology that still passes
        every file-existence check.
        """
        expected = self.expected_hydrogen_count()
        itp = os.path.join(lig_dir, "LIG.itp")
        if not os.path.exists(itp):
            return
        in_atoms = False
        total = hydrogens = 0
        with open(itp, "r") as handle:
            for line in handle:
                stripped = line.strip()
                if stripped.startswith("["):
                    in_atoms = stripped.replace(" ", "") == "[atoms]"
                    continue
                if not in_atoms or not stripped or stripped.startswith(";"):
                    continue
                fields = stripped.split()
                if len(fields) >= 5 and fields[0].isdigit():
                    total += 1
                    if fields[1].lower().startswith("h") or fields[4].upper().startswith("H"):
                        hydrogens += 1
        if total and hydrogens == 0 and expected:
            raise ValueError(
                f"Generated ligand topology has {total} atoms and no hydrogens "
                f"where {expected} were expected "
                f"({itp}). The parameterization silently dropped them; the "
                "resulting system would have wrong masses and charges. "
                "Re-run with a ligand file that carries explicit hydrogens."
            )

    def check_required_files(self, output_dir):
        """
        Check if all required files exist in the output directory

        Parameters:
        -----------
        output_dir : str
            Directory to check

        Returns:
        --------
        bool : True if all required files exist
        """
        required_files = ["LIG.gro", "LIG.itp", "LIG.top", "atomtypes_LIG.itp", "posre_LIG.itp"]

        for filename in required_files:
            filepath = os.path.join(output_dir, filename)
            if not os.path.exists(filepath):
                return False

        return True
