"""Environment and dependency checking tools."""

import json
import importlib.util

from ._common import _StdoutToStderr, _ensure_prism_importable, logger


def register(mcp):

    @mcp.tool()
    def check_dependencies() -> str:
        """Check if all required dependencies for PRISM are installed.

        Checks:
        - GROMACS (molecular dynamics engine, required)
        - pdbfixer (protein structure preparation, recommended)
        - AmberTools / antechamber (for GAFF/GAFF2 ligand force fields)
        - OpenFF toolkit (for OpenFF ligand force field)
        - PROPKA (for pKa-based protonation prediction of ionizable residues)
        - Gaussian g16 (for high-precision RESP charge calculation, optional)
        - MDTraj (for trajectory analysis)
        - RDKit (for chemical informatics)

        Call this tool first to verify the user's environment before building.

        Returns a JSON object with each dependency name mapped to true/false.
        """
        logger.info("Checking PRISM dependencies...")
        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                import shutil

                deps = {
                    "gromacs": False,
                    "pdbfixer": False,
                    "antechamber": False,
                    "openff": False,
                    "propka": False,
                    "gaussian_g16": False,
                    "mdtraj": False,
                    "rdkit": False,
                }

                # GROMACS
                try:
                    from prism.utils.environment import GromacsEnvironment
                    GromacsEnvironment()
                    deps["gromacs"] = True
                except Exception:
                    pass

                # pdbfixer
                deps["pdbfixer"] = shutil.which("pdbfixer") is not None

                # AmberTools (antechamber)
                deps["antechamber"] = shutil.which("antechamber") is not None

                # OpenFF
                deps["openff"] = importlib.util.find_spec("openff.toolkit") is not None

                # PROPKA
                deps["propka"] = importlib.util.find_spec("propka") is not None

                # Gaussian g16
                deps["gaussian_g16"] = shutil.which("g16") is not None

                # MDTraj
                deps["mdtraj"] = importlib.util.find_spec("mdtraj") is not None

                # RDKit
                deps["rdkit"] = importlib.util.find_spec("rdkit") is not None

            logger.info(f"Dependencies: {deps}")
            return json.dumps(deps, indent=2)

        except Exception as e:
            logger.error(f"check_dependencies failed: {e}")
            return json.dumps({"error": str(e)})

    @mcp.tool()
    def list_forcefields() -> str:
        """List all available force fields for protein-ligand system building.

        Returns two categories:
        1. Protein force fields: detected from the user's GROMACS installation.
           Recommended: amber14sb (well-validated, modern AMBER force field).
        2. Ligand force fields: supported by PRISM for small molecule parameterization.
           Recommended: gaff2 (improved GAFF with better torsion parameters).

        Common protein + ligand combinations:
        - amber14sb + gaff2  (recommended default)
        - amber99sb + gaff   (classic, widely used)
        - amber19sb + gaff2  (latest AMBER, requires OPC water: --water opc)
        - charmm36  + cgenff (CHARMM ecosystem)

        Use this tool to help the user choose appropriate force fields.
        """
        logger.info("Listing available force fields...")
        try:
            with _StdoutToStderr():
                _ensure_prism_importable()

                # Protein force fields from GROMACS
                protein_ffs = []
                try:
                    from prism.utils.environment import GromacsEnvironment
                    env = GromacsEnvironment()
                    protein_ffs = env.list_force_fields()
                except Exception:
                    protein_ffs = ["(GROMACS not detected - cannot list protein force fields)"]

                # Ligand force fields - hardcoded from PRISM's actual support
                ligand_ffs = {
                    "gaff": "GAFF force field (AmberTools). Classic, widely used.",
                    "gaff2": "GAFF2 force field (AmberTools). Improved torsion parameters. RECOMMENDED.",
                    "openff": "Open Force Field (openff-toolkit). Modern, data-driven.",
                    "cgenff": "CGenFF (CHARMM General FF). Requires web-downloaded files from cgenff.com. CLI only (use --forcefield-path).",
                    "opls": "OPLS-AA (via LigParGen server). Requires internet.",
                    "mmff": "MMFF-based (via SwissParam server). Requires internet.",
                    "match": "MATCH (via SwissParam server). Requires internet.",
                    "hybrid": "Hybrid MMFF-MATCH (via SwissParam server). Requires internet.",
                }

            result = {
                "protein_forcefields": protein_ffs,
                "ligand_forcefields": ligand_ffs,
                "recommended": {
                    "protein": "amber14sb",
                    "ligand": "gaff2",
                    "water": "tip3p",
                },
            }
            return json.dumps(result, indent=2, ensure_ascii=False)

        except Exception as e:
            logger.error(f"list_forcefields failed: {e}")
            return json.dumps({"error": str(e)})
