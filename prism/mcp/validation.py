"""Input/output validation and topology inspection tools."""

import os
import re
import json

from ._common import _StdoutToStderr, _ensure_prism_importable, logger


def _validate_ligand(ligand_path: str) -> dict:
    """Validate a ligand file (MOL2, SDF, or PDB)."""
    lig = {
        "path": ligand_path,
        "valid": False,
        "warnings": [],
    }

    if not os.path.exists(ligand_path):
        lig["warnings"].append(f"File not found: {ligand_path}")
        return lig

    if os.path.getsize(ligand_path) == 0:
        lig["warnings"].append("File is empty")
        return lig

    ext = os.path.splitext(ligand_path)[1].lower()
    lig["format"] = ext.lstrip(".")

    try:
        with open(ligand_path, "r") as f:
            content = f.read()

        if ext == ".mol2":
            if "@<TRIPOS>ATOM" not in content:
                lig["warnings"].append("Missing @<TRIPOS>ATOM section in MOL2 file")
            else:
                atom_section = content.split("@<TRIPOS>ATOM")[1]
                next_section = atom_section.split("@<TRIPOS>")[0] if "@<TRIPOS>" in atom_section else atom_section
                atom_lines = [l for l in next_section.strip().splitlines() if l.strip()]
                lig["atoms"] = len(atom_lines)
                lig["valid"] = len(atom_lines) > 0

        elif ext == ".sdf":
            lines = content.splitlines()
            if len(lines) >= 4:
                counts_line = lines[3].strip()
                match = re.match(r"(\d+)\s+(\d+)", counts_line)
                if match:
                    lig["atoms"] = int(match.group(1))
                    lig["valid"] = True
                else:
                    lig["warnings"].append(f"Cannot parse SDF counts line: '{counts_line}'")
            else:
                lig["warnings"].append("SDF file too short (fewer than 4 lines)")
            if "$$$$" not in content:
                lig["warnings"].append("Missing $$$$ end marker in SDF file")

        elif ext == ".pdb":
            atom_count = sum(1 for line in content.splitlines() if line.startswith(("ATOM", "HETATM")))
            lig["atoms"] = atom_count
            lig["valid"] = atom_count > 0
            if atom_count == 0:
                lig["warnings"].append("No ATOM/HETATM records in ligand PDB")

        else:
            lig["warnings"].append(f"Unrecognized ligand format: '{ext}'. Supported: .mol2, .sdf, .pdb")

    except Exception as e:
        lig["warnings"].append(f"Error reading ligand file: {e}")

    return lig


def _parse_topology_molecules(topol_path: str) -> list:
    """Parse the [ molecules ] section from a GROMACS topology file."""
    molecules = []
    in_molecules = False
    try:
        with open(topol_path, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith(";"):
                    continue
                if stripped.startswith("["):
                    section = stripped.strip("[] \t").lower()
                    in_molecules = section == "molecules"
                    continue
                if in_molecules:
                    if stripped.startswith("#"):
                        continue
                    parts = stripped.split()
                    if len(parts) >= 2:
                        try:
                            molecules.append({"name": parts[0], "count": int(parts[1])})
                        except ValueError:
                            pass
    except Exception:
        pass
    return molecules


def register(mcp):

    @mcp.tool()
    def validate_input_files(protein_path: str, ligand_path: str) -> str:
        """Pre-build validation of protein PDB and ligand files.

        Checks both files for common issues before building, saving time by
        catching problems early. Validates file existence, format, and content.

        For proteins: counts atoms, residues, chains, histidines, checks END record.
        For ligands: detects format (MOL2/SDF/PDB), counts atoms, checks structure.

        Call this BEFORE build_system() to verify inputs are correct.

        Args:
            protein_path: Absolute path to the protein PDB file.
            ligand_path: Absolute path to the ligand file (MOL2, SDF, or PDB).
        """
        logger.info(f"validate_input_files: {protein_path}, {ligand_path}")

        result = {"protein": {}, "ligand": {}}

        # --- Protein PDB validation ---
        prot = result["protein"]
        prot["path"] = protein_path
        prot["valid"] = False
        prot["warnings"] = []

        if not os.path.exists(protein_path):
            prot["warnings"].append(f"File not found: {protein_path}")
            result["ligand"] = _validate_ligand(ligand_path)
            return json.dumps(result, indent=2)

        if os.path.getsize(protein_path) == 0:
            prot["warnings"].append("File is empty")
            result["ligand"] = _validate_ligand(ligand_path)
            return json.dumps(result, indent=2)

        try:
            atom_count = 0
            chains = set()
            residue_ids = set()
            residue_names = set()
            has_end = False

            standard_aa = {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                "THR", "TRP", "TYR", "VAL",
                # Common protonation variants
                "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "ASH", "GLH", "LYN",
            }

            with open(protein_path, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        atom_count += 1
                        chain_id = line[21] if len(line) > 21 else " "
                        chains.add(chain_id.strip() or " ")
                        resname = line[17:20].strip()
                        resseq = line[22:26].strip()
                        residue_ids.add((chain_id, resseq, resname))
                        residue_names.add(resname)
                    elif line.startswith(("END", "ENDMDL")):
                        has_end = True

            if atom_count == 0:
                prot["warnings"].append("No ATOM/HETATM records found in PDB file")
            else:
                prot["valid"] = True

            his_residues = {
                (c, s) for c, s, r in residue_ids
                if r in ("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP")
            }

            ignore_res = {"HOH", "WAT", "SOL", "NA", "CL", "K", "MG", "CA", "ZN"}
            non_standard = residue_names - standard_aa - ignore_res
            if non_standard:
                prot["warnings"].append(
                    f"Non-standard residues found: {sorted(non_standard)}. "
                    f"These may be ligands, cofactors, or modified residues."
                )

            if not has_end:
                prot["warnings"].append("No END/ENDMDL record found (may be truncated)")

            prot["atoms"] = atom_count
            prot["residues"] = len(residue_ids)
            prot["chains"] = sorted(chains)
            prot["his_count"] = len(his_residues)
            prot["has_end_record"] = has_end

        except Exception as e:
            prot["warnings"].append(f"Error reading PDB: {e}")

        # --- Ligand validation ---
        result["ligand"] = _validate_ligand(ligand_path)

        return json.dumps(result, indent=2)

    @mcp.tool()
    def validate_build_output(output_dir: str, build_mode: str = "md") -> str:
        """Post-build validation — verify that PRISM produced all expected output files.

        Checks for required files, parses topology for molecule summary, and reads
        system atom count from the coordinate file.

        Call this AFTER a build completes to confirm everything is in order.

        Args:
            output_dir: The PRISM output directory path (same as passed to build_system).
            build_mode: The build type to validate. Default: "md".
                Options: "md", "pmf", "rest2", "mmpbsa".
        """
        logger.info(f"validate_build_output: {output_dir} (mode={build_mode})")

        result = {
            "valid": False,
            "build_mode": build_mode,
            "files_found": [],
            "files_missing": [],
            "topology_summary": {},
            "system_atoms": None,
            "warnings": [],
        }

        if not os.path.isdir(output_dir):
            result["warnings"].append(f"Output directory not found: {output_dir}")
            return json.dumps(result, indent=2)

        # --- Define expected files per build mode ---
        expected_files = {}

        if build_mode == "md":
            expected_files = {
                "solv_ions.gro": "GMX_PROLIG_MD/solv_ions.gro",
                "topol.top": "GMX_PROLIG_MD/topol.top",
                "localrun.sh": "GMX_PROLIG_MD/localrun.sh",
                "em.mdp": "mdps/em.mdp",
                "nvt.mdp": "mdps/nvt.mdp",
                "npt.mdp": "mdps/npt.mdp",
                "md.mdp": "mdps/md.mdp",
            }
        elif build_mode == "pmf":
            expected_files = {
                "solv_ions.gro": "GMX_PROLIG_PMF/solv_ions.gro",
                "topol.top": "GMX_PROLIG_PMF/topol.top",
                "smd_run.sh": "GMX_PROLIG_PMF/smd_run.sh",
                "umbrella_run.sh": "GMX_PROLIG_PMF/umbrella_run.sh",
                "smd.mdp": "mdps/smd.mdp",
                "umbrella.mdp": "mdps/umbrella.mdp",
            }
        elif build_mode == "rest2":
            expected_files = {
                "rest2_run.sh": "GMX_PROLIG_REST2/rest2_run.sh",
            }
        elif build_mode == "mmpbsa":
            expected_files = {
                "mmpbsa_run.sh": "GMX_PROLIG_MMPBSA/mmpbsa_run.sh",
                "mmpbsa.in": "GMX_PROLIG_MMPBSA/mmpbsa.in",
                "topol.top": "GMX_PROLIG_MMPBSA/topol.top",
            }
        else:
            result["warnings"].append(f"Unknown build_mode '{build_mode}'. Options: md, pmf, rest2, mmpbsa")
            return json.dumps(result, indent=2)

        # --- Check each expected file ---
        for label, rel_path in expected_files.items():
            full_path = os.path.join(output_dir, rel_path)
            if os.path.exists(full_path) and os.path.getsize(full_path) > 0:
                result["files_found"].append(label)
            else:
                result["files_missing"].append(label)

        # --- REST2: also check for replica topologies ---
        if build_mode == "rest2":
            rest2_dir = os.path.join(output_dir, "GMX_PROLIG_REST2")
            if os.path.isdir(rest2_dir):
                replica_tops = sorted(f for f in os.listdir(rest2_dir) if re.match(r"topol_r\d+\.top", f))
                if replica_tops:
                    result["files_found"].append(f"replica_topologies ({len(replica_tops)} files)")
                else:
                    result["files_missing"].append("topol_r*.top (replica topologies)")

        # --- Parse topology for molecule summary ---
        topol_path = None
        for rel in expected_files.values():
            if rel.endswith("topol.top"):
                topol_path = os.path.join(output_dir, rel)
                break

        if topol_path and os.path.exists(topol_path):
            molecules = _parse_topology_molecules(topol_path)
            if molecules:
                result["topology_summary"]["molecules"] = molecules

        # --- Read atom count from .gro file ---
        gro_path = None
        for rel in expected_files.values():
            if rel.endswith(".gro"):
                gro_path = os.path.join(output_dir, rel)
                break

        if gro_path and os.path.exists(gro_path):
            try:
                with open(gro_path, "r") as f:
                    f.readline()  # title
                    atom_line = f.readline().strip()
                    result["system_atoms"] = int(atom_line)
            except (ValueError, StopIteration):
                result["warnings"].append("Could not parse atom count from .gro file")

        # --- Check prism_config.yaml ---
        config_path = os.path.join(output_dir, "prism_config.yaml")
        if os.path.exists(config_path):
            result["files_found"].append("prism_config.yaml")
        else:
            result["warnings"].append("prism_config.yaml not found (non-critical)")

        # --- Final validity ---
        result["valid"] = len(result["files_missing"]) == 0

        return json.dumps(result, indent=2)

    @mcp.tool()
    def check_topology(topology_path: str) -> str:
        """Deep inspection of a GROMACS topology file for troubleshooting.

        Parses the topology to extract force field, water model, molecule list,
        and all #include directives. Checks whether included files exist on disk.
        Flags common issues like missing solvent, ions, or ligand entries.

        Use this tool when a build produces unexpected results or when the user
        reports simulation errors related to topology.

        Args:
            topology_path: Absolute path to the GROMACS topology file (topol.top).
        """
        logger.info(f"check_topology: {topology_path}")

        result = {
            "path": topology_path,
            "forcefield": None,
            "water_model": None,
            "molecules": [],
            "includes": [],
            "total_molecule_entries": 0,
            "warnings": [],
        }

        if not os.path.exists(topology_path):
            result["warnings"].append(f"Topology file not found: {topology_path}")
            return json.dumps(result, indent=2)

        topol_dir = os.path.dirname(os.path.abspath(topology_path))

        # Resolve GMXLIB for include path checking
        gmxlib = os.environ.get("GMXLIB", "")
        gmxdata = os.environ.get("GMXDATA", "")
        gmx_share = ""
        if gmxdata:
            gmx_share = os.path.join(gmxdata, "top")

        def _include_exists(inc_file: str) -> bool:
            if os.path.exists(os.path.join(topol_dir, inc_file)):
                return True
            if gmxlib and os.path.exists(os.path.join(gmxlib, inc_file)):
                return True
            if gmx_share and os.path.exists(os.path.join(gmx_share, inc_file)):
                return True
            return False

        try:
            with open(topology_path, "r") as f:
                content = f.read()

            # --- Extract #include directives ---
            include_pattern = re.compile(r'#include\s+"([^"]+)"')
            includes = include_pattern.findall(content)

            for inc in includes:
                exists = _include_exists(inc)
                result["includes"].append({"file": inc, "exists": exists})
                if not exists:
                    result["warnings"].append(f"Include file not found: {inc}")

            # --- Detect force field from first include ---
            if includes:
                first_inc = includes[0]
                ff_match = re.match(r"([^/]+)\.ff/", first_inc)
                if ff_match:
                    result["forcefield"] = ff_match.group(1)

            # --- Detect water model ---
            water_models = {"tip3p", "tip4p", "tip4pew", "tip5p", "spc", "spce", "opc"}
            for inc in includes:
                basename = os.path.splitext(os.path.basename(inc))[0].lower()
                if basename in water_models:
                    result["water_model"] = basename
                    break

            # --- Parse [ molecules ] section ---
            result["molecules"] = _parse_topology_molecules(topology_path)
            result["total_molecule_entries"] = len(result["molecules"])

            # --- Check for common issues ---
            mol_names = {m["name"] for m in result["molecules"]}

            if not any(n in mol_names for n in ("SOL", "HOH", "WAT", "TIP3", "SPC")):
                result["warnings"].append("System has no solvent (no SOL/HOH/WAT entry in [ molecules ])")

            if not any(n in mol_names for n in ("NA", "CL", "K", "Na+", "Cl-", "SOD", "CLA")):
                result["warnings"].append("System has no ions (no NA/CL entry in [ molecules ])")

            if not any(
                n for n in mol_names
                if n not in (
                    "Protein", "Protein_chain_A", "SOL", "HOH", "WAT",
                    "NA", "CL", "K", "Na+", "Cl-", "SOD", "CLA", "TIP3", "SPC",
                ) and not n.startswith("Protein")
            ):
                result["warnings"].append("No ligand found in topology (no non-protein/solvent/ion molecule)")

        except Exception as e:
            result["warnings"].append(f"Error parsing topology: {e}")

        return json.dumps(result, indent=2)
