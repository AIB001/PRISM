#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Topology generation and fixing utilities for SystemBuilder.
"""

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    from ..colors import print_success, print_warning, print_info, path, number
except ImportError:
    from prism.utils.colors import print_success, print_warning, print_info, path, number

try:
    from prism.forcefield.adapters import CGenFFAdapter
except ImportError:
    # Fallback if adapters not available
    class CGenFFAdapter:
        @classmethod
        def detect(cls, lig_ff_path):
            return (lig_ff_path / "charmm36.ff").exists()

        @classmethod
        def find_charmm_ff_dir(cls, lig_ff_path):
            charmm_ff = lig_ff_path / "charmm36.ff"
            return charmm_ff if charmm_ff.exists() else None


_PARAM_SECTIONS = {"bondtypes", "pairtypes", "angletypes", "dihedraltypes", "atomtypes"}


def _strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def _section_name(line: str) -> Optional[str]:
    stripped = line.strip().lower()
    if stripped.startswith("[") and stripped.endswith("]"):
        return stripped.strip("[] ").lower()
    return None


def _read_itp_sections(file_path: Path) -> Dict[str, List[str]]:
    sections: Dict[str, List[str]] = {}
    current: Optional[str] = None
    for raw_line in file_path.read_text().splitlines():
        section = _section_name(raw_line)
        if section:
            current = section
            sections.setdefault(section, [])
            continue
        if current:
            sections[current].append(raw_line)
    return sections


def _parse_atomtypes(file_path: Path) -> set[str]:
    atomtypes = set()
    for raw_line in _read_itp_sections(file_path).get("atomtypes", []):
        stripped = _strip_comment(raw_line)
        if not stripped:
            continue
        tokens = stripped.split()
        if tokens:
            atomtypes.add(tokens[0])
    return atomtypes


def _make_param_key(section: str, tokens: List[str]) -> Optional[Tuple[str, ...]]:
    if section == "bondtypes":
        if len(tokens) < 3:
            return None
        return (tokens[0], tokens[1], tokens[2])
    if section == "pairtypes":
        if len(tokens) < 3:
            return None
        return (tokens[0], tokens[1], tokens[2])
    if section == "angletypes":
        if len(tokens) < 4:
            return None
        return (tokens[0], tokens[1], tokens[2], tokens[3])
    if section == "dihedraltypes":
        if len(tokens) < 5:
            return None
        return (tokens[0], tokens[1], tokens[2], tokens[3], tokens[4])
    return None


def _parse_param_entries(file_paths: Iterable[Path], section: str) -> List[Dict[str, object]]:
    entries: List[Dict[str, object]] = []
    for file_path in file_paths:
        if not file_path.exists():
            continue
        for raw_line in _read_itp_sections(file_path).get(section, []):
            stripped = _strip_comment(raw_line)
            if not stripped:
                continue
            tokens = stripped.split()
            key = _make_param_key(section, tokens)
            if key is None:
                continue
            entries.append({"key": key, "line": raw_line.rstrip() + "\n", "source": str(file_path)})
    return entries


def _match_param_key(section: str, required: Tuple[str, ...], candidate: Tuple[str, ...]) -> bool:
    if section in {"bondtypes", "pairtypes"}:
        return required == candidate or required == (candidate[1], candidate[0], candidate[2])
    if section == "angletypes":
        rev = (candidate[2], candidate[1], candidate[0], candidate[3])
        return required == candidate or required == rev
    if section == "dihedraltypes":
        if required[-1] != candidate[-1]:
            return False

        def _matches(req: Tuple[str, ...], cand: Tuple[str, ...]) -> bool:
            return all(c == "X" or r == c for r, c in zip(req[:4], cand[:4]))

        rev = (candidate[3], candidate[2], candidate[1], candidate[0], candidate[4])
        return _matches(required, candidate) or _matches(required, rev)
    return False


def _entry_specificity(entry_key: Tuple[str, ...]) -> int:
    return sum(1 for token in entry_key[:-1] if token != "X")


def _extract_required_param_keys(lig_itp_path: Path) -> Dict[str, List[Tuple[str, ...]]]:
    sections = _read_itp_sections(lig_itp_path)
    atom_types: Dict[str, str] = {}
    for raw_line in sections.get("atoms", []):
        stripped = _strip_comment(raw_line)
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) >= 2:
            atom_types[tokens[0]] = tokens[1]

    required: Dict[str, List[Tuple[str, ...]]] = {
        section: [] for section in ("bondtypes", "pairtypes", "angletypes", "dihedraltypes")
    }

    for raw_line in sections.get("bonds", []):
        stripped = _strip_comment(raw_line)
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) < 3:
            continue
        key = (atom_types[tokens[0]], atom_types[tokens[1]], tokens[2])
        required["bondtypes"].append(key)

    for raw_line in sections.get("pairs", []):
        stripped = _strip_comment(raw_line)
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) < 3:
            continue
        key = (atom_types[tokens[0]], atom_types[tokens[1]], tokens[2])
        required["pairtypes"].append(key)

    for raw_line in sections.get("angles", []):
        stripped = _strip_comment(raw_line)
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) < 4:
            continue
        key = (atom_types[tokens[0]], atom_types[tokens[1]], atom_types[tokens[2]], tokens[3])
        required["angletypes"].append(key)

    for raw_line in sections.get("dihedrals", []):
        stripped = _strip_comment(raw_line)
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) < 5:
            continue
        key = (
            atom_types[tokens[0]],
            atom_types[tokens[1]],
            atom_types[tokens[2]],
            atom_types[tokens[3]],
            tokens[4],
        )
        required["dihedraltypes"].append(key)

    return required


class TopologyProcessorMixin:
    """Mixin for topology-related operations."""

    def _generate_topology(
        self,
        fixed_pdb: str,
        ff_idx: int,
        water_idx: int,
        ff_info: dict = None,
        water_info: dict = None,
        nter: str = None,
        cter: str = None,
    ) -> Tuple[str, str]:
        """Generate topology using pdb2gmx, handling histidines."""
        print("\n  Step 2: Generating topology with pdb2gmx...")

        protein_gro = self.model_dir / "pro.gro"
        topol_top = self.model_dir / "topol.top"

        if protein_gro.exists() and not self.overwrite:
            print(f"Protein topology already generated at {protein_gro}, skipping pdb2gmx.")
            return str(protein_gro), str(topol_top)

        # Rename HIS → HISE for force fields with CMAP (e.g. amber19sb)
        # This must happen before pdb2gmx to ensure correct residue names in topology
        his_renamed = self._rename_his_for_cmap(fixed_pdb, ff_info)

        # Check protonation method
        protonation_method = self.config.get("protonation", {}).get("method", "gromacs")

        # Build pdb2gmx command
        command = [
            self.gmx_command,
            "pdb2gmx",
            "-f",
            fixed_pdb,
            "-o",
            str(protein_gro),
            "-p",
            str(topol_top),
            "-ignh",  # Always use -ignh to regenerate standardized hydrogens
        ]

        # Print information about protonation method
        if protonation_method == "gromacs":
            print_info("  Using GROMACS protonation (pdb2gmx -ignh: ignore and regenerate all hydrogens)")
        elif protonation_method == "propka":
            print_info("  Using PROPKA-predicted protonation states (residues renamed, pdb2gmx -ignh regenerates H)")
            print("  -> PROPKA's role: pKa-based prediction of histidine protonation states (HID/HIE/HIP)")
            print("  -> GROMACS pdb2gmx: Regenerates standardized hydrogen atoms with correct naming")
        else:
            # Default to gromacs method
            print_warning(f"  Unknown protonation method '{protonation_method}', defaulting to GROMACS")

        # Use -ff and -water flags if force field info is provided
        # This is more reliable than interactive menu indexing
        if ff_info and "dir" in ff_info:
            # Copy the explicitly selected force field into the working directory
            # so pdb2gmx resolves exactly one match for -ff.
            import shutil

            ff_name = ff_info["dir"].replace(".ff", "")
            ff_basename = ff_info["dir"]
            local_ff_dir = self.model_dir / ff_basename

            ff_path = Path(ff_info.get("path", ""))
            if not ff_path.exists() or not ff_path.is_dir():
                print(f"Warning: Force field path not found: {ff_path}")
                print(f"Using force field from system search paths: {ff_name}")
            else:
                # Validate force field completeness before copying
                required_files = ["forcefield.itp", "ffbonded.itp", "ffnonbonded.itp"]

                # Check for ion files (support both old and new force field structures)
                # New: ions_tip3p.itp, ions_spce.itp (water-model-specific)
                # Old: ions.itp (generic)
                has_ions = False
                if water_info and "name" in water_info:
                    water_name = water_info["name"]
                    water_specific_ions = f"ions_{water_name}.itp"
                    if (ff_path / water_specific_ions).exists():
                        required_files.append(water_specific_ions)
                        has_ions = True
                    elif (ff_path / "ions.itp").exists():
                        # Old-style force field with generic ions.itp
                        required_files.append("ions.itp")
                        has_ions = True
                        print(f"  Using old-style ion file: ions.itp (for {water_name})")
                elif (ff_path / "ions.itp").exists():
                    # No water model specified, but ions.itp exists
                    required_files.append("ions.itp")
                    has_ions = True

                if not has_ions:
                    print(f"  Warning: No ion file found for water model {water_info.get('name', 'unknown')}")

                missing = [f for f in required_files if not (ff_path / f).exists()]
                if missing:
                    raise RuntimeError(
                        f"Force field is incomplete: {ff_info['name']}\n"
                        f"Missing files: {', '.join(missing)}\n"
                        f"Force field path: {ff_path}\n"
                        f"Source: {ff_info.get('source', 'unknown')}\n"
                        f"Please check your force field installation."
                    )

                # Copy if not exists, or validate and replace if incomplete
                need_copy = not local_ff_dir.exists()
                if not need_copy and self.overwrite:
                    # Check if existing copy is complete
                    missing_local = [f for f in required_files if not (local_ff_dir / f).exists()]
                    if missing_local:
                        print(f"Existing force field copy is incomplete, replacing...")
                        shutil.rmtree(local_ff_dir)
                        need_copy = True

                if need_copy:
                    shutil.copytree(ff_path, local_ff_dir)
                    source_str = ff_info.get("source", ff_path)
                    print(f"Copied force field to working directory: {local_ff_dir}")
                    print(f"  Source: {source_str}")
                else:
                    print(f"Using existing force field copy: {local_ff_dir}")

            command.extend(["-ff", ff_name])
            print(f"Using force field: {ff_name}")
        else:
            print(f"DEBUG: Force field index = {ff_idx} (interactive menu)")

        # pdb2gmx -water flag only accepts these hardcoded model names.
        # Non-standard models (e.g. OPC, OPC3) must use interactive selection.
        PDB2GMX_WATER_MODELS = {"spc", "spce", "tip3p", "tip4p", "tip4pew", "tip5p", "tips3p", "none"}

        use_water_flag = False
        if water_info and "name" in water_info:
            water_name = water_info["name"].lower()
            if water_name in PDB2GMX_WATER_MODELS:
                command.extend(["-water", water_name])
                use_water_flag = True
                print(f"Using water model: {water_name}")
            else:
                # Non-standard water model — use interactive menu selection
                command.extend(["-water", "select"])
                print(f"Using water model: {water_name} (via interactive selection, index {water_idx})")
        else:
            print(f"DEBUG: Water model index = {water_idx} (interactive menu)")

        # Handle terminal type selection
        # If nter or cter are specified, use interactive terminal selection mode
        use_terminal_selection = nter is not None or cter is not None

        if use_terminal_selection:
            command.append("-ter")  # Enable terminal selection
            print_info(f"  Using custom terminal types: N-term={nter or 'auto'}, C-term={cter or 'auto'}")

        # Handle histidine protonation states
        histidine_count = self._count_histidines(fixed_pdb)
        input_lines = []

        # Build interactive input: ff selection (if no -ff flag), then water selection (if no -water flag)
        use_ff_flag = bool(ff_info and "dir" in ff_info)
        if not use_ff_flag:
            input_lines.append(str(ff_idx))
        if not use_water_flag:
            input_lines.append(str(water_idx))

        # Add terminal type selections if specified
        if use_terminal_selection:
            # Need to provide terminal type menu selections
            # For now, we'll provide the terminal type names directly to the interactive prompts
            # This requires knowing the menu indices, which is force field dependent
            # As a workaround, we'll map common terminal names to expected indices
            if nter:
                # Map terminal name to menu selection number
                # This is a simplified approach - actual indices vary by force field
                nter_input = self._get_terminal_menu_index(nter, "nter", ff_info)
                if nter_input:
                    input_lines.append(nter_input)
            else:
                input_lines.append("0")  # Default selection

            if cter:
                cter_input = self._get_terminal_menu_index(cter, "cter", ff_info)
                if cter_input:
                    input_lines.append(cter_input)
            else:
                input_lines.append("0")  # Default selection

        if histidine_count > 0:
            print(f"Found {histidine_count} histidine residue(s), defaulting to HIE protonation state.")
            # HIE is typically the first and most common choice for neutral pH.
            input_lines.extend(["1"] * histidine_count)

        input_text = "\n".join(input_lines) + "\n" if input_lines else None

        self._run_command(command, str(self.model_dir), input_str=input_text)

        # Check if pdb2gmx already handled metals (any ion topology exists)
        # With improved PRISM cleaner, metals are converted from HETATM to ATOM
        # This allows pdb2gmx to directly recognize and include them
        #
        # pdb2gmx can generate ion topologies with different naming patterns:
        # - topol_Ion.itp (single ion chain)
        # - topol_Ion_chain_*.itp (multiple ion chains)
        ion_files = list(self.model_dir.glob("topol_Ion*.itp"))

        if ion_files:
            print_success(
                f"Metals processed by pdb2gmx ({number(len(ion_files))} ion topology file(s) found)", prefix="  ✓"
            )
            for itp_file in ion_files:
                print(f"  - {itp_file.name}")
        else:
            # Fallback: pdb2gmx didn't process metals, add them manually
            # This shouldn't happen with the new HETATM->ATOM conversion in cleaner.py
            print("\n⚠ Warning: No ion topologies detected by pdb2gmx")
            print("  Attempting manual metal addition as fallback...")
            self._add_metals_to_topology(fixed_pdb, str(protein_gro), str(topol_top))

        return str(protein_gro), str(topol_top)

    def _fix_topology_multi_ligands(self, topol_path: str, lig_ff_dirs: list):
        """Includes multiple ligand parameters into the main topology file."""
        print(f"\n  Step 3: Including {number(len(lig_ff_dirs))} ligand parameters in topology...")

        with open(topol_path, "r") as f:
            lines = f.readlines()

        all_atomtype_lines = []  # Store individual atomtype lines
        seen_atomtypes = set()  # Track seen atomtype names to avoid duplicates
        all_includes = []
        all_bonded_includes = []  # Store CGenFF *_ffbonded.itp includes
        all_molecules = []
        is_cgenff = False
        topol_dir = Path(topol_path).parent
        main_ff_dir = self._resolve_main_forcefield_dir(lines, topol_dir)
        main_bonded_files = []
        main_nonbonded_files = []
        main_atomtypes = set()
        if main_ff_dir:
            main_bonded_files = [
                p for p in [main_ff_dir / "ffbonded.itp", main_ff_dir / "ffmissingdihedrals.itp"] if p.exists()
            ]
            main_nonbonded_files = [p for p in [main_ff_dir / "ffnonbonded.itp"] if p.exists()]
            for file_path in main_nonbonded_files:
                main_atomtypes.update(_parse_atomtypes(file_path))

        # Process each ligand force field directory
        for idx, lig_ff_dir in enumerate(lig_ff_dirs, 1):
            lig_ff_path = Path(lig_ff_dir)
            lig_itp_path = lig_ff_path / "LIG.itp"
            lig_top_path = lig_ff_path / "LIG.top"
            atomtypes_itp_path = lig_ff_path / "atomtypes_LIG.itp"

            if not lig_itp_path.exists():
                raise FileNotFoundError(f"Ligand {idx} ITP file not found: {lig_itp_path}")

            # Check if this is CGenFF (has charmm*.ff subdirectory)
            charmm_ff_dir = CGenFFAdapter.find_charmm_ff_dir(lig_ff_path)
            if charmm_ff_dir:
                is_cgenff = True
                lig_include_path = lig_itp_path

                supplement_path = self._build_cgenff_parameter_supplement(
                    lig_ff_path=lig_ff_path,
                    lig_itp_path=lig_itp_path,
                    charmm_ff_dir=charmm_ff_dir,
                    main_bonded_files=main_bonded_files,
                    main_nonbonded_files=main_nonbonded_files,
                )
                if supplement_path:
                    if len(lig_ff_dirs) == 1:
                        relative_bonded_path = f"../{lig_ff_path.name}/{supplement_path.name}"
                    else:
                        relative_bonded_path = f"../Ligand_Forcefield/{lig_ff_path.name}/{supplement_path.name}"
                    all_bonded_includes.append(f'#include "{relative_bonded_path}"\n')

            # Collect atomtypes - with deduplication
            # For CGenFF, also collect atomtypes from atomtypes_LIG.itp
            if not atomtypes_itp_path.exists():
                raise FileNotFoundError(f"Ligand {idx} atomtypes file not found: {atomtypes_itp_path}")
            with open(atomtypes_itp_path, "r") as f:
                for line in f:
                    # Skip headers and comments
                    if line.strip().startswith("[") or line.strip().startswith(";") or not line.strip():
                        continue
                    # Extract atomtype name (first column)
                    atomtype_name = line.split()[0] if line.split() else None
                    if atomtype_name in main_atomtypes:
                        continue
                    if atomtype_name and atomtype_name not in seen_atomtypes:
                        seen_atomtypes.add(atomtype_name)
                        all_atomtype_lines.append(line)

            # Collect include statement - use relative path from GMX_PROLIG_MD
            # Single ligand: ../LIG.xxx2gmx/LIG.itp (no Ligand_Forcefield parent)
            # Multi-ligand: ../Ligand_Forcefield/LIG.xxx2gmx_N/LIG.itp
            lig_include_filename = lig_include_path.name
            if len(lig_ff_dirs) == 1:
                relative_path = f"../{lig_ff_path.name}/{lig_include_filename}"
            else:
                relative_path = f"../Ligand_Forcefield/{lig_ff_path.name}/{lig_include_filename}"
            all_includes.append(f'#include "{relative_path}"\n')

            # Collect molecule entry
            # Single ligand: use "LIG" for backward compatibility
            # Multiple ligands: use "LIG_N" to distinguish them
            if len(lig_ff_dirs) == 1:
                mol_name = "LIG"
            else:
                mol_name = f"LIG_{idx}"
            all_molecules.append(f"{mol_name:<20} 1  ; Ligand {idx}\n")

            print_success(f"  Processed ligand {number(idx)}: {path(lig_ff_path.name)}")

        # Check if already included to prevent duplication
        # Single ligand: check for LIG.xxx2gmx pattern
        # Multi-ligand: check for Ligand_Forcefield
        if len(lig_ff_dirs) == 1:
            if any("LIG.itp" in line and "#include" in line for line in lines):
                print("Topology already includes ligand parameters.")
                return
        else:
            if any("Ligand_Forcefield" in line for line in lines):
                print("Topology already includes multi-ligand parameters.")
                return

        # Build new topology
        new_lines = []
        molecules_added = False

        for line in lines:
            new_lines.append(line)

            # Insert atomtypes and ligand ITPs after the main forcefield include
            if "#include" in line and ".ff/forcefield.itp" in line:
                if is_cgenff:
                    # CGenFF: Include bonded parameters first, then LIG.itp files
                    print_info("  Detected CGenFF force field - using integrated CHARMM atomtypes")

                    # CRITICAL FIX: Include atomtypes for CGenFF ligands
                    # Even though CHARMM36 has atomtypes, ligand-specific atomtypes
                    # from atomtypes_LIG.itp must be included in the main topology
                    new_lines.append("\n; Include CGenFF ligand atomtypes\n")
                    new_lines.append("[ atomtypes ]\n")
                    new_lines.append(";type, bondingtype, atomic_number, mass, charge, ptype, sigma, epsilon\n")
                    # Write all unique atomtype lines
                    for atomtype_line in all_atomtype_lines:
                        new_lines.append(atomtype_line)
                    new_lines.append("\n")

                    # Include ligand-specific bonded parameters BEFORE LIG.itp files
                    if all_bonded_includes:
                        new_lines.append("\n; Include CGenFF ligand-specific bonded parameters\n")
                        for bonded_include in all_bonded_includes:
                            new_lines.append(bonded_include)

                    new_lines.append("\n; Include CGenFF ligand topologies (with CHARMM36 parameters)\n")
                else:
                    # Other force fields: Include merged atomtypes (deduplicated)
                    new_lines.append("\n; Include ligand atomtypes (merged and deduplicated from all ligands)\n")
                    new_lines.append("[ atomtypes ]\n")
                    # Add comment header
                    new_lines.append(";type, bondingtype, atomic_number, mass, charge, ptype, sigma, epsilon\n")
                    # Write all unique atomtype lines
                    for atomtype_line in all_atomtype_lines:
                        new_lines.append(atomtype_line)
                    new_lines.append("\n; Include ligand topologies\n")

                # Add all ligand includes
                for include in all_includes:
                    new_lines.append(include)

            # Add ligands to the molecules section
            if line.strip() == "[ molecules ]":
                molecules_added = True
            elif molecules_added and line.strip() and not line.strip().startswith(";"):
                # Insert all LIG entries before the first molecule (e.g., Protein_chain_A)
                for mol_entry in all_molecules:
                    new_lines.insert(-1, mol_entry)
                molecules_added = False  # Prevent adding them again

        with open(topol_path, "w") as f:
            f.writelines(new_lines)

        print_success(f"Topology file updated with {number(len(lig_ff_dirs))} ligand parameters")

    def _resolve_main_forcefield_dir(self, topol_lines: List[str], topol_dir: Path) -> Optional[Path]:
        for line in topol_lines:
            if "#include" not in line or ".ff/forcefield.itp" not in line:
                continue
            match = re.search(r'"([^"]+\.ff)/forcefield\.itp"', line)
            if not match:
                continue
            ff_dir = (topol_dir / match.group(1)).resolve()
            if ff_dir.exists():
                return ff_dir
        return None

    def _build_cgenff_parameter_supplement(
        self,
        lig_ff_path: Path,
        lig_itp_path: Path,
        charmm_ff_dir: Path,
        main_bonded_files: List[Path],
        main_nonbonded_files: List[Path],
    ) -> Optional[Path]:
        source_files: List[Path] = []
        charmm36_itp = charmm_ff_dir / "charmm36.itp"
        if charmm36_itp.exists():
            source_files.append(charmm36_itp)
        full_ffbonded = charmm_ff_dir / "ffbonded.itp"
        if full_ffbonded.exists():
            source_files.append(full_ffbonded)
        for ffbonded_file in sorted(charmm_ff_dir.glob("*_ffbonded.itp")):
            if ffbonded_file.exists() and ffbonded_file.stat().st_size > 0:
                source_files.append(ffbonded_file)

        if not source_files:
            print_warning(f"  No CGenFF bonded parameter source files found in {charmm_ff_dir}")
            return None

        required = _extract_required_param_keys(lig_itp_path)
        supplement_sections: Dict[str, List[str]] = {
            section: [] for section in ("bondtypes", "pairtypes", "angletypes", "dihedraltypes")
        }

        for section, main_files, source_pool in (
            ("bondtypes", main_bonded_files, source_files),
            ("pairtypes", main_nonbonded_files, source_files),
            ("angletypes", main_bonded_files, source_files),
            ("dihedraltypes", main_bonded_files, source_files),
        ):
            main_entries = _parse_param_entries(main_files, section)
            source_entries = _parse_param_entries(source_pool, section)
            seen_lines = set()
            for req_key in required[section]:
                if any(_match_param_key(section, req_key, entry["key"]) for entry in main_entries):
                    continue
                matches = [entry for entry in source_entries if _match_param_key(section, req_key, entry["key"])]
                if not matches:
                    continue
                best_entry = max(matches, key=lambda entry: _entry_specificity(entry["key"]))
                if best_entry["line"] not in seen_lines:
                    supplement_sections[section].append(best_entry["line"])
                    seen_lines.add(best_entry["line"])

        supplement_parts: List[str] = [
            "; Auto-generated CGenFF supplement with ligand-specific parameters\n",
            "; Generated by PRISM to avoid including a full secondary CHARMM force field.\n",
            "\n",
        ]
        wrote_any = False
        for section in ("bondtypes", "pairtypes", "angletypes", "dihedraltypes"):
            lines = supplement_sections[section]
            if not lines:
                continue
            wrote_any = True
            supplement_parts.append(f"[ {section} ]\n")
            supplement_parts.extend(lines)
            if lines and not lines[-1].endswith("\n"):
                supplement_parts.append("\n")
            supplement_parts.append("\n")

        if not wrote_any:
            return None

        supplement_path = lig_ff_path / "cgenff_bonded_LIG.itp"
        supplement_path.write_text("".join(supplement_parts))
        print_info(f"  Built CGenFF bonded supplement: {supplement_path.name}")
        return supplement_path

    def _update_topology_molecules(self, topol_path: str, genion_stdout: str):
        """
        Parses the stdout of `gmx genion` and updates the `[ molecules ]`
        section of the topology file to ensure consistency.

        Modern GROMACS versions (2024+) automatically update the topology when using -p flag,
        so we verify the update was successful rather than parsing stdout.
        """
        print("Checking topology update from genion...")

        # First try to parse traditional output format for backward compatibility
        mol_counts = {}
        pattern = re.compile(r"^\s*([A-Z0-9_]+)\s+:\s+(\d+)")
        parsing = False

        for line in genion_stdout.splitlines():
            if "Number of molecules:" in line:
                parsing = True
                continue
            if parsing:
                match = pattern.match(line)
                if match:
                    mol_name, count = match.groups()
                    mol_counts[mol_name] = int(count)
                elif not line.strip():
                    break

        # If parsing failed, check if GROMACS already updated the topology (modern versions)
        if not mol_counts:
            print("Modern GROMACS detected - topology should be auto-updated by genion -p flag.")
            # Verify the topology has ions by checking the [molecules] section
            with open(topol_path, "r") as f:
                content = f.read()

            if any(ion in content for ion in ["NA ", "CL ", "K ", "BR ", "CA "]):
                print_success("Topology updated with ions", prefix="  ✓")
                return
            else:
                print("Warning: Could not verify ion addition in topology. Manual check recommended.")
                return

        with open(topol_path, "r") as f:
            top_lines = f.readlines()

        # Find the [ molecules ] section
        start_index, end_index = -1, -1
        for i, line in enumerate(top_lines):
            if line.strip() == "[ molecules ]":
                start_index = i
            elif start_index != -1 and line.strip().startswith("["):
                end_index = i
                break
        if start_index != -1 and end_index == -1:
            end_index = len(top_lines)

        if start_index == -1:
            raise RuntimeError("Could not find [ molecules ] section in topology file.")

        # Reconstruct the section
        new_molecules_section = [top_lines[start_index]]

        # Preserve non-ion/water molecules from the original file
        original_mols = {}
        for i in range(start_index + 1, end_index):
            line = top_lines[i].strip()
            if not line or line.startswith(";"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                mol_name = parts[0]
                if mol_name not in mol_counts:
                    original_mols[mol_name] = parts[1]

        # Write preserved molecules first
        for mol, count in original_mols.items():
            new_molecules_section.append(f"{mol:<18} {count}\n")

        # Then write the updated counts from genion
        for mol, count in mol_counts.items():
            if count > 0:
                new_molecules_section.append(f"{mol:<18} {count}\n")

        # Replace the old section with the new one
        final_lines = top_lines[:start_index] + new_molecules_section + top_lines[end_index:]

        with open(topol_path, "w") as f:
            f.writelines(final_lines)

        print("Topology `[ molecules ]` section has been successfully updated.")

    def _fix_improper_from_grompp_errors(self, stderr: str) -> int:
        """
        Parse grompp error output and remove specific failing improper dihedrals.

        Uses GROMACS's own error messages to identify exactly which topology lines
        have undefined improper dihedral parameters. This is 100% correct because:
        - GROMACS uses order-dependent matching for improper dihedrals
        - Only lines that ACTUALLY fail parameter lookup are removed
        - Modified force fields that add improper parameters won't trigger errors,
          so their dihedrals are never touched

        Args:
            stderr: The stderr output from a failed grompp run

        Returns:
            Number of improper dihedrals removed
        """
        # Parse GROMACS error blocks:
        #   ERROR N [file FILENAME, line LINENUM]:
        #     No default Per. Imp. Dih. types
        error_pattern = re.compile(
            r"ERROR\s+\d+\s+\[file\s+(.+?),\s+line\s+(\d+)\]:\s*\n\s*No default Per\. Imp\. Dih\. types", re.MULTILINE
        )

        # Group errors by file
        errors_by_file = {}
        for match in error_pattern.finditer(stderr):
            filename = match.group(1)
            line_num = int(match.group(2))
            errors_by_file.setdefault(filename, set()).add(line_num)

        if not errors_by_file:
            return 0

        total_fixed = 0
        for filename, line_nums in errors_by_file.items():
            # Resolve file path (grompp reports paths relative to working dir)
            filepath = self.model_dir / filename
            if not filepath.exists():
                filepath = Path(filename)
            if not filepath.exists():
                print(f"  Warning: Cannot find {filename}, skipping")
                continue

            with open(filepath, "r") as f:
                lines = f.readlines()

            # Backup before first modification
            backup = str(filepath) + ".backup"
            if not Path(backup).exists():
                with open(backup, "w") as f:
                    f.writelines(lines)

            # Comment out the specific failing lines
            for line_num in sorted(line_nums):
                idx = line_num - 1  # Convert to 0-based
                if 0 <= idx < len(lines):
                    lines[idx] = f"; REMOVED by PRISM (no FF params): {lines[idx]}"
                    total_fixed += 1

            with open(filepath, "w") as f:
                f.writelines(lines)

        return total_fixed
