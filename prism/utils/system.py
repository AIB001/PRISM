# !filepath: prism/utils/system.py

# !/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM System Builder - Handles GROMACS system setup, solvation, and ion addition.
"""

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple

# Import color utilities
try:
    from .colors import print_success, print_warning, print_info, success, warning, path, number
except ImportError:
    try:
        from prism.utils.colors import print_success, print_warning, print_info, success, warning, path, number
    except ImportError:
        # Fallback
        def print_success(x, **kwargs): print(f"✓ {x}")
        def print_warning(x, **kwargs): print(f"⚠ {x}")
        def print_info(x, **kwargs): print(f"ℹ {x}")
        def success(x): return f"✓ {x}"
        def warning(x): return f"⚠ {x}"
        def path(x): return x
        def number(x): return x


class SystemBuilder:
    """
    Builds a complete protein-ligand system for GROMACS simulations.

    This class encapsulates the GROMACS commands needed to process a protein,
    combine it with a ligand, create a simulation box, solvate it, and add ions.
    It is designed to be robust, especially in handling the output of GROMACS
    tools to ensure consistency between coordinate and topology files.
    """

    def __init__(self, config: Dict, output_dir: str, overwrite: bool = False,
                 pmf_mode: bool = False, box_extension: Optional[Tuple[float, float, float]] = None):
        """
        Initializes the SystemBuilder.

        Args:
            config: A dictionary containing the simulation configuration.
            output_dir: The root directory for all output files.
            overwrite: If True, overwrite existing files.
            pmf_mode: If True, build system for PMF calculations (uses GMX_PROLIG_PMF).
            box_extension: Tuple of (x, y, z) extension values in nm for PMF box.
                          Only used when pmf_mode=True. Default: (0, 0, 2.0)
        """
        self.config = config
        self.output_dir = Path(output_dir)
        self.overwrite = overwrite
        self.gmx_command = self.config.get("general", {}).get("gmx_command", "gmx")

        # PMF mode settings
        self.pmf_mode = pmf_mode
        self.box_extension = box_extension if box_extension else (0.0, 0.0, 2.0)

        # Set model directory based on mode
        if pmf_mode:
            self.model_dir = self.output_dir / "GMX_PROLIG_PMF"
        else:
            self.model_dir = self.output_dir / "GMX_PROLIG_MD"
        self.model_dir.mkdir(exist_ok=True)

    def _run_command(self, command: list, work_dir: str, input_str: Optional[str] = None) -> Tuple[str, str]:
        """
        Executes a shell command and handles errors.

        Args:
            command: The command to execute as a list of strings.
            work_dir: The directory in which to run the command.
            input_str: A string to be passed to the command's stdin.

        Returns:
            A tuple containing the stdout and stderr of the command.

        Raises:
            RuntimeError: If the command returns a non-zero exit code.
        """
        cmd_str = ' '.join(map(str, command))
        print(f"Executing in {work_dir}: {cmd_str}")

        process = subprocess.Popen(
            command,
            cwd=work_dir,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate(input=input_str)

        if process.returncode != 0:
            print("--- STDOUT ---")
            print(stdout)
            print("--- STDERR ---")
            print(stderr)
            raise RuntimeError(f"Command failed with exit code {process.returncode}: {cmd_str}")

        return stdout, stderr

    def build(self, cleaned_protein_path: str, lig_ff_dirs, ff_idx: int, water_idx: int,
              ff_info: dict = None, water_info: dict = None, nter: str = None, cter: str = None) -> str:
        """
        Runs the full system building workflow with support for multiple ligands.

        Args:
            cleaned_protein_path: Path to the cleaned protein PDB file.
            lig_ff_dirs: Path to directory containing ligand force field files (str),
                        or list of paths for multiple ligands, or None for protein-only.
            ff_idx: The index of the protein force field to use.
            water_idx: The index of the water model to use.
            ff_info: Dictionary containing force field info (name, dir, path)
            nter: N-terminal type for pdb2gmx (e.g., 'NH3+', 'NH2'). Default: auto-select
            cter: C-terminal type for pdb2gmx (e.g., 'COO-', 'COOH'). Default: auto-select
            water_info: Dictionary containing water model info (name, label, description)

        Returns:
            The path to the directory containing the final model files.
        """
        try:
            print_info("Building GROMACS Model")

            if (self.model_dir / "solv_ions.gro").exists() and not self.overwrite:
                print(f"Final model file solv_ions.gro already exists in {self.model_dir}, skipping model building.")
                return str(self.model_dir)

            # Store water model name for use in solvation
            self.water_model_name = water_info.get('name', 'tip3p') if water_info else 'tip3p'
            # Store force field info for terminal atom fixing
            self.ff_info = ff_info

            # Normalize lig_ff_dirs to always be a list (or None for protein-only)
            if lig_ff_dirs is None:
                lig_ff_list = None
            elif isinstance(lig_ff_dirs, str):
                # Legacy single ligand support
                lig_ff_list = [lig_ff_dirs]
            elif isinstance(lig_ff_dirs, list):
                lig_ff_list = lig_ff_dirs
            else:
                raise TypeError(f"lig_ff_dirs must be None, str, or list, got {type(lig_ff_dirs)}")

            fixed_pdb = self._fix_protein(cleaned_protein_path)

            # Apply PROPKA renaming AFTER pdbfixer (pdbfixer may revert HIS names)
            protonation_method = self.config.get('protonation', {}).get('method', 'gromacs')
            if protonation_method == 'propka':
                self._apply_propka_renaming(fixed_pdb)

            protein_gro, topol_top = self._generate_topology(fixed_pdb, ff_idx, water_idx, ff_info, water_info, nter, cter)

            # Process ligands if provided
            if lig_ff_list is not None and len(lig_ff_list) > 0:
                self._fix_topology_multi_ligands(topol_top, lig_ff_list)
                complex_gro = self._combine_protein_multi_ligands(protein_gro, lig_ff_list)
            else:
                # Protein-only system
                print("\nSkipping ligand processing (protein-only system)")
                complex_gro = protein_gro

            boxed_gro = self._create_box(complex_gro)
            solvated_gro = self._solvate(boxed_gro, topol_top)
            self._add_ions(solvated_gro, topol_top)

            print_success("System building completed")
            print(f"System files are in {self.model_dir}")
            return str(self.model_dir)

        except (RuntimeError, FileNotFoundError) as e:
            print(f"\nError during system build: {e}")
            import traceback
            traceback.print_exc()
            return ""

    def _fix_protein(self, cleaned_protein: str) -> str:
        """Fix missing atoms in protein using pdbfixer."""
        print("\n  Step 1: Fixing protein with pdbfixer...")
        fixed_pdb = self.model_dir / "fixed_clean_protein.pdb"

        if fixed_pdb.exists() and not self.overwrite:
            print(f"Fixed protein PDB already exists at {fixed_pdb}, skipping pdbfixer.")
            return str(fixed_pdb)

        # Check if pdbfixer is available
        if not shutil.which("pdbfixer"):
            print("Warning: pdbfixer command not found. Skipping protein fixing step.")
            print("RECOMMENDED INSTALLATION:")
            print("  mamba install -c conda-forge pdbfixer")
            print("  # OR: conda install -c conda-forge pdbfixer")
            shutil.copy(cleaned_protein, fixed_pdb)
        else:
            # Run pdbfixer to fix missing atoms
            # Note: pdbfixer with --keep-heterogens=all will move metals to chain 'M'
            # This is OK - we extract them from chain 'M' after pdb2gmx
            command = [
                "pdbfixer",
                cleaned_protein,
                "--add-residues",
                "--output", str(fixed_pdb),
                "--add-atoms=heavy",
                "--keep-heterogens=all",
                f"--ph={self.config['simulation']['pH']}"
            ]
            self._run_command(command, str(self.model_dir.parent))

        # Fix terminal atoms for GROMACS compatibility (OXT → O or OC1/OC2)
        print("Fixing terminal atoms for GROMACS compatibility...")
        from prism.utils.cleaner import fix_terminal_atoms
        # Get force field name if available
        ff_name = getattr(self, 'ff_info', {}).get('name') if hasattr(self, 'ff_info') and self.ff_info else None
        fix_terminal_atoms(str(fixed_pdb), str(fixed_pdb), force_field=ff_name, verbose=True)

        return str(fixed_pdb)

    def _apply_propka_renaming(self, pdb_file: str):
        """Apply PROPKA pKa-based histidine renaming to fixed PDB.

        This must run AFTER pdbfixer (which may revert HIS names back to HIS)
        and BEFORE _generate_topology / _rename_his_for_cmap.
        """
        from .protonation import PropkaProtonator

        ph = self.config.get('protonation', {}).get(
            'ph', self.config.get('simulation', {}).get('pH', 7.0))

        print(f"\n  Applying PROPKA pKa-based HIS renaming (pH {ph})...")
        protonator = PropkaProtonator(ph=ph, verbose=True)
        stats = protonator.rename_histidines(pdb_file, pdb_file)

        if stats['renamed']:
            print(f"  PROPKA: Renamed {len(stats['renamed'])} HIS residue(s):")
            for (chain, resnum), state in stats['renamed'].items():
                print(f"    Chain {chain} Residue {resnum}: HIS -> {state}")
        else:
            print("  PROPKA: No histidine residues found or needed renaming")

    def _rename_his_for_cmap(self, pdb_file: str, ff_info: dict) -> bool:
        """Rename HIS → HISE in PDB for force fields with CMAP corrections.

        GROMACS pdb2gmx keeps the original PDB residue name in the topology
        even after selecting a specific protonation state (e.g. HIE).  For
        force fields that use CMAP corrections (like amber19sb), the CMAP
        type lookup is residue-name-dependent (e.g. XC-HIE, N-HIE), so
        leaving the name as HIS causes 'Unknown cmap torsion' errors.

        Renaming HIS → HISE in the PDB triggers the r2b mapping (HISE → HIE),
        so pdb2gmx writes the correct residue name and skips the interactive
        histidine protonation prompt.

        Returns True if any residues were renamed.
        """
        if not ff_info or 'path' not in ff_info:
            return False

        ff_dir = Path(ff_info['path'])
        if not (ff_dir / "cmap.itp").exists():
            return False

        # Read and rename HIS → HIE (default HIE protonation at neutral pH)
        # PDB format (0-indexed): positions 17-19 = residue name (3 chars).
        # Replace "HIS" with "HIE" (same width) so pdb2gmx uses the [ HIE ]
        # rtp block directly and writes "HIE" as the residue name in topology.
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        renamed_residues = set()
        new_lines = []
        for line in lines:
            if (line.startswith('ATOM') or line.startswith('HETATM')) \
                    and line[17:20] == 'HIS' and line[20] == ' ':
                resnum = line[22:26].strip()
                chain = line[21].strip()
                res_id = f"{chain}:{resnum}" if chain else resnum
                renamed_residues.add(res_id)
                line = line[:17] + 'HIE' + line[20:]
            new_lines.append(line)

        if renamed_residues:
            with open(pdb_file, 'w') as f:
                f.writelines(new_lines)
            print(f"  Renamed {len(renamed_residues)} HIS residue(s) → HIE for CMAP compatibility")
            return True

        return False

    def _generate_topology(self, fixed_pdb: str, ff_idx: int, water_idx: int,
                          ff_info: dict = None, water_info: dict = None,
                          nter: str = None, cter: str = None) -> Tuple[str, str]:
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
        protonation_method = self.config.get('protonation', {}).get('method', 'gromacs')

        # Build pdb2gmx command
        command = [
            self.gmx_command, "pdb2gmx",
            "-f", fixed_pdb,
            "-o", str(protein_gro),
            "-p", str(topol_top),
            "-ignh"  # Always use -ignh to regenerate standardized hydrogens
        ]

        # Print information about protonation method
        if protonation_method == 'gromacs':
            print_info("  Using GROMACS protonation (pdb2gmx -ignh: ignore and regenerate all hydrogens)")
        elif protonation_method == 'propka':
            print_info("  Using PROPKA-predicted protonation states (residues renamed, pdb2gmx -ignh regenerates H)")
            print("  -> PROPKA's role: pKa-based prediction of histidine protonation states (HID/HIE/HIP)")
            print("  -> GROMACS pdb2gmx: Regenerates standardized hydrogen atoms with correct naming")
        else:
            # Default to gromacs method
            print_warning(f"  Unknown protonation method '{protonation_method}', defaulting to GROMACS")

        # Use -ff and -water flags if force field info is provided
        # This is more reliable than interactive menu indexing
        if ff_info and 'dir' in ff_info:
            # Remove .ff extension from directory name for -ff flag
            ff_name = ff_info['dir'].replace('.ff', '')
            command.extend(["-ff", ff_name])
            print(f"Using force field: {ff_name} (from {ff_info.get('path', 'N/A')})")
        else:
            print(f"DEBUG: Force field index = {ff_idx} (interactive menu)")

        # pdb2gmx -water flag only accepts these hardcoded model names.
        # Non-standard models (e.g. OPC, OPC3) must use interactive selection.
        PDB2GMX_WATER_MODELS = {'spc', 'spce', 'tip3p', 'tip4p', 'tip4pew', 'tip5p', 'tips3p', 'none'}

        use_water_flag = False
        if water_info and 'name' in water_info:
            water_name = water_info['name'].lower()
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
        use_terminal_selection = (nter is not None or cter is not None)

        if use_terminal_selection:
            command.append("-ter")  # Enable terminal selection
            print_info(f"  Using custom terminal types: N-term={nter or 'auto'}, C-term={cter or 'auto'}")

        # Handle histidine protonation states
        histidine_count = self._count_histidines(fixed_pdb)
        input_lines = []

        # Build interactive input: ff selection (if no -ff flag), then water selection (if no -water flag)
        use_ff_flag = bool(ff_info and 'dir' in ff_info)
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
                nter_input = self._get_terminal_menu_index(nter, 'nter', ff_info)
                if nter_input:
                    input_lines.append(nter_input)
            else:
                input_lines.append("0")  # Default selection

            if cter:
                cter_input = self._get_terminal_menu_index(cter, 'cter', ff_info)
                if cter_input:
                    input_lines.append(cter_input)
            else:
                input_lines.append("0")  # Default selection

        if histidine_count > 0:
            print(f"Found {histidine_count} histidine residue(s), defaulting to HIE protonation state.")
            # HIE is typically the first and most common choice for neutral pH.
            input_lines.extend(["1"] * histidine_count)

        input_text = '\n'.join(input_lines)+'\n' if input_lines else None

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
            print_success(f"Metals processed by pdb2gmx ({number(len(ion_files))} ion topology file(s) found)", prefix="  ✓")
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

        all_atomtype_lines = []  # Store individual atomtype lines
        seen_atomtypes = set()   # Track seen atomtype names to avoid duplicates
        all_includes = []
        all_bonded_includes = []  # Store CGenFF *_ffbonded.itp includes
        all_molecules = []
        is_cgenff = False

        # Process each ligand force field directory
        for idx, lig_ff_dir in enumerate(lig_ff_dirs, 1):
            lig_ff_path = Path(lig_ff_dir)
            lig_itp_path = lig_ff_path / "LIG.itp"
            atomtypes_itp_path = lig_ff_path / "atomtypes_LIG.itp"

            if not lig_itp_path.exists():
                raise FileNotFoundError(f"Ligand {idx} ITP file not found: {lig_itp_path}")

            # Check if this is CGenFF (has charmm36.ff subdirectory)
            charmm_ff_dir = lig_ff_path / "charmm36.ff"
            if charmm_ff_dir.exists():
                is_cgenff = True

                # Collect ligand-specific bonded parameter files (*_ffbonded.itp)
                ffbonded_files = list(charmm_ff_dir.glob("*_ffbonded.itp"))
                for ffbonded_file in ffbonded_files:
                    # Check if file has content (not empty)
                    if ffbonded_file.stat().st_size > 200:  # More than just header
                        # Build relative path from GMX_PROLIG_MD directory
                        if len(lig_ff_dirs) == 1:
                            relative_bonded_path = f"../{lig_ff_path.name}/charmm36.ff/{ffbonded_file.name}"
                        else:
                            relative_bonded_path = f"../Ligand_Forcefield/{lig_ff_path.name}/charmm36.ff/{ffbonded_file.name}"
                        all_bonded_includes.append(f'#include "{relative_bonded_path}"\n')

            # Collect atomtypes (if not CGenFF) - with deduplication
            if not is_cgenff:
                if not atomtypes_itp_path.exists():
                    raise FileNotFoundError(f"Ligand {idx} atomtypes file not found: {atomtypes_itp_path}")
                with open(atomtypes_itp_path, 'r') as f:
                    for line in f:
                        # Skip headers and comments
                        if line.strip().startswith('[') or line.strip().startswith(';') or not line.strip():
                            continue
                        # Extract atomtype name (first column)
                        atomtype_name = line.split()[0] if line.split() else None
                        if atomtype_name and atomtype_name not in seen_atomtypes:
                            seen_atomtypes.add(atomtype_name)
                            all_atomtype_lines.append(line)

            # Collect include statement - use relative path from GMX_PROLIG_MD
            # Single ligand: ../LIG.xxx2gmx/LIG.itp (no Ligand_Forcefield parent)
            # Multi-ligand: ../Ligand_Forcefield/LIG.xxx2gmx_N/LIG.itp
            if len(lig_ff_dirs) == 1:
                relative_path = f"../{lig_ff_path.name}/LIG.itp"
            else:
                relative_path = f"../Ligand_Forcefield/{lig_ff_path.name}/LIG.itp"
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

        # Read current topology
        with open(topol_path, 'r') as f:
            lines = f.readlines()

        # Check if already included to prevent duplication
        # Single ligand: check for LIG.xxx2gmx pattern
        # Multi-ligand: check for Ligand_Forcefield
        if len(lig_ff_dirs) == 1:
            if any('LIG.itp' in line and '#include' in line for line in lines):
                print("Topology already includes ligand parameters.")
                return
        else:
            if any('Ligand_Forcefield' in line for line in lines):
                print("Topology already includes multi-ligand parameters.")
                return

        # Build new topology
        new_lines = []
        molecules_added = False

        for line in lines:
            new_lines.append(line)

            # Insert atomtypes and ligand ITPs after the main forcefield include
            if '#include' in line and '.ff/forcefield.itp' in line:
                if is_cgenff:
                    # CGenFF: Include bonded parameters first, then LIG.itp files
                    print_info("  Detected CGenFF force field - using integrated CHARMM atomtypes")

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
            elif molecules_added and line.strip() and not line.strip().startswith(';'):
                # Insert all LIG entries before the first molecule (e.g., Protein_chain_A)
                for mol_entry in all_molecules:
                    new_lines.insert(-1, mol_entry)
                molecules_added = False  # Prevent adding them again

        with open(topol_path, 'w') as f:
            f.writelines(new_lines)

        print_success(f"Topology file updated with {number(len(lig_ff_dirs))} ligand parameters")

    def _combine_protein_multi_ligands(self, protein_gro: str, lig_ff_dirs: list) -> str:
        """Combines protein and multiple ligand GRO files into a complex."""
        print(f"\n  Step 4: Combining protein with {number(len(lig_ff_dirs))} ligand(s)...")
        complex_gro = self.model_dir / "pro_lig.gro"

        if complex_gro.exists() and not self.overwrite:
            print("Using existing complex.gro file.")
            return str(complex_gro)

        # Read protein coordinates
        with open(protein_gro, "r") as f_prot:
            prot_lines = f_prot.readlines()

        # Read all ligand coordinates
        all_lig_lines = []
        total_lig_atoms = 0

        for idx, lig_ff_dir in enumerate(lig_ff_dirs, 1):
            lig_gro_path = Path(lig_ff_dir) / "LIG.gro"
            if not lig_gro_path.exists():
                raise FileNotFoundError(f"Ligand {idx} GRO file not found: {lig_gro_path}")

            with open(lig_gro_path, "r") as f_lig:
                lig_lines = f_lig.readlines()

            # Extract coordinates (skip header, atom count, and box vectors)
            lig_coords = lig_lines[2:-1]

            # Renumber residue numbers and atom numbers for each ligand
            renumbered_coords = []
            for line in lig_coords:
                if len(line) >= 20:
                    # Update residue number (columns 0-4): each ligand gets its own residue number
                    new_res_num = idx
                    # Update atom number (columns 15-19): sequential numbering
                    current_atom_num = int(line[15:20].strip())
                    new_atom_num = total_lig_atoms + current_atom_num
                    # Reconstruct the line
                    new_line = f"{new_res_num:5d}" + line[5:15] + f"{new_atom_num:5d}" + line[20:]
                    renumbered_coords.append(new_line)
                else:
                    renumbered_coords.append(line)

            total_lig_atoms += len(lig_coords)
            all_lig_lines.extend(renumbered_coords)
            print(f"  Ligand {number(idx)}: {number(len(lig_coords))} atoms")

        # Calculate total atoms
        prot_atom_count = int(prot_lines[1].strip())
        total_atoms = prot_atom_count + total_lig_atoms

        # Write combined GRO file
        with open(complex_gro, "w") as f_out:
            f_out.write("Complex Protein-Ligands\n")
            f_out.write(f"{total_atoms:5d}\n")

            # Write all ligand coordinates first (to match topology molecule order)
            f_out.writelines(all_lig_lines)

            # Write protein coordinates, renumbering both residue and atom numbers
            prot_coords = prot_lines[2:-1]
            num_ligands = len(lig_ff_dirs)

            for line in prot_coords:
                if len(line) >= 20:
                    # Update residue number (columns 0-4): shift by number of ligands
                    current_res_num = int(line[:5].strip())
                    new_res_num = current_res_num + num_ligands
                    # Update atom number (columns 15-19): shift by total ligand atoms
                    current_atom_num = int(line[15:20].strip())
                    new_atom_num = current_atom_num + total_lig_atoms
                    # Reconstruct the line
                    new_line = f"{new_res_num:5d}" + line[5:15] + f"{new_atom_num:5d}" + line[20:]
                    f_out.write(new_line)
                else:
                    f_out.write(line)

            # Write box vectors from protein file
            f_out.write(prot_lines[-1])

        print_success(f"Combined system: {number(prot_atom_count)} protein atoms + {number(total_lig_atoms)} ligand atoms = {number(total_atoms)} total")
        return str(complex_gro)

    def _create_box(self, complex_gro: str) -> str:
        """Creates the simulation box."""
        print("\n  Step 5: Creating simulation box...")
        boxed_gro = self.model_dir / "pro_lig_newbox.gro"

        if boxed_gro.exists() and not self.overwrite:
            print("Using existing boxed.gro file.")
            return str(boxed_gro)

        box_cfg = self.config['box']

        if self.pmf_mode:
            # PMF mode: Create rectangular box, then extend ONLY Z direction
            print(f"  PMF Mode: Creating rectangular box with Z-axis extension for pulling")
            print(f"  Z extension: {self.box_extension[2]:.1f} nm")

            # First create rectangular (triclinic) box - NOT cubic!
            # This ensures X and Y dimensions fit the protein tightly
            temp_boxed = self.model_dir / "temp_boxed.gro"
            command = [
                self.gmx_command, "editconf",
                "-f", complex_gro,
                "-o", str(temp_boxed),
                "-bt", "triclinic",  # Use triclinic (rectangular) for PMF, not cubic
                "-d", str(box_cfg['distance'])
            ]
            if box_cfg.get('center', True):
                command.append("-c")

            self._run_command(command, str(self.model_dir))

            # Read box vectors and extend ONLY Z direction (keep X, Y unchanged)
            self._extend_box_z(str(temp_boxed), str(boxed_gro))

            # Clean up temp file
            if temp_boxed.exists():
                temp_boxed.unlink()
        else:
            # Normal mode: Standard box creation
            command = [
                self.gmx_command, "editconf",
                "-f", complex_gro,
                "-o", str(boxed_gro),
                "-bt", box_cfg['shape'],
                "-d", str(box_cfg['distance'])
            ]
            if box_cfg.get('center', True):
                command.append("-c")

            self._run_command(command, str(self.model_dir))

        return str(boxed_gro)

    def _extend_box_z(self, input_gro: str, output_gro: str):
        """
        Extend box in Z direction for PMF calculations.

        For SMD pulling:
        - After alignment, ligand is above protein in +Z direction
        - Ligand will be pulled further in +Z direction
        - Extra space is added at the +Z end of the box (no atom translation)
        - This ensures the pulling space is in the correct direction
        """
        with open(input_gro, 'r') as f:
            lines = f.readlines()

        if len(lines) < 3:
            raise ValueError(f"Invalid GRO file: {input_gro}")

        # Parse box vectors from last line
        box_line = lines[-1].strip()
        box_parts = box_line.split()

        if len(box_parts) >= 3:
            # Box vectors: v1x v2y v3z [v1y v1z v2x v2z v3x v3y]
            v1x = float(box_parts[0])
            v2y = float(box_parts[1])
            v3z = float(box_parts[2])

            # For PMF mode: ONLY extend Z direction, keep X and Y unchanged
            new_v1x = v1x  # X unchanged
            new_v2y = v2y  # Y unchanged
            new_v3z = v3z + self.box_extension[2]  # Only Z extended

            print(f"    Original box: {v1x:.3f} x {v2y:.3f} x {v3z:.3f} nm")
            print(f"    Extended box: {new_v1x:.3f} x {new_v2y:.3f} x {new_v3z:.3f} nm (Z +{self.box_extension[2]:.1f} nm)")
            print(f"    Extra space added at +Z end for pulling direction")

            # For SMD: Do NOT translate atoms - keep them at original positions
            # The extra space will be at the +Z end where the ligand will be pulled
            new_lines = []
            new_lines.append(lines[0])  # Title
            new_lines.append(lines[1])  # Atom count

            # Copy all atom lines without modification
            for line in lines[2:-1]:
                new_lines.append(line)

            # Write new box vectors (only Z changed)
            if len(box_parts) > 3:
                # Preserve off-diagonal elements if present
                new_box_line = f"   {new_v1x:.5f}   {new_v2y:.5f}   {new_v3z:.5f}"
                for i in range(3, len(box_parts)):
                    new_box_line += f"   {box_parts[i]}"
                new_box_line += "\n"
            else:
                new_box_line = f"   {new_v1x:.5f}   {new_v2y:.5f}   {new_v3z:.5f}\n"

            new_lines.append(new_box_line)

            # Write output file
            with open(output_gro, 'w') as f:
                f.writelines(new_lines)

            print(f"    Box extended successfully")
        else:
            raise ValueError(f"Cannot parse box vectors from: {box_line}")

    def _solvate(self, boxed_gro: str, topol_top: str) -> str:
        """Solvates the system."""
        print("\n  Step 6: Solvating the system...")
        solvated_gro = self.model_dir / "solv.gro"

        if solvated_gro.exists() and not self.overwrite:
            print("Using existing solv.gro file.")
            return str(solvated_gro)

        # Select appropriate water model coordinate file
        # Map water model names to their coordinate files
        water_coord_files = {
            'tip3p': 'spc216.gro',  # TIP3P uses SPC geometry
            'tip3p_original': 'spc216.gro',
            'spc': 'spc216.gro',
            'spce': 'spc216.gro',
            'opc3': 'spc216.gro',  # OPC3 is 3-point, uses SPC geometry
            'tip4p': 'tip4p.gro',  # TIP4P has its own file with virtual sites
            'tip4pew': 'tip4p.gro',  # TIP4P/Ew uses same geometry
            'opc': 'tip4p.gro',   # OPC is 4-point with virtual site
            'tip5p': 'tip5p.gro'  # TIP5P has its own file
        }

        water_model = getattr(self, 'water_model_name', 'tip3p').lower()
        water_file = water_coord_files.get(water_model, 'spc216.gro')

        print(f"Using water model: {water_model} (coordinate file: {water_file})")

        command = [
            self.gmx_command, "solvate",
            "-cp", boxed_gro,
            "-cs", water_file,
            "-o", str(solvated_gro),
            "-p", topol_top
        ]
        self._run_command(command, str(self.model_dir))
        return str(solvated_gro)

    def _add_ions(self, solvated_gro: str, topol_top: str):
        """
        Adds ions to neutralize the system and achieve a target salt concentration.
        This method uses a single, robust call to `gmx genion`.
        """
        print("\n  Step 7: Adding ions...")
        ions_gro = self.model_dir / "solv_ions.gro"

        if ions_gro.exists() and not self.overwrite:
            print("Using existing solv_ions.gro file.")
            return

        ions_cfg = self.config['ions']
        temp_tpr = self.model_dir / "ions.tpr"

        # A minimal mdp is needed for grompp
        ions_mdp_path = self.model_dir / "ions.mdp"
        with open(ions_mdp_path, "w") as f:
            f.write("integrator=steep\nnsteps=0")

        grompp_cmd = [
            self.gmx_command, "grompp",
            "-f", str(ions_mdp_path),
            "-c", solvated_gro,
            "-p", topol_top,
            "-o", str(temp_tpr),
            "-maxwarn", "5"  # Allow some warnings
        ]

        # Run grompp with improper dihedral error recovery.
        # Some AMBER force fields (e.g. amber14sb) are missing improper dihedral
        # type definitions that pdb2gmx generates. Rather than guessing which
        # impropers to remove, we let GROMACS tell us exactly which lines fail.
        max_retries = 3
        grompp_success = False
        for attempt in range(max_retries + 1):
            cmd_str = ' '.join(map(str, grompp_cmd))
            print(f"Executing in {self.model_dir}: {cmd_str}")

            result = subprocess.Popen(
                grompp_cmd,
                cwd=str(self.model_dir),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            stdout, stderr = result.communicate()

            if result.returncode == 0:
                grompp_success = True
                break

            # Check for fixable improper dihedral errors
            if "No default Per. Imp. Dih. types" in stderr and attempt < max_retries:
                fixed = self._fix_improper_from_grompp_errors(stderr)
                if fixed > 0:
                    print(f"\n  Fixed {fixed} improper dihedral(s) without FF parameters, retrying grompp (attempt {attempt + 2}/{max_retries + 1})...")
                    # Clean up failed tpr before retry
                    if temp_tpr.exists():
                        temp_tpr.unlink()
                    continue

            # Non-recoverable error or nothing left to fix
            print("--- STDOUT ---")
            print(stdout)
            print("--- STDERR ---")
            print(stderr)
            raise RuntimeError(f"Command failed with exit code {result.returncode}: {cmd_str}")

        if not grompp_success:
            raise RuntimeError("grompp failed after all retry attempts")

        genion_cmd = [
            self.gmx_command, "genion",
            "-s", str(temp_tpr),
            "-o", str(ions_gro),
            "-p", topol_top,
            "-pname", ions_cfg['positive_ion'],
            "-nname", ions_cfg['negative_ion'],
        ]
        if ions_cfg.get('neutral', True):
            genion_cmd.append("-neutral")
        if ions_cfg.get('concentration', 0) > 0:
            genion_cmd.extend(["-conc", str(ions_cfg['concentration'])])

        # Provide "SOL" (the typical name for water) as input to stdin
        stdout, _ = self._run_command(genion_cmd, str(self.model_dir), input_str="SOL")

        # --- CRITICAL STEP: Parse genion output and fix topology ---
        self._update_topology_molecules(topol_top, stdout)

        # Clean up temporary files
        temp_tpr.unlink()
        ions_mdp_path.unlink()

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
        pattern = re.compile(r'^\s*([A-Z0-9_]+)\s+:\s+(\d+)')
        parsing = False
        
        for line in genion_stdout.splitlines():
            if 'Number of molecules:' in line:
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
            
            if any(ion in content for ion in ['NA ', 'CL ', 'K ', 'BR ', 'CA ']):
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
            elif start_index != -1 and line.strip().startswith('['):
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
        for i in range(start_index+1, end_index):
            line = top_lines[i].strip()
            if not line or line.startswith(';'): continue
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
        final_lines = top_lines[:start_index]+new_molecules_section+top_lines[end_index:]

        with open(topol_path, "w") as f:
            f.writelines(final_lines)

        print("Topology `[ molecules ]` section has been successfully updated.")

    def _get_terminal_menu_index(self, terminal_type: str, terminus: str, ff_info: dict = None) -> str:
        """
        Map terminal type name to menu selection index for pdb2gmx.

        This is a simplified mapping for common terminal types.
        Different force fields may have different menus.

        Args:
            terminal_type: Terminal type name (e.g., 'NH3+', 'COO-')
            terminus: 'nter' or 'cter'
            ff_info: Force field information dict

        Returns:
            Menu index as string
        """
        # Common N-terminal types (CHARMM36 and AMBER force fields)
        nter_map = {
            'None': '0',
            'NH3+': '1',  # Common protonated N-terminus
            'NH3': '1',   # Alias
            'NH2': '2',   # Neutral N-terminus
            'GLY-NH3+': '3',
            'PRO-NH2+': '4',
        }

        # Common C-terminal types
        cter_map = {
            'None': '0',
            'COO-': '1',  # Common deprotonated C-terminus (charged)
            'COO': '1',   # Alias
            'COOH': '2',  # Neutral C-terminus
        }

        if terminus == 'nter':
            return nter_map.get(terminal_type, '0')
        elif terminus == 'cter':
            return cter_map.get(terminal_type, '0')
        else:
            return '0'

    def _count_histidines(self, pdb_file: str) -> int:
        """Counts unique histidine residues that need interactive protonation selection.

        Only counts generic 'HIS' residues (3-char name with space at position 20).
        Residues already named HIE/HID/HIP (by PROPKA) or HISE/HISD/HISH (by
        _rename_his_for_cmap) don't trigger pdb2gmx prompts and are excluded.
        """
        his_residues = set()
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[17:21] == 'HIS ':
                    resnum = line[22:26].strip()
                    chain = line[21].strip()
                    res_id = f"{chain}:{resnum}" if chain else resnum
                    his_residues.add(res_id)
        return len(his_residues)

    def _extract_metals_from_pdb(self, pdb_file: str) -> list:
        """
        Extract metal ion information from PDB file.

        IMPORTANT: With the new cleaner.py implementation, metals are converted
        from HETATM to ATOM records, so we need to check BOTH record types.

        However, we must EXCLUDE protein residues like CYM, HID, HIE which are
        also ATOM records but are amino acids, not metals.

        Returns:
            List of dicts with keys: residue_name, atom_name, chain, resnum, coords
        """
        metals = []
        # True metal ions (must match cleaner.py definitions)
        metal_names = {'ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'CO', 'NI', 'CD', 'HG',
                      'ZN2', 'MG2', 'CA2', 'FE2', 'FE3', 'CU2', 'MN2'}

        # Amino acid residues that should NEVER be treated as metals
        # These include all standard + special protonation states
        protein_residues = {
            # Standard amino acids
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            # Special protonation states
            'HID', 'HIE', 'HIP',  # Histidine variants
            'CYM', 'CYX',  # Cysteine variants (deprotonated, disulfide)
            'LYN',  # Lysine neutral
            'ASH', 'GLH',  # Protonated acidic residues
            # Terminal variants
            'ACE', 'NME', 'NH2'
        }

        with open(pdb_file, 'r') as f:
            for line in f:
                # Check both ATOM and HETATM (metals may be converted to ATOM by cleaner)
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    residue_name = line[17:20].strip().upper()
                    atom_name = line[12:16].strip().upper()

                    # CRITICAL: Skip if this is a protein residue
                    if residue_name in protein_residues:
                        continue

                    # Check if this is a metal (by residue name or atom name)
                    if residue_name in metal_names or atom_name in metal_names:
                        chain = line[21].strip()
                        resnum = line[22:26].strip()
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        metals.append({
                            'residue_name': residue_name,
                            'atom_name': atom_name,
                            'chain': chain,
                            'resnum': resnum,
                            'coords': (x, y, z)
                        })

        return metals

    def _add_metals_to_topology(self, pdb_file: str, gro_file: str, topol_file: str):
        """
        Add metal ions to topology and coordinates after pdb2gmx.

        pdb2gmx only processes ATOM records (protein), so we need to manually
        add HETATM records (metals) to both topology and coordinates.
        """
        metals = self._extract_metals_from_pdb(pdb_file)

        if not metals:
            return  # No metals to add

        print(f"\nAdding {len(metals)} metal ion(s) to topology and coordinates...")

        # Add metals to topology
        with open(topol_file, 'r') as f:
            topol_lines = f.readlines()

        # Find [ molecules ] section
        molecules_idx = -1
        for i, line in enumerate(topol_lines):
            if line.strip() == '[ molecules ]':
                molecules_idx = i
                break

        if molecules_idx == -1:
            print("Warning: Could not find [ molecules ] section in topology")
            return

        # Count metals by residue name
        from collections import Counter
        metal_counts = Counter(m['residue_name'] for m in metals)

        # Insert metal molecule entries after [ molecules ] header
        insert_idx = molecules_idx + 1
        # Skip comment lines
        while insert_idx < len(topol_lines) and topol_lines[insert_idx].strip().startswith(';'):
            insert_idx += 1

        for metal_name, count in metal_counts.items():
            topol_lines.insert(insert_idx, f"{metal_name:<20} {count}\n")
            insert_idx += 1
            print(f"  Added {metal_name}: {count}")

        with open(topol_file, 'w') as f:
            f.writelines(topol_lines)

        # Add metals to GRO file
        with open(gro_file, 'r') as f:
            gro_lines = f.readlines()

        # Parse current atom count
        atom_count = int(gro_lines[1].strip())

        # Convert metal coordinates to GRO format and append before box line
        metal_gro_lines = []
        for i, metal in enumerate(metals, start=1):
            # GRO format: resnr resname atomname atomnr x y z
            resnum = metal['resnum']
            resname = metal['residue_name']
            atomname = metal['atom_name']
            atom_num = atom_count + i
            x, y, z = metal['coords']
            # Convert Angstroms to nm
            x_nm, y_nm, z_nm = x/10.0, y/10.0, z/10.0

            # GRO format line
            line = f"{int(resnum):5d}{resname:5s}{atomname:>5s}{atom_num:5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n"
            metal_gro_lines.append(line)

        # Update atom count
        gro_lines[1] = f"{atom_count + len(metals):5d}\n"

        # Insert metal lines before box line (last line)
        gro_lines = gro_lines[:-1] + metal_gro_lines + [gro_lines[-1]]

        with open(gro_file, 'w') as f:
            f.writelines(gro_lines)

        print(f"Successfully added {len(metals)} metal ion(s) to system")

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
            r'ERROR\s+\d+\s+\[file\s+(.+?),\s+line\s+(\d+)\]:\s*\n\s*No default Per\. Imp\. Dih\. types',
            re.MULTILINE
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

            with open(filepath, 'r') as f:
                lines = f.readlines()

            # Backup before first modification
            backup = str(filepath) + '.backup'
            if not Path(backup).exists():
                with open(backup, 'w') as f:
                    f.writelines(lines)

            # Comment out the specific failing lines
            for line_num in sorted(line_nums):
                idx = line_num - 1  # Convert to 0-based
                if 0 <= idx < len(lines):
                    lines[idx] = f"; REMOVED by PRISM (no FF params): {lines[idx]}"
                    total_fixed += 1

            with open(filepath, 'w') as f:
                f.writelines(lines)

        return total_fixed


if __name__ == '__main__':
    print("This script is intended to be used as a module within the PRISM package.")
