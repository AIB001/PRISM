#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GAFF2 force field generator wrapper for PRISM

GAFF2 is an improved version of GAFF with:
- Better torsion parameters
- Improved charge models
- Enhanced coverage of chemical space
- Better agreement with QM calculations
"""

import os

# Import the GAFF generator as base
try:
    from .gaff import GAFFForceFieldGenerator
except ImportError:
    from gaff import GAFFForceFieldGenerator

# Import color utilities
try:
    from ..utils.colors import print_success, print_info, print_warning, success, path, number
except ImportError:
    try:
        from prism.utils.colors import print_success, print_info, print_warning, success, path, number
    except ImportError:
        # Fallback if colors not available
        def print_success(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}✓ {x}")

        def print_info(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}ℹ {x}")

        def print_warning(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}⚠ {x}")

        def success(x):
            return f"✓ {x}"

        def path(x):
            return x

        def number(x):
            return x


class GAFF2ForceFieldGenerator(GAFFForceFieldGenerator):
    """GAFF2 force field generator wrapper - extends GAFF with gaff2 atom types"""

    def __init__(self, ligand_path, output_dir, overwrite=False, charge_mode="bcc"):
        """
        Initialize GAFF2 force field generator.

        Parameters
        ----------
        ligand_path : str
            Path to ligand file (MOL2 or SDF)
        output_dir : str
            Output directory for generated files
        overwrite : bool
            Whether to overwrite existing files
        charge_mode : str
            Charge calculation method: 'bcc' for AM1-BCC (default), 'gas' for fast
            gas-phase charges (useful when RESP charges will be applied later)
        """
        super().__init__(ligand_path, output_dir, overwrite, charge_mode)

        # Override the print message to show GAFF2
        print(f"\n{'='*60}")
        print_info("Using GAFF2 Force Field (improved version)")
        print("GAFF2 features:")
        print(f"  - {success('Better torsion parameters')}")
        print(f"  - {success('Improved charge models')}")
        print(f"  - {success('Enhanced coverage of chemical space')}")
        print(f"{'='*60}")

    def run(self):
        """Run the GAFF2 force field generation workflow"""
        print(f"\n{'='*60}")
        print("Starting GAFF2 Force Field Generation")
        print(f"{'='*60}")

        try:
            # Check if output already exists BEFORE trying to generate parameters
            lig_dir = os.path.join(self.output_dir, "LIG.amb2gmx")
            if os.path.exists(lig_dir) and not self.overwrite:
                if self.check_required_files(lig_dir):
                    print_info(f"Using cached GAFF2 parameters from: {path(lig_dir)}")
                    print("All required files found:")
                    for f in sorted(os.listdir(lig_dir)):
                        print_success(f, prefix="  -")
                    print_warning("Use --overwrite to regenerate")
                    return lig_dir

            # Generate AMBER parameters (using GAFF2)
            print("  Generating AMBER parameters...")
            ff_dir = self.generate_amber_parameters()
            print_success("  AMBER parameters generated")

            # Find acpype output
            print("  Converting to GROMACS format...")
            ff_files = self.find_acpype_output(ff_dir)

            # Standardize to LIG naming
            lig_dir = self.standardize_to_LIG(ff_files)
            print_success("  Converted to GROMACS format")

            # Cleanup
            self.cleanup_temp_files()

            print(f"\n{'='*60}")
            print_success("GAFF2 force field generation completed!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {path(lig_dir)}")
            print("\nGenerated files:")
            for f in sorted(os.listdir(lig_dir)):
                print(f"  - {f}")

            return lig_dir

        except Exception as e:
            print(f"\nError during GAFF2 force field generation: {e}")
            raise

    def _get_atom_type_flag(self):
        """Get the atom type flag for antechamber - GAFF2 version"""
        return "gaff2"

    def _get_leaprc(self):
        """Get the appropriate leaprc file - GAFF2 version"""
        return "leaprc.gaff2"

    def _process_mol2_format(self, amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file):
        """Process MOL2 format files using GAFF2"""
        print("Processing MOL2 format with GAFF2 atom types...")

        # Get the ff_dir from amber_mol2 path
        ff_dir = os.path.dirname(amber_mol2)

        # Clean up any existing antechamber temp files before starting
        self._clean_antechamber_temp_files(ff_dir)

        # Copy input file to working directory to avoid path issues
        import shutil

        input_mol2 = os.path.join(ff_dir, f"input_{self.ligand_name}.mol2")
        shutil.copy2(self.ligand_path, input_mol2)

        print(f"Generating charges using {self.charge_mode.upper()} method with GAFF2 atom types...")
        # Use -dr no to skip acdoctor validation, because mol2 files with mixed
        # aromatic + double bond types can cause false "Weird atomic valence" errors
        cmd = [
            "antechamber",
            "-i",
            os.path.basename(input_mol2),
            "-fi",
            "mol2",
            "-o",
            os.path.basename(amber_mol2),
            "-fo",
            "mol2",
            "-c",
            self.charge_mode,  # Use configured charge mode
            "-s",
            "2",
            "-at",
            self._get_atom_type_flag(),  # Use GAFF2 atom types
            "-dr",
            "no",
        ]
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])

        try:
            self.run_command(cmd, cwd=ff_dir)
        except Exception as e:
            if self.charge_mode == "bcc":
                print(f"\nAM1-BCC charge generation failed. Attempting fallback with gas phase charges...")
                # Try with simpler gas phase charges as fallback
                cmd_fallback = [
                    "antechamber",
                    "-i",
                    os.path.basename(input_mol2),
                    "-fi",
                    "mol2",
                    "-o",
                    os.path.basename(amber_mol2),
                    "-fo",
                    "mol2",
                    "-c",
                    "gas",
                    "-s",
                    "2",
                    "-at",
                    self._get_atom_type_flag(),  # Use GAFF2 atom types
                    "-dr",
                    "no",
                ]
                if self.net_charge != 0:
                    cmd_fallback.extend(["-nc", str(self.net_charge)])
                self.run_command(cmd_fallback, cwd=ff_dir)
            else:
                raise  # Re-raise if not using bcc mode

        print("Generating prep file with GAFF2...")
        # Generate prep file from the charged mol2 file
        # Use -dr no to skip acdoctor validation (same reason as above)
        cmd = [
            "antechamber",
            "-i",
            os.path.basename(amber_mol2),
            "-fi",
            "mol2",
            "-o",
            os.path.basename(prep_file),
            "-fo",
            "prepi",
            "-dr",
            "no",
            "-at",
            self._get_atom_type_flag(),  # Use GAFF2 atom types
        ]
        # Add net charge if specified
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])

        try:
            self.run_command(cmd, cwd=ff_dir)
        except Exception as e:
            print("\nPrep file generation failed. Trying with explicit charge retention...")
            # Alternative: explicitly tell antechamber to read charges from mol2
            cmd_alt = [
                "antechamber",
                "-i",
                os.path.basename(amber_mol2),
                "-fi",
                "mol2",
                "-o",
                os.path.basename(prep_file),
                "-fo",
                "prepi",
                "-dr",
                "no",
                "-c",
                "rc",  # Read charges from input file
                "-at",
                self._get_atom_type_flag(),  # Use GAFF2 atom types
            ]
            if self.net_charge != 0:
                cmd_alt.extend(["-nc", str(self.net_charge)])
            self.run_command(cmd_alt, cwd=ff_dir)

        print("Generating GAFF2 force field parameters...")
        self.run_command(
            [
                "parmchk2",
                "-i",
                os.path.basename(prep_file),
                "-f",
                "prepi",
                "-o",
                os.path.basename(frcmod_file),
                "-a",
                "Y",  # Print all force field parameters
                "-s",
                "2",  # GAFF2 parameter set
            ],
            cwd=ff_dir,
        )

        print("Creating AMBER topology with GAFF2...")
        self._create_topology_with_tleap(amber_mol2, frcmod_file, prmtop_file, rst7_file)

        # Clean up intermediate input file
        try:
            os.remove(input_mol2)
        except:
            pass

    def _create_topology_with_tleap(self, amber_mol2, frcmod_file, prmtop_file, rst7_file):
        """Create AMBER topology using tleap with GAFF2"""
        ff_dir = os.path.dirname(amber_mol2)
        tleap_input = os.path.join(ff_dir, f"tleap_{self.ligand_name}.in")

        with open(tleap_input, "w") as f:
            f.write(f"""source leaprc.protein.ff14SB
source {self._get_leaprc()}

LIG = loadmol2 {os.path.basename(amber_mol2)}
loadamberparams {os.path.basename(frcmod_file)}

saveamberparm LIG {os.path.basename(prmtop_file)} {os.path.basename(rst7_file)}

quit
""")

        self.run_command(["tleap", "-f", os.path.basename(tleap_input)], cwd=ff_dir)

    def _process_sdf_format_direct(self, amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file):
        """Process SDF format files with direct antechamber conversion using GAFF2"""
        print(f"Processing SDF format (direct method) with {self.charge_mode.upper()} charges and GAFF2 atom types...")

        ff_dir = os.path.dirname(amber_mol2)

        cmd = [
            "antechamber",
            "-i",
            self.ligand_path,
            "-fi",
            "sdf",
            "-o",
            amber_mol2,
            "-fo",
            "mol2",
            "-c",
            self.charge_mode,  # Use configured charge mode
            "-s",
            "2",
            "-at",
            self._get_atom_type_flag(),  # Use GAFF2 atom types
        ]
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])
        self.run_command(cmd, cwd=ff_dir)

        # Continue with prep file generation
        self._process_mol2_intermediate(amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file)

    def _process_mol2_intermediate(
        self, amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file, charge_method=None
    ):
        """Process intermediate MOL2 file with GAFF2

        Note: When processing MOL2 from SDF (already has BCC charges), we don't specify
        a charge method (-c flag) and let antechamber handle it automatically. This matches
        the behavior of _process_mol2_format() which also omits the -c flag.
        """
        print("Generating prep file from charged MOL2 with GAFF2...")
        if charge_method:
            print_warning(f"Ignoring charge_method='{charge_method}' for GAFF2 MOL2 processing (charges already set).")

        ff_dir = os.path.dirname(amber_mol2)

        # Verify input file exists
        if not os.path.exists(amber_mol2):
            raise FileNotFoundError(f"Input MOL2 file not found: {amber_mol2}")

        # Build command WITHOUT -c flag (let antechamber handle charges)
        cmd = [
            "antechamber",
            "-i",
            os.path.basename(amber_mol2),
            "-fi",
            "mol2",
            "-o",
            os.path.basename(prep_file),
            "-fo",
            "prepi",
            "-at",
            self._get_atom_type_flag(),  # Use GAFF2 atom types
        ]

        # Add net charge if specified
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])

        try:
            print(f"Attempting prep file generation without explicit charge method...")
            self.run_command(cmd, cwd=ff_dir)
        except Exception as e:
            print(f"\nPrep file generation failed: {e}")
            print("Trying with explicit charge retention (-c rc)...")
            # Fallback: explicitly tell antechamber to read charges from mol2
            cmd_alt = [
                "antechamber",
                "-i",
                os.path.basename(amber_mol2),
                "-fi",
                "mol2",
                "-o",
                os.path.basename(prep_file),
                "-fo",
                "prepi",
                "-c",
                "rc",  # Read charges from input file
                "-at",
                self._get_atom_type_flag(),  # Use GAFF2 atom types
            ]
            if self.net_charge != 0:
                cmd_alt.extend(["-nc", str(self.net_charge)])
            self.run_command(cmd_alt, cwd=ff_dir)

        print("Generating GAFF2 force field parameters...")
        self.run_command(
            [
                "parmchk2",
                "-i",
                os.path.basename(prep_file),
                "-f",
                "prepi",
                "-o",
                os.path.basename(frcmod_file),
                "-s",
                "2",  # GAFF2 parameter set
            ],
            cwd=ff_dir,
        )

        print("Creating AMBER topology with GAFF2...")
        self._create_topology_with_tleap(amber_mol2, frcmod_file, prmtop_file, rst7_file)
