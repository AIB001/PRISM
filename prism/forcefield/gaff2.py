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


class GAFF2ForceFieldGenerator(GAFFForceFieldGenerator):
    """GAFF2 force field generator wrapper - extends GAFF with gaff2 atom types"""

    def __init__(self, ligand_path, output_dir, overwrite=False):
        """Initialize GAFF2 force field generator"""
        super().__init__(ligand_path, output_dir, overwrite)

        # Override the print message to show GAFF2
        print(f"\n{'='*60}")
        print("NOTE: Using GAFF2 Force Field (improved version)")
        print("GAFF2 features:")
        print("  - Better torsion parameters")
        print("  - Improved charge models")
        print("  - Enhanced coverage of chemical space")
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
                    print(f"\nUsing cached GAFF2 force field parameters from: {lig_dir}")
                    print("All required files found:")
                    for f in sorted(os.listdir(lig_dir)):
                        print(f"  - {f}")
                    print("\n(Use --overwrite to regenerate)")
                    return lig_dir

            # Generate AMBER parameters (using GAFF2)
            ff_dir = self.generate_amber_parameters()

            # Find acpype output
            ff_files = self.find_acpype_output(ff_dir)

            # Standardize to LIG naming
            lig_dir = self.standardize_to_LIG(ff_files)

            # Cleanup
            self.cleanup_temp_files()

            print(f"\n{'='*60}")
            print("GAFF2 force field generation completed successfully!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {lig_dir}")
            print("\nGenerated files:")
            for f in os.listdir(lig_dir):
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

        print("Generating AM1-BCC charges with GAFF2 atom types...")
        cmd = [
            "antechamber",
            "-i", os.path.basename(input_mol2),
            "-fi", "mol2",
            "-o", os.path.basename(amber_mol2),
            "-fo", "mol2",
            "-c", "bcc",
            "-s", "2",
            "-at", self._get_atom_type_flag()  # Use GAFF2 atom types
        ]
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])

        try:
            self.run_command(cmd, cwd=ff_dir)
        except Exception as e:
            print(f"\nAM1-BCC charge generation failed. Attempting fallback with gas phase charges...")
            # Try with simpler gas phase charges as fallback
            cmd_fallback = [
                "antechamber",
                "-i", os.path.basename(input_mol2),
                "-fi", "mol2",
                "-o", os.path.basename(amber_mol2),
                "-fo", "mol2",
                "-c", "gas",
                "-s", "2",
                "-at", self._get_atom_type_flag()  # Use GAFF2 atom types
            ]
            if self.net_charge != 0:
                cmd_fallback.extend(["-nc", str(self.net_charge)])
            self.run_command(cmd_fallback, cwd=ff_dir)

        print("Generating prep file with GAFF2...")
        # Generate prep file from the charged mol2 file
        cmd = [
            "antechamber",
            "-i", os.path.basename(amber_mol2),
            "-fi", "mol2",
            "-o", os.path.basename(prep_file),
            "-fo", "prepi",
            "-at", self._get_atom_type_flag()  # Use GAFF2 atom types
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
                "-i", os.path.basename(amber_mol2),
                "-fi", "mol2",
                "-o", os.path.basename(prep_file),
                "-fo", "prepi",
                "-c", "rc",  # Read charges from input file
                "-at", self._get_atom_type_flag()  # Use GAFF2 atom types
            ]
            if self.net_charge != 0:
                cmd_alt.extend(["-nc", str(self.net_charge)])
            self.run_command(cmd_alt, cwd=ff_dir)

        print("Generating GAFF2 force field parameters...")
        self.run_command([
            "parmchk2",
            "-i", os.path.basename(prep_file),
            "-f", "prepi",
            "-o", os.path.basename(frcmod_file),
            "-a", "Y",  # Print all force field parameters
            "-s", "2"   # GAFF2 parameter set
        ], cwd=ff_dir)

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

        with open(tleap_input, 'w') as f:
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
        print("Processing SDF format (direct method) with GAFF2 atom types...")

        ff_dir = os.path.dirname(amber_mol2)

        cmd = [
            "antechamber",
            "-i", self.ligand_path,
            "-fi", "sdf",
            "-o", amber_mol2,
            "-fo", "mol2",
            "-c", "bcc",
            "-s", "2",
            "-at", self._get_atom_type_flag()  # Use GAFF2 atom types
        ]
        if self.net_charge != 0:
            cmd.extend(["-nc", str(self.net_charge)])
        self.run_command(cmd, cwd=ff_dir)

        # Continue with prep file generation
        self._process_mol2_intermediate(amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file)

    def _process_mol2_intermediate(self, amber_mol2, prep_file, frcmod_file, prmtop_file, rst7_file, charge_method=None):
        """Process intermediate MOL2 file with GAFF2

        Note: When processing MOL2 from SDF (already has BCC charges), we don't specify
        a charge method (-c flag) and let antechamber handle it automatically. This matches
        the behavior of _process_mol2_format() which also omits the -c flag.
        """
        print("Generating prep file from charged MOL2 with GAFF2...")

        ff_dir = os.path.dirname(amber_mol2)

        # Verify input file exists
        if not os.path.exists(amber_mol2):
            raise FileNotFoundError(f"Input MOL2 file not found: {amber_mol2}")

        # Build command WITHOUT -c flag (let antechamber handle charges)
        cmd = [
            "antechamber",
            "-i", os.path.basename(amber_mol2),
            "-fi", "mol2",
            "-o", os.path.basename(prep_file),
            "-fo", "prepi",
            "-at", self._get_atom_type_flag()  # Use GAFF2 atom types
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
                "-i", os.path.basename(amber_mol2),
                "-fi", "mol2",
                "-o", os.path.basename(prep_file),
                "-fo", "prepi",
                "-c", "rc",  # Read charges from input file
                "-at", self._get_atom_type_flag()  # Use GAFF2 atom types
            ]
            if self.net_charge != 0:
                cmd_alt.extend(["-nc", str(self.net_charge)])
            self.run_command(cmd_alt, cwd=ff_dir)

        print("Generating GAFF2 force field parameters...")
        self.run_command([
            "parmchk2",
            "-i", os.path.basename(prep_file),
            "-f", "prepi",
            "-o", os.path.basename(frcmod_file),
            "-s", "2"  # GAFF2 parameter set
        ], cwd=ff_dir)

        print("Creating AMBER topology with GAFF2...")
        self._create_topology_with_tleap(amber_mol2, frcmod_file, prmtop_file, rst7_file)
