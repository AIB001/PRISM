#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Solvation and ion addition utilities for SystemBuilder.
"""


class SolvationProcessorMixin:
    """Mixin for solvation and ion addition operations."""

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
            "tip3p": "spc216.gro",  # TIP3P uses SPC geometry
            "tip3p_original": "spc216.gro",
            "spc": "spc216.gro",
            "spce": "spc216.gro",
            "opc3": "spc216.gro",  # OPC3 is 3-point, uses SPC geometry
            "tip4p": "tip4p.gro",  # TIP4P has its own file with virtual sites
            "tip4pew": "tip4p.gro",  # TIP4P/Ew uses same geometry
            "opc": "tip4p.gro",  # OPC is 4-point with virtual site
            "tip5p": "tip5p.gro",  # TIP5P has its own file
        }

        water_model = getattr(self, "water_model_name", "tip3p").lower()
        water_file = water_coord_files.get(water_model, "spc216.gro")

        print(f"Using water model: {water_model} (coordinate file: {water_file})")

        command = [
            self.gmx_command,
            "solvate",
            "-cp",
            boxed_gro,
            "-cs",
            water_file,
            "-o",
            str(solvated_gro),
            "-p",
            topol_top,
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

        ions_cfg = self.config["ions"]
        temp_tpr = self.model_dir / "ions.tpr"

        # A minimal mdp is needed for grompp
        ions_mdp_path = self.model_dir / "ions.mdp"
        with open(ions_mdp_path, "w") as f:
            f.write("integrator=steep\nnsteps=0")

        grompp_cmd = [
            self.gmx_command,
            "grompp",
            "-f",
            str(ions_mdp_path),
            "-c",
            solvated_gro,
            "-p",
            topol_top,
            "-o",
            str(temp_tpr),
            "-maxwarn",
            "5",  # Allow some warnings
        ]

        # Run grompp with improper dihedral error recovery.
        # Some AMBER force fields (e.g. amber14sb) are missing improper dihedral
        # type definitions that pdb2gmx generates. Rather than guessing which
        # impropers to remove, we let GROMACS tell us exactly which lines fail.
        max_retries = 3
        grompp_success = False
        for attempt in range(max_retries + 1):
            cmd_str = " ".join(map(str, grompp_cmd))
            print(f"Executing in {self.model_dir}: {cmd_str}")

            result = subprocess.Popen(
                grompp_cmd,
                cwd=str(self.model_dir),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            stdout, stderr = result.communicate()

            if result.returncode == 0:
                grompp_success = True
                break

            # Check for fixable improper dihedral errors
            if "No default Per. Imp. Dih. types" in stderr and attempt < max_retries:
                fixed = self._fix_improper_from_grompp_errors(stderr)
                if fixed > 0:
                    print(
                        f"\n  Fixed {fixed} improper dihedral(s) without FF parameters, retrying grompp (attempt {attempt + 2}/{max_retries + 1})..."
                    )
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
            self.gmx_command,
            "genion",
            "-s",
            str(temp_tpr),
            "-o",
            str(ions_gro),
            "-p",
            topol_top,
            "-pname",
            ions_cfg["positive_ion"],
            "-nname",
            ions_cfg["negative_ion"],
        ]
        if ions_cfg.get("neutral", True):
            genion_cmd.append("-neutral")
        if ions_cfg.get("concentration", 0) > 0:
            genion_cmd.extend(["-conc", str(ions_cfg["concentration"])])

        # Provide "SOL" (the typical name for water) as input to stdin
        stdout, _ = self._run_command(genion_cmd, str(self.model_dir), input_str="SOL")

        # --- CRITICAL STEP: Parse genion output and fix topology ---
        self._update_topology_molecules(topol_top, stdout)

        # Clean up temporary files
        temp_tpr.unlink()
        ions_mdp_path.unlink()
