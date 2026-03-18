#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Base class for SystemBuilder with common utilities.
"""

import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple

# Import color utilities
try:
    from ..colors import print_success, print_warning, print_info, success, warning, path, number
except ImportError:
    try:
        from prism.utils.colors import print_success, print_warning, print_info, success, warning, path, number
    except ImportError:
        # Fallback
        def print_success(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}✓ {x}")

        def print_warning(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}⚠ {x}")

        def print_info(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}ℹ {x}")

        def success(x):
            return f"✓ {x}"

        def warning(x):
            return f"⚠ {x}"

        def path(x):
            return x

        def number(x):
            return x


class SystemBuilderBase:
    """Base class for GROMACS system building operations."""

    def __init__(
        self,
        config: Dict,
        output_dir: str,
        overwrite: bool = False,
        pmf_mode: bool = False,
        box_extension: Optional[Tuple[float, float, float]] = None,
    ):
        """
        Initializes the SystemBuilder base.

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
        cmd_str = " ".join(map(str, command))
        print(f"Executing in {work_dir}: {cmd_str}")

        process = subprocess.Popen(
            command, cwd=work_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        stdout, stderr = process.communicate(input=input_str)

        if process.returncode != 0:
            print("--- STDOUT ---")
            print(stdout)
            print("--- STDERR ---")
            print(stderr)
            raise RuntimeError(f"Command failed with exit code {process.returncode}: {cmd_str}")

        return stdout, stderr

    def _get_terminal_menu_index(self, terminal_type: str, terminus: str, ff_info: dict = None) -> str:
        """
        Get the menu index for N/C-terminal selection in pdb2gmx.

        Args:
            terminal_type: Either 'nter' or 'cter'
            terminus: The terminal type (e.g., 'NH3+', 'COO-', 'None')
            ff_info: Force field info dict containing 'name' key

        Returns:
            Menu index as string, or empty string for auto-selection
        """
        if not terminus or terminus.lower() == "none":
            return ""

        # Get force field name
        ff_name = ff_info.get("name", "").lower() if ff_info else ""

        # Terminal options mapping for different force fields
        terminal_options = {
            "nter": {
                "amber": {"NH3+": "0", "NH2": "2", "None": "1"},
                "charmm": {"NH3+": "0", "NH2": "1", "None": "2"},
                "default": {"NH3+": "0", "NH2": "1", "None": "2"},
            },
            "cter": {
                "amber": {"COO-": "0", "COOH": "1", "None": "2"},
                "charmm": {"COO-": "0", "COOH": "1", "None": "2"},
                "default": {"COO-": "0", "COOH": "1", "None": "2"},
            },
        }

        # Determine which mapping to use
        if "amber" in ff_name:
            mapping = terminal_options[terminal_type]["amber"]
        elif "charmm" in ff_name:
            mapping = terminal_options[terminal_type]["charmm"]
        else:
            mapping = terminal_options[terminal_type]["default"]

        return mapping.get(terminus, "")

    def _count_histidines(self, pdb_file: str) -> int:
        """Counts unique histidine residues that need interactive protonation selection.

        Only counts generic 'HIS' residues (3-char name with space at position 20).
        Residues already named HIE/HID/HIP (by PROPKA) or HISE/HISD/HISH (by
        _rename_his_for_cmap) don't trigger pdb2gmx prompts and are excluded.
        """
        his_residues = set()
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") and line[17:21] == "HIS ":
                    resnum = line[22:26].strip()
                    chain = line[21].strip()
                    res_id = f"{chain}:{resnum}" if chain else resnum
                    his_residues.add(res_id)
        return len(his_residues)
