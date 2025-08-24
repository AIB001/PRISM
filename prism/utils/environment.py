#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Environment - GROMACS environment detection and configuration
"""

import os
import subprocess
import platform
from pathlib import Path


class GromacsEnvironment:
    """Handle GROMACS environment detection and configuration"""

    def __init__(self):
        self.gmx_command = None
        self.gmx_bin_path = None
        self.gmx_lib = None
        self.force_fields = {}
        self.force_field_names = {}  # Map names to indices
        self.water_models = {}
        self._detect_gromacs()

    def _detect_gromacs(self):
        """Detect GROMACS installation and environment"""
        # Try to find gmx command
        gmx_commands = ['gmx', 'gmx_mpi', 'gmx_d', 'gmx_mpi_d']

        for cmd in gmx_commands:
            gmx_path = self._which(cmd)
            if gmx_path:
                self.gmx_command = cmd
                self.gmx_bin_path = gmx_path
                break

        if not self.gmx_command:
            raise RuntimeError("GROMACS not found. Please install GROMACS or add it to PATH.")

        print(f"Found GROMACS command: {self.gmx_command} at {self.gmx_bin_path}")

        # Get GMXLIB - check environment variable first
        self.gmx_lib = os.environ.get('GMXLIB')

        if not self.gmx_lib:
            # Deduce from binary path
            # e.g., /opt/gromacs/2024.2/bin/gmx -> /opt/gromacs/2024.2/share/gromacs/top
            bin_dir = os.path.dirname(self.gmx_bin_path)
            gromacs_root = os.path.dirname(bin_dir)  # Go up from bin to root

            # Try common relative paths from root
            potential_paths = [
                os.path.join(gromacs_root, 'share', 'gromacs', 'top'),
                os.path.join(gromacs_root, 'share', 'top'),
                os.path.join(gromacs_root, 'top')
            ]

            for path in potential_paths:
                if os.path.exists(path) and os.path.isdir(path):
                    self.gmx_lib = path
                    break

        if not self.gmx_lib:
            # Try to detect from gmx version output as fallback
            try:
                result = subprocess.run([self.gmx_command, '-version'],
                                        capture_output=True, text=True)
                for line in result.stdout.split('\n'):
                    if 'Data prefix' in line:
                        data_prefix = line.split(':')[1].strip()
                        potential_lib = os.path.join(data_prefix, 'share', 'gromacs', 'top')
                        if os.path.exists(potential_lib):
                            self.gmx_lib = potential_lib
                            break
            except:
                pass

        if not self.gmx_lib:
            print("Warning: Could not detect GMXLIB. Using default force field list.")
        else:
            print(f"GROMACS library path: {self.gmx_lib}")
            self._scan_force_fields()

    def _which(self, command):
        """Find the full path of a command"""
        if platform.system() == 'Windows':
            command = command + '.exe'

        # Check if command contains path separator
        if os.path.sep in command:
            if os.path.isfile(command) and os.access(command, os.X_OK):
                return os.path.abspath(command)
            return None

        # Search in PATH
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, command)
            if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                return os.path.abspath(exe_file)

        return None

    def _scan_force_fields(self):
        """Scan available force fields from GROMACS installation"""
        if not self.gmx_lib or not os.path.exists(self.gmx_lib):
            return

        # Reset force fields
        self.force_fields = {}
        self.force_field_names = {}
        ff_index = 1

        # Scan for force field directories
        for item in sorted(os.listdir(self.gmx_lib)):
            item_path = os.path.join(self.gmx_lib, item)
            if os.path.isdir(item_path) and item.endswith('.ff'):
                ff_name = item[:-3]  # Remove .ff extension
                self.force_fields[ff_index] = {
                    "name": ff_name,
                    "dir": item
                }
                self.force_field_names[ff_name.lower()] = ff_index

                # Check for water models in this force field
                water_model_file = os.path.join(item_path, 'watermodels.dat')
                if os.path.exists(water_model_file):
                    self._parse_water_models(water_model_file, ff_index)

                ff_index += 1

        print(f"Found {len(self.force_fields)} force fields")

    def _parse_water_models(self, water_model_file, ff_index):
        """Parse water models from watermodels.dat file"""
        if ff_index not in self.water_models:
            self.water_models[ff_index] = {}

        wm_index = 1
        with open(water_model_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split(None, 2)  # Split into max 3 parts
                    if len(parts) >= 2:
                        model_name = parts[0]
                        model_label = parts[1]
                        model_desc = parts[2] if len(parts) > 2 else ""
                        self.water_models[ff_index][wm_index] = {
                            "name": model_name,
                            "label": model_label,
                            "description": model_desc
                        }
                        wm_index += 1

    def get_force_field_index(self, ff_name):
        """Get force field index from name"""
        ff_name_lower = ff_name.lower()

        # Direct match
        if ff_name_lower in self.force_field_names:
            return self.force_field_names[ff_name_lower]

        # Try with common variations
        variations = [
            ff_name_lower,
            ff_name_lower.replace('-', ''),
            ff_name_lower.replace('_', ''),
            ff_name_lower.replace('-', '_')
        ]

        for var in variations:
            if var in self.force_field_names:
                return self.force_field_names[var]

        # Partial match
        for stored_name, idx in self.force_field_names.items():
            if ff_name_lower in stored_name or stored_name in ff_name_lower:
                return idx

        return None

    def get_water_model_index(self, ff_index, water_name):
        """Get water model index from name for a specific force field"""
        if ff_index not in self.water_models:
            return None

        water_name_lower = water_name.lower()

        for wm_idx, wm_info in self.water_models[ff_index].items():
            if wm_info['name'].lower() == water_name_lower:
                return wm_idx

        return None

    def list_force_fields(self):
        """List all available force fields"""
        if not self.force_fields:
            print("No force fields detected. Using defaults.")
            return []

        ff_list = []
        for idx, ff_info in sorted(self.force_fields.items()):
            ff_list.append(f"{ff_info['name']} (index: {idx})")
        return ff_list

    def list_water_models(self, ff_index):
        """List water models for a specific force field"""
        if ff_index not in self.water_models:
            return ["tip3p", "tip4p", "spc", "spce", "none"]

        wm_list = []
        for idx, wm_info in sorted(self.water_models[ff_index].items()):
            desc = f" - {wm_info['description']}" if wm_info.get('description') else ""
            wm_list.append(f"{wm_info['name']}{desc}")
        return wm_list

    def get_force_field_config(self):
        """Get force field configuration for PRISM config"""
        if self.force_fields:
            return self.force_fields
        else:
            # Return default if detection failed
            return {
                1: {"name": "amber99sb", "dir": "amber99sb.ff"},
                2: {"name": "amber99sb-ildn", "dir": "amber99sb-ildn.ff"},
                3: {"name": "amber03", "dir": "amber03.ff"},
                4: {"name": "amber14sb", "dir": "amber14sb.ff"},
                5: {"name": "charmm27", "dir": "charmm27.ff"},
                6: {"name": "oplsaa", "dir": "oplsaa.ff"},
                7: {"name": "gromos54a7", "dir": "gromos54a7.ff"}
            }

    def get_water_models_for_forcefield(self, ff_index):
        """Get water models for a specific force field"""
        if ff_index in self.water_models:
            return self.water_models[ff_index]
        else:
            # Return default water models
            return {
                1: {"name": "tip3p", "description": "TIP 3-point, recommended"},
                2: {"name": "tip4p", "description": "TIP 4-point"},
                3: {"name": "spc", "description": "simple point charge"},
                4: {"name": "spce", "description": "extended simple point charge"},
                5: {"name": "none", "description": "no water model"}
            }