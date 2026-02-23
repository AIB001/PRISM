#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Force Field Installer - Install additional force fields to GROMACS
"""

import os
import subprocess
from pathlib import Path


class ForceFieldInstaller:
    """Install additional force fields to GROMACS installation"""

    def __init__(self, gromacs_top_dir=None, prism_ff_dir=None):
        """
        Initialize force field installer

        Args:
            gromacs_top_dir: Path to GROMACS top directory (auto-detected if None)
            prism_ff_dir: Path to PRISM force field directory (auto-detected if None)
        """
        # Auto-detect GROMACS top directory
        if gromacs_top_dir is None:
            gromacs_top_dir = self._detect_gromacs_top_dir()
        self.gromacs_top_dir = Path(gromacs_top_dir)

        # Auto-detect PRISM force field directory
        if prism_ff_dir is None:
            prism_ff_dir = self._detect_prism_ff_dir()
        self.prism_ff_dir = Path(prism_ff_dir)

        # Scan available force fields
        self.available_ffs = self._scan_available_forcefields()
        self.installed_ffs = self._scan_installed_forcefields()

    def _detect_gromacs_top_dir(self):
        """Detect GROMACS top directory"""
        # Try to get from gmx command
        try:
            result = subprocess.run(["gmx", "--version"], capture_output=True, text=True, timeout=5)
            # Parse installation directory from version output
            for line in result.stdout.split("\n"):
                if "Data prefix:" in line:
                    prefix = line.split(":")[1].strip()
                    top_dir = os.path.join(prefix, "share", "gromacs", "top")
                    if os.path.exists(top_dir):
                        return top_dir
        except:
            pass

        # Common installation paths
        common_paths = [
            "/usr/local/gromacs/share/gromacs/top",
            "/usr/local/gromacs-2025.1/share/gromacs/top",
            "/usr/local/gromacs-2024/share/gromacs/top",
            "/opt/gromacs/share/gromacs/top",
        ]

        for path in common_paths:
            if os.path.exists(path):
                return path

        raise RuntimeError("Could not detect GROMACS installation directory")

    def _detect_prism_ff_dir(self):
        """Detect PRISM force field directory"""
        # Get PRISM package location
        try:
            import prism

            prism_dir = Path(prism.__file__).parent
            ff_dir = prism_dir / "configs" / "forcefield"
            if ff_dir.exists():
                return str(ff_dir)
        except:
            pass

        # Try relative to current file
        current_file = Path(__file__).resolve()
        ff_dir = current_file.parent.parent / "configs" / "forcefield"
        if ff_dir.exists():
            return str(ff_dir)

        raise RuntimeError("Could not find PRISM force field directory")

    def _scan_available_forcefields(self):
        """Scan available force fields in PRISM configs"""
        if not self.prism_ff_dir.exists():
            return {}

        forcefields = {}
        for item in sorted(self.prism_ff_dir.iterdir()):
            if item.is_dir() and item.name.endswith(".ff"):
                # Get description from forcefield.doc if exists
                desc_file = item / "forcefield.doc"
                if desc_file.exists():
                    try:
                        with open(desc_file, "r") as f:
                            description = f.readline().strip()
                    except:
                        description = self._get_ff_description(item.name)
                else:
                    description = self._get_ff_description(item.name)

                forcefields[item.name] = {"name": item.name, "path": item, "description": description}

        return forcefields

    def _get_ff_description(self, ff_name):
        """Get human-readable description for force field"""
        descriptions = {
            "a99SBdisp.ff": "Amber99SB-disp (improved dispersion interactions)",
            "amber14sb.ff": "Amber14SB (recommended for proteins)",
            "amber14sb_OL15.ff": "Amber14SB + OL15 (proteins + nucleic acids)",
            "charmm36-jul2022.ff": "CHARMM36 (July 2022 update)",
        }
        return descriptions.get(ff_name, ff_name.replace(".ff", ""))

    def _scan_installed_forcefields(self):
        """Scan installed force fields in GROMACS"""
        if not self.gromacs_top_dir.exists():
            return set()

        installed = set()
        for item in self.gromacs_top_dir.iterdir():
            if item.is_dir() and item.name.endswith(".ff"):
                installed.add(item.name)

        return installed

    def list_available_forcefields(self):
        """List available force fields for installation"""
        if not self.available_ffs:
            print("No additional force fields found in PRISM configs.")
            return

        print("\n" + "=" * 80)
        print("PRISM Additional Force Fields")
        print("=" * 80)
        print(f"\nGROMACS top directory: {self.gromacs_top_dir}")
        print(f"PRISM force field directory: {self.prism_ff_dir}\n")

        # List available force fields
        print("Available force fields for installation:\n")
        for idx, (ff_name, ff_info) in enumerate(self.available_ffs.items(), 1):
            status = "✓ INSTALLED" if ff_name in self.installed_ffs else "  NOT INSTALLED"
            print(f"  {idx}. {ff_name:<25} - {ff_info['description']}")
            print(f"     Status: {status}")

        print("\n" + "=" * 80)

    def install_forcefield(self, ff_name):
        """Install a single force field"""
        if ff_name not in self.available_ffs:
            print(f"Error: Force field '{ff_name}' not found in PRISM configs")
            return False

        if ff_name in self.installed_ffs:
            print(f"Force field '{ff_name}' is already installed")
            return True

        source = self.available_ffs[ff_name]["path"]
        target = self.gromacs_top_dir / ff_name

        print(f"\nInstalling {ff_name}...")
        print(f"  Source: {source}")
        print(f"  Target: {target}")

        # Check if we need sudo
        if not os.access(self.gromacs_top_dir, os.W_OK):
            print("\n  GROMACS directory requires sudo privileges.")
            print("  Please enter your password when prompted.\n")
            cmd = ["sudo", "cp", "-r", str(source), str(target)]
        else:
            cmd = ["cp", "-r", str(source), str(target)]

        try:
            result = subprocess.run(cmd, check=True)
            print(f"✓ Successfully installed {ff_name}")
            self.installed_ffs.add(ff_name)
            return True
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to install {ff_name}: {e}")
            return False
        except Exception as e:
            print(f"✗ Error installing {ff_name}: {e}")
            return False

    def install_forcefields_interactive(self):
        """Interactive installation of force fields"""
        self.list_available_forcefields()

        if not self.available_ffs:
            return

        # Check if all are already installed
        to_install = [ff for ff in self.available_ffs if ff not in self.installed_ffs]
        if not to_install:
            print("\n✓ All available force fields are already installed!")
            return

        print("\nOptions:")
        print("  - Enter numbers (e.g., 1,2,3) to install specific force fields")
        print("  - Enter 'all' to install all missing force fields")
        print("  - Press Enter to cancel")

        choice = input("\nYour choice: ").strip()

        if not choice:
            print("Installation cancelled.")
            return

        # Parse choice
        if choice.lower() == "all":
            ffs_to_install = list(self.available_ffs.keys())
        else:
            try:
                indices = [int(x.strip()) for x in choice.split(",")]
                ff_list = list(self.available_ffs.keys())
                ffs_to_install = [ff_list[i - 1] for i in indices if 1 <= i <= len(ff_list)]
            except (ValueError, IndexError) as e:
                print(f"Invalid input: {e}")
                return

        if not ffs_to_install:
            print("No force fields selected.")
            return

        # Install selected force fields
        print(f"\nInstalling {len(ffs_to_install)} force field(s)...")
        print("=" * 80)

        success_count = 0
        for ff_name in ffs_to_install:
            if self.install_forcefield(ff_name):
                success_count += 1

        print("\n" + "=" * 80)
        print(f"Installation complete: {success_count}/{len(ffs_to_install)} successful")

        if success_count > 0:
            print("\n✓ The following force fields are now available:")
            for ff in ffs_to_install:
                if ff in self.installed_ffs:
                    print(f"  - {ff}")

    def install_all_forcefields(self, silent=False):
        """Install all available force fields (for --install-forcefields --all)"""
        if not self.available_ffs:
            if not silent:
                print("No additional force fields found.")
            return

        to_install = [ff for ff in self.available_ffs if ff not in self.installed_ffs]

        if not to_install:
            if not silent:
                print("All force fields are already installed.")
            return

        if not silent:
            print(f"\nInstalling {len(to_install)} force field(s)...")
            print("=" * 80)

        success_count = 0
        for ff_name in to_install:
            if self.install_forcefield(ff_name):
                success_count += 1

        if not silent:
            print("\n" + "=" * 80)
            print(f"Installation complete: {success_count}/{len(to_install)} successful")


def main():
    """Main function for standalone usage"""
    installer = ForceFieldInstaller()
    installer.install_forcefields_interactive()


if __name__ == "__main__":
    main()
