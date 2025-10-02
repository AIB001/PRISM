#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Package Manager

Comprehensive package management and installation system for PRISM dependencies
and optional components.
"""

import sys
import subprocess
import importlib
import pkg_resources
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from enum import Enum
import json
import warnings

try:
    from .logging_system import PrismLogger
except ImportError:
    import logging
    PrismLogger = logging.getLogger


class DependencyType(Enum):
    """Types of dependencies"""
    REQUIRED = "required"
    OPTIONAL = "optional"
    DEVELOPMENT = "development"
    SYSTEM = "system"


class InstallationStatus(Enum):
    """Installation status"""
    INSTALLED = "installed"
    MISSING = "missing"
    OUTDATED = "outdated"
    INCOMPATIBLE = "incompatible"
    UNKNOWN = "unknown"


@dataclass
class DependencyInfo:
    """Information about a dependency"""
    name: str
    type: DependencyType
    version_required: Optional[str] = None
    version_installed: Optional[str] = None
    status: InstallationStatus = InstallationStatus.UNKNOWN
    install_command: Optional[str] = None
    description: Optional[str] = None
    import_name: Optional[str] = None
    check_function: Optional[str] = None
    alternatives: List[str] = None


class PrismPackageManager:
    """Package manager for PRISM dependencies and components"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger("package_manager")
        self.dependencies = self._define_dependencies()
        self._cache = {}
    
    def _define_dependencies(self) -> Dict[str, DependencyInfo]:
        """Define all PRISM dependencies"""
        deps = {}
        
        # Required Python packages
        deps["numpy"] = DependencyInfo(
            name="numpy",
            type=DependencyType.REQUIRED,
            version_required=">=1.19.0",
            install_command="pip install numpy",
            description="Numerical computing library",
            import_name="numpy"
        )
        
        deps["scipy"] = DependencyInfo(
            name="scipy", 
            type=DependencyType.REQUIRED,
            version_required=">=1.7.0",
            install_command="pip install scipy",
            description="Scientific computing library",
            import_name="scipy"
        )
        
        deps["matplotlib"] = DependencyInfo(
            name="matplotlib",
            type=DependencyType.REQUIRED,
            version_required=">=3.3.0",
            install_command="pip install matplotlib",
            description="Plotting and visualization library",
            import_name="matplotlib"
        )
        
        deps["pandas"] = DependencyInfo(
            name="pandas",
            type=DependencyType.REQUIRED,
            version_required=">=1.2.0",
            install_command="pip install pandas",
            description="Data analysis and manipulation library",
            import_name="pandas"
        )
        
        deps["pyyaml"] = DependencyInfo(
            name="PyYAML",
            type=DependencyType.REQUIRED,
            version_required=">=5.4.0",
            install_command="pip install PyYAML",
            description="YAML parser and emitter",
            import_name="yaml"
        )
        
        # Optional molecular dynamics packages
        deps["mdtraj"] = DependencyInfo(
            name="mdtraj",
            type=DependencyType.OPTIONAL,
            version_required=">=1.9.0",
            install_command="pip install mdtraj",
            description="Molecular dynamics trajectory analysis",
            import_name="mdtraj"
        )
        
        deps["openmm"] = DependencyInfo(
            name="openmm",
            type=DependencyType.OPTIONAL,
            version_required=">=7.5.0",
            install_command="conda install -c conda-forge openmm",
            description="High-performance molecular simulation toolkit",
            import_name="openmm"
        )
        
        deps["rdkit"] = DependencyInfo(
            name="rdkit-pypi",
            type=DependencyType.OPTIONAL,
            version_required=">=2021.03.0",
            install_command="pip install rdkit-pypi",
            description="Cheminformatics and machine learning toolkit",
            import_name="rdkit",
            alternatives=["conda install -c conda-forge rdkit"]
        )
        
        deps["openff"] = DependencyInfo(
            name="openff-toolkit",
            type=DependencyType.OPTIONAL,
            version_required=">=0.10.0",
            install_command="pip install openff-toolkit",
            description="Open Force Field toolkit for molecular simulations",
            import_name="openff.toolkit"
        )
        
        # PDB handling
        deps["pdbfixer"] = DependencyInfo(
            name="pdbfixer",
            type=DependencyType.OPTIONAL,
            version_required=">=1.7.0",
            install_command="conda install -c conda-forge pdbfixer",
            description="PDB structure preparation and fixing",
            import_name="pdbfixer"
        )
        
        # Web and networking
        deps["requests"] = DependencyInfo(
            name="requests",
            type=DependencyType.OPTIONAL,
            version_required=">=2.25.0",
            install_command="pip install requests",
            description="HTTP library for API access",
            import_name="requests"
        )
        
        # Development dependencies
        deps["pytest"] = DependencyInfo(
            name="pytest",
            type=DependencyType.DEVELOPMENT,
            version_required=">=6.0.0",
            install_command="pip install pytest",
            description="Testing framework",
            import_name="pytest"
        )
        
        deps["black"] = DependencyInfo(
            name="black",
            type=DependencyType.DEVELOPMENT,
            version_required=">=21.0.0",
            install_command="pip install black",
            description="Code formatter",
            import_name="black"
        )
        
        deps["flake8"] = DependencyInfo(
            name="flake8",
            type=DependencyType.DEVELOPMENT,
            version_required=">=3.9.0",
            install_command="pip install flake8",
            description="Code linter",
            import_name="flake8"
        )
        
        # System dependencies
        deps["gromacs"] = DependencyInfo(
            name="gromacs",
            type=DependencyType.SYSTEM,
            description="Molecular dynamics simulation package",
            check_function="gmx --version",
            install_command="# Install GROMACS from https://gromacs.org or package manager"
        )
        
        deps["amber"] = DependencyInfo(
            name="amber",
            type=DependencyType.SYSTEM,
            description="Molecular dynamics simulation suite",
            check_function="antechamber -h",
            install_command="# Install AmberTools from https://ambermd.org"
        )
        
        return deps
    
    def check_dependency(self, name: str, force_refresh: bool = False) -> DependencyInfo:
        """
        Check status of a specific dependency
        
        Args:
            name: Dependency name
            force_refresh: Force refresh of cached status
            
        Returns:
            Updated DependencyInfo object
        """
        if name not in self.dependencies:
            self.logger.warning(f"Unknown dependency: {name}")
            return DependencyInfo(name=name, type=DependencyType.OPTIONAL, 
                                status=InstallationStatus.UNKNOWN)
        
        if name in self._cache and not force_refresh:
            return self._cache[name]
        
        dep = self.dependencies[name]
        
        # Check system dependencies
        if dep.type == DependencyType.SYSTEM:
            dep.status = self._check_system_dependency(dep)
        else:
            # Check Python packages
            dep.status, dep.version_installed = self._check_python_package(dep)
        
        self._cache[name] = dep
        return dep
    
    def _check_python_package(self, dep: DependencyInfo) -> Tuple[InstallationStatus, Optional[str]]:
        """Check Python package installation status"""
        try:
            # Try to import the package
            if dep.import_name:
                module = importlib.import_module(dep.import_name)
                
                # Get version if available
                version = None
                for attr in ['__version__', 'version', 'VERSION']:
                    if hasattr(module, attr):
                        version = getattr(module, attr)
                        break
                
                # If no version found, try pkg_resources
                if not version:
                    try:
                        version = pkg_resources.get_distribution(dep.name).version
                    except:
                        pass
                
                # Check version compatibility if required
                if dep.version_required and version:
                    if self._check_version_compatibility(version, dep.version_required):
                        return InstallationStatus.INSTALLED, version
                    else:
                        return InstallationStatus.OUTDATED, version
                else:
                    return InstallationStatus.INSTALLED, version
            
            else:
                # Check using pkg_resources only
                pkg_resources.get_distribution(dep.name)
                return InstallationStatus.INSTALLED, None
                
        except ImportError:
            return InstallationStatus.MISSING, None
        except pkg_resources.DistributionNotFound:
            return InstallationStatus.MISSING, None
        except Exception as e:
            self.logger.debug(f"Error checking {dep.name}: {e}")
            return InstallationStatus.UNKNOWN, None
    
    def _check_system_dependency(self, dep: DependencyInfo) -> InstallationStatus:
        """Check system dependency availability"""
        if not dep.check_function:
            return InstallationStatus.UNKNOWN
        
        try:
            # Run the check command
            result = subprocess.run(
                dep.check_function.split(),
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                return InstallationStatus.INSTALLED
            else:
                return InstallationStatus.MISSING
                
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            return InstallationStatus.MISSING
        except Exception as e:
            self.logger.debug(f"Error checking system dependency {dep.name}: {e}")
            return InstallationStatus.UNKNOWN
    
    def _check_version_compatibility(self, installed: str, required: str) -> bool:
        """Check if installed version meets requirements"""
        try:
            from packaging import version
            
            # Parse requirement (e.g., ">=1.19.0")
            if required.startswith(">="):
                min_version = required[2:].strip()
                return version.parse(installed) >= version.parse(min_version)
            elif required.startswith(">"):
                min_version = required[1:].strip()
                return version.parse(installed) > version.parse(min_version)
            elif required.startswith("=="):
                exact_version = required[2:].strip()
                return version.parse(installed) == version.parse(exact_version)
            else:
                # Assume exact match
                return installed == required
                
        except ImportError:
            # Fallback to string comparison
            self.logger.warning("packaging library not available, using string comparison")
            return True  # Assume compatible
        except Exception as e:
            self.logger.debug(f"Version comparison error: {e}")
            return True  # Assume compatible
    
    def check_all_dependencies(self, include_types: List[DependencyType] = None) -> Dict[str, DependencyInfo]:
        """
        Check all dependencies of specified types
        
        Args:
            include_types: List of dependency types to check (default: all)
            
        Returns:
            Dictionary of dependency names to their info
        """
        if include_types is None:
            include_types = list(DependencyType)
        
        results = {}
        
        for name, dep in self.dependencies.items():
            if dep.type in include_types:
                results[name] = self.check_dependency(name)
        
        return results
    
    def get_missing_dependencies(self, include_types: List[DependencyType] = None) -> List[DependencyInfo]:
        """Get list of missing dependencies"""
        all_deps = self.check_all_dependencies(include_types)
        return [dep for dep in all_deps.values() 
                if dep.status in [InstallationStatus.MISSING, InstallationStatus.OUTDATED]]
    
    def get_installation_commands(self, missing_only: bool = True, 
                                include_types: List[DependencyType] = None) -> List[str]:
        """
        Get installation commands for dependencies
        
        Args:
            missing_only: Only return commands for missing/outdated dependencies
            include_types: Dependency types to include
            
        Returns:
            List of installation commands
        """
        commands = []
        
        if missing_only:
            deps = self.get_missing_dependencies(include_types)
        else:
            deps = list(self.check_all_dependencies(include_types).values())
        
        for dep in deps:
            if dep.install_command and not dep.install_command.startswith("#"):
                commands.append(dep.install_command)
            elif dep.alternatives:
                commands.extend([cmd for cmd in dep.alternatives if not cmd.startswith("#")])
        
        return commands
    
    def install_dependency(self, name: str, force: bool = False) -> bool:
        """
        Install a specific dependency
        
        Args:
            name: Dependency name
            force: Force reinstallation even if already installed
            
        Returns:
            True if installation successful
        """
        dep = self.check_dependency(name, force_refresh=True)
        
        if not force and dep.status == InstallationStatus.INSTALLED:
            self.logger.info(f"{name} is already installed")
            return True
        
        if not dep.install_command or dep.install_command.startswith("#"):
            self.logger.warning(f"No installation command available for {name}")
            return False
        
        self.logger.info(f"Installing {name}...")
        
        try:
            # Execute installation command
            result = subprocess.run(
                dep.install_command.split(),
                capture_output=True,
                text=True,
                timeout=300  # 5 minutes timeout
            )
            
            if result.returncode == 0:
                self.logger.info(f"Successfully installed {name}")
                # Clear cache to force recheck
                if name in self._cache:
                    del self._cache[name]
                return True
            else:
                self.logger.error(f"Failed to install {name}: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            self.logger.error(f"Installation of {name} timed out")
            return False
        except Exception as e:
            self.logger.error(f"Error installing {name}: {e}")
            return False
    
    def install_missing_required(self) -> Tuple[int, int]:
        """
        Install all missing required dependencies
        
        Returns:
            Tuple of (successful_installs, failed_installs)
        """
        missing = self.get_missing_dependencies([DependencyType.REQUIRED])
        successful = 0
        failed = 0
        
        for dep in missing:
            if self.install_dependency(dep.name):
                successful += 1
            else:
                failed += 1
        
        return successful, failed
    
    def generate_requirements_txt(self, output_file: str = "requirements.txt",
                                include_types: List[DependencyType] = None) -> str:
        """
        Generate requirements.txt file
        
        Args:
            output_file: Output file path
            include_types: Dependency types to include
            
        Returns:
            Content of requirements.txt file
        """
        if include_types is None:
            include_types = [DependencyType.REQUIRED, DependencyType.OPTIONAL]
        
        lines = ["# PRISM PMF System Requirements", "# Generated automatically", ""]
        
        for dep_type in include_types:
            type_deps = [dep for dep in self.dependencies.values() if dep.type == dep_type]
            if type_deps:
                lines.append(f"# {dep_type.value.title()} dependencies")
                
                for dep in type_deps:
                    if dep.name and not dep.type == DependencyType.SYSTEM:
                        if dep.version_required:
                            lines.append(f"{dep.name}{dep.version_required}")
                        else:
                            lines.append(dep.name)
                
                lines.append("")
        
        content = "\\n".join(lines)
        
        # Write to file
        try:
            Path(output_file).write_text(content)
            self.logger.info(f"Requirements file written to {output_file}")
        except Exception as e:
            self.logger.error(f"Failed to write requirements file: {e}")
        
        return content
    
    def generate_install_script(self, output_file: str = "install_prism.sh",
                              include_system: bool = True) -> str:
        """
        Generate installation script
        
        Args:
            output_file: Output script file path
            include_system: Include system dependencies
            
        Returns:
            Content of installation script
        """
        lines = [
            "#!/bin/bash",
            "# PRISM PMF System Installation Script",
            "# Generated automatically",
            "",
            "set -e  # Exit on error",
            "",
            "echo 'Installing PRISM PMF System dependencies...'",
            ""
        ]
        
        # Python dependencies
        python_commands = self.get_installation_commands(
            missing_only=False,
            include_types=[DependencyType.REQUIRED, DependencyType.OPTIONAL]
        )
        
        if python_commands:
            lines.extend([
                "echo 'Installing Python dependencies...'",
                ""
            ])
            
            for cmd in python_commands:
                lines.append(f"{cmd}")
            
            lines.append("")
        
        # System dependencies
        if include_system:
            system_deps = [dep for dep in self.dependencies.values() 
                          if dep.type == DependencyType.SYSTEM]
            
            if system_deps:
                lines.extend([
                    "echo 'System dependencies (install manually):'",
                    ""
                ])
                
                for dep in system_deps:
                    lines.append(f"# {dep.name}: {dep.description}")
                    lines.append(f"# {dep.install_command}")
                    lines.append("")
        
        lines.extend([
            "echo 'PRISM installation completed!'",
            "echo 'Run: python -c \"import prism; prism.check_dependencies()\" to verify'"
        ])
        
        content = "\\n".join(lines)
        
        # Write to file
        try:
            script_path = Path(output_file)
            script_path.write_text(content)
            script_path.chmod(0o755)  # Make executable
            self.logger.info(f"Installation script written to {output_file}")
        except Exception as e:
            self.logger.error(f"Failed to write installation script: {e}")
        
        return content
    
    def get_dependency_report(self) -> Dict[str, Any]:
        """Generate comprehensive dependency report"""
        all_deps = self.check_all_dependencies()
        
        report = {
            "timestamp": str(Path(__file__).stat().st_mtime),
            "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            "total_dependencies": len(all_deps),
            "by_type": {},
            "by_status": {},
            "dependencies": {}
        }
        
        # Count by type
        for dep_type in DependencyType:
            type_deps = [dep for dep in all_deps.values() if dep.type == dep_type]
            report["by_type"][dep_type.value] = {
                "total": len(type_deps),
                "installed": len([d for d in type_deps if d.status == InstallationStatus.INSTALLED]),
                "missing": len([d for d in type_deps if d.status == InstallationStatus.MISSING])
            }
        
        # Count by status
        for status in InstallationStatus:
            status_deps = [dep for dep in all_deps.values() if dep.status == status]
            report["by_status"][status.value] = len(status_deps)
        
        # Detailed dependency info
        for name, dep in all_deps.items():
            report["dependencies"][name] = {
                "type": dep.type.value,
                "status": dep.status.value,
                "version_required": dep.version_required,
                "version_installed": dep.version_installed,
                "description": dep.description
            }
        
        return report
    
    def print_dependency_status(self, include_types: List[DependencyType] = None):
        """Print formatted dependency status"""
        all_deps = self.check_all_dependencies(include_types)
        
        print("\\n🔍 PRISM Dependency Status")
        print("=" * 50)
        
        for dep_type in DependencyType:
            if include_types and dep_type not in include_types:
                continue
                
            type_deps = [dep for dep in all_deps.values() if dep.type == dep_type]
            if not type_deps:
                continue
            
            print(f"\\n📦 {dep_type.value.title()} Dependencies:")
            
            for dep in sorted(type_deps, key=lambda x: x.name):
                status_icon = {
                    InstallationStatus.INSTALLED: "✅",
                    InstallationStatus.MISSING: "❌", 
                    InstallationStatus.OUTDATED: "⚠️",
                    InstallationStatus.INCOMPATIBLE: "❌",
                    InstallationStatus.UNKNOWN: "❓"
                }.get(dep.status, "❓")
                
                version_info = ""
                if dep.version_installed:
                    version_info = f" (v{dep.version_installed})"
                elif dep.version_required:
                    version_info = f" (requires {dep.version_required})"
                
                print(f"  {status_icon} {dep.name}{version_info}")
                if dep.description:
                    print(f"      {dep.description}")
                
                if dep.status in [InstallationStatus.MISSING, InstallationStatus.OUTDATED]:
                    if dep.install_command and not dep.install_command.startswith("#"):
                        print(f"      Install: {dep.install_command}")


def main():
    """Main function for standalone usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description="PRISM Package Manager")
    parser.add_argument("--check", action="store_true", help="Check all dependencies")
    parser.add_argument("--install-missing", action="store_true", help="Install missing required dependencies")
    parser.add_argument("--generate-requirements", help="Generate requirements.txt file")
    parser.add_argument("--generate-script", help="Generate installation script")
    parser.add_argument("--report", help="Generate JSON dependency report")
    
    args = parser.parse_args()
    
    manager = PrismPackageManager()
    
    if args.check:
        manager.print_dependency_status()
    
    if args.install_missing:
        print("\\n🔧 Installing missing required dependencies...")
        successful, failed = manager.install_missing_required()
        print(f"Installation complete: {successful} successful, {failed} failed")
    
    if args.generate_requirements:
        content = manager.generate_requirements_txt(args.generate_requirements)
        print(f"Requirements file generated: {args.generate_requirements}")
    
    if args.generate_script:
        content = manager.generate_install_script(args.generate_script)
        print(f"Installation script generated: {args.generate_script}")
    
    if args.report:
        report = manager.get_dependency_report()
        with open(args.report, 'w') as f:
            json.dump(report, f, indent=2)
        print(f"Dependency report generated: {args.report}")


if __name__ == "__main__":
    main()