#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM exception hierarchy.

This module defines a hierarchical exception system for PRISM that provides better
error messages and easier error handling.
"""


class PRISMError(Exception):
    """Base exception class for all PRISM errors.

    All exceptions in PRISM inherit from this class, allowing users to catch
    any PRISM-specific error with a single except clause.
    """

    pass


class ConfigurationError(PRISMError):
    """Configuration error.

    Raised when:
    - Configuration file is missing or invalid
    - Required configuration options are missing
    - Configuration values are out of range or invalid
    """

    pass


class ForceFieldError(PRISMError):
    """Force field generation or processing error.

    Raised when:
    - Force field generator fails
    - Required dependencies are missing (e.g. AmberTools, openff-toolkit)
    - Force field parameters are invalid or incomplete
    - Output files are missing after generation
    """

    pass


class DependencyError(ForceFieldError):
    """Dependency missing error.

    Raised when an external dependency required for force field generation
    is not installed or not found in PATH.
    """

    pass


class ParameterError(ForceFieldError):
    """Force field parameter error.

    Raised when:
    - Parameter file is missing or invalid
    - Required parameters are missing
    - Parameter values are invalid
    """

    pass


class TopologyError(PRISMError):
    """Topology processing error.

    Raised when:
    - Topology file is invalid
    - Required sections are missing
    - Parameter inclusion fails
    - Atom mapping doesn't match expected
    """

    pass


class SystemBuildError(PRISMError):
    """System building error.

    Raised when:
    - GROMACS command fails (grompp, genion, etc.)
    - Solvation fails
    - Ion addition fails
    - Output files are missing
    """

    pass


class FEPError(PRISMError):
    """Free Energy Perturbation error.

    Raised when:
    - Hybrid topology building fails
    - Mapping between ligands fails
    - Window setup fails
    - Lambda schedule configuration is invalid
    """

    pass


class PMFError(PRISMError):
    """Potential of Mean Force error.

    Raised when:
    - System setup for PMF fails
    - Umbrella window setup fails
    - Pulling coordinate setup fails
    """

    pass


class FileError(PRISMError):
    """File system error.

    Raised when:
    - Required file is missing
    - File has invalid format
    - Permission error when reading/writing
    """

    pass


class GROMACSError(PRISMError):
    """GROMACS execution error.

    Raised when a GROMACS command returns non-zero exit code.
    Stores the command that failed and the exit code.
    """

    def __init__(self, command: str, exit_code: int, stderr: str = ""):
        self.command = command
        self.exit_code = exit_code
        self.stderr = stderr
        message = f"GROMACS command failed with exit code {exit_code}:\n" f"Command: {command}\n" f"{stderr}"
        super().__init__(message)
