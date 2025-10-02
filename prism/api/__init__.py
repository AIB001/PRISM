#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM API Module

Provides programmatic access to PRISM PMF calculation functionality through
Python library APIs, internal module communication, and data exchange interfaces.
"""

from .core import PMFCalculator, WorkflowManager, DataAnalyzer
from .interfaces import PMFInterface, MDInterface, SystemInterface
from .data import FormatConverter, DataExporter, DataImporter
from .client import PrismAPIClient
from .exceptions import APIError, ValidationError, CalculationError

__version__ = "1.0.0"
__author__ = "PRISM Development Team"

__all__ = [
    # Core API classes
    'PMFCalculator',
    'WorkflowManager',
    'DataAnalyzer',

    # Interface classes
    'PMFInterface',
    'MDInterface',
    'SystemInterface',

    # Data handling
    'FormatConverter',
    'DataExporter',
    'DataImporter',

    # Client
    'PrismAPIClient',
    'create_client',

    # Exceptions
    'APIError',
    'ValidationError',
    'CalculationError',

    # Utility functions
    'get_api_info'
]

# API version and compatibility
API_VERSION = "v1"
MIN_PYTHON_VERSION = (3, 7)
SUPPORTED_FORMATS = ['gromacs', 'amber', 'charmm', 'openmm']


def create_client(base_url=None, api_key=None, timeout=30):
    """
    Create a PRISM API client instance.

    Parameters:
    -----------
    base_url : str, optional
        Base URL for API server (default: local instance)
    api_key : str, optional
        API authentication key
    timeout : int, optional
        Request timeout in seconds (default: 30)

    Returns:
    --------
    PrismAPIClient
        Configured API client instance

    Examples:
    ---------
    >>> import prism.api as api
    >>> client = api.create_client()
    >>> result = client.run_pmf_calculation(system_dir="./my_system")

    >>> # With remote server
    >>> client = api.create_client(
    ...     base_url="https://api.prism-pmf.org",
    ...     api_key="your_api_key"
    ... )
    """
    return PrismAPIClient(base_url=base_url, api_key=api_key, timeout=timeout)


def get_api_info():
    """
    Get API module information.

    Returns:
    --------
    dict
        API module information including version and capabilities
    """
    return {
        'version': __version__,
        'api_version': API_VERSION,
        'min_python_version': '.'.join(map(str, MIN_PYTHON_VERSION)),
        'supported_formats': SUPPORTED_FORMATS,
        'available_interfaces': ['PMF', 'MD', 'System'],
        'data_formats': ['JSON', 'YAML', 'HDF5', 'CSV']
    }