#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Utils - Utility modules for PRISM Builder
"""

from .environment import GromacsEnvironment
from .config import ConfigurationManager
from .mdp import MDPGenerator
from .system import SystemBuilder

# Error handling (optional import to avoid breaking existing code)
try:
    from .error_handler import (
        ErrorHandler, get_error_handler, handle_error,
        log_warning, log_info, log_debug,
        PRISMError, PRISMSystemError, PRISMConfigurationError,
        PRISMSimulationError, PRISMValidationError, ErrorSeverity
    )
    _has_error_handler = True
except ImportError:
    _has_error_handler = False

__all__ = [
    'GromacsEnvironment',
    'ConfigurationManager',
    'MDPGenerator',
    'SystemBuilder',
]

# Add error handler exports if available
if _has_error_handler:
    __all__.extend([
        'ErrorHandler', 'get_error_handler', 'handle_error',
        'log_warning', 'log_info', 'log_debug',
        'PRISMError', 'PRISMSystemError', 'PRISMConfigurationError',
        'PRISMSimulationError', 'PRISMValidationError', 'ErrorSeverity'
    ])