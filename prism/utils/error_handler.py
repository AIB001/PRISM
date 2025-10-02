#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Error Handler - Centralized error handling and logging

This module provides centralized error handling, logging, and exception management
for the PRISM system.
"""

import sys
import traceback
import logging
from typing import Optional, Any, Dict
from pathlib import Path


class PRISMError(Exception):
    """Base exception for PRISM-related errors"""
    pass


class PRISMSystemError(PRISMError):
    """System-level errors (file I/O, environment)"""
    pass


class PRISMConfigurationError(PRISMError):
    """Configuration-related errors"""
    pass


class PRISMSimulationError(PRISMError):
    """Simulation execution errors"""
    pass


class PRISMValidationError(PRISMError):
    """Input validation errors"""
    pass


class ErrorSeverity:
    """Error severity levels"""
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


class ErrorHandler:
    """Centralized error handling and logging"""

    def __init__(self, logger_name: str = "prism"):
        self.logger = logging.getLogger(logger_name)
        self._setup_logging()

    def _setup_logging(self):
        """Setup basic logging configuration"""
        if not self.logger.handlers:
            handler = logging.StreamHandler(sys.stdout)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            self.logger.setLevel(logging.INFO)

    def handle_error(self, error: Exception, context: str = "",
                    reraise: bool = True) -> Optional[Dict[str, Any]]:
        """Handle and log errors with context"""
        error_info = {
            'type': type(error).__name__,
            'message': str(error),
            'context': context,
            'traceback': traceback.format_exc()
        }

        self.logger.error(f"Error in {context}: {error}")
        if self.logger.level <= logging.DEBUG:
            self.logger.debug(f"Traceback: {error_info['traceback']}")

        if reraise:
            raise error
        else:
            return error_info

    def log_warning(self, message: str, context: str = ""):
        """Log warning message"""
        full_message = f"{context}: {message}" if context else message
        self.logger.warning(full_message)

    def log_info(self, message: str, context: str = ""):
        """Log info message"""
        full_message = f"{context}: {message}" if context else message
        self.logger.info(full_message)

    def log_debug(self, message: str, context: str = ""):
        """Log debug message"""
        full_message = f"{context}: {message}" if context else message
        self.logger.debug(full_message)


# Global error handler instance
_global_error_handler = None


def get_error_handler() -> ErrorHandler:
    """Get global error handler instance"""
    global _global_error_handler
    if _global_error_handler is None:
        _global_error_handler = ErrorHandler()
    return _global_error_handler


def handle_error(error: Exception, context: str = "", reraise: bool = True):
    """Convenience function for error handling"""
    return get_error_handler().handle_error(error, context, reraise)


def log_warning(message: str, context: str = ""):
    """Convenience function for warning logging"""
    get_error_handler().log_warning(message, context)


def log_info(message: str, context: str = ""):
    """Convenience function for info logging"""
    get_error_handler().log_info(message, context)


def log_debug(message: str, context: str = ""):
    """Convenience function for debug logging"""
    get_error_handler().log_debug(message, context)


def validate_file_exists(file_path: str, description: str = "File") -> Path:
    """Validate that a file exists"""
    path = Path(file_path)
    if not path.exists():
        raise PRISMValidationError(f"{description} not found: {file_path}")
    if not path.is_file():
        raise PRISMValidationError(f"{description} is not a file: {file_path}")
    return path


def validate_directory_exists(dir_path: str, description: str = "Directory") -> Path:
    """Validate that a directory exists"""
    path = Path(dir_path)
    if not path.exists():
        raise PRISMValidationError(f"{description} not found: {dir_path}")
    if not path.is_dir():
        raise PRISMValidationError(f"{description} is not a directory: {dir_path}")
    return path


def safe_import(module_name: str, description: str = "") -> Optional[Any]:
    """Safely import a module with error handling"""
    try:
        module = __import__(module_name, fromlist=[''])
        return module
    except ImportError as e:
        context = f"importing {description}" if description else f"importing {module_name}"
        log_warning(f"Failed to import {module_name}: {e}", context)
        return None


def create_error_context(operation: str, **kwargs) -> str:
    """Create error context string"""
    context_parts = [operation]
    for key, value in kwargs.items():
        context_parts.append(f"{key}={value}")
    return " ".join(context_parts)


# Exception context manager
class ErrorContext:
    """Context manager for error handling"""

    def __init__(self, operation: str, **kwargs):
        self.context = create_error_context(operation, **kwargs)
        self.error_handler = get_error_handler()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            self.error_handler.handle_error(exc_val, self.context, reraise=False)
            return True  # Suppress the exception
        return False


# Decorator for function error handling
def handle_errors(operation: str = "", reraise: bool = True):
    """Decorator for automatic error handling"""
    def decorator(func):
        def wrapper(*args, **kwargs):
            context = operation or f"function {func.__name__}"
            try:
                return func(*args, **kwargs)
            except Exception as e:
                return handle_error(e, context, reraise)
        return wrapper
    return decorator


__all__ = [
    'PRISMError',
    'PRISMSystemError',
    'PRISMConfigurationError',
    'PRISMSimulationError',
    'PRISMValidationError',
    'ErrorSeverity',
    'ErrorHandler',
    'get_error_handler',
    'handle_error',
    'log_warning',
    'log_info',
    'log_debug',
    'validate_file_exists',
    'validate_directory_exists',
    'safe_import',
    'create_error_context',
    'ErrorContext',
    'handle_errors'
]