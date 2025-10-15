#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Error Handling Utilities

This module provides decorators, context managers, and utilities for
robust error handling and recovery in PMF calculations.
"""

import time
import logging
import functools
from typing import Callable, Optional, List, Type, Union, Any, Dict
from contextlib import contextmanager
from pathlib import Path

from .exceptions import (
    PMFError, PMFErrorCode, create_error_from_exception,
    ResourceError, ExternalToolError, PMFFileError
)

logger = logging.getLogger(__name__)


def with_retry(
    max_attempts: int = 3,
    backoff_factor: float = 2.0,
    retry_on: Union[Type[Exception], List[Type[Exception]]] = PMFError,
    retry_condition: Optional[Callable[[Exception], bool]] = None
):
    """
    Decorator for automatic retry of operations with exponential backoff.
    
    Parameters:
    -----------
    max_attempts : int
        Maximum number of retry attempts
    backoff_factor : float
        Exponential backoff factor (delay = backoff_factor ** attempt)
    retry_on : Exception type or list of types
        Exception types that should trigger retry
    retry_condition : callable, optional
        Custom function to determine if exception should trigger retry
        
    Example:
    --------
    @with_retry(max_attempts=3, retry_on=[PMFFileError, ExternalToolError])
    def copy_files(src, dst):
        # File operation that might fail
        pass
    """
    if not isinstance(retry_on, (list, tuple)):
        retry_on = [retry_on]
    
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_exception = None
            
            for attempt in range(max_attempts):
                try:
                    return func(*args, **kwargs)
                    
                except Exception as exc:
                    # Convert to PMFError if needed
                    if not isinstance(exc, PMFError):
                        exc = create_error_from_exception(exc)
                    
                    last_exception = exc
                    
                    # Check if we should retry
                    should_retry = False
                    
                    # Check exception type
                    if any(isinstance(exc, retry_type) for retry_type in retry_on):
                        should_retry = True
                    
                    # Check custom condition
                    if retry_condition and retry_condition(exc):
                        should_retry = True
                    
                    # Check if error is recoverable
                    if isinstance(exc, PMFError) and not exc.recoverable:
                        should_retry = False
                    
                    # Don't retry on last attempt
                    if attempt == max_attempts - 1:
                        should_retry = False
                    
                    if not should_retry:
                        raise exc
                    
                    # Calculate delay and wait
                    delay = backoff_factor ** attempt
                    logger.warning(
                        f"Attempt {attempt + 1}/{max_attempts} failed: {exc}. "
                        f"Retrying in {delay:.1f}s..."
                    )
                    time.sleep(delay)
            
            # This should never be reached, but just in case
            raise last_exception
            
        return wrapper
    return decorator


def with_timeout(timeout_seconds: int, operation_name: str = None):
    """
    Decorator to add timeout to operations.
    
    Parameters:
    -----------
    timeout_seconds : int
        Maximum time to allow for operation
    operation_name : str, optional
        Name of operation for error messages
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            import signal
            
            def timeout_handler(signum, frame):
                raise TimeoutError(
                    operation=operation_name or func.__name__,
                    timeout_seconds=timeout_seconds
                )
            
            # Set timeout
            old_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(timeout_seconds)
            
            try:
                result = func(*args, **kwargs)
                signal.alarm(0)  # Cancel timeout
                return result
            finally:
                signal.signal(signal.SIGALRM, old_handler)
                
        return wrapper
    return decorator


@contextmanager
def error_context(operation: str, context_data: Optional[Dict[str, Any]] = None):
    """
    Context manager to provide additional context for errors.
    
    Parameters:
    -----------
    operation : str
        Description of the operation being performed
    context_data : dict, optional
        Additional context data to include in errors
        
    Example:
    --------
    with error_context("SMD preparation", {"system_dir": "/path/to/system"}):
        prepare_smd_files()
    """
    try:
        yield
    except Exception as exc:
        # Enhance exception with context
        if isinstance(exc, PMFError):
            # Add context to existing PMFError
            if context_data:
                exc.context.update(context_data)
            exc.context['operation'] = operation
            raise exc
        else:
            # Convert to PMFError with context
            enhanced_context = {'operation': operation}
            if context_data:
                enhanced_context.update(context_data)
            
            pmf_error = create_error_from_exception(exc, enhanced_context)
            raise pmf_error from exc


@contextmanager
def safe_file_operation(description: str, cleanup_files: Optional[List[Path]] = None):
    """
    Context manager for safe file operations with cleanup on failure.
    
    Parameters:
    -----------
    description : str
        Description of the file operation
    cleanup_files : List[Path], optional
        Files to clean up if operation fails
        
    Example:
    --------
    with safe_file_operation("copying system files", cleanup_files=[temp_file]):
        shutil.copy2(source, destination)
    """
    try:
        yield
    except Exception as exc:
        # Clean up files on failure
        if cleanup_files:
            for file_path in cleanup_files:
                try:
                    if file_path.exists():
                        if file_path.is_file():
                            file_path.unlink()
                        elif file_path.is_dir():
                            import shutil
                            shutil.rmtree(file_path)
                        logger.debug(f"Cleaned up: {file_path}")
                except Exception as cleanup_exc:
                    logger.warning(f"Failed to clean up {file_path}: {cleanup_exc}")
        
        # Re-raise original exception
        if isinstance(exc, PMFError):
            raise exc
        else:
            raise PMFFileError(
                message=f"File operation failed: {description}",
                error_code=PMFErrorCode.FILE_CORRUPTED,
                recoverable=True,
                context={'description': description, 'cleanup_files': [str(f) for f in cleanup_files or []]},
                cause=exc
            ) from exc


class ErrorCollector:
    """
    Utility class to collect and manage multiple errors during batch operations.
    
    Example:
    --------
    collector = ErrorCollector()
    for window in windows:
        try:
            process_window(window)
        except Exception as e:
            collector.add_error(f"window_{window}", e)
    
    if collector.has_errors():
        collector.raise_summary()
    """
    
    def __init__(self):
        self.errors: Dict[str, Exception] = {}
        self.warnings: Dict[str, str] = {}
    
    def add_error(self, identifier: str, error: Exception):
        """Add an error for a specific item"""
        self.errors[identifier] = error
        logger.error(f"Error in {identifier}: {error}")
    
    def add_warning(self, identifier: str, message: str):
        """Add a warning for a specific item"""
        self.warnings[identifier] = message
        logger.warning(f"Warning in {identifier}: {message}")
    
    def has_errors(self) -> bool:
        """Check if any errors were collected"""
        return len(self.errors) > 0
    
    def has_warnings(self) -> bool:
        """Check if any warnings were collected"""
        return len(self.warnings) > 0
    
    def get_error_summary(self) -> str:
        """Get formatted summary of all errors"""
        lines = []
        
        if self.errors:
            lines.append(f"Errors ({len(self.errors)}):")
            for identifier, error in self.errors.items():
                lines.append(f"  • {identifier}: {error}")
        
        if self.warnings:
            lines.append(f"Warnings ({len(self.warnings)}):")
            for identifier, warning in self.warnings.items():
                lines.append(f"  • {identifier}: {warning}")
        
        return "\n".join(lines)
    
    def raise_summary(self):
        """Raise a summary error if any errors were collected"""
        if not self.has_errors():
            return
        
        error_count = len(self.errors)
        warning_count = len(self.warnings)
        
        message = f"Batch operation failed: {error_count} errors"
        if warning_count > 0:
            message += f", {warning_count} warnings"
        
        context = {
            'error_count': error_count,
            'warning_count': warning_count,
            'failed_items': list(self.errors.keys()),
            'warnings': dict(self.warnings)
        }
        
        # Determine recovery suggestions based on error types
        recovery_suggestions = ["Check individual item error details"]
        
        if any(isinstance(e, ExternalToolError) for e in self.errors.values()):
            recovery_suggestions.append("Check GROMACS installation and environment")
        
        if any(isinstance(e, PMFFileError) for e in self.errors.values()):
            recovery_suggestions.append("Check file permissions and disk space")
        
        raise PMFError(
            message=message,
            error_code=PMFErrorCode.WORKFLOW_STATE_INCONSISTENT,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


def validate_prerequisites(requirements: Dict[str, Callable[[], bool]], step_name: str):
    """
    Validate that all prerequisites are met before executing a step.
    
    Parameters:
    -----------
    requirements : Dict[str, Callable]
        Dictionary mapping requirement names to validation functions
    step_name : str
        Name of the step being validated
        
    Raises:
    -------
    PrerequisiteNotMetError
        If any requirements are not met
        
    Example:
    --------
    validate_prerequisites({
        'smd_completed': lambda: smd_dir.exists(),
        'files_present': lambda: all(f.exists() for f in required_files)
    }, 'umbrella_sampling')
    """
    from .exceptions import PrerequisiteNotMetError
    
    failed_requirements = []
    
    for req_name, check_func in requirements.items():
        try:
            if not check_func():
                failed_requirements.append(req_name)
        except Exception as exc:
            logger.error(f"Error checking requirement '{req_name}': {exc}")
            failed_requirements.append(req_name)
    
    if failed_requirements:
        raise PrerequisiteNotMetError(
            current_step=step_name,
            required_step=", ".join(failed_requirements),
            required_files=[]
        )


def handle_external_tool_error(command: str, exit_code: int, stderr: str, tool_name: str) -> None:
    """
    Handle errors from external tool execution (like GROMACS).
    
    Parameters:
    -----------
    command : str
        The command that was executed
    exit_code : int
        Exit code from the command
    stderr : str
        Standard error output
    tool_name : str
        Name of the external tool
        
    Raises:
    -------
    GromacsExecutionError or ExternalToolError
        Appropriate error based on the tool and error
    """
    from .exceptions import GromacsExecutionError
    
    # Parse common error patterns and provide specific suggestions
    recovery_suggestions = []
    
    if "command not found" in stderr.lower():
        recovery_suggestions.extend([
            f"Install {tool_name} and ensure it's in PATH",
            f"Source {tool_name} environment setup script"
        ])
    elif "permission denied" in stderr.lower():
        recovery_suggestions.extend([
            "Check file and directory permissions",
            "Ensure write access to output directories"
        ])
    elif "no space left" in stderr.lower():
        recovery_suggestions.extend([
            "Free up disk space",
            "Clean up temporary files"
        ])
    elif "memory" in stderr.lower() or "out of memory" in stderr.lower():
        recovery_suggestions.extend([
            "Reduce system size or increase available memory",
            "Use more conservative simulation parameters"
        ])
    else:
        recovery_suggestions.extend([
            f"Check {tool_name} documentation for error details",
            "Verify input file formats and parameters"
        ])
    
    if tool_name.lower() == "gromacs" or "gmx" in command:
        raise GromacsExecutionError(
            tool=tool_name,
            command=command,
            exit_code=exit_code,
            stderr_output=stderr
        )
    else:
        raise ExternalToolError(
            message=f"{tool_name} execution failed",
            error_code=PMFErrorCode.GROMACS_EXECUTION_FAILED,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context={
                'tool': tool_name,
                'command': command,
                'exit_code': exit_code,
                'stderr': stderr
            }
        )