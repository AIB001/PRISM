#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Exception Hierarchy

This module defines a comprehensive exception hierarchy for PMF calculations,
providing structured error handling with recovery capabilities and detailed
error information.
"""

from typing import Optional, Dict, Any, List
from enum import Enum


class PMFErrorCode(Enum):
    """Standard error codes for PMF operations"""
    
    # System and validation errors
    SYSTEM_NOT_FOUND = "PMF_001"
    SYSTEM_INVALID = "PMF_002"
    SYSTEM_INCOMPLETE = "PMF_003"
    
    # File and I/O errors
    FILE_NOT_FOUND = "PMF_101"
    FILE_CORRUPTED = "PMF_102"
    FILE_PERMISSION = "PMF_103"
    DIRECTORY_NOT_WRITABLE = "PMF_104"
    
    # Configuration errors
    CONFIG_INVALID = "PMF_201"
    CONFIG_MISSING_REQUIRED = "PMF_202"
    CONFIG_VALUE_OUT_OF_RANGE = "PMF_203"
    
    # Workflow and dependency errors
    PREREQUISITE_NOT_MET = "PMF_301"
    STEP_ALREADY_COMPLETED = "PMF_302"
    WORKFLOW_STATE_INCONSISTENT = "PMF_303"
    
    # SMD specific errors
    SMD_PREPARATION_FAILED = "PMF_401"
    SMD_EXECUTION_FAILED = "PMF_402"
    SMD_RESULTS_INVALID = "PMF_403"
    
    # Umbrella sampling errors
    UMBRELLA_WINDOW_FAILED = "PMF_501"
    UMBRELLA_INCOMPLETE = "PMF_502"
    UMBRELLA_CONVERGENCE_FAILED = "PMF_503"
    
    # Analysis errors
    WHAM_EXECUTION_FAILED = "PMF_601"
    WHAM_CONVERGENCE_FAILED = "PMF_602"
    ANALYSIS_DATA_INSUFFICIENT = "PMF_603"
    
    # External tool errors
    GROMACS_NOT_FOUND = "PMF_701"
    GROMACS_EXECUTION_FAILED = "PMF_702"
    GROMACS_VERSION_INCOMPATIBLE = "PMF_703"
    
    # Resource errors
    INSUFFICIENT_DISK_SPACE = "PMF_801"
    INSUFFICIENT_MEMORY = "PMF_802"
    TIMEOUT_EXCEEDED = "PMF_803"


class PMFError(Exception):
    """
    Base exception class for all PMF-related errors.
    
    Provides structured error information including error codes,
    recovery suggestions, and context data.
    """
    
    def __init__(
        self,
        message: str,
        error_code: PMFErrorCode = None,
        recoverable: bool = False,
        recovery_suggestions: Optional[List[str]] = None,
        context: Optional[Dict[str, Any]] = None,
        cause: Optional[Exception] = None
    ):
        """
        Initialize PMF error.
        
        Parameters:
        -----------
        message : str
            Human-readable error message
        error_code : PMFErrorCode, optional
            Standardized error code for programmatic handling
        recoverable : bool, optional
            Whether this error can potentially be recovered from
        recovery_suggestions : List[str], optional
            List of suggested recovery actions
        context : Dict[str, Any], optional
            Additional context information
        cause : Exception, optional
            Original exception that caused this error
        """
        self.message = message
        self.error_code = error_code
        self.recoverable = recoverable
        self.recovery_suggestions = recovery_suggestions or []
        self.context = context or {}
        self.cause = cause
        
        super().__init__(message)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert error to dictionary for logging/serialization"""
        return {
            'message': self.message,
            'error_code': self.error_code.value if self.error_code else None,
            'recoverable': self.recoverable,
            'recovery_suggestions': self.recovery_suggestions,
            'context': self.context,
            'cause': str(self.cause) if self.cause else None
        }
    
    def __str__(self) -> str:
        """String representation with error code"""
        code_str = f"[{self.error_code.value}] " if self.error_code else ""
        return f"{code_str}{self.message}"


# System and Validation Errors
class PMFSystemError(PMFError):
    """Errors related to PMF system setup and validation"""
    pass


class SystemNotFoundException(PMFSystemError):
    """PMF system directory not found"""
    
    def __init__(self, system_path: str, searched_paths: Optional[List[str]] = None):
        message = f"PMF system not found at: {system_path}"
        context = {
            'system_path': system_path,
            'searched_paths': searched_paths or []
        }
        recovery_suggestions = [
            "Verify the system directory path is correct",
            "Check if MD simulation has completed successfully",
            "Run PMF system builder first if using PMF-optimized workflow"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.SYSTEM_NOT_FOUND,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


class SystemValidationError(PMFSystemError):
    """PMF system failed validation checks"""
    
    def __init__(self, validation_score: float, missing_files: List[str], 
                 errors: List[str], warnings: List[str]):
        message = f"System validation failed (score: {validation_score:.1f}%)"
        context = {
            'validation_score': validation_score,
            'missing_files': missing_files,
            'errors': errors,
            'warnings': warnings
        }
        recovery_suggestions = [
            "Check MD simulation completion status",
            "Verify all required files are present",
            "Re-run failed simulation steps if necessary"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.SYSTEM_INVALID,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# File and I/O Errors
class PMFFileError(PMFError):
    """File and I/O related errors"""
    pass


class RequiredFileNotFoundError(PMFFileError):
    """Required file missing for PMF operation"""
    
    def __init__(self, file_path: str, step: str, alternative_files: Optional[List[str]] = None):
        message = f"Required file not found for {step}: {file_path}"
        context = {
            'file_path': file_path,
            'step': step,
            'alternative_files': alternative_files or []
        }
        recovery_suggestions = [
            f"Ensure {step} step has completed successfully",
            "Check file paths and permissions",
            "Re-run the previous step if necessary"
        ]
        
        if alternative_files:
            recovery_suggestions.append(f"Alternative files: {', '.join(alternative_files)}")
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.FILE_NOT_FOUND,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


class FileCorruptionError(PMFFileError):
    """File appears to be corrupted or invalid"""
    
    def __init__(self, file_path: str, expected_format: str, actual_issue: str):
        message = f"File corruption detected in {file_path}: {actual_issue}"
        context = {
            'file_path': file_path,
            'expected_format': expected_format,
            'actual_issue': actual_issue
        }
        recovery_suggestions = [
            "Re-generate the corrupted file",
            "Check disk space and file system integrity",
            "Verify the generating process completed successfully"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.FILE_CORRUPTED,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# Configuration Errors
class PMFConfigurationError(PMFError):
    """Configuration-related errors"""
    pass


class ConfigurationValidationError(PMFConfigurationError):
    """Configuration validation failed"""
    
    def __init__(self, invalid_params: Dict[str, str], config_file: Optional[str] = None):
        message = f"Configuration validation failed: {len(invalid_params)} invalid parameters"
        context = {
            'invalid_params': invalid_params,
            'config_file': config_file
        }
        recovery_suggestions = [
            "Check configuration parameter values and types",
            "Refer to documentation for valid parameter ranges",
            "Use configuration validation tools"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.CONFIG_INVALID,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# Workflow and Dependency Errors
class PMFWorkflowError(PMFError):
    """Workflow execution and dependency errors"""
    pass


class PrerequisiteNotMetError(PMFWorkflowError):
    """Required prerequisite step not completed"""
    
    def __init__(self, current_step: str, required_step: str, 
                 required_files: Optional[List[str]] = None):
        message = f"Cannot execute {current_step}: {required_step} not completed"
        context = {
            'current_step': current_step,
            'required_step': required_step,
            'required_files': required_files or []
        }
        recovery_suggestions = [
            f"Complete {required_step} step first",
            "Check workflow status and file integrity",
            "Use workflow checker to diagnose issues"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.PREREQUISITE_NOT_MET,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


class WorkflowStateError(PMFWorkflowError):
    """Workflow state is inconsistent or corrupted"""
    
    def __init__(self, expected_state: str, actual_state: str, state_file: Optional[str] = None):
        message = f"Workflow state inconsistent: expected {expected_state}, got {actual_state}"
        context = {
            'expected_state': expected_state,
            'actual_state': actual_state,
            'state_file': state_file
        }
        recovery_suggestions = [
            "Reset workflow state and restart from last known good step",
            "Check state file integrity",
            "Use workflow status checker for diagnosis"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.WORKFLOW_STATE_INCONSISTENT,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# SMD Specific Errors
class SMDError(PMFError):
    """SMD-specific errors"""
    pass


class SMDPreparationError(SMDError):
    """SMD preparation failed"""
    
    def __init__(self, reason: str, missing_components: Optional[List[str]] = None):
        message = f"SMD preparation failed: {reason}"
        context = {
            'reason': reason,
            'missing_components': missing_components or []
        }
        recovery_suggestions = [
            "Check PMF system building completion",
            "Verify system equilibration status",
            "Check GROMACS installation and tools"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.SMD_PREPARATION_FAILED,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


class SMDExecutionError(SMDError):
    """SMD execution failed"""
    
    def __init__(self, exit_code: int, stderr_output: str, 
                 simulation_time: Optional[float] = None):
        message = f"SMD execution failed with exit code {exit_code}"
        context = {
            'exit_code': exit_code,
            'stderr_output': stderr_output,
            'simulation_time': simulation_time
        }
        recovery_suggestions = [
            "Check GROMACS error messages in stderr",
            "Verify system stability and parameters",
            "Check for hardware/resource issues",
            "Try reducing pull rate or adjusting parameters"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.SMD_EXECUTION_FAILED,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# Umbrella Sampling Errors
class UmbrellaError(PMFError):
    """Umbrella sampling errors"""
    pass


class UmbrellaIncompleteError(UmbrellaError):
    """Umbrella sampling incomplete or insufficient"""
    
    def __init__(self, completed_windows: int, total_windows: int, 
                 failed_windows: Optional[List[str]] = None):
        completion_rate = (completed_windows / total_windows * 100) if total_windows > 0 else 0
        message = f"Umbrella sampling incomplete: {completed_windows}/{total_windows} windows ({completion_rate:.1f}%)"
        
        context = {
            'completed_windows': completed_windows,
            'total_windows': total_windows,
            'completion_rate': completion_rate,
            'failed_windows': failed_windows or []
        }
        recovery_suggestions = [
            "Re-run failed umbrella windows",
            "Check individual window error logs",
            "Consider adjusting umbrella parameters",
            "Ensure sufficient computational resources"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.UMBRELLA_INCOMPLETE,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# Analysis Errors
class PMFAnalysisError(PMFError):
    """PMF analysis errors"""
    pass


class WHAMExecutionError(PMFAnalysisError):
    """WHAM analysis execution failed"""
    
    def __init__(self, command: str, exit_code: int, stderr_output: str):
        message = f"WHAM analysis failed with exit code {exit_code}"
        context = {
            'command': command,
            'exit_code': exit_code,
            'stderr_output': stderr_output
        }
        recovery_suggestions = [
            "Check WHAM error messages",
            "Verify umbrella sampling data quality",
            "Check WHAM input file format",
            "Try adjusting WHAM parameters"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.WHAM_EXECUTION_FAILED,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


class WHAMConvergenceError(PMFAnalysisError):
    """WHAM analysis failed to converge"""
    
    def __init__(self, iterations: int, tolerance: float, final_error: float):
        message = f"WHAM failed to converge after {iterations} iterations (error: {final_error:.2e}, tolerance: {tolerance:.2e})"
        context = {
            'iterations': iterations,
            'tolerance': tolerance,
            'final_error': final_error
        }
        recovery_suggestions = [
            "Increase maximum iterations",
            "Adjust convergence tolerance",
            "Check umbrella sampling overlap",
            "Verify data quality and completeness"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.WHAM_CONVERGENCE_FAILED,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# External Tool Errors
class ExternalToolError(PMFError):
    """External tool execution errors"""
    pass


class GromacsNotFoundError(ExternalToolError):
    """GROMACS tools not found or not accessible"""
    
    def __init__(self, missing_tools: List[str], search_paths: Optional[List[str]] = None):
        message = f"GROMACS tools not found: {', '.join(missing_tools)}"
        context = {
            'missing_tools': missing_tools,
            'search_paths': search_paths or []
        }
        recovery_suggestions = [
            "Install GROMACS and ensure it's in PATH",
            "Source GROMACS environment setup script",
            "Check GROMACS installation directory",
            "Verify executable permissions"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.GROMACS_NOT_FOUND,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


class GromacsExecutionError(ExternalToolError):
    """GROMACS tool execution failed"""
    
    def __init__(self, tool: str, command: str, exit_code: int, stderr_output: str):
        message = f"GROMACS {tool} failed with exit code {exit_code}"
        context = {
            'tool': tool,
            'command': command,
            'exit_code': exit_code,
            'stderr_output': stderr_output
        }
        recovery_suggestions = [
            f"Check {tool} error messages",
            "Verify input file formats and parameters",
            "Check system resources and permissions",
            "Try alternative GROMACS parameters"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.GROMACS_EXECUTION_FAILED,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# Resource Errors
class ResourceError(PMFError):
    """Resource-related errors (disk, memory, etc.)"""
    pass


class InsufficientDiskSpaceError(ResourceError):
    """Insufficient disk space for operation"""
    
    def __init__(self, required_space_gb: float, available_space_gb: float, operation: str):
        message = f"Insufficient disk space for {operation}: need {required_space_gb:.1f}GB, have {available_space_gb:.1f}GB"
        context = {
            'required_space_gb': required_space_gb,
            'available_space_gb': available_space_gb,
            'operation': operation
        }
        recovery_suggestions = [
            "Free up disk space",
            "Clean up temporary files",
            "Move operation to different disk",
            "Compress or archive old data"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.INSUFFICIENT_DISK_SPACE,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


class TimeoutError(ResourceError):
    """Operation exceeded timeout limit"""
    
    def __init__(self, operation: str, timeout_seconds: int, elapsed_seconds: Optional[int] = None):
        message = f"Operation '{operation}' timed out after {timeout_seconds}s"
        context = {
            'operation': operation,
            'timeout_seconds': timeout_seconds,
            'elapsed_seconds': elapsed_seconds
        }
        recovery_suggestions = [
            "Increase timeout limit",
            "Check system performance",
            "Consider running on more powerful hardware",
            "Split operation into smaller chunks"
        ]
        
        super().__init__(
            message=message,
            error_code=PMFErrorCode.TIMEOUT_EXCEEDED,
            recoverable=True,
            recovery_suggestions=recovery_suggestions,
            context=context
        )


# Utility functions for error handling
def create_error_from_exception(exc: Exception, context: Optional[Dict[str, Any]] = None) -> PMFError:
    """Convert generic exception to PMFError with context"""
    if isinstance(exc, PMFError):
        return exc
    
    # Map common exceptions to PMF errors
    if isinstance(exc, FileNotFoundError):
        return RequiredFileNotFoundError(
            file_path=str(exc.filename) if exc.filename else "unknown",
            step="unknown"
        )
    elif isinstance(exc, PermissionError):
        return PMFFileError(
            message=f"Permission denied: {exc}",
            error_code=PMFErrorCode.FILE_PERMISSION,
            recoverable=True,
            context=context
        )
    elif isinstance(exc, TimeoutError):
        return TimeoutError(
            operation="unknown",
            timeout_seconds=0
        )
    else:
        return PMFError(
            message=f"Unexpected error: {exc}",
            recoverable=False,
            context=context,
            cause=exc
        )


def format_error_report(error: PMFError) -> str:
    """Format PMFError as detailed report"""
    lines = []
    lines.append(f"ERROR: {error}")
    
    if error.error_code:
        lines.append(f"Code: {error.error_code.value}")
    
    lines.append(f"Recoverable: {'Yes' if error.recoverable else 'No'}")
    
    if error.recovery_suggestions:
        lines.append("\nSuggested Actions:")
        for suggestion in error.recovery_suggestions:
            lines.append(f"  • {suggestion}")
    
    if error.context:
        lines.append("\nContext:")
        for key, value in error.context.items():
            lines.append(f"  {key}: {value}")
    
    if error.cause:
        lines.append(f"\nOriginal Error: {error.cause}")
    
    return "\n".join(lines)