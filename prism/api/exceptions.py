#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM API Exceptions

Custom exception classes for PRISM API operations.
"""


class APIError(Exception):
    """Base exception for PRISM API operations"""
    
    def __init__(self, message: str, error_code: str = None, details: dict = None):
        super().__init__(message)
        self.message = message
        self.error_code = error_code or "API_ERROR"
        self.details = details or {}
    
    def to_dict(self) -> dict:
        """Convert exception to dictionary format"""
        return {
            'error': self.error_code,
            'message': self.message,
            'details': self.details
        }


class ValidationError(APIError):
    """Exception for input validation errors"""
    
    def __init__(self, message: str, field: str = None, value = None):
        super().__init__(message, "VALIDATION_ERROR")
        self.field = field
        self.value = value
        if field:
            self.details['field'] = field
        if value is not None:
            self.details['value'] = str(value)


class CalculationError(APIError):
    """Exception for calculation and simulation errors"""
    
    def __init__(self, message: str, stage: str = None, exit_code: int = None):
        super().__init__(message, "CALCULATION_ERROR")
        self.stage = stage
        self.exit_code = exit_code
        if stage:
            self.details['stage'] = stage
        if exit_code is not None:
            self.details['exit_code'] = exit_code


class ConfigurationError(APIError):
    """Exception for configuration-related errors"""
    
    def __init__(self, message: str, config_key: str = None):
        super().__init__(message, "CONFIG_ERROR")
        self.config_key = config_key
        if config_key:
            self.details['config_key'] = config_key


class ResourceError(APIError):
    """Exception for resource allocation and management errors"""
    
    def __init__(self, message: str, resource_type: str = None, required_amount = None):
        super().__init__(message, "RESOURCE_ERROR")
        self.resource_type = resource_type
        self.required_amount = required_amount
        if resource_type:
            self.details['resource_type'] = resource_type
        if required_amount is not None:
            self.details['required_amount'] = required_amount


class WorkflowError(APIError):
    """Exception for workflow execution errors"""
    
    def __init__(self, message: str, workflow_id: str = None, task_id: str = None):
        super().__init__(message, "WORKFLOW_ERROR")
        self.workflow_id = workflow_id
        self.task_id = task_id
        if workflow_id:
            self.details['workflow_id'] = workflow_id
        if task_id:
            self.details['task_id'] = task_id


class DataError(APIError):
    """Exception for data handling and format errors"""
    
    def __init__(self, message: str, data_format: str = None, file_path: str = None):
        super().__init__(message, "DATA_ERROR")
        self.data_format = data_format
        self.file_path = file_path
        if data_format:
            self.details['data_format'] = data_format
        if file_path:
            self.details['file_path'] = file_path


class TimeoutError(APIError):
    """Exception for operation timeouts"""
    
    def __init__(self, message: str, timeout_seconds: int = None, operation: str = None):
        super().__init__(message, "TIMEOUT_ERROR")
        self.timeout_seconds = timeout_seconds
        self.operation = operation
        if timeout_seconds:
            self.details['timeout_seconds'] = timeout_seconds
        if operation:
            self.details['operation'] = operation