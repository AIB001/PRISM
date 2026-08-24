"""Generation-specific exceptions with stable machine-readable codes."""


class GenerationError(Exception):
    code = "GENERATION_ERROR"

    def __init__(self, message: str, code: str = None):
        super().__init__(message)
        if code is not None:
            self.code = code


class InputError(GenerationError):
    code = "INVALID_INPUT"


class ConfigurationError(GenerationError):
    code = "MODEL_NOT_CONFIGURED"


class CapabilityError(GenerationError):
    code = "UNSUPPORTED_POCKET_KIND"


class ExecutionFailure(GenerationError):
    code = "MODEL_EXECUTION_FAILED"
