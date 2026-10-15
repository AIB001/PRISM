"""Actionable failures; an unknown scientific result is never a pass."""


class DesignError(ValueError):
    def __init__(self, message, code="INVALID_DESIGN"):
        super().__init__(message)
        self.code = code


def require(condition, message, code="INVALID_DESIGN"):
    if not condition:
        raise DesignError(message, code)
