#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Logging Utilities - Unified logging with colored terminal output

Provides a custom logging handler that uses ColorPrinter for terminal output,
preserving the existing colored output while gaining logging infrastructure
(levels, filters, file output).

Usage:
    from prism.utils.log import get_logger

    logger = get_logger(__name__)
    logger.info("Processing ligand...")
    logger.warning("File already exists")
    logger.error("Failed to generate topology")
"""

import logging

from .colors import get_color_printer


class ColoredHandler(logging.StreamHandler):
    """Logging handler that uses ColorPrinter for colored terminal output."""

    # Map logging levels to ColorPrinter method names
    LEVEL_STYLE = {
        logging.DEBUG: "dim",
        logging.INFO: None,  # no color wrapping
        logging.WARNING: "warning",
        logging.ERROR: "error",
        logging.CRITICAL: "error",
    }

    def __init__(self, stream=None):
        super().__init__(stream)
        self.color_printer = get_color_printer()

    def emit(self, record):
        try:
            msg = self.format(record)
            style = self.LEVEL_STYLE.get(record.levelno)
            if style and hasattr(self.color_printer, style):
                # Apply color to the prefix only for warnings/errors
                if record.levelno == logging.WARNING:
                    msg = f"{self.color_printer.warning('⚠')} {msg}"
                elif record.levelno >= logging.ERROR:
                    msg = f"{self.color_printer.error('✗')} {msg}"
                else:
                    msg = getattr(self.color_printer, style)(msg)
            self.stream.write(msg + self.terminator)
            self.flush()
        except Exception:
            self.handleError(record)


def get_logger(name: str) -> logging.Logger:
    """
    Get a PRISM logger with colored output.

    Parameters
    ----------
    name : str
        Logger name (typically __name__)

    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = ColoredHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False  # Prevent duplicate output from root logger
    return logger


# Convenience functions for formatted output via logging
def log_header(logger, msg, char="=", width=80):
    """Log a header with separator line"""
    cp = get_color_printer()
    separator = cp.header(char * width)
    header_text = cp.header(msg)
    logger.info(f"\n{separator}\n{header_text}\n{separator}\n")


def log_subheader(logger, msg, prefix="==="):
    """Log a subheader"""
    cp = get_color_printer()
    logger.info(f"\n{cp.subheader(prefix)} {cp.subheader(msg)}\n")


def log_step(logger, step_num, total_steps, description):
    """Log a workflow step"""
    cp = get_color_printer()
    step_text = cp.highlight(f"Step {step_num}/{total_steps}:")
    logger.info(f"{step_text} {description}")


def log_success(logger, msg, prefix="✓"):
    """Log a success message"""
    cp = get_color_printer()
    logger.info(f"{cp.success(prefix)} {msg}")
