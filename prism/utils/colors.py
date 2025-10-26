#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Color Output Utilities - Colored terminal output for better UX
"""

import sys
import os


class Colors:
    """ANSI color codes for terminal output"""

    # Basic colors
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'

    # Bright/Bold colors
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'

    # Styles
    BOLD = '\033[1m'
    DIM = '\033[2m'
    ITALIC = '\033[3m'
    UNDERLINE = '\033[4m'
    BLINK = '\033[5m'
    REVERSE = '\033[7m'

    # Reset
    RESET = '\033[0m'

    # Backgrounds
    BG_BLACK = '\033[40m'
    BG_RED = '\033[41m'
    BG_GREEN = '\033[42m'
    BG_YELLOW = '\033[43m'
    BG_BLUE = '\033[44m'
    BG_MAGENTA = '\033[45m'
    BG_CYAN = '\033[46m'
    BG_WHITE = '\033[47m'

    @staticmethod
    def is_color_supported():
        """
        Check if the terminal supports color output

        Returns:
            bool: True if colors are supported
        """
        # Check for NO_COLOR environment variable
        if os.getenv('NO_COLOR'):
            return False

        # Check for FORCE_COLOR environment variable
        if os.getenv('FORCE_COLOR'):
            return True

        # Check if output is a terminal
        if not hasattr(sys.stdout, 'isatty'):
            return False

        if not sys.stdout.isatty():
            return False

        # Check platform
        if sys.platform == 'win32':
            # Windows 10+ supports ANSI colors
            try:
                import platform
                version = platform.version()
                # Windows 10 build 10586 and later support ANSI
                if 'Windows-10' in version:
                    return True
            except:
                pass

            # Try to enable ANSI on Windows
            try:
                import ctypes
                kernel32 = ctypes.windll.kernel32
                kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)
                return True
            except:
                return False

        # Unix-like systems generally support colors
        return True


class ColorPrinter:
    """Utility class for printing colored text"""

    def __init__(self, use_colors=None):
        """
        Initialize color printer

        Args:
            use_colors: Force enable/disable colors (None = auto-detect)
        """
        if use_colors is None:
            self.use_colors = Colors.is_color_supported()
        else:
            self.use_colors = use_colors

    def _colorize(self, text, *color_codes):
        """Apply color codes to text"""
        if not self.use_colors:
            return text

        color_str = ''.join(color_codes)
        return f"{color_str}{text}{Colors.RESET}"

    # Semantic color methods
    def success(self, text):
        """Green text for success messages"""
        return self._colorize(text, Colors.BRIGHT_GREEN, Colors.BOLD)

    def error(self, text):
        """Red text for error messages"""
        return self._colorize(text, Colors.BRIGHT_RED, Colors.BOLD)

    def warning(self, text):
        """Yellow text for warning messages"""
        return self._colorize(text, Colors.BRIGHT_YELLOW, Colors.BOLD)

    def info(self, text):
        """Cyan text for info messages"""
        return self._colorize(text, Colors.CYAN)

    def header(self, text):
        """Bold blue text for headers"""
        return self._colorize(text, Colors.BRIGHT_BLUE, Colors.BOLD)

    def subheader(self, text):
        """Blue text for subheaders"""
        return self._colorize(text, Colors.BLUE, Colors.BOLD)

    def dim(self, text):
        """Dimmed text for less important info"""
        return self._colorize(text, Colors.DIM)

    def highlight(self, text):
        """Highlighted text (magenta)"""
        return self._colorize(text, Colors.BRIGHT_MAGENTA, Colors.BOLD)

    def path(self, text):
        """Path/filename highlighting"""
        return self._colorize(text, Colors.CYAN, Colors.UNDERLINE)

    def number(self, text):
        """Number highlighting"""
        return self._colorize(text, Colors.BRIGHT_YELLOW)

    # Print methods
    def print_success(self, text, prefix="✓"):
        """Print success message"""
        print(f"{self.success(prefix)} {text}")

    def print_error(self, text, prefix="✗"):
        """Print error message"""
        print(f"{self.error(prefix)} {text}")

    def print_warning(self, text, prefix="⚠"):
        """Print warning message"""
        print(f"{self.warning(prefix)} {text}")

    def print_info(self, text, prefix="ℹ"):
        """Print info message"""
        print(f"{self.info(prefix)} {text}")

    def print_header(self, text, char="=", width=80):
        """Print a header with separator line"""
        separator = self.header(char * width)
        header_text = self.header(text)
        print(f"\n{separator}")
        print(header_text)
        print(f"{separator}\n")

    def print_subheader(self, text, prefix="==="):
        """Print a subheader"""
        print(f"\n{self.subheader(prefix)} {self.subheader(text)}\n")

    def print_step(self, step_num, total_steps, description):
        """Print a step in a workflow"""
        step_text = self.highlight(f"Step {step_num}/{total_steps}:")
        print(f"{step_text} {description}")

    def print_progress(self, current, total, description=""):
        """Print progress information"""
        percentage = (current / total) * 100 if total > 0 else 0
        progress_text = self.number(f"[{current}/{total}]")
        percent_text = self.number(f"{percentage:.1f}%")
        print(f"{progress_text} {percent_text} {description}")

    def print_file_info(self, label, path):
        """Print file information"""
        print(f"  {label}: {self.path(path)}")

    def print_param_info(self, label, value):
        """Print parameter information"""
        print(f"  {label}: {self.number(str(value))}")


# Global color printer instance
_color_printer = None

def get_color_printer():
    """Get the global color printer instance"""
    global _color_printer
    if _color_printer is None:
        _color_printer = ColorPrinter()
    return _color_printer


def set_color_enabled(enabled):
    """Enable or disable colored output globally"""
    global _color_printer
    _color_printer = ColorPrinter(use_colors=enabled)


# Convenience functions
def success(text):
    """Return success-colored text"""
    return get_color_printer().success(text)

def error(text):
    """Return error-colored text"""
    return get_color_printer().error(text)

def warning(text):
    """Return warning-colored text"""
    return get_color_printer().warning(text)

def info(text):
    """Return info-colored text"""
    return get_color_printer().info(text)

def header(text):
    """Return header-colored text"""
    return get_color_printer().header(text)

def subheader(text):
    """Return subheader-colored text"""
    return get_color_printer().subheader(text)

def dim(text):
    """Return dimmed text"""
    return get_color_printer().dim(text)

def highlight(text):
    """Return highlighted text"""
    return get_color_printer().highlight(text)

def path(text):
    """Return path-colored text"""
    return get_color_printer().path(text)

def number(text):
    """Return number-colored text"""
    return get_color_printer().number(text)


# Print functions
def print_success(text, prefix="✓"):
    """Print success message"""
    get_color_printer().print_success(text, prefix)

def print_error(text, prefix="✗"):
    """Print error message"""
    get_color_printer().print_error(text, prefix)

def print_warning(text, prefix="⚠"):
    """Print warning message"""
    get_color_printer().print_warning(text, prefix)

def print_info(text, prefix="ℹ"):
    """Print info message"""
    get_color_printer().print_info(text, prefix)

def print_header(text, char="=", width=80):
    """Print header with separator"""
    get_color_printer().print_header(text, char, width)

def print_subheader(text, prefix="==="):
    """Print subheader"""
    get_color_printer().print_subheader(text, prefix)

def print_step(step_num, total_steps, description):
    """Print workflow step"""
    get_color_printer().print_step(step_num, total_steps, description)


if __name__ == '__main__':
    # Simple test
    cp = ColorPrinter()
    print("\nPRISM Color Module")
    cp.print_success("Colors are working!")
    cp.print_warning("This is a warning")
    cp.print_error("This is an error")
    cp.print_info("This is info")
