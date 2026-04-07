#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
JavaScript code loader for contact visualization.

Loads JavaScript code from separate asset files.
"""

from pathlib import Path


# Get the directory where this file is located
_ASSETS_DIR = Path(__file__).parent / "assets"


def _load_js_file(filename: str) -> str:
    """Load JavaScript code from asset file.

    Args:
        filename: Name of the JavaScript file

    Returns:
        JavaScript code as string
    """
    js_path = _ASSETS_DIR / filename
    return js_path.read_text(encoding="utf-8")


def get_contact_map_class_part1() -> str:
    """First part of ContactMap class - initialization and events"""
    return _load_js_file("contact_map_part1.js")


def get_contact_map_class_part2() -> str:
    """Second part of ContactMap class - core functionality"""
    return _load_js_file("contact_map_part2.js")


def get_drawing_methods() -> str:
    """Drawing methods for visualization"""
    return _load_js_file("drawing_methods.js")


def get_utility_functions() -> str:
    """Utility functions for the visualization"""
    return _load_js_file("utility_functions.js")


def get_javascript_code() -> str:
    """Return the complete JavaScript code as a string.

    Combines all parts into a single JavaScript code block.

    Returns:
        Complete JavaScript code for HTML embedding
    """
    return (
        get_contact_map_class_part1()
        + get_contact_map_class_part2()
        + get_drawing_methods()
        + get_utility_functions()
        + r"""
    </script>
</body>
</html>"""
    )
