#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
JavaScript Code Module for Contact Visualization.

This module provides JavaScript code for interactive HTML contact visualization.
The JavaScript code is loaded from separate asset files for better maintainability.

Usage:
    from prism.analysis.contact.javascript_code import get_javascript_code
    js_code = get_javascript_code()
"""

from .loader import (
    get_javascript_code,
    get_contact_map_class_part1,
    get_contact_map_class_part2,
    get_drawing_methods,
    get_utility_functions,
)

__all__ = [
    "get_javascript_code",
    "get_contact_map_class_part1",
    "get_contact_map_class_part2",
    "get_drawing_methods",
    "get_utility_functions",
]
