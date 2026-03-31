#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Backward-compatible exports for FEP XVG parsing.

This wrapper preserves the older import path:
`prism.fep.analysis.xvg_parser`.
"""

from .core.xvg_parser import XVGParser, find_xvg_file, parse_leg_data

__all__ = [
    "XVGParser",
    "find_xvg_file",
    "parse_leg_data",
]
