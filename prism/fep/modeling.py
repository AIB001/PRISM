#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM-FEP modeling - Backward compatibility wrapper.

This module has been reorganized into a package structure.
All classes are re-exported for backward compatibility.
"""

from .modeling import FEPScaffoldBuilder, FEPScaffoldLayout

__all__ = ["FEPScaffoldBuilder", "FEPScaffoldLayout"]
