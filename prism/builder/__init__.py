#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Builder Package - Main builder class for protein-ligand systems

This package provides the PRISMBuilder class and CLI entry point.
Split from the original monolithic builder.py for maintainability.
"""

from .core import PRISMBuilder
from .cli import main

__all__ = ['PRISMBuilder', 'main']
