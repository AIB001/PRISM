#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Trajectory Processing Module.

This module provides trajectory file processing, format conversion,
and chain selection utilities for PRISM analysis workflows.
"""

from .processor import TrajectoryProcessor
from .manager import TrajectoryManager

__all__ = ["TrajectoryProcessor", "TrajectoryManager"]
