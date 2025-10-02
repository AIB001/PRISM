#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Data Management and Storage System

This module provides comprehensive data management capabilities for PRISM
molecular dynamics simulations, including efficient storage, retrieval,
compression, and processing of large-scale simulation data.
"""

from .storage import (
    DataStorageManager,
    HierarchicalStorage,
    CompressionManager
)

from .processing import (
    DataProcessor,
    StreamProcessor,
    BatchProcessor
)

from .formats import (
    TrajectoryHandler,
    EnergyHandler,
    PMFDataHandler
)

from .cache import (
    CacheManager,
    MemoryCache,
    DiskCache
)

__all__ = [
    # Storage management
    'DataStorageManager',
    'HierarchicalStorage', 
    'CompressionManager',
    
    # Data processing
    'DataProcessor',
    'StreamProcessor',
    'BatchProcessor',
    
    # Format handlers
    'TrajectoryHandler',
    'EnergyHandler',
    'PMFDataHandler',
    
    # Caching
    'CacheManager',
    'MemoryCache',
    'DiskCache'
]