#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel processing utilities for MDTraj-based analysis.

Provides centralized management of OpenMP threading for MDTraj computations
and parallel processing utilities for analysis modules.
"""

import os
import logging
from typing import Optional, Callable, Any
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class ParallelProcessor:
    """Centralized parallel processing manager for MDTraj computations"""

    def __init__(self, n_threads: Optional[int] = None):
        """
        Initialize parallel processor.

        Parameters
        ----------
        n_threads : int, optional
            Number of threads to use. If None, automatically determine optimal count.
        """
        self.original_omp_threads = os.environ.get('OMP_NUM_THREADS')
        self.n_threads = n_threads or self._get_optimal_thread_count()

    def _get_optimal_thread_count(self) -> int:
        """Determine optimal thread count based on CPU availability"""
        cpu_count = os.cpu_count()
        if cpu_count is None:
            return 1
        # Leave one core for system tasks, but ensure at least 1 thread
        return max(1, cpu_count - 1)

    @contextmanager
    def configure_omp_threads(self):
        """
        Context manager to temporarily configure OpenMP threads.

        Sets OMP_NUM_THREADS environment variable for MDTraj parallel processing
        and restores original setting when done.
        """
        # Save current setting
        original_setting = os.environ.get('OMP_NUM_THREADS')

        # Configure threads
        os.environ['OMP_NUM_THREADS'] = str(self.n_threads)
        logger.debug(f"Configured OpenMP with {self.n_threads} threads")

        try:
            yield self.n_threads
        finally:
            # Restore original setting
            if original_setting is not None:
                os.environ['OMP_NUM_THREADS'] = original_setting
            elif 'OMP_NUM_THREADS' in os.environ:
                del os.environ['OMP_NUM_THREADS']
            logger.debug("Restored original OpenMP thread configuration")

    def get_thread_info(self) -> dict:
        """
        Get information about current thread configuration.

        Returns
        -------
        dict
            Thread configuration information
        """
        return {
            'configured_threads': self.n_threads,
            'cpu_count': os.cpu_count(),
            'original_omp_threads': self.original_omp_threads,
            'current_omp_threads': os.environ.get('OMP_NUM_THREADS')
        }


def parallel_compute(func: Callable, *args, **kwargs) -> Any:
    """
    Execute a computation with automatic parallel processing configuration.

    Parameters
    ----------
    func : callable
        Function to execute
    *args, **kwargs
        Arguments to pass to function

    Returns
    -------
    Any
        Result of function execution
    """
    processor = ParallelProcessor()
    with processor.configure_omp_threads():
        return func(*args, **kwargs)


# Global instance for convenience
default_processor = ParallelProcessor()