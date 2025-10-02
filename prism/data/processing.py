#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Data Processing System

High-performance data processing capabilities for molecular dynamics simulations,
including stream processing, batch processing, and parallel data operations.
"""

import os
import sys
import threading
import multiprocessing as mp
import queue
import time
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Callable, Iterator, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
from abc import ABC, abstractmethod
import numpy as np
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from contextlib import contextmanager

from ..utils.logging_system import PrismLogger


class ProcessingMode(Enum):
    """Data processing modes"""
    SEQUENTIAL = "sequential"
    THREADED = "threaded" 
    MULTIPROCESS = "multiprocess"
    HYBRID = "hybrid"


class DataFormat(Enum):
    """Supported data formats"""
    NUMPY = "numpy"
    PANDAS = "pandas"
    JSON = "json"
    BINARY = "binary"
    TRAJECTORY = "trajectory"
    ENERGY = "energy"
    PMF = "pmf"


@dataclass
class ProcessingJob:
    """Data processing job specification"""
    job_id: str
    input_data: Any
    processor_func: Callable
    processor_args: Tuple = ()
    processor_kwargs: Dict[str, Any] = None
    priority: int = 0
    timeout: Optional[float] = None
    retry_count: int = 0
    max_retries: int = 3
    
    def __post_init__(self):
        if self.processor_kwargs is None:
            self.processor_kwargs = {}


@dataclass
class ProcessingResult:
    """Data processing result"""
    job_id: str
    status: str  # "completed", "failed", "timeout"
    output_data: Any = None
    error_message: Optional[str] = None
    processing_time: float = 0.0
    memory_peak_mb: float = 0.0
    retry_count: int = 0


class DataProcessor(ABC):
    """Abstract base class for data processors"""
    
    def __init__(self, name: str):
        self.name = name
        self.logger = PrismLogger(f"processor.{name}")
        self.processing_stats = {
            'jobs_completed': 0,
            'jobs_failed': 0,
            'total_processing_time': 0.0,
            'average_processing_time': 0.0
        }
    
    @abstractmethod
    def process(self, data: Any, *args, **kwargs) -> Any:
        """Process data and return result"""
        pass
    
    def validate_input(self, data: Any) -> bool:
        """Validate input data format"""
        return True
    
    def get_stats(self) -> Dict[str, Any]:
        """Get processing statistics"""
        return dict(self.processing_stats)
    
    def _update_stats(self, processing_time: float, success: bool):
        """Update processing statistics"""
        if success:
            self.processing_stats['jobs_completed'] += 1
        else:
            self.processing_stats['jobs_failed'] += 1
        
        self.processing_stats['total_processing_time'] += processing_time
        total_jobs = self.processing_stats['jobs_completed'] + self.processing_stats['jobs_failed']
        
        if total_jobs > 0:
            self.processing_stats['average_processing_time'] = \
                self.processing_stats['total_processing_time'] / total_jobs


class StreamProcessor:
    """High-performance stream data processor"""
    
    def __init__(self, buffer_size: int = 1000, max_workers: int = None):
        self.buffer_size = buffer_size
        self.max_workers = max_workers or min(32, (os.cpu_count() or 1) + 4)
        self.logger = PrismLogger("stream_processor")
        
        # Processing queues
        self.input_queue = queue.Queue(maxsize=buffer_size)
        self.output_queue = queue.Queue(maxsize=buffer_size)
        self.error_queue = queue.Queue()
        
        # Worker management
        self.workers = []
        self.shutdown_event = threading.Event()
        self.processing_stats = {
            'items_processed': 0,
            'items_failed': 0,
            'processing_rate': 0.0,
            'queue_size': 0
        }
        
        # Rate tracking
        self._start_time = time.time()
        self._last_stats_update = time.time()
    
    def start_processing(self, processor: DataProcessor, mode: ProcessingMode = ProcessingMode.THREADED):
        """Start stream processing with specified mode"""
        self.logger.info(f"Starting stream processing with {self.max_workers} workers ({mode.value} mode)")
        
        if mode == ProcessingMode.THREADED:
            self._start_threaded_workers(processor)
        elif mode == ProcessingMode.MULTIPROCESS:
            self._start_multiprocess_workers(processor)
        else:
            raise ValueError(f"Unsupported processing mode: {mode}")
    
    def _start_threaded_workers(self, processor: DataProcessor):
        """Start threaded workers"""
        def worker():
            while not self.shutdown_event.is_set():
                try:
                    # Get job from queue
                    job = self.input_queue.get(timeout=1.0)
                    if job is None:  # Poison pill
                        break
                    
                    # Process data
                    start_time = time.time()
                    try:
                        result = processor.process(job.input_data, *job.processor_args, **job.processor_kwargs)
                        processing_time = time.time() - start_time
                        
                        # Create result
                        processing_result = ProcessingResult(
                            job_id=job.job_id,
                            status="completed",
                            output_data=result,
                            processing_time=processing_time
                        )
                        
                        self.output_queue.put(processing_result)
                        self.processing_stats['items_processed'] += 1
                        
                        processor._update_stats(processing_time, True)
                        
                    except Exception as e:
                        processing_time = time.time() - start_time
                        
                        processing_result = ProcessingResult(
                            job_id=job.job_id,
                            status="failed",
                            error_message=str(e),
                            processing_time=processing_time
                        )
                        
                        self.error_queue.put(processing_result)
                        self.processing_stats['items_failed'] += 1
                        
                        processor._update_stats(processing_time, False)
                        self.logger.error(f"Processing failed for job {job.job_id}: {e}")
                    
                    finally:
                        self.input_queue.task_done()
                
                except queue.Empty:
                    continue
                except Exception as e:
                    self.logger.error(f"Worker error: {e}")
        
        # Start worker threads
        for i in range(self.max_workers):
            worker_thread = threading.Thread(target=worker, name=f"StreamWorker-{i}")
            worker_thread.daemon = True
            worker_thread.start()
            self.workers.append(worker_thread)
    
    def submit_job(self, job: ProcessingJob, timeout: Optional[float] = None) -> bool:
        """Submit job for processing"""
        try:
            self.input_queue.put(job, timeout=timeout)
            return True
        except queue.Full:
            self.logger.warning(f"Input queue full, could not submit job {job.job_id}")
            return False
    
    def get_result(self, timeout: Optional[float] = None) -> Optional[ProcessingResult]:
        """Get processing result"""
        try:
            return self.output_queue.get(timeout=timeout)
        except queue.Empty:
            return None
    
    def get_error(self, timeout: Optional[float] = None) -> Optional[ProcessingResult]:
        """Get processing error"""
        try:
            return self.error_queue.get(timeout=timeout)
        except queue.Empty:
            return None
    
    def get_stats(self) -> Dict[str, Any]:
        """Get processing statistics"""
        current_time = time.time()
        elapsed_time = current_time - self._start_time
        
        # Calculate processing rate
        if elapsed_time > 0:
            self.processing_stats['processing_rate'] = \
                self.processing_stats['items_processed'] / elapsed_time
        
        self.processing_stats['queue_size'] = self.input_queue.qsize()
        
        return dict(self.processing_stats)
    
    def shutdown(self, timeout: float = 10.0):
        """Shutdown stream processor"""
        self.logger.info("Shutting down stream processor")
        
        # Signal shutdown
        self.shutdown_event.set()
        
        # Send poison pills to workers
        for _ in self.workers:
            try:
                self.input_queue.put(None, timeout=1.0)
            except queue.Full:
                pass
        
        # Wait for workers to finish
        for worker in self.workers:
            worker.join(timeout=timeout / len(self.workers))
        
        self.logger.info("Stream processor shutdown complete")


class BatchProcessor:
    """Batch data processor with parallel execution"""
    
    def __init__(self, max_workers: int = None, chunk_size: int = 1000):
        self.max_workers = max_workers or os.cpu_count()
        self.chunk_size = chunk_size
        self.logger = PrismLogger("batch_processor")
    
    def process_batch(self, data_items: List[Any], processor: DataProcessor,
                     mode: ProcessingMode = ProcessingMode.MULTIPROCESS,
                     progress_callback: Optional[Callable] = None) -> List[ProcessingResult]:
        """Process batch of data items"""
        
        total_items = len(data_items)
        self.logger.info(f"Processing batch of {total_items} items using {mode.value} mode")
        
        if mode == ProcessingMode.SEQUENTIAL:
            return self._process_sequential(data_items, processor, progress_callback)
        elif mode == ProcessingMode.THREADED:
            return self._process_threaded(data_items, processor, progress_callback)
        elif mode == ProcessingMode.MULTIPROCESS:
            return self._process_multiprocess(data_items, processor, progress_callback)
        else:
            raise ValueError(f"Unsupported processing mode: {mode}")
    
    def _process_sequential(self, data_items: List[Any], processor: DataProcessor,
                          progress_callback: Optional[Callable] = None) -> List[ProcessingResult]:
        """Process data sequentially"""
        results = []
        
        for i, item in enumerate(data_items):
            start_time = time.time()
            
            try:
                output = processor.process(item)
                processing_time = time.time() - start_time
                
                result = ProcessingResult(
                    job_id=f"batch_item_{i}",
                    status="completed",
                    output_data=output,
                    processing_time=processing_time
                )
                
                processor._update_stats(processing_time, True)
                
            except Exception as e:
                processing_time = time.time() - start_time
                
                result = ProcessingResult(
                    job_id=f"batch_item_{i}",
                    status="failed",
                    error_message=str(e),
                    processing_time=processing_time
                )
                
                processor._update_stats(processing_time, False)
            
            results.append(result)
            
            if progress_callback:
                progress_callback(i + 1, len(data_items))
        
        return results
    
    def _process_threaded(self, data_items: List[Any], processor: DataProcessor,
                        progress_callback: Optional[Callable] = None) -> List[ProcessingResult]:
        """Process data using thread pool"""
        results = [None] * len(data_items)
        completed_count = 0
        
        def process_item(item_index, item):
            start_time = time.time()
            
            try:
                output = processor.process(item)
                processing_time = time.time() - start_time
                
                result = ProcessingResult(
                    job_id=f"batch_item_{item_index}",
                    status="completed",
                    output_data=output,
                    processing_time=processing_time
                )
                
                processor._update_stats(processing_time, True)
                return item_index, result
                
            except Exception as e:
                processing_time = time.time() - start_time
                
                result = ProcessingResult(
                    job_id=f"batch_item_{item_index}",
                    status="failed",
                    error_message=str(e),
                    processing_time=processing_time
                )
                
                processor._update_stats(processing_time, False)
                return item_index, result
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            futures = {
                executor.submit(process_item, i, item): i 
                for i, item in enumerate(data_items)
            }
            
            # Collect results
            for future in as_completed(futures):
                item_index, result = future.result()
                results[item_index] = result
                completed_count += 1
                
                if progress_callback:
                    progress_callback(completed_count, len(data_items))
        
        return results
    
    def _process_multiprocess(self, data_items: List[Any], processor: DataProcessor,
                            progress_callback: Optional[Callable] = None) -> List[ProcessingResult]:
        """Process data using process pool"""
        results = [None] * len(data_items)
        completed_count = 0
        
        # Create worker function that can be pickled
        def worker_function(item_data):
            item_index, item = item_data
            start_time = time.time()
            
            try:
                # Create new processor instance in worker process
                worker_processor = type(processor)(processor.name)
                output = worker_processor.process(item)
                processing_time = time.time() - start_time
                
                return item_index, ProcessingResult(
                    job_id=f"batch_item_{item_index}",
                    status="completed",
                    output_data=output,
                    processing_time=processing_time
                )
                
            except Exception as e:
                processing_time = time.time() - start_time
                
                return item_index, ProcessingResult(
                    job_id=f"batch_item_{item_index}",
                    status="failed",
                    error_message=str(e),
                    processing_time=processing_time
                )
        
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            futures = {
                executor.submit(worker_function, (i, item)): i
                for i, item in enumerate(data_items)
            }
            
            # Collect results
            for future in as_completed(futures):
                item_index, result = future.result()
                results[item_index] = result
                completed_count += 1
                
                # Update processor stats in main process
                if result.status == "completed":
                    processor._update_stats(result.processing_time, True)
                else:
                    processor._update_stats(result.processing_time, False)
                
                if progress_callback:
                    progress_callback(completed_count, len(data_items))
        
        return results
    
    def process_large_dataset(self, data_iterator: Iterator[Any], 
                            processor: DataProcessor,
                            output_handler: Callable[[List[ProcessingResult]], None],
                            mode: ProcessingMode = ProcessingMode.MULTIPROCESS) -> Dict[str, Any]:
        """Process large dataset in chunks"""
        
        total_processed = 0
        total_failed = 0
        start_time = time.time()
        
        chunk = []
        
        for item in data_iterator:
            chunk.append(item)
            
            if len(chunk) >= self.chunk_size:
                # Process chunk
                results = self.process_batch(chunk, processor, mode)
                
                # Handle results
                output_handler(results)
                
                # Update counters
                completed = sum(1 for r in results if r.status == "completed")
                failed = len(results) - completed
                
                total_processed += completed
                total_failed += failed
                
                self.logger.info(f"Processed chunk: {completed} completed, {failed} failed")
                
                # Reset chunk
                chunk = []
        
        # Process remaining items
        if chunk:
            results = self.process_batch(chunk, processor, mode)
            output_handler(results)
            
            completed = sum(1 for r in results if r.status == "completed")
            failed = len(results) - completed
            
            total_processed += completed
            total_failed += failed
        
        total_time = time.time() - start_time
        
        return {
            'total_processed': total_processed,
            'total_failed': total_failed,
            'total_time': total_time,
            'processing_rate': total_processed / total_time if total_time > 0 else 0
        }


class TrajectoryProcessor(DataProcessor):
    """Specialized processor for trajectory data"""
    
    def __init__(self):
        super().__init__("trajectory_processor")
    
    def process(self, trajectory_data: np.ndarray, *args, **kwargs) -> Dict[str, Any]:
        """Process trajectory data"""
        if not isinstance(trajectory_data, np.ndarray):
            raise ValueError("Expected numpy array for trajectory data")
        
        # Calculate basic properties
        n_frames, n_atoms, n_dims = trajectory_data.shape
        
        # Calculate center of mass trajectory
        masses = kwargs.get('masses', np.ones(n_atoms))
        com_trajectory = np.average(trajectory_data, axis=1, weights=masses)
        
        # Calculate RMSD from first frame
        ref_frame = trajectory_data[0]
        rmsd_values = np.sqrt(np.mean((trajectory_data - ref_frame)**2, axis=(1, 2)))
        
        # Calculate radius of gyration
        rg_values = []
        for frame in trajectory_data:
            com_frame = np.average(frame, axis=0, weights=masses)
            distances_sq = np.sum((frame - com_frame)**2, axis=1)
            rg = np.sqrt(np.average(distances_sq, weights=masses))
            rg_values.append(rg)
        
        return {
            'n_frames': n_frames,
            'n_atoms': n_atoms,
            'center_of_mass': com_trajectory,
            'rmsd': rmsd_values,
            'radius_of_gyration': rg_values,
            'trajectory_length_ns': n_frames * kwargs.get('dt', 0.002) / 1000  # Convert ps to ns
        }
    
    def validate_input(self, data: Any) -> bool:
        """Validate trajectory input data"""
        if not isinstance(data, np.ndarray):
            return False
        
        if len(data.shape) != 3:
            return False
        
        if data.shape[2] != 3:  # Should be 3D coordinates
            return False
        
        return True


class EnergyProcessor(DataProcessor):
    """Specialized processor for energy data"""
    
    def __init__(self):
        super().__init__("energy_processor")
    
    def process(self, energy_data: Dict[str, np.ndarray], *args, **kwargs) -> Dict[str, Any]:
        """Process energy data"""
        results = {}
        
        for energy_term, values in energy_data.items():
            if not isinstance(values, np.ndarray):
                continue
            
            results[energy_term] = {
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'min': float(np.min(values)),
                'max': float(np.max(values)),
                'drift': self._calculate_drift(values),
                'autocorr_time': self._calculate_autocorrelation_time(values)
            }
        
        return results
    
    def _calculate_drift(self, values: np.ndarray) -> float:
        """Calculate energy drift over time"""
        if len(values) < 2:
            return 0.0
        
        time_points = np.arange(len(values))
        slope, _ = np.polyfit(time_points, values, 1)
        return float(slope)
    
    def _calculate_autocorrelation_time(self, values: np.ndarray) -> float:
        """Calculate autocorrelation time"""
        if len(values) < 10:
            return 0.0
        
        # Simple autocorrelation calculation
        values_normalized = values - np.mean(values)
        autocorr = np.correlate(values_normalized, values_normalized, mode='full')
        autocorr = autocorr[autocorr.size // 2:]
        autocorr = autocorr / autocorr[0]
        
        # Find where autocorrelation drops to 1/e
        threshold = 1.0 / np.e
        crossing_points = np.where(autocorr < threshold)[0]
        
        if len(crossing_points) > 0:
            return float(crossing_points[0])
        else:
            return float(len(autocorr))


# Convenience functions
def create_stream_processor(buffer_size: int = 1000, max_workers: int = None) -> StreamProcessor:
    """Create stream processor with default settings"""
    return StreamProcessor(buffer_size=buffer_size, max_workers=max_workers)


def create_batch_processor(max_workers: int = None, chunk_size: int = 1000) -> BatchProcessor:
    """Create batch processor with default settings"""
    return BatchProcessor(max_workers=max_workers, chunk_size=chunk_size)


@contextmanager
def processing_pipeline(processors: List[DataProcessor], 
                       mode: ProcessingMode = ProcessingMode.THREADED):
    """Context manager for processing pipeline"""
    stream_processors = []
    
    try:
        # Setup stream processors for each stage
        for processor in processors:
            stream_proc = create_stream_processor()
            stream_proc.start_processing(processor, mode)
            stream_processors.append(stream_proc)
        
        yield stream_processors
        
    finally:
        # Shutdown all stream processors
        for stream_proc in stream_processors:
            stream_proc.shutdown()


def parallel_map(func: Callable, data: List[Any], 
                max_workers: int = None, 
                mode: ProcessingMode = ProcessingMode.THREADED) -> List[Any]:
    """Parallel map function with different execution modes"""
    
    # Create a simple processor wrapper
    class MapProcessor(DataProcessor):
        def __init__(self, map_func):
            super().__init__("map_processor")
            self.map_func = map_func
        
        def process(self, data, *args, **kwargs):
            return self.map_func(data)
    
    processor = MapProcessor(func)
    batch_processor = BatchProcessor(max_workers=max_workers)
    
    results = batch_processor.process_batch(data, processor, mode)
    
    # Extract output data from results
    return [r.output_data for r in results if r.status == "completed"]