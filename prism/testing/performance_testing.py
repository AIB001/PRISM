#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Performance Testing and Benchmarking

Comprehensive performance testing framework for PRISM components with
benchmarking, profiling, and performance regression detection.
"""

import time
import statistics
import tracemalloc
from pathlib import Path
from typing import Dict, Any, List, Optional, Callable, Union
from dataclasses import dataclass, asdict
from enum import Enum
import json
import psutil
import threading
from contextlib import contextmanager

from ..utils.logging_system import PrismLogger


class BenchmarkType(Enum):
    """Types of performance benchmarks"""
    EXECUTION_TIME = "execution_time"
    MEMORY_USAGE = "memory_usage"
    CPU_UTILIZATION = "cpu_utilization"
    THROUGHPUT = "throughput"
    SCALABILITY = "scalability"
    REGRESSION = "regression"


@dataclass
class PerformanceMetrics:
    """Performance metrics collection"""
    execution_time: float
    peak_memory_mb: float
    average_memory_mb: float
    cpu_percent: float
    disk_io_mb: float
    network_io_kb: Optional[float] = None
    throughput: Optional[float] = None
    custom_metrics: Dict[str, float] = None
    
    def __post_init__(self):
        if self.custom_metrics is None:
            self.custom_metrics = {}
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass 
class BenchmarkResult:
    """Result of a performance benchmark"""
    benchmark_name: str
    benchmark_type: BenchmarkType
    metrics: PerformanceMetrics
    baseline_metrics: Optional[PerformanceMetrics] = None
    performance_ratio: Optional[float] = None
    status: str = "completed"
    error_message: Optional[str] = None
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}
        
        # Calculate performance ratio if baseline is available
        if self.baseline_metrics and self.baseline_metrics.execution_time > 0:
            self.performance_ratio = self.baseline_metrics.execution_time / self.metrics.execution_time
    
    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data['benchmark_type'] = self.benchmark_type.value
        return data
    
    def is_regression(self, threshold: float = 0.1) -> bool:
        """Check if this result represents a performance regression"""
        if self.performance_ratio is None:
            return False
        
        # Performance regression if new execution time is >10% slower
        return self.performance_ratio < (1 - threshold)


class PerformanceProfiler:
    """Performance profiler for detailed analysis"""
    
    def __init__(self):
        self.logger = PrismLogger("performance_profiler")
        self.monitoring_active = False
        self.monitoring_thread = None
        self.metrics_history: List[Dict[str, float]] = []
        
    @contextmanager
    def profile_execution(self, track_memory: bool = True):
        """Context manager for profiling code execution"""
        # Start memory tracking
        if track_memory:
            tracemalloc.start()
        
        # Get process for monitoring
        process = psutil.Process()
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Start monitoring
        self.start_monitoring(process)
        
        start_time = time.time()
        
        try:
            yield self
        finally:
            execution_time = time.time() - start_time
            
            # Stop monitoring
            self.stop_monitoring()
            
            # Get final memory info
            final_memory = process.memory_info().rss / 1024 / 1024  # MB
            
            # Get memory peak if tracking
            peak_memory = final_memory
            if track_memory:
                current, peak = tracemalloc.get_traced_memory()
                peak_memory = max(peak_memory, peak / 1024 / 1024)  # MB
                tracemalloc.stop()
            
            # Calculate averages from monitoring
            avg_cpu = statistics.mean([m['cpu_percent'] for m in self.metrics_history]) if self.metrics_history else 0
            avg_memory = statistics.mean([m['memory_mb'] for m in self.metrics_history]) if self.metrics_history else final_memory
            
            # Create metrics
            self.last_metrics = PerformanceMetrics(
                execution_time=execution_time,
                peak_memory_mb=peak_memory,
                average_memory_mb=avg_memory,
                cpu_percent=avg_cpu,
                disk_io_mb=0,  # Would need more detailed tracking
                custom_metrics={
                    'initial_memory_mb': initial_memory,
                    'final_memory_mb': final_memory,
                    'memory_delta_mb': final_memory - initial_memory
                }
            )
    
    def start_monitoring(self, process: psutil.Process, interval: float = 0.1):
        """Start background monitoring of system metrics"""
        self.monitoring_active = True
        self.metrics_history.clear()
        
        def monitor():
            while self.monitoring_active:
                try:
                    cpu_percent = process.cpu_percent()
                    memory_info = process.memory_info()
                    memory_mb = memory_info.rss / 1024 / 1024
                    
                    self.metrics_history.append({
                        'timestamp': time.time(),
                        'cpu_percent': cpu_percent,
                        'memory_mb': memory_mb
                    })
                    
                    time.sleep(interval)
                except Exception as e:
                    self.logger.warning(f"Monitoring error: {e}")
                    break
        
        self.monitoring_thread = threading.Thread(target=monitor, daemon=True)
        self.monitoring_thread.start()
    
    def stop_monitoring(self):
        """Stop background monitoring"""
        self.monitoring_active = False
        if self.monitoring_thread:
            self.monitoring_thread.join(timeout=1.0)


class PerformanceBenchmark:
    """Individual performance benchmark"""
    
    def __init__(self, name: str, benchmark_type: BenchmarkType,
                 target_function: Callable, benchmark_args: tuple = (),
                 benchmark_kwargs: Dict[str, Any] = None):
        self.name = name
        self.benchmark_type = benchmark_type
        self.target_function = target_function
        self.benchmark_args = benchmark_args
        self.benchmark_kwargs = benchmark_kwargs or {}
        self.logger = PrismLogger(f"benchmark.{name}")
        
    def run(self, iterations: int = 1, warmup_iterations: int = 0) -> BenchmarkResult:
        """Run the benchmark"""
        self.logger.info(f"Running benchmark: {self.name} ({iterations} iterations)")
        
        # Warmup runs
        for i in range(warmup_iterations):
            try:
                self.target_function(*self.benchmark_args, **self.benchmark_kwargs)
            except Exception as e:
                self.logger.warning(f"Warmup iteration {i} failed: {e}")
        
        # Actual benchmark runs
        profiler = PerformanceProfiler()
        execution_times = []
        memory_usages = []
        
        try:
            for i in range(iterations):
                with profiler.profile_execution():
                    result = self.target_function(*self.benchmark_args, **self.benchmark_kwargs)
                
                execution_times.append(profiler.last_metrics.execution_time)
                memory_usages.append(profiler.last_metrics.peak_memory_mb)
            
            # Calculate aggregate metrics
            avg_execution_time = statistics.mean(execution_times)
            avg_memory = statistics.mean(memory_usages)
            peak_memory = max(memory_usages)
            
            # Get the last profiler metrics as base
            final_metrics = PerformanceMetrics(
                execution_time=avg_execution_time,
                peak_memory_mb=peak_memory,
                average_memory_mb=avg_memory,
                cpu_percent=profiler.last_metrics.cpu_percent,
                disk_io_mb=profiler.last_metrics.disk_io_mb,
                custom_metrics={
                    'iterations': iterations,
                    'execution_time_std': statistics.stdev(execution_times) if len(execution_times) > 1 else 0,
                    'min_execution_time': min(execution_times),
                    'max_execution_time': max(execution_times)
                }
            )
            
            return BenchmarkResult(
                benchmark_name=self.name,
                benchmark_type=self.benchmark_type,
                metrics=final_metrics,
                status="completed",
                metadata={
                    'iterations': iterations,
                    'warmup_iterations': warmup_iterations,
                    'target_function': self.target_function.__name__
                }
            )
            
        except Exception as e:
            self.logger.error(f"Benchmark failed: {e}")
            return BenchmarkResult(
                benchmark_name=self.name,
                benchmark_type=self.benchmark_type,
                metrics=PerformanceMetrics(0, 0, 0, 0, 0),
                status="failed",
                error_message=str(e)
            )


class BenchmarkSuite:
    """Collection of performance benchmarks"""
    
    def __init__(self, name: str):
        self.name = name
        self.benchmarks: List[PerformanceBenchmark] = []
        self.baseline_results: Dict[str, BenchmarkResult] = {}
        self.logger = PrismLogger(f"benchmark_suite.{name}")
        
    def add_benchmark(self, benchmark: PerformanceBenchmark):
        """Add a benchmark to the suite"""
        self.benchmarks.append(benchmark)
        self.logger.info(f"Added benchmark: {benchmark.name}")
    
    def add_execution_time_benchmark(self, name: str, target_function: Callable,
                                   args: tuple = (), kwargs: Dict[str, Any] = None):
        """Add execution time benchmark"""
        benchmark = PerformanceBenchmark(
            name=name,
            benchmark_type=BenchmarkType.EXECUTION_TIME,
            target_function=target_function,
            benchmark_args=args,
            benchmark_kwargs=kwargs or {}
        )
        self.add_benchmark(benchmark)
    
    def add_memory_benchmark(self, name: str, target_function: Callable,
                           args: tuple = (), kwargs: Dict[str, Any] = None):
        """Add memory usage benchmark"""
        benchmark = PerformanceBenchmark(
            name=name,
            benchmark_type=BenchmarkType.MEMORY_USAGE,
            target_function=target_function,
            benchmark_args=args,
            benchmark_kwargs=kwargs or {}
        )
        self.add_benchmark(benchmark)
    
    def add_throughput_benchmark(self, name: str, target_function: Callable,
                               args: tuple = (), kwargs: Dict[str, Any] = None):
        """Add throughput benchmark"""
        benchmark = PerformanceBenchmark(
            name=name,
            benchmark_type=BenchmarkType.THROUGHPUT,
            target_function=target_function,
            benchmark_args=args,
            benchmark_kwargs=kwargs or {}
        )
        self.add_benchmark(benchmark)
    
    def load_baseline(self, baseline_file: Path):
        """Load baseline results from file"""
        try:
            with open(baseline_file, 'r') as f:
                baseline_data = json.load(f)
            
            for benchmark_name, result_data in baseline_data.items():
                result_data['benchmark_type'] = BenchmarkType(result_data['benchmark_type'])
                # Reconstruct metrics
                metrics_data = result_data['metrics']
                metrics = PerformanceMetrics(**metrics_data)
                
                # Reconstruct result
                self.baseline_results[benchmark_name] = BenchmarkResult(
                    benchmark_name=benchmark_name,
                    benchmark_type=result_data['benchmark_type'],
                    metrics=metrics,
                    status=result_data['status'],
                    metadata=result_data.get('metadata', {})
                )
            
            self.logger.info(f"Loaded {len(self.baseline_results)} baseline results")
            
        except Exception as e:
            self.logger.error(f"Failed to load baseline: {e}")
    
    def save_baseline(self, baseline_file: Path, results: List[BenchmarkResult]):
        """Save results as baseline"""
        baseline_data = {}
        for result in results:
            baseline_data[result.benchmark_name] = result.to_dict()
        
        baseline_file.parent.mkdir(parents=True, exist_ok=True)
        with open(baseline_file, 'w') as f:
            json.dump(baseline_data, f, indent=2)
        
        self.logger.info(f"Saved baseline with {len(results)} results to {baseline_file}")
    
    def run_all(self, iterations: int = 3, warmup_iterations: int = 1) -> List[BenchmarkResult]:
        """Run all benchmarks in the suite"""
        self.logger.info(f"Running benchmark suite: {self.name} ({len(self.benchmarks)} benchmarks)")
        
        results = []
        for benchmark in self.benchmarks:
            result = benchmark.run(iterations, warmup_iterations)
            
            # Add baseline comparison if available
            if benchmark.name in self.baseline_results:
                result.baseline_metrics = self.baseline_results[benchmark.name].metrics
                # Recalculate performance ratio
                if result.baseline_metrics.execution_time > 0:
                    result.performance_ratio = result.baseline_metrics.execution_time / result.metrics.execution_time
            
            results.append(result)
        
        return results


class BenchmarkRunner:
    """High-level benchmark runner with reporting"""
    
    def __init__(self, output_dir: Optional[Path] = None):
        self.output_dir = output_dir or Path("benchmark_results")
        self.output_dir.mkdir(exist_ok=True)
        self.logger = PrismLogger("benchmark_runner")
        
    def run_suite(self, suite: BenchmarkSuite, iterations: int = 3,
                 save_baseline: bool = False) -> List[BenchmarkResult]:
        """Run benchmark suite and generate reports"""
        self.logger.info(f"Running benchmark suite: {suite.name}")
        
        results = suite.run_all(iterations)
        
        # Generate reports
        self._generate_reports(suite, results)
        
        # Save as baseline if requested
        if save_baseline:
            baseline_file = self.output_dir / f"{suite.name}_baseline.json"
            suite.save_baseline(baseline_file, results)
        
        return results
    
    def compare_with_baseline(self, suite: BenchmarkSuite, baseline_file: Path,
                            regression_threshold: float = 0.1) -> Dict[str, Any]:
        """Compare current results with baseline"""
        # Load baseline
        suite.load_baseline(baseline_file)
        
        # Run current benchmarks
        results = suite.run_all()
        
        # Analyze regressions
        regressions = []
        improvements = []
        
        for result in results:
            if result.baseline_metrics:
                if result.is_regression(regression_threshold):
                    regressions.append(result)
                elif result.performance_ratio and result.performance_ratio > 1.1:  # 10% improvement
                    improvements.append(result)
        
        comparison_report = {
            'total_benchmarks': len(results),
            'regressions': len(regressions),
            'improvements': len(improvements),
            'regression_details': [r.to_dict() for r in regressions],
            'improvement_details': [i.to_dict() for i in improvements],
            'threshold': regression_threshold
        }
        
        # Save comparison report
        report_file = self.output_dir / f"{suite.name}_comparison.json"
        with open(report_file, 'w') as f:
            json.dump(comparison_report, f, indent=2)
        
        return comparison_report
    
    def _generate_reports(self, suite: BenchmarkSuite, results: List[BenchmarkResult]):
        """Generate benchmark reports"""
        # JSON report
        json_report = {
            'suite_name': suite.name,
            'total_benchmarks': len(results),
            'timestamp': time.time(),
            'results': [result.to_dict() for result in results],
            'summary': self._calculate_summary(results)
        }
        
        json_file = self.output_dir / f"{suite.name}_results.json"
        with open(json_file, 'w') as f:
            json.dump(json_report, f, indent=2)
        
        # HTML report
        html_report = self._generate_html_report(suite, results)
        html_file = self.output_dir / f"{suite.name}_results.html"
        with open(html_file, 'w') as f:
            f.write(html_report)
        
        self.logger.info(f"Generated benchmark reports: {json_file}, {html_file}")
    
    def _calculate_summary(self, results: List[BenchmarkResult]) -> Dict[str, Any]:
        """Calculate summary statistics"""
        if not results:
            return {}
        
        successful_results = [r for r in results if r.status == "completed"]
        
        if not successful_results:
            return {'message': 'No successful benchmark results'}
        
        execution_times = [r.metrics.execution_time for r in successful_results]
        memory_usages = [r.metrics.peak_memory_mb for r in successful_results]
        
        return {
            'successful_benchmarks': len(successful_results),
            'failed_benchmarks': len(results) - len(successful_results),
            'total_execution_time': sum(execution_times),
            'average_execution_time': statistics.mean(execution_times),
            'max_execution_time': max(execution_times),
            'min_execution_time': min(execution_times),
            'average_memory_usage': statistics.mean(memory_usages),
            'max_memory_usage': max(memory_usages),
            'regressions_detected': len([r for r in successful_results if r.is_regression()])
        }
    
    def _generate_html_report(self, suite: BenchmarkSuite, results: List[BenchmarkResult]) -> str:
        """Generate HTML benchmark report"""
        summary = self._calculate_summary(results)
        
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>PRISM Benchmark Report - {suite.name}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 10px; border-radius: 5px; }}
                .summary {{ margin: 20px 0; }}
                .benchmark {{ margin: 10px 0; padding: 10px; border-left: 4px solid #ccc; }}
                .success {{ border-color: #4CAF50; background-color: #f1f8e9; }}
                .regression {{ border-color: #F44336; background-color: #ffebee; }}
                .improvement {{ border-color: #2196F3; background-color: #e3f2fd; }}
                .metrics {{ font-family: monospace; background-color: #f5f5f5; padding: 5px; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>PRISM Benchmark Report</h1>
                <h2>{suite.name}</h2>
            </div>
            
            <div class="summary">
                <h3>Summary</h3>
                <p><strong>Total Benchmarks:</strong> {summary.get('successful_benchmarks', 0) + summary.get('failed_benchmarks', 0)}</p>
                <p><strong>Successful:</strong> {summary.get('successful_benchmarks', 0)}</p>
                <p><strong>Failed:</strong> {summary.get('failed_benchmarks', 0)}</p>
                <p><strong>Total Execution Time:</strong> {summary.get('total_execution_time', 0):.3f}s</p>
                <p><strong>Average Execution Time:</strong> {summary.get('average_execution_time', 0):.3f}s</p>
                <p><strong>Peak Memory Usage:</strong> {summary.get('max_memory_usage', 0):.1f} MB</p>
            </div>
            
            <div class="results">
                <h3>Benchmark Results</h3>
                <table>
                    <tr>
                        <th>Benchmark</th>
                        <th>Type</th>
                        <th>Status</th>
                        <th>Execution Time (s)</th>
                        <th>Peak Memory (MB)</th>
                        <th>Performance Ratio</th>
                    </tr>
        """
        
        for result in results:
            status_class = "success"
            if result.is_regression():
                status_class = "regression"
            elif result.performance_ratio and result.performance_ratio > 1.1:
                status_class = "improvement"
            
            ratio_text = f"{result.performance_ratio:.2f}" if result.performance_ratio else "N/A"
            
            html += f"""
                    <tr class="{status_class}">
                        <td>{result.benchmark_name}</td>
                        <td>{result.benchmark_type.value}</td>
                        <td>{result.status}</td>
                        <td>{result.metrics.execution_time:.3f}</td>
                        <td>{result.metrics.peak_memory_mb:.1f}</td>
                        <td>{ratio_text}</td>
                    </tr>
            """
        
        html += """
                </table>
            </div>
        </body>
        </html>
        """
        
        return html


# Convenience functions
def create_benchmark_suite(name: str) -> BenchmarkSuite:
    """Create a new benchmark suite"""
    return BenchmarkSuite(name)


def create_execution_time_benchmark(name: str, target_function: Callable,
                                  args: tuple = (), kwargs: Dict[str, Any] = None) -> PerformanceBenchmark:
    """Create execution time benchmark"""
    return PerformanceBenchmark(
        name=name,
        benchmark_type=BenchmarkType.EXECUTION_TIME,
        target_function=target_function,
        benchmark_args=args,
        benchmark_kwargs=kwargs or {}
    )