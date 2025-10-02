#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Log Analysis and Reporting Tools

Provides comprehensive log analysis, pattern detection, and report generation
for PRISM molecular dynamics simulations with performance analysis and
troubleshooting assistance.
"""

import re
import json
import statistics
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, Any, List, Optional, Tuple, Union, Callable
from dataclasses import dataclass, asdict
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.figure import Figure


@dataclass
class LogPattern:
    """Log pattern for analysis"""
    pattern: str
    description: str
    severity: str
    category: str
    regex: Optional[re.Pattern] = None
    
    def __post_init__(self):
        self.regex = re.compile(self.pattern, re.IGNORECASE)


@dataclass
class AnalysisResult:
    """Result of log analysis"""
    total_entries: int
    time_range: Tuple[datetime, datetime]
    error_count: int
    warning_count: int
    performance_metrics: Dict[str, Any]
    detected_patterns: List[Dict[str, Any]]
    recommendations: List[str]
    
    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data['time_range'] = [self.time_range[0].isoformat(), self.time_range[1].isoformat()]
        return data


class PrismLogAnalyzer:
    """Comprehensive log analyzer for PRISM simulations"""
    
    def __init__(self):
        self.patterns = self._load_analysis_patterns()
        self.performance_extractors = self._setup_performance_extractors()
    
    def _load_analysis_patterns(self) -> List[LogPattern]:
        """Load predefined analysis patterns"""
        return [
            # Error patterns
            LogPattern(
                r"RuntimeError|Exception|Error:|Failed|failed",
                "General error or exception",
                "error",
                "errors"
            ),
            LogPattern(
                r"GROMACS.*failed|gmx.*failed|mdrun.*failed",
                "GROMACS execution failure",
                "error",
                "gromacs"
            ),
            LogPattern(
                r"File not found|No such file|Permission denied",
                "File system related error",
                "error",
                "filesystem"
            ),
            LogPattern(
                r"Memory|memory.*error|out of memory|allocation failed",
                "Memory related issue",
                "error",
                "memory"
            ),
            
            # Warning patterns
            LogPattern(
                r"Warning|warning|WARN",
                "General warning",
                "warning",
                "warnings"
            ),
            LogPattern(
                r"pull rate.*high|pulling.*fast|rate.*artifact",
                "SMD pull rate warning",
                "warning",
                "scientific"
            ),
            LogPattern(
                r"temperature.*inconsist|temp.*differ",
                "Temperature inconsistency",
                "warning",
                "scientific"
            ),
            LogPattern(
                r"convergence|converge.*fail|not converged",
                "Convergence issues",
                "warning",
                "scientific"
            ),
            
            # Performance patterns
            LogPattern(
                r"(\d+\.?\d*)\s*ns/day",
                "Performance metric: ns/day",
                "info",
                "performance"
            ),
            LogPattern(
                r"(\d+\.?\d*)\s*hour/ns",
                "Performance metric: hour/ns",
                "info",
                "performance"
            ),
            LogPattern(
                r"Step\s+(\d+)",
                "Simulation step progress",
                "info",
                "progress"
            ),
            
            # Success patterns
            LogPattern(
                r"completed successfully|finished|done|success",
                "Successful completion",
                "info",
                "success"
            ),
            
            # Resource usage patterns
            LogPattern(
                r"CPU.*(\d+\.?\d*)%|cpu.*usage.*(\d+\.?\d*)",
                "CPU usage information",
                "info",
                "resources"
            ),
            LogPattern(
                r"Memory.*(\d+\.?\d*)\s*(MB|GB)|memory.*usage.*(\d+\.?\d*)",
                "Memory usage information",
                "info",
                "resources"
            )
        ]
    
    def _setup_performance_extractors(self) -> Dict[str, Callable]:
        """Setup performance data extractors"""
        return {
            'ns_per_day': lambda line: self._extract_numeric(line, r'(\d+\.?\d*)\s*ns/day'),
            'hours_per_ns': lambda line: self._extract_numeric(line, r'(\d+\.?\d*)\s*hour/ns'),
            'cpu_percent': lambda line: self._extract_numeric(line, r'CPU.*(\d+\.?\d*)%'),
            'memory_mb': lambda line: self._extract_numeric(line, r'Memory.*(\d+\.?\d*)\s*MB'),
            'step_number': lambda line: self._extract_numeric(line, r'Step\s+(\d+)', int)
        }
    
    def _extract_numeric(self, text: str, pattern: str, dtype=float) -> Optional[float]:
        """Extract numeric value from text using regex"""
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            try:
                return dtype(match.group(1))
            except (ValueError, IndexError):
                pass
        return None
    
    def analyze_log_file(self, log_file: Union[str, Path]) -> AnalysisResult:
        """Analyze a single log file"""
        log_file = Path(log_file)
        
        if not log_file.exists():
            raise FileNotFoundError(f"Log file not found: {log_file}")
        
        entries = []
        detected_patterns = []
        performance_data = defaultdict(list)
        
        # Read and parse log file
        with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line_num, line in enumerate(f, 1):
                entry = self._parse_log_entry(line.strip(), line_num)
                if entry:
                    entries.append(entry)
                
                # Check for patterns
                for pattern in self.patterns:
                    if pattern.regex.search(line):
                        detected_patterns.append({
                            'line_number': line_num,
                            'pattern': pattern.description,
                            'severity': pattern.severity,
                            'category': pattern.category,
                            'text': line.strip()[:200]  # Truncate long lines
                        })
                
                # Extract performance data
                for metric, extractor in self.performance_extractors.items():
                    value = extractor(line)
                    if value is not None:
                        performance_data[metric].append({
                            'line_number': line_num,
                            'value': value,
                            'timestamp': entry.get('timestamp') if entry else None
                        })
        
        # Calculate analysis results
        if not entries:
            return AnalysisResult(
                total_entries=0,
                time_range=(datetime.now(), datetime.now()),
                error_count=0,
                warning_count=0,
                performance_metrics={},
                detected_patterns=[],
                recommendations=["No valid log entries found"]
            )
        
        # Extract timestamps
        timestamps = [entry['timestamp'] for entry in entries if entry.get('timestamp')]
        time_range = (min(timestamps), max(timestamps)) if timestamps else (datetime.now(), datetime.now())
        
        # Count issues
        error_count = len([p for p in detected_patterns if p['severity'] == 'error'])
        warning_count = len([p for p in detected_patterns if p['severity'] == 'warning'])
        
        # Analyze performance metrics
        performance_metrics = self._analyze_performance_data(performance_data)
        
        # Generate recommendations
        recommendations = self._generate_recommendations(detected_patterns, performance_metrics)
        
        return AnalysisResult(
            total_entries=len(entries),
            time_range=time_range,
            error_count=error_count,
            warning_count=warning_count,
            performance_metrics=performance_metrics,
            detected_patterns=detected_patterns,
            recommendations=recommendations
        )
    
    def _parse_log_entry(self, line: str, line_number: int) -> Optional[Dict[str, Any]]:
        """Parse a single log entry"""
        if not line.strip():
            return None
        
        entry = {
            'line_number': line_number,
            'text': line,
            'timestamp': None,
            'level': None,
            'message': line
        }
        
        # Try to parse structured JSON log
        if line.startswith('{'):
            try:
                json_data = json.loads(line)
                entry.update(json_data)
                if 'timestamp' in json_data:
                    entry['timestamp'] = datetime.fromisoformat(json_data['timestamp'].replace('Z', '+00:00'))
                return entry
            except json.JSONDecodeError:
                pass
        
        # Try to parse standard log format
        timestamp_match = re.match(r'(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})', line)
        if timestamp_match:
            try:
                entry['timestamp'] = datetime.strptime(timestamp_match.group(1), '%Y-%m-%d %H:%M:%S')
            except ValueError:
                pass
        
        # Extract log level
        level_match = re.search(r'(DEBUG|INFO|WARNING|ERROR|CRITICAL)', line)
        if level_match:
            entry['level'] = level_match.group(1)
        
        return entry
    
    def _analyze_performance_data(self, performance_data: Dict[str, List]) -> Dict[str, Any]:
        """Analyze performance data"""
        metrics = {}
        
        for metric, data_points in performance_data.items():
            if not data_points:
                continue
            
            values = [point['value'] for point in data_points]
            
            metrics[metric] = {
                'count': len(values),
                'min': min(values),
                'max': max(values),
                'mean': statistics.mean(values),
                'median': statistics.median(values),
                'std_dev': statistics.stdev(values) if len(values) > 1 else 0
            }
            
            # Add performance assessment
            if metric == 'ns_per_day':
                avg_ns_day = metrics[metric]['mean']
                if avg_ns_day > 10:
                    metrics[metric]['assessment'] = 'Excellent performance'
                elif avg_ns_day > 5:
                    metrics[metric]['assessment'] = 'Good performance'
                elif avg_ns_day > 1:
                    metrics[metric]['assessment'] = 'Moderate performance'
                else:
                    metrics[metric]['assessment'] = 'Poor performance'
            
            elif metric == 'hours_per_ns':
                avg_hours_ns = metrics[metric]['mean']
                if avg_hours_ns < 0.1:
                    metrics[metric]['assessment'] = 'Excellent performance'
                elif avg_hours_ns < 0.5:
                    metrics[metric]['assessment'] = 'Good performance'
                elif avg_hours_ns < 2:
                    metrics[metric]['assessment'] = 'Moderate performance'
                else:
                    metrics[metric]['assessment'] = 'Poor performance'
        
        return metrics
    
    def _generate_recommendations(self, patterns: List[Dict], performance: Dict[str, Any]) -> List[str]:
        """Generate recommendations based on analysis"""
        recommendations = []
        
        # Error-based recommendations
        error_patterns = [p for p in patterns if p['severity'] == 'error']
        if error_patterns:
            error_categories = Counter(p['category'] for p in error_patterns)
            
            if 'gromacs' in error_categories:
                recommendations.append("Check GROMACS installation and ensure proper environment setup")
            
            if 'filesystem' in error_categories:
                recommendations.append("Verify file paths and permissions for input/output files")
            
            if 'memory' in error_categories:
                recommendations.append("Consider reducing system size or increasing available memory")
        
        # Warning-based recommendations
        warning_patterns = [p for p in patterns if p['severity'] == 'warning']
        warning_categories = Counter(p['category'] for p in warning_patterns)
        
        if 'scientific' in warning_categories:
            recommendations.append("Review simulation parameters for scientific accuracy")
        
        # Performance-based recommendations
        if 'ns_per_day' in performance:
            ns_day = performance['ns_per_day'].get('mean', 0)
            if ns_day < 1:
                recommendations.append("Consider optimizing simulation performance: use GPU acceleration, adjust parallelization")
        
        if 'hours_per_ns' in performance:
            hours_ns = performance['hours_per_ns'].get('mean', 0)
            if hours_ns > 2:
                recommendations.append("Simulation performance is slow: check hardware utilization and simulation settings")
        
        # General recommendations
        if not recommendations:
            recommendations.append("No major issues detected in the log analysis")
        
        return recommendations
    
    def analyze_multiple_logs(self, log_files: List[Union[str, Path]]) -> Dict[str, AnalysisResult]:
        """Analyze multiple log files"""
        results = {}
        
        for log_file in log_files:
            log_file = Path(log_file)
            try:
                results[str(log_file)] = self.analyze_log_file(log_file)
            except Exception as e:
                results[str(log_file)] = f"Analysis failed: {e}"
        
        return results
    
    def generate_summary_report(self, analysis_results: Union[AnalysisResult, Dict[str, AnalysisResult]]) -> str:
        """Generate human-readable summary report"""
        if isinstance(analysis_results, AnalysisResult):
            return self._generate_single_report(analysis_results)
        else:
            return self._generate_multi_report(analysis_results)
    
    def _generate_single_report(self, result: AnalysisResult) -> str:
        """Generate report for single analysis result"""
        lines = []
        lines.append("=" * 70)
        lines.append("PRISM LOG ANALYSIS REPORT")
        lines.append("=" * 70)
        
        # Summary statistics
        lines.append(f"Total log entries: {result.total_entries:,}")
        lines.append(f"Time range: {result.time_range[0].strftime('%Y-%m-%d %H:%M:%S')} to {result.time_range[1].strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Duration: {(result.time_range[1] - result.time_range[0]).total_seconds() / 3600:.2f} hours")
        lines.append(f"Errors detected: {result.error_count}")
        lines.append(f"Warnings detected: {result.warning_count}")
        
        # Performance metrics
        if result.performance_metrics:
            lines.append("\nPERFORMANCE METRICS:")
            lines.append("-" * 40)
            
            for metric, data in result.performance_metrics.items():
                metric_name = metric.replace('_', ' ').title()
                lines.append(f"{metric_name}:")
                lines.append(f"  Average: {data['mean']:.2f}")
                lines.append(f"  Range: {data['min']:.2f} - {data['max']:.2f}")
                if 'assessment' in data:
                    lines.append(f"  Assessment: {data['assessment']}")
                lines.append("")
        
        # Pattern analysis
        if result.detected_patterns:
            pattern_categories = Counter(p['category'] for p in result.detected_patterns)
            lines.append("DETECTED PATTERNS:")
            lines.append("-" * 40)
            
            for category, count in pattern_categories.most_common():
                lines.append(f"{category.title()}: {count} occurrences")
            
            # Show most frequent issues
            lines.append("\nMOST FREQUENT ISSUES:")
            pattern_descriptions = Counter(p['pattern'] for p in result.detected_patterns 
                                         if p['severity'] in ['error', 'warning'])
            for description, count in pattern_descriptions.most_common(5):
                lines.append(f"  • {description} ({count} times)")
        
        # Recommendations
        if result.recommendations:
            lines.append("\nRECOMMENDATIONS:")
            lines.append("-" * 40)
            for i, rec in enumerate(result.recommendations, 1):
                lines.append(f"{i}. {rec}")
        
        lines.append("\n" + "=" * 70)
        return "\n".join(lines)
    
    def _generate_multi_report(self, results: Dict[str, AnalysisResult]) -> str:
        """Generate report for multiple analysis results"""
        lines = []
        lines.append("=" * 70)
        lines.append("PRISM MULTI-LOG ANALYSIS REPORT")
        lines.append("=" * 70)
        
        total_entries = sum(r.total_entries for r in results.values() if isinstance(r, AnalysisResult))
        total_errors = sum(r.error_count for r in results.values() if isinstance(r, AnalysisResult))
        total_warnings = sum(r.warning_count for r in results.values() if isinstance(r, AnalysisResult))
        
        lines.append(f"Files analyzed: {len(results)}")
        lines.append(f"Total log entries: {total_entries:,}")
        lines.append(f"Total errors: {total_errors}")
        lines.append(f"Total warnings: {total_warnings}")
        lines.append("")
        
        # Individual file summaries
        lines.append("INDIVIDUAL FILE SUMMARIES:")
        lines.append("-" * 50)
        
        for log_file, result in results.items():
            if isinstance(result, AnalysisResult):
                lines.append(f"\n{Path(log_file).name}:")
                lines.append(f"  Entries: {result.total_entries:,}")
                lines.append(f"  Errors: {result.error_count}")
                lines.append(f"  Warnings: {result.warning_count}")
            else:
                lines.append(f"\n{Path(log_file).name}: {result}")
        
        lines.append("\n" + "=" * 70)
        return "\n".join(lines)
    
    def create_performance_plot(self, performance_data: Dict[str, List], 
                              output_file: Optional[Path] = None) -> Optional[Figure]:
        """Create performance visualization"""
        if not performance_data:
            return None
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle('PRISM Simulation Performance Analysis', fontsize=16)
        
        plot_configs = [
            ('ns_per_day', 'Performance (ns/day)', axes[0, 0]),
            ('hours_per_ns', 'Performance (hours/ns)', axes[0, 1]),
            ('cpu_percent', 'CPU Usage (%)', axes[1, 0]),
            ('memory_mb', 'Memory Usage (MB)', axes[1, 1])
        ]
        
        for metric, title, ax in plot_configs:
            if metric in performance_data:
                data = performance_data[metric]
                values = [point['value'] for point in data]
                
                if values:
                    ax.plot(range(len(values)), values, marker='o', markersize=3)
                    ax.set_title(title)
                    ax.set_xlabel('Data Point')
                    ax.set_ylabel('Value')
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(title)
            else:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(title)
        
        plt.tight_layout()
        
        if output_file:
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig


# Convenience functions
def analyze_prism_log(log_file: Union[str, Path]) -> AnalysisResult:
    """Analyze a PRISM log file"""
    analyzer = PrismLogAnalyzer()
    return analyzer.analyze_log_file(log_file)


def analyze_simulation_logs(log_directory: Union[str, Path], pattern: str = "*.log") -> Dict[str, AnalysisResult]:
    """Analyze all log files in a directory"""
    log_dir = Path(log_directory)
    log_files = list(log_dir.glob(pattern))
    
    analyzer = PrismLogAnalyzer()
    return analyzer.analyze_multiple_logs(log_files)


def generate_log_report(log_file: Union[str, Path], output_file: Optional[Path] = None) -> str:
    """Generate log analysis report"""
    result = analyze_prism_log(log_file)
    analyzer = PrismLogAnalyzer()
    report = analyzer.generate_summary_report(result)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
    
    return report