#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Code Quality Assurance and Static Analysis

Comprehensive code quality checking, static analysis, and coverage analysis
tools for maintaining high-quality PRISM codebase.
"""

import ast
import os
import sys
import subprocess
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Set, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
import importlib.util
import re

from ..utils.logging_system import PrismLogger


class QualityLevel(Enum):
    """Code quality severity levels"""
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


class IssueCategory(Enum):
    """Categories of code quality issues"""
    STYLE = "style"
    COMPLEXITY = "complexity"
    MAINTAINABILITY = "maintainability"
    PERFORMANCE = "performance"
    SECURITY = "security"
    DOCUMENTATION = "documentation"
    TESTING = "testing"
    IMPORTS = "imports"
    NAMING = "naming"
    DESIGN = "design"


@dataclass
class QualityIssue:
    """Individual code quality issue"""
    file_path: str
    line_number: int
    column_number: int
    issue_type: str
    category: IssueCategory
    severity: QualityLevel
    message: str
    rule_id: Optional[str] = None
    suggestion: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'file_path': self.file_path,
            'line_number': self.line_number,
            'column_number': self.column_number,
            'issue_type': self.issue_type,
            'category': self.category.value,
            'severity': self.severity.value,
            'message': self.message,
            'rule_id': self.rule_id,
            'suggestion': self.suggestion
        }


@dataclass
class QualityMetrics:
    """Code quality metrics"""
    total_lines: int
    code_lines: int
    comment_lines: int
    blank_lines: int
    functions_count: int
    classes_count: int
    complexity_score: float
    maintainability_index: float
    documentation_coverage: float
    test_coverage: float
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class QualityReport:
    """Comprehensive quality report"""
    total_issues: int
    issues_by_severity: Dict[str, int]
    issues_by_category: Dict[str, int]
    quality_score: float
    metrics: QualityMetrics
    issues: List[QualityIssue]
    files_analyzed: int
    analysis_time: float
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'total_issues': self.total_issues,
            'issues_by_severity': self.issues_by_severity,
            'issues_by_category': self.issues_by_category,
            'quality_score': self.quality_score,
            'metrics': self.metrics.to_dict(),
            'issues': [issue.to_dict() for issue in self.issues],
            'files_analyzed': self.files_analyzed,
            'analysis_time': self.analysis_time
        }


class ASTAnalyzer:
    """AST-based code analysis"""
    
    def __init__(self):
        self.logger = PrismLogger("ast_analyzer")
        
    def analyze_file(self, file_path: Path) -> List[QualityIssue]:
        """Analyze a Python file using AST"""
        issues = []
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            tree = ast.parse(content, filename=str(file_path))
            
            # Run various AST checks
            issues.extend(self._check_complexity(tree, file_path))
            issues.extend(self._check_naming_conventions(tree, file_path))
            issues.extend(self._check_function_length(tree, file_path))
            issues.extend(self._check_imports(tree, file_path))
            issues.extend(self._check_docstrings(tree, file_path))
            
        except SyntaxError as e:
            issues.append(QualityIssue(
                file_path=str(file_path),
                line_number=e.lineno or 0,
                column_number=e.offset or 0,
                issue_type="syntax_error",
                category=IssueCategory.STYLE,
                severity=QualityLevel.ERROR,
                message=f"Syntax error: {e.msg}",
                rule_id="E999"
            ))
        except Exception as e:
            self.logger.warning(f"Error analyzing {file_path}: {e}")
        
        return issues
    
    def _check_complexity(self, tree: ast.AST, file_path: Path) -> List[QualityIssue]:
        """Check cyclomatic complexity"""
        issues = []
        
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                complexity = self._calculate_complexity(node)
                
                if complexity > 10:
                    severity = QualityLevel.ERROR if complexity > 15 else QualityLevel.WARNING
                    issues.append(QualityIssue(
                        file_path=str(file_path),
                        line_number=node.lineno,
                        column_number=node.col_offset,
                        issue_type="high_complexity",
                        category=IssueCategory.COMPLEXITY,
                        severity=severity,
                        message=f"Function '{node.name}' has high complexity ({complexity})",
                        rule_id="C901",
                        suggestion="Consider breaking this function into smaller functions"
                    ))
        
        return issues
    
    def _calculate_complexity(self, node: ast.AST) -> int:
        """Calculate cyclomatic complexity"""
        complexity = 1  # Base complexity
        
        for child in ast.walk(node):
            if isinstance(child, (ast.If, ast.While, ast.For, ast.AsyncFor)):
                complexity += 1
            elif isinstance(child, ast.ExceptHandler):
                complexity += 1
            elif isinstance(child, (ast.And, ast.Or)):
                complexity += 1
        
        return complexity
    
    def _check_naming_conventions(self, tree: ast.AST, file_path: Path) -> List[QualityIssue]:
        """Check naming conventions"""
        issues = []
        
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                if not self._is_snake_case(node.name) and not node.name.startswith('_'):
                    issues.append(QualityIssue(
                        file_path=str(file_path),
                        line_number=node.lineno,
                        column_number=node.col_offset,
                        issue_type="naming_convention",
                        category=IssueCategory.NAMING,
                        severity=QualityLevel.WARNING,
                        message=f"Function '{node.name}' should use snake_case naming",
                        rule_id="N802",
                        suggestion="Use snake_case for function names"
                    ))
            
            elif isinstance(node, ast.ClassDef):
                if not self._is_pascal_case(node.name):
                    issues.append(QualityIssue(
                        file_path=str(file_path),
                        line_number=node.lineno,
                        column_number=node.col_offset,
                        issue_type="naming_convention",
                        category=IssueCategory.NAMING,
                        severity=QualityLevel.WARNING,
                        message=f"Class '{node.name}' should use PascalCase naming",
                        rule_id="N801",
                        suggestion="Use PascalCase for class names"
                    ))
        
        return issues
    
    def _check_function_length(self, tree: ast.AST, file_path: Path) -> List[QualityIssue]:
        """Check function length"""
        issues = []
        
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                # Calculate function length
                if hasattr(node, 'end_lineno') and node.end_lineno:
                    length = node.end_lineno - node.lineno
                    
                    if length > 50:
                        severity = QualityLevel.ERROR if length > 100 else QualityLevel.WARNING
                        issues.append(QualityIssue(
                            file_path=str(file_path),
                            line_number=node.lineno,
                            column_number=node.col_offset,
                            issue_type="long_function",
                            category=IssueCategory.MAINTAINABILITY,
                            severity=severity,
                            message=f"Function '{node.name}' is too long ({length} lines)",
                            rule_id="R701",
                            suggestion="Consider breaking this function into smaller functions"
                        ))
        
        return issues
    
    def _check_imports(self, tree: ast.AST, file_path: Path) -> List[QualityIssue]:
        """Check import style and organization"""
        issues = []
        imports = []
        
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                imports.append((node.lineno, node))
        
        # Check for imports not at top of file (ignoring docstrings)
        first_non_import_line = None
        for line_no, node in imports:
            if first_non_import_line is None:
                # Find first non-import statement
                for check_node in ast.walk(tree):
                    if (hasattr(check_node, 'lineno') and 
                        check_node.lineno > line_no and
                        not isinstance(check_node, (ast.Import, ast.ImportFrom, ast.Expr))):
                        first_non_import_line = check_node.lineno
                        break
            
            if first_non_import_line and line_no > first_non_import_line:
                issues.append(QualityIssue(
                    file_path=str(file_path),
                    line_number=line_no,
                    column_number=0,
                    issue_type="import_order",
                    category=IssueCategory.IMPORTS,
                    severity=QualityLevel.WARNING,
                    message="Import not at top of file",
                    rule_id="E402",
                    suggestion="Move imports to top of file"
                ))
        
        return issues
    
    def _check_docstrings(self, tree: ast.AST, file_path: Path) -> List[QualityIssue]:
        """Check for missing docstrings"""
        issues = []
        
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                docstring = ast.get_docstring(node)
                
                if not docstring:
                    node_type = "Function" if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) else "Class"
                    
                    # Skip private methods and special methods for now
                    if not node.name.startswith('_'):
                        issues.append(QualityIssue(
                            file_path=str(file_path),
                            line_number=node.lineno,
                            column_number=node.col_offset,
                            issue_type="missing_docstring",
                            category=IssueCategory.DOCUMENTATION,
                            severity=QualityLevel.INFO,
                            message=f"{node_type} '{node.name}' is missing a docstring",
                            rule_id="D100",
                            suggestion=f"Add a docstring to {node_type.lower()} '{node.name}'"
                        ))
        
        return issues
    
    def _is_snake_case(self, name: str) -> bool:
        """Check if name follows snake_case convention"""
        return re.match(r'^[a-z][a-z0-9_]*$', name) is not None
    
    def _is_pascal_case(self, name: str) -> bool:
        """Check if name follows PascalCase convention"""
        return re.match(r'^[A-Z][a-zA-Z0-9]*$', name) is not None


class StaticAnalyzer:
    """Static code analysis using external tools"""
    
    def __init__(self):
        self.logger = PrismLogger("static_analyzer")
        
    def run_flake8(self, source_dir: Path) -> List[QualityIssue]:
        """Run flake8 static analysis"""
        issues = []
        
        try:
            cmd = ["flake8", "--format=json", str(source_dir)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.stdout:
                flake8_results = json.loads(result.stdout)
                
                for item in flake8_results:
                    issues.append(QualityIssue(
                        file_path=item['filename'],
                        line_number=item['line_number'],
                        column_number=item['column_number'],
                        issue_type="flake8_violation",
                        category=self._categorize_flake8_code(item['code']),
                        severity=self._severity_from_flake8_code(item['code']),
                        message=item['text'],
                        rule_id=item['code']
                    ))
                    
        except subprocess.CalledProcessError:
            self.logger.warning("flake8 not available or failed")
        except json.JSONDecodeError:
            self.logger.warning("Could not parse flake8 output")
        except Exception as e:
            self.logger.warning(f"Error running flake8: {e}")
        
        return issues
    
    def run_pylint(self, source_dir: Path) -> List[QualityIssue]:
        """Run pylint static analysis"""
        issues = []
        
        try:
            cmd = ["pylint", "--output-format=json", str(source_dir)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.stdout:
                pylint_results = json.loads(result.stdout)
                
                for item in pylint_results:
                    issues.append(QualityIssue(
                        file_path=item['path'],
                        line_number=item['line'],
                        column_number=item['column'],
                        issue_type=item['type'],
                        category=self._categorize_pylint_type(item['type']),
                        severity=self._severity_from_pylint_type(item['type']),
                        message=item['message'],
                        rule_id=item['message-id']
                    ))
                    
        except subprocess.CalledProcessError:
            self.logger.warning("pylint not available or failed")
        except json.JSONDecodeError:
            self.logger.warning("Could not parse pylint output")
        except Exception as e:
            self.logger.warning(f"Error running pylint: {e}")
        
        return issues
    
    def run_mypy(self, source_dir: Path) -> List[QualityIssue]:
        """Run mypy type checking"""
        issues = []
        
        try:
            cmd = ["mypy", "--show-error-codes", "--json-report", "/tmp/mypy_report", str(source_dir)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            # Parse mypy output from stderr (mypy writes to stderr)
            for line in result.stderr.split('\n'):
                if ':' in line and 'error:' in line:
                    parts = line.split(':')
                    if len(parts) >= 4:
                        file_path = parts[0]
                        line_number = int(parts[1]) if parts[1].isdigit() else 0
                        column_number = int(parts[2]) if parts[2].isdigit() else 0
                        message = ':'.join(parts[3:]).strip()
                        
                        issues.append(QualityIssue(
                            file_path=file_path,
                            line_number=line_number,
                            column_number=column_number,
                            issue_type="type_error",
                            category=IssueCategory.STYLE,
                            severity=QualityLevel.ERROR,
                            message=message,
                            rule_id="mypy"
                        ))
                        
        except subprocess.CalledProcessError:
            self.logger.warning("mypy not available or failed")
        except Exception as e:
            self.logger.warning(f"Error running mypy: {e}")
        
        return issues
    
    def _categorize_flake8_code(self, code: str) -> IssueCategory:
        """Categorize flake8 error codes"""
        if code.startswith('E'):
            return IssueCategory.STYLE
        elif code.startswith('W'):
            return IssueCategory.STYLE
        elif code.startswith('F'):
            return IssueCategory.DESIGN
        elif code.startswith('C'):
            return IssueCategory.COMPLEXITY
        elif code.startswith('N'):
            return IssueCategory.NAMING
        else:
            return IssueCategory.STYLE
    
    def _severity_from_flake8_code(self, code: str) -> QualityLevel:
        """Determine severity from flake8 code"""
        if code.startswith('E'):
            return QualityLevel.ERROR
        elif code.startswith('W'):
            return QualityLevel.WARNING
        elif code.startswith('F'):
            return QualityLevel.ERROR
        else:
            return QualityLevel.WARNING
    
    def _categorize_pylint_type(self, msg_type: str) -> IssueCategory:
        """Categorize pylint message types"""
        mapping = {
            'error': IssueCategory.DESIGN,
            'warning': IssueCategory.MAINTAINABILITY,
            'refactor': IssueCategory.MAINTAINABILITY,
            'convention': IssueCategory.STYLE,
            'info': IssueCategory.STYLE
        }
        return mapping.get(msg_type, IssueCategory.STYLE)
    
    def _severity_from_pylint_type(self, msg_type: str) -> QualityLevel:
        """Determine severity from pylint type"""
        mapping = {
            'error': QualityLevel.ERROR,
            'warning': QualityLevel.WARNING,
            'refactor': QualityLevel.INFO,
            'convention': QualityLevel.INFO,
            'info': QualityLevel.INFO
        }
        return mapping.get(msg_type, QualityLevel.INFO)


class CoverageAnalyzer:
    """Test coverage analysis"""
    
    def __init__(self):
        self.logger = PrismLogger("coverage_analyzer")
    
    def analyze_coverage(self, source_dir: Path, test_dir: Path) -> Dict[str, Any]:
        """Analyze test coverage"""
        coverage_data = {
            'line_coverage': 0.0,
            'branch_coverage': 0.0,
            'function_coverage': 0.0,
            'files_covered': 0,
            'total_files': 0,
            'uncovered_files': [],
            'detailed_coverage': {}
        }
        
        try:
            # Try to run coverage analysis
            cmd = ["coverage", "run", "--source", str(source_dir), "-m", "pytest", str(test_dir)]
            subprocess.run(cmd, capture_output=True, text=True)
            
            # Get coverage report
            cmd = ["coverage", "report", "--format=json"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.stdout:
                coverage_json = json.loads(result.stdout)
                coverage_data.update({
                    'line_coverage': coverage_json['totals']['percent_covered'],
                    'files_covered': coverage_json['totals']['num_statements'],
                    'total_files': len(coverage_json['files'])
                })
                
                # Detailed file coverage
                for filename, file_data in coverage_json['files'].items():
                    coverage_data['detailed_coverage'][filename] = {
                        'lines_covered': file_data['summary']['num_statements'],
                        'lines_missing': file_data['summary']['missing_lines'],
                        'percent_covered': file_data['summary']['percent_covered']
                    }
                    
                    if file_data['summary']['percent_covered'] < 50:
                        coverage_data['uncovered_files'].append(filename)
            
        except subprocess.CalledProcessError:
            self.logger.warning("coverage tool not available")
        except json.JSONDecodeError:
            self.logger.warning("Could not parse coverage output")
        except Exception as e:
            self.logger.warning(f"Error running coverage analysis: {e}")
        
        return coverage_data


class CodeQualityChecker:
    """Main code quality checker"""
    
    def __init__(self, source_dir: Path):
        self.source_dir = source_dir
        self.logger = PrismLogger("code_quality_checker")
        self.ast_analyzer = ASTAnalyzer()
        self.static_analyzer = StaticAnalyzer()
        self.coverage_analyzer = CoverageAnalyzer()
        
    def run_full_analysis(self, test_dir: Optional[Path] = None) -> QualityReport:
        """Run comprehensive quality analysis"""
        import time
        start_time = time.time()
        
        self.logger.info(f"Starting quality analysis for {self.source_dir}")
        
        # Collect all Python files
        python_files = list(self.source_dir.rglob("*.py"))
        python_files = [f for f in python_files if not f.name.startswith('.')]
        
        all_issues = []
        
        # AST analysis
        self.logger.info("Running AST analysis...")
        for file_path in python_files:
            issues = self.ast_analyzer.analyze_file(file_path)
            all_issues.extend(issues)
        
        # Static analysis
        self.logger.info("Running static analysis...")
        all_issues.extend(self.static_analyzer.run_flake8(self.source_dir))
        
        # Calculate metrics
        metrics = self._calculate_metrics(python_files)
        
        # Coverage analysis
        if test_dir:
            self.logger.info("Running coverage analysis...")
            coverage_data = self.coverage_analyzer.analyze_coverage(self.source_dir, test_dir)
            metrics.test_coverage = coverage_data['line_coverage']
        
        # Calculate quality score
        quality_score = self._calculate_quality_score(all_issues, metrics)
        
        # Categorize issues
        issues_by_severity = {}
        issues_by_category = {}
        
        for issue in all_issues:
            severity = issue.severity.value
            category = issue.category.value
            
            issues_by_severity[severity] = issues_by_severity.get(severity, 0) + 1
            issues_by_category[category] = issues_by_category.get(category, 0) + 1
        
        analysis_time = time.time() - start_time
        
        report = QualityReport(
            total_issues=len(all_issues),
            issues_by_severity=issues_by_severity,
            issues_by_category=issues_by_category,
            quality_score=quality_score,
            metrics=metrics,
            issues=all_issues,
            files_analyzed=len(python_files),
            analysis_time=analysis_time
        )
        
        self.logger.info(f"Quality analysis completed: {len(all_issues)} issues found in {analysis_time:.2f}s")
        
        return report
    
    def _calculate_metrics(self, python_files: List[Path]) -> QualityMetrics:
        """Calculate code metrics"""
        total_lines = 0
        code_lines = 0
        comment_lines = 0
        blank_lines = 0
        functions_count = 0
        classes_count = 0
        complexity_scores = []
        
        for file_path in python_files:
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                
                total_lines += len(lines)
                
                # Count line types
                for line in lines:
                    stripped = line.strip()
                    if not stripped:
                        blank_lines += 1
                    elif stripped.startswith('#'):
                        comment_lines += 1
                    else:
                        code_lines += 1
                
                # AST analysis for functions and classes
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                    
                    tree = ast.parse(content)
                    
                    for node in ast.walk(tree):
                        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                            functions_count += 1
                            complexity = self.ast_analyzer._calculate_complexity(node)
                            complexity_scores.append(complexity)
                        elif isinstance(node, ast.ClassDef):
                            classes_count += 1
                
                except:
                    pass  # Skip files with syntax errors
                    
            except Exception as e:
                self.logger.warning(f"Error analyzing {file_path}: {e}")
        
        # Calculate derived metrics
        avg_complexity = sum(complexity_scores) / len(complexity_scores) if complexity_scores else 0
        comment_ratio = comment_lines / total_lines if total_lines > 0 else 0
        
        # Simple maintainability index calculation
        maintainability_index = max(0, 171 - 5.2 * avg_complexity - 0.23 * (functions_count or 1) + 16.2 * comment_ratio)
        
        return QualityMetrics(
            total_lines=total_lines,
            code_lines=code_lines,
            comment_lines=comment_lines,
            blank_lines=blank_lines,
            functions_count=functions_count,
            classes_count=classes_count,
            complexity_score=avg_complexity,
            maintainability_index=maintainability_index,
            documentation_coverage=comment_ratio * 100,
            test_coverage=0.0  # Will be updated by coverage analysis
        )
    
    def _calculate_quality_score(self, issues: List[QualityIssue], metrics: QualityMetrics) -> float:
        """Calculate overall quality score (0-100)"""
        base_score = 100.0
        
        # Deduct points for issues
        for issue in issues:
            if issue.severity == QualityLevel.CRITICAL:
                base_score -= 10
            elif issue.severity == QualityLevel.ERROR:
                base_score -= 5
            elif issue.severity == QualityLevel.WARNING:
                base_score -= 2
            elif issue.severity == QualityLevel.INFO:
                base_score -= 0.5
        
        # Adjust for complexity
        if metrics.complexity_score > 10:
            base_score -= (metrics.complexity_score - 10) * 2
        
        # Adjust for maintainability
        if metrics.maintainability_index < 50:
            base_score -= (50 - metrics.maintainability_index) * 0.5
        
        # Bonus for good documentation
        if metrics.documentation_coverage > 20:
            base_score += min(5, (metrics.documentation_coverage - 20) * 0.2)
        
        return max(0, min(100, base_score))
    
    def generate_html_report(self, report: QualityReport, output_file: Path):
        """Generate HTML quality report"""
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>PRISM Code Quality Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 10px; border-radius: 5px; }}
                .score {{ font-size: 2em; font-weight: bold; color: {'green' if report.quality_score > 80 else 'orange' if report.quality_score > 60 else 'red'}; }}
                .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }}
                .metric {{ background-color: #f9f9f9; padding: 10px; border-radius: 5px; }}
                .issues {{ margin: 20px 0; }}
                .issue {{ margin: 5px 0; padding: 5px; border-left: 4px solid #ccc; }}
                .critical {{ border-color: #d32f2f; background-color: #ffebee; }}
                .error {{ border-color: #f57c00; background-color: #fff3e0; }}
                .warning {{ border-color: #fbc02d; background-color: #fffde7; }}
                .info {{ border-color: #1976d2; background-color: #e3f2fd; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>PRISM Code Quality Report</h1>
                <div class="score">Quality Score: {report.quality_score:.1f}/100</div>
            </div>
            
            <div class="metrics">
                <div class="metric">
                    <h3>Overview</h3>
                    <p>Files Analyzed: {report.files_analyzed}</p>
                    <p>Total Issues: {report.total_issues}</p>
                    <p>Analysis Time: {report.analysis_time:.2f}s</p>
                </div>
                
                <div class="metric">
                    <h3>Code Metrics</h3>
                    <p>Total Lines: {report.metrics.total_lines:,}</p>
                    <p>Code Lines: {report.metrics.code_lines:,}</p>
                    <p>Functions: {report.metrics.functions_count}</p>
                    <p>Classes: {report.metrics.classes_count}</p>
                </div>
                
                <div class="metric">
                    <h3>Quality Metrics</h3>
                    <p>Complexity: {report.metrics.complexity_score:.1f}</p>
                    <p>Maintainability: {report.metrics.maintainability_index:.1f}</p>
                    <p>Documentation: {report.metrics.documentation_coverage:.1f}%</p>
                    <p>Test Coverage: {report.metrics.test_coverage:.1f}%</p>
                </div>
            </div>
            
            <div class="issues">
                <h3>Issues by Severity</h3>
        """
        
        for severity, count in report.issues_by_severity.items():
            html += f"<p>{severity.title()}: {count}</p>"
        
        html += """
                <h3>Recent Issues</h3>
        """
        
        # Show first 20 issues
        for issue in report.issues[:20]:
            severity_class = issue.severity.value
            html += f"""
                <div class="issue {severity_class}">
                    <strong>{issue.file_path}:{issue.line_number}</strong> - {issue.message}
                    <br><small>{issue.category.value} | {issue.rule_id or 'N/A'}</small>
                </div>
            """
        
        html += """
            </div>
        </body>
        </html>
        """
        
        with open(output_file, 'w') as f:
            f.write(html)


# Convenience functions
def analyze_code_quality(source_dir: Path, test_dir: Optional[Path] = None,
                        output_dir: Optional[Path] = None) -> QualityReport:
    """Analyze code quality for a directory"""
    checker = CodeQualityChecker(source_dir)
    report = checker.run_full_analysis(test_dir)
    
    if output_dir:
        output_dir.mkdir(exist_ok=True)
        
        # Save JSON report
        json_file = output_dir / "quality_report.json"
        with open(json_file, 'w') as f:
            json.dump(report.to_dict(), f, indent=2)
        
        # Save HTML report
        html_file = output_dir / "quality_report.html"
        checker.generate_html_report(report, html_file)
    
    return report