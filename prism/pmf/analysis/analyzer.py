#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF Analysis module - WHAM analysis and visualization
"""

import os
import subprocess
from pathlib import Path
import numpy as np
from typing import Dict, List, Optional, Tuple
import logging

from ..utils.exceptions import (
    PMFAnalysisError, WHAMExecutionError, WHAMConvergenceError,
    RequiredFileNotFoundError, UmbrellaIncompleteError, PMFErrorCode
)
from ..utils.error_handling import (
    with_retry, error_context, validate_prerequisites,
    ErrorCollector, handle_external_tool_error
)
from ..utils.plotting import PMFPlotter

logger = logging.getLogger(__name__)


class PMFAnalyzer:
    """PMF analysis and WHAM runner"""
    
    def __init__(self, pmf_system):
        self.pmf_system = pmf_system
        self.analysis_dir = pmf_system.analysis_dir
        self.config = pmf_system.config['analysis']
    
    def run_wham(self, begin_time=None, bootstrap=None, energy_unit='kCal', **kwargs):
        """
        Run WHAM analysis on umbrella sampling results
        
        Parameters:
        -----------
        begin_time : int, optional
            Start time for analysis (ps)
        bootstrap : int, optional
            Number of bootstrap iterations
        energy_unit : str
            Energy unit for output
        
        Returns:
        --------
        Dict : WHAM analysis results
        
        Raises:
        -------
        PMFAnalysisError
            If analysis setup or execution fails
        WHAMExecutionError
            If WHAM command execution fails
        UmbrellaIncompleteError
            If umbrella sampling data is insufficient
        """
        with error_context("WHAM analysis", {
            "begin_time": begin_time, 
            "bootstrap": bootstrap, 
            "energy_unit": energy_unit
        }):
            # Validate prerequisites
            self._validate_wham_prerequisites()
            
            # Create analysis directory safely
            try:
                self.analysis_dir.mkdir(parents=True, exist_ok=True)
            except Exception as exc:
                raise PMFAnalysisError(
                    message=f"Failed to create analysis directory: {exc}",
                    error_code=PMFErrorCode.ANALYSIS_DATA_INSUFFICIENT,
                    recoverable=True,
                    context={"analysis_dir": str(self.analysis_dir)}
                ) from exc
            
            # Get parameters with validation
            if begin_time is None:
                begin_time = self.config.get('begin_time_ps', 2000)
            if bootstrap is None:
                bootstrap = self.config.get('bootstrap', 50)
            
            # Validate parameters
            self._validate_wham_parameters(begin_time, bootstrap, energy_unit)
            
            # Prepare file lists with validation
            wham_files = self._prepare_wham_files()
            
            # Run WHAM with error handling
            return self._run_wham_command(wham_files, begin_time, bootstrap, energy_unit)
    
    def visualize(self, pmf_file=None, error_file=None, **kwargs):
        """Generate PMF visualization plots using unified plotting utilities"""
        if pmf_file is None:
            pmf_file = self.analysis_dir / "pmf.xvg"
        if error_file is None:
            error_file = self.analysis_dir / "pmferror.xvg"

        # Use unified PMFPlotter
        plotter = PMFPlotter(self.analysis_dir)
        plots = plotter.plot_all_pmf(pmf_file, error_file)

        return plots
    
    def _prepare_wham_files(self):
        """Prepare TPR and pullf file lists for WHAM with validation"""
        umbrella_dir = self.pmf_system.umbrella_dir
        
        if not umbrella_dir.exists():
            raise RequiredFileNotFoundError(
                file_path=str(umbrella_dir),
                step="WHAM analysis preparation",
                alternative_files=[]
            )
        
        tpr_files = []
        pullf_files = []
        failed_windows = []
        
        # Collect window directories
        window_dirs = list(umbrella_dir.glob("window_*"))
        if not window_dirs:
            raise UmbrellaIncompleteError(
                completed_windows=0,
                total_windows=0,
                failed_windows=["No window directories found"]
            )
        
        # Validate each umbrella window
        error_collector = ErrorCollector()
        
        for window_dir in sorted(window_dirs):
            window_name = window_dir.name
            tpr_file = window_dir / "results" / "umbrella.tpr"
            pullf_file = window_dir / "results" / "umbrella_pullf.xvg"
            
            try:
                # Check file existence
                if not tpr_file.exists():
                    error_collector.add_error(window_name, f"Missing TPR file: {tpr_file}")
                    continue
                
                if not pullf_file.exists():
                    error_collector.add_error(window_name, f"Missing pullf file: {pullf_file}")
                    continue
                
                # Check file integrity
                if tpr_file.stat().st_size == 0:
                    error_collector.add_error(window_name, f"Empty TPR file: {tpr_file}")
                    continue
                
                if pullf_file.stat().st_size < 100:  # Minimum reasonable size
                    error_collector.add_error(window_name, f"Suspiciously small pullf file: {pullf_file}")
                    continue
                
                # Files are valid
                tpr_files.append(tpr_file)
                pullf_files.append(pullf_file)
                
            except Exception as exc:
                error_collector.add_error(window_name, f"Validation error: {exc}")
        
        # Check if we have sufficient data
        total_windows = len(window_dirs)
        completed_windows = len(tpr_files)
        completion_rate = completed_windows / total_windows if total_windows > 0 else 0
        
        if completed_windows == 0:
            raise UmbrellaIncompleteError(
                completed_windows=0,
                total_windows=total_windows,
                failed_windows=[error for error in error_collector.errors.keys()]
            )
        
        if completion_rate < 0.8:  # Require at least 80% completion
            logger.warning(f"Only {completion_rate:.1%} of umbrella windows completed")
            if error_collector.has_errors():
                logger.warning("Window validation errors:")
                logger.warning(error_collector.get_error_summary())
        
        if completed_windows < 5:  # Minimum windows for meaningful PMF
            raise UmbrellaIncompleteError(
                completed_windows=completed_windows,
                total_windows=total_windows,
                failed_windows=[f"Insufficient windows: need at least 5, have {completed_windows}"]
            )
        
        logger.info(f"Found {completed_windows}/{total_windows} completed umbrella windows")
        
        # Write file lists
        tpr_list = self.analysis_dir / "tpr-files.dat"
        pullf_list = self.analysis_dir / "pullf-files.dat"
        
        with open(tpr_list, 'w') as f:
            for tpr_file in tpr_files:
                f.write(f"{tpr_file.absolute()}\n")
        
        with open(pullf_list, 'w') as f:
            for pullf_file in pullf_files:
                f.write(f"{pullf_file.absolute()}\n")
        
        return {
            'tpr_list': tpr_list,
            'pullf_list': pullf_list,
            'n_windows': len(tpr_files)
        }
    
    @with_retry(max_attempts=2, retry_on=[WHAMExecutionError])
    def _run_wham_command(self, wham_files, begin_time, bootstrap, energy_unit):
        """Execute GROMACS WHAM command with enhanced error handling"""
        
        wham_cmd = [
            "gmx", "wham",
            "-it", str(wham_files['tpr_list']),
            "-if", str(wham_files['pullf_list']),
            "-b", str(begin_time),
            "-o", "pmf.xvg",
            "-unit", energy_unit,
            "-bsres", "pmferror.xvg",
            "-bsprof", "bsprofile.xvg",
            "-nBootstrap", str(bootstrap),
            "-bs-method", "b-hist",
            "-v"
        ]
        
        logger.info(f"Running WHAM analysis with {wham_files['n_windows']} windows")
        logger.debug(f"WHAM command: {' '.join(wham_cmd)}")
        
        try:
            with error_context("WHAM execution", {
                "command": " ".join(wham_cmd),
                "n_windows": wham_files['n_windows'],
                "begin_time": begin_time,
                "bootstrap": bootstrap
            }):
                result = subprocess.run(
                    wham_cmd,
                    cwd=self.analysis_dir,
                    capture_output=True,
                    text=True,
                    check=False,  # We'll handle return codes manually
                    timeout=3600
                )
                
                # Check for errors
                if result.returncode != 0:
                    handle_external_tool_error(
                        command=" ".join(wham_cmd),
                        exit_code=result.returncode,
                        stderr=result.stderr,
                        tool_name="GROMACS WHAM"
                    )
                
                # Check for convergence issues in output
                self._check_wham_convergence(result.stdout, result.stderr)
                
                logger.info("WHAM analysis completed successfully")
            
            # Calculate binding energy
            pmf_file = self.analysis_dir / "pmf.xvg"
            binding_energy = self._calculate_binding_energy(pmf_file)
            
            return {
                'pmf_file': str(pmf_file),
                'error_file': str(self.analysis_dir / "pmferror.xvg"),
                'profile_file': str(self.analysis_dir / "bsprofile.xvg"),
                'binding_energy': binding_energy,
                'n_windows': wham_files['n_windows'],
                'begin_time': begin_time,
                'bootstrap': bootstrap,
                'energy_unit': energy_unit
            }
            
        except subprocess.CalledProcessError as e:
            logger.error(f"WHAM analysis failed: {e.stderr}")
            raise WHAMExecutionError(
                command=" ".join(wham_cmd),
                exit_code=e.returncode,
                stderr_output=e.stderr
            )
        except subprocess.TimeoutExpired:
            logger.error("WHAM analysis timed out")
            raise PMFAnalysisError(
                message="WHAM analysis timed out after 3600 seconds",
                error_code=PMFErrorCode.TIMEOUT_EXCEEDED,
                recoverable=True,
                recovery_suggestions=[
                    "Increase timeout limit",
                    "Check system performance",
                    "Consider reducing data size"
                ]
            )
    
    def _validate_wham_prerequisites(self):
        """Validate prerequisites for WHAM analysis"""
        errors = []
        
        # Check umbrella directory exists
        umbrella_dir = self.pmf_system.umbrella_dir
        if not umbrella_dir.exists():
            raise RequiredFileNotFoundError(
                file_path=str(umbrella_dir),
                step="WHAM analysis",
                alternative_files=[]
            )
        
        # Check for completed umbrella windows
        window_dirs = list(umbrella_dir.glob("window_*"))
        if not window_dirs:
            raise UmbrellaIncompleteError(
                completed_windows=0,
                total_windows=0,
                failed_windows=["No umbrella windows found"]
            )
        
        # Check each window has required files
        valid_windows = 0
        for window_dir in window_dirs:
            tpr_file = window_dir / "results" / "umbrella.tpr"
            pullf_file = window_dir / "results" / "umbrella_pullf.xvg"
            
            if tpr_file.exists() and pullf_file.exists():
                valid_windows += 1
        
        if valid_windows == 0:
            raise UmbrellaIncompleteError(
                completed_windows=0,
                total_windows=len(window_dirs),
                failed_windows=["No valid umbrella windows with complete results"]
            )
        
        logger.info(f"WHAM prerequisites validated: {valid_windows}/{len(window_dirs)} windows ready")
    
    def _validate_wham_parameters(self, begin_time, bootstrap, energy_unit):
        """Validate WHAM analysis parameters"""
        errors = []
        
        # Validate begin_time
        if begin_time < 0:
            errors.append("begin_time must be non-negative")
        elif begin_time > 50000:  # Reasonable upper limit
            errors.append("begin_time appears too large (>50ns)")
        
        # Validate bootstrap
        if bootstrap < 1:
            errors.append("bootstrap iterations must be positive")
        elif bootstrap > 1000:  # Reasonable upper limit
            errors.append("bootstrap iterations appears too large (>1000)")
        
        # Validate energy unit
        valid_units = ['kJ', 'kcal', 'kCal', 'kT']
        if energy_unit not in valid_units:
            errors.append(f"energy_unit must be one of: {valid_units}")
        
        if errors:
            raise PMFAnalysisError(
                message=f"Invalid WHAM parameters: {'; '.join(errors)}",
                error_code=PMFErrorCode.CONFIG_INVALID,
                recoverable=True,
                recovery_suggestions=[
                    "Check parameter values and ranges",
                    "Refer to GROMACS WHAM documentation"
                ],
                context={
                    "begin_time": begin_time,
                    "bootstrap": bootstrap,
                    "energy_unit": energy_unit,
                    "validation_errors": errors
                }
            )
    
    def _check_wham_convergence(self, stdout, stderr):
        """Check WHAM output for convergence issues"""
        # Check for convergence warnings in output
        convergence_issues = []
        
        if stdout:
            stdout_lower = stdout.lower()
            if "not converged" in stdout_lower:
                convergence_issues.append("WHAM did not converge")
            if "warning" in stdout_lower and "iteration" in stdout_lower:
                convergence_issues.append("Convergence warnings detected")
        
        if stderr:
            stderr_lower = stderr.lower()
            if "error" in stderr_lower:
                convergence_issues.append("Errors detected in WHAM output")
            if "failed" in stderr_lower:
                convergence_issues.append("WHAM execution failures detected")
        
        # If we detected issues, raise convergence error
        if convergence_issues:
            raise WHAMConvergenceError(
                iterations=0,  # We don't parse exact iteration count
                tolerance=0.0,
                final_error=0.0
            )
        
        logger.debug("WHAM convergence check passed")
    
    def _calculate_binding_energy(self, pmf_file):
        """Calculate binding energy from PMF data with error handling"""
        try:
            with error_context("binding energy calculation", {"pmf_file": str(pmf_file)}):
                # Check file exists and is readable
                if not pmf_file.exists():
                    raise RequiredFileNotFoundError(
                        file_path=str(pmf_file),
                        step="binding energy calculation",
                        alternative_files=[]
                    )
                
                # Load PMF data
                try:
                    data = np.loadtxt(pmf_file, comments=['#', '@'])
                except Exception as exc:
                    raise PMFAnalysisError(
                        message=f"Failed to read PMF file: {exc}",
                        error_code=PMFErrorCode.FILE_CORRUPTED,
                        recoverable=False,
                        context={"pmf_file": str(pmf_file)}
                    ) from exc
                
                # Validate data structure
                if data.ndim != 2 or data.shape[1] < 2:
                    raise PMFAnalysisError(
                        message=f"Invalid PMF file format: expected 2 columns, got {data.shape}",
                        error_code=PMFErrorCode.ANALYSIS_DATA_INSUFFICIENT,
                        recoverable=False,
                        context={"data_shape": data.shape}
                    )
                
                distances = data[:, 0]
                pmf_values = data[:, 1]
                
                # Check for reasonable data
                if len(pmf_values) < 2:
                    raise PMFAnalysisError(
                        message="Insufficient PMF data points for binding energy calculation",
                        error_code=PMFErrorCode.ANALYSIS_DATA_INSUFFICIENT,
                        recoverable=False,
                        context={"data_points": len(pmf_values)}
                    )
                
                # Calculate binding energy
                min_pmf = np.min(pmf_values)
                max_pmf = np.max(pmf_values)
                binding_energy = max_pmf - min_pmf
                
                min_distance = distances[np.argmin(pmf_values)]
                max_distance = distances[np.argmax(pmf_values)]
                
                logger.info(f"Binding energy calculated: {binding_energy:.2f} kcal/mol")
                
                return {
                    'value': binding_energy,
                    'min_pmf': min_pmf,
                    'max_pmf': max_pmf,
                    'min_distance': min_distance,
                    'max_distance': max_distance
                }
                
        except PMFAnalysisError:
            raise  # Re-raise PMF errors as-is
        except Exception as exc:
            logger.warning(f"Could not calculate binding energy: {exc}")
            return None


class WhamRunner:
    """Standalone WHAM runner for advanced analysis"""
    
    def __init__(self, tpr_files, pullf_files, output_dir):
        self.tpr_files = tpr_files
        self.pullf_files = pullf_files
        self.output_dir = Path(output_dir)
    
    def run(self, **kwargs):
        """Run WHAM with custom parameters"""
        # Implementation for advanced WHAM analysis
        pass
    
    def bootstrap_analysis(self, n_iterations=100):
        """Run bootstrap error analysis"""
        # Implementation for bootstrap analysis
        pass