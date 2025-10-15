#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Workflow Manager

This module manages the complete PMF calculation workflow, integrating SMD, 
umbrella sampling, and WHAM analysis following PRISM architectural patterns.
"""

import os
import shutil
import subprocess
from pathlib import Path
import yaml
import logging
from typing import Dict, List, Optional, Union

from ..methods.smd import SMDManager
from ..methods.umbrella import UmbrellaManager
from ..analysis.analyzer import PMFAnalyzer
from ..analysis.validator import MDResultsValidator, ValidationLevel, validate_md_results
from ..utils.exceptions import (
    PMFWorkflowError, PrerequisiteNotMetError, SystemNotFoundException,
    WorkflowStateError, PMFAnalysisError, PMFErrorCode
)
from ..utils.error_handling import (
    with_retry, error_context, validate_prerequisites, ErrorCollector
)
from ..utils.plotting import PMFPlotter

logger = logging.getLogger(__name__)


class PMFWorkflow:
    """
    PRISM PMF Workflow Manager
    
    Manages the complete PMF calculation workflow with support for both
    automated and step-by-step execution modes, following PRISM patterns.
    Enhanced with PMF-specific box remodeling capabilities.
    """
    
    def __init__(self, system_dir: Union[str, Path], output_dir: Union[str, Path], 
                 config: Dict):
        """
        Initialize PMF workflow manager
        
        Parameters:
        -----------
        system_dir : str or Path
            Path to system directory containing GMX_PMF_SYSTEM
        output_dir : str or Path
            Output directory for PMF calculations
        config : Dict
            PMF configuration dictionary
        """
        self.system_dir = Path(system_dir)
        self.output_dir = Path(output_dir)
        
        # Inherit force field information from original system
        self.config = self._inherit_system_config(config)
        
        # Locate MD results directory
        self.md_results_dir = self._locate_md_results()
        
        # Create output structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up working directories
        self.smd_dir = self.output_dir / "smd"
        self.umbrella_dir = self.output_dir / "umbrella"
        self.analysis_dir = self.output_dir / "analysis"
        
        # Initialize managers
        self.smd_manager = SMDManager(self)
        self.umbrella_manager = UmbrellaManager(self)
        self.analyzer = PMFAnalyzer(self)
        
        # Load existing workflow state if available
        self._load_workflow_state()
        
        # Track workflow state (enhanced with PMF system building)
        self.state = {
            'pmf_system_built': False,
            'pmf_system_equilibrated': False,
            'smd_prepared': False,
            'smd_completed': False,
            'umbrella_prepared': False, 
            'umbrella_completed': False,
            'analysis_completed': False
        }
        
        # PMF system builder (created when needed)
        self.pmf_builder = None
        self.rebuilt_system_dir = self.output_dir / "GMX_PMF_SYSTEM"
        
        logger.info(f"PMF workflow initialized")
        logger.info(f"  System: {self.md_results_dir}")
        logger.info(f"  Output: {self.output_dir}")
        logger.info(f"  Enhanced: PMF box remodeling enabled")
    
    def _locate_md_results(self) -> Path:
        """Locate MD results directory using enhanced validator"""
        logger.debug(f"Locating MD results in: {self.system_dir}")
        
        # Use the enhanced MD validator to locate and validate results
        validator = MDResultsValidator(ValidationLevel.BASIC)
        validation_result = validator.validate(self.system_dir)
        
        if validation_result.is_valid:
            # Extract the MD directory from validation result
            md_dir = self._extract_md_dir_from_validation(validation_result)
            logger.info(f"MD results validation: {validation_result.summary}")
            return md_dir
        else:
            # Enhanced error message with validation details
            error_details = []
            if validation_result.errors:
                error_details.extend(validation_result.errors)
            if validation_result.files_missing:
                error_details.append(f"Missing files: {', '.join(validation_result.files_missing[:3])}")
            
            error_msg = f"Invalid MD results directory: {self.system_dir}"
            if error_details:
                error_msg += f"\nIssues found: {'; '.join(error_details[:2])}"
            
            raise FileNotFoundError(error_msg)
    
    def _extract_md_dir_from_validation(self, validation_result) -> Path:
        """Extract MD directory path from validation result"""
        # Get the parent directory of the first found file
        if validation_result.files_found:
            first_file = validation_result.files_found[0]
            # Find the directory containing the topology file
            for file_info in validation_result.files_found:
                if "topol.top" in str(file_info.path):
                    return file_info.path.parent
            
            # Fallback to first file's parent
            return first_file.path.parent
        
        # Last resort - try the old logic
        return self._locate_md_results_fallback()
    
    def _locate_md_results_fallback(self) -> Path:
        """Locate PMF system directory - only looks for GMX_PMF_SYSTEM"""
        # Check if system_dir itself contains GMX files (direct PMF system)
        if (self.system_dir / "solv_ions.gro").exists() and (self.system_dir / "topol.top").exists():
            return self.system_dir
        
        # Look for PMF rebuilt system (GMX_PMF_SYSTEM)
        pmf_system_dir = self.system_dir / "GMX_PMF_SYSTEM"
        if pmf_system_dir.exists() and (pmf_system_dir / "solv_ions.gro").exists():
            return pmf_system_dir
            
        # Search in subdirectories for PMF system only
        for subdir in self.system_dir.iterdir():
            if subdir.is_dir():
                pmf_candidate = subdir / "GMX_PMF_SYSTEM"
                if pmf_candidate.exists() and (pmf_candidate / "solv_ions.gro").exists():
                    return pmf_candidate
        
        raise FileNotFoundError(
            f"Could not locate PMF system (GMX_PMF_SYSTEM) in {self.system_dir}. "
            "Please run PMF system remodeling first using pmf_remodel_only.py"
        )
    
    def build_pmf_system(self, equilibrate: bool = True, **kwargs) -> Dict:
        """
        Build PMF-optimized system with Z-axis alignment and extended box (Step 0)
        
        Parameters:
        -----------
        equilibrate : bool, optional
            Whether to run equilibration after building (default: True)
        **kwargs : optional
            Additional PMF builder parameters (pulling_distance, box_distance, etc.)
            
        Returns:
        --------
        Dict : PMF system build results
        """
        logger.info("=== Building PMF-Optimized System ===")
        logger.info("This step extracts Protein+LIG, aligns to Z-axis, and extends box for pulling")
        
        # Import PMF builder here to avoid circular imports
        from ..builders.pmf_builder import PMFBuilder
        
        # Create PMF builder if not exists
        if not self.pmf_builder:
            self.pmf_builder = PMFBuilder(
                md_results_dir=self.md_results_dir,
                output_dir=self.output_dir,
                config=self.config,
                **kwargs
            )
        
        # Build PMF-optimized system
        build_results = self.pmf_builder.build(equilibrate=equilibrate)
        
        # Update state
        self.state['pmf_system_built'] = True
        if equilibrate:
            self.state['pmf_system_equilibrated'] = True
        
        # Update system references
        self.rebuilt_system_dir = Path(build_results['system_dir'])
        
        # Save state
        self._save_workflow_state()
        
        results = {
            **build_results,
            'workflow_step': 'pmf_system_build',
            'next_step': 'prepare_smd() - Use PMF-optimized system for SMD'
        }
        
        logger.info("PMF system build completed")
        logger.info(f"   PMF system: {self.rebuilt_system_dir}")
        logger.info(f"   Z-axis aligned: {build_results.get('alignment', {}).get('z_axis_aligned', True)}")
        logger.info(f"   Equilibrated: {'Yes' if equilibrate else 'No'}")
        
        return results
    
    def prepare_smd(self, use_pmf_system: bool = True) -> Dict:
        """
        Prepare SMD simulation (Step 1)
        
        Parameters:
        -----------
        use_pmf_system : bool, optional
            Whether to use PMF-rebuilt system (default: True)
            
        Returns:
        --------
        Dict : SMD preparation results with execution instructions
        
        Raises:
        -------
        PrerequisiteNotMetError
            If required prerequisites are not met
        PMFWorkflowError
            If SMD preparation fails
        """
        with error_context("SMD preparation", {"use_pmf_system": use_pmf_system}):
            logger.info("=== Preparing SMD Simulation ===")
            
            # Validate workflow prerequisites
            self._validate_smd_workflow_prerequisites(use_pmf_system)
            
            # Check if PMF system should be used but doesn't exist
            if use_pmf_system and not self.state['pmf_system_built']:
                logger.warning("PMF system not built yet. Consider running build_pmf_system() first.")
                logger.warning("Proceeding with PMF system build...")
                self.build_pmf_system(equilibrate=True)
            
            try:
                result = self.smd_manager.prepare(use_pmf_rebuilt_system=use_pmf_system)
                self.state['smd_prepared'] = True
                
                # Save state
                self._save_workflow_state()
                
                logger.info("SMD preparation completed")
                logger.info(f"System type: {'PMF-rebuilt (Z-optimized)' if use_pmf_system else 'Original MD'}")
                logger.info(f"Run SMD: cd {result['smd_dir']} && bash run_smd.sh")
                
                return result
                
            except Exception as exc:
                logger.error(f"SMD preparation failed: {exc}")
                raise PMFWorkflowError(
                    message=f"SMD preparation failed: {exc}",
                    error_code=PMFErrorCode.SMD_PREPARATION_FAILED,
                    recoverable=True,
                    recovery_suggestions=[
                        "Check system validation results",
                        "Verify PMF system build completion",
                        "Check disk space and permissions"
                    ],
                    context={"use_pmf_system": use_pmf_system}
                ) from exc
    
    def prepare_umbrella(self) -> Dict:
        """
        Prepare umbrella sampling (Step 2 - requires SMD completion)
        
        Returns:
        --------
        Dict : Umbrella preparation results with execution instructions
        """
        logger.info("=== Preparing Umbrella Sampling ===")
        
        # Check SMD completion
        if not self._check_smd_completion():
            raise RuntimeError(
                "SMD simulation must be completed before umbrella preparation. "
                f"Run: cd {self.smd_dir} && bash run_smd.sh"
            )
        
        result = self.umbrella_manager.setup_windows()
        self.state['umbrella_prepared'] = True
        
        # Save state
        self._save_workflow_state()
        
        logger.info("Umbrella sampling preparation completed")
        logger.info(f"Run umbrella: cd {result['umbrella_dir']} && bash run_all_umbrella.sh parallel")
        
        return result
    
    def prepare_analysis(self) -> Dict:
        """
        Prepare WHAM analysis (Step 3 - requires umbrella completion)
        
        Returns:
        --------
        Dict : Analysis preparation results
        """
        logger.info("=== Preparing WHAM Analysis ===")
        
        # Check umbrella completion
        if not self._check_umbrella_completion():
            raise RuntimeError(
                "Umbrella sampling must be completed before analysis. "
                f"Run: cd {self.umbrella_dir} && bash run_all_umbrella.sh parallel"
            )
        
        # Analysis preparation (create directories, validate files)
        self.analysis_dir.mkdir(exist_ok=True)
        
        result = {
            'analysis_dir': str(self.analysis_dir),
            'status': 'ready_for_analysis',
            'next_action': 'Call run_analysis() to perform WHAM analysis'
        }
        
        logger.info("Analysis preparation completed")
        return result
    
    def run_analysis(self) -> Dict:
        """
        Run WHAM analysis and generate results
        
        Returns:
        --------
        Dict : Complete analysis results with binding energy and plots
        """
        logger.info("=== Running WHAM Analysis ===")
        
        # Check prerequisites
        if not self._check_umbrella_completion():
            raise RuntimeError("Umbrella sampling must be completed before analysis")
        
        # Run WHAM analysis
        wham_results = self.analyzer.run_wham()
        
        # Generate visualizations
        plots = self.analyzer.visualize()
        
        # Combine results
        results = {
            **wham_results,
            'plots': plots,
            'analysis_dir': str(self.analysis_dir)
        }
        
        self.state['analysis_completed'] = True
        self._save_workflow_state()
        
        # Generate summary report
        self._generate_summary_report(results)
        
        logger.info("PMF analysis completed successfully!")
        if 'binding_energy' in results:
            be = results['binding_energy']
            logger.info(f"Binding energy: {be['value']:.2f} {results.get('energy_unit', 'kcal/mol')}")
        
        return results
    
    def run_complete(self, build_pmf_system: bool = True, equilibrate_pmf_system: bool = True) -> Dict:
        """
        Run complete PMF workflow (automated mode)
        
        Parameters:
        -----------
        build_pmf_system : bool, optional
            Whether to build PMF-optimized system first (default: True)
        equilibrate_pmf_system : bool, optional
            Whether to equilibrate PMF system after building (default: True)
            
        Returns:
        --------
        Dict : Complete workflow results
        """
        logger.info("=== Running Complete PMF Workflow ===")
        logger.info(f"Enhanced workflow: {'With PMF system building' if build_pmf_system else 'Using existing system'}")
        
        try:
            # Step 0: Build PMF-optimized system (optional but recommended)
            if build_pmf_system and not self.state['pmf_system_built']:
                logger.info("Step 0: Building PMF-optimized system...")
                pmf_build_results = self.build_pmf_system(equilibrate=equilibrate_pmf_system)
                logger.info("PMF system ready for calculations")
            
            # Step 1: Prepare and run SMD
            if not self.state['smd_completed']:
                smd_prep = self.prepare_smd(use_pmf_system=build_pmf_system) 
                logger.info("Running SMD simulation...")
                self._run_smd_automated()
                self.state['smd_completed'] = True
            
            # Step 2: Prepare and run umbrella sampling  
            if not self.state['umbrella_completed']:
                umbrella_prep = self.prepare_umbrella()
                logger.info("Running umbrella sampling...")
                self._run_umbrella_automated()
                self.state['umbrella_completed'] = True
            
            # Step 3: Run analysis
            if not self.state['analysis_completed']:
                results = self.run_analysis()
            else:
                results = self._load_existing_results()
            
            # Add workflow summary
            results['workflow_summary'] = {
                'pmf_system_used': build_pmf_system,
                'system_equilibrated': equilibrate_pmf_system if build_pmf_system else False,
                'z_axis_optimized': build_pmf_system,
                'total_stages': 4 if build_pmf_system else 3
            }
            
            logger.info("Complete PMF workflow finished successfully!")
            if build_pmf_system:
                logger.info("   Used PMF-optimized system with Z-axis alignment")
            
            return results
            
        except Exception as e:
            logger.error(f"PMF workflow failed: {e}")
            status = self.get_status()
            logger.info(f"Workflow stopped at: {status['current_stage']}")
            raise
    
    def get_status(self) -> Dict:
        """Get comprehensive workflow status"""
        # Determine current stage (enhanced with PMF system building)
        if self.state['analysis_completed']:
            current_stage = "completed"
        elif self.state['umbrella_completed'] or self._check_umbrella_completion():
            current_stage = "ready_for_analysis"
        elif self.state['umbrella_prepared']:
            current_stage = "running_umbrella"
        elif self.state['smd_completed'] or self._check_smd_completion():
            current_stage = "ready_for_umbrella"
        elif self.state['smd_prepared']:
            current_stage = "running_smd"
        elif self.state['pmf_system_built']:
            current_stage = "pmf_system_ready"
        else:
            current_stage = "ready_to_start"
        
        # Check file status
        file_status = self._check_file_status()
        
        # Get next action (enhanced with PMF system building)
        next_actions = {
            "ready_to_start": "Call build_pmf_system() [recommended] or prepare_smd(use_pmf_system=False)",
            "pmf_system_ready": "Call prepare_smd() to use PMF-optimized system",
            "running_smd": f"Run SMD: cd {self.smd_dir} && bash run_smd.sh",
            "ready_for_umbrella": "Call prepare_umbrella()",
            "running_umbrella": f"Run umbrella: cd {self.umbrella_dir} && bash run_all_umbrella.sh parallel",
            "ready_for_analysis": "Call run_analysis()",
            "completed": "Workflow complete. Results available in analysis directory."
        }
        
        return {
            'current_stage': current_stage,
            'state': self.state.copy(),
            'files': file_status,
            'pmf_system_info': {
                'built': self.state['pmf_system_built'],
                'equilibrated': self.state.get('pmf_system_equilibrated', False),
                'path': str(self.rebuilt_system_dir) if self.state['pmf_system_built'] else None
            },
            'next_action': next_actions.get(current_stage, "Unknown")
        }
    
    def analyze_smd(self) -> Dict:
        """Analyze SMD results and generate plots"""
        if not self._check_smd_completion():
            raise RuntimeError("SMD simulation not completed")
        
        logger.info("Analyzing SMD results...")
        
        # Read SMD data
        smd_results_dir = self.smd_dir / "results"
        pullf_file = smd_results_dir / "smd_pullf.xvg"
        pullx_file = smd_results_dir / "smd_pullx.xvg"
        
        # Generate SMD analysis plots
        plots = self._generate_smd_plots(pullf_file, pullx_file)
        
        # Calculate SMD statistics
        stats = self._calculate_smd_stats(pullf_file, pullx_file)
        
        results = {
            'smd_analysis': stats,
            'plots': plots,
            'analysis_dir': str(self.smd_dir / "analysis")
        }
        
        logger.info("SMD analysis completed")
        return results
    
    def analyze_pmf(self) -> Dict:
        """Analyze final PMF results with detailed plots and statistics"""
        if not self.state['analysis_completed']:
            raise RuntimeError("WHAM analysis not completed")
        
        logger.info("Performing detailed PMF analysis...")
        
        pmf_file = self.analysis_dir / "pmf.xvg"
        if not pmf_file.exists():
            raise FileNotFoundError(f"PMF file not found: {pmf_file}")
        
        # Generate detailed analysis
        detailed_plots = self._generate_detailed_pmf_plots(pmf_file)
        pmf_statistics = self._calculate_pmf_statistics(pmf_file)
        
        results = {
            'pmf_statistics': pmf_statistics,
            'detailed_plots': detailed_plots,
            'analysis_dir': str(self.analysis_dir)
        }
        
        logger.info("Detailed PMF analysis completed")
        return results
    
    def clean(self, components: Optional[List[str]] = None) -> None:
        """Clean PMF calculation results"""
        if components is None:
            components = ['smd', 'umbrella', 'analysis']
        
        for component in components:
            if component == 'smd' and self.smd_dir.exists():
                shutil.rmtree(self.smd_dir)
                logger.info("SMD results cleaned")
                self.state['smd_prepared'] = False
                self.state['smd_completed'] = False
                
            elif component == 'umbrella' and self.umbrella_dir.exists():
                shutil.rmtree(self.umbrella_dir)  
                logger.info("Umbrella results cleaned")
                self.state['umbrella_prepared'] = False
                self.state['umbrella_completed'] = False
                
            elif component == 'analysis' and self.analysis_dir.exists():
                shutil.rmtree(self.analysis_dir)
                logger.info("Analysis results cleaned")
                self.state['analysis_completed'] = False
        
        self._save_workflow_state()
    
    def _check_smd_completion(self) -> bool:
        """Check if SMD simulation is completed"""
        # Remove smd.gro as SMD trajectories are difficult to extract final structures
        required_files = [
            self.smd_dir / "results" / "smd_pullf.xvg",
            self.smd_dir / "results" / "smd_pullx.xvg",
            self.smd_dir / "results" / "smd.xtc"  # Use trajectory file instead
        ]
        return all(f.exists() for f in required_files)
    
    def _check_umbrella_completion(self) -> bool:
        """Check if umbrella sampling is completed"""
        if not self.umbrella_dir.exists():
            return False
        
        # Count completed windows
        completed = len(list(self.umbrella_dir.glob("window_*/results/umbrella_pullf.xvg")))
        total = len(list(self.umbrella_dir.glob("window_*")))
        
        return completed >= max(5, total * 0.8)  # At least 5 windows or 80% completion
    
    def _check_file_status(self) -> Dict:
        """Check status of key workflow files using enhanced validator"""
        # Use enhanced validator for system files
        validation_result = validate_md_results(self.md_results_dir, ValidationLevel.STANDARD)
        
        # Extract system file status from validation
        system_files = {
            'validation_score': validation_result.score,
            'validation_grade': validation_result.grade,
            'is_valid': validation_result.is_valid,
            'structure_files': len([f for f in validation_result.files_found 
                                   if f.file_type.value == 'structure']),
            'topology_files': len([f for f in validation_result.files_found 
                                  if f.file_type.value == 'topology']),
            'trajectory_files': len([f for f in validation_result.files_found 
                                    if f.file_type.value == 'trajectory']),
            'warnings': len(validation_result.warnings),
            'errors': len(validation_result.errors)
        }
        
        # Check PMF system files
        pmf_system_files = {
            'built': self.state['pmf_system_built'],
            'structure_exists': (self.rebuilt_system_dir / "solv_ions.gro").exists() if self.rebuilt_system_dir.exists() else False,
            'topology_exists': (self.rebuilt_system_dir / "topol.top").exists() if self.rebuilt_system_dir.exists() else False,
            'equilibrated': self.state.get('pmf_system_equilibrated', False),
            'equilibration_structure': (self.output_dir / "equilibration" / "npt" / "npt.gro").exists(),
            'equilibration_checkpoint': (self.output_dir / "equilibration" / "npt" / "npt.cpt").exists()
        }
        
        return {
            'system_files': system_files,
            'pmf_system_files': pmf_system_files,
            'smd_files': {
                'prepared': (self.smd_dir / "smd.mdp").exists() if self.smd_dir.exists() else False,
                'completed': self._check_smd_completion()
            },
            'umbrella_files': {
                'prepared': len(list(self.umbrella_dir.glob("window_*"))) if self.umbrella_dir.exists() else 0,
                'completed': self._check_umbrella_completion()
            },
            'analysis_files': {
                'pmf_data': (self.analysis_dir / "pmf.xvg").exists() if self.analysis_dir.exists() else False,
                'plots': (self.analysis_dir / "pmf_curve.png").exists() if self.analysis_dir.exists() else False
            }
        }
    
    def _save_workflow_state(self) -> None:
        """Save workflow state to file"""
        state_file = self.output_dir / "workflow_state.yaml"
        with open(state_file, 'w') as f:
            yaml.dump(self.state, f)
    
    def _validate_smd_workflow_prerequisites(self, use_pmf_system: bool) -> None:
        """Validate prerequisites for SMD workflow step"""
        
        requirements = {
            'system_directory_exists': lambda: self.system_dir.exists(),
            'output_directory_writable': lambda: self._check_directory_writable(self.output_dir)
        }
        
        if use_pmf_system:
            # Additional checks for PMF system workflow
            requirements.update({
                'md_system_validated': lambda: self._check_md_system_validation()
            })
        
        validate_prerequisites(requirements, 'SMD workflow')
    
    def _validate_umbrella_workflow_prerequisites(self) -> None:
        """Validate prerequisites for umbrella sampling workflow step"""
        
        requirements = {
            'smd_completed': lambda: self._check_smd_completion(),
            'smd_results_valid': lambda: self._check_smd_results_integrity()
        }
        
        validate_prerequisites(requirements, 'umbrella sampling workflow')
    
    def _validate_analysis_workflow_prerequisites(self) -> None:
        """Validate prerequisites for analysis workflow step"""
        
        requirements = {
            'umbrella_completed': lambda: self._check_umbrella_completion(),
            'umbrella_results_sufficient': lambda: self._check_umbrella_completeness()
        }
        
        validate_prerequisites(requirements, 'analysis workflow')
    
    def _check_directory_writable(self, directory: Path) -> bool:
        """Check if directory is writable"""
        try:
            directory.mkdir(parents=True, exist_ok=True)
            test_file = directory / ".write_test"
            test_file.write_text("test")
            test_file.unlink()
            return True
        except (PermissionError, OSError):
            return False
    
    def _check_md_system_validation(self) -> bool:
        """Check if MD system passes validation"""
        try:
            validation_result = validate_md_results(self.system_dir, ValidationLevel.STANDARD)
            if not validation_result.is_valid:
                from ..utils.exceptions import SystemValidationError
                raise SystemValidationError(
                    validation_score=validation_result.score,
                    missing_files=validation_result.files_missing,
                    errors=validation_result.errors,
                    warnings=validation_result.warnings
                )
            return True
        except Exception as exc:
            logger.error(f"MD system validation failed: {exc}")
            return False
    
    def _check_smd_results_integrity(self) -> bool:
        """Check SMD results file integrity"""
        smd_dir = self.smd_dir
        
        # Check file sizes and basic integrity
        critical_files = [
            smd_dir / "results" / "smd_pullf.xvg",
            smd_dir / "results" / "smd_pullx.xvg"
        ]
        
        for file_path in critical_files:
            if not file_path.exists():
                return False
            
            # Check if file has reasonable content (not empty and has multiple lines)
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    # Should have header + data lines
                    if len(lines) < 10:  # Minimum reasonable number of lines
                        logger.warning(f"SMD result file {file_path} appears too short")
                        return False
            except Exception as exc:
                logger.error(f"Failed to read SMD result file {file_path}: {exc}")
                return False
        
        return True
    
    def _check_umbrella_completeness(self) -> bool:
        """Check if umbrella sampling has sufficient completion"""
        if not self.umbrella_dir.exists():
            return False
        
        # Count completed windows with reasonable data
        completed_windows = 0
        total_windows = len(list(self.umbrella_dir.glob("window_*")))
        
        if total_windows == 0:
            return False
        
        for window_dir in self.umbrella_dir.glob("window_*"):
            result_file = window_dir / "results" / "umbrella_pullf.xvg"
            if result_file.exists() and result_file.stat().st_size > 1000:  # Reasonable size
                completed_windows += 1
        
        completion_rate = completed_windows / total_windows
        
        if completion_rate < 0.8:  # Require at least 80% completion
            raise UmbrellaIncompleteError(
                completed_windows=completed_windows,
                total_windows=total_windows,
                failed_windows=[d.name for d in self.umbrella_dir.glob("window_*") 
                              if not (d / "results" / "umbrella_pullf.xvg").exists()]
            )
        
        return True
    
    def _load_workflow_state(self) -> None:
        """Load workflow state from file"""
        state_file = self.output_dir / "workflow_state.yaml"
        if state_file.exists():
            with open(state_file, 'r') as f:
                saved_state = yaml.safe_load(f)
                self.state.update(saved_state)
    
    def _run_smd_automated(self) -> None:
        """Run SMD simulation in automated mode"""
        os.chdir(self.smd_dir)
        subprocess.run(["bash", "run_smd.sh"], check=True)
        self.state['smd_completed'] = True
    
    def _run_umbrella_automated(self) -> None:
        """Run umbrella sampling in automated mode"""  
        os.chdir(self.umbrella_dir)
        subprocess.run(["bash", "run_all_umbrella.sh", "sequential"], check=True)
        self.state['umbrella_completed'] = True
    
    def _generate_summary_report(self, results: Dict) -> None:
        """Generate comprehensive summary report"""
        report_file = self.analysis_dir / "pmf_analysis_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("PRISM PMF Analysis Report\n")
            f.write("=" * 50 + "\n\n")
            
            if 'binding_energy' in results:
                be = results['binding_energy']
                f.write(f"Binding Energy: {be['value']:.2f} {results.get('energy_unit', 'kcal/mol')}\n")
                f.write(f"Minimum PMF: {be['min_pmf']:.2f} at {be['min_distance']:.3f} nm\n")
                f.write(f"Maximum PMF: {be['max_pmf']:.2f} at {be['max_distance']:.3f} nm\n\n")
            
            f.write(f"Analysis Directory: {self.analysis_dir}\n")
            f.write(f"PMF Data: {self.analysis_dir}/pmf.xvg\n") 
            f.write(f"Error Data: {self.analysis_dir}/pmferror.xvg\n")
            
            if 'plots' in results:
                f.write("\nGenerated Plots:\n")
                for plot_name, plot_path in results['plots'].items():
                    f.write(f"  {plot_name}: {plot_path}\n")
        
        logger.info(f"Summary report saved: {report_file}")
    
    def _generate_smd_plots(self, pullf_file: Path, pullx_file: Path) -> Dict:
        """Generate SMD analysis plots using unified plotting utilities"""
        analysis_dir = pullf_file.parent / 'analysis'
        analysis_dir.mkdir(exist_ok=True)

        plotter = PMFPlotter(analysis_dir)
        plots = plotter.plot_all_smd(pullf_file, pullx_file)

        return plots
    
    def _calculate_smd_stats(self, pullf_file: Path, pullx_file: Path) -> Dict:
        """Calculate SMD statistics"""
        import numpy as np

        stats = {}

        try:
            # Read force data
            pullf_data = np.loadtxt(pullf_file, comments=['#', '@'])
            force = pullf_data[:, 1]

            # Read distance data
            pullx_data = np.loadtxt(pullx_file, comments=['#', '@'])
            distance = pullx_data[:, 1]

            # Calculate statistics
            stats['max_force'] = float(np.max(force))
            stats['avg_force'] = float(np.mean(force))
            stats['min_force'] = float(np.min(force))
            stats['std_force'] = float(np.std(force))
            stats['total_distance'] = float(distance[-1] - distance[0])
            stats['final_distance'] = float(distance[-1])
            stats['initial_distance'] = float(distance[0])

            logger.info(f"SMD statistics calculated: max_force={stats['max_force']:.2f} kJ/mol/nm")

        except Exception as e:
            logger.warning(f"Could not calculate SMD stats: {e}")

        return stats
    
    def _generate_detailed_pmf_plots(self, pmf_file: Path) -> Dict:
        """Generate detailed PMF plots"""
        # Implementation would go here
        return {}
    
    def _calculate_pmf_statistics(self, pmf_file: Path) -> Dict:
        """Calculate PMF statistics"""
        # Implementation would go here
        return {}
    
    def _load_existing_results(self) -> Dict:
        """Load existing analysis results"""
        pmf_file = self.analysis_dir / "pmf.xvg"
        if pmf_file.exists():
            return self.analyzer.run_wham()
        else:
            raise FileNotFoundError("No existing PMF results found")

    def _inherit_system_config(self, user_config: Dict) -> Dict:
        """
        Inherit force field and system settings from original PRISM system

        Parameters:
        -----------
        user_config : Dict
            User-provided PMF configuration

        Returns:
        --------
        Dict : Combined configuration with inherited system settings
        """
        logger.info("Inheriting system configuration from original PRISM build")

        # Try to load original prism_config.yaml
        original_config_path = self.system_dir / "prism_config.yaml"
        inherited_config = {}

        if original_config_path.exists():
            try:
                with open(original_config_path, 'r') as f:
                    original_config = yaml.safe_load(f)

                # Inherit key sections from original config
                sections_to_inherit = [
                    'forcefield', 'water_model', 'ligand_forcefield',
                    'simulation', 'ions', 'constraints', 'energy_minimization',
                    'electrostatics', 'vdw', 'temperature_coupling', 'pressure_coupling'
                ]

                for section in sections_to_inherit:
                    if section in original_config:
                        inherited_config[section] = original_config[section]
                        logger.debug(f"Inherited {section} configuration from original system")

                logger.info("Successfully inherited configuration from original PRISM system")

            except Exception as e:
                logger.warning(f"Could not load original config: {e}")
                logger.info("Using PMF template defaults")
        else:
            logger.warning(f"Original config not found at: {original_config_path}")
            logger.info("Using PMF template defaults")

        # Merge user config with inherited config (user config takes precedence)
        merged_config = inherited_config.copy()
        merged_config.update(user_config)

        # Ensure PMF-specific sections exist
        pmf_defaults = self._get_pmf_defaults()
        for section, defaults in pmf_defaults.items():
            if section not in merged_config:
                merged_config[section] = defaults
                logger.debug(f"Added default PMF section: {section}")

        return merged_config

    def _get_pmf_defaults(self) -> Dict:
        """
        Get default PMF-specific configuration sections

        Returns:
        --------
        Dict : Default PMF configuration sections
        """
        return {
            'box': {
                'pulling_distance': 3.0
            },
            'pmf_system': {
                'reference_group': 'Protein',
                'moving_group': 'LIG',
                'align_axis': 'z',
                'center_system': True
            },
            'extraction': {
                'frame_number': -1,
                'output_structure': True
            },
            'equilibration': {
                'mode': 'auto',  # 'auto' or 'manual'
                'hardware': {
                    'gpu_id': 0,
                    'ntmpi': 1,
                    'ntomp': 10,
                    'gpu_acceleration': True
                },
                'em': {
                    'integrator': 'steep',
                    'emtol': 10.0,
                    'nsteps': 50000,
                    'maxwarn': 3
                },
                'nvt': {
                    'nsteps': 50000,
                    'dt': 0.002,
                    'maxwarn': 3
                },
                'npt': {
                    'nsteps': 100000,
                    'dt': 0.002,
                    'maxwarn': 3
                },
                'script': {
                    'filename': 'equilibration_run.sh',
                    'gpu_optimization': True,
                    'checkpoint_resume': True,
                    'verbose_output': True,
                    'create_directories': True
                }
            },
            'smd': {
                'pull_rate': 0.005,
                'pull_k': 1000.0,
                'nsteps': 600000,
                'dt': 0.002,
                'output_interval': 2500
            },
            'distance': {
                'start': 0.3,
                'end': 3.5
            },
            'umbrella': {
                'window_spacing': 0.1,
                'force_constant': 1000.0,
                'nsteps': 250000,
                'dt': 0.002,
                'output_interval': 2500,
                'equilibration_time': 50000
            },
            'analysis': {
                'begin_time_ps': 5000,
                'bootstrap': 100,
                'energy_unit': 'kCal',
                'tolerance': 1e-6,
                'max_iterations': 10000,
                'bin_width': 0.01
            }
        }


def cli_main():
    """CLI entry point for prism-pmf command"""
    import sys
    import argparse
    from pathlib import Path
    
    parser = argparse.ArgumentParser(
        prog='prism-pmf',
        description='PRISM PMF Workflow Direct Command',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  prism-pmf create ./md_results ./pmf_output --auto
  prism-pmf status ./pmf_output
  prism-pmf run-complete ./md_results ./pmf_output
        """
    )
    
    subparsers = parser.add_subparsers(dest='action', help='PMF workflow actions')
    
    # Create workflow
    create_parser = subparsers.add_parser('create', help='Create new PMF workflow')
    create_parser.add_argument('system_dir', type=Path, help='System directory path')
    create_parser.add_argument('output_dir', type=Path, help='Output directory path')
    create_parser.add_argument('--auto', action='store_true', help='Run complete workflow automatically')
    create_parser.add_argument('--build-pmf-system', action='store_true', default=True, help='Build PMF-optimized system')
    
    # Status check
    status_parser = subparsers.add_parser('status', help='Check workflow status')
    status_parser.add_argument('output_dir', type=Path, help='Workflow output directory')
    
    # Run complete workflow
    run_parser = subparsers.add_parser('run-complete', help='Run complete PMF workflow')
    run_parser.add_argument('system_dir', type=Path, help='System directory path')
    run_parser.add_argument('output_dir', type=Path, help='Output directory path')
    run_parser.add_argument('--build-pmf-system', action='store_true', default=True, help='Build PMF-optimized system')
    
    args = parser.parse_args()
    
    if not args.action:
        parser.print_help()
        return 1
    
    try:
        if args.action == 'create':
            config = {}  # Use default configuration
            workflow = PMFWorkflow(args.system_dir, args.output_dir, config)
            
            if args.auto:
                print("Running complete PMF workflow...")
                results = workflow.run_complete(build_pmf_system=args.build_pmf_system)
                print(f"PMF workflow completed successfully!")
                if 'binding_energy' in results:
                    be = results['binding_energy']
                    print(f"Binding energy: {be['value']:.2f} {results.get('energy_unit', 'kcal/mol')}")
            else:
                print("PMF workflow created. Use 'run-complete' to execute.")
            
            return 0
            
        elif args.action == 'status':
            config = {}
            workflow = PMFWorkflow(args.system_dir, args.output_dir, config)
            status = workflow.get_status()
            
            print(f"Workflow status: {status['current_stage']}")
            print(f"Next action: {status['next_action']}")
            
            return 0
            
        elif args.action == 'run-complete':
            config = {}
            workflow = PMFWorkflow(args.system_dir, args.output_dir, config)
            
            print("Starting complete PMF workflow...")
            results = workflow.run_complete(build_pmf_system=args.build_pmf_system)
            
            print("PMF workflow completed successfully!")
            if 'binding_energy' in results:
                be = results['binding_energy']
                print(f"Binding energy: {be['value']:.2f} {results.get('energy_unit', 'kcal/mol')}")
            
            return 0
            
    except Exception as e:
        print(f"Error: {e}")
        return 1