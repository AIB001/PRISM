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

from .smd import SMDManager
from .umbrella import UmbrellaManager  
from .analyzer import PMFAnalyzer
from .md_validator import MDResultsValidator, ValidationLevel, validate_md_results

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
            Path to system directory containing GMX_PROLIG_MD
        output_dir : str or Path
            Output directory for PMF calculations
        config : Dict
            PMF configuration dictionary
        """
        self.system_dir = Path(system_dir)
        self.output_dir = Path(output_dir)
        self.config = config
        
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
        """Fallback MD results location logic"""
        # Check if system_dir itself contains GMX files
        if (self.system_dir / "solv_ions.gro").exists():
            return self.system_dir
        
        # Look for GMX_PROLIG_MD subdirectory
        gmx_dir = self.system_dir / "GMX_PROLIG_MD"
        if gmx_dir.exists():
            return gmx_dir
            
        # Search in subdirectories
        for subdir in self.system_dir.iterdir():
            if subdir.is_dir():
                candidate = subdir / "GMX_PROLIG_MD"
                if candidate.exists():
                    return candidate
        
        raise FileNotFoundError(
            f"Could not locate valid MD results in {self.system_dir}. "
            "Please ensure the system has been built with PRISM."
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
        from .pmf_builder import PMFBuilder
        
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
        """
        logger.info("=== Preparing SMD Simulation ===")
        
        # Check if PMF system should be used but doesn't exist
        if use_pmf_system and not self.state['pmf_system_built']:
            logger.warning("PMF system not built yet. Consider running build_pmf_system() first.")
            logger.warning("Proceeding with PMF system build...")
            self.build_pmf_system(equilibrate=True)
        
        result = self.smd_manager.prepare(use_pmf_rebuilt_system=use_pmf_system)
        self.state['smd_prepared'] = True
        
        # Save state
        self._save_workflow_state()
        
        logger.info("SMD preparation completed")
        logger.info(f"System type: {'PMF-rebuilt (Z-optimized)' if use_pmf_system else 'Original MD'}")
        logger.info(f"Run SMD: cd {result['smd_dir']} && bash run_smd.sh")
        
        return result
    
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
        required_files = [
            self.smd_dir / "results" / "smd_pullf.xvg",
            self.smd_dir / "results" / "smd_pullx.xvg",
            self.smd_dir / "results" / "smd.gro"
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
            'equilibration_structure': (self.output_dir / "equilibration" / "npt" / "npt_final.gro").exists(),
            'equilibration_checkpoint': (self.output_dir / "equilibration" / "npt" / "npt_final.cpt").exists()
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
        """Generate SMD analysis plots"""
        # Implementation would go here
        return {}
    
    def _calculate_smd_stats(self, pullf_file: Path, pullx_file: Path) -> Dict:
        """Calculate SMD statistics"""
        # Implementation would go here
        return {}
    
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