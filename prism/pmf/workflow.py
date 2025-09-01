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

logger = logging.getLogger(__name__)


class PMFWorkflow:
    """
    PRISM PMF Workflow Manager
    
    Manages the complete PMF calculation workflow with support for both
    automated and step-by-step execution modes, following PRISM patterns.
    """
    
    def __init__(self, system_dir: Union[str, Path], output_dir: Union[str, Path], 
                 config: Dict, work_in_place: bool = True):
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
        work_in_place : bool, optional
            Whether to work directly in MD results directory (default: True)
        """
        self.system_dir = Path(system_dir)
        self.output_dir = Path(output_dir)
        self.config = config
        self.work_in_place = work_in_place
        
        # Locate MD results directory
        self.md_results_dir = self._locate_md_results()
        
        # Create output structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if work_in_place:
            # Work directly in MD results directory, but output results to specified path
            # This avoids file copying and potential missing files
            self.smd_dir = self.md_results_dir.parent / "pmf_smd"
            self.umbrella_dir = self.md_results_dir.parent / "pmf_umbrella"
            self.analysis_dir = self.output_dir / "analysis"
        else:
            # Traditional approach: separate working directories
            self.smd_dir = self.output_dir / "smd"
            self.umbrella_dir = self.output_dir / "umbrella"  
            self.analysis_dir = self.output_dir / "analysis"
        
        # Initialize managers
        self.smd_manager = SMDManager(self)
        self.umbrella_manager = UmbrellaManager(self)
        self.analyzer = PMFAnalyzer(self)
        
        # Track workflow state
        self.state = {
            'smd_prepared': False,
            'smd_completed': False,
            'umbrella_prepared': False, 
            'umbrella_completed': False,
            'analysis_completed': False
        }
        
        logger.info(f"PMF workflow initialized")
        logger.info(f"  System: {self.md_results_dir}")
        logger.info(f"  Output: {self.output_dir}")
    
    def _locate_md_results(self) -> Path:
        """Locate MD results directory following PRISM patterns"""
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
            f"Could not locate GMX_PROLIG_MD in {self.system_dir}. "
            "Please ensure the system has been built with PRISM."
        )
    
    def prepare_smd(self) -> Dict:
        """
        Prepare SMD simulation (Step 1)
        
        Returns:
        --------
        Dict : SMD preparation results with execution instructions
        """
        logger.info("=== Preparing SMD Simulation ===")
        
        result = self.smd_manager.prepare()
        self.state['smd_prepared'] = True
        
        # Save state
        self._save_workflow_state()
        
        logger.info("SMD preparation completed")
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
    
    def run_complete(self) -> Dict:
        """
        Run complete PMF workflow (automated mode)
        
        Returns:
        --------
        Dict : Complete workflow results
        """
        logger.info("=== Running Complete PMF Workflow ===")
        
        try:
            # Step 1: Prepare and run SMD
            if not self.state['smd_completed']:
                smd_prep = self.prepare_smd() 
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
            
            logger.info("Complete PMF workflow finished successfully!")
            return results
            
        except Exception as e:
            logger.error(f"PMF workflow failed: {e}")
            status = self.get_status()
            logger.info(f"Workflow stopped at: {status['current_stage']}")
            raise
    
    def get_status(self) -> Dict:
        """Get comprehensive workflow status"""
        # Determine current stage
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
        else:
            current_stage = "ready_to_start"
        
        # Check file status
        file_status = self._check_file_status()
        
        # Get next action
        next_actions = {
            "ready_to_start": "Call prepare_smd()",
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
        """Check status of key workflow files"""
        return {
            'system_files': {
                'structure': (self.md_results_dir / "solv_ions.gro").exists(),
                'topology': (self.md_results_dir / "topol.top").exists()
            },
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