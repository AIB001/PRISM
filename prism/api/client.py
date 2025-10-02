#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM API Client

Client interface for programmatic access to PRISM functionality through
command-line and programmatic APIs.
"""

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
import time
import uuid

from ..utils.logging_system import PrismLogger
from .exceptions import *
from .core import PMFCalculator, WorkflowManager, DataAnalyzer


class PrismAPIClient:
    """Client for accessing PRISM PMF functionality programmatically"""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None, logger: Optional[PrismLogger] = None):
        """
        Initialize PRISM API client
        
        Args:
            config: Configuration dictionary
            logger: Optional logger instance
        """
        self.logger = logger or PrismLogger("prism_api_client")
        self.config = config or {}
        
        # Initialize core components
        self.pmf_calculator = PMFCalculator(logger)
        self.workflow_manager = WorkflowManager(logger)
        self.data_analyzer = DataAnalyzer(logger)
        
        # Client settings
        self.timeout = self.config.get('timeout', 3600)  # 1 hour default
        self.working_directory = Path(self.config.get('working_directory', '.'))
        self.temp_directory = Path(self.config.get('temp_directory', tempfile.gettempdir()))
        
        self.logger.info("PRISM API client initialized")
    
    # High-level convenience methods
    def run_pmf_calculation(self, system_file: str, ligand: str, protein: str, 
                           method: str = "umbrella_sampling", **kwargs) -> Dict[str, Any]:
        """
        Run complete PMF calculation with single method call
        
        Args:
            system_file: Path to system structure file
            ligand: Ligand residue name
            protein: Protein residue name
            method: PMF calculation method
            **kwargs: Additional configuration parameters
            
        Returns:
            Dictionary containing calculation results
        """
        try:
            self.logger.info(f"Starting PMF calculation: {ligand}-{protein}")
            
            # Run calculation using PMF calculator
            result = self.pmf_calculator.run_pmf_calculation(
                system_file=system_file,
                ligand=ligand, 
                protein=protein,
                method=method,
                **kwargs
            )
            
            self.logger.info(f"PMF calculation completed: {result.calculation_id}")
            return result.to_dict()
            
        except Exception as e:
            self.logger.error(f"PMF calculation failed: {e}")
            raise CalculationError(f"PMF calculation failed: {e}")
    
    def submit_workflow(self, config: Dict[str, Any], priority: str = "normal") -> Dict[str, Any]:
        """
        Submit workflow for execution
        
        Args:
            config: Workflow configuration dictionary
            priority: Workflow priority level
            
        Returns:
            Dictionary containing workflow information
        """
        try:
            from .core import WorkflowPriority
            
            # Convert priority string to enum
            priority_enum = WorkflowPriority(priority.lower())
            
            # Create and submit workflow
            workflow_id = self.workflow_manager.create_workflow(config, priority_enum)
            self.workflow_manager.submit_workflow(workflow_id)
            
            # Get workflow status
            status = self.workflow_manager.get_status(workflow_id)
            
            self.logger.info(f"Workflow submitted: {workflow_id}")
            return status
            
        except Exception as e:
            self.logger.error(f"Workflow submission failed: {e}")
            raise WorkflowError(f"Workflow submission failed: {e}")
    
    def get_workflow_status(self, workflow_id: str) -> Dict[str, Any]:
        """
        Get workflow execution status
        
        Args:
            workflow_id: Workflow identifier
            
        Returns:
            Dictionary containing workflow status
        """
        try:
            return self.workflow_manager.get_status(workflow_id)
        except Exception as e:
            self.logger.error(f"Failed to get workflow status: {e}")
            raise WorkflowError(f"Failed to get workflow status: {e}", workflow_id)
    
    def analyze_pmf_data(self, pmf_file: str, output_format: str = "dict") -> Union[Dict[str, Any], str]:
        """
        Analyze PMF profile data
        
        Args:
            pmf_file: Path to PMF data file
            output_format: Output format ('dict', 'json')
            
        Returns:
            Analysis results in specified format
        """
        try:
            # Load PMF data (simplified implementation)
            pmf_data = self._load_pmf_file(pmf_file)
            
            # Analyze data
            analysis = self.data_analyzer.analyze_pmf_profile(pmf_data)
            
            if output_format == "json":
                return json.dumps(analysis, indent=2)
            else:
                return analysis
                
        except Exception as e:
            self.logger.error(f"PMF data analysis failed: {e}")
            raise DataError(f"PMF analysis failed: {e}", file_path=pmf_file)
    
    # CLI integration methods
    def call_prism_cli(self, command: List[str], input_data: Optional[str] = None,
                       capture_output: bool = True, timeout: Optional[int] = None) -> Dict[str, Any]:
        """
        Call PRISM CLI command programmatically
        
        Args:
            command: CLI command as list of arguments
            input_data: Optional input data to pass to command
            capture_output: Whether to capture stdout/stderr
            timeout: Command timeout in seconds
            
        Returns:
            Dictionary containing command results
        """
        try:
            timeout = timeout or self.timeout
            
            # Construct full command
            full_command = ['prism'] + command
            
            self.logger.info(f"Executing CLI command: {' '.join(full_command)}")
            
            # Execute command
            result = subprocess.run(
                full_command,
                input=input_data,
                capture_output=capture_output,
                text=True,
                timeout=timeout,
                cwd=str(self.working_directory)
            )
            
            # Parse result
            output = {
                'returncode': result.returncode,
                'stdout': result.stdout if capture_output else None,
                'stderr': result.stderr if capture_output else None,
                'success': result.returncode == 0
            }
            
            if result.returncode != 0:
                self.logger.error(f"CLI command failed: {result.stderr}")
                raise CalculationError(f"CLI command failed: {result.stderr}", 
                                     stage="cli_execution", exit_code=result.returncode)
            
            self.logger.info("CLI command completed successfully")
            return output
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"CLI command timed out after {timeout} seconds")
            raise TimeoutError(f"CLI command timed out", timeout, "cli_command")
        except Exception as e:
            self.logger.error(f"CLI command execution failed: {e}")
            raise APIError(f"CLI command failed: {e}")
    
    def calculate_pmf_cli(self, system_file: str, ligand: str, protein: str,
                         output_dir: Optional[str] = None, **kwargs) -> Dict[str, Any]:
        """
        Run PMF calculation using CLI interface
        
        Args:
            system_file: Path to system file
            ligand: Ligand residue name
            protein: Protein residue name
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Dictionary containing calculation results
        """
        # Build CLI command
        command = [
            'api', 'calculate-pmf',
            '--input', system_file,
            '--ligand', ligand,
            '--protein', protein,
            '--format', 'json'
        ]
        
        if output_dir:
            command.extend(['--output', output_dir])
        
        # Add additional parameters
        for key, value in kwargs.items():
            command.extend([f'--{key.replace("_", "-")}', str(value)])
        
        # Execute command
        result = self.call_prism_cli(command)
        
        # Parse JSON output
        if result['success'] and result['stdout']:
            try:
                return json.loads(result['stdout'])
            except json.JSONDecodeError:
                return {'raw_output': result['stdout']}
        else:
            raise CalculationError(f"PMF calculation failed: {result['stderr']}")
    
    def get_status_cli(self, workflow_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Get system or workflow status using CLI interface
        
        Args:
            workflow_id: Optional workflow ID to get specific status
            
        Returns:
            Dictionary containing status information
        """
        command = ['api', 'get-status', '--format', 'json']
        
        if workflow_id:
            command.extend(['--workflow-id', workflow_id])
        
        result = self.call_prism_cli(command)
        
        if result['success'] and result['stdout']:
            try:
                return json.loads(result['stdout'])
            except json.JSONDecodeError:
                return {'raw_output': result['stdout']}
        else:
            raise APIError(f"Status request failed: {result['stderr']}")
    
    def submit_workflow_cli(self, config_file: str, async_mode: bool = False) -> Dict[str, Any]:
        """
        Submit workflow using CLI interface
        
        Args:
            config_file: Path to workflow configuration file
            async_mode: Whether to run asynchronously
            
        Returns:
            Dictionary containing submission results
        """
        command = ['api', 'submit-workflow', '--config', config_file, '--format', 'json']
        
        if async_mode:
            command.append('--async')
        
        result = self.call_prism_cli(command)
        
        if result['success'] and result['stdout']:
            try:
                return json.loads(result['stdout'])
            except json.JSONDecodeError:
                return {'raw_output': result['stdout']}
        else:
            raise WorkflowError(f"Workflow submission failed: {result['stderr']}")
    
    # Batch operations
    def batch_pmf_calculations(self, calculation_configs: List[Dict[str, Any]], 
                              max_concurrent: int = 5) -> List[Dict[str, Any]]:
        """
        Run multiple PMF calculations in batch
        
        Args:
            calculation_configs: List of calculation configuration dictionaries
            max_concurrent: Maximum concurrent calculations
            
        Returns:
            List of calculation results
        """
        results = []
        
        self.logger.info(f"Starting batch PMF calculations: {len(calculation_configs)} jobs")
        
        # Simple sequential execution for demo (real implementation would use threading/async)
        for i, config in enumerate(calculation_configs):
            try:
                self.logger.info(f"Processing calculation {i+1}/{len(calculation_configs)}")
                
                result = self.run_pmf_calculation(**config)
                result['batch_index'] = i
                results.append(result)
                
            except Exception as e:
                self.logger.error(f"Batch calculation {i+1} failed: {e}")
                results.append({
                    'batch_index': i,
                    'error': str(e),
                    'status': 'failed'
                })
        
        self.logger.info(f"Batch calculations completed: {len(results)} results")
        return results
    
    def export_results(self, results: Union[Dict, List[Dict]], 
                      output_file: str, format: str = "json") -> bool:
        """
        Export calculation results to file
        
        Args:
            results: Results data to export
            output_file: Output file path
            format: Export format ('json', 'yaml', 'csv')
            
        Returns:
            True if export successful
        """
        try:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            if format == "json":
                with open(output_path, 'w') as f:
                    json.dump(results, f, indent=2)
            elif format == "yaml":
                import yaml
                with open(output_path, 'w') as f:
                    yaml.dump(results, f, default_flow_style=False)
            elif format == "csv":
                self._export_to_csv(results, output_path)
            else:
                raise DataError(f"Unsupported export format: {format}", data_format=format)
            
            self.logger.info(f"Results exported to: {output_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"Export failed: {e}")
            raise DataError(f"Export failed: {e}", file_path=output_file)
    
    # Helper methods
    def _load_pmf_file(self, pmf_file: str) -> List[Dict[str, float]]:
        """Load PMF data from file (simplified implementation)"""
        pmf_path = Path(pmf_file)
        
        if not pmf_path.exists():
            raise DataError(f"PMF file not found: {pmf_file}", file_path=pmf_file)
        
        # Mock data loading - in real implementation would parse actual file format
        mock_data = []
        for i in range(30):  # 30 data points
            x = i * 0.1
            energy = 20 * ((x - 1.5) ** 2) - 8.5  # Gaussian-like potential
            error = 0.3 + 0.2 * abs(x - 1.5)
            
            mock_data.append({
                'reaction_coordinate': x,
                'pmf_energy': energy,
                'error_estimate': error,
                'sampling_count': 10000 - i * 100
            })
        
        return mock_data
    
    def _export_to_csv(self, results: Union[Dict, List[Dict]], output_path: Path):
        """Export results to CSV format"""
        import csv
        
        # Convert single dict to list
        if isinstance(results, dict):
            results = [results]
        
        if not results:
            return
        
        # Get all unique keys from all result dictionaries
        all_keys = set()
        for result in results:
            all_keys.update(self._flatten_dict(result).keys())
        
        # Write CSV
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=sorted(all_keys))
            writer.writeheader()
            
            for result in results:
                flattened = self._flatten_dict(result)
                writer.writerow(flattened)
    
    def _flatten_dict(self, d: Dict[str, Any], parent_key: str = '', sep: str = '.') -> Dict[str, Any]:
        """Flatten nested dictionary for CSV export"""
        items = []
        
        for k, v in d.items():
            new_key = f"{parent_key}{sep}{k}" if parent_key else k
            
            if isinstance(v, dict):
                items.extend(self._flatten_dict(v, new_key, sep=sep).items())
            elif isinstance(v, list):
                # Convert list to string representation
                items.append((new_key, str(v)))
            else:
                items.append((new_key, v))
        
        return dict(items)
    
    def __enter__(self):
        """Context manager entry"""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.logger.info("PRISM API client session ended")


# Convenience function for quick API access
def create_client(config: Optional[Dict[str, Any]] = None) -> PrismAPIClient:
    """
    Create PRISM API client with default configuration
    
    Args:
        config: Optional configuration dictionary
        
    Returns:
        Configured PrismAPIClient instance
    """
    return PrismAPIClient(config)