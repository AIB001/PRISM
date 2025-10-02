#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM CLI Main Entry Point

Main command-line interface for PRISM PMF calculations with comprehensive
subcommand system for workflow management, monitoring, and analysis.
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

from ..utils.logging_system import PrismLogger
from .utils import CLIFormatter, ProgressBar
from .commands import (
    WorkflowCommands, MonitorCommands, ConfigCommands,
    AnalysisCommands, SystemCommands
)


@dataclass
class CLIContext:
    """Context object passed to all CLI commands"""
    config_path: Optional[Path] = None
    verbose: bool = False
    quiet: bool = False
    output_format: str = "table"
    working_directory: Path = Path.cwd()
    log_level: str = "INFO"


class PrismCLI:
    """Main PRISM CLI application"""
    
    def __init__(self):
        self.logger = PrismLogger("prism_cli")
        self.formatter = CLIFormatter()
        self.context = CLIContext()
        
        # Command handlers
        self.workflow_commands = WorkflowCommands(self.context)
        self.monitor_commands = MonitorCommands(self.context)
        self.config_commands = ConfigCommands(self.context)
        self.analysis_commands = AnalysisCommands(self.context)
        self.system_commands = SystemCommands(self.context)
        
        # Setup argument parser
        self.parser = self._create_parser()
    
    def _create_parser(self) -> argparse.ArgumentParser:
        """Create the main argument parser with subcommands"""
        parser = argparse.ArgumentParser(
            prog='prism',
            description='PRISM PMF Calculation System - Advanced molecular dynamics workflow management',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  prism workflow create pmf-system.gro --ligand LIG --protein Protein
  prism workflow status workflow_001
  prism monitor system --real-time
  prism analysis plot-pmf results/pmf_profile.dat
  prism config validate config.yaml
  
For more information on each command, use:
  prism <command> --help
            """
        )
        
        # Global options
        parser.add_argument(
            '--version', action='version',
            version='PRISM PMF System v1.0.0'
        )
        
        parser.add_argument(
            '-v', '--verbose', action='store_true',
            help='Enable verbose output'
        )
        
        parser.add_argument(
            '-q', '--quiet', action='store_true',
            help='Suppress non-essential output'
        )
        
        parser.add_argument(
            '--config', type=Path,
            help='Path to configuration file'
        )
        
        parser.add_argument(
            '--output-format', choices=['table', 'json', 'yaml'],
            default='table', help='Output format for structured data'
        )
        
        parser.add_argument(
            '--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
            default='INFO', help='Logging level'
        )
        
        parser.add_argument(
            '--working-dir', type=Path, default=Path.cwd(),
            help='Working directory for operations'
        )
        
        # Subcommands
        subparsers = parser.add_subparsers(
            dest='command', help='Available commands',
            metavar='COMMAND'
        )
        
        # Workflow management commands
        workflow_parser = subparsers.add_parser(
            'workflow', help='Workflow management commands',
            description='Create, manage, and monitor PMF calculation workflows'
        )
        self.workflow_commands.setup_parser(workflow_parser)
        
        # System monitoring commands
        monitor_parser = subparsers.add_parser(
            'monitor', help='System monitoring commands',
            description='Monitor system resources, workflows, and performance'
        )
        self.monitor_commands.setup_parser(monitor_parser)
        
        # Configuration management commands
        config_parser = subparsers.add_parser(
            'config', help='Configuration management commands',
            description='Manage system configuration and settings'
        )
        self.config_commands.setup_parser(config_parser)
        
        # Analysis and visualization commands
        analysis_parser = subparsers.add_parser(
            'analysis', help='Analysis and reporting commands',
            description='Analyze PMF results and generate reports'
        )
        self.analysis_commands.setup_parser(analysis_parser)
        
        # System administration commands
        system_parser = subparsers.add_parser(
            'system', help='System administration commands',
            description='System administration and maintenance tools'
        )
        self.system_commands.setup_parser(system_parser)
        
        return parser
    
    def _setup_logging(self):
        """Setup logging based on CLI arguments"""
        log_level = getattr(logging, self.context.log_level.upper())
        
        if self.context.quiet:
            log_level = logging.ERROR
        elif self.context.verbose:
            log_level = logging.DEBUG
        
        # Configure the logger
        logger = logging.getLogger('prism')
        logger.setLevel(log_level)
        
        # Create console handler if not exists
        if not logger.handlers:
            handler = logging.StreamHandler(sys.stderr)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
    
    def _update_context(self, args):
        """Update CLI context from parsed arguments"""
        self.context.config_path = args.config
        self.context.verbose = args.verbose
        self.context.quiet = args.quiet
        self.context.output_format = args.output_format
        self.context.working_directory = args.working_dir
        self.context.log_level = args.log_level
    
    def run(self, argv: Optional[List[str]] = None) -> int:
        """Run the CLI application"""
        try:
            # Parse arguments
            args = self.parser.parse_args(argv)
            
            # Update context
            self._update_context(args)
            
            # Setup logging
            self._setup_logging()
            
            # Handle no command case
            if not hasattr(args, 'command') or args.command is None:
                self.parser.print_help()
                return 1
            
            # Route to appropriate command handler
            if args.command == 'workflow':
                return self.workflow_commands.handle(args)
            elif args.command == 'monitor':
                return self.monitor_commands.handle(args)
            elif args.command == 'config':
                return self.config_commands.handle(args)
            elif args.command == 'analysis':
                return self.analysis_commands.handle(args)
            elif args.command == 'system':
                return self.system_commands.handle(args)
            else:
                self.formatter.error(f"Unknown command: {args.command}")
                return 1
                
        except KeyboardInterrupt:
            self.formatter.info("Operation cancelled by user")
            return 130
        except Exception as e:
            if self.context.verbose:
                self.logger.exception("CLI error occurred")
            else:
                self.formatter.error(f"Error: {e}")
            return 1
    
    def print_banner(self):
        """Print the application banner"""
        banner = """
╔══════════════════════════════════════════════════════════════╗
║                    PRISM PMF System                          ║
║              Advanced Molecular Dynamics                     ║
║            Potential of Mean Force Calculator                ║
║                                                              ║
║  🧪 Intelligent Workflow Management                          ║
║  ⚡ High-Performance Computing                              ║  
║  📊 Real-time Monitoring & Analysis                         ║
║  🔧 Automated Configuration & Optimization                  ║
╚══════════════════════════════════════════════════════════════╝
        """
        print(banner)
    
    def print_quick_start(self):
        """Print quick start guide"""
        quick_start = """
🚀 Quick Start Guide:

1. Create a new PMF workflow:
   prism workflow create system.gro --ligand LIG --protein Protein

2. Monitor system status:
   prism monitor system --dashboard

3. Check workflow progress:
   prism workflow list --status running

4. Analyze results:
   prism analysis plot-pmf results/pmf_profile.dat

5. Get help for any command:
   prism <command> --help

For detailed documentation: prism help
        """
        print(quick_start)


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point for the PRISM CLI"""
    cli = PrismCLI()
    
    # Print banner if no arguments provided
    if not argv and len(sys.argv) == 1:
        cli.print_banner()
        cli.print_quick_start()
        return 0
    
    return cli.run(argv)


if __name__ == '__main__':
    sys.exit(main())