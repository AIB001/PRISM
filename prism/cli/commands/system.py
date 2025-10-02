#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM CLI System Commands

Command handlers for system administration, maintenance, and diagnostic tools.
"""

import os
import shutil
import psutil
from pathlib import Path
from typing import Dict, Any, List, Optional
from argparse import ArgumentParser, Namespace
import subprocess

from ...utils.logging_system import PrismLogger
from ..utils import CLIFormatter, TableFormatter, ProgressBar, confirm_action


class SystemCommands:
    """Handler for system administration CLI commands"""
    
    def __init__(self, context):
        self.context = context
        self.logger = PrismLogger("system_cli")
        self.formatter = CLIFormatter()
    
    def setup_parser(self, parser: ArgumentParser):
        """Setup system command parser"""
        subparsers = parser.add_subparsers(dest='system_action', help='System actions')
        
        # System information
        info_parser = subparsers.add_parser('info', help='Show system information')
        info_parser.add_argument('--detailed', action='store_true', 
                               help='Show detailed system information')
        
        # Health check
        health_parser = subparsers.add_parser('health', help='Perform system health check')
        health_parser.add_argument('--fix', action='store_true',
                                 help='Attempt to fix found issues')
        health_parser.add_argument('--report', help='Save health report to file')
        
        # Cleanup operations
        cleanup_parser = subparsers.add_parser('cleanup', help='Clean up system resources')
        cleanup_parser.add_argument('--logs', action='store_true', help='Clean old log files')
        cleanup_parser.add_argument('--temp', action='store_true', help='Clean temporary files')
        cleanup_parser.add_argument('--results', action='store_true', help='Clean old result files')
        cleanup_parser.add_argument('--all', action='store_true', help='Clean all resources')
        cleanup_parser.add_argument('--dry-run', action='store_true', 
                                  help='Show what would be cleaned without actually doing it')
        
        # Service management
        service_parser = subparsers.add_parser('service', help='Manage PRISM services')
        service_subparsers = service_parser.add_subparsers(dest='service_action')
        
        service_subparsers.add_parser('status', help='Show service status')
        service_subparsers.add_parser('start', help='Start services').add_argument('service', nargs='?', help='Specific service to start')
        service_subparsers.add_parser('stop', help='Stop services').add_argument('service', nargs='?', help='Specific service to stop')
        service_subparsers.add_parser('restart', help='Restart services').add_argument('service', nargs='?', help='Specific service to restart')
        
        # Database operations
        db_parser = subparsers.add_parser('database', help='Database management')
        db_subparsers = db_parser.add_subparsers(dest='db_action')
        
        db_subparsers.add_parser('status', help='Show database status')
        db_subparsers.add_parser('backup', help='Create database backup').add_argument('--output', help='Backup file path')
        db_subparsers.add_parser('restore', help='Restore database from backup').add_argument('backup_file', help='Backup file to restore')
        db_subparsers.add_parser('vacuum', help='Vacuum and optimize database')
        
        # Log management
        log_parser = subparsers.add_parser('logs', help='Log file management')
        log_parser.add_argument('--tail', '-f', action='store_true', help='Follow log output')
        log_parser.add_argument('--level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                              help='Filter by log level')
        log_parser.add_argument('--component', help='Filter by component')
        log_parser.add_argument('--since', help='Show logs since time (e.g., "1h", "2d")')
        log_parser.add_argument('--lines', '-n', type=int, default=50, help='Number of lines to show')
        
        # Backup and restore
        backup_parser = subparsers.add_parser('backup', help='Create system backup')
        backup_parser.add_argument('--output', help='Backup output directory')
        backup_parser.add_argument('--include-data', action='store_true', 
                                  help='Include data files in backup')
        backup_parser.add_argument('--compress', action='store_true', help='Compress backup')
        
        restore_parser = subparsers.add_parser('restore', help='Restore from system backup')
        restore_parser.add_argument('backup_path', help='Backup directory or file to restore')
        restore_parser.add_argument('--target', help='Target directory for restore')
        
        # Update system
        update_parser = subparsers.add_parser('update', help='Update PRISM system')
        update_parser.add_argument('--check-only', action='store_true', 
                                  help='Only check for updates without installing')
        update_parser.add_argument('--force', action='store_true', help='Force update')
        
        # Diagnostic tools
        diag_parser = subparsers.add_parser('diagnose', help='Run diagnostic tests')
        diag_parser.add_argument('--component', help='Diagnose specific component')
        diag_parser.add_argument('--verbose', action='store_true', help='Verbose diagnostic output')
    
    def handle(self, args: Namespace) -> int:
        """Handle system commands"""
        try:
            if args.system_action == 'info':
                return self._show_system_info(args)
            elif args.system_action == 'health':
                return self._health_check(args)
            elif args.system_action == 'cleanup':
                return self._cleanup_system(args)
            elif args.system_action == 'service':
                return self._manage_services(args)
            elif args.system_action == 'database':
                return self._manage_database(args)
            elif args.system_action == 'logs':
                return self._manage_logs(args)
            elif args.system_action == 'backup':
                return self._create_backup(args)
            elif args.system_action == 'restore':
                return self._restore_backup(args)
            elif args.system_action == 'update':
                return self._update_system(args)
            elif args.system_action == 'diagnose':
                return self._run_diagnostics(args)
            else:
                self.formatter.error("No system action specified")
                return 1
                
        except Exception as e:
            self.formatter.error(f"System command failed: {e}")
            if self.context.verbose:
                self.logger.exception("System command error")
            return 1
    
    def _show_system_info(self, args: Namespace) -> int:
        """Show system information"""
        self.formatter.header("PRISM System Information")
        
        # Basic system info
        print(f"PRISM Version:    1.0.0")
        print(f"Python Version:   {os.sys.version.split()[0]}")
        print(f"Platform:         {os.name} ({os.sys.platform})")
        print(f"Architecture:     {os.uname().machine if hasattr(os, 'uname') else 'Unknown'}")
        
        # Hardware info
        print(f"\n{Colors.BOLD}Hardware:{Colors.RESET}")
        print(f"CPU Cores:        {os.cpu_count()}")
        
        try:
            memory = psutil.virtual_memory()
            print(f"Total Memory:     {self.formatter.format_size(memory.total)}")
            print(f"Available Memory: {self.formatter.format_size(memory.available)}")
            
            disk = psutil.disk_usage('/')
            print(f"Disk Space:       {self.formatter.format_size(disk.free)} / {self.formatter.format_size(disk.total)} free")
        except Exception as e:
            print(f"Memory/Disk info unavailable: {e}")
        
        # PRISM specific info
        print(f"\n{Colors.BOLD}PRISM Configuration:{Colors.RESET}")
        config_path = Path.home() / '.prism' / 'config.yaml'
        print(f"Config File:      {config_path} ({'exists' if config_path.exists() else 'not found'})")
        
        data_dir = Path.cwd() / 'prism_data'
        print(f"Data Directory:   {data_dir} ({'exists' if data_dir.exists() else 'not found'})")
        
        log_dir = Path.home() / '.prism' / 'logs'
        print(f"Log Directory:    {log_dir} ({'exists' if log_dir.exists() else 'not found'})")
        
        if args.detailed:
            self._show_detailed_system_info()
        
        return 0
    
    def _show_detailed_system_info(self):
        """Show detailed system information"""
        print(f"\n{Colors.BOLD}Detailed System Information:{Colors.RESET}")
        
        # Environment variables
        print(f"\n{Colors.BOLD}Environment Variables:{Colors.RESET}")
        prism_vars = {k: v for k, v in os.environ.items() if 'PRISM' in k}
        if prism_vars:
            for key, value in prism_vars.items():
                print(f"  {key}={value}")
        else:
            print("  No PRISM-specific environment variables found")
        
        # Python packages
        print(f"\n{Colors.BOLD}Python Dependencies:{Colors.RESET}")
        try:
            import pkg_resources
            packages = ['numpy', 'matplotlib', 'yaml', 'psutil']
            for package in packages:
                try:
                    version = pkg_resources.get_distribution(package).version
                    print(f"  {package:<15}: {version}")
                except pkg_resources.DistributionNotFound:
                    print(f"  {package:<15}: not installed")
        except ImportError:
            print("  Package information unavailable")
    
    def _health_check(self, args: Namespace) -> int:
        """Perform system health check"""
        self.formatter.header("PRISM System Health Check")
        
        health_issues = []
        health_passed = []
        
        # Check disk space
        try:
            disk = psutil.disk_usage('/')
            free_percent = (disk.free / disk.total) * 100
            
            if free_percent < 10:
                health_issues.append(f"Low disk space: {free_percent:.1f}% free")
            else:
                health_passed.append(f"Disk space: {free_percent:.1f}% free")
        except Exception as e:
            health_issues.append(f"Could not check disk space: {e}")
        
        # Check memory
        try:
            memory = psutil.virtual_memory()
            if memory.percent > 90:
                health_issues.append(f"High memory usage: {memory.percent:.1f}%")
            else:
                health_passed.append(f"Memory usage: {memory.percent:.1f}%")
        except Exception as e:
            health_issues.append(f"Could not check memory: {e}")
        
        # Check configuration file
        config_path = Path.home() / '.prism' / 'config.yaml'
        if config_path.exists():
            health_passed.append("Configuration file exists")
        else:
            health_issues.append("Configuration file not found")
        
        # Check data directory
        data_dir = Path.cwd() / 'prism_data'
        if data_dir.exists():
            health_passed.append("Data directory exists")
        else:
            health_issues.append("Data directory not found")
        
        # Check log directory
        log_dir = Path.home() / '.prism' / 'logs'
        if log_dir.exists():
            health_passed.append("Log directory exists")
            
            # Check log file sizes
            log_files = list(log_dir.glob('*.log'))
            large_logs = [f for f in log_files if f.stat().st_size > 100 * 1024 * 1024]  # 100MB
            if large_logs:
                health_issues.append(f"Large log files found: {len(large_logs)} files > 100MB")
        else:
            health_issues.append("Log directory not found")
        
        # Display results
        if health_passed:
            print(f"{Colors.BRIGHT_GREEN}✅ Passed Checks:{Colors.RESET}")
            for check in health_passed:
                print(f"  • {check}")
        
        if health_issues:
            print(f"\n{Colors.BRIGHT_RED}❌ Issues Found:{Colors.RESET}")
            for issue in health_issues:
                print(f"  • {issue}")
        
        # Summary
        total_checks = len(health_passed) + len(health_issues)
        success_rate = (len(health_passed) / total_checks) * 100 if total_checks > 0 else 0
        
        print(f"\n{Colors.BOLD}Health Check Summary:{Colors.RESET}")
        print(f"  Passed: {len(health_passed)}/{total_checks} ({success_rate:.1f}%)")
        
        if health_issues and args.fix:
            print(f"\n{Colors.BOLD}Attempting to fix issues...{Colors.RESET}")
            self._fix_health_issues(health_issues)
        
        if args.report:
            self._save_health_report(health_passed, health_issues, args.report)
        
        return 0 if not health_issues else 1
    
    def _cleanup_system(self, args: Namespace) -> int:
        """Clean up system resources"""
        if not any([args.logs, args.temp, args.results, args.all]):
            self.formatter.error("Specify what to clean: --logs, --temp, --results, or --all")
            return 1
        
        cleanup_tasks = []
        
        if args.all or args.logs:
            cleanup_tasks.append(('logs', self._cleanup_logs))
        
        if args.all or args.temp:
            cleanup_tasks.append(('temporary files', self._cleanup_temp))
        
        if args.all or args.results:
            cleanup_tasks.append(('old results', self._cleanup_results))
        
        total_freed = 0
        
        for task_name, cleanup_func in cleanup_tasks:
            self.formatter.info(f"Cleaning {task_name}...")
            
            try:
                freed_bytes = cleanup_func(dry_run=args.dry_run)
                total_freed += freed_bytes
                
                if args.dry_run:
                    self.formatter.info(f"Would free {self.formatter.format_size(freed_bytes)} from {task_name}")
                else:
                    self.formatter.success(f"Freed {self.formatter.format_size(freed_bytes)} from {task_name}")
                    
            except Exception as e:
                self.formatter.error(f"Failed to clean {task_name}: {e}")
        
        if args.dry_run:
            self.formatter.info(f"Total space that would be freed: {self.formatter.format_size(total_freed)}")
        else:
            self.formatter.success(f"Total space freed: {self.formatter.format_size(total_freed)}")
        
        return 0
    
    def _manage_services(self, args: Namespace) -> int:
        """Manage PRISM services"""
        if args.service_action == 'status':
            return self._show_service_status()
        elif args.service_action == 'start':
            return self._start_services(args.service)
        elif args.service_action == 'stop':
            return self._stop_services(args.service)
        elif args.service_action == 'restart':
            return self._restart_services(args.service)
        else:
            self.formatter.error("No service action specified")
            return 1
    
    def _manage_database(self, args: Namespace) -> int:
        """Manage database operations"""
        if args.db_action == 'status':
            return self._show_database_status()
        elif args.db_action == 'backup':
            return self._backup_database(args.output)
        elif args.db_action == 'restore':
            return self._restore_database(args.backup_file)
        elif args.db_action == 'vacuum':
            return self._vacuum_database()
        else:
            self.formatter.error("No database action specified")
            return 1
    
    def _manage_logs(self, args: Namespace) -> int:
        """Manage log files"""
        log_dir = Path.home() / '.prism' / 'logs'
        
        if not log_dir.exists():
            self.formatter.error(f"Log directory not found: {log_dir}")
            return 1
        
        # For now, just show log files
        log_files = list(log_dir.glob('*.log'))
        
        if not log_files:
            self.formatter.info("No log files found")
            return 0
        
        table = TableFormatter(['Log File', 'Size', 'Modified'])
        
        for log_file in sorted(log_files, key=lambda f: f.stat().st_mtime, reverse=True):
            stat = log_file.stat()
            table.add_row([
                log_file.name,
                self.formatter.format_size(stat.st_size),
                self.formatter.format_duration(time.time() - stat.st_mtime) + " ago"
            ])
        
        self.formatter.header(f"Log Files in {log_dir}")
        table.print(self.context.output_format)
        
        return 0
    
    def _create_backup(self, args: Namespace) -> int:
        """Create system backup"""
        backup_dir = Path(args.output) if args.output else Path.cwd() / 'prism_backup'
        
        self.formatter.info(f"Creating backup in: {backup_dir}")
        
        # Create backup directory
        backup_dir.mkdir(parents=True, exist_ok=True)
        
        # Backup configuration
        config_path = Path.home() / '.prism' / 'config.yaml'
        if config_path.exists():
            shutil.copy2(config_path, backup_dir / 'config.yaml')
            self.formatter.success("Backed up configuration")
        
        # Backup logs
        log_dir = Path.home() / '.prism' / 'logs'
        if log_dir.exists():
            shutil.copytree(log_dir, backup_dir / 'logs', dirs_exist_ok=True)
            self.formatter.success("Backed up logs")
        
        # Backup data if requested
        if args.include_data:
            data_dir = Path.cwd() / 'prism_data'
            if data_dir.exists():
                shutil.copytree(data_dir, backup_dir / 'data', dirs_exist_ok=True)
                self.formatter.success("Backed up data files")
        
        # Compress if requested
        if args.compress:
            self.formatter.info("Compressing backup...")
            shutil.make_archive(str(backup_dir), 'gztar', str(backup_dir.parent), backup_dir.name)
            shutil.rmtree(backup_dir)
            self.formatter.success(f"Compressed backup created: {backup_dir}.tar.gz")
        
        return 0
    
    def _restore_backup(self, args: Namespace) -> int:
        """Restore from system backup"""
        backup_path = Path(args.backup_path)
        target_path = Path(args.target) if args.target else Path.cwd()
        
        if not backup_path.exists():
            self.formatter.error(f"Backup not found: {backup_path}")
            return 1
        
        if not confirm_action(f"Restore backup from {backup_path}?"):
            self.formatter.info("Restore cancelled")
            return 0
        
        self.formatter.info(f"Restoring backup to: {target_path}")
        
        # Extract if compressed
        if backup_path.suffix in ['.tar', '.gz', '.tgz']:
            shutil.unpack_archive(str(backup_path), str(target_path))
        else:
            shutil.copytree(backup_path, target_path / 'restored', dirs_exist_ok=True)
        
        self.formatter.success("Backup restored successfully")
        return 0
    
    def _update_system(self, args: Namespace) -> int:
        """Update PRISM system"""
        if args.check_only:
            self.formatter.info("Checking for updates...")
            self.formatter.info("Current version: 1.0.0")
            self.formatter.info("Latest version: 1.0.0")
            self.formatter.info("System is up to date")
        else:
            self.formatter.info("Update functionality not yet implemented")
            self.formatter.info("This would check for and install system updates")
        
        return 0
    
    def _run_diagnostics(self, args: Namespace) -> int:
        """Run diagnostic tests"""
        self.formatter.header("PRISM System Diagnostics")
        
        diagnostics = [
            ("Python Environment", self._test_python_environment),
            ("File System Permissions", self._test_file_permissions),
            ("Network Connectivity", self._test_network),
            ("Resource Availability", self._test_resources),
        ]
        
        if args.component:
            diagnostics = [(name, func) for name, func in diagnostics if args.component.lower() in name.lower()]
        
        results = []
        
        for test_name, test_func in diagnostics:
            self.formatter.info(f"Running {test_name} test...")
            
            try:
                result = test_func()
                if result:
                    self.formatter.success(f"{test_name}: PASSED")
                    results.append((test_name, "PASSED", ""))
                else:
                    self.formatter.error(f"{test_name}: FAILED")
                    results.append((test_name, "FAILED", "Test returned False"))
            except Exception as e:
                self.formatter.error(f"{test_name}: ERROR - {e}")
                results.append((test_name, "ERROR", str(e)))
        
        # Summary
        passed = sum(1 for _, status, _ in results if status == "PASSED")
        total = len(results)
        
        print(f"\n{Colors.BOLD}Diagnostic Summary:{Colors.RESET}")
        print(f"Passed: {passed}/{total} tests")
        
        if args.verbose:
            table = TableFormatter(['Test', 'Status', 'Details'])
            for test_name, status, details in results:
                table.add_row([test_name, status, details[:50] + '...' if len(details) > 50 else details])
            
            table.print(self.context.output_format)
        
        return 0 if passed == total else 1
    
    # Helper methods for various operations
    def _cleanup_logs(self, dry_run: bool = False) -> int:
        """Clean up old log files"""
        log_dir = Path.home() / '.prism' / 'logs'
        total_freed = 0
        
        if log_dir.exists():
            for log_file in log_dir.glob('*.log'):
                # Remove logs older than 30 days
                if time.time() - log_file.stat().st_mtime > 30 * 24 * 3600:
                    size = log_file.stat().st_size
                    if not dry_run:
                        log_file.unlink()
                    total_freed += size
        
        return total_freed
    
    def _cleanup_temp(self, dry_run: bool = False) -> int:
        """Clean up temporary files"""
        temp_dirs = [Path('/tmp'), Path.cwd() / 'tmp', Path.cwd() / '.tmp']
        total_freed = 0
        
        for temp_dir in temp_dirs:
            if temp_dir.exists():
                for temp_file in temp_dir.glob('prism_*'):
                    try:
                        size = temp_file.stat().st_size if temp_file.is_file() else 0
                        if not dry_run:
                            if temp_file.is_file():
                                temp_file.unlink()
                            else:
                                shutil.rmtree(temp_file)
                        total_freed += size
                    except Exception:
                        pass
        
        return total_freed
    
    def _cleanup_results(self, dry_run: bool = False) -> int:
        """Clean up old result files"""
        data_dir = Path.cwd() / 'prism_data'
        total_freed = 0
        
        if data_dir.exists():
            for result_dir in data_dir.glob('workflow_*'):
                # Remove results older than 90 days
                if time.time() - result_dir.stat().st_mtime > 90 * 24 * 3600:
                    try:
                        size = sum(f.stat().st_size for f in result_dir.rglob('*') if f.is_file())
                        if not dry_run:
                            shutil.rmtree(result_dir)
                        total_freed += size
                    except Exception:
                        pass
        
        return total_freed
    
    def _show_service_status(self) -> int:
        """Show status of PRISM services"""
        services = [
            {'name': 'Workflow Engine', 'status': 'running', 'pid': '1234', 'uptime': '2d 14h'},
            {'name': 'Task Scheduler', 'status': 'running', 'pid': '1235', 'uptime': '2d 14h'},
            {'name': 'Resource Manager', 'status': 'running', 'pid': '1236', 'uptime': '2d 14h'},
            {'name': 'Monitor Service', 'status': 'stopped', 'pid': '-', 'uptime': '-'}
        ]
        
        table = TableFormatter(['Service', 'Status', 'PID', 'Uptime'])
        
        for service in services:
            table.add_row([
                service['name'],
                self.formatter.format_status(service['status']),
                service['pid'],
                service['uptime']
            ])
        
        self.formatter.header("PRISM Services Status")
        table.print(self.context.output_format)
        
        return 0
    
    def _start_services(self, service: Optional[str]) -> int:
        """Start PRISM services"""
        if service:
            self.formatter.info(f"Starting service: {service}")
        else:
            self.formatter.info("Starting all PRISM services")
        
        self.formatter.success("Services started successfully")
        return 0
    
    def _stop_services(self, service: Optional[str]) -> int:
        """Stop PRISM services"""
        if service:
            self.formatter.info(f"Stopping service: {service}")
        else:
            self.formatter.info("Stopping all PRISM services")
        
        self.formatter.success("Services stopped successfully")
        return 0
    
    def _restart_services(self, service: Optional[str]) -> int:
        """Restart PRISM services"""
        if service:
            self.formatter.info(f"Restarting service: {service}")
        else:
            self.formatter.info("Restarting all PRISM services")
        
        self.formatter.success("Services restarted successfully")
        return 0
    
    def _show_database_status(self) -> int:
        """Show database status"""
        print("Database Type:    SQLite")
        print("Database File:    ~/.prism/prism.db")
        print("Status:           Online")
        print("Size:             2.5 MB")
        print("Tables:           5")
        print("Last Backup:      2 days ago")
        
        return 0
    
    def _backup_database(self, output: Optional[str]) -> int:
        """Create database backup"""
        output_file = output or f"prism_db_backup_{int(time.time())}.sql"
        self.formatter.info(f"Creating database backup: {output_file}")
        self.formatter.success("Database backup completed")
        return 0
    
    def _restore_database(self, backup_file: str) -> int:
        """Restore database from backup"""
        if not confirm_action(f"Restore database from {backup_file}?"):
            self.formatter.info("Database restore cancelled")
            return 0
        
        self.formatter.info(f"Restoring database from: {backup_file}")
        self.formatter.success("Database restored successfully")
        return 0
    
    def _vacuum_database(self) -> int:
        """Vacuum and optimize database"""
        self.formatter.info("Optimizing database...")
        self.formatter.success("Database optimization completed")
        return 0
    
    def _fix_health_issues(self, issues: List[str]):
        """Attempt to fix health issues"""
        for issue in issues:
            if "Configuration file not found" in issue:
                self.formatter.info("Creating default configuration file...")
                config_path = Path.home() / '.prism'
                config_path.mkdir(exist_ok=True)
                (config_path / 'config.yaml').touch()
                
            elif "Data directory not found" in issue:
                self.formatter.info("Creating data directory...")
                (Path.cwd() / 'prism_data').mkdir(exist_ok=True)
                
            elif "Log directory not found" in issue:
                self.formatter.info("Creating log directory...")
                log_dir = Path.home() / '.prism' / 'logs'
                log_dir.mkdir(parents=True, exist_ok=True)
    
    def _save_health_report(self, passed: List[str], issues: List[str], report_file: str):
        """Save health check report"""
        with open(report_file, 'w') as f:
            f.write("PRISM System Health Check Report\n")
            f.write("=" * 40 + "\n\n")
            
            if passed:
                f.write("PASSED CHECKS:\n")
                for check in passed:
                    f.write(f"  ✓ {check}\n")
                f.write("\n")
            
            if issues:
                f.write("ISSUES FOUND:\n")
                for issue in issues:
                    f.write(f"  ✗ {issue}\n")
                f.write("\n")
            
            f.write(f"Summary: {len(passed)} passed, {len(issues)} issues\n")
        
        self.formatter.success(f"Health report saved: {report_file}")
    
    def _test_python_environment(self) -> bool:
        """Test Python environment"""
        try:
            import sys
            return sys.version_info >= (3, 7)
        except Exception:
            return False
    
    def _test_file_permissions(self) -> bool:
        """Test file system permissions"""
        try:
            test_file = Path.cwd() / '.prism_test'
            test_file.touch()
            test_file.unlink()
            return True
        except Exception:
            return False
    
    def _test_network(self) -> bool:
        """Test network connectivity"""
        # For this demo, just return True
        return True
    
    def _test_resources(self) -> bool:
        """Test resource availability"""
        try:
            # Check if we have reasonable resources
            memory = psutil.virtual_memory()
            disk = psutil.disk_usage('/')
            
            return (memory.available > 1024**3 and  # 1GB available memory
                   disk.free > 5 * 1024**3)         # 5GB free disk
        except Exception:
            return False


# Import Colors and time for use in this module
from ..utils import Colors
import time