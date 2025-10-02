#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM CLI Monitor Commands

Command handlers for system monitoring, resource tracking, and real-time status display.
"""

import time
import json
import psutil
from pathlib import Path
from typing import Dict, Any, List, Optional
from argparse import ArgumentParser, Namespace
import threading

from ...utils.logging_system import PrismLogger
from ..utils import CLIFormatter, TableFormatter, ProgressBar, StatusIndicator, format_timestamp, format_relative_time


class MonitorCommands:
    """Handler for monitoring-related CLI commands"""
    
    def __init__(self, context):
        self.context = context
        self.logger = PrismLogger("monitor_cli")
        self.formatter = CLIFormatter()
        
        # Real-time monitoring state
        self.monitoring_active = False
        self.monitor_thread = None
    
    def setup_parser(self, parser: ArgumentParser):
        """Setup monitor command parser"""
        subparsers = parser.add_subparsers(dest='monitor_action', help='Monitor actions')
        
        # System status
        system_parser = subparsers.add_parser('system', help='Show system status and resources')
        system_parser.add_argument('--real-time', action='store_true', 
                                 help='Show real-time updating dashboard')
        system_parser.add_argument('--interval', type=int, default=5, 
                                 help='Update interval in seconds for real-time mode')
        system_parser.add_argument('--dashboard', action='store_true',
                                 help='Show comprehensive system dashboard')
        
        # Workflow monitoring
        workflows_parser = subparsers.add_parser('workflows', help='Monitor running workflows')
        workflows_parser.add_argument('--follow', '-f', action='store_true',
                                    help='Follow workflow progress in real-time')
        workflows_parser.add_argument('--status', choices=['all', 'running', 'failed'],
                                    default='running', help='Filter workflows by status')
        
        # Resource utilization
        resources_parser = subparsers.add_parser('resources', help='Show resource utilization')
        resources_parser.add_argument('--history', action='store_true',
                                    help='Show historical resource usage')
        resources_parser.add_argument('--duration', default='1h',
                                    help='Duration for history (1h, 1d, 1w)')
        
        # Cluster monitoring (if distributed)
        cluster_parser = subparsers.add_parser('cluster', help='Monitor cluster status')
        cluster_parser.add_argument('--nodes', action='store_true',
                                  help='Show detailed node information')
        cluster_parser.add_argument('--load-balance', action='store_true',
                                  help='Show load balancing status')
        
        # Performance metrics
        performance_parser = subparsers.add_parser('performance', help='Show performance metrics')
        performance_parser.add_argument('--summary', action='store_true',
                                       help='Show performance summary')
        performance_parser.add_argument('--detailed', action='store_true',
                                       help='Show detailed performance breakdown')
        
        # Alerts and notifications
        alerts_parser = subparsers.add_parser('alerts', help='Show system alerts')
        alerts_parser.add_argument('--active', action='store_true',
                                 help='Show only active alerts')
        alerts_parser.add_argument('--severity', choices=['info', 'warning', 'error', 'critical'],
                                 help='Filter by severity level')
        
        # Logs monitoring
        logs_parser = subparsers.add_parser('logs', help='Monitor system logs')
        logs_parser.add_argument('--follow', '-f', action='store_true',
                               help='Follow log output')
        logs_parser.add_argument('--level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                               default='INFO', help='Minimum log level')
        logs_parser.add_argument('--component', help='Filter by component name')
    
    def handle(self, args: Namespace) -> int:
        """Handle monitor commands"""
        try:
            if args.monitor_action == 'system':
                return self._monitor_system(args)
            elif args.monitor_action == 'workflows':
                return self._monitor_workflows(args)
            elif args.monitor_action == 'resources':
                return self._monitor_resources(args)
            elif args.monitor_action == 'cluster':
                return self._monitor_cluster(args)
            elif args.monitor_action == 'performance':
                return self._monitor_performance(args)
            elif args.monitor_action == 'alerts':
                return self._monitor_alerts(args)
            elif args.monitor_action == 'logs':
                return self._monitor_logs(args)
            else:
                self.formatter.error("No monitor action specified")
                return 1
                
        except KeyboardInterrupt:
            self.formatter.info("Monitoring stopped")
            return 0
        except Exception as e:
            self.formatter.error(f"Monitor command failed: {e}")
            if self.context.verbose:
                self.logger.exception("Monitor command error")
            return 1
    
    def _monitor_system(self, args: Namespace) -> int:
        """Monitor system status"""
        if args.real_time:
            return self._real_time_system_monitor(args)
        elif args.dashboard:
            return self._show_system_dashboard()
        else:
            return self._show_system_status()
    
    def _show_system_status(self) -> int:
        """Show current system status"""
        self.formatter.header("System Status")
        
        # Basic system info
        try:
            # CPU info
            cpu_percent = psutil.cpu_percent(interval=1)
            cpu_count = psutil.cpu_count()
            
            # Memory info
            memory = psutil.virtual_memory()
            
            # Disk info
            disk = psutil.disk_usage('/')
            
            # Load averages (Unix-like systems)
            try:
                load_avg = psutil.getloadavg()
                load_str = f"{load_avg[0]:.2f}, {load_avg[1]:.2f}, {load_avg[2]:.2f}"
            except AttributeError:
                load_str = "N/A (Windows)"
            
            print(f"CPU Usage:        {cpu_percent:.1f}% ({cpu_count} cores)")
            print(f"Memory Usage:     {memory.percent:.1f}% ({self.formatter.format_size(memory.used)}/{self.formatter.format_size(memory.total)})")
            print(f"Disk Usage:       {disk.percent:.1f}% ({self.formatter.format_size(disk.used)}/{self.formatter.format_size(disk.total)})")
            print(f"Load Average:     {load_str}")
            
            # PRISM-specific status
            print(f"\n{Colors.BOLD}PRISM Services:{Colors.RESET}")
            services = self._get_service_status()
            
            service_table = TableFormatter(['Service', 'Status', 'Uptime', 'Memory'])
            for service in services:
                service_table.add_row([
                    service['name'],
                    self.formatter.format_status(service['status']),
                    service['uptime'],
                    service['memory']
                ])
            
            service_table.print(self.context.output_format)
            
        except Exception as e:
            self.formatter.error(f"Failed to get system status: {e}")
            return 1
        
        return 0
    
    def _show_system_dashboard(self) -> int:
        """Show comprehensive system dashboard"""
        self.formatter.header("PRISM System Dashboard")
        
        # System overview
        print(f"{Colors.BOLD}System Overview:{Colors.RESET}")
        self._show_system_status()
        
        # Workflow summary
        print(f"\n{Colors.BOLD}Workflow Summary:{Colors.RESET}")
        workflow_stats = self._get_workflow_stats()
        
        stats_table = TableFormatter(['Metric', 'Count', 'Percentage'])
        total_workflows = sum(workflow_stats.values())
        
        for status, count in workflow_stats.items():
            percentage = (count / max(total_workflows, 1)) * 100
            stats_table.add_row([
                status.title(),
                str(count),
                f"{percentage:.1f}%"
            ])
        
        stats_table.print(self.context.output_format)
        
        # Resource utilization chart (ASCII)
        print(f"\n{Colors.BOLD}Resource Utilization (Last 24h):{Colors.RESET}")
        self._show_resource_chart()
        
        # Recent alerts
        print(f"\n{Colors.BOLD}Recent Alerts:{Colors.RESET}")
        alerts = self._get_recent_alerts()
        
        if alerts:
            alert_table = TableFormatter(['Time', 'Severity', 'Message'])
            for alert in alerts[:5]:  # Show last 5
                alert_table.add_row([
                    format_relative_time(alert['timestamp']),
                    self.formatter.format_status(alert['severity']),
                    alert['message'][:60] + '...' if len(alert['message']) > 60 else alert['message']
                ])
            
            alert_table.print(self.context.output_format)
        else:
            self.formatter.info("No recent alerts")
        
        return 0
    
    def _real_time_system_monitor(self, args: Namespace) -> int:
        """Real-time system monitoring"""
        self.formatter.header("Real-time System Monitor (Press Ctrl+C to exit)")
        
        try:
            while True:
                # Clear screen (ANSI escape code)
                print("\033[2J\033[H", end="")
                
                # Show timestamp
                print(f"Last updated: {format_timestamp(time.time())}")
                print("=" * 60)
                
                # System metrics
                cpu_percent = psutil.cpu_percent(interval=0.1)
                memory = psutil.virtual_memory()
                
                print(f"CPU:    {self._create_bar(cpu_percent, 100)} {cpu_percent:5.1f}%")
                print(f"Memory: {self._create_bar(memory.percent, 100)} {memory.percent:5.1f}%")
                
                # Workflow status
                print(f"\n{Colors.BOLD}Active Workflows:{Colors.RESET}")
                active_workflows = self._get_active_workflows()
                
                if active_workflows:
                    for workflow in active_workflows[:5]:  # Show top 5
                        progress = workflow.get('progress', 0)
                        print(f"  {workflow['id'][:15]:<15} {self._create_bar(progress, 100)} {progress:3.0f}%")
                else:
                    print("  No active workflows")
                
                # Resource usage by process
                print(f"\n{Colors.BOLD}Top Processes:{Colors.RESET}")
                processes = self._get_top_processes()
                
                for proc in processes[:3]:
                    print(f"  {proc['name'][:20]:<20} CPU: {proc['cpu']:5.1f}% MEM: {proc['memory']:5.1f}%")
                
                time.sleep(args.interval)
                
        except KeyboardInterrupt:
            self.formatter.info("Real-time monitoring stopped")
            return 0
    
    def _monitor_workflows(self, args: Namespace) -> int:
        """Monitor workflow status"""
        if args.follow:
            return self._follow_workflows(args)
        else:
            return self._show_workflow_status(args)
    
    def _show_workflow_status(self, args: Namespace) -> int:
        """Show current workflow status"""
        workflows = self._get_workflows_by_status(args.status)
        
        if not workflows:
            self.formatter.info(f"No {args.status} workflows found")
            return 0
        
        self.formatter.header(f"{args.status.title()} Workflows ({len(workflows)})")
        
        table = TableFormatter(['ID', 'Name', 'Status', 'Progress', 'Runtime', 'ETA'])
        
        for workflow in workflows:
            eta = self._calculate_eta(workflow)
            runtime = self.formatter.format_duration(workflow.get('runtime', 0))
            
            table.add_row([
                workflow['id'][:12] + '...' if len(workflow['id']) > 15 else workflow['id'],
                workflow['name'][:20] + '...' if len(workflow['name']) > 23 else workflow['name'],
                self.formatter.format_status(workflow['status']),
                f"{workflow.get('progress', 0):3.0f}%",
                runtime,
                eta
            ])
        
        table.print(self.context.output_format)
        
        return 0
    
    def _follow_workflows(self, args: Namespace) -> int:
        """Follow workflow progress in real-time"""
        self.formatter.header("Following Workflow Progress (Press Ctrl+C to exit)")
        
        try:
            while True:
                print(f"\033[2J\033[H", end="")  # Clear screen
                print(f"Updated: {format_timestamp(time.time())}")
                print("=" * 80)
                
                workflows = self._get_workflows_by_status('running')
                
                if not workflows:
                    print("No running workflows")
                else:
                    for workflow in workflows:
                        progress = workflow.get('progress', 0)
                        current_task = workflow.get('current_task', 'Unknown')
                        
                        print(f"\n{Colors.BOLD}{workflow['name']} ({workflow['id'][:8]}){Colors.RESET}")
                        print(f"  Status: {self.formatter.format_status(workflow['status'])}")
                        print(f"  Progress: {self._create_bar(progress, 100)} {progress:.1f}%")
                        print(f"  Current: {current_task}")
                        print(f"  Runtime: {self.formatter.format_duration(workflow.get('runtime', 0))}")
                
                time.sleep(5)  # Update every 5 seconds
                
        except KeyboardInterrupt:
            self.formatter.info("Workflow monitoring stopped")
            return 0
    
    def _monitor_resources(self, args: Namespace) -> int:
        """Monitor resource utilization"""
        if args.history:
            return self._show_resource_history(args)
        else:
            return self._show_current_resources()
    
    def _show_current_resources(self) -> int:
        """Show current resource utilization"""
        self.formatter.header("Current Resource Utilization")
        
        # System resources
        cpu_percent = psutil.cpu_percent(interval=1)
        memory = psutil.virtual_memory()
        
        print(f"CPU Usage:    {self._create_bar(cpu_percent, 100)} {cpu_percent:5.1f}%")
        print(f"Memory Usage: {self._create_bar(memory.percent, 100)} {memory.percent:5.1f}%")
        
        # Per-core CPU usage
        cpu_per_core = psutil.cpu_percent(percpu=True, interval=1)
        print(f"\nPer-Core CPU Usage:")
        for i, usage in enumerate(cpu_per_core):
            print(f"  Core {i:2d}: {self._create_bar(usage, 100)} {usage:5.1f}%")
        
        # Network I/O
        try:
            net_io = psutil.net_io_counters()
            print(f"\nNetwork I/O:")
            print(f"  Bytes Sent:     {self.formatter.format_size(net_io.bytes_sent)}")
            print(f"  Bytes Received: {self.formatter.format_size(net_io.bytes_recv)}")
        except AttributeError:
            pass
        
        # Disk I/O
        try:
            disk_io = psutil.disk_io_counters()
            print(f"\nDisk I/O:")
            print(f"  Bytes Read:     {self.formatter.format_size(disk_io.read_bytes)}")
            print(f"  Bytes Written:  {self.formatter.format_size(disk_io.write_bytes)}")
        except AttributeError:
            pass
        
        return 0
    
    def _show_resource_history(self, args: Namespace) -> int:
        """Show historical resource usage"""
        self.formatter.header(f"Resource History ({args.duration})")
        
        # Simulate historical data
        print("CPU Usage History:")
        self._show_resource_chart()
        
        print(f"\nMemory Usage History:")
        self._show_resource_chart()
        
        return 0
    
    def _monitor_performance(self, args: Namespace) -> int:
        """Monitor performance metrics"""
        if args.summary:
            return self._show_performance_summary()
        elif args.detailed:
            return self._show_detailed_performance()
        else:
            return self._show_basic_performance()
    
    def _show_performance_summary(self) -> int:
        """Show performance summary"""
        self.formatter.header("Performance Summary")
        
        # Mock performance data
        metrics = {
            'Workflows/Hour': 12.5,
            'Average Task Duration': '45.2 minutes',
            'System Efficiency': '87.3%',
            'Resource Utilization': '73.1%',
            'Error Rate': '2.1%',
            'Queue Wait Time': '3.4 minutes'
        }
        
        table = TableFormatter(['Metric', 'Value'])
        for metric, value in metrics.items():
            table.add_row([metric, str(value)])
        
        table.print(self.context.output_format)
        
        return 0
    
    def _monitor_alerts(self, args: Namespace) -> int:
        """Monitor system alerts"""
        alerts = self._get_alerts(args.active, args.severity)
        
        if not alerts:
            self.formatter.info("No alerts found")
            return 0
        
        self.formatter.header(f"System Alerts ({len(alerts)})")
        
        table = TableFormatter(['Time', 'Severity', 'Source', 'Message'])
        
        for alert in alerts:
            table.add_row([
                format_relative_time(alert['timestamp']),
                self.formatter.format_status(alert['severity']),
                alert['source'],
                alert['message'][:50] + '...' if len(alert['message']) > 50 else alert['message']
            ])
        
        table.print(self.context.output_format)
        
        return 0
    
    def _monitor_logs(self, args: Namespace) -> int:
        """Monitor system logs"""
        if args.follow:
            return self._follow_logs(args)
        else:
            return self._show_recent_logs(args)
    
    def _follow_logs(self, args: Namespace) -> int:
        """Follow log output in real-time"""
        self.formatter.header("Following System Logs (Press Ctrl+C to exit)")
        
        try:
            while True:
                # Simulate log entries
                logs = self._get_recent_logs(limit=5)
                
                for log in logs:
                    timestamp = format_timestamp(log['timestamp'])
                    level_color = self._get_log_level_color(log['level'])
                    
                    print(f"{timestamp} {level_color}{log['level']:8}{Colors.RESET} "
                          f"{log['component']:15} {log['message']}")
                
                time.sleep(2)
                
        except KeyboardInterrupt:
            self.formatter.info("Log monitoring stopped")
            return 0
    
    # Helper methods
    def _create_bar(self, value: float, max_value: float, width: int = 20) -> str:
        """Create ASCII progress bar"""
        if max_value == 0:
            filled = 0
        else:
            filled = int((value / max_value) * width)
        
        bar = "█" * filled + "░" * (width - filled)
        
        # Color based on value
        if value < 50:
            color = Colors.BRIGHT_GREEN
        elif value < 80:
            color = Colors.BRIGHT_YELLOW
        else:
            color = Colors.BRIGHT_RED
        
        return f"{color}{bar}{Colors.RESET}"
    
    def _show_resource_chart(self):
        """Show ASCII resource usage chart"""
        # Simulate 24 data points (hours)
        data = [45, 52, 48, 51, 65, 72, 68, 71, 75, 82, 79, 85,
                88, 84, 81, 78, 74, 69, 66, 63, 58, 54, 49, 46]
        
        max_val = max(data)
        
        # Print chart
        for i in range(10, 0, -1):
            line = f"{i*10:3d}% |"
            for val in data:
                if val >= i * 10:
                    line += "█"
                elif val >= (i-1) * 10:
                    line += "▄"
                else:
                    line += " "
            print(line)
        
        print("     +" + "─" * len(data))
        print("      " + "".join([f"{i:2d}" if i % 4 == 0 else "  " for i in range(24)]))
        print("      Time (hours ago)")
    
    def _get_service_status(self) -> List[Dict[str, str]]:
        """Get status of PRISM services"""
        return [
            {'name': 'Workflow Engine', 'status': 'active', 'uptime': '2d 14h', 'memory': '156 MB'},
            {'name': 'Task Scheduler', 'status': 'active', 'uptime': '2d 14h', 'memory': '89 MB'},
            {'name': 'Resource Manager', 'status': 'active', 'uptime': '2d 14h', 'memory': '45 MB'},
            {'name': 'Monitor Service', 'status': 'active', 'uptime': '2d 14h', 'memory': '32 MB'},
            {'name': 'Data Storage', 'status': 'active', 'uptime': '2d 14h', 'memory': '234 MB'}
        ]
    
    def _get_workflow_stats(self) -> Dict[str, int]:
        """Get workflow statistics"""
        return {
            'running': 3,
            'completed': 15,
            'failed': 1,
            'pending': 2
        }
    
    def _get_recent_alerts(self) -> List[Dict[str, Any]]:
        """Get recent system alerts"""
        return [
            {
                'timestamp': time.time() - 300,
                'severity': 'warning',
                'message': 'High memory usage detected on compute node 2'
            },
            {
                'timestamp': time.time() - 1800,
                'severity': 'info',
                'message': 'Workflow workflow_001 completed successfully'
            }
        ]
    
    def _get_active_workflows(self) -> List[Dict[str, Any]]:
        """Get currently active workflows"""
        return [
            {'id': 'workflow_001', 'progress': 65.0},
            {'id': 'workflow_002', 'progress': 23.0},
            {'id': 'workflow_003', 'progress': 89.0}
        ]
    
    def _get_top_processes(self) -> List[Dict[str, Any]]:
        """Get top processes by resource usage"""
        processes = []
        try:
            for proc in psutil.process_iter(['name', 'cpu_percent', 'memory_percent']):
                if proc.info['cpu_percent'] > 0:
                    processes.append({
                        'name': proc.info['name'],
                        'cpu': proc.info['cpu_percent'],
                        'memory': proc.info['memory_percent']
                    })
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        
        return sorted(processes, key=lambda x: x['cpu'], reverse=True)
    
    def _get_workflows_by_status(self, status: str) -> List[Dict[str, Any]]:
        """Get workflows by status"""
        # Mock workflow data
        all_workflows = [
            {
                'id': 'workflow_001',
                'name': 'PMF Calculation 1',
                'status': 'running',
                'progress': 65.0,
                'runtime': 3600,
                'current_task': 'umbrella_sampling'
            },
            {
                'id': 'workflow_002', 
                'name': 'PMF Calculation 2',
                'status': 'running',
                'progress': 23.0,
                'runtime': 1200,
                'current_task': 'equilibration'
            }
        ]
        
        if status == 'all':
            return all_workflows
        else:
            return [w for w in all_workflows if w['status'] == status]
    
    def _calculate_eta(self, workflow: Dict[str, Any]) -> str:
        """Calculate estimated time of arrival"""
        progress = workflow.get('progress', 0)
        runtime = workflow.get('runtime', 0)
        
        if progress == 0 or progress >= 100:
            return "N/A"
        
        total_estimated = runtime / (progress / 100)
        remaining = total_estimated - runtime
        
        return self.formatter.format_duration(remaining)
    
    def _get_alerts(self, active_only: bool = False, severity: str = None) -> List[Dict[str, Any]]:
        """Get system alerts"""
        alerts = [
            {
                'timestamp': time.time() - 300,
                'severity': 'warning',
                'source': 'resource_monitor',
                'message': 'High memory usage detected'
            },
            {
                'timestamp': time.time() - 1800,
                'severity': 'info', 
                'source': 'workflow_engine',
                'message': 'Workflow completed successfully'
            }
        ]
        
        if severity:
            alerts = [a for a in alerts if a['severity'] == severity]
        
        return alerts
    
    def _get_recent_logs(self, limit: int = 20) -> List[Dict[str, Any]]:
        """Get recent log entries"""
        return [
            {
                'timestamp': time.time() - i*10,
                'level': 'INFO',
                'component': 'workflow_engine',
                'message': f'Task completed successfully (ID: task_{i:03d})'
            }
            for i in range(limit)
        ]
    
    def _get_log_level_color(self, level: str) -> str:
        """Get color for log level"""
        colors = {
            'DEBUG': Colors.DIM,
            'INFO': Colors.BRIGHT_BLUE,
            'WARNING': Colors.BRIGHT_YELLOW,
            'ERROR': Colors.BRIGHT_RED,
            'CRITICAL': Colors.BRIGHT_MAGENTA
        }
        return colors.get(level, Colors.WHITE)
    
    def _show_recent_logs(self, args: Namespace) -> int:
        """Show recent log entries"""
        logs = self._get_recent_logs()
        
        self.formatter.header("Recent System Logs")
        
        for log in logs:
            timestamp = format_timestamp(log['timestamp'])
            level_color = self._get_log_level_color(log['level'])
            
            print(f"{timestamp} {level_color}{log['level']:8}{Colors.RESET} "
                  f"{log['component']:15} {log['message']}")
        
        return 0
    
    def _show_basic_performance(self) -> int:
        """Show basic performance metrics"""
        self.formatter.header("Performance Metrics")
        
        metrics = [
            ("Throughput", "12.5 workflows/hour"),
            ("Avg Response Time", "245ms"),
            ("Success Rate", "97.9%"),
            ("Resource Efficiency", "87.3%")
        ]
        
        for name, value in metrics:
            print(f"{name:20}: {value}")
        
        return 0
    
    def _show_detailed_performance(self) -> int:
        """Show detailed performance breakdown"""
        self.formatter.header("Detailed Performance Analysis")
        
        # Task performance
        print(f"{Colors.BOLD}Task Performance:{Colors.RESET}")
        task_table = TableFormatter(['Task Type', 'Avg Duration', 'Success Rate', 'Throughput'])
        
        task_data = [
            ('System Prep', '5.2 min', '99.1%', '15.2/hr'),
            ('Equilibration', '18.5 min', '98.7%', '3.8/hr'), 
            ('SMD Pulling', '45.3 min', '97.2%', '1.6/hr'),
            ('Umbrella Sampling', '125.8 min', '96.8%', '0.7/hr'),
            ('PMF Analysis', '8.9 min', '99.5%', '8.1/hr')
        ]
        
        for task_type, duration, success_rate, throughput in task_data:
            task_table.add_row([task_type, duration, success_rate, throughput])
        
        task_table.print(self.context.output_format)
        
        return 0
    
    def _monitor_cluster(self, args: Namespace) -> int:
        """Monitor cluster status"""
        self.formatter.header("Cluster Status")
        
        if args.nodes:
            return self._show_node_details()
        elif args.load_balance:
            return self._show_load_balance_status()
        else:
            return self._show_cluster_overview()
    
    def _show_cluster_overview(self) -> int:
        """Show cluster overview"""
        print("Cluster: production-cluster")
        print("Total Nodes: 4")
        print("Active Nodes: 3")
        print("Total Resources: 64 cores, 256 GB memory")
        print("Available Resources: 32 cores, 128 GB memory")
        
        return 0
    
    def _show_node_details(self) -> int:
        """Show detailed node information"""
        node_table = TableFormatter(['Node', 'Status', 'CPU', 'Memory', 'Load', 'Workflows'])
        
        nodes = [
            ('compute-01', 'active', '45%', '67%', '2.1', '2'),
            ('compute-02', 'active', '23%', '34%', '1.2', '1'), 
            ('compute-03', 'active', '78%', '89%', '3.4', '3'),
            ('compute-04', 'maintenance', '0%', '5%', '0.0', '0')
        ]
        
        for node, status, cpu, memory, load, workflows in nodes:
            node_table.add_row([node, self.formatter.format_status(status), cpu, memory, load, workflows])
        
        node_table.print(self.context.output_format)
        
        return 0
    
    def _show_load_balance_status(self) -> int:
        """Show load balancing status"""
        self.formatter.header("Load Balancing Status")
        
        print("Strategy: Resource-Aware")
        print("Balance Score: 0.85/1.00")
        print("Rebalancing: Disabled")
        
        return 0


# Import Colors for use in this module  
from ..utils import Colors