#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM CLI Workflow Commands

Command handlers for workflow management, creation, execution, and monitoring.
"""

import json
import time
from pathlib import Path
from typing import Dict, Any, List, Optional
from argparse import ArgumentParser, Namespace

from ...utils.logging_system import PrismLogger
from ..utils import CLIFormatter, TableFormatter, ProgressBar, format_timestamp, format_relative_time


class WorkflowCommands:
    """Handler for workflow-related CLI commands"""
    
    def __init__(self, context):
        self.context = context
        self.logger = PrismLogger("workflow_cli")
        self.formatter = CLIFormatter()
        
        # Mock workflow storage (in real implementation, this would be a database)
        self.workflows = {}
    
    def setup_parser(self, parser: ArgumentParser):
        """Setup workflow command parser"""
        subparsers = parser.add_subparsers(dest='workflow_action', help='Workflow actions')
        
        # Create workflow
        create_parser = subparsers.add_parser('create', help='Create a new PMF workflow')
        create_parser.add_argument('system_file', help='Input system file (GRO/PDB)')
        create_parser.add_argument('--ligand', default='LIG', help='Ligand residue name')
        create_parser.add_argument('--protein', default='Protein', help='Protein selection')
        create_parser.add_argument('--name', help='Workflow name')
        create_parser.add_argument('--template', choices=['pmf', 'umbrella', 'smd'], 
                                 default='pmf', help='Workflow template')
        create_parser.add_argument('--submit', action='store_true', 
                                 help='Submit workflow immediately after creation')
        
        # List workflows
        list_parser = subparsers.add_parser('list', help='List workflows')
        list_parser.add_argument('--status', choices=['all', 'running', 'completed', 'failed'],
                               default='all', help='Filter by status')
        list_parser.add_argument('--limit', type=int, default=20, help='Maximum number to show')
        
        # Show workflow details
        show_parser = subparsers.add_parser('show', help='Show workflow details')
        show_parser.add_argument('workflow_id', help='Workflow ID')
        show_parser.add_argument('--tasks', action='store_true', help='Show task details')
        show_parser.add_argument('--logs', action='store_true', help='Show execution logs')
        
        # Submit workflow
        submit_parser = subparsers.add_parser('submit', help='Submit workflow for execution')
        submit_parser.add_argument('workflow_id', help='Workflow ID')
        submit_parser.add_argument('--priority', choices=['low', 'normal', 'high'], 
                                 default='normal', help='Execution priority')
        
        # Cancel workflow
        cancel_parser = subparsers.add_parser('cancel', help='Cancel running workflow')
        cancel_parser.add_argument('workflow_id', help='Workflow ID')
        cancel_parser.add_argument('--force', action='store_true', 
                                 help='Force cancellation without confirmation')
        
        # Resume workflow
        resume_parser = subparsers.add_parser('resume', help='Resume failed/cancelled workflow')
        resume_parser.add_argument('workflow_id', help='Workflow ID')
        resume_parser.add_argument('--from-task', help='Resume from specific task')
        
        # Delete workflow
        delete_parser = subparsers.add_parser('delete', help='Delete workflow')
        delete_parser.add_argument('workflow_id', help='Workflow ID')
        delete_parser.add_argument('--force', action='store_true',
                                 help='Delete without confirmation')
        
        # Clone workflow
        clone_parser = subparsers.add_parser('clone', help='Clone existing workflow')
        clone_parser.add_argument('workflow_id', help='Source workflow ID')
        clone_parser.add_argument('--name', help='New workflow name')
        
        # Template management
        template_parser = subparsers.add_parser('template', help='Workflow template management')
        template_subparsers = template_parser.add_subparsers(dest='template_action')
        
        template_list = template_subparsers.add_parser('list', help='List available templates')
        template_show = template_subparsers.add_parser('show', help='Show template details')
        template_show.add_argument('template_name', help='Template name')
    
    def handle(self, args: Namespace) -> int:
        """Handle workflow commands"""
        try:
            if args.workflow_action == 'create':
                return self._create_workflow(args)
            elif args.workflow_action == 'list':
                return self._list_workflows(args)
            elif args.workflow_action == 'show':
                return self._show_workflow(args)
            elif args.workflow_action == 'submit':
                return self._submit_workflow(args)
            elif args.workflow_action == 'cancel':
                return self._cancel_workflow(args)
            elif args.workflow_action == 'resume':
                return self._resume_workflow(args)
            elif args.workflow_action == 'delete':
                return self._delete_workflow(args)
            elif args.workflow_action == 'clone':
                return self._clone_workflow(args)
            elif args.workflow_action == 'template':
                return self._handle_templates(args)
            else:
                self.formatter.error("No workflow action specified")
                return 1
                
        except Exception as e:
            self.formatter.error(f"Workflow command failed: {e}")
            if self.context.verbose:
                self.logger.exception("Workflow command error")
            return 1
    
    def _create_workflow(self, args: Namespace) -> int:
        """Create a new PMF workflow"""
        system_file = Path(args.system_file)
        
        if not system_file.exists():
            self.formatter.error(f"System file not found: {system_file}")
            return 1
        
        # Generate workflow ID
        workflow_id = f"workflow_{int(time.time())}"
        if args.name:
            workflow_id = f"{args.name}_{int(time.time())}"
        
        self.formatter.info(f"Creating {args.template} workflow: {workflow_id}")
        
        # Progress simulation
        progress = ProgressBar(100, "Creating workflow")
        
        # Simulate workflow creation steps
        tasks = [
            "Validating system file",
            "Analyzing molecular structure", 
            "Setting up force field parameters",
            "Creating simulation topology",
            "Generating task dependencies",
            "Validating workflow configuration"
        ]
        
        for i, task in enumerate(tasks):
            progress.update((i + 1) * 100 // len(tasks), task)
            time.sleep(0.3)  # Simulate work
        
        progress.finish()
        
        # Store workflow (mock)
        workflow_data = {
            'id': workflow_id,
            'name': args.name or f"PMF Calculation {workflow_id}",
            'template': args.template,
            'system_file': str(system_file),
            'ligand': args.ligand,
            'protein': args.protein,
            'status': 'created',
            'created_time': time.time(),
            'tasks': self._generate_task_list(args.template),
            'estimated_runtime': 480,  # 8 hours
            'priority': 'normal'
        }
        
        self.workflows[workflow_id] = workflow_data
        
        self.formatter.success(f"Workflow created successfully: {workflow_id}")
        
        # Show workflow summary
        self._print_workflow_summary(workflow_data)
        
        # Auto-submit if requested
        if args.submit:
            self.formatter.info("Submitting workflow for execution...")
            args.workflow_id = workflow_id
            return self._submit_workflow(args)
        
        self.formatter.info(f"To submit for execution: prism workflow submit {workflow_id}")
        return 0
    
    def _list_workflows(self, args: Namespace) -> int:
        """List workflows"""
        workflows = list(self.workflows.values())
        
        # Add some mock workflows for demonstration
        if not workflows:
            workflows = self._generate_mock_workflows()
        
        # Filter by status
        if args.status != 'all':
            workflows = [w for w in workflows if w['status'] == args.status]
        
        # Sort by creation time (newest first)
        workflows.sort(key=lambda w: w.get('created_time', 0), reverse=True)
        
        # Limit results
        workflows = workflows[:args.limit]
        
        if not workflows:
            self.formatter.info(f"No workflows found with status: {args.status}")
            return 0
        
        # Create table
        table = TableFormatter([
            'ID', 'Name', 'Template', 'Status', 'Created', 'Runtime', 'Progress'
        ])
        
        for workflow in workflows:
            table.add_row([
                workflow['id'][:12] + '...' if len(workflow['id']) > 15 else workflow['id'],
                workflow['name'][:25] + '...' if len(workflow['name']) > 28 else workflow['name'],
                workflow['template'],
                self.formatter.format_status(workflow['status']),
                format_relative_time(workflow.get('created_time', time.time())),
                f"{workflow.get('estimated_runtime', 0)}m",
                self._get_progress_bar(workflow)
            ])
        
        self.formatter.header(f"Workflows ({len(workflows)} of {len(self.workflows)})")
        table.print(self.context.output_format)
        
        return 0
    
    def _show_workflow(self, args: Namespace) -> int:
        """Show detailed workflow information"""
        workflow_id = args.workflow_id
        
        if workflow_id not in self.workflows:
            # Check if it's a partial ID match
            matches = [wid for wid in self.workflows.keys() if wid.startswith(workflow_id)]
            if len(matches) == 1:
                workflow_id = matches[0]
            elif len(matches) > 1:
                self.formatter.error(f"Ambiguous workflow ID. Matches: {', '.join(matches[:3])}")
                return 1
            else:
                self.formatter.error(f"Workflow not found: {workflow_id}")
                return 1
        
        workflow = self.workflows[workflow_id]
        
        # Basic information
        self.formatter.header(f"Workflow: {workflow['name']}")
        print(f"ID:               {workflow['id']}")
        print(f"Template:         {workflow['template']}")
        print(f"Status:           {self.formatter.format_status(workflow['status'])}")
        print(f"System File:      {workflow['system_file']}")
        print(f"Ligand:           {workflow['ligand']}")
        print(f"Protein:          {workflow['protein']}")
        print(f"Created:          {format_timestamp(workflow.get('created_time', time.time()))}")
        print(f"Estimated Time:   {workflow.get('estimated_runtime', 0)} minutes")
        
        if workflow['status'] in ['running', 'completed']:
            start_time = workflow.get('start_time', workflow.get('created_time', time.time()))
            print(f"Started:          {format_timestamp(start_time)}")
            
            if workflow['status'] == 'completed':
                end_time = workflow.get('end_time', time.time())
                duration = end_time - start_time
                print(f"Completed:        {format_timestamp(end_time)}")
                print(f"Actual Runtime:   {self.formatter.format_duration(duration)}")
        
        # Task details
        if args.tasks:
            print(f"\n{Colors.BOLD}Tasks:{Colors.RESET}")
            task_table = TableFormatter(['Task', 'Status', 'Duration', 'Resources'])
            
            for task in workflow.get('tasks', []):
                task_table.add_row([
                    task['name'],
                    self.formatter.format_status(task.get('status', 'pending')),
                    self.formatter.format_duration(task.get('duration', 0)),
                    f"{task.get('cpu_cores', 1)} cores, {task.get('memory_gb', 4):.1f} GB"
                ])
            
            task_table.print(self.context.output_format)
        
        # Execution logs
        if args.logs:
            print(f"\n{Colors.BOLD}Recent Logs:{Colors.RESET}")
            logs = workflow.get('logs', [
                "Workflow created and validated",
                "Task dependencies resolved", 
                "Resource requirements calculated",
                "Workflow queued for execution"
            ])
            
            for i, log in enumerate(logs[-10:], 1):  # Last 10 logs
                timestamp = time.time() - (len(logs) - i) * 60
                print(f"  {format_timestamp(timestamp)} | {log}")
        
        return 0
    
    def _submit_workflow(self, args: Namespace) -> int:
        """Submit workflow for execution"""
        workflow_id = args.workflow_id
        
        if workflow_id not in self.workflows:
            self.formatter.error(f"Workflow not found: {workflow_id}")
            return 1
        
        workflow = self.workflows[workflow_id]
        
        if workflow['status'] in ['running', 'completed']:
            self.formatter.error(f"Workflow is already {workflow['status']}")
            return 1
        
        # Update workflow status
        workflow['status'] = 'submitted'
        workflow['priority'] = getattr(args, 'priority', 'normal')
        workflow['submit_time'] = time.time()
        
        self.formatter.success(f"Workflow submitted: {workflow_id}")
        self.formatter.info(f"Priority: {workflow['priority']}")
        self.formatter.info(f"Estimated completion time: {workflow['estimated_runtime']} minutes")
        
        # Simulate starting execution
        workflow['status'] = 'running'
        workflow['start_time'] = time.time()
        
        self.formatter.info("Workflow execution started")
        
        return 0
    
    def _cancel_workflow(self, args: Namespace) -> int:
        """Cancel running workflow"""
        workflow_id = args.workflow_id
        
        if workflow_id not in self.workflows:
            self.formatter.error(f"Workflow not found: {workflow_id}")
            return 1
        
        workflow = self.workflows[workflow_id]
        
        if workflow['status'] != 'running':
            self.formatter.error(f"Workflow is not running (status: {workflow['status']})")
            return 1
        
        # Confirmation
        if not args.force:
            from ..utils import confirm_action
            if not confirm_action(f"Cancel workflow {workflow_id}?"):
                self.formatter.info("Cancellation aborted")
                return 0
        
        # Cancel workflow
        workflow['status'] = 'cancelled'
        workflow['end_time'] = time.time()
        
        self.formatter.success(f"Workflow cancelled: {workflow_id}")
        return 0
    
    def _resume_workflow(self, args: Namespace) -> int:
        """Resume failed or cancelled workflow"""
        workflow_id = args.workflow_id
        
        if workflow_id not in self.workflows:
            self.formatter.error(f"Workflow not found: {workflow_id}")
            return 1
        
        workflow = self.workflows[workflow_id]
        
        if workflow['status'] not in ['failed', 'cancelled']:
            self.formatter.error(f"Cannot resume workflow with status: {workflow['status']}")
            return 1
        
        # Resume workflow
        workflow['status'] = 'running'
        workflow['resume_time'] = time.time()
        
        if args.from_task:
            self.formatter.info(f"Resuming from task: {args.from_task}")
        
        self.formatter.success(f"Workflow resumed: {workflow_id}")
        return 0
    
    def _delete_workflow(self, args: Namespace) -> int:
        """Delete workflow"""
        workflow_id = args.workflow_id
        
        if workflow_id not in self.workflows:
            self.formatter.error(f"Workflow not found: {workflow_id}")
            return 1
        
        workflow = self.workflows[workflow_id]
        
        # Confirmation
        if not args.force:
            from ..utils import confirm_action
            if not confirm_action(f"Delete workflow {workflow_id}? This cannot be undone."):
                self.formatter.info("Deletion aborted")
                return 0
        
        # Check if running
        if workflow['status'] == 'running':
            self.formatter.error("Cannot delete running workflow. Cancel it first.")
            return 1
        
        # Delete workflow
        del self.workflows[workflow_id]
        
        self.formatter.success(f"Workflow deleted: {workflow_id}")
        return 0
    
    def _clone_workflow(self, args: Namespace) -> int:
        """Clone existing workflow"""
        source_id = args.workflow_id
        
        if source_id not in self.workflows:
            self.formatter.error(f"Source workflow not found: {source_id}")
            return 1
        
        source = self.workflows[source_id]
        
        # Create new workflow ID
        new_id = f"workflow_{int(time.time())}"
        if args.name:
            new_id = f"{args.name}_{int(time.time())}"
        
        # Clone workflow
        new_workflow = source.copy()
        new_workflow['id'] = new_id
        new_workflow['name'] = args.name or f"{source['name']} (Clone)"
        new_workflow['status'] = 'created'
        new_workflow['created_time'] = time.time()
        
        # Remove execution-specific fields
        for field in ['start_time', 'end_time', 'submit_time', 'resume_time']:
            new_workflow.pop(field, None)
        
        self.workflows[new_id] = new_workflow
        
        self.formatter.success(f"Workflow cloned: {new_id}")
        return 0
    
    def _handle_templates(self, args: Namespace) -> int:
        """Handle template management"""
        if args.template_action == 'list':
            return self._list_templates()
        elif args.template_action == 'show':
            return self._show_template(args.template_name)
        else:
            self.formatter.error("No template action specified")
            return 1
    
    def _list_templates(self) -> int:
        """List available workflow templates"""
        templates = [
            {
                'name': 'pmf',
                'description': 'Complete PMF calculation with umbrella sampling',
                'tasks': 6,
                'estimated_time': '8-12 hours'
            },
            {
                'name': 'umbrella',
                'description': 'Umbrella sampling with pre-generated windows',
                'tasks': 4,
                'estimated_time': '4-8 hours'
            },
            {
                'name': 'smd',
                'description': 'Steered molecular dynamics only',
                'tasks': 3,
                'estimated_time': '2-4 hours'
            }
        ]
        
        table = TableFormatter(['Template', 'Description', 'Tasks', 'Est. Time'])
        
        for template in templates:
            table.add_row([
                template['name'],
                template['description'],
                str(template['tasks']),
                template['estimated_time']
            ])
        
        self.formatter.header("Available Workflow Templates")
        table.print(self.context.output_format)
        
        return 0
    
    def _show_template(self, template_name: str) -> int:
        """Show template details"""
        templates = {
            'pmf': {
                'name': 'PMF Complete',
                'description': 'Full potential of mean force calculation',
                'tasks': [
                    'System preparation and validation',
                    'Energy minimization', 
                    'Equilibration simulation',
                    'Steered molecular dynamics',
                    'Umbrella sampling windows',
                    'PMF calculation and analysis'
                ],
                'resources': '8-16 CPU cores, 16-32 GB memory',
                'estimated_time': '8-12 hours'
            }
        }
        
        if template_name not in templates:
            self.formatter.error(f"Template not found: {template_name}")
            return 1
        
        template = templates[template_name]
        
        self.formatter.header(f"Template: {template['name']}")
        print(f"Description:      {template['description']}")
        print(f"Estimated Time:   {template['estimated_time']}")
        print(f"Resources:        {template['resources']}")
        print(f"\nTasks:")
        
        for i, task in enumerate(template['tasks'], 1):
            print(f"  {i}. {task}")
        
        return 0
    
    def _generate_task_list(self, template: str) -> List[Dict[str, Any]]:
        """Generate task list for workflow template"""
        if template == 'pmf':
            return [
                {'name': 'system_prep', 'status': 'pending', 'cpu_cores': 2, 'memory_gb': 4.0},
                {'name': 'minimization', 'status': 'pending', 'cpu_cores': 4, 'memory_gb': 8.0},
                {'name': 'equilibration', 'status': 'pending', 'cpu_cores': 8, 'memory_gb': 16.0},
                {'name': 'smd_pulling', 'status': 'pending', 'cpu_cores': 16, 'memory_gb': 32.0},
                {'name': 'umbrella_sampling', 'status': 'pending', 'cpu_cores': 32, 'memory_gb': 64.0},
                {'name': 'pmf_analysis', 'status': 'pending', 'cpu_cores': 4, 'memory_gb': 8.0}
            ]
        else:
            return [
                {'name': 'system_prep', 'status': 'pending', 'cpu_cores': 2, 'memory_gb': 4.0},
                {'name': 'simulation', 'status': 'pending', 'cpu_cores': 8, 'memory_gb': 16.0},
                {'name': 'analysis', 'status': 'pending', 'cpu_cores': 4, 'memory_gb': 8.0}
            ]
    
    def _generate_mock_workflows(self) -> List[Dict[str, Any]]:
        """Generate mock workflows for demonstration"""
        workflows = []
        statuses = ['running', 'completed', 'failed', 'pending']
        
        for i in range(8):
            workflow_id = f"workflow_{1000 + i:03d}"
            status = statuses[i % len(statuses)]
            
            workflow = {
                'id': workflow_id,
                'name': f'PMF Calculation {i+1}',
                'template': 'pmf',
                'system_file': f'/data/systems/system_{i+1}.gro',
                'ligand': 'LIG',
                'protein': 'Protein',
                'status': status,
                'created_time': time.time() - (i * 3600),  # 1 hour intervals
                'estimated_runtime': 480 + (i * 30),
                'tasks': self._generate_task_list('pmf')
            }
            
            if status in ['running', 'completed']:
                workflow['start_time'] = workflow['created_time'] + 300
            if status == 'completed':
                workflow['end_time'] = workflow['start_time'] + workflow['estimated_runtime'] * 60
            
            workflows.append(workflow)
            self.workflows[workflow_id] = workflow
        
        return workflows
    
    def _print_workflow_summary(self, workflow: Dict[str, Any]):
        """Print workflow creation summary"""
        print(f"\n{Colors.BOLD}Workflow Summary:{Colors.RESET}")
        print(f"  Name:             {workflow['name']}")
        print(f"  Template:         {workflow['template']}")
        print(f"  System File:      {workflow['system_file']}")
        print(f"  Tasks:            {len(workflow['tasks'])}")
        print(f"  Estimated Time:   {workflow['estimated_runtime']} minutes")
        print(f"  Status:           {self.formatter.format_status(workflow['status'])}")
    
    def _get_progress_bar(self, workflow: Dict[str, Any]) -> str:
        """Get progress bar representation"""
        if workflow['status'] == 'completed':
            return "████████████ 100%"
        elif workflow['status'] == 'running':
            # Simulate progress based on time
            elapsed = time.time() - workflow.get('start_time', time.time())
            estimated = workflow.get('estimated_runtime', 480) * 60
            progress = min(elapsed / estimated, 0.95)  # Cap at 95% while running
            
            filled = int(progress * 12)
            bar = "█" * filled + "░" * (12 - filled)
            return f"{bar} {progress*100:.0f}%"
        elif workflow['status'] == 'failed':
            return "████░░░░░░░░ ERR"
        else:
            return "░░░░░░░░░░░░   0%"


# Import Colors for use in this module
from ..utils import Colors