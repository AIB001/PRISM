#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM CLI Utilities

Utility classes for command-line interface formatting, progress tracking,
and data presentation.
"""

import sys
import time
import json
import yaml
from typing import Dict, Any, List, Optional, Union
from pathlib import Path
from dataclasses import dataclass
from datetime import datetime
import threading


class Colors:
    """Terminal color codes"""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    
    # Regular colors
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    
    # Bright colors
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    
    @classmethod
    def disable(cls):
        """Disable colors for non-terminal output"""
        cls.RESET = cls.BOLD = cls.DIM = ''
        cls.BLACK = cls.RED = cls.GREEN = cls.YELLOW = ''
        cls.BLUE = cls.MAGENTA = cls.CYAN = cls.WHITE = ''
        cls.BRIGHT_RED = cls.BRIGHT_GREEN = cls.BRIGHT_YELLOW = ''
        cls.BRIGHT_BLUE = cls.BRIGHT_MAGENTA = cls.BRIGHT_CYAN = ''


class CLIFormatter:
    """Formatter for CLI output with colors and styling"""
    
    def __init__(self, use_colors: bool = None):
        if use_colors is None:
            use_colors = sys.stdout.isatty()
        
        if not use_colors:
            Colors.disable()
    
    def success(self, message: str):
        """Print success message"""
        print(f"{Colors.BRIGHT_GREEN}✅ {message}{Colors.RESET}")
    
    def error(self, message: str):
        """Print error message"""
        print(f"{Colors.BRIGHT_RED}❌ {message}{Colors.RESET}", file=sys.stderr)
    
    def warning(self, message: str):
        """Print warning message"""
        print(f"{Colors.BRIGHT_YELLOW}⚠️  {message}{Colors.RESET}")
    
    def info(self, message: str):
        """Print info message"""
        print(f"{Colors.BRIGHT_BLUE}ℹ️  {message}{Colors.RESET}")
    
    def header(self, message: str):
        """Print header message"""
        print(f"{Colors.BOLD}{Colors.BRIGHT_CYAN}{message}{Colors.RESET}")
    
    def subheader(self, message: str):
        """Print subheader message"""
        print(f"{Colors.BOLD}{message}{Colors.RESET}")
    
    def dim(self, message: str):
        """Print dimmed message"""
        print(f"{Colors.DIM}{message}{Colors.RESET}")
    
    def highlight(self, message: str):
        """Print highlighted message"""
        print(f"{Colors.BRIGHT_WHITE}{Colors.BOLD}{message}{Colors.RESET}")
    
    def format_status(self, status: str) -> str:
        """Format status with appropriate color"""
        status_colors = {
            'running': Colors.BRIGHT_BLUE,
            'completed': Colors.BRIGHT_GREEN,
            'failed': Colors.BRIGHT_RED,
            'cancelled': Colors.YELLOW,
            'pending': Colors.CYAN,
            'active': Colors.BRIGHT_GREEN,
            'inactive': Colors.DIM,
            'warning': Colors.BRIGHT_YELLOW,
            'error': Colors.BRIGHT_RED
        }
        
        color = status_colors.get(status.lower(), Colors.WHITE)
        return f"{color}{status.upper()}{Colors.RESET}"
    
    def format_duration(self, seconds: float) -> str:
        """Format duration in human-readable format"""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            return f"{seconds/60:.1f}m"
        else:
            return f"{seconds/3600:.1f}h"
    
    def format_size(self, bytes_size: int) -> str:
        """Format file size in human-readable format"""
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if bytes_size < 1024.0:
                return f"{bytes_size:.1f}{unit}"
            bytes_size /= 1024.0
        return f"{bytes_size:.1f}PB"


class ProgressBar:
    """Progress bar for long-running operations"""
    
    def __init__(self, total: int = 100, description: str = "Progress", 
                 width: int = 50, show_percent: bool = True,
                 show_time: bool = True):
        self.total = total
        self.description = description
        self.width = width
        self.show_percent = show_percent
        self.show_time = show_time
        
        self.current = 0
        self.start_time = time.time()
        self.last_update = 0
        self.finished = False
        
        # For animation
        self.spinner_chars = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
        self.spinner_pos = 0
    
    def update(self, progress: int = None, description: str = None):
        """Update progress bar"""
        if self.finished:
            return
        
        if progress is not None:
            self.current = min(progress, self.total)
        else:
            self.current += 1
        
        if description is not None:
            self.description = description
        
        # Throttle updates
        now = time.time()
        if now - self.last_update < 0.1 and self.current < self.total:
            return
        self.last_update = now
        
        self._render()
        
        if self.current >= self.total:
            self.finish()
    
    def _render(self):
        """Render the progress bar"""
        if self.total == 0:
            # Indeterminate progress (spinner)
            spinner = self.spinner_chars[self.spinner_pos % len(self.spinner_chars)]
            self.spinner_pos += 1
            
            elapsed = time.time() - self.start_time
            time_str = f" [{self._format_time(elapsed)}]" if self.show_time else ""
            
            print(f"\r{Colors.BRIGHT_BLUE}{spinner}{Colors.RESET} {self.description}{time_str}", 
                  end="", flush=True)
        else:
            # Determinate progress (bar)
            percent = (self.current / self.total) * 100
            filled = int((self.current / self.total) * self.width)
            
            bar = "█" * filled + "░" * (self.width - filled)
            
            percent_str = f" {percent:5.1f}%" if self.show_percent else ""
            
            elapsed = time.time() - self.start_time
            if self.current > 0:
                eta = (elapsed / self.current) * (self.total - self.current)
                time_str = f" [{self._format_time(elapsed)}<{self._format_time(eta)}]" if self.show_time else ""
            else:
                time_str = f" [{self._format_time(elapsed)}]" if self.show_time else ""
            
            print(f"\r{self.description}: {Colors.BRIGHT_GREEN}{bar}{Colors.RESET}{percent_str}{time_str}",
                  end="", flush=True)
    
    def _format_time(self, seconds: float) -> str:
        """Format time duration"""
        if seconds < 60:
            return f"{int(seconds)}s"
        elif seconds < 3600:
            return f"{int(seconds//60)}m{int(seconds%60):02d}s"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            return f"{hours}h{minutes:02d}m"
    
    def finish(self):
        """Finish the progress bar"""
        if self.finished:
            return
        
        self.finished = True
        elapsed = time.time() - self.start_time
        
        if self.total == 0:
            print(f"\r{Colors.BRIGHT_GREEN}✅{Colors.RESET} {self.description} completed in {self._format_time(elapsed)}")
        else:
            print(f"\r{self.description}: {Colors.BRIGHT_GREEN}{'█' * self.width}{Colors.RESET} 100.0% [{self._format_time(elapsed)}]")
        print()


class TableFormatter:
    """Formatter for tabular data"""
    
    def __init__(self, headers: List[str], max_width: int = 120):
        self.headers = headers
        self.max_width = max_width
        self.rows = []
        self.column_widths = [len(h) for h in headers]
    
    def add_row(self, row: List[str]):
        """Add a row to the table"""
        str_row = [str(cell) for cell in row]
        self.rows.append(str_row)
        
        # Update column widths
        for i, cell in enumerate(str_row):
            if i < len(self.column_widths):
                self.column_widths[i] = max(self.column_widths[i], len(cell))
    
    def _truncate_cell(self, text: str, width: int) -> str:
        """Truncate cell content if too long"""
        if len(text) <= width:
            return text
        return text[:width-3] + "..."
    
    def render(self, output_format: str = "table") -> str:
        """Render the table in specified format"""
        if output_format == "json":
            return self._render_json()
        elif output_format == "yaml":
            return self._render_yaml()
        else:
            return self._render_table()
    
    def _render_table(self) -> str:
        """Render as ASCII table"""
        if not self.rows:
            return "No data available"
        
        # Adjust column widths to fit max_width
        available_width = self.max_width - len(self.headers) * 3 - 1
        total_width = sum(self.column_widths)
        
        if total_width > available_width:
            # Scale down proportionally
            scale_factor = available_width / total_width
            self.column_widths = [max(8, int(w * scale_factor)) for w in self.column_widths]
        
        output = []
        
        # Header
        header_line = "│ " + " │ ".join(
            self._truncate_cell(header, width).ljust(width)
            for header, width in zip(self.headers, self.column_widths)
        ) + " │"
        
        separator = "├" + "┼".join("─" * (width + 2) for width in self.column_widths) + "┤"
        top_border = "┌" + "┬".join("─" * (width + 2) for width in self.column_widths) + "┐"
        bottom_border = "└" + "┴".join("─" * (width + 2) for width in self.column_widths) + "┘"
        
        output.append(top_border)
        output.append(header_line)
        output.append(separator)
        
        # Rows
        for row in self.rows:
            row_line = "│ " + " │ ".join(
                self._truncate_cell(cell, width).ljust(width)
                for cell, width in zip(row, self.column_widths)
            ) + " │"
            output.append(row_line)
        
        output.append(bottom_border)
        return "\n".join(output)
    
    def _render_json(self) -> str:
        """Render as JSON"""
        data = []
        for row in self.rows:
            row_dict = {}
            for i, header in enumerate(self.headers):
                value = row[i] if i < len(row) else ""
                row_dict[header] = value
            data.append(row_dict)
        
        return json.dumps(data, indent=2)
    
    def _render_yaml(self) -> str:
        """Render as YAML"""
        data = []
        for row in self.rows:
            row_dict = {}
            for i, header in enumerate(self.headers):
                value = row[i] if i < len(row) else ""
                row_dict[header] = value
            data.append(row_dict)
        
        return yaml.dump(data, default_flow_style=False)
    
    def print(self, output_format: str = "table"):
        """Print the table"""
        print(self.render(output_format))


class StatusIndicator:
    """Animated status indicator for long operations"""
    
    def __init__(self, message: str = "Working"):
        self.message = message
        self.running = False
        self.thread = None
        self.spinner_chars = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
        self.position = 0
    
    def start(self):
        """Start the status indicator"""
        if self.running:
            return
        
        self.running = True
        self.thread = threading.Thread(target=self._animate, daemon=True)
        self.thread.start()
    
    def stop(self, final_message: str = None):
        """Stop the status indicator"""
        self.running = False
        if self.thread:
            self.thread.join()
        
        # Clear the line and show final message
        print(f"\r{' ' * (len(self.message) + 10)}\r", end="", flush=True)
        if final_message:
            print(final_message)
    
    def _animate(self):
        """Animation loop"""
        while self.running:
            char = self.spinner_chars[self.position % len(self.spinner_chars)]
            print(f"\r{Colors.BRIGHT_BLUE}{char}{Colors.RESET} {self.message}...", 
                  end="", flush=True)
            self.position += 1
            time.sleep(0.1)


def format_timestamp(timestamp: float) -> str:
    """Format timestamp for display"""
    dt = datetime.fromtimestamp(timestamp)
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def format_relative_time(timestamp: float) -> str:
    """Format timestamp relative to now"""
    now = time.time()
    diff = now - timestamp
    
    if diff < 60:
        return f"{int(diff)}s ago"
    elif diff < 3600:
        return f"{int(diff/60)}m ago"
    elif diff < 86400:
        return f"{int(diff/3600)}h ago"
    else:
        return f"{int(diff/86400)}d ago"


def confirm_action(message: str, default: bool = False) -> bool:
    """Ask for user confirmation"""
    suffix = " [Y/n]" if default else " [y/N]"
    
    while True:
        try:
            response = input(f"{message}{suffix}: ").strip().lower()
            
            if not response:
                return default
            elif response in ('y', 'yes'):
                return True
            elif response in ('n', 'no'):
                return False
            else:
                print("Please answer 'y' or 'n'")
        except (KeyboardInterrupt, EOFError):
            print("\nCancelled.")
            return False