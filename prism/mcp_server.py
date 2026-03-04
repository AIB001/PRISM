#!/usr/bin/env python3
"""
PRISM MCP Server - AI Agent interface for protein-ligand system building.

This MCP server exposes PRISM's molecular dynamics system building capabilities
as tools that AI agents (Claude, etc.) can call through natural language.

Usage:
    # Test with MCP Inspector (browser-based debug UI)
    mcp dev prism/mcp_server.py

    # Connect to Claude Code
    claude mcp add --transport stdio prism -- python prism/mcp_server.py
"""

import sys
import os

# Fix: running this script directly sets sys.path[0] to its parent dir (prism/),
# which causes prism/mcp/ to shadow the external 'mcp' SDK package.
# Replace it with the project root so both 'prism' and 'mcp' resolve correctly.
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.dirname(_script_dir)
if sys.path and os.path.abspath(sys.path[0]) == _script_dir:
    sys.path[0] = _project_root

from prism.mcp import mcp

if __name__ == "__main__":
    mcp.run(transport="stdio")
