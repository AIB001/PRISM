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

from prism.mcp import mcp

if __name__ == "__main__":
    mcp.run(transport="stdio")
