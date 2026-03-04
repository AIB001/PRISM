"""PRISM MCP Server package — modular AI agent interface for MD system building."""

from mcp.server.fastmcp import FastMCP

mcp = FastMCP("prism")

# Register all submodule tools/resources onto the shared mcp instance
from . import environment, build, validation, analysis, pmf_analysis, trajectory, resources

environment.register(mcp)
build.register(mcp)
validation.register(mcp)
analysis.register(mcp)
pmf_analysis.register(mcp)
trajectory.register(mcp)
resources.register(mcp)
