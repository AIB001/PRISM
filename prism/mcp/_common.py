"""Shared utilities for PRISM MCP server submodules."""

import os
import sys
import logging

# ==========================================================================
#  CRITICAL: Redirect all logging to stderr.
#
#  In stdio transport mode, stdout is exclusively used for MCP JSON-RPC
#  messages. ANY non-MCP output on stdout will corrupt the protocol.
#  PRISM internally uses print() extensively for progress output, so we
#  must redirect stdout to stderr before calling any PRISM code.
# ==========================================================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("prism-mcp")


class _StdoutToStderr:
    """Context manager that redirects stdout to stderr.

    PRISM code uses print() everywhere for colored progress output.
    When running as an MCP stdio server, those prints would go to stdout
    and corrupt the JSON-RPC communication channel.
    """

    def __enter__(self):
        self._original = sys.stdout
        sys.stdout = sys.stderr
        return self

    def __exit__(self, *args):
        sys.stdout = self._original


def _ensure_prism_importable():
    """Add PRISM root to sys.path if not already present."""
    # prism/mcp/_common.py -> prism/mcp -> prism -> project root
    prism_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    if prism_root not in sys.path:
        sys.path.insert(0, prism_root)
