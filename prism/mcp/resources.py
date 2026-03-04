"""MCP resources and prompts."""

import os


def register(mcp):

    @mcp.resource("prism://config/default")
    def get_default_config() -> str:
        """Return the default PRISM configuration YAML as a reference template.

        This shows all configurable parameters with their default values.
        Users can export this to a file with: prism --export-config my_config.yaml
        """
        config_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "configs", "default_config.yaml",
        )
        if os.path.exists(config_path):
            with open(config_path, "r") as f:
                return f.read()
        return "# default_config.yaml not found"

    def _load_agent_prompt() -> str:
        """Load the agent prompt markdown from disk."""
        prompt_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "prompts", "agent_prompt.md",
        )
        with open(prompt_path, "r") as f:
            return f.read()

    @mcp.resource("prism://prompts/agent")
    def get_agent_prompt() -> str:
        """System prompt for AI agents using PRISM tools."""
        return _load_agent_prompt()

    @mcp.prompt()
    def prism_agent() -> str:
        """PRISM assistant system prompt.

        Injects the full PRISM workflow guide into the conversation, covering:
        environment setup (GROMACS 2026.0, conda, AmberTools), recommended
        workflow (validate → build → verify), force field choices, build modes
        (MD/PMF/REST2/MMPBSA), and troubleshooting guidance.

        Use this prompt to turn the AI agent into a PRISM-aware MD simulation
        assistant that knows exactly how to use every PRISM tool.
        """
        return _load_agent_prompt()
