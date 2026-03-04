"""MCP resources and prompts."""

import os


def register(mcp):

    _prompts_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "prompts",
    )

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

    def _load_prompt(filename: str) -> str:
        """Load a prompt markdown file from the prompts directory."""
        path = os.path.join(_prompts_dir, filename)
        with open(path, "r") as f:
            return f.read()

    @mcp.resource("prism://prompts/agent")
    def get_agent_prompt() -> str:
        """System prompt for AI agents using PRISM tools."""
        return _load_prompt("agent_prompt.md")

    @mcp.resource("prism://prompts/cadd_workflow")
    def get_cadd_workflow() -> str:
        """CADD-Agent workflow guide for multi-software drug design pipeline.

        Describes how to chain chemblfind → moljam → docking → PRISM tools
        for a complete computer-aided drug design workflow.
        """
        return _load_prompt("skills/cadd_workflow.md")

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
        return _load_prompt("agent_prompt.md")

    @mcp.prompt()
    def cadd_agent() -> str:
        """CADD-Agent system prompt for full drug design pipeline.

        Injects the CADD workflow guide that orchestrates multiple MCP servers:
        chemblfind (ChEMBL search) → moljam (data quality) → docking →
        PRISM (MD simulation & analysis). Use this when the user wants to run
        the complete computer-aided drug design pipeline.
        """
        return _load_prompt("agent_prompt.md") + "\n\n" + _load_prompt("skills/cadd_workflow.md")
