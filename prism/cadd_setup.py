#!/usr/bin/env python3
"""CADD-Agent setup: auto-detect MCP servers and configure Claude Code."""

import json
import os
import shutil
import subprocess
import sys


GITHUB_URLS = {
    "chemblfind": "https://github.com/AIB001/ChEMBLFind",
    "molscope": "https://github.com/AIB001/MolScope",
    "autodock_mcp": "https://github.com/AIB001/AutodockVina_MCP",
}

CLAUDE_MD_CONTENT = """\
# Global CADD-Agent Instructions

When the user asks about drug design, inhibitor screening, molecule search, target-based drug discovery, or any CADD-related task, you MUST:

1. Read the CADD workflow guide by calling `ReadMcpResourceTool` with URI `prism://prompts/cadd_workflow`
2. Follow the workflow exactly as described in that guide

## Key Rules

- **Confirm before every action** -- present parameters to the user and wait for approval before calling any tool.
- **Directory structure** -- create `docking/` before docking, create `MD/` before building. All outputs go under these directories.
- **MD building defaults** -- use `prism protein.pdb ligand.mol2 -lff gaff2 -ff amber14sb -o <ID> --gaussian hf --isopt false --protonation propka`. Do NOT override box, salt, or temperature parameters -- always use PRISM's built-in defaults for those.
- **All file paths must be absolute.**
- **Save state** -- after MD building, write `cadd_state.json` for cross-session resumption.
"""


def _find_python() -> str:
    """Return the absolute path to the current Python interpreter."""
    return sys.executable


def _find_command(name: str) -> str | None:
    """Find a command on PATH."""
    return shutil.which(name)


def _find_package_location(package: str) -> str | None:
    """Find a pip-installed package's location."""
    try:
        result = subprocess.run(
            [_find_python(), "-m", "pip", "show", package],
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode != 0:
            return None
        for line in result.stdout.splitlines():
            if line.startswith("Location:"):
                return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return None


def _find_module_file(package: str, module_path: str) -> str | None:
    """Find the absolute path to a module file within a package."""
    loc = _find_package_location(package)
    if not loc:
        return None
    full = os.path.join(loc, package, module_path)
    if os.path.isfile(full):
        return full
    # Try editable install (source directory)
    try:
        result = subprocess.run(
            [_find_python(), "-c", f"import {package}; import os; print(os.path.dirname({package}.__file__))"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            pkg_dir = result.stdout.strip()
            full = os.path.join(pkg_dir, module_path)
            if os.path.isfile(full):
                return full
    except Exception:
        pass
    return None


def detect_servers() -> dict:
    """Detect installed CADD MCP servers and return config dict.

    Returns dict of {server_name: config_dict} for found servers,
    and prints status for each.
    """
    python = _find_python()
    servers = {}
    missing = []

    # 1. chemblfind
    chemblfind_cmd = _find_command("chemblfind-mcp")
    if chemblfind_cmd:
        servers["chemblfind"] = {
            "type": "stdio",
            "command": chemblfind_cmd,
            "args": [],
            "timeout": 300,
        }
        print(f"  [OK] chemblfind: {chemblfind_cmd}")
    else:
        missing.append("chemblfind")
        print(f"  [MISSING] chemblfind: not found")
        print(f"            Install: pip install git+{GITHUB_URLS['chemblfind']}")

    # 2. MolScope
    molscope_server = _find_module_file("molscope", "server.py")
    if molscope_server:
        servers["molscope"] = {
            "type": "stdio",
            "command": python,
            "args": [molscope_server],
            "timeout": 120,
        }
        print(f"  [OK] molscope: {molscope_server}")
    else:
        missing.append("molscope")
        print(f"  [MISSING] molscope: not found")
        print(f"            Install: pip install git+{GITHUB_URLS['molscope']}")

    # 3. autodock
    autodock_loc = _find_package_location("autodock-mcp")
    if autodock_loc:
        servers["autodock"] = {
            "type": "stdio",
            "command": python,
            "args": ["-m", "autodock_mcp.mcp_server"],
            "timeout": 600,
        }
        print(f"  [OK] autodock: autodock_mcp.mcp_server")
    else:
        missing.append("autodock_mcp")
        print(f"  [MISSING] autodock: not found")
        print(f"            Install: pip install git+{GITHUB_URLS['autodock_mcp']}")

    # 4. PRISM (always available since we're running from it)
    prism_server = _find_module_file("prism", "mcp_server.py")
    if not prism_server:
        # Fallback: look relative to this file
        prism_server = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "mcp_server.py"
        )
    if os.path.isfile(prism_server):
        servers["prism"] = {
            "type": "stdio",
            "command": python,
            "args": [prism_server],
            "timeout": 600,
        }
        print(f"  [OK] prism: {prism_server}")
    else:
        print(f"  [WARNING] prism: mcp_server.py not found")

    return servers, missing


def setup_cadd_agent():
    """Auto-detect CADD tools, configure ~/.claude/settings.json and CLAUDE.md,
    and copy resource templates to the current directory."""

    print("\n=== CADD-Agent Setup ===\n")
    print("Detecting installed MCP servers...\n")

    servers, missing = detect_servers()

    if missing:
        print(f"\n[WARNING] {len(missing)} package(s) not found.")
        print("The CADD pipeline requires all 4 servers. Install missing packages:")
        for pkg in missing:
            print(f"  pip install git+{GITHUB_URLS.get(pkg, 'N/A')}")
        print()
        response = input("Continue with available servers? [y/N]: ").strip().lower()
        if response not in ("y", "yes"):
            print("Setup cancelled.")
            return False

    if not servers:
        print("\n[ERROR] No MCP servers found. Please install at least one package.")
        return False

    # --- Update ~/.claude/settings.json ---
    claude_dir = os.path.expanduser("~/.claude")
    os.makedirs(claude_dir, exist_ok=True)
    settings_path = os.path.join(claude_dir, "settings.json")

    if os.path.isfile(settings_path):
        with open(settings_path, "r") as f:
            settings = json.load(f)
    else:
        settings = {}

    existing_mcp = settings.get("mcpServers", {})
    existing_mcp.update(servers)
    settings["mcpServers"] = existing_mcp

    with open(settings_path, "w") as f:
        json.dump(settings, f, indent=2)
    print(f"\n[OK] Updated {settings_path}")
    print(f"     Added {len(servers)} MCP server(s): {', '.join(servers.keys())}")

    # --- Write ~/.claude/CLAUDE.md ---
    claude_md_path = os.path.join(claude_dir, "CLAUDE.md")
    if os.path.isfile(claude_md_path):
        with open(claude_md_path, "r") as f:
            existing = f.read()
        if "CADD-Agent" in existing:
            print(f"[OK] {claude_md_path} already contains CADD-Agent instructions (skipped)")
        else:
            with open(claude_md_path, "a") as f:
                f.write("\n\n" + CLAUDE_MD_CONTENT)
            print(f"[OK] Appended CADD-Agent instructions to {claude_md_path}")
    else:
        with open(claude_md_path, "w") as f:
            f.write(CLAUDE_MD_CONTENT)
        print(f"[OK] Created {claude_md_path}")

    # --- Copy resource templates to current directory ---
    resources_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "prompts", "resources",
    )
    cwd = os.getcwd()

    # Copy CLAUDE.md template (skip if exists)
    claude_md_src = os.path.join(resources_dir, "CLAUDE.md")
    claude_md_dst = os.path.join(cwd, "CLAUDE.md")
    if os.path.isfile(claude_md_src):
        if os.path.isfile(claude_md_dst):
            print(f"[SKIP] {claude_md_dst} already exists")
        else:
            shutil.copy2(claude_md_src, claude_md_dst)
            print(f"[OK] Copied CLAUDE.md to {cwd}/")

    # Write .mcp.json with actual detected paths (only if not exists)
    mcp_json_dst = os.path.join(cwd, ".mcp.json")
    if os.path.isfile(mcp_json_dst):
        print(f"[SKIP] {mcp_json_dst} already exists")
    else:
        mcp_config = {"mcpServers": servers}
        with open(mcp_json_dst, "w") as f:
            json.dump(mcp_config, f, indent=2)
        print(f"[OK] Wrote {mcp_json_dst} with detected server paths")

    print("\n=== Setup Complete ===\n")
    print("You can now start Claude Code in any directory and use the CADD pipeline.")
    print("Try: claude  ->  'I want to screen inhibitors for CDK7'\n")

    return True


if __name__ == "__main__":
    setup_cadd_agent()
