# CADD-Agent Configuration Resources

Template files for setting up a CADD-Agent workspace with Claude Code.

## Quick Setup

### Option A: Global (recommended)

Add MCP servers to `~/.claude/settings.json` and copy `CLAUDE.md` to `~/.claude/CLAUDE.md`. This makes the CADD pipeline available in **any** directory.

### Option B: Per-project

Copy both files to your project directory:

```bash
cp CLAUDE.md /path/to/your/project/
cp .mcp.json /path/to/your/project/
```

Then edit `.mcp.json` to set the correct absolute paths for `molscope` and `prism` server scripts.

## Files

| File | Purpose |
|------|---------|
| `CLAUDE.md` | Instructs Claude Code to auto-load the CADD workflow when users ask about drug design tasks |
| `.mcp.json` | MCP server configuration template — **edit paths before use** |

## Requirements

All 4 MCP servers must be installed:
- `chemblfind` — `pip install chemblfind`
- `MolScope` — `pip install molscope`
- `autodock_mcp` — `pip install autodock-mcp`
- `PRISM` — `pip install -e .` (this repository)
