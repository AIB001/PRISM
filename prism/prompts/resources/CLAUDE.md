# CADD-Agent Workspace

This is a CADD (Computer-Aided Drug Design) workspace with 4 MCP servers configured:
- **chemblfind** — ChEMBL molecule search
- **MolScope** — Chemical space coverage selection
- **autodock** — AutoDock Vina blind docking
- **PRISM** — MD system building & analysis

## Auto-load Workflow

When the user asks about drug design, inhibitor screening, molecule search, target-based drug discovery, or any CADD-related task, you MUST:

1. Read the CADD workflow guide by calling: `ReadMcpResourceTool` with URI `prism://prompts/cadd_workflow`
2. Follow the workflow exactly as described in that guide

## Key Rules

- **Confirm before every action** — present parameters to the user and wait for approval before calling any tool.
- **Directory structure** — create `docking/` before docking, create `MD/` before building. All outputs go under these directories.
- **MD building defaults** — use `prism protein.pdb ligand.mol2 -lff gaff2 -ff amber14sb -o <ID> --gaussian hf --isopt false --protonation propka`. Do NOT override box, salt, or temperature parameters — always use PRISM's built-in defaults for those.
- **All file paths must be absolute.**
- **Save state** — after MD building, write `cadd_state.json` for cross-session resumption.
