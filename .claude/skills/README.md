# Skills Documentation Guide

## Overview

The `.claude/skills/` directory contains skill documentation for the PRISM project, used to guide Claude Code in completing specific tasks. These files are included in the git repository and can be accessed at any time.

## Skills File List

### FEP Related
- **`fep-system-debug.md`** - FEP system setup and runtime debugging
- **`fep-analysis.md`** - FEP result analysis and visualization (multi-estimator)
- **`fep-forcefield-integration.md`** - FEP force field integration testing
- **`fep-mapping-viz.md`** - FEP atom mapping visualization

### General Debugging
- **`gaff-debug.md`** - GAFF force field generation debugging

## Usage Instructions

### How to Use Skills Files

1. **In Claude Code**:
   - Skills files are automatically loaded
   - When starting a conversation, Claude will reference relevant skills files
   - You can directly ask: "Use the fep-system-debug skill to..."

2. **From command line**:
   ```bash
   # View available skills
   ls .claude/skills/

   # Read a specific skill
   cat .claude/skills/fep-system-debug.md
   ```

### Test Data References

Code examples in skills files reference the `examples/` directory, but this directory is not in git. To use these examples:

#### Option 1: Use Your Own Data

```bash
# Replace example paths
mol2_a = "examples/ligands/oMeEtPh.mol2"  # Example path
mol2_a = "my_data/ligand.mol2"           # Your data
```

#### Option 2: Create Example Data Directory (Optional)

If you want to keep the complete examples, you can create an `examples/` directory locally:

```bash
# Copy example data from test directory (not committed to git)
mkdir -p examples/ligands
cp tests/gxf/FEP/unit_test/oMeEtPh-EtPh/*.mol2 examples/ligands/ 2>/dev/null || true

# Add to .gitignore (if not already there)
echo "examples/" >> .gitignore
```

## LSP Tool Notes

All skills files include LSP tool usage instructions:

- **LSP Tools Preferred**: Use `mcp__cclsp__` tools for code navigation
- **Availability**: Only available in Claude Code tool context
- **Fallback**: Use `Grep` tool when LSP is unavailable

## Code Search Tool Comparison

| Tool | Advantages | Use Cases |
|------|-----------|-----------|
| `mcp__cclsp__find_definition` | Precise symbol location | Find class/function definitions |
| `mcp__cclsp__find_references` | Find all references | Refactoring impact analysis |
| `mcp__cclsp__find_workspace_symbols` | Global symbol search | Discover related features |
| `Grep` | Text search, widely available | Quick keyword searches |

## Updating Skills Documentation

When updating skills files:

1. **Avoid hardcoded test paths**: Use `examples/` instead of `tests/`
2. **Keep code examples concise**: Show only core logic
3. **Add necessary explanations**: Explain how to replace with user's own data
4. **Update LSP notes**: Ensure all documents have LSP tool sections

## Git Tracking Status

Skills files are now tracked by git:

```bash
# Check tracking status
git status .claude/skills/

# View tracked files
git ls-files .claude/skills/
```

## Related Documentation

- **Project Guide**: `CLAUDE.md`
- **FEP Tutorial**: `/home/gxf1212/data/work/PRISM-Tutorial`
- **Test Data**: `tests/gxf/FEP/unit_test/` (local, not in git)
