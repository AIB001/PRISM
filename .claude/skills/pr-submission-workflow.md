# PR Submission Workflow

Document the complete process and requirements for PR submission and description file management.

## PR Description File Management

### Fixed Filename Rule
**Sole PR Description File**: `tests/gxf/FEP/pr_body_expanded.md`

**Strict Rules**:
- ✅ Keep: `tests/gxf/FEP/pr_body_expanded.md` (fixed name, never rename)
- ❌ Delete: All other PR-related .md files (pr_body_update.md, pr_description.md, etc.)

**Important Notes**:
- This file contains the **complete description** for PR #4 (entire FEP workflow)
- When making new commits, update this file to include **all accumulated changes**
- Use this file to update the remote PR so colleagues can see the complete change history during review
- It's fine to write a long description; the important thing is **completeness and detail**

### PR Description Cleanup
```bash
# Delete all other PR files
rm -f /tmp/pr_body_update.md /tmp/pr_description.md /tmp/prism_format_support.md

# Verify only one file remains
ls -la /tmp/pr*.md
```

## PR Submission Process

### 1. Check Uncommitted Changes
```bash
cd /data2/gxf1212/work/PRISM
git status --short
```

### 2. Stage All Changes
```bash
git add -A
git status --short
```

### 3. Commit Changes
```bash
git commit -m "$(cat <<'EOF'
[module] Brief description of changes

- Main change 1
- Main change 2
- Main change 3

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

**Commit Message Format**:
- Use `[module]` prefix (e.g., [fep], [builder], [analysis])
- Briefly describe main changes
- List key changes
- Include Co-Authored-By tag

### 4. Push to Remote
```bash
git push origin gxf
```

### 5. Update PR Description (Critical Step)
```bash
gh api repos/AIB001/PRISM/pulls/4 -X PATCH -F body=@tests/gxf/FEP/pr_body_expanded.md
```

**Important**:
- After each new push, update the PR description using `tests/gxf/FEP/pr_body_expanded.md`
- PR description should contain the **complete FEP workflow description**, not just the latest changes
- Reflect **all changes compared to main branch** from the beginning
- Help colleagues understand the complete work content during review

### 6. Verify PR Update
```bash
gh pr view 4 --json title,state,updatedAt --jq '{title, state, updated_at}'
```

## PR Description Content Requirements

### Required Sections

```markdown
## Summary
Brief overview (1-2 sentences)

## Main Features
### FEP workflow
- Complete FEP workflow features
- Atom mapping, hybrid topology, etc.

### FEP analysis
- Analysis methods, Bootstrap, etc.

### Builder and force-field integration
- Force field integration, RTF support, etc.

### Analysis and architecture refactoring
- Architecture refactoring, module splitting, etc.

## Command-line usage
Actual CLI usage examples (not test file paths)

## Python API usage
Python interface usage examples

## Current validation status
### Validated end-to-end
- Fully validated test systems

### Validated for modeling / scaffold generation
- Modeling validated test systems

### Additional validated paths
- Other validated paths

## Review guidance
Review order and recommendations

## Notes
Important notes and limitations

## Files Changed
All changed files organized by category
```

### Content Rules
1. **Complete and detailed**: Include all main features, validation results, usage methods
2. **Review-friendly**: Colleagues can see complete change content, not just latest changes
3. **Show comparison**: Display complete change history from main to now
4. **Keep updated**: Update PR description with each new commit, including all accumulated changes
5. **Real code**: Describe based on actual code, don't fabricate content
6. **Highlight focus**: FEP module is the main content, describe all FEP features in detail
7. **List features**: Detail feature points, expand on each feature category
8. **Test files**: Don't tell colleagues where test files are, just mention comprehensive test systems exist

### Must-Keep Content
- ✅ Main Features (all main feature categories)
- ✅ Command-line usage (actual CLI interface)
- ✅ Python API usage (Python interface)
- ✅ Current validation status (complete validation results)
- ✅ Review guidance (review guide)
- ✅ Files Changed (file list organized by category)
- ✅ Notes (important notes)

## Current PR Status

### PR #4: "Integrate FEP workflow with analysis and builder refactors"
- **Status**: OPEN
- **Branch**: gxf → main
- **URL**: https://github.com/AIB001/PRISM/pull/4
- **Main content**: FEP workflow (atom mapping, hybrid topology, scaffold generation, analysis, force-field integration)
- **Latest additions**: RTF force field support, SwissParam-based force fields (MMFF/MATCH/Both)

### Common Commands
```bash
# View PR status
gh pr view 4

# View PR commit count
gh pr view 4 --json commits --jq '.commits | length'

# View PR file changes
gh pr view 4 --json changedFiles --jq '.changedFiles'

# Check PR comments
gh pr view 4 --json comments --jq '.comments | length'

# Update PR description
gh api repos/AIB001/PRISM/pulls/4 -X PATCH -F body=@tests/gxf/FEP/pr_body_expanded.md
```

## Important Reminders

1. **Fixed filename**: Always use `tests/gxf/FEP/pr_body_expanded.md`, do not change filename
2. **Cleanup old files**: Delete other PR description files before each PR update
3. **Commit format**: Use `[module] description` format, don't make it too long
4. **Complete description**: PR description should be complete and detailed, include all accumulated changes
5. **Update every time**: Update PR description after each new push
6. **Show comparison**: PR description reflects complete change comparison from main to now

## Git Workflow

### Branch Management
- **Development branch**: gxf
- **Main branch**: main
- **PR target**: main

### Submission Process
1. Develop and test on gxf branch
2. Commit changes to gxf branch
3. Push to origin/gxf
4. **Update PR #4 description** (including all accumulated changes)
5. Wait for review and merge

### Don't Do
- ❌ Don't push directly to main branch
- ❌ Don't create new PRs (use PR #4)
- ❌ Don't arbitrarily change PR description filename
- ❌ Don't write too long in commit messages
- ❌ Don't describe only latest changes, include all accumulated changes
