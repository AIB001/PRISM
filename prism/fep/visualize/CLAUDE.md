# PRISM-FEP Visualization Module

This file provides guidance for working with FEP visualization (HTML/PNG generation and rendering).

## Directory Naming Rule (CRITICAL)

- Use simple case names for force-field combinations: `<protein_ff>-mut_<ligand_ff>` (example: `charmm36m-mut_mmff`).
- Never create nested duplicate system directories like `.../GMX_PROLIG_FEP/GMX_PROLIG_FEP/`.
- Avoid ad-hoc suffixes (`_pkgfix*`, `_final*`, `_new*`) in directory names.
- Default output directory for FEP cases should follow `<protein_ff>-mut_<ligand_ff>` when users do not pass an explicit output path.
- For 42-38 maintenance, move legacy non-canonical outputs into `tests/gxf/FEP/unit_test/42-38/Archive/`.

## Module Overview

The visualization module generates HTML and PNG visualizations of atom mapping results for quality control:
- **HTML Canvas**: Interactive 2D visualization with zoom/pan
- **Bond order rendering**: Professional aromatic ring display
- **Warning banner**: Alerts for unclassified atoms
- **Multi-force field support**: GAFF, CGenFF, OpenFF, OPLS-AA

## Critical Quality Constraints

### HTML Generation Quality Standards

**Mandatory Requirements**:
- ✅ **Zero gray atoms**: All atoms must be classified (no `rgb(200, 200, 200)`)
- ✅ **Correct atom labels**: Show real names (not "Atom15" placeholders)
- ✅ **Consistent common counts**: Left and right molecules must have same common atom count
- ✅ **Element validation**: Common atoms must have same element on both sides
- ✅ **Chemical sanity**: No carbons at chain ends, no multi-bond hydrogens

**Verification Commands**:
```bash
# Check for unclassified atoms
grep '"classification": "unknown"' mapping.html | wc -l  # Should be 0

# Check for gray atoms
grep "rgb(200, 200, 200)" mapping.html | wc -l  # Should be 0

# Check for placeholder atom names
grep "Atom[0-9]+" mapping.html | wc -l  # Should be 0
```

**Forbidden Practices**:
- ❌ Using purple or flashy backgrounds
- ❌ Deviating from mapping module style
- ❌ Debugging via screenshots (must use log analysis)
- ❌ Displaying Chinese punctuation (use English only)

### Frontend Layout Requirements

**Control Panel Layout**:
- ✅ Coloring Mode buttons (FEP Classification / Element)
- ✅ Legend (mode-specific, not duplicated)
- ✅ Controls section (Drag to pan, Scroll to zoom, Hover for atom info)
- ✅ Working Directory display (not highlighted)
- ✅ Total Charge display (4 decimal places precision)

**YAML Parameters Display**:
- ✅ Collapsible section showing FEP-specific parameters
- ✅ Non-default values highlighted
- ✅ Tooltips for important parameters
- ✅ Parameters: dist_cutoff, charge_cutoff, force field type

**Layout Constraints**:
- Coloring Mode: Single-select radio buttons (not two independent button sets)
- Legend width: Consistent with other panels, not content-stretched
- Legend: Centered or aligned with controls
- Avoid extra whitespace in Legend area
- No奇怪换行 in English text

### Bond Order Rendering Requirements

**MOL2 File Priority**:
1. **Primary**: Use input MOL2/SDF files (have correct bond orders)
2. **Fallback**: RDKit sanitization if MOL2 unavailable
3. **Never**: Trust GRO files for bond orders (no bond info)

**Implementation**:
```python
# ✅ CORRECT: Prioritize MOL2 for bond orders
if mol2_file.exists():
    assign_bond_orders_from_mol2(mol2_file)
else:
    # Fallback to RDKit (may fail for complex molecules)
    Chem.SanitizeMol(mol)
```

**Aromatic Ring Rendering**:
- ✅ Must show alternating single/double bonds
- ✅ Preserve bond orders from MOL2 template
- ✅ RDKit MCS matching for template alignment
- ❌ Don't use RDKit's Kekulize if it changes original bond orders

### Coordinate and Mapping Validation

**Coordinate Source Priority** (CRITICAL):
1. **MOL2** (highest priority - has bond orders)
2. **PDB** (secondary - no bond orders but full atom names)
3. **GRO** (fallback only - 5-char atom name limit)

**Atom Name Matching**:
- ✅ Track atom name changes during force field conversion
- ✅ Validate common atoms have same element on both sides
- ❌ Don't use GRO-only coordinates if MOL2/PDB available

**Mapping Validation Checks**:
```python
# Sanity checks for mapping results
assert common_count_left == common_count_right, "Common counts must match"

for atom_left, atom_right in common_pairs:
    assert atom_left.element == atom_right.element, "Elements must match"
    assert atom_left.name != "C" or not is_chain_end(atom_left), "Carbon at chain end?"
    assert atom_right.name != "H" or not has_multiple_bonds(atom_right), "H with multiple bonds?"
```

### Visualization-Specific Constraints

**Coloring Mode Implementation**:
- ✅ FEP Classification: common (gray), transformed A (red), transformed B (blue), surrounding (green)
- ✅ Element Colors: standard CPK coloring (C=gray, N=blue, O=red, H=white, etc.)
- ✅ Mode switch must be mutually exclusive (radio buttons, not checkboxes)

**Statistics Panel Requirements**:
- ✅ Show atom counts: common, transformed A, transformed B, surrounding
- ✅ Show total charges for both molecules
- ✅ Show mapping quality score (e.g., "✅ 31 common atoms")
- ✅ Real-time calculation from atoms (not copied from files)

**Build Log Display**:
- ✅ Collapsible section at bottom of HTML
- ✅ Show dist_cutoff, charge_cutoff, and other YAML parameters
- ✅ Include source index information if available
- ❌ Don't mix build log with dist_cutoff explanation

### Error Handling and User Feedback

**Warning Banner Triggers**:
- ⚠️ Gray atoms present (unclassified)
- ⚠️ Total charge ≠ 0 for neutral molecules
- ⚠️ Common count mismatch between molecules
- ⚠️ Placeholder atom names detected

**Error Messages**:
- ✅ Clear, actionable error messages
- ✅ Suggestions for fixing common issues
- ❌ No generic "error occurred" messages

**Debugging Support**:
- ✅ Include mapping log in HTML for troubleshooting
- ✅ Show atom names, types, charges in hover tooltips
- ✅ Display source_a_index and source_b_index when available
- ❌ Don't expose internal Python errors to end users

## Key Functions

### visualize_mapping_html()
Main entry point for HTML visualization generation.

**Parameters**:
- `mapping`: AtomMapping result from DistanceAtomMapper
- `pdb_a`, `pdb_b`: PDB coordinate files (receptor-aligned)
- `mol2_a`, `mol2_b`: Optional MOL2 files for bond order template
- `atoms_a`, `atoms_b`: Optional atom lists (filters dummy atoms automatically)
- `total_charge_a`, `total_charge_b**: Total charges for display
- `config`: FEP configuration dict (shows parameters in UI)
- `build_log`: Optional build log string for debugging

**Output**: Interactive HTML file with:
- Side-by-side molecule comparison
- FEP classification coloring (common/transformed/surrounding)
- Element coloring mode
- Zoom and pan controls
- Statistics panel (atom counts, total charges)
- Warning banner for unclassified atoms

## Architecture

```
prism/fep/visualize/
├── html.py           # HTML generation and template handling
├── molecule.py       # RDKit Mol preparation and bond order assignment
├── mapping.py        # 2D alignment and coordinate matching
├── highlight.py      # Classification color schemes
└── templates/
    ├── index.html    # Main HTML template
    ├── script.js     # Canvas rendering and interaction
    └── styles.css    # Styling and layout
```

## Molecule Preparation Pipeline

### prepare_mol_with_charges_and_labels()
**Universal strategy for ALL force fields** (GAFF2/OPLS/OpenFF/CGenFF):

1. **Coordinates**: Always use PDB (receptor-aligned)
2. **Bond orders**: Use MOL2 as template (if available)
3. **Atom names**: Match RDKit atoms to hybrid topology atoms
4. **Charges**: Annotate with mapping-modified charges

**Coordinate matching** (0.6 Å threshold):
```python
# Step 1: Match by 3D coordinates
for rdkit_atom in mol.GetAtoms():
    pos = rdkit_mol.GetConformer().GetAtomPosition(rdkit_atom.GetIdx())
    for topo_atom in topology_atoms:
        if distance(pos, topo_atom.coord) < 0.6:
            # Match found
```

**Fallback**: Base name matching (H07 → H07A)

### assign_bond_orders_from_mol2()
**Bond order correction for aromatic rings**:

1. Load MOL2 file with RDKit
2. Find maximum common substructure (MCS) between PDB and MOL2
3. Transfer bond orders from MOL2 to PDB
4. Sanitize molecule (infer aromaticity, hybridization)

**Fallback**: If MOL2 unavailable, use RDKit sanitization only.

## Quality Control Features

### Unclassified Atom Detection
**Gray atoms** (`rgb(200, 200, 200)`) indicate mapping failures.

**Causes**:
- RDKit `Chem.AddHs()` added extra hydrogens (FIXED: removed)
- Coordinate distance > 0.6 Å threshold
- Element mismatch
- Missing atom in hybrid topology

**Verification**:
```bash
grep "rgb(200, 200, 200)" output/mapping.html  # Should return 0
grep '"classification": "unknown"' output/mapping.html  # Should return 0
```

### Warning Banner
**Auto-generated** when `unclassified_count > 0`:
```html
<div class="warning-banner">
    <div class="warning-icon">⚠️</div>
    <div class="warning-content">
        <strong>Warning: N unclassified atom(s) detected</strong>
        <span>Some atoms could not be mapped...</span>
    </div>
</div>
```

### Total Charge Display
**Must use mapping-modified charges**, not original GAFF charges:
```python
# ✅ Correct: After mapping
total_charge_a = sum(atom.charge for atom in atoms_a)

# ❌ Wrong: Original charges
total_charge_a = real_total_charge_a
```

## Multi-Force Field Support

### OPLS-AA Special Handling
**Issue**: OPLS uses sequential atom types (opls_800, opls_801, ...) without chemical environment info

**Solution**: Match atoms by distance + element + charge (ignore type)

**Code**: `molecule.py::prepare_mol_with_charges_and_labels()`

### GAFF/CGenFF
**Full atom type matching** because types encode chemical environment:
- GAFF: ca (aromatic C), c3 (sp3 C), ha (aromatic H), hc (sp3 H)
- CGenFF: CA, CB, CG (different carbon environments)

## Canvas Rendering

### Transform Pipeline
**Critical**: Avoid double scaling
```javascript
// ✅ Correct: Single transform
ctx.translate(centerX, centerY);
ctx.scale(zoom, zoom);
ctx.translate(-centerX, -centerY);

// ❌ Wrong: Double scaling (causes zoom artifacts)
ctx.scale(zoom, zoom);
ctx.translate(x * zoom, y * zoom);  // Don't scale twice!
```

### Bond Rendering
**Bond order visualization**:
- Single bonds: Single line
- Double bonds: Double line (offset by 0.15 * scale)
- Aromatic bonds: Special coloring

**Code**: `templates/script.js::_drawBond()`

### Atom Labels
**Display priority**:
1. Real atom name from hybrid topology (e.g., "CD2", "H09A")
2. Fallback to element symbol (e.g., "C", "H")
3. Never show "Atom15" placeholder

## Common Issues

### PNG Generation Fails
**Symptom**: `ValueError: Depict error: Substructure match with reference not found`

**Cause**: Multi-point mutations (e.g., 39-8) too different for RDKit alignment

**Solution**: HTML works fine. PNG generates with warning using default coordinates.

### Atoms Appear as "Atom15"
**Problem**: RDKit atoms couldn't match to topology atoms

**Debugging**:
1. Check coordinate distance (< 0.6 Å)
2. Verify element symbols match
3. Check base name fallback table

**Code**: `molecule.py::prepare_mol_with_charges_and_labels()`, lines 256-336

### Button Styling Issues
**Symptom**: Only text highlights, not entire button

**Solution**: Set `margin: 0` in `.toggle-buttons` class

**Code**: `templates/styles.css:447`

## Testing

### Verification Checklist
- [ ] No gray atoms (`rgb(200, 200, 200)`)
- [ ] No unknown classifications
- [ ] Total charges ≈ 0 (mapping-modified)
- [ ] Common counts equal on both sides
- [ ] All atoms show real names (not "Atom15")
- [ ] Buttons highlight entire area
- [ ] Aromatic rings show bond order

### Test Scripts
```bash
cd tests/gxf/FEP/unit_test/oMeEtPh-EtPh

# Regenerate all modes
python tests/gxf/FEP/unit_test/regenerate_gaff_html.py --all-modes

# Check with Playwright
mcp__plugin_playwright_playwright__browser_navigate \
  file:///path/to/oMeEtPh-EtPh_gaff_ref.html
```

## Code Locations

- **HTML generation**: `prism/fep/visualize/html.py`
- **Molecule preparation**: `prism/fep/visualize/molecule.py`
- **2D alignment**: `prism/fep/visualize/mapping.py`
- **Color schemes**: `prism/fep/visualize/highlight.py`
- **Templates**: `prism/fep/visualize/templates/`

## Related Documentation

- `CLAUDE.md` (parent) - FEP module overview
- `../core/CLAUDE.md` - Atom mapping details
- `.claude/skills/fep-visualization.md` - Visualization verification guide
