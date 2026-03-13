#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML Canvas visualization for FEP atom mapping - FIXED VERSION

Key fixes:
1. Proper Canvas transform scaling (no double scaling)
2. Correct molecule separation and centering
3. All text in English
4. Bonds and atoms stay connected during zoom
"""

import json
from pathlib import Path
from typing import Optional, Dict, List

from rdkit import Chem

from prism.fep.core.mapping import AtomMapping
from .molecule import prepare_mol_with_charges_and_labels


def _load_template(filename: str) -> str:
    """Load template file from templates directory."""
    template_dir = Path(__file__).parent / "templates"
    template_path = template_dir / filename
    with open(template_path, "r", encoding="utf-8") as f:
        return f.read()


def visualize_mapping_html(
    mapping: AtomMapping,
    pdb_a: str,
    pdb_b: str,
    mol2_a: Optional[str] = None,
    mol2_b: Optional[str] = None,
    atoms_a: Optional[List] = None,
    atoms_b: Optional[List] = None,
    output_path: Optional[str] = None,  # type: ignore[assignment]
    title: str = "FEP Mapping Visualization",
    ligand_a_name: str = "Ligand A",
    ligand_b_name: str = "Ligand B",
    config: Optional[dict] = None,
) -> None:
    """Generate interactive HTML visualization of atom mapping using Canvas.

    Parameters
    ----------
    config : dict, optional
        FEP configuration dictionary. If provided, displays parameters in collapsible panel.
    """
    # Prepare molecules
    mol_a = prepare_mol_with_charges_and_labels(pdb_a, mol2_a, atoms_a)
    mol_b = prepare_mol_with_charges_and_labels(pdb_b, mol2_b, atoms_b)

    # Prepare canvas data
    canvas_data_a = _prepare_canvas_data(mol_a, atoms_a, mapping, "a")
    canvas_data_b = _prepare_canvas_data(mol_b, atoms_b, mapping, "b")

    # Build correspondence map
    correspondence = _build_correspondence_map(mapping, canvas_data_a, canvas_data_b)

    # Generate HTML
    html = _generate_canvas_html(
        canvas_data_a, canvas_data_b, correspondence, mapping, title, ligand_a_name, ligand_b_name, config
    )

    # Save HTML
    if output_path is None:
        output_path = "mapping_visualization.html"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(str(output_path), "w", encoding="utf-8") as f:
        f.write(html)

    print(f"  ✓ HTML saved to {output_path}")


def _prepare_canvas_data(mol: Chem.Mol, atoms: List, mapping: AtomMapping, mol_id: str) -> Dict:
    """Prepare data for Canvas rendering - NO SCALING here."""
    # Build classification sets
    common_atoms_a = {a[0].name for a in mapping.common}
    common_atoms_b = {a[1].name for a in mapping.common}
    transformed_a = {a.name for a in mapping.transformed_a}
    transformed_b = {a.name for a in mapping.transformed_b}
    surrounding_a = {a.name for a in mapping.surrounding_a}
    surrounding_b = {a.name for a in mapping.surrounding_b}

    # Build charge and type lookup
    charge_dict = {atom.name: atom.charge for atom in atoms}
    type_dict = {atom.name: atom.atom_type for atom in atoms}

    # Get 2D coordinates
    conf = mol.GetConformer()
    atoms_data = []

    # Collect coordinates and calculate center
    raw_coords = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        raw_coords.append((pos.x, pos.y))

    # Calculate molecule center
    if raw_coords:
        avg_x = sum(c[0] for c in raw_coords) / len(raw_coords)
        avg_y = sum(c[1] for c in raw_coords) / len(raw_coords)
    else:
        avg_x, avg_y = 0.0, 0.0

    # Create atom data with centered coordinates (NO SCALING)
    for atom in mol.GetAtoms():
        name = atom.GetProp("name") if atom.HasProp("name") else f"Atom{atom.GetIdx()}"
        pos = conf.GetAtomPosition(atom.GetIdx())
        element = atom.GetSymbol()

        # Determine classification
        if mol_id == "a":
            if name in common_atoms_a:
                classification = "common"
            elif name in transformed_a:
                classification = "transformed"
            elif name in surrounding_a:
                classification = "surrounding"
            else:
                classification = "unknown"
        else:
            if name in common_atoms_b:
                classification = "common"
            elif name in transformed_b:
                classification = "transformed"
            elif name in surrounding_b:
                classification = "surrounding"
            else:
                classification = "unknown"

        # Get colors and radius
        fep_color = _get_classification_color(classification)
        element_color = _get_element_color(element)
        radius = _get_atom_radius(element)

        # Center and scale coordinates (RDKit 2D coords are ~1-2 Å units, scale to pixels)
        scale = 30.0
        x = (pos.x - avg_x) * scale
        y = (pos.y - avg_y) * scale

        atoms_data.append(
            {
                "id": f"{mol_id.upper()}{atom.GetIdx() + 1}",
                "name": name,
                "element": element,
                "x": float(x),
                "y": float(y),
                "radius": radius,
                "charge": charge_dict.get(name, 0.0),
                "type": type_dict.get(name, ""),
                "classification": classification,
                "fepColor": fep_color,
                "elementColor": element_color,
            }
        )

    # Build bonds list
    bonds = []
    for bond in mol.GetBonds():
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        bonds.append([f"{mol_id.upper()}{atom1_idx + 1}", f"{mol_id.upper()}{atom2_idx + 1}"])

    return {"atoms": atoms_data, "bonds": bonds}


def _get_classification_color(classification: str) -> str:
    """Get FEP classification color."""
    colors = {
        "common": "rgb(204, 229, 77)",
        "transformed": "rgb(255, 77, 77)",
        "surrounding": "rgb(77, 153, 255)",
        "unknown": "rgb(200, 200, 200)",
    }
    return colors.get(classification, "rgb(200, 200, 200)")


def _get_element_color(element: str) -> str:
    """Get element-based color (optimized CPK scheme for better visibility)."""
    colors = {
        "C": "#909090",  # Gray (better than black for visibility)
        "N": "#3050F8",  # Blue
        "O": "#FF0D0D",  # Red
        "S": "#FFC832",  # Gold/Orange (better than pure yellow)
        "P": "#FF8000",  # Orange
        "F": "#90E050",  # Light green
        "Cl": "#1FF01F",  # Green
        "Br": "#A62929",  # Dark red/brown
        "I": "#940094",  # Purple
        "H": "#E8E8E8",  # Light gray (visible on white background)
    }
    return colors.get(element, "#909090")


def _get_atom_radius(element: str) -> float:
    """Get atom radius for rendering."""
    radii = {"C": 15, "N": 14, "O": 13, "S": 18, "P": 17, "F": 12, "Cl": 16, "Br": 18, "I": 20, "H": 8}
    return radii.get(element, 15)


def _build_correspondence_map(mapping: AtomMapping, canvas_data_a: Dict, canvas_data_b: Dict) -> Dict:
    """Build correspondence map between atoms in A and B (including common and surrounding)."""
    correspondence = {}

    # Add common atom pairs
    for atom_a, atom_b in mapping.common:
        idx_a = None
        idx_b = None
        for i, atom in enumerate(canvas_data_a["atoms"]):
            if atom["name"] == atom_a.name:
                idx_a = i
                break
        for i, atom in enumerate(canvas_data_b["atoms"]):
            if atom["name"] == atom_b.name:
                idx_b = i
                break

        if idx_a is not None and idx_b is not None:
            correspondence[f"a_{idx_a}"] = f"b_{idx_b}"
            correspondence[f"b_{idx_b}"] = f"a_{idx_a}"

    # Add surrounding atom pairs (they also have correspondence)
    # Surrounding atoms are position-matched but have different charges/types
    for atom_a in mapping.surrounding_a:
        for atom_b in mapping.surrounding_b:
            # Check if they are at similar positions (distance-based matching)
            if atom_a.element == atom_b.element:
                import numpy as np

                dist = np.linalg.norm(atom_a.coord - atom_b.coord)
                if dist < 0.6:  # Same distance threshold as mapping
                    idx_a = None
                    idx_b = None
                    for i, atom in enumerate(canvas_data_a["atoms"]):
                        if atom["name"] == atom_a.name:
                            idx_a = i
                            break
                    for i, atom in enumerate(canvas_data_b["atoms"]):
                        if atom["name"] == atom_b.name:
                            idx_b = i
                            break

                    if idx_a is not None and idx_b is not None:
                        correspondence[f"a_{idx_a}"] = f"b_{idx_b}"
                        correspondence[f"b_{idx_b}"] = f"a_{idx_a}"
                    break

    return correspondence


def _generate_canvas_html(
    canvas_data_a: Dict,
    canvas_data_b: Dict,
    correspondence: Dict,
    mapping: AtomMapping,
    title: str,
    ligand_a_name: str,
    ligand_b_name: str,
    config: Optional[dict] = None,
) -> str:
    """Generate complete HTML with proper Canvas transforms."""

    # Build configuration panel HTML
    config_panel_html = ""
    if config:
        import os

        defaults = {"dist_cutoff": 0.6, "charge_cutoff": 0.05, "charge_common": "mean", "charge_reception": "pert"}
        tooltips = {
            "dist_cutoff": "Maximum distance (Å) between atoms to be considered as common",
            "charge_cutoff": "Maximum charge difference for common atoms",
            "charge_common": "Charge assignment for common atoms (ref/mut/mean)",
            "charge_reception": "Charge reception mode (pert/unique/surround)",
            "working_dir": "Working directory for this FEP calculation",
        }

        fep_cfg = config.get("fep", {}).get("mapping", {})
        case_info = config.get("case", {})

        def make_item(label, value, key, unit=""):
            is_default = key in defaults and value == defaults[key]
            cls = "" if is_default else " non-default"
            tip = tooltips.get(key, "")
            val = f"{value} {unit}".strip() if value is not None else "N/A"
            return f'<div class="config-item{cls}" title="{tip}"><span class="param-label">{label}:</span><span class="param-value">{val}</span></div>'

        rows = []

        # Working directory (never highlight as non-default) - First
        wd = os.getcwd()
        if case_info.get("path"):
            wd = os.path.join(wd, case_info["path"])
        # Force default style for working directory
        rows.append(
            f'<div class="config-row full-width"><div class="config-item" title="{tooltips["working_dir"]}"><span class="param-label">Working Directory:</span><span class="param-value">{wd}</span></div></div>'
        )

        # Force field information and FEP parameters (all in one row)
        all_items = []

        # Force field type and params
        ff_info = config.get("forcefield", {})
        if ff_info:
            if ff_info.get("type"):
                ff_type = ff_info["type"].upper()
                all_items.append(
                    f'<div class="config-item" title="Ligand force field type"><span class="param-label">Force Field:</span><span class="param-value">{ff_type}</span></div>'
                )
            ff_params = ff_info.get("params", {})
            if ff_params:
                # Display key force field parameters
                for key, value in ff_params.items():
                    param_label = key.replace("_", " ").title()
                    all_items.append(
                        f'<div class="config-item" title="Force field parameter: {key}"><span class="param-label">{param_label}:</span><span class="param-value">{value}</span></div>'
                    )

        # FEP parameters
        if fep_cfg:
            all_items.extend(
                [
                    make_item("Distance Cutoff", fep_cfg.get("dist_cutoff"), "dist_cutoff", "Å"),
                    make_item("Charge Cutoff", fep_cfg.get("charge_cutoff"), "charge_cutoff"),
                    make_item("Charge Common", fep_cfg.get("charge_common"), "charge_common"),
                    make_item("Charge Reception", fep_cfg.get("charge_reception"), "charge_reception"),
                ]
            )

        # Add all items in one row
        if all_items:
            rows.append(f'<div class="config-row">{"".join(all_items)}</div>')

        if rows:
            config_panel_html = f"""
    <div class="config-panel">
        <div class="config-header" onclick="toggleConfig()">
            <div><h3 style="display:inline;">Configuration Parameters</h3><span class="config-note">(Non-default values highlighted)</span></div>
            <span class="toggle-icon" id="config-toggle">▼</span>
        </div>
        <div class="config-content" id="config-content">
            {"".join(rows)}
        </div>
    </div>
"""

    # Convert data to JSON
    atoms_a_json = json.dumps(canvas_data_a["atoms"])
    atoms_b_json = json.dumps(canvas_data_b["atoms"])
    bonds_a_json = json.dumps(canvas_data_a["bonds"])
    bonds_b_json = json.dumps(canvas_data_b["bonds"])
    correspondence_json = json.dumps(correspondence)

    # Load templates
    css_content = _load_template("styles.css")
    js_content = _load_template("script.js")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
{css_content}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <h1>{title}</h1>
    </div>

    <div class="toolbar">
        <h3 class="toolbar-title">Display Options</h3>
        <div class="toolbar-content">
            <div class="toolbar-left">
                <div class="toolbar-row">
                    <strong>Coloring Mode:</strong>
                    <label><input type="radio" name="colorMode" value="fep" checked> FEP Classification</label>
                    <label><input type="radio" name="colorMode" value="element"> Element</label>
                    <label><input type="checkbox" id="toggle-charges"> Show Charges</label>
                    <label><input type="checkbox" id="toggle-labels" checked> Show Labels</label>
                </div>
                <div class="toolbar-row">
                    <strong>Controls:</strong> Drag to pan | Scroll to zoom (independent control for each molecule) | Hover for atom info
                </div>
            </div>
            <div class="toolbar-right">
                <div class="toolbar-buttons">
                    <button onclick="resetView()">Reset View</button>
                    <button onclick="exportPNG()">Export PNG</button>
                    <button class="secondary" onclick="scrollToAtomDetails()">View Atom Details</button>
                    <button class="secondary" onclick="window.print()">Print</button>
                </div>
            </div>
        </div>
    </div>

    <div class="legend">
        <h3>Legend</h3>
        <div class="legend-content">
            <div class="legend-left">
                <div class="legend-section">
                    <h4>FEP Classification</h4>
                    <div class="legend-items">
                        <div class="legend-item">
                            <div class="color-box common"></div>
                            <span><strong>Common:</strong> {len(mapping.common)}</span>
                        </div>
                        <div class="legend-item">
                            <div class="color-box transformed"></div>
                            <span><strong>Transformed A:</strong> {len(mapping.transformed_a)}</span>
                        </div>
                        <div class="legend-item">
                            <div class="color-box transformed"></div>
                            <span><strong>Transformed B:</strong> {len(mapping.transformed_b)}</span>
                        </div>
                        <div class="legend-item">
                            <div class="color-box surrounding"></div>
                            <span><strong>Surrounding A:</strong> {len(mapping.surrounding_a)}</span>
                        </div>
                        <div class="legend-item">
                            <div class="color-box surrounding"></div>
                            <span><strong>Surrounding B:</strong> {len(mapping.surrounding_b)}</span>
                        </div>
                    </div>
                </div>
            </div>
            <div class="legend-right">
                <div class="legend-section">
                    <h4>Element Colors</h4>
                    <div class="legend-items">
                        <div class="legend-item">
                            <div class="element-box dark-bg" style="background: #909090; border-color: #666;">C</div>
                            <span>Carbon</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box dark-bg" style="background: #3050F8; border-color: #1030D0;">N</div>
                            <span>Nitrogen</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box dark-bg" style="background: #FF0D0D; border-color: #CC0000;">O</div>
                            <span>Oxygen</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box light-bg" style="background: #FFC832; border-color: #CC9900; color: #333;">S</div>
                            <span>Sulfur</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box light-bg" style="background: #FF8000; border-color: #CC6600; color: #333;">P</div>
                            <span>Phosphorus</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box light-bg" style="background: #90E050; border-color: #60B020; color: #333;">F</div>
                            <span>Fluorine</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box light-bg" style="background: #1FF01F; border-color: #00C000; color: #333;">Cl</div>
                            <span>Chlorine</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box dark-bg" style="background: #A62929; border-color: #801010;">Br</div>
                            <span>Bromine</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box dark-bg" style="background: #940094; border-color: #600060;">I</div>
                            <span>Iodine</span>
                        </div>
                        <div class="legend-item">
                            <div class="element-box light-bg" style="background: #E8E8E8; border-color: #999; color: #333;">H</div>
                            <span>Hydrogen</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    {config_panel_html}

    <div class="canvas-container">
        <div class="molecule-label molecule-label-a">{ligand_a_name}</div>
        <div class="molecule-label molecule-label-b">{ligand_b_name}</div>
        <canvas id="canvas" width="1400" height="700"></canvas>
        <div class="tooltip" id="tooltip"></div>
    </div>

    <div class="atom-list" id="atom-list">
        <h3>Atom Details & Statistics</h3>

        <div class="stats-summary">
            <div class="stat-item">
                <div class="stat-value">{len(mapping.common)}</div>
                <div class="stat-label">Common Atoms</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">{len(mapping.transformed_a) + len(mapping.transformed_b)}</div>
                <div class="stat-label">Transformed Atoms</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">{len(mapping.surrounding_a) + len(mapping.surrounding_b)}</div>
                <div class="stat-label">Surrounding Atoms</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">{len(mapping.common) * 2 + len(mapping.transformed_a) + len(mapping.transformed_b) + len(mapping.surrounding_a) + len(mapping.surrounding_b)}</div>
                <div class="stat-label">Total Atoms</div>
            </div>
        </div>

        <table class="atom-table" id="atom-table">
            <thead>
                <tr>
                    <th>Ligand A: {ligand_a_name}</th>
                    <th>Element</th>
                    <th>Charge</th>
                    <th>Classification</th>
                    <th style="border-left: 2px solid #ddd;">Ligand B: {ligand_b_name}</th>
                    <th>Element</th>
                    <th>Charge</th>
                    <th>Classification</th>
                </tr>
            </thead>
            <tbody id="atom-table-body">
                <!-- Table content will be populated by JavaScript -->
            </tbody>
        </table>
    </div>

    <script>
        // Data
        const ATOMS_A = {atoms_a_json};
        const ATOMS_B = {atoms_b_json};
        const BONDS_A = {bonds_a_json};
        const BONDS_B = {bonds_b_json};
        const CORRESPONDENCE = {correspondence_json};

{js_content}
    </script>
</div>
</body>
</html>"""
