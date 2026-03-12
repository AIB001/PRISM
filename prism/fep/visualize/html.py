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


def visualize_mapping_html(
    mapping: AtomMapping,
    pdb_a: str,
    pdb_b: str,
    mol2_a: str,
    mol2_b: str,
    atoms_a: List,
    atoms_b: List,
    output_path: Optional[str] = None,  # type: ignore[assignment]
    title: str = "FEP Atom Mapping Visualization",
    ligand_a_name: str = "Ligand A",
    ligand_b_name: str = "Ligand B",
) -> None:
    """Generate interactive HTML visualization of atom mapping using Canvas."""
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
        canvas_data_a, canvas_data_b, correspondence, mapping, title, ligand_a_name, ligand_b_name
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
) -> str:
    """Generate complete HTML with proper Canvas transforms."""

    # Convert data to JSON
    atoms_a_json = json.dumps(canvas_data_a["atoms"])
    atoms_b_json = json.dumps(canvas_data_b["atoms"])
    bonds_a_json = json.dumps(canvas_data_a["bonds"])
    bonds_b_json = json.dumps(canvas_data_b["bonds"])
    correspondence_json = json.dumps(correspondence)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .header {{
            text-align: center;
            margin-bottom: 20px;
        }}
        .toolbar {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            text-align: center;
        }}
        .toolbar-group {{
            display: inline-block;
            margin: 0 20px;
            padding: 0 20px;
            border-right: 1px solid #ddd;
        }}
        .toolbar-group:last-child {{
            border-right: none;
        }}
        .toolbar label {{
            margin-right: 15px;
            cursor: pointer;
            font-size: 14px;
        }}
        .toolbar button {{
            padding: 8px 16px;
            background: #4CAF50;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            margin: 0 3px;
        }}
        .toolbar button:hover {{
            background: #45a049;
        }}
        .toolbar button.secondary {{
            background: #2196F3;
        }}
        .toolbar button.secondary:hover {{
            background: #1976D2;
        }}
        .canvas-container {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            position: relative;
        }}
        canvas {{
            border: 1px solid #ddd;
            border-radius: 4px;
            cursor: grab;
            display: block;
            width: 100%;
        }}
        canvas:active {{
            cursor: grabbing;
        }}
        .info-panel {{
            background: white;
            padding: 15px 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            text-align: center;
            color: #666;
            font-size: 13px;
        }}
        .molecule-label {{
            position: absolute;
            top: 30px;
            font-size: 18px;
            font-weight: bold;
            color: #333;
            background: rgba(255, 255, 255, 0.9);
            padding: 8px 16px;
            border-radius: 6px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.15);
        }}
        .molecule-label-a {{
            left: 25%;
            transform: translateX(-50%);
        }}
        .molecule-label-b {{
            right: 25%;
            transform: translateX(50%);
        }}
        .legend {{
            background: white;
            padding: 15px 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            width: 100%;
            max-width: none;
            margin: 0 0 20px 0;
        }}
        .legend h3 {{
            margin: 0 0 15px 0;
            color: #333;
            text-align: center;
            font-size: 18px;
        }}
        .legend-content {{
            display: flex;
            gap: 30px;
            align-items: flex-start;
            flex-wrap: wrap;
        }}
        .legend-left {{
            flex: 1;
            min-width: 300px;
        }}
        .legend-right {{
            min-width: 200px;
            flex-shrink: 0;
        }}
        @media (max-width: 768px) {{
            .legend-content {{
                flex-direction: column;
            }}
            .legend-right {{
                width: 100%;
            }}
        }}
        .legend-section {{
            margin-bottom: 15px;
        }}
        .legend-section:last-child {{
            margin-bottom: 0;
        }}
        .legend-section h4 {{
            margin: 0 0 8px 0;
            color: #555;
            font-size: 14px;
            font-weight: 600;
        }}
        .legend-items {{
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            align-items: center;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            font-size: 12px;
            white-space: nowrap;
        }}
        .color-box {{
            width: 28px;
            height: 18px;
            margin-right: 8px;
            border-radius: 3px;
            border: 1px solid #ddd;
            flex-shrink: 0;
        }}
        .common {{ background: rgb(204, 229, 77); }}
        .transformed {{ background: rgb(255, 77, 77); }}
        .surrounding {{ background: rgb(77, 153, 255); }}
        .element-box {{
            width: 20px;
            height: 20px;
            margin-right: 8px;
            border-radius: 50%;
            border: 2px solid #333;
            flex-shrink: 0;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 10px;
            font-weight: bold;
        }}
        .element-box.light-bg {{
            color: #333;
        }}
        .element-box.dark-bg {{
            color: white;
        }}
        .stats {{
            text-align: center;
            color: #666;
            font-size: 13px;
            background: white;
            padding: 12px 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            width: 100%;
            max-width: none;
            margin: 0;
        }}
        .tooltip {{
            position: absolute;
            background: rgba(0, 0, 0, 0.92);
            color: white;
            padding: 12px;
            border-radius: 8px;
            font-size: 13px;
            pointer-events: none;
            z-index: 1000;
            display: none;
            min-width: 220px;
            box-shadow: 0 6px 16px rgba(0,0,0,0.4);
            border: 1px solid rgba(255,255,255,0.15);
        }}
        .tooltip-title {{
            font-weight: bold;
            margin-bottom: 8px;
            color: #ffd700;
            font-size: 16px;
        }}
        .tooltip-row {{
            margin: 5px 0;
            font-size: 12px;
        }}
        .tooltip-label {{
            color: #aaa;
            margin-right: 10px;
        }}
        .correspondence {{
            margin-top: 10px;
            padding-top: 8px;
            border-top: 1px solid rgba(255,255,255,0.2);
            color: #4CAF50;
            font-size: 12px;
        }}
        .atom-list {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin: 20px 0;
            display: block;
        }}
        .atom-list h3 {{
            margin: 0 0 15px 0;
            color: #333;
            font-size: 18px;
        }}
        .stats-summary {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            margin-bottom: 20px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
        }}
        .stat-item {{
            text-align: center;
        }}
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            color: #333;
        }}
        .stat-label {{
            font-size: 12px;
            color: #666;
            margin-top: 4px;
        }}
        .atom-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
        }}
        .atom-table th {{
            background: #f8f9fa;
            padding: 10px 8px;
            text-align: left;
            border-bottom: 2px solid #ddd;
            font-weight: 600;
            color: #555;
            font-size: 12px;
        }}
        .atom-table td {{
            padding: 6px 8px;
            border-bottom: 1px solid #eee;
        }}
        .atom-table tr:hover {{
            background: #f8f9fa;
        }}
        .atom-name {{
            font-weight: bold;
            font-family: monospace;
        }}
        .atom-charge {{
            font-family: monospace;
            text-align: right;
        }}
        .classification-badge {{
            padding: 2px 8px;
            border-radius: 3px;
            font-size: 10px;
            font-weight: bold;
            display: inline-block;
        }}
        .badge-common {{ background: rgb(204, 229, 77); color: #333; }}
        .badge-transformed {{ background: rgb(255, 77, 77); color: white; }}
        .badge-surrounding {{ background: rgb(77, 153, 255); color: white; }}
            border-radius: 4px;
            font-size: 12px;
        }}
        .atom-item:hover {{
            background: #f0f0f0;
        }}
        .badge {{
            padding: 3px 6px;
            border-radius: 3px;
            font-size: 10px;
            color: white;
            font-weight: bold;
        }}
        .badge-common {{ background: rgb(154, 205, 50); }}
        .badge-transformed {{ background: rgb(220, 53, 69); }}
        .badge-surrounding {{ background: rgb(13, 110, 253); }}
        .info-panel {{
            background: white;
            padding: 10px 15px;
            border-radius: 4px;
            position: absolute;
            bottom: 10px;
            left: 10px;
            font-size: 12px;
            color: #666;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>{title}</h1>
    </div>

    <div class="toolbar">
        <div class="toolbar-group">
            <strong>Coloring Mode:</strong>
            <label><input type="radio" name="colorMode" value="fep" checked> FEP Classification</label>
            <label><input type="radio" name="colorMode" value="element"> Element</label>
        </div>
        <div class="toolbar-group">
            <label><input type="checkbox" id="toggle-charges"> Show Charges</label>
            <label><input type="checkbox" id="toggle-labels" checked> Show Labels</label>
        </div>
        <div class="toolbar-group">
            <button onclick="resetView()">Reset View</button>
            <button onclick="exportPNG()">Export PNG</button>
            <button class="secondary" onclick="scrollToAtomDetails()">View Atom Details</button>
            <button class="secondary" onclick="window.print()">Print</button>
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

    <div class="canvas-container">
        <div class="molecule-label molecule-label-a">Ligand A: {ligand_a_name}</div>
        <div class="molecule-label molecule-label-b">Ligand B: {ligand_b_name}</div>
        <canvas id="canvas" width="1400" height="700"></canvas>
        <div class="tooltip" id="tooltip"></div>
    </div>

    <div class="info-panel">
        <strong>Controls:</strong> Drag to pan | Scroll to zoom (independent control for each molecule) | Hover for atom info
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

        // Canvas setup
        const canvas = document.getElementById('canvas');
        const ctx = canvas.getContext('2d');
        const tooltip = document.getElementById('tooltip');
        const infoPanel = document.getElementById('infoPanel');

        // Canvas dimensions
        const canvasWidth = 1400;
        const canvasHeight = 700;

        // Layout: each molecule gets half the canvas
        const halfWidth = canvasWidth / 2;
        const centerX_A = halfWidth / 2;     // 350 - center of left half
        const centerX_B = halfWidth + halfWidth / 2;  // 1050 - center of right half
        const centerY = canvasHeight / 2;   // 350

        // View state - INDEPENDENT for each molecule
        let molAState = {{
            viewOffset: {{ x: 0, y: 0 }},
            zoom: 1.0,
            isDragging: false,
            dragStart: {{ x: 0, y: 0 }}
        }};
        let molBState = {{
            viewOffset: {{ x: 0, y: 0 }},
            zoom: 1.0,
            isDragging: false,
            dragStart: {{ x: 0, y: 0 }}
        }};
        let showLabels = true;
        let showCharges = false;
        let colorMode = 'fep';
        let hoveredAtom = null;
        let hoveredMolecule = null;
        let activeMolecule = null;  // Which molecule is being interacted with

        // Calculate initial zoom for a single molecule
        function calculateInitialZoomForMolecule(atoms) {{
            let maxX = 0, maxY = 0;
            atoms.forEach(atom => {{
                maxX = Math.max(maxX, Math.abs(atom.x));
                maxY = Math.max(maxY, Math.abs(atom.y));
            }});

            const padding = 50;
            const availableWidth = halfWidth - padding * 2;
            const availableHeight = canvasHeight - padding * 2;

            const zoomX = availableWidth / (maxX * 2 + 0.1);
            const zoomY = availableHeight / (maxY * 2 + 0.1);

            return Math.min(zoomX, zoomY, 1.5);
        }}

        // Initialize zoom for each molecule independently
        molAState.zoom = calculateInitialZoomForMolecule(ATOMS_A);
        molBState.zoom = calculateInitialZoomForMolecule(ATOMS_B);

        // Get color based on current mode
        function getAtomFillColor(atom) {{
            if (colorMode === 'fep') {{
                return atom.fepColor;
            }} else {{
                return atom.elementColor;
            }}
        }}

        // Draw bonds
        function drawBonds(atoms, bonds, offsetX, offsetY) {{
            bonds.forEach(([id1, id2]) => {{
                const atom1 = atoms.find(a => a.id === id1);
                const atom2 = atoms.find(a => a.id === id2);
                if (atom1 && atom2) {{
                    ctx.beginPath();
                    ctx.moveTo(atom1.x, atom1.y);
                    ctx.lineTo(atom2.x, atom2.y);
                    ctx.strokeStyle = '#666';
                    ctx.lineWidth = 2;
                    ctx.stroke();
                }}
            }});
        }}

        // Draw atoms
        function drawAtoms(atoms, offsetX, offsetY, moleculeId) {{
            atoms.forEach((atom, index) => {{
                const isHovered = hoveredAtom === index && hoveredMolecule === moleculeId;

                // Draw atom circle
                ctx.beginPath();
                ctx.arc(atom.x, atom.y, atom.radius, 0, 2 * Math.PI);

                // Fill color based on mode
                ctx.fillStyle = getAtomFillColor(atom);
                ctx.fill();

                // Highlight border if hovered
                if (isHovered) {{
                    ctx.strokeStyle = '#FFD700';
                    ctx.lineWidth = 3;
                    ctx.stroke();
                }} else {{
                    // Border always uses element color
                    ctx.strokeStyle = atom.elementColor;
                    ctx.lineWidth = 1.5;
                    ctx.stroke();
                }}

                // Draw label
                if (showLabels) {{
                    ctx.fillStyle = (atom.element === 'H' || colorMode === 'fep') ? '#333' : '#fff';
                    ctx.font = 'bold 11px Arial';
                    ctx.textAlign = 'center';
                    ctx.textBaseline = 'middle';
                    ctx.fillText(atom.name, atom.x, atom.y);
                }}

                // Draw charge if requested
                if (showCharges) {{
                    ctx.fillStyle = '#333';
                    ctx.font = '9px Arial';
                    ctx.fillText(atom.charge.toFixed(4), atom.x, atom.y + atom.radius + 12);
                }}
            }});
        }}

        // Main draw function - uses Canvas transforms
        function draw() {{
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            // Draw divider line
            ctx.save();
            ctx.setTransform(1, 0, 0, 1, 0, 0); // Reset transform for absolute positioning
            ctx.strokeStyle = '#ddd';
            ctx.lineWidth = 1;
            ctx.setLineDash([5, 5]);
            ctx.beginPath();
            ctx.moveTo(halfWidth, 50);
            ctx.lineTo(halfWidth, canvasHeight - 20);
            ctx.stroke();
            ctx.setLineDash([]);
            ctx.restore();

            // Draw molecule A on left side - INDEPENDENT TRANSFORM
            ctx.save();
            ctx.translate(molAState.viewOffset.x + centerX_A, molAState.viewOffset.y + centerY);
            ctx.scale(molAState.zoom, molAState.zoom);
            drawBonds(ATOMS_A, BONDS_A, 0, 0);
            drawAtoms(ATOMS_A, 0, 0, 'a');
            ctx.restore();

            // Draw molecule B on right side - INDEPENDENT TRANSFORM
            ctx.save();
            ctx.translate(molBState.viewOffset.x + centerX_B, molBState.viewOffset.y + centerY);
            ctx.scale(molBState.zoom, molBState.zoom);
            drawBonds(ATOMS_B, BONDS_B, 0, 0);
            drawAtoms(ATOMS_B, 0, 0, 'b');
            ctx.restore();
        }}

        // Find atom at position - converts mouse to world coordinates for each molecule
        function findAtomAtPosition(mouseX, mouseY) {{
            // Check molecule A with its independent transform
            const worldX_A = (mouseX - molAState.viewOffset.x - centerX_A) / molAState.zoom;
            const worldY_A = (mouseY - molAState.viewOffset.y - centerY) / molAState.zoom;

            for (let i = 0; i < ATOMS_A.length; i++) {{
                const atom = ATOMS_A[i];
                const dx = worldX_A - atom.x;
                const dy = worldY_A - atom.y;
                if (Math.sqrt(dx * dx + dy * dy) < atom.radius) {{
                    return {{ atom, index: i, molecule: 'a' }};
                }}
            }}

            // Check molecule B with its independent transform
            const worldX_B = (mouseX - molBState.viewOffset.x - centerX_B) / molBState.zoom;
            const worldY_B = (mouseY - molBState.viewOffset.y - centerY) / molBState.zoom;

            for (let i = 0; i < ATOMS_B.length; i++) {{
                const atom = ATOMS_B[i];
                const dx = worldX_B - atom.x;
                const dy = worldY_B - atom.y;
                if (Math.sqrt(dx * dx + dy * dy) < atom.radius) {{
                    return {{ atom, index: i, molecule: 'b' }};
                }}
            }}

            return null;
        }}

        // Show tooltip
        function showTooltip(atom, clientX, clientY, molecule, index) {{
            const classificationLabels = {{
                'common': 'Common',
                'transformed': 'Transformed',
                'surrounding': 'Surrounding'
            }};

            let content = `
                <div class="tooltip-title">${{atom.name}}</div>
                <div class="tooltip-row">
                    <span class="tooltip-label">Element:</span>
                    <span>${{atom.element}}</span>
                </div>
                <div class="tooltip-row">
                    <span class="tooltip-label">Charge:</span>
                    <span>${{atom.charge.toFixed(4)}}</span>
                </div>
                <div class="tooltip-row">
                    <span class="tooltip-label">Classification:</span>
                    <span>${{classificationLabels[atom.classification] || atom.classification}}</span>
                </div>
            `;

            // Add correspondence info
            const key = `${{molecule}}_${{index}}`;
            if (CORRESPONDENCE[key]) {{
                const parts = CORRESPONDENCE[key].split('_');
                const targetMol = parts[0];
                const targetIdx = parseInt(parts[1]);
                const targetAtoms = targetMol === 'a' ? ATOMS_A : ATOMS_B;
                const targetAtom = targetAtoms[targetIdx];

                if (targetAtom) {{
                    content += `
                        <div class="correspondence">
                            ↔ Corresponds to: ${{targetAtom.name}} (${{targetAtom.element}}, charge: ${{targetAtom.charge.toFixed(4)}})
                        </div>
                    `;
                }}
            }}

            tooltip.innerHTML = content;
            tooltip.style.display = 'block';
            tooltip.style.left = (clientX + 15) + 'px';
            tooltip.style.top = (clientY - 10) + 'px';
        }}

        // Hide tooltip
        function hideTooltip() {{
            tooltip.style.display = 'none';
        }}

        // Mouse events - INDEPENDENT for each molecule
        canvas.addEventListener('mousedown', (e) => {{
            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;

            // Determine which molecule was clicked based on X position
            if (mouseX < halfWidth) {{
                activeMolecule = 'a';
                molAState.isDragging = true;
                molAState.dragStart = {{ x: e.clientX - molAState.viewOffset.x, y: e.clientY - molAState.viewOffset.y }};
            }} else {{
                activeMolecule = 'b';
                molBState.isDragging = true;
                molBState.dragStart = {{ x: e.clientX - molBState.viewOffset.x, y: e.clientY - molBState.viewOffset.y }};
            }}
            canvas.style.cursor = 'grabbing';
        }});

        canvas.addEventListener('mousemove', (e) => {{
            // Handle dragging for active molecule
            if (molAState.isDragging && activeMolecule === 'a') {{
                molAState.viewOffset.x = e.clientX - molAState.dragStart.x;
                molAState.viewOffset.y = e.clientY - molAState.dragStart.y;
                draw();
            }} else if (molBState.isDragging && activeMolecule === 'b') {{
                molBState.viewOffset.x = e.clientX - molBState.dragStart.x;
                molBState.viewOffset.y = e.clientY - molBState.dragStart.y;
                draw();
            }} else {{
                // Handle hover for both molecules
                const rect = canvas.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;

                const found = findAtomAtPosition(x, y);

                if (found) {{
                    if (hoveredAtom !== found.index || hoveredMolecule !== found.molecule) {{
                        hoveredAtom = found.index;
                        hoveredMolecule = found.molecule;
                        showTooltip(found.atom, e.clientX, e.clientY, found.molecule, found.index);
                        draw();
                    }}
                }} else {{
                    if (hoveredAtom !== null) {{
                        hoveredAtom = null;
                        hoveredMolecule = null;
                        hideTooltip();
                        draw();
                    }}
                }}
            }}
        }});

        canvas.addEventListener('mouseup', () => {{
            molAState.isDragging = false;
            molBState.isDragging = false;
            activeMolecule = null;
            canvas.style.cursor = 'grab';
        }});

        // Zoom - INDEPENDENT for each molecule (like prism/analysis/contact)
        canvas.addEventListener('wheel', (e) => {{
            e.preventDefault();
            const delta = e.deltaY > 0 ? 0.9 : 1.1;

            // Determine which molecule to zoom based on mouse X position
            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;

            let targetState;
            let centerX;
            if (mouseX < halfWidth) {{
                targetState = molAState;
                centerX = centerX_A;
            }} else {{
                targetState = molBState;
                centerX = centerX_B;
            }}

            const newZoom = targetState.zoom * delta;

            if (newZoom >= 0.3 && newZoom <= 3) {{
                const mouseY = e.clientY - rect.top;

                // Convert to world coordinates before zoom
                const worldX = (mouseX - targetState.viewOffset.x - centerX) / targetState.zoom;
                const worldY = (mouseY - targetState.viewOffset.y - centerY) / targetState.zoom;

                targetState.zoom = newZoom;

                // Adjust view offset to zoom towards mouse position
                targetState.viewOffset.x = mouseX - worldX * targetState.zoom - centerX;
                targetState.viewOffset.y = mouseY - worldY * targetState.zoom - centerY;

                draw();
            }}
        }});
        // Toggle handlers
        document.getElementById('toggle-charges').addEventListener('change', (e) => {{
            showCharges = e.target.checked;
            draw();  // Redraw to show/hide charge labels on canvas
            draw();
        }});

        document.getElementById('toggle-labels').addEventListener('change', (e) => {{
            showLabels = e.target.checked;
            draw();
        }});

        // Color mode handlers
        document.querySelectorAll('input[name="colorMode"]').forEach(radio => {{
            radio.addEventListener('change', (e) => {{
                colorMode = e.target.value;
                draw();
            }});
        }});

        // Render atom list
        function renderAtomList() {{
            renderAtomTable();
        }}

        function renderAtomTable() {{
            const tbody = document.getElementById('atom-table-body');
            const classificationLabels = {{
                'common': 'Common',
                'transformed': 'Transformed',
                'surrounding': 'Surrounding'
            }};

            let html = '';

            // Create a map of all atoms with their pairs
            const atomPairs = [];

            // Add common atom pairs
            ATOMS_A.forEach((atomA, indexA) => {{
                if (atomA.classification === 'common') {{
                    const keyA = `a_${{indexA}}`;
                    const correspondingKey = CORRESPONDENCE[keyA];
                    if (correspondingKey) {{
                        const parts = correspondingKey.split('_');
                        const indexB = parseInt(parts[1]);
                        const atomB = ATOMS_B[indexB];
                        if (atomB) {{
                            atomPairs.push({{ atomA, atomB, type: 'common' }});
                        }}
                    }}
                }}
            }});

            // Add transformed A atoms (no pair in B)
            ATOMS_A.forEach(atomA => {{
                if (atomA.classification === 'transformed') {{
                    atomPairs.push({{ atomA, atomB: null, type: 'transformed_a' }});
                }}
            }});

            // Add transformed B atoms (no pair in A)
            ATOMS_B.forEach(atomB => {{
                if (atomB.classification === 'transformed') {{
                    atomPairs.push({{ atomA: null, atomB, type: 'transformed_b' }});
                }}
            }});

            // Add surrounding atom pairs
            ATOMS_A.forEach((atomA, indexA) => {{
                if (atomA.classification === 'surrounding') {{
                    const keyA = `a_${{indexA}}`;
                    const correspondingKey = CORRESPONDENCE[keyA];
                    if (correspondingKey) {{
                        const parts = correspondingKey.split('_');
                        const indexB = parseInt(parts[1]);
                        const atomB = ATOMS_B[indexB];
                        if (atomB) {{
                            atomPairs.push({{ atomA, atomB, type: 'surrounding' }});
                        }}
                    }}
                }}
            }});

            // Render table rows
            atomPairs.forEach(pair => {{
                const {{ atomA, atomB, type }} = pair;

                html += '<tr>';

                // Ligand A columns
                if (atomA) {{
                    const badgeClass = 'badge-' + atomA.classification;
                    html += `
                        <td class="atom-name">${{atomA.name}}</td>
                        <td>${{atomA.element}}</td>
                        <td class="atom-charge">${{atomA.charge.toFixed(4)}}</td>
                        <td><span class="classification-badge ${{badgeClass}}">${{classificationLabels[atomA.classification]}}</span></td>
                    `;
                }} else {{
                    html += '<td>—</td><td>—</td><td>—</td><td><span class="classification-badge badge-transformed">—</span></td>';
                }}

                // Ligand B columns
                if (atomB) {{
                    const badgeClass = 'badge-' + atomB.classification;
                    html += `
                        <td class="atom-name" style="border-left: 2px solid #ddd;">${{atomB.name}}</td>
                        <td>${{atomB.element}}</td>
                        <td class="atom-charge">${{atomB.charge.toFixed(4)}}</td>
                        <td><span class="classification-badge ${{badgeClass}}">${{classificationLabels[atomB.classification]}}</span></td>
                    `;
                }} else {{
                    html += '<td style="border-left: 2px solid #ddd;">—</td><td>—</td><td>—</td><td><span class="classification-badge badge-transformed">—</span></td>';
                }}

                html += '</tr>';
            }});

            tbody.innerHTML = html;
        }}

        function renderAtomListColumn(containerId, atoms, molId) {{
            // Legacy function - no longer used
            const container = document.getElementById(containerId);
            const classificationLabels = {{
                'common': 'Common',
                'transformed': 'Transformed',
                'surrounding': 'Surrounding'
            }};

            let html = '';
            atoms.forEach((atom, index) => {{
                const badgeClass = 'badge-' + atom.classification;
                const key = `${{molId.toLowerCase()}}_${{index}}`;
                let corresponding = '';

                if (CORRESPONDENCE[key]) {{
                    const parts = CORRESPONDENCE[key].split('_');
                    const targetMol = parts[0];
                    const targetIdx = parseInt(parts[1]);
                    const targetAtoms = targetMol === 'a' ? ATOMS_A : ATOMS_B;
                    const targetAtom = targetAtoms[targetIdx];
                    if (targetAtom) {{
                        corresponding = `<span style="color: #4CAF50; margin-left: 10px;">↔ ${{targetAtom.name}}</span>`;
                    }}
                }}

                html += `
                    <div class="atom-item">
                        <div>
                            <span style="font-weight: bold;">${{atom.name}}</span>
                            <span class="badge ${{badgeClass}}">${{classificationLabels[atom.classification]}}</span>
                            <span style="margin-left: 8px;">${{atom.element}}</span>
                            <span style="margin-left: 8px;">charge: ${{atom.charge.toFixed(4)}}</span>
                            ${{corresponding}}
                        </div>
                    </div>
                `;
            }});
            container.innerHTML = html;
        }}

        // Reset view - reset BOTH molecules independently
        function resetView() {{
            molAState.viewOffset = {{ x: 0, y: 0 }};
            molAState.zoom = calculateInitialZoomForMolecule(ATOMS_A);
            molBState.viewOffset = {{ x: 0, y: 0 }};
            molBState.zoom = calculateInitialZoomForMolecule(ATOMS_B);
            hoveredAtom = null;
            hoveredMolecule = null;
            activeMolecule = null;
            hideTooltip();
            draw();
        }}

        // Export PNG
        function exportPNG() {{
            const link = document.createElement('a');
            link.download = 'fep_mapping.png';
            link.href = canvas.toDataURL();
            link.click();
        }}

        function scrollToAtomDetails() {{
            const atomList = document.getElementById('atom-list');
            if (atomList) {{
                atomList.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
            }}
        }}

        // Initial draw
        draw();
        renderAtomTable();

        console.log('Interactive Canvas visualization loaded');
        console.log('Ligand A:', ATOMS_A.length, 'atoms, initial zoom:', molAState.zoom.toFixed(3));
        console.log('Ligand B:', ATOMS_B.length, 'atoms, initial zoom:', molBState.zoom.toFixed(3));
        console.log('Correspondence map:', CORRESPONDENCE);
        console.log('Instructions: Each molecule has INDEPENDENT pan/zoom control');
        console.log('  - Drag on left/right side to pan that molecule');
        console.log('  - Scroll on left/right side to zoom that molecule');
    </script>
</body>
</html>"""
