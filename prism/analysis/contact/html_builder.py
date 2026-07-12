#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML Builder Module
Generate a complete, self-contained interactive HTML report with embedded
JavaScript visualization.

The report renders an interaction diagram (HTML5 canvas "interaction wheel")
where edges are coloured by *interaction type* (hydrogen bond, salt bridge,
pi-stacking, pi-cation, halogen, hydrophobic) and weighted by trajectory
occupancy, residue nodes are coloured by chemistry class, and the whole shell
supports a light/dark colour-blind-safe theme. It stays a single offline file
with no external/CDN dependencies.
"""

import json

from .utils import convert_numpy_types
from .interactions import INTERACTION_STYLE, INTERACTION_TYPES, RESIDUE_CLASS_COLORS


class HTMLBuilder:
    """Build complete HTML output with embedded JavaScript"""

    def generate_html(
        self,
        trajectory_file,
        topology_file,
        ligand_file,
        ligand_name,
        total_frames,
        ligand_data,
        contacts,
        stats,
        interaction_summary=None,
    ):
        """Generate the complete HTML file"""

        # Ensure all data is JSON serializable
        ligand_data = convert_numpy_types(ligand_data)
        contacts = convert_numpy_types(contacts)
        interaction_summary = convert_numpy_types(interaction_summary or {})

        # Import JavaScript code
        from .javascript_code import get_javascript_code

        # Create the HTML using string replacement to avoid format string issues
        html_template = self._get_html_template()

        # Replace placeholders one by one to avoid curly brace issues
        html_content = html_template
        html_content = html_content.replace("[[trajectory_file]]", str(trajectory_file))
        html_content = html_content.replace("[[topology_file]]", str(topology_file))
        html_content = html_content.replace("[[ligand_file]]", str(ligand_file))
        html_content = html_content.replace("[[ligand_name]]", str(ligand_name))
        html_content = html_content.replace("[[total_frames]]", str(total_frames))
        html_content = html_content.replace("[[total_contacts]]", str(stats["total_contacts"]))
        html_content = html_content.replace("[[high_freq_contacts]]", str(stats["high_freq_contacts"]))
        html_content = html_content.replace("[[max_freq_percent]]", str(stats["max_freq_percent"]))
        html_content = html_content.replace("[[n_ligand_atoms]]", str(len(ligand_data["atoms"])))
        html_content = html_content.replace("[[ligand_atoms_json]]", json.dumps(ligand_data["atoms"], indent=2))
        html_content = html_content.replace("[[contacts_json]]", json.dumps(contacts, indent=2))
        html_content = html_content.replace("[[ligand_bonds_json]]", json.dumps(ligand_data["bonds"], indent=2))

        # Interaction-typing metadata injected as JS globals.
        html_content = html_content.replace(
            "[[interaction_style_json]]", json.dumps(INTERACTION_STYLE, indent=2)
        )
        html_content = html_content.replace(
            "[[interaction_order_json]]", json.dumps(list(INTERACTION_TYPES))
        )
        html_content = html_content.replace(
            "[[residue_class_colors_json]]", json.dumps(RESIDUE_CLASS_COLORS, indent=2)
        )
        html_content = html_content.replace(
            "[[interaction_summary_json]]", json.dumps(interaction_summary, indent=2)
        )

        # Add JavaScript code
        html_content += get_javascript_code()

        return html_content

    def _get_html_template(self):
        """Return the HTML template with placeholders using [[ ]] instead of { }"""
        return r"""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PRISM · Protein–Ligand Interaction Report</title>
    <style>
        :root {
            --bg: #eef1f7;
            --bg-grad-1: #6a8dff;
            --bg-grad-2: #8b5cf6;
            --surface: #ffffff;
            --surface-2: #f6f8fc;
            --surface-3: #eef1f7;
            --border: rgba(20, 30, 60, 0.10);
            --text: #1f2733;
            --text-muted: #6b7689;
            --accent: #5b6ef5;
            --accent-2: #8b5cf6;
            --shadow: 0 10px 40px rgba(20, 30, 80, 0.12);
            --shadow-sm: 0 2px 10px rgba(20, 30, 80, 0.06);
            --radius: 16px;
            --radius-sm: 10px;
            --canvas-bg-1: #ffffff;
            --canvas-bg-2: #f4f7fc;
            --grid: rgba(91, 110, 245, 0.08);
        }
        html[data-theme="dark"] {
            --bg: #0e1117;
            --bg-grad-1: #1b2540;
            --bg-grad-2: #2a1b46;
            --surface: #161b25;
            --surface-2: #1c222e;
            --surface-3: #232a38;
            --border: rgba(255, 255, 255, 0.10);
            --text: #e6e9ef;
            --text-muted: #97a1b3;
            --accent: #7c8cff;
            --accent-2: #a78bfa;
            --shadow: 0 10px 40px rgba(0, 0, 0, 0.5);
            --shadow-sm: 0 2px 10px rgba(0, 0, 0, 0.4);
            --canvas-bg-1: #11151d;
            --canvas-bg-2: #161b25;
            --grid: rgba(124, 140, 255, 0.10);
        }
        * { box-sizing: border-box; }
        body {
            margin: 0;
            padding: 24px;
            font-family: 'Inter', 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background:
                radial-gradient(1200px 600px at 12% -10%, var(--bg-grad-1) 0%, transparent 55%),
                radial-gradient(1000px 500px at 110% 0%, var(--bg-grad-2) 0%, transparent 50%),
                var(--bg);
            color: var(--text);
            min-height: 100vh;
            transition: background 0.4s ease, color 0.3s ease;
        }
        .container {
            max-width: 1500px;
            margin: 0 auto;
            background: var(--surface);
            border-radius: 22px;
            box-shadow: var(--shadow);
            padding: 32px 34px 38px;
            border: 1px solid var(--border);
        }
        header.app-head {
            display: flex;
            align-items: flex-start;
            justify-content: space-between;
            gap: 16px;
            margin-bottom: 10px;
        }
        .title-wrap h1 {
            margin: 0;
            font-size: 28px;
            font-weight: 800;
            letter-spacing: -0.5px;
            background: linear-gradient(120deg, var(--accent) 0%, var(--accent-2) 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }
        .title-wrap .subtitle {
            color: var(--text-muted);
            font-size: 14px;
            margin-top: 4px;
            font-weight: 400;
        }
        .theme-toggle {
            border: 1px solid var(--border);
            background: var(--surface-2);
            color: var(--text);
            border-radius: 999px;
            padding: 9px 16px;
            cursor: pointer;
            font-size: 13px;
            font-weight: 600;
            display: inline-flex;
            align-items: center;
            gap: 8px;
            transition: all 0.2s ease;
            box-shadow: var(--shadow-sm);
        }
        .theme-toggle:hover { transform: translateY(-1px); border-color: var(--accent); }
        .file-info {
            background: var(--surface-2);
            border: 1px solid var(--border);
            border-radius: var(--radius-sm);
            padding: 12px 16px;
            margin: 14px 0 18px;
            font-size: 12.5px;
            color: var(--text-muted);
            display: flex;
            flex-wrap: wrap;
            gap: 8px 22px;
        }
        .file-info b { color: var(--text); font-weight: 600; }
        .toolbar {
            display: flex;
            justify-content: center;
            gap: 9px;
            margin-bottom: 14px;
            flex-wrap: wrap;
            align-items: center;
        }
        .control-btn {
            padding: 9px 16px;
            border: 1px solid transparent;
            border-radius: var(--radius-sm);
            background: linear-gradient(135deg, var(--accent) 0%, var(--accent-2) 100%);
            color: #fff;
            font-size: 13px;
            cursor: pointer;
            transition: all 0.22s cubic-bezier(0.4, 0, 0.2, 1);
            font-weight: 600;
            box-shadow: var(--shadow-sm);
        }
        .control-btn:hover { transform: translateY(-2px); box-shadow: 0 6px 18px rgba(91,110,245,0.35); }
        .control-btn.secondary {
            background: var(--surface-2);
            color: var(--text);
            border: 1px solid var(--border);
        }
        .control-btn.secondary:hover { border-color: var(--accent); }
        .control-btn.active { background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%); color: #fff; }
        .zoom-controls { display: flex; justify-content: center; gap: 10px; margin-bottom: 14px; }
        .zoom-btn {
            width: 42px; height: 42px; border-radius: 50%;
            border: 1px solid var(--border); background: var(--surface-2);
            color: var(--accent); font-size: 19px; cursor: pointer; font-weight: 700;
            transition: all 0.2s ease; box-shadow: var(--shadow-sm);
        }
        .zoom-btn:hover { background: var(--accent); color: #fff; transform: scale(1.08); }
        .export-controls {
            display: inline-flex; gap: 8px; align-items: center;
            background: var(--surface-2); padding: 6px 10px; border-radius: var(--radius-sm);
            border: 1px solid var(--border);
        }
        .export-controls select, .export-controls input {
            color: var(--text); background: var(--surface); padding: 8px;
            border: 1px solid var(--border); border-radius: 8px; font-size: 13px;
        }
        .export-controls input { width: 130px; }
        #canvas {
            border: 1px solid var(--border);
            border-radius: var(--radius);
            background: linear-gradient(to bottom, var(--canvas-bg-1), var(--canvas-bg-2));
            cursor: grab; display: block; margin: 0 auto; max-width: 100%;
            box-shadow: var(--shadow-sm);
        }
        #canvas:active { cursor: grabbing; }
        .canvas-card {
            background: var(--surface-2);
            border: 1px solid var(--border);
            border-radius: var(--radius);
            padding: 14px;
            overflow: auto;
        }
        .help-text {
            text-align: center; color: var(--text-muted); font-size: 12.5px;
            margin: 12px 0; padding: 9px; background: var(--surface-2);
            border-radius: var(--radius-sm); border: 1px solid var(--border);
        }
        /* occupancy + filters bar */
        .filter-bar {
            display: flex; flex-wrap: wrap; align-items: center; gap: 18px;
            background: var(--surface-2); border: 1px solid var(--border);
            border-radius: var(--radius-sm); padding: 12px 16px; margin: 12px 0;
        }
        .filter-bar label { font-size: 12.5px; color: var(--text-muted); font-weight: 600; }
        .occ-control { display: flex; align-items: center; gap: 10px; }
        .occ-control input[type=range] { width: 180px; accent-color: var(--accent); }
        .occ-value {
            font-weight: 700; color: var(--accent); min-width: 42px; text-align: right;
            font-variant-numeric: tabular-nums;
        }
        #typeFilters { display: flex; flex-wrap: wrap; gap: 8px; }
        .type-chip {
            display: inline-flex; align-items: center; gap: 7px; cursor: pointer;
            padding: 5px 11px; border-radius: 999px; font-size: 12px; font-weight: 600;
            background: var(--surface); border: 1px solid var(--border); color: var(--text);
            user-select: none; transition: all 0.15s ease;
        }
        .type-chip.off { opacity: 0.4; }
        .type-chip .swatch { width: 22px; height: 0; border-top-width: 3px; border-top-style: solid; }
        /* info grid */
        .panel-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
            gap: 16px; margin-top: 20px;
        }
        .card {
            background: var(--surface-2); border: 1px solid var(--border);
            border-radius: var(--radius); padding: 18px 20px; box-shadow: var(--shadow-sm);
        }
        .card h3 {
            margin: 0 0 14px; font-size: 14px; font-weight: 700; letter-spacing: 0.2px;
            color: var(--text); display: flex; align-items: center; gap: 8px;
        }
        .stats { display: grid; grid-template-columns: repeat(2, 1fr); gap: 12px; }
        .stat-item {
            background: var(--surface); border: 1px solid var(--border);
            border-radius: 12px; padding: 16px; text-align: center;
        }
        .stat-value {
            font-size: 26px; font-weight: 800;
            background: linear-gradient(135deg, var(--accent) 0%, var(--accent-2) 100%);
            -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;
            font-variant-numeric: tabular-nums;
        }
        .stat-label { font-size: 11px; color: var(--text-muted); margin-top: 6px; text-transform: uppercase; letter-spacing: 0.8px; }
        .legend-list { display: flex; flex-direction: column; gap: 9px; }
        .legend-item { display: flex; align-items: center; gap: 11px; font-size: 13px; color: var(--text); }
        .legend-item .count { margin-left: auto; color: var(--text-muted); font-size: 12px; font-weight: 600; font-variant-numeric: tabular-nums; }
        .legend-line { width: 30px; height: 0; border-top-width: 3px; border-top-style: solid; flex: none; }
        .legend-dot { width: 16px; height: 16px; border-radius: 50%; flex: none; box-shadow: 0 1px 3px rgba(0,0,0,0.25); }
        .color-bar-grad {
            height: 14px; width: 100%; border-radius: 8px; margin: 6px 0 4px;
            background: linear-gradient(to right, #95a5a6 0%, #3498db 30%, #8e44ad 60%, #e74c3c 100%);
        }
        .color-labels { display: flex; justify-content: space-between; font-size: 11px; color: var(--text-muted); font-weight: 600; }
        table.itable { width: 100%; border-collapse: collapse; font-size: 12.5px; }
        table.itable th, table.itable td { text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--border); }
        table.itable th { color: var(--text-muted); font-weight: 700; text-transform: uppercase; font-size: 10.5px; letter-spacing: 0.6px; }
        table.itable td .pill {
            display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11px;
            font-weight: 700; color: #fff;
        }
        table.itable .occbar { height: 7px; border-radius: 4px; background: var(--accent); }
        .occ-cell { display: flex; align-items: center; gap: 8px; }
        .table-wrap { max-height: 320px; overflow: auto; }
        .tooltip {
            position: absolute; background: rgba(20,26,38,0.96); color: #fff;
            padding: 10px 13px; border-radius: 8px; font-size: 12.5px; pointer-events: none;
            z-index: 1000; display: none; box-shadow: 0 8px 25px rgba(0,0,0,0.4);
            border: 1px solid rgba(255,255,255,0.12); max-width: 280px;
        }
        .muted { color: var(--text-muted); font-size: 12px; }
        @media (max-width: 720px) {
            body { padding: 12px; }
            .container { padding: 20px 16px 28px; }
            .stats { grid-template-columns: 1fr 1fr; }
        }
    </style>
</head>
<body>
    <div class="container">
        <header class="app-head">
            <div class="title-wrap">
                <h1>Protein–Ligand Interaction Report</h1>
                <div class="subtitle">Trajectory-resolved interaction diagram · edges coloured by interaction type, weighted by occupancy</div>
            </div>
            <button class="theme-toggle" id="themeToggle" onclick="toggleTheme()">🌙 Dark</button>
        </header>

        <div class="file-info">
            <span><b>Trajectory:</b> [[trajectory_file]]</span>
            <span><b>Topology:</b> [[topology_file]]</span>
            <span><b>Ligand:</b> [[ligand_file]] ([[ligand_name]])</span>
            <span><b>Frames:</b> [[total_frames]]</span>
        </div>

        <div class="zoom-controls">
            <button class="zoom-btn" onclick="zoomIn()">+</button>
            <button class="zoom-btn" onclick="zoomOut()">−</button>
            <button class="zoom-btn" onclick="centerView()" style="width: auto; padding: 0 15px; font-size:14px;">Center</button>
        </div>

        <div class="toolbar">
            <button class="control-btn" onclick="resetPositions()">Reset Layout</button>
            <button class="control-btn" onclick="rotateWheel()">Rotate Wheel</button>
            <button class="control-btn" onclick="toggleDistanceLock()" id="distanceLockBtn">🔒 Distance Locked</button>
            <button class="control-btn" onclick="toggle3DMode()" id="mode3DBtn">2D Mode</button>
            <button class="control-btn secondary" onclick="toggleConnections()" id="connectBtn">Hide Connections</button>
            <button class="control-btn secondary" onclick="toggleHydrogens()" id="hydrogenBtn">Hide H atoms</button>
            <button class="control-btn secondary" onclick="toggleGrid()" id="gridBtn">Hide Grid</button>
            <button class="control-btn secondary" onclick="rotateCanvas180()">Rotate 180°</button>
            <div class="export-controls">
                <select id="exportQuality" class="control-btn secondary" style="padding: 8px;">
                    <option value="1">Standard (1x)</option>
                    <option value="2">High (2x)</option>
                    <option value="4">Ultra HD (4x)</option>
                    <option value="8">8K (8x)</option>
                </select>
                <input type="text" id="exportFilename" placeholder="filename" value="interaction_report">
                <button class="control-btn secondary" onclick="exportImage()">Export PNG</button>
            </div>
        </div>

        <div class="filter-bar">
            <div class="occ-control">
                <label for="occSlider">Min. occupancy</label>
                <input type="range" id="occSlider" min="0" max="100" step="5" value="0" oninput="onOccupancyChange(this.value)">
                <span class="occ-value" id="occValue">0%</span>
            </div>
            <div id="typeFilters" aria-label="Interaction type filters"></div>
        </div>

        <p class="help-text">🖱️ 2D: drag residues to rearrange · 3D: drag to rotate the ligand · Shift/middle-drag to pan · scroll to zoom · hover an edge for details</p>

        <div class="canvas-card">
            <canvas id="canvas" width="1200" height="800"></canvas>
        </div>

        <div class="panel-grid">
            <div class="card">
                <h3>📊 Summary</h3>
                <div class="stats">
                    <div class="stat-item"><div class="stat-value">[[total_contacts]]</div><div class="stat-label">Contacting residues</div></div>
                    <div class="stat-item"><div class="stat-value">[[high_freq_contacts]]</div><div class="stat-label">High occupancy</div></div>
                    <div class="stat-item"><div class="stat-value">[[max_freq_percent]]%</div><div class="stat-label">Max occupancy</div></div>
                    <div class="stat-item"><div class="stat-value">[[n_ligand_atoms]]</div><div class="stat-label">Ligand atoms</div></div>
                </div>
                <div style="margin-top:16px;">
                    <div class="stat-label" style="margin-bottom:6px;">Contact frequency scale</div>
                    <div class="color-bar-grad"></div>
                    <div class="color-labels"><span>0%</span><span>25%</span><span>50%</span><span>75%</span><span>100%</span></div>
                </div>
            </div>

            <div class="card">
                <h3>🔗 Interaction types</h3>
                <div class="legend-list" id="interactionLegend"></div>
                <div class="muted" id="interactionLegendEmpty" style="display:none; margin-top:6px;">
                    No interaction typing available — edges are coloured by contact frequency.
                </div>
            </div>

            <div class="card">
                <h3>🧬 Residue chemistry</h3>
                <div class="legend-list" id="residueLegend"></div>
            </div>
        </div>

        <div class="card" style="margin-top:16px;">
            <h3>📋 Interaction details</h3>
            <div class="table-wrap">
                <table class="itable" id="interactionTable">
                    <thead>
                        <tr><th>Residue</th><th>Class</th><th>Interactions</th><th>Top occupancy</th></tr>
                    </thead>
                    <tbody id="interactionTableBody"></tbody>
                </table>
            </div>
        </div>
    </div>

    <div class="tooltip" id="tooltip"></div>

    <script>
        const LIGAND_ATOMS = [[ligand_atoms_json]];
        const CONTACTS = [[contacts_json]];
        const LIGAND_BONDS = [[ligand_bonds_json]];
        const INTERACTION_STYLE = [[interaction_style_json]];
        const INTERACTION_ORDER = [[interaction_order_json]];
        const RESIDUE_CLASS_COLORS = [[residue_class_colors_json]];
        const INTERACTION_SUMMARY = [[interaction_summary_json]];"""
