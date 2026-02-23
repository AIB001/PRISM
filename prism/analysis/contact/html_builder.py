#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML Builder Module
Generate complete HTML file with JavaScript visualization
"""

import json
from .utils import convert_numpy_types


class HTMLBuilder:
    """Build complete HTML output with embedded JavaScript"""

    def generate_html(
        self, trajectory_file, topology_file, ligand_file, ligand_name, total_frames, ligand_data, contacts, stats
    ):
        """Generate the complete HTML file"""

        # Ensure all data is JSON serializable
        ligand_data = convert_numpy_types(ligand_data)
        contacts = convert_numpy_types(contacts)

        # Import JavaScript code
        from .javascript_code import get_javascript_code

        # Create the HTML using string replacement to avoid format string issues
        html_template = self._get_html_template()

        # Replace placeholders one by one to avoid curly brace issues
        html_content = html_template
        html_content = html_content.replace("[[trajectory_file]]", trajectory_file)
        html_content = html_content.replace("[[topology_file]]", topology_file)
        html_content = html_content.replace("[[ligand_file]]", ligand_file)
        html_content = html_content.replace("[[ligand_name]]", ligand_name)
        html_content = html_content.replace("[[total_frames]]", str(total_frames))
        html_content = html_content.replace("[[total_contacts]]", str(stats["total_contacts"]))
        html_content = html_content.replace("[[high_freq_contacts]]", str(stats["high_freq_contacts"]))
        html_content = html_content.replace("[[max_freq_percent]]", str(stats["max_freq_percent"]))
        html_content = html_content.replace("[[n_ligand_atoms]]", str(len(ligand_data["atoms"])))
        html_content = html_content.replace("[[ligand_atoms_json]]", json.dumps(ligand_data["atoms"], indent=2))
        html_content = html_content.replace("[[contacts_json]]", json.dumps(contacts, indent=2))
        html_content = html_content.replace("[[ligand_bonds_json]]", json.dumps(ligand_data["bonds"], indent=2))

        # Add JavaScript code
        html_content += get_javascript_code()

        return html_content

    def _get_html_template(self):
        """Return the HTML template with placeholders using [[ ]] instead of { }"""
        return r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Enhanced Protein-Ligand Contact Analysis</title>
    <style>
        body {
            margin: 0; 
            padding: 20px; 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }
        .container {
            max-width: 1400px; 
            margin: 0 auto; 
            background: rgba(255, 255, 255, 0.98); 
            border-radius: 20px; 
            box-shadow: 0 20px 60px rgba(0,0,0,0.3), 0 0 100px rgba(102, 126, 234, 0.1);
            padding: 35px;
            backdrop-filter: blur(10px);
        }
        h1 {
            text-align: center; 
            color: #2C3E50; 
            font-size: 32px; 
            margin-bottom: 10px; 
            font-weight: 700;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            text-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .subtitle {
            text-align: center; 
            color: #7F8C8D; 
            font-size: 16px; 
            margin-bottom: 30px;
            font-weight: 300;
            letter-spacing: 0.5px;
        }
        #canvas {
            border: 1px solid rgba(102, 126, 234, 0.2);
            border-radius: 15px; 
            background: linear-gradient(to bottom, #ffffff, #f8f9fa);
            cursor: grab; 
            display: block; 
            margin: 0 auto;
            box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        }
        #canvas:active { cursor: grabbing; }
        .controls {
            display: flex; 
            justify-content: center; 
            gap: 12px; 
            margin-bottom: 20px; 
            flex-wrap: wrap;
            align-items: center;
        }
        .control-btn {
            padding: 12px 24px; 
            border: none; 
            border-radius: 10px; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; 
            font-size: 14px; 
            cursor: pointer; 
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            font-weight: 500;
            letter-spacing: 0.3px;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);
        }
        .control-btn:hover {
            transform: translateY(-2px) scale(1.02);
            box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4);
        }
        .control-btn.secondary {
            background: linear-gradient(135deg, #95a5a6 0%, #7f8c8d 100%);
            box-shadow: 0 4px 15px rgba(149, 165, 166, 0.3);
        }
        .control-btn.secondary:hover {
            box-shadow: 0 6px 20px rgba(149, 165, 166, 0.4);
        }
        .control-btn.active {
            background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%);
            box-shadow: 0 4px 15px rgba(231, 76, 60, 0.3);
        }
        .zoom-controls {
            display: flex; 
            justify-content: center; 
            gap: 10px; 
            margin-bottom: 15px; 
        }
        .zoom-btn {
            width: 45px; 
            height: 45px; 
            border-radius: 50%; 
            border: 2px solid transparent;
            background: linear-gradient(white, white) padding-box,
                        linear-gradient(135deg, #667eea 0%, #764ba2 100%) border-box;
            color: #667eea; 
            font-size: 20px; 
            cursor: pointer; 
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            font-weight: 600;
            box-shadow: 0 2px 10px rgba(102, 126, 234, 0.2);
        }
        .zoom-btn:hover {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) padding-box,
                        linear-gradient(135deg, #667eea 0%, #764ba2 100%) border-box;
            color: white;
            transform: scale(1.1);
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
        }
        .info-panel {
            background: linear-gradient(135deg, #f5f7fa 0%, #ecf0f1 100%);
            border-radius: 15px; 
            padding: 25px; 
            margin-top: 25px;
            box-shadow: inset 0 2px 10px rgba(0,0,0,0.05);
        }
        .file-info {
            background: linear-gradient(135deg, #dfe6e9 0%, #b2bec3 100%);
            border-radius: 10px; 
            padding: 15px; 
            margin-bottom: 15px; 
            font-size: 13px;
            color: #2c3e50;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }
        .stats {
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); 
            gap: 15px; 
            margin: 20px 0; 
        }
        .stat-item {
            background: white; 
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
            transition: all 0.3s ease;
            border: 1px solid rgba(102, 126, 234, 0.1);
        }
        .stat-item:hover {
            transform: translateY(-3px);
            box-shadow: 0 8px 25px rgba(102, 126, 234, 0.15);
        }
        .stat-value {
            font-size: 28px; 
            font-weight: bold; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .stat-label {
            font-size: 12px; 
            color: #7F8C8D; 
            margin-top: 8px;
            text-transform: uppercase;
            letter-spacing: 1px;
        }
        .legend {
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); 
            gap: 15px; 
            margin-top: 20px; 
        }
        .legend-item {
            display: flex; 
            align-items: center; 
            gap: 12px; 
            font-size: 14px; 
            color: #2C3E50;
            padding: 8px;
            border-radius: 8px;
            transition: background 0.2s;
        }
        .legend-item:hover {
            background: rgba(102, 126, 234, 0.05);
        }
        .legend-color {
            width: 24px; 
            height: 24px; 
            border-radius: 50%;
            box-shadow: 0 2px 5px rgba(0,0,0,0.2);
        }
        .color-bar {
            margin: 20px 0;
            padding: 20px;
            background: white;
            border-radius: 12px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.08);
        }
        .color-gradient-container {
            position: relative;
            width: 100%;
            margin: 15px 0 25px 0;  /* 底部增加空间 */
        }
        .color-gradient {
            height: 40px;
            width: 100%;
            background: linear-gradient(to right, 
                #95a5a6 0%, 
                #3498db 25%, 
                #9b59b6 50%, 
                #e74c3c 75%, 
                #c0392b 100%);
            border-radius: 20px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.15), inset 0 2px 4px rgba(255,255,255,0.3);
        }
        .color-labels {
            display: flex;
            justify-content: space-between;
            width: 100%;
            margin-top: 8px;
            font-size: 12px;
            color: #7F8C8D;
            font-weight: 600;
        }
        .color-bar-title {
            font-size: 16px;
            font-weight: bold;
            color: #2C3E50;
            margin-bottom: 10px;
        }
        .tooltip {
            position: absolute; 
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
            color: white; 
            padding: 12px 16px; 
            border-radius: 8px; 
            font-size: 13px; 
            pointer-events: none; 
            z-index: 1000; 
            display: none;
            box-shadow: 0 8px 25px rgba(0,0,0,0.3);
            border: 1px solid rgba(255,255,255,0.1);
        }
        .help-text {
            text-align: center; 
            color: #7F8C8D; 
            font-size: 13px; 
            margin: 15px 0;
            padding: 10px;
            background: rgba(102, 126, 234, 0.05);
            border-radius: 8px;
        }
        .export-controls {
            display: inline-flex;
            gap: 10px;
            align-items: center;
            background: rgba(255,255,255,0.5);
            padding: 8px 12px;
            border-radius: 10px;
        }
        .export-controls select, .export-controls input {
            color: #000000; 
            padding: 10px;
            border: 1px solid rgba(102, 126, 234, 0.3);
            border-radius: 8px;
            font-size: 14px;
            transition: all 0.3s;
        }
        .export-controls select:focus, .export-controls input:focus {
            outline: none;
            border-color: #667eea;
            box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Enhanced Protein-Ligand Contact Analysis</h1>
        <p class="subtitle">Interactive visualization with improved 3D mode and alignment</p>
        
        <div class="file-info">
            <strong>Analysis Details:</strong> Trajectory: [[trajectory_file]] | Topology: [[topology_file]] | Ligand: [[ligand_file]] ([[ligand_name]]) | Frames: [[total_frames]]
        </div>
        
        <div class="zoom-controls">
            <button class="zoom-btn" onclick="zoomIn()">+</button>
            <button class="zoom-btn" onclick="zoomOut()">−</button>
            <button class="zoom-btn" onclick="centerView()" style="width: auto; padding: 0 15px;">Center</button>
        </div>

        <div class="controls">
            <button class="control-btn" onclick="resetPositions()">Reset Layout</button>
            <button class="control-btn" onclick="rotateWheel()">Rotate Wheel</button>
            <button class="control-btn" onclick="toggleDistanceLock()" id="distanceLockBtn">🔒 Distance Locked</button>
            <button class="control-btn" onclick="toggle3DMode()" id="mode3DBtn">2D Mode</button>
            <button class="control-btn secondary" onclick="toggleConnections()" id="connectBtn">Hide Connections</button>
            <button class="control-btn secondary" onclick="toggleHydrogens()" id="hydrogenBtn">Hide H atoms</button>
            <button class="control-btn secondary" onclick="toggleGrid()" id="gridBtn">Hide Grid</button>
            <button class="control-btn" onclick="rotateCanvas180()">Rotate 180°</button>
            <div class="export-controls" style="display: inline-flex; gap: 10px; align-items: center;">
                <select id="exportQuality" class="control-btn secondary" style="padding: 10px;">
                    <option value="1">Standard (1x)</option>
                    <option value="2">High (2x)</option>
                    <option value="4">Ultra HD (4x)</option>
                    <option value="8">8K (8x)</option>
                </select>
                <input type="text" id="exportFilename" placeholder="filename" value="contact_analysis" style="padding: 8px; border: 1px solid #95A5A6; border-radius: 5px;">
                <button class="control-btn secondary" onclick="exportImage()">Export PNG</button>
            </div>
        </div>

        <p class="help-text">🖱️ In 2D: Drag residues to rotate | In 3D: Drag to rotate molecule | Middle-click/Shift+Drag to pan | Scroll to zoom</p>

        <canvas id="canvas" width="1200" height="800"></canvas>
        
        <div class="info-panel">
            <h3 style="margin-top: 0;">Contact Analysis Results</h3>
            <div class="stats">
                <div class="stat-item"><div class="stat-value">[[total_contacts]]</div><div class="stat-label">Total Contacts</div></div>
                <div class="stat-item"><div class="stat-value">[[high_freq_contacts]]</div><div class="stat-label">High Frequency</div></div>
                <div class="stat-item"><div class="stat-value">[[max_freq_percent]]%</div><div class="stat-label">Max Frequency</div></div>
                <div class="stat-item"><div class="stat-value">[[n_ligand_atoms]]</div><div class="stat-label">Ligand Atoms</div></div>
            </div>
            
            <div class="color-bar">
                <div class="color-bar-title">Contact Frequency</div>
                <div class="color-gradient-container">
                    <div class="color-gradient"></div>
                </div>
                <div class="color-labels">
                    <span>0%</span>
                    <span>25%</span>
                    <span>50%</span>
                    <span>75%</span>
                    <span>100%</span>
                </div>
            </div>
            <div class="legend">
                <div class="legend-item"><div class="legend-color" style="background: #27ae60;"></div><span>Carbon atoms</span></div>
                <div class="legend-item"><div class="legend-color" style="background: #e91e63;"></div><span>Oxygen atoms</span></div>
                <div class="legend-item"><div class="legend-color" style="background: #2196f3;"></div><span>Nitrogen atoms</span></div>
                <div class="legend-item"><div class="legend-color" style="background: #ffc107;"></div><span>Sulfur atoms</span></div>
            </div>
        </div>
    </div>

    <div class="tooltip" id="tooltip"></div>

    <script>
        const LIGAND_ATOMS = [[ligand_atoms_json]];
        const CONTACTS = [[contacts_json]];
        const LIGAND_BONDS = [[ligand_bonds_json]];"""
