#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
JavaScript Code Module
Contains the complete JavaScript visualization code for the HTML
Due to size constraints, this module returns the JavaScript in parts
"""


def get_javascript_code():
    """Return the complete JavaScript code as a string"""
    return (
        get_contact_map_class_part1()
        + get_contact_map_class_part2()
        + get_drawing_methods()
        + get_utility_functions()
        + r"""
    </script>
</body>
</html>"""
    )


def get_contact_map_class_part1():
    """First part of ContactMap class - initialization and events"""
    return r"""
        class ContactMap {
            constructor() {
                try {
                    this.canvas = document.getElementById('canvas');
                    this.ctx = this.canvas.getContext('2d');
                    this.tooltip = document.getElementById('tooltip');
                    this.showConnections = true;
                    this.showHydrogens = true;
                    this.showGrid = true;
                    this.distanceLocked = true;
                    this.is3DMode = false;
                    this.isDragging = false;
                    this.isPanning = false;
                    this.dragTarget = null;
                    this.wheelRotation = 0;
                    this.offset = {x: 0, y: 0};
                    this.panStart = {x: 0, y: 0};
                    this.viewOffset = {x: 0, y: 0};
                    this.zoom = 1.0;
                    this.centerX = 600;
                    this.centerY = 400;
                    
                    // 3D rotation angles
                    this.rotationX = 0;
                    this.rotationY = 0;
                    this.rotationZ = 0;
                    
                    // Initialize data
                    this.initializeData();
                    
                    // Store initial state for reset
                    this.saveInitialState();
                    
                    this.initEvents();
                    this.animate();
                } catch(e) {
                    console.error('Error initializing ContactMap:', e);
                }
            }
            
            initializeData() {
                // Calculate ligand center
                let ligandCenterX = 0;
                let ligandCenterY = 0;
                if (LIGAND_ATOMS.length > 0) {
                    LIGAND_ATOMS.forEach(atom => {
                        ligandCenterX += atom.x;
                        ligandCenterY += atom.y;
                    });
                    ligandCenterX /= LIGAND_ATOMS.length;
                    ligandCenterY /= LIGAND_ATOMS.length;
                }
                
                // Initialize ligand atoms with proper 3D coordinates
                this.ligandAtoms = LIGAND_ATOMS.map(atom => {
                    // Use provided 3D coordinates or generate from 2D
                    let x3d = atom.x3d;
                    let y3d = atom.y3d;
                    let z3d = atom.z3d;
                    
                    // If 3D coordinates are missing or invalid, use 2D coordinates
                    if (x3d === undefined || y3d === undefined || z3d === undefined) {
                        x3d = atom.x;
                        y3d = atom.y;
                        z3d = 0;
                    }
                    
                    return {
                        ...atom,
                        x: atom.x + this.centerX,
                        y: atom.y + this.centerY,
                        x3d: x3d * 3,  // Scale up 3D coordinates for better visibility
                        y3d: y3d * 3,
                        z3d: z3d * 3,
                        fixed: true
                    };
                });
                
                // Group contacts by ligand atom to handle overlaps
                const contactsByLigandAtom = {};
                CONTACTS.forEach(contact => {
                    const atomId = contact.ligandAtom || 'L0';
                    if (!contactsByLigandAtom[atomId]) {
                        contactsByLigandAtom[atomId] = [];
                    }
                    contactsByLigandAtom[atomId].push(contact);
                });
                
                // Initialize contacts with proper alignment and spacing
                this.contacts = [];
                
                // Process each group of contacts
                Object.keys(contactsByLigandAtom).forEach(ligandAtomId => {
                    const contactGroup = contactsByLigandAtom[ligandAtomId];
                    const ligandAtom = this.ligandAtoms.find(a => a.id === ligandAtomId);
                    
                    if (ligandAtom) {
                        // Calculate vector from ligand center to this atom
                        const atomRelX = ligandAtom.x - this.centerX;
                        const atomRelY = ligandAtom.y - this.centerY;
                        
                        // Calculate base angle from center through atom
                        let baseAngle = Math.atan2(atomRelY, atomRelX);
                        const baseDirLength = Math.sqrt(atomRelX * atomRelX + atomRelY * atomRelY);
                        
                        // If atom is at center, distribute evenly
                        if (baseDirLength < 1) {
                            baseAngle = 0;
                        }
                        
                        // For multiple residues on same atom, create larger angular offset for TOP3
                        const spreadAngle = Math.PI / 6; // 30 degrees total spread
                        
                        contactGroup.forEach((contact, groupIdx) => {
                            // Calculate offset angle for multiple contacts
                            let angle = baseAngle;
                            if (contactGroup.length > 1) {
                                // Give TOP3 contacts more space
                                const spreadFactor = contact.isTop3 ? 1.5 : 1.0;
                                const offset = (groupIdx - (contactGroup.length - 1) / 2) * 
                                             (spreadAngle * spreadFactor / Math.max(1, contactGroup.length - 1));
                                angle = baseAngle + offset;
                            }
                            
                            // Adjust distance based on whether it's TOP3
                            const baseDistance = contact.pixelDistance || 250;
                            const distance = contact.isTop3 ? baseDistance * 0.7 : baseDistance * 0.5;
                            
                            // Calculate position along the extension line
                            const dirX = Math.cos(angle);
                            const dirY = Math.sin(angle);
                            
                            // Position residue along the line from ligand center through ligand atom
                            const x = ligandAtom.x + dirX * distance;
                            const y = ligandAtom.y + dirY * distance;
                            
                            this.contacts.push({
                                ...contact,
                                x: x,
                                y: y,
                                angle: angle,
                                radius: distance,
                                fixed: false,
                                initialX: x,
                                initialY: y,
                                pixelDistance: distance
                            });
                        });
                    } else {
                        // Fallback positioning
                        contactGroup.forEach((contact, idx) => {
                            const angle = (this.contacts.length * 2 * Math.PI / CONTACTS.length) - Math.PI/2;
                            const baseDistance = contact.pixelDistance || 250;
                            const distance = contact.isTop3 ? baseDistance * 0.7 : baseDistance * 0.5;
                            const x = this.centerX + Math.cos(angle) * distance;
                            const y = this.centerY + Math.sin(angle) * distance;
                            
                            this.contacts.push({
                                ...contact,
                                x: x,
                                y: y,
                                angle: angle,
                                radius: distance,
                                fixed: false,
                                initialX: x,
                                initialY: y,
                                pixelDistance: distance
                            });
                        });
                    }
                });
                
                this.ligandBonds = LIGAND_BONDS;
            }
            
            saveInitialState() {
                this.initialState = {
                    ligandAtoms: this.ligandAtoms.map(atom => ({...atom})),
                    contacts: this.contacts.map(contact => ({...contact})),
                    zoom: 1.0,
                    viewOffset: {x: 0, y: 0},
                    rotationX: 0,
                    rotationY: 0,
                    rotationZ: 0
                };
            }
            
            restoreInitialState() {
                if (!this.initialState) return;
                
                this.ligandAtoms = this.initialState.ligandAtoms.map(atom => ({...atom}));
                this.contacts = this.initialState.contacts.map(contact => ({...contact}));
                this.zoom = this.initialState.zoom;
                this.viewOffset = {...this.initialState.viewOffset};
                this.rotationX = this.initialState.rotationX;
                this.rotationY = this.initialState.rotationY;
                this.rotationZ = this.initialState.rotationZ;
                this.wheelRotation = 0;
            }"""


def get_contact_map_class_part2():
    """Second part of ContactMap class - event handlers"""
    return r"""
            
            initEvents() {
                this.canvas.addEventListener('mousedown', e => this.onMouseDown(e));
                this.canvas.addEventListener('mousemove', e => this.onMouseMove(e));
                this.canvas.addEventListener('mouseup', e => this.onMouseUp(e));
                this.canvas.addEventListener('mouseleave', e => this.onMouseLeave(e));
                this.canvas.addEventListener('wheel', e => this.onWheel(e));
                this.canvas.addEventListener('contextmenu', e => e.preventDefault());
            }
            
            getMousePos(e) {
                const rect = this.canvas.getBoundingClientRect();
                const x = (e.clientX - rect.left - this.viewOffset.x) / this.zoom;
                const y = (e.clientY - rect.top - this.viewOffset.y) / this.zoom;
                return { x, y };
            }
            
            onMouseDown(e) {
                const pos = this.getMousePos(e);
                
                if (e.button === 1 || (e.button === 0 && e.shiftKey)) {
                    // Middle click or Shift+click for panning
                    this.isPanning = true;
                    this.panStart = { x: e.clientX - this.viewOffset.x, y: e.clientY - this.viewOffset.y };
                    this.canvas.style.cursor = 'move';
                    e.preventDefault();
                } else if (e.button === 0) {
                    if (this.is3DMode) {
                        // 3D mode rotation
                        this.isDragging = true;
                        this.dragStart = { x: e.clientX, y: e.clientY };
                        this.startRotationX = this.rotationX;
                        this.startRotationY = this.rotationY;
                        this.canvas.style.cursor = 'grabbing';
                    } else {
                        // 2D mode
                        this.dragTarget = this.findElementAt(pos.x, pos.y);
                        
                        if (this.dragTarget && !this.dragTarget.fixed && this.dragTarget.radius) {
                            this.isDragging = true;
                            const ligandAtom = this.ligandAtoms.find(a => a.id === this.dragTarget.ligandAtom);
                            this.rotationCenter = ligandAtom ? { x: ligandAtom.x, y: ligandAtom.y } : { x: this.centerX, y: this.centerY };
                            
                            if (this.distanceLocked) {
                                const dx = pos.x - this.rotationCenter.x;
                                const dy = pos.y - this.rotationCenter.y;
                                this.dragStartAngle = Math.atan2(dy, dx);
                                this.targetStartAngle = Math.atan2(this.dragTarget.y - this.rotationCenter.y, 
                                                                 this.dragTarget.x - this.rotationCenter.x);
                            } else {
                                this.offset = { x: pos.x - this.dragTarget.x, y: pos.y - this.dragTarget.y };
                            }
                            this.canvas.style.cursor = 'grabbing';
                        } else if (!this.dragTarget) {
                            this.isPanning = true;
                            this.panStart = { x: e.clientX - this.viewOffset.x, y: e.clientY - this.viewOffset.y };
                            this.canvas.style.cursor = 'move';
                        }
                    }
                }
            }
            
            onMouseMove(e) {
                if (this.isPanning) {
                    this.viewOffset.x = e.clientX - this.panStart.x;
                    this.viewOffset.y = e.clientY - this.panStart.y;
                } else if (this.isDragging) {
                    if (this.is3DMode) {
                        const dx = e.clientX - this.dragStart.x;
                        const dy = e.clientY - this.dragStart.y;
                        this.rotationY = this.startRotationY + dx * 0.01;
                        this.rotationX = this.startRotationX + dy * 0.01;
                        this.rotationX = Math.max(-Math.PI, Math.min(Math.PI, this.rotationX));
                    } else if (this.dragTarget) {
                        const pos = this.getMousePos(e);
                        
                        if (this.distanceLocked) {
                            const dx = pos.x - this.rotationCenter.x;
                            const dy = pos.y - this.rotationCenter.y;
                            const currentAngle = Math.atan2(dy, dx);
                            const angleDiff = currentAngle - this.dragStartAngle;
                            
                            const dist = Math.sqrt(Math.pow(this.dragTarget.x - this.rotationCenter.x, 2) + 
                                                 Math.pow(this.dragTarget.y - this.rotationCenter.y, 2));
                            
                            const newAngle = this.targetStartAngle + angleDiff;
                            this.dragTarget.x = this.rotationCenter.x + dist * Math.cos(newAngle);
                            this.dragTarget.y = this.rotationCenter.y + dist * Math.sin(newAngle);
                            this.dragTarget.angle = newAngle;
                        } else {
                            this.dragTarget.x = pos.x - this.offset.x;
                            this.dragTarget.y = pos.y - this.offset.y;
                            
                            const dx = this.dragTarget.x - this.rotationCenter.x;
                            const dy = this.dragTarget.y - this.rotationCenter.y;
                            this.dragTarget.angle = Math.atan2(dy, dx);
                            this.dragTarget.radius = Math.sqrt(dx * dx + dy * dy);
                        }
                    }
                } else {
                    const pos = this.getMousePos(e);
                    this.handleHover(pos.x, pos.y, e);
                }
            }
            
            onMouseUp(e) {
                this.isDragging = false;
                this.isPanning = false;
                this.dragTarget = null;
                this.canvas.style.cursor = 'grab';
            }
            
            onMouseLeave(e) {
                this.onMouseUp(e);
                this.hideTooltip();
            }
            
            onWheel(e) {
                e.preventDefault();
                const delta = e.deltaY > 0 ? 0.9 : 1.1;
                const newZoom = this.zoom * delta;
                
                if (newZoom >= 0.3 && newZoom <= 3) {
                    const rect = this.canvas.getBoundingClientRect();
                    const mouseX = e.clientX - rect.left;
                    const mouseY = e.clientY - rect.top;
                    
                    const worldX = (mouseX - this.viewOffset.x) / this.zoom;
                    const worldY = (mouseY - this.viewOffset.y) / this.zoom;
                    
                    this.zoom = newZoom;
                    
                    this.viewOffset.x = mouseX - worldX * this.zoom;
                    this.viewOffset.y = mouseY - worldY * this.zoom;
                }
            }
            
            findElementAt(x, y) {
                // Check contacts first (smaller hit areas for TOP3)
                for (let contact of this.contacts) {
                    const hitRadius = contact.isTop3 ? 50 : 40;
                    const dist = Math.sqrt((x - contact.x) ** 2 + (y - contact.y) ** 2);
                    if (dist < hitRadius) return contact;
                }
                for (let atom of this.ligandAtoms) {
                    const dist = Math.sqrt((x - atom.x) ** 2 + (y - atom.y) ** 2);
                    if (dist < atom.radius) return atom;
                }
                return null;
            }
            
            handleHover(x, y, e) {
                const element = this.findElementAt(x, y);
                if (element) {
                    this.showTooltip(element, e.clientX, e.clientY);
                    this.canvas.style.cursor = element.fixed ? 'default' : 'grab';
                } else {
                    this.hideTooltip();
                    this.canvas.style.cursor = 'grab';
                }
            }
            
            showTooltip(element, x, y) {
                let content = '';
                if (element.frequency !== undefined) {
                    const distanceText = element.avgDistance ? ` | Distance: ${element.avgDistance.toFixed(2)}nm` : '';
                    const top3Text = element.isTop3 ? ' <span style="color: gold;">★ TOP 3</span>' : '';
                    content = `${element.residue}${top3Text}<br>Frequency: ${(element.frequency * 100).toFixed(1)}%${distanceText}`;
                } else {
                    content = `Ligand Atom: ${element.element}${element.id.substring(1)}<br>Contacts: ${element.contacts || 0}`;
                }
                this.tooltip.innerHTML = content;
                this.tooltip.style.display = 'block';
                this.tooltip.style.left = (x + 10) + 'px';
                this.tooltip.style.top = (y - 10) + 'px';
            }
            
            hideTooltip() {
                this.tooltip.style.display = 'none';
            }"""


def get_drawing_methods():
    """Return complete drawing methods for ContactMap class"""
    return r"""
            
            rotate3D(x, y, z, rotX, rotY, rotZ) {
                let x1 = x, y1 = y, z1 = z;
                const cosX = Math.cos(rotX);
                const sinX = Math.sin(rotX);
                const y2 = y1 * cosX - z1 * sinX;
                const z2 = y1 * sinX + z1 * cosX;
                const cosY = Math.cos(rotY);
                const sinY = Math.sin(rotY);
                const x3 = x1 * cosY + z2 * sinY;
                const z3 = -x1 * sinY + z2 * cosY;
                const cosZ = Math.cos(rotZ);
                const sinZ = Math.sin(rotZ);
                const x4 = x3 * cosZ - y2 * sinZ;
                const y4 = x3 * sinZ + y2 * cosZ;
                return { x: x4, y: y4, z: z3 };
            }
            
            project3D(x, y, z) {
                const focalLength = 800;
                const scale = focalLength / (focalLength + z);
                return {
                    x: x * scale,
                    y: y * scale,
                    scale: scale,
                    depth: z
                };
            }
            
            draw3DScene() {
                const elements3D = [];
                this.ligandAtoms.forEach(atom => {
                    const rotated = this.rotate3D(
                        atom.x3d * 2, atom.y3d * 2, atom.z3d * 2,
                        this.rotationX, this.rotationY, this.rotationZ
                    );
                    const projected = this.project3D(rotated.x, rotated.y, rotated.z);
                    elements3D.push({
                        type: 'atom',
                        data: atom,
                        x: this.centerX + projected.x,
                        y: this.centerY + projected.y,
                        z: projected.depth,
                        scale: projected.scale
                    });
                });
                
                this.ligandBonds.forEach(([id1, id2]) => {
                    const atom1 = this.ligandAtoms.find(a => a.id === id1);
                    const atom2 = this.ligandAtoms.find(a => a.id === id2);
                    if (atom1 && atom2) {
                        const rot1 = this.rotate3D(
                            atom1.x3d * 2, atom1.y3d * 2, atom1.z3d * 2,
                            this.rotationX, this.rotationY, this.rotationZ
                        );
                        const rot2 = this.rotate3D(
                            atom2.x3d * 2, atom2.y3d * 2, atom2.z3d * 2,
                            this.rotationX, this.rotationY, this.rotationZ
                        );
                        const proj1 = this.project3D(rot1.x, rot1.y, rot1.z);
                        const proj2 = this.project3D(rot2.x, rot2.y, rot2.z);
                        elements3D.push({
                            type: 'bond',
                            x1: this.centerX + proj1.x,
                            y1: this.centerY + proj1.y,
                            x2: this.centerX + proj2.x,
                            y2: this.centerY + proj2.y,
                            z: (proj1.depth + proj2.depth) / 2,
                            scale: (proj1.scale + proj2.scale) / 2,
                            atom1: atom1,
                            atom2: atom2
                        });
                    }
                });
                
                elements3D.sort((a, b) => b.z - a.z);
                
                elements3D.forEach(elem => {
                    if (elem.type === 'bond') {
                        this.draw3DBond(elem);
                    } else if (elem.type === 'atom') {
                        this.draw3DAtom(elem);
                    }
                });
                
                this.draw3DAxisIndicator();
            }
            
            draw3DBond(bond) {
                const width = 3 * bond.scale;
                const opacity = 0.3 + bond.scale * 0.5;
                this.ctx.globalAlpha = opacity;
                this.ctx.strokeStyle = '#34495e';
                this.ctx.lineWidth = width;
                if (bond.atom1.element === 'H' || bond.atom2.element === 'H') {
                    this.ctx.lineWidth = width * 0.5;
                }
                this.ctx.beginPath();
                this.ctx.moveTo(bond.x1, bond.y1);
                this.ctx.lineTo(bond.x2, bond.y2);
                this.ctx.stroke();
                this.ctx.globalAlpha = 1.0;
            }
            
            draw3DAtom(elem) {
                const atom = elem.data;
                const radius = atom.radius * elem.scale;
                const opacity = 0.7 + elem.scale * 0.3;
                
                this.ctx.globalAlpha = opacity;
                this.ctx.fillStyle = this.getAtomColor(atom.element);
                this.ctx.beginPath();
                this.ctx.arc(elem.x, elem.y, radius, 0, 2 * Math.PI);
                this.ctx.fill();
                
                const gradient = this.ctx.createRadialGradient(
                    elem.x - radius * 0.3, elem.y - radius * 0.3, 0,
                    elem.x, elem.y, radius
                );
                gradient.addColorStop(0, 'rgba(255, 255, 255, 0.5)');
                gradient.addColorStop(0.5, 'rgba(255, 255, 255, 0.1)');
                gradient.addColorStop(1, 'rgba(0, 0, 0, 0.4)');
                this.ctx.fillStyle = gradient;
                this.ctx.fill();
                
                this.ctx.strokeStyle = elem.z > 0 ? 'rgba(255,255,255,0.8)' : 'rgba(0,0,0,0.3)';
                this.ctx.lineWidth = atom.element === 'H' ? 1 : 2;
                this.ctx.stroke();
                
                if (atom.element !== 'H' || elem.scale > 0.7) {
                    this.ctx.globalAlpha = opacity;
                    this.ctx.fillStyle = atom.element === 'H' ? '#333333' : '#ffffff';
                    this.ctx.font = `bold ${Math.floor(12 * elem.scale)}px Arial`;
                    this.ctx.textAlign = 'center';
                    this.ctx.textBaseline = 'middle';
                    const label = atom.element + atom.id.substring(1);
                    this.ctx.fillText(label, elem.x, elem.y);
                }
                this.ctx.globalAlpha = 1.0;
            }
            
            draw3DAxisIndicator() {
                const x = 80;
                const y = 80;
                const length = 40;
                const axes = [
                    { x: 1, y: 0, z: 0, color: '#e74c3c', label: 'X' },
                    { x: 0, y: 1, z: 0, color: '#2ecc71', label: 'Y' },
                    { x: 0, y: 0, z: 1, color: '#3498db', label: 'Z' }
                ];
                
                axes.forEach(axis => {
                    const rotated = this.rotate3D(
                        axis.x * length, axis.y * length, axis.z * length,
                        this.rotationX, this.rotationY, this.rotationZ
                    );
                    this.ctx.strokeStyle = axis.color;
                    this.ctx.lineWidth = 2;
                    this.ctx.beginPath();
                    this.ctx.moveTo(x, y);
                    this.ctx.lineTo(x + rotated.x, y + rotated.y);
                    this.ctx.stroke();
                    this.ctx.fillStyle = axis.color;
                    this.ctx.font = 'bold 12px Arial';
                    this.ctx.textAlign = 'center';
                    this.ctx.textBaseline = 'middle';
                    this.ctx.fillText(axis.label, x + rotated.x * 1.2, y + rotated.y * 1.2);
                });
                
                this.ctx.fillStyle = '#34495e';
                this.ctx.beginPath();
                this.ctx.arc(x, y, 3, 0, 2 * Math.PI);
                this.ctx.fill();
            }
            
            getAtomColor(element) {
                const colors = {
                    'C': '#27ae60',
                    'N': '#2196f3',
                    'O': '#e91e63',
                    'S': '#ffc107',
                    'H': '#ffffff',
                    'P': '#ff5722',
                    'F': '#9c27b0',
                    'Cl': '#4caf50',
                    'Br': '#795548',
                    'I': '#607d8b'
                };
                return colors[element] || '#7f8c8d';
            }
            
            getFrequencyColor(freq) {
                if (freq <= 0.3) {
                    const t = freq / 0.3;
                    const r = Math.floor(149 + (52 - 149) * t);
                    const g = Math.floor(165 + (152 - 165) * t);
                    const b = Math.floor(166 + (219 - 166) * t);
                    return `rgb(${r}, ${g}, ${b})`;
                } else if (freq <= 0.6) {
                    const t = (freq - 0.3) / 0.3;
                    const r = Math.floor(52 + (142 - 52) * t);
                    const g = Math.floor(152 + (68 - 152) * t);
                    const b = Math.floor(219 + (173 - 219) * t);
                    return `rgb(${r}, ${g}, ${b})`;
                } else {
                    const t = (freq - 0.6) / 0.4;
                    const r = Math.floor(142 + (231 - 142) * t);
                    const g = Math.floor(68 + (76 - 68) * t);
                    const b = Math.floor(173 + (60 - 173) * t);
                    return `rgb(${r}, ${g}, ${b})`;
                }
            }
            
            drawGrid() {
                this.ctx.strokeStyle = 'rgba(236, 240, 241, 0.3)';
                this.ctx.lineWidth = 1;
                const step = 50;
                const startX = -this.viewOffset.x / this.zoom;
                const startY = -this.viewOffset.y / this.zoom;
                const endX = (this.canvas.width - this.viewOffset.x) / this.zoom;
                const endY = (this.canvas.height - this.viewOffset.y) / this.zoom;
                
                for (let x = Math.floor(startX / step) * step; x <= endX; x += step) {
                    this.ctx.beginPath();
                    this.ctx.moveTo(x, startY);
                    this.ctx.lineTo(x, endY);
                    this.ctx.stroke();
                }
                for (let y = Math.floor(startY / step) * step; y <= endY; y += step) {
                    this.ctx.beginPath();
                    this.ctx.moveTo(startX, y);
                    this.ctx.lineTo(endX, y);
                    this.ctx.stroke();
                }
            }
            
            drawConnections() {
                this.contacts.forEach((contact, index) => {
                    const ligandAtom = this.ligandAtoms.find(a => a.id === contact.ligandAtom);
                    if (!ligandAtom) return;
                    
                    const color = this.getFrequencyColor(contact.frequency);
                    
                    let targetX = contact.x;
                    let targetY = contact.y;
                    
                    // Adjust target for TOP3 contacts
                    if (contact.isTop3) {
                        if (contact.contactAtomX !== undefined && contact.contactAtomY !== undefined) {
                            targetX = contact.contactAtomX;
                            targetY = contact.contactAtomY;
                        } else {
                            const angle = contact.angle || Math.atan2(ligandAtom.y - contact.y, ligandAtom.x - contact.x);
                            targetX = contact.x + Math.cos(angle) * 25;  // Reduced from 30
                            targetY = contact.y + Math.sin(angle) * 25;
                        }
                    }
                    
                    // Draw enhanced connections for TOP3
                    if (contact.isTop3) {
                        // Draw glow effect
                        this.ctx.save();
                        this.ctx.shadowBlur = 15;
                        this.ctx.shadowColor = color;
                        this.ctx.strokeStyle = color;
                        this.ctx.lineWidth = 4 + contact.frequency * 4;
                        this.ctx.globalAlpha = 0.3;
                        this.ctx.beginPath();
                        this.ctx.moveTo(ligandAtom.x, ligandAtom.y);
                        this.ctx.lineTo(targetX, targetY);
                        this.ctx.stroke();
                        this.ctx.restore();
                    }
                    
                    // Draw main connection
                    this.ctx.strokeStyle = color;
                    this.ctx.lineWidth = contact.isTop3 ? (3 + contact.frequency * 3) : (1.5 + contact.frequency * 2);
                    this.ctx.globalAlpha = contact.isTop3 ? 0.8 : (0.5 + contact.frequency * 0.3);
                    this.ctx.beginPath();
                    this.ctx.moveTo(ligandAtom.x, ligandAtom.y);
                    this.ctx.lineTo(targetX, targetY);
                    this.ctx.stroke();
                    this.ctx.globalAlpha = 1.0;
                });
            }
            
            drawLigand() {
                this.ctx.strokeStyle = '#34495e';
                this.ctx.lineWidth = 3;
                
                // 绘制键
                this.ligandBonds.forEach(([id1, id2]) => {
                    const atom1 = this.ligandAtoms.find(a => a.id === id1);
                    const atom2 = this.ligandAtoms.find(a => a.id === id2);
                    if (atom1 && atom2) {
                        // 如果隐藏氢原子模式开启，跳过涉及氢的键
                        if (!this.showHydrogens && (atom1.element === 'H' || atom2.element === 'H')) {
                            return;
                        }
                        
                        if (atom1.element === 'H' || atom2.element === 'H') {
                            this.ctx.lineWidth = 1.5;
                        } else {
                            this.ctx.lineWidth = 3;
                        }
                        this.ctx.beginPath();
                        this.ctx.moveTo(atom1.x, atom1.y);
                        this.ctx.lineTo(atom2.x, atom2.y);
                        this.ctx.stroke();
                    }
                });
                
                // 绘制原子
                this.ligandAtoms.forEach(atom => {
                    // 如果隐藏氢原子模式开启，跳过氢原子
                    if (!this.showHydrogens && atom.element === 'H') {
                        return;
                    }
                    
                    const color = this.getAtomColor(atom.element);
                    
                    this.ctx.fillStyle = color;
                    this.ctx.beginPath();
                    this.ctx.arc(atom.x, atom.y, atom.radius, 0, 2 * Math.PI);
                    this.ctx.fill();
                    
                    this.ctx.strokeStyle = atom.element === 'H' ? '#cccccc' : '#ffffff';
                    this.ctx.lineWidth = atom.element === 'H' ? 1 : 2;
                    this.ctx.stroke();
                    
                    const atomNumber = atom.id.substring(1);
                    let label = atom.element;
                    
                    if (atom.element === 'H') {
                        label = 'H' + atomNumber;
                        this.ctx.fillStyle = '#333333';
                        this.ctx.font = 'bold 9px Arial';
                    } else {
                        label = atom.element + atomNumber;
                        this.ctx.fillStyle = '#ffffff';
                        this.ctx.font = 'bold 12px Arial';
                    }
                    
                    this.ctx.textAlign = 'center';
                    this.ctx.textBaseline = 'middle';
                    this.ctx.fillText(label, atom.x, atom.y);
                });
            }
            
            drawContacts() {
                this.contacts.forEach((contact, i) => {
                    const color = this.getFrequencyColor(contact.frequency);
                    this.drawOrientedStructure(contact, color, i + 1);
                });
            }
            
            drawOrientedStructure(contact, color, index) {
                const x = contact.x, y = contact.y;
                const residueType = contact.residueType;
                
                const ligandAtom = this.ligandAtoms.find(a => a.id === contact.ligandAtom);
                if (!ligandAtom) {
                    this.drawSimpleLabel(x, y, residueType, contact.residue, color, index);
                    return;
                }
                
                const dx = ligandAtom.x - x;
                const dy = ligandAtom.y - y;
                const angle = Math.atan2(dy, dx);
                
                // Use contact.isTop3 instead of index
                if (contact.isTop3) {
                    this.drawOrientedAAStructure(x, y, residueType, contact.residue, color, index, angle, contact.frequency);
                } else {
                    this.drawSimpleLabel(x, y, residueType, contact.residue, color, index);
                }
            }
            
            drawOrientedAAStructure(x, y, residueType, residueName, color, index, angle, frequency) {
                // Calculate contact atom position (closer to structure)
                const contactAtomPos = {
                    x: x + Math.cos(angle) * 25,  // Reduced from 30
                    y: y + Math.sin(angle) * 25
                };
                
                const contactIdx = this.contacts.findIndex(c => c.residue === residueName);
                if (contactIdx >= 0) {
                    this.contacts[contactIdx].contactAtomX = contactAtomPos.x;
                    this.contacts[contactIdx].contactAtomY = contactAtomPos.y;
                }
                
                this.ctx.save();
                this.ctx.translate(x, y);
                this.ctx.rotate(angle);
                
                // Smaller backbone structure
                const backbone = [
                    {name: 'N', x: 22, y: 0, color: '#4169E1', radius: 11},    // Reduced sizes
                    {name: 'Ca', x: 7, y: 0, color: '#808080', radius: 11},
                    {name: 'C', x: -7, y: 0, color: '#808080', radius: 11},
                    {name: 'O1', x: -22, y: -7, color: '#DC143C', radius: 9},
                    {name: 'O2', x: -22, y: 7, color: '#DC143C', radius: 9}
                ];
                
                // Draw bonds with thinner lines
                this.ctx.strokeStyle = '#000000';
                this.ctx.lineWidth = 4;  // Reduced from 5
                
                this.ctx.beginPath();
                this.ctx.moveTo(backbone[0].x - backbone[0].radius * 0.8, backbone[0].y);
                this.ctx.lineTo(backbone[1].x + backbone[1].radius * 0.8, backbone[1].y);
                this.ctx.stroke();
                
                this.ctx.beginPath();
                this.ctx.moveTo(backbone[1].x - backbone[1].radius * 0.8, backbone[1].y);
                this.ctx.lineTo(backbone[2].x + backbone[2].radius * 0.8, backbone[2].y);
                this.ctx.stroke();
                
                this.ctx.lineWidth = 3;  // Reduced from 4
                this.ctx.beginPath();
                this.ctx.moveTo(backbone[2].x - backbone[2].radius * 0.7, backbone[2].y - 4);
                this.ctx.lineTo(backbone[3].x + backbone[3].radius * 0.7, backbone[3].y);
                this.ctx.stroke();
                
                // Draw side chain (smaller)
                this.drawOrientedSideChain(backbone[1].x, backbone[1].y, residueType);
                
                // Draw atoms
                backbone.forEach(atom => {
                    // Add subtle glow for high frequency contacts
                    if (frequency > 0.7) {
                        this.ctx.save();
                        this.ctx.shadowBlur = 10;
                        this.ctx.shadowColor = atom.color;
                        this.ctx.fillStyle = atom.color;
                        this.ctx.beginPath();
                        this.ctx.arc(atom.x, atom.y, atom.radius, 0, 2 * Math.PI);
                        this.ctx.fill();
                        this.ctx.restore();
                    }
                    
                    this.ctx.fillStyle = atom.color;
                    this.ctx.strokeStyle = 'black';
                    this.ctx.lineWidth = 1.5;  // Reduced from 2
                    
                    this.ctx.beginPath();
                    this.ctx.arc(atom.x, atom.y, atom.radius, 0, 2 * Math.PI);
                    this.ctx.fill();
                    this.ctx.stroke();
                    
                    if (atom.name === 'N' || atom.name.startsWith('O')) {
                        this.ctx.fillStyle = 'white';
                        this.ctx.font = 'bold 13px Arial';  // Reduced from 16px
                        this.ctx.textAlign = 'center';
                        this.ctx.textBaseline = 'middle';
                        const label = atom.name.startsWith('O') ? 'O' : atom.name;
                        this.ctx.fillText(label, atom.x, atom.y);
                    }
                });
                
                this.ctx.restore();
                
                // Draw residue name
                this.ctx.fillStyle = 'darkred';
                this.ctx.font = 'bold 20px Arial';  // Reduced from 26px
                this.ctx.textAlign = 'center';
                this.ctx.textBaseline = 'bottom';
                this.ctx.fillText(residueName, x, y - 38);  // Adjusted position
                
                // Draw rank badge
                this.ctx.fillStyle = '#FFD700';
                this.ctx.strokeStyle = '#000000';
                this.ctx.lineWidth = 1.5;
                this.ctx.beginPath();
                this.ctx.arc(x + 40, y - 25, 14, 0, 2 * Math.PI);  // Smaller badge
                this.ctx.fill();
                this.ctx.stroke();
                this.ctx.fillStyle = 'black';
                this.ctx.font = 'bold 14px Arial';  // Reduced from 18px
                this.ctx.textAlign = 'center';
                this.ctx.textBaseline = 'middle';
                // 计算TOP3排名 - 使用residueName参数
                const top3Contacts = this.contacts.filter(c => c.isTop3).sort((a, b) => b.frequency - a.frequency);
                const rankIndex = top3Contacts.findIndex(c => c.residue === residueName) + 1;
                this.ctx.fillText('#' + rankIndex, x + 40, y - 25);
            }
            
            drawOrientedSideChain(x, y, residueType) {
                this.ctx.strokeStyle = '#000000';
                this.ctx.lineWidth = 3;  // Reduced from 4
                
                if (residueType === 'GLY') {
                    this.ctx.fillStyle = '#808080';
                    this.ctx.font = 'bold 12px Arial';  // Reduced from 14px
                    this.ctx.textAlign = 'center';
                    this.ctx.fillText('H', x, y + 20);  // Adjusted position
                } else {
                    this.ctx.beginPath();
                    this.ctx.moveTo(x, y + 11);  // Adjusted positions
                    this.ctx.lineTo(x, y + 28);
                    this.ctx.stroke();
                    
                    this.ctx.fillStyle = '#808080';
                    this.ctx.strokeStyle = 'black';
                    this.ctx.lineWidth = 1.5;  // Reduced from 2
                    this.ctx.beginPath();
                    this.ctx.arc(x, y + 28, 8, 0, 2 * Math.PI);  // Smaller R group
                    this.ctx.fill();
                    this.ctx.stroke();
                    
                    this.ctx.fillStyle = 'white';
                    this.ctx.font = 'bold 10px Arial';  // Reduced from 12px
                    this.ctx.textAlign = 'center';
                    this.ctx.fillText('R', x, y + 28);
                }
            }
            
            drawSimpleLabel(x, y, residueType, residueName, color, index) {
                const label = residueName;
                const radius = 28;  // Slightly smaller
                
                // Add subtle shadow for depth
                this.ctx.save();
                this.ctx.shadowBlur = 5;
                this.ctx.shadowColor = 'rgba(0,0,0,0.2)';
                this.ctx.shadowOffsetX = 2;
                this.ctx.shadowOffsetY = 2;
                
                this.ctx.fillStyle = 'white';
                this.ctx.strokeStyle = color;  // Use frequency color instead of fixed color
                this.ctx.lineWidth = 3;
                this.ctx.beginPath();
                this.ctx.arc(x, y, radius, 0, 2 * Math.PI);
                this.ctx.fill();
                this.ctx.stroke();
                
                this.ctx.restore();
                
                this.ctx.fillStyle = 'darkblue';
                this.ctx.font = 'bold 13px Arial';  // Slightly smaller
                this.ctx.textAlign = 'center';
                this.ctx.textBaseline = 'middle';
                this.ctx.fillText(label, x, y);
            }
            
            draw() {
                this.ctx.save();
                this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

                this.ctx.translate(this.viewOffset.x, this.viewOffset.y);
                this.ctx.scale(this.zoom, this.zoom);

                if (this.is3DMode) {
                    if (this.showGrid) this.drawGrid();
                    this.draw3DScene();
                } else {
                    if (this.showGrid) this.drawGrid();
                    if (this.showConnections) this.drawConnections();
                    this.drawLigand();
                    this.drawContacts();
                }

                this.ctx.restore();
            }
            
            animate() {
                this.draw();
                requestAnimationFrame(() => this.animate());
            }
            
            resetPositions() {
                this.restoreInitialState();
            }
            
            rotateWheel() {
                const rotationStep = Math.PI / 6;
                const steps = 20;
                let currentStep = 0;
                
                const animate = () => {
                    if (currentStep < steps) {
                        const delta = rotationStep / steps;
                        this.wheelRotation += delta;
                        
                        this.contacts.forEach(contact => {
                            const ligandAtom = this.ligandAtoms.find(a => a.id === contact.ligandAtom);
                            const centerX = ligandAtom ? ligandAtom.x : this.centerX;
                            const centerY = ligandAtom ? ligandAtom.y : this.centerY;
                            
                            const dx = contact.x - centerX;
                            const dy = contact.y - centerY;
                            const dist = Math.sqrt(dx * dx + dy * dy);
                            const currentAngle = Math.atan2(dy, dx);
                            const newAngle = currentAngle + delta;
                            
                            contact.x = centerX + dist * Math.cos(newAngle);
                            contact.y = centerY + dist * Math.sin(newAngle);
                            contact.angle = newAngle;
                        });
                        
                        currentStep++;
                        setTimeout(animate, 20);
                    }
                };
                animate();
            }
        }"""


def get_utility_functions():
    """Return utility functions JavaScript code"""
    return r"""
        
        let contactMap;
        window.addEventListener('DOMContentLoaded', () => { 
            contactMap = new ContactMap(); 
        });
        
        function toggle3DMode() {
            if (!contactMap.is3DMode) {
                contactMap.is3DMode = true;
                const btn = document.getElementById('mode3DBtn');
                btn.textContent = '3D Mode';
                btn.classList.add('active');
                contactMap.rotationX = 0;
                contactMap.rotationY = 0;
                contactMap.rotationZ = 0;
            } else {
                contactMap.is3DMode = false;
                const btn = document.getElementById('mode3DBtn');
                btn.textContent = '2D Mode';
                btn.classList.remove('active');
                contactMap.restoreInitialState();
            }
        }
        
        function toggleDistanceLock() {
            contactMap.distanceLocked = !contactMap.distanceLocked;
            const btn = document.getElementById('distanceLockBtn');
            if (contactMap.distanceLocked) {
                btn.textContent = '🔒 Distance Locked';
                btn.classList.remove('secondary');
            } else {
                btn.textContent = '🔓 Distance Free';
                btn.classList.add('secondary');
            }
        }
        
        function resetPositions() { contactMap.resetPositions(); }
        function rotateWheel() { contactMap.rotateWheel(); }
        function toggleConnections() { 
            contactMap.showConnections = !contactMap.showConnections; 
            const btn = document.getElementById('connectBtn');
            btn.textContent = contactMap.showConnections ? 'Hide Connections' : 'Show Connections';
        }
        function toggleHydrogens() {
            contactMap.showHydrogens = !contactMap.showHydrogens;
            const btn = document.getElementById('hydrogenBtn');
            btn.textContent = contactMap.showHydrogens ? 'Hide H atoms' : 'Show H atoms';
        }

        function toggleGrid() {
            contactMap.showGrid = !contactMap.showGrid;
            const btn = document.getElementById('gridBtn');
            btn.textContent = contactMap.showGrid ? 'Hide Grid' : 'Show Grid';
        }

        function rotateCanvas180() {
            // Animate 180-degree rotation
            const steps = 30;
            let currentStep = 0;
            const rotationPerStep = Math.PI / steps;

            const animate = () => {
                if (currentStep < steps) {
                    // Rotate all contacts around the center
                    contactMap.contacts.forEach(contact => {
                        const dx = contact.x - contactMap.centerX;
                        const dy = contact.y - contactMap.centerY;
                        const dist = Math.sqrt(dx * dx + dy * dy);
                        const currentAngle = Math.atan2(dy, dx);
                        const newAngle = currentAngle + rotationPerStep;

                        contact.x = contactMap.centerX + dist * Math.cos(newAngle);
                        contact.y = contactMap.centerY + dist * Math.sin(newAngle);
                        contact.angle = newAngle;
                    });

                    // Rotate ligand atoms around the center
                    contactMap.ligandAtoms.forEach(atom => {
                        const dx = atom.x - contactMap.centerX;
                        const dy = atom.y - contactMap.centerY;
                        const dist = Math.sqrt(dx * dx + dy * dy);
                        const currentAngle = Math.atan2(dy, dx);
                        const newAngle = currentAngle + rotationPerStep;

                        atom.x = contactMap.centerX + dist * Math.cos(newAngle);
                        atom.y = contactMap.centerY + dist * Math.sin(newAngle);
                    });

                    currentStep++;
                    setTimeout(animate, 16); // ~60fps
                }
            };
            animate();
        }

        function centerView() { 
            // Calculate the bounding box of all visible elements
            let minX = Infinity, maxX = -Infinity;
            let minY = Infinity, maxY = -Infinity;
            
            // Include ligand atoms
            contactMap.ligandAtoms.forEach(atom => {
                minX = Math.min(minX, atom.x);
                maxX = Math.max(maxX, atom.x);
                minY = Math.min(minY, atom.y);
                maxY = Math.max(maxY, atom.y);
            });
            
            // Include contacts (only if visible)
            if (!contactMap.is3DMode) {
                contactMap.contacts.forEach(contact => {
                    minX = Math.min(minX, contact.x);
                    maxX = Math.max(maxX, contact.x);
                    minY = Math.min(minY, contact.y);
                    maxY = Math.max(maxY, contact.y);
                });
            }
            
            // Calculate center and size of content
            const contentWidth = maxX - minX;
            const contentHeight = maxY - minY;
            const contentCenterX = (minX + maxX) / 2;
            const contentCenterY = (minY + maxY) / 2;
            
            // Calculate zoom to fit content
            const canvasWidth = contactMap.canvas.width;
            const canvasHeight = contactMap.canvas.height;
            const padding = 50; // Padding around content
            
            const zoomX = (canvasWidth - padding * 2) / contentWidth;
            const zoomY = (canvasHeight - padding * 2) / contentHeight;
            contactMap.zoom = Math.min(zoomX, zoomY, 2); // Cap at 2x zoom
            
            // Center the view
            contactMap.viewOffset.x = canvasWidth / 2 - contentCenterX * contactMap.zoom;
            contactMap.viewOffset.y = canvasHeight / 2 - contentCenterY * contactMap.zoom;
        }
        
        function exportImage() {
            const quality = parseInt(document.getElementById('exportQuality').value);
            const filename = document.getElementById('exportFilename').value || 'contact_analysis';
            
            // Create high-resolution canvas
            const hdCanvas = document.createElement('canvas');
            hdCanvas.width = contactMap.canvas.width * quality;
            hdCanvas.height = contactMap.canvas.height * quality;
            const hdCtx = hdCanvas.getContext('2d');
            
            // Save current state
            const originalZoom = contactMap.zoom;
            const originalOffset = {...contactMap.viewOffset};
            
            // Scale up for high quality
            hdCtx.scale(quality, quality);
            
            // Temporarily replace context
            const originalCtx = contactMap.ctx;
            contactMap.ctx = hdCtx;
            
            // Enable high quality rendering
            contactMap.ctx.imageSmoothingEnabled = true;
            contactMap.ctx.imageSmoothingQuality = 'high';
            
            // Draw the scene
            contactMap.draw();
            
            // Restore original context
            contactMap.ctx = originalCtx;
            
            // Create download link
            const link = document.createElement('a');
            link.download = `${filename}_${quality}x.png`;
            link.href = hdCanvas.toDataURL('image/png', 1.0);
            link.click();
            
            // Clean up
            hdCanvas.remove();
        }
        
        function zoomIn() { 
            const newZoom = Math.min(contactMap.zoom * 1.2, 3);
            const rect = contactMap.canvas.getBoundingClientRect();
            const centerX = rect.width / 2;
            const centerY = rect.height / 2;
            
            const worldX = (centerX - contactMap.viewOffset.x) / contactMap.zoom;
            const worldY = (centerY - contactMap.viewOffset.y) / contactMap.zoom;
            
            contactMap.zoom = newZoom;
            
            contactMap.viewOffset.x = centerX - worldX * contactMap.zoom;
            contactMap.viewOffset.y = centerY - worldY * contactMap.zoom;
        }
        
        function zoomOut() { 
            const newZoom = Math.max(contactMap.zoom * 0.8, 0.3);
            const rect = contactMap.canvas.getBoundingClientRect();
            const centerX = rect.width / 2;
            const centerY = rect.height / 2;
            
            const worldX = (centerX - contactMap.viewOffset.x) / contactMap.zoom;
            const worldY = (centerY - contactMap.viewOffset.y) / contactMap.zoom;
            
            contactMap.zoom = newZoom;
            
            contactMap.viewOffset.x = centerX - worldX * contactMap.zoom;
            contactMap.viewOffset.y = centerY - worldY * contactMap.zoom;
        }"""
