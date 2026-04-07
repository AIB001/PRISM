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
        }