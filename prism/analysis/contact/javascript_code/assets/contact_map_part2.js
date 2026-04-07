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
            }