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
        }