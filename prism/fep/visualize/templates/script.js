        // Canvas setup
        const canvas = document.getElementById('canvas');
        const ctx = canvas.getContext('2d');
        const tooltip = document.getElementById('tooltip');

        // Canvas dimensions
        const canvasWidth = 1400;
        const canvasHeight = 700;

        // Layout: each molecule gets half the canvas
        const halfWidth = canvasWidth / 2;
        const centerX_A = halfWidth / 2;     // 350 - center of left half
        const centerX_B = halfWidth + halfWidth / 2;  // 1050 - center of right half
        const centerY = canvasHeight / 2;   // 350

        // View state - INDEPENDENT for each molecule
        let molAState = {
            viewOffset: { x: 0, y: 0 },
            zoom: 1.0,
            isDragging: false,
            dragStart: { x: 0, y: 0 }
        };
        let molBState = {
            viewOffset: { x: 0, y: 0 },
            zoom: 1.0,
            isDragging: false,
            dragStart: { x: 0, y: 0 }
        };
        let showLabels = true;
        let showCharges = false;
        let colorMode = 'fep';
        let hoveredAtom = null;
        let hoveredMolecule = null;
        let activeMolecule = null;  // Which molecule is being interacted with

        // Calculate initial zoom for a single molecule
        function calculateInitialZoomForMolecule(atoms) {
            let maxX = 0, maxY = 0;
            atoms.forEach(atom => {
                maxX = Math.max(maxX, Math.abs(atom.x));
                maxY = Math.max(maxY, Math.abs(atom.y));
            });

            const padding = 50;
            const availableWidth = halfWidth - padding * 2;
            const availableHeight = canvasHeight - padding * 2;

            const zoomX = availableWidth / (maxX * 2 + 0.1);
            const zoomY = availableHeight / (maxY * 2 + 0.1);

            return Math.min(zoomX, zoomY, 1.5);
        }

        // Initialize zoom for each molecule independently
        molAState.zoom = calculateInitialZoomForMolecule(ATOMS_A);
        molBState.zoom = calculateInitialZoomForMolecule(ATOMS_B);

        // Get color based on current mode
        function getAtomFillColor(atom) {
            if (colorMode === 'fep') {
                return atom.fepColor;
            } else {
                return atom.elementColor;
            }
        }

        // Draw bonds
        function drawBonds(atoms, bonds, offsetX, offsetY) {
            bonds.forEach(bond => {
                // 兼容旧格式：数组 [id1, id2]
                const atom1Id = bond.atom1 || bond[0];
                const atom2Id = bond.atom2 || bond[1];
                const bondType = bond.type || 1;  // 默认单键

                const atom1 = atoms.find(a => a.id === atom1Id);
                const atom2 = atoms.find(a => a.id === atom2Id);

                if (atom1 && atom2) {
                    ctx.beginPath();
                    ctx.moveTo(atom1.x, atom1.y);
                    ctx.lineTo(atom2.x, atom2.y);

                    // 根据键序设置样式
                    if (bondType === 2) {
                        ctx.lineWidth = 4;  // DOUBLE: 粗线
                    } else if (bondType === 3) {
                        ctx.lineWidth = 6;  // TRIPLE: 更粗
                    } else if (bondType === 12) {
                        ctx.lineWidth = 2;
                        ctx.setLineDash([4, 4]);  // AROMATIC: 虚线
                    } else {
                        ctx.lineWidth = 2;  // SINGLE: 默认
                    }

                    ctx.strokeStyle = '#666';
                    ctx.stroke();
                    ctx.setLineDash([]);  // 重置
                }
            });
        }

        // Draw atoms
        function drawAtoms(atoms, offsetX, offsetY, moleculeId) {
            atoms.forEach((atom, index) => {
                const isHovered = hoveredAtom === index && hoveredMolecule === moleculeId;

                // Draw atom circle
                ctx.beginPath();
                ctx.arc(atom.x, atom.y, atom.radius, 0, 2 * Math.PI);

                // Fill color based on mode
                ctx.fillStyle = getAtomFillColor(atom);
                ctx.fill();

                // Highlight border if hovered
                if (isHovered) {
                    ctx.strokeStyle = '#FFD700';
                    ctx.lineWidth = 3;
                    ctx.stroke();
                } else {
                    // Border always uses element color
                    ctx.strokeStyle = atom.elementColor;
                    ctx.lineWidth = 1.5;
                    ctx.stroke();
                }

                // Draw label
                if (showLabels) {
                    ctx.fillStyle = (atom.element === 'H' || colorMode === 'fep') ? '#333' : '#fff';
                    ctx.font = 'bold 11px Arial';
                    ctx.textAlign = 'center';
                    ctx.textBaseline = 'middle';
                    ctx.fillText(atom.name, atom.x, atom.y);
                }

                // Draw charge if requested
                if (showCharges) {
                    ctx.fillStyle = '#333';
                    ctx.font = '9px Arial';
                    ctx.fillText(atom.charge.toFixed(4), atom.x, atom.y + atom.radius + 12);
                }
            });
        }

        // Main draw function - uses Canvas transforms
        function draw() {
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
        }

        // Find atom at position - converts mouse to world coordinates for each molecule
        function findAtomAtPosition(mouseX, mouseY) {
            // Check molecule A with its independent transform
            const worldX_A = (mouseX - molAState.viewOffset.x - centerX_A) / molAState.zoom;
            const worldY_A = (mouseY - molAState.viewOffset.y - centerY) / molAState.zoom;

            for (let i = 0; i < ATOMS_A.length; i++) {
                const atom = ATOMS_A[i];
                const dx = worldX_A - atom.x;
                const dy = worldY_A - atom.y;
                if (Math.sqrt(dx * dx + dy * dy) < atom.radius) {
                    return { atom, index: i, molecule: 'a' };
                }
            }

            // Check molecule B with its independent transform
            const worldX_B = (mouseX - molBState.viewOffset.x - centerX_B) / molBState.zoom;
            const worldY_B = (mouseY - molBState.viewOffset.y - centerY) / molBState.zoom;

            for (let i = 0; i < ATOMS_B.length; i++) {
                const atom = ATOMS_B[i];
                const dx = worldX_B - atom.x;
                const dy = worldY_B - atom.y;
                if (Math.sqrt(dx * dx + dy * dy) < atom.radius) {
                    return { atom, index: i, molecule: 'b' };
                }
            }

            return null;
        }

        // Show tooltip
        function showTooltip(atom, clientX, clientY, molecule, index) {
            const classificationLabels = {
                'common': 'Common',
                'transformed': 'Transformed',
                'surrounding': 'Surrounding'
            };

            let content = `
                <div class="tooltip-title">${atom.name}</div>
                <div class="tooltip-row">
                    <span class="tooltip-label">Element:</span>
                    <span>${atom.element}</span>
                </div>
                <div class="tooltip-row">
                    <span class="tooltip-label">Charge:</span>
                    <span>${atom.charge.toFixed(4)}</span>
                </div>
                <div class="tooltip-row">
                    <span class="tooltip-label">Classification:</span>
                    <span>${classificationLabels[atom.classification] || atom.classification}</span>
                </div>
            `;

            // Add correspondence info
            const key = `${molecule}_${index}`;
            if (CORRESPONDENCE[key]) {
                const parts = CORRESPONDENCE[key].split('_');
                const targetMol = parts[0];
                const targetIdx = parseInt(parts[1]);
                const targetAtoms = targetMol === 'a' ? ATOMS_A : ATOMS_B;
                const targetAtom = targetAtoms[targetIdx];

                if (targetAtom) {
                    content += `
                        <div class="correspondence">
                            ↔ Corresponds to: ${targetAtom.name} (${targetAtom.element}, charge: ${targetAtom.charge.toFixed(4)})
                        </div>
                    `;
                }
            }

            tooltip.innerHTML = content;
            tooltip.style.display = 'block';
            tooltip.style.left = (clientX + 15) + 'px';
            tooltip.style.top = (clientY - 10) + 'px';
        }

        // Hide tooltip
        function hideTooltip() {
            tooltip.style.display = 'none';
        }

        // Mouse events - INDEPENDENT for each molecule
        canvas.addEventListener('mousedown', (e) => {
            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;

            // Determine which molecule was clicked based on X position
            if (mouseX < halfWidth) {
                activeMolecule = 'a';
                molAState.isDragging = true;
                molAState.dragStart = { x: e.clientX - molAState.viewOffset.x, y: e.clientY - molAState.viewOffset.y };
            } else {
                activeMolecule = 'b';
                molBState.isDragging = true;
                molBState.dragStart = { x: e.clientX - molBState.viewOffset.x, y: e.clientY - molBState.viewOffset.y };
            }
            canvas.style.cursor = 'grabbing';
        });

        canvas.addEventListener('mousemove', (e) => {
            // Handle dragging for active molecule
            if (molAState.isDragging && activeMolecule === 'a') {
                molAState.viewOffset.x = e.clientX - molAState.dragStart.x;
                molAState.viewOffset.y = e.clientY - molAState.dragStart.y;
                draw();
            } else if (molBState.isDragging && activeMolecule === 'b') {
                molBState.viewOffset.x = e.clientX - molBState.dragStart.x;
                molBState.viewOffset.y = e.clientY - molBState.dragStart.y;
                draw();
            } else {
                // Handle hover for both molecules
                const rect = canvas.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;

                const found = findAtomAtPosition(x, y);

                if (found) {
                    if (hoveredAtom !== found.index || hoveredMolecule !== found.molecule) {
                        hoveredAtom = found.index;
                        hoveredMolecule = found.molecule;
                        showTooltip(found.atom, e.clientX, e.clientY, found.molecule, found.index);
                        draw();
                    }
                } else {
                    if (hoveredAtom !== null) {
                        hoveredAtom = null;
                        hoveredMolecule = null;
                        hideTooltip();
                        draw();
                    }
                }
            }
        });

        canvas.addEventListener('mouseup', () => {
            molAState.isDragging = false;
            molBState.isDragging = false;
            activeMolecule = null;
            canvas.style.cursor = 'grab';
        });

        // Zoom - INDEPENDENT for each molecule (like prism/analysis/contact)
        canvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            const delta = e.deltaY > 0 ? 0.9 : 1.1;

            // Determine which molecule to zoom based on mouse X position
            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;

            let targetState;
            let centerX;
            if (mouseX < halfWidth) {
                targetState = molAState;
                centerX = centerX_A;
            } else {
                targetState = molBState;
                centerX = centerX_B;
            }

            const newZoom = targetState.zoom * delta;

            if (newZoom >= 0.3 && newZoom <= 3) {
                const mouseY = e.clientY - rect.top;

                // Convert to world coordinates before zoom
                const worldX = (mouseX - targetState.viewOffset.x - centerX) / targetState.zoom;
                const worldY = (mouseY - targetState.viewOffset.y - centerY) / targetState.zoom;

                targetState.zoom = newZoom;

                // Adjust view offset to zoom towards mouse position
                targetState.viewOffset.x = mouseX - worldX * targetState.zoom - centerX;
                targetState.viewOffset.y = mouseY - worldY * targetState.zoom - centerY;

                draw();
            }
        });
        // Toggle handlers
        document.getElementById('toggle-charges').addEventListener('change', (e) => {
            showCharges = e.target.checked;
            draw();  // Redraw to show/hide charge labels on canvas
            draw();
        });

        document.getElementById('toggle-labels').addEventListener('change', (e) => {
            showLabels = e.target.checked;
            draw();
        });

        // Color mode handlers
        document.querySelectorAll('input[name="colorMode"]').forEach(radio => {
            radio.addEventListener('change', (e) => {
                colorMode = e.target.value;
                switchLegendTab(colorMode);  // Auto-switch legend to match coloring mode
                draw();
            });
        });

        // Render atom list
        function renderAtomList() {
            renderAtomTable();
        }

        function renderAtomTable() {
            const tbody = document.getElementById('atom-table-body');
            const classificationLabels = {
                'common': 'Common',
                'transformed': 'Transformed',
                'surrounding': 'Surrounding'
            };

            let html = '';

            // Create a map of all atoms with their pairs
            const atomPairs = [];

            // Add transformed A atoms (no pair in B)
            ATOMS_A.forEach(atomA => {
                if (atomA.classification === 'transformed') {
                    atomPairs.push({ atomA, atomB: null, type: 'transformed_a' });
                }
            });

            // Add transformed B atoms (no pair in A)
            ATOMS_B.forEach(atomB => {
                if (atomB.classification === 'transformed') {
                    atomPairs.push({ atomA: null, atomB, type: 'transformed_b' });
                }
            });

            // Add surrounding atom pairs
            ATOMS_A.forEach((atomA, indexA) => {
                if (atomA.classification === 'surrounding') {
                    const keyA = `a_${indexA}`;
                    const correspondingKey = CORRESPONDENCE[keyA];
                    if (correspondingKey) {
                        const parts = correspondingKey.split('_');
                        const indexB = parseInt(parts[1]);
                        const atomB = ATOMS_B[indexB];
                        if (atomB) {
                            atomPairs.push({ atomA, atomB, type: 'surrounding' });
                        }
                    }
                }
            });

            // Add common atom pairs
            ATOMS_A.forEach((atomA, indexA) => {
                if (atomA.classification === 'common') {
                    const keyA = `a_${indexA}`;
                    const correspondingKey = CORRESPONDENCE[keyA];
                    if (correspondingKey) {
                        const parts = correspondingKey.split('_');
                        const indexB = parseInt(parts[1]);
                        const atomB = ATOMS_B[indexB];
                        if (atomB) {
                            atomPairs.push({ atomA, atomB, type: 'common' });
                        }
                    }
                }
            });

            // Render table rows
            atomPairs.forEach(pair => {
                const { atomA, atomB, type } = pair;

                html += '<tr>';

                // Ligand A columns
                if (atomA) {
                    const badgeClass = 'badge-' + atomA.classification;
                    html += `
                        <td class="atom-name">${atomA.name}</td>
                        <td>${atomA.element}</td>
                        <td>${atomA.type}</td>
                        <td class="atom-charge">${atomA.charge.toFixed(4)}</td>
                        <td><span class="classification-badge ${badgeClass}">${classificationLabels[atomA.classification]}</span></td>
                    `;
                } else {
                    html += '<td>—</td><td>—</td><td>—</td><td>—</td><td><span class="classification-badge badge-transformed">—</span></td>';
                }

                // Ligand B columns
                if (atomB) {
                    const badgeClass = 'badge-' + atomB.classification;
                    html += `
                        <td class="atom-name" style="border-left: 2px solid #ddd;">${atomB.name}</td>
                        <td>${atomB.element}</td>
                        <td>${atomB.type}</td>
                        <td class="atom-charge">${atomB.charge.toFixed(4)}</td>
                        <td><span class="classification-badge ${badgeClass}">${classificationLabels[atomB.classification]}</span></td>
                    `;
                } else {
                    html += '<td style="border-left: 2px solid #ddd;">—</td><td>—</td><td>—</td><td>—</td><td><span class="classification-badge badge-transformed">—</span></td>';
                }

                html += '</tr>';
            });

            tbody.innerHTML = html;
        }

        // Reset view - reset BOTH molecules independently
        function resetView() {
            molAState.viewOffset = { x: 0, y: 0 };
            molAState.zoom = calculateInitialZoomForMolecule(ATOMS_A);
            molBState.viewOffset = { x: 0, y: 0 };
            molBState.zoom = calculateInitialZoomForMolecule(ATOMS_B);
            hoveredAtom = null;
            hoveredMolecule = null;
            activeMolecule = null;
            hideTooltip();
            draw();
        }

        // Export PNG
        function exportPNG() {
            const link = document.createElement('a');
            link.download = 'fep_mapping.png';
            link.href = canvas.toDataURL();
            link.click();
        }

        function toggleConfig() {
            const content = document.getElementById('config-content');
            const icon = document.getElementById('config-toggle');
            if (content && icon) {
                content.classList.toggle('collapsed');
                icon.classList.toggle('collapsed');
            }
        }

        function scrollToAtomDetails() {
            const atomList = document.getElementById('atom-list');
            if (atomList) {
                atomList.scrollIntoView({ behavior: 'smooth', block: 'start' });
            }
        }

        function switchLegendTab(tabName) {
            // Hide all tab contents
            const contents = document.querySelectorAll('.legend-tab-content');
            contents.forEach(content => content.classList.remove('active'));

            // Show selected tab content
            const selectedContent = document.getElementById('legend-' + tabName);
            if (selectedContent) {
                selectedContent.classList.add('active');
            }
        }

        // Toggle Mapping/Build Log panel
        const logToggle = document.getElementById('log-toggle');
        const logContent = document.getElementById('log-content');
        const logToggleIcon = document.getElementById('log-toggle-icon');
        if (logToggle && logContent && logToggleIcon) {
            logToggle.addEventListener('click', () => {
                logContent.classList.toggle('collapsed');
                logToggleIcon.textContent = logContent.classList.contains('collapsed') ? '▶' : '▼';
            });
        }

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
