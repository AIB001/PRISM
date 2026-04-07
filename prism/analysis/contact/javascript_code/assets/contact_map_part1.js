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
            }