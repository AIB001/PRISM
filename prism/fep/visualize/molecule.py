#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Molecule processing utilities for FEP visualization.

Handles conversion from PDB files to RDKit Mol objects with proper
atom labeling, bond order correction, and charge annotation.
"""

from typing import List, Optional
from pathlib import Path
import re

# Constants for coordinate matching
COORDINATE_MATCH_TOLERANCE_ANGSTROM = 1.0  # Max distance (Å) for atom coordinate matching


def _normalize_atom_name(name: str) -> str:
    """Return a normalized atom name for cross-force-field matching."""
    cleaned = re.sub(r"[^A-Za-z0-9]", "", name.strip())
    return re.sub(r"[AB]$", "", cleaned)


def _assign_atom_metadata(rdkit_atom: "Chem.Atom", atom) -> None:
    """Copy name and charge metadata from a topology atom onto an RDKit atom."""
    rdkit_atom.SetProp("atomLabel", atom.name)
    rdkit_atom.SetProp("name", atom.name)
    rdkit_atom.SetProp("atomNote", f"{atom.charge:+.4f}")


def _create_mol_from_atoms(atoms: List) -> "Chem.Mol":
    """
    Create RDKit Mol object from a list of Atom objects.

    This function creates an RDKit molecule directly from atoms, bypassing
    PDB file parsing. This is necessary for hybrid topology atoms whose names
    (N00A, C01A, etc.) are not in standard PDB format.

    Parameters
    ----------
    atoms : List[Atom]
        List of Atom objects with name, element, coord, etc.

    Returns
    -------
    Chem.Mol
        RDKit Mol object with atoms and coordinates
    """
    from rdkit import Chem

    # Create empty molecule
    mol = Chem.RWMol()

    # Add atoms
    for atom in atoms:
        rdkit_atom = Chem.Atom(atom.element)  # Create atom from element symbol
        mol.AddAtom(rdkit_atom)

    # Set conformer (coordinates)
    conf = Chem.Conformer(len(atoms))
    for i, atom in enumerate(atoms):
        conf.SetAtomPosition(i, (atom.coord[0], atom.coord[1], atom.coord[2]))
    mol.AddConformer(conf)

    # Add bonds based on distance (simple heuristic)
    # This will be overwritten by assign_bond_orders_from_mol2 if MOL2 is available
    _add_bonds_by_distance(mol, atoms)

    # Return as regular Mol
    return mol.GetMol()


def _add_bonds_by_distance(mol: "Chem.RWMol", atoms: List, tolerance: float = 1.7) -> None:
    """
    Add bonds to mol based on distance between atoms.

    Simple heuristic: two atoms are bonded if their distance is less than tolerance.
    This is a fallback for when MOL2 file is not available.

    Parameters
    ----------
    mol : Chem.RWMol
        RDKit molecule to add bonds to
    atoms : List[Atom]
        List of Atom objects with coordinates
    tolerance : float
        Maximum distance for bond (default 1.7 Å)
    """
    from rdkit import Chem

    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            coord_i = atoms[i].coord
            coord_j = atoms[j].coord
            dist = (
                (coord_i[0] - coord_j[0]) ** 2 + (coord_i[1] - coord_j[1]) ** 2 + (coord_i[2] - coord_j[2]) ** 2
            ) ** 0.5

            if dist < tolerance:
                mol.AddBond(i, j, Chem.BondType.SINGLE)


def _assign_bond_orders_from_mol2_to_atoms_mol(
    target_mol: "Chem.Mol",
    mol2_file: str,
    coord_threshold: Optional[float] = None,
) -> "Chem.Mol":
    """
    Transfer bond orders from a MOL2 template onto a molecule created from the
    mapped atom list.

    The target molecule already has the correct atom names and receptor-aligned
    coordinates. We only need the bond topology from the MOL2 file. Matching is
    done by element + 3D coordinates, which is robust across force fields where
    atom names may differ.
    """
    from rdkit import Chem
    import numpy as np

    if coord_threshold is None:
        coord_threshold = COORDINATE_MATCH_TOLERANCE_ANGSTROM

    mol2_mol = Chem.MolFromMol2File(mol2_file, sanitize=False, removeHs=False)
    if mol2_mol is None:
        raise ValueError(f"RDKit could not parse MOL2: {mol2_file}")

    target_rw = Chem.RWMol(target_mol)
    # Remove any heuristic bonds before adding template bonds.
    while target_rw.GetNumBonds():
        bond = target_rw.GetBondWithIdx(0)
        target_rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    target_conf = target_rw.GetConformer()
    mol2_conf = mol2_mol.GetConformer()

    template_to_target = {}
    used_target = set()
    for mol2_idx in range(mol2_mol.GetNumAtoms()):
        mol2_atom = mol2_mol.GetAtomWithIdx(mol2_idx)
        mol2_elem = mol2_atom.GetSymbol()
        mol2_pos = np.array(mol2_conf.GetAtomPosition(mol2_idx))

        best_idx = None
        best_dist = float("inf")
        for target_idx in range(target_rw.GetNumAtoms()):
            if target_idx in used_target:
                continue
            target_atom = target_rw.GetAtomWithIdx(target_idx)
            if target_atom.GetSymbol() != mol2_elem:
                continue
            target_pos = np.array(target_conf.GetAtomPosition(target_idx))
            dist = np.linalg.norm(mol2_pos - target_pos)
            if dist < coord_threshold and dist < best_dist:
                best_idx = target_idx
                best_dist = dist

        if best_idx is not None:
            template_to_target[mol2_idx] = best_idx
            used_target.add(best_idx)

    for bond in mol2_mol.GetBonds():
        begin = template_to_target.get(bond.GetBeginAtomIdx())
        end = template_to_target.get(bond.GetEndAtomIdx())
        if begin is None or end is None or begin == end:
            continue
        if target_rw.GetBondBetweenAtoms(begin, end) is None:
            target_rw.AddBond(begin, end, bond.GetBondType())

    return target_rw.GetMol()


def pdb_to_mol(pdb_file: str) -> "Chem.Mol":
    """
    Convert PDB file to RDKit Mol object with proper labeling.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file

    Returns
    -------
    Chem.Mol
        RDKit Mol object with:
        - Atom name labels from PDB
        - Inferred bond orders and aromaticity
        - 2D coordinates for depiction

    Examples
    --------
    >>> mol = pdb_to_mol('ligand.pdb')
    >>> from rdkit import Chem
    >>> Chem.MolToSmiles(mol)  # Get SMILES
    """
    try:
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise ImportError(
            "RDKit is required for FEP visualization. " "Install with: conda install -c conda-forge rdkit"
        ) from e

    if not Path(pdb_file).exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    # Read PDB file
    # IMPORTANT: Use removeHs=False to keep all atoms (including explicit hydrogens)
    # CHARMM-GUI PDB files already have all hydrogens explicitly defined
    # We do NOT call Chem.AddHs() because it adds extra implicit hydrogens
    # that create gray atoms in visualization (names like "Atom32", "Atom33")
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if mol is None:
        raise ValueError(f"Failed to read PDB file: {pdb_file}")

    # DO NOT add hydrogens - CHARMM-GUI files already have all hydrogens
    # Chem.AddHs() adds implicit hydrogens that don't match the mapping data
    # if not remove_hs:
    #     mol = Chem.AddHs(mol)  # <-- REMOVED to fix gray atom bug

    # Sanitize: infer bond orders, aromaticity, hybridization
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # If sanitization fails, try kekulization only
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            pass  # Proceed with unsanitized mol

    # Add atom name labels from PDB
    for atom in mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info:
            name = pdb_info.GetName().strip()
            atom.SetProp("atomLabel", name)
            atom.SetProp("name", name)

    # Generate 2D coordinates for depiction
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        # If 2D coord generation fails, keep 3D coords
        pass

    return mol


def assign_bond_orders_from_mol2(
    pdb_file: str,
    mol2_file: str,
) -> "Chem.Mol":
    """
    Assign bond orders to PDB molecule using mol2 file as template.

    Uses RDKit's MolFromMol2File to read correct bond orders from mol2,
    then assigns those bond orders to the PDB structure.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file with coordinates
    mol2_file : str
        Path to mol2 file with correct bond orders

    Returns
    -------
    Chem.Mol
        RDKit Mol object with corrected bond orders and stereochemistry

    Examples
    --------
    >>> mol = assign_bond_orders_from_mol2('ligand.pdb', 'ligand.mol2')
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdmolops
        from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
    except ImportError as e:
        raise ImportError("RDKit is required. Install with: conda install -c conda-forge rdkit") from e

    if not Path(mol2_file).exists():
        raise FileNotFoundError(f"Mol2 file not found: {mol2_file}")
    if not Path(pdb_file).exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    # Read mol2 file to get bond order template
    try:
        mol2_template = Chem.MolFromMol2File(mol2_file, sanitize=True, removeHs=False)
    except Exception as e:
        raise ValueError(f"Failed to read mol2 file {mol2_file}: {e}")

    # Read PDB file with coordinates
    # Use sanitize=False for hybrid topologies with unusual valence states
    pdb_mol = Chem.MolFromPDBFile(pdb_file, sanitize=False, removeHs=False)
    if pdb_mol is None:
        raise ValueError(f"Failed to read PDB file: {pdb_file}")

    # Assign bond orders from mol2 template to PDB structure
    try:
        pdb_mol = AssignBondOrdersFromTemplate(mol2_template, pdb_mol)
    except Exception as e:
        raise ValueError(f"Failed to assign bond orders from template: {e}")

    # Assign stereochemistry from 3D coordinates
    try:
        rdmolops.AssignStereochemistryFrom3D(pdb_mol)
    except Exception:
        pass  # Stereochemistry assignment may fail, continue

    return pdb_mol


def assign_bond_orders_from_mol2_by_coords(
    pdb_file: str,
    mol2_file: str,
    coord_threshold: Optional[float] = None,
) -> "Chem.Mol":
    """
    Assign bond orders to PDB molecule using mol2 file as template with coordinate-based matching.

    This function is designed for cases where atom names differ between PDB and MOL2 files
    (e.g., OPLS-AA force field uses different naming than original MOL2).

    Strategy:
    1. Read the entire molecule from MOL2 - this already has correct bond orders
    2. Match atoms from MOL2 to PDB by 3D coordinates (coordinates are identical or very close)
    3. Transfer the coordinates from PDB to MOL2 molecule
    4. This preserves all bond order information from MOL2

    This is the most reliable approach when PDB coordinates are receptor-aligned
    but the original MOL2 already contains correct bond order information.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file with (receptor-aligned) coordinates
    mol2_file : str
        Path to mol2 file with correct bond orders
    coord_threshold : float, optional
        Distance threshold for coordinate matching.
        If not provided, uses the global `COORDINATE_MATCH_TOLERANCE_ANGSTROM`.

    Returns
    -------
    Chem.Mol
        RDKit Mol object with corrected bond orders and stereochemistry

    Examples
    --------
    >>> mol = assign_bond_orders_from_mol2_by_coords('ligand.pdb', 'ligand.mol2')
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdmolops
        import numpy as np
    except ImportError as e:
        raise ImportError("RDKit is required. Install with: conda install -c conda-forge rdkit") from e

    if not Path(mol2_file).exists():
        raise FileNotFoundError(f"Mol2 file not found: {mol2_file}")
    if not Path(pdb_file).exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    # Use default tolerance if not specified
    if coord_threshold is None:
        coord_threshold = COORDINATE_MATCH_TOLERANCE_ANGSTROM

    # Read MOL2 - this already has correct bond orders
    # Use sanitize=False because we do it later anyway, and sanitize can remove bonds
    # that RDKit considers invalid but are actually valid for the hybrid topology
    try:
        mol = Chem.MolFromMol2File(mol2_file, sanitize=False, removeHs=False)
    except Exception as e:
        raise ValueError(f"Failed to read mol2 file {mol2_file}: {e}")

    if mol is None:
        raise ValueError(f"RDKit could not parse MOL2: {mol2_file}")

    # Read PDB to get receptor-aligned coordinates
    pdb_mol = Chem.MolFromPDBFile(pdb_file, sanitize=False, removeHs=False)
    if pdb_mol is None:
        raise ValueError(f"Failed to read PDB file: {pdb_file}")

    mol_count = mol.GetNumAtoms()
    pdb_count = pdb_mol.GetNumAtoms()

    print(f"  MOL2 atoms: {mol_count}, PDB atoms: {pdb_count}")

    mol_conf = mol.GetConformer()
    pdb_conf = pdb_mol.GetConformer()

    matched = 0
    unmatched = 0

    # Always use coordinate-based matching
    # Even when atom counts match, the order may differ (PDB from hybrid topology,
    # MOL2 from original ligand). Coordinate matching is more robust.
    print(f"  Using coordinate matching with threshold {coord_threshold} Å")
    for mol_idx in range(mol.GetNumAtoms()):
        mol_atom = mol.GetAtomWithIdx(mol_idx)
        mol_elem = mol_atom.GetSymbol()
        mol_pos = np.array(mol_conf.GetAtomPosition(mol_idx))

        # Find matching atom in PDB by element and coordinate
        best_match = None
        best_dist = float("inf")

        for pdb_idx in range(pdb_mol.GetNumAtoms()):
            pdb_atom = pdb_mol.GetAtomWithIdx(pdb_idx)
            pdb_elem = pdb_atom.GetSymbol()

            if mol_elem != pdb_elem:
                continue

            pdb_pos = np.array(pdb_conf.GetAtomPosition(pdb_idx))
            dist = np.linalg.norm(mol_pos - pdb_pos)

            if dist < best_dist and dist < coord_threshold:
                best_dist = dist
                best_match = pdb_idx

        if best_match is not None:
            # Update MOL2 atom position with PDB (receptor-aligned) coordinates
            pdb_pos = pdb_conf.GetAtomPosition(best_match)
            mol_conf.SetAtomPosition(mol_idx, pdb_pos)
            matched += 1
        else:
            unmatched += 1

    if unmatched > 0:
        print(f"  Warning: {unmatched} MOL2 atoms could not be matched to PDB coordinates")

    print(f"  ✓ Matched {matched}/{mol_count} atoms, transferred PDB coordinates to MOL2")

    # The MOL2 already has correct bond orders - we just transferred PDB coordinates.
    # But we need to copy the PDB atom names (which match the hybrid topology naming)
    # back to the molecule so that downstream label matching works correctly.
    # PDB atom names are like C0A, C02, N00A (for hybrid topology), while MOL2 has original names.
    for mol_idx in range(mol.GetNumAtoms()):
        if mol_idx < pdb_mol.GetNumAtoms():
            pdb_atom = pdb_mol.GetAtomWithIdx(mol_idx)
            pdb_info = pdb_atom.GetPDBResidueInfo()
            if pdb_info:
                mol_atom = mol.GetAtomWithIdx(mol_idx)
                # Copy the entire PDB residue info object - it already contains
                # the correct atom name that matches the hybrid topology atom list
                mol_atom.SetPDBResidueInfo(pdb_info)

    # The MOL2 already has correct bond orders - just assign stereochemistry
    try:
        rdmolops.AssignStereochemistryFrom3D(mol)
    except Exception:
        pass  # Stereochemistry assignment may fail, continue

    return mol


def auto_detect_mol2(pdb_file: str, search_dir: Optional[str] = None, charmm_gui: bool = False) -> Optional[str]:
    """
    Automatically detect a MOL2 file corresponding to the given PDB file.

    Search strategy:
    1. A same-stem ``.mol2`` in the current directory
    2. ``basename + '_3D.mol2'`` in the current directory (only when
       ``charmm_gui=True``)
    3. Special basename handling such as ``39_matched -> 39`` (only when
       ``charmm_gui=True``)
    4. A same-stem file in the parent directory
    5. A matching file in the optional search directory

    Strategy 2/3 (`_3D.mol2` naming) is a CHARMM-GUI convention. For GAFF/OpenFF/OPLS
    workflows the user should name the MOL2 file with the same stem as the PDB (Strategy 1)
    so that bond order templates are found reliably without picking up unrelated files.

    Parameters
    ----------
    pdb_file : str
        PDB file path
    search_dir : str, optional
        Optional extra search directory
    charmm_gui : bool, optional
        Enable CHARMM-GUI/CGenFF search strategies (``_3D.mol2`` pattern).
        Default False.

    Returns
    -------
    str or None
        Path to the detected MOL2 file, or ``None`` if no match is found

    Examples
    --------
    >>> auto_detect_mol2("/path/to/42.pdb")           # finds 42.mol2 only
    >>> auto_detect_mol2("/path/to/39.pdb", charmm_gui=True)  # also finds 39_3D.mol2
    >>> auto_detect_mol2("/path/to/output/39_matched.pdb", charmm_gui=True)
    "/path/to/39_3D.mol2"
    """
    from pathlib import Path

    pdb_path = Path(pdb_file)
    search_paths = []
    parent_dir = pdb_path.parent

    # Strategy 1: same-stem .mol2 in the current directory.
    search_paths.append(pdb_path.with_suffix(".mol2"))

    if charmm_gui:
        # Strategy 2: basename_3D.mol2 in the current directory (CHARMM-GUI convention).
        search_paths.append(parent_dir / f"{pdb_path.stem}_3D.mol2")

        # Strategy 3: strip common suffixes such as "_matched" and search again.
        basename = pdb_path.stem
        for suffix in ["_matched", "_aligned", "_processed", "_generated"]:
            if basename.endswith(suffix):
                base_stem = basename[: -len(suffix)]
                # Search both the current directory and its parent.
                search_paths.append(parent_dir / f"{base_stem}.mol2")
                search_paths.append(parent_dir / f"{base_stem}_3D.mol2")
                search_paths.append(parent_dir.parent / f"{base_stem}.mol2")
                search_paths.append(parent_dir.parent / f"{base_stem}_3D.mol2")
                break

    # Strategy 4: search for a same-stem file in the parent directory.
    if parent_dir != parent_dir.parent:
        search_paths.append(parent_dir.parent / f"{pdb_path.stem}.mol2")
        if charmm_gui:
            search_paths.append(parent_dir.parent / f"{pdb_path.stem}_3D.mol2")

    # Strategy 5: search the optional extra directory.
    if search_dir:
        search_dir = Path(search_dir)
        search_paths.append(search_dir / f"{pdb_path.stem}.mol2")
        if charmm_gui:
            search_paths.append(search_dir / f"{pdb_path.stem}_3D.mol2")

    # Return the first existing candidate.
    for path in search_paths:
        if path.exists():
            return str(path)

    return None


def prepare_mol_with_charges_and_labels(
    pdb_file: str,
    mol2_file: Optional[str],
    atoms: List,  # List of Atom objects from prism.fep.core.mapping
    charmm_gui: bool = False,
) -> "Chem.Mol":
    """
    Prepare RDKit Mol object with correct bond orders, atom names and charge labels.

    Universal strategy for ALL force fields (GAFF2/OPLS/OpenFF/CGenFF):
    - Always use PDB for coordinates (receptor-aligned)
    - Use MOL2 as bond order template only (via RDKit's MCS matching)
    - Distance matching works because PDB coords == atoms coords

    This solves the coordinate mismatch problem where MOL2 files contain
    parameterization-generated coordinates that differ from receptor-aligned
    PDB coordinates.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file with receptor-aligned coordinates
    mol2_file : str, optional
        Path to mol2 file with correct bond orders (template only)
    atoms : List[Atom]
        List of Atom objects with receptor-aligned coordinates
    charmm_gui : bool, optional
        Enable CHARMM-GUI/CGenFF MOL2 auto-detection strategies (``_3D.mol2``).
        Default False.

    Returns
    -------
    Chem.Mol
        RDKit Mol object with atom labels and charge notes
    """
    try:
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise ImportError("RDKit is required. Install with: conda install -c conda-forge rdkit") from e

    # Step 0: Auto-detect MOL2 file if not provided
    if mol2_file is None:
        detected_mol2 = auto_detect_mol2(pdb_file, charmm_gui=charmm_gui)
        if detected_mol2:
            mol2_file = detected_mol2

    # Step 1: Create RDKit Mol object from atoms list (NOT from PDB file)
    # This is necessary because hybrid topology atom names (N00A, C01A, etc.)
    # are not in standard PDB format and RDKit cannot parse them correctly
    #
    # We create a mol object manually from the atoms list, which already has
    # the correct atom names, elements, and coordinates
    try:
        mol = _create_mol_from_atoms(atoms)
        print(f"✓ Created RDKit Mol from {len(atoms)} atoms")
    except Exception as e:
        raise ValueError(f"Failed to create RDKit Mol from atoms list: {e}")

    # Step 2: Assign bond orders from MOL2 template (if available).
    # This uses the visualization PDB plus a MOL2 template and preserves the
    # topology-derived atom ordering/labels from the PDB.
    if mol2_file is not None and Path(mol2_file).exists():
        try:
            mol = _assign_bond_orders_from_mol2_to_atoms_mol(mol, mol2_file)
            print(f"✓ Assigned bond orders from MOL2 template: {mol2_file}")
        except Exception as e:
            error_msg = str(e)
            # If the error is "No matching found", try coordinate-based matching
            if "No matching found" in error_msg:
                print(f"Warning: Standard bond order assignment failed: {e}")
                print(f"  Trying coordinate-based bond order matching...")
                try:
                    mol = _assign_bond_orders_from_mol2_to_atoms_mol(mol, mol2_file)
                    print(f"✓ Assigned bond orders from MOL2 using coordinate matching: {mol2_file}")
                except Exception as e2:
                    print(f"Warning: Coordinate-based bond order assignment also failed: {e2}")
                    print(f"  Falling back to heuristic bond inference")
                    # Fall through - PDB mol already has coordinates
                    pass
            else:
                print(f"Warning: Failed to assign bond orders from MOL2: {e}")
                print(f"  Falling back to heuristic bond inference")
                # Fall through - PDB mol already has coordinates
                pass

    # Step 3: Sanitize molecule (with timeout for hybrid topologies)
    # Hybrid topologies may have unusual valence states that cause infinite loops
    # Skip sanitization for hybrid topologies to avoid hanging
    # We only need the molecule for visualization, not for chemical calculations

    # Skip sanitization entirely for hybrid topologies
    # The PDB already has coordinates and bond orders from the topology
    # Sanitization can hang on molecules with unusual valence states

    # DON'T call Kekulize - it clears aromatic flags and converts aromatic bonds
    # to single/double alternating pattern. We want to preserve AROMATIC bond type.
    # The MOL2 already has correct bond orders including aromaticity.

    print("Note: Skipping full sanitization for hybrid topology (visualization only)")

    # Step 4: Assign atom metadata directly (mol created from atoms, so order matches)
    # Since we created the mol object directly from the atoms list, the atom order
    # should match exactly. We can assign metadata directly without complex matching.
    for mol_idx, atom in enumerate(atoms):
        if mol_idx >= mol.GetNumAtoms():
            print(f"Warning: Atom index {mol_idx} out of range in RDKit mol")
            break

        rdkit_atom = mol.GetAtomWithIdx(mol_idx)
        _assign_atom_metadata(rdkit_atom, atom)

    # Step 5: Generate 2D coordinates for depiction
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass  # Keep 3D coords if 2D generation fails

    return mol
