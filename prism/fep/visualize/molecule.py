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
        from rdkit import Chem
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
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise ImportError("RDKit is required. Install with: conda install -c conda-forge rdkit") from e

    # Step 0: Auto-detect MOL2 file if not provided
    if mol2_file is None:
        detected_mol2 = auto_detect_mol2(pdb_file, charmm_gui=charmm_gui)
        if detected_mol2:
            mol2_file = detected_mol2

    # Step 1: Read PDB for coordinates (ALWAYS - receptor-aligned)
    # This is the KEY fix: we always use PDB coordinates, not MOL2 coordinates
    # sanitize=False is critical for hybrid topologies with unusual valence states
    # flavor=0 ensures RDKit preserves PDB atom names exactly as written (4-char names like N00A)
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=False, removeHs=False, flavor=0)
    if mol is None:
        raise ValueError(f"Failed to read PDB: {pdb_file}")

    # Step 2: Assign bond orders from MOL2 template (if available)
    # Uses RDKit's AssignBondOrdersFromTemplate which does MCS matching
    # This is coordinate-independent and works for all force fields
    if mol2_file is not None and Path(mol2_file).exists():
        try:
            mol = assign_bond_orders_from_mol2(pdb_file, mol2_file)
            print(f"✓ Assigned bond orders from MOL2 template: {mol2_file}")
        except Exception as e:
            error_msg = str(e)
            # If the error is "No matching found", try coordinate-based matching
            if "No matching found" in error_msg:
                print(f"Warning: Standard bond order assignment failed: {e}")
                print(f"  Trying coordinate-based bond order matching...")
                try:
                    mol = assign_bond_orders_from_mol2_by_coords(pdb_file, mol2_file)
                    print(f"✓ Assigned bond orders from MOL2 using coordinate matching: {mol2_file}")
                except Exception as e2:
                    print(f"Warning: Coordinate-based bond order assignment also failed: {e2}")
                    print(f"  Falling back to RDKit-inferred bond orders from PDB")
                    # Fall through - PDB mol already has coordinates
                    pass
            else:
                print(f"Warning: Failed to assign bond orders from MOL2: {e}")
                print(f"  Falling back to RDKit-inferred bond orders from PDB")
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

    # Step 4: Match atoms and assign labels
    charge_dict = {atom.name: atom.charge for atom in atoms}
    atom_coords = [atom.coord for atom in atoms]
    atoms_by_name = {atom.name: idx for idx, atom in enumerate(atoms)}
    atoms_by_normalized_name = {}
    for idx, atom in enumerate(atoms):
        atoms_by_normalized_name.setdefault(_normalize_atom_name(atom.name), []).append(idx)

    conf = mol.GetConformer()
    used_indices = set()

    for mol_idx in range(mol.GetNumAtoms()):
        rdkit_atom = mol.GetAtomWithIdx(mol_idx)
        pdb_info = rdkit_atom.GetPDBResidueInfo()
        pdb_name = pdb_info.GetName().strip() if pdb_info else ""
        normalized_pdb_name = _normalize_atom_name(pdb_name) if pdb_name else ""
        mol_element = rdkit_atom.GetSymbol()

        matched_idx = None

        # Priority 1: exact atom name match
        if pdb_name:
            candidate_idx = atoms_by_name.get(pdb_name)
            if candidate_idx is not None and candidate_idx not in used_indices:
                candidate_atom = atoms[candidate_idx]
                if candidate_atom.element == mol_element:
                    matched_idx = candidate_idx

        # Priority 2: normalized atom name match (strip A/B suffixes)
        if matched_idx is None and normalized_pdb_name:
            for candidate_idx in atoms_by_normalized_name.get(normalized_pdb_name, []):
                if candidate_idx in used_indices:
                    continue
                candidate_atom = atoms[candidate_idx]
                if candidate_atom.element == mol_element:
                    matched_idx = candidate_idx
                    break

        # Priority 3: nearest coordinate match with same element
        if matched_idx is None:
            pos = conf.GetAtomPosition(mol_idx)
            mol_coord = [pos.x, pos.y, pos.z]
            min_dist = float("inf")
            best_idx = None

            for atom_idx, atom in enumerate(atoms):
                if atom_idx in used_indices or atom.element != mol_element:
                    continue

                dist = sum((a - b) ** 2 for a, b in zip(mol_coord, atom_coords[atom_idx])) ** 0.5
                if dist < min_dist:
                    min_dist = dist
                    best_idx = atom_idx

            if best_idx is not None and min_dist < COORDINATE_MATCH_TOLERANCE_ANGSTROM:
                matched_idx = best_idx

        # Priority 4: nearest coordinate match regardless of element
        # Fallback for force fields (e.g., OPLS-AA LigParGen) that rename atoms
        # but preserve coordinates and order
        if matched_idx is None:
            pos = conf.GetAtomPosition(mol_idx)
            mol_coord = [pos.x, pos.y, pos.z]
            min_dist = float("inf")
            best_idx = None

            for atom_idx, atom in enumerate(atoms):
                if atom_idx in used_indices:
                    continue

                dist = sum((a - b) ** 2 for a, b in zip(mol_coord, atom_coords[atom_idx])) ** 0.5
                if dist < min_dist:
                    min_dist = dist
                    best_idx = atom_idx

            if best_idx is not None and min_dist < COORDINATE_MATCH_TOLERANCE_ANGSTROM:
                matched_idx = best_idx

        if matched_idx is not None:
            _assign_atom_metadata(rdkit_atom, atoms[matched_idx])
            used_indices.add(matched_idx)
            continue

        # Fallback: keep PDB name if present
        if pdb_name:
            rdkit_atom.SetProp("atomLabel", pdb_name)
            rdkit_atom.SetProp("name", pdb_name)
            if pdb_name in charge_dict:
                rdkit_atom.SetProp("atomNote", f"{charge_dict[pdb_name]:+.4f}")

    # Step 5: Generate 2D coordinates for depiction
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass  # Keep 3D coords if 2D generation fails

    return mol
