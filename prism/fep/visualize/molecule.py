#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Molecule processing utilities for FEP visualization.

Handles conversion from PDB files to RDKit Mol objects with proper
atom labeling, bond order correction, and charge annotation.
"""

from typing import List, Optional
from pathlib import Path


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
    pdb_mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
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


def auto_detect_mol2(pdb_file: str, search_dir: Optional[str] = None) -> Optional[str]:
    """
    自动检测与PDB文件对应的MOL2文件。

    搜索策略：
    1. 相同目录下，同名但扩展名为.mol2的文件
    2. 相同目录下，basename + "_3D.mol2" 的文件
    3. 处理特殊basename（如 "39_matched" -> "39"）
    4. 父目录中搜索同名文件
    5. 搜索指定目录中的匹配文件

    Parameters
    ----------
    pdb_file : str
        PDB文件路径
    search_dir : str, optional
        额外的搜索目录

    Returns
    -------
    str or None
        找到的MOL2文件路径，或None

    Examples
    --------
    >>> auto_detect_mol2("/path/to/39.pdb")
    "/path/to/39_3D.mol2"
    >>> auto_detect_mol2("/path/to/output/39_matched.pdb")
    "/path/to/39_3D.mol2"
    """
    from pathlib import Path

    pdb_path = Path(pdb_file)
    search_paths = []
    parent_dir = pdb_path.parent

    # 策略1: 同目录同名.mol2
    search_paths.append(pdb_path.with_suffix(".mol2"))

    # 策略2: 同目录 basename_3D.mol2
    search_paths.append(parent_dir / f"{pdb_path.stem}_3D.mol2")

    # 策略3: 处理特殊basename（如 "39_matched" -> "39"）
    # 移除常见后缀：_matched, _aligned, _processed 等
    basename = pdb_path.stem
    for suffix in ["_matched", "_aligned", "_processed", "_generated"]:
        if basename.endswith(suffix):
            base_stem = basename[: -len(suffix)]
            # 在当前目录和父目录中搜索
            search_paths.append(parent_dir / f"{base_stem}.mol2")
            search_paths.append(parent_dir / f"{base_stem}_3D.mol2")
            search_paths.append(parent_dir.parent / f"{base_stem}.mol2")
            search_paths.append(parent_dir.parent / f"{base_stem}_3D.mol2")
            break

    # 策略4: 父目录中搜索同名文件
    if parent_dir != parent_dir.parent:
        search_paths.append(parent_dir.parent / f"{pdb_path.stem}.mol2")
        search_paths.append(parent_dir.parent / f"{pdb_path.stem}_3D.mol2")

    # 策略5: 额外搜索目录
    if search_dir:
        search_dir = Path(search_dir)
        search_paths.append(search_dir / f"{pdb_path.stem}.mol2")
        search_paths.append(search_dir / f"{pdb_path.stem}_3D.mol2")

    # 返回第一个存在的文件
    for path in search_paths:
        if path.exists():
            return str(path)

    return None


def prepare_mol_with_charges_and_labels(
    pdb_file: str,
    mol2_file: Optional[str],
    atoms: List,  # List of Atom objects from prism.fep.core.mapping
) -> "Chem.Mol":
    """
    Prepare RDKit Mol object with correct bond orders, atom names and charge labels.

    Strategy:
    - For OpenFF/GAFF: Use MOL2/SDF as primary coordinate source (user input)
    - For CGenFF: Use PDB (from CHARMM-GUI) with MOL2 for bond orders if needed
    - Avoid SMILES round-trip to prevent atom reordering

    Parameters
    ----------
    pdb_file : str
        Path to PDB file (used for CGenFF or as fallback)
    mol2_file : str, optional
        Path to mol2 file with correct bond orders (PRIMARY for OpenFF/GAFF)
    atoms : List[Atom]
        List of Atom objects with name and charge properties

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
        detected_mol2 = auto_detect_mol2(pdb_file)
        if detected_mol2:
            mol2_file = detected_mol2

    # Step 1: Read molecule from best available source
    if mol2_file is not None and Path(mol2_file).exists():
        # For OpenFF/GAFF: read MOL2 directly (has coords + bonds)
        # This avoids atom reordering issues
        try:
            mol = Chem.MolFromMol2File(mol2_file, removeHs=False, sanitize=False)
            if mol is None:
                raise ValueError(f"Failed to read MOL2 file: {mol2_file}")

            # IMPORTANT: Match atoms BEFORE computing 2D coords
            # because atoms have 3D coordinates, and 2D coords will change positions
            charge_dict = {atom.name: atom.charge for atom in atoms}

            # Create a mapping from base names (without A/B suffix) to full names
            # This handles cases where hybrid topology atoms have A/B suffixes (e.g., H22A)
            # but MOL2 files have base names (e.g., H22)
            base_name_to_full = {}
            for atom in atoms:
                base_name = atom.name
                # Remove trailing A/B suffix (e.g., H22A -> H22, H02B -> H02)
                if base_name and base_name[-1] in "AB" and base_name[-2].isdigit():
                    base_name = base_name[:-1]
                base_name_to_full[base_name] = atom.name

            # Read MOL2 charges for all atoms (including those not in hybrid topology)
            # This is needed for atoms marked as dummy in hybrid topology but real in MOL2
            mol2_charges = {}
            try:
                with open(mol2_file, "r") as f:
                    in_atoms = False
                    for line in f:
                        if line.startswith("@<TRIPOS>ATOM"):
                            in_atoms = True
                            continue
                        if line.startswith("@"):
                            in_atoms = False
                            continue
                        if in_atoms:
                            parts = line.split()
                            if len(parts) >= 9:
                                atom_name = parts[1]
                                try:
                                    charge = float(parts[8])
                                    mol2_charges[atom_name] = charge
                                except (ValueError, IndexError):
                                    mol2_charges[atom_name] = 0.0
            except Exception as e:
                print(f"Warning: Could not read MOL2 charges: {e}")
                mol2_charges = {}

            # Build coordinate lookup from atoms parameter
            # atoms have coordinates in Angstroms from GRO (converted from nm)
            atom_coords = {i: atom.coord for i, atom in enumerate(atoms)}

            # Get RDKit mol conformer coordinates (3D, in Angstroms)
            conf = mol.GetConformer()
            mol_coords = []
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                mol_coords.append([pos.x, pos.y, pos.z])

            # Match by closest coordinates WITH ELEMENT CHECK
            # IMPORTANT: do not assume MOL2 atom order matches hybrid-topology order.
            # Hybrid topology may reorder atoms (common/transformed/dummy layout), and
            # using direct index matching can assign hydrogen names to carbon atoms and
            # vice versa, which corrupts both labels and displayed element identities.
            # CRITICAL: Must check element symbol match to prevent C/H mismatches!
            used_atom_indices = set()
            direct_match = False

            if not direct_match or len(used_atom_indices) == 0:
                # Original coordinate-based matching with element validation
                for mol_idx in range(mol.GetNumAtoms()):
                    mol_coord = mol_coords[mol_idx]
                    rdkit_atom = mol.GetAtomWithIdx(mol_idx)
                    mol_element = rdkit_atom.GetSymbol()

                    # Find closest atom in atoms list with matching element
                    min_dist = float("inf")
                    best_atom_idx = None

                    for atom_idx, atom_coord in atom_coords.items():
                        if atom_idx in used_atom_indices:
                            continue

                        # CRITICAL: Check element match first!
                        # This prevents carbon positions from being labeled with hydrogen names
                        target_atom = atoms[atom_idx]
                        if target_atom.element != mol_element:
                            continue

                        # Calculate Euclidean distance
                        dist = sum((a - b) ** 2 for a, b in zip(mol_coord, atom_coord)) ** 0.5

                        if dist < min_dist:
                            min_dist = dist
                            best_atom_idx = atom_idx

                    # If match found within reasonable tolerance (0.6 Å)
                    if best_atom_idx is not None and min_dist < 0.6:
                        name = atoms[best_atom_idx].name
                        rdkit_atom.SetProp("atomLabel", name)
                        rdkit_atom.SetProp("name", name)

                        # Add charge as note
                        if name in charge_dict:
                            charge = charge_dict[name]
                            rdkit_atom.SetProp("atomNote", f"{charge:+.4f}")

                        used_atom_indices.add(best_atom_idx)
                    else:
                        # Try to match using base name (without A/B suffix)
                        rdkit_atom = mol.GetAtomWithIdx(mol_idx)
                        rdkit_name = rdkit_atom.GetProp("name") if rdkit_atom.HasProp("name") else f"Atom{mol_idx}"

                        # Extract base name from RDKit atom name
                        base_name = rdkit_name
                        if base_name and base_name[-1] in "AB" and base_name[-2].isdigit():
                            base_name = base_name[:-1]

                        # Check if base name matches any atom in our list
                        if base_name in base_name_to_full:
                            full_name = base_name_to_full[base_name]
                            # Find the atom with this full name
                            for atom_idx, atom in enumerate(atoms):
                                if atom.name == full_name and atom_idx not in used_atom_indices:
                                    # Verify distance is reasonable
                                    atom_coord = atom_coords[atom_idx]
                                    dist = sum((a - b) ** 2 for a, b in zip(mol_coord, atom_coord)) ** 0.5

                                    if dist < 1.0:  # Use 1.0 Å tolerance for base name matching
                                        rdkit_atom.SetProp("atomLabel", full_name)
                                        rdkit_atom.SetProp("name", full_name)

                                        # Add charge as note
                                        if full_name in charge_dict:
                                            charge = charge_dict[full_name]
                                            rdkit_atom.SetProp("atomNote", f"{charge:+.4f}")

                                        used_atom_indices.add(atom_idx)
                                        break
                            else:
                                # Found base name match but no unused atom with that full name
                                # This shouldn't happen, but handle gracefully
                                print(
                                    f"Warning: Base name match found for {rdkit_name} -> {full_name}, but all atoms with that name are already used"
                                )
                        else:
                            # No match found even with base name
                            print(
                                f"Warning: Could not match RDKit atom {mol_idx} ({rdkit_name}) to any atom in list (min_dist={min_dist:.3f})"
                            )

            # For atoms that have no match in the hybrid topology:
            # Check if they should be preserved (real hydrogens connected to matched carbons)
            # Some atoms may be marked as dummy in hybrid topology but are real in MOL2
            unmatched = []
            to_preserve = set()

            for atom in mol.GetAtoms():
                if not atom.HasProp("name"):
                    # This atom didn't match any hybrid topology atom
                    # Check if it's a hydrogen connected to a matched atom
                    if atom.GetSymbol() == "H":
                        # Check neighbors
                        for neighbor in atom.GetNeighbors():
                            if neighbor.HasProp("name"):
                                # This hydrogen is connected to a matched atom
                                # Preserve it with MOL2 name (try _TriposAtomName first, then _MolFileAtomName)
                                if atom.HasProp("_TriposAtomName"):
                                    mol_name = atom.GetProp("_TriposAtomName")
                                elif atom.HasProp("_MolFileAtomName"):
                                    mol_name = atom.GetProp("_MolFileAtomName")
                                else:
                                    mol_name = f"H{atom.GetIdx()}"
                                atom.SetProp("name", mol_name)
                                atom.SetProp("atomLabel", mol_name)
                                # Estimate charge from similar atoms (0.03 for hc, 0.13 for ha)
                                # Use a default charge
                                # Use actual charge from MOL2 file if available
                                if mol_name in mol2_charges:
                                    charge = mol2_charges[mol_name]
                                else:
                                    charge = 0.03  # Default charge for hc hydrogens
                                atom.SetProp("atomNote", f"{charge:+.4f}")
                                to_preserve.add(atom.GetIdx())
                                break

                    if atom.GetIdx() not in to_preserve:
                        unmatched.append(atom.GetIdx())

            # Remove only atoms that should be removed (excluding preserved ones)
            # But only remove if we successfully matched some atoms
            if unmatched and len(unmatched) < mol.GetNumAtoms():
                from rdkit.Chem import RWMol

                rwmol = RWMol(mol)
                for idx in sorted(unmatched, reverse=True):
                    rwmol.RemoveAtom(idx)
                mol = rwmol.GetMol()
            elif unmatched:
                # All atoms are unmatched - distance matching completely failed
                # This means coordinates don't match between MOL2 and atoms list
                # Don't use the MOL2 molecule directly - instead raise exception to fall through
                # to PDB + MOL2 bond order approach
                raise ValueError("Distance matching failed - all atoms unmatched")

            # Try to sanitize for proper bond orders (after matching)
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                # If sanitization fails, try kekulization only
                try:
                    Chem.Kekulize(mol, clearAromaticFlags=True)
                except Exception:
                    pass  # Proceed with unsanitized mol

            # Generate 2D coordinates for depiction (after matching!)
            try:
                AllChem.Compute2DCoords(mol)
            except Exception:
                pass  # Keep 3D coords if 2D generation fails

            return mol

        except Exception as e:
            print(f"Warning: Failed to read MOL2, falling back to PDB: {e}")
            # Fall through to PDB-based approach

    # Fallback: Use PDB + MOL2 for bond orders (original CGenFF workflow)
    if mol2_file is not None:
        mol = assign_bond_orders_from_mol2(pdb_file, mol2_file)
    else:
        mol = pdb_to_mol(pdb_file)

    # IMPORTANT: For CHARMM-GUI files, skip the SMILES round-trip
    # The PDB file already has all atoms with correct names and coordinates
    # SMILES round-trip is only needed for OpenFF/GAFF which use SMILES as input
    # For CHARMM-GUI, use the PDB molecule directly to preserve atom names

    # Build charge lookup and coordinate list from Atom list
    charge_dict = {atom.name: atom.charge for atom in atoms}
    atom_coords = [(atom.coord[0], atom.coord[1], atom.coord[2]) for atom in atoms]

    # Match atoms by distance and assign correct names from atoms list
    # This handles cases where PDB and ITP use different atom naming conventions (e.g., OPLS)
    from rdkit import Chem
    from rdkit.Chem import AllChem

    conf = mol.GetConformer()
    used_atom_indices = set()

    for mol_idx in range(mol.GetNumAtoms()):
        rdkit_atom = mol.GetAtomWithIdx(mol_idx)
        mol_coord = conf.GetAtomPosition(mol_idx)
        mol_coord_tuple = (mol_coord.x, mol_coord.y, mol_coord.z)

        # Find closest atom from the atoms list
        min_dist = float("inf")
        best_atom_idx = None
        for atom_idx, atom_coord in enumerate(atom_coords):
            if atom_idx in used_atom_indices:
                continue
            dist = sum((a - b) ** 2 for a, b in zip(mol_coord_tuple, atom_coord)) ** 0.5
            if dist < min_dist:
                min_dist = dist
                best_atom_idx = atom_idx

        # Match if distance is reasonable (< 1.0 Å)
        if best_atom_idx is not None and min_dist < 1.0:
            atom = atoms[best_atom_idx]
            name = atom.name
            rdkit_atom.SetProp("atomLabel", name)
            rdkit_atom.SetProp("name", name)

            # Add charge as note
            if name in charge_dict:
                charge = charge_dict[name]
                rdkit_atom.SetProp("atomNote", f"{charge:+.4f}")

            used_atom_indices.add(best_atom_idx)
        else:
            # No match found - use PDB name as fallback
            pdb_info = rdkit_atom.GetPDBResidueInfo()
            if pdb_info:
                name = pdb_info.GetName().strip()
                rdkit_atom.SetProp("atomLabel", name)
                rdkit_atom.SetProp("name", name)

                # Try to look up charge with PDB name
                if name in charge_dict:
                    charge = charge_dict[name]
                    rdkit_atom.SetProp("atomNote", f"{charge:+.4f}")

    # Generate 2D coordinates
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass  # Keep 3D coords if 2D generation fails

    return mol
