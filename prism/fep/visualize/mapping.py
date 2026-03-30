#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PNG visualization for FEP atom mapping.

Generates side-by-side 2D structure images with color-coded atom
classification highlighting. Replicates FEbuilder's visualization style.
"""

from pathlib import Path
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS

from prism.fep.core.mapping import AtomMapping
from .molecule import prepare_mol_with_charges_and_labels
from .highlight import create_highlight_info


def _align_mols_2d(mol_a: Chem.Mol, mol_b: Chem.Mol) -> bool:
    """
    Align 2D depictions to their maximum common structure.

    Replicates FEbuilder's align_mols_2d function.

    Returns
    -------
    bool
        True if alignment succeeded, False if it failed (but coordinates were generated)
    """
    try:
        mcs = Chem.rdFMCS.FindMCS(
            [mol_a, mol_b],
            atomCompare=rdFMCS.AtomCompare.CompareAny,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=True,
        )
        core = Chem.MolFromSmarts(mcs.smartsString)
        _ = AllChem.Compute2DCoords(core)

        for mol in [mol_a, mol_b]:
            _ = AllChem.Compute2DCoords(mol)
            _ = AllChem.GenerateDepictionMatching2DStructure(mol, core)
            _ = AllChem.NormalizeDepiction(mol)

        return True

    except ValueError as e:
        # 2D alignment failed - use default coordinates
        import warnings

        warnings.warn(f"Warning: 2D alignment failed, continuing without alignment: {e}", UserWarning)

        # Fallback: generate default 2D coordinates without alignment
        for mol in [mol_a, mol_b]:
            try:
                _ = AllChem.Compute2DCoords(mol)
            except Exception:
                # If even basic coordinate generation fails, skip this molecule
                pass

        return False


def visualize_mapping_png(
    mapping: AtomMapping,
    pdb_a: str,
    pdb_b: str,
    mol2_a: Optional[str],
    mol2_b: Optional[str],
    output_path: Optional[str] = None,
    edge: int = 600,
    legends: Tuple[str, str] = ("Reference", "Mutant"),
) -> None:
    """
    Generate PNG visualization of atom mapping.

    Replicates FEbuilder's see_common() function with three-color
    highlighting:
    - Green: Common atoms
    - Red: Transformed atoms
    - Blue: Surrounding atoms

    Parameters
    ----------
    mapping : AtomMapping
        Atom mapping result from DistanceAtomMapper
    pdb_a : str
        Path to reference ligand PDB file
    pdb_b : str
        Path to mutant ligand PDB file
    output_path : str, optional
        Output PNG file path (default: mapping_visualization.png)
    edge : int, optional
        Height of single ligand figure in pixels (default: 600)
    legends : Tuple[str, str], optional
        Figure legends for (reference, mutant) (default: ('Reference', 'Mutant'))

    Examples
    --------
    >>> from prism.fep.core.mapping import DistanceAtomMapper
    >>> from prism.fep.io import read_rtf_for_fep
    >>> from prism.fep.visualize import visualize_mapping_png
    >>> lig_a = read_rtf_for_fep('25.rtf', '25.pdb')
    >>> lig_b = read_rtf_for_fep('36.rtf', '36.pdb')
    >>> mapper = DistanceAtomMapper()
    >>> mapping = mapper.map(lig_a, lig_b)
    >>> visualize_mapping_png(mapping, '25.pdb', '36.pdb', '25_3D.mol2', '36_3D.mol2', '25-36_mapping.png')
    """
    # Collect all atoms from both ligands for charge labeling
    atoms_a = [a for a, _ in mapping.common] + mapping.transformed_a + mapping.surrounding_a
    atoms_b = [b for _, b in mapping.common] + mapping.transformed_b + mapping.surrounding_b

    # Generate RDKit Mol objects with bond order correction and charge labels (FEbuilder style)
    mol_a = prepare_mol_with_charges_and_labels(pdb_a, mol2_a, atoms_a)
    mol_b = prepare_mol_with_charges_and_labels(pdb_b, mol2_b, atoms_b)

    # Extract atom names by classification
    # IMPORTANT: Use ligand-specific names for each molecule
    common_names_a = [a.name for a, _ in mapping.common]
    common_names_b = [b.name for _, b in mapping.common]
    transformed_a_names = [a.name for a in mapping.transformed_a]
    surrounding_a_names = [a.name for a in mapping.surrounding_a]

    transformed_b_names = [a.name for a in mapping.transformed_b]
    surrounding_b_names = [a.name for a in mapping.surrounding_b]

    # Generate highlight information
    highlight_a = create_highlight_info(mol_a, common_names_a, transformed_a_names, surrounding_a_names)
    highlight_b = create_highlight_info(mol_b, common_names_b, transformed_b_names, surrounding_b_names)

    # Align 2D depictions (FEbuilder style)
    _align_mols_2d(highlight_a["mol"], highlight_b["mol"])

    # Configure drawing options (FEbuilder style)
    opts = Draw.MolDrawOptions()
    opts.centreMoleculesBeforeDrawing = True
    opts.baseFontSize = 0.45
    opts.annotationFontScale = 0.6
    opts.legendFraction = 0.3
    opts.legendFontSize = int(edge / 15)

    # Generate side-by-side image
    img = Draw.MolsToGridImage(
        [highlight_a["mol"], highlight_b["mol"]],
        subImgSize=(edge, edge),
        molsPerRow=2,
        legends=list(legends),
        highlightAtomLists=[highlight_a["hlist"], highlight_b["hlist"]],
        highlightAtomColors=[highlight_a["hmap"], highlight_b["hmap"]],
        drawOptions=opts,
    )

    # Save image
    if output_path is None:
        output_path = "mapping_visualization.png"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    img.save(str(output_path))
    print(f"  ✓ PNG saved to {output_path}")
