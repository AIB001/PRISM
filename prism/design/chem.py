"""Chemical identity excludes atom-map identifiers; coordinates are never aligned."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors

from .errors import DesignError, require


def chemistry_copy(mol, keep_maps=False):
    result = Chem.Mol(mol)
    maps = [a.GetAtomMapNum() for a in result.GetAtoms()]
    for atom in result.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(result)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    if keep_maps:
        for atom, amap in zip(result.GetAtoms(), maps):
            atom.SetAtomMapNum(amap)
    return result


def load_reference(path):
    try:
        supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False, strictParsing=True)
        require(len(supplier) == 1, "Reference SDF must contain exactly one molecule", "REFERENCE_COUNT")
        mol = supplier[0]
        require(mol is not None, "Reference SDF cannot be parsed", "REFERENCE_PARSE_FAILED")
        mol = chemistry_copy(mol, keep_maps=True)
        require(mol.GetNumConformers() == 1, "Reference must have one conformer", "REFERENCE_POSE")
        require(
            np.isfinite(mol.GetConformer().GetPositions()).all(), "Non-finite reference coordinates", "REFERENCE_POSE"
        )
        require(len(Chem.GetMolFrags(mol)) == 1, "Multi-component reference is unsupported", "MULTICOMPONENT")
        maps = [a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        require(
            all(m > 0 for m in maps) and len(set(maps)) == len(maps),
            "Heavy atom maps must be unique positive integers",
            "ATOM_MAP_INVALID",
        )
        return mol
    except DesignError:
        raise
    except Exception as exc:
        raise DesignError(f"Invalid reference molecule: {exc}", "REFERENCE_PARSE_FAILED") from exc


def heavy_copy(mol, keep_maps=True):
    result = Chem.Mol(mol)
    for atom in result.GetAtoms():
        atom.SetIntProp("_prism_atom_map", atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
    result = chemistry_copy(result)
    result = Chem.RemoveHs(result)
    require(
        all(a.GetAtomicNum() > 1 for a in result.GetAtoms()),
        "Stereo-defining explicit H requires a dedicated backend mapping protocol",
        "EXPLICIT_H_UNSUPPORTED",
    )
    if keep_maps:
        for atom in result.GetAtoms():
            atom.SetAtomMapNum(atom.GetIntProp("_prism_atom_map"))
    return result


def canonical_smiles(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(chemistry_copy(mol)), isomericSmiles=True, canonical=True)


def descriptors(mol):
    m = Chem.RemoveHs(chemistry_copy(mol))
    return {
        "crippen_logp": float(Crippen.MolLogP(m)),
        "rotatable_bonds_strict": int(
            rdMolDescriptors.CalcNumRotatableBonds(m, rdMolDescriptors.NumRotatableBondsOptions.Strict)
        ),
        "formal_charge": int(Chem.GetFormalCharge(m)),
        "heavy_atoms": m.GetNumHeavyAtoms(),
    }


def map_to_index(mol):
    return {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}


def connected_subset(mol, indices):
    indices = set(indices)
    if not indices:
        return False
    seen = {next(iter(indices))}
    todo = list(seen)
    while todo:
        idx = todo.pop()
        for atom in mol.GetAtomWithIdx(idx).GetNeighbors():
            other = atom.GetIdx()
            if other in indices and other not in seen:
                seen.add(other)
                todo.append(other)
    return seen == indices
