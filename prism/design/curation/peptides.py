"""Test peptide connectivity by virtual amide hydrolysis, not amide count.

The returned hydrolysate is a graph diagnostic, not a simulated chemical reaction.
Amide-sidechain hydrolysis can change Asn/Gln to acidic residues; its edges are
not counted as peptide links unless BOTH resulting fragments contain an amino
acid backbone. Non-alpha backbones and inseparable crosslinks are held for review.
"""

from collections import defaultdict

from rdkit import Chem

from ..chem import chemistry_copy

PROTOCOL = "amide-hydrolysis-amino-acid-connectivity/1"
# N protection (Boc/Cbz/sulfonyl) does not erase an amino-acid skeleton.
# The hydrolysate can retain these protecting groups. Connectivity between
# separate amino-acid fragments, rather than a free NH2, establishes a peptide.
ALPHA = Chem.MolFromSmarts("[N;!a]-[C;X4]-[C;X3](=[O])-[O;H1,-1]")
ACID = Chem.MolFromSmarts("[C;X3](=[O])-[O;H1,-1]")


def _amide_bonds(mol):
    found = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() == 7:
            a, b = b, a
        if a.GetAtomicNum() != 6 or b.GetAtomicNum() != 7 or a.GetIsAromatic() or b.GetIsAromatic():
            continue
        oxygen = [
            n
            for n in a.GetNeighbors()
            if n.GetAtomicNum() == 8
            and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
        ]
        if not oxygen:
            continue
        # Excludes urea/carbamate, whose carbonyl carbon is not attached to a C
        # skeleton and therefore is not an amino-carboxylic acid carboxyl group.
        if not any(n.GetAtomicNum() == 6 for n in a.GetNeighbors()):
            continue
        found.append((bond.GetIdx(), a.GetIdx(), b.GetIdx()))
    return found


def classify_peptide(molecule):
    result = {
        "protocol": PROTOCOL,
        "status": "unknown",
        "exclude_from_small_molecule_training": True,
        "amide_bonds_tested": [],
        "peptide_bonds": [],
        "hydrolysis_fragments": [],
        "amino_acid_fragment_count": 0,
        "largest_connected_amino_acid_group": 0,
        "original_heavy_atom_coverage_in_amino_acids": 0.0,
        "uncertainty_flags": [],
    }
    if molecule is None:
        result["uncertainty_flags"] = ["NO_PARSEABLE_GRAPH"]
        return result
    mol = chemistry_copy(molecule)
    for atom in mol.GetAtoms():
        atom.SetIntProp("_source_index_1based", atom.GetIdx() + 1)
    mol = Chem.RemoveHs(mol)
    amides = _amide_bonds(mol)
    source_id = {a.GetIdx(): a.GetIntProp("_source_index_1based") for a in mol.GetAtoms()}
    result["amide_bonds_tested"] = [
        {
            "carbonyl_atom_map": source_id[c],
            "nitrogen_atom_map": source_id[n],
            "in_ring": mol.GetBondWithIdx(b).IsInRing(),
        }
        for b, c, n in amides
    ]
    if not amides:
        result.update(status="nonpeptide", exclude_from_small_molecule_training=False)
        return result
    rw = Chem.RWMol(mol)
    for _, carbon, nitrogen in amides:
        rw.RemoveBond(carbon, nitrogen)
        oxygen = rw.AddAtom(Chem.Atom(8))
        rw.AddBond(carbon, oxygen, Chem.BondType.SINGLE)
        # Lost acyl valence is completed by implicit H; no original atom type,
        # bond order or formal charge elsewhere is 'repaired'.
        rw.GetAtomWithIdx(nitrogen).SetNoImplicit(False)
    product = rw.GetMol()
    try:
        Chem.SanitizeMol(product)
        atom_groups = []
        fragments = Chem.GetMolFrags(product, asMols=True, sanitizeFrags=True, fragsMolAtomMapping=atom_groups)
    except Exception as exc:
        result["uncertainty_flags"] = ["VIRTUAL_HYDROLYSIS_FAILED"]
        result["error"] = str(exc)
        return result
    fragment_for_atom = {idx: group for group, atoms in enumerate(atom_groups) for idx in atoms}
    aa_groups, nonalpha_groups = set(), set()
    original_covered = set()
    for group, (frag, indices) in enumerate(zip(fragments, atom_groups)):
        alpha = frag.GetSubstructMatches(ALPHA)
        sources = sorted(source_id[i] for i in indices if i in source_id)
        has_extended = False
        if not alpha:
            acids = [m[0] for m in frag.GetSubstructMatches(ACID)]
            nitrogens = [
                a.GetIdx()
                for a in frag.GetAtoms()
                if a.GetAtomicNum() == 7
                and not a.GetIsAromatic()
                and not any(
                    n.GetAtomicNum() == 6
                    and any(
                        b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(n).GetAtomicNum() == 8
                        for b in n.GetBonds()
                    )
                    for n in a.GetNeighbors()
                )
            ]
            for nitrogen in nitrogens:
                for carbon in acids:
                    path = Chem.GetShortestPath(frag, nitrogen, carbon)
                    if 4 <= len(path) <= 5 and all(
                        frag.GetAtomWithIdx(i).GetAtomicNum() == 6 and not frag.GetAtomWithIdx(i).GetIsAromatic()
                        for i in path[1:-1]
                    ):
                        has_extended = True
        if alpha:
            aa_groups.add(group)
            original_covered.update(sources)
        elif has_extended:
            nonalpha_groups.add(group)
        result["hydrolysis_fragments"].append(
            {
                "fragment_index": group,
                "smiles": Chem.MolToSmiles(frag),
                "source_atom_maps": sources,
                "alpha_amino_acid_backbones": len(alpha),
                "amino_acid_core_atom_maps": [[source_id.get(indices[i]) for i in match] for match in alpha],
                "nonalpha_amino_carboxylic_acid_candidate": has_extended,
            }
        )
    adjacency = defaultdict(set)
    for _, carbon, nitrogen in amides:
        a, b = fragment_for_atom[carbon], fragment_for_atom[nitrogen]
        if a in aa_groups and b in aa_groups and a != b:
            adjacency[a].add(b)
            adjacency[b].add(a)
            result["peptide_bonds"].append(
                {
                    "carbonyl_atom_map": source_id[carbon],
                    "nitrogen_atom_map": source_id[nitrogen],
                    "amino_acid_fragments": [a, b],
                }
            )
    largest = 0
    visited = set()
    for group in aa_groups:
        pending, component = [group], set()
        while pending:
            node = pending.pop()
            if node in component:
                continue
            component.add(node)
            pending.extend(adjacency[node] - component)
        visited.update(component)
        largest = max(largest, len(component))
    result.update(
        amino_acid_fragment_count=len(aa_groups),
        largest_connected_amino_acid_group=largest,
        original_heavy_atom_coverage_in_amino_acids=len(original_covered) / mol.GetNumHeavyAtoms(),
    )
    complex_residues = any(f["alpha_amino_acid_backbones"] > 1 for f in result["hydrolysis_fragments"])
    if largest >= 2:
        result["status"] = (
            "peptide" if result["original_heavy_atom_coverage_in_amino_acids"] >= 0.6 else "peptide_conjugate"
        )
    elif nonalpha_groups or complex_residues or (len(aa_groups) >= 2 and not result["peptide_bonds"]):
        result["status"] = "peptide_like_review"
        result["uncertainty_flags"].append("NONALPHA_OR_CROSSLINKED_BACKBONE_REQUIRES_REVIEW")
    else:
        result.update(status="nonpeptide", exclude_from_small_molecule_training=False)
    return result
