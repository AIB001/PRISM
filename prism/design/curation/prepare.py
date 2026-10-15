"""Prepare inspectable structural evidence; NEVER mark automatic drafts reviewed."""

import argparse
from collections import Counter, defaultdict
import csv
import io
import json
from pathlib import Path
import re

from rdkit import Chem

from ..chem import chemistry_copy
from ..io import file_ref, read_json, write_json
from .peptides import classify_peptide

PDB_ID = re.compile(r"^[0-9][a-z0-9]{3}$")
PATTERNS = {
    "carboxylic_acid_or_carboxylate": "[CX3](=O)[O;H1,-1]",
    "amide": "[CX3](=O)[NX3,NX4+]",
    "ester": "[CX3](=O)[OX2][#6]",
    "urea": "[NX3][CX3](=O)[NX3]",
    "carbamate": "[NX3][CX3](=O)[OX2][#6]",
    "sulfonamide": "[SX4](=O)(=O)[NX3]",
    "nitrile": "[#6]#[N]",
    "phosphate_or_phosphonate": "[P](=O)([O])[O]",
    "hydroxyl": "[O;H1;!$([O]-[C,P,S]=O)]",
    "thiol": "[S;H1]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "nitro": "[N+](=O)[O-]",
    "azo": "[N]=[N]",
    "guanidine": "[NX3][CX3](=[NX2,NX3+])[NX3]",
    "beta_lactam": "[C;R](=O)1[N;R][C;R][C;R]1",
    "epoxide": "[O;r3]1[C;r3][C;r3]1",
    "boronic_acid": "[B]([O;H1])[O;H1]",
    "aromatic_six_ring": "[a;r6]1[a;r6][a;r6][a;r6][a;r6][a;r6]1",
}
QUERIES = {key: Chem.MolFromSmarts(value) for key, value in PATTERNS.items()}
# Count functional centres, not SMARTS permutations (one phosphate yields
# three P/O/O matches). Preserve the union of atom IDs as the evidence set.
CENTRE = {"ketone": 1, "urea": 1, "carbamate": 1, "guanidine": 1}


def functional_groups(mol):
    groups = {}
    for name, query in QUERIES.items():
        centres = {}
        for match in mol.GetSubstructMatches(query):
            key = tuple(sorted(match)) if name in {"azo", "aromatic_six_ring"} else match[CENTRE.get(name, 0)]
            centres.setdefault(key, set()).update(i + 1 for i in match)
        if centres:
            groups[name] = [sorted(atoms) for atoms in centres.values()]
    return groups


CSV_FIELDS = [
    "PDB code",
    "Representative_PDB",
    "Protein Name",
    "Ligand Name",
    "Complex Type",
    "Formula",
    "Canonical SMILES",
    "Affinity Data",
    "pKd pKi pIC50",
    "UniProt AC",
    "Reference",
    "Pubmed ID",
    "Duplicate_PDBs",
]


def csv_evidence(directory, ids):
    profiles, by_id = [], defaultdict(list)
    for path in sorted(Path(directory).glob("*.csv")):
        raw = path.read_bytes()
        for encoding in ("utf-8-sig", "gb18030"):
            try:
                text = raw.decode(encoding)
                break
            except UnicodeDecodeError:
                continue
        else:
            profiles.append({"file": file_ref(path), "status": "UNREADABLE_ENCODING"})
            continue
        reader = csv.DictReader(io.StringIO(text))
        counts, keys = Counter(), Counter()
        for row_number, row in enumerate(reader, 2):
            counts["rows"] += 1
            key = (row.get("PDB code") or row.get("Representative_PDB") or "").strip().lower()
            if not PDB_ID.fullmatch(key):
                counts["invalid_or_blank_primary_keys"] += 1
                continue
            keys[key] += 1
            matches = [(key, "primary")]
            aliases = row.get("Duplicate_PDBs") or row.get("All_PDBs_for_SMILES") or ""
            for alias in re.split(r"[;,\s]+", aliases):
                alias = alias.strip().lower()
                if PDB_ID.fullmatch(alias) and alias != key:
                    matches.append((alias, "duplicate_alias_not_independent_metadata"))
            for ident, relation in matches:
                if ident not in ids:
                    continue
                counts["selected_case_links"] += 1
                by_id[ident].append(
                    {
                        "file": path.name,
                        "csv_record_number_1based_including_header": row_number,
                        "relation": relation,
                        "fields": {k: row[k] for k in CSV_FIELDS if row.get(k)},
                    }
                )
        profiles.append(
            {
                "file": file_ref(path),
                "status": "READ",
                "encoding": encoding,
                "columns": reader.fieldnames,
                "counts": dict(counts),
                "unique_valid_primary_keys": len(keys),
                "repeated_primary_keys": {k: n for k, n in keys.items() if n > 1},
            }
        )
    return profiles, by_id


def prepare(root, csv_directory):
    root = Path(root).resolve()
    selection = read_json(root / "selection_manifest.json")
    profiles, metadata = csv_evidence(csv_directory, set(selection["selected_ids"]))
    write_json(root / "csv_source_profile.json", profiles)
    rows, counts = [], Counter()
    for ident in selection["selected_ids"]:
        folder = root / "complexes" / ident
        evidence = read_json(folder / "evidence.json")
        mol = None
        if evidence.get("reference"):
            mol = Chem.SDMolSupplier(str(root / evidence["reference"]["mapped_sdf"]), sanitize=False, removeHs=False)[0]
            # classify_peptide must retain original source row IDs, not invent
            # heavy-only numbering. The mapped reference has that original order.
            mol = chemistry_copy(mol)
        peptide = classify_peptide(mol)
        csv_rows = metadata.get(ident, [])
        csv_peptide = [r for r in csv_rows if r["file"] == "pdbbind_2020_peptides.csv"]
        if csv_peptide and peptide["status"] == "nonpeptide":
            peptide["csv_disagreement"] = (
                "CSV peptide label lacks confirmed amino-acid connectivity; individual review required"
            )
        if not csv_peptide and peptide["status"] in {"peptide", "peptide_conjugate"}:
            peptide["csv_disagreement"] = "Structural peptide absent from historical peptide CSV"
        groups = functional_groups(mol) if mol is not None else {}
        descriptors = evidence.get("ligand", {}).get("descriptors", {})
        flags = list(evidence["quality_flags"])
        # Old flag used amide count as a peptide proxy. Preserve its lineage,
        # while never using that flag to establish peptide identity.
        flags = [f for f in flags if f != "LARGE_OR_PEPTIDE_LIKE_OUTSIDE_INITIAL_MVP"]
        if descriptors.get("heavy_atoms", 0) > 80:
            flags.append("LARGE_LIGAND_GT80_HEAVY_ATOMS")
        if groups.get("beta_lactam"):
            flags.append("BETA_LACTAM_CLASSIFICATION_REVIEW")
        if csv_peptide:
            flags.append("HISTORICAL_CSV_PEPTIDE_LABEL_CHECKED")
        primary = [r for r in csv_rows if r["file"] == "merged_pdb_data.csv" and r["relation"] == "primary"]
        names = sorted({r["fields"].get("Ligand Name", "") for r in primary})
        protein_names = sorted({r["fields"].get("Protein Name", "") for r in primary})
        caption = (
            f"结构图含 {descriptors.get('heavy_atoms', '未知')} 个重原子，"
            f"{descriptors.get('aromatic_rings', '未知')} 个芳香环，"
            f"形式电荷 {descriptors.get('formal_charge', '未知')}。"
        )
        record = {
            "schema": "prism-design-curation/1",
            "complex_id": ident,
            "selection_rank": evidence["selection_rank"],
            "evidence": file_ref(folder / "evidence.json"),
            "reference": evidence.get("reference"),
            "measurement_status": evidence["measurement_status"],
            "smiles": evidence.get("ligand", {}).get("canonical_isomeric_smiles"),
            "descriptors": descriptors,
            "functional_groups": groups,
            "peptide_analysis": peptide,
            "source_quality_flags": flags,
            "csv_ligand_names": names,
            "csv_protein_names": protein_names,
            "csv_evidence": csv_rows,
            "automatic_draft": {
                "structural_caption": caption,
                "warning": "Automatic draft only; not an individual semantic review",
            },
            "individual_review": {
                "status": "pending",
                "decision": None,
                "structure_reasoning_zh": None,
                "corrections_zh": None,
            },
            "training_eligible": False,
            "training_label_kind": None,
            "observed_optimization_pair": False,
            "affinity_improvement_label": None,
        }
        write_json(folder / "curation_record.json", record)
        (folder / "annotation.auto_draft.txt").write_text(
            f"{ident} | 自动草稿，尚未逐条语义复核\n{caption}\nSMILES: {record['smiles']}\n"
            f"肽分类（断键图检验）: {peptide['status']}；连通氨基酸单元: {peptide['largest_connected_amino_acid_group']}\n"
            f"官能团: {', '.join(groups)}\n质量标记: {', '.join(flags)}\n"
            "氢键、溶剂暴露、活性改善均不能由这一草稿直接确认。\n",
            encoding="utf-8",
        )
        cuts = evidence.get("cut_candidates", [])
        growth = evidence.get("growth_candidates", [])
        compact = {
            "rank": record["selection_rank"],
            "id": ident,
            "smiles": record["smiles"],
            "name": names,
            "target": protein_names,
            "peptide": peptide["status"],
            "aa": peptide["largest_connected_amino_acid_group"],
            "hydrolysis": [f["smiles"] for f in peptide["hydrolysis_fragments"] if f["alpha_amino_acid_backbones"]],
            "csv_peptide": bool(csv_peptide),
            "groups": {k: len(v) for k, v in groups.items()},
            "charge": descriptors.get("formal_charge"),
            "rotors": descriptors.get("rotatable_bonds_strict"),
            "flags": flags,
            "cuts": [
                {
                    k: c[k]
                    for k in ("candidate_id", "fragment_smiles", "editable_heavy_atoms", "mean_protein_neighbors_4A")
                }
                for c in cuts[:3]
            ],
            "growth": growth[:2],
            "nearest_protein_A": evidence.get("geometry", {}).get("nearest_protein_heavy_atom_A"),
        }
        rows.append(compact)
        counts[peptide["status"]] += 1
    with (root / "review_queue.jsonl").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")
    write_json(
        root / "automated_profile.json",
        {
            "records": len(rows),
            "peptide_status_counts": dict(counts),
            "individually_reviewed": 0,
            "note": "Profile is not a completed semantic annotation dataset",
        },
    )
    print(dict(counts), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--csv-directory", required=True)
    args = parser.parse_args()
    prepare(args.root, args.csv_directory)


if __name__ == "__main__":
    main()
