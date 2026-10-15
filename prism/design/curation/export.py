"""Export reviewed labels, structure drawings and portable M1 records."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from ..chem import chemistry_copy, heavy_copy, load_reference
from ..io import read_json, sha256, verify_ref, write_json
from ..m1.data import SCHEMA, chemical_keys, prepare_manifest, validate_cut, validate_record
from .prepare import QUERIES, functional_groups
from .review import apply_reviews


def protein_groups(path, metadata):
    # Include exact sequence identities even when companion UniProt labels are
    # absent or inconsistent. These are NOT 30%-identity homology clusters.
    from Bio.SeqUtils import seq1

    chains = {}
    residues = set()
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM  "):
            continue
        key = (line[21:22], line[22:27])
        if key in residues:
            continue
        residues.add(key)
        chains.setdefault(key[0], []).append(seq1(line[17:20].strip(), custom_map={"MSE": "M"}))
    groups = [
        "sequence:" + hashlib.sha256("".join(seq).encode()).hexdigest() for seq in chains.values() if len(seq) >= 20
    ]
    for row in metadata:
        if row["relation"] != "primary":
            continue
        for acc in re.split(r"[,;\s]+", row["fields"].get("UniProt AC", "")):
            if re.fullmatch(r"[A-Z0-9]{6,10}(?:-[0-9]+)?", acc):
                groups.append("uniprot:" + acc.split("-")[0])
    if not groups:
        raise ValueError("No usable protein sequence/group evidence")
    return sorted(set(groups))


def export_reviewed(root):
    root = Path(root).resolve()
    # Recompute group-count evidence after centre-deduplication, independently
    # of the explicit semantic decisions. Nothing changes in original sources.
    corrections = []
    for ident in read_json(root / "selection_manifest.json")["selected_ids"]:
        path = root / "complexes" / ident / "curation_record.json"
        r = read_json(path)
        if r["reference"]:
            mol = chemistry_copy(load_reference(root / r["reference"]["mapped_sdf"]))
            groups = functional_groups(mol)
            old_counts = {k: len(mol.GetSubstructMatches(q)) for k, q in QUERIES.items() if mol.HasSubstructMatch(q)}
            new_counts = {k: len(v) for k, v in groups.items()}
            if old_counts != new_counts:
                corrections.append({"id": ident, "before": old_counts, "after": new_counts})
            r["functional_groups"] = groups
            r["functional_group_count_method"] = "unique chemical centre; atom-map evidence is union of matching atoms"
            write_json(path, r)
            # These are newly generated scientific structure drawings, not
            # claimed manual 3D inspections.
            drawer = rdMolDraw2D.MolDraw2DSVG(680, 480)
            drawmol = Chem.RemoveHs(mol)
            drawmol.RemoveAllConformers()
            rdMolDraw2D.PrepareAndDrawMolecule(drawer, drawmol, legend=ident)
            drawer.FinishDrawing()
            (path.parent / "structure_2d.svg").write_text(drawer.GetDrawingText(), encoding="utf-8")
    write_json(root / "functional_group_count_corrections.json", corrections)
    records = apply_reviews(root, require_complete=True)
    fields = ["rank", "id", "peptide_class", "training_eligible", "smiles", "note", "caption", "exclusions"]
    with (root / "annotations_reviewed.csv").open("w", encoding="utf-8-sig", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        for r in records:
            writer.writerow(
                dict(
                    rank=r["selection_rank"],
                    id=r["complex_id"],
                    peptide_class=r["peptide_class"],
                    training_eligible=r["training_eligible"],
                    smiles=r["smiles"],
                    note=r["individual_review"]["structure_reasoning_zh"],
                    caption=r["structural_caption_zh"],
                    exclusions=";".join(r["training_exclusion_reasons"]),
                )
            )
    directory = root / "m1_dataset"
    directory.mkdir(exist_ok=True)
    output, excluded, rejected_masks = [], [], []
    for r in records:
        ident = r["complex_id"]
        if not r["training_eligible"]:
            excluded.append({"id": ident, "reason": r["training_exclusion_reasons"]})
            continue
        e = read_json(verify_ref(root, r["evidence"]))
        try:
            mol = heavy_copy(load_reference(root / r["reference"]["mapped_sdf"]))
        except ValueError as exc:
            excluded.append({"id": ident, "reason": ["HEAVY_ATOM_PROTOCOL: " + str(exc)]})
            continue
        maps = [a.GetAtomMapNum() for a in mol.GetAtoms()]
        masks = []
        for cut in e["cut_candidates"]:
            mask = {
                "fixed": cut["fixed_atom_maps"],
                "variable": cut["editable_atom_maps"],
                "attachment_fixed": cut["retained_attachment_atom_map"],
                "attachment_variable": cut["replaced_attachment_atom_map"],
            }
            try:
                validate_cut(mol, maps, mask)
            except ValueError as exc:
                rejected_masks.append({"id": ident, "candidate": cut["candidate_id"], "reason": str(exc)})
                continue
            masks.append(mask)
        if not masks:
            excluded.append({"id": ident, "reason": ["no_legal_fragment_reconstruction_mask"]})
            continue
        dest = directory / "inputs" / ident
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / r["reference"]["mapped_sdf"], dest / "reference.sdf")
        for key in ("pocket", "protein"):
            source = verify_ref(root, e["files"][key])
            shutil.copyfile(source, dest / f"{key}.pdb")
        files = {
            key: {"path": str((dest / name).relative_to(directory)), "sha256": sha256(dest / name)}
            for key, name in (("reference", "reference.sdf"), ("pocket", "pocket.pdb"), ("protein", "protein.pdb"))
        }
        keys = chemical_keys(mol)
        keys["protein_groups"] = protein_groups(dest / "protein.pdb", r["csv_evidence"])
        record = {
            "schema": SCHEMA,
            "record_id": ident,
            "complex_group": ident,
            "task": "fragment_reconstruction",
            "text": r["structural_caption_zh"],
            "label_kind": "observed_structure_caption",
            "review": {
                "status": "completed",
                "decision": "nonpeptide",
                "quality_pass": True,
                "reviewer_kind": "assistant",
                "human_reviewed": False,
            },
            "files": files,
            "atom_maps": maps,
            "mask_candidates": masks,
            "split_keys": keys,
            "provenance": {
                "curation_record_sha256": sha256(root / "complexes" / ident / "curation_record.json"),
                "optimization_pair": False,
                "source": "local PDBbind 2020; user-provided license applies",
                "protein_grouping": "exact sequence and primary companion UniProt; homologs not fully controlled",
            },
        }
        try:
            validate_record(record, directory)
        except (ValueError, RuntimeError) as exc:
            excluded.append({"id": ident, "reason": ["M1_RECORD_VALIDATION: " + str(exc)]})
            continue
        output.append(record)
    manifest = directory / "records.jsonl"
    with manifest.open("w", encoding="utf-8") as handle:
        for r in output:
            handle.write(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n")
    write_json(directory / "excluded.json", excluded)
    write_json(directory / "rejected_masks.json", rejected_masks)
    split = prepare_manifest(manifest, directory / "split.json")
    report = {
        "curation_decisions": len(records),
        "caption_eligible": sum(r["training_eligible"] for r in records),
        "fragment_reconstruction_records": len(output),
        "excluded_from_reconstruction": len(excluded),
        "split_counts": split["counts"],
        "split_components": split["n_components"],
        "functional_group_count_corrected_cases": len(corrections),
        "manifest_sha256": sha256(manifest),
        "individual_corrections": len(read_json(root / "individual_corrections.json")),
    }
    write_json(root / "export_summary.json", report)
    print(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    export_reviewed(parser.parse_args().root)
